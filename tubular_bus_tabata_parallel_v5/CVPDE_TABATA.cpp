#include "CVPDE_TABATA.h"
#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <unordered_map> // [优化] 引入哈希表
#include <algorithm>
#include <memory>
#include <cstdint>

using namespace mfem;

#define TICK(t) t -= MPI_Wtime()
#define TOCK(t) t += MPI_Wtime()

void ExchangeUpwindData(ParMesh &mesh, const ParFiniteElementSpace &fespace,
                        const std::vector<UpwindInfo> &Upwindelements,
                        std::unordered_map<uint64_t, UpwindData> &upwind_data, // [优化] 哈希表
                        const std::vector<mfem::DenseMatrix> &global_dshapes) 
{
    int my_rank, num_procs;
    MPI_Comm comm = fespace.GetComm();
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    double t_identify = 0.0, t_size_ex = 0.0, t_req_ex = 0.0, t_pack = 0.0, t_comm_data = 0.0, t_unpack = 0.0;
    
    TICK(t_identify);
    std::map<int, std::vector<int>> send_requests;
    for (size_t i = 0; i < Upwindelements.size(); ++i) {
        if (Upwindelements[i].rank != my_rank && Upwindelements[i].rank >= 0 && Upwindelements[i].rank < num_procs) {
            send_requests[Upwindelements[i].rank].push_back(Upwindelements[i].elem_id);
        }
    }
    TOCK(t_identify);

    TICK(t_size_ex);
    std::vector<int> send_counts(num_procs, 0), recv_counts(num_procs, 0);
    for (int r = 0; r < num_procs; ++r) send_counts[r] = send_requests[r].size();
    MPI_Alltoall(&send_counts[0], 1, MPI_INT, &recv_counts[0], 1, MPI_INT, comm);
    TOCK(t_size_ex);

    TICK(t_req_ex);
    std::map<int, std::vector<int>> sorted_requests;
    for (auto &req : send_requests) {
        sorted_requests[req.first] = req.second;
        std::sort(sorted_requests[req.first].begin(), sorted_requests[req.first].end());
    }
    std::vector<MPI_Request> requests;
    std::vector<std::vector<int>> recv_elem_ids(num_procs);
    for (int r = 0; r < num_procs; ++r) {
        if (send_counts[r] > 0) {
            MPI_Request req;
            MPI_Isend(sorted_requests[r].data(), send_counts[r], MPI_INT, r, 0, comm, &req);
            requests.push_back(req);
        }
        if (recv_counts[r] > 0) {
            recv_elem_ids[r].resize(recv_counts[r]);
            MPI_Request req;
            MPI_Irecv(recv_elem_ids[r].data(), recv_counts[r], MPI_INT, r, 0, comm, &req);
            requests.push_back(req);
        }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
    TOCK(t_req_ex);

    TICK(t_pack);
    const int dim = mesh.Dimension();
    std::vector<std::vector<int>> send_dofs(num_procs);
    std::vector<std::vector<double>> send_dshape(num_procs);
    std::vector<int> send_dof_counts(num_procs, 0), send_dshape_counts(num_procs, 0);
    
    for (int r = 0; r < num_procs; ++r) {
        if (recv_counts[r] == 0) continue;
        const int DOFS_PER_ELEM = 4;
        
        for (int elem_id : recv_elem_ids[r]) {
            Array<int> dofs;
            fespace.GetElementDofs(elem_id, dofs);
            for (int i = 0; i < dofs.Size(); ++i) {
                send_dofs[r].push_back(fespace.GetGlobalTDofNumber(dofs[i]));
            }
            
            const DenseMatrix &dshape_e = global_dshapes[elem_id];
            
            for (int i = 0; i < DOFS_PER_ELEM; ++i) {
                for (int d = 0; d < dim; ++d) {
                    send_dshape[r].push_back(dshape_e(i, d));
                }
            }
        }
        send_dof_counts[r] = send_dofs[r].size();
        send_dshape_counts[r] = send_dshape[r].size();
    }
    TOCK(t_pack);

    TICK(t_comm_data);
    std::vector<int> recv_dof_counts(num_procs), recv_dshape_counts(num_procs);
    for (int r = 0; r < num_procs; ++r) {
        recv_dof_counts[r] = send_counts[r] * 4;
        recv_dshape_counts[r] = send_counts[r] * 4 * dim;
    }
    
    std::vector<std::vector<int>> recv_dofs(num_procs);
    std::vector<std::vector<double>> recv_dshape(num_procs);

    for (int r = 0; r < num_procs; ++r) {
        if (send_dof_counts[r] > 0) {
            MPI_Request req; MPI_Isend(send_dofs[r].data(), send_dof_counts[r], MPI_INT, r, 10, comm, &req); requests.push_back(req);
        }
        if (recv_dof_counts[r] > 0) {
            recv_dofs[r].resize(recv_dof_counts[r]);
            MPI_Request req; MPI_Irecv(recv_dofs[r].data(), recv_dof_counts[r], MPI_INT, r, 10, comm, &req); requests.push_back(req);
        }
    }
    for (int r = 0; r < num_procs; ++r) {
        if (send_dshape_counts[r] > 0) {
            MPI_Request req; MPI_Isend(send_dshape[r].data(), send_dshape_counts[r], MPI_DOUBLE, r, 11, comm, &req); requests.push_back(req);
        }
        if (recv_dshape_counts[r] > 0) {
            recv_dshape[r].resize(recv_dshape_counts[r]);
            MPI_Request req; MPI_Irecv(recv_dshape[r].data(), recv_dshape_counts[r], MPI_DOUBLE, r, 11, comm, &req); requests.push_back(req);
        }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
    TOCK(t_comm_data);

    TICK(t_unpack);
    for (int r = 0; r < num_procs; ++r) {
        if (recv_dof_counts[r] == 0) continue;

        auto &elem_ids = sorted_requests[r];
        const int DOFS_PER_ELEM = 4; 
        
        std::vector<std::vector<int>> dofs_per_elem(elem_ids.size(), std::vector<int>(DOFS_PER_ELEM));
        int dof_offset = 0;
        for (size_t e = 0; e < elem_ids.size(); ++e) {
            for (int i = 0; i < DOFS_PER_ELEM; ++i)
                dofs_per_elem[e][i] = recv_dofs[r][dof_offset++];
        }
        
        std::vector<DenseMatrix> dshape_per_elem(elem_ids.size(), DenseMatrix(DOFS_PER_ELEM, dim));
        int dshape_offset = 0;
        for (size_t e = 0; e < elem_ids.size(); ++e) {
            for (int i = 0; i < DOFS_PER_ELEM; ++i)
                for (int d = 0; d < dim; ++d)
                    dshape_per_elem[e](i, d) = recv_dshape[r][dshape_offset++];
        }
        
        for (size_t e = 0; e < elem_ids.size(); ++e) {
            UpwindData data;
            data.dofs = dofs_per_elem[e];
            data.dshape = dshape_per_elem[e];
            // [极致优化] O(1) 哈希构建
            upwind_data[PackRankElemKey(r, elem_ids[e])] = data;
        }
    }
    TOCK(t_unpack);
}

void SolveContinuityTabata(
    ParMesh &mesh,
    mfem::ParFiniteElementSpace &fespace,
    mfem::ParFiniteElementSpace &vector_fes,
    std::vector<double> Volume,
    const Vector &rho_previous,
    const ParGridFunction &nodeE_poisson,
    const std::vector<UpwindInfo> &Upwindelements,
    Vector &rho_new,
    const std::vector<int> &corona_dofs,
    const std::vector<mfem::DenseMatrix> &global_dshapes) 
{
    int my_rank;
    MPI_Comm comm = fespace.GetComm();
    MPI_Comm_rank(comm, &my_rank);
    
    double t_func_start = MPI_Wtime();
    double t_precalc_local = 0.0, t_exchange_func = 0.0, t_hypre_init = 0.0, t_pre_noupwind = 0.0;
    double t_assembly_loop = 0.0, t_hypre_assemble = 0.0, t_rhs_assemble = 0.0, t_solve_setup = 0.0, t_solve_solve = 0.0;

    int numTet = mesh.GetNE();
    int dim = mesh.Dimension();

    TICK(t_precalc_local);
    TOCK(t_precalc_local);

    TICK(t_exchange_func);
    // [优化] O(1) 哈希表接收 Ghost
    std::unordered_map<uint64_t, UpwindData> upwind_data;
    ExchangeUpwindData(mesh, fespace, Upwindelements, upwind_data, global_dshapes); 
    TOCK(t_exchange_func);

    TICK(t_hypre_init);
    int num_procs;
    MPI_Comm_size(comm, &num_procs);

    HYPRE_IJMatrix A;
    HYPRE_BigInt ilower, iupper, jlower, jupper;
    const HYPRE_BigInt *tdof_offsets = fespace.GetTrueDofOffsets();
    ilower = tdof_offsets[0]; iupper = tdof_offsets[0 + 1] - 1;
    jlower = tdof_offsets[0]; jupper = tdof_offsets[0 + 1] - 1;

    HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    int local_rows = iupper - ilower + 1;
    if (local_rows > 0) {
        std::vector<int> row_sizes(local_rows, 50); 
        HYPRE_IJMatrixSetRowSizes(A, row_sizes.data());
    }
    HYPRE_IJMatrixInitialize(A);
    TOCK(t_hypre_init);

    TICK(t_pre_noupwind);
    int local_ndofs = fespace.GetNDofs();
    std::vector<char> is_noupwind_local(local_ndofs, 0);
    for (int ld = 0; ld < local_ndofs; ++ld) {
        if (Upwindelements[ld].elem_id < 0) is_noupwind_local[ld] = 1;
    }
    const double BIG = 1e12;
    for (int ld = 0; ld < local_ndofs; ++ld) {
        if (!is_noupwind_local[ld]) continue;
        HYPRE_BigInt row = fespace.GetGlobalTDofNumber(ld);
        if (row >= ilower && row <= iupper) {
            HYPRE_Int ncols = 1; HYPRE_BigInt row_arr[1] = {row}, col_arr[1] = {row}; HYPRE_Complex val_arr[1] = {BIG};
            HYPRE_IJMatrixSetValues(A, 1, &ncols, row_arr, col_arr, val_arr);
        }
    }
    TOCK(t_pre_noupwind);

    Array<int> inner_dofs(corona_dofs.size());
    for (size_t i = 0; i < corona_dofs.size(); i++) inner_dofs[i] = fespace.GetGlobalTDofNumber(corona_dofs[i]);
    std::vector<bool> isBC(fespace.GlobalTrueVSize(), false);
    for (int td : inner_dofs) if(td >= 0) isBC[td] = true;

    TICK(t_assembly_loop);
    std::vector<HYPRE_BigInt> batch_cols(4);
    std::vector<double> batch_vals(4);
    Array<int> upwind_local_dofs; 
    HYPRE_BigInt row_single[1]; 
    HYPRE_Int ncols_single = 1, ncols_batch = 4;
    HYPRE_Real val_single[1];

    for (int e = 0; e < numTet; e++) {
        Array<int> dofs;
        fespace.GetElementDofs(e, dofs);
        double volume = Volume[e];

        for (int i = 0; i < dofs.Size(); i++) {
            int ldof_i = dofs[i];
            HYPRE_BigInt gdof_i = fespace.GetGlobalTDofNumber(ldof_i);
            row_single[0] = gdof_i;

            if (isBC[gdof_i]) {
                val_single[0] = 1.0;
                HYPRE_IJMatrixSetValues(A, 1, &ncols_single, row_single, row_single, val_single);
                continue;
            }

            val_single[0] = rho_previous[ldof_i] * volume / 4.0;
            HYPRE_IJMatrixAddToValues(A, 1, &ncols_single, row_single, row_single, val_single);

            int up_elem = Upwindelements[ldof_i].elem_id;
            int up_rk = Upwindelements[ldof_i].rank;
            if (up_elem < 0) continue;

            const DenseMatrix* p_upwind_dshape = nullptr;
            
            if (up_rk == my_rank) {
                if (up_elem >= numTet) continue;
                fespace.GetElementDofs(up_elem, upwind_local_dofs);
                if (upwind_local_dofs.Size() != 4) continue;

                for(int k = 0; k < 4; ++k) batch_cols[k] = fespace.GetGlobalTDofNumber(upwind_local_dofs[k]);
                
                if (up_elem < (int)global_dshapes.size()) {
                     p_upwind_dshape = &global_dshapes[up_elem];
                }
            } else {
                // [极致优化] O(1) 哈希查询
                uint64_t key = PackRankElemKey(up_rk, up_elem);
                auto it = upwind_data.find(key);
                if (it == upwind_data.end()) continue; 
                
                const auto &ud = it->second;
                if (ud.dofs.size() != 4) continue;

                for(int k = 0; k < 4; ++k) batch_cols[k] = ud.dofs[k];
                p_upwind_dshape = &ud.dshape;
            }

            if (!p_upwind_dshape) continue;

            for (int j = 0; j < 4; j++) {
                double conv = 0.0;
                for (int d = 0; d < dim; d++) {
                    int vldof = vector_fes.DofToVDof(ldof_i, d);
                    double velo_d = nodeE_poisson(vldof);
                    conv += velo_d * (*p_upwind_dshape)(j, d);
                }
                batch_vals[j] = conv * volume / 4.0;
            }
            HYPRE_IJMatrixAddToValues(A, 1, &ncols_batch, row_single, batch_cols.data(), batch_vals.data());
        }
    }
    TOCK(t_assembly_loop);

    TICK(t_hypre_assemble);
    HYPRE_IJMatrixAssemble(A);
    hypre_ParCSRMatrix *par_a;
    HYPRE_IJMatrixGetObject(A, (void **)&par_a);
    TOCK(t_hypre_assemble);

    TICK(t_rhs_assemble);
    HYPRE_IJVector B;
    HYPRE_IJVectorCreate(comm, ilower, iupper, &B);
    HYPRE_IJVectorSetObjectType(B, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(B);

    double SMALL_VAL = 1e-12 * BIG;
    for (int ld = 0; ld < local_ndofs; ++ld) {
        if (!is_noupwind_local[ld]) continue;
        HYPRE_BigInt row = fespace.GetGlobalTDofNumber(ld);
        if (row >= ilower && row <= iupper) HYPRE_IJVectorSetValues(B, 1, &row, &SMALL_VAL);
    }

    for (int local_ldof : corona_dofs) {
        HYPRE_BigInt gi = fespace.GetGlobalTDofNumber(local_ldof);
        if (gi >= ilower && gi <= iupper) {
            HYPRE_Real gval = 0.0;
            if (local_ldof >= 0 && local_ldof < rho_previous.Size()) gval = rho_previous[local_ldof];
            HYPRE_IJVectorSetValues(B, 1, &gi, &gval);
        }
    }

    HYPRE_IJVectorAssemble(B);
    hypre_ParVector *par_b;
    HYPRE_IJVectorGetObject(B, (void **)&par_b);
    TOCK(t_rhs_assemble);

    TICK(t_solve_setup);
    HYPRE_IJVector X;
    HYPRE_IJVectorCreate(comm, ilower, iupper, &X);
    HYPRE_IJVectorSetObjectType(X, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(X);
    HYPRE_IJVectorAssemble(X);
    hypre_ParVector *par_x;
    HYPRE_IJVectorGetObject(X, (void **)&par_x);

    HYPRE_Solver solver, precond;
    HYPRE_ParCSRGMRESCreate(comm, &solver);
    HYPRE_GMRESSetKDim(solver, 30);        
    HYPRE_GMRESSetMaxIter(solver, 2000);   
    HYPRE_GMRESSetTol(solver, 1e-6);      
    HYPRE_GMRESSetPrintLevel(solver, 0);   

    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1); 
    HYPRE_BoomerAMGSetTol(precond, 0.0);   
    HYPRE_BoomerAMGSetInterpType(precond, 6);        
    HYPRE_BoomerAMGSetRelaxType(precond, 3);         
    HYPRE_BoomerAMGSetCoarsenType(precond, 8);       
    HYPRE_BoomerAMGSetCycleType(precond, 1);         
    HYPRE_BoomerAMGSetStrongThreshold(precond, 0.85);
    HYPRE_BoomerAMGSetRestriction(precond, 1);

    HYPRE_GMRESSetPrecond(solver,
                          (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup,
                          precond);

    HYPRE_ParCSRGMRESSetup(solver, par_a, par_b, par_x);
    TOCK(t_solve_setup);

    TICK(t_solve_solve);
    HYPRE_ParCSRGMRESSolve(solver, par_a, par_b, par_x);
    TOCK(t_solve_solve);

    mfem::HypreParVector parx;
    parx.WrapHypreParVector(par_x, false);
    ParGridFunction rho_new_gf(&fespace);
    rho_new_gf.SetFromTrueDofs(parx);
    rho_new = rho_new_gf;

    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(B);
    HYPRE_IJVectorDestroy(X);
}