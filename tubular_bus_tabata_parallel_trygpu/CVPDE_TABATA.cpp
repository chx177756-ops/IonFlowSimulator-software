#include "CVPDE_TABATA.h"
#include "TabataAssemblyGPU.h"
#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

using namespace mfem;

// 辅助宏：用于累计时间
#define TICK(t) t -= MPI_Wtime()
#define TOCK(t) t += MPI_Wtime()

// 非本地进城上流单元的MPI通信
void ExchangeUpwindData(ParMesh &mesh, const ParFiniteElementSpace &fespace,
                        const std::vector<UpwindInfo> &Upwindelements,
                        std::map<std::pair<int, int>, UpwindData> &upwind_data,
                        const std::vector<mfem::DenseMatrix> &global_dshapes)
{
    int my_rank, num_procs;
    MPI_Comm comm = fespace.GetComm();
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    // 计时变量
    double t_identify = 0.0;
    double t_size_ex = 0.0;
    double t_req_ex = 0.0;
    double t_pack = 0.0;
    double t_comm_data = 0.0;
    double t_unpack = 0.0;
    double t_total_start = MPI_Wtime();

    // 步骤 1：识别非本地上流元
    TICK(t_identify);
    std::map<int, std::vector<int>> send_requests;
    for (size_t i = 0; i < Upwindelements.size(); ++i)
    {
        if (Upwindelements[i].rank != my_rank && Upwindelements[i].rank >= 0 && Upwindelements[i].rank < num_procs)
        {
            send_requests[Upwindelements[i].rank].push_back(Upwindelements[i].elem_id);
        }
    }
    TOCK(t_identify);

    // 步骤 2：交换请求大小
    TICK(t_size_ex);
    std::vector<int> send_counts(num_procs, 0);
    std::vector<int> recv_counts(num_procs, 0);
    for (int r = 0; r < num_procs; ++r)
    {
        send_counts[r] = send_requests[r].size();
    }
    MPI_Alltoall(&send_counts[0], 1, MPI_INT, &recv_counts[0], 1, MPI_INT, comm);
    TOCK(t_size_ex);

    // 步骤 3：发送请求的 elem_id
    TICK(t_req_ex);
    std::map<int, std::vector<int>> sorted_requests;
    for (auto &req : send_requests)
    {
        sorted_requests[req.first] = req.second;
        std::sort(sorted_requests[req.first].begin(), sorted_requests[req.first].end());
    }
    std::vector<MPI_Request> requests;
    std::vector<std::vector<int>> recv_elem_ids(num_procs);
    for (int r = 0; r < num_procs; ++r)
    {
        if (send_counts[r] > 0)
        {
            MPI_Request req;
            MPI_Isend(sorted_requests[r].data(), send_counts[r], MPI_INT, r, 0, comm, &req);
            requests.push_back(req);
        }
        if (recv_counts[r] > 0)
        {
            recv_elem_ids[r].resize(recv_counts[r]);
            MPI_Request req;
            MPI_Irecv(recv_elem_ids[r].data(), recv_counts[r], MPI_INT, r, 0, comm, &req);
            requests.push_back(req);
        }
    }
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
    TOCK(t_req_ex);

    // 步骤 4：处理接收到的请求，准备数据 (打包)
    TICK(t_pack);
    const int dim = mesh.Dimension();
    std::vector<std::vector<int>> send_dofs(num_procs);
    std::vector<std::vector<double>> send_dshape(num_procs);
    std::vector<int> send_dof_counts(num_procs, 0);
    std::vector<int> send_dshape_counts(num_procs, 0);
    
    // [优化] 变量提升，避免循环内分配

    for (int r = 0; r < num_procs; ++r)
    {
        if (recv_counts[r] == 0) continue;

        const int DOFS_PER_ELEM = 4;

        for (int elem_id : recv_elem_ids[r])
        {
            Array<int> dofs;
            fespace.GetElementDofs(elem_id, dofs);

            for (int i = 0; i < dofs.Size(); ++i)
            {
                int gdof = fespace.GetGlobalTDofNumber(dofs[i]);
                send_dofs[r].push_back(gdof);
            }

            // [优化] 直接读取预计算的形函数梯度
            const DenseMatrix &dshape_e = global_dshapes[elem_id];

            for (int i = 0; i < DOFS_PER_ELEM; ++i)
            {
                for (int d = 0; d < dim; ++d)
                {
                    send_dshape[r].push_back(dshape_e(i, d));
                }
            }
        }
        send_dof_counts[r] = send_dofs[r].size();
        send_dshape_counts[r] = send_dshape[r].size();
    }
    TOCK(t_pack);

    // 步骤 5：交换数据
    TICK(t_comm_data);
    std::vector<int> recv_dof_counts(num_procs);
    std::vector<int> recv_dshape_counts(num_procs);
    
    for (int r = 0; r < num_procs; ++r) {
        recv_dof_counts[r] = send_counts[r] * 4;
        recv_dshape_counts[r] = send_counts[r] * 4 * dim;
    }
    
    std::vector<std::vector<int>> recv_dofs(num_procs);
    std::vector<std::vector<double>> recv_dshape(num_procs);

    // 5A: DOFs
    for (int r = 0; r < num_procs; ++r)
    {
        if (send_dof_counts[r] > 0) {
            MPI_Request req;
            MPI_Isend(send_dofs[r].data(), send_dof_counts[r], MPI_INT, r, 10, comm, &req);
            requests.push_back(req);
        }
        if (recv_dof_counts[r] > 0) {
            recv_dofs[r].resize(recv_dof_counts[r]);
            MPI_Request req;
            MPI_Irecv(recv_dofs[r].data(), recv_dof_counts[r], MPI_INT, r, 10, comm, &req);
            requests.push_back(req);
        }
    }
    // 5B: dShapes
    for (int r = 0; r < num_procs; ++r)
    {
        if (send_dshape_counts[r] > 0) {
            MPI_Request req;
            MPI_Isend(send_dshape[r].data(), send_dshape_counts[r], MPI_DOUBLE, r, 11, comm, &req);
            requests.push_back(req);
        }
        if (recv_dshape_counts[r] > 0) {
            recv_dshape[r].resize(recv_dshape_counts[r]);
            MPI_Request req;
            MPI_Irecv(recv_dshape[r].data(), recv_dshape_counts[r], MPI_DOUBLE, r, 11, comm, &req);
            requests.push_back(req);
        }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
    TOCK(t_comm_data);

    // 步骤 6：解析接收到的数据
    TICK(t_unpack);
    for (int r = 0; r < num_procs; ++r)
    {
        if (recv_dof_counts[r] == 0) continue;

        auto &elem_ids = sorted_requests[r];
        const int DOFS_PER_ELEM = 4; 
        
        std::vector<std::vector<int>> dofs_per_elem(elem_ids.size(), std::vector<int>(DOFS_PER_ELEM));
        int dof_offset = 0;
        for (size_t e = 0; e < elem_ids.size(); ++e)
        {
            for (int i = 0; i < DOFS_PER_ELEM; ++i)
                dofs_per_elem[e][i] = recv_dofs[r][dof_offset++];
        }
        
        std::vector<DenseMatrix> dshape_per_elem(elem_ids.size(), DenseMatrix(DOFS_PER_ELEM, dim));
        int dshape_offset = 0;
        for (size_t e = 0; e < elem_ids.size(); ++e)
        {
            for (int i = 0; i < DOFS_PER_ELEM; ++i)
                for (int d = 0; d < dim; ++d)
                    dshape_per_elem[e](i, d) = recv_dshape[r][dshape_offset++];
        }
        
        for (size_t e = 0; e < elem_ids.size(); ++e)
        {
            UpwindData data;
            data.dofs = dofs_per_elem[e];
            data.dshape = dshape_per_elem[e];
            upwind_data[{r, elem_ids[e]}] = data;
        }
    }
    TOCK(t_unpack);

    // --- 输出通信计时信息 (汇总到 Rank 0) ---
    double local_times[6] = {t_identify, t_size_ex, t_req_ex, t_pack, t_comm_data, t_unpack};
    double max_times[6];
    MPI_Reduce(local_times, max_times, 6, MPI_DOUBLE, MPI_MAX, 0, comm);

    // if (my_rank == 0) {
    //     printf("  [Timer-Exchange] Identify: %.4f s\n", max_times[0]);
    //     printf("  [Timer-Exchange] SizeExch: %.4f s\n", max_times[1]);
    //     printf("  [Timer-Exchange] ReqExch : %.4f s\n", max_times[2]);
    //     printf("  [Timer-Exchange] PackData: %.4f s\n", max_times[3]);
    //     printf("  [Timer-Exchange] CommData: %.4f s\n", max_times[4]);
    //     printf("  [Timer-Exchange] Unpack  : %.4f s\n", max_times[5]);
    //     printf("  [Timer-Exchange] TOTAL   : %.4f s\n", MPI_Wtime() - t_total_start);
    // }
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
    
    // 总计时开始
    double t_func_start = MPI_Wtime();
    
    // 分段计时器
    double t_precalc_local = 0.0;
    double t_exchange_func = 0.0;
    double t_hypre_init = 0.0;
    double t_pre_noupwind = 0.0;
    double t_assembly_loop = 0.0;
    double t_hypre_assemble = 0.0;
    double t_rhs_assemble = 0.0;
    double t_solve_setup = 0.0;
    double t_solve_solve = 0.0;

    int numTet = mesh.GetNE();
    int dim = mesh.Dimension();

    // 1. 预处理本地各单元 (优化版：变量提升 + 避免循环内重复构造)
    TICK(t_precalc_local);
    // [优化] 形函数梯度已在 main 中预计算，通过 global_dshapes 传入，此处零开销
    TOCK(t_precalc_local);

    // 2. 预处理非本进程单元信息 (MPI)
    TICK(t_exchange_func);
    std::map<std::pair<int, int>, UpwindData> upwind_data;
    ExchangeUpwindData(mesh, fespace, Upwindelements, upwind_data, global_dshapes);
    TOCK(t_exchange_func);

    // 3. 创建 hypre 矩阵 (CPU 主机路径, 装配后切 GPU)
    TICK(t_hypre_init);

    int num_procs;
    MPI_Comm_size(comm, &num_procs);

    HYPRE_IJMatrix A;
    HYPRE_BigInt ilower, iupper;
    HYPRE_BigInt jlower, jupper;
    const HYPRE_BigInt *tdof_offsets = fespace.GetTrueDofOffsets();

    ilower = tdof_offsets[0];
    iupper = tdof_offsets[0 + 1] - 1;
    jlower = tdof_offsets[0];
    jupper = tdof_offsets[0 + 1] - 1;

    HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    // [Step 6 优化核心]: 预分配内存
    // 对于 P1 四面体网格，设置 50 是一个非常安全且高效的上限。
    int local_rows = iupper - ilower + 1;
    if (local_rows > 0) {
        std::vector<int> row_sizes(local_rows, 50); 
        HYPRE_IJMatrixSetRowSizes(A, row_sizes.data());
    }

    HYPRE_IJMatrixInitialize(A);
    TOCK(t_hypre_init);

    // 4. Pre-handle no-upwind nodes
    TICK(t_pre_noupwind);
    int local_ndofs = fespace.GetNDofs();
    std::vector<char> is_noupwind_local(local_ndofs, 0);
    int local_noupwind_count = 0;
    for (int ld = 0; ld < local_ndofs; ++ld)
    {
        if (Upwindelements[ld].elem_id < 0)
        {
            is_noupwind_local[ld] = 1;
            ++local_noupwind_count;
        }
    }
    const double BIG = 1e12;
    TOCK(t_pre_noupwind);

    Array<int> inner_dofs;
    inner_dofs.SetSize(corona_dofs.size());
    for (size_t i = 0; i < corona_dofs.size(); i++)
    {
        inner_dofs[i] = fespace.GetGlobalTDofNumber(corona_dofs[i]);
    }
    std::vector<bool> isBC(fespace.GlobalTrueVSize(), false);
    for (int td : inner_dofs)
    {
        if(td >= 0) isBC[td] = true;
    }

    // 5. 批量收集 + H2D + 一次 AddToValues (官方方式)
    // 5. GPU Kernel + batch AddToValues + corrective
    TICK(t_assembly_loop);

    // static cached elem_dofs+dshape_flat (Step1 dshape 已缓存)
    static std::vector<int> se_dofs; static std::vector<double> sd_flat;
    if(se_dofs.empty()){se_dofs.resize(numTet*4);sd_flat.resize(numTet*4*dim);
        for(int e=0;e<numTet;e++){Array<int>d;fespace.GetElementDofs(e,d);
            for(int i=0;i<4;i++){se_dofs[e*4+i]=d[i];
                for(int dd=0;dd<dim;dd++)sd_flat[(e*4+i)*dim+dd]=global_dshapes[e](i,dd);}}}

    std::vector<GIdx> h_ug(local_ndofs*4,-1);std::vector<double> h_ud(local_ndofs*4*dim,0.0);
    {Array<int> ul;for(int ld=0;ld<local_ndofs;ld++){int ue=Upwindelements[ld].elem_id,ur=Upwindelements[ld].rank;if(ue<0)continue;const DenseMatrix*p=nullptr;
        if(ur==my_rank){if(ue>=numTet)continue;fespace.GetElementDofs(ue,ul);if(ul.Size()!=4)continue;
            for(int k=0;k<4;k++)h_ug[ld*4+k]=fespace.GetGlobalTDofNumber(ul[k]);if(ue<(int)global_dshapes.size())p=&global_dshapes[ue];}
        else{auto it=upwind_data.find({ur,ue});if(it==upwind_data.end())continue;const auto&ud=it->second;if(ud.dofs.size()!=4)continue;
            for(int k=0;k<4;k++)h_ug[ld*4+k]=ud.dofs[k];p=&ud.dshape;}
        if(p)for(int j=0;j<4;j++)for(int dd=0;dd<dim;dd++)h_ud[(ld*4+j)*dim+dd]=(*p)(j,dd);}}

    std::vector<GIdx> h_gd(local_ndofs);std::vector<int> h_bc_v(local_ndofs);std::vector<double> h_vl(local_ndofs*dim);
    for(int ld=0;ld<local_ndofs;ld++){h_gd[ld]=fespace.GetGlobalTDofNumber(ld);HYPRE_BigInt gd=h_gd[ld];h_bc_v[ld]=(gd>=0&&gd<(HYPRE_BigInt)isBC.size())?(isBC[gd]?1:0):0;for(int dd=0;dd<dim;dd++)h_vl[ld*dim+dd]=nodeE_poisson(vector_fes.DofToVDof(ld,dd));}
    std::vector<double> h_rp(rho_previous.GetData(),rho_previous.GetData()+local_ndofs);
    std::vector<int> h_rptr(local_ndofs+1);std::vector<GIdx> h_ci;
    for(int ld=0;ld<local_ndofs;ld++){std::set<GIdx> cs;cs.insert(h_gd[ld]);int ub=ld*4;if(h_ug[ub]>=0)for(int j=0;j<4;j++){GIdx gj=h_ug[ub+j];if(gj>=0)cs.insert(gj);}for(auto c:cs)h_ci.push_back(c);h_rptr[ld+1]=(int)h_ci.size();}
    int nnz=(int)h_ci.size();

    double*d_vals;int nnz_out;
    TabataGPUComputeValues(numTet,dim,local_ndofs,my_rank,se_dofs,sd_flat,Volume,h_ug,h_ud,h_gd,h_bc_v,h_rp,h_vl,h_rptr,h_ci,&d_vals,&nnz_out);

    std::vector<HYPRE_Int> s5_nc(local_ndofs);std::vector<HYPRE_BigInt> s5_rw(local_ndofs);
    for(int ld=0;ld<local_ndofs;ld++){s5_nc[ld]=h_rptr[ld+1]-h_rptr[ld];s5_rw[ld]=h_gd[ld];}
    // Static GPU buffers for AddToValues (sizes stable after first call)
    static HYPRE_Int    *d_nc5=nullptr; static int nc5_cap=0;
    static HYPRE_BigInt *d_rw5=nullptr; static int rw5_cap=0;
    static HYPRE_BigInt *d_cl5=nullptr; static int cl5_cap=0;
    if(!d_nc5||nc5_cap<local_ndofs){if(d_nc5)cudaFree(d_nc5);cudaMalloc(&d_nc5,local_ndofs*sizeof(HYPRE_Int));nc5_cap=local_ndofs;}
    if(!d_rw5||rw5_cap<local_ndofs){if(d_rw5)cudaFree(d_rw5);cudaMalloc(&d_rw5,local_ndofs*sizeof(HYPRE_BigInt));rw5_cap=local_ndofs;}
    int cl5_bytes=nnz*sizeof(HYPRE_BigInt);
    if(!d_cl5||cl5_cap<cl5_bytes){if(d_cl5)cudaFree(d_cl5);cudaMalloc(&d_cl5,cl5_bytes);cl5_cap=cl5_bytes;}
    cudaMemcpy(d_nc5,s5_nc.data(),local_ndofs*sizeof(HYPRE_Int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_rw5,s5_rw.data(),local_ndofs*sizeof(HYPRE_BigInt),cudaMemcpyHostToDevice);
    cudaMemcpy(d_cl5,h_ci.data(),cl5_bytes,cudaMemcpyHostToDevice);
    HYPRE_IJMatrixAddToValues(A,local_ndofs,d_nc5,d_rw5,d_cl5,d_vals);
    cudaFree(d_vals);

    // Batch corrective — static GPU buffers
    {std::vector<HYPRE_BigInt> nw_r,bc_r;
        for(int ld=0;ld<local_ndofs;ld++){if(!is_noupwind_local[ld])continue;HYPRE_BigInt r=fespace.GetGlobalTDofNumber(ld);if(r>=ilower&&r<=iupper)nw_r.push_back(r);}
        for(int lc:corona_dofs){HYPRE_BigInt g=fespace.GetGlobalTDofNumber(lc);if(g>=ilower&&g<=iupper)bc_r.push_back(g);}
        static HYPRE_BigInt *d_nw_dr=nullptr,*d_bc_dr=nullptr;
        static HYPRE_Real   *d_nw_dv=nullptr,*d_bc_dv=nullptr;
        static HYPRE_Int    *d_nw_dnc=nullptr,*d_bc_dnc=nullptr;
        static int nw_cap=0,bc_cap=0;
        auto bs=[&](std::vector<HYPRE_BigInt>&rows,HYPRE_Real val,HYPRE_BigInt**d_dr,HYPRE_Real**d_dv,HYPRE_Int**d_dnc,int*cap){
            if(rows.empty())return;int n=(int)rows.size();
            if(!*d_dr||*cap<n){
                if(*d_dr)cudaFree(*d_dr);if(*d_dv)cudaFree(*d_dv);if(*d_dnc)cudaFree(*d_dnc);
                cudaMalloc(d_dr,n*sizeof(HYPRE_BigInt));cudaMalloc(d_dv,n*sizeof(HYPRE_Real));cudaMalloc(d_dnc,n*sizeof(HYPRE_Int));
                *cap=n;
            }
            cudaMemcpy(*d_dr,rows.data(),n*sizeof(HYPRE_BigInt),cudaMemcpyHostToDevice);
            std::vector<HYPRE_Real>vv(n,val);cudaMemcpy(*d_dv,vv.data(),n*sizeof(HYPRE_Real),cudaMemcpyHostToDevice);
            std::vector<HYPRE_Int>nc(n,1);cudaMemcpy(*d_dnc,nc.data(),n*sizeof(HYPRE_Int),cudaMemcpyHostToDevice);
            HYPRE_IJMatrixSetValues(A,n,*d_dnc,*d_dr,*d_dr,*d_dv);
        };
        bs(nw_r,BIG,&d_nw_dr,&d_nw_dv,&d_nw_dnc,&nw_cap);
        bs(bc_r,1.0,&d_bc_dr,&d_bc_dv,&d_bc_dnc,&bc_cap);
    }

    // 6. 矩阵组装 Finalize (CPU 路径, 最后统一切 GPU)
    TICK(t_hypre_assemble);
    
    // 执行最终组装
    TOCK(t_assembly_loop);
    // 由于在 Step 3 进行了预分配，这里不再需要昂贵的 realloc 操作
    HYPRE_IJMatrixAssemble(A);
    
    hypre_ParCSRMatrix *par_a;
    HYPRE_IJMatrixGetObject(A, (void **)&par_a);

    // #if 0 // --- DUMP DIAG (verify) ---
    // {
    //     hypre_CSRMatrix *diag = hypre_ParCSRMatrixDiag(par_a);
    //     ...
    //     fclose(fp);
    // }
    //
    // // --- DUMP FULL MATRIX (GPU) ---
    // {
    //     static int first=1;
    //     if(first){first=0; ... fclose(fp); }
    // }

    TOCK(t_hypre_assemble);

    // 7. RHS 组装
    TICK(t_rhs_assemble);
    HYPRE_IJVector B;
    HYPRE_IJVectorCreate(comm, ilower, iupper, &B);
    HYPRE_IJVectorSetObjectType(B, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(B);

    std::vector<HYPRE_BigInt> h_r7; std::vector<HYPRE_Real> h_v7;
    double SV = 1e-12 * BIG;
    for (int ld = 0; ld < local_ndofs; ++ld) {
        if (!is_noupwind_local[ld]) continue;
        HYPRE_BigInt row = fespace.GetGlobalTDofNumber(ld);
        if (row >= ilower && row <= iupper) { h_r7.push_back(row); h_v7.push_back(SV); }
    }
    for (int local_ldof : corona_dofs) {
        HYPRE_BigInt gi = fespace.GetGlobalTDofNumber(local_ldof);
        if (gi >= ilower && gi <= iupper) {
            h_r7.push_back(gi);
            h_v7.push_back((local_ldof>=0&&local_ldof<(int)rho_previous.Size())?rho_previous[local_ldof]:0.0);
        }
    }
    int NR=(int)h_r7.size();
    if(NR>0){
        static HYPRE_BigInt *d_r7=nullptr; static HYPRE_Real *d_v7=nullptr; static int r7_cap=0;
        if(!d_r7||r7_cap<NR){
            if(d_r7)cudaFree(d_r7); if(d_v7)cudaFree(d_v7);
            cudaMalloc(&d_r7,NR*sizeof(HYPRE_BigInt)); cudaMalloc(&d_v7,NR*sizeof(HYPRE_Real)); r7_cap=NR;
        }
        cudaMemcpy(d_r7,h_r7.data(),NR*sizeof(HYPRE_BigInt),cudaMemcpyHostToDevice);
        cudaMemcpy(d_v7,h_v7.data(),NR*sizeof(HYPRE_Real),cudaMemcpyHostToDevice);
        HYPRE_IJVectorSetValues(B,NR,d_r7,d_v7);
    }

    HYPRE_IJVectorAssemble(B);
    hypre_ParVector *par_b;
    HYPRE_IJVectorGetObject(B, (void **)&par_b);

    // #if 0 // --- DUMP RHS (GPU) ---
    // {
    //     static int first=1;
    //     if(first){first=0;
    //         ...
    //         fclose(fp);
    //     }
    // }
    // #endif

    TOCK(t_rhs_assemble);

    // 8. 创建解向量 X & Solver Setup (PCG + AMG Preconditioner)
    TICK(t_solve_setup);
    HYPRE_IJVector X;
    HYPRE_IJVectorCreate(comm, ilower, iupper, &X);
    HYPRE_IJVectorSetObjectType(X, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(X);
    // X 全行设零，static GPU buffer
    {
        static HYPRE_BigInt *d_xr=nullptr; static HYPRE_Real *d_xv=nullptr; static int xr_cap=0;
        static std::vector<HYPRE_BigInt> h_xr_static; static std::vector<HYPRE_Real> h_xv_static;
        HYPRE_BigInt lo=(HYPRE_BigInt)ilower, hi=(HYPRE_BigInt)iupper;
        int NX=(int)(hi-lo+1);
        if((int)h_xr_static.size()!=NX){
            h_xr_static.clear();h_xv_static.clear();
            for(HYPRE_BigInt g=lo;g<=hi;g++){h_xr_static.push_back(g);h_xv_static.push_back(0.0);}
            if(d_xr)cudaFree(d_xr);if(d_xv)cudaFree(d_xv);
            cudaMalloc(&d_xr,NX*sizeof(HYPRE_BigInt));cudaMalloc(&d_xv,NX*sizeof(HYPRE_Real));
            cudaMemcpy(d_xr,h_xr_static.data(),NX*sizeof(HYPRE_BigInt),cudaMemcpyHostToDevice);
            cudaMemcpy(d_xv,h_xv_static.data(),NX*sizeof(HYPRE_Real),cudaMemcpyHostToDevice);
            xr_cap=NX;
        }
        HYPRE_IJVectorSetValues(X,NX,d_xr,d_xv);
    }
    HYPRE_IJVectorAssemble(X);
    hypre_ParVector *par_x;
    HYPRE_IJVectorGetObject(X, (void **)&par_x);

    HYPRE_Solver solver, precond;

    // // --- A. 创建 PCG 求解器 ---
    // HYPRE_ParCSRPCGCreate(comm, &solver);
    // HYPRE_PCGSetMaxIter(solver, 2000);
    // HYPRE_PCGSetTol(solver, 1e-6);
    // HYPRE_PCGSetPrintLevel(solver, 0);

    // --- A. 创建 GMRES 求解器 ---
    HYPRE_ParCSRGMRESCreate(comm, &solver);
    HYPRE_ParCSRGMRESSetKDim(solver, 30);
    HYPRE_ParCSRGMRESSetMaxIter(solver, 2000);
    HYPRE_ParCSRGMRESSetTol(solver, 1e-6);
    HYPRE_ParCSRGMRESSetPrintLevel(solver, 0);

    // --- B. 创建 AMG 预条件子 ---
    HYPRE_BoomerAMGCreate(&precond);
    HYPRE_BoomerAMGSetPrintLevel(precond, 0);
    HYPRE_BoomerAMGSetMaxIter(precond, 1);
    HYPRE_BoomerAMGSetTol(precond, 0.0);

    HYPRE_BoomerAMGSetInterpType(precond, 6);
    HYPRE_BoomerAMGSetRelaxType(precond, 3);
    HYPRE_BoomerAMGSetCoarsenType(precond, 8);
    HYPRE_BoomerAMGSetCycleType(precond, 1);
    HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
    HYPRE_BoomerAMGSetRelaxOrder(precond, 0);
    HYPRE_BoomerAMGSetKeepTranspose(precond, 1);
    HYPRE_SetSpGemmUseVendor(0);
    HYPRE_SetUseGpuRand(1);

    // --- C. 将 AMG 挂载为 GMRES 的预条件子 ---
    HYPRE_ParCSRGMRESSetPrecond(solver,
                          (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSolve,
                          (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSetup,
                          precond);

    // --- D. Setup GMRES ---
    HYPRE_ParCSRGMRESSetup(solver, par_a, par_b, par_x);
    TOCK(t_solve_setup);

    // 9. Solver Solve
    TICK(t_solve_solve);
    HYPRE_ParCSRGMRESSolve(solver, par_a, par_b, par_x);
    TOCK(t_solve_solve);

    mfem::HypreParVector parx;
    parx.WrapHypreParVector(par_x, false);
    ParGridFunction rho_new_gf(&fespace);
    rho_new_gf.SetFromTrueDofs(parx);
    rho_new = rho_new_gf;

    // GMRES 迭代残差可能产生微小负值，裁剪到非负
    for (int i = 0; i < rho_new.Size(); i++)
        if (rho_new(i) < 0.0) rho_new(i) = 0.0;

    // 销毁顺序
    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(precond);
    
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(B);
    HYPRE_IJVectorDestroy(X);

    // --- 汇总并输出主函数时间 ---
    double local_times[9] = {
        t_precalc_local, t_exchange_func, t_hypre_init, t_pre_noupwind, 
        t_assembly_loop, t_hypre_assemble, t_rhs_assemble, t_solve_setup, t_solve_solve
    };
    double max_times[9];
    
    MPI_Reduce(local_times, max_times, 9, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (my_rank == 0) {
        printf("\n=== Tabata Solver (GMRES+AMG) Profile (MAX across ranks) ===\n");
        printf("1. Precalc Local Shapes : %.4f s\n", max_times[0]);
        printf("2. Exchange Upwind Info : %.4f s\n", max_times[1]);
        printf("3. Hypre Matrix Init    : %.4f s\n", max_times[2]);
        printf("4. Pre-handle NoUpwind  : %.4f s\n", max_times[3]);
        printf("5. Matrix Assembly Loop : %.4f s\n", max_times[4]);
        printf("6. Matrix Final Assemble: %.4f s\n", max_times[5]);
        printf("7. RHS Assembly         : %.4f s\n", max_times[6]);
        printf("8. GMRES+AMG Setup      : %.4f s\n", max_times[7]);
        printf("9. GMRES+AMG Solve      : %.4f s\n", max_times[8]);
        printf("----------------------------------------------------------\n");
        printf("   TOTAL Function Time  : %.4f s\n", MPI_Wtime() - t_func_start);
        printf("==========================================================\n\n");
    }
}