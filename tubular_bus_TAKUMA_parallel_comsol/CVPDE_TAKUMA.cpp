#include "CVPDE_TAKUMA.h"
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace mfem;

// 辅助函数：解二次方程 (内联优化)
inline double SolveQuadraticRho(double A, double B, double C) {
    const double eps = 1e-14;
    if (std::abs(A) < eps) {
        if (std::abs(B) > eps) return std::max(0.0, -C / B);
        return 0.0;
    }
    double delta = B*B - 4.0*A*C;
    if (delta < 0) return std::max(0.0, -B / (2.0 * A));
    double sqrt_delta = std::sqrt(delta);
    double r1 = (-B + sqrt_delta) / (2.0 * A);
    double r2 = (-B - sqrt_delta) / (2.0 * A);
    return std::max(0.0, std::max(r1, r2));
}

void ExchangeUpwindData(ParMesh &mesh, const ParFiniteElementSpace &fespace,
                        const Vector& rho_local,
                        const std::vector<UpwindInfo> &Upwindelements,
                        std::unordered_map<uint64_t, UpwindData> &upwind_data,
                        const std::vector<mfem::DenseMatrix> &global_dshapes)
{
    int my_rank, num_procs;
    MPI_Comm comm = fespace.GetComm();
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    // 1. 识别请求 (略去中间不变的 MPI 收集逻辑，直接跳至末尾解析)
    std::map<int, std::vector<int>> send_requests;
    for (const auto& info : Upwindelements) {
        if (info.rank != my_rank && info.rank >= 0 && info.elem_id >= 0) {
            send_requests[info.rank].push_back(info.elem_id);
        }
    }
    for(auto& pair : send_requests) {
        std::sort(pair.second.begin(), pair.second.end());
        pair.second.erase(std::unique(pair.second.begin(), pair.second.end()), pair.second.end());
    }

    std::vector<int> send_counts(num_procs, 0), recv_counts(num_procs, 0);
    for(int r=0; r<num_procs; ++r) send_counts[r] = send_requests[r].size();
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

    std::vector<MPI_Request> requests;
    std::vector<std::vector<int>> recv_elem_ids(num_procs);
    for (int r = 0; r < num_procs; ++r) {
        if (send_counts[r] > 0) {
            MPI_Request req;
            MPI_Isend(send_requests[r].data(), send_counts[r], MPI_INT, r, 0, comm, &req);
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

    const int dim = mesh.Dimension();
    const int DOFS_PER_ELEM = 4; 
    
    std::vector<std::vector<int>> send_dofs(num_procs);
    std::vector<std::vector<double>> send_dshape(num_procs);
    std::vector<std::vector<double>> send_rhos(num_procs);
    std::vector<int> cnt_dofs(num_procs, 0), cnt_dshape(num_procs, 0), cnt_rhos(num_procs, 0);

    for (int r = 0; r < num_procs; ++r) {
        for (int elem_id : recv_elem_ids[r]) {
            Array<int> dofs;
            fespace.GetElementDofs(elem_id, dofs); 
            const DenseMatrix &dshape_e = global_dshapes[elem_id];
            for(int i=0; i<dofs.Size(); ++i) {
                send_dofs[r].push_back(fespace.GetGlobalTDofNumber(dofs[i]));
                double val = (dofs[i] >= 0 && dofs[i] < rho_local.Size()) ? rho_local(dofs[i]) : 0.0;
                send_rhos[r].push_back(val);
                for(int d=0; d<dim; ++d) send_dshape[r].push_back(dshape_e(i, d));
            }
        }
        cnt_dofs[r] = send_dofs[r].size();
        cnt_rhos[r] = send_rhos[r].size();
        cnt_dshape[r] = send_dshape[r].size();
    }

    std::vector<std::vector<int>> recv_dofs(num_procs);
    std::vector<std::vector<double>> recv_rhos(num_procs);
    std::vector<std::vector<double>> recv_dshape(num_procs);

    auto Communicate = [&](auto& send_buf, auto& cnt_send, auto& recv_buf, MPI_Datatype type, int tag_base, int multiplier) {
        std::vector<int> cnt_recv(num_procs);
        for(int r=0; r<num_procs; ++r) cnt_recv[r] = send_counts[r] * multiplier;
        
        for(int r=0; r<num_procs; ++r) {
            if(cnt_send[r] > 0) {
                MPI_Request req;
                MPI_Isend(send_buf[r].data(), cnt_send[r], type, r, tag_base, comm, &req);
                requests.push_back(req);
            }
            if(cnt_recv[r] > 0) {
                recv_buf[r].resize(cnt_recv[r]);
                MPI_Request req;
                MPI_Irecv(recv_buf[r].data(), cnt_recv[r], type, r, tag_base, comm, &req);
                requests.push_back(req);
            }
        }
    };

    Communicate(send_dofs, cnt_dofs, recv_dofs, MPI_INT, 100, DOFS_PER_ELEM);
    Communicate(send_rhos, cnt_rhos, recv_rhos, MPI_DOUBLE, 101, DOFS_PER_ELEM);
    Communicate(send_dshape, cnt_dshape, recv_dshape, MPI_DOUBLE, 102, DOFS_PER_ELEM * dim);

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    for(int r=0; r<num_procs; ++r) {
        if(recv_dofs[r].empty()) continue;
        auto& ids = send_requests[r];
        int ptr_dof = 0, ptr_shape = 0;
        for(size_t i=0; i<ids.size(); ++i) {
            int elem_id = ids[i];
            UpwindData data;
            data.dshape.SetSize(DOFS_PER_ELEM, dim);
            for(int k=0; k<DOFS_PER_ELEM; ++k) {
                data.dofs.push_back(recv_dofs[r][ptr_dof]);
                data.rho_vals.push_back(recv_rhos[r][ptr_dof]);
                ptr_dof++;
                for(int d=0; d<dim; ++d) {
                    data.dshape(k, d) = recv_dshape[r][ptr_shape++];
                }
            }
            // [极致优化] 填入 unordered_map，O(1) 索引
            upwind_data[PackRankElemKey(r, elem_id)] = data;
        }
    }
}

void SolveContinuityTakuma_3D(
    mfem::ParMesh &mesh,
    mfem::ParFiniteElementSpace &fespace,
    mfem::ParFiniteElementSpace &vector_fes,
    const mfem::Vector &rho_prev,
    const mfem::ParGridFunction &phi,
    const mfem::ParGridFunction &E,
    const std::vector<UpwindInfo> &Upwindelements,
    mfem::Vector &rho_new,
    const std::vector<int> &corona_dofs,
    double K, mfem::Vector w_vec, double epsilon,
    const std::vector<mfem::DenseMatrix> &global_dshapes)
{
    int numLocalDofs = fespace.GetNDofs();
    int my_rank;
    MPI_Comm_rank(mesh.GetComm(), &my_rank);

    // [优化] 使用 O(1) 的 hash map
    std::unordered_map<uint64_t, UpwindData> remote_data;
    ExchangeUpwindData(mesh, fespace, rho_prev, Upwindelements, remote_data, global_dshapes);

    std::vector<std::pair<double, int>> sorted_nodes(numLocalDofs);
    for(int i=0; i<numLocalDofs; ++i) {
        sorted_nodes[i] = {phi(i), i};
    }
    std::sort(sorted_nodes.begin(), sorted_nodes.end(), [](const auto& a, const auto& b){
        return a.first > b.first;
    });

    std::vector<char> is_fixed(numLocalDofs, 0);
    for(int idx : corona_dofs) if(idx >=0 && idx < numLocalDofs) is_fixed[idx] = 1;

    rho_new = rho_prev; 

    // [极致优化] 将动态数组移至热循环外部，彻底消除内存分配！
    std::vector<double> neighbor_rhos; 
    std::vector<int> neighbor_local_dofs; 
    neighbor_rhos.reserve(16);       // 预留足够空间避免扩张
    neighbor_local_dofs.reserve(16);

    Array<int> dofs_cache; // 缓存dofs

    for(const auto& item : sorted_nodes) {
        int i = item.second; 
        if(is_fixed[i]) continue; 

        const auto& up_info = Upwindelements[i];
        if(up_info.elem_id < 0) continue; 

        int up_elem = up_info.elem_id;
        int up_rank = up_info.rank;

        double Ex = E(vector_fes.DofToVDof(i, 0));
        double Ey = E(vector_fes.DofToVDof(i, 1));
        double Ez = E(vector_fes.DofToVDof(i, 2));
        double Vx = w_vec(0) + K * Ex;
        double Vy = w_vec(1) + K * Ey;
        double Vz = w_vec(2) + K * Ez;

        double A_coef = K * epsilon;
        double B_coef = 0.0;
        double C_coef = 0.0;

        const DenseMatrix* dshape_ptr = nullptr;
        
        // [极致优化] 重置 size 为 0，不释放内存容量
        neighbor_rhos.clear(); 
        neighbor_local_dofs.clear(); 

        if(up_rank == my_rank) {
            dshape_ptr = &global_dshapes[up_elem];
            fespace.GetElementDofs(up_elem, dofs_cache);
            for(int k=0; k<dofs_cache.Size(); ++k) {
                neighbor_rhos.push_back(rho_new(dofs_cache[k])); 
                neighbor_local_dofs.push_back(dofs_cache[k]);
            }
        } else {
            // [极致优化] O(1) 哈希极速匹配
            uint64_t key = PackRankElemKey(up_rank, up_elem);
            auto it = remote_data.find(key);
            if(it == remote_data.end()) continue; 
            
            dshape_ptr = &it->second.dshape;
            neighbor_rhos = it->second.rho_vals; // vector 浅拷贝或赋值，容量不变极快
        }

        int my_global_dof = fespace.GetGlobalTDofNumber(i); 
        int self_k = -1;

        for(size_t k=0; k<neighbor_rhos.size(); ++k) {
            bool is_self = false;
            if(up_rank == my_rank) {
                if(neighbor_local_dofs[k] == i) is_self = true;
            } else {
                uint64_t key = PackRankElemKey(up_rank, up_elem);
                if(remote_data[key].dofs[k] == my_global_dof) is_self = true;
            }

            double v_dot_dN = Vx * (*dshape_ptr)(k, 0) + Vy * (*dshape_ptr)(k, 1) + Vz * (*dshape_ptr)(k, 2);

            if(is_self) {
                self_k = k;
                B_coef = v_dot_dN * epsilon;
            } else {
                C_coef += (v_dot_dN * neighbor_rhos[k]) * epsilon;
            }
        }

        if (self_k == -1) continue; 
        rho_new(i) = SolveQuadraticRho(A_coef, B_coef, C_coef);
    }
}