#include "CVPDE_TAKUMA.h"
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map> // [优化] 引入哈希表
#include <stack>
#include <iostream>
#include <memory>
#include <cstdint>       // [优化] 引入 uint64_t

using namespace mfem;

// =============================================================================
// 辅助函数
// =============================================================================
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

// =============================================================================
// 数据交换函数 (Ghost Cells)
// =============================================================================
void ExchangeUpwindData(ParMesh &mesh, const ParFiniteElementSpace &fespace,
                        const Vector& rho_local,
                        const std::vector<UpwindInfo> &Upwindelements,
                        std::unordered_map<uint64_t, UpwindData> &upwind_data, // [优化] 使用哈希表
                        const std::vector<mfem::DenseMatrix> &global_dshapes)  // [优化] 传入全局梯度
{
    int my_rank, num_procs;
    MPI_Comm comm = fespace.GetComm();
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    // 1. 识别请求
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

    // 2. 交换计数
    std::vector<int> send_counts(num_procs, 0), recv_counts(num_procs, 0);
    for(int r=0; r<num_procs; ++r) send_counts[r] = send_requests[r].size();
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, comm);

    // 3. 发送 elem_id
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

    // 4. 打包数据
    const int dim = mesh.Dimension();
    const int DOFS_PER_ELEM = 4;
    
    std::vector<std::vector<int>> send_dofs(num_procs);
    std::vector<std::vector<double>> send_dshape(num_procs);
    std::vector<std::vector<double>> send_rhos(num_procs); 

    std::vector<int> cnt_dofs(num_procs, 0), cnt_dshape(num_procs, 0), cnt_rhos(num_procs, 0);

    // [极致优化] 删除了局部的 Jacobian 求逆代码，直接读取 global_dshapes
    for (int r = 0; r < num_procs; ++r) {
        for (int elem_id : recv_elem_ids[r]) {
            Array<int> dofs;
            fespace.GetElementDofs(elem_id, dofs); 
            
            const DenseMatrix &dshape_e = global_dshapes[elem_id];

            for(int i=0; i<dofs.Size(); ++i) {
                send_dofs[r].push_back(fespace.GetGlobalTDofNumber(dofs[i]));
                double val = 0.0;
                if(dofs[i] >= 0 && dofs[i] < rho_local.Size()) val = rho_local(dofs[i]);
                send_rhos[r].push_back(val);
                for(int d=0; d<dim; ++d) send_dshape[r].push_back(dshape_e(i, d));
            }
        }
        cnt_dofs[r] = send_dofs[r].size();
        cnt_rhos[r] = send_rhos[r].size();
        cnt_dshape[r] = send_dshape[r].size();
    }

    // 5. 传输
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
    requests.clear();

    // 6. 解包
    for(int r=0; r<num_procs; ++r) {
        if(recv_dofs[r].empty()) continue;
        auto& ids = send_requests[r];
        int ptr_dof = 0;
        int ptr_shape = 0;
        
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
            // [极致优化] 填入 O(1) 哈希表
            upwind_data[PackRankElemKey(r, elem_id)] = data;
        }
    }
}

// =============================================================================
// 主求解函数：双栈波前推进 (Double Stack Wavefront)
// =============================================================================

void SolveContinuityTakuma_3D(
    mfem::ParMesh &mesh,
    mfem::ParFiniteElementSpace &fespace,
    mfem::ParFiniteElementSpace &vector_fes,
    const mfem::Vector &rho_prev,
    const mfem::ParGridFunction &E,
    const std::vector<UpwindInfo> &Upwindelements,
    mfem::Vector &rho_new,
    const std::vector<int> &corona_dofs,
    double K, mfem::Vector w_vec, double epsilon,
    const std::vector<mfem::DenseMatrix> &global_dshapes) // [优化] 确保签名传入了预计算阵列
{
    int numLocalDofs = fespace.GetNDofs();
    int my_rank;
    MPI_Comm_rank(mesh.GetComm(), &my_rank);
    double t_func_start = MPI_Wtime();

    // [极致优化] 删除了原有的 local_dshapes 预处理！直接使用传进来的 global_dshapes！

    // 2. 数据交换 (Ghost Cells)
    std::unordered_map<uint64_t, UpwindData> remote_data; // [优化] 使用 O(1) 哈希表
    ExchangeUpwindData(mesh, fespace, rho_prev, Upwindelements, remote_data, global_dshapes);

    // 3. 初始化状态
    std::vector<bool> is_fixed(numLocalDofs, false);
    for(int idx : corona_dofs) if(idx >=0 && idx < numLocalDofs) is_fixed[idx] = true;

    std::vector<bool> is_computed(numLocalDofs, false);
    for(int i=0; i<numLocalDofs; ++i) {
        if(is_fixed[i]) is_computed[i] = true;
    }

    rho_new = rho_prev; 

    // 4. 构建拓扑依赖图 (Dependency Graph)
    std::vector<std::vector<int>> forward_graph(numLocalDofs);
    std::vector<int> in_degree(numLocalDofs, 0);

    Array<int> graph_dofs_cache; // [优化] 提升循环内变量
    for (int i = 0; i < numLocalDofs; ++i) {
        if (is_fixed[i]) continue; 

        const auto& up_info = Upwindelements[i];
        if (up_info.elem_id < 0) continue; 

        int up_elem = up_info.elem_id;
        int up_rank = up_info.rank;

        if (up_rank == my_rank) {
            fespace.GetElementDofs(up_elem, graph_dofs_cache);
            for (int k = 0; k < graph_dofs_cache.Size(); ++k) {
                int dependency_node = graph_dofs_cache[k];
                if (dependency_node != i) {
                    forward_graph[dependency_node].push_back(i);
                    in_degree[i]++;
                }
            }
        }
    }

    // ---------------------------------------------------------
    // 5. 物理计算核心 (Lambda)
    // ---------------------------------------------------------
    
    // [极致优化] 将 Lambda 内部频繁分配的数组提取到外部！
    std::vector<double> neighbor_rhos; 
    neighbor_rhos.reserve(16);
    Array<int> lambda_dofs_cache; 

    auto ComputeNodePhysics = [&](int i) {
        const auto& up_info = Upwindelements[i];
        if(up_info.elem_id < 0) return; 

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
        int my_global_dof = fespace.GetGlobalTDofNumber(i);
        int self_index_in_element = -1;

        // [极致优化] 复用内存，只清空计数
        neighbor_rhos.clear(); 

        if(up_rank == my_rank) {
            // [极致优化] 使用全局静态缓存，告别 local_dshapes
            dshape_ptr = &global_dshapes[up_elem];
            
            fespace.GetElementDofs(up_elem, lambda_dofs_cache);
            for(int k=0; k<lambda_dofs_cache.Size(); ++k) {
                neighbor_rhos.push_back(rho_new(lambda_dofs_cache[k])); 
                if(lambda_dofs_cache[k] == i) self_index_in_element = k;
            }
        } else {
            // [极致优化] O(1) 哈希查询
            uint64_t key = PackRankElemKey(up_rank, up_elem);
            auto it = remote_data.find(key);
            if(it == remote_data.end()) return; 
            
            dshape_ptr = &it->second.dshape;
            neighbor_rhos = it->second.rho_vals;
            for(size_t k=0; k<it->second.dofs.size(); ++k) {
                if(it->second.dofs[k] == my_global_dof) {
                    self_index_in_element = k;
                    break;
                }
            }
        }

        if (self_index_in_element == -1) return;

        for(size_t k=0; k<neighbor_rhos.size(); ++k) {
            double v_dot_dN = Vx * (*dshape_ptr)(k, 0) + 
                              Vy * (*dshape_ptr)(k, 1) + 
                              Vz * (*dshape_ptr)(k, 2);
            if((int)k == self_index_in_element) {
                B_coef = v_dot_dN * epsilon;
            } else {
                C_coef += (v_dot_dN * neighbor_rhos[k]) * epsilon;
            }
        }
        rho_new(i) = SolveQuadraticRho(A_coef, B_coef, C_coef);
    };

    // ---------------------------------------------------------
    // 6. 救援检测 (Rescue Check)
    // ---------------------------------------------------------
    
    // [极致优化] 提领救援函数内部的动态数组
    Array<int> rescue_dofs_cache;
    
    auto CanRescue = [&](int i) -> bool {
        const auto& up_info = Upwindelements[i];
        if (up_info.elem_id < 0) return false;
        
        int up_elem = up_info.elem_id;
        int up_rank = up_info.rank;
        
        if (up_rank != my_rank) return true; // 远程边界默认可信

        fespace.GetElementDofs(up_elem, rescue_dofs_cache);
        double rho_threshold = 1e-9; 

        for (int k = 0; k < rescue_dofs_cache.Size(); ++k) {
            int neighbor = rescue_dofs_cache[k];
            if (neighbor == i) continue;
            
            // 只要有一个邻居是已计算且有值的，就可以救援
            if ((is_fixed[neighbor] || is_computed[neighbor]) && std::abs(rho_new(neighbor)) > rho_threshold) {
                return true;
            }
        }
        return false;
    };

    // ---------------------------------------------------------
    // 7. 算法主循环 (Wavefront Propagation)
    // ---------------------------------------------------------
    std::vector<int> stack_A; 
    std::vector<int> stack_B; 

    // 初始化栈：Corona 点 + 远程入流点
    for (int idx : corona_dofs) {
        if(idx >= 0 && idx < numLocalDofs) stack_A.push_back(idx);
    }
    
    for (int i = 0; i < numLocalDofs; ++i) {
        if (is_fixed[i]) continue;
        const auto& up_info = Upwindelements[i];
        // 远程驱动点视为源头
        if (up_info.rank != my_rank && up_info.rank >= 0) {
            stack_A.push_back(i);
        }
    }

    bool active = true;
    
    while (active) {
        active = false;

        // Stage I: 拓扑波前推进
        while (!stack_A.empty()) {
            while (!stack_A.empty()) {
                int u = stack_A.back();
                stack_A.pop_back();

                if (!is_fixed[u]) {
                    if (!is_computed[u]) {
                        ComputeNodePhysics(u);
                        is_computed[u] = true;
                    }
                }

                for (int v : forward_graph[u]) {
                    if (is_computed[v]) continue;
                    in_degree[v]--;
                    if (in_degree[v] == 0) {
                        stack_B.push_back(v);
                    }
                }
            }
            stack_A = stack_B;
            stack_B.clear();
        }

        // Stage II: 救援机制 (处理循环依赖和死区边缘)
        std::vector<int> rescued_nodes;
        for (int i = 0; i < numLocalDofs; ++i) {
            if (!is_fixed[i] && !is_computed[i]) {
                if (CanRescue(i)) {
                    ComputeNodePhysics(i); // 强制计算一次
                    is_computed[i] = true;
                    rescued_nodes.push_back(i);
                }
            }
        }

        // Stage III: 重启推进
        if (!rescued_nodes.empty()) {
            for (int node : rescued_nodes) {
                stack_A.push_back(node);
            }
            active = true;
        }
    }

    // 8. 最终清理 (死区置零)
    for (int i = 0; i < numLocalDofs; ++i) {
        if (!is_fixed[i] && !is_computed[i]) {
            rho_new(i) = 0.0;
            is_computed[i] = true;
        }
    }
    if (my_rank == 0) printf("   TOTAL TAKUMA Time  : %.6f s\n", MPI_Wtime() - t_func_start);
}