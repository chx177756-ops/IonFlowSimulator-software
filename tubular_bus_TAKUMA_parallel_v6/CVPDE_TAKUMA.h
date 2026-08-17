#ifndef CVPDE_TAKUMA_H
#define CVPDE_TAKUMA_H

#include "mfem.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>
#include "ParallelUpwindElementFinder.h"

// [极致优化] 采用位运算将 rank 和 elem_id 压缩为一个 64 位整数作为 Hash Key，
// 实现 O(1) 的无冲撞极速查找，替代原本极慢的 std::map<pair> 红黑树
inline uint64_t PackRankElemKey(int rank, int elem_id) {
    return (static_cast<uint64_t>(rank) << 32) | static_cast<uint32_t>(elem_id);
}

struct UpwindData {
    std::vector<int> dofs;         // 全局自由度索引
    mfem::DenseMatrix dshape;      // 导数形状函数
    std::vector<double> rho_vals;  // 邻居进程的 rho 值 (Ghost values)
};

void ExchangeUpwindData(mfem::ParMesh &mesh, 
                        const mfem::ParFiniteElementSpace &fespace,
                        const mfem::Vector& rho_local,
                        const std::vector<UpwindInfo> &Upwindelements,
                        std::unordered_map<uint64_t, UpwindData> &upwind_data,
                        const std::vector<mfem::DenseMatrix> &global_dshapes);

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
    const std::vector<mfem::DenseMatrix> &global_dshapes);

#endif // CVPDE_TAKUMA_H