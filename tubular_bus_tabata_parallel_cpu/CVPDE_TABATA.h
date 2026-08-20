#ifndef CVPDE_TABATA_H
#define CVPDE_TABATA_H

#include <unordered_map> // [优化] 使用 unordered_map
#include <vector>
#include <cstdint>       // [优化] 引入 uint64_t
#include "mfem.hpp"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include "ParallelUpwindElementFinder.h" 

// [极致优化] 采用位运算将 rank 和 elem_id 压缩为一个 64 位整数作为 Hash Key
inline uint64_t PackRankElemKey(int rank, int elem_id) {
    return (static_cast<uint64_t>(rank) << 32) | static_cast<uint32_t>(elem_id);
}

struct UpwindData {
    std::vector<int> dofs; 
    mfem::DenseMatrix dshape; 
};

void ExchangeUpwindData(
    mfem::ParMesh& mesh,
    const mfem::ParFiniteElementSpace& fespace,
    const std::vector<UpwindInfo>& Upwindelements,
    std::unordered_map<uint64_t, UpwindData>& upwind_data, // [优化] 替换为 unordered_map
    const std::vector<mfem::DenseMatrix> &global_dshapes); 

void SolveContinuityTabata(
    mfem::ParMesh &mesh,
    mfem::ParFiniteElementSpace &fespace,
    mfem::ParFiniteElementSpace &vector_fes,
    std::vector<double> Area,
    const mfem::Vector &rho_previous,
    const mfem::ParGridFunction &nodeE_poisson,
    const std::vector<UpwindInfo> &Upwindelements,
    mfem::Vector &rho_new,
    const std::vector<int> &corona_dofs,
    double K_mobility,
    const mfem::Vector &w_vec,
    const std::vector<mfem::DenseMatrix> &global_dshapes);

#endif // CVPDE_TABATA_H