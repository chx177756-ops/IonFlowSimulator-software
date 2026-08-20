#ifndef CVPDE_TAKUMA_H
#define CVPDE_TAKUMA_H

#include "mfem.hpp"
#include <vector>
#include <unordered_map> // 替换 map
#include <set>
#include <stack>
#include <cstdint>       // 引入 uint64_t
#include "ParallelUpwindElementFinder.h" 

// [极致优化] 位运算打包 Key
inline uint64_t PackRankElemKey(int rank, int elem_id) {
    return (static_cast<uint64_t>(rank) << 32) | static_cast<uint32_t>(elem_id);
}

struct UpwindData {
    std::vector<int> dofs;         
    mfem::DenseMatrix dshape;      
    std::vector<double> rho_vals;  
};

void ExchangeUpwindData(mfem::ParMesh &mesh, 
                        const mfem::ParFiniteElementSpace &fespace,
                        const mfem::Vector& rho_local,
                        const std::vector<UpwindInfo> &Upwindelements,
                        std::unordered_map<uint64_t, UpwindData> &upwind_data, // 改为 unordered_map
                        const std::vector<mfem::DenseMatrix> &global_dshapes); 

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
    const std::vector<mfem::DenseMatrix> &global_dshapes);

#endif // CVPDE_TAKUMA_H