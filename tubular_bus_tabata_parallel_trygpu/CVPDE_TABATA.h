#ifndef CVPDE_TABATA_H
#define CVPDE_TABATA_H

#include <map>
#include <vector>
#include "mfem.hpp"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include "ParallelUpwindElementFinder.h" // 包含 UpwindInfo 定义

// 定义存储非本地数据的结构
struct UpwindData {
    std::vector<int> dofs; // 全局自由度索引
    mfem::DenseMatrix dshape; // 导数形状函数 (3 x dim)
};

// 交换非本地上流元数据
void ExchangeUpwindData(
    mfem::ParMesh& mesh,
    const mfem::ParFiniteElementSpace& fespace,
    const std::vector<UpwindInfo>& Upwindelements,
    std::map<std::pair<int, int>, UpwindData>& upwind_data,
    const std::vector<mfem::DenseMatrix>& global_dshapes);

// 求解连续性方程（Tabata 上风方案）
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
    const std::vector<mfem::DenseMatrix>& global_dshapes);

#endif // CVPDE_TABATA_H