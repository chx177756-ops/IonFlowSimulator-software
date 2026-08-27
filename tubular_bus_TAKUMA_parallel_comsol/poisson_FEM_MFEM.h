#ifndef POISSON_FEM_MFEM_H
#define POISSON_FEM_MFEM_H

#include "mfem.hpp"
#include <vector>
#include <memory>
#include <utility>

class PoissonSolver
{
public:
    // [极致优化] 使用键值对数组灵活接收多个不同电位的边界条件组，适应任意复杂的换流站模型
    PoissonSolver(mfem::ParMesh &mesh, 
                  mfem::ParFiniteElementSpace &fespace, 
                  const std::vector<std::pair<std::vector<int>, double>> &bdr_conditions);

    ~PoissonSolver();

    // 求解函数：轻量级，只组装 RHS 并求解
    void Solve(const mfem::Vector &rho, mfem::ParGridFunction &phi);

private:
    mfem::ParFiniteElementSpace &fespace;
    
    // 边界条件相关
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;
    std::vector<std::pair<std::vector<int>, double>> bdr_conditions_list;
    
    // --- 持久化对象 (避免重复分配) ---
    mfem::ParBilinearForm *a;        // 刚度矩阵形式
    mfem::ParBilinearForm *m;        // 质量矩阵形式
    mfem::HypreParMatrix *A;         // 组装并消除边界后的 Hypre 矩阵
    mfem::ParLinearForm *b;          // 线性形式 (直接接收 M*rho 的结果)

    mfem::CGSolver *cg;
    mfem::HypreSmoother *smoother; 

    // 临时向量 (避免循环内 malloc)
    mfem::Vector B, X;
};

#endif // POISSON_FEM_MFEM_H