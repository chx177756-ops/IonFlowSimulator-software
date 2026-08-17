#ifndef POISSON_FEM_MFEM_H
#define POISSON_FEM_MFEM_H

#include "mfem.hpp"
#include <vector>
#include <memory>

class PoissonSolver
{
public:
    PoissonSolver(mfem::ParMesh &mesh, 
                  mfem::ParFiniteElementSpace &fespace, 
                  const std::vector<int> &inner_attrs,
                  const std::vector<int> &outer_attrs,
                  const std::vector<int> &artificial_attrs,
                  double phi_inner_val);

    ~PoissonSolver();

    void Solve(const mfem::Vector &rho, mfem::ParGridFunction &phi);

private:
    mfem::ParFiniteElementSpace &fespace;
    
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;
    mfem::ConstantCoefficient inner_coeff;
    mfem::ConstantCoefficient outer_coeff;
    mfem::Array<int> inner_bdr_marker;
    mfem::Array<int> outer_bdr_marker;
    
    // 1. 矩阵相关
    mfem::ParBilinearForm *a;        
    mfem::ParBilinearForm *m;        // [核心优化] 质量矩阵形式
    mfem::HypreParMatrix *A;         
    
    // 2. RHS 相关
    mfem::ParLinearForm *b;          

    // 3. 求解器组件
    mfem::CGSolver *cg;
    mfem::HypreSmoother *smoother; 

    // 4. 临时向量
    mfem::Vector B, X;
};

#endif // POISSON_FEM_MFEM_H