#ifndef POISSON_FEM_MFEM_H
#define POISSON_FEM_MFEM_H

#include "mfem.hpp"
#include "param_parser.h"
#include <vector>
#include <map>
#include <memory>
#include <string>

// 支持常数和函数表达式的边界条件系数 (与 v6 CPU 版本接口一致)
class ParamCoefficient : public mfem::Coefficient {
    double m_constant;
    std::string m_expr;
    bool m_useFunc;
public:
    ParamCoefficient(double c, const std::string& e = "", bool f = false)
        : mfem::Coefficient(), m_constant(c), m_expr(e), m_useFunc(f) {}
    double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override;
};

class PoissonSolver
{
public:
    PoissonSolver(mfem::ParMesh &mesh,
                  mfem::ParFiniteElementSpace &fespace,
                  const std::map<int, double>& bdr_voltages,
                  const std::map<int, std::string>& bdr_expressions = {});

    ~PoissonSolver();

    void Solve(const mfem::Vector &rho, mfem::ParGridFunction &phi);

private:
    mfem::ParFiniteElementSpace &fespace;

    std::map<int, double> m_bdr_voltages;
    std::map<int, std::string> m_bdr_expressions;

    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    mfem::ParBilinearForm *a;        // 刚度矩阵
    mfem::ParBilinearForm *m;        // 质量矩阵 (RHS 组装)
    mfem::HypreParMatrix *A;         // 消除边界后的 Hypre 矩阵

    mfem::ParLinearForm *b;

    mfem::CGSolver *cg;
    mfem::HypreSmoother *smoother;

    mfem::Vector B, X;
};

#endif // POISSON_FEM_MFEM_H
