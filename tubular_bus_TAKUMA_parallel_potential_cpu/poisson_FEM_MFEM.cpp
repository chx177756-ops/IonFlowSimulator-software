#include "poisson_FEM_MFEM.h"

using namespace mfem;

double ParamCoefficient::Eval(ElementTransformation& T, const IntegrationPoint& ip) {
    if (!m_useFunc) return m_constant;
    Vector pt(3);
    T.Transform(ip, pt);
    std::map<std::string, double> vars = {{"x",pt[0]}, {"y",pt[1]}, {"z",pt[2]}};
    bool ok;
    double v = ParamParser::evaluate(m_expr, vars, "", &ok);
    return ok ? v : m_constant;
}

PoissonSolver::PoissonSolver(ParMesh &mesh,
                             ParFiniteElementSpace &fespace,
                             const std::map<int, double>& bdr_voltages,
                             const std::map<int, std::string>& bdr_expressions)
    : fespace(fespace),
      m_bdr_voltages(bdr_voltages),
      m_bdr_expressions(bdr_expressions),
      a(nullptr), m(nullptr), A(nullptr),
      b(nullptr), 
      cg(nullptr), 
      smoother(nullptr)
{
    // ==========================================
    // 1. 动态生成 Dirichlet 边界列表
    // ==========================================
    int max_attr = mesh.bdr_attributes.Max();
    ess_bdr.SetSize(max_attr); 
    ess_bdr = 0;

    for (auto const& [tag, voltage] : m_bdr_voltages) {
        if (tag >= 1 && tag <= max_attr) {
            ess_bdr[tag - 1] = 1;
        }
    }

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // 2. 刚度矩阵 A 组装与边界消除
    a = new ParBilinearForm(&fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator);
    a->Assemble();
    a->Finalize();
    
    A = new HypreParMatrix();
    a->FormSystemMatrix(ess_tdof_list, *A);

    // 3. 质量矩阵 M 组装
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator());
    m->Assemble();
    m->Finalize();

    b = new ParLinearForm(&fespace);

    // 4. 求解器配置
    cg = new CGSolver(fespace.GetComm());
    cg->SetRelTol(1e-6); 
    cg->SetMaxIter(2000);
    cg->SetPrintLevel(0); 
    
    smoother = new HypreSmoother(*A, 
                                 mfem::HypreSmoother::Type::l1GS,
                                 2, 1.0, 1.0, 2, 0.3);
    cg->SetPreconditioner(*smoother);
    cg->SetOperator(*A);
}

PoissonSolver::~PoissonSolver()
{
    if (cg) delete cg;
    if (smoother) delete smoother;
    if (A) delete A; 
    if (a) delete a;
    if (m) delete m; 
    if (b) delete b;
}

void PoissonSolver::Solve(const Vector &rho, ParGridFunction &phi)
{
    // ==========================================
    // 1. 动态遍历施加各个面的电压
    // ==========================================
    phi = 0.0;
    int max_attr = ess_bdr.Size();

    for (auto const& [tag, voltage] : m_bdr_voltages) {
        if (tag >= 1 && tag <= max_attr) {
            Array<int> marker(max_attr);
            marker = 0;
            marker[tag - 1] = 1;

            auto itExpr = m_bdr_expressions.find(tag);
            bool hasExpr = (itExpr != m_bdr_expressions.end());
            ParamCoefficient volt_coeff(voltage, hasExpr ? itExpr->second : "", hasExpr);
            phi.ProjectBdrCoefficient(volt_coeff, marker);
        }
    }

    // 2. 极速组装 RHS
    m->Mult(rho, *b);

    // 3. 生成线性系统 RHS
    a->FormLinearSystem(ess_tdof_list, phi, *b, *A, X, B);

    // 4. 求解
    X = 0.0; 
    cg->Mult(B, X);

    // 5. 恢复解
    a->RecoverFEMSolution(X, *b, phi);
}