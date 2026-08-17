#include "poisson_FEM_MFEM.h"

using namespace mfem;

PoissonSolver::PoissonSolver(ParMesh &mesh, 
                             ParFiniteElementSpace &fespace, 
                             const std::vector<int> &inner_attrs,
                             const std::vector<int> &outer_attrs,
                             const std::vector<int> &artificial_attrs,
                             double phi_inner_val)
    : fespace(fespace), 
      inner_coeff(phi_inner_val), 
      outer_coeff(0.0),
      a(nullptr), m(nullptr), A(nullptr), 
      b(nullptr), 
      cg(nullptr), 
      smoother(nullptr)
{
    int max_attr = mesh.bdr_attributes.Max();
    inner_bdr_marker.SetSize(max_attr); inner_bdr_marker = 0;
    outer_bdr_marker.SetSize(max_attr); outer_bdr_marker = 0;
    ess_bdr.SetSize(max_attr); ess_bdr = 0;

    for (int a : inner_attrs) if (a >= 1 && a <= max_attr) inner_bdr_marker[a-1] = 1;
    for (int a : outer_attrs) if (a >= 1 && a <= max_attr) outer_bdr_marker[a-1] = 1;
    for (int a : artificial_attrs) if (a >= 1 && a <= max_attr) outer_bdr_marker[a-1] = 1;
    
    for (int i = 0; i < max_attr; ++i) 
        if (inner_bdr_marker[i] || outer_bdr_marker[i]) ess_bdr[i] = 1;

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    a = new ParBilinearForm(&fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator);
    a->Assemble();
    a->Finalize();
    
    A = new HypreParMatrix();
    a->FormSystemMatrix(ess_tdof_list, *A);

    // [核心优化] 质量矩阵 M 组装 (一次性操作)
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator());
    m->Assemble();
    m->Finalize();

    b = new ParLinearForm(&fespace);

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
    if (m) delete m; // 释放质量矩阵
    if (b) delete b;
}

void PoissonSolver::Solve(const Vector &rho, ParGridFunction &phi)
{
    phi = 0.0;
    phi.ProjectBdrCoefficient(inner_coeff, inner_bdr_marker);
    phi.ProjectBdrCoefficient(outer_coeff, outer_bdr_marker);

    // [核心优化] 极速组装 RHS：稀疏矩阵向量乘法 b = M * rho
    m->Mult(rho, *b);

    a->FormLinearSystem(ess_tdof_list, phi, *b, *A, X, B);

    X = 0.0; 
    cg->Mult(B, X);

    a->RecoverFEMSolution(X, *b, phi);
}