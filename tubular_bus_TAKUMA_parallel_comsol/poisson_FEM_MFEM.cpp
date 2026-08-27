#include "poisson_FEM_MFEM.h"

using namespace mfem;

PoissonSolver::PoissonSolver(ParMesh &mesh, 
                             ParFiniteElementSpace &fespace, 
                             const std::vector<std::pair<std::vector<int>, double>> &bdr_conditions)
    : fespace(fespace), 
      a(nullptr), m(nullptr), A(nullptr), 
      b(nullptr), cg(nullptr), smoother(nullptr),
      bdr_conditions_list(bdr_conditions)
{
    // ==========================================
    // 1. 边界条件准备与组合提取
    // ==========================================
    int max_attr = mesh.bdr_attributes.Max();
    ess_bdr.SetSize(max_attr); ess_bdr = 0;

    // 提取所有需要被设定为 Essential (Dirichlet) 的边界
    for (const auto& pair : bdr_conditions) {
        for (int attr : pair.first) {
            if (attr >= 1 && attr <= max_attr) {
                ess_bdr[attr-1] = 1;
            }
        }
    }

    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // ==========================================
    // 2. 刚度矩阵 A 组装与边界消除 (一次性操作)
    // ==========================================
    a = new ParBilinearForm(&fespace);
    a->AddDomainIntegrator(new DiffusionIntegrator);
    a->Assemble();
    a->Finalize();
    
    A = new HypreParMatrix();
    a->FormSystemMatrix(ess_tdof_list, *A);

    // ==========================================
    // 3. 质量矩阵 M 组装 (一次性操作)
    // ==========================================
    m = new ParBilinearForm(&fespace);
    m->AddDomainIntegrator(new MassIntegrator());
    m->Assemble();
    m->Finalize();

    // ==========================================
    // 4. RHS 对象初始化
    // ==========================================
    b = new ParLinearForm(&fespace);

    // ==========================================
    // 5. 求解器配置
    // ==========================================
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
    // 1. 设置边界条件值 (基于传入的各组电位配置)
    // ==========================================
    phi = 0.0;
    int max_attr = ess_bdr.Size();

    for (const auto& pair : bdr_conditions_list) {
        if (pair.first.empty()) continue; // 防御性检查
        mfem::Array<int> bdr_marker(max_attr);
        bdr_marker = 0;
        for (int attr : pair.first) {
            if (attr >= 1 && attr <= max_attr) {
                bdr_marker[attr-1] = 1;
            }
        }
        mfem::ConstantCoefficient voltage_coeff(pair.second);
        
        // 【核心修复】直接传入对象本身（即传引用），千万不要加 & 取地址！
        phi.ProjectBdrCoefficient(voltage_coeff, bdr_marker);
    }

    // ==========================================
    // 2. 极速组装 RHS : b = M * rho
    // ==========================================
    m->Mult(rho, *b);

    // ==========================================
    // 3. 生成线性系统 RHS
    // ==========================================
    // 此步骤会处理 Dirichlet 边界条件对真正右端项 B 的反向代入影响
    a->FormLinearSystem(ess_tdof_list, phi, *b, *A, X, B);

    // ==========================================
    // 4. 求解与恢复
    // ==========================================
    X = 0.0; // 初始猜测
    cg->Mult(B, X);
    a->RecoverFEMSolution(X, *b, phi);
}