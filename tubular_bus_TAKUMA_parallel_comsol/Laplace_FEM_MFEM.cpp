#include "Laplace_FEM_MFEM.h"

mfem::ParGridFunction Laplace_FEM_MFEM(const mfem::ParMesh& mesh, double phi_0, 
                                    mfem::ParFiniteElementSpace& fespace) {
    // 创建 dir_bdr 数组
    mfem::Array<int> dir_bdr(mesh.bdr_attributes.Max());
    dir_bdr = 0; // 初始化为0（无边界条件）

    // 标记需要设置边界条件的属性
    dir_bdr[2] = 1; // 属性 3
    dir_bdr[3] = 1; // 属性 4

    // 定义两个边界属性数组
    mfem::Array<int> inner_bdr(mesh.bdr_attributes.Max()), outer_bdr(mesh.bdr_attributes.Max());
    inner_bdr = 0;
    outer_bdr = 0;
    inner_bdr[2] = 1; // 标记内边界（属性 3）
    outer_bdr[3] = 1; // 标记外边界（属性 4）

    // 创建 GridFunction 并初始化
    mfem::ParGridFunction phi(&fespace);
    phi = 0.0; // 初始化内部区域为0

    // 设置内边界（属性3）= phi_0
    mfem::ConstantCoefficient inner_coeff(phi_0);
    phi.ProjectBdrCoefficient(inner_coeff, inner_bdr);

    // 设置外边界（属性4）= 0   
    mfem::ConstantCoefficient outer_coeff(0.0);
    phi.ProjectBdrCoefficient(outer_coeff, outer_bdr);

    // 组装双线性形式和线性形式
    mfem::ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator); // -∇²φ
    a.Assemble();

    mfem::ParLinearForm b(&fespace);
    mfem::ConstantCoefficient zero_coeff(0.0);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(zero_coeff)); // 0
    b.Assemble();

    // 获取 essential true DOFs
    mfem::Array<int> ess_tdof_list;
    fespace.GetEssentialTrueDofs(dir_bdr, ess_tdof_list);

    // 形成线性系统
    mfem::HypreParMatrix A;
    mfem::Vector B, X;
    a.FormLinearSystem(ess_tdof_list, phi, b, A, X, B);

    // 求解线性系统
    mfem::HypreBoomerAMG M(A);
    M.SetPrintLevel(0);
    mfem::CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(0);
    cg.SetPreconditioner(M);
    cg.SetOperator(A);
    cg.Mult(B, X);

    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::ofstream filename("laplace_rank_" + std::to_string(rank) + "_.csv");
    // for (int i = 0 ; i < X.Size() ; i ++){
    //     filename << X[ i ] <<std::endl;
    // }
    // filename.close();

    // 恢复解
    a.RecoverFEMSolution(X, b, phi);

    return phi; // 返回设置好边界条件并求解后的 phi
}