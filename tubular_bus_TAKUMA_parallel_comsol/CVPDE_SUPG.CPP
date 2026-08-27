#include "mfem.hpp"
using namespace mfem;
using namespace std;

class CustomSUPGIntegrator : public BilinearFormIntegrator
{
private:
    const ParGridFunction *nodeE_poisson_gf_ptr;
    const GridFunction *rho_previous_gf_ptr;

public:
    CustomSUPGIntegrator(const ParGridFunction *nodeE_poisson_gf,
                         const GridFunction *rho_prev_gf)
        : nodeE_poisson_gf_ptr(nodeE_poisson_gf), rho_previous_gf_ptr(rho_prev_gf) {}

    void AssembleElementMatrix(const FiniteElement &el,
                               ElementTransformation &Trans,
                               DenseMatrix &elmat) override;
};

void CustomSUPGIntegrator::AssembleElementMatrix(const FiniteElement &el,
                                                 ElementTransformation &Trans,
                                                 DenseMatrix &elmat)
{
    int num_nodes = 4;
    int dim = el.GetDim();
    elmat.SetSize(num_nodes, num_nodes);
    elmat = 0.0;

    DenseMatrix grad_N_phys_transposed(num_nodes, dim);
    el.CalcPhysDShape(Trans, grad_N_phys_transposed);
    grad_N_phys_transposed.Transpose();

    double Ve = Trans.Weight() / 6.0;

    int vdim = nodeE_poisson_gf_ptr->ParFESpace()->GetVDim();
    MFEM_VERIFY(vdim == dim, "Error: VDim mismatch.");

    Array<int> dofs_E;
    int element_index = Trans.ElementNo;
    nodeE_poisson_gf_ptr->ParFESpace()->GetElementVDofs(element_index, dofs_E);
    MFEM_VERIFY(dofs_E.Size() == num_nodes * vdim, "Error: DOF size mismatch!");

    Vector E_vec_avg(vdim); E_vec_avg = 0.0;
    for (int i = 0; i < num_nodes; ++i)
    {
        E_vec_avg(0) += (*nodeE_poisson_gf_ptr)(dofs_E[i + 0 * num_nodes]);
        E_vec_avg(1) += (*nodeE_poisson_gf_ptr)(dofs_E[i + 1 * num_nodes]);
        E_vec_avg(2) += (*nodeE_poisson_gf_ptr)(dofs_E[i + 2 * num_nodes]);
    }
    E_vec_avg /= num_nodes;

    Vector rho_prev_local_vals(num_nodes);
    rho_previous_gf_ptr->GetElementDofValues(element_index, rho_prev_local_vals);

    double he = 0.0;
    DenseMatrix coord_nodes_phys(dim, num_nodes);
    const IntegrationRule &ref_nodes_rule = el.GetNodes();
    for (int j = 0; j < num_nodes; ++j)
    {
        Vector trans_coords(dim);
        Trans.Transform(ref_nodes_rule[j], trans_coords);
        for (int d = 0; d < dim; ++d)
        {
            coord_nodes_phys(d, j) = trans_coords(d);
        }
    }
    for (int p = 0; p < num_nodes; ++p)
    {
        for (int q = p + 1; q < num_nodes; ++q)
        {
            double dist_sq = 0.0;
            for (int d = 0; d < dim; ++d)
                dist_sq += pow(coord_nodes_phys(d, p) - coord_nodes_phys(d, q), 2);
            he = max(he, sqrt(dist_sq));
        }
    }

    double E_mag = E_vec_avg.Norml2();
    double tau_supg = (E_mag > 1e-12) ? 0.5 * he / E_mag : 0.0;

    DenseMatrix ae_convection(num_nodes, num_nodes);
    for (int i = 0; i < num_nodes; i++)
        for (int j = 0; j < num_nodes; j++)
            ae_convection(i, j) = (E_vec_avg * grad_N_phys_transposed.GetColumn(j)) * (Ve / 4.0);

    DenseMatrix ae_reaction(num_nodes, num_nodes);
    for (int i = 0; i < num_nodes; i++)
        for (int j = 0; j < num_nodes; j++)
        {
            double val = 0.0;
            for (int k = 0; k < num_nodes; k++)
            {
                double coeff = (i == j && j == k) ? Ve / 20.0 :
                               ((i == j || i == k || j == k) ? Ve / 60.0 : Ve / 120.0);
                val += rho_prev_local_vals(k) * coeff;
            }
            ae_reaction(i, j) = val;
        }

    DenseMatrix ae_convection_supg(num_nodes, num_nodes);
    for (int i = 0; i < num_nodes; i++)
        for (int j = 0; j < num_nodes; j++)
            ae_convection_supg(i, j) = tau_supg * (E_vec_avg * grad_N_phys_transposed.GetColumn(j)) *
                                                  (E_vec_avg * grad_N_phys_transposed.GetColumn(i)) * Ve;

    DenseMatrix ae_reaction_supg(num_nodes, num_nodes);
    Vector E_dot_grad_N_tau(num_nodes);
    for (int i = 0; i < num_nodes; i++)
    {
        E_dot_grad_N_tau(i) = tau_supg * (E_vec_avg * grad_N_phys_transposed.GetColumn(i));
        for (int j = 0; j < num_nodes; j++)
        {
            double sum_coeff = 0.0;
            for (int m = 0; m < num_nodes; m++)
            {
                double coeff = (m == j) ? Ve / 10.0 : Ve / 20.0;
                sum_coeff += rho_prev_local_vals(m) * coeff;
            }
            ae_reaction_supg(i, j) = E_dot_grad_N_tau(i) * sum_coeff;
        }
    }

    for (int i = 0; i < num_nodes; i++)
        for (int j = 0; j < num_nodes; j++)
            elmat(i, j) = ae_convection(i,j) + ae_reaction(i,j)
                        + ae_convection_supg(i,j) + ae_reaction_supg(i,j);
}


// 封装为函数：用于主程序中调用
void SolveConvectionSUPG(ParFiniteElementSpace &scalar_fes,
                         const ParGridFunction &nodeE_poisson,
                         const GridFunction &rho_previous_gf,
                         double rho_0,
                         ParGridFunction &rho_new,
                         Array<int> &ess_bdr)
{
    rho_new = rho_0;

    ParBilinearForm a_con(&scalar_fes);
    a_con.AddDomainIntegrator(new CustomSUPGIntegrator(&nodeE_poisson, &rho_previous_gf));
    a_con.Assemble();

    ParLinearForm b_con(&scalar_fes);
    b_con = 0.0;
    b_con.Assemble();

    Array<int> ess_tdof_list;
    scalar_fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    HypreParMatrix A_con;
    Vector B_con, X_con;
    a_con.FormLinearSystem(ess_tdof_list, rho_new, b_con, A_con, X_con, B_con);

    HypreBoomerAMG M(A_con);
    M.SetPrintLevel(0);
    GMRESSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(0);
    cg.SetPreconditioner(M);
    cg.SetOperator(A_con);
    cg.Mult(B_con, X_con);
    
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // ofstream rho_new_output("rho_new" + std::to_string(rank) + ".txt");
    // rho_new_output << std::scientific << std::setprecision(12); // 例如：1.23456789e+02
    // for(int i = 0; i < X_con.Size(); i++)
    // {
    //     rho_new_output << X_con[i] <<endl;
    // }
    // rho_new_output.close();

    a_con.RecoverFEMSolution(X_con, b_con, rho_new);
}
