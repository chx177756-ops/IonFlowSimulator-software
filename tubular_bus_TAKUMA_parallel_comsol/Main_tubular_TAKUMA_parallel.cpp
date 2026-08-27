#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <mpi.h>
#include <algorithm> 
#include <sys/stat.h> 
#include <iomanip>   

#include "mfem.hpp"

#include "ComputeElementVolume.hpp"
#include "GetNodeE.h"
#include "poisson_FEM_MFEM.h" 
#include "ParallelUpwindElementFinder.h"
#include "CVPDE_TAKUMA.h" 
#include "ReadComsolMesh.hpp" // [引入] COMSOL 网格解析器

#include <fstream>
#include <string>

double GetPeakMemoryMB() {
    std::ifstream status_file("/proc/self/status");
    std::string line;
    while (std::getline(status_file, line)) {
        if (line.substr(0, 6) == "VmHWM:") {
            long kb = 0;
            sscanf(line.c_str(), "VmHWM: %ld kB", &kb);
            return kb / 1024.0; 
        }
    }
    return 0.0;
}

using namespace std;
using namespace mfem;

struct PerfTimer {
    double total_time = 0.0;
    int call_count = 0;
    void Start() { start_t = MPI_Wtime(); }
    void Stop() { 
        total_time += (MPI_Wtime() - start_t); 
        call_count++; 
    }
private:
    double start_t = 0.0;
};

int main(int argc, char *argv[])
{
    Mpi::Init(argc, argv);
    Hypre::Init();
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank = 0, num_procs = 1;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    double program_start_time = MPI_Wtime();

    PerfTimer t_mesh_read;      
    PerfTimer t_init_setup;     
    PerfTimer t_poisson;        
    PerfTimer t_grad_calc;      
    PerfTimer t_upwind_find;    
    PerfTimer t_takuma;         
    PerfTimer t_output;         
    PerfTimer t_loop_overhead;  

    double rela_error_rho = 1.0, max_rela_error_rho = 1.0, rela_error_E = 1.0;
    int count_iteration = 0;
    
    std::vector<double> rela_error_rho_history;
    std::vector<double> rela_error_E_history;
    std::vector<double> max_abs_error_rho_history;
    std::vector<double> rate_convergence_history;
    std::vector<double> rate_update_history;
    std::vector<int> iterations;

    double w_rho = 1;
    double rate_convergence = 0.0;
    double Goal_convergence = 0.95;

    double E0 = 6e5;                    
    double rho_surface = 1e4;
    const double epsilon_0 = 8.8542e-12; 
    
    const double K_mobility = 1; 
    Vector w_vec(3); w_vec = 0.0;     

    string mesh_file = "MESH/pipe_I_ybj_neg2_clear.mphtxt";
    string output_folder = "results";
    int order = 1;
    
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&output_folder, "-otpt", "--output", "results saving folder.");
    args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
    args.AddOption(&E0, "-e0", "--e-onset", "Corona onset field strength (V/m)");
    args.AddOption(&rho_surface, "-rho", "--rho-surface", "Initial surface charge density"); 
    
    args.ParseCheck();

    if (my_rank == 0) {
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Current E0 (Onset Field): " << E0 << " V/m" << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Reading COMSOL mesh..." << std::endl;
    }
    
    // =====================================================================
    // 调用 COMSOL 解析器
    // =====================================================================
    t_mesh_read.Start();
    Mesh *serial_mesh_ptr = nullptr;
    std::map<std::string, std::vector<int>> comsol_selections; 

    if (my_rank == 0) {
        serial_mesh_ptr = ReadComsolMesh(mesh_file, comsol_selections);
        
        std::cout << "\n--- COMSOL Selections Auto-Detected ---" << std::endl;
        for (const auto& pair : comsol_selections) {
            std::cout << "Label: [" << pair.first << "] -> Entities: { ";
            for (int val : pair.second) std::cout << val << " ";
            std::cout << "}\n";
        }
        std::cout << "---------------------------------------\n\n";
    }
    
    if (serial_mesh_ptr == nullptr) {
        serial_mesh_ptr = ReadComsolMesh(mesh_file, comsol_selections);
    }

    ParMesh mesh(comm, *serial_mesh_ptr);
    
    if (my_rank != 0) {
        delete serial_mesh_ptr;
        serial_mesh_ptr = nullptr;
    }
    t_mesh_read.Stop();

    // === 计时：初始化设置 ===
    t_init_setup.Start();

    H1_FECollection scalar_fec(order, mesh.Dimension());
    ParFiniteElementSpace scalar_fes(&mesh, &scalar_fec); 
    H1_FECollection vector_fec(order, mesh.Dimension());
    ParFiniteElementSpace vector_fes(&mesh, &vector_fec, 3); 

    std::vector<double> Volume = ComputeElementVolume(mesh);

    if (my_rank == 0) std::cout << "Precomputing physical shape function gradients..." << std::endl;
    int numTet = mesh.GetNE();
    int dim = mesh.Dimension();
    std::vector<mfem::DenseMatrix> global_dshapes(numTet);
    mfem::DenseMatrix dshape_ref;
    for (int e = 0; e < numTet; ++e) {
        int ndof_e = scalar_fes.GetFE(e)->GetDof();
        global_dshapes[e].SetSize(ndof_e, dim);
        mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
        const mfem::IntegrationRule &ir = mfem::IntRules.Get(mesh.GetElementType(e), 0);
        T->SetIntPoint(&ir.IntPoint(0));
        const mfem::DenseMatrix &Jinv = T->InverseJacobian();
        scalar_fes.GetFE(e)->CalcDShape(ir.IntPoint(0), dshape_ref);
        mfem::Mult(dshape_ref, Jinv, global_dshapes[e]);
    }

    // =====================================================================
    // [新模型重构] 自动映射 COMSOL 物理组到代码内部边界变量
    // =====================================================================
    auto get_selection = [&](const std::string& name) -> std::vector<int> {
        if (comsol_selections.find(name) != comsol_selections.end()) {
            return comsol_selections[name];
        } else {
            if(my_rank == 0) std::cout << "WARNING: Label '" << name << "' not found in COMSOL file!\n";
            return {};
        }
    };

    // 提取新模型对应的所有部件 (COMSOL selection label)
    const std::vector<int> ATTR_T1_TOWER    = get_selection("T1 Tower");
    const std::vector<int> ATTR_T2_TOWER    = get_selection("T2 Tower");
    const std::vector<int> ATTR_FLOATING    = get_selection("Floating Conductor");
    const std::vector<int> ATTR_ARTIFICIAL  = get_selection("Artificial_boundry");

    const std::vector<int> ATTR_TUBULAR     = get_selection("Tubular");
    const std::vector<int> ATTR_SPHERES     = get_selection("spheres rings connection");
    
    auto collect_boundary_dofs = [&](const std::vector<int>& target_attrs) -> std::vector<int> {
        std::set<int> s;
        if (target_attrs.empty()) return {}; 
        for (int i = 0; i < mesh.GetNBE(); ++i) {
            Element *bel = mesh.GetBdrElement(i);
            if (!bel) continue;
            
            int attr = bel->GetAttribute();
            if (std::find(target_attrs.begin(), target_attrs.end(), attr) != target_attrs.end()) {
                Array<int> dofs;
                scalar_fes.GetBdrElementDofs(i, dofs); 
                for (int j = 0; j < dofs.Size(); ++j) s.insert(dofs[j]);
            }
        }
        return std::vector<int>(s.begin(), s.end());
    };

    // 定义电位边界配置向量：支持任意多种不同电位
    std::vector<std::pair<std::vector<int>, double>> bdr_conditions;
    
    // 配置 0V 边界 (T1 Tower、T2 Tower、Floating Conductor、Artificial_boundry)
    std::vector<int> bdr_0V;
    bdr_0V.insert(bdr_0V.end(), ATTR_T1_TOWER.begin(), ATTR_T1_TOWER.end());
    bdr_0V.insert(bdr_0V.end(), ATTR_T2_TOWER.begin(), ATTR_T2_TOWER.end());
    bdr_0V.insert(bdr_0V.end(), ATTR_FLOATING.begin(), ATTR_FLOATING.end());
    bdr_0V.insert(bdr_0V.end(), ATTR_ARTIFICIAL.begin(), ATTR_ARTIFICIAL.end());
    bdr_conditions.push_back({bdr_0V, 0.0});

    // 配置 1000kV 边界 (Tubular、spheres rings connection)
    std::vector<int> bdr_1000kV;
    bdr_1000kV.insert(bdr_1000kV.end(), ATTR_TUBULAR.begin(), ATTR_TUBULAR.end());
    bdr_1000kV.insert(bdr_1000kV.end(), ATTR_SPHERES.begin(), ATTR_SPHERES.end());
    bdr_conditions.push_back({bdr_1000kV, 1000e3});

    int numLocalDofs = scalar_fes.GetNDofs();
    Vector rho_local(numLocalDofs);
    rho_local = 0.0;

    // 使用更新后的 PoissonSolver
    PoissonSolver poisson_solver(mesh, scalar_fes, bdr_conditions);

    t_init_setup.Stop(); 

    // 初始化无空间电荷时的泊松场电势分布
    t_poisson.Start();
    ParGridFunction phi0(&scalar_fes);
    poisson_solver.Solve(rho_local, phi0);
    t_poisson.Stop();

    t_init_setup.Start(); 

    DenseMatrix E_elem(mesh.GetNE(), 3);
    for (int e = 0; e < mesh.GetNE(); ++e) {
        ElementTransformation *tr = mesh.GetElementTransformation(e);
        Vector grad;
        phi0.GetGradient(*tr, grad);
        E_elem(e,0) = -grad(0); E_elem(e,1) = -grad(1); E_elem(e,2) = -grad(2);
    }
    
    ParGridFunction nodeE(&vector_fes);
    ParGridFunction nodeEn(&scalar_fes);
    GetNodeE(E_elem, mesh, Volume, vector_fes, nodeE, nodeEn);

    std::vector<double> nodeEn_initial(numLocalDofs, 0.0);
    for (int ld = 0; ld < numLocalDofs; ++ld) nodeEn_initial[ld] = (ld < nodeEn.Size()) ? nodeEn(ld) : 0.0;

    // ================================================================
    // [起晕判断] 在 Tubular、T1 Tower、T2 Tower 上检测表面电场是否超过起晕场强
    // ================================================================
    std::vector<int> corona_check_attrs;
    corona_check_attrs.insert(corona_check_attrs.end(), ATTR_TUBULAR.begin(), ATTR_TUBULAR.end());
    corona_check_attrs.insert(corona_check_attrs.end(), ATTR_T1_TOWER.begin(), ATTR_T1_TOWER.end());
    corona_check_attrs.insert(corona_check_attrs.end(), ATTR_T2_TOWER.begin(), ATTR_T2_TOWER.end());
    std::vector<int> vertex_corona_check = collect_boundary_dofs(corona_check_attrs);
    
    std::vector<char> is_corona_local(numLocalDofs, 0);
    for (int idx : vertex_corona_check) {
        if (idx >= 0 && idx < numLocalDofs) {
            if (nodeEn_initial[idx] > E0) is_corona_local[idx] = 1;
        }
    }
    
    GroupCommunicator &gcomm = scalar_fes.GroupComm();
    Vector is_corona_vec(scalar_fes.GetVSize());
    for (int i = 0; i < scalar_fes.GetVSize(); ++i) is_corona_vec(i) = is_corona_local[i] ? 1.0 : 0.0;
    gcomm.Reduce<double>(is_corona_vec.HostReadWrite(), GroupCommunicator::Max);
    gcomm.Bcast<double>(is_corona_vec.HostReadWrite());

    std::vector<int> vertex_corona;
    for (int i = 0; i < numLocalDofs; ++i) {
        if (is_corona_vec(i) > 0.5) vertex_corona.push_back(i);
    }

    long long global_corona_count = 0;
    long long local_corona_count = vertex_corona.size();
    MPI_Allreduce(&local_corona_count, &global_corona_count, 1, MPI_LONG_LONG, MPI_SUM, comm);

    // ================================================================
    // [初值设定] 准备所有金属边界，以防电荷在金属内部异常累积
    // (Artificial_boundry 是人工计算域边界, 非金属, 不纳入)
    // ================================================================
    std::vector<int> all_inner_metals;
    all_inner_metals.insert(all_inner_metals.end(), ATTR_T1_TOWER.begin(), ATTR_T1_TOWER.end());
    all_inner_metals.insert(all_inner_metals.end(), ATTR_T2_TOWER.begin(), ATTR_T2_TOWER.end());
    all_inner_metals.insert(all_inner_metals.end(), ATTR_FLOATING.begin(), ATTR_FLOATING.end());
    all_inner_metals.insert(all_inner_metals.end(), ATTR_TUBULAR.begin(), ATTR_TUBULAR.end());
    all_inner_metals.insert(all_inner_metals.end(), ATTR_SPHERES.begin(), ATTR_SPHERES.end());
    std::vector<int> vertex_inner_metals = collect_boundary_dofs(all_inner_metals);

    bool run_main_loop = true;
    if (global_corona_count == 0)
    {
        if (my_rank == 0) {
            std::cout << "==========================================================" << std::endl;
            std::cout << "NO CORONA INCEPTION DETECTED (Max Surf E < E0)." << std::endl;
            std::cout << "Skipping iterative calculation loop." << std::endl;
            std::cout << "==========================================================" << std::endl;
        }
        run_main_loop = false;
    }
    else
    {
        for (int i = 0; i < numLocalDofs; ++i) rho_local(i) = rho_surface / 100.0;
        // 所有内部金属件赋微小本底值
        for (int idx : vertex_inner_metals) if (idx >=0 && idx < numLocalDofs) rho_local(idx) = 1e-12;
        // 唯独覆盖满足起晕条件的节点 (Tubular + T1/T2 Tower)
        for (int idx : vertex_corona) if (idx >=0 && idx < numLocalDofs) rho_local(idx) = rho_surface;
    }
    t_init_setup.Stop(); 

    int set_updatetimes = 100;
    double tolerance_E = 0.01;
    double tolerance_rho = 0.01;
    int count_update = 0;

    ParGridFunction v_total(&vector_fes);

    if (run_main_loop) 
    {
        if (my_rank == 0) std::cout << "Start iterative loop (TAKUMA Method)..." << std::endl;

        while (rela_error_E >= tolerance_E)
        {
            ++count_iteration;
            iterations.push_back(count_iteration);
            if (my_rank == 0) std::cout << "Iteration: " << count_iteration << std::endl;

            Vector rho_prev = rho_local;

            t_poisson.Start();
            ParGridFunction phi(&scalar_fes);
            poisson_solver.Solve(rho_local, phi);
            t_poisson.Stop();

            t_grad_calc.Start();
            // [极致优化] 利用预计算全局 dShape 直接乘加
            Array<int> dofs_cache;
            for (int e = 0; e < mesh.GetNE(); ++e) {
                scalar_fes.GetElementDofs(e, dofs_cache);
                const mfem::DenseMatrix &dshape_e = global_dshapes[e];
                int ndof_e = dofs_cache.Size();
                
                double gx = 0.0, gy = 0.0, gz = 0.0;
                for (int i = 0; i < ndof_e; ++i) {
                    double phi_val = phi(dofs_cache[i]); 
                    gx += phi_val * dshape_e(i, 0);
                    gy += phi_val * dshape_e(i, 1);
                    gz += phi_val * dshape_e(i, 2);
                }
                E_elem(e,0) = -gx; E_elem(e,1) = -gy; E_elem(e,2) = -gz;
            }
            GetNodeE(E_elem, mesh, Volume, vector_fes, nodeE, nodeEn);

            v_total = nodeE;       
            v_total *= K_mobility; 
            VectorConstantCoefficient w_coeff(w_vec);
            ParGridFunction w_gf(&vector_fes);
            w_gf.ProjectCoefficient(w_coeff);
            v_total += w_gf;       
            t_grad_calc.Stop();

            t_upwind_find.Start();
            UpwindElementFinder finder(mesh, v_total, vector_fes, global_dshapes); 
            auto Upwindelements = finder.ComputeUpwindElements();
            t_upwind_find.Stop();

            t_takuma.Start();
            Vector rho_calc_local = rho_local; 
            
            SolveContinuityTakuma_3D(
                mesh, scalar_fes, vector_fes, 
                rho_local, phi, nodeE, 
                Upwindelements, rho_calc_local, vertex_corona, 
                K_mobility, w_vec, epsilon_0,
                global_dshapes
            );
            t_takuma.Stop();

            t_loop_overhead.Start();
            
            for(int i=0; i<numLocalDofs; ++i) {
                rho_local(i) = w_rho * rho_calc_local(i) + (1.0 - w_rho) * rho_prev(i);
            }
            scalar_fes.GroupComm().Bcast<double>(rho_local.GetData());

            std::vector<double> rel_rho_local_vec(numLocalDofs, 0.0);
            double max_rel_local = 0.0;
            double max_abs_local = 0.0;
            const double eps_denom = 1e-16; 
            const double ABS_THRESHOLD = 1e-3; 

            for (int i = 0; i < numLocalDofs; ++i)
            {
                double val_calc = std::abs(rho_local(i));
                double val_prev = std::abs(rho_prev(i));
                double diff = std::abs(rho_local(i) - rho_prev(i));
                
                double rel = 0.0;
                if (val_calc < ABS_THRESHOLD && val_prev < ABS_THRESHOLD) {
                    rel = 0.0;
                }
                else {
                    rel = diff / (val_prev + eps_denom);
                }

                rel_rho_local_vec[i] = rel;
                if (rel > max_rel_local) max_rel_local = rel;
                if (diff > max_abs_local) max_abs_local = diff;
            }

            MPI_Allreduce(&max_rel_local, &max_rela_error_rho, 1, MPI_DOUBLE, MPI_MAX, comm);
            double max_abs_global = 0.0;
            MPI_Allreduce(&max_abs_local, &max_abs_global, 1, MPI_DOUBLE, MPI_MAX, comm);

            const HYPRE_BigInt *trueOffsets = scalar_fes.GetTrueDofOffsets();
            HYPRE_BigInt my_true_lo = trueOffsets[my_rank];
            HYPRE_BigInt my_true_hi = trueOffsets[my_rank+1] - 1;
            long long local_count_below = 0;
            long long local_total_true = 0;
            
            for (int ld = 0; ld < numLocalDofs; ++ld)
            {
                HYPRE_BigInt gd = scalar_fes.GetGlobalTDofNumber(ld);
                if (gd >= my_true_lo && gd <= my_true_hi)
                {
                    ++local_total_true; 
                    if (rel_rho_local_vec[ld] < tolerance_rho) ++local_count_below;
                }
            }
            long long global_count_below = 0;
            long long global_total_true = 0;
            MPI_Allreduce(&local_count_below, &global_count_below, 1, MPI_LONG_LONG, MPI_SUM, comm);
            MPI_Allreduce(&local_total_true, &global_total_true, 1, MPI_LONG_LONG, MPI_SUM, comm);

            if (global_total_true > 0)
                rate_convergence = double(global_count_below) / double(global_total_true);
            else
                rate_convergence = 0.0;

            rela_error_rho_history.push_back(max_rela_error_rho);
            max_abs_error_rho_history.push_back(max_abs_global);
            rate_convergence_history.push_back(rate_convergence);

            if (my_rank == 0) {
                printf("after relaxation: max_rela_error_rho=%.3e, max_abs_rho=%.3e, convergence_rate=%.3f\n",
                    max_rela_error_rho, max_abs_global, rate_convergence);
            }

            if (rate_convergence >= Goal_convergence)
            {
                t_loop_overhead.Stop(); 
                t_poisson.Start();      
                
                ParGridFunction phi_new(&scalar_fes);
                poisson_solver.Solve(rho_local, phi_new);
                
                t_poisson.Stop();       
                t_grad_calc.Start();    

                for (int e = 0; e < mesh.GetNE(); ++e) {
                    ElementTransformation *tr = mesh.GetElementTransformation(e);
                    Vector grad;
                    phi_new.GetGradient(*tr, grad);
                    E_elem(e,0) = -grad(0); E_elem(e,1) = -grad(1); E_elem(e,2) = -grad(2);
                }
                GetNodeE(E_elem, mesh, Volume, vector_fes, nodeE, nodeEn);
                
                t_grad_calc.Stop();     
                t_loop_overhead.Start(); 

                double max_E_local = 0.0;
                double min_E_local = 1e300;
                double sum_E_local = 0.0;
                int ncor_local = 0;

                for (int idx : vertex_corona)
                {
                    if (idx < 0 || idx >= numLocalDofs) continue;
                    double val = nodeEn(idx);
                    if (val > max_E_local) max_E_local = val;
                    if (val < min_E_local) min_E_local = val;
                    sum_E_local += val;
                    ++ncor_local;
                }

                double max_E_global = 0.0;
                double min_E_global = 0.0;
                double sum_E_global = 0.0;
                int ncor_global = 0;

                MPI_Allreduce(&max_E_local, &max_E_global, 1, MPI_DOUBLE, MPI_MAX, comm);
                MPI_Allreduce(&min_E_local, &min_E_global, 1, MPI_DOUBLE, MPI_MIN, comm);
                MPI_Allreduce(&sum_E_local, &sum_E_global, 1, MPI_DOUBLE, MPI_SUM, comm);
                MPI_Allreduce(&ncor_local, &ncor_global, 1, MPI_INT, MPI_SUM, comm);

                double mean_E_global = (ncor_global > 0) ? (sum_E_global / ncor_global) : 0.0;

                if (my_rank == 0) {
                printf("corona E: max=%.6e, min=%.6e, mean=%.6e\n", max_E_global, min_E_global, mean_E_global);
                }

                double E_standard = max_E_global; 

                if (E0 != 0.0) rela_error_E = std::abs(E_standard - E0) / std::abs(E0);
                else rela_error_E = std::abs(E_standard - E0);

                rela_error_E_history.push_back(rela_error_E);
                if (my_rank == 0) printf("rela_error_E (using max_E) = %.6e\n", rela_error_E);

                if (rela_error_E < tolerance_E)
                {
                    if (my_rank == 0) printf("Converged after %d iterations.\n", count_iteration);
                    t_loop_overhead.Stop();
                    break;
                }
                else
                {
                    if (my_rank == 0) printf("Updating corona rho (count -> %d).\n", count_update+1);
                    ++count_update;

                    int miu = 1;
                    if (rela_error_E > 0.2) miu = 8;
                    else if (rela_error_E > 0.1) miu = 4;
                    else if (rela_error_E > 0.01) miu = 2;
                    
                    double rate_update = 0.0;
                    if (E_standard + E0 != 0.0) rate_update = miu * (E_standard - E0) / (E_standard + E0);
                    rate_update_history.push_back(rate_update);

                    for (int idx : vertex_corona) {
                        if (idx >= 0 && idx < numLocalDofs) {
                            rho_local(idx) = rho_local(idx) * (1.0 + rate_update);
                        }
                    }
                }
            } 

            if (my_rank == 0) printf("End Iter %d\n", count_iteration);
            
            t_loop_overhead.Stop();

            if (count_update >= set_updatetimes) break;
        }

        if (my_rank == 0) std::cout << "Iteration loop finished." << std::endl;
    } 

t_output.Start();
    {
        int myid = my_rank;
        int nprocs = num_procs;
        
        // 1. 计算总场电势 phi_final
        ParGridFunction phi_final(&scalar_fes);
        poisson_solver.Solve(rho_local, phi_final);

        // 2. 计算总场强
        DenseMatrix E_elem_final(mesh.GetNE(), 3);
        for (int e = 0; e < mesh.GetNE(); ++e) {
            ElementTransformation *tr = mesh.GetElementTransformation(e);
            Vector grad;
            phi_final.GetGradient(*tr, grad);
            E_elem_final(e,0) = -grad(0); E_elem_final(e,1) = -grad(1); E_elem_final(e,2) = -grad(2);
        }

        ParGridFunction nodeE_final(&vector_fes);
        ParGridFunction nodeEn_final(&scalar_fes);
        GetNodeE(E_elem_final, mesh, Volume, vector_fes, nodeE_final, nodeEn_final);

        // 3. 建立 顶点到自由度 的映射
        int localNV = mesh.GetNV();
        std::vector<int> vertex_to_dof(localNV, -1);

        for (int e = 0; e < mesh.GetNE(); ++e) {
            Array<int> el_dofs;
            scalar_fes.GetElementDofs(e, el_dofs);
            Element *el = mesh.GetElement(e);
            const int *verts = el->GetVertices();
            int nv = el->GetNVertices();
            int nd = el_dofs.Size();
            int min_nd_nv = std::min(nd, nv);
            for (int j = 0; j < min_nd_nv; ++j) {
                int v_local = verts[j];
                int dof_local = el_dofs[j];
                if (v_local >= 0 && v_local < localNV) vertex_to_dof[v_local] = dof_local;
            }
        }

        // 4. 提取本地顶点数据 (新增了 send_V0)
        int my_localNV = localNV;
        std::vector<double> send_coords(3 * my_localNV);
        std::vector<double> send_rho(my_localNV);
        std::vector<double> send_E(my_localNV);
        std::vector<double> send_E0(my_localNV);      
        std::vector<double> send_V0(my_localNV); // [新增] 标称场电压 V0

        for (int i = 0; i < my_localNV; ++i) {
            const double *c = mesh.GetVertex(i);
            send_coords[3*i + 0] = c[0]; send_coords[3*i + 1] = c[1]; send_coords[3*i + 2] = c[2];

            int dof = vertex_to_dof[i];
            if (dof >= 0 && dof < rho_local.Size()) send_rho[i] = rho_local(dof); else send_rho[i] = 0.0;
            if (dof >= 0 && dof < nodeEn_final.Size()) send_E[i] = nodeEn_final(dof); else send_E[i] = 0.0;
            if (dof >= 0 && dof < (int)nodeEn_initial.size()) send_E0[i] = nodeEn_initial[dof]; else send_E0[i] = 0.0;
            // [新增] 提取 phi0(标称场电压)
            if (dof >= 0 && dof < phi0.Size()) send_V0[i] = phi0(dof); else send_V0[i] = 0.0;
        }

        // 5. MPI 收集信息准备
        std::vector<int> all_NV_counts(nprocs, 0);
        MPI_Allgather(&my_localNV, 1, MPI_INT, all_NV_counts.data(), 1, MPI_INT, comm);

        std::vector<int> displs(nprocs, 0);
        int total_entries = 0;
        for (int r = 0; r < nprocs; ++r) {
            displs[r] = total_entries;
            total_entries += all_NV_counts[r];
        }

        std::vector<int> recv_counts_coords(nprocs), recv_displs_coords(nprocs);
        for (int r = 0; r < nprocs; ++r) {
            recv_counts_coords[r] = all_NV_counts[r] * 3;
            recv_displs_coords[r] = displs[r] * 3;
        }
        
        // 6. 执行 MPI 全局收集
        std::vector<double> recv_coords(total_entries * 3);
        MPI_Allgatherv(send_coords.data(), 3*my_localNV, MPI_DOUBLE, recv_coords.data(), recv_counts_coords.data(), recv_displs_coords.data(), MPI_DOUBLE, comm);

        std::vector<int> recv_counts_scalar = all_NV_counts;
        std::vector<int> recv_displs_scalar = displs;
        
        std::vector<double> recv_rho(total_entries);
        MPI_Allgatherv(send_rho.data(), my_localNV, MPI_DOUBLE, recv_rho.data(), recv_counts_scalar.data(), recv_displs_scalar.data(), MPI_DOUBLE, comm);

        std::vector<double> recv_E(total_entries);
        MPI_Allgatherv(send_E.data(), my_localNV, MPI_DOUBLE, recv_E.data(), recv_counts_scalar.data(), recv_displs_scalar.data(), MPI_DOUBLE, comm);

        std::vector<double> recv_E0(total_entries);
        MPI_Allgatherv(send_E0.data(), my_localNV, MPI_DOUBLE, recv_E0.data(), recv_counts_scalar.data(), recv_displs_scalar.data(), MPI_DOUBLE, comm);

        // [新增] 收集标称场电压 V0
        std::vector<double> recv_V0(total_entries);
        MPI_Allgatherv(send_V0.data(), my_localNV, MPI_DOUBLE, recv_V0.data(), recv_counts_scalar.data(), recv_displs_scalar.data(), MPI_DOUBLE, comm);

        // 7. Rank 0 节点负责数据匹配与保存
        if (myid == 0 && serial_mesh_ptr != nullptr)
        {
            std::cout << "[Output] Gathering data on Rank 0..." << std::endl;
            int NV_serial = serial_mesh_ptr->GetNV(); 
            
            // [修改] 在结构体中增加 V0
            struct PointData {
                int id;         
                double x, y, z; 
                double rho, E, E0, V0; 
            };

            std::vector<PointData> parallel_points(total_entries);
            for(int i = 0; i < total_entries; ++i) {
                parallel_points[i].id = -1;
                parallel_points[i].x = recv_coords[3*i+0];
                parallel_points[i].y = recv_coords[3*i+1];
                parallel_points[i].z = recv_coords[3*i+2];
                parallel_points[i].rho = recv_rho[i];
                parallel_points[i].E   = recv_E[i];
                parallel_points[i].E0  = recv_E0[i];
                parallel_points[i].V0  = recv_V0[i]; // [新增]
            }

            std::vector<PointData> serial_points(NV_serial);
            for(int i = 0; i < NV_serial; ++i) {
                const double* v = serial_mesh_ptr->GetVertex(i); 
                serial_points[i].id = i; 
                serial_points[i].x = v[0];
                serial_points[i].y = v[1];
                serial_points[i].z = v[2];
                serial_points[i].rho = 0.0; 
                serial_points[i].V0  = 0.0; // [新增] 初始化
            }

            auto comp_strict = [](const PointData& a, const PointData& b) {
                if (std::isnan(a.x) || std::isnan(b.x)) return false;
                if (std::abs(a.x - b.x) > 1e-14) return a.x < b.x;
                if (std::abs(a.y - b.y) > 1e-14) return a.y < b.y;
                return a.z < b.z;
            };

            std::cout << "[Output] Sorting..." << std::endl;
            std::sort(parallel_points.begin(), parallel_points.end(), comp_strict);
            std::sort(serial_points.begin(), serial_points.end(), comp_strict);

            const double search_tol = 1e-5; 
            const double match_tol_sq = 1e-12; 
            size_t p_start = 0; 
            int matched_count = 0;

            for (size_t s_idx = 0; s_idx < serial_points.size(); ++s_idx) {
                double sx = serial_points[s_idx].x;
                
                while (p_start < parallel_points.size() && parallel_points[p_start].x < sx - search_tol) {
                    p_start++;
                }

                for (size_t k = p_start; k < parallel_points.size(); ++k) {
                    if (parallel_points[k].x > sx + search_tol) break;

                    if (std::abs(parallel_points[k].y - serial_points[s_idx].y) > search_tol) continue;
                    if (std::abs(parallel_points[k].z - serial_points[s_idx].z) > search_tol) continue;

                    double dx = parallel_points[k].x - sx;
                    double dy = parallel_points[k].y - serial_points[s_idx].y;
                    double dz = parallel_points[k].z - serial_points[s_idx].z;
                    
                    if (dx*dx + dy*dy + dz*dz < match_tol_sq) {
                        serial_points[s_idx].rho = parallel_points[k].rho;
                        serial_points[s_idx].E   = parallel_points[k].E;
                        serial_points[s_idx].E0  = parallel_points[k].E0;
                        serial_points[s_idx].V0  = parallel_points[k].V0; // [新增]
                        matched_count++;
                        break; 
                    }
                }
            }
            std::cout << "[Output] Matched " << matched_count << " / " << NV_serial << " vertices." << std::endl;

            std::vector<double> final_rho(NV_serial, 0.0);
            std::vector<double> final_E(NV_serial, 0.0);
            std::vector<double> final_E0(NV_serial, 0.0);
            std::vector<double> final_V0(NV_serial, 0.0); // [新增]
            std::vector<double> final_x(NV_serial, 0.0);
            std::vector<double> final_y(NV_serial, 0.0);
            std::vector<double> final_z(NV_serial, 0.0);
            
            for(const auto& p : serial_points) {
                if (p.id >= 0 && p.id < NV_serial) {
                    final_rho[p.id] = p.rho;
                    final_E[p.id]   = p.E;
                    final_E0[p.id]  = p.E0;
                    final_V0[p.id]  = p.V0; // [新增]
                    final_x[p.id]   = p.x;
                    final_y[p.id]   = p.y;
                    final_z[p.id]   = p.z;
                }
            }

            std::cout << "[Output] Writing files to: " << output_folder << std::endl;
            #ifdef _WIN32
                _mkdir(output_folder.c_str());
            #else
                mkdir(output_folder.c_str(), 0777);
            #endif

            auto safe_write = [&](string filename, const vector<double>& data) {
                string fullpath = output_folder + filename;
                ofstream fout(fullpath);
                if(fout.is_open()) {
                    fout.precision(14); 
                    fout << std::scientific;
                    for(auto v : data) fout << v << "\n";
                    fout.close();
                } else {
                    std::cerr << "Error opening file: " << fullpath << std::endl;
                }
            };
            
            safe_write("/rho.txt", final_rho);
            safe_write("/E.txt", final_E);
            safe_write("/E0.txt", final_E0);
            safe_write("/V0.txt", final_V0); // [新增] 将标称场电压导出为 V0.txt
            safe_write("/x.txt", final_x);
            safe_write("/y.txt", final_y);
            safe_write("/z.txt", final_z);
            
            cout << "Rank 0 wrote results and coordinates successfully." << endl;
        }
    }
    t_output.Stop();

    if (my_rank == 0 && serial_mesh_ptr != nullptr) {
        delete serial_mesh_ptr;
        serial_mesh_ptr = nullptr;
    }

    double program_end_time = MPI_Wtime();
    double total_elapsed = program_end_time - program_start_time;
    
    double local_peak_mem = GetPeakMemoryMB();
    double global_max_peak_mem = 0.0;
    double global_sum_peak_mem = 0.0;
    
    MPI_Reduce(&local_peak_mem, &global_max_peak_mem, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_peak_mem, &global_sum_peak_mem, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    
    if (my_rank == 0) {
        std::cout << "\n======================================================\n";
        std::cout << "             PERFORMANCE REPORT (Rank 0)              \n";
        std::cout << "======================================================\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Total Wall Time: " << total_elapsed << " s\n\n";

        auto print_metric = [&](const string& name, const PerfTimer& t) {
            double percent = (total_elapsed > 0) ? (t.total_time / total_elapsed * 100.0) : 0.0;
            std::cout << std::left << std::setw(25) << name 
                      << " | Calls: " << std::setw(5) << t.call_count 
                      << " | Time: " << std::setw(8) << t.total_time << " s"
                      << " | " << std::setw(5) << percent << " %\n";
        };

        print_metric("Mesh Read & Distribute", t_mesh_read);
        print_metric("Init Setup & Bdr", t_init_setup);
        print_metric("Poisson Solver", t_poisson);
        print_metric("Gradient (GetNodeE)", t_grad_calc);
        print_metric("Upwind Search", t_upwind_find);
        print_metric("Takuma Solver", t_takuma);
        print_metric("Loop Overhead/Misc", t_loop_overhead);
        print_metric("Final Output", t_output);

        std::cout << "------------------------------------------------------\n";
        std::cout << "Max Process Peak Memory (VmHWM): " << std::setw(8) << global_max_peak_mem << " MB\n";
        std::cout << "Total Machine Peak Memory:       " << std::setw(8) << global_sum_peak_mem << " MB\n";
        std::cout << "======================================================\n";
        
    }

    Mpi::Finalize();
    return 0;
}