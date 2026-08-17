#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <mpi.h>
#include <algorithm>
#include <sys/stat.h>
#include <iomanip>
#include <fstream>
#include <string>

#include "mfem.hpp"
#include "nlohmann/json.hpp"

#include "ComputeElementVolume.hpp"
#include "GetNodeE.h"
#include "poisson_FEM_MFEM.h"
#include "ParallelUpwindElementFinder.h"
#include "CVPDE_TABATA.h"
#include "param_parser.h"

// --- 获取当前进程历史最高的物理内存占用 (VmHWM) 单位：MB ---
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
using json = nlohmann::json;

// --- 简单的计时器结构体 ---
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
    // ---------------- MPI/Hypre/GPU 初始化 ----------------
    Mpi::Init(argc, argv);
    cudaSetDevice(0);
    Hypre::Init();
    // MFEM Device = CPU: Poisson 求解器用 CPU(确保与纯CPU版结果一致)
    // Hypre = DEVICE: Tabata 求解器用 GPU 加速
    Device device("cpu");
    mfem::Hypre::configure_runtime_policy_from_mfem = false;
    HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
    HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
    HYPRE_DeviceInitialize();
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank = 0, num_procs = 1;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    // 全局计时开始
    double program_start_time = MPI_Wtime();

    // 定义各环节计时器
    PerfTimer t_mesh_read;      
    PerfTimer t_init_setup;     
    PerfTimer t_poisson;        
    PerfTimer t_grad_calc;      
    PerfTimer t_upwind_find;    
    PerfTimer t_tabata;         
    PerfTimer t_output;         
    PerfTimer t_loop_overhead;  

    // ---------------- 参数 ----------------
    double rela_error_rho = 1.0, max_rela_error_rho = 1.0, rela_error_E = 1.0;
    int count_iteration = 0;
    
    std::vector<double> rela_error_rho_history;
    std::vector<double> rela_error_E_history;
    std::vector<double> max_abs_error_rho_history;
    std::vector<double> rate_convergence_history;
    std::vector<double> rate_update_history;
    std::vector<int> iterations;

    double rate_convergence = 0.0;

    // ==========================================================
    // JSON 动态加载
    // ==========================================================
    std::string config_file = "config.json";
    if (argc > 2 && std::string(argv[1]) == "-c") config_file = argv[2];
    else if (argc == 2) config_file = argv[1];

    json config;
    try {
        std::ifstream f(config_file);
        if (!f.is_open()) { if (my_rank == 0) std::cerr << "[ERROR] 无法打开配置文件" << std::endl; MPI_Abort(comm, 1); }
        config = json::parse(f);
    } catch (...) { if (my_rank == 0) std::cerr << "[ERROR] JSON 解析失败" << std::endl; MPI_Abort(comm, 1); }

    string mesh_file = config.value("mesh_file", "mesh.msh");
    string output_folder = config.value("output_folder", "results");
    int order = config.value("order", 1);
    bool export_txt = config.value("export_txt", true);

    double E0 = config["physics"]["E_onset"];
    double rho_surface = config["physics"]["rho_surface"];
    double K_mobility = config["physics"]["K_mobility"];
    double V0 = 0;
    std::vector<double> w_arr = config["physics"]["wind_velocity"];
    Vector w_vec(3); w_vec(0)=w_arr[0]; w_vec(1)=w_arr[1]; w_vec(2)=(w_arr.size()>2)?w_arr[2]:0.0;

    double w_rho = config["solver"].value("w_rho", 0.25);
    double Goal_convergence = config["solver"]["goal_convergence"];
    int set_updatetimes = config["solver"]["max_update_times"];
    double tolerance_E = config["solver"]["tolerance_E"];
    double tolerance_rho = config["solver"]["tolerance_rho"];

    // 动态边界参数
    std::map<int, double> bdr_voltages;
    std::map<int, std::string> bdr_expressions;  // 函数表达式(可选, 与 v6 接口一致)
    std::vector<int> corona_tags, all_dirichlet_tags, inner_tags, outer_tags, artificial_tags;
    for (const auto& bdr : config["boundaries"]) {
        int tag = bdr["tag"]; double voltage = bdr["voltage"];
        bool is_corona = bdr["is_corona"];
        bool useFunc = bdr.value("use_function", false);
        std::string expr = bdr.value("voltage_expression", "");
        bdr_voltages[tag] = voltage;
        if (useFunc && !expr.empty()) bdr_expressions[tag] = expr;
        all_dirichlet_tags.push_back(tag);
        if (is_corona) corona_tags.push_back(tag);
        if (voltage >= 1.0) { inner_tags.push_back(tag); V0 = voltage; }
        else if (voltage <= 0.0) outer_tags.push_back(tag);
        else artificial_tags.push_back(tag);
    }

    if (my_rank == 0) {
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Dynamic Config from: " << config_file << std::endl;
        std::cout << "E_onset: " << E0 << " V/m, V0: " << V0 << std::endl;
        std::cout << "Boundary Conditions Applied: " << bdr_voltages.size() << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        std::cout << "Reading mesh..." << std::endl;
    }

    t_mesh_read.Start();

    // 检测预分区并行网格文件: <mesh_file>.part.000000
    std::ostringstream part0_name;
    part0_name << mesh_file << ".part." << std::setfill('0') << std::setw(6) << 0;
    std::ifstream test_part(part0_name.str());
    bool use_parallel_mesh = test_part.good() && num_procs > 1;
    test_part.close();

    // 【新增】验证分区数与 MPI 进程数一致, 不一致则回退串行加载
    // (分区文件声明了固定分区拓扑, 若与启动的进程数不符, ParMesh 构造时
    //  会因 rank 数不匹配触发 MPI_ERR_RANK 崩溃)
    if (use_parallel_mesh) {
        int num_parts = 0;
        for (int i = 0; i < 4096; i++) {
            std::ostringstream oss;
            oss << mesh_file << ".part." << std::setfill('0') << std::setw(6) << i;
            std::ifstream f(oss.str());
            if (f.good()) num_parts = i + 1;
            else break;
        }
        if (num_parts != num_procs) {
            if (my_rank == 0)
                std::cout << "[Mesh] 分区数 (" << num_parts << ") != MPI 进程数 ("
                          << num_procs << "), 回退串行加载..." << std::endl;
            use_parallel_mesh = false;
        }
    }

    Mesh serial_mesh;
    ParMesh *pmesh = nullptr;

    if (use_parallel_mesh) {
        if (my_rank == 0)
            std::cout << "[Mesh] Loading pre-partitioned parallel mesh..." << std::endl;
        std::ostringstream my_name;
        my_name << mesh_file << ".part." << std::setfill('0') << std::setw(6) << my_rank;
        std::ifstream ifs(my_name.str());
        pmesh = new ParMesh(comm, ifs, /*refine=*/false, /*generate_edges=*/1, /*fix_orientation=*/true);
    } else {
        serial_mesh = Mesh(mesh_file.c_str(), 1, 1);
        pmesh = new ParMesh(comm, serial_mesh);
    }
    ParMesh &mesh = *pmesh;

    if (my_rank != 0) serial_mesh.Clear();
    t_mesh_read.Stop();

    // [bench] MESH_BENCH=1 时仅计时网格加载后退出
    if (std::getenv("MESH_BENCH")) {
        double tmax;
        MPI_Reduce(&t_mesh_read.total_time, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        if (my_rank == 0) printf("MESH_LOAD_TIME: %.4f\n", tmax);
        MPI_Finalize();
        return 0;
    }

    // .pgroups.json 加载
    std::map<int, std::vector<int>> pg_nodes;
    {
        std::string pg_file = mesh_file + ".pgroups.json";
        std::ifstream pgf(pg_file);
        if (pgf.is_open()) {
            json pgj = json::parse(pgf);
            for (auto& g : pgj["groups"]) {
                int tag = g["tag"];
                for (auto& n : g["nodes"])
                    pg_nodes[tag].push_back(n.get<int>() - 1);
                if (my_rank == 0) std::cout << "[PG] group '" << g["name"] << "' tag=" << tag
                    << " nodes=" << g["nodes"].size() << std::endl;
            }
            if (my_rank == 0) std::cout << "[PG] Loaded " << pg_nodes.size() << " groups" << std::endl;
        }
    }

    // === 计时：初始化设置 ===
    t_init_setup.Start();

    H1_FECollection scalar_fec(order, mesh.Dimension());
    ParFiniteElementSpace scalar_fes(&mesh, &scalar_fec); 
    H1_FECollection vector_fec(order, mesh.Dimension());
    ParFiniteElementSpace vector_fes(&mesh, &vector_fec, 3); 

    std::vector<double> Volume = ComputeElementVolume(mesh);

    // [优化] 预计算所有单元的形函数梯度 (P1单元为常数, 全程序只算1次)
    // UpwindElementFinder 和 ExchangeUpwindData 均复用此缓存
    int numTet = mesh.GetNE();
    int dim = mesh.Dimension();
    std::vector<mfem::DenseMatrix> global_dshapes(numTet);
    {
        mfem::DenseMatrix dshape_ref(4, dim);
        mfem::DenseMatrix Jinv(dim);
        for (int e = 0; e < numTet; ++e) {
            global_dshapes[e].SetSize(4, dim);
            mfem::ElementTransformation *T = mesh.GetElementTransformation(e);
            const mfem::IntegrationRule &ir = mfem::IntRules.Get(mesh.GetElementType(e), 2);
            T->SetIntPoint(&ir.IntPoint(0));
            CalcInverse(T->Jacobian(), Jinv);
            scalar_fes.GetFE(e)->CalcDShape(ir.IntPoint(0), dshape_ref);
            Mult(dshape_ref, Jinv, global_dshapes[e]);
        }
    }
    if (my_rank == 0) std::cout << "[优化] 形函数梯度预计算完成 (" << numTet << " 个单元)" << std::endl;

    // 动态边界提取 (优先 .pgroups.json)
    auto collect_boundary_dofs = [&](int target_attr) -> std::vector<int> {
        auto it = pg_nodes.find(target_attr);
        if (it != pg_nodes.end() && !it->second.empty()) {
            std::set<int> s;
            for (int vi : it->second)
                if (vi >= 0 && vi < mesh.GetNV()) {
                    Array<int> vdofs; scalar_fes.GetVertexDofs(vi, vdofs);
                    for (int j = 0; j < vdofs.Size(); j++) s.insert(vdofs[j]);
                }
            return std::vector<int>(s.begin(), s.end());
        }
        std::set<int> s;
        for (int i = 0; i < mesh.GetNBE(); ++i) {
            Element *bel = mesh.GetBdrElement(i);
            if (!bel) continue;
            if (bel->GetAttribute() == target_attr) {
                Array<int> dofs; scalar_fes.GetBdrElementDofs(i, dofs);
                for (int j = 0; j < dofs.Size(); ++j) s.insert(dofs[j]);
            }
        }
        return std::vector<int>(s.begin(), s.end());
    };

    // 合并 inner/outer 顶点
    std::vector<int> vertex_inner;
    for (int tag : inner_tags) { auto d = collect_boundary_dofs(tag); vertex_inner.insert(vertex_inner.end(), d.begin(), d.end()); }
    std::sort(vertex_inner.begin(), vertex_inner.end());
    vertex_inner.erase(std::unique(vertex_inner.begin(), vertex_inner.end()), vertex_inner.end());

    std::vector<int> vertex_outer;
    for (int tag : outer_tags) { auto d = collect_boundary_dofs(tag); vertex_outer.insert(vertex_outer.end(), d.begin(), d.end()); }
    std::sort(vertex_outer.begin(), vertex_outer.end());
    vertex_outer.erase(std::unique(vertex_outer.begin(), vertex_outer.end()), vertex_outer.end());

    // 起晕候选节点: 仅来自 is_corona=true 的边界
    std::vector<int> vertex_corona_candidates;
    for (int tag : corona_tags) { auto d = collect_boundary_dofs(tag); vertex_corona_candidates.insert(vertex_corona_candidates.end(), d.begin(), d.end()); }
    std::sort(vertex_corona_candidates.begin(), vertex_corona_candidates.end());
    vertex_corona_candidates.erase(std::unique(vertex_corona_candidates.begin(), vertex_corona_candidates.end()), vertex_corona_candidates.end());

    // 初始 rho
    int numLocalDofs = scalar_fes.GetNDofs();
    Vector rho_local(numLocalDofs);
    rho_local = 0.0;

    PoissonSolver poisson_solver(mesh, scalar_fes, bdr_voltages, bdr_expressions);

    t_init_setup.Stop(); 

    // 初始 Poisson (此时rho_local=0，解出来的是拉普拉斯场，用于寻找起晕点)
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

    // 判定起晕点
    std::vector<double> nodeEn_initial(numLocalDofs, 0.0);
    for (int ld = 0; ld < numLocalDofs; ++ld) nodeEn_initial[ld] = (ld < nodeEn.Size()) ? nodeEn(ld) : 0.0;

    std::vector<char> is_corona_local(numLocalDofs, 0);
    for (int idx : vertex_corona_candidates) {
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
        // 赋予初始空间电荷
        for (int i = 0; i < numLocalDofs; ++i) rho_local(i) = rho_surface / 100.0;
        for (int idx : vertex_inner) if (idx >=0 && idx < numLocalDofs) rho_local(idx) = 1e-12;
        for (int idx : vertex_corona) if (idx >=0 && idx < numLocalDofs) rho_local(idx) = rho_surface;
    }
    t_init_setup.Stop(); // 初始化结束

    // ---------------- 迭代循环 ----------------
    int count_update = 0;

    ParGridFunction v_total(&vector_fes);
    std::vector<UpwindInfo> Upwindelements;

    if (run_main_loop) 
    {
        if (my_rank == 0) std::cout << "Start iterative loop (TABATA Method)..." << std::endl;
        while (rela_error_E >= tolerance_E)
        {
            ++count_iteration;
            iterations.push_back(count_iteration);
            if (my_rank == 0) std::cout << "Iteration: " << count_iteration << std::endl;

            Vector rho_prev = rho_local;

            // === 计时：Poisson 求解 ===
            // (此时的 rho_local 包含了你最初赋予的大量电荷)
            t_poisson.Start();
            ParGridFunction phi(&scalar_fes);
            poisson_solver.Solve(rho_local, phi);
            t_poisson.Stop();

            // === 计时：电场与速度更新 ===
            t_grad_calc.Start();
            for (int e = 0; e < mesh.GetNE(); ++e) {
                ElementTransformation *tr = mesh.GetElementTransformation(e);
                Vector grad;
                phi.GetGradient(*tr, grad);
                E_elem(e,0) = -grad(0); E_elem(e,1) = -grad(1); E_elem(e,2) = -grad(2);
            }
            GetNodeE(E_elem, mesh, Volume, vector_fes, nodeE, nodeEn);

            v_total = nodeE;       
            v_total *= K_mobility; 
            VectorConstantCoefficient w_coeff(w_vec);
            ParGridFunction w_gf(&vector_fes);
            w_gf.ProjectCoefficient(w_coeff);
            v_total += w_gf;       
            t_grad_calc.Stop();

            // === 计时：上流元搜索 ===
            t_upwind_find.Start();
                UpwindElementFinder finder(mesh, v_total, vector_fes, global_dshapes);
                Upwindelements = finder.ComputeUpwindElements(); 
                t_upwind_find.Stop();

            // === 计时：Tabata 求解 ===
            t_tabata.Start();
            Vector rho_calc_local = rho_local; 
            
            SolveContinuityTabata(mesh, scalar_fes, vector_fes, Volume, rho_local, nodeE, Upwindelements, rho_calc_local, vertex_corona, global_dshapes);

            t_tabata.Stop();

            // === 计时：循环杂项 ===
            t_loop_overhead.Start();
            
            // 松弛更新
            for(int i=0; i<numLocalDofs; ++i) {
                rho_local(i) = w_rho * rho_calc_local(i) + (1.0 - w_rho) * rho_prev(i);
            }
            scalar_fes.GroupComm().Bcast<double>(rho_local.GetData());

            // 误差计算
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

            // 计算收敛率
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

            // Corona 更新 (检测激波事件)
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

    // ============================
    // 输出部分
    // ============================
    t_output.Start();
    {
        int myid = my_rank;
        int nprocs = num_procs;

        ParGridFunction phi_final(&scalar_fes);
        poisson_solver.Solve(rho_local, phi_final);

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

        // VTK/ParaView 输出 (后处理依赖此格式)
        {
            ParGridFunction rho_gf(&scalar_fes); rho_gf = rho_local;
            ParGridFunction nodeE0_final(&scalar_fes);
            for (int i = 0; i < scalar_fes.GetVSize(); ++i)
                nodeE0_final(i) = (i < (int)nodeEn_initial.size()) ? nodeEn_initial[i] : 0.0;

            if (my_rank == 0) {
                #ifdef _WIN32
                    _mkdir(output_folder.c_str());
                #else
                    mkdir(output_folder.c_str(), 0777);
                #endif
            }
            MPI_Barrier(comm);

            ParaViewDataCollection paraview_dc("Takuma_Results", &mesh);
            paraview_dc.SetPrefixPath(output_folder);
            paraview_dc.SetLevelsOfDetail(order);
            paraview_dc.SetHighOrderOutput(false);
            paraview_dc.RegisterField("rho", &rho_gf);
            paraview_dc.RegisterField("E_scalar", &nodeEn_final);
            paraview_dc.RegisterField("E_vector", &nodeE_final);
            paraview_dc.RegisterField("phi", &phi_final);
            paraview_dc.RegisterField("E0_initial", &nodeE0_final);
            paraview_dc.Save();

            if (my_rank == 0)
                std::cout << "[Output] VTK results saved to " << output_folder << "/Takuma_Results/" << std::endl;
        }

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

        int my_localNV = localNV;
        std::vector<double> send_coords(3 * my_localNV);
        std::vector<double> send_rho(my_localNV);
        std::vector<double> send_E(my_localNV);
        std::vector<double> send_E0(my_localNV);      

        for (int i = 0; i < my_localNV; ++i) {
            const double *c = mesh.GetVertex(i);
            send_coords[3*i + 0] = c[0]; send_coords[3*i + 1] = c[1]; send_coords[3*i + 2] = c[2];

            int dof = vertex_to_dof[i];
            if (dof >= 0 && dof < rho_local.Size()) send_rho[i] = rho_local(dof); else send_rho[i] = 0.0;
            if (dof >= 0 && dof < nodeEn_final.Size()) send_E[i] = nodeEn_final(dof); else send_E[i] = 0.0;
            if (dof >= 0 && dof < (int)nodeEn_initial.size()) send_E0[i] = nodeEn_initial[dof]; else send_E0[i] = 0.0;
        }

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

        if (myid == 0)
        {
            std::cout << "[Output] Gathering data on Rank 0..." << std::endl;

            struct PointData {
                int id;
                double x, y, z;
                double rho, E, E0;
            };

            std::vector<PointData> parallel_points(total_entries);
            for(int i = 0; i < total_entries; ++i) {
                parallel_points[i].id = -1;
                parallel_points[i].x = recv_coords[3*i+0];
                parallel_points[i].y = recv_coords[3*i+1];
                parallel_points[i].z = recv_coords[3*i+2];
                parallel_points[i].rho = recv_rho[i];
                parallel_points[i].E = recv_E[i];
                parallel_points[i].E0 = recv_E0[i];
            }

            auto comp_strict = [](const PointData& a, const PointData& b) {
                if (std::isnan(a.x) || std::isnan(b.x)) return false;
                if (std::abs(a.x - b.x) > 1e-14) return a.x < b.x;
                if (std::abs(a.y - b.y) > 1e-14) return a.y < b.y;
                return a.z < b.z;
            };

            const double match_tol_sq = 1e-12;

            std::vector<PointData> serial_points;
            if (serial_mesh.GetNV() > 0) {
                // 串行加载: 用 serial_mesh 顶点做匹配
                int NV_serial = serial_mesh.GetNV();
                serial_points.resize(NV_serial);
                for(int i = 0; i < NV_serial; ++i) {
                    const double* v = serial_mesh.GetVertex(i);
                    serial_points[i].id = i;
                    serial_points[i].x = v[0];
                    serial_points[i].y = v[1];
                    serial_points[i].z = v[2];
                    serial_points[i].rho = 0.0;
                }
            } else {
                // 并行加载: 从预保存的顶点坐标文件恢复串行顺序
                std::string vert_file = mesh_file + ".part.vertices";
                std::ifstream vf(vert_file);
                if (vf.is_open()) {
                    int nv, sd;
                    vf >> nv >> sd;
                    serial_points.resize(nv);
                    for (int i = 0; i < nv; i++) {
                        double vx, vy, vz = 0.0;
                        vf >> vx >> vy;
                        if (sd >= 3) vf >> vz;
                        serial_points[i].id = i;
                        serial_points[i].x = vx;
                        serial_points[i].y = vy;
                        serial_points[i].z = vz;
                        serial_points[i].rho = 0.0;
                    }
                    vf.close();
                    std::cout << "[Output] Loaded " << nv << " vertex coords from " << vert_file << std::endl;
                } else {
                    // 回退: 从收集数据中去重
                    std::cout << "[Output] Building global vertex list from gathered data..." << std::endl;
                    std::vector<PointData> tmp(parallel_points);
                    std::sort(tmp.begin(), tmp.end(), comp_strict);
                    for (size_t i = 0; i < tmp.size(); ++i) {
                        if (i > 0) {
                            double dx = tmp[i].x - tmp[i-1].x;
                            double dy = tmp[i].y - tmp[i-1].y;
                            double dz = tmp[i].z - tmp[i-1].z;
                            if (dx*dx + dy*dy + dz*dz < match_tol_sq) continue;
                        }
                        PointData p = tmp[i];
                        p.id = (int)serial_points.size();
                        p.rho = 0.0; p.E = 0.0; p.E0 = 0.0;
                        serial_points.push_back(p);
                    }
                    std::cout << "[Output] Deduplicated " << tmp.size() << " -> "
                              << serial_points.size() << " unique vertices" << std::endl;
                }
            }

            int NV_serial = (int)serial_points.size();

            std::cout << "[Output] Sorting..." << std::endl;
            std::sort(parallel_points.begin(), parallel_points.end(), comp_strict);
            std::sort(serial_points.begin(), serial_points.end(), comp_strict);

            const double search_tol = 1e-5;
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
                        matched_count++;
                        break; 
                    }
                }
            }
            std::cout << "[Output] Matched " << matched_count << " / " << NV_serial << " vertices." << std::endl;

            std::vector<double> final_rho(NV_serial, 0.0);
            std::vector<double> final_E(NV_serial, 0.0);
            std::vector<double> final_E0(NV_serial, 0.0);
            
            for(const auto& p : serial_points) {
                if (p.id >= 0 && p.id < NV_serial) {
                    final_rho[p.id] = p.rho;
                    final_E[p.id]   = p.E;
                    final_E0[p.id]  = p.E0;
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
            
            cout << "Rank 0 wrote results successfully." << endl;
        }
    }
    t_output.Stop();

    // ============================
    // 性能报告输出 (仅 Rank 0)
    // ============================
    double program_end_time = MPI_Wtime();
    double total_elapsed = program_end_time - program_start_time;
    
    double local_peak_mem = GetPeakMemoryMB();
    double global_max_peak_mem = 0.0;
    double global_sum_peak_mem = 0.0;
    
    MPI_Reduce(&local_peak_mem, &global_max_peak_mem, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_peak_mem, &global_sum_peak_mem, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    

    if (my_rank == 0) {
        std::cout << "\n======================================================\n";
        std::cout << "              PERFORMANCE REPORT (Rank 0)             \n";
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
        print_metric("Tabata Solver", t_tabata);
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