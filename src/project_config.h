#pragma once
#include <string>
#include <vector>

struct BoundarySetup {
    int tag;
    std::string name;
    bool useFunction = false;
    double voltage = 0.0;
    std::string expression;
    bool is_corona = false;
};

// 物理场全局设置
struct PhysicsSetup {
    double E_onset = 600000.0;
    double rho_surface = 10000.0;
    double K_mobility = 1.0;
    double wind_x = 0.0, wind_y = 0.0, wind_z = 0.0;
};

// 求解器与并行控制设置
struct SolverSetup {
    std::string solver_type = "CPU";  // "CPU"=Takuma电位, "TABATA_CPU"=Tabata CPU, "GPU"=Tabata GPU, "DOUBLESTACK"=Takuma双层堆叠
    double w_rho = 1.0;
    double goal_convergence = 0.95;
    int max_update_times = 100;
    double tolerance_E = 0.01;
    double tolerance_rho = 0.01;

    int order = 1;
    int num_cores = 4;
};

// 参数系统
struct ParamEntry {
    std::string name;
    std::string expression;
    double value = 0.0;
    std::string unit;
    std::string description;
};

struct ParameterPage {
    std::string name;
    std::vector<ParamEntry> entries;
};

// 案例库
struct CaseInfo {
    std::string name;
    std::string category;
    std::string description;
    std::string folderName;
};

// 全局配置汇总
struct ProjectConfig {
    PhysicsSetup physics;
    SolverSetup solver;
    std::vector<BoundarySetup> boundaries;
    std::vector<ParameterPage> paramPages;
};