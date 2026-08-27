#include "ReadComsolMesh.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

using namespace mfem;

Mesh* ReadComsolMesh(const std::string& filename, std::map<std::string, std::vector<int>>& selection_dict) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        MFEM_ABORT("Failed to open COMSOL mesh file: " + filename);
    }

    std::string line;
    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<int>> tet_nodes;
    std::vector<int> tet_attrs;
    std::vector<std::vector<int>> tri_nodes;
    std::vector<int> tri_attrs;

    int current_nodes_per_elem = 0;
    
    // 用于 Selection 块解析的临时变量
    std::string current_label = "";
    int num_entities_to_read = 0;
    std::vector<int> current_entities;

    // 引入状态机
    enum State { NONE, READING_VERTICES, READING_ELEMENTS, READING_ATTRS, PARSING_SELECTION, READING_SELECTION_ENTITIES };
    State state = NONE;

    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n")); // 清除行首空格
        if (line.empty()) continue;

        // === 1. 侦测关键字，切换流水线状态 ===
        if (line.find("Mesh vertex coordinates") != std::string::npos) {
            state = READING_VERTICES; continue;
        }
        if (line.find("number of nodes per element") != std::string::npos ||
            line.find("number of vertices per element") != std::string::npos) {
            std::istringstream iss(line);
            iss >> current_nodes_per_elem; continue;
        }
        if (line.find("# Elements") != std::string::npos && line.find("type") == std::string::npos) {
            state = READING_ELEMENTS; continue;
        }
        if (line.find("Geometric entity indices") != std::string::npos) {
            state = READING_ATTRS; continue;
        }
        // [新增] 侦测到物理组 Selection 块
        if (line.find("Selection # class") != std::string::npos) {
            state = PARSING_SELECTION;
            current_label = "";
            num_entities_to_read = 0;
            current_entities.clear();
            continue;
        }

        // 遇到注释或分隔符中止普通读取状态，但 Selection 中的注释带有关键信息，需特别处理
        if ((line[0] == '#' || line.find("---") != std::string::npos) && 
             state != PARSING_SELECTION && state != READING_SELECTION_ENTITIES) {
            state = NONE; continue;
        }

        // === 2. 按状态提取数据 ===
        if (state == READING_VERTICES) {
            std::istringstream iss(line);
            double x, y, z;
            if (iss >> x >> y >> z) vertices.push_back({x, y, z});
        }
        else if (state == READING_ELEMENTS) {
            std::istringstream iss(line);
            std::vector<int> nodes(current_nodes_per_elem);
            bool valid = true;
            for (int i = 0; i < current_nodes_per_elem; ++i) {
                if (!(iss >> nodes[i])) { valid = false; break; }
            }
            if (valid) {
                if (current_nodes_per_elem == 4) tet_nodes.push_back(nodes);
                else if (current_nodes_per_elem == 3) tri_nodes.push_back(nodes);
            }
        }
        else if (state == READING_ATTRS) {
            std::istringstream iss(line);
            int attr;
            if (iss >> attr) {
                attr = attr + 1; // [MFEM 必须为正整数]
                if (current_nodes_per_elem == 4) tet_attrs.push_back(attr);
                else if (current_nodes_per_elem == 3) tri_attrs.push_back(attr);
            }
        }
        // [新增] 解析 Selection 头部信息
        else if (state == PARSING_SELECTION) {
            if (line.find("# Label") != std::string::npos) {
                // 提取标签字符串 (例如从 "4 ring # Label" 提取出 "ring")
                size_t first_space = line.find(' ');
                size_t hash_pos = line.find(" # Label");
                if (first_space != std::string::npos && hash_pos != std::string::npos && hash_pos > first_space) {
                    current_label = line.substr(first_space + 1, hash_pos - first_space - 1);
                }
            }
            else if (line.find("# Number of entities") != std::string::npos) {
                std::istringstream iss(line);
                iss >> num_entities_to_read;
            }
            else if (line.find("# Entities") != std::string::npos) {
                if (num_entities_to_read > 0) state = READING_SELECTION_ENTITIES;
                else state = NONE;
            }
        }
        // [新增] 提取 Selection 具体编号
        else if (state == READING_SELECTION_ENTITIES) {
            std::istringstream iss(line);
            int ent;
            while (iss >> ent) {
                // ⚠️ 极其关键：必须 +1 才能与前面网格属性 +1 的偏移对齐！
                current_entities.push_back(ent + 1); 
                num_entities_to_read--;
            }
            if (num_entities_to_read <= 0) {
                if (!current_label.empty()) {
                    selection_dict[current_label] = current_entities;
                }
                state = NONE;
            }
        }
    }
    file.close();

    if (vertices.empty() || tet_nodes.empty()) {
        MFEM_ABORT("COMSOL Parser Error: Failed to extract vertices or tetrahedra.");
    }

    Mesh *mesh = new Mesh(3, vertices.size(), tet_nodes.size(), tri_nodes.size(), 3);
    for (const auto& v : vertices) mesh->AddVertex(v[0], v[1], v[2]);
    for (size_t i = 0; i < tet_nodes.size(); ++i) mesh->AddTet(tet_nodes[i].data(), tet_attrs[i]);
    for (size_t i = 0; i < tri_nodes.size(); ++i) mesh->AddBdrTriangle(tri_nodes[i].data(), tri_attrs[i]);

    mesh->FinalizeTopology();
    mesh->SetAttributes(); 
    return mesh;
}