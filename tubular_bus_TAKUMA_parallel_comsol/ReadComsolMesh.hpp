#ifndef READ_COMSOL_MESH_HPP
#define READ_COMSOL_MESH_HPP

#include "mfem.hpp"
#include <string>
#include <map>
#include <vector>

// 解析 COMSOL .mphtxt 并返回 MFEM 网格，同时提取物理组的选择标签 (Selection Labels)
mfem::Mesh* ReadComsolMesh(const std::string& filename, 
                           std::map<std::string, std::vector<int>>& selection_dict);

#endif // READ_COMSOL_MESH_HPP