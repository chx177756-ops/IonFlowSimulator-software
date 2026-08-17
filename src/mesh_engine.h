#pragma once
#include <string>
#include <map>
#include <set>
#include <vector>
#include "mesh_data.h"

// 用于封装从外部导入的网格物理组信息
struct ImportedGroup {
    int tag;
    std::string name;
    int dimension;
};

class MeshEngine {
public:
    // 触发耗时的网格生成，立即落盘并释放内存，仅返回表面渲染数据 (Rendering Proxy)
    static FaceMeshMap generateMesh(
        const std::string& stepFilePath,
        const std::map<std::string, std::pair<std::set<int>, MeshGroupParams>>& physicalGroupData,
        int firstFaceId,
        int order,
        bool optimize,
        int numThreads,
        const std::string& cacheMshPath, // 持久化缓存路径
        std::string& outLog
    );

    // 导入外部 MSH 文件，提取物理组结构与渲染表面
    static FaceMeshMap importMesh(
        const std::string& mshFilePath, 
        std::vector<ImportedGroup>& outGroups
    );

    // 从已打开的Gmsh模型中提取表面网格
    static FaceMeshMap extractSurfaceMeshForOpenGL();
    // 二进制缓存: 保存/加载 FaceMeshMap
    static void saveMeshCache(const std::string& path, const FaceMeshMap& data);
    static FaceMeshMap loadMeshCache(const std::string& path);

    // 清理 Gmsh 引擎状态
    static void cleanup();
    // 确保 Gmsh 已初始化
    static void ensureInit();
    // 从当前 Gmsh 模型提取物理组→网格映射, 保存为 .pgroups.json
    static void savePhysicalGroupsJson(const std::string& mshPath);

private:
    static bool s_isGmshInit;
};