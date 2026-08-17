#include "mesh_engine.h"
#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cstdio>
#include <thread>
#include <algorithm>
#include <unordered_map>

bool MeshEngine::s_isGmshInit = false;

void MeshEngine::cleanup() {
    if (s_isGmshInit) {
        try { gmsh::finalize(); } catch (...) {}
        s_isGmshInit = false;
    }
}

void MeshEngine::ensureInit() {
    if (!s_isGmshInit) { gmsh::initialize(); s_isGmshInit = true; }
}

void MeshEngine::saveMeshCache(const std::string& path, const FaceMeshMap& data) {
    std::ofstream f(path, std::ios::binary);
    int32_t n = data.size(); f.write((char*)&n, 4);
    for (auto& [faceId, rm] : data) {
        f.write((char*)&faceId, 4);
        int32_t vc = rm.vertices.size(), ic = rm.indices.size();
        f.write((char*)&vc, 4); f.write((char*)&ic, 4);
        f.write((char*)rm.vertices.data(), vc * 4);
        f.write((char*)rm.indices.data(), ic * 4);
    }
}

FaceMeshMap MeshEngine::loadMeshCache(const std::string& path) {
    FaceMeshMap result;
    std::ifstream f(path, std::ios::binary);
    if (!f) return result;
    int32_t n; f.read((char*)&n, 4);
    for (int32_t i = 0; i < n; i++) {
        int32_t faceId, vc, ic;
        f.read((char*)&faceId, 4); f.read((char*)&vc, 4); f.read((char*)&ic, 4);
        RenderableMesh rm;
        rm.vertices.resize(vc); rm.indices.resize(ic);
        f.read((char*)rm.vertices.data(), vc * 4);
        f.read((char*)rm.indices.data(), ic * 4);
        result[faceId] = rm;
    }
    return result;
}

// ====================================================================
// 【架构核心】：降维提取，仅返回用于 OpenGL 渲染的表面网格
// ====================================================================
FaceMeshMap MeshEngine::extractSurfaceMeshForOpenGL() {
    FaceMeshMap faceMeshes;

    // 获取全局节点 (所有维度的所有节点, 避免面边界处遗漏)
    std::vector<std::size_t> globalNodeTags;
    std::vector<double> globalCoords, globalParamCoord;
    try {
        gmsh::model::mesh::getNodes(globalNodeTags, globalCoords, globalParamCoord);
    } catch (...) { return faceMeshes; }
    if (globalNodeTags.empty()) return faceMeshes;

    std::unordered_map<std::size_t, unsigned int> globalTagToIdx;
    for (size_t i = 0; i < globalNodeTags.size(); i++)
        globalTagToIdx[globalNodeTags[i]] = static_cast<unsigned int>(i);

    // 获取所有 2D 单元 (不限定几何面, 避免面边界处单元丢失)
    std::vector<int> elemTypes;
    std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    try {
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, 2);
    } catch (...) { return faceMeshes; }

    // 获取几何面 → 物理组映射: 每个 2D 几何面属于哪个物理组
    std::map<int, std::vector<int>> faceToGroups; // faceId → [groupTag, ...]
    try {
        std::vector<std::pair<int,int>> physGroups;
        gmsh::model::getPhysicalGroups(physGroups);
        for (auto& pg : physGroups) {
            if (pg.first != 2) continue;
            std::vector<int> entityTags;
            gmsh::model::getEntitiesForPhysicalGroup(pg.first, pg.second, entityTags);
            for (int et : entityTags)
                faceToGroups[et].push_back(pg.second);
        }
    } catch (...) {}

    // 获取单元 → 几何面映射: 每个 2D 单元属于哪个几何面
    // Gmsh: getElements 返回的 elemTags 是单元 tag, 需用 getElements 的 dim+tag 参数...
    // 简化: 对所有 2D 几何面分别获取单元, 与全局列表对照
    std::vector<std::pair<int,int>> dimTags;
    try { gmsh::model::getEntities(dimTags, 2); } catch (...) {}

    // 全局 2D 单元按几何面分组
    std::map<int, std::vector<unsigned int>> faceToElemIndices;
    std::map<int, std::set<std::size_t>> faceToNodeTags;

    for (auto& dt : dimTags) {
        int faceId = dt.second;
        std::vector<int> fElemTypes;
        std::vector<std::vector<std::size_t>> fElemTags, fNodeTags;
        try {
            gmsh::model::mesh::getElements(fElemTypes, fElemTags, fNodeTags, 2, faceId);
        } catch (...) { continue; }

        for (size_t ei = 0; ei < fElemTypes.size(); ei++) {
            int etype = fElemTypes[ei];
            if (etype != 2 && etype != 3 && etype != 9 && etype != 10 && etype != 16) continue;
            const auto& fnodes = fNodeTags[ei];
            int npt = (etype == 2) ? 3 : (etype == 9) ? 6 :
                      (etype == 3) ? 4 : (etype == 10) ? 9 : 8;
            for (size_t j = 0; j < fnodes.size(); j += npt) {
                for (int k = 0; k < std::min(npt, 4); k++)
                    faceToNodeTags[faceId].insert(fnodes[j + k]);
            }
        }
    }

    // 为每个面构建渲染数据
    for (auto& dt : dimTags) {
        int faceId = dt.second;
        // 构建此面的节点子集
        const auto& faceNodes = faceToNodeTags[faceId];
        if (faceNodes.empty()) continue;

        RenderableMesh rMesh;
        std::unordered_map<std::size_t, unsigned int> localMap;
        for (std::size_t tag : faceNodes) {
            auto it = globalTagToIdx.find(tag);
            if (it == globalTagToIdx.end()) continue;
            localMap[tag] = static_cast<unsigned int>(rMesh.vertices.size() / 3);
            rMesh.vertices.push_back(static_cast<float>(globalCoords[3 * it->second]));
            rMesh.vertices.push_back(static_cast<float>(globalCoords[3 * it->second + 1]));
            rMesh.vertices.push_back(static_cast<float>(globalCoords[3 * it->second + 2]));
        }

        // 获取此面的单元
        std::vector<int> fElemTypes;
        std::vector<std::vector<std::size_t>> fElemTags, fNodeTags;
        try {
            gmsh::model::mesh::getElements(fElemTypes, fElemTags, fNodeTags, 2, faceId);
        } catch (...) { continue; }

        for (size_t ei = 0; ei < fElemTypes.size(); ei++) {
            int etype = fElemTypes[ei];
            const auto& fnodes = fNodeTags[ei];
            if (etype == 2 || etype == 9) {
                int npt = (etype == 2) ? 3 : 6;
                for (size_t j = 0; j + 2 < fnodes.size(); j += npt) {
                    auto it0 = localMap.find(fnodes[j]);
                    auto it1 = localMap.find(fnodes[j+1]);
                    auto it2 = localMap.find(fnodes[j+2]);
                    if (it0 != localMap.end() && it1 != localMap.end() && it2 != localMap.end()) {
                        rMesh.indices.push_back(it0->second);
                        rMesh.indices.push_back(it1->second);
                        rMesh.indices.push_back(it2->second);
                    }
                }
            } else if (etype == 3 || etype == 10 || etype == 16) {
                int npt = (etype == 3) ? 4 : (etype == 10 ? 9 : 8);
                for (size_t j = 0; j + 3 < fnodes.size(); j += npt) {
                    auto it0 = localMap.find(fnodes[j]);
                    auto it1 = localMap.find(fnodes[j+1]);
                    auto it2 = localMap.find(fnodes[j+2]);
                    auto it3 = localMap.find(fnodes[j+3]);
                    if (it0 != localMap.end() && it1 != localMap.end() &&
                        it2 != localMap.end() && it3 != localMap.end()) {
                        rMesh.indices.push_back(it0->second);
                        rMesh.indices.push_back(it1->second);
                        rMesh.indices.push_back(it2->second);
                        rMesh.indices.push_back(it0->second);
                        rMesh.indices.push_back(it2->second);
                        rMesh.indices.push_back(it3->second);
                    }
                }
            }
        }
        if (!rMesh.indices.empty()) faceMeshes[faceId] = rMesh;
    }
    return faceMeshes;
}

// ====================================================================
// 【引擎生成】：解析求距、Sigmoid平滑、持久化落盘与内存核爆清理
// ====================================================================
FaceMeshMap MeshEngine::generateMesh(
    const std::string& stepFilePath,
    const std::map<std::string, std::pair<std::set<int>, MeshGroupParams>>& physicalGroupData,
    int firstFaceId, int order, bool optimize, int numThreads, const std::string& cacheMshPath, std::string& outLog)
{
    auto safeSetOption = [](const std::string& newApi, const std::string& oldApi, double val) {
        try { gmsh::option::setNumber(newApi, val); } 
        catch (...) { try { if (!oldApi.empty()) gmsh::option::setNumber(oldApi, val); } catch (...) {} }
    };
    auto safeSetFieldNumbers = [](int id, const std::string& prop, const std::string& altProp, const std::vector<double>& val) {
        try { gmsh::model::mesh::field::setNumbers(id, prop, val); } 
        catch (...) { if (!altProp.empty()) { try { gmsh::model::mesh::field::setNumbers(id, altProp, val); } catch (...) {} } }
    };
    auto safeSetFieldNumber = [](int id, const std::string& prop, const std::string& altProp, double val) {
        try { gmsh::model::mesh::field::setNumber(id, prop, val); } 
        catch (...) { if (!altProp.empty()) { try { gmsh::model::mesh::field::setNumber(id, altProp, val); } catch (...) {} } }
    };

    try {
        if (!s_isGmshInit) { gmsh::initialize(); s_isGmshInit = true; }
        gmsh::clear(); 

        try { gmsh::option::setNumber("General.Cancel", 0); } catch (...) {}

        if (numThreads > 0) {
            safeSetOption("General.NumThreads", "", numThreads);
            safeSetOption("Geometry.NumThreads", "", numThreads);
            safeSetOption("Mesh.MaxNumThreads1D", "", numThreads);
            safeSetOption("Mesh.MaxNumThreads2D", "", numThreads);
            safeSetOption("Mesh.MaxNumThreads3D", "", numThreads);
        }

        gmsh::model::add("IonFlowModel");
        std::vector<std::pair<int, int>> importedShapes; 
        fprintf(stderr, "[DEBUG mesh] importing: %s\n", stepFilePath.c_str());
        gmsh::model::occ::importShapes(stepFilePath, importedShapes);
        fprintf(stderr, "[DEBUG mesh] imported shapes=%zu\n", importedShapes.size());
        gmsh::model::occ::synchronize();
        // 碎片化: 仅多实体时创建共享拓扑, 单体跳过避免冗余
        if (importedShapes.size() > 1) {
            fprintf(stderr, "[DEBUG mesh] fragmenting %zu shapes...\n", importedShapes.size());
            gmsh::vectorpair outTags;
            std::vector<gmsh::vectorpair> outTagMap;
            gmsh::model::occ::fragment(importedShapes, {}, outTags, outTagMap, -1, true);
            fprintf(stderr, "[DEBUG mesh] after fragment\n");
        } else {
            fprintf(stderr, "[DEBUG mesh] single body, skip fragment\n");
        }
        gmsh::model::occ::synchronize();

        std::vector<std::pair<int, int>> allFaces;
        gmsh::model::getEntities(allFaces, 2);

        // DEBUG: 列出所有Gmsh面
        fprintf(stderr, "[DEBUG mesh] allFaces count=%zu: ", allFaces.size());
        for (auto& dt : allFaces) fprintf(stderr, "%d ", dt.second);
        fprintf(stderr, "\n");

        double globalMaxSize = 1.0, globalMinSize = 1000.0;
        for (const auto& [gName, dataPair] : physicalGroupData) {
            if (dataPair.second.sizeMax > globalMaxSize) globalMaxSize = dataPair.second.sizeMax;
            if (dataPair.second.sizeMin < globalMinSize) globalMinSize = dataPair.second.sizeMin;
        }
        
        safeSetOption("Mesh.MeshSizeMin", "Mesh.CharacteristicLengthMin", 1e-7);
        safeSetOption("Mesh.MeshSizeMax", "Mesh.CharacteristicLengthMax", globalMaxSize);
        fprintf(stderr, "[DEBUG mesh] globalSize min=%.4f max=%.4f\n", globalMinSize, globalMaxSize);
        safeSetOption("Mesh.MeshSizeFromPoints", "Mesh.CharacteristicLengthFromPoints", 0);
        safeSetOption("Mesh.MeshSizeFromCurvature", "Mesh.CharacteristicLengthFromCurvature", 6);
        safeSetOption("Mesh.MeshSizeExtendFromBoundary", "Mesh.CharacteristicLengthExtendFromBoundary", 0);
        safeSetOption("Mesh.Algorithm", "", 6);     // 2D Frontal-Delaunay
        safeSetOption("Mesh.Algorithm3D", "", 1);   // 3D Delaunay 算法 (HXT=10有bug)
        safeSetOption("Mesh.Smoothing", "", 1);

        int tagCounter = 1, fieldTagCounter = 1; 
        std::vector<double> thresholdFieldTags; 

        // 内部几何: 用固定首面ID计算偏移(不受group选择变化影响)
        int globalOffset = (firstFaceId > 1000) ? (firstFaceId - 1) : 0;
        fprintf(stderr, "[DEBUG mesh] firstFaceId=%d globalOffset=%d\n", firstFaceId, globalOffset);

        for (const auto& [groupName, dataPair] : physicalGroupData) {
            const std::set<int>& faceIds = dataPair.first;
            const MeshGroupParams& params = dataPair.second;

            if (faceIds.empty()) continue;

            fprintf(stderr, "[DEBUG mesh] group='%s' faceIds={", groupName.c_str());
            for (int id : faceIds) fprintf(stderr, "%d ", id);
            fprintf(stderr, "} sizeMin=%.2f sizeMax=%.2f\n", params.sizeMin, params.sizeMax);

            std::vector<int> validIds; std::vector<double> doubleFaceIds;
            for (int id : faceIds) {
                auto it = std::find_if(allFaces.begin(), allFaces.end(), [id](const std::pair<int,int>& p){ return p.second == id; });
                if (it != allFaces.end()) { validIds.push_back(id); doubleFaceIds.push_back(static_cast<double>(id)); }
            }

            if (validIds.empty() && !faceIds.empty() && globalOffset > 0) {
                for (int id : faceIds) {
                    int gmshTag = id - globalOffset;
                    if (gmshTag > 0 && std::find_if(allFaces.begin(), allFaces.end(),
                            [gmshTag](const std::pair<int,int>& p){ return p.second == gmshTag; }) != allFaces.end()) {
                        validIds.push_back(gmshTag);
                        doubleFaceIds.push_back(static_cast<double>(gmshTag));
                    } else {
                        fprintf(stderr, "[DEBUG mesh] WARNING: gmshTag=%d (from faceId=%d - offset=%d) not in Gmsh model, skipping\n",
                                gmshTag, id, globalOffset);
                    }
                }
            }
            fprintf(stderr, "[DEBUG mesh] group='%s' → validIds={", groupName.c_str());
            for (int v : validIds) fprintf(stderr, "%d ", v);
            fprintf(stderr, "}\n");

            if (!validIds.empty()) {
                gmsh::model::addPhysicalGroup(2, validIds, tagCounter);
                gmsh::model::setPhysicalName(2, tagCounter, groupName); 
                tagCounter++;

                if (!params.enableLocalField) continue; 

                int distFieldId = fieldTagCounter++;
                gmsh::model::mesh::field::add("Distance", distFieldId);
                safeSetFieldNumbers(distFieldId, "FacesList", "SurfacesList", doubleFaceIds);

                int threshFieldId = fieldTagCounter++;
                gmsh::model::mesh::field::add("Threshold", threshFieldId);
                safeSetFieldNumber(threshFieldId, "InField", "IField", distFieldId);
                safeSetFieldNumber(threshFieldId, "SizeMin", "LcMin", params.sizeMin);
                safeSetFieldNumber(threshFieldId, "SizeMax", "LcMax", params.sizeMax);
                safeSetFieldNumber(threshFieldId, "DistMin", "", params.distMin);
                safeSetFieldNumber(threshFieldId, "DistMax", "", params.distMax);
                safeSetFieldNumber(threshFieldId, "Sigmoid", "", 1.0);

                thresholdFieldTags.push_back(static_cast<double>(threshFieldId));
            }
        }

        if (!thresholdFieldTags.empty()) {
            int minFieldId = fieldTagCounter++;
            gmsh::model::mesh::field::add("Min", minFieldId);
            safeSetFieldNumbers(minFieldId, "FieldsList", "", thresholdFieldTags);
            gmsh::model::mesh::field::setAsBackgroundMesh(minFieldId);
        }

        std::vector<std::pair<int, int>> volumes;
        gmsh::model::getEntities(volumes, 3); 
        std::vector<int> volIds;
        for (const auto& vol : volumes) volIds.push_back(vol.second);
        
        fprintf(stderr, "[DEBUG mesh] volumes=%zu, before generate...\n", volIds.size());
        if (!volIds.empty()) {
            gmsh::model::addPhysicalGroup(3, volIds, tagCounter);
            gmsh::model::setPhysicalName(3, tagCounter, "Air_Domain");
            gmsh::model::mesh::generate(3);
        } else {
            gmsh::model::mesh::generate(2);
        }
        fprintf(stderr, "[DEBUG mesh] after generate\n");

        if (order == 2) { gmsh::model::mesh::setOrder(2); fprintf(stderr, "[DEBUG mesh] setOrder done\n"); }
        if (optimize) { gmsh::model::mesh::optimize("Netgen"); fprintf(stderr, "[DEBUG mesh] optimize done\n"); }

        std::vector<std::size_t> nodeTags;
        std::vector<double> coord, paramCoord;
        gmsh::model::mesh::getNodes(nodeTags, coord, paramCoord);
        std::size_t numNodes = nodeTags.size();
        fprintf(stderr, "[DEBUG mesh] nodes=%zu\n", numNodes);

        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags, nodeTagsPerElem;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElem);
        std::size_t numElements = 0;
        for (const auto& tags : elementTags) {
            numElements += tags.size();
        }
        fprintf(stderr, "[DEBUG mesh] elements=%zu\n", numElements);

        outLog = "========== 网格生成报告 ==========\n";
        outLog += "状态: ✅ 成功\n";
        outLog += "节点数 (Nodes): " + std::to_string(numNodes) + "\n";
        outLog += "单元数 (Elements): " + std::to_string(numElements) + "\n";
        outLog += "==================================\n";

        fprintf(stderr, "[DEBUG mesh] extracting surface mesh...\n");
        FaceMeshMap resultData = extractSurfaceMeshForOpenGL();
        fprintf(stderr, "[DEBUG mesh] extract done, faces=%zu\n", resultData.size());

        safeSetOption("Mesh.MshFileVersion", "", 2.2);
        gmsh::write(cacheMshPath);
        fprintf(stderr, "[DEBUG mesh] .msh written\n"); fflush(stderr);

        savePhysicalGroupsJson(cacheMshPath);
        fprintf(stderr, "[DEBUG mesh] .pgroups.json written\n"); fflush(stderr);

        gmsh::clear();
        fprintf(stderr, "[DEBUG mesh] returning data\n");

        return resultData;

    } catch (std::exception &e) {
        std::cerr << "\n[致命错误] " << e.what() << std::endl;
        cleanup();
        return FaceMeshMap{}; 
    }
}

// ====================================================================
// 导入外部网格文件并提取物理组
// ====================================================================
void MeshEngine::savePhysicalGroupsJson(const std::string& mshPath) {
    fprintf(stderr, "[DEBUG PG] start, getting groups...\n"); fflush(stderr);
    std::vector<std::pair<int,int>> dimTags;
    gmsh::model::getPhysicalGroups(dimTags);
    fprintf(stderr, "[DEBUG PG] %zu groups found\n", dimTags.size());
    if (dimTags.empty()) return;

    // 1. 收集所有 PG 信息
    struct PGInfo { int dim, tag; std::string name; std::vector<int> entities;
                    std::vector<std::size_t> elements; std::set<std::size_t> nodes; };
    std::vector<PGInfo> pgInfos;
    for (auto& dt : dimTags) {
        PGInfo info; info.dim = dt.first; info.tag = dt.second;
        gmsh::model::getPhysicalName(info.dim, info.tag, info.name);
        gmsh::model::getEntitiesForPhysicalGroup(info.dim, info.tag, info.entities);
        fprintf(stderr, "[DEBUG PG] group '%s' dim=%d entities=%zu\n", info.name.c_str(), info.dim, info.entities.size()); fflush(stderr);
        for (int et : info.entities) {
            std::vector<int> eTypes; std::vector<std::vector<std::size_t>> eTags, nTags;
            gmsh::model::mesh::getElements(eTypes, eTags, nTags, info.dim, et);
            for (size_t ti = 0; ti < eTags.size(); ti++) {
                info.elements.insert(info.elements.end(), eTags[ti].begin(), eTags[ti].end());
                info.nodes.insert(nTags[ti].begin(), nTags[ti].end());
            }
        }
        pgInfos.push_back(info);
    }

    // 2. 写 .pgroups.json
    std::string pgPath = mshPath + ".pgroups.json";
    std::ofstream pgf(pgPath);
    pgf << "{\n  \"mesh_file\": \"" << mshPath << "\",\n  \"groups\": [\n";
    for (size_t gi = 0; gi < pgInfos.size(); gi++) {
        auto& info = pgInfos[gi];
        pgf << "    {\"name\":\"" << info.name << "\", \"dim\":" << info.dim
            << ", \"tag\":" << info.tag << ", \"entities\":[";
        for (size_t ei = 0; ei < info.entities.size(); ei++) { if (ei) pgf << ","; pgf << info.entities[ei]; }
        pgf << "], \"elements\":[";
        for (size_t ei = 0; ei < info.elements.size(); ei++) { if (ei) pgf << ","; pgf << info.elements[ei]; }
        pgf << "], \"nodes\":[";
        size_t ni = 0; for (auto n : info.nodes) { if (ni++) pgf << ","; pgf << n; }
        pgf << "]}";
        if (gi < pgInfos.size() - 1) pgf << ",";
        pgf << "\n";
    }
    pgf << "  ]\n}\n";
    fprintf(stderr, "[DEBUG] .pgroups.json written: %zu groups\n", pgInfos.size());
}

FaceMeshMap MeshEngine::importMesh(const std::string& filePath, std::vector<ImportedGroup>& outGroups) {
    outGroups.clear();
    try {
        if (!s_isGmshInit) { gmsh::initialize(); s_isGmshInit = true; }
        gmsh::clear();
        std::cout << "[INFO] 正在解析网格文件: " << filePath << std::endl;
        gmsh::open(filePath);

        std::vector<std::pair<int, int>> dimTags;
        gmsh::model::getPhysicalGroups(dimTags);
        if (dimTags.empty())
            std::cout << "[INFO] 文件中无物理组信息" << std::endl;

        for (const auto& dt : dimTags) {
            int dim = dt.first, tag = dt.second;
            std::string name;
            gmsh::model::getPhysicalName(dim, tag, name);
            if (name.empty()) name = "Imported_Boundary_" + std::to_string(tag);
            outGroups.push_back({tag, name, dim});
        }

        FaceMeshMap renderData = extractSurfaceMeshForOpenGL();
        gmsh::clear();
        return renderData;

    } catch (std::exception &e) {
        std::cerr << "[解析错误] " << e.what() << std::endl;
        return FaceMeshMap{};
    }
}