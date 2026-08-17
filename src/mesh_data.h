#pragma once
#include <vector>
#include <map>
#include <set>
#include <string>
#include <gmsh.h>
#include <unordered_map>

struct RenderableMesh {
    std::vector<float> vertices; 
    std::vector<unsigned int> indices; 
};

typedef std::map<int, RenderableMesh> FaceMeshMap;

// =========================================================
// 物理组专属网格参数结构体
// =========================================================
struct MeshGroupParams {
    bool enableLocalField = false; 
    
    // 【修改】：移除了已经废弃的 sampling 参数
    double sizeMin = 0.05;     
    double sizeMax = 8.0;      
    double distMin = 0.0;      
    double distMax = 30.0;     
};

inline FaceMeshMap extractGmshMeshToOpenGL() {
    FaceMeshMap faceMeshes;
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);

    if (nodeTags.empty()) return faceMeshes;

    std::unordered_map<std::size_t, unsigned int> tagToIndexMap;
    std::vector<float> globalVertices;
    globalVertices.reserve(coords.size()); 

    for (size_t i = 0; i < nodeTags.size(); ++i) {
        tagToIndexMap[nodeTags[i]] = static_cast<unsigned int>(i);
        globalVertices.push_back(static_cast<float>(coords[3 * i]));
        globalVertices.push_back(static_cast<float>(coords[3 * i + 1]));
        globalVertices.push_back(static_cast<float>(coords[3 * i + 2]));
    }

    std::vector<std::pair<int, int>> dimTags;
    gmsh::model::getEntities(dimTags, 2); 

    for (const auto& dt : dimTags) {
        int faceId = dt.second; 
        RenderableMesh rMesh;
        rMesh.vertices = globalVertices; 

        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags, nodeTagsPerElem;
        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElem, 2, faceId);

        for (size_t i = 0; i < elementTypes.size(); ++i) {
            int etype = elementTypes[i];
            const auto& nodes = nodeTagsPerElem[i];
            if (etype == 2 || etype == 9) { // 三角形 (1阶/2阶, 只取前3个角点)
                int npt = (etype == 2) ? 3 : 6;
                for (size_t j = 0; j + 2 < nodes.size(); j += npt) {
                    rMesh.indices.push_back(tagToIndexMap[nodes[j]]);
                    rMesh.indices.push_back(tagToIndexMap[nodes[j+1]]);
                    rMesh.indices.push_back(tagToIndexMap[nodes[j+2]]);
                }
            } else if (etype == 3 || etype == 10 || etype == 16) { // 四边形 → 两个三角形
                int npt = (etype == 3) ? 4 : (etype == 10 ? 9 : 8);
                for (size_t j = 0; j + 3 < nodes.size(); j += npt) {
                    unsigned int a = tagToIndexMap[nodes[j]];
                    unsigned int b = tagToIndexMap[nodes[j+1]];
                    unsigned int c = tagToIndexMap[nodes[j+2]];
                    unsigned int d = tagToIndexMap[nodes[j+3]];
                    rMesh.indices.push_back(a); rMesh.indices.push_back(b); rMesh.indices.push_back(c);
                    rMesh.indices.push_back(a); rMesh.indices.push_back(c); rMesh.indices.push_back(d);
                }
            }
        }
        if (!rMesh.indices.empty()) faceMeshes[faceId] = rMesh;
    }
    return faceMeshes;
}