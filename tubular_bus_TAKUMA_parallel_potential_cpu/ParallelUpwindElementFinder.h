#ifndef PARALLEL_UPWIND_ELEMENT_FINDER_H
#define PARALLEL_UPWIND_ELEMENT_FINDER_H

#include <mfem.hpp>
#include <memory>
#include <vector>

class UpwindInfo
{
public:
    int elem_id;         // 上流元的ID (>=0 表示有效)
    int rank;            // 拥有上流元的进程Rank
    double cosin_theta;  // 占位字段
    double centroid[3];  // 上流元的质心坐标

    UpwindInfo(int e = -1, int r = -1, double c = -1e12, double cx = 0.0, double cy = 0.0, double cz = 0.0)
        : elem_id(e), rank(r), cosin_theta(c), centroid{cx, cy, cz} {}
};

class UpwindElementFinder
{
private:
    mfem::ParMesh &mesh;
    mfem::ParGridFunction &nodeE_poisson;
    mfem::ParFiniteElementSpace &vector_fes;
    const std::vector<mfem::DenseMatrix> &global_dshapes; // [新增] 全局预计算梯度引用
    int numNod;
    int dim;
    int rank;
    std::unique_ptr<mfem::Table> vertex_to_elem;

public:
    // [修改] 构造函数参数新增 global_dshapes
    UpwindElementFinder(mfem::ParMesh &mesh,
                        mfem::ParGridFunction &nodeE_poisson,
                        mfem::ParFiniteElementSpace &vector_fes,
                        const std::vector<mfem::DenseMatrix> &global_dshapes);
    std::vector<UpwindInfo> ComputeUpwindElements();
    UpwindInfo FindUpwindElementForNode(int NodeIndex);
    void SynchronizeSharedNodes(std::vector<UpwindInfo> &upwindInfos);
};

#endif // PARALLEL_UPWIND_ELEMENT_FINDER_H