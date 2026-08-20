#include "GetNodeE.h"
#include <cmath>
#include <mfem.hpp>

void GetNodeE(const mfem::DenseMatrix &elementE,
              const mfem::ParMesh &mesh,
              const std::vector<double> &Volume,
              mfem::ParFiniteElementSpace &fespace,
              mfem::ParGridFunction &nodeE,
              mfem::ParGridFunction &nodeEn)
{
    // check
    MFEM_VERIFY(fespace.GetVDim() == 3, "FESpace must have vdim == 3");
    MFEM_VERIFY(nodeE.Size() == 3 * fespace.GetNDofs(), "Size mismatch");

    // [优化] 使用 static 变量避免在循环中重复分配内存
    // 注意：这在纯 MPI 环境下是安全的（每个进程独立）。
    // 如果使用了 OpenMP 混合并行，这里需要改为 thread_local 或者传入工作空间对象。
    static mfem::Vector weight;
    static mfem::Vector ex, ey, ez;
    static mfem::Vector temp_weight, temp_ex, temp_ey, temp_ez;

    int vSize = fespace.GetVSize();

    // 检查尺寸是否匹配（处理网格可能的动态变化，虽然后者很少见）
    if (weight.Size() != vSize) {
        weight.SetSize(vSize);
        ex.SetSize(vSize);
        ey.SetSize(vSize);
        ez.SetSize(vSize);
        
        // 同时也预分配临时通信向量
        temp_weight.SetSize(vSize);
        temp_ex.SetSize(vSize);
        temp_ey.SetSize(vSize);
        temp_ez.SetSize(vSize);
    }

    // 重置为 0
    nodeE = 0.0;
    nodeEn = 0.0;
    weight = 0.0; ex = 0.0; ey = 0.0; ez = 0.0;

    // 累加元素贡献到节点
    for (int i = 0; i < mesh.GetNE(); i++)
    {
        const mfem::Element *el = mesh.GetElement(i);
        const int *nodes = el->GetVertices();
        double V = Volume[i];
        
        // 预取值以减少数组访问
        double E_x = elementE(i, 0);
        double E_y = elementE(i, 1);
        double E_z = elementE(i, 2);

        for (int j = 0; j < el->GetNVertices(); j++)
        {
            int node = nodes[j];
            ex(node) += V * E_x;
            ey(node) += V * E_y;
            ez(node) += V * E_z;
            weight(node) += V;
        }
    }
  
    // 执行规约操作 (Reduce + Broadcast)
    mfem::GroupCommunicator &gcomm = fespace.GroupComm();

    // 复用预分配的临时向量
    temp_ey = ey;
    gcomm.Reduce<double>(temp_ey.HostReadWrite(), mfem::GroupCommunicator::Sum);
    gcomm.Bcast<double>(temp_ey.HostReadWrite());
    ey = temp_ey;

    temp_ex = ex;
    gcomm.Reduce<double>(temp_ex.HostReadWrite(), mfem::GroupCommunicator::Sum);
    gcomm.Bcast<double>(temp_ex.HostReadWrite());
    ex = temp_ex;

    temp_ez = ez;
    gcomm.Reduce<double>(temp_ez.HostReadWrite(), mfem::GroupCommunicator::Sum);
    gcomm.Bcast<double>(temp_ez.HostReadWrite());
    ez = temp_ez;

    temp_weight = weight;
    gcomm.Reduce<double>(temp_weight.HostReadWrite(), mfem::GroupCommunicator::Sum);
    gcomm.Bcast<double>(temp_weight.HostReadWrite());
    weight = temp_weight; 

    // 同步共享节点
    nodeE.Update();  
    nodeEn.Update(); 

    nodeE.ExchangeFaceNbrData();  
    nodeEn.ExchangeFaceNbrData(); 

    // 归一化处理
    int ndofs = fespace.GetNDofs();
    for (int i = 0; i < ndofs; i++)
    {
        double w = weight(i);
        if (w > 1e-16) // 避免除以零
        {
            double inv_w = 1.0 / w;
            double val_ex = ex(i) * inv_w;
            double val_ey = ey(i) * inv_w;
            double val_ez = ez(i) * inv_w;

            nodeE(i) = val_ex;                  
            nodeE(i + ndofs) = val_ey;    
            nodeE(i + 2 * ndofs) = val_ez; 
            nodeEn(i) = std::sqrt(val_ex*val_ex + val_ey*val_ey + val_ez*val_ez); // 使用 sqrt 替代 hypot 可能稍微快一点
        }
    }
}