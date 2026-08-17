#ifndef COMPUTE_ELEMENT_VOLUME_HPP
#define COMPUTE_ELEMENT_VOLUME_HPP

#include <vector>
#include "mfem.hpp" // 假设 mfem 库的头文件是 mfem.hpp

// 计算网格中所有元素的面积
std::vector<double> ComputeElementVolume(mfem::ParMesh &mesh);

#endif // COMPUTE_ELEMENT_VOLUME_HPP