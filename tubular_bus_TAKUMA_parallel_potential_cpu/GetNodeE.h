#ifndef GET_NODE_E_HPP
#define GET_NODE_E_HPP

#include <vector>
#include <mfem.hpp>

void GetNodeE(const mfem::DenseMatrix& elementE, 
              const mfem::ParMesh& mesh,
              const std::vector<double>& Area,
              mfem::ParFiniteElementSpace& fespace,
              mfem::ParGridFunction& nodeE,
              mfem::ParGridFunction& nodeEn);

#endif // GET_NODE_E_HPP
