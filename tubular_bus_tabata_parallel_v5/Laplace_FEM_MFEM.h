#ifndef LAPLACE_FEM_MFEM_H
#define LAPLACE_FEM_MFEM_H

#include "mfem.hpp"

mfem::ParGridFunction Laplace_FEM_MFEM(const mfem::ParMesh& mesh, double phi_0, 
                                    mfem::ParFiniteElementSpace& fespace);

#endif // LAPLACE_FEM_MFEM_H