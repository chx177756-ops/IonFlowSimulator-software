#pragma once
#include "mfem.hpp"

void SolveConvectionSUPG(mfem::ParFiniteElementSpace &scalar_fes,
                         const mfem::ParGridFunction &nodeE_poisson,
                         const mfem::GridFunction &rho_previous_gf,
                         double rho_0,
                         mfem::ParGridFunction &rho_new,
                         mfem::Array<int> &ess_bdr);
