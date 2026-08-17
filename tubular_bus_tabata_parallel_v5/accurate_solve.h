#ifndef ACCURATE_SOLVE_H
#define ACCURATE_SOLVE_H

#include <vector>
#include <cmath>
#include "mfem.hpp"

// 函数声明
void accurate_solve(
    const std::vector<double>& param,
    mfem::ParMesh& mesh,
    mfem::Vector& rho_values,
    std::vector<double>& r_coords 
);

#endif // ACCURATE_SOLVE_H