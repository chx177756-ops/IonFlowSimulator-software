#ifndef TABATA_ASSEMBLY_GPU_H
#define TABATA_ASSEMBLY_GPU_H
#include <vector>
using GIdx = int;
void TabataGPUComputeValues(int numTet, int dim, int local_ndofs, int my_rank,
    const std::vector<int>&, const std::vector<double>&, const std::vector<double>&,
    const std::vector<GIdx>&, const std::vector<double>&, const std::vector<GIdx>&,
    const std::vector<int>&, const std::vector<double>&, const std::vector<double>&,
    const std::vector<int>&, const std::vector<GIdx>&, double**, int*);
#endif
