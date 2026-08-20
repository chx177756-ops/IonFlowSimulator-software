#include "mfem.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>

using namespace mfem;

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <mesh_file> <num_parts> [output_prefix]" << std::endl;
        std::cout << "  mesh_file      : GMSH mesh file (e.g. mesh.msh)" << std::endl;
        std::cout << "  num_parts      : number of partitions" << std::endl;
        std::cout << "  output_prefix  : output file prefix (default: <mesh_file>.part)" << std::endl;
        return 1;
    }

    const char *mesh_file = argv[1];
    int num_parts = std::atoi(argv[2]);
    std::string prefix = (argc > 3) ? argv[3] : (std::string(mesh_file) + ".part");

    std::cout << "Loading mesh: " << mesh_file << " ..." << std::endl;
    Mesh serial_mesh(mesh_file, 1, 1);
    std::cout << "  Elements: " << serial_mesh.GetNE()
              << ", Vertices: " << serial_mesh.GetNV()
              << ", Dimension: " << serial_mesh.Dimension() << std::endl;

    std::cout << "Partitioning into " << num_parts << " parts ..." << std::endl;
    int *partitioning = serial_mesh.GeneratePartitioning(num_parts);
    MeshPartitioner partitioner(serial_mesh, num_parts, partitioning);

    for (int i = 0; i < num_parts; i++) {
        MeshPart mp;
        partitioner.ExtractPart(i, mp);

        std::ostringstream fname;
        fname << prefix << "." << std::setfill('0') << std::setw(6) << i;
        std::ofstream ofs(fname.str());
        ofs.precision(14);
        mp.Print(ofs);
        ofs.close();

        std::cout << "  Written: " << fname.str()
                  << " (" << mp.num_elements << " elements, "
                  << mp.num_vertices << " vertices)" << std::endl;
    }

    delete[] partitioning;  // MeshPartitioner uses MakeRef(..., false), 不接管内存

    // 保存串行网格顶点坐标 (用于并行加载时恢复输出顺序)
    {
        std::string vert_file = prefix + ".vertices";
        std::ofstream vf(vert_file);
        vf.precision(14);
        int nv = serial_mesh.GetNV();
        int sd = serial_mesh.SpaceDimension();
        vf << nv << " " << sd << "\n";
        for (int i = 0; i < nv; i++) {
            const double *v = serial_mesh.GetVertex(i);
            vf << v[0];
            for (int d = 1; d < sd; d++) vf << " " << v[d];
            vf << "\n";
        }
        vf.close();
        std::cout << "  Vertices reference: " << vert_file << " (" << nv << " vertices)" << std::endl;
    }

    std::cout << "Done. Parallel mesh files: " << prefix << ".*" << std::endl;
    return 0;
}
