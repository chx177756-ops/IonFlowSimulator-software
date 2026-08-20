#include "ParallelUpwindElementFinder.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstdlib> // srand, rand
#include <ctime>   // time
#include <limits>

using namespace mfem;

UpwindElementFinder::UpwindElementFinder(mfem::ParMesh &mesh,
                                         mfem::ParGridFunction &nodeE_poisson,
                                         mfem::ParFiniteElementSpace &vector_fes,
                                         const std::vector<mfem::DenseMatrix> &global_dshapes)
    : mesh(mesh),
      nodeE_poisson(nodeE_poisson),
      vector_fes(vector_fes),
      global_dshapes(global_dshapes), // [优化] 复用预计算
      numNod(mesh.GetNV()),
      dim(mesh.Dimension())
{
    MPI_Comm_rank(mesh.GetComm(), &rank);

    // build vertex -> element table
    vertex_to_elem.reset(mesh.GetVertexToElementTable());

    // seed random for fallback selection
    std::srand(static_cast<unsigned>(std::time(nullptr)) + rank);
}

std::vector<UpwindInfo> UpwindElementFinder::ComputeUpwindElements()
{
    std::vector<UpwindInfo> upwindInfos(numNod);
    for (int NodeIndex = 0; NodeIndex < numNod; NodeIndex++)
    {
        upwindInfos[NodeIndex] = FindUpwindElementForNode(NodeIndex);
    }
    SynchronizeSharedNodes(upwindInfos);
    return upwindInfos;
}

UpwindInfo UpwindElementFinder::FindUpwindElementForNode(int NodeIndex)
{
    // Get node position
    const double *targetPos = mesh.GetVertex(NodeIndex);

    // Get convection velocity at node
    mfem::Vector convectionVelocity(dim);
    for (int d = 0; d < dim; ++d)
    {
        int vdof = vector_fes.DofToVDof(NodeIndex, d);
        double val = nodeE_poisson(vdof);
        convectionVelocity(d) = val;
    }

    double Vx = convectionVelocity(0);
    double Vy = convectionVelocity(1);
    double Vz = (dim >= 3) ? convectionVelocity(2) : 0.0;

    double convNorm = std::sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
    const double EPS_V = 1e-14;

    // get neighbor elements
    int num_elems = vertex_to_elem->RowSize(NodeIndex);
    const int *elems = vertex_to_elem->GetRow(NodeIndex);

    // Velocity nearly zero logic (fallback)
    if (convNorm < EPS_V)
    {
        if (num_elems > 0)
        {
            int pick = std::rand() % num_elems;
            int chosen_elem = elems[pick];
            Element *elem = mesh.GetElement(chosen_elem);
            double centroid[3] = {0.0, 0.0, 0.0};
            int nv = elem->GetNVertices();
            for (int j = 0; j < nv; ++j)
            {
                const double *v = mesh.GetVertex(elem->GetVertices()[j]);
                centroid[0] += v[0]; centroid[1] += v[1]; centroid[2] += v[2];
            }
            centroid[0] /= nv; centroid[1] /= nv; centroid[2] /= nv;
            return UpwindInfo(chosen_elem, rank, 0.0, centroid[0], centroid[1], centroid[2]);
        }
        else
        {
            return UpwindInfo(-1, rank, -1e12, 0.0, 0.0, 0.0);
        }
    }

    // Normal case: use MATLAB's b,c,d criterion
    std::vector<int> candidates;
    candidates.reserve(num_elems);

    const double TOL = 1e-10;

    for (int ie = 0; ie < num_elems; ++ie)
    {
        int elem_id = elems[ie];
        Element *elem = mesh.GetElement(elem_id);

        int localIndex = -1;
        int nv = elem->GetNVertices();
        for (int k = 0; k < nv; ++k)
        {
            if (elem->GetVertices()[k] == NodeIndex)
            {
                localIndex = k;
                break;
            }
        }
        if (localIndex < 0) continue;

        // [优化] 直接读取预计算的形函数梯度，避免重复 Jacobian 求逆+CalcDShape
        const DenseMatrix &dshape_e = global_dshapes[elem_id];
        int ndof_e = dshape_e.Height();  // P1四面体 = 4

        // Apply MATLAB b,c,d test
        bool isUpwind = true;
        for (int loc = 0; loc < ndof_e; ++loc)
        {
            if (loc == localIndex) continue;
            double b_i = dshape_e(loc, 0); 
            double c_i = (dim >= 2) ? dshape_e(loc, 1) : 0.0;
            double d_i = (dim >= 3) ? dshape_e(loc, 2) : 0.0;

            double test_val = b_i * Vx + c_i * Vy + d_i * Vz;
            if (test_val > TOL)
            {
                isUpwind = false;
                break;
            }
        }

        if (isUpwind)
        {
            candidates.push_back(elem_id);
        }
    } 

    if (!candidates.empty())
    {
        int pick = 0;
        if ((int)candidates.size() > 1) pick = std::rand() % static_cast<int>(candidates.size());
        int chosen_elem = candidates[pick];

        Element *elemC = mesh.GetElement(chosen_elem);
        double centroid[3] = {0.0, 0.0, 0.0};
        int nvC = elemC->GetNVertices();
        for (int j = 0; j < nvC; ++j)
        {
            const double *v = mesh.GetVertex(elemC->GetVertices()[j]);
            centroid[0] += v[0]; centroid[1] += v[1]; centroid[2] += v[2];
        }
        centroid[0] /= nvC; centroid[1] /= nvC; centroid[2] /= nvC;

        return UpwindInfo(chosen_elem, rank, 1.0, centroid[0], centroid[1], centroid[2]);
    }

    return UpwindInfo(-1, rank, -1e12, 0.0, 0.0, 0.0);
}

void UpwindElementFinder::SynchronizeSharedNodes(std::vector<UpwindInfo> &upwindInfos)
{
    mfem::GroupCommunicator &group_comm = vector_fes.GroupComm();
    const int nverts = numNod;

    mfem::Array<int> candidate_flag(vector_fes.GetVSize());
    mfem::Array<int> elem_ids(vector_fes.GetVSize());
    mfem::Array<int> ranks_array(vector_fes.GetVSize());
    mfem::Vector centroid_x(vector_fes.GetVSize());
    mfem::Vector centroid_y(vector_fes.GetVSize());
    mfem::Vector centroid_z(vector_fes.GetVSize());

    for (int i = 0; i < nverts; i++)
    {
        candidate_flag[i] = (upwindInfos[i].elem_id >= 0) ? 1 : 0;
        elem_ids[i] = (upwindInfos[i].elem_id >= 0) ? upwindInfos[i].elem_id : -1000000000;
        ranks_array[i] = (upwindInfos[i].elem_id >= 0) ? upwindInfos[i].rank : -1000000000;
        centroid_x[i] = upwindInfos[i].centroid[0];
        centroid_y[i] = upwindInfos[i].centroid[1];
        centroid_z[i] = upwindInfos[i].centroid[2];
    }

    if (nverts > 0)
    {
        group_comm.Reduce<int>(candidate_flag.HostReadWrite(), mfem::GroupCommunicator::Max);
        group_comm.Bcast<int>(candidate_flag.HostReadWrite());

        group_comm.Reduce<int>(elem_ids.HostReadWrite(), mfem::GroupCommunicator::Max);
        group_comm.Bcast<int>(elem_ids.HostReadWrite());

        group_comm.Reduce<int>(ranks_array.HostReadWrite(), mfem::GroupCommunicator::Max);
        group_comm.Bcast<int>(ranks_array.HostReadWrite());

        group_comm.Reduce<double>(centroid_x.HostReadWrite(), mfem::GroupCommunicator::Max);
        group_comm.Bcast<double>(centroid_x.HostReadWrite());

        group_comm.Reduce<double>(centroid_y.HostReadWrite(), mfem::GroupCommunicator::Max);
        group_comm.Bcast<double>(centroid_y.HostReadWrite());

        group_comm.Reduce<double>(centroid_z.HostReadWrite(), mfem::GroupCommunicator::Max);
        group_comm.Bcast<double>(centroid_z.HostReadWrite());
    }

    for (int i = 0; i < nverts; i++)
    {
        if (candidate_flag[i] > 0) 
        {
            int eid = elem_ids[i];
            if (eid < 0 || eid == -1000000000)
            {
                upwindInfos[i].elem_id = -1;
                upwindInfos[i].rank = -1;
            }
            else
            {
                upwindInfos[i].elem_id = eid;
                upwindInfos[i].rank = ranks_array[i];
                upwindInfos[i].cosin_theta = 1.0; 
                upwindInfos[i].centroid[0] = centroid_x[i];
                upwindInfos[i].centroid[1] = centroid_y[i];
                upwindInfos[i].centroid[2] = centroid_z[i];
            }
        }
        else
        {
            upwindInfos[i].elem_id = -1;
            upwindInfos[i].rank = -1;
            upwindInfos[i].centroid[0] = 0.0;
        }
    }
}