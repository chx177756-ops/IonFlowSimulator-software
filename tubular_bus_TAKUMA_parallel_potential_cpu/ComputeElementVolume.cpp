#include "ComputeElementVolume.hpp"
#include "mfem.hpp"

// Function to compute the volume of each element in a parallel mesh.
// This version is specifically tailored for tetrahedral elements.
std::vector<double> ComputeElementVolume(mfem::ParMesh &mesh)
{
    // Initialize a vector to store the volume of each element.
    // The size of the vector is equal to the number of elements in the mesh.
    std::vector<double> volumes(mesh.GetNE());

    // Iterate over each element in the mesh.
    for (int i = 0; i < mesh.GetNE(); i++)
    {
        // Get the element transformation for the current element.
        // This object allows mapping from reference element to physical element.
        mfem::ElementTransformation *trans = mesh.GetElementTransformation(i);

        // Get the element itself to determine its type.
        mfem::Element *elem = mesh.GetElement(i);
        mfem::Geometry::Type geom_type;

        // Determine the geometry type based on the element type.
        switch (elem->GetType())
        {
            case mfem::Element::TETRAHEDRON:
                geom_type = mfem::Geometry::TETRAHEDRON;
                break;
            // If your mesh might contain other 3D element types (like HEXAHEDRON),
            // you should add cases for them here.
            // For example:
            // case mfem::Element::HEXAHEDRON:
            //     geom_type = mfem::Geometry::CUBE;
            //     break;
            default:
                // Abort if an unsupported element type is encountered.
                MFEM_ABORT("Unsupported element type for volume calculation. Only TETRAHEDRON is explicitly handled.");
        }

        // Get an integration rule for the determined geometry type.
        // For linear elements, order 0 (a single point at the centroid) or order 1
        // is typically sufficient because the Jacobian is constant.
        const mfem::IntegrationRule &ir = mfem::IntRules.Get(geom_type, 0); // Using order 0 for simplicity and efficiency

        // Set the integration point for the transformation.
        // For a single point, we just use the first (and only) point in the rule.
        trans->SetIntPoint(&ir.IntPoint(0));

        // Calculate the volume of the element.
        // trans->Weight() returns the Jacobian determinant multiplied by the weight
        // of the integration point in the reference element, effectively providing
        // the weighted volume contribution for that point. For linear elements and
        // a single integration point with weight 1, this directly gives the volume.
        volumes[i] = trans->Weight() * (1 / 6.0); // For linear tetrahedron, this is enough.
                                      // It's |det(J)| * weight_in_ref_space.
                                      // For a standard reference tetrahedron, its volume is 1/6.
                                      // MFEM's Weight() typically accounts for this reference scaling.
    }

    return volumes;
}