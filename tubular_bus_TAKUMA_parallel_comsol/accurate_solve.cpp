#include "accurate_solve.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <numeric> // For std::accumulate or similar if needed for integration

// Helper function for f(r)
static double calculate_f(double r, double A_val, double Eonset_val, double a_val);

void accurate_solve(
    const std::vector<double>& param,
    mfem::ParMesh& mesh,
    mfem::Vector& rho_values,
    std::vector<double>& r_coords 
) {
    // Parse parameters
    double a = param[0]; // Renamed r1 to 'a' for clarity matching MATLAB
    double b = param[1]; // Renamed r2 to 'b' for clarity matching MATLAB
    double V0 = param[2] * 1e3;  // V, Renamed phi_0 to V0 for clarity matching MATLAB
    double epsilon = param[3]; // Renamed epsilon_0 to 'epsilon' for clarity matching MATLAB
    double Eonset = param[4] * 1e3;  // V/m, Renamed E0 to 'Eonset' for clarity matching MATLAB
    // double c = 0.099935938490612;
    double c = 0.001984530974949; // This 'c' seems to be a constant, not from param

    // Calculate A
    double A = ( std::pow(a * b * V0 / (b - a), 2) - std::pow(a, 4) * std::pow(Eonset, 2) ) / ( std::pow(c, 3) - std::pow(a, 3) );

    // Get mesh coordinates and calculate r_coords using GetVertex
    int numVertices = mesh.GetNV(); // Get the number of vertices in the mesh
    r_coords.resize(numVertices);

    for (int i = 0; i < numVertices; ++i) {
        const double* v = mesh.GetVertex(i); // Get a pointer to the vertex coordinates
        double x = v[0];
        double y = v[1];
        double z = v[2]; 
        r_coords[i] = std::sqrt(x*x + y*y + z*z);
    }

    // Resize output vectors
    rho_values.SetSize(r_coords.size()); // mfem::Vector uses SetSize

    // Loop through each radial coordinate to calculate values
    for (size_t i = 0; i < r_coords.size(); ++i) {
        double r_val = r_coords[i];

        // Ensure f(r) is not zero to avoid division by zero in rho and E
        // Also handle the sqrt argument being non-negative
        double f_r = calculate_f(r_val, A, Eonset, a);
        if (f_r < 1e-12) { // Check for very small f_r to prevent division by zero
            // Handle error or set values to a safe default (e.g., 0)
            std::cerr << "Warning: f(r) is too small at r = " << r_val << ". Potential division by zero." << std::endl;
            rho_values[i] = 0.0;
            continue;
        }
        rho_values[i] = 1.5 * epsilon * A / f_r * 1e-3 * (1.0 / epsilon); // Apply 1e-3 scaling here
    }


}

// Helper function definitions

static double calculate_f(double r, double A_val, double Eonset_val, double a_val) {
    // f = @(r) sqrt( A * r.^3 + Eonset^2 * a^4 - A * a^3 );
    double arg_sqrt = A_val * std::pow(r, 3) + std::pow(Eonset_val, 2) * std::pow(a_val, 4) - A_val * std::pow(a_val, 3);
    if (arg_sqrt < 0) {
        std::cerr << "Error: Argument to sqrt in calculate_f is negative: " << arg_sqrt << " at r = " << r << std::endl;
        return 0.0;
    }
    return std::sqrt(arg_sqrt);
}

