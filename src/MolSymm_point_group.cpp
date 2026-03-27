#include "MolSymm.hpp"

std::string Molecule::detect_point_group() const {
    if (natoms == 1) return "Kh";
    if (natoms == 2) return elements[0] == elements[1] ? "Dinfh" : "Cinfv";

    Eigen::MatrixXd coords_centerized = coordinates.colwise() - coordinates.rowwise().mean();

    Eigen::Matrix3d moments_of_inertia_tensor(Eigen::Matrix3d::Zero());

    for (int iatom = 0; iatom < natoms; ++ iatom) {
        moments_of_inertia_tensor += atomic_weights[iatom] * (
            coords_centerized.col(iatom).transpose() * coords_centerized.col(iatom) * Eigen::Matrix3d::Identity() 
          - coords_centerized.col(iatom) * coords_centerized.col(iatom).transpose());
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(moments_of_inertia_tensor);
    Eigen::Vector3d moments_of_inertia = eigen_solver.eigenvalues();
    // Eigen::Matrix3d principal_axes = eigen_solver.eigenvectors();

    std::cout << moments_of_inertia.transpose() << std::endl;

    return "Not detected.";
}
