#include "MolSymm.hpp"

std::string Molecule::detect_point_group(double tol) const {
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

    Eigen::MatrixXi atomic_numbers_to_compare = atomic_numbers.replicate(1, natoms);

    Eigen::MatrixXd distance_to_compare;
    {
        Eigen::MatrixXd coords_expanded_row = coordinates.replicate(natoms, 1);
        coords_expanded_row.resize(ncoords, natoms * natoms);
        Eigen::MatrixXd coords_expanded_col = coordinates.replicate(1, natoms);
        
        Eigen::MatrixXd diff = coords_expanded_col - coords_expanded_row;
        distance_to_compare = diff.colwise().norm().reshaped(natoms, natoms);
    }

    Eigen::VectorXd sort_method(natoms);
    for (int iatom = 0; iatom < natoms; ++ iatom) {
        std::iota(sort_method.data(), sort_method.data() + natoms, 0);
        std::sort(sort_method.data(), sort_method.data() + natoms,
            [iatom, &distance_to_compare](int i, int j) 
            { return distance_to_compare(i, iatom) < distance_to_compare(j, iatom); });
        distance_to_compare.col(iatom) = distance_to_compare.col(iatom)(sort_method).eval();
        atomic_numbers_to_compare.col(iatom) = atomic_numbers_to_compare.col(iatom)(sort_method).eval();
    }

    std::vector<bool> touched(natoms, false);
    std::vector<std::vector<int> > SEAs; // symmetry equavalent atoms
    for (int iatom = 0; iatom < natoms; ++ iatom) {
        if (touched[iatom]) continue;
        touched[iatom] = true;
        SEAs.push_back(std::vector<int>(1, iatom));
        for (int jatom = 0; jatom < natoms; ++ jatom) {
            if (touched[jatom]) continue;
            if (elements[iatom] != elements[jatom]) continue;
            if (!((atomic_numbers_to_compare.col(iatom).array() == atomic_numbers_to_compare.col(jatom).array()).all())) continue;
            if (((distance_to_compare.col(iatom).tail(natoms - 1) - distance_to_compare.col(jatom).tail(natoms - 1)).array().abs() > tol).any()) continue;

            touched[jatom] = true;
            SEAs.back().push_back(jatom);
        }
    }

    /*
    for (const std::vector<int> &SEA_group : SEAs) {
        for (int iatom : SEA_group) {
            fmt::print(" {:d}", iatom + 1);
        }
        fmt::print("\n");
    }
    */

    return "Not detected.";
}
