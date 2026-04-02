#include "MolSymm.hpp"

std::string Molecule::detect_point_group(double tol) const {
    // quick return
    if (natoms == 0) throw std::runtime_error("Error: you should load a molecule first.");
    if (natoms == 1) return "Kh";
    if (natoms == 2) return elements[0] == elements[1] ? "Dinfh" : "Cinfv";

    Eigen::MatrixXd coords_centered = coordinates.colwise() - coordinates.rowwise().mean();

    Eigen::Matrix3d moments_of_inertia_tensor(Eigen::Matrix3d::Zero());

    // get moments of inertia
    for (int iatom = 0; iatom < natoms; ++ iatom) {
        moments_of_inertia_tensor += atomic_weights[iatom] * (
            coords_centered.col(iatom).transpose() * coords_centered.col(iatom) * Eigen::Matrix3d::Identity() 
          - coords_centered.col(iatom) * coords_centered.col(iatom).transpose());
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(moments_of_inertia_tensor);
    Eigen::Vector3d moments_of_inertia = eigen_solver.eigenvalues();
    Eigen::Matrix3d principal_axes = eigen_solver.eigenvectors();
    if (principal_axes.determinant() < 0.) principal_axes.col(0) = - principal_axes.col(0);
    coords_centered = principal_axes.transpose() * coords_centered; // rotate principal axes to x y z

    // detect SEA
    Eigen::MatrixXi atomic_numbers_to_compare = atomic_numbers.replicate(1, natoms);

    Eigen::MatrixXd distance_to_compare;
    {
        Eigen::MatrixXd coords_expanded_row = coordinates.replicate(natoms, 1);
        coords_expanded_row.resize(ncoords, natoms * natoms);
        Eigen::MatrixXd coords_expanded_col = coordinates.replicate(1, natoms);
        
        Eigen::MatrixXd diff = coords_expanded_col - coords_expanded_row;
        distance_to_compare = diff.colwise().norm().reshaped(natoms, natoms);
    }

    // argsort
    Eigen::VectorXd sort_method(natoms);
    for (int iatom = 0; iatom < natoms; ++ iatom) {
        std::iota(sort_method.data(), sort_method.data() + natoms, 0);
        std::sort(sort_method.data(), sort_method.data() + natoms,
            [iatom, &distance_to_compare](int i, int j) 
            { return distance_to_compare(i, iatom) < distance_to_compare(j, iatom); });
        distance_to_compare.col(iatom) = distance_to_compare.col(iatom)(sort_method).eval(); // must be a .eval() or there is aliasing issue.
        atomic_numbers_to_compare.col(iatom) = atomic_numbers_to_compare.col(iatom)(sort_method).eval(); // must be a .eval() or there is aliasing issue.
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

    Eigen::MatrixXd coords_operated(coords_centered);

    std::function<bool()> is_sym_okay = [&touched, &SEAs, &coords_centered, &coords_operated, tol]() {
        std::fill(touched.begin(), touched.end(), false);
        bool sym_okay = true;
        for (const std::vector<int>& SEA_group : SEAs) {
            for (int iatom : SEA_group) {
                for (int jatom : SEA_group) {
                    if (((coords_operated.col(iatom) - coords_centered.col(jatom)).array().abs() <= tol).all()) {
                        touched[iatom] = true;
                        break;
                    }
                }
                if (!touched[iatom]) {
                    sym_okay = false;
                    break;
                }
            }
            if (!sym_okay) break;
        }
        return sym_okay;
    };

    const double inertia_tol = tol * 2. * coords_centered.colwise().norm() * atomic_weights;
    if (moments_of_inertia[0] <= inertia_tol && moments_of_inertia[2] - moments_of_inertia[1] <= inertia_tol) {
        // linear, I_A = 0, I_B = I_C
        // fmt::print("{:s}\n", "linear");
        // Cinfv or Dinfh
        // only need to check symmetry center

        // coords_operated = - coords_centered;
        // return is_sym_okay() ? "Dinfh" : "Cinfv";

        // the molecule is on axis x

        bool sym_okay = true;
        for (const std::vector<int>& SEA_group : SEAs) {
            if (SEA_group.size() == 2) {
                if (std::abs(coords_centered(coord_x, SEA_group[0]) + coords_centered(coord_x, SEA_group[1])) / 2. > tol) {
                    sym_okay = false;
                    break;
                }
            } else if (SEA_group.size() == 1) {
                if (std::abs(coords_centered(coord_x, SEA_group[0])) > tol) {
                    sym_okay = false;
                    break;
                }
            } else {
                throw std::runtime_error("Error: this should never happen.");
            }
        }
        return sym_okay ? "Dinfh" : "Cinfv";

    } else if (moments_of_inertia[1] - moments_of_inertia[0] <= inertia_tol && moments_of_inertia[2] - moments_of_inertia[1] <= inertia_tol) {
        // more than one main-axes where n > 2, I_A = I_B = I_C, a.k.a. "spherial-like"
        // T, Td, Th, O, Oh, I, Ih
        fmt::print("{:s}\n", "spherial-like");

    } else if (moments_of_inertia[1] - moments_of_inertia[0] <= inertia_tol || moments_of_inertia[2] - moments_of_inertia[1] <= inertia_tol) {
        // symmetric, I_A = I_B \ne I_C or I_A \ne I_B = I_C
        // Dnd for n>= 2, (Cn, Cnh, Cnv, Dn, Dnh) for n > 2, Cni for (n > 1 and n is odd, a.k.a. S(2n)) and Sn for 4 | n
        fmt::print("{:s}\n", "symmetric");

        // rotate the principal axis corresponding to the unequivalent moment of inertia to x axis
        if (moments_of_inertia[1] - moments_of_inertia[0] <= inertia_tol) {
            Eigen::VectorXd swp = - coords_centered.row(coord_x);
            coords_centered.row(coord_x) = coords_centered.row(coord_z);
            coords_centered.row(coord_z) = swp;
        }

        // detect the main axis n>2
        // first find the maximum possible Cn axis
        int max_Cn_try = 2;
        for (const std::vector<int>& SEA_group : SEAs) {
            if (SEA_group.size() <= max_Cn_try) continue;
            std::vector<double> compare_x(SEA_group.size());
            for (int i = 0; i < SEA_group.size(); ++ i) {
                compare_x[i] = coords_centered(coord_x, SEA_group[i]);
            }
            std::sort(compare_x.begin(), compare_x.end());
            int max_Cn_try_current = 1;
            std::vector<int> max_Cn_try_list;
            for (int i = 1; i < compare_x.size(); ++ i) {
                if (compare_x[i] - compare_x[i - 1] <= tol) {
                    ++ max_Cn_try_current;
                } else {
                    max_Cn_try_list.push_back(max_Cn_try_current);
                    max_Cn_try_current = 1;
                }
            }
            max_Cn_try_list.push_back(max_Cn_try_current);
            for (int max_Cn_try_current : max_Cn_try_list) {
                if (max_Cn_try_current > max_Cn_try) max_Cn_try = max_Cn_try_current;
            }
        }

        /*
        max_Cn_try = 2
        for SEA_group in SEAs:
            if len(SEA_group) <= max_Cn_try: continue
            # at least 2 elements in compare_x
            compare_x = np.empty((len(SEA_group),), dtype=np.double)
            for i, iatom in enumerate(SEA_group):
                compare_x[i] = coords_centered[iatom, coord_x]
            compare_x.sort()
            max_Cn_try_current: int = 1
            max_Cn_try_list: list = []
            for i in range(1, len(compare_x)):
                if compare_x[i] - compare_x[i - 1] <= tol:
                    max_Cn_try_current += 1
                else:
                    max_Cn_try_list.append(max_Cn_try_current)
                    max_Cn_try_current = 1
            max_Cn_try_list.append(max_Cn_try_current)
            for max_Cn_try_current in max_Cn_try_list:
                if max_Cn_try_current > max_Cn_try: max_Cn_try = max_Cn_try_current
        */

    } else {
        // asymmetric, I_A \ne I_B \ne I_C
        // D2, D2h, C2, C2h, C2v, C1, Ci, Cs
        // fmt::print("{:s}\n", "asymmetric");
        bool has_x_C2, has_y_C2, has_z_C2, has_xOy_mirror, has_yOz_mirror, has_zOx_mirror, has_sym_center;

        coords_operated.row(coord_x) =   coords_centered.row(coord_x);
        coords_operated.row(coord_y) = - coords_centered.row(coord_y);
        coords_operated.row(coord_z) = - coords_centered.row(coord_z);
        has_x_C2 = is_sym_okay();
        coords_operated.row(coord_x) = - coords_centered.row(coord_x);
        coords_operated.row(coord_y) =   coords_centered.row(coord_y);
        has_y_C2 = is_sym_okay();
        coords_operated.row(coord_y) = - coords_centered.row(coord_y);
        coords_operated.row(coord_z) =   coords_centered.row(coord_z);
        has_z_C2 = is_sym_okay();

        // there can only be 0, 1 or 3 C2 in this situation
        if (has_x_C2 && has_y_C2 && has_z_C2) {
            // D2, D2h
            coords_operated.row(coord_x) =   coords_centered.row(coord_x);
            has_zOx_mirror = is_sym_okay();
            coords_operated.row(coord_x) = - coords_centered.row(coord_x);
            coords_operated.row(coord_y) =   coords_centered.row(coord_y);
            has_yOz_mirror = is_sym_okay();
            coords_operated.row(coord_x) =   coords_centered.row(coord_x);
            coords_operated.row(coord_z) = - coords_centered.row(coord_z);
            has_xOy_mirror = is_sym_okay();
            return has_xOy_mirror && has_yOz_mirror && has_zOx_mirror ? "D2h" : "D2";
        }

        // C2, C2h, C2v, C1, Ci, Cs
        // rotate the C2 axis to x axis
        if (has_z_C2) {
            Eigen::RowVectorXd swp = - coords_centered.row(coord_x);
            coords_centered.row(coord_x) = coords_centered.row(coord_z);
            coords_centered.row(coord_z) = swp;
            has_z_C2 = false;
            has_x_C2 = true;
        } else if (has_y_C2) {
            Eigen::RowVectorXd swp = - coords_centered.row(coord_y);
            coords_centered.row(coord_y) = coords_centered.row(coord_x);
            coords_centered.row(coord_x) = swp;
            has_y_C2 = false;
            has_x_C2 = true;
        }

        if (has_x_C2) {
            // C2, C2h, C2v
            coords_operated.row(coord_x) = - coords_centered.row(coord_x);
            coords_operated.row(coord_y) =   coords_centered.row(coord_y);
            coords_operated.row(coord_z) =   coords_centered.row(coord_z);
            has_yOz_mirror = is_sym_okay();
            if (has_yOz_mirror) return "C2h";
            // C2, C2v
            coords_operated.row(coord_x) =   coords_centered.row(coord_x);
            coords_operated.row(coord_z) = - coords_centered.row(coord_z);
            has_xOy_mirror = is_sym_okay();
            return has_xOy_mirror ? "C2v" : "C2";
        }

        // C1, Ci, Cs
        coords_operated.row(coord_x) =   coords_centered.row(coord_x);
        coords_operated.row(coord_y) = - coords_centered.row(coord_y);
        coords_operated.row(coord_z) =   coords_centered.row(coord_z);
        has_zOx_mirror = is_sym_okay();
        coords_operated.row(coord_x) = - coords_centered.row(coord_x);
        coords_operated.row(coord_y) =   coords_centered.row(coord_y);
        has_yOz_mirror = is_sym_okay();
        coords_operated.row(coord_x) =   coords_centered.row(coord_x);
        coords_operated.row(coord_z) = - coords_centered.row(coord_z);
        has_xOy_mirror = is_sym_okay();
        if (has_xOy_mirror || has_yOz_mirror || has_zOx_mirror) return "Cs";

        // C1, Ci
        coords_operated.row(coord_x) = - coords_centered.row(coord_x);
        coords_operated.row(coord_y) = - coords_centered.row(coord_y);
        has_sym_center = is_sym_okay();
        return has_sym_center ? "Ci" : "C1";
    }

    return "Not detected.";
}
