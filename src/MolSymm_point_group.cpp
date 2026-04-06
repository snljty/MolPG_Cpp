#include "MolSymm.hpp"

std::string Molecule::detect_point_group(double tol) const {
    // quick return
    if (natoms == 0) throw std::runtime_error("Error: you should load a molecule first.");
    if (natoms == 1) {
        new_coordinates.fill(0.);
        return "Kh";
    }
    if (natoms == 2) {
        double bond_length = (coordinates.col(1) - coordinates.col(0)).norm();
        new_coordinates.fill(0.);
        new_coordinates(coord_x, 0) = - bond_length / 2.;
        new_coordinates(coord_x, 1) = bond_length / 2.;
        return elements[0] == elements[1] ? "Dinfh" : "Cinfv";
    }

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
    new_coordinates = coords_centered;

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
        // more than one main-axes where n > 2, I_A = I_B = I_C, a.k.a. "spherical-like"
        // T, Td, Th, O, Oh, I, Ih
        // fmt::print("{:s}\n", "spherical-like");

        std::function<void(const Eigen::Vector3d&)> rotate_half_around_axis = 
        [&coords_centered, &coords_operated](const Eigen::Vector3d& axis) {
            // assume axis is already normalized
            // \cos\pi = -1, \sin\pi = 0
            // R = \begin{bmatrix}
            // \cos\theta + (1-\cos\theta){n_x}^2 & (1-\cos\theta)n_xn_y-\sin\theta n_z & (1-\cos\theta)n_xn_z+\sin\theta n_y \\
            // (1-\cos\theta)n_yn_x+\sin\theta n_z & \cos\theta + (1-\cos\theta){n_y}^2 & (1-\cos\theta)n_yn_z-\sin\theta n_x \\
            // (1-\cos\theta)n_zn_x-\sin\theta n_y & (1-\cos\theta)n_zn_y+\sin\theta n_x & \cos\theta + (1-\cos\theta){n_z}^2 \\
            // \end{bmatrix}
            Eigen::Matrix3d rot_mat = 2. * axis * axis.transpose() - Eigen::Matrix3d::Identity();
            coords_operated = rot_mat * coords_centered;
        };

        std::function<void(const Eigen::Vector3d&)> rotate_quarter_around_axis = 
        [&coords_centered, &coords_operated](const Eigen::Vector3d& axis) {
            // assume axis is already normalized
            // \cos\frac\pi2 = 0, \sin\frac\pi2 = 1
            Eigen::Matrix3d cross_mat;
            cross_mat << 0., - axis[coord_z], axis[coord_y], 
                         axis[coord_z], 0., - axis[coord_x], 
                         - axis[coord_y], axis[coord_x], 0.;
            Eigen::Matrix3d rot_mat = axis * axis.transpose() + cross_mat;
            coords_operated = rot_mat * coords_centered;
        };

        std::function<void(const Eigen::Vector3d&)> flip_against_plane = 
        [&coords_centered, &coords_operated](const Eigen::Vector3d& normal_axis) {
            // assume normal_axis is already normalized
            Eigen::MatrixXd projection = normal_axis * (normal_axis.transpose() * coords_centered);
            coords_operated = coords_centered - 2. * projection;
        };

        // T series have 3 C2, O series have 6 individual C2 and 3 C4, I series have 15 C2

        Eigen::MatrixXd all_C2(ncoords, 15);
        int num_C2_found = 0;
        Eigen::Vector3d axis_point;
        double axis_point_norm;

        // first, check C2 through center of two SEAs
        for (const std::vector<int>& SEA_group : SEAs) {
            if (SEA_group.size() < 2) continue;
            for (int iatom : SEA_group) {
                for (int jatom : SEA_group) {
                    if (iatom >= jatom) continue;
                    axis_point = (coords_centered.col(iatom) + coords_centered.col(jatom)) / 2.;
                    axis_point_norm = axis_point.norm();
                    if (axis_point_norm <= tol) continue;
                    axis_point /= axis_point_norm;
                    bool already_found = false;
                    for (int C2_index = 0; C2_index < num_C2_found; ++ C2_index) {
                        if ((axis_point + all_C2.col(C2_index)).norm() <= tol || 
                            (axis_point - all_C2.col(C2_index)).norm() <= tol) {
                            already_found = true;
                            break;
                        }
                    }
                    if (!already_found) {
                        rotate_half_around_axis(axis_point);
                        if (is_sym_okay()) {
                            if (num_C2_found >= 15) throw std::runtime_error("Error: this should never happen.");
                            all_C2.col(num_C2_found) = axis_point;
                            ++ num_C2_found;
                        }
                    }
                }
            }
        }

        // second, check C2 through each atom
        for (int iatom = 0; iatom > natoms; ++ iatom) {
            axis_point_norm = coords_centered.col(iatom).norm();
            if (axis_point_norm <= tol) continue;
            axis_point = coords_centered.col(iatom) / axis_point_norm;
            bool already_found = false;
            for (int C2_index = 0; C2_index < num_C2_found; ++ C2_index) {
                if ((axis_point + all_C2.col(C2_index)).norm() <= tol || 
                    (axis_point - all_C2.col(C2_index)).norm() <= tol) {
                    already_found = true;
                    break;
                }
            }
            if (!already_found) {
                rotate_half_around_axis(axis_point);
                if (is_sym_okay()) {
                    if (num_C2_found >= 15) throw std::runtime_error("Error: this should never happen.");
                    all_C2.col(num_C2_found) = axis_point;
                    ++ num_C2_found;
                }
            }
        }

        if (num_C2_found == 3) {
            // T, Td, Th
            flip_against_plane(all_C2.col(0));
            bool has_sigma_h = is_sym_okay();
            Eigen::Vector3d mirror_point = (all_C2.col(0) + all_C2.col(1)) / M_SQRT2;
            flip_against_plane(mirror_point);
            bool has_sigma_d = is_sym_okay();
            Eigen::Matrix3d rot_mat = all_C2.leftCols(ncoords).transpose();
            if (rot_mat.determinant() < 0.) rot_mat.row(coord_x) = - rot_mat.row(coord_x);
            coords_centered = rot_mat * coords_centered;
            new_coordinates = coords_centered;
            return has_sigma_h ? "Th" : has_sigma_d ? "Td" : "T";
        } else if (num_C2_found == 9) {
            // O, Oh
            Eigen::Vector3i C4_index;
            int num_C4_found = 0;
            for (int C2_index = 0; C2_index < 9; ++ C2_index) {
                rotate_quarter_around_axis(all_C2.col(C2_index));
                if (is_sym_okay()) {
                    C4_index[num_C4_found] = C2_index;
                    ++ num_C4_found;
                }
            }
            if (num_C4_found != 3) throw std::runtime_error("Error: this should never happen.");
            flip_against_plane(all_C2.col(C4_index[0]));
            bool has_sigma_h = is_sym_okay();
            Eigen::Matrix3d rot_mat;
            for (int i = 0; i < ncoords; ++ i) {
                rot_mat.row(i) = all_C2.col(C4_index[i]).transpose();
            }
            if (rot_mat.determinant() < 0.) rot_mat.row(coord_x) = - rot_mat.row(coord_x);
            coords_centered = rot_mat * coords_centered;
            new_coordinates = coords_centered;
            return has_sigma_h ? "Oh" : "O";
        } else if (num_C2_found == 15) {
            // I, Ih
            flip_against_plane(all_C2.col(0));
            bool has_sigma_h = is_sym_okay();
            Eigen::Vector3i C2_use_index(Eigen::Vector3i::Zero());
            bool found_second_orthogonal_C2 = false;
            for (C2_use_index[1] = 1; C2_use_index[1] < 15; ++ C2_use_index[1]) {
                if (std::abs(all_C2.col(C2_use_index[1]).transpose() * all_C2.col(0)) <= tol) {
                    found_second_orthogonal_C2 = true;
                    break;
                }
            }
            if (!found_second_orthogonal_C2) throw std::runtime_error("Error: this should never happen.");
            bool found_third_orthogonal_C2 = false;
            for (C2_use_index[2] = 1; C2_use_index[2] < 15; ++ C2_use_index[2]) {
                if (C2_use_index[2] == C2_use_index[1]) continue;
                if (std::abs(all_C2.col(C2_use_index[2]).transpose() * all_C2.col(0)) <= tol && 
                    std::abs(all_C2.col(C2_use_index[2]).transpose() * all_C2.col(C2_use_index[1])) <= tol) {
                    found_third_orthogonal_C2 = true;
                    break;
                }
            }
            if (!found_third_orthogonal_C2) throw std::runtime_error("Error: this should never happen.");
            Eigen::Matrix3d rot_mat;
            for (int i = 0; i < ncoords; ++ i) {
                rot_mat.row(i) = all_C2.col(C2_use_index[i]);
            }
            if (rot_mat.determinant() < 0.) rot_mat.row(coord_x) = - rot_mat.row(coord_x);
            coords_centered = rot_mat * coords_centered;
            new_coordinates = coords_centered;
            return has_sigma_h ? "Ih" : "I";
        } else {
            throw std::runtime_error("Error: this should never happen.");
        }

    } else if (moments_of_inertia[1] - moments_of_inertia[0] <= inertia_tol || moments_of_inertia[2] - moments_of_inertia[1] <= inertia_tol) {
        // symmetric, I_A = I_B \ne I_C or I_A \ne I_B = I_C
        // Dnd for n >= 2, (Cn, Cnh, Cnv, Dn, Dnh) for n > 2, Cni for (n > 1 and n is odd, a.k.a. S(2n)) and Sn for 4 | n
        // fmt::print("{:s}\n", "symmetric");

        // rotate the principal axis corresponding to the unequivalent moment of inertia to x axis
        if (moments_of_inertia[1] - moments_of_inertia[0] <= inertia_tol) {
            Eigen::VectorXd swp = - coords_centered.row(coord_x);
            coords_centered.row(coord_x) = coords_centered.row(coord_z);
            coords_centered.row(coord_z) = swp;
        }

        // detect the main axis n > 2
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

        std::function<void(int)> rotate_around_x_by_n = [&coords_centered, &coords_operated](int n) {
            const double angle = 2. * M_PI / static_cast<double>(n);
            const double cos_theta = std::cos(angle);
            const double sin_theta = std::sin(angle);
            Eigen::Matrix2d rot_mat;
            rot_mat << cos_theta, - sin_theta, sin_theta, cos_theta;
            coords_operated.row(coord_x) = coords_centered.row(coord_x);
            coords_operated.bottomRows(2).noalias() = rot_mat * coords_centered.bottomRows(2);
        };

        // find the major Cn
        int major_Cn;
        for (major_Cn = max_Cn_try; major_Cn > 1; -- major_Cn) {
            rotate_around_x_by_n(major_Cn);
            if (is_sym_okay()) break;
        }
        if (major_Cn == 1) throw std::runtime_error("Error: this should never happen.");

        // find available C2
        bool has_minor_C2 = false;
        Eigen::Vector2d axis_point;
        double axis_point_norm;
        Eigen::MatrixXd projection;
        coords_operated.row(coord_x) = - coords_centered.row(coord_x);
        // first, check centers of two SEAs
        for (const std::vector<int>& SEA_group : SEAs) {
            if (SEA_group.size() < 2) continue;
            for (int iatom : SEA_group) {
                for (int jatom : SEA_group) {
                    if (iatom >= jatom) continue;
                    axis_point = (coords_centered.col(iatom).tail(2) + coords_centered.col(jatom).tail(2)) / 2.;
                    axis_point_norm = axis_point.norm();
                    if (axis_point_norm <= tol) continue;
                    axis_point /= axis_point_norm;
                    // test a C2 through origin point and center of iatom and jatom
                    projection = axis_point * (axis_point.transpose() * coords_centered.bottomRows(2));
                    coords_operated.bottomRows(2) = 2. * projection - coords_centered.bottomRows(2);
                    if (is_sym_okay()) {
                        has_minor_C2 = true;
                        break;
                    }
                }
                if (has_minor_C2) break;
            }
            if (has_minor_C2) break;
        }

        if (!has_minor_C2) {
            // second, check C2 through each atom
            for (int iatom = 0; iatom < natoms; ++ iatom) {
                axis_point_norm = coords_centered.col(iatom).tail(2).norm();
                if (axis_point_norm <= tol) continue;
                axis_point = coords_centered.col(iatom).tail(2) / axis_point_norm;
                // test a C2 through origin point and iatom
                projection = axis_point * (axis_point.transpose() * coords_centered.bottomRows(2));
                coords_operated.bottomRows(2) = 2. * projection - coords_centered.bottomRows(2);
                if (is_sym_okay()) {
                    has_minor_C2 = true;
                    break;
                }
            }
        }

        if (has_minor_C2) {
            Eigen::Matrix2d rot_mat;
            rot_mat << axis_point[0], axis_point[1], - axis_point[1], axis_point[0];
            coords_centered.bottomRows(2) = rot_mat * coords_centered.bottomRows(2);
            new_coordinates = coords_centered;
        }

        // find sigma_h
        coords_operated.row(coord_x) = - coords_centered.row(coord_x);
        coords_operated.bottomRows(2) = coords_centered.bottomRows(2);
        bool has_sigma_h = is_sym_okay();

        if (has_minor_C2) {
            // Dn (n > 2), Dnh (n > 2), Dnd (n >= 2)
            if (major_Cn == 2) return "D2d";
            if (has_sigma_h) return fmt::format("D{:d}h", major_Cn);
            // Dn, Dnd for n > 2
            // if exists sigma_d, it divides two minor C2. one minor C2 is already on axis y.
            coords_operated.row(coord_x) = coords_centered.row(coord_x);
            double angle = M_PI / static_cast<double>(major_Cn) / 2.;
            axis_point << std::cos(angle), std::sin(angle);
            projection = axis_point * (axis_point.transpose() * coords_centered.bottomRows(2));
            coords_operated.bottomRows(2) = 2. * projection - coords_centered.bottomRows(2);
            bool has_sigma_d = is_sym_okay();
            return has_sigma_d ? fmt::format("D{:d}d", major_Cn) : fmt::format("D{:d}", major_Cn);
        }

        // (Cn, Cnv, Cnh) for n > 3, Cni for odd i > 1, S4n for positive n
        coords_operated = - coords_centered;
        bool has_sym_center = is_sym_okay();
        // S{4n+2} = C{2n+1} + i (C{2n+1}i), S{2n+1} = C{2n+1} + sigma_h (C{2n+1}h)
        // if there is S4n, there must be C2n, and S4n does not contain i or sigma_h
        if (major_Cn % 2 != 0) {
            if (has_sym_center and has_sigma_h) throw std::runtime_error("Error: this should never happen.");
            if (has_sym_center) return fmt::format("C{:d}i", major_Cn);
            if (has_sigma_h) return fmt::format("C{:d}h", major_Cn); // only odd Cnh here
        } else {
            if (has_sigma_h) {
                if (!has_sym_center) throw std::runtime_error("Error: this should never happen.");
                return fmt::format("C{:d}h", major_Cn); // only even Cnh here
            }
            // if major_Cn is n, then the maximum S, if any, must be S{2n} or Sn.
            // if major_Cn is even, then if the maximum S is Sn, then:
            // 1. if n = 4k, then S{4k} only has C{2k}, this is a paradox.
            // 2. if n = 4k+2, then S{4k+2} = C{2k+1} + i, however C{4k+2} has C2, and C2 + i generates sigma_h, 
            // but C{even} + sigma_h has been discussed above, hence here we only need to check S{2n}, and since 
            // n is even here, that is a S{4m} point group if there is S{2n}.
            const int S_order = major_Cn * 2;
            rotate_around_x_by_n(S_order);
            coords_operated.row(coord_x) = - coords_centered.row(coord_x);
            bool has_Sn = is_sym_okay();
            if (has_Sn) {
                if (S_order % 4 != 0) throw std::runtime_error("Error: this should never happen.");
                return fmt::format("S{:d}", S_order);
            }
        }

        // Cn or Cnv for n > 3
        // find available sigma v
        bool has_sigma_v = false;
        Eigen::Vector2d mirror_point;
        double mirror_point_norm;
        coords_operated.row(coord_x) = coords_centered.row(coord_x);
        // first, check centers of two SEAs
        for (const std::vector<int>& SEA_group : SEAs) {
            if (SEA_group.size() < 2) continue;
            for (int iatom : SEA_group) {
                for (int jatom : SEA_group) {
                    if (iatom >= jatom) continue;
                    mirror_point = (coords_centered.col(iatom).tail(2) + coords_centered.col(jatom).tail(2)) / 2.;
                    mirror_point_norm = mirror_point.norm();
                    if (mirror_point_norm <= tol) continue;
                    mirror_point /= mirror_point_norm;
                    // test a C2 through origin point and center of iatom and jatom
                    projection = mirror_point * (mirror_point.transpose() * coords_centered.bottomRows(2));
                    coords_operated.bottomRows(2) = 2. * projection - coords_centered.bottomRows(2);
                    if (is_sym_okay()) {
                        has_sigma_v = true;
                        break;
                    }
                }
                if (has_sigma_v) break;
            }
            if (has_sigma_v) break;
        }

        if (!has_sigma_v) {
            // second, check C2 through each atom
            for (int iatom = 0; iatom < natoms; ++ iatom) {
                mirror_point_norm = coords_centered.col(iatom).tail(2).norm();
                if (mirror_point_norm <= tol) continue;
                mirror_point = coords_centered.col(iatom).tail(2) / mirror_point_norm;
                // test a C2 through origin point and iatom
                projection = mirror_point * (mirror_point.transpose() * coords_centered.bottomRows(2));
                coords_operated.bottomRows(2) = 2. * projection - coords_centered.bottomRows(2);
                if (is_sym_okay()) {
                    has_sigma_v = true;
                    break;
                }
            }
        }

        if (has_sigma_v) {
            // rotate the found sigma_v to y axis
            Eigen::Matrix2d rot_mat;
            rot_mat << mirror_point[0], mirror_point[1], - mirror_point[1], mirror_point[0];
            coords_centered.bottomRows(2) = rot_mat * coords_centered.bottomRows(2);
            new_coordinates = coords_centered;
            return fmt::format("C{:d}v", major_Cn);
        } else {
            return fmt::format("C{:d}", major_Cn);
        }

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
