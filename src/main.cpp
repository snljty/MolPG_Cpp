#include "MolSymm.hpp"

int main(int argc, const char *argv[]) {
    if (argc - 1 != 1 && argc -1 != 2) throw std::invalid_argument(fmt::format("Usage: {:s} xxx.xyz [tolerance]", argv[0]));
    std::string molname(argv[1]);
    Molecule mol(molname);
    std::string result;
    if (argc - 1 == 2) {
        double tol = std::stod(std::string(argv[2]));
        result = mol.detect_point_group(tol);
    } else {
        result = mol.detect_point_group();
    }
    fmt::print("{:s}\n", result);
    // mol.use_new_coordinates();
    // mol.write_gjf("new.gjf");

    return 0;
}

