#include "MolSymm.hpp"

int main(int argc, const char *argv[]) {
    if (argc - 1 != 1) throw std::invalid_argument(fmt::format("Usage: {:s} xxx.xyz", argv[0]));
    Molecule mol(argv[1]);
    fmt::print("{:s}\n", mol.detect_point_group());
    // mol.use_new_coordinates();
    // mol.write_gjf("new.gjf");

    return 0;
}

