#include "MolSymm.hpp"

Molecule::Molecule(const std::string& ifilename) : periodic_table(PeriodicTable::get_instance()) {
    if (!ifilename.empty()) read(ifilename);
}

void Molecule::read(const std::string& ifilename) {
    std::string::size_type pos = ifilename.rfind('.');
    if (pos == std::string::npos) throw std::invalid_argument("Error: the suffix cannot be determined.");
    if (ifilename.substr(pos) == ".xyz") {
        return read_xyz(ifilename);
    } else {
        throw std::invalid_argument(fmt::format("Error: cannot understand the suffix \"{:s}\".", ifilename.substr(pos)));
    }
}

void Molecule::resize() {
    // note: this method does not change this->natoms .
    elements.resize(natoms);
    coordinates.resize(ncoords, natoms);
    atomic_numbers.resize(natoms);
    atomic_weights.resize(natoms);
}

void Molecule::read_xyz(const std::string& ifilename) {
    std::string::size_type pos = ifilename.rfind('.');
    if (pos == std::string::npos || ifilename.substr(pos) != ".xyz") throw std::invalid_argument("Error: the suffix must be \".xyz\".");

    std::ifstream ifile(ifilename);
    if (!ifile.is_open()) throw std::ios_base::failure(fmt::format("Error: cannot open file \"{:s}\" for reading.", ifilename));
    std::string line;

    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read amount of atoms.");
    natoms = std::stoi(line);
    if (natoms <= 0) throw std::invalid_argument("Error: amount of atoms must be positive.");
    resize();

    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read comment line.");

    for (int iatom = 0; iatom < natoms; ++ iatom) {
        if (!std::getline(ifile, line)) throw std::ios_base::failure(fmt::format("Error: cannot read atom {:d}.", iatom + 1));
        std::istringstream iss(line);
        if (!(iss >> elements[iatom] 
                >> coordinates(coord_x, iatom) 
                >> coordinates(coord_y, iatom) 
                >> coordinates(coord_z, iatom))) {
            throw std::ios_base::failure(fmt::format("Error: cannot read atom {:d}.", iatom + 1));
        }
        atomic_numbers[iatom] = periodic_table.get_atomic_number(elements[iatom]);
        atomic_weights[iatom] = periodic_table.get_atomic_weight(atomic_numbers[iatom]);
    }

    ifile.close();
}

