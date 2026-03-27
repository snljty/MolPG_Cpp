#include "MolSymm.hpp"

PeriodicTable& PeriodicTable::get_instance() {
    static PeriodicTable instance;
    return instance;
}

PeriodicTable::PeriodicTable() {
    for (int i = 1; i <= max_elements_num; ++ i) {
        elements_table.insert({elements_list[i], i});
    }
}

int PeriodicTable::get_atomic_number(std::string_view element_name) const {
    std::unordered_map<std::string_view, int>::const_iterator it = elements_table.find(element_name);
    return (it != elements_table.end()) ? it->second : 0;
}

double PeriodicTable::get_atomic_weight(int atomic_index) const {
    if (atomic_index > max_elements_num) throw std::invalid_argument(fmt::format("Error: no element with atomic index {:d}.", atomic_index));
    return elements_average_weight[atomic_index];
}
