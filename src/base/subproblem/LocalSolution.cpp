#include "LocalSolution.hpp"
#include "Logger.hpp"

LocalSolution::LocalSolution(std::vector<double>& x, std::vector<double>& bound_multipliers, std::vector<double>& constraint_multipliers) :
phase_1_required(false), x(x), bound_multipliers(bound_multipliers), constraint_multipliers(constraint_multipliers) {
}

std::ostream& operator<<(std::ostream &stream, LocalSolution& solution) {
    unsigned int max_size = 50;

    if (solution.status == OPTIMAL) {
        stream << GREEN "Status: optimal\n" RESET;
    }
    else if (solution.status == UNBOUNDED_PROBLEM) {
        stream << GREEN "Status: unbounded\n" RESET;
    }
    else {
        stream << RED "Status " << solution.status << ": Beware peasant, something went wrong\n" RESET;
    }

    stream << MAGENTA;
    stream << "d^* = ";
    print_vector(stream, solution.x, 0, max_size);

    stream << "objective = " << solution.objective << "\n";
    stream << "norm = " << solution.norm << "\n";

    stream << "active set at upper bound =";
    for (unsigned int k = 0; k < solution.active_set.at_upper_bound.size(); k++) {
        unsigned int index = solution.active_set.at_upper_bound[k];
        if (index < solution.x.size()) {
            stream << " x" << index;
        }
        else {
            stream << " c" << (index - solution.x.size());
        }
    }
    stream << "\n";

    stream << "active set at lower bound =";
    for (unsigned int k = 0; k < solution.active_set.at_lower_bound.size(); k++) {
        unsigned int index = solution.active_set.at_lower_bound[k];
        if (index < solution.x.size()) {
            stream << " x" << index;
        }
        else {
            stream << " c" << (index - solution.x.size());
        }
    }
    stream << "\n";

    stream << "general feasible =";
    for (unsigned int k = 0; k < solution.constraint_partition.feasible_set.size(); k++) {
        int index = solution.constraint_partition.feasible_set[k];
        stream << " c" << index;
    }
    stream << "\n";

    stream << "general infeasible =";
    for (unsigned int k = 0; k < solution.constraint_partition.infeasible_set.size(); k++) {
        int index = solution.constraint_partition.infeasible_set[k];
        stream << " c" << index;
        if (solution.constraint_partition.constraint_status[index] == INFEASIBLE_LOWER) {
            stream << " (lower)";
        }
        else if (solution.constraint_partition.constraint_status[index] == INFEASIBLE_UPPER) {
            stream << " (upper)";
        }
    }
    stream << "\n";

    stream << "bound multipliers = ";
    print_vector(stream, solution.bound_multipliers);

    stream << "constraint multipliers = ";
    print_vector(stream, solution.constraint_multipliers);

    stream << RESET;

    return stream;
}
