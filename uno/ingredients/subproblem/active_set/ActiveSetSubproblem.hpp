// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_ACTIVESETSUBPROBLEM_H
#define UNO_ACTIVESETSUBPROBLEM_H

#include "ingredients/subproblem/Subproblem.hpp"

class ActiveSetSubproblem : public Subproblem {
public:
   ActiveSetSubproblem(size_t max_number_variables, size_t max_number_constraints);
   ~ActiveSetSubproblem() override = default;

   void initialize(Statistics& statistics, const NonlinearProblem& problem, Iterate& first_iterate) override;
   void set_initial_point(const std::vector<double>& initial_point) override;
   void initialize_feasibility_problem() override;
   void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;

   void set_unscaled_optimality_measure(const NonlinearProblem& problem, Iterate& iterate) override;
   [[nodiscard]] double generate_predicted_unscaled_optimality_reduction_model(const NonlinearProblem& problem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const override;

   void postprocess_accepted_iterate(const NonlinearProblem& model, Iterate& iterate) override;

protected:
   std::vector<double> initial_point{};
   std::vector<Interval> variable_displacement_bounds{};
   std::vector<Interval> linearized_constraint_bounds{};

   void set_variable_displacement_bounds(const NonlinearProblem& problem, const Iterate& current_iterate);
   void set_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& current_constraints);
   static void compute_dual_displacements(const NonlinearProblem& problem, const Iterate& current_iterate, Direction& direction);
   void shift_linearized_constraint_bounds(const NonlinearProblem& problem, const std::vector<double>& trial_constraints);
};

#endif // UNO_ACTIVESETSUBPROBLEM_H