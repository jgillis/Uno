// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SUBPROBLEM_H
#define UNO_SUBPROBLEM_H

#include <vector>
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "reformulation/l1RelaxedProblem.hpp"
#include "Direction.hpp"
#include "linear_algebra/Vector.hpp"
#include "tools/Statistics.hpp"

/*! \class Subproblem
 * \brief Subproblem
 */
class Subproblem {
public:
   Subproblem(size_t max_number_variables, size_t max_number_constraints);
   virtual ~Subproblem() = default;

   // virtual methods implemented by subclasses
   virtual void generate_initial_iterate(const NonlinearProblem& problem, Iterate& initial_iterate) = 0;
   virtual Direction solve(Statistics& statistics, const NonlinearProblem& problem, Iterate& current_iterate,
         const WarmstartInformation& warmstart_information) = 0;

   void set_trust_region_radius(double new_trust_region_radius);
   virtual void initialize_feasibility_problem() = 0;
   virtual void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) = 0;
   virtual void exit_feasibility_problem(const NonlinearProblem& problem, Iterate& trial_iterate) = 0;

   // progress measures
   [[nodiscard]] virtual std::function<double(double)> compute_predicted_optimality_reduction_model(const NonlinearProblem& problem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const = 0;
   virtual void set_auxiliary_measure(const NonlinearProblem& problem, Iterate& iterate) = 0;
   [[nodiscard]] virtual double compute_predicted_auxiliary_reduction_model(const NonlinearProblem& problem,
         const Iterate& current_iterate, const Direction& direction, double step_length) const = 0;

   virtual void postprocess_iterate(const NonlinearProblem& problem, Iterate& iterate) = 0;

   [[nodiscard]] virtual size_t get_hessian_evaluation_count() const = 0;
   virtual void set_initial_point(const std::vector<double>& initial_point) = 0;

   size_t number_subproblems_solved{0};
   // when the parameterization of the subproblem (e.g. penalty or barrier parameter) is updated, signal it
   bool subproblem_definition_changed{false};

protected:
   Direction direction;
   Evaluations evaluations;
   double trust_region_radius{INF<double>};
};

#endif // UNO_SUBPROBLEM_H
