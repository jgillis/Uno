// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_GLOBALIZATIONMECHANISM_H
#define UNO_GLOBALIZATIONMECHANISM_H

#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategy.hpp"
#include "optimization/Model.hpp"
#include "preprocessing/Scaling.hpp"
#include "tools/Statistics.hpp"

class GlobalizationMechanism {
public:
   GlobalizationMechanism(ConstraintRelaxationStrategy& constraint_relaxation_strategy, const Options& options);
   virtual ~GlobalizationMechanism() = default;

   virtual void initialize(Iterate& initial_iterate) = 0;
   virtual Iterate compute_next_iterate(Statistics& statistics, const Model& model, Iterate& current_iterate) = 0;

   [[nodiscard]] size_t get_hessian_evaluation_count() const;
   [[nodiscard]] size_t get_number_subproblems_solved() const;

protected:
   // reference to allow polymorphism
   ConstraintRelaxationStrategy& constraint_relaxation_strategy; /*!< Constraint relaxation strategy */
   const double tight_tolerance; /*!< Tight tolerance of the termination criteria */
   const double loose_tolerance; /*!< Loose tolerance of the termination criteria */
   size_t loose_tolerance_consecutive_iterations{0};
   const size_t loose_tolerance_consecutive_iteration_threshold;
   const Norm progress_norm;
   const double unbounded_objective_threshold;

   static Iterate assemble_trial_iterate(Iterate& current_iterate, const Direction& direction, double primal_step_length,
         double dual_step_length, double bound_dual_step_length);
   bool check_termination_with_small_step(const Model& model, const Direction& direction, Iterate& trial_iterate) const;
   [[nodiscard]] TerminationStatus check_convergence(const Model& model, Iterate& current_iterate);
   [[nodiscard]] TerminationStatus check_convergence(const Model& model, Iterate& current_iterate, double tolerance) const;
};

#endif // UNO_GLOBALIZATIONMECHANISM_H
