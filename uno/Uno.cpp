// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "Uno.hpp"
#include "ingredients/constraint_relaxation_strategy/ConstraintRelaxationStrategyFactory.hpp"
#include "ingredients/globalization_mechanism/GlobalizationMechanismFactory.hpp"
#include "ingredients/globalization_strategy/GlobalizationStrategyFactory.hpp"
#include "ingredients/subproblem/SubproblemFactory.hpp"
#include "optimization/Iterate.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"
#include "tools/Timer.hpp"

Uno::Uno(GlobalizationMechanism& globalization_mechanism, const Options& options) :
      globalization_mechanism(globalization_mechanism),
      max_iterations(options.get_unsigned_int("max_iterations")),
      time_limit(options.get_double("time_limit")) {
}

Result Uno::solve(Statistics& statistics, const Model& model, Iterate& current_iterate) {
   Timer timer{};
   size_t major_iterations = 0;

   std::cout << "\nProblem " << model.name << '\n';
   std::cout << model.number_variables << " variables, " << model.number_constraints << " constraints\n\n";

   // use the current point to initialize the strategies and generate the initial iterate
   try {
      this->globalization_mechanism.initialize(current_iterate);
   }
   catch (const std::exception& e) {
      ERROR << RED << "An error occurred at the initial iterate: " << e.what() << RESET;
      throw;
   }

   bool termination = false;
   try {
      // check for termination
      while (not termination) {
         statistics.new_line();
         major_iterations++;
         DEBUG << "### Outer iteration " << major_iterations << '\n';

         // compute an acceptable iterate by solving a subproblem at the current point
         current_iterate = this->globalization_mechanism.compute_next_iterate(statistics, model, current_iterate);

         // compute the status of the next iterate
         Uno::add_statistics(statistics, current_iterate, major_iterations);
         if (Logger::level == INFO) statistics.print_current_line();

         statistics.add_iteration();
         // add iteration to map, if planned to serialize
         // if (statistics.serialize_iterations == true){
         // }

         termination = this->termination_criteria(current_iterate.status, major_iterations, timer.get_duration());
      }
   }
   catch (const std::runtime_error& e) {
      ERROR << RED << e.what() << RESET;
      throw;
   }
   catch (std::exception& exception) {
      ERROR << RED << exception.what() << RESET;
   }
   statistics.serialize();
   Uno::postprocess_iterate(model, current_iterate, current_iterate.status);

   if (Logger::level == INFO) statistics.print_footer();

   const size_t number_subproblems_solved = this->globalization_mechanism.get_number_subproblems_solved();
   const size_t hessian_evaluation_count = this->globalization_mechanism.get_hessian_evaluation_count();
   Result result = {std::move(current_iterate), model.number_variables, model.number_constraints, major_iterations, timer.get_duration(),
         Iterate::number_eval_objective, Iterate::number_eval_constraints, Iterate::number_eval_objective_gradient,
         Iterate::number_eval_jacobian, hessian_evaluation_count, number_subproblems_solved};
   return result;
}

void Uno::add_statistics(Statistics& statistics, const Iterate& iterate, size_t major_iterations) {
   statistics.add_statistic(std::string("iters"), major_iterations);
   if (iterate.is_objective_computed) {
      statistics.add_statistic("objective", iterate.evaluations.objective);
   }
   else {
      statistics.add_statistic("objective", "-");
   }
}

bool Uno::termination_criteria(TerminationStatus current_status, size_t iteration, double current_time) const {
   return current_status != TerminationStatus::NOT_OPTIMAL || this->max_iterations <= iteration || this->time_limit <= current_time;
}

void Uno::postprocess_iterate(const Model& model, Iterate& iterate, TerminationStatus termination_status) {
   // in case the objective was not yet evaluated, evaluate it
   iterate.evaluate_objective(model);
   model.postprocess_solution(iterate, termination_status);
   DEBUG2 << "Final iterate:\n" << iterate;
}

void join(const std::vector<std::string>& vector, char separator) {
   if (not vector.empty()) {
      std::cout << vector[0];
      for (size_t i: Range(1, vector.size())) {
         std::cout << separator << ' ' << vector[i];
      }
   }
}

void Uno::print_available_strategies() {
   std::cout << "Available strategies:\n";
   std::cout << "Constraint relaxation strategies: ";
   join(ConstraintRelaxationStrategyFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Globalization mechanisms: ";
   join(GlobalizationMechanismFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Globalization strategies: ";
   join(GlobalizationStrategyFactory::available_strategies(), ',');
   std::cout << '\n';
   std::cout << "Subproblems: ";
   join(SubproblemFactory::available_strategies(), ',');
   std::cout << '\n';
}