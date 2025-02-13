// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "BarrierParameterUpdateStrategy.hpp"
#include "linear_algebra/VectorExpression.hpp"
#include "tools/Logger.hpp"

BarrierParameterUpdateStrategy::BarrierParameterUpdateStrategy(const Options& options):
   barrier_parameter(options.get_double("barrier_initial_parameter")),
   tolerance(options.get_double("tolerance")),
   parameters({
      options.get_double("barrier_k_mu"),
      options.get_double("barrier_theta_mu"),
      options.get_double("barrier_k_epsilon"),
      options.get_double("barrier_update_fraction")
   }) {
}

double BarrierParameterUpdateStrategy::get_barrier_parameter() const {
   return this->barrier_parameter;
}

void BarrierParameterUpdateStrategy::set_barrier_parameter(double new_barrier_parameter) {
   assert(0. <= new_barrier_parameter && "The barrier parameter should be positive.");
   this->barrier_parameter = new_barrier_parameter;
}

bool BarrierParameterUpdateStrategy::update_barrier_parameter(const NonlinearProblem& problem, const Iterate& current_iterate) {
   // primal-dual errors
   const double scaled_stationarity = current_iterate.residuals.optimality_stationarity/current_iterate.residuals.stationarity_scaling;
   double primal_dual_error = std::max({
      scaled_stationarity,
      current_iterate.residuals.infeasibility,
      current_iterate.residuals.optimality_complementarity / current_iterate.residuals.complementarity_scaling
   });
   DEBUG << "Max scaled primal-dual error for barrier subproblem is " << primal_dual_error << '\n';

   // update the barrier parameter (Eq. 7 in IPOPT paper)
   const double tolerance_fraction = this->tolerance / this->parameters.update_fraction;
   bool parameter_updated = false;
   while (primal_dual_error <= this->parameters.k_epsilon * this->barrier_parameter && tolerance_fraction < this->barrier_parameter) {
      this->barrier_parameter = std::max(tolerance_fraction, std::min(this->parameters.k_mu * this->barrier_parameter,
            std::pow(this->barrier_parameter, this->parameters.theta_mu)));
      DEBUG << "Barrier parameter mu updated to " << this->barrier_parameter << '\n';
      // update complementarity error
      double scaled_complementarity_error = BarrierParameterUpdateStrategy::compute_shifted_complementarity_error(problem, current_iterate,
            this->barrier_parameter) / current_iterate.residuals.complementarity_scaling;
      primal_dual_error = std::max({
         scaled_stationarity,
         current_iterate.residuals.infeasibility,
         scaled_complementarity_error
      });
      DEBUG << "Max scaled primal-dual error for barrier subproblem is " << primal_dual_error << '\n';
      parameter_updated = true;
   }
   return parameter_updated;
}

double BarrierParameterUpdateStrategy::compute_shifted_complementarity_error(const NonlinearProblem& problem, const Iterate& iterate,
      double shift_value) {
   VectorExpression<double> shifted_bound_complementarity(problem.number_variables, [&](size_t i) {
      double result = 0.;
      if (0. < iterate.multipliers.lower_bounds[i]) { // lower bound
         result = std::max(result, std::abs(iterate.multipliers.lower_bounds[i] * (iterate.primals[i] - problem.get_variable_lower_bound(i)) - shift_value));
      }
      if (iterate.multipliers.upper_bounds[i] < 0.) { // upper bound
         result = std::max(result, std::abs(iterate.multipliers.upper_bounds[i] * (iterate.primals[i] - problem.get_variable_upper_bound(i)) - shift_value));
      }
      return result;
   });
   return norm_inf(shifted_bound_complementarity); // TODO use a generic norm
}
