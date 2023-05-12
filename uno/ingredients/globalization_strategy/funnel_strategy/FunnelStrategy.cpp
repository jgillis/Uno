// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "FunnelStrategy.hpp"
#include "funnel/FunnelFactory.hpp"

FunnelStrategy::FunnelStrategy(Statistics& statistics, const Options& options) :
      GlobalizationStrategy(options),
      funnel(FunnelFactory::create(options)),
      parameters({
         options.get_double("funnel_kappa_initial_upper_bound"),
         options.get_double("funnel_kappa_initial_multiplication"),
         options.get_double("funnel_delta"),
         options.get_double("funnel_ubd"),
         options.get_double("funnel_fact"),
         options.get_double("funnel_switching_infeasibility_exponent")
      }) {
   statistics.add_column("funnel width", Statistics::double_width, options.get_int("statistics_funnel_size_column_order"));

}

void FunnelStrategy::initialize(const Iterate& initial_iterate) {
   // set the funnel upper bound
   double upper_bound = std::max(this->parameters.kappa_initial_upper_bound,
                                 this->parameters.kappa_initial_multiplication * initial_iterate.progress.infeasibility);
   this->funnel->initial_upper_bound = upper_bound;
   this->funnel->initialize();
   this->current_phase = 2;
   this->initial_funnel_upper_bound = upper_bound;

}

void FunnelStrategy::reset() {
   // re-initialize the restoration funnel
   this->funnel->reset();
   // TODO: we should set the ub of the optimality funnel. But now, our 2 funnels live independently...
   this->funnel->initial_upper_bound = this->initial_funnel_upper_bound;
}

void FunnelStrategy::register_current_progress(const ProgressMeasures& current_progress_measures) {
   const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   this->funnel->add(current_progress_measures.infeasibility, current_optimality_measure);
}

bool FunnelStrategy::is_infeasibility_acceptable(double infeasibility_measure) const {
   return (this->funnel->acceptable(infeasibility_measure));//not used anywhere
}

bool FunnelStrategy::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
   // return predicted_reduction > 0;
}

/* check acceptability of step(s) (funnel & sufficient reduction)
 * funnel methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool FunnelStrategy::is_iterate_acceptable(Statistics& statistics, const Iterate& /*trial_iterate*/,
      const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction,
      double /*objective_multiplier*/) {
   const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   
   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   DEBUG << "\t\tCurrent: η = " << current_progress_measures.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   DEBUG << "\t\tTrial:   η = " << trial_progress_measures.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   DEBUG << "\t\tUnconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
         " = " <<  unconstrained_predicted_reduction << '\n';

   GlobalizationStrategy::check_finiteness(current_progress_measures, 1.);
   GlobalizationStrategy::check_finiteness(trial_progress_measures, 1.);
   statistics.add_statistic("funnel width", this->funnel->get_funnel_size());
   
   DEBUG << "\t\t" <<*this->funnel << '\n';

   bool accept = false;

   bool funnel_reduction_mechanism = false;
   
   // std::cout << "Phase in funnel mechanism" << static_cast<int>(this->get_phase()) << std::endl;
   // this->current_phase = Phase::FEASIBILITY_RESTORATION;
   DEBUG << "\t\tPhase in funnel mechanism" << static_cast<int>(this->current_phase) << "\n";

   // if (this->current_phase == Phase::OPTIMALITY){
   if (this->current_phase == 2){
      DEBUG  << "\t\tCurrent phase OPTIMALITY\n";
   } else if (this->current_phase == 1) {
      DEBUG  << "\t\tCurrent phase RESTORATION\n";
   } else {
      DEBUG << "WARNING!!: Phase is not in {1,2}!!\n";
   }

   // if (this->current_phase == Phase::OPTIMALITY){
   // if (this->current_phase == 2){
   
      // check acceptance
      const bool funnel_acceptable = this->funnel->acceptable(trial_progress_measures.infeasibility);
      
      if (funnel_acceptable) {
         DEBUG << "\t\tFunnel condition acceptable\n";

         // check acceptance wrt current point
         // const bool improves_current_iterate = this->funnel->acceptable_wrt_current_iterate(current_progress_measures.infeasibility,
         //       current_optimality_measure, trial_progress_measures.infeasibility, trial_optimality_measure);
         // if (improves_current_iterate) {
         //    DEBUG << "Acceptable with respect to current point\n";

            const double actual_reduction = this->funnel->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
                  trial_optimality_measure);
            DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';

            DEBUG << "\t\tCheck the switching condition ....\n";
            // switching condition: the unconstrained predicted reduction is sufficiently positive
            if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
               // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
               DEBUG << "\t\tCheck the armijo condition for descent ....\n";
               if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
                  DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
                  accept = true;
               }
               else { // switching condition holds, but not Armijo condition
                  DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
               }
            }
            else {
               DEBUG << "\t\tTrial iterate violates switching condition ...\n";
               funnel_reduction_mechanism = true;
            }
            // else { // switching condition violated: predicted reduction is not promising
            //    this->funnel->add(current_progress_measures.infeasibility, current_optimality_measure);
            //    DEBUG << "Trial iterate was accepted by violating switching condition\n";
            //    DEBUG << "Current iterate was added to the funnel\n";
            //    accept = true;
            // }

         // }
         // else {
         //    DEBUG << "Not acceptable with respect to current point\n";
         // }
      } else {
         funnel_reduction_mechanism = true;
         DEBUG << "\t\tFunnel condition NOT acceptable\n";
      }

   // }
   // else {
   //    funnel_reduction_mechanism = true;
   // }

   if (funnel_reduction_mechanism){
       DEBUG << "\t\tEntering funnel reduction mechanism\n";

      // Feasibility measures
      const double current_infeasibility_measure = current_progress_measures.infeasibility;
      const double trial_infeasibility_measure = trial_progress_measures.infeasibility;
      const double predicted_infeasibility_reduction = predicted_reduction.infeasibility;
      
      if (predicted_infeasibility_reduction < 0){
         std::cout << "WARNING: predicted infeasibility reduction smaller 0!!" << std::endl;
      }

      const double actual_feasibility_reduction = this->funnel->compute_actual_reduction(current_infeasibility_measure, current_progress_measures.infeasibility,
                  trial_infeasibility_measure);
      
      if (this->armijo_sufficient_decrease(predicted_infeasibility_reduction, actual_feasibility_reduction)) {
         accept = true;
         DEBUG << "\t\tChanging funnel radius\n";
         this->funnel->update_funnel_parameter(current_infeasibility_measure, 
                                               trial_infeasibility_measure);
      } else {
         DEBUG << "\t\tFunnel is not changed and iterate not accepted\n";
      }
   }

   return accept;
}
