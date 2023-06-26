// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "FunnelRestorationMethod.hpp"
// #include "funnel/FunnelFactory.hpp"

FunnelRestorationMethod::FunnelRestorationMethod(Statistics& statistics, const Options& options) :
      FunnelMethod(statistics, options)
   {}


void FunnelRestorationMethod::update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure) {

   if (trial_infeasibility_measure <= this->funnel_width){
      if (current_infeasibility_measure > this->funnel_width){
         this->funnel_width = std::min(this->parameters.kappa_infeasibility_1 *this->funnel_width,
         trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (funnel_width - trial_infeasibility_measure));
         // std::cout << "Funnel Update Restoration: Current iterate outside funnel" << std::endl;
      } else {
         this->funnel_width = std::max(this->parameters.kappa_infeasibility_1 *this->funnel_width, 
         trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure));
         // std::cout << "Funnel Update Restoration: Current iterate inside funnel" << std::endl;
      }
   } //else: do not reduce the funnel

   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
   
}


/* check acceptability of step(s) (funnel & sufficient reduction)
 * funnel methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool FunnelRestorationMethod::is_iterate_acceptable(Statistics& statistics, const Iterate& /*trial_iterate*/,
      const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction,
      double /*objective_multiplier*/) {

   const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   
   // // // Feasibility measures
   // const double current_infeasibility_measure = current_progress_measures.infeasibility;
   // const double trial_infeasibility_measure = trial_progress_measures.infeasibility;

   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   DEBUG << "\t\tCurrent: η = " << current_progress_measures.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   DEBUG << "\t\tTrial:   η = " << trial_progress_measures.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   DEBUG << "\t\tUnconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
         " = " <<  unconstrained_predicted_reduction << '\n';

   statistics.add_statistic("funnel width", this->get_funnel_width());
   
   //Output funnel parameter
   DEBUG << "\t\t" << *this << '\n';

   bool accept = false;
   bool funnel_acceptable = false;
   
   const double actual_reduction = this->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
         trial_optimality_measure);
   DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';

   // switching condition: the unconstrained predicted reduction is sufficiently positive
   if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
      DEBUG << "\t\tTrial iterate satisfies switching condition ....\n";
      // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
      if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
         DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
         accept = true;
      }
      else { // switching condition holds, but not Armijo condition
         DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
      }
   }
   else {
      DEBUG << "\t\tTrial iterate violates switching condition ...THIS SHOULD NOT HAPPEN!\n";
   }

   // We are in restoration phase, i.e., optimality measure contains the infeasibility!
   if (accept){
      if (this->is_infeasibility_acceptable_to_funnel(trial_optimality_measure)){
         this->current_iterate_acceptable_to_funnel = true;
      }
      else {
         this->current_iterate_acceptable_to_funnel = false;
      }



      funnel_acceptable = this->is_infeasibility_acceptable_to_funnel(trial_optimality_measure);
      // check acceptance   
      if (funnel_acceptable) {
         DEBUG << "\t\tFunnel condition satisfied\n";
         
         // do update 1
         this->update_funnel_width(current_optimality_measure, 
                                    trial_optimality_measure);
         this->funnel_width = this->get_funnel_width(); 
         this->current_iterate_acceptable_to_funnel = true;
         DEBUG << "\t\t New Iterate inside of FUNNEL\n";
      }
      else {
         this->current_iterate_acceptable_to_funnel = false;
      }
      // } else {

      //    // do update 2
      //    DEBUG << "\t\tFunnel condition NOT satisfied\n";
      //    // this->funnel->update_funnel_parameter(current_optimality_measure, 
      //    //                                        trial_optimality_measure);
      //    this->funnel_width = this->funnel->get_funnel_size();   
      // }
      // Funnel has changed, check, if still acceptable
   }

   return accept;
}
