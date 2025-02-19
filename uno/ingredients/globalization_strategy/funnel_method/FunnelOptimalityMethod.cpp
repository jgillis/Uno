// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "FunnelOptimalityMethod.hpp"

FunnelOptimalityMethod::FunnelOptimalityMethod(Statistics& statistics, const Options& options) :
      FunnelMethod(statistics, options)
     { }

/* check acceptability of step(s) (funnel & sufficient reduction)
 * funnel methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */

void FunnelOptimalityMethod::update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure) {

   this->funnel_width = std::max(this->parameters.kappa_infeasibility_1 *this->funnel_width,
      trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (funnel_width - trial_infeasibility_measure));

   DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
   
}

//! check acceptability wrt current point
bool FunnelOptimalityMethod::acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
      double trial_optimality_measure) {
   return (trial_optimality_measure <= current_optimality_measure - this->parameters.gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->parameters.beta * current_infeasibility_measure);
}

bool FunnelOptimalityMethod::is_iterate_acceptable(Statistics& statistics, const Iterate& /*trial_iterate*/,
      const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction,
      double /*objective_multiplier*/) {
   const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   

   const double current_infeasibility_measure = current_progress_measures.infeasibility;
   const double trial_infeasibility_measure = trial_progress_measures.infeasibility;


   // unconstrained predicted reduction:
   // - ignore the predicted infeasibility reduction
   // - scale the scaled optimality measure with 1
   const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   const double unconstrained_predicted_reduction_infeasibility = predicted_reduction.infeasibility;

   DEBUG << "\t\tCurrent: η = " << current_progress_measures.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   DEBUG << "\t\tTrial:   η = " << trial_progress_measures.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   DEBUG << "\t\tUnconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
         " = " <<  unconstrained_predicted_reduction << '\n';
   DEBUG << "\t\tUnconstrained predicted infeasibility reduction: " << unconstrained_predicted_reduction_infeasibility << '\n';

   statistics.add_statistic("funnel width", this->get_funnel_width());
   
   DEBUG << "\t\t" <<*this << '\n';

   bool accept = false;
   bool funnel_reduction_mechanism = false;
   bool funnel_acceptable = false;

   funnel_acceptable = this->is_infeasibility_acceptable_to_funnel(trial_progress_measures.infeasibility);
   if (funnel_acceptable) {
      if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
         DEBUG << "\t\tTrial iterate satisfies switching condition ....\n";
         // check acceptance 
         const double actual_reduction = this->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
               trial_optimality_measure);
         DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';
         
         // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
         if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
            DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
            accept = true;
            // decrease funnel here??? ......
         }
         else { // switching condition holds, but not Armijo condition
            DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
         }
      } else {
         DEBUG << "\t\tTrial iterate ACCEPTED by violating the switching condition ...\n";
         funnel_reduction_mechanism = true;
         accept = true; // accept the step and reduce the tr-radius
      }
   } else {
      DEBUG << "\t\tTrial iterate REJECTED by violating Funnel condition\n";
      // acceptable = false;
      // funnel_reduction_mechanism = false;
   }

   // if (this->current_iterate_acceptable_to_funnel){
   //    DEBUG << "\t\tCurrent iterate IN Funnel\n";
   //    // switching condition: the unconstrained predicted reduction is sufficiently positive
   //    if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
   //       DEBUG << "\t\tTrial iterate satisfies switching condition ....\n";
   //       // check acceptance 
   //       const double actual_reduction = this->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
   //             trial_optimality_measure);
   //       DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';
         
   //       // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
   //       if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
   //          DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
   //          accept = true;
   //          // decrease funnel here??? ......
   //       }
   //       else { // switching condition holds, but not Armijo condition
   //          DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
   //       }
   //    } else {
   //       DEBUG << "\t\tTrial iterate violates switching condition ...\n";
   //       funnel_acceptable = this->is_infeasibility_acceptable_to_funnel(trial_progress_measures.infeasibility);
   //       if (funnel_acceptable) {
   //          DEBUG << "\t\tFunnel condition acceptable\n";
   //          funnel_reduction_mechanism = true;
   //          accept = true; // accept the step and reduce the tr-radius
   //       } else {
   //          DEBUG << "\t\tFunnel condition NOT acceptable\n";
   //       }
   //    }
   // } else {
   //    // step is rejected --> this is very much like to tolerance-tube approach
   //    DEBUG << "\t\tCurrent iterate NOT in Funnel\n";

   //    const double actual_reduction = this->compute_actual_reduction(current_infeasibility_measure, current_progress_measures.infeasibility,
   //       trial_infeasibility_measure);
   //    DEBUG << "\t\tActual reduction feasibility: " << actual_reduction << '\n';
   //    if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction_infeasibility, actual_reduction)) {
   //       DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition for INFEASIBILITY\n";
   //       accept = true;
   //    } else {
   //       DEBUG << "\t\tTrial iterate was NOT ACCEPTED by violating Armijo condition for INFEASIBILITY\n";

   //    }
   // }
   // -------------------------------------------------------------------------

   if (funnel_reduction_mechanism){ // steps needs to be accepted for this...
       DEBUG << "\t\tEntering funnel reduction mechanism\n";

      // // Feasibility measures
      this->update_funnel_width(current_infeasibility_measure, 
                                             trial_infeasibility_measure);
      // this->funnel_width = this->get_funnel_width(); //?                                            ;
   }

   // if (accept){
   //    if (this->is_infeasibility_acceptable_to_funnel(trial_infeasibility_measure)){
   //       this->current_iterate_acceptable_to_funnel = true;
   //       DEBUG << "New Iterate inside of FUNNEL\n";
   //    }
   //    else {
   //       this->current_iterate_acceptable_to_funnel = false;
   //    }
   // }

   return accept;
}
