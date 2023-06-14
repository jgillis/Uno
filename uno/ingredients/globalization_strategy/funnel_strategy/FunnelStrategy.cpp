// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include "FunnelStrategy.hpp"
// #include "funnel/FunnelFactory.hpp"

FunnelStrategy::FunnelStrategy(Statistics& statistics, const Options& options) :
      GlobalizationStrategy(options),
      parameters({
         options.get_double("funnel_kappa_initial_upper_bound"),
         options.get_double("funnel_kappa_initial_multiplication"),
         options.get_double("funnel_delta"),
         options.get_double("funnel_ubd"),
         options.get_double("funnel_fact"),
         options.get_double("funnel_switching_infeasibility_exponent"),
         options.get_double("funnel_kappa_infeasibility_1"),
         options.get_double("funnel_kappa_infeasibility_2"),
         options.get_double("funnel_beta"),
         options.get_double("funnel_gamma")
      }) {
   statistics.add_column("funnel width", Statistics::double_width, options.get_int("statistics_funnel_size_column_order"));

}

void FunnelStrategy::initialize(const Iterate& initial_iterate) {
   // set the funnel upper bound
   double upper_bound = std::max(this->parameters.kappa_initial_upper_bound,
                                 this->parameters.kappa_initial_multiplication * initial_iterate.progress.infeasibility);

   // std::cout << "Initial kappa upper bound: " << upper_bound << std::endl;
   // std::cout << "Initial infeasibility: " << initial_iterate.progress.infeasibility << std::endl;
   // std::cout << "Initial funnel upper bound: " << upper_bound << std::endl;
   this->initial_funnel_upper_bound = upper_bound;
   this->funnel_width = initial_funnel_upper_bound;
   this->current_iterate_acceptable_to_funnel = true;
   this->initial_funnel_upper_bound = upper_bound;

}

void FunnelStrategy::reset() {
   // re-initialize the restoration funnel
   // this->funnel->reset();
   // this->funnel->initial_upper_bound = this->initial_funnel_upper_bound;
}

void FunnelStrategy::register_current_progress(const ProgressMeasures& /*current_progress_measures*/) {
   // const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   // this->funnel->add(current_progress_measures.infeasibility, current_optimality_measure);
}

bool FunnelStrategy::is_infeasibility_acceptable_to_funnel(double infeasibility_measure) const {
   if (infeasibility_measure <= this->funnel_width){
      return true;
   }
   else {
      DEBUG << "\t\tNot acceptable to funnel.\n";
      return false;
   }
}

bool FunnelStrategy::is_infeasibility_acceptable(double infeasibility_measure) const {
   return this->is_infeasibility_acceptable_to_funnel(infeasibility_measure);
}

bool FunnelStrategy::switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const {
   return predicted_reduction > switching_fraction * std::pow(current_infeasibility, this->parameters.switching_infeasibility_exponent);
}

void FunnelStrategy::update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure) {

   // this->funnel_width = std::max(this->parameters.kappa_infeasibility_1 *this->funnel_width, 
   //    trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure));

   // DEBUG << "\t\tNew funnel parameter is: " << this->funnel_width << "\n"; 
   
}

//! check acceptability wrt current point, Maybe useful later
// bool Funnel::acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
//       double trial_optimality_measure) {
//    return (trial_optimality_measure <= current_optimality_measure - this->parameters.gamma * trial_infeasibility_measure) ||
//           (trial_infeasibility_measure < this->parameters.beta * current_infeasibility_measure);
// }


double FunnelStrategy::get_funnel_width(){
   return this->funnel_width;
}

double FunnelStrategy::compute_actual_reduction(double current_optimality_measure, double /*current_infeasibility_measure*/, double trial_optimality_measure) {
   return current_optimality_measure - trial_optimality_measure;
}

//! print: print the current funnel parameter
std::ostream& operator<<(std::ostream& stream, FunnelStrategy& funnel) {
   stream << "************\n";
   stream << "\t\t  Current funnel width:\n";
   stream << "\t\t\t" << funnel.funnel_width << '\n';
   stream << "\t\t************\n";
   return stream;
}














/* check acceptability of step(s) (funnel & sufficient reduction)
 * funnel methods enforce an *unconstrained* sufficient decrease condition
 * precondition: feasible step
 * */
bool FunnelStrategy::is_iterate_acceptable(Statistics& statistics, const Iterate& /*trial_iterate*/,
      const ProgressMeasures& current_progress_measures, const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction,
      double /*objective_multiplier*/) {
   // const double current_optimality_measure = current_progress_measures.optimality(1.) + current_progress_measures.auxiliary_terms;
   // const double trial_optimality_measure = trial_progress_measures.optimality(1.) + trial_progress_measures.auxiliary_terms;
   
   // // unconstrained predicted reduction:
   // // - ignore the predicted infeasibility reduction
   // // - scale the scaled optimality measure with 1
   // const double unconstrained_predicted_reduction = predicted_reduction.optimality(1.) + predicted_reduction.auxiliary_terms;
   // DEBUG << "\t\tCurrent: η = " << current_progress_measures.infeasibility << ",\t ω = " << current_optimality_measure << '\n';
   // DEBUG << "\t\tTrial:   η = " << trial_progress_measures.infeasibility << ",\t ω = " << trial_optimality_measure << '\n';
   // DEBUG << "\t\tUnconstrained predicted reduction: " << predicted_reduction.optimality(1.) << " + " << predicted_reduction.auxiliary_terms <<
   //       " = " <<  unconstrained_predicted_reduction << '\n';

   // GlobalizationStrategy::check_finiteness(current_progress_measures, 1.);
   // GlobalizationStrategy::check_finiteness(trial_progress_measures, 1.);
   // statistics.add_statistic("funnel width", this->funnel->get_funnel_size());
   
   // DEBUG << "\t\t" <<*this->funnel << '\n';

   // bool accept = false;
   // bool funnel_reduction_mechanism = false;
   // bool funnel_acceptable = false;

   // funnel_acceptable = this->is_infeasibility_acceptable(trial_progress_measures.infeasibility);
   
   // // check acceptance   
   // if (funnel_acceptable) {
   //    DEBUG << "\t\tFunnel condition acceptable or in feasibility restoration phase\n";

   //       const double actual_reduction = this->funnel->compute_actual_reduction(current_optimality_measure, current_progress_measures.infeasibility,
   //             trial_optimality_measure);
   //       DEBUG << "\t\tActual reduction: " << actual_reduction << '\n';

   //       DEBUG << "\t\tCheck the switching condition ....\n";
   //       // switching condition: the unconstrained predicted reduction is sufficiently positive
   //       if (this->switching_condition(unconstrained_predicted_reduction, current_progress_measures.infeasibility, this->parameters.delta)) {
   //          // unconstrained Armijo sufficient decrease condition (predicted reduction should be positive)
   //          DEBUG << "\t\tCheck the armijo condition for descent ....\n";
   //          if (this->armijo_sufficient_decrease(unconstrained_predicted_reduction, actual_reduction)) {
   //             DEBUG << "\t\tTrial iterate was ACCEPTED by satisfying Armijo condition\n";
   //             accept = true;
   //             // decrease funnel here ......
   //          }
   //          else { // switching condition holds, but not Armijo condition
   //             DEBUG << "\t\tArmijo condition not satisfied, trial iterate REJECTED\n";
   //          }
   //       }
   //       else {
   //          DEBUG << "\t\tTrial iterate violates switching condition ...\n";
   //          funnel_reduction_mechanism = true;
   //          accept = true; // accept the step and reduce the tr-radius
   //       }

   // } else {
   //    funnel_reduction_mechanism = false;
   //    // step is rejected
   //    DEBUG << "\t\tFunnel condition NOT acceptable\n";
   // }

   // if (funnel_reduction_mechanism){
   //     DEBUG << "\t\tEntering funnel reduction mechanism\n";

   //    // // Feasibility measures
   //    const double current_infeasibility_measure = current_progress_measures.infeasibility;
   //    const double trial_infeasibility_measure = trial_progress_measures.infeasibility;
   //    this->funnel->update_funnel_parameter(current_infeasibility_measure, 
   //                                           trial_infeasibility_measure);
   //    this->funnel_width = this->funnel->get_funnel_size()                                             ;
   // }

   // if (accept){
   //    if (trial_progress_measures.infeasibility < this->funnel->get_funnel_size()){
   //       this->current_iterate_acceptable_to_funnel = true;
   //    }
   //    if (this->current_phase == 1){
   //       //update funnel ....
   //    }
   // }

   // return accept;
   return true;
}
