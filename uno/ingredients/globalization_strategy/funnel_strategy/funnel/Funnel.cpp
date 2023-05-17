// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <algorithm>
#include "Funnel.hpp"
#include "tools/Logger.hpp"
#include "tools/Range.hpp"

Funnel::Funnel(const Options& options) :
      capacity(options.get_unsigned_int("max_iterations")),
      funnel_bounds(this->capacity), //length of funnel parameter vector
      infeasibility(this->capacity), //length of infeasibility vector
      optimality(this->capacity), //length of optimality vector
      parameters({
         options.get_double("funnel_kappa_infeasibility_1"),
         options.get_double("funnel_kappa_infeasibility_2"),
         options.get_double("funnel_beta"),
         options.get_double("funnel_gamma")
      }) {
   this->reset();
}

void Funnel::initialize(){
   this->current_upper_bound = this->initial_upper_bound;
   DEBUG << "Initial funnel parameter is: " << this->current_upper_bound << "\n"; 
}

void Funnel::reset() {
   // this->current_upper_bound = INF<double>;
   this->current_upper_bound = this->initial_upper_bound;
   this->number_entries = 0;
}

void Funnel::update_funnel_parameter(double current_infeasibility_measure, double trial_infeasibility_measure) {

   this->current_upper_bound = std::max(this->parameters.kappa_infeasibility_1 *this->current_upper_bound, 
      trial_infeasibility_measure + this->parameters.kappa_infeasibility_2 * (current_infeasibility_measure - trial_infeasibility_measure));

   DEBUG << "\t\tNew funnel parameter is: " << this->current_upper_bound << "\n"; 
   
}

double Funnel::get_funnel_size(){
   return this->current_upper_bound;
}

// bool Fu

//  add (infeasibility_measure, optimality_measure) to the funnel
void Funnel::add(double infeasibility_measure, double optimality_measure) {
   
   this->number_entries++;

   this->funnel_bounds[this->number_entries] = this->current_upper_bound;
   this->infeasibility[this->number_entries] = infeasibility_measure;
   this->optimality[this->number_entries] = optimality_measure;
}

// bool Funnel::acceptable_wrt_upper_bound(double infeasibility_measure) const {
//    return (infeasibility_measure < this->parameters.beta * this->upper_bound);
// }

// // return true if (infeasibility_measure, optimality_measure) acceptable, false otherwise
// bool Funnel::acceptable(double infeasibility_measure, double optimality_measure) {
//    // check upper bound first
//    if (not this->acceptable_wrt_upper_bound(infeasibility_measure)) {
//       DEBUG << "Rejected because of funnel upper bound\n";
//       return false;
//    }

//    // TODO: use binary search
//    size_t position = 0;
//    while (position < this->number_entries && infeasibility_measure >= this->parameters.beta * this->infeasibility[position]) {
//       position++;
//    }

//    // check acceptability
//    if (position == 0) {
//       return true; // acceptable as left-most entry
//    }
//    // until here, the optimality measure was not required
//    else if (optimality_measure <= this->optimality[position - 1] - this->parameters.gamma * infeasibility_measure) {
//       return true; // point acceptable
//    }
//    DEBUG << "Rejected because of funnel domination\n";
//    return false;
// }

bool Funnel::acceptable(double infeasibility_measure) {

   if (infeasibility_measure <= this->current_upper_bound){
      return true;
   }
   else {
      DEBUG << "\t\tREJECTED because of funnel condition.\n";
      return false;
   }
}

//! check acceptability wrt current point
bool Funnel::acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
      double trial_optimality_measure) {
   return (trial_optimality_measure <= current_optimality_measure - this->parameters.gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->parameters.beta * current_infeasibility_measure);
}

double Funnel::compute_actual_reduction(double current_optimality_measure, double /*current_infeasibility_measure*/, double trial_optimality_measure) {
   return current_optimality_measure - trial_optimality_measure;
}

//! print: print the current funnel parameter
std::ostream& operator<<(std::ostream& stream, Funnel& funnel) {
   stream << "************\n";
   stream << "\t\t  Current funnel parameter:\n";
   stream << "\t\t\t" << funnel.current_upper_bound << '\n';
   stream << "\t\t************\n";
   return stream;
}