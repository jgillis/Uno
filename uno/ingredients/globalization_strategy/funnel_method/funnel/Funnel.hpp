// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FUNNEL_H
#define UNO_FUNNEL_H

#include <vector>
#include <memory>
#include "tools/Options.hpp"
#include "tools/Infinity.hpp"

struct FunnelParameters {
   double kappa_infeasibility_1;
   double kappa_infeasibility_2;
   double beta; /*!< Margin around funnel */
   double gamma; /*!< Margin around funnel (sloping margin) */
};

class Funnel {
public:
   double initial_upper_bound{INF<double>}; /*!< Upper bound on constraint violation */

   explicit Funnel(const Options& options);
   virtual ~Funnel() = default;

   void reset();
   // [[nodiscard]] double get_smallest_infeasibility() const;
   virtual void add(double infeasibility_measure, double optimality_measure);
   virtual void initialize();
   
   virtual void update_funnel_parameter(double current_infeasibility_measure, double trial_infeasibility_measure);

   virtual bool acceptable(double infeasibility_measure);
   virtual bool acceptable_wrt_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure);

   virtual double compute_actual_reduction(double current_optimality_measure, double current_infeasibility_measure, double trial_optimality_measure);

   friend std::ostream& operator<<(std::ostream& stream, Funnel& funnel);

   double get_funnel_size();

protected:
   double current_upper_bound; // funnel parameter
   const size_t capacity; /*!< Max funnel size */

   std::vector<double> funnel_bounds{}; // all Funnel parameters
   std::vector<double> infeasibility{}; // of all iterates
   std::vector<double> optimality{}; // of all iterates

   size_t number_entries{0};
   const FunnelParameters parameters; /*!< Set of parameters */

   // [[nodiscard]] bool is_empty() const;
   // [[nodiscard]] bool acceptable_wrt_upper_bound(double infeasibility_measure) const;
   // void left_shift(size_t start, size_t shift_size);
   // void right_shift(size_t start, size_t shift_size);
};

#endif // UNO_FILTER_H
