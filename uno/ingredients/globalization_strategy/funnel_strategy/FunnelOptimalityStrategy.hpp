// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FUNNELOPTIMALITYSTRATEGY_H
#define UNO_FUNNELOPTIMALITYSTRATEGY_H

#include "../GlobalizationStrategy.hpp"
#include "funnel/Funnel.hpp"
#include "tools/Options.hpp"
#include "tools/Infinity.hpp"
#include "FunnelStrategy.hpp"

/*! \class TwoPhaseConstants
 * \brief Constants for funnel strategy
 *
 *  Set of constants to control the funnel strategy
 */
// struct FunnelStrategyParameters {
//    double kappa_initial_upper_bound;
//    double kappa_initial_multiplication;
//    double delta; /*!< Switching constant */
//    double upper_bound;
//    double infeasibility_fraction;
//    double switching_infeasibility_exponent;
// };

// enum class Phase {FEASIBILITY_RESTORATION = 1, OPTIMALITY = 2};

/*! \class FunnelStrategy
 * \brief Step acceptance strategy based on a funnel
 *
 *  Strategy that accepts or declines a trial step
 */
class FunnelOptimalityStrategy : public FunnelStrategy {
public:
   explicit FunnelOptimalityStrategy(Statistics& statistics, const Options& options);

   // void initialize(const Iterate& initial_iterate) override;
   // [[nodiscard]] bool is_infeasibility_acceptable(double infeasibility_measure) const override;
   [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, const Iterate& trial_iterate, const ProgressMeasures& current_progress_measures,
         const ProgressMeasures& trial_progress_measures, const ProgressMeasures& predicted_reduction, double objective_multiplier) override;
   
   void update_funnel_width(double current_infeasibility_measure, double trial_infeasibility_measure);
   // void reset() override;
   // void register_current_progress(const ProgressMeasures& current_progress_measures) override;

// protected:
//    // pointer to allow polymorphism
//    const std::unique_ptr<Funnel> funnel;
//    double initial_funnel_upper_bound{INF<double>};
//    // const FunnelStrategyParameters parameters; /*!< Set of constants */


//    [[nodiscard]] bool switching_condition(double predicted_reduction, double current_infeasibility, double switching_fraction) const;
};

#endif // UNO_FUNNELSTRATEGY_H
