// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "filter_method/LeyfferFilterMethod.hpp"
#include "filter_method/WaechterFilterMethod.hpp"
#include "funnel_method/FunnelMethod.hpp"
#include "funnel_method/FunnelRestorationMethod.hpp"
#include "funnel_method/FunnelOptimalityMethod.hpp"

std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(Statistics& statistics, const std::string& strategy_type,
      bool accept_when_switching_violated, const Options& options) {
   if (strategy_type == "l1_merit") {
      return std::make_unique<l1MeritFunction>(statistics, options);
   }
   else if (strategy_type == "leyffer_filter_method") {
      return std::make_unique<LeyfferFilterMethod>(statistics, accept_when_switching_violated, options);
   }
   else if (strategy_type == "waechter_filter_method") {
      return std::make_unique<WaechterFilterMethod>(statistics, options);
   }
   else if (strategy_type == "funnel_method") {
      return std::make_unique<FunnelMethod>(statistics, options);
   }
   else if (strategy_type == "funnel_restoration_method") {
      return std::make_unique<FunnelRestorationMethod>(statistics, options);
   }
   else if (strategy_type == "funnel_optimality_method") {
      return std::make_unique<FunnelOptimalityMethod>(statistics, options);
   }
   throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
}

std::vector<std::string> GlobalizationStrategyFactory::available_strategies() {
   return {"l1_merit", "leyffer_filter_strategy", "waechter_filter_strategy", "funnel_strategy"};
}