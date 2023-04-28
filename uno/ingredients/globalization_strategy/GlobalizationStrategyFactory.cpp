// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "GlobalizationStrategyFactory.hpp"
#include "l1MeritFunction.hpp"
#include "filter_strategy/LeyfferFilterStrategy.hpp"
#include "filter_strategy/WaechterFilterStrategy.hpp"
#include "funnel_strategy/FunnelStrategy.hpp"

std::unique_ptr <GlobalizationStrategy> GlobalizationStrategyFactory::create(Statistics& statistics, const std::string& strategy_type,
      const Options& options) {
   if (strategy_type == "l1_merit") {
      return std::make_unique<l1MeritFunction>(statistics, options);
   }
   else if (strategy_type == "leyffer_filter_strategy") {
      return std::make_unique<LeyfferFilterStrategy>(statistics, options);
   }
   else if (strategy_type == "waechter_filter_strategy") {
      return std::make_unique<WaechterFilterStrategy>(statistics, options);
   }
   else if (strategy_type == "funnel_strategy") {
      return std::make_unique<LeyfferFilterStrategy>(statistics, options);
   }
   throw std::invalid_argument("GlobalizationStrategy " + strategy_type + " is not supported");
}

std::vector<std::string> GlobalizationStrategyFactory::available_strategies() {
   return {"l1_merit", "leyffer_filter_strategy", "waechter_filter_strategy"};
}