// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "FunnelFactory.hpp"

// FunnelFactory class
std::unique_ptr<Funnel> FunnelFactory::create(const Options& options) {
   return std::make_unique<Funnel>(options);
}