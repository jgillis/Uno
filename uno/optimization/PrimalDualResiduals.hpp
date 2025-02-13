// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALRESIDUALS_H
#define UNO_PRIMALDUALRESIDUALS_H

#include "tools/Infinity.hpp"

class PrimalDualResiduals {
public:
   double optimality_stationarity{INF<double>};
   double feasibility_stationarity{INF<double>};
   double infeasibility{INF<double>};
   double optimality_complementarity{INF<double>};
   double feasibility_complementarity{INF<double>};
   double stationarity_scaling{INF<double>};
   double complementarity_scaling{INF<double>};
};

#endif // UNO_PRIMALDUALRESIDUALS_H