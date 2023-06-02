// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CASADISOLVER_H
#define UNO_CASADISOLVER_H

#include <vector>
#include "QPSolver.hpp"
#include "solvers/LP/LPSolver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "tools/Options.hpp"

#include <casadi/casadi.hpp>

class CASADISolver : public QPSolver {
public:
   CASADISolver(size_t max_number_variables, size_t number_constraints, size_t number_hessian_nonzeros, bool quadratic_programming, const Options& options);

   Direction solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) override;

   Direction solve_QP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
         const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
         const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point,
         const WarmstartInformation& warmstart_information) override;

private:
      const int fortran_shift{1};

      size_t number_calls{0};
      const bool print_subproblem;

};

#endif // UNO_CASADISOLVER_H
