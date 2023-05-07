// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <algorithm>
#include "CasadiSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "tools/Logger.hpp"
#include "tools/Infinity.hpp"

using namespace casadi;

typedef std::map<std::string, Sparsity> SparsityDict;


// preallocate a bunch of stuff
CASADISolver::CASADISolver(size_t max_number_variables, size_t number_constraints, size_t number_hessian_nonzeros, bool quadratic_programming,
         const Options& options):
      QPSolver() {
}

Direction CASADISolver::solve_QP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point) {

   SparsityDict qp_struct;
   Function solver = conic("solver", "qrqp", qp_struct);
      
 
}

Direction CASADISolver::solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point) {

}
