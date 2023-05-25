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
      const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {

   std::cout << "constraint_jacobian: " << std::endl;

   std::vector<casadi_int> row, col;
   std::vector<double> values;

   casadi_int i=0; //row
   for (const auto& cj_row : constraint_jacobian) {
      cj_row.for_each([&](size_t c, double val) {
         row.push_back(i);
         col.push_back(c);
         values.push_back(val);
      });
      i++;
   }
   DM A = DM::triplet(row, col, values, number_constraints, number_variables);

   row.clear();
   col.clear();
   values.clear();

   hessian.for_each([&](size_t i, size_t j,double entry) {
      row.push_back(i);
      col.push_back(j);
      values.push_back(entry);
   });
   DM H = DM::triplet(row, col, values, number_variables, number_variables);


   SparsityDict qp_struct = {{"a", A.sparsity()}, {"h", H.sparsity()}};
   Function solver = conic("solver", "qrqp", qp_struct,{{"print_problem", true}});

   DMDict args;
   args["x0"] = DM(std::vector<double>(initial_point.begin(), initial_point.begin()+number_variables));
   std::cout << "test" << args << std::endl;
   args["a"] = A;
   args["h"] = H;
   std::vector<double> g(number_variables);
   linear_objective.for_each([&](size_t i, double entry) {
      g[i] = entry;
   });
   args["g"] = DM(g);
   std::vector<double> lbx(number_variables);
   std::vector<double> ubx(number_variables);

   i = 0;
   for (const auto & interval: variables_bounds) {
      lbx[i] = interval.lb;
      ubx[i] = interval.ub;
      uout() << interval.lb << ":" << interval.ub << std::endl;
      i++;
      if (i==number_variables) break;
   }
   args["lbx"] = DM(lbx);
   args["ubx"] = DM(ubx);
   std::vector<double> lba(number_constraints);
   std::vector<double> uba(number_constraints);
   casadi_assert_dev(constraint_bounds.size()==number_constraints);
   i = 0;
   for (const auto & interval: constraint_bounds) {
      std::cout << "i" << i << number_constraints << std::endl;
      lba[i] = interval.lb;
      uba[i] = interval.ub;
      i++;
   }
   args["lba"] = DM(lba);
   args["uba"] = DM(uba);

   solver(args);

   // Populate Direction

}

Direction CASADISolver::solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {

      casadi_error("Not implemented yet");
}
