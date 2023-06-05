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
CASADISolver::CASADISolver(size_t max_number_variables,
         size_t number_constraints,
         size_t number_hessian_nonzeros,
         bool quadratic_programming,
         const Options& options):
      QPSolver(),
      print_subproblem(options.get_bool("BQPD_print_subproblem")) {
}

Direction CASADISolver::solve_QP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {

   // Taken from BQPD
   if (this->print_subproblem) {
      DEBUG << "QP:\n";
      DEBUG << "Hessian: " << hessian;
   }

   // Taken from BQPD
   if (this->print_subproblem) {
      DEBUG << "objective gradient: " << linear_objective;
      for (size_t j: Range(number_constraints)) {
         DEBUG << "gradient c" << j << ": " << constraint_jacobian[j];
      }
      for (size_t i: Range(number_variables)) {
         DEBUG << "d_x" << i << " in [" << variables_bounds[i].lb << ", " << variables_bounds[i].ub << "]\n";
      }
      for (size_t j: Range(number_constraints)) {
         DEBUG << "linearized c" << j << " in [" << constraint_bounds[j].lb << ", " << constraint_bounds[j].ub << "]\n";
      }
   }

   // ---------------------------------------------------
   // Prepare argument dictionary for Casadi solver
   // ---------------------------------------------------

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

   // It is a symmetric matrix, but Uno just saves the upper triangular part
   hessian.for_each([&](size_t i, size_t j,double entry) {
      row.push_back(i);
      col.push_back(j);
      values.push_back(entry);

      if (i != j){
         // Lower diagonal part also needed
         row.push_back(j);
         col.push_back(i);
         values.push_back(entry);
      }
   });
   DM H = DM::triplet(row, col, values, number_variables, number_variables);
   
   DMDict args;
   args["x0"] = DM(std::vector<double>(initial_point.begin(), initial_point.begin()+number_variables));
   args["a"] = A;
   args["h"] = H;

   DEBUG << "direction initial point: \n";
   for (size_t i: Range(number_variables)) {
         DEBUG <<  initial_point[i] << "\n";
      }

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
      // uout() << interval.lb << ":" << interval.ub << std::endl;
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
      // std::cout << "i" << i << number_constraints << std::endl;
      lba[i] = interval.lb;
      uba[i] = interval.ub;
      i++;
   }
   args["lba"] = DM(lba);
   args["uba"] = DM(uba);

   std::cout << "Casadi input to the solver: " << args << std::endl;

   // ---------------------------------------------------
   // Create the Casadi QP solver and solve the problem
   // ---------------------------------------------------

   SparsityDict qp_struct = {{"a", A.sparsity()}, {"h", H.sparsity()}};

   // Dict opts_osqp;
   // opts_osqp["verbose"] = false;
   // opts_osqp["eps_abs"] = 1e-8;
   // opts_osqp["eps_rel"] = 1e-8;
   Dict opts_highs;
   opts_highs["output_flag"] = false;

   Dict opts_conic;
   opts_conic["highs"] = opts_highs;

   // opts_conic["printLevel"] = "none";
   opts_conic["verbose"] = true;
   opts_conic["print_problem"] = false;
   opts_conic["error_on_fail"] = false;
   Function solver = conic("solver", "highs", qp_struct, opts_conic);

   DMDict res = solver(args);
   Dict memory_solver = solver.stats();


   // ---------------------------------------------------
   // Postprocess the direction
   // ---------------------------------------------------
   Direction direction(number_variables, number_constraints);

   // Solver status
   // ----------------
   direction.status = CasadiSolver::status_from_bqpd_status(memory_solver["success"],
                                                            memory_solver["unified_return_status"]);
   std::cout << "QP success: " << memory_solver["success"] << std::endl;
   std::cout << "Unified Return Status: " << memory_solver["unified_return_status"] << std::endl;
   // Primal variables
   // ----------------
   copy_from(direction.primals, res["x"].nonzeros());
   //direction.status = ..... tbd
   this->number_calls++;

   // project solution into bounds
   for (size_t i: Range(number_variables)) {
      direction.primals[i] = std::min(std::max(direction.primals[i], variables_bounds[i].lb), variables_bounds[i].ub);
   }
   
   // Dual variables
   // ----------------
   // Analyze the constraints here ..................
   direction.subproblem_objective = res["cost"].nonzeros().front();
   ConstraintPartition constraint_partition(number_constraints);

   // TODO: check signs (validate with BQPSolver answer)
   //       do we need to construct activate set? see BQPDSolver::analyze_constraints

   for (size_t i: Range(number_variables)) {
         
         if (res["lam_x"].nonzeros()[i] < 0){ //lower bounds active
            direction.multipliers.lower_bounds[i] = -res["lam_x"].nonzeros()[i];
            direction.multipliers.upper_bounds[i] = 0.0;

            direction.active_set.bounds.at_lower_bound.push_back(i);
         } else if (res["lam_x"].nonzeros()[i] > 0) { // upper bound active
            direction.multipliers.lower_bounds[i] = 0.0;
            direction.multipliers.upper_bounds[i] = -res["lam_x"].nonzeros()[i];

            direction.active_set.bounds.at_upper_bound.push_back(i);
         } else {
            direction.multipliers.lower_bounds[i] = 0.0;
            direction.multipliers.upper_bounds[i] = 0.0;
         }
   }

   // Sign convention of Casadi and Uno is apparently different
   for (size_t j: Range(number_constraints)) {
         // so far, we just deal with feasible subproblems, so:
         constraint_partition.feasible.push_back(j);
         if (res["lam_a"].nonzeros()[j] != 0){
            direction.multipliers.constraints[j] = -res["lam_a"].nonzeros()[j];

            if (res["lam_a"].nonzeros()[j] < 0){ // lower bound active
               direction.active_set.constraints.at_lower_bound.push_back(j);
            } else {
               direction.active_set.constraints.at_upper_bound.push_back(j);
            }
         }
      }

   direction.constraint_partition = constraint_partition;

   // Analysis not over yet ......
   DEBUG << "direction multipliers ub: \n";
   for (size_t i: Range(number_variables)) {
         DEBUG <<  direction.multipliers.upper_bounds[i] << "\n";
      }
   DEBUG << "direction multipliers lb: \n";
   for (size_t i: Range(number_variables)) {
         DEBUG << direction.multipliers.lower_bounds[i] << "\n";
      }
   DEBUG << "direction constraints multipliers: \n";
   for (size_t j: Range(number_constraints)) {
         DEBUG << direction.multipliers.constraints[j]<< "\n";
      }
   return direction;

}

SubproblemStatus CasadiSolver::status_from_casadi_status(bool success, std::string casadi_status) {
   
   if (success == true){
      return SubproblemStatus::OPTIMAL;
   } else {
      switch (casadi_status) {
      case "Infeasible":
         return SubproblemStatus::INFEASIBLE;
      default:
         DEBUG << "The return status is " << casadi_status;
         WARNING << YELLOW << " error: ...\n" << RESET;
         return SubproblemStatus::ERROR;
      }
   throw std::invalid_argument("The Casadi solver ifail is not consistent with the Uno status values");
   }
}

Direction CASADISolver::solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {

      casadi_error("Not implemented yet");
}
