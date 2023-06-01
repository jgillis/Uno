// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include <algorithm>
#include "BQPDSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "optimization/EvaluationErrors.hpp"
#include "tools/Logger.hpp"
#include "tools/Infinity.hpp"

#define BIG 1e30

extern "C" {
// fortran common block used in bqpd/bqpd.f
extern struct {
   int kk, ll, kkk, lll, mxws, mxlws;
} wsc_;

// fortran common for inertia correction in wdotd
extern struct {
   double alpha;
} kktalphac_;

extern void
bqpd_(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
      double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
      int* info, int* iprint, int* nout);
}

// preallocate a bunch of stuff
BQPDSolver::BQPDSolver(size_t max_number_variables, size_t number_constraints, size_t number_hessian_nonzeros, bool quadratic_programming,
         const Options& options):
      QPSolver(), number_hessian_nonzeros(number_hessian_nonzeros),
      lb(max_number_variables + number_constraints),
      ub(max_number_variables + number_constraints), jacobian(max_number_variables * (number_constraints + 1)),
      jacobian_sparsity(max_number_variables * (number_constraints + 1) + number_constraints + 3),
      kmax(quadratic_programming ? options.get_int("BQPD_kmax") : 0), alp(this->mlp), lp(this->mlp), active_set(max_number_variables + number_constraints),
      w(max_number_variables + number_constraints), gradient_solution(max_number_variables), residuals(max_number_variables + number_constraints),
      e(max_number_variables + number_constraints),
      size_hessian_sparsity(quadratic_programming ? number_hessian_nonzeros + max_number_variables + 3 : 0),
      size_hessian_workspace(number_hessian_nonzeros + this->kmax * (this->kmax + 9) / 2 + 2 * max_number_variables + number_constraints + this->mxwk0),
      size_hessian_sparsity_workspace(this->size_hessian_sparsity + this->kmax + this->mxiwk0),
      hessian_values(this->size_hessian_workspace), hessian_sparsity(this->size_hessian_sparsity_workspace),
      print_subproblem(options.get_bool("BQPD_print_subproblem")) {
   // default active set
   for (size_t i: Range(max_number_variables + number_constraints)) {
      this->active_set[i] = static_cast<int>(i) + this->fortran_shift;
   }
}

Direction BQPDSolver::solve_QP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const SymmetricMatrix<double>& hessian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {
   if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
      this->save_lagrangian_hessian_to_local_format(hessian);
   }
   if (this->print_subproblem) {
      DEBUG << "QP:\n";
      DEBUG << "Hessian: " << hessian;
   }
   return this->solve_subproblem(number_variables, number_constraints, variables_bounds, constraint_bounds, linear_objective, constraint_jacobian,
         initial_point, warmstart_information);
}

Direction BQPDSolver::solve_LP(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {
   if (this->print_subproblem) {
      DEBUG << "LP:\n";
   }
   return this->solve_subproblem(number_variables, number_constraints, variables_bounds, constraint_bounds, linear_objective, constraint_jacobian,
         initial_point, warmstart_information);
}

Direction BQPDSolver::solve_subproblem(size_t number_variables, size_t number_constraints, const std::vector<Interval>& variables_bounds,
      const std::vector<Interval>& constraint_bounds, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian, const std::vector<double>& initial_point,
      const WarmstartInformation& warmstart_information) {
   // initialize wsc_ common block (Hessian & workspace for BQPD)
   // setting the common block here ensures that several instances of BQPD can run simultaneously
   wsc_.kk = static_cast<int>(this->number_hessian_nonzeros);
   wsc_.ll = static_cast<int>(this->size_hessian_sparsity);
   wsc_.mxws = static_cast<int>(this->size_hessian_workspace);
   wsc_.mxlws = static_cast<int>(this->size_hessian_sparsity_workspace);
   kktalphac_.alpha = 0; // inertia control

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

   // Jacobian (objective and constraints)
   if (warmstart_information.objective_changed || warmstart_information.constraints_changed) {
      this->save_gradients_to_local_format(number_constraints, linear_objective, constraint_jacobian);
   }

   // set variable bounds
   if (warmstart_information.variable_bounds_changed) {
      for (size_t i: Range(number_variables)) {
         this->lb[i] = (variables_bounds[i].lb == -INF<double>) ? -BIG : variables_bounds[i].lb;
         this->ub[i] = (variables_bounds[i].ub == INF<double>) ? BIG : variables_bounds[i].ub;
      }
   }
   // set constraint bounds
   if (warmstart_information.constraint_bounds_changed) {
      for (size_t j: Range(number_constraints)) {
         this->lb[number_variables + j] = (constraint_bounds[j].lb == -INF<double>) ? -BIG : constraint_bounds[j].lb;
         this->ub[number_variables + j] = (constraint_bounds[j].ub == INF<double>) ? BIG : constraint_bounds[j].ub;
      }
   }

   Direction direction(number_variables, number_constraints);
   copy_from(direction.primals, initial_point);
   const int n = static_cast<int>(number_variables);
   const int m = static_cast<int>(number_constraints);

   BQPDMode mode = this->determine_mode(warmstart_information);
   const int mode_integer = static_cast<int>(mode);

   // solve the LP/QP
   bqpd_(&n, &m, &this->k, &this->kmax, this->jacobian.data(), this->jacobian_sparsity.data(), direction.primals.data(), this->lb.data(),
         this->ub.data(), &direction.subproblem_objective, &this->fmin, this->gradient_solution.data(), this->residuals.data(), this->w.data(),
         this->e.data(), this->active_set.data(), this->alp.data(), this->lp.data(), &this->mlp, &this->peq_solution, this->hessian_values.data(),
         this->hessian_sparsity.data(), &mode_integer, &this->ifail, this->info.data(), &this->iprint, &this->nout);
   BQPDStatus bqpd_status = BQPDSolver::bqpd_status_from_int(this->ifail);
   direction.status = BQPDSolver::status_from_bqpd_status(bqpd_status);
   this->number_calls++;

   // project solution into bounds
   for (size_t i: Range(number_variables)) {
      direction.primals[i] = std::min(std::max(direction.primals[i], variables_bounds[i].lb), variables_bounds[i].ub);
   }
   this->analyze_constraints(number_variables, number_constraints, direction);
   
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

BQPDMode BQPDSolver::determine_mode(const WarmstartInformation& warmstart_information) const {
   BQPDMode mode = (this->number_calls == 0) ? BQPDMode::ACTIVE_SET_EQUALITIES : BQPDMode::USER_DEFINED;
   // if problem changed, use cold start
   if (warmstart_information.problem_changed) {
      mode = BQPDMode::ACTIVE_SET_EQUALITIES;
   }
   // if only the variable bounds changed, reuse the active set estimate and the Jacobian information
   else if (warmstart_information.variable_bounds_changed && not warmstart_information.objective_changed &&
         not warmstart_information.constraints_changed && not warmstart_information.constraint_bounds_changed) {
      mode = BQPDMode::UNCHANGED_ACTIVE_SET_AND_JACOBIAN;
   }
   return mode;
}

// save Hessian (in arbitrary format) to a "weak" CSC format: compressed columns but row indices are not sorted, nor unique
void BQPDSolver::save_lagrangian_hessian_to_local_format(const SymmetricMatrix<double>& hessian) {
   const size_t header_size = 1;
   // pointers withing the single array
   int* row_indices = &this->hessian_sparsity[header_size];
   int* column_starts = &this->hessian_sparsity[header_size + hessian.number_nonzeros];
   // header
   this->hessian_sparsity[0] = static_cast<int>(hessian.number_nonzeros + 1);
   // count the elements in each column
   for (size_t j: Range(hessian.dimension + 1)) {
      column_starts[j] = 0;
   }
   hessian.for_each([&](size_t /*i*/, size_t j, double /*entry*/) {
      column_starts[j + 1]++;
   });
   // carry over the column starts
   for (size_t j: Range(1, hessian.dimension + 1)) {
      column_starts[j] += column_starts[j-1];
      column_starts[j-1] += this->fortran_shift;
   }
   column_starts[hessian.dimension] += this->fortran_shift;
   // copy the entries
   std::vector<int> current_indices(hessian.dimension);
   hessian.for_each([&](size_t i, size_t j, double entry) {
      const size_t index = static_cast<size_t>(column_starts[j] + current_indices[j] - this->fortran_shift);
      assert(index <= static_cast<size_t>(column_starts[j+1]) &&
         "BQPD: error in converting the Hessian matrix to the local format. Try setting the sparse format to CSC");
      this->hessian_values[index] = entry;
      row_indices[index] = static_cast<int>(i) + this->fortran_shift;
      current_indices[j]++;
   });
}

void BQPDSolver::save_gradients_to_local_format(size_t number_constraints, const SparseVector<double>& linear_objective,
      const RectangularMatrix<double>& constraint_jacobian) {
   size_t current_index = 0;
   linear_objective.for_each([&](size_t i, double derivative) {
      this->jacobian[current_index] = derivative;
      this->jacobian_sparsity[current_index + 1] = static_cast<int>(i) + this->fortran_shift;
      current_index++;
   });
   for (size_t j: Range(number_constraints)) {
      constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->jacobian[current_index] = derivative;
         this->jacobian_sparsity[current_index + 1] = static_cast<int>(i) + this->fortran_shift;
         current_index++;
      });
   }
   current_index++;
   this->jacobian_sparsity[0] = static_cast<int>(current_index);
   // header
   size_t size = 1;
   this->jacobian_sparsity[current_index] = static_cast<int>(size);
   current_index++;
   size += linear_objective.size();
   this->jacobian_sparsity[current_index] = static_cast<int>(size);
   current_index++;
   for (size_t j: Range(number_constraints)) {
      size += constraint_jacobian[j].size();
      this->jacobian_sparsity[current_index] = static_cast<int>(size);
      current_index++;
   }
}

void BQPDSolver::analyze_constraints(size_t number_variables, size_t number_constraints, Direction& direction) {
   ConstraintPartition constraint_partition(number_constraints);

   // active constraints
   for (size_t j: Range(number_variables - static_cast<size_t>(this->k))) {
      const size_t index = static_cast<size_t>(std::abs(this->active_set[j]) - this->fortran_shift);

      if (index < number_variables) {
         // bound constraint
         if (0 <= this->active_set[j]) { // lower bound active
            direction.multipliers.lower_bounds[index] = this->residuals[index];
            direction.active_set.bounds.at_lower_bound.push_back(index);
         }
         else { // upper bound active */
            direction.multipliers.upper_bounds[index] = -this->residuals[index];
            direction.active_set.bounds.at_upper_bound.push_back(index);
         }
      }
      else {
         // general constraint
         size_t constraint_index = index - number_variables;
         constraint_partition.feasible.push_back(constraint_index);
         if (0 <= this->active_set[j]) { // lower bound active
            direction.multipliers.constraints[constraint_index] = this->residuals[index];
            direction.active_set.constraints.at_lower_bound.push_back(constraint_index);
         }
         else { // upper bound active
            direction.multipliers.constraints[constraint_index] = -this->residuals[index];
            direction.active_set.constraints.at_upper_bound.push_back(constraint_index);
         }
      }
   }

   // inactive constraints
   for (size_t j: Range(number_variables - static_cast<size_t>(this->k), number_variables + number_constraints)) {
      size_t index = static_cast<size_t>(std::abs(this->active_set[j]) - this->fortran_shift);

      if (number_variables <= index) { // general constraints
         size_t constraint_index = index - number_variables;
         if (this->residuals[index] < 0.) { // infeasible constraint
            constraint_partition.infeasible.push_back(constraint_index);
            if (this->active_set[j] < 0) { // upper bound violated
               constraint_partition.upper_bound_infeasible.push_back(constraint_index);
            }
            else { // lower bound violated
               constraint_partition.lower_bound_infeasible.push_back(constraint_index);
            }
         }
         else { // feasible constraint
            constraint_partition.feasible.push_back(constraint_index);
         }
      }
   }
   direction.constraint_partition = constraint_partition;
}

BQPDStatus BQPDSolver::bqpd_status_from_int(int ifail) {
   assert(0 <= ifail && ifail <= 9 && "BQPDSolver.bqpd_status_from_int: ifail does not belong to [0, 9]");
   return static_cast<BQPDStatus>(ifail);
}

SubproblemStatus BQPDSolver::status_from_bqpd_status(BQPDStatus bqpd_status) {
   switch (bqpd_status) {
      case BQPDStatus::OPTIMAL:
         return SubproblemStatus::OPTIMAL;
      case BQPDStatus::UNBOUNDED_PROBLEM:
         return SubproblemStatus::UNBOUNDED_PROBLEM;
      case BQPDStatus::BOUND_INCONSISTENCY:
         WARNING << YELLOW << "BQPD error: bound inconsistency\n" << RESET;
         return SubproblemStatus::INFEASIBLE;
      case BQPDStatus::INFEASIBLE:
         return SubproblemStatus::INFEASIBLE;
      // errors
      case BQPDStatus::INCORRECT_PARAMETER:
         WARNING << YELLOW << "BQPD error: incorrect parameter\n" << RESET;
         return SubproblemStatus::ERROR;
      case BQPDStatus::LP_INSUFFICIENT_SPACE:
         WARNING << YELLOW << "BQPD error: LP insufficient space\n" << RESET;
         return SubproblemStatus::ERROR;
      case BQPDStatus::HESSIAN_INSUFFICIENT_SPACE:
         WARNING << YELLOW << "BQPD kmax too small, continue anyway\n" << RESET;
         return SubproblemStatus::ERROR;
      case BQPDStatus::SPARSE_INSUFFICIENT_SPACE:
         WARNING << YELLOW << "BQPD error: sparse insufficient space\n" << RESET;
         return SubproblemStatus::ERROR;
      case BQPDStatus::MAX_RESTARTS_REACHED:
         WARNING << YELLOW << "BQPD max restarts reached\n" << RESET;
         return SubproblemStatus::ERROR;
      case BQPDStatus::UNDEFINED:
         WARNING << YELLOW << "BQPD error: undefined\n" << RESET;
         return SubproblemStatus::ERROR;
   }
   throw std::invalid_argument("The BQPD ifail is not consistent with the Uno status values");
}