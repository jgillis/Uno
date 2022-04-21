#include "HessianModel.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "solvers/linear/LinearSolverFactory.hpp"

HessianModel::HessianModel(size_t dimension, size_t maximum_number_nonzeros, const std::string& sparse_format) :
      hessian(SymmetricMatrixFactory::create(sparse_format, dimension, maximum_number_nonzeros)) {
}

// Exact Hessian
ExactHessian::ExactHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options) :
   HessianModel(dimension, maximum_number_nonzeros, options.at("sparse_format")) {
}

void ExactHessian::evaluate(const NonlinearReformulation& problem, const std::vector<double>& primal_variables, const std::vector<double>& constraint_multipliers) {
   // evaluate Lagrangian Hessian
   problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, *this->hessian);
   this->evaluation_count++;
}

// Convexified Hessian
ConvexifiedHessian::ConvexifiedHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options):
      HessianModel(dimension, maximum_number_nonzeros, options.at("sparse_format")),
      linear_solver(LinearSolverFactory::create(options.at("linear_solver"), dimension, maximum_number_nonzeros)),
      regularization_initial_value(stod(options.at("regularization_initial_value"))),
      regularization_increase_factor(stod(options.at("regularization_increase_factor"))) {
}

void ConvexifiedHessian::evaluate(const NonlinearReformulation& problem, const std::vector<double>& primal_variables, const std::vector<double>& constraint_multipliers) {
   // evaluate Lagrangian Hessian
   problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, *this->hessian);
   this->evaluation_count++;
   // regularize (only on the original variables) to make the problem strictly convex
   DEBUG << "hessian before convexification: " << *this->hessian;
   this->regularize(*this->hessian, problem.get_number_original_variables());
}

// Nocedal and Wright, p51
void ConvexifiedHessian::regularize(SymmetricMatrix& matrix, size_t number_original_variables) {
   //assert(size_block_to_regularize <= matrix.dimension && "The block to regularize is larger than the matrix");

   const double smallest_diagonal_entry = matrix.smallest_diagonal_entry();
   DEBUG << "The minimal diagonal entry of the matrix is " << matrix.smallest_diagonal_entry() << "\n";

   double regularization = (smallest_diagonal_entry <= 0.) ? this->regularization_initial_value - smallest_diagonal_entry : 0.;
   /*
   double regularization = (smallest_diagonal_entry <= 0.) ? this->regularization_initial_value - smallest_diagonal_entry : this->regularization_initial_value;
   */
   bool good_inertia = false;
   bool regularized = false;
   while (!good_inertia) {
      DEBUG << "Testing factorization with regularization factor " << regularization << "\n";
      if (0. < regularization) {
         // if regularization was already perform, delete the previous entries
         if (regularized) {
            for (size_t i = 0; i < number_original_variables; i++) {
               matrix.pop();
            }
         }
         matrix.add_identity_multiple(regularization, number_original_variables);
         regularized = true;
      }
      this->linear_solver->do_symbolic_factorization(matrix);
      this->linear_solver->do_numerical_factorization(matrix);

      if (this->linear_solver->rank() == number_original_variables && this->linear_solver->number_negative_eigenvalues() == 0) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n";
      }
      else {
         DEBUG << "rank: " << this->linear_solver->rank() << ", negative eigenvalues: " << this->linear_solver->number_negative_eigenvalues() << "\n";
         regularization = (regularization == 0.) ? this->regularization_initial_value : this->regularization_increase_factor * regularization;
         assert(regularization < std::numeric_limits<double>::infinity() && "The regularization coefficient diverged");
      }
   }
}

// Factory
std::unique_ptr<HessianModel> HessianModelFactory::create(const std::string& hessian_model, size_t dimension, size_t maximum_number_nonzeros,
      bool convexify, const Options& options) {
   if (hessian_model == "exact") {
      if (convexify) {
         return std::make_unique<ConvexifiedHessian>(dimension, maximum_number_nonzeros, options);
      }
      else {
         return std::make_unique<ExactHessian>(dimension, maximum_number_nonzeros, options);
      }
   }
   throw std::invalid_argument("Hessian model " + hessian_model + " does not exist");
}