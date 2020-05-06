#include <cmath>
#include "BQPDSolver.hpp"

#define BIG 1e30

extern "C" {

    /* fortran common block used in bqpd/bqpd.f */
    extern struct {
        int kk, ll, kkk, lll, mxws, mxlws;
    } wsc_;

    /* fortran common for inertia correction in wdotd */
    extern struct {
        double alpha;
    } kktalphac_;

    extern void bqpd_(int *n, int *m, int *k, int *kmax, double *a, int *la, double *x, double *bl, double *bu,
            double *f, double *fmin, double *g, double *r, double *w, double *e, int *ls, double *alp, int *lp,
            int *mlp, int *peq, double *ws, int *lws, int *mode, int *ifail, int *info, int *iprint, int *nout);
}

/* preallocate a bunch of stuff */
BQPDSolver::BQPDSolver(int number_variables, int number_constraints, int maximum_number_nonzeros):
QPSolver(), n_(number_variables), m_(number_constraints), maximum_number_nonzeros(maximum_number_nonzeros), lb(n_ + m_), ub(n_ + m_), use_fortran(1), kmax_(500), mlp_(1000), mxwk0_(2000000), mxiwk0_(500000), info_(100), alp_(mlp_), lp_(mlp_), ls_(n_ + m_), w_(n_ + m_), gradient_solution_(n_), residuals_(n_ + m_), e_(n_ + m_), nhr_(maximum_number_nonzeros), nhi_(maximum_number_nonzeros + n_ + 3), mxws_(nhr_ + kmax_ * (kmax_ + 9) / 2 + 2 * n_ + m_ + mxwk0_), mxlws_(nhi_ + kmax_ + mxiwk0_), ws_(mxws_), lws_(mxlws_), k_(0), mode_(COLD_START), iprint_(0), nout_(6), fmin_(-1e20) {
    kktalphac_.alpha = 0;
    /* initialize wsc_ common block (Hessian & workspace for bqpd) */
    wsc_.kk = this->nhr_;
    wsc_.ll = this->nhi_;
    wsc_.mxws = this->mxws_;
    wsc_.mxlws = this->mxlws_;

    // active set
    for (int i = 0; i < this->n_ + this->m_; i++) {
        this->ls_[i] = i + 1;
    }
}

SubproblemSolution BQPDSolver::solve_QP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, CSCMatrix& hessian, std::vector<double>& x) {
    /* Hessian */
    for (int i = 0; i < hessian.number_nonzeros(); i++) {
        this->ws_[i] = hessian.matrix[i];
    }
    /* Hessian sparsity */
    this->lws_[0] = hessian.number_nonzeros() + 1;
    for (int i = 0; i < hessian.number_nonzeros(); i++) {
        this->lws_[i + 1] = hessian.row_number[i] + (hessian.fortran_indexing ? 0 : this->use_fortran);
    }
    for (unsigned int i = 0; i < hessian.column_start.size(); i++) {
        this->lws_[hessian.number_nonzeros() + i + 1] = hessian.column_start[i] + (hessian.fortran_indexing ? 0 : this->use_fortran);
    }

    // if extra variables have been introduced, correct hessian.column_start
    int i = hessian.number_nonzeros() + hessian.column_start.size() + 1;
    int last_value = hessian.column_start[hessian.column_start.size() - 1];
    for (int j = hessian.dimension; j < this->n_; j++) {
        this->lws_[i] = last_value + (hessian.fortran_indexing ? 0 : this->use_fortran);
        i++;
    }

    DEBUG << "hessian: " << hessian;
    return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, x, this->kmax_);
}

SubproblemSolution BQPDSolver::solve_LP(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x) {
    return this->solve_subproblem(variables_bounds, constraints_bounds, linear_objective, constraints_jacobian, x, 0);
}

SubproblemSolution BQPDSolver::solve_subproblem(std::vector<Range>& variables_bounds, std::vector<Range>& constraints_bounds, std::map<int, double>& linear_objective, std::vector<std::map<int, double> >& constraints_jacobian, std::vector<double>& x, int kmax) {

    DEBUG << "gradient obj: "; print_vector(DEBUG, linear_objective);
    for (unsigned int j = 0; j < constraints_jacobian.size(); j++) {
        DEBUG << "gradient c" << j << ": ";
        print_vector(DEBUG, constraints_jacobian[j]);
    }
    for (unsigned int i = 0; i < variables_bounds.size(); i++) {
        DEBUG << "Δx" << i << " in [" << variables_bounds[i].lb << ", " << variables_bounds[i].ub << "]\n";
    }
    for (unsigned int j = 0; j < constraints_bounds.size(); j++) {
        DEBUG << "linearized c" << j << " in [" << constraints_bounds[j].lb << ", " << constraints_bounds[j].ub << "]\n";
    }

    /* Jacobian */
    // TODO preallocate
    std::vector<double> jacobian;
    std::vector<int> jacobian_sparsity;

    jacobian_sparsity.push_back(0); // header-related, to be modified later on
    build_jacobian(jacobian, jacobian_sparsity, linear_objective);
    for (unsigned int j = 0; j < constraints_jacobian.size(); j++) {
        build_jacobian(jacobian, jacobian_sparsity, constraints_jacobian[j]);
    }
    /* Jacobian header */
    jacobian_sparsity[0] = jacobian_sparsity.size();
    unsigned int total_size = 1;
    jacobian_sparsity.push_back(total_size);
    total_size += linear_objective.size();
    jacobian_sparsity.push_back(total_size);
    for (unsigned int j = 0; j < constraints_jacobian.size(); j++) {
        total_size += constraints_jacobian[j].size();
        jacobian_sparsity.push_back(total_size);
    }

    /* bounds */
    for (int i = 0; i < this->n_; i++) {
        this->lb[i] = (variables_bounds[i].lb == -INFINITY) ? -BIG : variables_bounds[i].lb;
        this->ub[i] = (variables_bounds[i].ub == INFINITY) ? BIG : variables_bounds[i].ub;
    }
    for (int j = 0; j < this->m_; j++) {
        this->lb[this->n_ + j] = constraints_bounds[j].lb;
        this->ub[this->n_ + j] = constraints_bounds[j].ub;
    }

    /* call BQPD */
    int mode = (int) this->mode_;
    bqpd_(&this->n_, &this->m_, &this->k_, &kmax, jacobian.data(), jacobian_sparsity.data(), x.data(),
            this->lb.data(), this->ub.data(), &this->f_solution_, &this->fmin_, this->gradient_solution_.data(),
            this->residuals_.data(), this->w_.data(), this->e_.data(), this->ls_.data(), this->alp_.data(),
            this->lp_.data(), &this->mlp_, &this->peq_solution_, this->ws_.data(), this->lws_.data(), &mode,
            &this->ifail_, this->info_.data(), &this->iprint_, &this->nout_);

    /* project solution into bounds: it's a ray! */
    for (unsigned int i = 0; i < x.size(); i++) {
        if (x[i] < variables_bounds[i].lb) {
            x[i] = variables_bounds[i].lb;
        }
        else if (variables_bounds[i].ub < x[i]) {
            x[i] = variables_bounds[i].ub;
        }
    }

    SubproblemSolution solution = this->generate_solution(x);
    return solution;
}

SubproblemSolution BQPDSolver::generate_solution(std::vector<double>& x) {
    Multipliers multipliers(this->n_, this->m_);
    SubproblemSolution solution(x, multipliers);

    /* active constraints */
    for (int j = 0; j < this->n_ - this->k_; j++) {
        int index = std::abs(this->ls_[j]) - this->use_fortran;

        if (this->ls_[j] < 0) { /* upper bound active */
            solution.active_set.at_upper_bound.insert(index);
        }
        else { /* lower bound active */
            solution.active_set.at_lower_bound.insert(index);
        }

        if (index < this->n_) {
            if (this->ls_[j] < 0) { /* upper bound active */
                solution.multipliers.upper_bounds[index] = -this->residuals_[index];
            }
            else { /* lower bound active */
                solution.multipliers.lower_bounds[index] = this->residuals_[index];
            }
        }
        else {
            int constraint_index = index - this->n_;
            solution.constraint_partition.feasible.insert(constraint_index);
            solution.constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
            solution.multipliers.constraints[constraint_index] = (this->ls_[j] < 0) ? - this->residuals_[index] : this->residuals_[index];
        }
    }

    /* inactive constraints */
    for (int j = this->n_ - this->k_; j < this->n_ + this->m_; j++) {
        int index = std::abs(this->ls_[j]) - this->use_fortran;

        if (this->n_ <= index) { // general constraints
            int constraint_index = index - this->n_;
            if (this->residuals_[index] < 0.) { // infeasible constraint
                solution.constraint_partition.infeasible.insert(constraint_index);
                if (this->ls_[j] < 0) { // upper bound violated
                    solution.constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_UPPER;
                }
                else { // lower bound violated
                    solution.constraint_partition.constraint_feasibility[constraint_index] = INFEASIBLE_LOWER;
                }
            }
            else { // feasible constraint
                solution.constraint_partition.feasible.insert(constraint_index);
                solution.constraint_partition.constraint_feasibility[constraint_index] = FEASIBLE;
            }
        }
    }
    solution.status = this->int_to_status(this->ifail_);
    // phase
    // phase_1_required
    solution.norm = norm_inf(x);
    solution.objective = this->f_solution_;
    return solution;
}

Status BQPDSolver::int_to_status(int ifail) {
    if (ifail < 0 || 10 <= ifail) {
        throw std::length_error("BQPDSolver.int_to_status: ifail does not belong to [0, 9]");
    }
    Status status = static_cast<Status>(ifail);
    return status;
}

void BQPDSolver::build_jacobian(std::vector<double>& full_jacobian, std::vector<int>& full_jacobian_sparsity, std::map<int, double>& jacobian) {
    for (std::map<int, double>::iterator it = jacobian.begin(); it != jacobian.end(); it++) {
        int i = it->first;
        double derivative = it->second;

        full_jacobian.push_back(derivative);
        full_jacobian_sparsity.push_back(i + this->use_fortran);
    }
    return;
}
