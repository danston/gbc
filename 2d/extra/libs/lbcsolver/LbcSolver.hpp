// Authors: Juyong Zhang, Bailin Deng, Zishun Liu, Giuseppe Patane, Sofien Bouaziz, Kai Hormann, and Ligang Liu, https://github.com/bldeng/LBC

#ifndef GBC_LBCSOLVER_HPP
#define GBC_LBCSOLVER_HPP

// STL includes.
#include <vector>
#include <iostream>
#include <algorithm>

// OMP includes.
// #include <omp.h>

// Local includes.
#include "DataTypes.hpp"
#include "DataSetup.hpp"

namespace LBC {

    struct Param {
        double penalty_weight;

        // The relaxation coefficient for ADMM.
        double relaxation_alpha;

        // Relative convergence threshold for primal and dual residuals.
        double rel_primal_eps, rel_dual_eps;

        // Absolute threshold for primal and dual residuals.
        double abs_eps;

        // Max iterations.
        int max_iterations;

        // How often should we check the convergence of the solvers?
        int convergence_check_frequency;

        // How often should we output the progress of the solver?
        // This variable represents the ratio between the output frequency and the convergence check frequency.
        int output_frequency_ratio;

        // Timer.
        bool use_timer;

        // Constructor that sets the default values.
        Param() :
                penalty_weight(10), relaxation_alpha(1.0), rel_primal_eps(1e-6), rel_dual_eps(1e-6),
                abs_eps(1e-8), max_iterations(10000), convergence_check_frequency(10),
                output_frequency_ratio(10),
                use_timer(true) { }
    };

    class LBCSolver {
    public:
        LBCSolver(const Param &param, const DataSetup &data_setup) :
                param_(param), control_points_(data_setup.get_LBC_solver_control_points()),
                data_points_(data_setup.get_LBC_solver_data_points()),
                grad_weights_(data_setup.get_LBC_solver_grad_weights()),
                G_(data_setup.get_LBC_solver_grad_operator()), E_(data_setup.get_LBC_solver_grad_const()),
                has_init_coord_(false), start_time_(0.0), end_time_(0.0), verbose_(false) {

            init();
        }

        virtual ~LBCSolver() { }

        void initialize_coordinates(const DenseMatrix &coord) {

            init_coord_ = coord;
            if (init_coord_.rows() != n_data_points_ || init_coord_.cols() != n_control_points_) {
                std::cerr << "Invalid size of initial coordinates" << std::endl;
                valid_init_data_ = false;
            }

            has_init_coord_ = true;
        }

        const DenseMatrix &get_coordinates() const {
            return Z_;
        }

        void solve() {

            if (!valid_init_data_) {
                std::cerr << "Invalid data, unable to solve...." << std::endl;
                return;
            }

            initialize_variables();

            start_timer();
            int iter = 0;
            optimization_end_ = false;

            while (!optimization_end_) {

                iter++;
                check_convergence_ = (iter % convergence_check_frequency_ == 0 || iter >= param_.max_iterations);
                output_progress_ = (iter % output_frequency_ == 0);

#pragma omp parallel
                {
                    this->update_x();

                    this->update_w();

                    this->update_y();

                    this->update_dual_variables(iter);
                }
            }

            end_timer();
            show_elapsed_time();
        }

        inline void showOutput(const bool verbose) {
            verbose_ = verbose;
        }

    protected:

        // Solver parameters.
        Param param_;

        // Positions of control points and sample data points.
        // The number of rows is the same as the problem dimension.
        DenseMatrix control_points_, data_points_;

        // Gradient weights on the cells, of size N_f x N_c, where N_f is the number of cells,
        // N_c is the number of control points.
        DenseMatrix grad_weights_;

        // G: Gradient operator for cells, with size (d N_f) x N_v, where N_f is the number of cells,
        // N_v is the number of non-boundary vertices, and d is the problem dimension.
        // Gt is the transpose of G.
        SparseMatrix G_, Gt_;

        // Constant factor in the face gradient computation.
        DenseMatrix E_;

        // Other variables used in the solver (see the LBC paper for their definitions).
        DenseMatrix W_, X_, Y_, M_, Mt_, H_, J_;

        // Relaxation variables for W and X.
        DenseMatrix W_h_, X_h_;

        // Scaled dual variables.
        DenseMatrix d1_, d2_;

        // The coordinate values computed from Y: Z = Y * M + H.
        DenseMatrix Z_, prev_Z_;

        // Temporary storage for face gradient values.
        DenseMatrix Grad_, prev_Grad_;

        // Initial values for the coordinate values.
        DenseMatrix init_coord_;
        bool has_init_coord_;

        // Temporary storage for computing the rhs of the linear system for Y.
        DenseMatrix pre_rhs_;

        // Boolean variables for the progress of the solver.
        bool optimization_converge_, optimization_end_;

        // Timer variables.
        double start_time_, end_time_;

        // Verbose.
        bool verbose_;

        // Are the initial data valid?
        bool valid_init_data_;

        // Problem dimension.
        int dim_;

        int n_control_points_;
        int n_data_points_;
        int n_cells_;
        int n_total_coef_;
        int n_grad_elements_;
        int n_Y_col_;

        bool check_convergence_;
        bool output_progress_;
        int output_frequency_, convergence_check_frequency_;

        // Variables for primal and dual residuals.
        DenseMatrix primal_residual_X_, primal_residual_W_;
        double primal_residual_X_sqr_norm_, primal_residual_W_sqr_norm_;
        double dual_residual_X_sqr_norm_, dual_residual_W_sqr_norm_;
        double primal_residual_sqr_norm_, dual_residual_sqr_norm_;
        double primal_residual_sqr_norm_threshold_, dual_residual_sqr_norm_threshold_;

        // Cholesky solver for the linear system.
#ifdef USE_CHOLMOD
        Eigen::CholmodDecomposition<SparseMatrix> solver_;
#else
        Eigen::SimplicialLDLT<SparseMatrix> solver_;
#endif

        // Update steps for variables W, X, Y.
        void update_w() {

#pragma omp for
            for (int i = 0; i < n_total_coef_; ++i) {

                const int v_index = i % n_data_points_;
                const int c_index = i / n_data_points_;

                W_(v_index, c_index) = std::min(1.0, std::max(0.0, W_(v_index, c_index)));
            }
        }

        void update_x() {

#pragma omp for
            for (int i = 0; i < n_grad_elements_; ++i) {

                const int f_index = i % n_cells_;
                const int r_index = dim_ * f_index;
                const int c_index = i / n_cells_;

                double a = param_.penalty_weight * X_.block(r_index, c_index, dim_, 1).norm();
                double current_weight = grad_weights_(f_index, c_index);

                if (a <= current_weight) {
                    X_.block(r_index, c_index, dim_, 1).setZero();
                } else {
                    X_.block(r_index, c_index, dim_, 1) *= (1.0 - current_weight / a);
                }
            }
        }

        void update_y() {

#pragma omp for
            for (int i = 0; i < n_control_points_; ++i) {
                X_h_.col(i) = X_.col(i) * param_.relaxation_alpha - Grad_.col(i) * (param_.relaxation_alpha - 1);
                W_h_.col(i) = W_.col(i) * param_.relaxation_alpha - Z_.col(i) * (param_.relaxation_alpha - 1);
                pre_rhs_.col(i) = Gt_ * (X_h_.col(i) - J_.col(i) + d1_.col(i)) + W_h_.col(i) - H_.col(i) + d2_.col(i);
            }

#pragma omp for
            for (int i = 0; i < n_Y_col_; ++i) {
                Y_.col(i) = solver_.solve(pre_rhs_ * Mt_.col(i));
            }
        }

        // Initialize primal and dual thresholds for solver convergence.
        void initialize_thresholds() {

            const double primal_threshold = std::max(param_.abs_eps,
                                               dim_ * n_cells_ * n_control_points_ * param_.rel_primal_eps);

            const double dual_threshold = std::max(param_.abs_eps,
                                             n_data_points_ * n_control_points_ * param_.rel_dual_eps);

            primal_residual_sqr_norm_threshold_ = primal_threshold * primal_threshold;
            dual_residual_sqr_norm_threshold_ = dual_threshold * dual_threshold;
        }

        // Pre-process the linear equality constraints, to compute matrices M and H.
        void initialize_linear_constraint_data() {

            DenseMatrix Kt, Bt;
            Kt.resize(dim_ + 1, n_control_points_);
            Bt.resize(dim_ + 1, n_data_points_);

            Kt.block(0, 0, dim_, n_control_points_) = control_points_;
            Kt.row(dim_).fill(1.0);
            Bt.block(0, 0, dim_, n_data_points_) = data_points_;
            Bt.row(dim_).fill(1.0);

            Eigen::JacobiSVD<DenseMatrix, Eigen::FullPivHouseholderQRPreconditioner>
                    jsvd(Kt, Eigen::ComputeFullU | Eigen::ComputeFullV);

            const int nrank = (int) jsvd.nonzeroSingularValues();

            Mt_ = jsvd.matrixV().block(0, nrank, n_control_points_, n_control_points_ - nrank);
            M_ = Mt_.transpose();
            n_Y_col_ = (int) M_.rows();
            H_ = jsvd.solve(Bt).transpose(); // least squares solving
        }

        // Pre-factorize the linear system matrix for Y.
        void initialize_solver() {

            // Construct identity matrix for data point size.
            SparseMatrix Id(n_data_points_, n_data_points_);

            std::vector<TripletD> triplets_I;
            for (int i = 0; i < n_data_points_; i++) {
                triplets_I.push_back(TripletD(i, i, 1.0));
            }
            Id.setFromTriplets(triplets_I.begin(), triplets_I.end());

            SparseMatrix M = Gt_ * G_ + Id;
            solver_.compute(M);
            if (solver_.info() != Eigen::Success) {
                std::cerr << "Linear system matrix factorization failed!" << std::endl;
                valid_init_data_ = false;
            }
        }

        // Initialization for variables used in the solver.
        void initialize_variables() {

            if (has_init_coord_) {
                Z_ = init_coord_;
            } else {
                Z_.resize(n_data_points_, n_control_points_);
                Z_.fill(1.0 / n_control_points_);
            }

            Y_ = (Z_ - H_) * Mt_;

            W_ = Z_;
            Grad_ = G_ * W_ + E_;
            J_ = G_ * H_ + E_;
            X_ = Grad_;
            X_h_ = X_;
            W_h_ = W_;

            primal_residual_W_.setZero(W_.rows(), W_.cols());
            primal_residual_X_.setZero(X_.rows(), X_.cols());

            pre_rhs_.setZero(H_.rows(), H_.cols());

            d1_.setZero(X_.rows(), X_.cols());
            d2_.setZero(W_.rows(), W_.cols());
        }

        void init() {

            valid_init_data_ = true;

            dim_ = (int) control_points_.rows();
            if (data_points_.rows() != dim_) {
                std::cerr << "Dimension mismatch between control points and data points!" << std::endl;
                valid_init_data_ = false;
            }

            n_control_points_ = (int) control_points_.cols();
            n_data_points_ = (int) data_points_.cols();
            n_total_coef_ = n_control_points_ * n_data_points_;
            if (n_control_points_ == 0 || n_data_points_ == 0) {
                std::cerr << "Invalid number of control points or data poins!" << std::endl;
                valid_init_data_ = false;
            }

            Gt_ = G_.transpose();
            n_cells_ = G_.rows() / dim_;
            n_grad_elements_ = n_cells_ * n_control_points_;
            if (G_.rows() % dim_) {
                std::cerr << "Invalid dimension for the gradient operator!" << std::endl;
                valid_init_data_ = false;
            }
            if (n_cells_ == 0) {
                std::cerr << "Invalid number of cells!" << std::endl;
                valid_init_data_ = false;
            }

            if (grad_weights_.rows() != n_cells_ || grad_weights_.cols() != n_control_points_) {
                std::cerr << "Invalid dimension of gradient weights!" << std::endl;
                valid_init_data_ = false;
            }

            convergence_check_frequency_ = std::max(1, param_.convergence_check_frequency);
            output_frequency_ = std::max(1, param_.output_frequency_ratio * convergence_check_frequency_);

            initialize_solver();
            initialize_linear_constraint_data();
            initialize_thresholds();
        }

        // Update the dual variables and check if the solver converges.
        void update_dual_variables(int iter_num) {

#pragma omp sections
            {
#pragma omp section
                {
                    if (check_convergence_) {
                        prev_Z_ = Z_;
                    }
                }

#pragma omp section
                {
                    if (check_convergence_) {
                        prev_Grad_ = Grad_;
                    }
                }

#pragma omp section
                {
                    primal_residual_X_ = X_;
                }

#pragma omp section
                {
                    primal_residual_W_ = W_;
                }
            }

#pragma omp for
            for (int i = 0; i < n_control_points_; i++) {
                Z_.col(i) = Y_ * M_.col(i) + H_.col(i);
            }

#pragma omp for
            for (int i = 0; i < n_control_points_; i++) {
                Grad_.col(i) = G_ * Z_.col(i) + E_.col(i);
            }

#pragma omp sections
            {
#pragma omp section
                {
                    if (check_convergence_) {
                        dual_residual_X_sqr_norm_ = (Grad_ - prev_Grad_).squaredNorm();
                    }
                }

#pragma omp section
                {
                    if (check_convergence_) {
                        dual_residual_W_sqr_norm_ = (Z_ - prev_Z_).squaredNorm();
                    }
                }

#pragma omp section
                {
                    if (check_convergence_) {
                        primal_residual_X_ -= Grad_;
                        primal_residual_X_sqr_norm_ = primal_residual_X_.squaredNorm();
                    }
                }

#pragma omp section
                {
                    if (check_convergence_) {
                        primal_residual_W_ -= Z_;
                        primal_residual_W_sqr_norm_ = primal_residual_W_.squaredNorm();
                    }
                }

#pragma omp section
                {
                    d1_ += (X_h_ - Grad_);
                    X_ = Grad_ - d1_;
                }

#pragma omp section
                {
                    d2_ += (W_h_ - Z_);
                    W_ = Z_ - d2_;
                }
            }

#pragma omp single
            {
                if (check_convergence_) {

                    primal_residual_sqr_norm_ = primal_residual_X_sqr_norm_ + primal_residual_W_sqr_norm_;
                    dual_residual_sqr_norm_   = (dual_residual_X_sqr_norm_ + dual_residual_W_sqr_norm_)
                                              * param_.penalty_weight * param_.penalty_weight;

                    optimization_converge_ = (primal_residual_sqr_norm_  <= primal_residual_sqr_norm_threshold_
                                              && dual_residual_sqr_norm_ <= dual_residual_sqr_norm_threshold_);

                    optimization_end_ = optimization_converge_ || iter_num >= param_.max_iterations;

                    if (optimization_converge_ && verbose_) {

                        std::cout << "Solver converged!\n" << std::endl;

                    } else if (optimization_end_ && verbose_) {

                        std::cout << "Maximum iteration reached!\n" << std::endl;
                    }

                    if ((output_progress_ || optimization_end_) && verbose_) {

                        std::cout << "Iteration " << iter_num << ":" << std::endl;

                        std::cout << "Primal residual squared norm: " << primal_residual_sqr_norm_ <<
                        ",  primal threshold:" << primal_residual_sqr_norm_threshold_ << std::endl;

                        std::cout << "Dual residual squared norm: " << dual_residual_sqr_norm_ <<
                        ",  dual threshold:" << dual_residual_sqr_norm_threshold_ << std::endl;

                        std::cout << std::endl;
                    }
                }
            }
        }

        // Timer methods.
        void start_timer() {
            if (param_.use_timer) {
                start_time_ = 0.0; // omp_get_wtime();
            }
        }

        void end_timer() {
            if (param_.use_timer) {
                end_time_ = 0.0; // omp_get_wtime();
            }
        }

        void show_elapsed_time() {
            if (param_.use_timer) {
                std::cout << "Solving time: " << std::max(0.0, end_time_ - start_time_) << " seconds." << std::endl;
            }
        }
    };

} // namespace LBC

#endif // GBC_LBCSOLVER_HPP
