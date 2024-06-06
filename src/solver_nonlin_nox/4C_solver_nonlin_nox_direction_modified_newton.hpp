/*----------------------------------------------------------------------------*/
/*! \file

\brief Implementation of a modified Newton approach



\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_MODIFIED_NEWTON_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_DIRECTION_MODIFIED_NEWTON_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_direction_newton.hpp"
#include "4C_solver_nonlin_nox_floating_point_exception.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    enum class CorrectionType : int;
    namespace Direction
    {
      namespace Test
      {
        class Generic;
      }  // namespace Test

      /** \brief Modified Newton method
       *
       *  If the system matrix shows bad properties and does not pass the applied
       *  default step tests, it will be modified by adding a scaled vector to
       *  the diagonal of the primal-primal block. Note that this improves
       *  in general simultaneously the linear solver performance, since the matrix
       *  becomes more diagonal dominant.
       *
       *  After a certain number of iterates which all asked for a decreasing
       *  regularization of the system matrix block, the here proposed algorithm will
       *  switch automatically back to the unmodified system. In this way
       *  quadratic convergence is maintained whenever possible.
       *
       *  \author hiermeier */
      class ModifiedNewton : public Newton
      {
       public:
        /// constructor
        ModifiedNewton(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& p);


        bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
            const ::NOX::Solver::Generic& solver) override;

       private:
        /** \brief compute the direction in case of a correction step
         *
         *  In difference to a default attempt the correction parameter is not
         *  changed but instead used as it is. If the solver fails, the correction
         *  fails and the line search method takes care. */
        bool compute_correction_direction(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
            const ::NOX::Solver::Generic& solver, NOX::Nln::CorrectionType corr_type);

        /** \brief update the successive correction counter
         *
         *  Increase the successive reduction counter if only one correction has been
         *  necessary. This is equivalent to the case that the primal correction factor
         *  has been reduced. If this happens N-times in a row, we try to use the
         *  unmodified system matrix again. */
        void update_successive_reduction_counter();

        /// return true if the system shall be modified
        bool use_unmodified_system() const;

        /// fill the default step test set
        void fill_default_step_tests(Teuchos::ParameterList& pmodnewton);

        /// test the default step, i.e. the (modified) Newton direction of length 1.0
        bool test_default_step_quality(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
            Teuchos::RCP<Epetra_Vector>& diagonal, bool first_test = false);

        /// compute the modified Newton direction
        bool compute_modified_newton(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
            const ::NOX::Solver::Generic& solver, Epetra_Vector* diagonal = nullptr);

        /// solve the unmodified system
        bool solve_unmodified_system(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
            const ::NOX::Solver::Generic& solver);

        /// solve the modified system
        bool solve_modified_system(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
            const ::NOX::Solver::Generic& solver);

        /// detect a stagnation, update the counter and inform the user if necessary
        void set_stagnation_counter(::NOX::Abstract::Vector& dir);

        /// modify the system matrix
        bool modify_system(::NOX::Abstract::Group& grp, Epetra_Vector* diagonal);
        bool modify_system(
            ::NOX::Abstract::Group& grp, Epetra_Vector* diagonal, const double primal_diag_corr);

        /// return the primal diagonal correction
        double get_primal_diag_correction(const bool first) const;

        /** return the very first primal diagonal correction, i.e. this is the first iterate asking
         *  for a modification */
        double get_first_primal_diag_correction() const;

        /// return default primal diagonal correction
        double get_primal_diag_correction() const;

        /// store correction factor as soon as a successful modification could be achieved
        void store_correction_factor();

        /// get diagonal vector (currently unused)
        Teuchos::RCP<Epetra_Vector> get_diagonal(const ::NOX::Abstract::Group& grp) const;

        /// print info to stream
        void print(std::ostream& os, const NOX::Nln::CorrectionType* corr_type = nullptr) const;

        /// throw error which can be caught
        void throw_error(
            const int line, const std::string& functionName, const std::string& errorMsg) const;

       private:
        /// object to handle floating point exceptions
        FloatingPointException fp_except_;

        /// correction factor for the diagonal of the primal field block
        double primal_diag_corr_ = 0.0;

        /** \brief lastly used and accepeted diagonal correction factor
         *
         *  i.e. stemming from a previous Newton iterate which asks for a correction as well) */
        double primal_diag_corr_last_ = 0.0;

        /// number of corrections (pointer to the value in the parameter list)
        int* corr_counter_ = nullptr;

        /// count successive reduction
        int successive_red_counter_ = 1000;

        /// count stagnations (shouldn't happen)
        unsigned stagnation_counter_ = 0;

        /// initial correction factor for the diagonal of the primal block
        const double init_primal_diag_corr_;

        /// minimal correction factor for the diagonal of the primal block
        const double min_primal_diag_corr_;

        /// maximal correction factor for the diagonal of the primal block
        const double max_primal_diag_corr_;

        /// reduction factor for the adaption of the primal diagonal correction
        const double primal_red_fac_;

        /// accretion factor for the adaption of the primal diagonal correction
        const double primal_acc_fac_;

        /// high accretion factor for the adaption of the primal diagonal correction
        const double primal_high_acc_fac_;

        /// see base class for more info
        bool use_adjustable_forcing_term_ = false;

        /// Unmodified diagonal of the primal system matrix block
        Teuchos::RCP<const Epetra_Vector> original_diag_ptr_ = Teuchos::null;

        //! NOX_Utils pointer
        Teuchos::RCP<::NOX::Utils> utils_;

        /// pointer to the parameter list
        Teuchos::ParameterList* params_;

        /// set of direction tests
        std::vector<Teuchos::RCP<NOX::Nln::Direction::Test::Generic>> dstests_;
      };
    }  // namespace Direction
  }    // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
