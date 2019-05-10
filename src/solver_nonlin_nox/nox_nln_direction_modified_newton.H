/*----------------------------------------------------------------------------*/
/*!
\file nox_nln_direction_modified_newton.H

\brief Implementation of a modified Newton approach

\maintainer Michael Hiermeier

\date Aug 15, 2018

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_DIRECTION_MODIFIED_NEWTON_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_DIRECTION_MODIFIED_NEWTON_H_

#include "nox_nln_direction_newton.H"
#include "nox_nln_floating_point_exception.H"

class Epetra_Vector;
namespace NOX
{
  class Utils;
  namespace NLN
  {
    enum class CorrectionType : int;
    namespace Direction
    {
      namespace Test
      {
        class Generic;
        ;
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
        ModifiedNewton(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& p);

        /// destructor
        virtual ~ModifiedNewton(){};

        // derived
        bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            const NOX::Solver::Generic& solver) override;

       private:
        /** \brief compute the direction in case of a correction step
         *
         *  In difference to a default attempt the correction parameter is not
         *  changed but instead used as it is. If the solver fails, the correction
         *  fails and the line search method takes care. */
        bool computeCorrectionDirection(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            const NOX::Solver::Generic& solver, NOX::NLN::CorrectionType corr_type);

        /** \brief update the successive correction counter
         *
         *  Increase the successive reduction counter if only one correction has been
         *  necessary. This is equivalent to the case that the primal correction factor
         *  has been reduced. If this happens N-times in a row, we try to use the
         *  unmodified system matrix again. */
        void updateSuccessiveReductionCounter();

        /// return true if the system shall be modified
        bool useUnmodifiedSystem() const;

        /// fill the default step test set
        void fillDefaultStepTests(Teuchos::ParameterList& pmodnewton);

        /// test the default step, i.e. the (modified) Newton direction of length 1.0
        bool testDefaultStepQuality(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            Teuchos::RCP<Epetra_Vector>& diagonal, bool first_test = false);

        /// compute the modified Newton direction
        bool computeModifiedNewton(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            const NOX::Solver::Generic& solver, Epetra_Vector* diagonal = NULL);

        /// solve the unmodified system
        bool solveUnmodifiedSystem(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            const NOX::Solver::Generic& solver);

        /// solve the modified system
        bool solveModifiedSystem(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
            const NOX::Solver::Generic& solver);

        /// detect a stagnation, update the counter and inform the user if necessary
        void setStagnationCounter(NOX::Abstract::Vector& dir);

        /// modify the system matrix
        bool modifySystem(NOX::Abstract::Group& grp, Epetra_Vector* diagonal);
        bool modifySystem(
            NOX::Abstract::Group& grp, Epetra_Vector* diagonal, const double primal_diag_corr);

        /// return the primal diagonal correction
        double getPrimalDiagCorrection(const bool first) const;

        /** return the very first primal diagonal correction, i.e. this is the first iterate asking
         *  for a modification */
        double getFirstPrimalDiagCorrection() const;

        /// return default primal diagonal correction
        double getPrimalDiagCorrection() const;

        /// store correction factor as soon as a successful modification could be achieved
        void storeCorrectionFactor();

        /// get diagonal vector (currently unused)
        Teuchos::RCP<Epetra_Vector> getDiagonal(const NOX::Abstract::Group& grp) const;

        /// print info to stream
        void print(std::ostream& os, const NOX::NLN::CorrectionType* corr_type = NULL) const;

        /// throw error which can be caught
        void throwError(
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
        int* corr_counter_ = NULL;

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
        bool useAdjustableForcingTerm_ = false;

        /// Unmodified diagonal of the primal system matrix block
        Teuchos::RCP<const Epetra_Vector> original_diag_ptr_ = Teuchos::null;

        //! NOX_Utils pointer
        Teuchos::RCP<NOX::Utils> utils_;

        /// pointer to the parameter list
        Teuchos::ParameterList* params_;

        /// set of direction tests
        std::vector<Teuchos::RCP<NOX::NLN::Direction::Test::Generic>> dstests_;
      };
    }  // namespace Direction
  }    // namespace NLN
}  // namespace NOX



#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_DIRECTION_MODIFIED_NEWTON_H_ */
