/*-----------------------------------------------------------*/
/*!

\brief %NOX::NLN implementation of a pseudo transient non-linear
       solver.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_SOLVER_PTC_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_SOLVER_PTC_H_


#include "nox_nln_solver_linesearchbased.H"  // base class of NOX::NLN::Solver::PseudoTransient
#include "nox_nln_abstract_prepostoperator.H"  // base class of NOX::NLN::LinSystem::PrePostOp::PseudoTransient
#include <boost/algorithm/string/predicate.hpp>  // case insensitive string compare

#include <NOX_Abstract_Vector.H>

namespace LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace LINALG
namespace NOX
{
  namespace Epetra
  {
    class Vector;
  }  // namespace  Epetra
  namespace NLN
  {
    class LinearSystem;
    namespace Solver
    {
      /*! \brief Pseudo Transient Continuation (PTC) non-linear solver
       *
       * This implementation is based on the NOX::Solver::PseudoTransient class
       * and the former BACI internal implementation of the PTC method. In contrast
       * to the NOX internal implementation we try to avoid the need of the Thyra
       * interface package.
       *
       * During the line search routine we use always the transient residual [3] for the
       * current iteration point \f$x_{n}\f$
       * \f[
       *   r_{t}(x) = \delta_{n}^{-1} V_{n} (x - x_{n}) + F(x)
       * \f]
       * and a desired merit function (e.g. the sum of squares merit function).
       *
       * See also:
       *
       * [1] C. T. Kelley, D. E. Keyes, "Convergence analysis of pseudo-transient continuation",
       *     SIAM J. Numer. Anal., Vol. 35, No. 2, pp. 508-523, 1998.
       *
       * [2] M. W. Gee, C. T. Kelley, R. B. Lehoucq, "Pseudo-transient continuation for nonlinear
       * transient elasticity", Int. J. Numer. Meth. Engng., Vol. 78, pp. 1209-1219, 2009.
       *
       * [3] M. Ceze, K. J. Fidkowski, "Constrained pseudo-transient continuation",
       *     Int. J. Numer. Meth. Engng., Vol. 102, pp. 1683-1703, 2015.
       *
       * \author Michael Hiermeier */
      class PseudoTransient : public NOX::NLN::Solver::LineSearchBased
      {
       public:
        //! Different pseudo time step control types
        enum TSCType
        {
          tsc_ser,  //!< switched evolution relaxation
          tsc_tte,  //!< temporal truncation error
          tsc_mrr   //!< model reduction ratio
        };

        //! Map pseudo time step control type stl_string to enum
        inline enum TSCType String2TSCType(const std::string& name)
        {
          TSCType type = tsc_ser;
          if (boost::iequals(name, "ser") or boost::iequals(name, "switched evolution relaxation"))
            type = tsc_ser;
          else if (boost::iequals(name, "tte") or boost::iequals(name, "temporal truncation error"))
            type = tsc_tte;
          else if (boost::iequals(name, "mrr") or boost::iequals(name, "model reduction ratio"))
            type = tsc_mrr;
          else
          {
            std::ostringstream msg;
            msg << "Unknown conversion from STL_STRING to TSCType enum for " << name << "."
                << std::endl;
            throwError("String2TSCType", msg.str());
          }
          return type;
        };

        //! Different types of scaling operator V
        enum ScaleOpType
        {
          scale_op_identity,      //!< Use the identity matrix (see [2])
          scale_op_cfl_diagonal,  //!< Use a diagonal matrix based on the local
                                  //!< Courant-Friedrichs-Lewy (CFL) number (see [1]).
          scale_op_lumped_mass,   //!< Use the lumped mass matrix as diagonal scaling operator (not
                                  //!< well tested!)
          scale_op_element_based,
          scale_op_element_based_constant
        };

        enum BuildOpType
        {
          build_op_everyiter,
          build_op_everytimestep
        };

        //! Map build operator type stl_string to enum
        inline enum BuildOpType String2BuildOpType(const std::string& name)
        {
          BuildOpType type = build_op_everyiter;
          if (boost::iequals(name, "every iter"))
            type = build_op_everyiter;
          else if (boost::iequals(name, "every timestep"))
            type = build_op_everytimestep;
          else
          {
            std::ostringstream msg;
            msg << "Unknown conversion from STL_STRING to BuildOperatorType enum for " << name
                << "." << std::endl;
            throwError("String2BuildOpType", msg.str());
          }

          return type;
        };

        //! Map scaling operator type stl_string to enum
        inline enum ScaleOpType String2ScaleOpType(const std::string& name)
        {
          ScaleOpType type = scale_op_identity;
          if (boost::iequals(name, "identity"))
            type = scale_op_identity;
          else if (boost::iequals(name, "cfl diagonal"))
            type = scale_op_cfl_diagonal;
          else if (boost::iequals(name, "lumped mass"))
            type = scale_op_lumped_mass;
          else if (boost::iequals(name, "element based"))
            type = scale_op_element_based;
          else
          {
            std::ostringstream msg;
            msg << "Unknown conversion from STL_STRING to ScaleOperatorType enum for " << name
                << "." << std::endl;
            throwError("String2ScaleOpType", msg.str());
          }

          return type;
        };

       public:
        //! constructor
        PseudoTransient(const Teuchos::RCP<NOX::Abstract::Group>& grp,
            const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
            const Teuchos::RCP<Teuchos::ParameterList>& params);

        //! destructor
        virtual ~PseudoTransient(){};

        //! reset the non-linear solver
        void reset(const NOX::Abstract::Vector& initialGuess,
            const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests);
        virtual void reset(const NOX::Abstract::Vector& initialGuess,
            const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests);

        void createScalingOperator();

        void createLinSystemPrePostOperator();

        void createGroupPrePostOperator();

        //! Evaluated the current nonlinear solver step. (derived)
        virtual NOX::StatusTest::StatusType step();

        //! Solve the stated non-linear problem of the current load/time step (derived)
        virtual NOX::StatusTest::StatusType solve();

        //! Returns the inverse pseudo time step
        const double& getInversePseudoTimeStep() const;

        //! Returns the scaling factor
        const double& getScalingFactor() const;

        //! Returns the scaling operator type
        const enum ScaleOpType& getScalingOperatorType() const;

        //! Returns the scaling diagonal operator
        const Epetra_Vector& getScalingDiagOperator() const;

        //! Returns the usePseudoTransientResidual status
        bool usePseudoTransientResidual() const;

        //! Returns the pseudo transient continuation status
        bool isPtcSolve() const;

       protected:
        //! Initialize the PTC specific variables and call the init function of the base class
        virtual void init();

        //! Print the non-linear solver update
        virtual void printUpdate();

        //! Compute the nodal pseudo velocity for the CFL update
        void computePseudoVelocity();

        /*! \brief Evaluate the model reduction ratio
         *
         *  We use the transient residual and the corresponding (quadratic) model
         * \f$m(\eta_{n}^{*})\f$. If we choose the Sum of Squares merit function, we get
         *
         *  \f[
         *    \rho = \frac{0.5 r_{t}^{T}(x_{n}) r_{t}^{T}(x_{n}) - 0.5 r_{t}^{T}(x_{n} +
         * \eta_{n}^{*} d) r_{t}^{T}(x_{n} + \eta_{n}^{*} d)}{0.5 r_{t}^{T}(x_{n}) r_{t}^{T}(x_{n})
         * - m(\eta_{n}^{*})}, \f] where \f{eqnarray*}{
         *    m(\eta) &=& 0.5 r_{t}^{T}(x_{n}) r_{t}(x_{n}) + \eta d_{n} \nabla_{x}
         * r_{t}(x)|_{x=x_{n}} r_{t}(x_{n}) + 0.5 \eta^{2} d_{n}^{T} \nabla_{x} r_{t}(x)|_{x=x_{n}}
         * \nabla_{x} r_{t}(x)|_{x=x_{n}}^{T} d_{n} \\
         *            &=& 0.5 F_{n}^{T} F_{n} + \eta d_{n} J_{t}^{T} F_{n} + 0.5 \eta^{2} d_{n}^{T}
         * J_{t}^{T} J_{t} d_{n}, \\ \nabla_{x} r_{t}(x)|_{x=x_{n}} &=& J_{t}^{T} = \delta_{n}^{-1}
         * V_{n} + \nabla_{x} F(x)|_{x=x_{n}}. \f}
         *
         *  Please note, that the transient or modified Jacobian is necessary for the evaluation.
         * Therefore, the order of the function calls has to be considered. We have to call the
         * evalModelReductionRatio() before we call the adjustPseudoTimeStep() function! */
        void evalModelReductionRatio();

        /*! \brief Adjust the pseudo time step if line search was active and the step length has
         * been changed
         *
         *  We use a least squares approximation based on the following idea.
         *  The pseudo time step is modified in such a way, that it represents in a least squares
         *  sense the actual used step, i.e. the direction vector scaled by the step length
         *  parameter \f$\eta_{n}^{*}\f$. The corresponding pseudo time step is called
         *  \f$\delta_{n}^{*}\f$:
         *
         *  \f[
         *    [(\delta_{n}^{*})^{-1} V_{n} + \nabla_{x} F_{n}^{T}] (x(t_{n}+\delta^{*}_{n}) - x_{n})
         * = -F_{n}, \f]
         *
         *  where we replace \f$(x(t_{n}+\delta^{*}_{n})\f$ by \f$\eta_{n}^{*} d_{n}\f$ and
         * reformulate the equation
         *
         *  \f[
         *    V_{n} d_{n} = - \delta_{n}^{*} ((\eta_{n}^{*})^{-1} F_{n} + \nabla_{x} F_{n}^{T}
         * d_{n}). \f]
         *
         *  Since we cannot find a suitable scalar \f$\delta_{n}^{*}\f$, which solves the problem
         * exactly, we formulate a least squares problem to minimize the error:
         *
         *  \f[
         *    \min\limits_{\delta_{n}} 0.5 \| V_{n} d_{n} - \delta_{n} ((\eta_{n}^{*})^{-1} F_{n} +
         * \nabla_{x} F_{n}^{T} d_{n}) \|^{2} \f]
         *
         *  The approximated pseudo time step is
         *  \f[
         *    \delta^{*}_{n} \approx - \frac{d_{n}^{T} V_{n} ((\eta_{n}^{*})^{-1} F_{n} + \nabla_{x}
         * F_{n}^{T} d_{n})}{\|(\eta_{n}^{*})^{-1} F-{n} + \nabla F_{n}^{T} d_{n} \|^{2}}. \f]
         *
         *  Please note, that we undo the changes to the jacobian, if the step length is not
         * equal 1.0! */
        void adjustPseudoTimeStep();

        /*! \brief Update the pseudo time step
         *
         *  There are three different options:
         *  - The switched evolution relaxation\n
         *    (in the most simple case)
         *    \f[
         *       \delta_{n} = \delta_{0}
         * \frac{\|\underline{F}(\underline{x}_{0})\|}{\|\underline{F}(\underline{x}_{n})\|} \f]
         *  - The temporal truncation error\n
         *
         *       MISSING AT THE MOMENT
         *
         *  - The model reduction ratio.\n
         *    See the evalModelReductionRatio() routine for more detail on how to evaluate the
         * ratio.
         *
         *    If \f$\rho\f$ is smaller than 0.2, the pseudo time step is decreased,
         *    else if \f$\rho\f$ is larger than 0.9, the pseudo time step is increased,
         *    otherwise we keep the current pseudo time step unchanged. */
        void updatePseudoTimeStep();

        //! Returns true if the scaling operator has been evaluated.
        bool IsScalingOperator() const { return isScalingOperator_; };

       private:
        //! Throw class specific error
        void throwError(const std::string& functionName, const std::string& errorMsg) const;

       protected:
        //! Inner Stopping test
        Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> iTestPtr_;

        /*! Pointer to a linear system pre/post operator, which is used
         *  to modify the jacobian directly from the non-linear solver. */
        Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prePostLinSysPtr_;

        /*! Pointer to a nln group pre/post operator, which is used
         *  to modify the right-hand-side directly from the non-linear solver. */
        Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prePostGroupPtr_;

        //! Pointer to the scaling Operator (Identity scaling)
        Teuchos::RCP<Epetra_Vector> scalingDiagOpPtr_;
        //! Pointer to the scaling Operator (Element based scaling)
        Teuchos::RCP<LINALG::SparseMatrix> scalingMatrixOpPtr_;

        //! @name Special pseudo transient continuation parameters
        //@{
        //! Pseudo step size for pseudo-transient stepping.
        double delta_;
        //! Inverse pseudo time step size for pseudo-transient stepping.
        double invDelta_;
        //! Initial pseudo time step size.
        double deltaInit_;
        //! Maximum pseudo time step size.
        double deltaMax_;
        //! Minimum pseudo time step size.
        double deltaMin_;
        //! Pseudo time step size from previous iteration.
        double deltaOld_;
        //! Pseudo-transient time.
        double pseudoTime_;
        //! Exponent for SER.
        double SER_alpha_;
        //! Scaling factor for modified Jacobian, i.e. the scaling operator.
        double scaleFactor_;

        //! solution time derivative used for scaling operator V in pseudo-transient paper.
        Teuchos::RCP<NOX::Abstract::Vector> xDot_;

        //! Are we solving the pseudo transient continuation problem at the moment?
        bool isPtcSolve_;

        /** \brief If set to true, the computeF function is used to modify the residual during the
         * line search call.
         *
         *  The pseudo transient residual is defined as
         *  \f[
         *    r_{t}(\eta) = V_{n} \delta_{n}^{-1} [x_{n} + \eta d_{n} - x_{n}] + F(x_{n} + \eta
         * d_{n}) \f]
         */
        bool usePseudoTransientResidual_;

        //! if it's true, we compute automatically an initial pseudo time step.
        bool calcDeltaInit_;

        //! Has the scaling operator been evaluated, yet?
        bool isScalingOperator_;

        //! Maximum number of iterations before pseudo-transient is disabled and the algorithm
        //! switches to a line search-based direct to steady state solve.
        int maxPseudoTransientIterations_;

        //! time step control type
        enum TSCType tscType_;

        //! scaling operator type
        enum ScaleOpType scaleOpType_;

        //! scaling operator type
        enum BuildOpType build_scaling_op_;

        //! vector norm type (necessary for the time step control)
        enum NOX::Abstract::Vector::NormType normType_;

        //! The current model reduction ratio.
        double modelReductionRatio_;
        //@}

      };  // class PseudoTransient
    }     // namespace Solver
    namespace LinSystem
    {
      namespace PrePostOp
      {
        /*! \brief Pseudo Transient Continuation (PTC) non-linear solver helper class
         *
         * This class implementation is an example for the usage of the
         * NOX::NLN::Abstract::PrePostOperator and its wrapper the
         * NOX::NLN::LinSystem::PrePostOperator classes. The PTC solver uses them to modify the
         * linear system or to be even more precise the jacobian.
         *
         * \author Michael Hiermeier */
        class PseudoTransient : public NOX::NLN::Abstract::PrePostOperator
        {
         public:
          //! constructor
          PseudoTransient(Teuchos::RCP<Epetra_Vector>& scalingDiagOp,
              Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr,
              const NOX::NLN::Solver::PseudoTransient& ptcsolver);

          //! destructor
          virtual ~PseudoTransient(){};

          //! derived
          virtual void runPostComputeJacobian(LINALG::SparseOperator& jac, const Epetra_Vector& x,
              const NOX::NLN::LinearSystem& linsys);

          //! derived
          virtual void runPostComputeFandJacobian(Epetra_Vector& rhs, LINALG::SparseOperator& jac,
              const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys);

         protected:
          //! modify the jacobian as defined by the scaling operator type
          virtual void modifyJacobian(LINALG::SparseMatrix& jac);

         private:
          //! read-only access
          const NOX::NLN::Solver::PseudoTransient& ptcsolver_;

          //! reference to the scaling diagonal operator
          Teuchos::RCP<Epetra_Vector>& scalingDiagOpPtr_;

          //! reference to the scaling diagonal operator matrix
          Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr_;

        };  // class PseudoTransient
      }     // namespace PrePostOp
    }       // namespace LinSystem
    namespace GROUP
    {
      namespace PrePostOp
      {
        /*! \brief Pseudo Transient Continuation (PTC) non-linear solver helper class
         *
         * The PTC solver uses this helper class to modify the right hand side.
         *
         * \author Michael Hiermeier */
        class PseudoTransient : public NOX::NLN::Abstract::PrePostOperator
        {
         public:
          //! constructor
          PseudoTransient(Teuchos::RCP<Epetra_Vector>& scalingDiagOpPtr,
              Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr,
              const NOX::NLN::Solver::PseudoTransient& ptcsolver);

          //! destructor
          virtual ~PseudoTransient(){};

          //! derived
          virtual void runPreComputeF(Epetra_Vector& F, const NOX::NLN::Group& grp);

          //! derived
          virtual void runPostComputeF(Epetra_Vector& F, const NOX::NLN::Group& grp);

         protected:
          Teuchos::RCP<NOX::Epetra::Vector> evalPseudoTransientFUpdate(const NOX::NLN::Group& grp);

         private:
          //! read-only access
          const NOX::NLN::Solver::PseudoTransient& ptcsolver_;

          //! reference to the scaling diagonal operator
          Teuchos::RCP<Epetra_Vector>& scalingDiagOpPtr_;

          //! reference to the scaling diagonal operator matrix
          Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr_;

          //! reference to the flag which indicates, if the current right-hand-side is the transient
          //! right-hand-side
          bool isPseudoTransientResidual_;

        };  // class PseudoTransient
      }     // namespace PrePostOp
    }       // namespace GROUP
  }         // namespace NLN
}  // namespace NOX

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_SOLVER_PTC_H_ */
