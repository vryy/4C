/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_ptc.H

\maintainer Michael Hiermeier

\date Oct 6, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_SOLVER_PTC_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_SOLVER_PTC_H_


#include "nox_nln_solver_linesearchbased.H"   // base class of NOX::NLN::Solver::PseudoTransient
#include "nox_nln_linearsystem_prepostoperator_generic.H"   // base class of NOX::NLN::LinSystem::PrePostOp::PseudoTransient
#include <boost/algorithm/string/predicate.hpp>   // case insensitive string compare

#include <NOX_Abstract_Vector.H>

namespace LINALG {
  class SparseOperator;
  class SparseMatrix;
} // namespace LINALG
namespace NOX {
namespace NLN {
  class LinearSystem;
namespace LinSystem {
namespace PrePostOp {
  class Generic;
} // namespace PrePostOp
} // namespace LinSystem
namespace Solver {
/*!
  \brief Pseudo Transient Continuation (PTC) non-linear solver

  This implementation is based on the NOX::Solver::PseudoTransient class
  and the former BACI internal implementation of the PTC method. In contrast
  to the NOX internal implementation we try to avoid the need of the Thyra
  interface package.
  A further modification is the usage of the internal line search call: Instead of
  using the new trial point after the line search method has already been applied,
  we use the trial point with a default step length of one for the update of the
  pseudo time step of the PTC method.

  See also:

  [1] C. T. Kelley, D. E. Keyes, "Convergence analysis of pseudo-transient continuation",
      SIAM J. Numer. Anal., Vol. 35, No. 2, pp. 508-523, 1998.

  [2] M. W. Gee, C. T. Kelley, R. B. Lehoucq, "Pseudo-transient continuation for nonlinear transient elasticity",
      Int. J. Numer. Meth. Engng., Vol. 78, pp. 1209-1219, 2009.

  \author Michael Hiermeier
 */
  class PseudoTransient : public NOX::NLN::Solver::LineSearchBased
  {
  public:
    //! Different pseudo time step control types
    enum TSCType
    {
      tsc_ser,      //!< switched evolution relaxation
      tsc_tte       //!< temporal truncation error
    };

    //! Map pseudo time step control type stl_string to enum
    inline enum TSCType String2TSCType(const std::string& name)
    {
      TSCType type = tsc_ser;
      if (boost::iequals(name,"ser"))
        type = tsc_ser;
      else if (boost::iequals(name,"tte"))
        type = tsc_tte;
      else
      {
        std::ostringstream msg;
        msg << "Unknown conversion from STL_STRING to TSCType enum for "
            << name << "." << std::endl;
        throwError("String2TSCType",msg.str());
      }
      return type;
    };

    //! Different types of scaling operator V
    enum ScaleOpType
    {
      scale_op_identity,       //!< Use the identity matrix (see [2])
      scale_op_cfl_diagonal    //!< Use a diagonal matrix based on the local Courant-Friedrichs-Lewy (CFL) number (see [1]).
    };

    //! Map scaling operator type stl_string to enum
    inline enum ScaleOpType String2ScaleOpType(const std::string& name)
    {
      ScaleOpType type = scale_op_identity;
      if (boost::iequals(name,"identity"))
        type = scale_op_identity;
      else if (boost::iequals(name,"cfl diagonal"))
        type = scale_op_cfl_diagonal;
      else
      {
        std::ostringstream msg;
        msg << "Unknown conversion from STL_STRING to ScaleOperatorType enum for "
            << name << "." << std::endl;
        throwError("String2ScaleOpType",msg.str());
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
    virtual ~PseudoTransient() {};

    //! reset the non-linear solver
    void reset(const NOX::Abstract::Vector& initialGuess,
            const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
            const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests);
    virtual void reset(const NOX::Abstract::Vector& initialGuess,
        const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests);
    virtual void reset(const NOX::Abstract::Vector& initialGuess);

    //! derived
    virtual NOX::StatusTest::StatusType step();

    //! derived
    virtual NOX::StatusTest::StatusType solve();

    //! returns the inverse pseudo time step
    const double& getInversePseudoTimeStep() const;

    //! returns the scaling operator type
    const enum ScaleOpType& getScalingOperatorType() const;

  protected:
    //! initialize the PTC specific variables and call the init function of the base class
    virtual void init();

    //! print the non-linear solver update
    virtual void printUpdate();

    //! compute the nodal pseudo velocity for the CFL update
    void computePseudoVelocity();

    /*!
      \brief Update the pseudo time step

      There are two different options:
      - The switched evolution relaxation
        In the most simple case

            \delta_{n} = \delta_{0} \frac{\|\underline{F}(\underline{x}_{0})\|}{\|\underline{F}(\underline{x}_{n})\|}

      - The temporal truncation error

            MISSING AT THE MOMENT
     */
    void updatePseudoTimeStep();

  private:
    //! throw class specific error
    void throwError(
        const std::string& functionName,
        const std::string& errorMsg) const;

  protected:
    //! Inner Stopping test
    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> iTestPtr_;

    //! Pointer to a linear system pre/post operator, which is used
    //! to modify the jacobian directly from the non-linear solver.
    Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic> prePostLinSysPtr_;

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
    //! solution time derivative used for scaling operator V in pseudo-transient paper.
    Teuchos::RCP<NOX::Abstract::Vector> xDot_;

//    //! If set to true, the candidate direction will use the transient residual instead of the steady-state residual.  This is a modification of the Kelley-Keyes paper.
//    bool use_transient_residual_;

    //! Maximum number of iterations before pseudo-transient is disabled and the algorithm switches to a line search-based direct to steady state solve.
    int maxPseudoTransientIterations_;

    //! time step control type
    enum TSCType tscType_;

    //! scaling operator type
    enum ScaleOpType scaleOpType_;

    //! vector norm type (necessary for the time step control)
    enum NOX::Abstract::Vector::NormType normType_;
    //@}

  }; // class PseudoTransient

} // namespace Solver
namespace LinSystem {
namespace PrePostOp {
/*!
  \brief Pseudo Transient Continuation (PTC) non-linear solver helper class

  This class implementation is an example for the use of the NOX::NLN::LinSystem::PrePostOp::Generic
  and its wrapper the NOX::NLN::LinSystem::PrePostOperator classes. The PTC solver uses them to modify
  the linear system or to be even more precise the jacobian. The pre/post operator gives you access
  to the linearsystem directly before and after the system is solved. At this moment the jacobian
  evaluation and the right-hand-side evaluation should be done. You can modify the system in such a
  way as your non-linear solver prescribes it. Here we use the specified scaling operator type for
  the decision how to modify the jacobian.

  \author Michael Hiermeier
 */
  class PseudoTransient : public Generic
  {
  public:
    //! constructor
    PseudoTransient(const NOX::NLN::Solver::PseudoTransient& ptcsolver);

    //! destructor
    virtual ~PseudoTransient() { };

    //! derived
    virtual void runPreApplyJacobianInverse(
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

    //! derived (does nothing at the moment)
    virtual void runPostApplyJacobianInverse(
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

  protected:
    //! modify the jacobian as defined by the scaling operator type
    virtual void modifyJacobian(LINALG::SparseMatrix& jac);

  private:
    //! read-only access
    const NOX::NLN::Solver::PseudoTransient& ptcsolver_;
  }; // class PseudoTransient
} // namespace PrePostOp
} // namespace LinSystem
} // namespace NLN
} // namespace NOX

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_SOLVER_PTC_H_ */
