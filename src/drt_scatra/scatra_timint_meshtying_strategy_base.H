/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_base.H

\brief Abstract interface for meshtying strategies in scalar transport problems (including standard
strategy without meshtying)

\level 2

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251


*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_BASE_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_BASE_H

#include "scatra_timint_convcheck_strategies.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

// forward declarations
namespace IO
{
  class InputControl;
}

namespace LINALG
{
  class KrylovProjector;
  class Solver;
  class SparseOperator;
}  // namespace LINALG

namespace SCATRA
{
  class ConvCheckStrategyBase;
  class ScaTraTimIntImpl;

  /*!
  \brief Abstract interface for meshtying strategies in scalar transport problems (including
  standard strategy without meshtying)

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. Every specific meshtying strategy, e.g. fluid-fluid meshtying, scatra-
  scatra interface coupling, and also the standard algorithm without meshtying, has to be
  implemented in a subclass derived from this abstract, purely virtual interface class.

  \author fang (fang@lnm.mw.tum.de)
  */

  class MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBase(SCATRA::ScaTraTimIntImpl* scatratimint)
        : convcheckstrategy_(Teuchos::null), scatratimint_(scatratimint)
    {
      return;
    }


    //! destructor
    virtual ~MeshtyingStrategyBase() { return; };

    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
        ) const
    {
      if (convcheckstrategy_ == Teuchos::null)
        dserror("Strategy for Newton-Raphson convergence check has not been instantiated!");

      return convcheckstrategy_->AbortNonlinIter(scatratimint, actresidual);
    }

    //! perform convergence check for outer iteration in partitioned coupling schemes
    bool AbortOuterIter(const ScaTraTimIntImpl& scatratimint  //!< scalar transport time integrator
        ) const
    {
      if (convcheckstrategy_ == Teuchos::null)
        dserror("Strategy for outer convergence check has not been instantiated!");

      return convcheckstrategy_->AbortOuterIter(scatratimint);
    }

    //! provide global state vectors for element evaluation
    virtual void AddTimeIntegrationSpecificVectors() const { return; };

    //! compute time step size
    virtual void ComputeTimeStepSize(double& dt) { return; };

    //! compute time derivatives of discrete state variables
    virtual void ComputeTimeDerivative() const { return; };

    //! condense global system of equations
    virtual void CondenseMatAndRHS(
        const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& residual,               //!< residual vector
        const bool calcinittimederiv = false  //!< flag for calculation of initial time derivative
        ) const
    {
      return;
    };

    //! return global map of degrees of freedom
    virtual const Epetra_Map& DofRowMap() const = 0;

    /*!
    \brief Evaluate a given condition

     Evaluate terms of your weak formulation on elements marked with a given condition.

    \param params (in):         List of parameters for use at element level
    \param systemmatrix1 (out): Sparse matrix that may be changed by
                                assembly of boundary element contributions.
                                May not be Teuchos::null.
                                Matrix must be systemmatrix->Filled()==false on input.
    \param systemmatrix2 (out): Sparse matrix that may be changed by
                                assembly of boundary element contributions.
                                May not be Teuchos::null.
                                Matrix must be systemmatrix->Filled()==false on input.
    \param systemvector1 (out): Vector to assemble BCs to.
                                The vector is NOT initialized to zero by this method.
    \param systemvector2 (out): Vector to assemble BCs to.
                                The vector is NOT initialized to zero by this method.
    \param systemvector3 (out): Vector to assemble BCs to.
                                The vector is NOT initialized to zero by this method.
    \param condstring (in):     Name of condition to be evaluated
    \param condid (in):         Condition ID

    \return void
    \date 08/16
    \author rauch
    */
    virtual void EvaluateCondition(Teuchos::ParameterList& params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring, const int condid)
    {
      dserror("EvaluateCondition(...) is not implemented in MeshtyingStrategyBase.");
    };

    //! compute meshtying residual terms and their linearizations
    virtual void EvaluateMeshtying() = 0;

    //! explicit predictor step to obtain better starting value for Newton-Raphson iteration
    virtual void ExplicitPredictor() const { return; };

    //! include Dirichlet conditions into condensation
    virtual void IncludeDirichletInCondensation() const { return; };

    //! init meshtying objects
    virtual void InitMeshtying() = 0;

    //! initialize system matrix
    virtual Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const = 0;

    //! output solution for post-processing
    virtual void Output() const { return; };

    //! output restart information
    virtual void OutputRestart() const { return; };

    //! read restart data
    virtual void ReadRestart(const int step,                  //!< restart step
        Teuchos::RCP<IO::InputControl> input = Teuchos::null  //!< control file manager
        ) const
    {
      return;
    };

    //! set general parameters for element evaluation
    virtual void SetElementGeneralParameters(Teuchos::ParameterList& parameters) const { return; };

    //! compute history vector, i.e., the history part of the right-hand side vector with all
    //! contributions from the previous time step
    virtual void SetOldPartOfRHS() const { return; };

    /*!
    \brief Set state on an discretization hidden in any MeshtyingStrategy.

     This method should encapsulate a standard SetState(...) of the discretizations
     you want to consider. See \ref HeterogeneousReactionStrategy for an example.

    \param nds (in): number of dofset
    \param name (in): Name of data
    \param state (in): vector of some data

    \return void
    \date 12/16
    \author rauch
    */
    virtual void SetState(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
    {
      dserror(
          "SetState(...) is not implemented in MeshtyingStrategyBase.\n"
          "SetState(...) allows to set states on discretizations held within any\n"
          "specific MeshtyingStrategy. See e.g. HeterogeneousReactionStrategy for\n"
          "an example implementation.");
    };

    //! setup meshtying objects
    virtual void SetupMeshtying() = 0;

    //! solve resulting linear system of equations
    virtual void Solve(const Teuchos::RCP<LINALG::Solver>& solver,  //!< solver
        const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,   //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& increment,               //!< increment vector
        const Teuchos::RCP<Epetra_Vector>& residual,                //!< residual vector
        const Teuchos::RCP<Epetra_Vector>& phinp,                   //!< state vector at time n+1
        const int& iteration,  //!< number of current Newton-Raphson iteration
        const Teuchos::RCP<LINALG::KrylovProjector>& projector  //!< Krylov projector
        ) const = 0;

    //! return linear solver for global system of linear equations
    virtual const LINALG::Solver& Solver() const = 0;

    //! update solution after convergence of the nonlinear Newton-Raphson iteration
    virtual void Update() const { return; };

   protected:
    //! instantiate strategy for Newton-Raphson convergence check
    virtual void InitConvCheckStrategy() = 0;

    //! strategy for Newton-Raphson convergence check called by scalar transport time integrator
    Teuchos::RCP<SCATRA::ConvCheckStrategyBase> convcheckstrategy_;

    //! scalar transport time integrator
    SCATRA::ScaTraTimIntImpl* scatratimint_;

   private:
    //! copy constructor
    MeshtyingStrategyBase(const MeshtyingStrategyBase& old);
  };  // class MeshtyingStrategyBase
}  // namespace SCATRA
#endif  // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_BASE_H
