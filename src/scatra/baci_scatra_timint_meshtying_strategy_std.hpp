/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for standard scalar transport problems (without meshtying)

\level 2



*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_STD_HPP
#define FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_STD_HPP

#include "baci_config.hpp"

#include "baci_scatra_timint_meshtying_strategy_base.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  struct SolverParams;
}

namespace SCATRA
{
  /*!
  \brief Standard solution strategy for standard scalar transport problems (without meshtying)

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the standard solution strategy for standard scalar
  transport problems without meshtying.

  */

  class MeshtyingStrategyStd : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyStd(SCATRA::ScaTraTimIntImpl* scatratimint);

    //! return global map of degrees of freedom
    const Epetra_Map& DofRowMap() const override;

    /*!
    \brief Evaluate a given condition

     Evaluate terms of your weak formulation on elements marked with a given condition.

    \return void
    \date 08/16
    \author rauch
    */
    void EvaluateCondition(Teuchos::ParameterList& params,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring,
        const int condid) override
    {
      dserror("EvaluateCondition(...) is not implemented in MeshtyingStrategyStd.");
    };

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying() override;

    //! init meshtying objects
    void InitMeshtying() override;

    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> InterfaceMaps() const override
    {
      dserror("InterfaceMaps() is not implemented in MeshtyingStrategyStd.");
      return Teuchos::null;
    }

    bool SystemMatrixInitializationNeeded() const override { return false; }

    Teuchos::RCP<CORE::LINALG::SparseOperator> InitSystemMatrix() const override
    {
      dserror(
          "This meshtying strategy does not need to initialize the system matrix, but relies "
          "instead on the initialization of the field. If this changes, you also need to change "
          "'SystemMatrixInitializationNeeded()' to return true");
      // dummy return
      return Teuchos::null;
    }

    //! setup meshtying objects
    void SetupMeshtying() override;

    //! solve resulting linear system of equations
    void Solve(const Teuchos::RCP<CORE::LINALG::Solver>& solver,         //!< solver
        const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& increment,                    //!< increment vector
        const Teuchos::RCP<Epetra_Vector>& residual,                     //!< residual vector
        const Teuchos::RCP<Epetra_Vector>& phinp,  //!< state vector at time n+1
        const int iteration,                       //!< number of current Newton-Raphson iteration
        CORE::LINALG::SolverParams& solver_params) const override;

    //! return linear solver for global system of linear equations
    const CORE::LINALG::Solver& Solver() const override;

   protected:
    //! instantiate strategy for Newton-Raphson convergence check
    void InitConvCheckStrategy() override;

   private:
    //! copy constructor
    MeshtyingStrategyStd(const MeshtyingStrategyStd& old);
  };  // class MeshtyingStrategyStd
}  // namespace SCATRA
BACI_NAMESPACE_CLOSE

#endif
