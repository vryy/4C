/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid-fluid meshtying strategy for standard scalar transport problems

\level 2


*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_H

#include "../drt_inpar/inpar_fluid.H"

#include "scatra_timint_meshtying_strategy_base.H"

// forward declaration
namespace FLD
{
  class Meshtying;
}

namespace SCATRA
{
  /*!
  \brief Fluid-fluid meshtying strategy for standard scalar transport problems

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the fluid-fluid meshtying strategy for standard
  scalar transport problems.

  */

  class MeshtyingStrategyFluid : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyFluid(SCATRA::ScaTraTimIntImpl* scatratimint);

    //! destructor
    virtual ~MeshtyingStrategyFluid() { return; };

    //! return global map of degrees of freedom
    const Epetra_Map& DofRowMap() const;

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying();

    //! include Dirichlet conditions into condensation
    void IncludeDirichletInCondensation() const;

    //! initialize meshtying objects
    virtual void InitMeshtying();

    bool SystemMatrixInitializationNeeded() const override { return true; }

    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const override;

    const Teuchos::RCP<LINALG::MultiMapExtractor> InterfaceMaps() const override
    {
      dserror("InterfaceMaps() is not implemented in MeshtyingStrategyFluid.");
      return Teuchos::null;
    }

    //! setup meshtying objects
    void SetupMeshtying() override;

    //! solve resulting linear system of equations
    void Solve(const Teuchos::RCP<LINALG::Solver>& solver,         //!< solver
        const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& increment,              //!< increment vector
        const Teuchos::RCP<Epetra_Vector>& residual,               //!< residual vector
        const Teuchos::RCP<Epetra_Vector>& phinp,                  //!< state vector at time n+1
        const int& iteration,  //!< number of current Newton-Raphson iteration
        const Teuchos::RCP<LINALG::KrylovProjector>& projector  //!< Krylov projector
        ) const;

    //! return linear solver for global system of linear equations
    const LINALG::Solver& Solver() const;

   protected:
    //! instantiate strategy for Newton-Raphson convergence check
    virtual void InitConvCheckStrategy();

    //! fluid-fluid meshtying algorithm for internal interface
    Teuchos::RCP<FLD::Meshtying> meshtying_;

    //! type of fluid-fluid meshtying
    enum INPAR::FLUID::MeshTying type_;

   private:
    //! copy constructor
    MeshtyingStrategyFluid(const MeshtyingStrategyFluid& old);
  };  // class MeshtyingStrategyFluid
}  // namespace SCATRA
#endif  // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_H
