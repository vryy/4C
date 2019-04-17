/*----------------------------------------------------------------------*/
/*!
 \file scatra_timint_meshtying_strategy_artery.H

 \brief routines for coupling between 1D arterial network and 2D/3D
        scatra-algorithm

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_SCATRA_TIMINT_MESHTYING_STRATEGY_ARTERY_H_
#define SRC_DRT_SCATRA_SCATRA_TIMINT_MESHTYING_STRATEGY_ARTERY_H_

#include "scatra_timint_meshtying_strategy_base.H"
#include "../linalg/linalg_utils.H"

// forward declaration
namespace DRT
{
  class Discretization;
}
namespace ADAPTER
{
  class Coupling;
  class ArtNet;
}  // namespace ADAPTER
namespace FSI
{
  class Monolithic;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScaTraArtCouplBase;
}

namespace SCATRA
{
  class MeshtyingStrategyArtery : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyArtery(
        SCATRA::ScaTraTimIntImpl* scatratimint  //!< scalar transport time integrator
    );

    //! destructor
    virtual ~MeshtyingStrategyArtery() { return; };

    //! return global map of degrees of freedom
    const Epetra_Map& DofRowMap() const;

    //! return global map of degrees of freedom
    Teuchos::RCP<const Epetra_Map> ArtScatraDofRowMap() const;

    //! evaluate
    void EvaluateMeshtying();

    //! init
    void InitMeshtying();

    //! init system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

    //! setup
    void SetupMeshtying();

    //! solver
    const LINALG::Solver& Solver() const;

    //! init the convergence check
    void InitConvCheckStrategy();

    //! solve resulting linear system of equations
    void Solve(const Teuchos::RCP<LINALG::Solver>& solver,         //!< solver
        const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& increment,              //!< increment vector
        const Teuchos::RCP<Epetra_Vector>& residual,               //!< residual vector
        const Teuchos::RCP<Epetra_Vector>& phinp,                  //!< state vector at time n+1
        const int& iteration,  //!< number of current Newton-Raphson iteration
        const Teuchos::RCP<LINALG::KrylovProjector>& projector  //!< Krylov projector
        ) const;

    void SetupSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
        const Teuchos::RCP<Epetra_Vector>& residual                             //!< residual vector
        ) const;

    //! init the convergence check
    void SetArteryScatraTimeIntegrator(Teuchos::RCP<SCATRA::ScaTraTimIntImpl> artscatratimint);

    //! set the artery time integrator
    void SetArteryTimeIntegrator(Teuchos::RCP<ADAPTER::ArtNet> arttimint);

    //! set the element pairs that are close as found by search algorithm
    void SetNearbyElePairs(const std::map<int, std::set<int>>* nearbyelepairs);

    //! prepare a time step
    void PrepareTimeStep() const;

    //! set the artery pressure
    void SetArteryPressure() const;

    //! apply mesh movement
    void ApplyMeshMovement();

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> CombinedSystemMatrix()
    {
      return comb_systemmatrix_;
    }

    //! get the combined rhs
    Teuchos::RCP<Epetra_Vector> CoupledRHS() const { return rhs_; }

    //! get the combined increment
    Teuchos::RCP<Epetra_Vector> CombinedIncrement() const { return comb_increment_; }

    //! access to time integrator
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ArtScatraField() { return artscatratimint_; }

    //! check if initial fields match
    void CheckInitialFields() const;

    void UpdateArtScatraIter(Teuchos::RCP<const Epetra_Vector> combined_inc);

   private:
    //! time integrators
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> artscatratimint_;
    Teuchos::RCP<ADAPTER::ArtNet> arttimint_;

    //! mesh tying object
    Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase> arttoscatracoupling_;

    //! the two discretizations
    Teuchos::RCP<DRT::Discretization> artscatradis_;
    Teuchos::RCP<DRT::Discretization> scatradis_;

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> comb_systemmatrix_;

    //! combined rhs
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! combined rhs
    Teuchos::RCP<Epetra_Vector> comb_increment_;

  };  // class MeshtyingStrategyArtery

}  // namespace SCATRA

#endif /* SRC_DRT_SCATRA_SCATRA_TIMINT_MESHTYING_STRATEGY_ARTERY_H_ */
