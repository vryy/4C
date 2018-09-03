/*!----------------------------------------------------------------------
\file porofluidmultiphase_meshtying_strategy_artery.H

\brief routines for coupling with artery network

\level 3

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/


*----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_ARTERY_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_ARTERY_H_

#include "porofluidmultiphase_meshtying_strategy_base.H"

// forward declaration
namespace POROMULTIPHASESCATRA
{
  class PoroMultiPhaseScaTraArtCouplBase;
}

namespace POROFLUIDMULTIPHASE
{
  class MeshtyingStrategyArtery : public MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyArtery(POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams);

    //! destructor
    virtual ~MeshtyingStrategyArtery();

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! prepare time step
    virtual void PrepareTimeStep();

    //! update
    virtual void Update();

    //! output
    virtual void Output();

    //! solve linear system of equations
    virtual void LinearSolve(Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<LINALG::SparseOperator> sysmat, Teuchos::RCP<Epetra_Vector> increment,
        Teuchos::RCP<Epetra_Vector> residual);

    //! calculate norms for convergence checks
    virtual void CalculateNorms(double& preresnorm, double& incprenorm, double& prenorm,
        const Teuchos::RCP<const Epetra_Vector> increment);

    //! create the field test
    virtual void CreateFieldTest();

    //! restart
    virtual void ReadRestart(const int step);

    //! evaluate mesh tying
    virtual void Evaluate();

    //! extract increments and update mesh tying
    virtual Teuchos::RCP<const Epetra_Vector> ExtractAndUpdateIter(
        const Teuchos::RCP<const Epetra_Vector> inc);

    // return arterial network time integrator
    virtual Teuchos::RCP<ADAPTER::ArtNet> ArtNetTimInt() { return artnettimint_; }

    //! access dof row map
    virtual Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const;

    //! right-hand side alias the dynamic force residual for coupled system
    virtual Teuchos::RCP<const Epetra_Vector> ArteryPorofluidRHS() const;

    //! access to block system matrix of artery poro problem
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ArteryPorofluidSysmat() const;

    //! get global (combined) increment of coupled problem
    virtual Teuchos::RCP<const Epetra_Vector> CombinedIncrement(
        const Teuchos::RCP<const Epetra_Vector> inc) const;

    // check if initial fields on coupled DOFs are equal
    void CheckInitialFields(Teuchos::RCP<const Epetra_Vector> vec_cont) const;

    //! apply the mesh movement
    void ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> disp) const;

   protected:
    //! artery time integration
    Teuchos::RCP<ADAPTER::ArtNet> artnettimint_;

    //! artery discretization
    Teuchos::RCP<DRT::Discretization> arterydis_;

    //! the mesh tying object
    Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplBase> arttoporofluidcoupling_;

    //! block systemmatrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> comb_systemmatrix_;

    //! global rhs
    Teuchos::RCP<Epetra_Vector> rhs_;

    //! global increment
    Teuchos::RCP<Epetra_Vector> comb_increment_;

    //! global solution at time n+1
    Teuchos::RCP<Epetra_Vector> comb_phinp_;
  };

}  // namespace POROFLUIDMULTIPHASE



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_ARTERY_H_ */
