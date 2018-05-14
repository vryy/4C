/*!----------------------------------------------------------------------
\file porofluidmultiphase_meshtying_strategy_std.H

\brief standard case without mesh tying

\level 3

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/


*----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_STD_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_STD_H_

#include "porofluidmultiphase_meshtying_strategy_base.H"

namespace POROFLUIDMULTIPHASE
{
class MeshtyingStrategyStd : public MeshtyingStrategyBase
{
  public:

    //! constructor
    explicit MeshtyingStrategyStd(
        POROFLUIDMULTIPHASE::TimIntImpl* porofluidmultitimint,
        const Teuchos::ParameterList&           probparams,
        const Teuchos::ParameterList&           poroparams
        );

    //! destructor
    virtual ~MeshtyingStrategyStd();

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! prepare time step
    virtual void PrepareTimeStep();

    //! update
    virtual void Update();

    //! output
    virtual void Output();

    //! solve linear system of equations
    virtual void LinearSolve(
        Teuchos::RCP<LINALG::Solver>         solver,
        Teuchos::RCP<LINALG::SparseOperator> sysmat,
        Teuchos::RCP<Epetra_Vector>          increment,
        Teuchos::RCP<Epetra_Vector>          residual
        );

    //! calculate norms for convergence checks
    virtual void CalculateNorms(
        double& preresnorm,
        double& incprenorm,
        double& prenorm,
        const Teuchos::RCP<const Epetra_Vector> increment
        );

    //! create the field test
    virtual void CreateFieldTest();

    //! restart
    virtual void ReadRestart(const int step);

    //! evaluate mesh tying
    virtual void Evaluate();

    //! extract increments and update mesh tying
    virtual Teuchos::RCP<const Epetra_Vector> ExtractAndUpdateIter(const Teuchos::RCP<const Epetra_Vector> inc);

    //! access to global (combined) increment of coupled problem
    virtual Teuchos::RCP<const Epetra_Vector> CombinedIncrement(const Teuchos::RCP<const Epetra_Vector> inc) const;

};

}



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_MESHTYING_STRATEGY_STD_H_ */
