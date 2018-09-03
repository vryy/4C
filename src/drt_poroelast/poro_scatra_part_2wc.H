/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_2wc.H

 \brief  partitioned two way coupled poroelasticity scalar transport interaction algorithms

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


#ifndef PORO_SCATRA_PART_2WC_H_
#define PORO_SCATRA_PART_2WC_H_

#include "poro_scatra_part.H"
#include <Epetra_Vector.h>

namespace POROELAST
{
  class PoroScatraPart2WC : public PoroScatraPart
  {
   public:
    //! explicit constructor
    explicit PoroScatraPart2WC(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor
    virtual ~PoroScatraPart2WC(){};

    //! full time loop
    virtual void Timeloop();

    //! read and set fields needed for restart
    void ReadRestart(int restart);

    //! prepare time step for single fields
    void PrepareTimeStep(bool printheader = true);

    //! perform iteration loop between fields
    void Solve();

    //! prepare output
    void PrepareOutput();

    //! update time step
    void Update();

    //! write output print to screen
    void Output();

   protected:
    //! perform iteration step of structure field
    void DoPoroStep();

    //! perform iteration step of scatra field
    void DoScatraStep();

    //! convergence check of outer loop
    bool ConvergenceCheck(int itnum);

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;
    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;
  };

}  // namespace POROELAST


#endif /* PORO_SCATRA_PART_2WC_H_ */
