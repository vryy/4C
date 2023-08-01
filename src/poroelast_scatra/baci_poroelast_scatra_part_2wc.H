/*----------------------------------------------------------------------*/
/*! \file

 \brief  partitioned two way coupled poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/


#ifndef BACI_POROELAST_SCATRA_PART_2WC_H
#define BACI_POROELAST_SCATRA_PART_2WC_H

#include "baci_poroelast_scatra_part.H"

#include <Epetra_Vector.h>

namespace POROELASTSCATRA
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

}  // namespace POROELASTSCATRA


#endif  // POROELAST_SCATRA_PART_2WC_H
