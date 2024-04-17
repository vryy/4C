/*----------------------------------------------------------------------*/
/*! \file

 \brief  partitioned two way coupled poroelasticity scalar transport interaction algorithms

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_POROELAST_SCATRA_PART_2WC_HPP
#define FOUR_C_POROELAST_SCATRA_PART_2WC_HPP

#include "baci_config.hpp"

#include "baci_poroelast_scatra_part.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

namespace POROELASTSCATRA
{
  class PoroScatraPart2WC : public PoroScatraPart
  {
   public:
    //! explicit constructor
    explicit PoroScatraPart2WC(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! full time loop
    void Timeloop() override;

    //! read and set fields needed for restart
    void ReadRestart(int restart) override;

    //! prepare time step for single fields
    void PrepareTimeStep(bool printheader = true) override;

    //! perform iteration loop between fields
    void Solve() override;

    //! prepare output
    void PrepareOutput() override;

    //! update time step
    void Update() override;

    //! write output print to screen
    void Output() override;

   protected:
    //! perform iteration step of structure field
    void DoPoroStep() override;

    //! perform iteration step of scatra field
    void DoScatraStep() override;

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


FOUR_C_NAMESPACE_CLOSE

#endif
