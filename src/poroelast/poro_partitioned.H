/*----------------------------------------------------------------------*/
/*! \file

 \brief  base class for partitioned poroelasticity algorithms

\level 2

 *-----------------------------------------------------------------------*/

#ifndef PORO_PARTITIONED_H_
#define PORO_PARTITIONED_H_

#include "inpar_poroelast.H"
#include "inpar_structure.H"
#include "poro_base.H"

namespace POROELAST
{
  //! base class of all monolithic Poroelasticity algorithms
  class Partitioned : public PoroBase
  {
   public:
    //! create using a Epetra_Comm
    explicit Partitioned(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! proceed one time step (prepare, solve, update)
    void DoTimeStep() override;

    //! initialise system
    void SetupSystem() override;

    //! dof row map of Structure field
    Teuchos::RCP<const Epetra_Map> DofRowMapStructure() override;

    //! dof row map of Fluid field
    Teuchos::RCP<const Epetra_Map> DofRowMapFluid() override;

   protected:
    //! prepare new time step
    void PrepareTimeStep() override;

    //! solve one time step of structural problem
    void DoStructStep();

    //! solve one time step of fluid problem
    void DoFluidStep();

    //! solve one time step (iteration between fields)
    void Solve() override;

    //! update and write output to screen and files after solved time step
    void UpdateAndOutput();

    //! convergence check of outer loop
    bool ConvergenceCheck(int itnum);

    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;

    Teuchos::RCP<Epetra_Vector> fluidveln_;  //!< global fluid velocities and pressures
  };

}  // namespace POROELAST

#endif /* PORO_PARTITIONED_H_ */
