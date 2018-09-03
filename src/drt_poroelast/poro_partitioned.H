/*----------------------------------------------------------------------*/
/*!
 \file poro_partitioned.H

 \brief  base class for partitioned poroelasticity algorithms

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *-----------------------------------------------------------------------*/

#ifndef PORO_PARTITIONED_H_
#define PORO_PARTITIONED_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "../drt_inpar/inpar_poroelast.H"
#include "../drt_inpar/inpar_structure.H"
#include "poro_base.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  /// base class of all monolithic Poroelasticity algorithms
  class Partitioned : public PoroBase
  {
   public:
    //! create using a Epetra_Comm
    explicit Partitioned(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! destructor
    virtual ~Partitioned(){};

    //! proceed one time step (prepare, solve, update)
    void DoTimeStep();

    //! initialise system
    virtual void SetupSystem();

    //! dof row map of Structure field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapStructure();

    //! dof row map of Fluid field
    virtual Teuchos::RCP<const Epetra_Map> DofRowMapFluid();

   protected:
    //! prepare new time step
    void PrepareTimeStep();

    //! solve one time step of structural problem
    void DoStructStep();

    //! solve one time step of fluid problem
    void DoFluidStep();

    //! solve one time step (iteration between fields)
    void Solve();

    //! update and write output to screen and files after solved time step
    void UpdateAndOutput();

    //! convergence check of outer loop
    bool ConvergenceCheck(int itnum);

    //! aitken relaxation
    void AitkenRelax();

    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;

    Teuchos::RCP<Epetra_Vector> fluidveln_;  //!< global fluid velocities and pressures

    //! @name Aitken relaxation

    //! difference of last two solutions
    // del = r^{i+1}_{n+1} = d^{i+1}_{n+1} - d^i_{n+1}
    Teuchos::RCP<Epetra_Vector> del_;
    //! difference of difference of last two pair of solutions
    // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
    Teuchos::RCP<Epetra_Vector> delhist_;
    //! relaxation factor of old iteration step i
    double omegan_;
    //! relaxation factor of current iteration step i+1
    double omeganp_;

    //@}

  };  // Partitioned

}  // namespace POROELAST

#endif /* PORO_PARTITIONED_H_ */
