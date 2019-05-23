/*!----------------------------------------------------------------------

\brief one way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#include "pasi_partitioned_onewaycoup.H"

#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

/*----------------------------------------------------------------------*
 | constructor                                           sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
PASI::PASI_PartOneWayCoup::PASI_PartOneWayCoup(const Epetra_Comm& comm,  //! communicator
    const Teuchos::ParameterList&
        pasi_params  //! input parameters for particle structure interaction
    )
    :  // instantiate pasi algorithm class
      PartitionedAlgo(comm, pasi_params)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.

}  // PASI::PASI_PartOneWayCoup::PASI_PartOneWayCoup()

/*----------------------------------------------------------------------*
 | partitioned one way coupled timeloop                  sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  while (NotFinished())
  {
    // redistribute load in parallel
    particles_->DynamicLoadBalancing();

    // increment time and step
    PrepareTimeStep(true);

    // solve structural time step
    StructStep();

    // set structural states in particle field
    SetStructStates();

    // solve particle time step
    ParticleStep();

    // output of structure field
    StructOutput();

    // output of particle field
    ParticleOutput();
  }

  return;
}  // PASI::PASI_PartOneWayCoup::Timeloop()
