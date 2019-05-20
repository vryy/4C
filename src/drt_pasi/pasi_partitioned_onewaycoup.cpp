/*---------------------------------------------------------------------------*/
/*!
\brief one way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned_onewaycoup.H"

#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
PASI::PASI_PartOneWayCoup::PASI_PartOneWayCoup(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PartitionedAlgo(comm, params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | partitioned one way coupled timeloop                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  while (NotFinished())
  {
    // redistribute load in parallel
    particles_->DynamicLoadBalancing();

    // counter and print header
    PrepareTimeStep();

    // structural time step
    StructStep();

    // set structural states
    SetStructStates();

    // particle time step
    ParticleStep();

    // output of structure field
    StructOutput();

    // output of particle field
    ParticleOutput();
  }
}
