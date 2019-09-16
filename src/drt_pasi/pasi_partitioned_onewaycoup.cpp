/*---------------------------------------------------------------------------*/
/*! \file
\brief one way coupled partitioned algorithm for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include "pasi_partitioned_onewaycoup.H"

#include "../drt_particle_algorithm/particle_algorithm.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

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
 | setup pasi algorithm                                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::Setup()
{
  // call base class setup
  PASI::PartitionedAlgo::Setup();

  // safety check
  {
    // get interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface =
        particlealgorithm_->GetParticleWallHandlerInterface();

    // get wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
        particlewallinterface->GetWallDataState();

    if (walldatastate->GetDispRow() == Teuchos::null or
        walldatastate->GetDispCol() == Teuchos::null)
      dserror("wall displacements not initialized!");
    if (walldatastate->GetVelCol() == Teuchos::null) dserror("wall velocities not initialized!");
    if (walldatastate->GetAccCol() == Teuchos::null) dserror("wall accelerations not initialized!");
  }
}

/*---------------------------------------------------------------------------*
 | partitioned one way coupled timeloop                       sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // time loop
  while (NotFinished())
  {
    // counter and print header
    PrepareTimeStep();

    // structural time step
    StructStep();

    // set structural states
    SetStructStates();

    // particle time step
    ParticleStep();

    // output of fields
    Output();
  }
}

/*---------------------------------------------------------------------------*
 | output of fields                                           sfuchs 09/2019 |
 *---------------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::Output()
{
  // output of structure field
  StructOutput();

  // output of particle field
  ParticleOutput();
}
