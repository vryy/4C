/*---------------------------------------------------------------------------*/
/*! \file
\brief one way coupled partitioned algorithm for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_pasi_partitioned_onewaycoup.hpp"

#include "baci_particle_algorithm.hpp"
#include "baci_particle_wall_datastate.hpp"
#include "baci_particle_wall_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PASI::PasiPartOneWayCoup::PasiPartOneWayCoup(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : PartitionedAlgo(comm, params)
{
  // empty constructor
}

void PASI::PasiPartOneWayCoup::Setup()
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

void PASI::PasiPartOneWayCoup::Timeloop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // time loop
  while (NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // pre evaluate time step
    PreEvaluateTimeStep();

    // structural time step
    StructStep();

    // extract interface states
    ExtractInterfaceStates();

    // set interface states
    SetInterfaceStates(intfdispnp_, intfvelnp_, intfaccnp_);

    // particle time step
    ParticleStep();

    // post evaluate time step
    PostEvaluateTimeStep();

    // output of fields
    Output();
  }
}

void PASI::PasiPartOneWayCoup::Output()
{
  // output of structure field
  StructOutput();

  // output of particle field
  ParticleOutput();
}

FOUR_C_NAMESPACE_CLOSE
