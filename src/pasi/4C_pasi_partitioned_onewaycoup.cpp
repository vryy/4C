/*---------------------------------------------------------------------------*/
/*! \file
\brief one way coupled partitioned algorithm for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_pasi_partitioned_onewaycoup.hpp"

#include "4C_particle_algorithm.hpp"
#include "4C_particle_wall_datastate.hpp"
#include "4C_particle_wall_interface.hpp"

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
        particlealgorithm_->get_particle_wall_handler_interface();

    // get wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
        particlewallinterface->GetWallDataState();

    if (walldatastate->GetDispRow() == Teuchos::null or
        walldatastate->GetDispCol() == Teuchos::null)
      FOUR_C_THROW("wall displacements not initialized!");
    if (walldatastate->GetVelCol() == Teuchos::null)
      FOUR_C_THROW("wall velocities not initialized!");
    if (walldatastate->GetAccCol() == Teuchos::null)
      FOUR_C_THROW("wall accelerations not initialized!");
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
    extract_interface_states();

    // set interface states
    SetInterfaceStates(intfdispnp_, intfvelnp_, intfaccnp_);

    // particle time step
    ParticleStep();

    // post evaluate time step
    post_evaluate_time_step();

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
