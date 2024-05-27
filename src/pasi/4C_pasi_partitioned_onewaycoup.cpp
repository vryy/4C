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
  check_is_init();
  check_is_setup();

  // time loop
  while (NotFinished())
  {
    // prepare time step
    prepare_time_step();

    // pre evaluate time step
    pre_evaluate_time_step();

    // structural time step
    struct_step();

    // extract interface states
    extract_interface_states();

    // set interface states
    set_interface_states(intfdispnp_, intfvelnp_, intfaccnp_);

    // particle time step
    particle_step();

    // post evaluate time step
    post_evaluate_time_step();

    // output of fields
    output();
  }
}

void PASI::PasiPartOneWayCoup::output()
{
  // output of structure field
  struct_output();

  // output of particle field
  particle_output();
}

FOUR_C_NAMESPACE_CLOSE
