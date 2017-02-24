/*!----------------------------------------------------------------------
\file pasi_partitioned_onewaycoup.cpp

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

#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_str_pasiwrapper.H"

#include "../drt_particle/particle_algorithm.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | constructor                                           sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
PASI::PASI_PartOneWayCoup::PASI_PartOneWayCoup(
    const Epetra_Comm&              comm,         //! communicator
    const Teuchos::ParameterList&   pasi_params   //! input parameters for particle structure interaction
    ) :
    // instantiate pasi algorithm class
    PartitionedAlgo(comm,pasi_params)
{

  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.

} // PASI::PASI_PartOneWayCoup::PASI_PartOneWayCoup()

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

    // update particle wall
    UpdateParticleWall();

    // output of structure field
    StructOutput();

    // solve particle time step
    ParticleStep();

    // output of particle field
    ParticleOutput();
  }

  return;
} // PASI::PASI_PartOneWayCoup::Timeloop()

/*----------------------------------------------------------------------*
 | structural time step                                  sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::StructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "                           STRUCTURE SOLVER                           " << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
  }

  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartOneWayCoup::StructStep");

  // Newton-Raphson iteration
  structure_->Solve();

  return;
} // PASI::PASI_PartOneWayCoup::StructStep()

/*----------------------------------------------------------------------*
 | structural output                                     sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::StructOutput()
{
  // calculate stresses, strains, energies
  structure_->PrepareOutput();

  // update all single field solvers
  structure_->Update();

  // write output to files
  structure_->Output();

  // write output to screen
  structure_->PrintStep();

  return;
} // PASI::PASI_PartOneWayCoup::StructOutput()

/*----------------------------------------------------------------------*
 | particle time step                                    sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::ParticleStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "                           PARTICLE SOLVER                            " << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
  }

  TEUCHOS_FUNC_TIME_MONITOR("PASI::PASI_PartOneWayCoup::ParticleStep");

  particles_->AdapterParticle()->UpdateExtActions();

  particles_->AdapterParticle()->IntegrateStep();

  particles_->NormDemThermoAdapt();

  return;
} // PASI::PASI_PartOneWayCoup::ParticleStep()

/*----------------------------------------------------------------------*
 | particle output                                       sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::ParticleOutput()
{
  particles_->AdapterParticle()->PrepareOutput();

  particles_->UpdateConnectivity();

  particles_->AdapterParticle()->Update();

  particles_->Output();

  return;
} // PASI::PASI_PartOneWayCoup::ParticleOutput()

/*----------------------------------------------------------------------*
 | update particle wall                                  sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
void PASI::PASI_PartOneWayCoup::UpdateParticleWall()
{
  // hand structural states to particle field
  particles_->SetStructStates(structure_->Dispn(),structure_->Dispnp(),structure_->Velnp());

  // set structural displacements to particle wall
  particles_->SetUpWallDiscret();

  return;
} // PASI::PASI_PartOneWayCoup::UpdateParticleWall()
