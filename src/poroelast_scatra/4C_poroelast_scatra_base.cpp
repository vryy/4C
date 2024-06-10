/*----------------------------------------------------------------------*/
/*! \file

 \brief  base class for all poroelasticity scalar transport interaction algorithms

\level 2


 *----------------------------------------------------------------------*/

#include "4C_poroelast_scatra_base.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_coupling_volmortar_utils.hpp"
#include "4C_fem_dofset_gidbased_wrapper.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_poroelast_base.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_scatra_timint_implicit.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
PoroElastScaTra::PoroScatraBase::PoroScatraBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : AlgorithmBase(comm, timeparams),
      matchinggrid_(Core::UTILS::IntegralValue<bool>(
          Global::Problem::Instance()->poro_scatra_control_params(), "MATCHINGGRID")),
      volcoupl_structurescatra_(Teuchos::null),
      volcoupl_fluidscatra_(Teuchos::null)
{
  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();

  // do some checks
  {
    Inpar::ScaTra::TimeIntegrationScheme timealgo =
        Core::UTILS::IntegralValue<Inpar::ScaTra::TimeIntegrationScheme>(scatradyn, "TIMEINTEGR");
    if (timealgo != Inpar::ScaTra::timeint_one_step_theta and
        timealgo != Inpar::ScaTra::timeint_stationary)
      FOUR_C_THROW(
          "scalar transport in porous media is limited in functionality (only one-step-theta "
          "scheme or stationary case possible)");

    //    Inpar::ScaTra::ConvForm convform
    //    = Core::UTILS::IntegralValue<Inpar::ScaTra::ConvForm>(scatradyn,"CONVFORM");
    //    if ( convform != Inpar::ScaTra::convform_convective )
    //      FOUR_C_THROW("The balance of mass is included in the formulation for scalar transport in
    //      porous media. "
    //          "Set 'CONVFORM' to 'convective' in the SCALAR TRANSPORT DYNAMIC section! ");

    Inpar::ScaTra::VelocityField velfield =
        Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(scatradyn, "VELOCITYFIELD");
    if (velfield != Inpar::ScaTra::velocity_Navier_Stokes)
      FOUR_C_THROW(
          "scalar transport is coupled with the porous medium. Set 'VELOCITYFIELD' to "
          "'Navier_Stokes' in the SCALAR TRANSPORT DYNAMIC section! ");

    //    bool skipinitder
    //    = Core::UTILS::IntegralValue<int>(scatradyn,"SKIPINITDER");
    //    if ( not skipinitder )
    //      FOUR_C_THROW("Calculation of initial time derivative not yet supported for scalar
    //      transport in porous media. Set 'SKIPINITDER' to 'yes' in the SCALAR TRANSPORT DYNAMIC
    //      section! ");
  }

  // the problem is two way coupled, thus each discretization must know the other discretization
  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->GetDis("porofluid");
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->GetDis("scatra");
  setup_coupling(structdis, fluiddis, scatradis);
  // Create the two uncoupled subproblems.
  // 1. poro problem
  poro_ = PoroElast::UTILS::CreatePoroAlgorithm(
      timeparams, comm, false, PoroElastScaTra::UTILS::BuildPoroScatraSplitter(structdis));

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
  // 2. scatra problem
  scatra_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
      timeparams, scatradyn, problem->SolverParams(linsolvernumber), "scatra", true));

  // now we can call Init() on the base algo.
  // time integrator is constructed and initialized inside.
  scatra_->Init();
  scatra_->ScaTraField()->set_number_of_dof_set_displacement(2);
  scatra_->ScaTraField()->set_number_of_dof_set_velocity(2);

  // only now we must call Setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls Setup() on the time integrator inside.
  scatra_->Setup();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::SetupSystem() { poro_->SetupSystem(); }

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::TestResults(const Epetra_Comm& comm)
{
  Global::Problem* problem = Global::Problem::Instance();

  problem->AddFieldTest(poro_->structure_field()->CreateFieldTest());
  problem->AddFieldTest(poro_->fluid_field()->CreateFieldTest());
  problem->AddFieldTest(scatra_->create_sca_tra_field_test());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::SetPoroSolution()
{
  set_mesh_disp();
  set_velocity_fields();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::SetScatraSolution()
{
  Teuchos::RCP<const Epetra_Vector> phinp_s = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> phin_s = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> phinp_f = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> phin_f = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> phidtnp = Teuchos::null;

  if (matchinggrid_)
  {
    phinp_s = scatra_->ScaTraField()->Phinp();
    phinp_f = phinp_s;
    phin_s = scatra_->ScaTraField()->Phin();
    phin_f = phin_s;
    phidtnp = scatra_->ScaTraField()->Phidtnp();
  }
  else
  {
    phinp_s = volcoupl_structurescatra_->apply_vector_mapping12(scatra_->ScaTraField()->Phinp());
    phinp_f = volcoupl_fluidscatra_->apply_vector_mapping12(scatra_->ScaTraField()->Phinp());
    phin_s = volcoupl_structurescatra_->apply_vector_mapping12(scatra_->ScaTraField()->Phin());
    phin_f = volcoupl_fluidscatra_->apply_vector_mapping12(scatra_->ScaTraField()->Phin());
    phidtnp = volcoupl_fluidscatra_->apply_vector_mapping12(scatra_->ScaTraField()->Phidtnp());
  }

  // porous structure
  poro_->structure_field()->discretization()->set_state(2, "scalar", phinp_s);
  poro_->structure_field()->discretization()->set_state(2, "scalarn", phin_s);

  // porous fluid
  poro_->fluid_field()->SetIterScalarFields(phinp_f, phin_f, phidtnp,
      // scatra_->ScaTraField()->discretization()
      poro_->fluid_field()->discretization(), 2);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::set_velocity_fields()
{
  Teuchos::RCP<const Epetra_Vector> convel = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> velnp = Teuchos::null;

  if (matchinggrid_)
  {
    convel = poro_->fluid_field()->ConvectiveVel();
    velnp = poro_->fluid_field()->Velnp();
  }
  else
  {
    convel = volcoupl_fluidscatra_->apply_vector_mapping21(poro_->fluid_field()->ConvectiveVel());
    velnp = volcoupl_fluidscatra_->apply_vector_mapping21(poro_->fluid_field()->Velnp());
  }

  scatra_->ScaTraField()->set_velocity_field(convel,  // convective vel.
      Teuchos::null,                                  // acceleration
      velnp,                                          // velocity
      Teuchos::null,                                  // fsvel
      true                                            // set pressure
  );
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::set_mesh_disp()
{
  Teuchos::RCP<const Epetra_Vector> dispnp = Teuchos::null;

  if (matchinggrid_)
  {
    dispnp = poro_->fluid_field()->Dispnp();
  }
  else
  {
    dispnp = volcoupl_fluidscatra_->apply_vector_mapping21(fluid_field()->Dispnp());
  }

  scatra_->ScaTraField()->ApplyMeshMovement(dispnp);

  Teuchos::RCP<const Epetra_Vector> sdispnp = Teuchos::null;

  if (matchinggrid_)
  {
    sdispnp = structure_field()->Dispnp();
  }
  else
  {
    sdispnp = volcoupl_structurescatra_->apply_vector_mapping21(structure_field()->Dispnp());
  }

  scatra_->ScaTraField()->discretization()->set_state(1, "displacement", sdispnp);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::replace_dof_sets(
    Teuchos::RCP<Core::FE::Discretization> structdis,
    Teuchos::RCP<Core::FE::Discretization> fluiddis,
    Teuchos::RCP<Core::FE::Discretization> scatradis)
{
  if (matchinggrid_)
  {
    // the problem is two way coupled, thus each discretization must know the other discretization

    if (poro_field()->HasSubmeshes())
    {
      Teuchos::RCP<Core::DOFSets::DofSetGIDBasedWrapper> structsubdofset = Teuchos::rcp(
          new Core::DOFSets::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));
      Teuchos::RCP<Core::DOFSets::DofSetGIDBasedWrapper> fluidsubdofset = Teuchos::rcp(
          new Core::DOFSets::DofSetGIDBasedWrapper(fluiddis, fluiddis->GetDofSetProxy()));
      Teuchos::RCP<Core::DOFSets::DofSetGIDBasedWrapper> scatrasubdofset = Teuchos::rcp(
          new Core::DOFSets::DofSetGIDBasedWrapper(scatradis, scatradis->GetDofSetProxy()));

      scatradis->ReplaceDofSet(1, structsubdofset);
      scatradis->ReplaceDofSet(2, fluidsubdofset);
      structdis->ReplaceDofSet(2, scatrasubdofset);
      fluiddis->ReplaceDofSet(2, scatrasubdofset);
    }
    else
    {
      // build a proxy of the structure discretization for the scatra field
      Teuchos::RCP<Core::DOFSets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
      // build a proxy of the fluid discretization for the scatra field
      Teuchos::RCP<Core::DOFSets::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
      // build a proxy of the fluid discretization for the structure/fluid field
      Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

      scatradis->ReplaceDofSet(1, structdofset);
      scatradis->ReplaceDofSet(2, fluiddofset);
      structdis->ReplaceDofSet(2, scatradofset);
      fluiddis->ReplaceDofSet(2, scatradofset);
    }

    fluiddis->fill_complete();
    scatradis->fill_complete();
    structdis->fill_complete();
  }
  else
  {
    FOUR_C_THROW("restart for non-matching poro-scatra not yet tested. Feel free to try");
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void PoroElastScaTra::PoroScatraBase::setup_coupling(
    Teuchos::RCP<Core::FE::Discretization> structdis,
    Teuchos::RCP<Core::FE::Discretization> fluiddis,
    Teuchos::RCP<Core::FE::Discretization> scatradis)
{
  if (not matchinggrid_)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_structurescatra_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());
    volcoupl_fluidscatra_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

    std::pair<int, int> dofsets12_structurescatra = std::pair<int, int>(2, 0);
    std::pair<int, int> dofsets21_structurescatra = std::pair<int, int>(1, 0);
    std::pair<int, int> dofsets12_fluidscatra = std::pair<int, int>(2, 0);
    std::pair<int, int> dofsets21_fluidscatra = std::pair<int, int>(2, 0);

    // setup projection matrices (use default material strategy)
    volcoupl_structurescatra_->Init(Global::Problem::Instance()->NDim(), structdis, scatradis,
        nullptr, nullptr, &dofsets12_structurescatra, &dofsets21_structurescatra, Teuchos::null);
    volcoupl_fluidscatra_->Init(Global::Problem::Instance()->NDim(), fluiddis, scatradis, nullptr,
        nullptr, &dofsets12_fluidscatra, &dofsets21_fluidscatra, Teuchos::null);

    volcoupl_structurescatra_->Setup(Global::Problem::Instance()->VolmortarParams());
    volcoupl_fluidscatra_->Setup(Global::Problem::Instance()->VolmortarParams());
  }
}

FOUR_C_NAMESPACE_CLOSE
