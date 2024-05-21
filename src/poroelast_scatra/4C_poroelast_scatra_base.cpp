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
#include "4C_discretization_dofset_gidbased_wrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_lib_utils_createdis.hpp"
#include "4C_poroelast_base.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_scatra_timint_implicit.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
POROELASTSCATRA::PoroScatraBase::PoroScatraBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : AlgorithmBase(comm, timeparams),
      matchinggrid_(CORE::UTILS::IntegralValue<bool>(
          GLOBAL::Problem::Instance()->PoroScatraControlParams(), "MATCHINGGRID")),
      volcoupl_structurescatra_(Teuchos::null),
      volcoupl_fluidscatra_(Teuchos::null)
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // do some checks
  {
    INPAR::SCATRA::TimeIntegrationScheme timealgo =
        CORE::UTILS::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn, "TIMEINTEGR");
    if (timealgo != INPAR::SCATRA::timeint_one_step_theta and
        timealgo != INPAR::SCATRA::timeint_stationary)
      FOUR_C_THROW(
          "scalar transport in porous media is limited in functionality (only one-step-theta "
          "scheme or stationary case possible)");

    //    INPAR::SCATRA::ConvForm convform
    //    = CORE::UTILS::IntegralValue<INPAR::SCATRA::ConvForm>(scatradyn,"CONVFORM");
    //    if ( convform != INPAR::SCATRA::convform_convective )
    //      FOUR_C_THROW("The balance of mass is included in the formulation for scalar transport in
    //      porous media. "
    //          "Set 'CONVFORM' to 'convective' in the SCALAR TRANSPORT DYNAMIC section! ");

    INPAR::SCATRA::VelocityField velfield =
        CORE::UTILS::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn, "VELOCITYFIELD");
    if (velfield != INPAR::SCATRA::velocity_Navier_Stokes)
      FOUR_C_THROW(
          "scalar transport is coupled with the porous medium. Set 'VELOCITYFIELD' to "
          "'Navier_Stokes' in the SCALAR TRANSPORT DYNAMIC section! ");

    //    bool skipinitder
    //    = CORE::UTILS::IntegralValue<int>(scatradyn,"SKIPINITDER");
    //    if ( not skipinitder )
    //      FOUR_C_THROW("Calculation of initial time derivative not yet supported for scalar
    //      transport in porous media. Set 'SKIPINITDER' to 'yes' in the SCALAR TRANSPORT DYNAMIC
    //      section! ");
  }

  // the problem is two way coupled, thus each discretization must know the other discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("porofluid");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  SetupCoupling(structdis, fluiddis, scatradis);
  // Create the two uncoupled subproblems.
  // 1. poro problem
  poro_ = POROELAST::UTILS::CreatePoroAlgorithm(
      timeparams, comm, false, POROELASTSCATRA::UTILS::BuildPoroScatraSplitter(structdis));

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
  // 2. scatra problem
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(
      timeparams, scatradyn, problem->SolverParams(linsolvernumber), "scatra", true));

  // now we can call Init() on the base algo.
  // time integrator is constructed and initialized inside.
  scatra_->Init();
  scatra_->ScaTraField()->SetNumberOfDofSetDisplacement(2);
  scatra_->ScaTraField()->SetNumberOfDofSetVelocity(2);

  // only now we must call Setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls Setup() on the time integrator inside.
  scatra_->Setup();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::SetupSystem() { poro_->SetupSystem(); }

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::TestResults(const Epetra_Comm& comm)
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  problem->AddFieldTest(poro_->StructureField()->CreateFieldTest());
  problem->AddFieldTest(poro_->FluidField()->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::SetPoroSolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::SetScatraSolution()
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
    phinp_s = volcoupl_structurescatra_->ApplyVectorMapping12(scatra_->ScaTraField()->Phinp());
    phinp_f = volcoupl_fluidscatra_->ApplyVectorMapping12(scatra_->ScaTraField()->Phinp());
    phin_s = volcoupl_structurescatra_->ApplyVectorMapping12(scatra_->ScaTraField()->Phin());
    phin_f = volcoupl_fluidscatra_->ApplyVectorMapping12(scatra_->ScaTraField()->Phin());
    phidtnp = volcoupl_fluidscatra_->ApplyVectorMapping12(scatra_->ScaTraField()->Phidtnp());
  }

  // porous structure
  poro_->StructureField()->Discretization()->SetState(2, "scalar", phinp_s);
  poro_->StructureField()->Discretization()->SetState(2, "scalarn", phin_s);

  // porous fluid
  poro_->FluidField()->SetIterScalarFields(phinp_f, phin_f, phidtnp,
      // scatra_->ScaTraField()->Discretization()
      poro_->FluidField()->Discretization(), 2);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::SetVelocityFields()
{
  Teuchos::RCP<const Epetra_Vector> convel = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> velnp = Teuchos::null;

  if (matchinggrid_)
  {
    convel = poro_->FluidField()->ConvectiveVel();
    velnp = poro_->FluidField()->Velnp();
  }
  else
  {
    convel = volcoupl_fluidscatra_->ApplyVectorMapping21(poro_->FluidField()->ConvectiveVel());
    velnp = volcoupl_fluidscatra_->ApplyVectorMapping21(poro_->FluidField()->Velnp());
  }

  scatra_->ScaTraField()->SetVelocityField(convel,  // convective vel.
      Teuchos::null,                                // acceleration
      velnp,                                        // velocity
      Teuchos::null,                                // fsvel
      true                                          // set pressure
  );
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::SetMeshDisp()
{
  Teuchos::RCP<const Epetra_Vector> dispnp = Teuchos::null;

  if (matchinggrid_)
  {
    dispnp = poro_->FluidField()->Dispnp();
  }
  else
  {
    dispnp = volcoupl_fluidscatra_->ApplyVectorMapping21(FluidField()->Dispnp());
  }

  scatra_->ScaTraField()->ApplyMeshMovement(dispnp);

  Teuchos::RCP<const Epetra_Vector> sdispnp = Teuchos::null;

  if (matchinggrid_)
  {
    sdispnp = StructureField()->Dispnp();
  }
  else
  {
    sdispnp = volcoupl_structurescatra_->ApplyVectorMapping21(StructureField()->Dispnp());
  }

  scatra_->ScaTraField()->Discretization()->SetState(1, "displacement", sdispnp);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::ReplaceDofSets(Teuchos::RCP<DRT::Discretization> structdis,
    Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> scatradis)
{
  if (matchinggrid_)
  {
    // the problem is two way coupled, thus each discretization must know the other discretization

    if (PoroField()->HasSubmeshes())
    {
      Teuchos::RCP<CORE::Dofsets::DofSetGIDBasedWrapper> structsubdofset = Teuchos::rcp(
          new CORE::Dofsets::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));
      Teuchos::RCP<CORE::Dofsets::DofSetGIDBasedWrapper> fluidsubdofset = Teuchos::rcp(
          new CORE::Dofsets::DofSetGIDBasedWrapper(fluiddis, fluiddis->GetDofSetProxy()));
      Teuchos::RCP<CORE::Dofsets::DofSetGIDBasedWrapper> scatrasubdofset = Teuchos::rcp(
          new CORE::Dofsets::DofSetGIDBasedWrapper(scatradis, scatradis->GetDofSetProxy()));

      scatradis->ReplaceDofSet(1, structsubdofset);
      scatradis->ReplaceDofSet(2, fluidsubdofset);
      structdis->ReplaceDofSet(2, scatrasubdofset);
      fluiddis->ReplaceDofSet(2, scatrasubdofset);
    }
    else
    {
      // build a proxy of the structure discretization for the scatra field
      Teuchos::RCP<CORE::Dofsets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
      // build a proxy of the fluid discretization for the scatra field
      Teuchos::RCP<CORE::Dofsets::DofSetInterface> fluiddofset = fluiddis->GetDofSetProxy();
      // build a proxy of the fluid discretization for the structure/fluid field
      Teuchos::RCP<CORE::Dofsets::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

      scatradis->ReplaceDofSet(1, structdofset);
      scatradis->ReplaceDofSet(2, fluiddofset);
      structdis->ReplaceDofSet(2, scatradofset);
      fluiddis->ReplaceDofSet(2, scatradofset);
    }

    fluiddis->FillComplete();
    scatradis->FillComplete();
    structdis->FillComplete();
  }
  else
  {
    FOUR_C_THROW("restart for non-matching poro-scatra not yet tested. Feel free to try");
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraBase::SetupCoupling(Teuchos::RCP<DRT::Discretization> structdis,
    Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> scatradis)
{
  if (not matchinggrid_)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_structurescatra_ = Teuchos::rcp(new CORE::ADAPTER::MortarVolCoupl());
    volcoupl_fluidscatra_ = Teuchos::rcp(new CORE::ADAPTER::MortarVolCoupl());

    std::pair<int, int> dofsets12_structurescatra = std::pair<int, int>(2, 0);
    std::pair<int, int> dofsets21_structurescatra = std::pair<int, int>(1, 0);
    std::pair<int, int> dofsets12_fluidscatra = std::pair<int, int>(2, 0);
    std::pair<int, int> dofsets21_fluidscatra = std::pair<int, int>(2, 0);

    // setup projection matrices (use default material strategy)
    volcoupl_structurescatra_->Init(GLOBAL::Problem::Instance()->NDim(), structdis, scatradis,
        nullptr, nullptr, &dofsets12_structurescatra, &dofsets21_structurescatra, Teuchos::null);
    volcoupl_fluidscatra_->Init(GLOBAL::Problem::Instance()->NDim(), fluiddis, scatradis, nullptr,
        nullptr, &dofsets12_fluidscatra, &dofsets21_fluidscatra, Teuchos::null);

    volcoupl_structurescatra_->Setup(GLOBAL::Problem::Instance()->VolmortarParams());
    volcoupl_fluidscatra_->Setup(GLOBAL::Problem::Instance()->VolmortarParams());
  }
}

FOUR_C_NAMESPACE_CLOSE
