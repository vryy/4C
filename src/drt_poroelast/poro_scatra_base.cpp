/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_base.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "poro_scatra_base.H"
#include "poro_base.H"

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_inpar/inpar_scatra.H"
#include "poroelast_utils.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Base::PORO_SCATRA_Base(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams):
    AlgorithmBase(comm, timeparams)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  const Teuchos::ParameterList& scatradyn  = problem->ScalarTransportDynamicParams();

  //do some checks
  {
    INPAR::SCATRA::TimeIntegrationScheme timealgo
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");
    if ( timealgo != INPAR::SCATRA::timeint_one_step_theta and timealgo != INPAR::SCATRA::timeint_stationary)
      dserror("scalar transport in porous media is limited in functionality (only one-step-theta scheme or stationary case possible)");

//    INPAR::SCATRA::ConvForm convform
//    = DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatradyn,"CONVFORM");
//    if ( convform != INPAR::SCATRA::convform_convective )
//      dserror("The balance of mass is included in the formulation for scalar transport in porous media. "
//          "Set 'CONVFORM' to 'convective' in the SCALAR TRANSPORT DYNAMIC section! ");

    INPAR::SCATRA::VelocityField velfield
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");
    if ( velfield != INPAR::SCATRA::velocity_Navier_Stokes )
      dserror("scalar transport is coupled with the porous medium. Set 'VELOCITYFIELD' to 'Navier_Stokes' in the SCALAR TRANSPORT DYNAMIC section! ");

    bool skipinitder
    = DRT::INPUT::IntegralValue<int>(scatradyn,"SKIPINITDER");
    if ( not skipinitder )
      dserror("Calculation of initial time derivative not yet supported for scalar transport in porous media. Set 'SKIPINITDER' to 'yes' in the SCALAR TRANSPORT DYNAMIC section! ");
  }

  //1.- Setup discretizations.
  SetupDiscretizations(comm);

  //2.- Create the two uncoupled subproblems.
  poro_ = POROELAST::UTILS::CreatePoroAlgorithm(timeparams, comm);
  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(timeparams,true,"scatra",problem->SolverParams(linsolvernumber)));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::SetupDiscretizations(const Epetra_Comm& comm)
{
  // Scheme    : the structure discretization is received from the input. Then, an ale-fluid disc.is cloned from the struct. one.
  //  After that, an ale-scatra disc. is cloned from the structure discretization.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  //1.2.-Set degrees of freedom in the str. discretization
  if (!structdis->Filled() or !structdis->HaveDofs())
    structdis->FillComplete();

  //2.- Access the fluid discretization, make sure it's empty, and fill it cloning the structural one.
  if (!fluiddis->Filled())
    fluiddis->FillComplete();

  if (structdis->NumGlobalNodes() == 0)
    dserror("Structure discretization is empty!");

  if (fluiddis->NumGlobalNodes()==0)
  {
    // create the fluid discretization
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,fluiddis);
  }
  else
    dserror("Structure AND Fluid discretization present. This is not supported.");

  //3.-Access the scatra discretization, make sure it's empty, and fill it by cloning the structural one.
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  if(!scatradis->Filled())
    scatradis->FillComplete();

  if (scatradis->NumGlobalNodes()==0)
  {
    // create the fluid scatra discretization
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroScatraCloneStrategy>(structdis,scatradis);
  }
  else
  dserror("Structure AND ScaTra discretization present. This is not supported.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::SetupSystem()
{
  poro_->SetupSystem();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(poro_->StructureField()->CreateFieldTest());
  problem->AddFieldTest(poro_->FluidField()->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::SetPoroSolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::SetScatraSolution()
{
  //porous structure
  poro_->StructureField()->ApplyCouplingState(scatra_->ScaTraField().Phinp(),"scalar",2);
  poro_->StructureField()->ApplyCouplingState(scatra_->ScaTraField().Phin(),"scalarn",2);

  //porous fluid
  poro_->FluidField()->SetIterLomaFields(scatra_->ScaTraField().Phinp(),
                                        scatra_->ScaTraField().Phin(),
                                        scatra_->ScaTraField().Phidtnp(),
                                        Teuchos::null,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        scatra_->ScaTraField().Discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::SetVelocityFields()
{
  scatra_->ScaTraField().SetVelocityField(
      poro_->FluidField()->ConvectiveVel(), //convective vel.
      Teuchos::null, //acceleration
      poro_->FluidField()->Velnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      poro_->FluidField()->Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::SetMeshDisp()
{
  scatra_->ScaTraField().ApplyMeshMovement(
      poro_->FluidField()->Dispnp(),
      poro_->FluidField()->Discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Base::AddDofSets(bool replace)
{
  // the problem is two way coupled, thus each discretization must know the other discretization
  Teuchos::RCP<DRT::DofSet> structdofset = Teuchos::null;
  Teuchos::RCP<DRT::DofSet> fluiddofset = Teuchos::null;
  Teuchos::RCP<DRT::DofSet> scatradofset = Teuchos::null;

  //get discretizations
  Teuchos::RCP<DRT::Discretization> structdis = PoroField()->StructureField()->Discretization();
  Teuchos::RCP<DRT::Discretization> fluiddis = PoroField()->FluidField()->Discretization();
  Teuchos::RCP<DRT::Discretization> scatradis = ScatraField().Discretization();

  if(PoroField()->HasSubmeshes())
  {
    // build a proxy of the structure discretization for the scatra field
    structdofset = structdis->GetDofSetProxy(structdis->NodeColMap(),structdis->ElementColMap());
    // build a proxy of the fluid discretization for the scatra field
    fluiddofset = fluiddis->GetDofSetProxy(fluiddis->NodeColMap(),fluiddis->ElementColMap());
    // build a proxy of the scatra discretization for the structure/fluid field
    scatradofset = scatradis->GetDofSetProxy(scatradis->NodeColMap(),scatradis->ElementColMap());
  }
  else
  {
    // build a proxy of the structure discretization for the scatra field
    structdofset = structdis->GetDofSetProxy();
    // build a proxy of the fluid discretization for the scatra field
    fluiddofset = fluiddis->GetDofSetProxy();
    // build a proxy of the fluid discretization for the structure/fluid field
    scatradofset = scatradis->GetDofSetProxy();
  }

  if(not replace)
  {
    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (scatradis->AddDofSet(structdofset) != 1)
      dserror("unexpected dof sets in scatra field");
    if (structdis->AddDofSet(scatradofset)!=2)
      dserror("unexpected dof sets in structure field");
    if (fluiddis->AddDofSet(scatradofset)!=2)
      dserror("unexpected dof sets in fluid field");
  }
  else
  {
    scatradis->ReplaceDofSet(1,structdofset);
    structdis->ReplaceDofSet(2,scatradofset);
    fluiddis->ReplaceDofSet(2,scatradofset);
  }
}
