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
#include "../drt_scatra/scatra_utils_clonestrategy.H"

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
    if ( timealgo != INPAR::SCATRA::timeint_one_step_theta )
      dserror("scalar transport in porous media is limited in functionality (only one-step-theta scheme possible)");

    INPAR::SCATRA::ConvForm convform
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatradyn,"CONVFORM");
    if ( convform != INPAR::SCATRA::convform_convective )
      dserror("The balance of mass is included in the formulation for scalar transport in porous media. "
          "Set 'CONVFORM' to 'convective' in the SCALAR TRANSPORT DYNAMIC section! ");

    INPAR::SCATRA::VelocityField velfield
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");
    if ( velfield != INPAR::SCATRA::velocity_Navier_Stokes )
      dserror("scalar transport is coupled with the porous medium. Set 'VELOCITYFIELD' to 'Navier_Stokes' in the SCALAR TRANSPORT DYNAMIC section! ");

    INPAR::SCATRA::ScaTraType scatratype
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::ScaTraType>(scatradyn,"SCATRATYPE");
    if ( scatratype != INPAR::SCATRA::scatratype_poro )
      dserror("Set 'SCATRATYPE' to 'Poroscatra' in the SCALAR TRANSPORT DYNAMIC section! ");
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
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,scatradis);
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
void POROELAST::PORO_SCATRA_Base::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    poro_->ReadRestart(restart);
    scatra_->ScaTraField().ReadRestart(restart);

    SetTimeStep(poro_->Time(), restart);
  }
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
