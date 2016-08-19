/*!------------------------------------------------------------------------------------------------*
 \file ssi_base.cpp

 \brief base class for all scalar structure algorithms

 \level 1

 \maintainer Anh-Tu Vuong
             vuong@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15264

 *------------------------------------------------------------------------------------------------*/

#include "ssi_base.H"

#include "ssi_partitioned.H"
#include "ssi_coupling.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
//for cloning
#include "../drt_lib/drt_utils_createdis.H"

#include"../drt_inpar/inpar_volmortar.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Base::SSI_Base(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    AlgorithmBase(comm, globaltimeparams),
    structure_(Teuchos::null),
    scatra_(Teuchos::null),
    zeros_(Teuchos::null),
    fieldcoupling_(DRT::INPUT::IntegralValue<INPAR::SSI::FieldCoupling>(DRT::Problem::Instance()->SSIControlParams(),"FIELDCOUPLING")),
    ssicoupling_(Teuchos::null)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
}


/*----------------------------------------------------------------------*
 | Setup this class                                         rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::Setup(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  //2.- Setup discretizations and coupling.
  SetupDiscretizationsAndFieldCoupling(comm,struct_disname, scatra_disname);

  //3.- Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis(struct_disname);

  // Set isale to false what should be the case in scatratosolid algorithm
  const INPAR::SSI::SolutionSchemeOverFields coupling
      = DRT::INPUT::IntegralValue<INPAR::SSI::SolutionSchemeOverFields>(problem->SSIControlParams(),"COUPALGO");

  bool isale = true;
  if(coupling == INPAR::SSI::ssi_OneWay_ScatraToSolid) isale = false;

  // determine which time params to use to build the single fields
  // in case of time stepping time params have to be read from single field sections
  // in case of equal timestep size for all fields the time params are controlled solely
  // by the problem section (e.g. ssi or cell dynamic)
  const Teuchos::ParameterList* structtimeparams = &globaltimeparams;
  const Teuchos::ParameterList* scatratimeparams = &globaltimeparams;
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->SSIControlParams(),"DIFFTIMESTEPSIZE"))
  {
    structtimeparams = &structparams;
    scatratimeparams = &scatraparams;
  }

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(*structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));

  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::Structure>(structure->StructureField());
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(*scatratimeparams,scatraparams,problem->SolverParams(linsolvernumber),scatra_disname,isale));
  zeros_ = LINALG::CreateVector(*structure_->DofRowMap(), true);

  return;
}


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 01/12  |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestart( int restart )
{
  if (restart)
  {
    structure_->ReadRestart(restart);

    const Teuchos::ParameterList& ssidyn = DRT::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure = DRT::INPUT::IntegralValue<int>(ssidyn,"RESTART_FROM_STRUCTURE");

    if (not restartfromstructure) //standard restart
    {
        scatra_->ScaTraField()->ReadRestart(restart);
    }
    else //restart from structure simulation
    {
      // Since there is no restart output for the scatra fiels available, we only have to fix the
      // time and step counter
      scatra_->ScaTraField()->SetTimeStep(structure_->TimeOld(),restart);
    }

    SetTimeStep(structure_->TimeOld(), restart);
  }

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  ssicoupling_->AssignMaterialPointers(structure_->Discretization(),scatra_->ScaTraField()->Discretization());

  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time (public)        AN, JH 10/14 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::ReadRestartfromTime( double restarttime )
{
  if ( restarttime > 0.0 )
  {
    const int restartstructure = SSI::Utils::CheckTimeStepping(structure_->Dt(), restarttime);
    const int restartscatra    = SSI::Utils::CheckTimeStepping(scatra_->ScaTraField()->Dt(), restarttime);

    structure_->ReadRestart(restartstructure);

    const Teuchos::ParameterList& ssidyn = DRT::Problem::Instance()->SSIControlParams();
    const bool restartfromstructure = DRT::INPUT::IntegralValue<int>(ssidyn,"RESTART_FROM_STRUCTURE");

    if (not restartfromstructure) //standard restart
    {
      scatra_->ScaTraField()->ReadRestart(restartscatra);
    }
    else //restart from structure simulation
    {
      // Since there is no restart output for the scatra fiels available, we only have to fix the
      // time and step counter
      scatra_->ScaTraField()->SetTimeStep(structure_->TimeOld(),restartscatra);
    }

    SetTimeStep(structure_->TimeOld(), restartstructure);
  }

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  ssicoupling_->AssignMaterialPointers(structure_->Discretization(),scatra_->ScaTraField()->Discretization());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& scatra_disname)
{
  // Scheme   : the structure discretization is received from the input. Then, an ale-scatra disc. is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);
  if(!structdis->Filled())
    structdis->FillComplete();
  if(!scatradis->Filled())
    scatradis->FillComplete();

  if (scatradis->NumGlobalNodes()==0)
  {
    // fill scatra discretization by cloning structure discretization
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,scatradis);

    // set implementation type
    for(int i=0; i<scatradis->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
      if(element == NULL)
        dserror("Invalid element type!");
      else
        element->SetImplType(DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(problem->SSIControlParams(),"SCATRATYPE"));
    }
  }
  else
  {
    if(fieldcoupling_==INPAR::SSI::coupling_volume_match)
      dserror("Reading a TRANSPORT discretization from the .dat file for the input parameter 'FIELDCOUPLING volume_matching' in the"
          "SSI CONTROL section is not supported! As this coupling relies on matching node (and sometimes element) IDs,"
          "the ScaTra discretization is cloned from the structure discretization. Delete the ScaTra discretization"
          "from your input file.");

    // copy conditions
    // this is actually only needed for copying TRANSPORT DIRICHLET/NEUMANN CONDITIONS
    // as standard DIRICHLET/NEUMANN CONDITIONS
    std::map<std::string,std::string> conditions_to_copy;
    SCATRA::ScatraFluidCloneStrategy clonestrategy;
    conditions_to_copy = clonestrategy.ConditionsToCopy();
    DRT::UTILS::DiscretizationCreatorBase creator;
    creator.CopyConditions(*scatradis,*scatradis,conditions_to_copy);

  }

  // setup the coupling adapters
  {
    DRT::Problem* problem = DRT::Problem::Instance();
    Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
    Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);

    //safety check
    {
      //check for ssi coupling condition
      std::vector<DRT::Condition*> ssicoupling;
      scatradis->GetCondition("SSICoupling",ssicoupling);
      const bool havessicoupling = (ssicoupling.size()>0);

      if(havessicoupling and
         (fieldcoupling_!= INPAR::SSI::coupling_boundary_nonmatch and fieldcoupling_!= INPAR::SSI::coupling_volumeboundary_match) )
        dserror("SSICoupling condition only valid in combination with FIELDCOUPLING set to 'boundary_nonmatching' "
            "or 'volumeboundary_matching' in SSI DYNAMIC section. ");

      if(fieldcoupling_==INPAR::SSI::coupling_volume_nonmatch)
      {
        const Teuchos::ParameterList& volmortarparams = DRT::Problem::Instance()->VolmortarParams();
        if (DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::CouplingType>(volmortarparams,"COUPLINGTYPE")!=
             INPAR::VOLMORTAR::couplingtype_coninter)
          dserror("Volmortar coupling only tested for consistent interpolation, "
              "i.e. 'COUPLINGTYPE consint' in VOLMORTAR COUPLING section. Try other couplings at own risk.");
      }
    }

    //build SSI coupling class
    switch(fieldcoupling_)
    {
    case INPAR::SSI::coupling_volume_match:
      ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolume());
      break;
    case INPAR::SSI::coupling_volume_nonmatch:
      ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingVolume());
      break;
    case INPAR::SSI::coupling_boundary_nonmatch:
      ssicoupling_ = Teuchos::rcp(new SSICouplingNonMatchingBoundary());
      break;
    case INPAR::SSI::coupling_volumeboundary_match:
      ssicoupling_ = Teuchos::rcp(new SSICouplingMatchingVolumeAndBoundary());
      break;
    default:
      dserror("unknown type of field coupling for SSI!");
      break;
    }

    // setup coupling objects including dof sets
    ssicoupling_->Setup(problem->NDim(),structdis,scatradis);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                                       Teuchos::RCP<const Epetra_Vector> vel )
{
  SetMeshDisp(disp);
  SetVelocityFields(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetScatraSolution( Teuchos::RCP<const Epetra_Vector> phi )
{
  ssicoupling_->SetScalarField(structure_,phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetVelocityFields( Teuchos::RCP<const Epetra_Vector> vel)
{
  ssicoupling_->SetVelocityFields(scatra_,zeros_,vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp )
{
  ssicoupling_->SetMeshDisp(scatra_,disp);
}


