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
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/adapter_coupling_mortar.H"

#include "../drt_lib/drt_globalproblem.H"
//for cloning
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_scatra_ele/scatra_ele.H"

//for coupling of nonmatching meshes
#include "../drt_adapter/adapter_coupling_volmortar.H"
#include "../drt_volmortar/volmortar_utils.H"
#include"../drt_inpar/inpar_volmortar.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"

#include "../drt_particle/binning_strategy.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Base::SSI_Base(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname):
    AlgorithmBase(comm, globaltimeparams),
    structure_(Teuchos::null),
    scatra_(Teuchos::null),
    zeros_(Teuchos::null),
    adapter_meshtying_(Teuchos::null),
    extractor_meshtying_(Teuchos::null),
    volcoupl_structurescatra_(Teuchos::null),
    fieldcoupling_(DRT::INPUT::IntegralValue<INPAR::SSI::FieldCoupling>(DRT::Problem::Instance()->SSIControlParams(),"FIELDCOUPLING"))
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  //2.- Setup discretizations and coupling.
  SetupDiscretizations(comm,struct_disname, scatra_disname);
  SetupFieldCoupling(struct_disname, scatra_disname);

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

  //some safety checks
  if (volcoupl_structurescatra_==Teuchos::null and fieldcoupling_==INPAR::SSI::coupling_volmortar )
    dserror("Something went terrible wrong. Sorry about this!");

  if (adapter_meshtying_==Teuchos::null and fieldcoupling_==INPAR::SSI::coupling_meshtying )
    dserror("Something went terrible wrong. Sorry about this!");
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
  switch(fieldcoupling_)
  {
  case INPAR::SSI::coupling_match:
    SetMaterialPointersMatchingGrid( structure_->Discretization(),scatra_->ScaTraField()->Discretization() );
    break;
  case INPAR::SSI::coupling_volmortar:
    volcoupl_structurescatra_->AssignMaterials(
                      structure_->Discretization(),
                      scatra_->ScaTraField()->Discretization());
    break;
  case INPAR::SSI::coupling_meshtying:
    //nothing to do in this case, since
    //transfering scalar state to structure discretization not implemented for
    //transport on structural boundary. Only SolidToScatra coupling available.
    break;
  default:
    dserror("unknown field coupling type in SetScatraSolution()");
    break;
  }

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
  switch(fieldcoupling_)
  {
  case INPAR::SSI::coupling_match:
    SetMaterialPointersMatchingGrid( structure_->Discretization(),scatra_->ScaTraField()->Discretization() );
    break;
  case INPAR::SSI::coupling_volmortar:
    volcoupl_structurescatra_->AssignMaterials(
                      structure_->Discretization(),
                      scatra_->ScaTraField()->Discretization());
    break;
  case INPAR::SSI::coupling_meshtying:
    //nothing to do in this case, since
    //transfering scalar state to structure discretization not implemented for
    //transport on structural boundary. Only SolidToScatra coupling available.
    break;
  default:
    dserror("unknown field coupling type in SetScatraSolution()");
    break;
  }

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
void SSI::SSI_Base::SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string scatra_disname)
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
    if (fieldcoupling_==INPAR::SSI::coupling_match)
      dserror("If you have different structure and scatra meshes use volmortar! If you have matching grids clone your structure mesh!");

    std::map<std::string,std::string> conditions_to_copy;
    SCATRA::ScatraFluidCloneStrategy clonestrategy;
    conditions_to_copy = clonestrategy.ConditionsToCopy();
    DRT::UTILS::DiscretizationCreatorBase creator;
    creator.CopyConditions(scatradis,scatradis,conditions_to_copy);
  }

  if(fieldcoupling_==INPAR::SSI::coupling_match)
  {
    structdis->FillComplete();
    scatradis->FillComplete();

    // build a proxy of the structure discretization for the scatra field
    Teuchos::RCP<DRT::DofSet> structdofset = structdis->GetDofSetProxy();
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<DRT::DofSet> scatradofset = scatradis->GetDofSetProxy();

    // check if scatra field has 2 discretizations, so that coupling is possible
    if (scatradis->AddDofSet(structdofset)!=1)
      dserror("unexpected dof sets in scatra field");
    if (structdis->AddDofSet(scatradofset)!=1)
      dserror("unexpected dof sets in structure field");

    SetMaterialPointersMatchingGrid(structdis,scatradis);
  }
  else
  {
    //first call FillComplete for single discretizations.
    //This way the physical dofs are numbered successively
    structdis->FillComplete();
    scatradis->FillComplete();

    //build auxiliary dofsets, i.e. pseudo dofs on each discretization
    const int ndofpernode_scatra = scatradis->NumDof(0,scatradis->lRowNode(0));
    const int ndofperelement_scatra  = 0;
    const int ndofpernode_struct = structdis->NumDof(0,structdis->lRowNode(0));
    const int ndofperelement_struct = 0;
    if (structdis->BuildDofSetAuxProxy(ndofpernode_scatra, ndofperelement_scatra, 0, true ) != 1)
      dserror("unexpected dof sets in structure field");
    if (scatradis->BuildDofSetAuxProxy(ndofpernode_struct, ndofperelement_struct, 0, true) != 1)
      dserror("unexpected dof sets in scatra field");

    //call AssignDegreesOfFreedom also for auxiliary dofsets
    //note: the order of FillComplete() calls determines the gid numbering!
    // 1. structure dofs
    // 2. scatra dofs
    // 3. structure auxiliary dofs
    // 4. scatra auxiliary dofs
    structdis->FillComplete(true, false,false);
    scatradis->FillComplete(true, false,false);
  }
}

/*----------------------------------------------------------------------*
 | Set material pointers for matching grids                  thon 01/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Base::SetMaterialPointersMatchingGrid(
    Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis
    )
{
  const int numelements = targetdis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* targetele = targetdis->lColElement(i);
    const int gid = targetele->Id();

    DRT::Element* sourceele = sourcedis->gElement(gid);

    //for coupling we add the source material to the target element and vice versa
    targetele->AddMaterial(sourceele->Material());
    sourceele->AddMaterial(targetele->Material());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                                       Teuchos::RCP<const Epetra_Vector> vel )
{
  SetMeshDisp(disp);
  SetVelocityField(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetScatraSolution( Teuchos::RCP<const Epetra_Vector> phi )
{
  structure_->Discretization()->SetState(1,"temperature",ScatraToStructure(phi));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetVelocityField( Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra_->ScaTraField()->SetVelocityField(
      StructureToScatra(zeros_), //convective vel.
      Teuchos::null, //acceleration
      StructureToScatra(vel), //velocity
      Teuchos::null, //fsvel
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp )
{
  scatra_->ScaTraField()->ApplyMeshMovement(StructureToScatra(disp),1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Base::SetupFieldCoupling(const std::string struct_disname, const std::string scatra_disname)
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

    if(havessicoupling and fieldcoupling_!=INPAR::SSI::coupling_meshtying)
      dserror("SSICoupling condition only valid in combination with FIELDCOUPLING 'meshtying' in SSI DYNAMIC section. "
          "If you want volume and surface coupling, FIELDCOUPLING 'volmortar' and "
          "a Mortar/S2I condition (and no SSICoupling condition) for the volume-surface-scatra coupling.");

    if(fieldcoupling_==INPAR::SSI::coupling_volmortar)
    {
      const Teuchos::ParameterList& volmortarparams = DRT::Problem::Instance()->VolmortarParams();
      if (DRT::INPUT::IntegralValue<INPAR::VOLMORTAR::CouplingType>(volmortarparams,"COUPLINGTYPE")!=
           INPAR::VOLMORTAR::couplingtype_coninter)
        dserror("Volmortar coupling only tested for consistent interpolation, "
            "i.e. 'COUPLINGTYPE consint' in VOLMORTAR COUPLING section. Try other couplings at own risk.");
    }
  }

  if(fieldcoupling_==INPAR::SSI::coupling_meshtying)
  {
    adapter_meshtying_ = Teuchos::rcp(new ADAPTER::CouplingMortar());

    std::vector<int> coupleddof(problem->NDim(), 1);
    // Setup of meshtying adapter
    adapter_meshtying_->Setup(structdis,
                            scatradis,
                            Teuchos::null,
                            coupleddof,
                            "SSICoupling",
                            structdis->Comm(),
                            false,
                            false,
                            0,
                            1
                            );

    extractor_meshtying_= Teuchos::rcp(new LINALG::MapExtractor(*structdis->DofRowMap(0),adapter_meshtying_->MasterDofMap(),true));
  }
  else if(fieldcoupling_==INPAR::SSI::coupling_volmortar)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_structurescatra_=Teuchos::rcp(new ADAPTER::MortarVolCoupl() );

    //setup projection matrices (use default material strategy)
    volcoupl_structurescatra_->Setup( structdis,
                                      scatradis);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> SSI::SSI_Base::StructureToScatra(const Teuchos::RCP<const Epetra_Vector> structurevector) const
{
  switch(fieldcoupling_)
  {
  case INPAR::SSI::coupling_match:
  {
    return structurevector;
    break;
  }
  case INPAR::SSI::coupling_volmortar:
  {
    return volcoupl_structurescatra_->ApplyVectorMapping21(structurevector);
    break;
  }
  case INPAR::SSI::coupling_meshtying:
  {
    return adapter_meshtying_->MasterToSlave(extractor_meshtying_->ExtractCondVector(structurevector));
    break;
  }
  default:
  {
    dserror("unknown field coupling type");
    return Teuchos::null;
    break;
  }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> SSI::SSI_Base::ScatraToStructure(const Teuchos::RCP<const Epetra_Vector> scatravector) const
{
  switch(fieldcoupling_)
  {
  case INPAR::SSI::coupling_match:
  {
    return scatravector;
    break;
  }
  case INPAR::SSI::coupling_volmortar:
  {
    return volcoupl_structurescatra_->ApplyVectorMapping12(scatravector);
    break;
  }
  case INPAR::SSI::coupling_meshtying:
  {
    return adapter_meshtying_->SlaveToMaster(extractor_meshtying_->ExtractCondVector(scatravector));
    break;
  }
  default:
  {
    dserror("unknown field coupling type");
    return Teuchos::null;
    break;
  }
  }
}
