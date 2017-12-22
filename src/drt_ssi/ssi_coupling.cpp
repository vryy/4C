/*----------------------------------------------------------------------*/
/*!
 \file ssi_coupling.cpp

 \brief helper classes for  scalar-structure coupling

 \level 3

   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "ssi_coupling.H"

//for coupling of nonmatching meshes
#include "../drt_adapter/adapter_coupling_volmortar.H"
#include "../drt_adapter/adapter_coupling_mortar.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_wrapper.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include"../drt_inpar/inpar_volmortar.H"
#include "../drt_volmortar/volmortar_utils.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_dofset_definedmapping_wrapper.H"
#include "../drt_lib/drt_dofset_gidbased_wrapper.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::Init(
    const int                         ndim,          /// dimension of the problem
    Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
    Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
    )
{
  SetIsSetup(false);

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (scatradis->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in scatra field");
  if (structdis->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");


  AssignMaterialPointers(structdis,scatradis);

  SetIsInit(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::Setup()
{
  CheckIsInit();

  SetIsSetup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        )
{
  const int numelements = scatradis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* scatratele = scatradis->lColElement(i);
    const int gid = scatratele->Id();

    DRT::Element* structele = structdis->gElement(gid);

    //for coupling we add the source material to the target element and vice versa
    scatratele->AddMaterial(structele->Material());
    structele->AddMaterial(scatratele->Material());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
    )
{
  scatra->ScaTraField()->ApplyMeshMovement(
      disp,
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
    Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
    )
{
  scatra->ScaTraField()->SetVelocityField(
        convvel, //convective vel.
        Teuchos::null, //acceleration
        vel, //velocity
        Teuchos::null, //fsvel
        1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetScalarField(
    Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
    Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
    )
{
  structure->Discretization()->SetState(1,"scalarfield",phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::Init(
    const int                             ndim,      /// dimension of the problem
    Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
    Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
    )
{
  SetIsSetup(false);

  // set pointers
  structdis_=structdis;
  scatradis_=scatradis;

  // set problem dimension
  problem_dimension_ = ndim;

  //first call FillComplete for single discretizations.
  //This way the physical dofs are numbered successively
  structdis_->FillComplete();
  scatradis_->FillComplete();

  //build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = scatradis->NumDof(0,scatradis->lRowNode(0));
  const int ndofperelement_scatra  = 0;
  const int ndofpernode_struct = structdis->NumDof(0,structdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber( ndofpernode_scatra, ndofperelement_scatra, 0, true ));
  if (structdis->AddDofSet(dofsetaux) != 1)
    dserror("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber( ndofpernode_struct, ndofperelement_struct, 0, true ));
  if (scatradis->AddDofSet(dofsetaux) != 1)
    dserror("unexpected dof sets in scatra field");

  //call AssignDegreesOfFreedom also for auxiliary dofsets
  //note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  structdis_->FillComplete(true, false,false);
  scatradis_->FillComplete(true, false,false);

  // setup mortar adapter for surface volume coupling
  adaptermeshtying_ = Teuchos::rcp(new ADAPTER::CouplingMortar());

  SetIsInit(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::Setup()
{
  CheckIsInit();

  std::vector<int> coupleddof(problem_dimension_, 1);
  // Setup of meshtying adapter
  adaptermeshtying_->Setup(
      structdis_,
      scatradis_,
      Teuchos::null,
      coupleddof,
      "SSICoupling",
      structdis_->Comm(),
      false,
      false,
      0,
      1
      );

  //extractor for coupled surface of structure discretization with surface scatra
  extractor_= Teuchos::rcp(new LINALG::MapExtractor(
      *structdis_->DofRowMap(0),
      adaptermeshtying_->MasterDofMap(),
      true));

  SetIsSetup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        )
{
  //nothing to do in this case, since
  //transferring scalar state to structure discretization not implemented for
  //transport on structural boundary. Only SolidToScatra coupling available.
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
    )
{
  scatra->ScaTraField()->ApplyMeshMovement(
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(disp)),
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
    Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
    )
{
  scatra->ScaTraField()->SetVelocityField(
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(convvel)), //convective vel.
      Teuchos::null, //acceleration
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(vel)), //velocity
      Teuchos::null, //fsvel
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetScalarField(
    Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
    Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
    )
{
  dserror("transferring scalar state to structure discretization not implemented for "
      "transport on structural boundary. Only SolidToScatra coupling available.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::Init(
    const int                         ndim,          /// dimension of the problem
    Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
    Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
    )
{
  SetIsSetup(false);

  //first call FillComplete for single discretizations.
  //This way the physical dofs are numbered successively
  structdis->FillComplete();
  scatradis->FillComplete();

  //build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = scatradis->NumDof(0,scatradis->lRowNode(0));
  const int ndofperelement_scatra  = 0;
  const int ndofpernode_struct = structdis->NumDof(0,structdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber( ndofpernode_scatra, ndofperelement_scatra, 0, true ));
  if (structdis->AddDofSet(dofsetaux) != 1)
    dserror("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber( ndofpernode_struct, ndofperelement_struct, 0, true ));
  if (scatradis->AddDofSet(dofsetaux) != 1)
    dserror("unexpected dof sets in scatra field");

  //call AssignDegreesOfFreedom also for auxiliary dofsets
  //note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  structdis->FillComplete(true, false,false);
  scatradis->FillComplete(true, false,false);

  // Scheme: non matching meshes --> volumetric mortar coupling...
  volcoupl_structurescatra_=Teuchos::rcp(new ADAPTER::MortarVolCoupl() );

  // init projection matrices (use default material strategy)
  volcoupl_structurescatra_->Init(
      structdis,
      scatradis);

  // parallel redistribution is performed in the global control
  // algorithm. We redistribute between Init(...) and Setup().
  //volcoupl_structurescatra_->Redistribute();

  SetIsInit(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::Setup()
{
  CheckIsInit();

  //setup projection matrices (use default material strategy)
  volcoupl_structurescatra_->Setup();

  SetIsSetup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        )
{
  volcoupl_structurescatra_->AssignMaterials(
      structdis,
      scatradis);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
    )
{
  scatra->ScaTraField()->ApplyMeshMovement(
      volcoupl_structurescatra_->ApplyVectorMapping21(disp),
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
    Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
    )
{
  scatra->ScaTraField()->SetVelocityField(
      volcoupl_structurescatra_->ApplyVectorMapping21(convvel), //convective vel.
      Teuchos::null, //acceleration
      volcoupl_structurescatra_->ApplyVectorMapping21(vel), //velocity
      Teuchos::null, //fsvel
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetScalarField(
    Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
    Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
    )
{
  structure->Discretization()->SetState(1,"scalarfield",volcoupl_structurescatra_->ApplyVectorMapping12(phi));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::Init(
    const int                         ndim,          /// dimension of the problem
    Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
    Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
    )
{
  SetIsSetup(false);

  // Note : We need to make sure that the parallel distribution of Volume and Boundary
  //        is the same externally! The best thing is if you do this in your *_dyn.cpp,
  //        i.e., your global control algorithm.

  {
    //get condition which defines the coupling on target discretization
    std::vector<DRT::Condition*> conds_struct;
    structdis->GetCondition("SSICouplingSolidToScatra", conds_struct);

    //get condition which defines the coupling on source discretization
    std::vector<DRT::Condition*> conds_scatra;
    scatradis->GetCondition("SSICouplingSolidToScatra", conds_scatra);

    // at least one condition needs to be defined on each discretization
    if(conds_struct.size() == 0 or conds_scatra.size() == 0)
      dserror("No coupling condition defined on one or both structure or scatra discretization!");

    std::set<int> couplingids;
    for (unsigned i=0; i<conds_struct.size(); ++i)
      couplingids.insert(conds_struct[i]->GetInt("coupling id"));

    Teuchos::RCP<DRT::DofSetGIDBasedWrapper> structgidmatchingdofset =
        Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(
            structdis,
            structdis->GetDofSetProxy()));

      Teuchos::RCP<DRT::DofSetDefinedMappingWrapper> newdofset_scatra =
          Teuchos::rcp(new DRT::DofSetDefinedMappingWrapper(
              structgidmatchingdofset,
              structdis,
              "SSICouplingSolidToScatra",
              couplingids));

      // add dofset and check if scatra field has 2 dofsets, so that coupling is possible
      if (scatradis->AddDofSet(newdofset_scatra)!=1)
        dserror("unexpected dof sets in scatra field");
  }

  {
    //get condition which defines the coupling on target discretization
    std::vector<DRT::Condition*> conds_struct;
    structdis->GetCondition("SSICouplingScatraToSolid", conds_struct);

    //get condition which defines the coupling on source discretization
    std::vector<DRT::Condition*> conds_scatra;
    scatradis->GetCondition("SSICouplingScatraToSolid", conds_scatra);

    // at least one condition needs to be defined on each discretization
    if(conds_struct.size() == 0 or conds_scatra.size() == 0)
      dserror("No coupling condition defined on one or both structure or scatra discretization!");

    std::set<int> couplingids;
    for (unsigned i=0; i<conds_struct.size(); ++i)
      couplingids.insert(conds_struct[i]->GetInt("coupling id"));

    Teuchos::RCP<DRT::DofSetGIDBasedWrapper> scatragidmatchingdofset =
        Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(
            scatradis,
            scatradis->GetDofSetProxy()));

    for (std::set<int>::iterator it=couplingids.begin(); it!=couplingids.end(); ++it)
    {
      std::set<int> tempset;
      tempset.insert(*it);

      Teuchos::RCP<DRT::DofSetDefinedMappingWrapper> newdofset_struct =
          Teuchos::rcp(new DRT::DofSetDefinedMappingWrapper(
              scatragidmatchingdofset,
              scatradis,
              "SSICouplingScatraToSolid",
              tempset));

      structdis->AddDofSet(newdofset_struct);
    }
  }

  // exchange material pointers for coupled material formulations
  AssignMaterialPointers(structdis,scatradis);

  SetIsInit(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::Setup()
{
  CheckIsInit();

  SetIsSetup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::AssignMaterialPointers(
        Teuchos::RCP<DRT::Discretization> structdis,     /// underlying structure discretization
        Teuchos::RCP<DRT::Discretization> scatradis      /// underlying scatra discretization
        )
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> disp             /// displacement field to set
    )
{
  scatra->ScaTraField()->ApplyMeshMovement(
      disp,
      1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, /// underlying scatra problem of the SSI problem
    Teuchos::RCP<const Epetra_Vector> convvel,         /// convective velocity field to set
    Teuchos::RCP<const Epetra_Vector> vel              /// velocity field to set
    )
{
  scatra->ScaTraField()->SetVelocityField(
        convvel, //convective vel.
        Teuchos::null, //acceleration
        vel, //velocity
        Teuchos::null, //fsvel
        1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetScalarField(
    Teuchos::RCP< ::ADAPTER::Structure> structure,     /// underlying structure of the SSI problem,
    Teuchos::RCP<const Epetra_Vector> phi              /// scalar field to set
    )
{
  structure->Discretization()->SetState(1,"scalarfield",phi);
}
