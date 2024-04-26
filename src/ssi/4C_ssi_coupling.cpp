/*----------------------------------------------------------------------*/
/*! \file
 \brief helper classes for  scalar-structure coupling

 \level 3


 *----------------------------------------------------------------------*/

#include "4C_ssi_coupling.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_coupling_volmortar_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_condition_utils.hpp"
#include "4C_lib_dofset_definedmapping_wrapper.hpp"
#include "4C_lib_dofset_gidbased_wrapper.hpp"
#include "4C_lib_dofset_predefineddofnumber.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::Init(const int ndim,
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  SetIsSetup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();
  auto scatradis = scatra_integrator->Discretization();
  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // add proxy dofssets of other fields to discretizations and check if number of dofsets is correct
  if (scatradis->AddDofSet(structdofset) != ++scatra_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  scatra_integrator->SetNumberOfDofSetDisplacement(scatra_dofset_counter);
  scatra_integrator->SetNumberOfDofSetVelocity(scatra_dofset_counter);
  if (structdis->AddDofSet(scatradofset) != ++structure_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in structure field");

  if (GLOBAL::Problem::Instance()->ELCHControlParams().get<int>("TEMPERATURE_FROM_FUNCT") != -1)
  {
    const int numDofsPerNodeTemp = 1;  // defined by temperature field

    Teuchos::RCP<DRT::DofSetInterface> dofsettemp =
        Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(numDofsPerNodeTemp, 0, 0, true));
    if (structdis->AddDofSet(dofsettemp) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
  }

  if (GLOBAL::Problem::Instance()->Materials()->FirstIdByType(
          CORE::Materials::m_scatra_multiscale) != -1 or
      GLOBAL::Problem::Instance()->Materials()->FirstIdByType(
          CORE::Materials::m_newman_multiscale) != -1)
  {
    auto dofsetmicro = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(1, 0, 0, true));
    if (scatradis->AddDofSet(dofsetmicro) != ++scatra_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra field");
    scatra_integrator->SetNumberOfDofSetMicroScale(scatra_dofset_counter);
    if (structdis->AddDofSet(dofsetmicro) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
  }

  if (ssi_base->IsS2IKineticsWithPseudoContact())
  {
    const int numDofsPerNodeStresses = 6;
    Teuchos::RCP<DRT::DofSetInterface> dofsetstresses =
        Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(numDofsPerNodeStresses, 0, 0, true));
    if (structdis->AddDofSet(dofsetstresses) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
    if (scatradis->AddDofSet(structdis->GetDofSetProxy(structure_dofset_counter)) !=
        ++scatra_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra field");
    scatra_integrator->SetNumberOfDofSetTwoTensorQuantity(scatra_dofset_counter);
  }

  AssignMaterialPointers(structdis, scatradis);

  SetIsInit(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::Setup()
{
  CheckIsInit();

  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::AssignMaterialPointers(
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<DRT::Discretization> scatradis)
{
  const int numelements = scatradis->NumMyColElements();

  for (int i = 0; i < numelements; ++i)
  {
    DRT::Element* scatratele = scatradis->lColElement(i);
    const int gid = scatratele->Id();

    DRT::Element* structele = structdis->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    scatratele->AddMaterial(structele->Material());
    structele->AddMaterial(scatratele->Material());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetMechanicalStressState(
    DRT::Discretization& scatradis, Teuchos::RCP<const Epetra_Vector> stress_state, unsigned nds)
{
  scatradis.SetState(nds, "mechanicalStressState", stress_state);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(disp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->SetVelocityField(convvel,  // convective vel.
      Teuchos::null,                                // acceleration
      vel,                                          // velocity
      Teuchos::null                                 // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetScalarField(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.SetState(nds, "scalarfield", phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetScalarFieldMicro(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.SetState(nds, "MicroCon", phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetTemperatureField(
    DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp)
{
  structdis.SetState(2, "tempfield", temp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetTemperatureField(
    DRT::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp)
{
  structdis.SetState(2, "tempfield", temp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::Init(const int ndim,
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  SetIsSetup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();

  // set pointers
  structdis_ = structdis;
  scatradis_ = scatra_integrator->Discretization();

  // set problem dimension
  problem_dimension_ = ndim;

  // first call FillComplete for single discretizations.
  // This way the physical dofs are numbered successively
  structdis_->FillComplete();
  scatradis_->FillComplete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = scatradis_->NumDof(0, scatradis_->lRowNode(0));
  const int ndofperelement_scatra = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(
      new DRT::DofSetPredefinedDoFNumber(ndofpernode_scatra, ndofperelement_scatra, 0, true));
  if (structdis->AddDofSet(dofsetaux) != ++structure_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(
      new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct, ndofperelement_struct, 0, true));
  if (scatradis_->AddDofSet(dofsetaux) != ++scatra_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  scatra_integrator->SetNumberOfDofSetDisplacement(scatra_dofset_counter);
  scatra_integrator->SetNumberOfDofSetVelocity(scatra_dofset_counter);

  // call AssignDegreesOfFreedom also for auxiliary dofsets
  // note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  structdis_->FillComplete(true, false, false);
  scatradis_->FillComplete(true, false, false);

  // setup mortar adapter for surface volume coupling
  adaptermeshtying_ = Teuchos::rcp(new CORE::ADAPTER::CouplingMortar(
      GLOBAL::Problem::Instance()->NDim(), GLOBAL::Problem::Instance()->MortarCouplingParams(),
      GLOBAL::Problem::Instance()->ContactDynamicParams(),
      GLOBAL::Problem::Instance()->SpatialApproximationType()));

  SetIsInit(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::Setup()
{
  CheckIsInit();

  std::vector<int> coupleddof(problem_dimension_, 1);
  // Setup of meshtying adapter
  adaptermeshtying_->Setup(structdis_, scatradis_, Teuchos::null, coupleddof, "SSICoupling",
      structdis_->Comm(), false, false, 0, 1);

  // extractor for coupled surface of structure discretization with surface scatra
  extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(
      *structdis_->DofRowMap(0), adaptermeshtying_->MasterDofMap(), true));

  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::AssignMaterialPointers(
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<DRT::Discretization> scatradis)
{
  // nothing to do in this case, since
  // transferring scalar state to structure discretization not implemented for
  // transport on structural boundary. Only SolidToScatra coupling available.
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(disp)));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->SetVelocityField(
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(convvel)),  // convective vel.
      Teuchos::null,                                                             // acceleration
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(vel)),      // velocity
      Teuchos::null                                                              // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetScalarField(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW(
      "transferring scalar state to structure discretization not implemented for "
      "transport on structural boundary. Only SolidToScatra coupling available.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetScalarFieldMicro(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW("transferring micro scalar state to structure discretization not implemented.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::Init(const int ndim,
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  SetIsSetup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();
  auto scatradis = scatra_integrator->Discretization();
  // first call FillComplete for single discretizations.
  // This way the physical dofs are numbered successively
  structdis->FillComplete();
  scatradis->FillComplete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = scatradis->NumDof(0, scatradis->lRowNode(0));
  const int ndofperelement_scatra = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(
      new DRT::DofSetPredefinedDoFNumber(ndofpernode_scatra, ndofperelement_scatra, 0, true));
  if (structdis->AddDofSet(dofsetaux) != ++structure_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(
      new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct, ndofperelement_struct, 0, true));
  if (scatradis->AddDofSet(dofsetaux) != ++scatra_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  scatra_integrator->SetNumberOfDofSetDisplacement(scatra_dofset_counter);
  scatra_integrator->SetNumberOfDofSetVelocity(scatra_dofset_counter);

  // call AssignDegreesOfFreedom also for auxiliary dofsets
  // note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  structdis->FillComplete(true, false, false);
  scatradis->FillComplete(true, false, false);

  // Scheme: non matching meshes --> volumetric mortar coupling...
  volcoupl_structurescatra_ = Teuchos::rcp(new CORE::ADAPTER::MortarVolCoupl());

  // init projection matrices (use default material strategy)
  volcoupl_structurescatra_->Init(ndim, structdis, scatradis);

  // parallel redistribution is performed in the global control
  // algorithm. We redistribute between Init(...) and Setup().
  // volcoupl_structurescatra_->Redistribute();

  SetIsInit(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::Setup()
{
  CheckIsInit();

  // setup projection matrices (use default material strategy)
  volcoupl_structurescatra_->Setup(GLOBAL::Problem::Instance()->VolmortarParams());

  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::AssignMaterialPointers(
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<DRT::Discretization> scatradis)
{
  volcoupl_structurescatra_->AssignMaterials(
      structdis, scatradis, GLOBAL::Problem::Instance()->VolmortarParams());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(volcoupl_structurescatra_->ApplyVectorMapping21(disp));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->SetVelocityField(
      volcoupl_structurescatra_->ApplyVectorMapping21(convvel),  // convective vel.
      Teuchos::null,                                             // acceleration
      volcoupl_structurescatra_->ApplyVectorMapping21(vel),      // velocity
      Teuchos::null                                              // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetScalarField(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.SetState(nds, "scalarfield", volcoupl_structurescatra_->ApplyVectorMapping12(phi));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetScalarFieldMicro(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW("transferring micro scalar state to structure discretization not implemented.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::Init(const int ndim,
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  SetIsSetup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();
  auto scatradis = scatra_integrator->Discretization();


  // Note : We need to make sure that the parallel distribution of Volume and Boundary
  //        is the same externally! The best thing is if you do this in your *_dyn.cpp,
  //        i.e., your global control algorithm.

  if (ssi_base->ScaTraManifoldBaseAlgorithm() == Teuchos::null)
  {
    {
      // get condition which defines the coupling on target discretization
      std::vector<DRT::Condition*> conds_struct;
      structdis->GetCondition("SSICouplingSolidToScatra", conds_struct);

      // get condition which defines the coupling on source discretization
      std::vector<DRT::Condition*> conds_scatra;
      scatradis->GetCondition("SSICouplingSolidToScatra", conds_scatra);

      // at least one condition needs to be defined on each discretization
      if (conds_struct.size() == 0 or conds_scatra.size() == 0)
        FOUR_C_THROW(
            "No coupling condition defined on one or both structure or scatra discretization!");

      std::set<int> couplingids;
      for (auto& cond_struct : conds_struct)
        couplingids.insert(cond_struct->Get<int>("coupling id"));

      Teuchos::RCP<DRT::DofSetGIDBasedWrapper> structgidmatchingdofset =
          Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));

      Teuchos::RCP<DRT::DofSetDefinedMappingWrapper> newdofset_scatra =
          Teuchos::rcp(new DRT::DofSetDefinedMappingWrapper(
              structgidmatchingdofset, structdis, "SSICouplingSolidToScatra", couplingids));

      // add dofset and check if scatra field has 2 dofsets, so that coupling is possible
      if (scatradis->AddDofSet(newdofset_scatra) != ++scatra_dofset_counter)
        FOUR_C_THROW("unexpected dof sets in scatra field");
      scatra_integrator->SetNumberOfDofSetDisplacement(scatra_dofset_counter);
      scatra_integrator->SetNumberOfDofSetVelocity(scatra_dofset_counter);
    }

    {
      // get condition which defines the coupling on target discretization
      std::vector<DRT::Condition*> conds_struct;
      structdis->GetCondition("SSICouplingScatraToSolid", conds_struct);

      // get condition which defines the coupling on source discretization
      std::vector<DRT::Condition*> conds_scatra;
      scatradis->GetCondition("SSICouplingScatraToSolid", conds_scatra);

      // at least one condition needs to be defined on each discretization
      if (conds_struct.size() == 0 or conds_scatra.size() == 0)
        FOUR_C_THROW(
            "No coupling condition defined on one or both structure or scatra discretization!");

      std::set<int> couplingids;
      for (auto& cond_struct : conds_struct)
        couplingids.insert(cond_struct->Get<int>("coupling id"));

      Teuchos::RCP<DRT::DofSetGIDBasedWrapper> scatragidmatchingdofset =
          Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(scatradis, scatradis->GetDofSetProxy()));

      for (int couplingid : couplingids)
      {
        std::set<int> tempset;
        tempset.insert(couplingid);

        Teuchos::RCP<DRT::DofSetDefinedMappingWrapper> newdofset_struct =
            Teuchos::rcp(new DRT::DofSetDefinedMappingWrapper(
                scatragidmatchingdofset, scatradis, "SSICouplingScatraToSolid", tempset));

        structdis->AddDofSet(newdofset_struct);
      }
    }
  }
  else
  {
    int scatra_manifold_dofset_counter(0);

    auto scatra_manifold_integrator = ssi_base->ScaTraManifold();
    auto scatra_manifold_dis = scatra_manifold_integrator->Discretization();

    // build a proxy of the structure discretization for the other fields
    auto structdofset = structdis->GetDofSetProxy();
    // build a proxy of the scatra discretization for the other fields
    auto scatradofset = scatradis->GetDofSetProxy();

    // add proxy dofssets of other fields to discretizations and check if number of dofsets is
    // correct
    if (scatradis->AddDofSet(structdofset) != ++scatra_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra field");
    scatra_integrator->SetNumberOfDofSetDisplacement(scatra_dofset_counter);
    scatra_integrator->SetNumberOfDofSetVelocity(scatra_dofset_counter);
    if (structdis->AddDofSet(scatradofset) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");

    // set dummy coupling id, as coupling between scatra_manifold dis and structdis/scatradis should
    // be setup for all conditions
    std::set<int> couplingids;
    couplingids.insert(0);

    auto structgidmatchingdofset =
        Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));

    auto proxy_structure_scatramanifold = Teuchos::rcp(new DRT::DofSetDefinedMappingWrapper(
        structgidmatchingdofset, scatra_manifold_dis, "SSISurfaceManifold", couplingids));

    if (scatra_manifold_dis->AddDofSet(proxy_structure_scatramanifold) !=
        ++scatra_manifold_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra manifold field");
    scatra_manifold_integrator->SetNumberOfDofSetDisplacement(scatra_manifold_dofset_counter);
    scatra_manifold_integrator->SetNumberOfDofSetVelocity(scatra_manifold_dofset_counter);
  }

  if (GLOBAL::Problem::Instance()->ELCHControlParams().get<int>("TEMPERATURE_FROM_FUNCT") != -1)
  {
    const int numDofsPerNodeTemp = 1;  // defined by temperature field

    Teuchos::RCP<DRT::DofSetInterface> dofsettemp =
        Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(numDofsPerNodeTemp, 0, 0, true));
    if (structdis->AddDofSet(dofsettemp) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
  }

  // exchange material pointers for coupled material formulations
  AssignMaterialPointers(structdis, scatradis);

  SetIsInit(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::Setup()
{
  CheckIsInit();

  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::AssignMaterialPointers(
    Teuchos::RCP<DRT::Discretization> structdis, Teuchos::RCP<DRT::Discretization> scatradis)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetMeshDisp(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(disp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetVelocityFields(
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->SetVelocityField(convvel,  // convective vel.
      Teuchos::null,                                // acceleration
      vel,                                          // velocity
      Teuchos::null                                 // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetScalarField(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.SetState(nds, "scalarfield", phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetScalarFieldMicro(
    DRT::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW("transferring micro scalar state to structure discretization not implemented.");
}

FOUR_C_NAMESPACE_CLOSE
