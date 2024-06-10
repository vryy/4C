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
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_dofset_definedmapping_wrapper.hpp"
#include "4C_fem_dofset_gidbased_wrapper.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::Init(const int ndim,
    Teuchos::RCP<Core::FE::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  set_is_setup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();
  auto scatradis = scatra_integrator->discretization();
  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<Core::DOFSets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<Core::DOFSets::DofSetInterface> scatradofset = scatradis->GetDofSetProxy();

  // add proxy dofssets of other fields to discretizations and check if number of dofsets is correct
  if (scatradis->AddDofSet(structdofset) != ++scatra_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  scatra_integrator->set_number_of_dof_set_displacement(scatra_dofset_counter);
  scatra_integrator->set_number_of_dof_set_velocity(scatra_dofset_counter);
  if (structdis->AddDofSet(scatradofset) != ++structure_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in structure field");

  if (Global::Problem::Instance()->ELCHControlParams().get<int>("TEMPERATURE_FROM_FUNCT") != -1)
  {
    const int numDofsPerNodeTemp = 1;  // defined by temperature field

    Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsettemp =
        Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(numDofsPerNodeTemp, 0, 0, true));
    if (structdis->AddDofSet(dofsettemp) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
  }

  if (Global::Problem::Instance()->Materials()->FirstIdByType(
          Core::Materials::m_scatra_multiscale) != -1 or
      Global::Problem::Instance()->Materials()->FirstIdByType(
          Core::Materials::m_newman_multiscale) != -1)
  {
    auto dofsetmicro = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(1, 0, 0, true));
    if (scatradis->AddDofSet(dofsetmicro) != ++scatra_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra field");
    scatra_integrator->set_number_of_dof_set_micro_scale(scatra_dofset_counter);
    if (structdis->AddDofSet(dofsetmicro) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
  }

  if (ssi_base->is_s2_i_kinetics_with_pseudo_contact())
  {
    const int numDofsPerNodeStresses = 6;
    Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetstresses = Teuchos::rcp(
        new Core::DOFSets::DofSetPredefinedDoFNumber(numDofsPerNodeStresses, 0, 0, true));
    if (structdis->AddDofSet(dofsetstresses) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
    if (scatradis->AddDofSet(structdis->GetDofSetProxy(structure_dofset_counter)) !=
        ++scatra_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra field");
    scatra_integrator->set_number_of_dof_set_two_tensor_quantity(scatra_dofset_counter);
  }

  assign_material_pointers(structdis, scatradis);

  set_is_init(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::Setup()
{
  check_is_init();

  set_is_setup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::assign_material_pointers(
    Teuchos::RCP<Core::FE::Discretization> structdis,
    Teuchos::RCP<Core::FE::Discretization> scatradis)
{
  const int numelements = scatradis->NumMyColElements();

  for (int i = 0; i < numelements; ++i)
  {
    Core::Elements::Element* scatratele = scatradis->lColElement(i);
    const int gid = scatratele->Id();

    Core::Elements::Element* structele = structdis->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    scatratele->AddMaterial(structele->Material());
    structele->AddMaterial(scatratele->Material());
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::set_mechanical_stress_state(
    Core::FE::Discretization& scatradis, Teuchos::RCP<const Epetra_Vector> stress_state,
    unsigned nds)
{
  scatradis.set_state(nds, "mechanicalStressState", stress_state);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::set_mesh_disp(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(disp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::set_velocity_fields(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->set_velocity_field(convvel,  // convective vel.
      Teuchos::null,                                  // acceleration
      vel,                                            // velocity
      Teuchos::null                                   // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetScalarField(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.set_state(nds, "scalarfield", phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetScalarFieldMicro(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.set_state(nds, "MicroCon", phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolume::SetTemperatureField(
    Core::FE::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp)
{
  structdis.set_state(2, "tempfield", temp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetTemperatureField(
    Core::FE::Discretization& structdis, Teuchos::RCP<const Epetra_Vector> temp)
{
  structdis.set_state(2, "tempfield", temp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::Init(const int ndim,
    Teuchos::RCP<Core::FE::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  set_is_setup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();

  // set pointers
  structdis_ = structdis;
  scatradis_ = scatra_integrator->discretization();

  // set problem dimension
  problem_dimension_ = ndim;

  // first call fill_complete for single discretizations.
  // This way the physical dofs are numbered successively
  structdis_->fill_complete();
  scatradis_->fill_complete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = scatradis_->NumDof(0, scatradis_->lRowNode(0));
  const int ndofperelement_scatra = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
      ndofpernode_scatra, ndofperelement_scatra, 0, true));
  if (structdis->AddDofSet(dofsetaux) != ++structure_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
      ndofpernode_struct, ndofperelement_struct, 0, true));
  if (scatradis_->AddDofSet(dofsetaux) != ++scatra_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  scatra_integrator->set_number_of_dof_set_displacement(scatra_dofset_counter);
  scatra_integrator->set_number_of_dof_set_velocity(scatra_dofset_counter);

  // call assign_degrees_of_freedom also for auxiliary dofsets
  // note: the order of fill_complete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  structdis_->fill_complete(true, false, false);
  scatradis_->fill_complete(true, false, false);

  // setup mortar adapter for surface volume coupling
  adaptermeshtying_ = Teuchos::rcp(new Core::Adapter::CouplingMortar(
      Global::Problem::Instance()->NDim(), Global::Problem::Instance()->mortar_coupling_params(),
      Global::Problem::Instance()->contact_dynamic_params(),
      Global::Problem::Instance()->spatial_approximation_type()));

  set_is_init(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::Setup()
{
  check_is_init();

  std::vector<int> coupleddof(problem_dimension_, 1);
  // Setup of meshtying adapter
  adaptermeshtying_->Setup(structdis_, scatradis_, Teuchos::null, coupleddof, "SSICoupling",
      structdis_->Comm(), Global::Problem::Instance()->FunctionManager(), false, false, 0, 1);

  // extractor for coupled surface of structure discretization with surface scatra
  extractor_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(
      *structdis_->dof_row_map(0), adaptermeshtying_->MasterDofMap(), true));

  set_is_setup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::assign_material_pointers(
    Teuchos::RCP<Core::FE::Discretization> structdis,
    Teuchos::RCP<Core::FE::Discretization> scatradis)
{
  // nothing to do in this case, since
  // transferring scalar state to structure discretization not implemented for
  // transport on structural boundary. Only SolidToScatra coupling available.
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::set_mesh_disp(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(disp)));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::set_velocity_fields(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->set_velocity_field(
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(convvel)),  // convective vel.
      Teuchos::null,                                                             // acceleration
      adaptermeshtying_->MasterToSlave(extractor_->ExtractCondVector(vel)),      // velocity
      Teuchos::null                                                              // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetScalarField(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW(
      "transferring scalar state to structure discretization not implemented for "
      "transport on structural boundary. Only SolidToScatra coupling available.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingBoundary::SetScalarFieldMicro(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW("transferring micro scalar state to structure discretization not implemented.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::Init(const int ndim,
    Teuchos::RCP<Core::FE::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  set_is_setup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();
  auto scatradis = scatra_integrator->discretization();
  // first call fill_complete for single discretizations.
  // This way the physical dofs are numbered successively
  structdis->fill_complete();
  scatradis->fill_complete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_scatra = scatradis->NumDof(0, scatradis->lRowNode(0));
  const int ndofperelement_scatra = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux;
  dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
      ndofpernode_scatra, ndofperelement_scatra, 0, true));
  if (structdis->AddDofSet(dofsetaux) != ++structure_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in structure field");
  dofsetaux = Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
      ndofpernode_struct, ndofperelement_struct, 0, true));
  if (scatradis->AddDofSet(dofsetaux) != ++scatra_dofset_counter)
    FOUR_C_THROW("unexpected dof sets in scatra field");
  scatra_integrator->set_number_of_dof_set_displacement(scatra_dofset_counter);
  scatra_integrator->set_number_of_dof_set_velocity(scatra_dofset_counter);

  // call assign_degrees_of_freedom also for auxiliary dofsets
  // note: the order of fill_complete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. scatra dofs
  // 3. structure auxiliary dofs
  // 4. scatra auxiliary dofs
  structdis->fill_complete(true, false, false);
  scatradis->fill_complete(true, false, false);

  // Scheme: non matching meshes --> volumetric mortar coupling...
  volcoupl_structurescatra_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

  // init projection matrices (use default material strategy)
  volcoupl_structurescatra_->Init(ndim, structdis, scatradis);

  // parallel redistribution is performed in the global control
  // algorithm. We redistribute between Init(...) and Setup().
  // volcoupl_structurescatra_->Redistribute();

  set_is_init(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::Setup()
{
  check_is_init();

  // setup projection matrices (use default material strategy)
  volcoupl_structurescatra_->Setup(Global::Problem::Instance()->VolmortarParams());

  set_is_setup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::assign_material_pointers(
    Teuchos::RCP<Core::FE::Discretization> structdis,
    Teuchos::RCP<Core::FE::Discretization> scatradis)
{
  volcoupl_structurescatra_->AssignMaterials(
      structdis, scatradis, Global::Problem::Instance()->VolmortarParams());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::set_mesh_disp(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(volcoupl_structurescatra_->apply_vector_mapping21(disp));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::set_velocity_fields(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->set_velocity_field(
      volcoupl_structurescatra_->apply_vector_mapping21(convvel),  // convective vel.
      Teuchos::null,                                               // acceleration
      volcoupl_structurescatra_->apply_vector_mapping21(vel),      // velocity
      Teuchos::null                                                // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetScalarField(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.set_state(nds, "scalarfield", volcoupl_structurescatra_->apply_vector_mapping12(phi));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingNonMatchingVolume::SetScalarFieldMicro(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW("transferring micro scalar state to structure discretization not implemented.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::Init(const int ndim,
    Teuchos::RCP<Core::FE::Discretization> structdis, Teuchos::RCP<SSI::SSIBase> ssi_base)
{
  set_is_setup(false);

  int scatra_dofset_counter = 0;
  int structure_dofset_counter = 0;

  auto scatra_integrator = ssi_base->ScaTraField();
  auto scatradis = scatra_integrator->discretization();


  // Note : We need to make sure that the parallel distribution of Volume and Boundary
  //        is the same externally! The best thing is if you do this in your *_dyn.cpp,
  //        i.e., your global control algorithm.

  if (ssi_base->sca_tra_manifold_base_algorithm() == Teuchos::null)
  {
    {
      // get condition which defines the coupling on target discretization
      std::vector<Core::Conditions::Condition*> conds_struct;
      structdis->GetCondition("SSICouplingSolidToScatra", conds_struct);

      // get condition which defines the coupling on source discretization
      std::vector<Core::Conditions::Condition*> conds_scatra;
      scatradis->GetCondition("SSICouplingSolidToScatra", conds_scatra);

      // at least one condition needs to be defined on each discretization
      if (conds_struct.size() == 0 or conds_scatra.size() == 0)
        FOUR_C_THROW(
            "No coupling condition defined on one or both structure or scatra discretization!");

      std::set<int> couplingids;
      for (auto& cond_struct : conds_struct)
        couplingids.insert(cond_struct->parameters().Get<int>("coupling id"));

      Teuchos::RCP<Core::DOFSets::DofSetGIDBasedWrapper> structgidmatchingdofset = Teuchos::rcp(
          new Core::DOFSets::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));

      Teuchos::RCP<Core::DOFSets::DofSetDefinedMappingWrapper> newdofset_scatra =
          Teuchos::rcp(new Core::DOFSets::DofSetDefinedMappingWrapper(
              structgidmatchingdofset, structdis, "SSICouplingSolidToScatra", couplingids));

      // add dofset and check if scatra field has 2 dofsets, so that coupling is possible
      if (scatradis->AddDofSet(newdofset_scatra) != ++scatra_dofset_counter)
        FOUR_C_THROW("unexpected dof sets in scatra field");
      scatra_integrator->set_number_of_dof_set_displacement(scatra_dofset_counter);
      scatra_integrator->set_number_of_dof_set_velocity(scatra_dofset_counter);
    }

    {
      // get condition which defines the coupling on target discretization
      std::vector<Core::Conditions::Condition*> conds_struct;
      structdis->GetCondition("SSICouplingScatraToSolid", conds_struct);

      // get condition which defines the coupling on source discretization
      std::vector<Core::Conditions::Condition*> conds_scatra;
      scatradis->GetCondition("SSICouplingScatraToSolid", conds_scatra);

      // at least one condition needs to be defined on each discretization
      if (conds_struct.size() == 0 or conds_scatra.size() == 0)
        FOUR_C_THROW(
            "No coupling condition defined on one or both structure or scatra discretization!");

      std::set<int> couplingids;
      for (auto& cond_struct : conds_struct)
        couplingids.insert(cond_struct->parameters().Get<int>("coupling id"));

      Teuchos::RCP<Core::DOFSets::DofSetGIDBasedWrapper> scatragidmatchingdofset = Teuchos::rcp(
          new Core::DOFSets::DofSetGIDBasedWrapper(scatradis, scatradis->GetDofSetProxy()));

      for (int couplingid : couplingids)
      {
        std::set<int> tempset;
        tempset.insert(couplingid);

        Teuchos::RCP<Core::DOFSets::DofSetDefinedMappingWrapper> newdofset_struct =
            Teuchos::rcp(new Core::DOFSets::DofSetDefinedMappingWrapper(
                scatragidmatchingdofset, scatradis, "SSICouplingScatraToSolid", tempset));

        structdis->AddDofSet(newdofset_struct);
      }
    }
  }
  else
  {
    int scatra_manifold_dofset_counter(0);

    auto scatra_manifold_integrator = ssi_base->ScaTraManifold();
    auto scatra_manifold_dis = scatra_manifold_integrator->discretization();

    // build a proxy of the structure discretization for the other fields
    auto structdofset = structdis->GetDofSetProxy();
    // build a proxy of the scatra discretization for the other fields
    auto scatradofset = scatradis->GetDofSetProxy();

    // add proxy dofssets of other fields to discretizations and check if number of dofsets is
    // correct
    if (scatradis->AddDofSet(structdofset) != ++scatra_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra field");
    scatra_integrator->set_number_of_dof_set_displacement(scatra_dofset_counter);
    scatra_integrator->set_number_of_dof_set_velocity(scatra_dofset_counter);
    if (structdis->AddDofSet(scatradofset) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");

    // set dummy coupling id, as coupling between scatra_manifold dis and structdis/scatradis should
    // be setup for all conditions
    std::set<int> couplingids;
    couplingids.insert(0);

    auto structgidmatchingdofset = Teuchos::rcp(
        new Core::DOFSets::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));

    auto proxy_structure_scatramanifold =
        Teuchos::rcp(new Core::DOFSets::DofSetDefinedMappingWrapper(
            structgidmatchingdofset, scatra_manifold_dis, "SSISurfaceManifold", couplingids));

    if (scatra_manifold_dis->AddDofSet(proxy_structure_scatramanifold) !=
        ++scatra_manifold_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in scatra manifold field");
    scatra_manifold_integrator->set_number_of_dof_set_displacement(scatra_manifold_dofset_counter);
    scatra_manifold_integrator->set_number_of_dof_set_velocity(scatra_manifold_dofset_counter);
  }

  if (Global::Problem::Instance()->ELCHControlParams().get<int>("TEMPERATURE_FROM_FUNCT") != -1)
  {
    const int numDofsPerNodeTemp = 1;  // defined by temperature field

    Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsettemp =
        Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(numDofsPerNodeTemp, 0, 0, true));
    if (structdis->AddDofSet(dofsettemp) != ++structure_dofset_counter)
      FOUR_C_THROW("unexpected dof sets in structure field");
  }

  // exchange material pointers for coupled material formulations
  assign_material_pointers(structdis, scatradis);

  set_is_init(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::Setup()
{
  check_is_init();

  set_is_setup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::assign_material_pointers(
    Teuchos::RCP<Core::FE::Discretization> structdis,
    Teuchos::RCP<Core::FE::Discretization> scatradis)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::set_mesh_disp(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> disp)
{
  scatra->ScaTraField()->ApplyMeshMovement(disp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::set_velocity_fields(
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra, Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> vel)
{
  scatra->ScaTraField()->set_velocity_field(convvel,  // convective vel.
      Teuchos::null,                                  // acceleration
      vel,                                            // velocity
      Teuchos::null                                   // fsvel
  );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetScalarField(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  dis.set_state(nds, "scalarfield", phi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSICouplingMatchingVolumeAndBoundary::SetScalarFieldMicro(
    Core::FE::Discretization& dis, Teuchos::RCP<const Epetra_Vector> phi, unsigned nds)
{
  FOUR_C_THROW("transferring micro scalar state to structure discretization not implemented.");
}

FOUR_C_NAMESPACE_CLOSE
