/*---------------------------------------------------------------------*/
/*! \file
\brief Augmented contact interface.

\level 2

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/
#include "contact_augmented_interface.H"
#include "contact_augmented_integrator.H"

#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_coupling3d.H"
#include "../drt_contact/contact_paramsinterface.H"

#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_binarytree.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/epetra_utils.H"

#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_io/io_pstream.H"

#include <Teuchos_Time.hpp>
#include <Epetra_Time.h>
#include <Epetra_IntVector.h>
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::InterfaceDataContainer::InterfaceDataContainer()
    : penBound_(-1.0),
      ct_(-1.0),
      slMaElementAreaRatio_(0),
      isTriangleOnMaster_(false),
      issetup_(false),
      assemble_strat_(INPAR::CONTACT::assemble_none),
      var_type_(INPAR::CONTACT::var_unknown),
      sndofrowmap_(Teuchos::null),
      stdofrowmap_(Teuchos::null)
{
  /* nothing to do */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Interface::Interface(
    const Teuchos::RCP<CONTACT::AUG::InterfaceDataContainer>& interfaceData_ptr)
    : ::CONTACT::CoInterface(interfaceData_ptr),
      interfaceData_ptr_(interfaceData_ptr),
      interfaceData_(*interfaceData_ptr_)
{
  /* do nothing */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Interface::Interface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr, int id,
    const Epetra_Comm& comm, int dim, const Teuchos::ParameterList& icontact, bool selfcontact)
    : CONTACT::CoInterface(interfaceData_ptr, id, comm, dim, icontact, selfcontact),
      interfaceData_ptr_(
          Teuchos::rcp_dynamic_cast<AUG::InterfaceDataContainer>(interfaceData_ptr, true)),
      interfaceData_(*interfaceData_ptr_)
{
  const Teuchos::ParameterList& p_aug = icontact.sublist("AUGMENTED");

  interfaceData_.SetAssembleStratType(
      DRT::INPUT::IntegralValue<INPAR::CONTACT::AssembleStrategy>(p_aug, "ASSEMBLE_STRATEGY"));

  interfaceData_.SetVariationalApproachType(
      DRT::INPUT::IntegralValue<INPAR::CONTACT::VariationalApproach>(
          p_aug, "VARIATIONAL_APPROACH"));

  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::Setup()
{
  if (interfaceData_.IsSetup()) return;

  // find smallest interface element edge length
  double myMinEdgeLength = 1.0e12;

  // find largest slave element area
  double myMaxAreaSl = 0.0;

  // find smallest master element area
  double myMinAreaMa = 1.0e12;

  int myTriangleOnMaster = 0;

  // ------------------------------------------------------------------------
  // loop over all slave elements
  for (int i = 0; i < selerowmap_->NumMyElements(); ++i)
  {
    int gid = selerowmap_->GID(i);
    DRT::Element* ele = idiscret_->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find slave element with gid %i", gid);

    MORTAR::MortarElement* sele = static_cast<MORTAR::MortarElement*>(ele);
    if (myMinEdgeLength > sele->MinEdgeSize()) myMinEdgeLength = sele->MinEdgeSize();

    if (myMaxAreaSl < sele->MoData().Area()) myMaxAreaSl = sele->MoData().Area();
  }

  // ------------------------------------------------------------------------
  // loop over all master elements
  for (int i = 0; i < melerowmap_->NumMyElements(); ++i)
  {
    int gid = melerowmap_->GID(i);
    DRT::Element* ele = idiscret_->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find master element with gid %i", gid);

    MORTAR::MortarElement* mele = static_cast<MORTAR::MortarElement*>(ele);
    if (myMinEdgeLength > mele->MinEdgeSize()) myMinEdgeLength = mele->MinEdgeSize();

    // no Mortar Data container on the master side
    double myarea = mele->ComputeArea();
    if (myMinAreaMa > myarea) myMinAreaMa = myarea;

    switch (ele->Shape())
    {
      case DRT::Element::tri3:
      case DRT::Element::tri6:
        myTriangleOnMaster = 1;
        break;
      default:
        // do nothing
        break;
    }
  }

  // ------------------------------------------------------------------------
  // communicate the minimal edge length and element area over all procs
  double lMins[2] = {myMinEdgeLength, myMinAreaMa};
  double gMins[2] = {1.0e12, 1.0e12};
  Comm().MinAll(&lMins[0], &gMins[0], 2);

  if (gMins[0] == 1.0e12 or gMins[0] < 0.0)
    dserror(
        "ERROR: Global minimal interface edge length couldn't"
        " be calculated! (value: %d)",
        gMins[0]);

  // set penetration bound to 50% of the minimal interface element edge length,
  // if no user defined value is provided.
  if (interfaceData_.PenBound() < 0.0) interfaceData_.PenBound() = 0.5 * gMins[0];

  if (gMins[1] == 1.0e12 or gMins[1] < 0.0)
    dserror(
        "ERROR: Global minimal master element area couldn't"
        " be calculated! (value: %d)",
        gMins[1]);

  const double gMinAreaMa = gMins[1];

  // ------------------------------------------------------------------------
  // communicate the maximal slave element area over all procs
  double lMaxs[2] = {myMaxAreaSl, static_cast<double>(myTriangleOnMaster)};
  double gMaxs[2] = {0.0, 0.0};
  Comm().MaxAll(&lMaxs[0], &gMaxs[0], 2);

  const double gMaxAreaSl = gMaxs[0];
  const bool isTriangleOnMaster = static_cast<bool>(gMaxs[1]);

  const int slMaElementAreaRatio = std::ceil(gMaxAreaSl / gMinAreaMa);

  interfaceData_.SetSlMaElementAreaRatio(slMaElementAreaRatio);
  interfaceData_.SetIsTriangleOnMaster(isTriangleOnMaster);

  interfaceData_.IsSetup() = true;

  // ------------------------------------------------------------------------
  // Setup assemble strategy
  SetupAssembleStrategy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::SetupAssembleStrategy()
{
  const enum INPAR::CONTACT::AssembleStrategy type = interfaceData_.AssembleStratType();
  switch (type)
  {
    case INPAR::CONTACT::assemble_node_based:
    {
      Teuchos::RCP<INTERFACE::AssembleStrategy> assemble_strategy =
          CreateNodeBasedAssembleStrategy();
      interfaceData_.SetAssembleStrategy(assemble_strategy);

      break;
    }
    default:
    {
      dserror("Unknown interface assemble strategy! (enum=%d)", type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::AUG::INTERFACE::AssembleStrategy>
CONTACT::AUG::Interface::CreateNodeBasedAssembleStrategy()
{
  const enum INPAR::CONTACT::VariationalApproach var_type =
      interfaceData_.VariationalApproachType();

  switch (var_type)
  {
    case INPAR::CONTACT::var_complete:
    {
      return Teuchos::rcp(
          new INTERFACE::NodeBasedAssembleStrategy<INTERFACE::CompleteAssemblePolicy>(this));
    }
    case INPAR::CONTACT::var_incomplete:
    {
      return Teuchos::rcp(
          new INTERFACE::NodeBasedAssembleStrategy<INTERFACE::IncompleteAssemblePolicy>(this));
    }
    default:
    {
      dserror("Unknown variational approach! (var_type= \"%s\" | %d)",
          INPAR::CONTACT::VariationalApproach2String(var_type).c_str(), var_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::Initialize()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // call initialization routine of the contact_interface
  CoInterface::Initialize();

  // setup member variables (has to be done only once)
  Setup();

  // *** Reset all augmented quantities *********************************
  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  const int nummynodes = SlaveColNodesBound()->NumMyElements();
  const int* mynodegids = SlaveColNodesBound()->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    DRT::Node* node = Discret().gNode(gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    if (not cnode) dserror("Dynamic cast for node gid %d failed!", gid);

    cnode->InitializeAugDataContainer(
        interfaceData_.SlMaElementAreaRatio(), interfaceData_.IsTriangleOnMaster());

    // reset nodal weighted gap
    cnode->AugData().GetWGap() = 1.0e12;

    // reset nodal scaling factor
    cnode->AugData().GetKappa() = 1.0e12;
    cnode->AugData().GetAugA() = 0.0;

    // reset kappaLin
    cnode->AugData().GetDeriv1st_Kappa().clear();
    cnode->AugData().GetDeriv2nd_Kappa().clear();
    cnode->AugData().GetDeriv1st_A().clear();
    cnode->AugData().GetDeriv2nd_A().clear();

    // reset variables of the complete variational approach
    cnode->AugData().GetDeriv1st_WGapSl().clear();
    cnode->AugData().GetDeriv1st_WGapMa().clear();

    // reset variables of the complete variational approach
    cnode->AugData().GetDeriv1st_WGapSl_Complete().clear();
    cnode->AugData().GetDeriv1st_WGapMa_Complete().clear();

    cnode->AugData().GetDeriv2nd_WGapSl().clear();

    cnode->AugData().GetDeriv2nd_WGapMa().clear();

    cnode->AugData().Get_Debug() = std::pair<int, double>();
    cnode->AugData().Get_DebugVec() = std::vector<std::pair<int, double>>(Dim());

    cnode->AugData().GetDeriv1st_Debug().clear();
    cnode->AugData().GetDeriv1st_DebugVec().clear();

    cnode->AugData().GetDeriv2nd_Debug().clear();
    cnode->AugData().GetDeriv2nd_DebugVec().clear();
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::UpdateMasterSlaveSets()
{
  MORTAR::MortarInterface::UpdateMasterSlaveSets();
  SplitSlaveDofs();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::RedEvaluate(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // interface needs to be complete
  if (!Filled() && Comm().MyPID() == 0)
    dserror("ERROR: FillComplete() not called on interface %", id_);

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  const int nummyelements = selecolmap_->NumMyElements();
  const int* myelementgids = selecolmap_->MyGlobalElements();

  for (int i = 0; i < nummyelements; ++i)
  {
    const int gid1 = myelementgids[i];

    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %", gid1);
    MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(ele1);

    if (selement->ZeroSized()) continue;

    /**************************************************************************
     *    Integrate all remaining quantities only over the slave interface    *
     **************************************************************************/
    // create a Augmented Lagrangian integrator instance with correct NumGP and Dim
    CONTACT::AUG::IntegrationWrapper augIntegrationWrapper(
        InterfaceParams(), selement->Shape(), Comm());
    switch (Dim())
    {
      case 2:
      case 3:
        augIntegrationWrapper.IntegrateDerivSlaveElement((*selement), Comm(), mparams_ptr);
        break;
      default:
        dserror("ERROR: RedEvaluate: Dim value has to be 2 or 3!");
        exit(EXIT_FAILURE);
    }
    /**************************************************************************
     *                       !!! end integration !!!                          *
     **************************************************************************/
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::EvaluateNodalNormals() const
{
  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  const unsigned nummyeles = snoderowmapbound_->NumMyElements();
  const int* mygids = snoderowmapbound_->MyGlobalElements();

  Deriv1stVecMap d_nodal_avg_normal;
  for (unsigned i = 0; i < nummyeles; ++i)
  {
    const int gid = mygids[i];
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::CoNode& cnode = dynamic_cast<CONTACT::CoNode&>(*node);

    /*------------------------------------------------------------------------*/
    std::vector<CONTACT::INTEGRATOR::ElementNormal> adj_ele_normals;
    const double avg_normal_length =
        CONTACT::INTEGRATOR::BuildAveragedNormalAtSlaveNode(adj_ele_normals, cnode);

    /*------------------------------------------------------------------------*/
    CONTACT::INTEGRATOR::Deriv1st_AveragedSlaveNormal(
        cnode, adj_ele_normals, avg_normal_length, d_nodal_avg_normal);

    /*------------------------------------------------------------------------*/
    CONTACT::INTEGRATOR::Deriv2nd_AveragedSlaveNormal(
        cnode, adj_ele_normals, avg_normal_length, d_nodal_avg_normal);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::ExportNodalNormals() const
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  ExportNodalNormalsOnly();

  ExportDeriv1stNodalNormals();

  ExportDeriv2ndNodalNormals();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::ExportNodalNormalsOnly() const
{
  std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>> normals;

  /*--------------------------------------------------------------------------*/
  // (0) Pack normals
  {
    const unsigned my_num_row_entries = snoderowmapbound_->NumMyElements();
    const int* my_rgids = snoderowmapbound_->MyGlobalElements();

    for (unsigned i = 0; i < my_num_row_entries; ++i)
    {
      int gid = my_rgids[i];
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      Teuchos::RCP<Epetra_SerialDenseMatrix>& normals_i = normals[gid];

      // fill nodal matrix
      normals_i = Teuchos::rcp(new Epetra_SerialDenseMatrix(3, 1));
      Epetra_SerialDenseMatrix& n_i = *normals_i;
      n_i(0, 0) = cnode->MoData().n()[0];
      n_i(1, 0) = cnode->MoData().n()[1];
      n_i(2, 0) = cnode->MoData().n()[2];
    }
  }

  /*--------------------------------------------------------------------------*/
  // (1) Export normals
  DRT::Exporter& ex = interfaceData_.Exporter();
  ex.Export(normals);

  /*--------------------------------------------------------------------------*/
  // (2) Unpack normals
  {
    const unsigned my_num_col_entries = snodecolmapbound_->NumMyElements();
    const int* my_cgids = snodecolmapbound_->MyGlobalElements();

    for (unsigned i = 0; i < my_num_col_entries; ++i)
    {
      int gid = my_cgids[i];
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      CoNode& cnode = dynamic_cast<CoNode&>(*node);

      Teuchos::RCP<Epetra_SerialDenseMatrix>& normals_i = normals[gid];

      if (normals_i.is_null()) dserror("GID %d could not be found!", gid);

      // fill nodal matrix
      const Epetra_SerialDenseMatrix& n_i = *normals_i;
      std::copy(n_i.A(), n_i.A() + 3, cnode.MoData().n());
    }
  }

  normals.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::ExportDeriv1stNodalNormals() const
{
  std::map<int, Deriv1stVecMap> export_d_normals;

  /*--------------------------------------------------------------------------*/
  // (0) Pack the 1-st order derivatives of the nodal normals
  {
    const unsigned my_row_num_entries = snoderowmapbound_->NumMyElements();
    const int* my_rgids = snoderowmapbound_->MyGlobalElements();

    for (unsigned i = 0; i < my_row_num_entries; ++i)
    {
      const int gid = my_rgids[i];
      CoNode& cnode = dynamic_cast<CoNode&>(*idiscret_->gNode(gid));

      const Deriv1stVecMap& d_normal = cnode.AugData().GetDeriv1st_N();
      Deriv1stVecMap& exp_d_n = export_d_normals[gid];

      GEN_DATA::copy(d_normal, exp_d_n, GEN_DATA::DeepCopy);
    }
  }

  /*--------------------------------------------------------------------------*/
  // (1) Export the 1-st order derivatives
#ifdef DEBUG
  const double t_start = Teuchos::Time::wallTime();
  const double send_size = CONTACT::INTEGRATOR::SizeOf(export_d_normals);
  std::cout << "Export 1-st order normal derivatives of proc " << Comm().MyPID() << " ("
            << std::scientific << send_size / 1.0e6 << " MByte)...";
#endif

  DRT::Exporter& ex = interfaceData_.Exporter();
  ex.Export(export_d_normals);

#ifdef DEBUG
  const double t_diff = Teuchos::Time::wallTime() - t_start;
  std::cout << " successfully! (" << t_diff << "secs)\n";
#endif

  /*--------------------------------------------------------------------------*/
  // (2) Unpack the 1-st order derivatives
  {
    const unsigned my_col_num_entries = snodecolmapbound_->NumMyElements();
    const int* my_cgids = snodecolmapbound_->MyGlobalElements();

    for (unsigned i = 0; i < my_col_num_entries; ++i)
    {
      const int gid = my_cgids[i];
      CoNode& cnode = dynamic_cast<CoNode&>(*idiscret_->gNode(gid));

      Deriv1stVecMap& d_normal = cnode.AugData().GetDeriv1st_N();
      const Deriv1stVecMap& exp_d_n = export_d_normals[gid];

      GEN_DATA::copy(exp_d_n, d_normal, GEN_DATA::DeepCopy);
    }
  }

  export_d_normals.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::ExportDeriv2ndNodalNormals() const
{
  std::map<int, Deriv2ndVecMap> export_dd_normals;

  /*--------------------------------------------------------------------------*/
  // (0) Pack the 2-nd order derivatives of the nodal normals
  {
    const unsigned my_num_entries = snoderowmapbound_->NumMyElements();
    const int* my_gids = snoderowmapbound_->MyGlobalElements();

    for (unsigned i = 0; i < my_num_entries; ++i)
    {
      const int gid = my_gids[i];
      CoNode& cnode = dynamic_cast<CoNode&>(*idiscret_->gNode(gid));

      const Deriv2ndVecMap& dd_normal = cnode.AugData().GetDeriv2nd_N();
      Deriv2ndVecMap& exp_dd_n = export_dd_normals[gid];

      GEN_DATA::copy(dd_normal, exp_dd_n, GEN_DATA::DeepCopy);
    }
  }

  /*--------------------------------------------------------------------------*/
  // (1) Export the 2-nd order derivatives
#ifdef DEBUG
  const double t_start = Teuchos::Time::wallTime();
  const double send_size = CONTACT::INTEGRATOR::SizeOf(export_dd_normals);
  std::cout << "Export 2-nd order normal derivatives of proc " << Comm().MyPID() << " ("
            << std::scientific << send_size / 1.0e6 << " MByte)...";
#endif

  DRT::Exporter& ex = interfaceData_.Exporter();
  ex.Export(export_dd_normals);

#ifdef DEBUG
  const double t_diff = Teuchos::Time::wallTime() - t_start;
  std::cout << " successfully! (" << t_diff << "secs)\n";
#endif

  /*--------------------------------------------------------------------------*/
  // (2) Unpack the 2-nd order derivatives
  {
    const unsigned my_col_num_entries = snodecolmapbound_->NumMyElements();
    const int* my_cgids = snodecolmapbound_->MyGlobalElements();

    for (unsigned i = 0; i < my_col_num_entries; ++i)
    {
      const int gid = my_cgids[i];
      CoNode& cnode = dynamic_cast<CoNode&>(*idiscret_->gNode(gid));

      Deriv2ndVecMap& dd_normal = cnode.AugData().GetDeriv2nd_N();
      const Deriv2ndVecMap& exp_dd_n = export_dd_normals[gid];

      GEN_DATA::copy(exp_dd_n, dd_normal, GEN_DATA::DeepCopy);
    }
  }

  export_dd_normals.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::EvalActiveContributions(
    const int rriter, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  // search algorithm
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("ERROR: Invalid search algorithm");

  // set active slave node map of this interface and start the
  // coupling evaluation
  cparams_ptr->Set<Epetra_Map>(activenodes_.get(), 1);
  EvaluateCoupling(*interfaceData_.SActiveEleColMap(), NULL, cparams_ptr);
  cparams_ptr->ClearEntry(GEN::AnyDataContainer::DataType::any, 1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleBMatrix(LINALG::SparseMatrix& BMatrix) const
{
  interfaceData_.AssembleStrategy().AssembleBMatrix(BMatrix);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDGLmLinMatrix(LINALG::SparseMatrix& dGLmLinMatrix) const
{
  interfaceData_.AssembleStrategy().AssembleDGLmLinMatrix(dGLmLinMatrix);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDGGLinMatrix(
    LINALG::SparseMatrix& dGGLinMatrix, const Epetra_Vector& cnVec) const
{
  interfaceData_.AssembleStrategy().AssembleDGGLinMatrix(dGGLinMatrix, cnVec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::Add_Var_A_GG(
    Epetra_Vector& sl_force_g, const Epetra_Vector& cnVec) const
{
  interfaceData_.AssembleStrategy().Add_Var_A_GG(sl_force_g, cnVec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::Assemble_SlForceLmInactive(Epetra_Vector& sl_force_lm_inactive,
    const Epetra_Vector& cnVec, const double inactive_scale) const
{
  interfaceData_.AssembleStrategy().Assemble_SlForceLmInactive(
      sl_force_lm_inactive, cnVec, inactive_scale);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleInactiveDDMatrix(LINALG::SparseMatrix& inactive_dd_matrix,
    const Epetra_Vector& cnVec, const double inactive_scale) const
{
  interfaceData_.AssembleStrategy().AssembleInactiveDDMatrix(
      inactive_dd_matrix, cnVec, inactive_scale);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDLmNWGapLinMatrix(
    LINALG::SparseMatrix& dLmNWGapLinMatrix, const enum MapType map_type) const
{
  interfaceData_.AssembleStrategy().AssembleDLmNWGapLinMatrix(dLmNWGapLinMatrix, map_type);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleActiveGapVectors(
    Epetra_Vector& aWGapVec, Epetra_Vector& wGapVec) const
{
  IO::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__ << IO::endl;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  int* myactiven_gids = activen_->MyGlobalElements();
  if (nummynodes != activen_->NumMyElements()) dserror("Dimension mismatch!");

  double* aWGapVec_vals = aWGapVec.Values();
  double* wGapVec_vals = wGapVec.Values();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    if (not cnode) dserror("Dynamic_cast to CoNode failed! [node-gid=%d]", gid);

    // calculate averaged weighted gap
    const double wGap = cnode->AugData().GetWGap();
    const double kappa = cnode->AugData().GetKappa();
    double aWGap = wGap;

    if (kappa == 1.0e12 or aWGap == 1.0e12)
      dserror(
          "ERROR: Kappa and/or the weighted gap should "
          "not be equal 1.0e12 for active nodes!");

    aWGap /= kappa;

    const int rgid = myactiven_gids[i];

    // --- averaged weighted gap vector
    int rlid = aWGapVec.Map().LID(rgid);
    if (rlid == -1) dserror("Sparse vector aWGapVec does not have global row %d", rgid);

    aWGapVec_vals[rlid] += aWGap;

    // --- weighted gap vector
    rlid = wGapVec.Map().LID(rgid);
    if (rlid == -1) dserror("Sparse vector wGapVec does not have global row %d", rgid);

    wGapVec_vals[rlid] += wGap;
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleGapVectorOfAllSlNodes(Epetra_Vector& wGapAllSlNodesVec) const
{
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = snoderowmap_->NumMyElements();
  const int* mynodegids = snoderowmap_->MyGlobalElements();

  int* myslndof_gids = interfaceData_.SNDofRowMap()->MyGlobalElements();
  if (nummynodes != interfaceData_.SNDofRowMap()->NumMyElements()) dserror("Dimension mismatch!");

  double* wGapAllVec_vals = wGapAllSlNodesVec.Values();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    if (not cnode) dserror("Dynamic_cast to CoNode failed! [node-gid=%d]", gid);
    const double wGap = cnode->AugData().GetWGap();

    const int rgid = myslndof_gids[i];

    // --- weighted gap vector
    const int rlid = wGapAllSlNodesVec.Map().LID(rgid);
    if (rlid == -1) dserror("Sparse vector wGapVec does not have global row %d", rgid);

    wGapAllVec_vals[rlid] += wGap;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleLmNVector(Epetra_Vector& lmNVec) const
{
  // get cn
  //  const double cn = InterfaceParams().get<double>("SEMI_SMOOTH_CN");

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  int* myactiven_gids = activen_->MyGlobalElements();
  if (nummynodes != activen_->NumMyElements()) dserror("Dimension mismatch!");

  double* lmNVec_vals = lmNVec.Values();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror(
          "ERROR: AssembleDGLmrhs: Cannot find slave"
          " node with gid %",
          gid);
    CoNode* cnode = static_cast<CoNode*>(node);

    double lmn = cnode->MoData().lm()[0];

    const int rgid = myactiven_gids[i];

    // --- normal lagrange multiplier vector
    int rlid = lmNVec.Map().LID(rgid);
    if (rlid == -1) dserror("Sparse vector lmNVec does not have global row %d", rgid);

    lmNVec_vals[rlid] += lmn;
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugAVector(
    Epetra_Vector& augAVec, Epetra_Vector& kappaVec) const
{
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = snoderowmap_->NumMyElements();
  int* mynodegids = snoderowmap_->MyGlobalElements();

  double* augA_values = augAVec.Values();
  double* kappa_values = kappaVec.Values();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);

    CoNode* cnode = static_cast<CoNode*>(node);

    // *** augmented Area ***
    const double augA = cnode->AugData().GetAugA();

    if (augA > 0.0)
    {
      const int lid = augAVec.Map().LID(gid);
      if (lid == -1) dserror("Sparse vector augAVec does not have global row %d", gid);

      augA_values[lid] += augA;
    }
    else
      dserror(
          "ERROR: The augmented nodal area shouldn't be equal/lower than zero! "
          "(value= %.2e)",
          augA);

    // *** kappa ***
    if (cnode->Active())
    {
      const double kappa = cnode->AugData().GetKappa();

      if (kappa > 0.0)
      {
        const int lid = kappaVec.Map().LID(gid);
        if (lid == -1) dserror("Sparse vector kappaVec does not have global row %d", gid);

        kappa_values[lid] += kappa;
      }
      else
        dserror(
            "ERROR: The weighted area kappa shouldn't be equal/lower than zero! "
            "(value= %.2e)",
            augA);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugInactiveRhs(
    Epetra_Vector& augInactiveRhs, Epetra_Vector& cnVec, const double inactive_scale) const
{
  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes = LINALG::SplitMap(*snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1 / ct_inv;

  const int nummynodes = augInactiveSlaveNodes->NumMyElements();
  const int* mynodegids = augInactiveSlaveNodes->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find inactive slave node with gid %", gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror(
          "ERROR: AugmentedInterface::AssembleInactiverhs: "
          "Node ownership inconsistency!");

    double cn_inv = 1 / cnVec[cnVec.Map().LID(gid)];
    double* lm = cnode->MoData().lm();
    double augA = cnode->AugData().GetAugA();

    std::vector<int> rGid(Dim());
    std::vector<int> rOwner(Dim(), cnode->Owner());
    Epetra_SerialDenseVector rVal(Dim());

    for (int rDof = 0; rDof < cnode->NumDof(); ++rDof)
    {
      rGid[rDof] = cnode->Dofs()[rDof];
      // normal direction
      if (rDof == 0) rVal[rDof] = 2.0 * inactive_scale * cn_inv * lm[rDof] * augA;
      // tangential direction
      else
        rVal[rDof] = ct_inv * lm[rDof] * augA;
    }

    LINALG::Assemble(augInactiveRhs, rVal, rGid, rOwner);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDLmTLmTRhs(Epetra_Vector& dLmTLmTRhs) const
{
  // get ct and invert it
  double ct_inv = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1 / ct_inv;

  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode)
      dserror(
          "ERROR: AssembleDLmTLmTRhs: Cannot find active"
          " slave node with gid %",
          gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // Get the Lagrange multiplier and txi of the current node
    double* lm = cnode->MoData().lm();
    // Get weighted element area
    const double kappa = cnode->AugData().GetKappa();

    std::vector<int> rGid(Dim() - 1);
    std::vector<int> rOwner(Dim() - 1, cnode->Owner());
    Epetra_SerialDenseVector rVal(Dim() - 1);

    /*------------------------------------------------------------------*
     |* 2-D case *******************************************************|
     | (Dim()-1) = 1 and j = 0                                          |
     |                                                                  |
     |        i                      (Dim()-1)*i+j                      |
     |==================================================================|
     |        0             ==>           0                             |
     |        1             ==>           1                             |
     |                       :                                          |
     |                       :                                          |
     |                                                                  |
     |* 3-D case *******************************************************|
     | (Dim()-1) = 2 and j=0,1                                          |
     |                                                                  |
     |        i                      (Dim()-1)*i+j                      |
     |==================================================================|
     |        0             ==>          0,1                            |
     |        1             ==>          2,3                            |
     |                       :                                          |
     |                       :                                          |
     *------------------------------------------------------------------*/
    for (int j = 0; j < (Dim() - 1); ++j)
    {
      rGid[j] = activet_->GID((Dim() - 1) * i + j);
      rVal[j] = ct_inv * lm[j + 1] * kappa;
    }

    // Assemble
    LINALG::Assemble(dLmTLmTRhs, rVal, rGid, rOwner);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDLmTLmTMatrix(LINALG::SparseMatrix& dLmTLmTMatrix) const
{
  // get ct and invert it
  double ct_inv = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1 / ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode)
      dserror(
          "ERROR: AssembleDLmTLmTrhs: Cannot find slave"
          " node with gid %",
          gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDLmTLmTrhs: Node ownership inconsistency!");

    // get weighted element area
    const double kappa = cnode->AugData().GetKappa();

    for (int j = 0; j < (Dim() - 1); ++j)
    {
      int rowId = activet_->GID((Dim() - 1) * i + j);
      int colId = rowId;
      /*-----------------------------------------------------------*
       | Attention:                                                |
       | D_{lm}[1/ct*(dlm_t * lm_t)] = D_{lm}[1/ct*[dlm_t*lm_t]]   |
       *-----------------------------------------------------------*/
      double val = ct_inv * kappa;

      // Assemble
      if (abs(val) > 1.0e-12) dLmTLmTMatrix.Assemble(val, rowId, colId);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleDLmTLmTLinMatrix(LINALG::SparseMatrix& dLmTLmTLinMatrix) const
{
  // get ct and invert it
  double ct_inv = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1 / ct_inv;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDGLmLin: Node ownership inconsistency!");

    double* lm = cnode->MoData().lm();
    Deriv1stMap& augALinMap = cnode->AugData().GetDeriv1st_Kappa();

    /*---------------------------------------------------------*
     | Attention:                                              |
     | D_{d}[-1/ct*(dlm_t * lm_t)] = D_{d}[-1/ct*[dlm_t*lm_t]] |
     *---------------------------------------------------------*/
    for (int j = 0; j < (Dim() - 1); ++j)
    {
      const int rowId = activet_->GID((Dim() - 1) * i + j);

      double tmp = ct_inv * lm[j + 1];

      AssembleMapIntoMatrix(rowId, tmp, augALinMap, dLmTLmTLinMatrix);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugInactiveDiagMatrix(Epetra_Vector& augInactiveDiagMatrix,
    const Epetra_Vector& cnVec, const double inactive_scale) const
{
  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes = LINALG::SplitMap(*snoderowmap_, *activenodes_);

  // get ct and invert it
  double ct_inv = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1 / ct_inv;

  const int nummynodes = augInactiveSlaveNodes->NumMyElements();
  const int* mynodegids = augInactiveSlaveNodes->MyGlobalElements();

  LINALG::SerialDenseVector vals;
  std::vector<int> rowIds(0);
  std::vector<int> rowner(0);

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode)
      dserror(
          "ERROR: AssembleAugInactiveMatrix: Cannot find"
          " inactive slave node with gid %",
          gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror(
          "ERROR: AugmentedInterface::AssembleAugInactiveMatrix:"
          " Node ownership inconsistency!");

    const double cn_inv_scale = 2.0 * inactive_scale / cnVec[cnVec.Map().LID(gid)];
    double augA = cnode->AugData().GetAugA();

    const int numdof = cnode->NumDof();

    if (vals.Length() != numdof)
    {
      vals.Resize(numdof);
      rowIds.resize(numdof, -1);
      rowner.resize(numdof, -1);
    }

    // normal direction
    vals(0) = cn_inv_scale * augA;

    // tangential directions
    std::fill(vals.A() + 1, vals.A() + numdof, ct_inv * augA);

    // copy dof ids
    std::copy(cnode->Dofs(), cnode->Dofs() + numdof, &rowIds[0]);

    // insert owner
    std::fill(&rowner[0], &rowner[0] + numdof, cnode->Owner());

    LINALG::Assemble(augInactiveDiagMatrix, vals, rowIds, rowner);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleAugInactiveLinMatrix(
    LINALG::SparseMatrix& augInactiveLinMatrix, const Epetra_Vector& cnVec,
    const double inactive_scale) const
{
  // get ct and invert it
  double ct_inv = InterfaceParams().get<double>("SEMI_SMOOTH_CT");
  ct_inv = 1 / ct_inv;

  const int nummynodes = inactivenodes_->NumMyElements();
  const int* mynodegids = inactivenodes_->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode)
      dserror(
          "ERROR: AssembleAugInactiveMatrix: Cannot find"
          " inactive slave node with gid %",
          gid);

    if (cnode->Owner() != Comm().MyPID())
      dserror(
          "ERROR: AugmentedInterface::AssembleAugInactiveMatrix:"
          " Node ownership inconsistency!");

    double cn_inv = 1 / cnVec[cnVec.Map().LID(gid)];
    double* lm = cnode->MoData().lm();
    Deriv1stMap& augALinMap = cnode->AugData().GetDeriv1st_A();

    for (int j = 0; j < Dim(); ++j)
    {
      int rowId = cnode->Dofs()[j];
      double tmp = 0.0;
      // normal direction
      if (j == 0) tmp = 2.0 * inactive_scale * cn_inv * lm[j];
      // tangential direction
      else
        tmp = ct_inv * lm[j];

      AssembleMapIntoMatrix(rowId, tmp, augALinMap, augInactiveLinMatrix);
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleContactPotentialTerms(
    const Epetra_Vector& cnVec, double& zn_gn, double& gn_gn, double& zn_zn, double& zt_zt) const
{
  double ct_inv = 1.0 / interfaceData_.Ct();

  // *** Active part *************************************************
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));

    if (cnode->Owner() != Comm().MyPID()) dserror("ERROR: Node ownership inconsistency!");

    const double cn = cnVec[cnVec.Map().LID(gid)];

    const double wgap = cnode->AugData().GetWGap();
    dsassert(wgap != 1.0e12, "ERROR: WGap is equal 1.e12 for a active node!");

    const double kappa = cnode->AugData().GetKappa();
    dsassert(kappa != 1.0e12, "ERROR: Kappa is equal 1.e12 for a active node! (node-id: %d)");

    // get the Lagrange multiplier
    const double* lm = cnode->MoData().lm();

    // *** ACTIVE - NORMAL DIRECTION ***
    if (wgap != 0.0)
    {
      // ** zn_i * awgap_i * A_i **
      zn_gn += wgap * lm[0];
      // ** cn/2 * awgap * awgap * A_i **
      gn_gn += 0.5 * cn * wgap * wgap / kappa;
    }

    // *** ACTIVE - TANGENTIAL DIRECTION ***
    for (int d = 1; d < Dim(); ++d)
    {
      // ** 1/(ct) * zt_i^T*zt_i * A_i **
      zt_zt += ct_inv * lm[d] * lm[d] * kappa;
    }
  }

  // *** Inactive part *************************************************
  Teuchos::RCP<Epetra_Map> augInactiveSlaveNodes = LINALG::SplitMap(*snoderowmap_, *activenodes_);
  // loop over proc's active slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < augInactiveSlaveNodes->NumMyElements(); ++i)
  {
    int gid = augInactiveSlaveNodes->GID(i);
    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));

    if (cnode->Owner() != Comm().MyPID()) dserror("ERROR: Node ownership inconsistency!");

    const int cn_lid = cnVec.Map().LID(gid);
    if (cn_lid < 0) dserror("Couldn't find the cn-LID for GID %d.", gid);
    const double cn_inv = 1.0 / cnVec[cn_lid];
    const double augA = cnode->AugData().GetAugA();

    // get the lagrange multiplier
    const double* lm = cnode->MoData().lm();

    // *** INACTIVE - NORMAL DIRECTION ***
    // ** 1/(cn) * zn_i * zn_i * A_i **
    zn_zn += cn_inv * lm[0] * lm[0] * augA;

    // *** INACTIVE - TANGENTIAL DIRECTION ***
    for (int d = 1; d < Dim(); ++d)
    {
      // ** 1/(ct) * zt_i^T*zt_i * A_i **
      zt_zt += ct_inv * lm[d] * lm[d] * augA;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::Interface::BuildActiveSet(bool init)
{
  // define local variables
  std::vector<int> myactivenodegids(0);
  std::vector<int> myactivedofgids(0);

  // loop over all slave nodes
  const int nummynodes = snoderowmap_->NumMyElements();
  const int* mynodegids = snoderowmap_->MyGlobalElements();

  for (int i = 0; i < nummynodes; ++i)
  {
    const int gid = mynodegids[i];

    CoNode* cnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find node with gid %i", gid);

    const int numdof = cnode->NumDof();

    // -------------------------------------------------------------------
    // RE-BUILDING OF THE ACTIVE SET
    // -------------------------------------------------------------------
    if (not init)
    {
      // check if node is active
      if (cnode->Active())
      {
        myactivenodegids.push_back(cnode->Id());

        for (int j = 0; j < numdof; ++j) myactivedofgids.push_back(cnode->Dofs()[j]);
      }
    }
  }  // end loop over all slave nodes

  // create interface local augmented active node map and augmented active dof map
  activenodes_ = Teuchos::rcp(
      new Epetra_Map(-1, (int)myactivenodegids.size(), &myactivenodegids[0], 0, Comm()));
  activedofs_ =
      Teuchos::rcp(new Epetra_Map(-1, (int)myactivedofgids.size(), &myactivedofgids[0], 0, Comm()));

  inactivenodes_ = LINALG::SplitMap(*snoderowmap_, *activenodes_);
  inactivedofs_ = LINALG::SplitMap(*sdofrowmap_, *activedofs_);

  SplitAugActiveDofs();

  BuildActiveColMaps();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::Interface::SetNodeInitiallyActive(
    const CONTACT::ParamsInterface& cparams, CONTACT::CoNode& cnode) const
{
  static const bool init_contact_by_gap =
      DRT::INPUT::IntegralValue<int>(InterfaceParams(), "INITCONTACTBYGAP");

  const bool node_init_active = cnode.IsInitActive();

  // Either init contact by definition or by gap
  if (node_init_active)
  {
    cnode.Active() = true;
    IO::cout(IO::debug) << "Node #" << std::setw(5) << cnode.Id()
                        << " is set initially active via the condition line.\n";
  }
  else if (init_contact_by_gap)
  {
    SetNodeInitiallyActiveByGap(cnode);
    if (cnode.Active())
      IO::cout(IO::debug) << "Node #" << std::setw(5) << cnode.Id()
                          << " is set initially active by gap.\n";
    else
      IO::cout(IO::debug) << "Node #" << std::setw(5) << cnode.Id()
                          << " is not set initially active by gap, "
                             "due to a too larger positive gap ["
                          << cnode.AugData().GetWGap() << "].\n";
  }

  return cnode.Active();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::SetNodeInitiallyActiveByGap(CONTACT::CoNode& cnode) const
{
  static const double initcontactval = InterfaceParams().get<double>("INITCONTACTGAPVALUE");

  if (cnode.AugData().GetWGap() < initcontactval) cnode.Active() = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AUG::Interface::SplitAugActiveDofs()
{
  // get out of here if augmented active set is empty
  if (activenodes_ == Teuchos::null or activenodes_->NumGlobalElements() == 0)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNGids(activenodes_->NumMyElements());
  std::vector<int> myTGids((Dim() - 1) * activenodes_->NumMyElements());

  // dimension check
  const double dimcheck = (activedofs_->NumGlobalElements()) / (activenodes_->NumGlobalElements());

  if (dimcheck != Dim()) dserror("ERROR: SplitAugActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all augmented active row nodes
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find slave node with gid %", gid);

    CoNode* cnode = static_cast<CoNode*>(node);

    // add first dof to nMap
    myNGids[countN] = cnode->Dofs()[0];
    ++countN;

    // add reamining dofs to tMap
    for (int j = 1; j < cnode->NumDof(); ++j)
    {
      myTGids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNGids.resize(countN);
  myTGids.resize(countT);

  // communicate countN and countT among procs
  int gCountN, gCountT;
  Comm().SumAll(&countN, &gCountN, 1);
  Comm().SumAll(&countT, &gCountT, 1);

  // check global dimensions
  if ((gCountN + gCountT) != activedofs_->NumGlobalElements())
    dserror("ERROR: SplitAugActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = Teuchos::rcp(new Epetra_Map(gCountN, countN, &myNGids[0], 0, Comm()));
  activet_ = Teuchos::rcp(new Epetra_Map(gCountT, countT, &myTGids[0], 0, Comm()));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::SplitSlaveDofs()
{
  // get out of here if augmented active set is empty
  if (snoderowmap_ == Teuchos::null or snoderowmap_->NumGlobalElements() == 0)
  {
    interfaceData_.SNDofRowMap() = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    interfaceData_.STDofRowMap() = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNGids(snoderowmap_->NumMyElements());
  std::vector<int> myTGids((Dim() - 1) * snoderowmap_->NumMyElements());

  // dimension check
  double dimcheck = (sdofrowmap_->NumGlobalElements()) / (snoderowmap_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitSlaveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all augmented active row nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (!cnode) dserror("ERROR: Cannot find slave node with gid %", gid);

    // add first dof to nMap
    myNGids[countN] = cnode->Dofs()[0];
    ++countN;

    // add reamining dofs to tMap
    for (int j = 1; j < cnode->NumDof(); ++j)
    {
      myTGids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNGids.resize(countN);
  myTGids.resize(countT);

  // communicate countN and countT among procs
  int lCount[2] = {countN, countT};
  int gCount[2] = {0, 0};
  Comm().SumAll(&lCount[0], &gCount[0], 2);

  const int gCountN = gCount[0];
  const int gCountT = gCount[1];

  // check global dimensions
  if ((gCountN + gCountT) != sdofrowmap_->NumGlobalElements())
    dserror("ERROR: SplitSlaveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  interfaceData_.SNDofRowMap() =
      Teuchos::rcp(new Epetra_Map(gCountN, countN, &myNGids[0], 0, Comm()));
  interfaceData_.STDofRowMap() =
      Teuchos::rcp(new Epetra_Map(gCountT, countT, &myTGids[0], 0, Comm()));

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::BuildActiveColMaps()
{
  Epetra_IntVector active_row_nodes(*snoderowmap_, true);
  Epetra_IntVector active_col_nodes(*snodecolmap_, true);

  {
    const unsigned my_num_row_entries = activenodes_->NumMyElements();
    const int* active_row_node_gids = activenodes_->MyGlobalElements();

    int* arow_node_vals = active_row_nodes.Values();

    for (unsigned i = 0; i < my_num_row_entries; ++i)
    {
      const int rgid = active_row_node_gids[i];
      const int rlid = active_row_nodes.Map().LID(rgid);
      if (rlid == -1)
        dserror("Active slave node GID #%d is no part of the slave row nodes?!", rgid);

      arow_node_vals[rlid] = 1;
    }
  }

  LINALG::Export(active_row_nodes, active_col_nodes);

  {
    std::vector<int> acol_gids;
    acol_gids.reserve(active_col_nodes.Map().NumMyElements());

    const unsigned my_num_col_entries = snodecolmap_->NumMyElements();
    const int* sl_col_node_gids = snodecolmap_->MyGlobalElements();

    int* acol_node_vals = active_col_nodes.Values();

    for (unsigned i = 0; i < my_num_col_entries; ++i)
    {
      if (acol_node_vals[i] == 1) acol_gids.push_back(sl_col_node_gids[i]);
    }

    interfaceData_.SActiveNodeColMap() = Teuchos::null;
    interfaceData_.SActiveNodeColMap() = Teuchos::rcp(
        new Epetra_Map(-1, (int)acol_gids.size(), acol_gids.data(), 0, idiscret_->Comm()));
  }

  BuildActiveSlaveElementColMap(*interfaceData_.SActiveNodeColMap());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::BuildActiveSlaveElementColMap(const Epetra_Map& sanode_col_map)
{
  const unsigned my_num_col_entries = sanode_col_map.NumMyElements();
  const int* active_col_node_gids = sanode_col_map.MyGlobalElements();

  std::set<int> col_ele_gid_set;

  for (unsigned i = 0; i < my_num_col_entries; ++i)
  {
    const int cngid = active_col_node_gids[i];
    DRT::Node* anode = idiscret_->gNode(cngid);

    DRT::Element** adj_eles = anode->Elements();
    const unsigned num_adj_eles = anode->NumElement();

    for (unsigned e = 0; e < num_adj_eles; ++e)
    {
      const int egid = adj_eles[e]->Id();
      col_ele_gid_set.insert(egid);
    }
  }

  std::vector<int> col_ele_gid_vec(col_ele_gid_set.begin(), col_ele_gid_set.end());
  interfaceData_.SActiveEleColMap() = Teuchos::null;
  interfaceData_.SActiveEleColMap() = Teuchos::rcp(new Epetra_Map(
      -1, (int)col_ele_gid_vec.size(), col_ele_gid_vec.data(), 0, idiscret_->Comm()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CONTACT::AUG::Interface::BuildActiveForceMap(
    const Epetra_Vector& force, const double threshold) const
{
  const int num_my_entries = force.Map().NumMyElements();
  const int* my_entry_gids = force.Map().MyGlobalElements();
  const double* fvalues = force.Values();

  std::vector<int> my_active_gids(0);
  my_active_gids.reserve(num_my_entries);

  for (int i = 0; i < num_my_entries; ++i)
  {
    if (std::abs(fvalues[i]) > threshold) my_active_gids.push_back(my_entry_gids[i]);
  }

  return Teuchos::rcp(new Epetra_Map(
      -1, static_cast<int>(my_active_gids.size()), &my_active_gids[0], 0, force.Comm()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleGradientBMatrixContribution(
    const Epetra_Vector& dincr, const Epetra_Vector& str_grad, Epetra_Vector& lmincr) const
{
  if (dincr.Map().NumMyElements() != str_grad.Map().NumMyElements())
    dserror("The number of local elements does not coincide!");

  const double* dincr_vals = dincr.Values();
  const double* str_grad_vals = str_grad.Values();

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  if (nummynodes != lmincr.Map().NumMyElements())
    dserror(
        "Number of elements does not coincide! Currently we"
        "support only Lagrange Multiplier vectors in normal direction"
        "only!");

  double* lmincr_vals = lmincr.Values();

  for (int lid = 0; lid < nummynodes; ++lid)
  {
    const int gid = mynodegids[lid];

    CoNode* cnode = static_cast<CoNode*>(idiscret_->gNode(gid));
    if (not cnode) dserror("ERROR: Cannot find slave node with gid %", gid);

    const int ndof = cnode->Dofs()[0];
    const int nlid_j = lmincr.Map().LID(ndof);
    if (nlid_j == -1)
    {
      dserror(
          "Couldn't find the normal-dof GID %d in the "
          "Lagrange multiplier increment vector on proc %d!",
          ndof, Comm().MyPID());
    }

    double& lmincr_j = lmincr_vals[nlid_j];

    // --- SLAVE SIDE
    {
      const Deriv2ndMap& varWGapLinSlMap = GetVarWGapLinOfSide<SideType::slave>(*cnode);

      AssembleGradientBMatrixContributionOfSide(
          dincr.Map(), varWGapLinSlMap, 1.0, str_grad_vals, dincr_vals, lmincr_j);
    }

    // --- MASTER SIDE
    {
      const Deriv2ndMap& varWGapLinMaMap = GetVarWGapLinOfSide<SideType::master>(*cnode);

      AssembleGradientBMatrixContributionOfSide(
          dincr.Map(), varWGapLinMaMap, -1.0, str_grad_vals, dincr_vals, lmincr_j);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleGradientBMatrixContributionOfSide(
    const Epetra_BlockMap& gslmadofrowmap, const Deriv2ndMap& varWGapLinSideMap,
    const double scalar, const double* const str_grad_vals, const double* const dincr_vals,
    double& lmincr_j) const
{
  for (Deriv2ndMap::const_iterator ci = varWGapLinSideMap.begin(); ci != varWGapLinSideMap.end();
       ++ci)
  {
    const int var_gid = ci->first;
    const int var_lid = gslmadofrowmap.LID(var_gid);
    if (var_lid == -1)
    {
      //      printf( "Couldn't find the variation dof GID %d in the "
      //          "structural gradient vector on proc %d!\n",
      //          var_gid, Comm()->MyPID() );
      continue;
    }

    const double str_grad_val = str_grad_vals[var_lid] * scalar;

    const Deriv1stMap& grad_varWGapSlMap = ci->second;
    for (Deriv1stMap::const_iterator cj = grad_varWGapSlMap.begin(); cj != grad_varWGapSlMap.end();
         ++cj)
    {
      const int lin_gid = cj->first;
      const int lin_lid = gslmadofrowmap.LID(lin_gid);
      if (lin_lid == -1)
      {
        //        printf( "Couldn't find the linearization dof GID %d in the "
        //            "displacement increment vector on proc %d!",
        //            lin_gid, Comm()->MyPID() );
        continue;
      }

      const double dincr_val = dincr_vals[lin_lid];

      lmincr_j += str_grad_val * dincr_val * cj->second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::AssembleGradientBBMatrixContribution(
    const Epetra_Vector& dincr, const Epetra_Vector& lm, Epetra_Vector& lmincr) const
{
  if (lmincr.Map().NumMyElements() != lm.Map().NumMyElements())
    dserror("The number of local elements does not coincide!");

  const double* dincr_vals = dincr.Values();
  const Epetra_BlockMap& dincr_block_map = dincr.Map();

  const double* lm_vals = lm.Values();
  const Epetra_BlockMap& lm_block_map = lm.Map();

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  const int nummynodes = activenodes_->NumMyElements();
  const int* mynodegids = activenodes_->MyGlobalElements();

  if (nummynodes != lmincr.Map().NumMyElements())
    dserror(
        "Number of elements does not coincide! Currently we"
        "support only Lagrange Multiplier vectors in normal direction!");

  double* lmincr_vals = lmincr.Values();

  for (int k = 0; k < nummynodes; ++k)
  {
    const int kgid = mynodegids[k];

    CoNode* cnode_k = static_cast<CoNode*>(idiscret_->gNode(kgid));
    dsassert(cnode_k, "ERROR: Cannot find slave node!");

    const int ndof_k = cnode_k->Dofs()[0];
    const int nlid_k = lm.Map().LID(ndof_k);
    if (nlid_k == -1)
    {
      dserror(
          "Couldn't find the normal-dof GID %d in the "
          "Lagrange multiplier vector on proc %d!",
          ndof_k, Comm().MyPID());
    }

    // 1-st summand
    const double lk = lm_vals[nlid_k];

    // 2-nd summand
    double& lmincr_k = lmincr_vals[nlid_k];

    // --- SLAVE SIDE
    {
      AssembleGradientBBMatrixContributionOfSide<SideType::slave>(nummynodes, mynodegids,
          dincr_block_map, lm_block_map, *cnode_k, 1.0, lk, dincr_vals, lm_vals, lmincr_vals,
          lmincr_k);
    }

    // --- MASTER SIDE
    {
      AssembleGradientBBMatrixContributionOfSide<SideType::master>(nummynodes, mynodegids,
          dincr_block_map, lm_block_map, *cnode_k, -1.0, lk, dincr_vals, lm_vals, lmincr_vals,
          lmincr_k);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <enum CONTACT::AUG::SideType side>
void CONTACT::AUG::Interface::AssembleGradientBBMatrixContributionOfSide(const int nummynodes,
    const int* const mynodegids, const Epetra_BlockMap& dincr_block_map,
    const Epetra_BlockMap& lm_block_map, const CoNode& cnode_k, const double scalar,
    const double lk, const double* const dincr_vals, const double* const lm_vals,
    double* const lmincr_vals, double& lmincr_k) const
{
  const Deriv1stMap& varWGapSideMap = GetVarWGapOfSide<side>(cnode_k);

  for (Deriv1stMap::const_iterator ci = varWGapSideMap.begin(); ci != varWGapSideMap.end(); ++ci)
  {
    const int var_kgid = ci->first;

    // 1-st summand
    const double bik_lk = ci->second * lk;

    // 2-nd summand
    const double bik = ci->second;

    for (int j = 0; j < nummynodes; ++j)
    {
      const int jgid = mynodegids[j];

      CoNode* cnode_j = static_cast<CoNode*>(idiscret_->gNode(jgid));
      dsassert(cnode_j, "ERROR: Cannot find slave node!");

      const Deriv2ndMap& varWGapLinSideMap = GetVarWGapLinOfSide<side>(*cnode_j);

      Deriv2ndMap::const_iterator cci = varWGapLinSideMap.find(var_kgid);
      if (varWGapLinSideMap.end() == cci) continue;

      const int ndof_j = cnode_j->Dofs()[0];
      const int nlid_j = lm_block_map.LID(ndof_j);
      if (nlid_j == -1)
      {
        dserror(
            "Couldn't find the normal-dof GID %d in the "
            "Lagrange multiplier increment vector on proc %d!",
            ndof_j, Comm().MyPID());
      }

      // 1-st summand
      double& lmincr_j = lmincr_vals[nlid_j];

      // 2-nd summand
      const double bik_lj = bik * lm_vals[nlid_j];

      for (Deriv1stMap::const_iterator cix = cci->second.begin(); cix != cci->second.end(); ++cix)
      {
        const int lin_gid = ci->first;
        const int lin_lid = dincr_block_map.LID(lin_gid);
        if (lin_lid == -1) continue;

        // 1-st summand: lmincr_{j} = B_{ij,x} B_{ik} lm^{k} dincr^{x}
        lmincr_j += cix->second * bik_lk * dincr_vals[lin_lid];

        // 2-nd summand: lmincr_{k} = B_{ij,x} B_{ik} lm^{j} dincr^{x}
        lmincr_k += cix->second * bik_lj * dincr_vals[lin_lid];
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
const CONTACT::AUG::Deriv1stMap&
CONTACT::AUG::Interface::GetVarWGapOfSide<CONTACT::AUG::SideType::master>(const CoNode& cnode) const
{
  return cnode.AugData().GetDeriv1st_WGapMa();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
const CONTACT::AUG::Deriv1stMap&
CONTACT::AUG::Interface::GetVarWGapOfSide<CONTACT::AUG::SideType::slave>(const CoNode& cnode) const
{
  return cnode.AugData().GetDeriv1st_WGapSl();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
const CONTACT::AUG::Deriv2ndMap&
CONTACT::AUG::Interface::GetVarWGapLinOfSide<CONTACT::AUG::SideType::master>(
    const CoNode& cnode) const
{
  return cnode.AugData().GetDeriv2nd_WGapMa();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
const CONTACT::AUG::Deriv2ndMap&
CONTACT::AUG::Interface::GetVarWGapLinOfSide<CONTACT::AUG::SideType::slave>(
    const CoNode& cnode) const
{
  return cnode.AugData().GetDeriv2nd_WGapSl();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class T>
void CONTACT::AUG::AssembleMapIntoMatrix(
    int row, double scal, const T& values, LINALG::SparseMatrix& mat, double threshold)
{
  dsassert(threshold >= 0.0, "The threshold value has to be positive!");

  for (typename T::const_iterator pp = values.begin(); pp != values.end(); ++pp)
  {
    const double val = scal * (pp->second);
    const int col = pp->first;

    if (std::abs(val) < threshold) continue;

    switch (mat.GetMatrixtype())
    {
      case LINALG::SparseMatrix::FE_MATRIX:
        mat.FEAssemble(val, row, col);
        break;
      case LINALG::SparseMatrix::CRS_MATRIX:
        mat.Assemble(val, row, col);
        break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Interface::GetMySquareOfWeightedGapGradientError() const
{
  const int num_my_anodes = activenodes_->NumMyElements();
  const int* my_anode_gids = activenodes_->MyGlobalElements();

  double my_square_error = 0.0;

  for (int i = 0; i < num_my_anodes; ++i)
  {
    const CoNode* cnode = dynamic_cast<const CoNode*>(idiscret_->gNode(my_anode_gids[i]));
    if (not cnode) dserror("Dynamic cast failed!");

    {
      const Deriv1stMap& d_wgap_sl = cnode->AugData().GetDeriv1st_WGapSl();
      const Deriv1stMap& d_wgap_sl_c = cnode->AugData().GetDeriv1st_WGapSl_Complete();

      unsigned sl_sanity_count = 0;

      for (auto& sl_complete : d_wgap_sl_c)
      {
        const int var_dof = sl_complete.first;

        double e = sl_complete.second;
        auto sl_incomplete = d_wgap_sl.find(var_dof);

        if (sl_incomplete != d_wgap_sl.end())
        {
          e -= sl_incomplete->second;
          ++sl_sanity_count;
        }

        std::cout << "GID #" << my_anode_gids[i] << " - sl-error = " << std::sqrt(e * e)
                  << std::endl;

        my_square_error += e * e;
      }

      if (sl_sanity_count != d_wgap_sl.size()) dserror("Size mismatch for the slave side!");
    }

    {
      const Deriv1stMap& d_wgap_ma = cnode->AugData().GetDeriv1st_WGapMa();
      const Deriv1stMap& d_wgap_ma_c = cnode->AugData().GetDeriv1st_WGapMa_Complete();

      unsigned ma_sanity_count = 0;

      for (auto& ma_complete : d_wgap_ma_c)
      {
        const int var_dof = ma_complete.first;

        double e = ma_complete.second;
        auto ma_incomplete = d_wgap_ma.find(var_dof);

        if (ma_incomplete != d_wgap_ma.end())
        {
          e -= ma_incomplete->second;
          ++ma_sanity_count;
        }

        std::cout << "GID #" << my_anode_gids[i] << " - ma-error = " << std::sqrt(e * e)
                  << std::endl;

        my_square_error += e * e;
      }

      if (ma_sanity_count != d_wgap_ma.size()) dserror("Size mismatch for the master side!");
    }
  }

  return my_square_error;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Interface::MyCharacteristicElementLength(
    const enum CONTACT::AUG::SideType stype) const
{
  Teuchos::RCP<const Epetra_Map> ele_map_ptr = interfaceData_.ElementRowMapPtr(stype);
  const int num_my_entries = ele_map_ptr->NumMyElements();
  const int* my_gids = ele_map_ptr->MyGlobalElements();

  double my_max_h_ele = -1.0;
  for (int i = 0; i < num_my_entries; ++i)
  {
    const DRT::Element* ele = idiscret_->gElement(my_gids[i]);
    if (not ele) dserror("Couldn't find the element! (ele-GID=%d)", my_gids[i]);

    const DRT::Node* const* nodes = ele->Nodes();
    switch (ele->Shape())
    {
      case DRT::Element::line2:
      {
        const CoNode& cnode0 = dynamic_cast<const CoNode&>(*nodes[0]);
        const CoNode& cnode1 = dynamic_cast<const CoNode&>(*nodes[1]);

        const LINALG::Matrix<3, 1> X0(cnode0.X(), true);
        LINALG::Matrix<3, 1> diffX(cnode1.X(), false);

        diffX.Update(-1.0, X0, 1.0);

        my_max_h_ele = std::max(my_max_h_ele, diffX.Norm2());

        break;
      }
      default:
        dserror(
            "You have to implement the characteristic element length"
            " calculation for the given element type. (ele-type=%s)",
            DRT::DistypeToString(ele->Shape()).c_str());
        exit(EXIT_FAILURE);
    }
  }

  return my_max_h_ele;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::StoreSeleEvalTimes(const Epetra_Vector& gseleevaltimes)
{
  Teuchos::RCP<Epetra_Vector>& et_selerow = interfaceData_.ElementEvalTimes().sele_row_;
  Teuchos::RCP<Epetra_Vector>& et_selecol = interfaceData_.ElementEvalTimes().sele_col_;

  et_selerow = Teuchos::rcp(new Epetra_Vector(*SlaveRowElements(), true));
  et_selecol = Teuchos::rcp(new Epetra_Vector(*SlaveColElements(), true));

  const int num_gselecol = gseleevaltimes.Map().NumMyElements();
  const int* gselecol_gids = gseleevaltimes.Map().MyGlobalElements();

  for (int i = 0; i < num_gselecol; ++i)
  {
    const int gsele_cgid = gselecol_gids[i];

    const int lsele_clid = SlaveColElements()->LID(gsele_cgid);
    if (lsele_clid == -1) continue;

    (*et_selecol)[lsele_clid] = gseleevaltimes[i];

    const int lsele_rlid = SlaveRowElements()->LID(gsele_cgid);
    if (lsele_rlid == -1) continue;

    (*et_selerow)[lsele_rlid] = gseleevaltimes[i];
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Interface::SplitIntoFarAndCloseSets(std::vector<int>& close_sele,
    std::vector<int>& far_sele, std::vector<int>& local_close_nodes,
    std::vector<int>& local_far_nodes) const
{
  const Teuchos::RCP<const Epetra_Vector>& et_selerow_ptr =
      interfaceData_.ElementEvalTimes().sele_row_;

  if (!et_selerow_ptr.is_null())
  {
    const Epetra_Vector& et_selerow = *et_selerow_ptr;

    const int num_selerow = SlaveRowElements()->NumMyElements();

    int num_zeros = 0;
    int num_non_zeros = 0;

    double max_non_zeros = 0.0;
    double min_non_zeros = std::numeric_limits<double>::max();
    double mean_non_zeros = 0.0;

    for (int i = 0; i < num_selerow; ++i)
    {
      if (et_selerow[i] == 0)
      {
        ++num_zeros;
        continue;
      }

      ++num_non_zeros;
      max_non_zeros = std::max(max_non_zeros, et_selerow[i]);
      min_non_zeros = std::min(min_non_zeros, et_selerow[i]);
      mean_non_zeros += et_selerow[i];
    }
    mean_non_zeros /= (num_non_zeros != 0 ? num_non_zeros : 1.0);

    int mypid = Comm().MyPID();
    std::cout << "(#" << mypid << ") Number of zero evaluate times:     " << num_zeros << "\n";
    std::cout << "(#" << mypid << ") Number of non-zero evaluate times: " << num_non_zeros << "\n";
    std::cout << "(#" << mypid << ") Maximal evaluate time:             " << max_non_zeros << "\n";
    std::cout << "(#" << mypid << ") Minimal non-zero evaluate time:    " << min_non_zeros << "\n";
    std::cout << "(#" << mypid << ") Mean over non-zero evaluate times: " << mean_non_zeros << "\n";

    double gmax_evaltime = 0.0;
    Comm().MaxAll(&max_non_zeros, &gmax_evaltime, 1);

    const double threshold = 0.1 * gmax_evaltime;

    // loop over all row elements to gather the local information
    for (int i = 0; i < SlaveRowElements()->NumMyElements(); ++i)
    {
      // get element
      int gid = SlaveRowElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find element with gid %", gid);

      const bool close = et_selerow[i] > threshold;
      if (close)
      {
        close_sele.push_back(gid);
        for (int k = 0; k < ele->NumNode(); ++k) local_close_nodes.push_back(ele->NodeIds()[k]);
      }
      else
      {
        far_sele.push_back(gid);
        for (int k = 0; k < ele->NumNode(); ++k) local_far_nodes.push_back(ele->NodeIds()[k]);
      }
    }
  }
  else
    CoInterface::SplitIntoFarAndCloseSets(close_sele, far_sele, local_close_nodes, local_far_nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::Interface::CollectRowNodeOwners(
    const DRT::Discretization& structure_dis) const
{
  Teuchos::RCP<Epetra_Map> smnodemap_ptr = LINALG::MergeMap(snoderowmap_, mnoderowmap_);

  Teuchos::RCP<Epetra_Vector> gsmnodeowners =
      Teuchos::rcp(new Epetra_Vector(*smnodemap_ptr, false));
  const int* myrnodegids = smnodemap_ptr->MyGlobalElements();
  const unsigned mynumrnodes = smnodemap_ptr->NumMyElements();

  for (unsigned i = 0; i < mynumrnodes; ++i)
  {
    const DRT::Node* node = idiscret_->gNode(myrnodegids[i]);
    (*gsmnodeowners)[i] = node->Owner();
  }

  return gsmnodeowners;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr(
    const enum SideType stype) const
{
  switch (stype)
  {
    case SideType::slave:
      return ElementRowMapPtr<SideType::slave>();
    case SideType::master:
      return ElementRowMapPtr<SideType::master>();
    case SideType::slave_master:
      return ElementRowMapPtr<SideType::slave_master>();
    default:
      dserror("Unknown sidetype. (enum=%d)", stype);
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr<CONTACT::AUG::SideType::master>() const
{
  if (MEleRowMap().is_null()) dserror("Master element row map ptr is NULL.");

  return MEleRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr<CONTACT::AUG::SideType::slave>() const
{
  if (SEleRowMap().is_null()) dserror("Slave element row map ptr is NULL.");

  return SEleRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <>
Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr<CONTACT::AUG::SideType::slave_master>() const
{
  Teuchos::RCP<const Epetra_Map> slelerowmap = ElementRowMapPtr<SideType::slave>();
  Teuchos::RCP<const Epetra_Map> maelerowmap = ElementRowMapPtr<SideType::master>();

  return LINALG::MergeMap(slelerowmap, maelerowmap);
}

/*----------------------------------------------------------------------------*/
// explicit template instantiations

template Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr<CONTACT::AUG::SideType::slave>() const;
template Teuchos::RCP<const Epetra_Map>
CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr<CONTACT::AUG::SideType::master>() const;
template Teuchos::RCP<const Epetra_Map> CONTACT::AUG::InterfaceDataContainer::ElementRowMapPtr<
    CONTACT::AUG::SideType::slave_master>() const;

template void CONTACT::AUG::AssembleMapIntoMatrix<CONTACT::AUG::Deriv1stMap>(
    int row, double scal, const Deriv1stMap& values, LINALG::SparseMatrix& mat, double threshold);
template void CONTACT::AUG::AssembleMapIntoMatrix<CONTACT::AUG::plain_double_map>(int row,
    double scal, const plain_double_map& values, LINALG::SparseMatrix& mat, double threshold);

template const CONTACT::AUG::Deriv1stMap&
CONTACT::AUG::Interface::GetVarWGapOfSide<CONTACT::AUG::SideType::master>(
    const CoNode& cnode) const;
template const CONTACT::AUG::Deriv1stMap&
CONTACT::AUG::Interface::GetVarWGapOfSide<CONTACT::AUG::SideType::slave>(const CoNode& cnode) const;

template const CONTACT::AUG::Deriv2ndMap&
CONTACT::AUG::Interface::GetVarWGapLinOfSide<CONTACT::AUG::SideType::master>(
    const CoNode& cnode) const;
template const CONTACT::AUG::Deriv2ndMap&
CONTACT::AUG::Interface::GetVarWGapLinOfSide<CONTACT::AUG::SideType::slave>(
    const CoNode& cnode) const;
