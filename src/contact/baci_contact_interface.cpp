/*---------------------------------------------------------------------*/
/*! \file
\brief One contact interface

\level 2


*/
/*---------------------------------------------------------------------*/
#include "baci_contact_interface.hpp"

#include "baci_binstrategy.hpp"
#include "baci_contact_coupling2d.hpp"
#include "baci_contact_coupling3d.hpp"
#include "baci_contact_element.hpp"
#include "baci_contact_friction_node.hpp"
#include "baci_contact_integrator.hpp"
#include "baci_contact_interpolator.hpp"
#include "baci_contact_line_coupling.hpp"
#include "baci_contact_nitsche_utils.hpp"
#include "baci_contact_node.hpp"
#include "baci_contact_selfcontact_binarytree_unbiased.hpp"
#include "baci_io.hpp"
#include "baci_linalg_utils_densematrix_communication.hpp"
#include "baci_linalg_utils_densematrix_multiply.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_mortar_binarytree.hpp"
#include "baci_mortar_defines.hpp"
#include "baci_mortar_dofset.hpp"
#include "baci_mortar_projector.hpp"
#include "baci_rebalance_graph_based.hpp"
#include "baci_scatra_ele_parameter_boundary.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::InterfaceDataContainer::InterfaceDataContainer()
    : selfcontact_(false),
      friction_(false),
      nonSmoothContact_(false),
      two_half_pass_(false),
      constr_direction_(INPAR::CONTACT::constr_vague),
      activenodes_(Teuchos::null),
      activedofs_(Teuchos::null),
      inactivenodes_(Teuchos::null),
      inactivedofs_(Teuchos::null),
      activen_(Teuchos::null),
      activet_(Teuchos::null),
      slipnodes_(Teuchos::null),
      slipdofs_(Teuchos::null),
      slipt_(Teuchos::null),
      nonsmoothnodes_(Teuchos::null),
      smoothnodes_(Teuchos::null),
      sdofVertexRowmap_(Teuchos::null),
      sdofVertexColmap_(Teuchos::null),
      sdofEdgeRowmap_(Teuchos::null),
      sdofEdgeColmap_(Teuchos::null),
      sdofSurfRowmap_(Teuchos::null),
      sdofSurfColmap_(Teuchos::null),
      nextendedghosting_(Teuchos::null),
      eextendedghosting_(Teuchos::null),
      binarytreeself_(Teuchos::null),
      cnValues_(Teuchos::null),
      ctValues_(Teuchos::null),
      smpairs_(0),
      smintpairs_(0),
      intcells_(0)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Interface> CONTACT::Interface::Create(const int id, const Epetra_Comm& comm,
    const int spatialDim, const Teuchos::ParameterList& icontact, const bool selfcontact)
{
  Teuchos::RCP<MORTAR::InterfaceDataContainer> interfaceData_ptr =
      Teuchos::rcp(new CONTACT::InterfaceDataContainer());
  return Teuchos::rcp(
      new Interface(interfaceData_ptr, id, comm, spatialDim, icontact, selfcontact));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Interface::Interface(const Teuchos::RCP<CONTACT::InterfaceDataContainer>& interfaceData)
    : MORTAR::Interface(interfaceData),
      interfaceData_(interfaceData),
      selfcontact_(interfaceData_->IsSelfContact()),
      friction_(interfaceData_->IsFriction()),
      nonSmoothContact_(interfaceData_->IsNonSmoothContact()),
      two_half_pass_(interfaceData_->IsTwoHalfPass()),
      constr_direction_(interfaceData_->ConstraintDirection()),
      activenodes_(interfaceData_->ActiveNodes()),
      activedofs_(interfaceData_->ActiveDofs()),
      inactivenodes_(interfaceData_->InActiveNodes()),
      inactivedofs_(interfaceData_->InActiveDofs()),
      activen_(interfaceData_->ActiveN()),
      activet_(interfaceData_->ActiveT()),
      slipnodes_(interfaceData_->SlipNodes()),
      slipdofs_(interfaceData_->SlipDofs()),
      slipt_(interfaceData_->SlipT()),
      nonsmoothnodes_(interfaceData_->NonSmoothNodes()),
      smoothnodes_(interfaceData_->SmoothNodes()),
      sdofVertexRowmap_(interfaceData_->SdofVertexRowmap()),
      sdofVertexColmap_(interfaceData_->SdofVertexColmap()),
      sdofEdgeRowmap_(interfaceData_->SdofEdgeRowmap()),
      sdofEdgeColmap_(interfaceData_->SdofEdgeColmap()),
      sdofSurfRowmap_(interfaceData_->SdofSurfRowmap()),
      sdofSurfColmap_(interfaceData_->SdofSurfColmap()),
      nextendedghosting_(interfaceData_->NExtendedGhosting()),
      eextendedghosting_(interfaceData_->EExtendedGhosting()),
      binarytreeself_(interfaceData_->BinaryTreeSelf()),
      cnValues_(interfaceData_->CnValues()),
      ctValues_(interfaceData_->CtValues()),
      smpairs_(interfaceData_->SMIntPairs()),
      smintpairs_(interfaceData_->SMIntPairs()),
      intcells_(interfaceData_->IntCells())
{
  /* do nothing */
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Interface::Interface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData,
    const int id, const Epetra_Comm& comm, const int spatialDim,
    const Teuchos::ParameterList& icontact, bool selfcontact)
    : MORTAR::Interface(interfaceData, id, comm, spatialDim, icontact),
      interfaceData_(
          Teuchos::rcp_dynamic_cast<CONTACT::InterfaceDataContainer>(interfaceData, true)),
      selfcontact_(interfaceData_->IsSelfContact()),
      friction_(interfaceData_->IsFriction()),
      nonSmoothContact_(interfaceData_->IsNonSmoothContact()),
      two_half_pass_(interfaceData_->IsTwoHalfPass()),
      constr_direction_(interfaceData_->ConstraintDirection()),
      activenodes_(interfaceData_->ActiveNodes()),
      activedofs_(interfaceData_->ActiveDofs()),
      inactivenodes_(interfaceData_->InActiveNodes()),
      inactivedofs_(interfaceData_->InActiveDofs()),
      activen_(interfaceData_->ActiveN()),
      activet_(interfaceData_->ActiveT()),
      slipnodes_(interfaceData_->SlipNodes()),
      slipdofs_(interfaceData_->SlipDofs()),
      slipt_(interfaceData_->SlipT()),
      nonsmoothnodes_(interfaceData_->NonSmoothNodes()),
      smoothnodes_(interfaceData_->SmoothNodes()),
      sdofVertexRowmap_(interfaceData_->SdofVertexRowmap()),
      sdofVertexColmap_(interfaceData_->SdofVertexColmap()),
      sdofEdgeRowmap_(interfaceData_->SdofEdgeRowmap()),
      sdofEdgeColmap_(interfaceData_->SdofEdgeColmap()),
      sdofSurfRowmap_(interfaceData_->SdofSurfRowmap()),
      sdofSurfColmap_(interfaceData_->SdofSurfColmap()),
      nextendedghosting_(interfaceData_->NExtendedGhosting()),
      eextendedghosting_(interfaceData_->EExtendedGhosting()),
      binarytreeself_(interfaceData_->BinaryTreeSelf()),
      cnValues_(interfaceData_->CnValues()),
      ctValues_(interfaceData_->CtValues()),
      smpairs_(interfaceData_->SMIntPairs()),
      smintpairs_(interfaceData_->SMIntPairs()),
      intcells_(interfaceData_->IntCells())
{
  selfcontact_ = selfcontact;
  nonSmoothContact_ = CORE::UTILS::IntegralValue<int>(icontact, "NONSMOOTH_GEOMETRIES");
  two_half_pass_ = icontact.get<bool>("Two_half_pass");
  constr_direction_ = CORE::UTILS::IntegralValue<INPAR::CONTACT::ConstraintDirection>(
      icontact, "CONSTRAINT_DIRECTIONS");
  smpairs_ = 0;
  smintpairs_ = 0;
  intcells_ = 0;

  // set frictional contact status
  INPAR::CONTACT::FrictionType ftype =
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(icontact, "FRICTION");
  if (ftype != INPAR::CONTACT::friction_none) friction_ = true;

  // set poro contact
  if (icontact.get<int>("PROBTYPE") == INPAR::CONTACT::poroelast ||
      icontact.get<int>("PROBTYPE") == INPAR::CONTACT::poroscatra ||
      icontact.get<int>("PROBTYPE") == INPAR::CONTACT::fpi)
  {
    SetPoroFlag(true);
    SetPoroType(INPAR::MORTAR::poroelast);
  }
  if (icontact.get<int>("PROBTYPE") == INPAR::CONTACT::poroscatra)
    SetPoroType(INPAR::MORTAR::poroscatra);

  // set ehl contact
  if (icontact.get<int>("PROBTYPE") == INPAR::CONTACT::ehl) SetEhlFlag(true);

  // check for redundant slave storage
  // needed for self contact but not wanted for general contact
  // for self contact this is ensured in BuildInterfaces in contact_strategy_factory.cpp
  // so we only print a warning here, as it is possible to have another contact interface with a
  // different ID that does not need to be a self contact interface
  if (!(selfcontact_ or nonSmoothContact_) &&
      interfaceData_->GetExtendGhosting() == INPAR::MORTAR::ExtendGhosting::redundant_all)
    if (Comm().MyPID() == 0)
      std::cout << "\n\nWARNING: We do not want redundant interface storage for contact where not "
                   "needed, as it is very expensive. But we need it e.g. for self contact."
                << std::endl;

  // initialize extended ghosting for RR loop
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                farah 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::UpdateMasterSlaveSets()
{
  // call mortar function
  MORTAR::Interface::UpdateMasterSlaveSets();

  //********************************************************************
  // DOFS
  //********************************************************************
  // do the same business for dofs
  // (get row and column maps of slave and master dofs seperately)
  if (nonSmoothContact_)
  {
    std::vector<int> sVc;  // slave column map
    std::vector<int> sVr;  // slave row map
    std::vector<int> sEc;  // master column map
    std::vector<int> sEr;  // master row map
    std::vector<int> sSc;  // master column map
    std::vector<int> sSr;  // master row map

    for (int i = 0; i < Discret().NodeColMap()->NumMyElements(); ++i)
    {
      int gid = Discret().NodeColMap()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);
      const bool isslave = mrtrnode->IsSlave();
      const int numdof = mrtrnode->NumDof();

      if (isslave)
      {
        // vertex
        if (mrtrnode->IsOnCorner())
        {
          for (int j = 0; j < numdof; ++j) sVc.push_back(mrtrnode->Dofs()[j]);

          if (Discret().NodeRowMap()->MyGID(gid))
            for (int j = 0; j < numdof; ++j) sVr.push_back(mrtrnode->Dofs()[j]);
        }
        // edge
        else if (mrtrnode->IsOnEdge())
        {
          for (int j = 0; j < numdof; ++j) sEc.push_back(mrtrnode->Dofs()[j]);

          if (Discret().NodeRowMap()->MyGID(gid))
            for (int j = 0; j < numdof; ++j) sEr.push_back(mrtrnode->Dofs()[j]);
        }
        // surface
        else if (!mrtrnode->IsOnCornerEdge())
        {
          for (int j = 0; j < numdof; ++j) sSc.push_back(mrtrnode->Dofs()[j]);

          if (Discret().NodeRowMap()->MyGID(gid))
            for (int j = 0; j < numdof; ++j) sSr.push_back(mrtrnode->Dofs()[j]);
        }
        else
        {
          dserror("unknown case!");
        }
      }
    }

    sdofVertexRowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sVr.size(), sVr.data(), 0, Comm()));
    sdofVertexColmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sVc.size(), sVc.data(), 0, Comm()));
    sdofEdgeRowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sEr.size(), sEr.data(), 0, Comm()));
    sdofEdgeColmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sEc.size(), sEc.data(), 0, Comm()));
    sdofSurfRowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sSr.size(), sSr.data(), 0, Comm()));
    sdofSurfColmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sSc.size(), sSc.data(), 0, Comm()));
  }
}

/*----------------------------------------------------------------------*
 |  create and fill cn vector                                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::SetCnCtValues(const int& iter)
{
  // get cn from the input file
  const double cn = InterfaceParams().get<double>("SEMI_SMOOTH_CN");
  const double ct = InterfaceParams().get<double>("SEMI_SMOOTH_CT");

  // set all nodal cn-values to the input value
  GetCn() = CORE::LINALG::CreateVector(*SlaveRowNodes(), true);
  int err = GetCn()->PutScalar(cn);
  if (err != 0) dserror("cn definition failed!");

  // set all nodal ct-values to the input value
  if (friction_)
  {
    GetCt() = CORE::LINALG::CreateVector(*SlaveRowNodes(), true);
    err = GetCt()->PutScalar(ct);
    if (err != 0) dserror("cn definition failed!");
  }

  // modification for edge/corner nodes
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %i", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // calculate characteristic edge length:
    DRT::Element* ele = cnode->Elements()[0];
    Element* cele = dynamic_cast<Element*>(ele);
    std::array<double, 3> pos1 = {0.0, 0.0, 0.0};
    std::array<double, 3> pos2 = {0.0, 0.0, 0.0};
    std::array<double, 3> vec = {0.0, 0.0, 0.0};

    pos1[0] = dynamic_cast<Node*>(cele->Nodes()[0])->X()[0];
    pos1[1] = dynamic_cast<Node*>(cele->Nodes()[0])->X()[1];
    pos1[2] = dynamic_cast<Node*>(cele->Nodes()[0])->X()[2];

    pos2[0] = dynamic_cast<Node*>(cele->Nodes()[1])->X()[0];
    pos2[1] = dynamic_cast<Node*>(cele->Nodes()[1])->X()[1];
    pos2[2] = dynamic_cast<Node*>(cele->Nodes()[1])->X()[2];

    vec[0] = pos1[0] - pos2[0];
    vec[1] = pos1[1] - pos2[1];
    vec[2] = pos1[2] - pos2[2];

    const double length = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (length < 1e-12)
    {
      std::cout << "*** WARNING *** element edge nearly zero" << std::endl;
      continue;
    }

    if (cnode->IsOnEdge())
    {
      GetCnRef()[GetCnRef().Map().LID(cnode->Id())] = cn * (length * length);
      if (friction_) GetCtRef()[GetCtRef().Map().LID(cnode->Id())] = ct * (length * length);
    }

    if (cnode->IsOnCorner())
    {
      GetCnRef()[GetCnRef().Map().LID(cnode->Id())] = cn * (length * length * length * length);
      if (friction_)
        GetCtRef()[GetCtRef().Map().LID(cnode->Id())] = ct * (length * length * length * length);
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::Interface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0) os << "Contact ";
  MORTAR::Interface::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  add contact node (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddNode(Teuchos::RCP<CONTACT::Node> cnode)
{
  idiscret_->AddNode(cnode);
  return;
}

/*----------------------------------------------------------------------*
 |  add contact element (public)                             mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddElement(Teuchos::RCP<CONTACT::Element> cele)
{
  // check for quadratic 2d slave elements to be modified
  if (cele->IsSlave() && (cele->Shape() == CORE::FE::CellType::line3)) quadslave_ = true;

  // check for quadratic 3d slave elements to be modified
  if (cele->IsSlave() &&
      (cele->Shape() == CORE::FE::CellType::quad8 || cele->Shape() == CORE::FE::CellType::tri6))
    quadslave_ = true;

  idiscret_->AddElement(cele);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::UpdateParallelLayoutAndDataStructures(const bool perform_rebalancing,
    const bool enforce_ghosting_update, const int maxdof, const double meanVelocity)
{
  if (perform_rebalancing)
  {
    Redistribute();
    FillCompleteNew(false, maxdof);
  }

  if (perform_rebalancing || enforce_ghosting_update)
  {
    // Assure that at least some maps are available
    if (!Filled()) FillCompleteNew(false, maxdof);

    // Finalize interface maps
    ExtendInterfaceGhostingSafely(meanVelocity);
    FillCompleteNew(true, maxdof);
  }

  // print new parallel distribution
  if (perform_rebalancing)
  {
    if (Comm().MyPID() == 0)
      std::cout << "Interface parallel distribution after rebalancing:" << std::endl;
    PrintParallelDistribution();
  }

  if (perform_rebalancing || enforce_ghosting_update) CreateSearchTree();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::FillCompleteNew(const bool isFinalParallelDistribution, const int maxdof)
{
  std::stringstream ss;
  ss << "CONTACT::Interface::FillCompleteNew of '" << Discret().Name() << "'";
  TEUCHOS_FUNC_TIME_MONITOR(ss.str());

  // store maximum global dof ID handed in
  // this ID is later needed when setting up the Lagrange multiplier
  // dof map, which of course must not overlap with existing dof ranges
  maxdofglobal_ = maxdof;

  /* We'd like to call idiscret_.FillComplete(true,false,false) but this will assign all nodes new
   * degrees of freedom which we don't want. We would like to use the degrees of freedom that were
   * stored in the mortar nodes. To do so, we have to create and set our own version of a
   * MORTAR::DofSet class before we call FillComplete on the interface discretization. The
   * specialized DofSet class will not assign new dofs but will assign the dofs stored in the nodes.
   */
  {
    Teuchos::RCP<MORTAR::DofSet> mrtrdofset = Teuchos::rcp(new MORTAR::DofSet());
    Discret().ReplaceDofSet(mrtrdofset);
  }

  // FillComplete the interface discretization
  Discret().FillComplete(isFinalParallelDistribution, false, false);

  // check whether crosspoints / edge nodes shall be considered or not
  InitializeCrossPoints();

  // check for const/linear interpolation of 2D/3D quadratic Lagrange multipliers
  InitializeLagMultConst();
  InitializeLagMultLin();

  // check/init corner/edge modification
  InitializeCornerEdge();

  // need row and column maps of slave and master nodes / elements / dofs
  // separately so we can easily address them
  UpdateMasterSlaveSets();

  // initialize node and element data container
  InitializeDataContainer();

  // Communicate quadslave status among ALL processors
  CommunicateQuadSlaveStatusAmongAllProcs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ExtendInterfaceGhostingSafely(const double meanVelocity)
{
  using Teuchos::RCP;

  if (Discret().NodeColMap() == nullptr) dserror("NodeColMap not set.");
  if (Discret().ElementColMap() == nullptr) dserror("ElementColMap not set.");

  // later we might export node and element column map to extended or even FULL overlap,
  // thus store the standard column maps first
  {
    // get standard nodal column map (overlap=1)
    oldnodecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().NodeColMap())));

    // get standard element column map (overlap=1)
    oldelecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().ElementColMap())));
  }

  switch (interfaceData_->GetExtendGhosting())
  {
    case INPAR::MORTAR::ExtendGhosting::redundant_all:
    {
      // to ease our search algorithms we'll afford the luxury to ghost all nodes
      // on all processors. To do so, we'll take the node row map and export it to
      // full overlap. Then we export the discretization to full overlap column map.
      // This way, also all mortar elements will be fully ghosted on all processors.

      // we want to do full ghosting on all procs
      std::vector<int> allproc(Comm().NumProc());
      for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

      // fill my own row node ids
      const Epetra_Map* noderowmap = Discret().NodeRowMap();
      std::vector<int> sdata(noderowmap->NumMyElements());
      for (int i = 0; i < noderowmap->NumMyElements(); ++i) sdata[i] = noderowmap->GID(i);

      // gather all gids of nodes redundantly
      std::vector<int> rdata;
      CORE::LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

      // build completely overlapping map of nodes (on ALL processors)
      Teuchos::RCP<Epetra_Map> newnodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
      sdata.clear();
      rdata.clear();

      // fill my own row element ids
      const Epetra_Map* elerowmap = Discret().ElementRowMap();
      sdata.resize(elerowmap->NumMyElements());
      for (int i = 0; i < elerowmap->NumMyElements(); ++i) sdata[i] = elerowmap->GID(i);

      // gather all gids of elements redundantly
      rdata.resize(0);
      CORE::LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

      // build complete overlapping map of elements (on ALL processors)
      Teuchos::RCP<Epetra_Map> newelecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
      sdata.clear();
      rdata.clear();
      allproc.clear();

      // redistribute the discretization of the interface according to the
      // new column layout
      Discret().ExportColumnNodes(*newnodecolmap);
      Discret().ExportColumnElements(*newelecolmap);

      break;
    }
    case INPAR::MORTAR::ExtendGhosting::redundant_master:
    {
      // to ease our search algorithms we'll afford the luxury to ghost all master
      // nodes on all processors. To do so, we'll take the master node row map and
      // export it to full overlap. Then we export the discretization to partially
      // full overlap column map. This way, also all master elements will be fully
      // ghosted on all processors.

      // at least for master, we want to do full ghosting on all procs
      std::vector<int> allproc(Comm().NumProc());
      for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

      // fill my own master row node ids
      const Epetra_Map* noderowmap = Discret().NodeRowMap();
      std::vector<int> sdata;
      for (int i = 0; i < noderowmap->NumMyElements(); ++i)
      {
        int gid = noderowmap->GID(i);
        DRT::Node* node = Discret().gNode(gid);
        if (!node) dserror("Cannot find node with gid %", gid);
        MORTAR::Node* mrtrnode = dynamic_cast<MORTAR::Node*>(node);
        if (!mrtrnode->IsSlave()) sdata.push_back(gid);
      }

      // gather all master row node gids redundantly
      std::vector<int> rdata;
      CORE::LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

      // add my own slave column node ids (non-redundant, standard overlap)
      const Epetra_Map* nodecolmap = Discret().NodeColMap();
      for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
      {
        int gid = nodecolmap->GID(i);
        DRT::Node* node = Discret().gNode(gid);
        if (!node) dserror("Cannot find node with gid %", gid);
        MORTAR::Node* mrtrnode = dynamic_cast<MORTAR::Node*>(node);
        if (mrtrnode->IsSlave()) rdata.push_back(gid);
      }

      // build new node column map (on ALL processors)
      Teuchos::RCP<Epetra_Map> newnodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
      sdata.clear();
      rdata.clear();

      // fill my own master row element ids
      const Epetra_Map* elerowmap = Discret().ElementRowMap();
      sdata.resize(0);
      for (int i = 0; i < elerowmap->NumMyElements(); ++i)
      {
        int gid = elerowmap->GID(i);
        DRT::Element* ele = Discret().gElement(gid);
        if (!ele) dserror("Cannot find element with gid %", gid);
        MORTAR::Element* mrtrele = dynamic_cast<MORTAR::Element*>(ele);
        if (!mrtrele->IsSlave()) sdata.push_back(gid);
      }

      // gather all gids of elements redundantly
      rdata.resize(0);
      CORE::LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

      // add my own slave column node ids (non-redundant, standard overlap)
      const Epetra_Map* elecolmap = Discret().ElementColMap();
      for (int i = 0; i < elecolmap->NumMyElements(); ++i)
      {
        int gid = elecolmap->GID(i);
        DRT::Element* ele = Discret().gElement(gid);
        if (!ele) dserror("Cannot find element with gid %", gid);
        MORTAR::Element* mrtrele = dynamic_cast<MORTAR::Element*>(ele);
        if (mrtrele->IsSlave()) rdata.push_back(gid);
      }

      // build new element column map (on ALL processors)
      Teuchos::RCP<Epetra_Map> newelecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
      sdata.clear();
      rdata.clear();
      allproc.clear();

      // redistribute the discretization of the interface according to the
      // new node / element column layout (i.e. master = full overlap)
      Discret().ExportColumnNodes(*newnodecolmap);
      Discret().ExportColumnElements(*newelecolmap);

      break;
    }
    case INPAR::MORTAR::ExtendGhosting::roundrobin:
    {
      // Nothing to do in case of Round-Robin
      break;
    }
    case INPAR::MORTAR::ExtendGhosting::binning:
    {
      // Extend master column map via binning

      // Create the binning strategy
      RCP<BINSTRATEGY::BinningStrategy> binningstrategy = SetupBinningStrategy(meanVelocity);

      // fill master and slave elements into bins
      std::map<int, std::set<int>> slavebinelemap;
      binningstrategy->DistributeElesToBins(Discret(), slavebinelemap, true);
      std::map<int, std::set<int>> masterbinelemap;
      binningstrategy->DistributeElesToBins(Discret(), masterbinelemap, false);

      // Extend ghosting of the master elements
      std::map<int, std::set<int>> ext_bin_to_ele_map;
      RCP<const Epetra_Map> extendedmastercolmap =
          binningstrategy->ExtendElementColMap(slavebinelemap, masterbinelemap, ext_bin_to_ele_map,
              Teuchos::null, Teuchos::null, Discret().ElementColMap());

      Discret().ExportColumnElements(*extendedmastercolmap);

      // get the node ids of the elements that are to be ghosted and create a proper node column
      // map for their export
      std::set<int> nodes;
      const int numMasterColElements = extendedmastercolmap->NumMyElements();
      for (int lid = 0; lid < numMasterColElements; ++lid)
      {
        DRT::Element* ele = Discret().gElement(extendedmastercolmap->GID(lid));
        const int* nodeids = ele->NodeIds();
        for (int inode = 0; inode < ele->NumNode(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), colnodes.data(), 0, Comm()));

      Discret().ExportColumnNodes(*nodecolmap);
      break;
    }
    default:
    {
      dserror("This case of redundant interface storage has not been covered, yet. Implement it!");
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 |  redistribute contact interface (public)                   popp 08/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Redistribute()
{
  const Teuchos::ParameterList& mortarParallelRedistParams =
      InterfaceParams().sublist("PARALLEL REDISTRIBUTION");

  // make sure we are supposed to be here
  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == INPAR::MORTAR::ParallelRedist::redist_none)
    dserror(
        "You are not supposed to be here since you did not enable PARALLEL_REDIST in the "
        "input file. ");

  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();
  Teuchos::Time time("", true);
  std::set<int>::const_iterator iter;

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  //**********************************************************************
  // (1) SLAVE splitting in close / non-close parts
  //**********************************************************************
  // perform contact search (still with non-optimal distribution)
  Initialize();
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("Invalid search algorithm");

  // split slave element row map and build redundant vector of
  // all close / non-close slave node ids on all procs
  std::vector<int> closeele, noncloseele;
  std::vector<int> localcns, localfns;

  SplitIntoFarAndCloseSets(closeele, noncloseele, localcns, localfns);

  // loop over all elements to reset candidates / search lists
  // (use standard slave column map)
  const int numMySlaveColElements = SlaveColElements()->NumMyElements();
  for (int i = 0; i < numMySlaveColElements; ++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("Cannot find ele with gid %i", gid);
    MORTAR::Element* mele = dynamic_cast<MORTAR::Element*>(ele);

    mele->MoData().SearchElements().resize(0);
  }

  // we need an arbitrary preliminary element row map
  Teuchos::RCP<Epetra_Map> slaveCloseRowEles =
      Teuchos::rcp(new Epetra_Map(-1, (int)closeele.size(), closeele.data(), 0, Comm()));
  Teuchos::RCP<Epetra_Map> slaveNonCloseRowEles =
      Teuchos::rcp(new Epetra_Map(-1, (int)noncloseele.size(), noncloseele.data(), 0, Comm()));
  Teuchos::RCP<Epetra_Map> masterRowEles = Teuchos::rcp(new Epetra_Map(*MasterRowElements()));

  // check for consistency
  if (slaveCloseRowEles->NumGlobalElements() == 0 && slaveNonCloseRowEles->NumGlobalElements() == 0)
    dserror("CONTACT Redistribute: Both slave sets (close/non-close) are empty");

  //**********************************************************************
  // (2) SPECIAL CASES and output to screen
  //**********************************************************************
  // print element overview
  if (!myrank)
  {
    int cl = slaveCloseRowEles->NumGlobalElements();
    int ncl = slaveNonCloseRowEles->NumGlobalElements();
    int ma = masterRowEles->NumGlobalElements();
    std::cout << "Element overview: " << cl << " / " << ncl << " / " << ma
              << "  (close-S / non-close-S / M)";
  }

  // print old parallel distribution
  if (myrank == 0)
    std::cout << "\nInterface parallel distribution before rebalancing:" << std::endl;
  PrintParallelDistribution();

  // use simple base class method if there are ONLY close or non-close elements
  // (return value TRUE, because redistribution performed)
  if (slaveCloseRowEles->NumGlobalElements() == 0 || slaveNonCloseRowEles->NumGlobalElements() == 0)
  {
    MORTAR::Interface::Redistribute();
    return;
  }

  //**********************************************************************
  // (3a) PREPARATIONS decide how many procs are used
  //**********************************************************************
  // first we assume that all procs will be used
  int scproc = numproc;
  int sncproc = numproc;
  int mproc = numproc;

  // minimum number of elements per proc
  const int minele = mortarParallelRedistParams.get<int>("MIN_ELEPROC");

  // Max. relative imbalance between subdomain sizes
  const double imbalance_tol = mortarParallelRedistParams.get<double>("IMBALANCE_TOL");

  // calculate real number of procs to be used
  if (minele > 0)
  {
    scproc = static_cast<int>((slaveCloseRowEles->NumGlobalElements()) / minele);
    sncproc = static_cast<int>((slaveNonCloseRowEles->NumGlobalElements()) / minele);
    mproc = static_cast<int>((masterRowEles->NumGlobalElements()) / minele);
    if (slaveCloseRowEles->NumGlobalElements() < 2 * minele) scproc = 1;
    if (slaveNonCloseRowEles->NumGlobalElements() < 2 * minele) sncproc = 1;
    if (masterRowEles->NumGlobalElements() < 2 * minele) mproc = 1;
    if (scproc > numproc) scproc = numproc;
    if (sncproc > numproc) sncproc = numproc;
    if (mproc > numproc) mproc = numproc;
  }

  // print message
  if (!myrank)
  {
    std::cout << "\nRedistributing interface '" << Discret().Name() << "' .........\n";
    std::cout << "Procs used for redistribution: " << scproc << " / " << sncproc << " / " << mproc
              << " (close-S / non-close-S / M)\n";
  }

  //**********************************************************************
  // (3b) PREPARATIONS build initial node graph
  //**********************************************************************
  // create graph object
  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *SlaveRowNodes(), 108, false));

  // loop over all row nodes to fill graph
  const int numMySlaveRowNodes = SlaveRowNodes()->NumMyElements();
  for (int k = 0; k < numMySlaveRowNodes; ++k)
  {
    int gid = SlaveRowNodes()->GID(k);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);

    // find adjacent elements first
    for (int k = 0; k < node->NumElement(); ++k)
    {
      // store adjacent nodes
      DRT::Element* ele = node->Elements()[k];
      int numnode = ele->NumNode();
      std::vector<int> nodeids(numnode);
      for (int n = 0; n < numnode; ++n) nodeids[n] = ele->NodeIds()[n];

      int err = graph->InsertGlobalIndices(gid, numnode, nodeids.data());
      if (err < 0) dserror("graph->InsertGlobalIndices returned %d", err);
      if (err == 1) dserror("graph->InsertGlobalIndices returned %d", err);
    }
  }

  // fill graph and optimize storage
  graph->FillComplete();
  graph->OptimizeStorage();

  //**********************************************************************
  // (4) CLOSE SLAVE redistribution
  //**********************************************************************

  // build redundant vector of all close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  std::vector<int> globalcns;
  std::set<int> setglobalcns;
  std::vector<int> scnids;
  CORE::LINALG::Gather<int>(localcns, globalcns, numproc, allproc.data(), Comm());
  for (int i = 0; i < (int)globalcns.size(); ++i) setglobalcns.insert(globalcns[i]);
  for (iter = setglobalcns.begin(); iter != setglobalcns.end(); ++iter) scnids.push_back(*iter);

  //**********************************************************************
  // call parallel redistribution
  Teuchos::RCP<const Epetra_CrsGraph> slaveCloseNodeGraph =
      CORE::REBALANCE::BuildGraph(idiscret_, slaveCloseRowEles);

  Teuchos::ParameterList slaveCloseRebalanceParams;
  slaveCloseRebalanceParams.set<std::string>("num parts", std::to_string(scproc));
  slaveCloseRebalanceParams.set<std::string>("imbalance tol", std::to_string(imbalance_tol));

  const auto& [slaveCloseRowNodes, slaveCloseColNodes] =
      CORE::REBALANCE::RebalanceNodeMaps(slaveCloseNodeGraph, slaveCloseRebalanceParams);
  //**********************************************************************

  //**********************************************************************
  // (5) NON-CLOSE SLAVE redistribution
  //**********************************************************************

  // build redundant vector of all non-close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  std::vector<int> globalfns;
  std::set<int> setglobalfns;
  std::vector<int> sncnids;
  CORE::LINALG::Gather<int>(localfns, globalfns, numproc, allproc.data(), Comm());
  for (int i = 0; i < (int)globalfns.size(); ++i) setglobalfns.insert(globalfns[i]);
  for (iter = setglobalfns.begin(); iter != setglobalfns.end(); ++iter) sncnids.push_back(*iter);

  //**********************************************************************
  // call parallel redistribution
  Teuchos::RCP<const Epetra_CrsGraph> slaveNonCloseNodeGraph =
      CORE::REBALANCE::BuildGraph(idiscret_, slaveNonCloseRowEles);

  Teuchos::ParameterList slaveNonCloseRebalanceParams;
  slaveNonCloseRebalanceParams.set<std::string>("num parts", std::to_string(sncproc));
  slaveNonCloseRebalanceParams.set<std::string>("imbalance tol", std::to_string(imbalance_tol));

  const auto& [slaveNonCloseRowNodes, snccolnodes] =
      CORE::REBALANCE::RebalanceNodeMaps(slaveNonCloseNodeGraph, slaveNonCloseRebalanceParams);
  //**********************************************************************

  //**********************************************************************
  // (6) MASTER redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> mrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> mcolnodes = Teuchos::null;

  RedistributeMasterSide(mrownodes, mcolnodes, masterRowEles, comm, mproc, imbalance_tol);

  //**********************************************************************
  // (7) Merge global interface node row and column map
  //**********************************************************************
  // merge slave node row map from close and non-close parts
  Teuchos::RCP<Epetra_Map> srownodes = Teuchos::null;

  //----------------------------------CASE 1: ONE OR BOTH SLAVE SETS EMPTY
  if (slaveCloseRowNodes == Teuchos::null || slaveNonCloseRowNodes == Teuchos::null)
  {
    dserror("CONTACT Redistribute: Both slave sets (close/non-close) are empty");
  }
  //-------------------------------------CASE 2: BOTH SLAVE SETS NON-EMPTY
  else
  {
    // find intersection set of close and non-close nodes
    std::set<int> intersec;
    for (iter = setglobalcns.begin(); iter != setglobalcns.end(); ++iter)
    {
      std::set<int>::const_iterator found = setglobalfns.find(*iter);
      if (found != setglobalfns.end()) intersec.insert(*found);
    }

    // build slave node row map
    const int numMySlaveCloseRowNodes = slaveCloseRowNodes->NumMyElements();
    const int numMySlaveNonCloseRowNodes = slaveNonCloseRowNodes->NumMyElements();
    std::vector<int> mygids(numMySlaveCloseRowNodes + numMySlaveNonCloseRowNodes);
    int count = slaveCloseRowNodes->NumMyElements();

    // first get GIDs of input slaveCloseRowNodes
    for (int i = 0; i < count; ++i) mygids[i] = slaveCloseRowNodes->GID(i);

    // then add GIDs of input slaveNonCloseRowNodes (only new ones)
    for (int i = 0; i < numMySlaveNonCloseRowNodes; ++i)
    {
      // check for intersection gid
      // don't do anything for intersection gids (slaveCloseRowNodes dominates!!!)
      std::set<int>::const_iterator found = intersec.find(slaveNonCloseRowNodes->GID(i));
      if (found != intersec.end()) continue;

      // check for overlap
      if (slaveCloseRowNodes->MyGID(slaveNonCloseRowNodes->GID(i)))
        dserror("CORE::LINALG::MergeMap: Result map is overlapping");

      // add new GIDs to mygids
      mygids[count] = slaveNonCloseRowNodes->GID(i);
      ++count;
    }
    mygids.resize(count);
    sort(mygids.begin(), mygids.end());
    srownodes = Teuchos::rcp(
        new Epetra_Map(-1, (int)mygids.size(), mygids.data(), 0, slaveCloseRowNodes->Comm()));
  }

  // merge interface node row map from slave and master parts
  Teuchos::RCP<Epetra_Map> rownodes = CORE::LINALG::MergeMap(srownodes, mrownodes, false);

  // IMPORTANT NOTE:
  // While merging from the two different slave parts of the discretization
  // (close slave, non-close slave) is feasible for the node row map,
  // this is not possible for the node column map. Some necessary
  // information on ghosting at the transition between close and non-close
  // slave region would always be missed! Thus, we reconstruct a
  // suitable slave node column map "by hand" here. This is quite simply
  // done by exporting the initial node graph to the new distribution
  // and by then asking for its column map.

  // create the output graph (with new slave node row map) and export to it
  Teuchos::RCP<Epetra_CrsGraph> outgraph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *srownodes, 108, false));
  Epetra_Export exporter(graph->RowMap(), *srownodes);
  int err = outgraph->Export(*graph, exporter, Add);
  if (err < 0) dserror("Graph export returned err=%d", err);

  // trash old graph
  graph = Teuchos::null;

  // call fill complete and optimize storage
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  // get column map from the graph -> build slave node column map
  // (do stupid conversion from Epetra_BlockMap to Epetra_Map)
  const Epetra_BlockMap& bcol = outgraph->ColMap();
  Teuchos::RCP<Epetra_Map> scolnodes = Teuchos::rcp(new Epetra_Map(
      bcol.NumGlobalElements(), bcol.NumMyElements(), bcol.MyGlobalElements(), 0, Comm()));

  // trash new graph
  outgraph = Teuchos::null;

  // merge interface node column map from slave and master parts
  Teuchos::RCP<Epetra_Map> colnodes = CORE::LINALG::MergeMap(scolnodes, mcolnodes, false);

  //**********************************************************************
  // (8) Get partitioning information into discretization
  //**********************************************************************
  // build reasonable element maps from the already valid and final node maps
  // (note that nothing is actually redistributed in here)
  const auto& [roweles, coleles] = Discret().BuildElementRowColumn(*rownodes, *colnodes);

  // export nodes and elements to the row map
  Discret().ExportRowNodes(*rownodes);
  Discret().ExportRowElements(*roweles);

  // export nodes and elements to the column map (create ghosting)
  Discret().ExportColumnNodes(*colnodes);
  Discret().ExportColumnElements(*coleles);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::SplitIntoFarAndCloseSets(std::vector<int>& closeele,
    std::vector<int>& noncloseele, std::vector<int>& localcns, std::vector<int>& localfns) const
{
  const bool performSplitting = CORE::UTILS::IntegralValue<bool>(
      InterfaceParams().sublist("PARALLEL REDISTRIBUTION"), "EXPLOIT_PROXIMITY");

  if (performSplitting)
  {
    // loop over all row elements to gather the local information
    for (int i = 0; i < SlaveRowElements()->NumMyElements(); ++i)
    {
      // get element
      int gid = SlaveRowElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("Cannot find element with gid %", gid);
      MORTAR::Element* cele = dynamic_cast<MORTAR::Element*>(ele);

      // store element id and adjacent node ids
      int close = cele->MoData().NumSearchElements();
      if (close > 0)
      {
        closeele.push_back(gid);
        for (int k = 0; k < cele->NumNode(); ++k) localcns.push_back(cele->NodeIds()[k]);
      }
      else
      {
        noncloseele.push_back(gid);
        for (int k = 0; k < cele->NumNode(); ++k) localfns.push_back(cele->NodeIds()[k]);
      }
    }
  }
  else
  {
    // loop over all row elements to gather the local information
    for (int i = 0; i < SlaveRowElements()->NumMyElements(); ++i)
    {
      // get element
      int gid = SlaveRowElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("Cannot find element with gid %", gid);
      MORTAR::Element* cele = dynamic_cast<MORTAR::Element*>(ele);

      // store element id and adjacent node ids
      noncloseele.push_back(gid);
      for (int k = 0; k < cele->NumNode(); ++k) localfns.push_back(cele->NodeIds()[k]);
    }
  }
}

/*----------------------------------------------------------------------*
 | collect distribution data (public)                         popp 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::CollectDistributionData(int& numColElements, int& numRowElements)
{
  // loop over proc's column slave elements of the interface
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* slaveElement = dynamic_cast<Element*>(ele1);

    // bool indicating coupling partners
    bool add = (slaveElement->MoData().NumSearchElements() > 0);

    // Check if this element has any coupling partners.
    // Increment element counter if so.
    if (add) ++numColElements;

    // check if - in addition - the active proc owns this element.
    // Increment input variable rowele if so.
    if (add && slaveElement->Owner() == Comm().MyPID()) ++numRowElements;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::CreateSearchTree()
{
  // warning
#ifdef MORTARGMSHCTN
  if (Dim() == 3 && Comm().MyPID() == 0)
  {
    std::cout << "\n******************************************************************\n";
    std::cout << "GMSH output of all contact tree nodes in 3D needs a lot of memory!\n";
    std::cout << "******************************************************************\n";
  }
#endif

  // binary tree search
  if (SearchAlg() == INPAR::MORTAR::search_binarytree)
  {
    //*****SELF CONTACT*****
    if (SelfContact())
    {
      // set state in interface to intialize all kinds of quantities
      Teuchos::RCP<Epetra_Vector> zero = Teuchos::rcp(new Epetra_Vector(*idiscret_->DofRowMap()));
      SetState(MORTAR::state_new_displacement, *zero);

      // create fully overlapping map of all contact elements
      Teuchos::RCP<Epetra_Map> elefullmap =
          CORE::LINALG::AllreduceEMap(*idiscret_->ElementRowMap());

      // create binary tree object for self contact search
      if (!TwoHalfPass())
      {
        // (NOTE THAT SELF CONTACT SEARCH IS NOT YET FULLY PARALLELIZED!)
        binarytreeself_ = Teuchos::rcp(new CONTACT::SelfBinaryTree(
            Discret(), InterfaceParams(), elefullmap, Dim(), SearchParam()));
      }
      else
      {
        // if we use the two half pass algorithm, we use the unbiased self binary tree
        // implementation
        binarytreeself_ = Teuchos::rcp(new CONTACT::UnbiasedSelfBinaryTree(
            Discret(), InterfaceParams(), elefullmap, Dim(), SearchParam()));
      }
      // initialize the self binary tree
      binarytreeself_->Init();
    }
    //*****TWO BODY CONTACT*****
    else
    {
      Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;
      switch (interfaceData_->GetExtendGhosting())
      {
        case INPAR::MORTAR::ExtendGhosting::roundrobin:
        case INPAR::MORTAR::ExtendGhosting::binning:
        {
          melefullmap = melecolmap_;
          break;
        }
        case INPAR::MORTAR::ExtendGhosting::redundant_all:
        case INPAR::MORTAR::ExtendGhosting::redundant_master:
        {
          melefullmap = CORE::LINALG::AllreduceEMap(*melerowmap_);
          break;
        }
        default:
        {
          dserror("Chosen parallel strategy not supported!");
          break;
        }
      }

      {
        // get update type of binary tree
        INPAR::MORTAR::BinaryTreeUpdateType updatetype =
            CORE::UTILS::IntegralValue<INPAR::MORTAR::BinaryTreeUpdateType>(
                InterfaceParams(), "BINARYTREE_UPDATETYPE");

        // create binary tree object for contact search and setup tree
        binarytree_ = Teuchos::rcp(new MORTAR::BinaryTree(Discret(), selecolmap_, melefullmap,
            Dim(), SearchParam(), updatetype, SearchUseAuxPos()));
        // initialize the binary tree
        binarytree_->Init();
      }
    }
  }

  // no binary tree search
  else
  {
    if (SelfContact()) dserror("Binarytree search needed for self contact");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Initialize Data Container for nodes and elements         farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::InitializeDataContainer()
{
  // call base class functionality
  MORTAR::Interface::InitializeDataContainer();

  // ==================
  // non-smooth contact:
  // we need this master node data container to create an averaged
  // nodal normal field on the master side for the smoothed cpp
  // normal field!
  if (CORE::UTILS::IntegralValue<int>(InterfaceParams(), "CPP_NORMALS") || nonSmoothContact_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = CORE::LINALG::AllreduceEMap(*(MasterRowNodes()));

    for (int i = 0; i < masternodes->NumMyElements(); ++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("Cannot find node with gid %i", gid);
      CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
      mnode->InitializeDataContainer();
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  initialize / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::Initialize()
{
  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
  {
    CONTACT::Node* node = dynamic_cast<CONTACT::Node*>(idiscret_->lColNode(i));

    // reset feasible projection and segmentation status
    node->HasProj() = false;
    node->HasSegment() = false;
  }

  // init normal data in master node data container for cpp calculation
  if (CORE::UTILS::IntegralValue<int>(InterfaceParams(), "CPP_NORMALS"))
  {
    for (int i = 0; i < MasterColNodes()->NumMyElements(); ++i)
    {
      int gid = MasterColNodes()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // reset derivative maps of normal vector
      for (int j = 0; j < (int)((cnode->Data().GetDerivN()).size()); ++j)
        (cnode->Data().GetDerivN())[j].clear();
      (cnode->Data().GetDerivN()).resize(0, 0);

      // reset derivative maps of tangent vectors
      for (int j = 0; j < (int)((cnode->Data().GetDerivTxi()).size()); ++j)
        (cnode->Data().GetDerivTxi())[j].clear();
      (cnode->Data().GetDerivTxi()).resize(0, 0);
      for (int j = 0; j < (int)((cnode->Data().GetDerivTeta()).size()); ++j)
        (cnode->Data().GetDerivTeta())[j].clear();
      (cnode->Data().GetDerivTeta()).resize(0, 0);

      for (int j = 0; j < (int)((cnode->Data().GetDerivTangent()).size()); ++j)
        (cnode->Data().GetDerivTangent())[j].clear();
      (cnode->Data().GetDerivTangent()).resize(0, 0);
    }
  }

  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i = 0; i < SlaveColNodesBound()->NumMyElements(); ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // reset nodal Mortar maps
    // for sts
    cnode->MoData().GetD().clear();
    cnode->MoData().GetM().clear();
    cnode->MoData().GetMmod().clear();
    // for nts
    cnode->MoData().GetDnts().clear();
    cnode->MoData().GetMnts().clear();
    // for lts
    cnode->MoData().GetDlts().clear();
    cnode->MoData().GetMlts().clear();
    // for ltl
    cnode->MoData().GetDltl().clear();
    cnode->MoData().GetMltl().clear();

    // reset derivative maps of normal vector
    for (int j = 0; j < (int)((cnode->Data().GetDerivN()).size()); ++j)
      (cnode->Data().GetDerivN())[j].clear();
    (cnode->Data().GetDerivN()).resize(0, 0);

    // reset derivative maps of tangent vectors
    for (int j = 0; j < (int)((cnode->Data().GetDerivTxi()).size()); ++j)
      (cnode->Data().GetDerivTxi())[j].clear();
    (cnode->Data().GetDerivTxi()).resize(0, 0);
    for (int j = 0; j < (int)((cnode->Data().GetDerivTeta()).size()); ++j)
      (cnode->Data().GetDerivTeta())[j].clear();
    (cnode->Data().GetDerivTeta()).resize(0, 0);

    // reset derivative map of Mortar matrices
    (cnode->Data().GetDerivD()).clear();
    (cnode->Data().GetDerivDlts()).clear();
    (cnode->Data().GetDerivDltl()).clear();
    (cnode->Data().GetDerivM()).clear();
    (cnode->Data().GetDerivMnts()).clear();
    (cnode->Data().GetDerivMlts()).clear();
    (cnode->Data().GetDerivMltl()).clear();

    // reset nodal weighted gap and derivative
    cnode->Data().Getg() = 1.0e12;
    cnode->Data().Getgnts() = 1.0e12;
    cnode->Data().Getglts() = 1.0e12;
    cnode->Data().Getgltl()[0] = 1.0e12;
    cnode->Data().Getgltl()[1] = 1.0e12;
    cnode->Data().Getgltl()[2] = 1.0e12;
    cnode->MoData().GetDscale() = 0.0;
    (cnode->Data().GetDerivG()).clear();
    (cnode->Data().GetDerivGnts()).clear();
    (cnode->Data().GetDerivGlts()).clear();
    for (int j = 0; j < (int)cnode->Data().GetDerivGltl().size(); ++j)
      cnode->Data().GetDerivGltl()[j].clear();
    for (int j = 0; j < (int)cnode->Data().GetDerivJumpltl().size(); ++j)
      cnode->Data().GetDerivJumpltl()[j].clear();
    //    (cnode->Data().GetDerivGltl()).resize(0);

    // reset nodal jump
    cnode->Data().Getjumpltl()[0] = 1.0e12;
    cnode->Data().Getjumpltl()[1] = 1.0e12;
    cnode->Data().Getjumpltl()[2] = 1.0e12;

    // hybrid formulation
    cnode->Data().GetAlphaN() = -1.0;
    cnode->Data().GetAlpha().clear();

    // reset derivative map of lagrange multipliers
    for (int j = 0; j < (int)((cnode->Data().GetDerivZ()).size()); ++j)
      (cnode->Data().GetDerivZ())[j].clear();
    (cnode->Data().GetDerivZ()).resize(0);

    if (friction_)
    {
      FriNode* frinode = dynamic_cast<FriNode*>(cnode);

      // reset SNodes and Mnodes
      frinode->FriData().GetSNodes().clear();
      frinode->FriData().GetMNodes().clear();

      // for gp slip
      if (CORE::UTILS::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
      {
        // reset jump deriv.
        for (int j = 0; j < (int)((frinode->FriData().GetDerivVarJump()).size()); ++j)
          (frinode->FriData().GetDerivVarJump())[j].clear();

        (frinode->FriData().GetDerivVarJump()).resize(2);

        // reset jumps
        frinode->FriData().jump_var()[0] = 0.0;
        frinode->FriData().jump_var()[1] = 0.0;
      }
    }

    // just do poro contact relevant stuff!
    if (InterfaceData().IsPoro())
    {
      cnode->PoroData().GetnCoup() = 0.0;
      cnode->PoroData().GetDerivnCoup().clear();
      cnode->PoroData().GetVelDerivnCoup().clear();
      cnode->PoroData().GetPresDerivnCoup().clear();
    }

    // just do ehl relevant stuff!
    if (ehl_) cnode->EhlData().Clear();
  }

  //**********************************************************************
  // In general, it is sufficient to reset search candidates only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need to reset the search candidates of
  // all slave elements in the fully overlapping column map there. This
  // is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact())
  {
    // loop over all elements to reset candidates / search lists
    // (use fully overlapping column map of S+M elements)
    for (int i = 0; i < idiscret_->NumMyColElements(); ++i)
    {
      DRT::Element* ele = idiscret_->lColElement(i);
      MORTAR::Element* mele = dynamic_cast<MORTAR::Element*>(ele);

      mele->MoData().SearchElements().resize(0);

      // dual shape function coefficient matrix
      mele->MoData().ResetDualShape();
      mele->MoData().ResetDerivDualShape();
    }
  }
  else
  {
    // loop over all elements to reset candidates / search lists
    // (use standard slave column map)
    for (int i = 0; i < SlaveColElements()->NumMyElements(); ++i)
    {
      int gid = SlaveColElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("Cannot find ele with gid %i", gid);
      MORTAR::Element* mele = dynamic_cast<MORTAR::Element*>(ele);

      mele->MoData().SearchElements().resize(0);

      // dual shape function coefficient matrix
      mele->MoData().ResetDualShape();
      mele->MoData().ResetDerivDualShape();
    }
  }

  // clear all Nitsche data
  if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM") ==
      INPAR::MORTAR::algorithm_gpts)
    for (int e = 0; e < Discret().ElementColMap()->NumMyElements(); ++e)
      dynamic_cast<MORTAR::Element*>(Discret().gElement(Discret().ElementColMap()->GID(e)))
          ->GetNitscheContainer()
          .Clear();

  // reset s/m pairs and intcell counters
  smpairs_ = 0;
  smintpairs_ = 0;
  intcells_ = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::SetElementAreas()
{
  //**********************************************************************
  // In general, it is sufficient to compute element areas only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need the element areas of all elements
  // (slave and master) in the fully overlapping column map there. At the
  // same time we initialize the element data containers for self contact.
  // This is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact() or CORE::UTILS::IntegralValue<int>(InterfaceParams(), "CPP_NORMALS") or
      nonSmoothContact_)
  {
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i = 0; i < idiscret_->NumMyColElements(); ++i)
    {
      MORTAR::Element* element = dynamic_cast<MORTAR::Element*>(idiscret_->lColElement(i));
      element->InitializeDataContainer();
      element->MoData().Area() = element->ComputeArea();
    }
  }
  else
  {
    // refer call back to base class version
    MORTAR::Interface::SetElementAreas();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  pre evaluate to calc normals                            farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::PreEvaluate(const int& step, const int& iter)
{
  TEUCHOS_FUNC_TIME_MONITOR("CONTACT::Interface::PreEvaluate");

  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("Invalid search algorithm");

    // TODO: maybe we can remove this debug functionality
#ifdef MORTARGMSHCELLS
  // reset integration cell GMSH files
  int proc = Comm().MyPID();
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";
  FILE* fp = fopen(filename.str().c_str(), "w");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \"Integration Cells Proc " << proc << "\" {" << std::endl;
  fprintf(fp, gmshfilecontent.str().c_str());
  fclose(fp);
#endif  // #ifdef MORTARGMSHCELLS

  // set global vector of cn values
  SetCnCtValues(iter);

  // cpp normals or averaged normal field?
  if (CORE::UTILS::IntegralValue<int>(InterfaceParams(), "CPP_NORMALS"))
  {
    // evaluate cpp nodal normals on slave side
    EvaluateCPPNormals();
  }
  else
  {
    // evaluate averaged nodal normals on slave side
    EvaluateNodalNormals();

    // export nodal normals to slave node column map
    // this call is very expensive and the computation
    // time scales directly with the proc number !
    ExportNodalNormals();
  }

  // set condition specific parameters needed for evaluation
  SetConditionSpecificParameters();

  // compute scaling between coupling types
  //  if(nonSmoothContact_)
  //    ComputeScaling();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  store nts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::StoreNTSvalues()
{
  // create iterators for data types
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all possibly non smooth nodes
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if integration is done
    if (cnode->MoData().GetDnts().size() < 1) continue;

    // if nonsmooth contact is activated and the node is no corner node continue
    // if non-smooth contact is not activated go on
    if (!cnode->IsOnCorner() and nonSmoothContact_) continue;

    //-------------------------------------------------------------------------------------
    // store D matrix entries
    // resize pairedvector to nts size
    if ((int)cnode->MoData().GetD().size() == 0)
      cnode->MoData().GetD().resize(cnode->MoData().GetDnts().size());

    for (CI p = cnode->MoData().GetDnts().begin(); p != cnode->MoData().GetDnts().end(); ++p)
      cnode->MoData().GetD()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store M matrix entries
    for (CImap p = cnode->MoData().GetMnts().begin(); p != cnode->MoData().GetMnts().end(); ++p)
      cnode->MoData().GetM()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store weighted gap
    cnode->Data().Getg() = cnode->Data().Getgnts();

    //-------------------------------------------------------------------------------------
    // store weighted gap linearization
    for (CImap p = cnode->Data().GetDerivGnts().begin(); p != cnode->Data().GetDerivGnts().end();
         ++p)
      cnode->Data().GetDerivG()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store D deriv
    // --> No D linearization!

    //-------------------------------------------------------------------------------------
    // store M deriv
    {
      // Mortar M derivatives
      std::map<int, std::map<int, double>>& mntsderiv = cnode->Data().GetDerivMnts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int, std::map<int, double>>::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int, double>& thismderivnts = cnode->Data().GetDerivMnts()[mgid];
        std::map<int, double>& thismderivmortar = cnode->Data().GetDerivM()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int, double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderivnts.end())
          dserror("StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }
  }  // end node loop


  return;
}


/*----------------------------------------------------------------------*
 |  store lts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::StoreLTSvalues()
{
  // create iterators for data types
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all possibly non smooth nodes
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    double msum = 0.0;
    double ssum = 0.0;

    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if this is an edge or a corner and nonsmooth contact is activated
    if (!cnode->IsOnEdge() and nonSmoothContact_) continue;
    if (cnode->IsOnCorner() and nonSmoothContact_) continue;

    // check if integration is done
    if (cnode->MoData().GetDlts().size() < 1) continue;

    //-------------------------------------------------------------------------------------
    // store D matrix entries
    // resize pairedvector to nts size
    if ((int)cnode->MoData().GetD().size() == 0)
      cnode->MoData().GetD().resize(
          cnode->MoData().GetDlts().size() + cnode->MoData().GetDltl().size());

    for (CI p = cnode->MoData().GetDlts().begin(); p != cnode->MoData().GetDlts().end(); ++p)
    {
      cnode->MoData().GetD()[p->first] += (p->second);
      ssum += (p->second);
    }

    //-------------------------------------------------------------------------------------
    // store M matrix entries
    for (CImap p = cnode->MoData().GetMlts().begin(); p != cnode->MoData().GetMlts().end(); ++p)
    {
      cnode->MoData().GetM()[p->first] += (p->second);
      msum += (p->second);
    }

    //-------------------------------------------------------------------------------------
    // store weighted gap
    cnode->Data().Getg() = cnode->Data().Getglts();

    //-------------------------------------------------------------------------------------
    // store weighted gap linearization
    for (CImap p = cnode->Data().GetDerivGlts().begin(); p != cnode->Data().GetDerivGlts().end();
         ++p)
      cnode->Data().GetDerivG()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store D deriv
    {
      // Mortar M derivatives
      std::map<int, std::map<int, double>>& mntsderiv = cnode->Data().GetDerivDlts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int, std::map<int, double>>::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int, double>& thismderivnts = cnode->Data().GetDerivDlts()[mgid];
        std::map<int, double>& thismderivmortar = cnode->Data().GetDerivD()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int, double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderivnts.end())
          dserror("StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }

    //-------------------------------------------------------------------------------------
    // store M deriv
    {
      // Mortar M derivatives
      std::map<int, std::map<int, double>>& mntsderiv = cnode->Data().GetDerivMlts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int, std::map<int, double>>::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int, double>& thismderivnts = cnode->Data().GetDerivMlts()[mgid];
        std::map<int, double>& thismderivmortar = cnode->Data().GetDerivM()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int, double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderivnts.end())
          dserror("StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }
    //    std::cout << "ssum = " << ssum << "  msum = " << msum << "  balance= " << ssum-msum <<
    //    std::endl;
    if (abs(ssum - msum) > 1e-12) dserror("no slave master balance!");

  }  // end node loop


  return;
}



/*----------------------------------------------------------------------*
 |  store lts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::StoreLTLvalues()
{
  dserror("StoreLTLvalues() is outdated!");
  return;
  //  // create iterators for data types
  //  typedef CORE::GEN::pairedvector<int,double>::const_iterator CI;
  //  typedef std::map<int,double>::const_iterator          CImap;
  //
  //  // loop over all possibly non smooth nodes
  //  for(int i=0; i<SlaveRowNodes()->NumMyElements();++i)
  //  {
  //    int gid = SlaveRowNodes()->GID(i);
  //    DRT::Node* node = idiscret_->gNode(gid);
  //    if (!node)
  //      dserror("Cannot find node with gid %",gid);
  //    Node* cnode = dynamic_cast<Node*>(node);
  //
  //    if(!cnode->IsOnEdge())
  //      continue;
  //    if(cnode->IsOnCorner())
  //      continue;
  //
  //    // check if integration is done
  //    if (cnode->MoData().GetDltl().size()<1)
  //      continue;
  //
  //    //-------------------------------------------------------------------------------------
  //    // store D matrix entries
  //    // resize pairedvector to nts size
  //    if ((int)cnode->MoData().GetD().size()==0)
  //      cnode->MoData().GetD().resize(cnode->MoData().GetDltl().size());
  //
  //    for (CI p=cnode->MoData().GetDltl().begin();p!=cnode->MoData().GetDltl().end();++p)
  //      cnode->MoData().GetD()[p->first] += (p->second);
  //
  //    //-------------------------------------------------------------------------------------
  //    // store M matrix entries
  //    for (CImap p=cnode->MoData().GetMltl().begin();p!=cnode->MoData().GetMltl().end();++p)
  //      cnode->MoData().GetM()[p->first] += (p->second);
  //
  //    //-------------------------------------------------------------------------------------
  //    // store weighted gap
  //    cnode->Data().Getg() = cnode->Data().Getgltl();
  //
  //    //-------------------------------------------------------------------------------------
  //    // store weighted gap linearization
  //    for (CImap
  //    p=cnode->Data().GetDerivGltl().begin();p!=cnode->Data().GetDerivGltl().end();++p)
  //      cnode->Data().GetDerivG()[p->first] += (p->second);
  //
  //    //-------------------------------------------------------------------------------------
  //    // store D deriv
  //    {
  //      // Mortar M derivatives
  //      std::map<int,std::map<int,double> >& mntsderiv = cnode->Data().GetDerivDltl();
  //
  //      // get sizes and iterator start
  //      int mastersize = (int)mntsderiv.size();
  //      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();
  //
  //      /********************************************** LinMMatrix **********/
  //      // loop over all master nodes in the DerivM-map of the current LM slave node
  //      for (int l=0;l<mastersize;++l)
  //      {
  //        int mgid = mntscurr->first;
  //        ++mntscurr;
  //
  //        // Mortar matrix M derivatives
  //        std::map<int,double>&thismderivnts    = cnode->Data().GetDerivDltl()[mgid];
  //        std::map<int,double>&thismderivmortar = cnode->Data().GetDerivD()[mgid];
  //
  //        int mapsize = (int)(thismderivnts.size());
  //
  //        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();
  //
  //        // loop over all directional derivative entries
  //        for (int c=0;c<mapsize;++c)
  //        {
  //          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
  //          ++mcolcurr;
  //        }
  //
  //        // check for completeness of DerivM-Derivatives-iteration
  //        if (mcolcurr!=thismderivnts.end())
  //          dserror("StoreNTS: Not all derivative entries of DerivM considered!");
  //      }
  //    }
  //
  //    //-------------------------------------------------------------------------------------
  //    // store M deriv
  //    {
  //      // Mortar M derivatives
  //      std::map<int,std::map<int,double> >& mntsderiv = cnode->Data().GetDerivMltl();
  //
  //      // get sizes and iterator start
  //      int mastersize = (int)mntsderiv.size();
  //      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();
  //
  //      /********************************************** LinMMatrix **********/
  //      // loop over all master nodes in the DerivM-map of the current LM slave node
  //      for (int l=0;l<mastersize;++l)
  //      {
  //        int mgid = mntscurr->first;
  //        ++mntscurr;
  //
  //        // Mortar matrix M derivatives
  //        std::map<int,double>&thismderivnts    = cnode->Data().GetDerivMltl()[mgid];
  //        std::map<int,double>&thismderivmortar = cnode->Data().GetDerivM()[mgid];
  //
  //        int mapsize = (int)(thismderivnts.size());
  //
  //        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();
  //
  //        // loop over all directional derivative entries
  //        for (int c=0;c<mapsize;++c)
  //        {
  //          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
  //          ++mcolcurr;
  //        }
  //
  //        // check for completeness of DerivM-Derivatives-iteration
  //        if (mcolcurr!=thismderivnts.end())
  //          dserror("StoreNTS: Not all derivative entries of DerivM considered!");
  //      }
  //    }
  //  }// end node loop
  //  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddLTLforcesFric(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");
  const double penaltytan = InterfaceParams().get<double>("PENALTYPARAMTAN");
  const double frcoeff = InterfaceParams().get<double>("FRCOEFF");

  std::array<double, 3> oldtraction = {0.0, 0.0, 0.0};

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    double x = cnode->FriData().tractionoldLTL()[0] * cnode->FriData().tractionoldLTL()[0];
    double y = cnode->FriData().tractionoldLTL()[1] * cnode->FriData().tractionoldLTL()[1];
    double z = cnode->FriData().tractionoldLTL()[2] * cnode->FriData().tractionoldLTL()[2];
    double tracvalue = sqrt(x + y + z);

    if (tracvalue > 1e-8)
    {
      oldtraction[0] = cnode->FriData().tractionoldLTL()[0];
      oldtraction[1] = cnode->FriData().tractionoldLTL()[1];
      oldtraction[2] = cnode->FriData().tractionoldLTL()[2];
      break;
    }
  }

  // maybe the old traction is here zero (first contact step)
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // check if this is active node
    if (cnode->Data().Getgltl()[0] < 1e8 and cnode->Data().Getgltl()[1] < 1e8 and
        cnode->Data().Getgltl()[2] < 1e8)
    {
      // normal force
      std::array<double, 3> fn = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < Dim(); ++dim) fn[dim] = -penalty * cnode->Data().Getgltl()[dim];

      // f trial tangential
      std::array<double, 3> ftrial = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < Dim(); ++dim)
        ftrial[dim] = oldtraction[dim] - penaltytan * cnode->Data().Getjumpltl()[dim];

      // trial norm
      double trialnorm =
          sqrt(ftrial[0] * ftrial[0] + ftrial[1] * ftrial[1] + ftrial[2] * ftrial[2]);

      // maxtrac
      double maxtrac = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);

      // real traction
      std::array<double, 3> ftan = {0.0, 0.0, 0.0};

      if (trialnorm - frcoeff * maxtrac <= 0.0)
      {
        for (int dim = 0; dim < Dim(); ++dim) ftan[dim] = ftrial[dim];
      }
      else
      {
        double coeff = frcoeff * maxtrac / trialnorm;
        for (int dim = 0; dim < Dim(); ++dim) ftan[dim] = coeff * ftrial[dim];
      }

      // store
      cnode->FriData().traction()[0] = ftan[0];
      cnode->FriData().traction()[1] = ftan[1];
      cnode->FriData().traction()[2] = ftan[2];

      // ASSEMBLE
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDltl()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDltl();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value = (p->second) * ftan[dim];
            const int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMltl()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMltl();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value = -(p->second) * ftan[dim];
            const int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }

      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddLTLstiffnessFric(Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");
  const double penaltytan = InterfaceParams().get<double>("PENALTYPARAMTAN");
  const double frcoeff = InterfaceParams().get<double>("FRCOEFF");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  std::array<double, 3> oldtraction = {0.0, 0.0, 0.0};

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    double x = cnode->FriData().tractionoldLTL()[0] * cnode->FriData().tractionoldLTL()[0];
    double y = cnode->FriData().tractionoldLTL()[1] * cnode->FriData().tractionoldLTL()[1];
    double z = cnode->FriData().tractionoldLTL()[2] * cnode->FriData().tractionoldLTL()[2];
    double tracvalue = sqrt(x + y + z);

    if (tracvalue > 1e-8)
    {
      oldtraction[0] = cnode->FriData().tractionoldLTL()[0];
      oldtraction[1] = cnode->FriData().tractionoldLTL()[1];
      oldtraction[2] = cnode->FriData().tractionoldLTL()[2];
      break;
    }
  }

  // maybe the old traction is here zero (first contact step)
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // check if this is active node
    if (cnode->Data().Getgltl()[0] < 1e8 and cnode->Data().Getgltl()[1] < 1e8 and
        cnode->Data().Getgltl()[2] < 1e8)
    {
      // state
      bool stick = true;
      CORE::GEN::pairedvector<int, double> coefflin(100);

      // normal force
      std::array<double, 3> fn = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < Dim(); ++dim) fn[dim] = -penalty * cnode->Data().Getgltl()[dim];

      // f trial tangential
      std::array<double, 3> ftrial = {0.0, 0.0, 0.0};
      for (int dim = 0; dim < Dim(); ++dim)
        ftrial[dim] = oldtraction[dim] - penaltytan * cnode->Data().Getjumpltl()[dim];

      double coeff = 0.0;

      // trial norm
      double trialnorm =
          sqrt(ftrial[0] * ftrial[0] + ftrial[1] * ftrial[1] + ftrial[2] * ftrial[2]);

      // maxtrac
      double maxtrac = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);

      // real traction
      std::array<double, 3> ftan = {0.0, 0.0, 0.0};

      if (trialnorm - frcoeff * maxtrac <= 0.0)
      {
        stick = true;
        for (int dim = 0; dim < Dim(); ++dim) ftan[dim] = ftrial[dim];
      }
      else
      {
        stick = false;
        coeff = frcoeff * maxtrac / trialnorm;
        for (int dim = 0; dim < Dim(); ++dim) ftan[dim] = coeff * ftrial[dim];

        CORE::GEN::pairedvector<int, double> fn_x(100);
        CORE::GEN::pairedvector<int, double> fn_y(100);
        CORE::GEN::pairedvector<int, double> fn_z(100);

        CORE::GEN::pairedvector<int, double> ft_x(100);
        CORE::GEN::pairedvector<int, double> ft_y(100);
        CORE::GEN::pairedvector<int, double> ft_z(100);

        for (CImap pp = cnode->Data().GetDerivGltl()[0].begin();
             pp != cnode->Data().GetDerivGltl()[0].end(); ++pp)
          fn_x[pp->first] -= penalty * (pp->second);
        for (CImap pp = cnode->Data().GetDerivGltl()[1].begin();
             pp != cnode->Data().GetDerivGltl()[1].end(); ++pp)
          fn_y[pp->first] -= penalty * (pp->second);
        for (CImap pp = cnode->Data().GetDerivGltl()[2].begin();
             pp != cnode->Data().GetDerivGltl()[2].end(); ++pp)
          fn_z[pp->first] -= penalty * (pp->second);

        for (CImap pp = cnode->Data().GetDerivJumpltl()[0].begin();
             pp != cnode->Data().GetDerivJumpltl()[0].end(); ++pp)
          ft_x[pp->first] -= penaltytan * (pp->second);
        for (CImap pp = cnode->Data().GetDerivJumpltl()[1].begin();
             pp != cnode->Data().GetDerivJumpltl()[1].end(); ++pp)
          ft_y[pp->first] -= penaltytan * (pp->second);
        for (CImap pp = cnode->Data().GetDerivJumpltl()[2].begin();
             pp != cnode->Data().GetDerivJumpltl()[2].end(); ++pp)
          ft_z[pp->first] -= penaltytan * (pp->second);

        CORE::GEN::pairedvector<int, double> maxtraclin(100);
        for (CI pp = fn_x.begin(); pp != fn_x.end(); ++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0 / maxtrac) * (pp->second) * 2.0 * fn[0] * pp->second;
        for (CI pp = fn_y.begin(); pp != fn_y.end(); ++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0 / maxtrac) * (pp->second) * 2.0 * fn[1] * pp->second;
        for (CI pp = fn_z.begin(); pp != fn_z.end(); ++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0 / maxtrac) * (pp->second) * 2.0 * fn[2] * pp->second;

        CORE::GEN::pairedvector<int, double> trialnormlin(100);
        for (CI pp = ft_x.begin(); pp != ft_x.end(); ++pp)
          trialnormlin[pp->first] -=
              0.5 * (1.0 / trialnorm) * (pp->second) * 2.0 * ftrial[0] * pp->second;
        for (CI pp = ft_y.begin(); pp != ft_y.end(); ++pp)
          trialnormlin[pp->first] -=
              0.5 * (1.0 / trialnorm) * (pp->second) * 2.0 * ftrial[1] * pp->second;
        for (CI pp = ft_z.begin(); pp != ft_z.end(); ++pp)
          trialnormlin[pp->first] -=
              0.5 * (1.0 / trialnorm) * (pp->second) * 2.0 * ftrial[2] * pp->second;

        for (CI pp = maxtraclin.begin(); pp != maxtraclin.end(); ++pp)
          coefflin[pp->first] += frcoeff * pp->second * (1.0 / trialnorm);

        for (CI pp = trialnormlin.begin(); pp != trialnormlin.end(); ++pp)
          coefflin[pp->first] -= frcoeff * pp->second * (1.0 / (trialnorm * trialnorm)) * maxtrac;
      }


      std::map<int, std::map<int, double>>& dderiv = cnode->Data().GetDerivDltl();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->Data().GetDerivDltl()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = csnode->Dofs()[prodj];
          std::map<int, double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = scolcurr->first;
            double val = ftan[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->FEAssemble(val, row, col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr != thisdderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->Data().GetDerivMltl();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->Data().GetDerivMltl()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = ftan[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(-val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      // ****************************************************************
      // **************************************************** stick state
      // ****************************************************************
      if (stick)
      {
        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetDltl()).size() > 0)
        {
          CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDltl();

          for (CI p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < Dim(); ++dim)
            {
              for (CImap pp = cnode->Data().GetDerivJumpltl()[dim].begin();
                   pp != cnode->Data().GetDerivJumpltl()[dim].end(); ++pp)
              {
                double value = penaltytan * (p->second) * (pp->second);
                kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          dserror("no d matrix entries available for ltlt contact");
        }
        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetMltl()).size() > 0)
        {
          std::map<int, double> map = cnode->MoData().GetMltl();

          for (CImap p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < Dim(); ++dim)
            {
              for (CImap pp = cnode->Data().GetDerivJumpltl()[dim].begin();
                   pp != cnode->Data().GetDerivJumpltl()[dim].end(); ++pp)
              {
                double value = -penaltytan * (p->second) * (pp->second);
                kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          dserror("no m matrix entries available for ltlt contact");
        }
      }
      // ****************************************************************
      // **************************************************** slip state
      // ****************************************************************
      else
      {
        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetDltl()).size() > 0)
        {
          CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDltl();

          for (CI p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < Dim(); ++dim)
            {
              for (CImap pp = cnode->Data().GetDerivJumpltl()[dim].begin();
                   pp != cnode->Data().GetDerivJumpltl()[dim].end(); ++pp)
              {
                double value = penaltytan * coeff * (p->second) * (pp->second);
                kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          dserror("no d matrix entries available for ltlt contact");
        }

        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetDltl()).size() > 0)
        {
          CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDltl();

          for (CI p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < Dim(); ++dim)
            {
              for (CI pp = coefflin.begin(); pp != coefflin.end(); ++pp)
              {
                double value = -penaltytan * ftan[dim] * (p->second) * (pp->second);
                kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          dserror("no d matrix entries available for ltlt contact");
        }

        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetMltl()).size() > 0)
        {
          std::map<int, double> map = cnode->MoData().GetMltl();

          for (CImap p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < Dim(); ++dim)
            {
              for (CImap pp = cnode->Data().GetDerivJumpltl()[dim].begin();
                   pp != cnode->Data().GetDerivJumpltl()[dim].end(); ++pp)
              {
                double value = -penaltytan * coeff * (p->second) * (pp->second);
                kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          dserror("no m matrix entries available for ltlt contact");
        }
        if ((cnode->MoData().GetMltl()).size() > 0)
        {
          std::map<int, double> map = cnode->MoData().GetMltl();

          for (CImap p = map.begin(); p != map.end(); ++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("Cannot find node with gid");
            Node* csnode = dynamic_cast<Node*>(snode);

            for (int dim = 0; dim < Dim(); ++dim)
            {
              for (CI pp = coefflin.begin(); pp != coefflin.end(); ++pp)
              {
                double value = penaltytan * ftan[dim] * (p->second) * (pp->second);
                kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
              }
            }
          }
        }
        else
        {
          dserror("no m matrix entries available for ltlt contact");
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add nts penalty forces master                           farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddNTSforcesMaster(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for corners
    if (!cnode->IsOnCorner()) continue;

    // is gap is in contact
    if (cnode->Data().Getgnts() < 1e-12)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDnts()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDnts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value =
                penalty * (p->second) * cnode->Data().Getgnts() * cnode->MoData().n()[dim];
            int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMnts()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMnts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value =
                -penalty * (p->second) * cnode->Data().Getgnts() * cnode->MoData().n()[dim];
            int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces master                  farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddLTSforcesMaster(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for edges (without corner)
    if (!cnode->IsOnCornerEdge()) continue;
    if (cnode->IsOnCorner()) continue;

    // scale penalty
    double penaltyLts = penalty * cnode->Data().Kappa();

    // is gap is in contact
    if (cnode->Data().Getglts() < 1e-12)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDlts()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDlts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value =
                penaltyLts * (p->second) * cnode->Data().Getglts() * cnode->MoData().n()[dim];
            int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMlts()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMlts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value =
                -penaltyLts * (p->second) * cnode->Data().Getglts() * cnode->MoData().n()[dim];
            int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddLTLforces(Teuchos::RCP<Epetra_FEVector> feff)
{
  // gap = g_n * n
  // D/M = sval/mval
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if this is valid node
    if (cnode->Data().Getgltl()[0] < 1e8 and cnode->Data().Getgltl()[1] < 1e8 and
        cnode->Data().Getgltl()[2] < 1e8)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDltl()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDltl();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value = penalty * (p->second) * cnode->Data().Getgltl()[dim];
            int ltlid = csnode->Dofs()[dim];
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMltl()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMltl();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            double value = -penalty * (p->second) * cnode->Data().Getgltl()[dim];
            int ltlid = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1, &ltlid, &value);
            if (err < 0) dserror("stop");
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddLTSstiffnessMaster(Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for edges (without corner)
    if (!cnode->IsOnCornerEdge()) continue;
    if (cnode->IsOnCorner()) continue;

    // scale penalty
    double penaltyLts = penalty * cnode->Data().Kappa();

    // is gap is in contact
    if (cnode->Data().Getglts() < 1e-12)
    {
      std::array<double, 3> lm = {0.0, 0.0, 0.0};
      lm[0] = penaltyLts * cnode->Data().Getglts() * cnode->MoData().n()[0];
      lm[1] = penaltyLts * cnode->Data().Getglts() * cnode->MoData().n()[1];
      lm[2] = penaltyLts * cnode->Data().Getglts() * cnode->MoData().n()[2];

      std::map<int, std::map<int, double>>& dderiv = cnode->Data().GetDerivDlts();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->Data().GetDerivDlts()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = csnode->Dofs()[prodj];
          std::map<int, double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = scolcurr->first;
            double val = lm[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->FEAssemble(-val, row, col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr != thisdderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->Data().GetDerivMlts();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->Data().GetDerivMlts()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDlts()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDlts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            // gap linearization
            for (CImap pp = cnode->Data().GetDerivGlts().begin();
                 pp != cnode->Data().GetDerivGlts().end(); ++pp)
            {
              double value = -penaltyLts * (p->second) * (pp->second) * cnode->MoData().n()[dim];
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->Data().GetDerivN()[dim].begin();
                 pp != cnode->Data().GetDerivN()[dim].end(); ++pp)
            {
              double value = -penaltyLts * (p->second) * (pp->second) * cnode->Data().Getglts();
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetMlts()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMlts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            for (CImap pp = cnode->Data().GetDerivGlts().begin();
                 pp != cnode->Data().GetDerivGlts().end(); ++pp)
            {
              double value = penaltyLts * (p->second) * (pp->second) * cnode->MoData().n()[dim];
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->Data().GetDerivN()[dim].begin();
                 pp != cnode->Data().GetDerivN()[dim].end(); ++pp)
            {
              double value = penaltyLts * (p->second) * (pp->second) * cnode->Data().Getglts();
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddNTSstiffnessMaster(Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < mnoderowmap_->NumMyElements(); ++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // only for corners
    if (!cnode->IsOnCorner()) continue;

    // is gap is in contact
    if (cnode->Data().Getgnts() < 1e-12)
    {
      std::array<double, 3> lm = {0.0, 0.0, 0.0};
      lm[0] = penalty * cnode->Data().Getgnts() * cnode->MoData().n()[0];
      lm[1] = penalty * cnode->Data().Getgnts() * cnode->MoData().n()[1];
      lm[2] = penalty * cnode->Data().Getgnts() * cnode->MoData().n()[2];

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->Data().GetDerivMnts();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->Data().GetDerivMnts()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDnts()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDnts();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            // gap linearization
            for (CImap pp = cnode->Data().GetDerivGnts().begin();
                 pp != cnode->Data().GetDerivGnts().end(); ++pp)
            {
              double value = -penalty * (p->second) * (pp->second) * cnode->MoData().n()[dim];
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->Data().GetDerivN()[dim].begin();
                 pp != cnode->Data().GetDerivN()[dim].end(); ++pp)
            {
              double value = -penalty * (p->second) * (pp->second) * cnode->Data().Getgnts();
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetMnts()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMnts();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            for (CImap pp = cnode->Data().GetDerivGnts().begin();
                 pp != cnode->Data().GetDerivGnts().end(); ++pp)
            {
              double value = penalty * (p->second) * (pp->second) * cnode->MoData().n()[dim];
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
            // normal linearization
            for (CI pp = cnode->Data().GetDerivN()[dim].begin();
                 pp != cnode->Data().GetDerivN()[dim].end(); ++pp)
            {
              double value = penalty * (p->second) * (pp->second) * cnode->Data().Getgnts();
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::AddLTLstiffness(Teuchos::RCP<CORE::LINALG::SparseMatrix> kteff)
{
  const double penalty = InterfaceParams().get<double>("PENALTYPARAM");

  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  typedef std::map<int, double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // check if this is valid node
    if (cnode->Data().Getgltl()[0] < 1e8 and cnode->Data().Getgltl()[1] < 1e8 and
        cnode->Data().Getgltl()[2] < 1e8)
    {
      std::array<double, 3> lm = {0.0, 0.0, 0.0};
      lm[0] = penalty * cnode->Data().Getgltl()[0];
      lm[1] = penalty * cnode->Data().Getgltl()[1];
      lm[2] = penalty * cnode->Data().Getgltl()[2];

      std::map<int, std::map<int, double>>& dderiv = cnode->Data().GetDerivDltl();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->Data().GetDerivDltl()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = csnode->Dofs()[prodj];
          std::map<int, double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = scolcurr->first;
            double val = lm[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->FEAssemble(-val, row, col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr != thisdderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int, std::map<int, double>>& mderiv = cnode->Data().GetDerivMltl();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->Data().GetDerivMltl()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj = 0; prodj < Dim(); ++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int, double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c = 0; c < mapsize; ++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(val, row, col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr != thismderiv.end())
            dserror("AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDltl()).size() > 0)
      {
        CORE::GEN::pairedvector<int, double> map = cnode->MoData().GetDltl();

        for (CI p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            for (CImap pp = cnode->Data().GetDerivGltl()[dim].begin();
                 pp != cnode->Data().GetDerivGltl()[dim].end(); ++pp)
            {
              double value = -penalty * (p->second) * (pp->second);
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        dserror("no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetMltl()).size() > 0)
      {
        std::map<int, double> map = cnode->MoData().GetMltl();

        for (CImap p = map.begin(); p != map.end(); ++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("Cannot find node with gid");
          Node* csnode = dynamic_cast<Node*>(snode);

          for (int dim = 0; dim < Dim(); ++dim)
          {
            for (CImap pp = cnode->Data().GetDerivGltl()[dim].begin();
                 pp != cnode->Data().GetDerivGltl()[dim].end(); ++pp)
            {
              double value = penalty * (p->second) * (pp->second);
              kteff->FEAssemble(value, csnode->Dofs()[dim], pp->first);
            }
          }
        }
      }
      else
      {
        dserror("no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  post evaluate to scale calculated terms                 farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::PostEvaluate(const int step, const int iter)
{
  // decide which type of coupling should be evaluated
  INPAR::MORTAR::AlgorithmType algo =
      CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  switch (algo)
  {
    //*********************************
    // Mortar Coupling (STS)    (2D/3D)
    //*********************************
    case INPAR::MORTAR::algorithm_mortar:
    {
      // non-smooth contact
      if (nonSmoothContact_)
      {
        // store lts into mortar data container
        StoreLTSvalues();

        // store nts into mortar data container
        StoreNTSvalues();
      }
      return;
      break;
    }
    //*********************************
    // Gauss-Point-To-Segment (GPTS)
    //*********************************
    case INPAR::MORTAR::algorithm_gpts:
    {
      // already stored
      return;
      break;
    }
    //*********************************
    // Line-to-Segment Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_lts:
    {
      // store lts into mortar data container
      StoreLTSvalues();
      break;
    }
    //*********************************
    // Node-to-Segment Coupling (2D/3D)
    //*********************************
    case INPAR::MORTAR::algorithm_nts:
    {
      // store nts into mortar data container
      StoreNTSvalues();
      break;
    }
    //*********************************
    // line-to-line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_ltl:
    {
      return;
      break;
    }
    //*********************************
    // Node-to-Line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_ntl:
    {
      dserror("not yet implemented!");
      break;
    }
    //*********************************
    // Segment-to-Line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_stl:
    {
      // store lts into mortar data container
      StoreLTSvalues();
      break;
    }
    //*********************************
    // Default case
    //*********************************
    default:
    {
      dserror("Unknown discr. type for constraints!");
      break;
    }
  }

#ifdef MORTARGMSHCELLS
  // finish integration cell GMSH files
  int proc = Comm().MyPID();
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";
  FILE* fp = fopen(filename.str().c_str(), "a");
  std::stringstream gmshfilecontent2;
  gmshfilecontent2 << "};" << std::endl;
  fprintf(fp, gmshfilecontent2.str().c_str());
  fclose(fp);

  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream newfilename;
  newfilename << "o/gmsh_output/cells_";
  if (step < 10)
    newfilename << 0 << 0 << 0 << 0;
  else if (step < 100)
    newfilename << 0 << 0 << 0;
  else if (step < 1000)
    newfilename << 0 << 0;
  else if (step < 10000)
    newfilename << 0;
  else if (step > 99999)
    dserror("Gmsh output implemented for a maximum of 99.999 time steps");
  newfilename << step;

  // second index = Newton iteration index
  newfilename << "_";
  if (iter < 10)
    newfilename << 0;
  else if (iter > 99)
    dserror("Gmsh output implemented for a maximum of 99 iterations");
  newfilename << iter << "_p" << proc << ".pos";

  // rename file
  rename(filename.str().c_str(), newfilename.str().c_str());
#endif  // #ifdef MORTARGMSHCELLS

  return;
}

/*----------------------------------------------------------------------*
 |  scale calculated entries for nts, mortar etc...         farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ScaleTerms()
{
  // only for 3d contact
  if (Dim() == 2) return;

  // scale ltl terms for point and line contact
  ScaleTermsLTL();

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  scale calculated entries for lts contact                farah 09/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ScaleTermsLTL()
{
  std::cout << "ScaleTerms LTL" << std::endl;
  dserror("ScaleTermsLTL() is outdated");

  return;
  //  // create iterators for data types
  //  typedef CORE::GEN::pairedvector<int,double>::const_iterator CI;
  //  typedef std::map<int,double>::const_iterator          CImap;
  //
  //  // loop over all possibly non smooth nodes
  //  for(int i=0; i<snoderowmap_->NumMyElements();++i)
  //  {
  //    int gid = snoderowmap_->GID(i);
  //    DRT::Node* node = idiscret_->gNode(gid);
  //    if (!node)
  //      dserror("Cannot find node with gid %",gid);
  //    Node* cnode = dynamic_cast<Node*>(node);
  //
  //    // only for edge nodes
  //    if(!cnode->IsOnEdge())
  //      continue;
  //
  //    // get scale factor alpha
  //    const double alpha = cnode->Data().GetAlphaN();
  //    if(alpha<-1e-12)
  //      continue;
  //
  //    std::cout << "PERFORM SCALING" << std::endl;
  //
  //    // check if integration is done
  ////    if (cnode->MoData().GetD().size()<1)
  ////      continue;
  //
  //    // ############################################################
  //    //                   create aux terms
  //    // ############################################################
  //    std::map<int, double> derivg_lts =
  //        cnode->Data().GetDerivGlts();
  //
  //    std::map<int, std::map<int, double> >derivm_lts =
  //        cnode->Data().GetDerivMlts();
  //
  //    std::map<int, std::map<int, double> > derivd_lts =
  //        cnode->Data().GetDerivDlts();
  //
  //    double g_lts = cnode->Data().Getglts();
  //
  //    std::map<int, double> m_lts =
  //        cnode->MoData().GetMlts();
  //
  //    CORE::GEN::pairedvector<int, double> d_lts =
  //        cnode->MoData().GetDlts();
  //
  //    // ############################################################
  //    //                   scale lts terms:
  //    // ############################################################
  //    {
  //      //-------------------------------------------------------------------------------------
  //      // scale M matrix entries for mortar data
  //      for (CImap p=m_lts.begin();p!=m_lts.end();++p)
  //        cnode->MoData().GetMlts()[p->first] = alpha*(p->second);
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale weighted gap
  //      cnode->Data().Getglts() = alpha * g_lts;
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale M deriv
  //      {
  //        // Mortar M derivatives
  //        std::map<int,std::map<int,double> >& mderiv = cnode->Data().GetDerivMlts();
  //
  //        // get sizes and iterator start
  //        int mastersize = (int)mderiv.size();
  //        std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();
  //
  //        /********************************************** LinMMatrix **********/
  //        // loop over all master nodes in the DerivM-map of the current LM slave node
  //        for (int l=0;l<mastersize;++l)
  //        {
  //          int mgid = mcurr->first;
  //          ++mcurr;
  //
  //          // Mortar matrix M derivatives
  //          std::map<int,double>&thismderiv = cnode->Data().GetDerivMlts()[mgid];
  //          std::map<int,double>&auxderiv   = derivm_lts[mgid];
  //
  //          int mapsize = (int)(thismderiv.size());
  //
  //          // inner product M_{jl,c} * z_j for index j
  //          for (int prodj=0;prodj<Dim();++prodj)
  //          {
  //            std::map<int,double>::iterator mcolcurr    = thismderiv.begin();
  //            std::map<int,double>::iterator mcolcurraux = auxderiv.begin();
  //
  //            // loop over all directional derivative entries
  //            for (int c=0;c<mapsize;++c)
  //            {
  //              mcolcurr->second = alpha * (mcolcurraux->second);
  //              ++mcolcurr;
  //            }
  //
  //            // check for completeness of DerivM-Derivatives-iteration
  //            if (mcolcurr!=thismderiv.end())
  //              dserror("ScaleTerms: Not all derivative entries of DerivM considered!");
  //          }
  //        }
  //      }
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale alpha lin
  //      {
  //        for (CI p=cnode->Data().GetAlpha().begin();p!=cnode->Data().GetAlpha().end();++p)
  //        {
  //          for (CImap pp=m_lts.begin();pp!=m_lts.end();++pp)
  //          {
  //            (cnode->Data().GetDerivMlts()[pp->first])[p->first] += (p->second) * (pp->second);
  //          }
  //        }
  //      }
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale gap derivative entries for mortar data
  //      for (CImap p=derivg_lts.begin();p!=derivg_lts.end();++p)
  //        cnode->Data().GetDerivGlts()[p->first] = (p->second)*alpha;
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale gap derivative with lin alpha entries for mortar data
  //      for (CI p=cnode->Data().GetAlpha().begin();p!=cnode->Data().GetAlpha().end();++p)
  //      {
  //        cnode->Data().GetDerivGlts()[p->first] += g_lts *(p->second);
  //      }
  //    }// end lts terms
  //
  //    // ############################################################
  //    // scale ltl terms and add to lts terms:
  //    // ############################################################
  //    {
  //      // get scale factor alpha
  //      const double alphaNew = 1.0 - cnode->Data().GetAlphaN();
  //
  //      if(cnode->MoData().GetDlts().size()>1)
  //        dserror("D matrix must be diagonal!");
  //
  //      //-------------------------------------------------------------------------------------
  //      const double dvalue = d_lts.begin()->second;
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale M matrix entries for mortar data --> Ntilde * D * (1-a)
  //      for (CImap p=cnode->MoData().GetMltl().begin();p!=cnode->MoData().GetMltl().end();++p)
  //        cnode->MoData().GetMlts()[p->first] += alphaNew*(p->second)*dvalue;
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale weighted gap
  //      cnode->Data().Getglts() += alphaNew *  cnode->Data().Getgltl() * dvalue;
  //
  //      //-------------------------------------------------------------------------------------
  //      // Mortar matrix D and M derivatives
  //      std::map<int,std::map<int,double> >& dderiv = cnode->Data().GetDerivDlts();
  //
  //      // get sizes and iterator start
  //      int slavesize = (int)dderiv.size();
  ////      if (slavesize>1)
  ////        dserror("wrong dlin dimension");
  //
  //      std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();
  //
  //      /********************************************** LinDMatrix **********/
  //      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
  //      for (int k=0;k<slavesize;++k)
  //      {
  //        int sgid = scurr->first;
  //        ++scurr;
  //
  //        DRT::Node* snode = idiscret_->gNode(sgid);
  //        if (!snode) dserror("Cannot find node with gid %",sgid);
  //
  //        // Mortar matrix D derivatives
  //        std::map<int,double>& thisdderivD        = derivd_lts[sgid];
  //        std::map<int,double>& thismderivmortarM  = cnode->Data().GetDerivMlts()[sgid];
  //
  //        int mapsize = (int)(thisdderivD.size());
  //
  //        // inner product D_{jk,c} * z_j for index j
  ////        for (int prodj=0;prodj<Dim();++prodj)
  ////        {
  //          std::map<int,double>::iterator scolcurrD = thisdderivD.begin();
  //
  //          // loop over all directional derivative entries
  //          for (int c=0;c<mapsize;++c)
  //          {
  //            // derivG
  //            cnode->Data().GetDerivGlts()[scolcurrD->first] +=
  //            alphaNew*(scolcurrD->second)*cnode->Data().Getgltl();
  //
  //            //derivM with delta D
  //            for (CImap
  //            p=cnode->MoData().GetMltl().begin();p!=cnode->MoData().GetMltl().end();++p)
  //              thismderivmortarM[scolcurrD->first] += alphaNew*(scolcurrD->second) * p->second;
  //            ++scolcurrD;
  //          }
  //
  //          // check for completeness of DerivD-Derivatives-iteration
  //          if (scolcurrD!=thisdderivD.end())
  //            dserror("ScaleTerms: Not all derivative entries of DerivD considered!");
  //        //}
  //      }
  //
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale M LTL deriv
  //      {
  //        // Mortar M derivatives
  //        std::map<int,std::map<int,double> >& mntsderiv = cnode->Data().GetDerivMltl();
  //
  //        // get sizes and iterator start
  //        int mastersize = (int)mntsderiv.size();
  //        std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();
  //
  //        /********************************************** LinMMatrix **********/
  //        // loop over all master nodes in the DerivM-map of the current LM slave node
  //        for (int l=0;l<mastersize;++l)
  //        {
  //          int mgid = mntscurr->first;
  //          ++mntscurr;
  //
  //          // Mortar matrix M derivatives
  //          std::map<int,double>&thismderivnts    = cnode->Data().GetDerivMltl()[mgid];
  //          std::map<int,double>&thismderivmortar = cnode->Data().GetDerivMlts()[mgid];
  //
  //          int mapsize = (int)(thismderivnts.size());
  //
  //          // inner product M_{jl,c} * z_j for index j
  ////          for (int prodj=0;prodj<Dim();++prodj)
  ////          {
  //            std::map<int,double>::iterator mcolcurr = thismderivnts.begin();
  //
  //            // loop over all directional derivative entries
  //            for (int c=0;c<mapsize;++c)
  //            {
  //              thismderivmortar[mcolcurr->first] += alphaNew * (mcolcurr->second) * dvalue;
  //              ++mcolcurr;
  //            }
  //
  //            // check for completeness of DerivM-Derivatives-iteration
  //            if (mcolcurr!=thismderivnts.end())
  //              dserror("ScaleTerms: Not all derivative entries of DerivM considered!");
  ////          }
  //        }
  //      }
  //
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale alpha lin TODO check this
  //      {
  //        for (CI p=cnode->Data().GetAlpha().begin();p!=cnode->Data().GetAlpha().end();++p)
  //        {
  //          for (CImap
  //          pp=cnode->MoData().GetMltl().begin();pp!=cnode->MoData().GetMltl().end();++pp)
  //          {
  //            (cnode->Data().GetDerivMlts()[pp->first])[p->first] -= (p->second) * (pp->second)
  //            * dvalue;
  //          }
  //        }
  //      }
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale gap derivative entries for mortar data
  //      for (CImap
  //      p=cnode->Data().GetDerivGltl().begin();p!=cnode->Data().GetDerivGltl().end();++p)
  //        cnode->Data().GetDerivGlts()[p->first] += alphaNew*(p->second)*dvalue;
  //
  //      //-------------------------------------------------------------------------------------
  //      // scale gap derivative with lin alpha entries for mortar data
  //      for (CI p=cnode->Data().GetAlpha().begin();p!=cnode->Data().GetAlpha().end();++p)
  //      {
  //        cnode->Data().GetDerivGlts()[p->first] -= cnode->Data().Getgltl()* dvalue
  //        *(p->second);
  //      }
  //    }// end scale nts
  //  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-segment coupl          farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateSTS(
    const Epetra_Map& selecolmap, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  MORTAR::Interface::EvaluateSTS(selecolmap, mparams_ptr);
  return;
  //  // loop over all slave col elements
  //  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  //  {
  //    int gid1 = selecolmap_->GID(i);
  //    DRT::Element* ele1 = idiscret_->gElement(gid1);
  //    if (!ele1)
  //      dserror("Cannot find slave element with gid %", gid1);
  //    MORTAR::Element* selement = dynamic_cast<MORTAR::Element*>(ele1);
  //
  //    // loop over all slave nodes of this element to check for active nodes
  //    bool eval = true;
  //    for(int k = 0; k<selement->NumNode(); ++k)
  //    {
  //      CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(selement->Nodes()[k]);
  //      if(cnode->Active())
  //        eval=false;
  //    }
  //    if(!eval)
  //      continue;
  //
  //    // skip zero-sized nurbs elements (slave)
  //    if (selement->ZeroSized())
  //      continue;
  //
  //    // empty vector of master element pointers
  //    std::vector<MORTAR::Element*> melements;
  //
  //    // loop over the candidate master elements of sele_
  //    // use slave element's candidate list SearchElements !!!
  //    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
  //    {
  //      int gid2 = selement->MoData().SearchElements()[j];
  //      DRT::Element* ele2 = idiscret_->gElement(gid2);
  //      if (!ele2)
  //        dserror("Cannot find master element with gid %", gid2);
  //      MORTAR::Element* melement = dynamic_cast<MORTAR::Element*>(ele2);
  //
  //      // skip zero-sized nurbs elements (master)
  //      if (melement->ZeroSized())
  //        continue;
  //
  //      melements.push_back(melement);
  //    }
  //
  //    // concrete coupling evaluation routine
  //    MortarCoupling(selement,melements,mparams_ptr);
  //  }

  return;
}

/*----------------------------------------------------------------------*
 |  protected evaluate routine                               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateCoupling(const Epetra_Map& selecolmap,
    const Epetra_Map* snoderowmap, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // ask if non-smooth contact is activated!
  if (nonSmoothContact_)
  {
    // 2D: only STS and nts has to be performed
    if (Dim() == 2)
    {
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateSTS(selecolmap, mparams_ptr);

      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateNTS();

      //********************************************************************
      // NTN is a special case of NTS and an additional implementation is
      // not required!
      //********************************************************************
    }
    else if (Dim() == 3)
    {
      //********************************************************************
      // TODO: remove this hack!
      // HACK: LTL is not yet included in nonsmooth contact framework!
      // However, we want to test the LTL code separately. Thus, the "if"-
      // statement is included:
      // decide which type of coupling should be evaluated
      //********************************************************************
      INPAR::MORTAR::AlgorithmType algo =
          CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");
      if (algo == INPAR::MORTAR::algorithm_ltl)
      {
        EvaluateLTL();
        return;
      }

      //********************************************************************
      // 1) try to project slave nodes onto master elements
      // 2) evaluate shape functions at projected positions
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      EvaluateNTS();

      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //********************************************************************
      EvaluateLTL();

      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //********************************************************************
      EvaluateLTS();

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      EvaluateSTS(selecolmap, mparams_ptr);

      //********************************************************************
      // perform LTS steps for master edges
      //********************************************************************
      //      EvaluateLTSMaster();

      //********************************************************************
      // perform NTS steps for master edges
      //********************************************************************
      //      EvaluateNTSMaster();

      //********************************************************************
      // NTN is a special case of NTS and an additional implementation is
      // not required!
      //********************************************************************
    }
    else
    {
      dserror("Wrong dimension!");
    }
  }
  else
  {
    //********************************************************************
    // call base routine for standard mortar/nts evaluation
    //********************************************************************
    MORTAR::Interface::EvaluateCoupling(selecolmap, snoderowmap, mparams_ptr);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Check and initialize corner/edge contact                 farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::InitializeCornerEdge()
{
  // return if nonsmooth contact is activated
  if (nonSmoothContact_) return;

  // call base function
  MORTAR::Interface::InitializeCornerEdge();

  return;
}


/*----------------------------------------------------------------------*
 |  stuff for non-smooth contact geometries                 farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::DetectNonSmoothGeometries()
{
  dserror("outdated!");

  std::vector<int> nonsmoothnodegids(0);
  std::vector<int> smoothnodegids(0);

  // loop over slave nodes to find nodes and eles
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("Node ownership inconsistency!");

    if (cnode->NumElement() < 2)
    {
      nonsmoothnodegids.push_back(cnode->Id());
      continue;
    }

    int loc = 0;
    std::array<double, 3> normalsele0 = {0.0, 0.0, 0.0};
    std::array<double, 3> normalsele1 = {0.0, 0.0, 0.0};

    // build normal at node for 1. ele
    Element* sele0 = dynamic_cast<Element*>(cnode->Elements()[0]);
    CORE::LINALG::SerialDenseMatrix elen0(6, 1);
    sele0->BuildNormalAtNode(cnode->Id(), loc, elen0);

    // create 1. unit normal
    normalsele0[0] = elen0(0, 0) / elen0(4, 0);
    normalsele0[1] = elen0(1, 0) / elen0(4, 0);
    normalsele0[2] = elen0(2, 0) / elen0(4, 0);

    // build normal at node for 2. ele
    loc = 0;
    Element* sele1 = dynamic_cast<Element*>(cnode->Elements()[1]);
    CORE::LINALG::SerialDenseMatrix elen1(6, 1);
    sele1->BuildNormalAtNode(cnode->Id(), loc, elen1);

    // create 2. unit normal
    normalsele1[0] = elen1(0, 0) / elen1(4, 0);
    normalsele1[1] = elen1(1, 0) / elen1(4, 0);
    normalsele1[2] = elen1(2, 0) / elen1(4, 0);

    double dot = normalsele0[0] * normalsele1[0] + normalsele0[1] * normalsele1[1] +
                 normalsele0[2] * normalsele1[2];
    double R = abs(dot);  // this is the abs ( cos(phi) ) --> 1.0 if parallel; 0.0 if orthogonal

    // critical node
    if (R < 0.8)  // circa 90 - 11.5 = 78.5
      nonsmoothnodegids.push_back(cnode->Id());
    else
      smoothnodegids.push_back(cnode->Id());
  }

  // create maps
  nonsmoothnodes_ = CORE::LINALG::CreateMap(nonsmoothnodegids, Comm());
  smoothnodes_ = CORE::LINALG::CreateMap(smoothnodegids, Comm());

  // debug output
  std::cout << "# nonsmooth nodes: " << nonsmoothnodes_->NumGlobalElements() << std::endl;
  std::cout << "#    smooth nodes: " << smoothnodes_->NumGlobalElements() << std::endl;

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  cpp to edge + Lin                                       farah 11/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::ComputeNormalNodeToEdge(MORTAR::Node& snode, MORTAR::Element& mele,
    double* normal, std::vector<CORE::GEN::pairedvector<int, double>>& normaltonodelin)
{
  // define tolerance
  const double tol = 1e-8;
  double dist = 1e12;
  int nrow = mele.NumNode();

  Node* node1 = dynamic_cast<Node*>(mele.Nodes()[0]);
  Node* node2 = dynamic_cast<Node*>(mele.Nodes()[1]);

  double length1 = sqrt(node1->MoData().EdgeTangent()[0] * node1->MoData().EdgeTangent()[0] +
                        node1->MoData().EdgeTangent()[1] * node1->MoData().EdgeTangent()[1] +
                        node1->MoData().EdgeTangent()[2] * node1->MoData().EdgeTangent()[2]);
  double length2 = sqrt(node2->MoData().EdgeTangent()[0] * node2->MoData().EdgeTangent()[0] +
                        node2->MoData().EdgeTangent()[1] * node2->MoData().EdgeTangent()[1] +
                        node2->MoData().EdgeTangent()[2] * node2->MoData().EdgeTangent()[2]);

  if (length1 < 1e-12 or length2 < 1e-12) return dist;

  // calc angle between tangents
  std::array<double, 3> t1 = {0.0, 0.0, 0.0};
  std::array<double, 3> t2 = {0.0, 0.0, 0.0};

  t1[0] = node1->MoData().EdgeTangent()[0];
  t1[1] = node1->MoData().EdgeTangent()[1];
  t1[2] = node1->MoData().EdgeTangent()[2];

  t2[0] = node2->MoData().EdgeTangent()[0];
  t2[1] = node2->MoData().EdgeTangent()[1];
  t2[2] = node2->MoData().EdgeTangent()[2];

  double test = t1[0] * t2[0] + t1[1] * t2[1] + t1[2] * t2[2];
  if (test < tol) dserror("tangents have wrong direction!");


  double f = 0.0;
  double df = 0.0;
  double xi = 0.0;

  // newton loop
  for (int k = 0; k < MORTARMAXITER; ++k)
  {
    //    std::cout << "k= " << k << std::endl;
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************
    CORE::LINALG::SerialDenseVector sval(nrow);
    CORE::LINALG::SerialDenseMatrix sderiv(nrow, 1);
    mele.EvaluateShape(&xi, sval, sderiv, nrow);

    // tangent part
    std::array<double, 3> tangent = {0.0, 0.0, 0.0};
    tangent[0] += sval[0] * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[0];
    tangent[1] += sval[0] * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[1];
    tangent[2] += sval[0] * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[2];

    tangent[0] += sval[1] * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[0];
    tangent[1] += sval[1] * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[1];
    tangent[2] += sval[1] * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[2];

    double tangentSlave = 0.0;
    tangentSlave = tangent[0] * snode.xspatial()[0] + tangent[1] * snode.xspatial()[1] +
                   tangent[2] * snode.xspatial()[2];

    // master part
    std::array<double, 3> master = {0.0, 0.0, 0.0};
    master[0] += sval[0] * node1->xspatial()[0];
    master[1] += sval[0] * node1->xspatial()[1];
    master[2] += sval[0] * node1->xspatial()[2];

    master[0] += sval[1] * node2->xspatial()[0];
    master[1] += sval[1] * node2->xspatial()[1];
    master[2] += sval[1] * node2->xspatial()[2];

    double tangentMaster = 0.0;
    tangentMaster = tangent[0] * master[0] + tangent[1] * master[1] + tangent[2] * master[2];

    f = tangentSlave - tangentMaster;
    if (abs(f) <= MORTARCONVTOL) break;
    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // lin tangent part
    std::array<double, 3> lintangent = {0.0, 0.0, 0.0};
    lintangent[0] += sderiv(0, 0) * node1->MoData().EdgeTangent()[0];
    lintangent[1] += sderiv(0, 0) * node1->MoData().EdgeTangent()[1];
    lintangent[2] += sderiv(0, 0) * node1->MoData().EdgeTangent()[2];

    lintangent[0] += sderiv(1, 0) * node2->MoData().EdgeTangent()[0];
    lintangent[1] += sderiv(1, 0) * node2->MoData().EdgeTangent()[1];
    lintangent[2] += sderiv(1, 0) * node2->MoData().EdgeTangent()[2];

    double lintangentSlave = 0.0;
    lintangentSlave = lintangent[0] * snode.xspatial()[0] + lintangent[1] * snode.xspatial()[1] +
                      lintangent[2] * snode.xspatial()[2];

    // lin master part
    std::array<double, 3> linmaster = {0.0, 0.0, 0.0};
    linmaster[0] += sderiv(0, 0) * node1->xspatial()[0];
    linmaster[1] += sderiv(0, 0) * node1->xspatial()[1];
    linmaster[2] += sderiv(0, 0) * node1->xspatial()[2];

    linmaster[0] += sderiv(1, 0) * node2->xspatial()[0];
    linmaster[1] += sderiv(1, 0) * node2->xspatial()[1];
    linmaster[2] += sderiv(1, 0) * node2->xspatial()[2];

    double lintangentMaster = 0.0;
    lintangentMaster =
        lintangent[0] * master[0] + lintangent[1] * master[1] + lintangent[2] * master[2];

    double tangentlinMaster = 0.0;
    tangentlinMaster =
        tangent[0] * linmaster[0] + tangent[1] * linmaster[1] + tangent[2] * linmaster[2];

    df = lintangentSlave - lintangentMaster - tangentlinMaster;
    if (abs(df) < 1e-12) dserror("df zero");
    xi += -f / df;
  }

  //**********************************************
  //   CHECK XI                                 //
  //**********************************************
  if (-1.0 - tol > xi or xi > 1.0 + tol) return dist;

  //**********************************************
  //   LINEARIZATION   df                       //
  //**********************************************

  CORE::LINALG::SerialDenseVector sval(nrow);
  CORE::LINALG::SerialDenseMatrix sderiv(nrow, 1);
  mele.EvaluateShape(&xi, sval, sderiv, nrow);

  // tangent part
  std::array<double, 3> tangent = {0.0, 0.0, 0.0};
  tangent[0] += sval[0] * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[0];
  tangent[1] += sval[0] * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[1];
  tangent[2] += sval[0] * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[2];

  tangent[0] += sval[1] * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[0];
  tangent[1] += sval[1] * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[1];
  tangent[2] += sval[1] * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[2];

  // master part
  std::array<double, 3> master = {0.0, 0.0, 0.0};
  master[0] += sval[0] * node1->xspatial()[0];
  master[1] += sval[0] * node1->xspatial()[1];
  master[2] += sval[0] * node1->xspatial()[2];

  master[0] += sval[1] * node2->xspatial()[0];
  master[1] += sval[1] * node2->xspatial()[1];
  master[2] += sval[1] * node2->xspatial()[2];

  // lin tangent part
  std::array<double, 3> lintangent = {0.0, 0.0, 0.0};
  lintangent[0] += sderiv(0, 0) * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[0];
  lintangent[1] += sderiv(0, 0) * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[1];
  lintangent[2] += sderiv(0, 0) * dynamic_cast<Node*>(mele.Nodes()[0])->MoData().EdgeTangent()[2];

  lintangent[0] += sderiv(1, 0) * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[0];
  lintangent[1] += sderiv(1, 0) * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[1];
  lintangent[2] += sderiv(1, 0) * dynamic_cast<Node*>(mele.Nodes()[1])->MoData().EdgeTangent()[2];

  // lin master part
  std::array<double, 3> linmaster = {0.0, 0.0, 0.0};
  linmaster[0] += sderiv(0, 0) * node1->xspatial()[0];
  linmaster[1] += sderiv(0, 0) * node1->xspatial()[1];
  linmaster[2] += sderiv(0, 0) * node1->xspatial()[2];

  linmaster[0] += sderiv(1, 0) * node2->xspatial()[0];
  linmaster[1] += sderiv(1, 0) * node2->xspatial()[1];
  linmaster[2] += sderiv(1, 0) * node2->xspatial()[2];


  //**********************************************
  //   LINEARIZATION    f                       //
  //**********************************************
  typedef CORE::GEN::pairedvector<int, double>::const_iterator _CI;

  std::vector<CORE::GEN::pairedvector<int, double>> linT(3, 100);  // added all sizes

  for (_CI p = node1->Data().GetDerivTangent()[0].begin();
       p != node1->Data().GetDerivTangent()[0].end(); ++p)
    linT[0][p->first] += sval[0] * p->second;
  for (_CI p = node1->Data().GetDerivTangent()[1].begin();
       p != node1->Data().GetDerivTangent()[1].end(); ++p)
    linT[1][p->first] += sval[0] * p->second;
  for (_CI p = node1->Data().GetDerivTangent()[2].begin();
       p != node1->Data().GetDerivTangent()[2].end(); ++p)
    linT[2][p->first] += sval[0] * p->second;

  for (_CI p = node2->Data().GetDerivTangent()[0].begin();
       p != node2->Data().GetDerivTangent()[0].end(); ++p)
    linT[0][p->first] += sval[1] * p->second;
  for (_CI p = node2->Data().GetDerivTangent()[1].begin();
       p != node2->Data().GetDerivTangent()[1].end(); ++p)
    linT[1][p->first] += sval[1] * p->second;
  for (_CI p = node2->Data().GetDerivTangent()[2].begin();
       p != node2->Data().GetDerivTangent()[2].end(); ++p)
    linT[2][p->first] += sval[1] * p->second;

  std::vector<CORE::GEN::pairedvector<int, double>> linXsl(3, 100);  // added all sizes
  linXsl[0][snode.Dofs()[0]] += 1.0;
  linXsl[1][snode.Dofs()[1]] += 1.0;
  linXsl[2][snode.Dofs()[2]] += 1.0;

  std::vector<CORE::GEN::pairedvector<int, double>> linXm(3, 100);  // added all sizes
  linXm[0][node1->Dofs()[0]] += sval[0];
  linXm[1][node1->Dofs()[1]] += sval[0];
  linXm[2][node1->Dofs()[2]] += sval[0];

  linXm[0][node2->Dofs()[0]] += sval[1];
  linXm[1][node2->Dofs()[1]] += sval[1];
  linXm[2][node2->Dofs()[2]] += sval[1];

  CORE::GEN::pairedvector<int, double> linf(100);  // added all sizes
  for (_CI p = linT[0].begin(); p != linT[0].end(); ++p)
    linf[p->first] += snode.xspatial()[0] * p->second;
  for (_CI p = linT[1].begin(); p != linT[1].end(); ++p)
    linf[p->first] += snode.xspatial()[1] * p->second;
  for (_CI p = linT[2].begin(); p != linT[2].end(); ++p)
    linf[p->first] += snode.xspatial()[2] * p->second;

  for (_CI p = linXsl[0].begin(); p != linXsl[0].end(); ++p)
    linf[p->first] += tangent[0] * p->second;
  for (_CI p = linXsl[1].begin(); p != linXsl[1].end(); ++p)
    linf[p->first] += tangent[1] * p->second;
  for (_CI p = linXsl[2].begin(); p != linXsl[2].end(); ++p)
    linf[p->first] += tangent[2] * p->second;

  for (_CI p = linT[0].begin(); p != linT[0].end(); ++p) linf[p->first] -= master[0] * p->second;
  for (_CI p = linT[1].begin(); p != linT[1].end(); ++p) linf[p->first] -= master[1] * p->second;
  for (_CI p = linT[2].begin(); p != linT[2].end(); ++p) linf[p->first] -= master[2] * p->second;

  for (_CI p = linXm[0].begin(); p != linXm[0].end(); ++p) linf[p->first] -= tangent[0] * p->second;
  for (_CI p = linXm[1].begin(); p != linXm[1].end(); ++p) linf[p->first] -= tangent[1] * p->second;
  for (_CI p = linXm[2].begin(); p != linXm[2].end(); ++p) linf[p->first] -= tangent[2] * p->second;

  CORE::GEN::pairedvector<int, double> linXi(100);  // added all sizes
  for (_CI p = linf.begin(); p != linf.end(); ++p) linXi[p->first] -= p->second / df;

  //**********************************************
  //   CALC NORMAL                              //
  //**********************************************
  std::array<double, 3> auxnormal = {0.0, 0.0, 0.0};
  auxnormal[0] = snode.xspatial()[0] - master[0];
  auxnormal[1] = snode.xspatial()[1] - master[1];
  auxnormal[2] = snode.xspatial()[2] - master[2];

  // remove numerical artifacts
  //  if(abs(auxnormal[0])<1e-12)
  //    auxnormal[0] = 0.0;
  //  if(abs(auxnormal[1])<1e-12)
  //    auxnormal[1] = 0.0;
  //  if(abs(auxnormal[2])<1e-12)
  //    auxnormal[2] = 0.0;

  // calc distance
  dist =
      sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] + auxnormal[2] * auxnormal[2]);

  if (abs(dist) < 1e-12) return 1e12;

  //*******************************************
  // Lin:
  std::vector<CORE::GEN::pairedvector<int, double>> auxlin(3, 100);  // added all sizes

  // xslave
  for (int k = 0; k < 3; ++k) (auxlin[k])[snode.Dofs()[k]] += 1.0;

  // xmaster n1
  for (int k = 0; k < 3; ++k) (auxlin[k])[node1->Dofs()[k]] -= sval[0];
  // xmaster n2
  for (int k = 0; k < 3; ++k) (auxlin[k])[node2->Dofs()[k]] -= sval[1];

  for (_CI p = linXi.begin(); p != linXi.end(); ++p)
  {
    for (int k = 0; k < 3; ++k)
    {
      (auxlin[k])[p->first] -= sderiv(0, 0) * node1->xspatial()[k] * p->second;
      (auxlin[k])[p->first] -= sderiv(1, 0) * node2->xspatial()[k] * p->second;
    }
  }

  normal[0] = auxnormal[0];
  normal[1] = auxnormal[1];
  normal[2] = auxnormal[2];

  //******************************
  // Orientation check:
  std::array<double, 3> slavebasednormal = {0.0, 0.0, 0.0};
  int nseg = snode.NumElement();
  DRT::Element** adjeles = snode.Elements();

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************
  CORE::LINALG::SerialDenseMatrix elens(6, nseg);
  MORTAR::Element* adjmrtrele = dynamic_cast<MORTAR::Element*>(adjeles[0]);

  // build element normal at current node
  // (we have to pass in the index i to be able to store the
  // normal and other information at the right place in elens)
  int i = 0;
  adjmrtrele->BuildNormalAtNode(snode.Id(), i, elens);

  // add (weighted) element normal to nodal normal n
  for (int j = 0; j < 3; ++j) slavebasednormal[j] += elens(j, 0) / elens(4, 0);

  // create unit normal vector
  const double length =
      sqrt(slavebasednormal[0] * slavebasednormal[0] + slavebasednormal[1] * slavebasednormal[1] +
           slavebasednormal[2] * slavebasednormal[2]);
  if (abs(length) < 1e-12)
  {
    dserror("Nodal normal length 0, node ID %i", snode.Id());
  }
  else
  {
    for (int j = 0; j < 3; ++j) slavebasednormal[j] /= length;
  }

  const double dotprod =
      -(normal[0] / dist * slavebasednormal[0] + normal[1] / dist * slavebasednormal[1] +
          normal[2] / dist * slavebasednormal[2]);

  if (dotprod < -1e-12)
  {
    // get the cpp normal
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];

    for (int j = 0; j < 3; ++j)
      for (_CI p = auxlin[j].begin(); p != auxlin[j].end(); ++p)
        (normaltonodelin[j])[p->first] -= (p->second);
  }
  else
  {
    // linearization
    for (int j = 0; j < 3; ++j)
      for (_CI p = auxlin[j].begin(); p != auxlin[j].end(); ++p)
        (normaltonodelin[j])[p->first] += (p->second);
  }

  return dist;
}

/*----------------------------------------------------------------------*
 |  cpp to node + Lin                                       farah 01/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::ComputeNormalNodeToNode(MORTAR::Node& snode, MORTAR::Node& mnode,
    double* normal, std::vector<CORE::GEN::pairedvector<int, double>>& normaltonodelin)
{
  const int dim = Dim();

  // distance between node and surface
  double gdist = 1e12;
  std::array<double, 3> gnormal = {0.0, 0.0, 0.0};
  std::vector<CORE::GEN::pairedvector<int, double>> glin(3, 1000);
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;

  double dist = 1e12;
  std::array<double, 3> auxnormal = {0.0, 0.0, 0.0};

  // loop over found master nodes
  std::vector<CORE::GEN::pairedvector<int, double>> auxlin(3, 1000);

  // calc vector
  auxnormal[0] = snode.xspatial()[0] - mnode.xspatial()[0];
  auxnormal[1] = snode.xspatial()[1] - mnode.xspatial()[1];
  auxnormal[2] = snode.xspatial()[2] - mnode.xspatial()[2];

  // remove numerical artifacts
  if (abs(auxnormal[0]) < 1e-12) auxnormal[0] = 0.0;
  if (abs(auxnormal[1]) < 1e-12) auxnormal[1] = 0.0;
  if (abs(auxnormal[2]) < 1e-12) auxnormal[2] = 0.0;

  // calc distance
  dist =
      sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] + auxnormal[2] * auxnormal[2]);

  // if nodes lying on each other: continue to next master node
  if (abs(dist) < 1e-12) return dist;

  //*******************************************
  // Lin:
  // xslave
  for (int k = 0; k < dim; ++k) (auxlin[k])[snode.Dofs()[k]] += 1.0;

  // xmaster
  for (int k = 0; k < dim; ++k) (auxlin[k])[mnode.Dofs()[k]] -= 1.0;

  // get normal
  gdist = dist;

  // normalize vector
  gnormal[0] = auxnormal[0];  /// dist;
  gnormal[1] = auxnormal[1];  /// dist;
  gnormal[2] = auxnormal[2];  /// dist;

  // linearization
  glin = auxlin;

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];

  //******************************
  // Orientation check:
  std::array<double, 3> slavebasednormal = {0.0, 0.0, 0.0};
  int nseg = snode.NumElement();
  DRT::Element** adjeles = snode.Elements();

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************
  CORE::LINALG::SerialDenseMatrix elens(6, nseg);
  MORTAR::Element* adjmrtrele = dynamic_cast<MORTAR::Element*>(adjeles[0]);

  // build element normal at current node
  // (we have to pass in the index i to be able to store the
  // normal and other information at the right place in elens)
  int i = 0;
  adjmrtrele->BuildNormalAtNode(snode.Id(), i, elens);

  // add (weighted) element normal to nodal normal n
  for (int j = 0; j < 3; ++j) slavebasednormal[j] += elens(j, 0) / elens(4, 0);

  // create unit normal vector
  const double length =
      sqrt(slavebasednormal[0] * slavebasednormal[0] + slavebasednormal[1] * slavebasednormal[1] +
           slavebasednormal[2] * slavebasednormal[2]);
  if (abs(length) < 1e-12)
  {
    dserror("Nodal normal length 0, node ID %i", snode.Id());
  }
  else
  {
    for (int j = 0; j < 3; ++j) slavebasednormal[j] /= length;
  }

  const double dotprod = -(normal[0] * slavebasednormal[0] + normal[1] * slavebasednormal[1] +
                           normal[2] * slavebasednormal[2]);

  if (dotprod < -1e-12)
  {
    // get the cpp normal
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];

    for (int j = 0; j < dim; ++j)
      for (CI p = glin[j].begin(); p != glin[j].end(); ++p)
        (normaltonodelin[j])[p->first] -= (p->second);
  }
  else
  {
    // linearization
    for (int j = 0; j < dim; ++j)
      for (CI p = glin[j].begin(); p != glin[j].end(); ++p)
        (normaltonodelin[j])[p->first] += (p->second);
  }


  return gdist;
}

/*----------------------------------------------------------------------*
 |  evaluate closest point normals                          farah 08/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateCPPNormals()
{
  // Build averaged normal field on physically smooth surface
  // loop over proc's master nodes of the interface
  // use row map and export to column map later
  for (int i = 0; i < MasterRowNodes()->NumMyElements(); ++i)
  {
    int gid = MasterRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* mrtrnode = dynamic_cast<Node*>(node);

    // build averaged normal at each master node
    mrtrnode->BuildAveragedNormal();

    // build tangent
    if (mrtrnode->IsOnEdge()) mrtrnode->BuildAveragedEdgeTangent();
  }

  // export nodal normals
  ExportMasterNodalNormals();

  // loop over slave nodes
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    MORTAR::Node* mrtrnode = dynamic_cast<MORTAR::Node*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) dserror("Node ownership inconsistency!");

    // vector with possible contacting master eles/nodes
    std::vector<MORTAR::Element*> meles;
    std::vector<MORTAR::Node*> mnodes;

    // fill vector with possibly contacting meles
    FindMEles(*mrtrnode, meles);

    // fallback solution if no mele is available
    if (meles.size() < 1)  // or !mrtrnode->IsOnCornerEdge())
    {
      Node* cnode = dynamic_cast<Node*>(mrtrnode);
      cnode->BuildAveragedNormal();
      continue;
    }


    // Here we have all found master elements for one slave node.
    // distance for cpp
    double normaltoline[3] = {0.0, 0.0, 0.0};
    std::vector<CORE::GEN::pairedvector<int, double>> normaltolineLin(3, 1);  // 1 dummy

    // Now, calculate distance between node and master line
    double dist = ComputeCPPNormal(*mrtrnode, meles, normaltoline, normaltolineLin);

    // if no projection was posible
    if (dist > 1e11)
    {
      Node* cnode = dynamic_cast<Node*>(mrtrnode);
      cnode->BuildAveragedNormal();
      continue;
    }

    // set the normal and its lineratization
    SetCPPNormal(*mrtrnode, normaltoline, normaltolineLin);
  }

  // export slave normals
  ExportNodalNormals();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  export master nodal normals (protected)                  farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ExportMasterNodalNormals()
{
  std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> triad;

  std::map<int, std::vector<int>> n_x_key;
  std::map<int, std::vector<int>> n_y_key;
  std::map<int, std::vector<int>> n_z_key;
  std::map<int, std::vector<int>> txi_x_key;
  std::map<int, std::vector<int>> txi_y_key;
  std::map<int, std::vector<int>> txi_z_key;
  std::map<int, std::vector<int>> teta_x_key;
  std::map<int, std::vector<int>> teta_y_key;
  std::map<int, std::vector<int>> teta_z_key;

  std::map<int, std::vector<double>> n_x_val;
  std::map<int, std::vector<double>> n_y_val;
  std::map<int, std::vector<double>> n_z_val;
  std::map<int, std::vector<double>> txi_x_val;
  std::map<int, std::vector<double>> txi_y_val;
  std::map<int, std::vector<double>> txi_z_val;
  std::map<int, std::vector<double>> teta_x_val;
  std::map<int, std::vector<double>> teta_y_val;
  std::map<int, std::vector<double>> teta_z_val;

  CORE::GEN::pairedvector<int, double>::iterator iter;

  const Teuchos::RCP<Epetra_Map> masternodes = CORE::LINALG::AllreduceEMap(*(mnoderowmap_));

  // build info on row map
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);

    // fill nodal matrix
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> loc =
        Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(3, 3));
    (*loc)(0, 0) = cnode->MoData().n()[0];
    (*loc)(1, 0) = cnode->MoData().n()[1];
    (*loc)(2, 0) = cnode->MoData().n()[2];
    (*loc)(0, 1) = cnode->Data().txi()[0];
    (*loc)(1, 1) = cnode->Data().txi()[1];
    (*loc)(2, 1) = cnode->Data().txi()[2];
    (*loc)(0, 2) = cnode->Data().teta()[0];
    (*loc)(1, 2) = cnode->Data().teta()[1];
    (*loc)(2, 2) = cnode->Data().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->Data().GetDerivN();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivtxi = cnode->Data().GetDerivTxi();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivteta = cnode->Data().GetDerivTeta();

    for (iter = derivn[0].begin(); iter != derivn[0].end(); ++iter)
    {
      n_x_key[gid].push_back(iter->first);
      n_x_val[gid].push_back(iter->second);
    }
    for (iter = derivn[1].begin(); iter != derivn[1].end(); ++iter)
    {
      n_y_key[gid].push_back(iter->first);
      n_y_val[gid].push_back(iter->second);
    }
    for (iter = derivn[2].begin(); iter != derivn[2].end(); ++iter)
    {
      n_z_key[gid].push_back(iter->first);
      n_z_val[gid].push_back(iter->second);
    }

    for (iter = derivtxi[0].begin(); iter != derivtxi[0].end(); ++iter)
    {
      txi_x_key[gid].push_back(iter->first);
      txi_x_val[gid].push_back(iter->second);
    }
    for (iter = derivtxi[1].begin(); iter != derivtxi[1].end(); ++iter)
    {
      txi_y_key[gid].push_back(iter->first);
      txi_y_val[gid].push_back(iter->second);
    }
    for (iter = derivtxi[2].begin(); iter != derivtxi[2].end(); ++iter)
    {
      txi_z_key[gid].push_back(iter->first);
      txi_z_val[gid].push_back(iter->second);
    }

    for (iter = derivteta[0].begin(); iter != derivteta[0].end(); ++iter)
    {
      teta_x_key[gid].push_back(iter->first);
      teta_x_val[gid].push_back(iter->second);
    }
    for (iter = derivteta[1].begin(); iter != derivteta[1].end(); ++iter)
    {
      teta_y_key[gid].push_back(iter->first);
      teta_y_val[gid].push_back(iter->second);
    }
    for (iter = derivteta[2].begin(); iter != derivteta[2].end(); ++iter)
    {
      teta_z_key[gid].push_back(iter->first);
      teta_z_val[gid].push_back(iter->second);
    }
  }

  // communicate from master node row to column map
  CORE::COMM::Exporter ex(*mnoderowmap_, *masternodes, Comm());
  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for (int i = 0; i < masternodes->NumMyElements(); ++i)
  {
    // only do something for ghosted nodes
    int gid = masternodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    int linsize = cnode->GetLinsize() + (int)(n_x_key[gid].size());

    if (cnode->Owner() == Comm().MyPID()) continue;

    // extract info
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> loc = triad[gid];
    cnode->MoData().n()[0] = (*loc)(0, 0);
    cnode->MoData().n()[1] = (*loc)(1, 0);
    cnode->MoData().n()[2] = (*loc)(2, 0);
    cnode->Data().txi()[0] = (*loc)(0, 1);
    cnode->Data().txi()[1] = (*loc)(1, 1);
    cnode->Data().txi()[2] = (*loc)(2, 1);
    cnode->Data().teta()[0] = (*loc)(0, 2);
    cnode->Data().teta()[1] = (*loc)(1, 2);
    cnode->Data().teta()[2] = (*loc)(2, 2);

    // extract derivative info
    std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->Data().GetDerivN();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivtxi = cnode->Data().GetDerivTxi();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivteta = cnode->Data().GetDerivTeta();

    for (int k = 0; k < (int)(derivn.size()); ++k) derivn[k].clear();
    derivn.resize(3, linsize);
    for (int k = 0; k < (int)(derivtxi.size()); ++k) derivtxi[k].clear();
    derivtxi.resize(3, linsize);
    for (int k = 0; k < (int)(derivteta.size()); ++k) derivteta[k].clear();
    derivteta.resize(3, linsize);

    cnode->Data().GetDerivN()[0].resize(linsize);
    cnode->Data().GetDerivN()[1].resize(linsize);
    cnode->Data().GetDerivN()[2].resize(linsize);

    cnode->Data().GetDerivTxi()[0].resize(linsize);
    cnode->Data().GetDerivTxi()[1].resize(linsize);
    cnode->Data().GetDerivTxi()[2].resize(linsize);

    cnode->Data().GetDerivTeta()[0].resize(linsize);
    cnode->Data().GetDerivTeta()[1].resize(linsize);
    cnode->Data().GetDerivTeta()[2].resize(linsize);

    for (int k = 0; k < (int)(n_x_key[gid].size()); ++k)
      (cnode->Data().GetDerivN()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k = 0; k < (int)(n_y_key[gid].size()); ++k)
      (cnode->Data().GetDerivN()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k = 0; k < (int)(n_z_key[gid].size()); ++k)
      (cnode->Data().GetDerivN()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k = 0; k < (int)(txi_x_key[gid].size()); ++k)
      (cnode->Data().GetDerivTxi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k = 0; k < (int)(txi_y_key[gid].size()); ++k)
      (cnode->Data().GetDerivTxi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k = 0; k < (int)(txi_z_key[gid].size()); ++k)
      (cnode->Data().GetDerivTxi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k = 0; k < (int)(teta_x_key[gid].size()); ++k)
      (cnode->Data().GetDerivTeta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k = 0; k < (int)(teta_y_key[gid].size()); ++k)
      (cnode->Data().GetDerivTeta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k = 0; k < (int)(teta_z_key[gid].size()); ++k)
      (cnode->Data().GetDerivTeta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate nodal normals (public)                          farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateAveragedNodalNormals()
{
  dserror("outdated function!");
  // safety
  if (smoothnodes_ == Teuchos::null) dserror("map of non smooth nodes is wrong!");

  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < smoothnodes_->NumMyElements(); ++i)
  {
    int gid = smoothnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* mrtrnode = dynamic_cast<Node*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calc scaling for hybrid formulation                     farah 09/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ComputeScaling()
{
  // ltl scaling only for 3D setting
  if (Dim() == 2) return;

  // compute ltl scaling for point and line contact
  ComputeScalingLTL();

  // bye
  return;
}


/*----------------------------------------------------------------------*
 |  calc scaling for hybrid formulation                     farah 01/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ComputeScalingLTL()
{
  std::cout << "ComputeScalingLTL" << std::endl;

  // define iterator for linerization
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;

  // angle in degree
  const double minAngle = InterfaceParams().get<double>("HYBRID_ANGLE_MIN");
  const double maxAngle = InterfaceParams().get<double>("HYBRID_ANGLE_MAX");

  // check
  if (minAngle < 0.0 or maxAngle < 0.0) dserror("invalid hybrid angle!");
  if (minAngle >= maxAngle) dserror("invalid hybrid angle!");

  // angle in rad
  const double alphaMin = minAngle * (M_PI / 180.0);  // min angle in rad
  const double alphaMax = maxAngle * (M_PI / 180.0);  // max angle in rad

  //======================================================
  // search slave edges and master edges and store them

  // guarantee uniquness of slave edges
  std::set<std::pair<int, int>> donebeforeS;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::Element>> lineElementsS;

    if (selement->Shape() == CORE::FE::CellType::quad4)
    {
      for (int j = 0; j < 4; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = selement->NodeIds()[3];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebeforeS.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebeforeS.find(actIDstw);

        // if not then create ele
        if (iter == donebeforeS.end() and itertw == donebeforeS.end())
        {
          // add to set of processed nodes
          donebeforeS.insert(actIDs);
          donebeforeS.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
              j, selement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<DRT::Node*, 2> nodes = {
              selement->Nodes()[nodeLIds[0]], selement->Nodes()[nodeLIds[1]]};
          lineEle->BuildNodalPointers(nodes.data());

          // init data container for dual shapes
          lineEle->InitializeDataContainer();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else
      dserror("LTL only for quad4!");

    // guarantee uniquness of master edges
    std::set<std::pair<int, int>> donebeforeM;

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::Element>> lineElementsM;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < selement->MoData().NumSearchElements(); ++k)
    {
      int gid2 = selement->MoData().SearchElements()[k];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      if (melement->Shape() == CORE::FE::CellType::quad4)
      {
        for (int j = 0; j < 4; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = melement->NodeIds()[3];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge = dynamic_cast<MORTAR::Node*>(melement->Nodes()[nodeLIds[0]])->IsOnEdge();
          bool node1Edge = dynamic_cast<MORTAR::Node*>(melement->Nodes()[nodeLIds[1]])->IsOnEdge();

          if (!node0Edge or !node1Edge) continue;

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebeforeM.find(actIDstw);

          // if not then create ele
          if (iter == donebeforeM.end() and itertw == donebeforeM.end())
          {
            // add to set of processed nodes
            donebeforeM.insert(actIDs);
            donebeforeM.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
                j, melement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<DRT::Node*, 2> nodes = {
                melement->Nodes()[nodeLIds[0]], melement->Nodes()[nodeLIds[1]]};
            lineEle->BuildNodalPointers(nodes.data());

            // init data container for dual shapes
            lineEle->InitializeDataContainer();

            // push back into vector
            lineElementsM.push_back(lineEle);
          }
        }  // end edge loop
      }
      else
        dserror("LTL only for quad4!");
    }  // end found mele loop

    // loop over slave edges
    for (int s = 0; s < (int)lineElementsS.size(); ++s)
    {
      double gR = 1e12;
      bool parallel = false;

      // loop over master edges
      for (int m = 0; m < (int)lineElementsM.size(); ++m)
      {
        // 1. search slave master pair with cpp point
        LineToLineCouplingPoint3d coup(
            *idiscret_, 3, InterfaceParams(), lineElementsS[s], lineElementsM[m]);

        parallel = coup.CheckParallelity();
        if (parallel)
        {
          continue;
          //          dserror("edges parallel");
          break;
        }
        // create empty points
        double sxi = 0.0;
        double mxi = 0.0;
        CORE::GEN::pairedvector<int, double> dsxi(
            3 * lineElementsM[m]->NumNode() + 3 * lineElementsS[s]->NumNode());
        CORE::GEN::pairedvector<int, double> dmxi(
            3 * lineElementsM[m]->NumNode() + 3 * lineElementsS[s]->NumNode());

        coup.LineIntersection(&sxi, &mxi, dsxi, dmxi);
        bool valid = coup.CheckIntersection(&sxi, &mxi);

        // no valid intersection
        if (!valid)
        {
          continue;
        }
        // found valid intersection
        else
        {
          std::cout << "VALID INTERSECTION" << std::endl;
          // 2. calc current angle for this element pair
          CORE::GEN::pairedvector<int, double> linAngle(
              100 + 3 * lineElementsM[m]->NumNode() + 3 * lineElementsS[s]->NumNode());
          gR = coup.CalcCurrentAngle(linAngle);

          // get vector of nodes
          std::vector<Node*> cnodes;
          for (int nn = 0; nn < lineElementsS[s]->NumNode(); ++nn)
            cnodes.push_back(dynamic_cast<Node*>(lineElementsS[s]->Nodes()[nn]));

          //======================================================
          // calc alpha:
          double alpha = 0.0;

          std::cout << "gR= " << gR << std::endl;

          // ltl line contact
          if (gR <= alphaMin)
          {
            std::cout << "LTS contact" << std::endl;
            alpha = 1.0;
          }
          // ltl hybrid
          else if (gR > alphaMin and gR < alphaMax)
          {
            std::cout << "current transition angle in degree = " << gR * (180.0 / M_PI)
                      << std::endl;
            alpha = 0.5 * (1.0 - cos(M_PI * (gR - alphaMax) / (alphaMin - alphaMax)));
            std::cout << "alpha = " << alpha << std::endl;
          }
          // ltl point contact
          else
          {
            std::cout << "LTL contact" << std::endl;
            alpha = 0.0;
          }

          // store to data container
          for (int nn = 0; nn < lineElementsS[s]->NumNode(); ++nn)
            cnodes[nn]->Data().GetAlphaN() = alpha;

          //======================================================
          // lin alpha:
          if (gR <= alphaMin)
          {
            // nothing: pure ltl line contact
          }
          // hybrid linearization:
          else if (gR > alphaMin and gR < alphaMax)
          {
            // clear old data
            for (int nn = 0; nn < lineElementsS[s]->NumNode(); ++nn)
              if (cnodes[nn]->Data().GetAlpha().size() == 0)
                cnodes[nn]->Data().GetAlpha().resize(1000);

            const double fac = (M_PI / (2.0 * (alphaMin - alphaMax))) *
                               sin(M_PI * (gR - alphaMax) / (alphaMin - alphaMax));
            for (CI p = linAngle.begin(); p != linAngle.end(); ++p)
            {
              for (int nn = 0; nn < lineElementsS[s]->NumNode(); ++nn)
                (cnodes[nn]->Data().GetAlpha())[p->first] = fac * (p->second);
            }
          }
          else
          {
            // nothing: pure ltl point contact
          }
          break;
        }  // valid projection
      }
    }
  }  // end slave loop



  //  //======================================================
  //  // loop over smooth slave nodes
  //  for (int i = 0; i < smoothnodes_->NumMyElements(); ++i)
  //  {
  //    int gid = smoothnodes_->GID(i);
  //    DRT::Node* node = idiscret_->gNode(gid);
  //    if (!node)
  //      dserror("Cannot find node with gid %", gid);
  //    Node* cnode = dynamic_cast<Node*>(node);
  //
  //    if (cnode->Owner() != Comm().MyPID())
  //      dserror("Node ownership inconsistency!");
  //
  //    // MORTAR
  //    cnode->Data().GetAlphaN() = -1.0;
  //  }
  //
  //  //======================================================
  //  // loop over non smooth slave nodes
  //  for (int i = 0; i < nonsmoothnodes_->NumMyElements(); ++i)
  //  {
  //    int gid = nonsmoothnodes_->GID(i);
  //    DRT::Node* node = idiscret_->gNode(gid);
  //    if (!node)
  //      dserror("Cannot find node with gid %", gid);
  //    Node* cnode = dynamic_cast<Node*>(node);
  //
  //    if (cnode->Owner() != Comm().MyPID())
  //      dserror("Node ownership inconsistency!");
  //
  //    // clear old data
  //    if(cnode->Data().GetAlpha().size() == 0)
  //      cnode->Data().GetAlpha().resize(1000);
  //
  //    // get cpp normal
  //    double normalcpp[3] = {0.0, 0.0, 0.0};
  //    normalcpp[0] = cnode->MoData().n()[0];
  //    normalcpp[1] = cnode->MoData().n()[1];
  //    normalcpp[2] = cnode->MoData().n()[2];
  //
  //    // NTS evtl...
  //    // calc element normal
  //    const int eles = cnode->NumElement();
  //    if(eles>2)
  //      dserror("element number to high!");
  //
  //    // calc smallest angle between cpp and nele
  //    double gR = 1e12;
  //    CORE::GEN::pairedvector<int,double> ddotcppfinal(1000);
  //
  //    for (int n = 0; n<eles; ++n)
  //    {
  //      // create element normal:
  //      int loc = 0;
  //      Element* seleaux = dynamic_cast<Element*>(cnode->Elements()[n]);
  //      CORE::LINALG::SerialDenseMatrix elens(6,1);
  //      seleaux->BuildNormalAtNode(cnode->Id(),loc,elens);
  //
  //      double normalsele[3] = {0.0, 0.0, 0.0};
  //      normalsele[0] = elens(0,0) / elens(4,0);
  //      normalsele[1] = elens(1,0) / elens(4,0);
  //      normalsele[2] = elens(2,0) / elens(4,0);
  //
  //      // calculate angle between cpp and elenormal
  //      double dotcpp = normalsele[0] * normalcpp[0] + normalsele[1] * normalcpp[1] +
  //      normalsele[2] * normalcpp[2]; double Rl = acos(dotcpp); // current angle in rad
  //
  //      if(Rl<0.0)
  //        dserror("angle less than 0.0");
  //
  //      //======================================================
  //      // lin angle
  //      // lin Rl
  //      double fac = (-1.0/(sqrt(1.0-dotcpp*dotcpp)));
  //
  //      // init lin
  //      CORE::GEN::pairedvector<int,double> ddotcpp(1000);
  //      std::vector<CORE::GEN::pairedvector<int,double> > dnele(3,1000);
  //
  //      // deriv unit normal!
  //      seleaux->DerivNormalAtNode(cnode->Id(),loc,elens,dnele);
  //
  //      // linearization of cpp normal
  //      for (CI
  //      p=cnode->Data().GetDerivN()[0].begin();p!=cnode->Data().GetDerivN()[0].end();++p)
  //        (ddotcpp)[p->first] += (p->second) * normalsele[0];
  //      for (CI
  //      p=cnode->Data().GetDerivN()[1].begin();p!=cnode->Data().GetDerivN()[1].end();++p)
  //        (ddotcpp)[p->first] += (p->second) * normalsele[1];
  //      for (CI
  //      p=cnode->Data().GetDerivN()[2].begin();p!=cnode->Data().GetDerivN()[2].end();++p)
  //        (ddotcpp)[p->first] += (p->second) * normalsele[2];
  //
  //      // linearization of ele normal
  //      for (CI p=dnele[0].begin();p!=dnele[0].end();++p)
  //        (ddotcpp)[p->first] += (p->second) * normalcpp[0];
  //      for (CI p=dnele[1].begin();p!=dnele[1].end();++p)
  //        (ddotcpp)[p->first] += (p->second) * normalcpp[1];
  //      for (CI p=dnele[2].begin();p!=dnele[2].end();++p)
  //        (ddotcpp)[p->first] += (p->second) * normalcpp[2];
  //
  //      // multiply with fac: -asin(dotcpp)
  //      for (CI p=ddotcpp.begin();p!=ddotcpp.end();++p)
  //        (ddotcpp)[p->first] = fac * (p->second);
  //
  //
  //      // safe smaller angle and linearization
  //      if(Rl<gR)
  //      {
  //        gR=Rl;
  //        ddotcppfinal = ddotcpp;
  //      }
  //    }// end ele loop
  //
  //
  //    //======================================================
  //    // calc alpha:
  //    double alpha = 0.0;
  //    if(gR<=alphaMin)
  //    {
  //      alpha = 1.0;
  //    }
  //    else if(gR>alphaMin and gR<alphaMax)
  //    {
  //      std::cout << "current transition angle in degree = " << gR*(180.0/3.141592) << std::endl;
  //      alpha=0.5*(1.0 - cos(3.141592*(gR-alphaMax)/(alphaMin-alphaMax)));
  //      std::cout << "alpha = " << alpha << std::endl;
  //    }
  //    else
  //    {
  //      alpha = 0.0;
  //    }
  //
  //    // store to data container
  //    cnode->Data().GetAlphaN() = alpha;
  //
  //    //======================================================
  //    // lin alpha:
  //    if(gR<=alphaMin)
  //    {
  //      // nothing: pure mortar
  //    }
  //    // hybrid linearization:
  //    else if (gR>alphaMin and gR<alphaMax)
  //    {
  //      const double fac = (3.141592/(2.0*(alphaMin-alphaMax))) *
  //      sin(3.141592*(gR-alphaMax)/(alphaMin-alphaMax)); for (CI
  //      p=ddotcppfinal.begin();p!=ddotcppfinal.end();++p)
  //        (cnode->Data().GetAlpha())[p->first] = fac * (p->second);
  //    }
  //    else
  //    {
  //      // nothing: pure nts
  //    }
  //  }// end node loop

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ScaleNormals()
{
  dserror("outdated!");

  if (Dim() == 2)
    ScaleNormals2D();
  else if (Dim() == 3)
    ScaleNormals3D();
  else
    dserror("Wrong dimension!");
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ScaleNormals2D()
{
  // define iterator for paired vector
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;

  //======================================================
  // loop over non smooth slave nodes
  for (int i = 0; i < nonsmoothnodes_->NumMyElements(); ++i)
  {
    int gid = nonsmoothnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    if (cnode->Owner() != Comm().MyPID()) dserror("Node ownership inconsistency!");

    // get cpp normal
    std::array<double, 3> normalcpp = {0.0, 0.0, 0.0};

    // get cpp normal lin.
    std::vector<CORE::GEN::pairedvector<int, double>> dcppaux(3, 1000);

    dcppaux = cnode->Data().GetDerivN();
    normalcpp[0] = cnode->MoData().n()[0];
    normalcpp[1] = cnode->MoData().n()[1];
    normalcpp[2] = cnode->MoData().n()[2];

    if (cnode->Active())
    {
      std::cout << "cnode->Data().GetAlphaN()] = " << cnode->Data().GetAlphaN() << std::endl;
      std::cout << "normalcpp[0] = " << normalcpp[0] << std::endl;
      std::cout << "normalcpp[1] = " << normalcpp[1] << std::endl;
      std::cout << "normalcpp[2] = " << normalcpp[2] << std::endl;
    }


    // calc element normal
    const int eles = cnode->NumElement();
    if (eles > 2) dserror("element number to high!");

    double gR = 1e12;
    int localid = -1;

    for (int n = 0; n < eles; ++n)
    {
      // create element normal:
      int loc = 0;
      Element* seleaux = dynamic_cast<Element*>(cnode->Elements()[n]);
      CORE::LINALG::SerialDenseMatrix elens(6, 1);
      seleaux->BuildNormalAtNode(cnode->Id(), loc, elens);

      std::array<double, 3> normalsele = {0.0, 0.0, 0.0};
      normalsele[0] = elens(0, 0) / elens(4, 0);
      normalsele[1] = elens(1, 0) / elens(4, 0);
      normalsele[2] = elens(2, 0) / elens(4, 0);

      // calculate angle between cpp and elenormal
      double dotcpp = normalsele[0] * normalcpp[0] + normalsele[1] * normalcpp[1] +
                      normalsele[2] * normalcpp[2];
      double Rl = acos(dotcpp);  // current angle in rad

      if (cnode->Active())
      {
        std::cout << "normalsele[0] = " << normalsele[0] << std::endl;
        std::cout << "normalsele[1] = " << normalsele[1] << std::endl;
        std::cout << "normalsele[2] = " << normalsele[2] << std::endl;

        std::cout << "dotcpp= " << dotcpp << std::endl;
        std::cout << "Rl= " << Rl << std::endl;
      }


      if (Rl < 0.0) dserror("angle less than 0.0");



      // safe smaller angle and linearization
      if (Rl < gR)
      {
        gR = Rl;
        localid = n;
      }
    }  // end ele loop

    // create new normals:
    const double alpha = cnode->Data().GetAlphaN();

    // create element normal:
    int loc = 0;
    Element* seleaux = dynamic_cast<Element*>(cnode->Elements()[localid]);
    CORE::LINALG::SerialDenseMatrix elens(6, 1);
    seleaux->BuildNormalAtNode(cnode->Id(), loc, elens);

    std::array<double, 3> normalsele = {0.0, 0.0, 0.0};
    normalsele[0] = elens(0, 0) / elens(4, 0);
    normalsele[1] = elens(1, 0) / elens(4, 0);
    normalsele[2] = elens(2, 0) / elens(4, 0);

    cnode->MoData().n()[0] = alpha * normalsele[0] + (1.0 - alpha) * normalcpp[0];
    cnode->MoData().n()[1] = alpha * normalsele[1] + (1.0 - alpha) * normalcpp[1];
    cnode->MoData().n()[2] = alpha * normalsele[2] + (1.0 - alpha) * normalcpp[2];


    // deriv unit normal!
    std::vector<CORE::GEN::pairedvector<int, double>> dnele(3, 1000);
    seleaux->DerivNormalAtNode(cnode->Id(), loc, elens, dnele);

    cnode->Data().GetDerivN()[0].clear();
    cnode->Data().GetDerivN()[1].clear();
    cnode->Data().GetDerivN()[2].clear();

    // lin ele normal * alpha
    for (CI p = dnele[0].begin(); p != dnele[0].end(); ++p)
      (cnode->Data().GetDerivN()[0])[p->first] += alpha * (p->second);
    for (CI p = dnele[1].begin(); p != dnele[1].end(); ++p)
      (cnode->Data().GetDerivN()[1])[p->first] += alpha * (p->second);
    for (CI p = dnele[2].begin(); p != dnele[2].end(); ++p)
      (cnode->Data().GetDerivN()[2])[p->first] += alpha * (p->second);

    // lin cpp normal * alpha
    for (CI p = dcppaux[0].begin(); p != dcppaux[0].end(); ++p)
      (cnode->Data().GetDerivN()[0])[p->first] += (1.0 - alpha) * (p->second);
    for (CI p = dcppaux[1].begin(); p != dcppaux[1].end(); ++p)
      (cnode->Data().GetDerivN()[1])[p->first] += (1.0 - alpha) * (p->second);
    for (CI p = dcppaux[2].begin(); p != dcppaux[2].end(); ++p)
      (cnode->Data().GetDerivN()[2])[p->first] += (1.0 - alpha) * (p->second);

    // lin alpha
    for (CI p = cnode->Data().GetAlpha().begin(); p != cnode->Data().GetAlpha().end(); ++p)
    {
      (cnode->Data().GetDerivN()[0])[p->first] += (p->second) * (normalsele[0] - normalcpp[0]);
      (cnode->Data().GetDerivN()[1])[p->first] += (p->second) * (normalsele[1] - normalcpp[1]);
      (cnode->Data().GetDerivN()[2])[p->first] += (p->second) * (normalsele[2] - normalcpp[2]);
    }

    // 2D Tangent!
    if (cnode->NumDof() == 2)
    {
      // simple definition for txi
      cnode->Data().txi()[0] = -cnode->MoData().n()[1];
      cnode->Data().txi()[1] = cnode->MoData().n()[0];
      cnode->Data().txi()[2] = 0.0;

      // teta is z-axis
      cnode->Data().teta()[0] = 0.0;
      cnode->Data().teta()[1] = 0.0;
      cnode->Data().teta()[2] = 1.0;

      cnode->Data().GetDerivTxi()[0].clear();
      cnode->Data().GetDerivTxi()[1].clear();
      cnode->Data().GetDerivTxi()[2].clear();


      for (CI p = cnode->Data().GetDerivN()[1].begin(); p != cnode->Data().GetDerivN()[1].end();
           ++p)
        (cnode->Data().GetDerivTxi()[0])[p->first] -= (p->second);
      for (CI p = cnode->Data().GetDerivN()[0].begin(); p != cnode->Data().GetDerivN()[0].end();
           ++p)
        (cnode->Data().GetDerivTxi()[1])[p->first] += (p->second);
    }
    else
      dserror("only 2D");



    double length = sqrt(cnode->MoData().n()[0] * cnode->MoData().n()[0] +
                         cnode->MoData().n()[1] * cnode->MoData().n()[1] +
                         cnode->MoData().n()[2] * cnode->MoData().n()[2]);

    if (cnode->Active())
    {
      std::cout << "cnode->MoData().n()[0]= " << cnode->MoData().n()[0] << std::endl;
      std::cout << "cnode->MoData().n()[1]= " << cnode->MoData().n()[1] << std::endl;
      std::cout << "cnode->MoData().n()[2]= " << cnode->MoData().n()[2] << std::endl;
      std::cout << "length= " << length << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ScaleNormals3D() { dserror("not yet implemented!"); }


/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 08/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::ComputeCPPNormal2D(MORTAR::Node& mrtrnode,
    std::vector<MORTAR::Element*> meles, double* normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& normaltolineLin)
{
  // define tolerance
  const double tol = 1e-8;
  const double validAngle = 20;

  // distance between node and surface
  double gdist = 1e12;  // distance
  std::array<double, 3> gnormal = {0.0, 0.0, 0.0};
  std::vector<CORE::GEN::pairedvector<int, double>> glin(3, 1);  // 1 dummy
  std::set<int> donebeforeMasterCorner;

  bool nodeOnNode = false;     // flag for node on node (corner on corner) setting
  bool pathdependent = false;  // flag if we have to check path from last converged check
  std::array<double, 3> vect = {0.0, 0.0, 0.0};  // patch

  // calc trajectory and node-to-node distance for corner nodes
  if (mrtrnode.IsOnCorner())
  {
    pathdependent = true;
    CONTACT::Node& coNode = dynamic_cast<CONTACT::Node&>(mrtrnode);
    if (coNode.Active()) pathdependent = false;

    // calculate path
    // check trajectory from considered node
    std::array<double, 3> Posn = {0.0, 0.0, 0.0};
    std::array<double, 3> Posnp = {0.0, 0.0, 0.0};
    Posn[0] = mrtrnode.X()[0] + mrtrnode.uold()[0];
    Posn[1] = mrtrnode.X()[1] + mrtrnode.uold()[1];
    Posn[2] = mrtrnode.X()[2] + mrtrnode.uold()[2];
    Posnp[0] = mrtrnode.xspatial()[0];
    Posnp[1] = mrtrnode.xspatial()[1];
    Posnp[2] = mrtrnode.xspatial()[2];

    vect[0] = Posnp[0] - Posn[0];
    vect[1] = Posnp[1] - Posn[1];
    vect[2] = Posnp[2] - Posn[2];

    double lvec = sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
    if (lvec < 1e-12)
    {
      pathdependent = false;
    }
    else
    {
      vect[0] /= lvec;
      vect[1] /= lvec;
      vect[2] /= lvec;
    }

    // Compute normal to node:
    // loop over found eles
    for (size_t ele = 0; ele < meles.size(); ++ele)
    {
      // get linsize
      int linsize = 0;
      for (int i = 0; i < meles[ele]->NumNode(); ++i)
      {
        DRT::Node* node = meles[ele]->Nodes()[i];
        if (!node) dserror("Cannot find master node");
        CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
        linsize += mnode->GetLinsize();

        // if master node is also corner node
        if (mnode->IsOnCorner())
        {
          std::set<int>::iterator iter = donebeforeMasterCorner.find(mnode->Id());

          // if not then create ele
          if (iter != donebeforeMasterCorner.end()) continue;

          // set master corner node id
          donebeforeMasterCorner.insert(mnode->Id());

          // aux variables
          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<CORE::GEN::pairedvector<int, double>> auxlin(
              3, linsize + 1 + meles[ele]->NumNode());

          // compute distance between corners
          dist = ComputeNormalNodeToNode(mrtrnode, *mnode, auxnormal, auxlin);

          // if nodes lying on each other
          if (abs(dist) < 1e-12)
          {
            nodeOnNode = true;
            continue;
          }

          // angle between trajectory and normal
          if (pathdependent)
          {
            double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1]);
            if (auxl < 1e-12) continue;
            double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl);

            angle = 180 * (angle / 3.14159265359);
            if (abs(angle) > validAngle) continue;
          }

          // get closest valid distance
          if (abs(dist) < abs(gdist))
          {
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }  // end node loop
    }    // end element loop
  }

  // loop over found eles
  for (size_t ele = 0; ele < meles.size(); ++ele)
  {
    bool cornerele = false;
    // get linsize
    int linsize = 0;
    for (int i = 0; i < meles[ele]->NumNode(); ++i)
    {
      DRT::Node* node = meles[ele]->Nodes()[i];
      if (!node) dserror("Cannot find master node");
      CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
      linsize += mnode->GetLinsize();

      if (mnode->IsOnCorner()) cornerele = true;
    }

    double xi[2] = {0.0, 0.0};
    double dist = 1e12;
    double auxnormal[3] = {0.0, 0.0, 0.0};
    std::vector<CORE::GEN::pairedvector<int, double>> auxlin(
        3, linsize + 1 + meles[ele]->NumNode());

    // check for nonsmooth mele
    if (cornerele and !nodeOnNode)
    {
      // perform CPP to find normals based on element normals
      MORTAR::Projector::Impl(*meles[ele])
          ->ProjectSNodeByMNormalLin(mrtrnode, *meles[ele], xi, auxnormal, dist, auxlin);
    }
    // compute normal with averaged nodal normal field from master surface
    else
    {
      // perform CPP to find normals based on averaged nodal normal field
      MORTAR::Projector::Impl(*meles[ele])
          ->ProjectSNodeByMNodalNormalLin(mrtrnode, *meles[ele], xi, auxnormal, dist, auxlin);
    }

    // check if found parameter space coordinate is within element domain
    if (meles[ele]->Shape() == CORE::FE::CellType::line2 or
        meles[ele]->Shape() == CORE::FE::CellType::line3)
    {
      if (-1.0 - tol > xi[0] or xi[0] > 1.0 + tol) continue;
    }
    else
    {
      dserror("Unknown ele type!");
    }

    // angle between trajectory and normal
    if (pathdependent)
    {
      double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1]);
      if (auxl < 1e-12) continue;
      double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl);

      angle = 180 * (angle / 3.14159265359);
      if (abs(angle) > validAngle) continue;
    }

    // get closest valid distance
    if (abs(dist) < abs(gdist))
    {
      gdist = dist;
      gnormal[0] = auxnormal[0];
      gnormal[1] = auxnormal[1];
      gnormal[2] = auxnormal[2];
      glin = auxlin;
    }
  }  // end mele loop

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];
  normaltolineLin = glin;

  return gdist;
}

/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 08/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::ComputeCPPNormal3D(MORTAR::Node& mrtrnode,
    std::vector<MORTAR::Element*> meles, double* normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& normaltolineLin)
{
  // define tolerance
  const double tol = 1e-8;

  // distance between node and surface
  bool pathdependent = true;
  const double validAngle = 5.0;
  double gdist = 1e12;  // distance
  std::array<double, 3> gnormal = {0.0, 0.0, 0.0};
  std::vector<CORE::GEN::pairedvector<int, double>> glin(3, 1);  // 1 dummy

  //******************************************************
  //             CALC TRAJECTORY
  //******************************************************
  std::array<double, 3> Posn = {0.0, 0.0, 0.0};
  std::array<double, 3> Posnp = {0.0, 0.0, 0.0};
  std::array<double, 3> vect = {0.0, 0.0, 0.0};

  Posn[0] = mrtrnode.X()[0] + mrtrnode.uold()[0];
  Posn[1] = mrtrnode.X()[1] + mrtrnode.uold()[1];
  Posn[2] = mrtrnode.X()[2] + mrtrnode.uold()[2];
  Posnp[0] = mrtrnode.xspatial()[0];
  Posnp[1] = mrtrnode.xspatial()[1];
  Posnp[2] = mrtrnode.xspatial()[2];

  vect[0] = Posnp[0] - Posn[0];
  vect[1] = Posnp[1] - Posn[1];
  vect[2] = Posnp[2] - Posn[2];

  double lvec = sqrt(vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]);
  if (lvec < 1e-12)
  {
    pathdependent = false;
  }
  else
  {
    vect[0] /= lvec;
    vect[1] /= lvec;
    vect[2] /= lvec;
  }

  if (dynamic_cast<CONTACT::Node&>(mrtrnode).Active()) pathdependent = false;

  //******************************************************
  //             COMPUTE NORMAL TO SURFACE
  //******************************************************
  // loop over found eles for all geometrical nodes
  for (size_t ele = 0; ele < meles.size(); ++ele)
  {
    double xi[2] = {0.0, 0.0};
    double dist = 1e12;
    double auxnormal[3] = {0.0, 0.0, 0.0};
    std::vector<CORE::GEN::pairedvector<int, double>> auxlin(3, 1000);

    // perform CPP to find normals
    bool success =
        MORTAR::Projector::Impl(*meles[ele])
            ->ProjectSNodeByMNodalNormalLin(mrtrnode, *meles[ele], xi, auxnormal, dist, auxlin);

    // newton not converged
    if (!success) continue;

    // check if found parameter space coordinate is within element domain
    if (meles[ele]->Shape() == CORE::FE::CellType::quad4 or
        meles[ele]->Shape() == CORE::FE::CellType::quad8 or
        meles[ele]->Shape() == CORE::FE::CellType::quad9)
    {
      if (-1.0 - tol > xi[0] or xi[0] > 1.0 + tol or -1.0 - tol > xi[1] or xi[1] > 1.0 + tol)
        continue;
    }
    else if (meles[ele]->Shape() == CORE::FE::CellType::tri3 or
             meles[ele]->Shape() == CORE::FE::CellType::tri6)
    {
      if (xi[0] < 0.0 - tol or xi[1] < 0.0 - tol or xi[0] > 1.0 + tol or xi[1] > 1.0 + tol or
          xi[0] + xi[1] > 1.0 + 2 * tol)
        continue;
    }
    else
    {
      dserror("Unknown ele type!");
    }

    // angle between trajectory and normal
    if (pathdependent)
    {
      double auxl = sqrt(
          auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] + auxnormal[2] * auxnormal[2]);
      if (auxl < 1e-12) continue;
      double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl -
                          vect[2] * auxnormal[2] / auxl);

      angle = 180 * (angle / 3.14159265359);
      if (abs(angle) > validAngle) continue;
    }

    if (dist < gdist)
    {
      gdist = dist;
      gnormal[0] = auxnormal[0];
      gnormal[1] = auxnormal[1];
      gnormal[2] = auxnormal[2];
      glin = auxlin;
    }
  }  // end mele loop
  //******************************************************
  //             COMPUTE NORMAL TO LINE
  //******************************************************
  if (mrtrnode.IsOnCornerEdge())  // only for edge or corner nodes possible
  {
    // guarantee uniquness
    std::set<std::pair<int, int>> donebefore;

    // calc
    for (size_t ele = 0; ele < meles.size(); ++ele)
    {
      // loop over master edges -> match node number for quad4
      for (int j = 0; j < meles[ele]->NumNode(); ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (meles[ele]->Shape() == CORE::FE::CellType::quad4)
        {
          if (j == 0)
          {
            nodeIds[0] = meles[ele]->NodeIds()[0];
            nodeIds[1] = meles[ele]->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = meles[ele]->NodeIds()[1];
            nodeIds[1] = meles[ele]->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = meles[ele]->NodeIds()[2];
            nodeIds[1] = meles[ele]->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = meles[ele]->NodeIds()[3];
            nodeIds[1] = meles[ele]->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            dserror("loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge =
            dynamic_cast<MORTAR::Node*>(meles[ele]->Nodes()[nodeLIds[0]])->IsOnCornerEdge();
        bool node1Edge =
            dynamic_cast<MORTAR::Node*>(meles[ele]->Nodes()[nodeLIds[1]])->IsOnCornerEdge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

        // if not then create ele
        if (iter == donebefore.end() and itertw == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(actIDs);
          donebefore.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
              j, meles[ele]->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          DRT::Node* nodes[2] = {
              meles[ele]->Nodes()[nodeLIds[0]], meles[ele]->Nodes()[nodeLIds[1]]};
          lineEle->BuildNodalPointers(nodes);

          // init data container for dual shapes
          lineEle->InitializeDataContainer();

          // call cpp function for edge to edge

          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<CORE::GEN::pairedvector<int, double>> auxlin(
              3, 100 + 1 + meles[ele]->NumNode());

          // compute distance between node and edge
          dist = ComputeNormalNodeToEdge(mrtrnode, *lineEle, auxnormal, auxlin);

          // angle between trajectory and normal
          if (pathdependent)
          {
            double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] +
                               auxnormal[2] * auxnormal[2]);
            if (auxl < 1e-12) continue;
            double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl -
                                vect[2] * auxnormal[2] / auxl);

            angle = 180 * (angle / 3.14159265359);
            if (abs(angle) > validAngle) continue;
          }

          if (dist <= gdist + tol)
          {
            std::cout << "CLOSE TO EDGE!!!" << std::endl;
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }  // end mele node loop
    }
  }

  //******************************************************
  //             COMPUTE NORMAL TO NODE
  //******************************************************
  if (mrtrnode.IsOnCorner())  // only for corner nodes possible
  {
    std::set<int> donebeforeMasterCorner;

    for (size_t ele = 0; ele < meles.size(); ++ele)
    {
      // get linsize
      int linsize = 0;
      for (int i = 0; i < meles[ele]->NumNode(); ++i)
      {
        DRT::Node* node = meles[ele]->Nodes()[i];
        if (!node) dserror("Cannot find master node");
        CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
        linsize += mnode->GetLinsize();

        // if master node is also corner node
        if (mnode->IsOnCorner())
        {
          std::set<int>::iterator iter = donebeforeMasterCorner.find(mnode->Id());

          // if not then create ele
          if (iter != donebeforeMasterCorner.end()) continue;

          // set master corner node id
          donebeforeMasterCorner.insert(mnode->Id());

          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<CORE::GEN::pairedvector<int, double>> auxlin(
              3, linsize + 1 + meles[ele]->NumNode());

          // compute distance between corners
          dist = ComputeNormalNodeToNode(mrtrnode, *mnode, auxnormal, auxlin);

          // angle between trajectory and normal
          if (pathdependent)
          {
            double auxl = sqrt(auxnormal[0] * auxnormal[0] + auxnormal[1] * auxnormal[1] +
                               auxnormal[2] * auxnormal[2]);
            if (auxl < 1e-12) continue;
            double angle = acos(-vect[0] * auxnormal[0] / auxl - vect[1] * auxnormal[1] / auxl -
                                vect[2] * auxnormal[2] / auxl);

            angle = 180 * (angle / 3.14159265359);
            if (abs(angle) > validAngle) continue;
          }

          if (dist < gdist)
          {
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }
    }
  }

  //******************************************************
  //             FINAL STORAGE
  //******************************************************

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];
  normaltolineLin = glin;

  // bye bye
  return gdist;
}

/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 05/16 |
 *----------------------------------------------------------------------*/
double CONTACT::Interface::ComputeCPPNormal(MORTAR::Node& mrtrnode,
    std::vector<MORTAR::Element*> meles, double* normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& normaltolineLin)
{
  // define distance
  double gdist = 1e12;

  //===================================================================
  //===================================================================
  //                           2D case
  //===================================================================
  //===================================================================
  if (Dim() == 2)
  {
    gdist = ComputeCPPNormal2D(mrtrnode, meles, normal, normaltolineLin);
  }
  //===================================================================
  //===================================================================
  //                           3D case
  //===================================================================
  //===================================================================
  else if (Dim() == 3)
  {
    gdist = ComputeCPPNormal3D(mrtrnode, meles, normal, normaltolineLin);
  }
  //===================================================================
  //===================================================================
  //                           Invalid
  //===================================================================
  //===================================================================
  else
  {
    dserror("invalid dimension!");
  }

  // return distance
  return gdist;
}

/*----------------------------------------------------------------------*
 |  set cpp normal                                           farah 01/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::SetCPPNormal(MORTAR::Node& snode, double* normal,
    std::vector<CORE::GEN::pairedvector<int, double>>& normallin)
{
  Node& cnode = dynamic_cast<Node&>(snode);

  const double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  if (length < 1e-12) dserror("normal length is zero!!!");

  // negative sign because it is a master normal!
  cnode.MoData().n()[0] = -normal[0] / length;
  cnode.MoData().n()[1] = -normal[1] / length;
  cnode.MoData().n()[2] = -normal[2] / length;

  //  if (cnode.IsOnEdge())
  //    std::cout << "normal =  " << cnode.MoData().n()[0] << "  "<< cnode.MoData().n()[1] << "  "<<
  //    cnode.MoData().n()[2] << std::endl;

  // prepare nodal storage maps for derivative
  if ((int)cnode.Data().GetDerivN().size() == 0)
    cnode.Data().GetDerivN().resize(3, normallin[0].size() * 3);
  if ((int)cnode.Data().GetDerivTxi().size() == 0)
    cnode.Data().GetDerivTxi().resize(3, normallin[0].size() * 3);
  if ((int)cnode.Data().GetDerivTeta().size() == 0)
    cnode.Data().GetDerivTeta().resize(3, normallin[0].size() * 3);

  // init tangent length
  double ltxi = -1.0;

  //------------------------------------------------------------------
  // 2D Tangent!
  if (cnode.NumDof() == 2)
  {
    // simple definition for txi
    cnode.Data().txi()[0] = -cnode.MoData().n()[1];
    cnode.Data().txi()[1] = cnode.MoData().n()[0];
    cnode.Data().txi()[2] = 0.0;

    // teta is z-axis
    cnode.Data().teta()[0] = 0.0;
    cnode.Data().teta()[1] = 0.0;
    cnode.Data().teta()[2] = 1.0;
  }
  // 3D Tangent!
  else
  {
    if (abs(cnode.MoData().n()[0]) > 1.0e-6 || abs(cnode.MoData().n()[1]) > 1.0e-6)
    {
      cnode.Data().txi()[0] = -cnode.MoData().n()[1];
      cnode.Data().txi()[1] = cnode.MoData().n()[0];
      cnode.Data().txi()[2] = 0.0;
    }
    else
    {
      cnode.Data().txi()[0] = 0.0;
      cnode.Data().txi()[1] = -cnode.MoData().n()[2];
      cnode.Data().txi()[2] = cnode.MoData().n()[1];
    }

    ltxi = sqrt(cnode.Data().txi()[0] * cnode.Data().txi()[0] +
                cnode.Data().txi()[1] * cnode.Data().txi()[1] +
                cnode.Data().txi()[2] * cnode.Data().txi()[2]);
    if (ltxi < 1e-12) dserror("tangent txi length is zero!!!");
    for (int j = 0; j < 3; ++j) cnode.Data().txi()[j] /= ltxi;

    // teta follows from corkscrew rule (teta = n x txi)
    cnode.Data().teta()[0] = cnode.MoData().n()[1] * cnode.Data().txi()[2] -
                             cnode.MoData().n()[2] * cnode.Data().txi()[1];
    cnode.Data().teta()[1] = cnode.MoData().n()[2] * cnode.Data().txi()[0] -
                             cnode.MoData().n()[0] * cnode.Data().txi()[2];
    cnode.Data().teta()[2] = cnode.MoData().n()[0] * cnode.Data().txi()[1] -
                             cnode.MoData().n()[1] * cnode.Data().txi()[0];
  }


  //------------------------------------------------------------------
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;

  for (CI p = normallin[0].begin(); p != normallin[0].end(); ++p)
    (cnode.Data().GetDerivN()[0])[p->first] -= (p->second);
  for (CI p = normallin[1].begin(); p != normallin[1].end(); ++p)
    (cnode.Data().GetDerivN()[1])[p->first] -= (p->second);
  for (CI p = normallin[2].begin(); p != normallin[2].end(); ++p)
    (cnode.Data().GetDerivN()[2])[p->first] -= (p->second);

  // normalize directional derivative
  // (length differs for weighted/unweighted case bot not the procedure!)
  // (be careful with reference / copy of derivative maps!)
  typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
  CORE::GEN::pairedvector<int, double>& derivnx = cnode.Data().GetDerivN()[0];
  CORE::GEN::pairedvector<int, double>& derivny = cnode.Data().GetDerivN()[1];
  CORE::GEN::pairedvector<int, double>& derivnz = cnode.Data().GetDerivN()[2];
  CORE::GEN::pairedvector<int, double> cderivnx = cnode.Data().GetDerivN()[0];
  CORE::GEN::pairedvector<int, double> cderivny = cnode.Data().GetDerivN()[1];
  CORE::GEN::pairedvector<int, double> cderivnz = cnode.Data().GetDerivN()[2];
  const double nxnx = cnode.MoData().n()[0] * cnode.MoData().n()[0];
  const double nxny = cnode.MoData().n()[0] * cnode.MoData().n()[1];
  const double nxnz = cnode.MoData().n()[0] * cnode.MoData().n()[2];
  const double nyny = cnode.MoData().n()[1] * cnode.MoData().n()[1];
  const double nynz = cnode.MoData().n()[1] * cnode.MoData().n()[2];
  const double nznz = cnode.MoData().n()[2] * cnode.MoData().n()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  std::vector<int> allkeysn;
  for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }
  for (CI p = derivny.begin(); p != derivny.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }
  for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }

  // normalize x-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivnx[allkeysn[j]];
    derivnx[allkeysn[j]] =
        (val - nxnx * val - nxny * cderivny[allkeysn[j]] - nxnz * cderivnz[allkeysn[j]]) / length;
  }

  // normalize y-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivny[allkeysn[j]];
    derivny[allkeysn[j]] =
        (val - nxny * cderivnx[allkeysn[j]] - nyny * val - nynz * cderivnz[allkeysn[j]]) / length;
  }

  // normalize z-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivnz[allkeysn[j]];
    derivnz[allkeysn[j]] =
        (val - nxnz * cderivnx[allkeysn[j]] - nynz * cderivny[allkeysn[j]] - nznz * val) / length;
  }

  //------------------------------------------------------------------
  // 2D Tangent!
  if (cnode.NumDof() == 2)
  {
    for (CI p = cnode.Data().GetDerivN()[1].begin(); p != cnode.Data().GetDerivN()[1].end(); ++p)
      (cnode.Data().GetDerivTxi()[0])[p->first] -= (p->second);
    for (CI p = cnode.Data().GetDerivN()[0].begin(); p != cnode.Data().GetDerivN()[0].end(); ++p)
      (cnode.Data().GetDerivTxi()[1])[p->first] += (p->second);
  }
  // 3D Tangent!
  else
  {
    // unnormalized tangent derivative txi
    // use definitions for txi from BuildAveragedNormal()
    if (abs(cnode.MoData().n()[0]) > 1.0e-6 || abs(cnode.MoData().n()[1]) > 1.0e-6)
    {
      CORE::GEN::pairedvector<int, double>& derivtxix = cnode.Data().GetDerivTxi()[0];
      CORE::GEN::pairedvector<int, double>& derivtxiy = cnode.Data().GetDerivTxi()[1];

      for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxix[p->first] -= (p->second);

      for (CI p = derivnx.begin(); p != derivnx.end(); ++p) derivtxiy[p->first] += (p->second);
    }
    else
    {
      CORE::GEN::pairedvector<int, double>& derivtxiy = cnode.Data().GetDerivTxi()[1];
      CORE::GEN::pairedvector<int, double>& derivtxiz = cnode.Data().GetDerivTxi()[2];

      for (CI p = derivnz.begin(); p != derivnz.end(); ++p) derivtxiy[p->first] -= (p->second);

      for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxiz[p->first] += (p->second);
    }

    if (ltxi < 1e-12) dserror("tangent txi length is zero!!!");

    // normalize txi directional derivative
    // (identical to normalization of normal derivative)
    typedef CORE::GEN::pairedvector<int, double>::const_iterator CI;
    CORE::GEN::pairedvector<int, double>& derivtxix = cnode.Data().GetDerivTxi()[0];
    CORE::GEN::pairedvector<int, double>& derivtxiy = cnode.Data().GetDerivTxi()[1];
    CORE::GEN::pairedvector<int, double>& derivtxiz = cnode.Data().GetDerivTxi()[2];
    CORE::GEN::pairedvector<int, double> cderivtxix = cnode.Data().GetDerivTxi()[0];
    CORE::GEN::pairedvector<int, double> cderivtxiy = cnode.Data().GetDerivTxi()[1];
    CORE::GEN::pairedvector<int, double> cderivtxiz = cnode.Data().GetDerivTxi()[2];
    const double txtx = cnode.Data().txi()[0] * cnode.Data().txi()[0];
    const double txty = cnode.Data().txi()[0] * cnode.Data().txi()[1];
    const double txtz = cnode.Data().txi()[0] * cnode.Data().txi()[2];
    const double tyty = cnode.Data().txi()[1] * cnode.Data().txi()[1];
    const double tytz = cnode.Data().txi()[1] * cnode.Data().txi()[2];
    const double tztz = cnode.Data().txi()[2] * cnode.Data().txi()[2];

    // build a vector with all keys from x,y,z maps
    // (we need this in order not to miss any entry!)
    std::vector<int> allkeyst;
    for (CI p = derivtxix.begin(); p != derivtxix.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }
    for (CI p = derivtxiy.begin(); p != derivtxiy.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }
    for (CI p = derivtxiz.begin(); p != derivtxiz.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }

    // normalize x-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxix[allkeyst[j]];
      derivtxix[allkeyst[j]] =
          (val - txtx * val - txty * cderivtxiy[allkeyst[j]] - txtz * cderivtxiz[allkeyst[j]]) /
          ltxi;
    }

    // normalize y-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxiy[allkeyst[j]];
      derivtxiy[allkeyst[j]] =
          (val - txty * cderivtxix[allkeyst[j]] - tyty * val - tytz * cderivtxiz[allkeyst[j]]) /
          ltxi;
    }

    // normalize z-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxiz[allkeyst[j]];
      derivtxiz[allkeyst[j]] =
          (val - txtz * cderivtxix[allkeyst[j]] - tytz * cderivtxiy[allkeyst[j]] - tztz * val) /
          ltxi;
    }

    // get normalized tangent derivative teta
    // use corkscrew rule from BuildAveragedNormal()
    CORE::GEN::pairedvector<int, double>& derivtetax = cnode.Data().GetDerivTeta()[0];
    CORE::GEN::pairedvector<int, double>& derivtetay = cnode.Data().GetDerivTeta()[1];
    CORE::GEN::pairedvector<int, double>& derivtetaz = cnode.Data().GetDerivTeta()[2];

    for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
    {
      derivtetay[p->first] -= cnode.Data().txi()[2] * (p->second);
      derivtetaz[p->first] += cnode.Data().txi()[1] * (p->second);
    }
    for (CI p = derivny.begin(); p != derivny.end(); ++p)
    {
      derivtetax[p->first] += cnode.Data().txi()[2] * (p->second);
      derivtetaz[p->first] -= cnode.Data().txi()[0] * (p->second);
    }
    for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
    {
      derivtetax[p->first] -= cnode.Data().txi()[1] * (p->second);
      derivtetay[p->first] += cnode.Data().txi()[0] * (p->second);
    }
    for (CI p = derivtxix.begin(); p != derivtxix.end(); ++p)
    {
      derivtetay[p->first] += cnode.MoData().n()[2] * (p->second);
      derivtetaz[p->first] -= cnode.MoData().n()[1] * (p->second);
    }
    for (CI p = derivtxiy.begin(); p != derivtxiy.end(); ++p)
    {
      derivtetax[p->first] -= cnode.MoData().n()[2] * (p->second);
      derivtetaz[p->first] += cnode.MoData().n()[0] * (p->second);
    }
    for (CI p = derivtxiz.begin(); p != derivtxiz.end(); ++p)
    {
      derivtetax[p->first] += cnode.MoData().n()[1] * (p->second);
      derivtetay[p->first] -= cnode.MoData().n()[0] * (p->second);
    }

    // OLD VERSION:

    //    if (abs(cnode.MoData().n()[0])>1.0e-6 || abs(cnode.MoData().n()[1])>1.0e-6 )
    //    {
    //      for (CI
    //      p=cnode.Data().GetDerivN()[1].begin();p!=cnode.Data().GetDerivN()[1].end();++p)
    //        (cnode.Data().GetDerivTxi()[0])[p->first] -= (p->second);
    //      for (CI
    //      p=cnode.Data().GetDerivN()[0].begin();p!=cnode.Data().GetDerivN()[0].end();++p)
    //        (cnode.Data().GetDerivTxi()[1])[p->first] += (p->second);
    //    }
    //    else
    //    {
    //      for (CI
    //      p=cnode.Data().GetDerivN()[2].begin();p!=cnode.Data().GetDerivN()[2].end();++p)
    //        (cnode.Data().GetDerivTxi()[1])[p->first] -= (p->second);
    //      for (CI
    //      p=cnode.Data().GetDerivN()[1].begin();p!=cnode.Data().GetDerivN()[1].end();++p)
    //        (cnode.Data().GetDerivTxi()[2])[p->first] += (p->second);
    //    }
    //
    //    double ltxi =
    //    sqrt(cnode.Data().txi()[0]*cnode.Data().txi()[0]+cnode.Data().txi()[1]*cnode.Data().txi()[1]+cnode.Data().txi()[2]*cnode.Data().txi()[2]);
    //    for (int j=0;j<3;++j) cnode.Data().txi()[j]/=ltxi;

    //     //teta follows from corkscrew rule (teta = n x txi)
    //    cnode.Data().teta()[0] =
    //    cnode.MoData().n()[1]*cnode.Data().txi()[2]-cnode.MoData().n()[2]*cnode.Data().txi()[1];
    //    cnode.Data().teta()[1] =
    //    cnode.MoData().n()[2]*cnode.Data().txi()[0]-cnode.MoData().n()[0]*cnode.Data().txi()[2];
    //    cnode.Data().teta()[2] =
    //    cnode.MoData().n()[0]*cnode.Data().txi()[1]-cnode.MoData().n()[1]*cnode.Data().txi()[0];
    //
    //    for (CI
    //    p=cnode.Data().GetDerivN()[1].begin();p!=cnode.Data().GetDerivN()[1].end();++p)
    //      (cnode.Data().GetDerivTeta()[0])[p->first] += (p->second) * cnode.Data().txi()[2];
    //    for (CI
    //    p=cnode.Data().GetDerivTxi()[2].begin();p!=cnode.Data().GetDerivTxi()[2].end();++p)
    //      (cnode.Data().GetDerivTeta()[0])[p->first] += (p->second) * cnode.MoData().n()[1];
    //    for (CI
    //    p=cnode.Data().GetDerivN()[2].begin();p!=cnode.Data().GetDerivN()[2].end();++p)
    //      (cnode.Data().GetDerivTeta()[0])[p->first] -= (p->second) * cnode.Data().txi()[1];
    //    for (CI
    //    p=cnode.Data().GetDerivTxi()[1].begin();p!=cnode.Data().GetDerivTxi()[1].end();++p)
    //      (cnode.Data().GetDerivTeta()[0])[p->first] -= (p->second) * cnode.MoData().n()[2];
    //
    //    for (CI
    //    p=cnode.Data().GetDerivN()[2].begin();p!=cnode.Data().GetDerivN()[2].end();++p)
    //      (cnode.Data().GetDerivTeta()[1])[p->first] += (p->second) * cnode.Data().txi()[0];
    //    for (CI
    //    p=cnode.Data().GetDerivTxi()[0].begin();p!=cnode.Data().GetDerivTxi()[0].end();++p)
    //      (cnode.Data().GetDerivTeta()[1])[p->first] += (p->second) * cnode.MoData().n()[2];
    //    for (CI
    //    p=cnode.Data().GetDerivN()[0].begin();p!=cnode.Data().GetDerivN()[0].end();++p)
    //      (cnode.Data().GetDerivTeta()[1])[p->first] -= (p->second) * cnode.Data().txi()[2];
    //    for (CI
    //    p=cnode.Data().GetDerivTxi()[2].begin();p!=cnode.Data().GetDerivTxi()[2].end();++p)
    //      (cnode.Data().GetDerivTeta()[1])[p->first] -= (p->second) * cnode.MoData().n()[0];
    //
    //    for (CI
    //    p=cnode.Data().GetDerivN()[0].begin();p!=cnode.Data().GetDerivN()[0].end();++p)
    //      (cnode.Data().GetDerivTeta()[2])[p->first] += (p->second) * cnode.Data().txi()[1];
    //    for (CI
    //    p=cnode.Data().GetDerivTxi()[1].begin();p!=cnode.Data().GetDerivTxi()[1].end();++p)
    //      (cnode.Data().GetDerivTeta()[2])[p->first] += (p->second) * cnode.MoData().n()[0];
    //    for (CI
    //    p=cnode.Data().GetDerivN()[1].begin();p!=cnode.Data().GetDerivN()[1].end();++p)
    //      (cnode.Data().GetDerivTeta()[2])[p->first] -= (p->second) * cnode.Data().txi()[0];
    //    for (CI
    //    p=cnode.Data().GetDerivTxi()[0].begin();p!=cnode.Data().GetDerivTxi()[0].end();++p)
    //      (cnode.Data().GetDerivTeta()[2])[p->first] -= (p->second) * cnode.MoData().n()[1];
  }

  return;
}


/*----------------------------------------------------------------------*
 |  export nodal normals (public)                             popp 11/10|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::ExportNodalNormals() const
{
  // create empty data objects
  std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> triad;

  std::map<int, std::vector<int>> n_x_key;
  std::map<int, std::vector<int>> n_y_key;
  std::map<int, std::vector<int>> n_z_key;
  std::map<int, std::vector<int>> txi_x_key;
  std::map<int, std::vector<int>> txi_y_key;
  std::map<int, std::vector<int>> txi_z_key;
  std::map<int, std::vector<int>> teta_x_key;
  std::map<int, std::vector<int>> teta_y_key;
  std::map<int, std::vector<int>> teta_z_key;

  std::map<int, std::vector<double>> n_x_val;
  std::map<int, std::vector<double>> n_y_val;
  std::map<int, std::vector<double>> n_z_val;
  std::map<int, std::vector<double>> txi_x_val;
  std::map<int, std::vector<double>> txi_y_val;
  std::map<int, std::vector<double>> txi_z_val;
  std::map<int, std::vector<double>> teta_x_val;
  std::map<int, std::vector<double>> teta_y_val;
  std::map<int, std::vector<double>> teta_z_val;

  std::map<int, double>::iterator iter;
  CORE::GEN::pairedvector<int, double>::iterator _iter;

  // build info on row map
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // fill nodal matrix
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> loc =
        Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(3, 3));
    (*loc)(0, 0) = cnode->MoData().n()[0];
    (*loc)(1, 0) = cnode->MoData().n()[1];
    (*loc)(2, 0) = cnode->MoData().n()[2];
    (*loc)(0, 1) = cnode->Data().txi()[0];
    (*loc)(1, 1) = cnode->Data().txi()[1];
    (*loc)(2, 1) = cnode->Data().txi()[2];
    (*loc)(0, 2) = cnode->Data().teta()[0];
    (*loc)(1, 2) = cnode->Data().teta()[1];
    (*loc)(2, 2) = cnode->Data().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->Data().GetDerivN();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivtxi = cnode->Data().GetDerivTxi();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivteta = cnode->Data().GetDerivTeta();

    for (_iter = derivn[0].begin(); _iter != derivn[0].end(); ++_iter)
    {
      n_x_key[gid].push_back(_iter->first);
      n_x_val[gid].push_back(_iter->second);
    }
    for (_iter = derivn[1].begin(); _iter != derivn[1].end(); ++_iter)
    {
      n_y_key[gid].push_back(_iter->first);
      n_y_val[gid].push_back(_iter->second);
    }
    for (_iter = derivn[2].begin(); _iter != derivn[2].end(); ++_iter)
    {
      n_z_key[gid].push_back(_iter->first);
      n_z_val[gid].push_back(_iter->second);
    }

    for (_iter = derivtxi[0].begin(); _iter != derivtxi[0].end(); ++_iter)
    {
      txi_x_key[gid].push_back(_iter->first);
      txi_x_val[gid].push_back(_iter->second);
    }
    for (_iter = derivtxi[1].begin(); _iter != derivtxi[1].end(); ++_iter)
    {
      txi_y_key[gid].push_back(_iter->first);
      txi_y_val[gid].push_back(_iter->second);
    }
    for (_iter = derivtxi[2].begin(); _iter != derivtxi[2].end(); ++_iter)
    {
      txi_z_key[gid].push_back(_iter->first);
      txi_z_val[gid].push_back(_iter->second);
    }

    for (_iter = derivteta[0].begin(); _iter != derivteta[0].end(); ++_iter)
    {
      teta_x_key[gid].push_back(_iter->first);
      teta_x_val[gid].push_back(_iter->second);
    }
    for (_iter = derivteta[1].begin(); _iter != derivteta[1].end(); ++_iter)
    {
      teta_y_key[gid].push_back(_iter->first);
      teta_y_val[gid].push_back(_iter->second);
    }
    for (_iter = derivteta[2].begin(); _iter != derivteta[2].end(); ++_iter)
    {
      teta_z_key[gid].push_back(_iter->first);
      teta_z_val[gid].push_back(_iter->second);
    }
  }


  // communicate from slave node row to column map
  CORE::COMM::Exporter& ex = interfaceData_->Exporter();

  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for (int i = 0; i < snodecolmapbound_->NumMyElements(); ++i)
  {
    // only do something for ghosted nodes
    int gid = snodecolmapbound_->GID(i);
    if (snoderowmapbound_->MyGID(gid)) continue;

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    int linsize = cnode->GetLinsize() + (int)(n_x_key[gid].size());

    // extract info
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> loc = triad[gid];
    cnode->MoData().n()[0] = (*loc)(0, 0);
    cnode->MoData().n()[1] = (*loc)(1, 0);
    cnode->MoData().n()[2] = (*loc)(2, 0);
    cnode->Data().txi()[0] = (*loc)(0, 1);
    cnode->Data().txi()[1] = (*loc)(1, 1);
    cnode->Data().txi()[2] = (*loc)(2, 1);
    cnode->Data().teta()[0] = (*loc)(0, 2);
    cnode->Data().teta()[1] = (*loc)(1, 2);
    cnode->Data().teta()[2] = (*loc)(2, 2);

    // extract derivative info
    std::vector<CORE::GEN::pairedvector<int, double>>& derivn = cnode->Data().GetDerivN();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivtxi = cnode->Data().GetDerivTxi();
    std::vector<CORE::GEN::pairedvector<int, double>>& derivteta = cnode->Data().GetDerivTeta();

    for (int k = 0; k < (int)(derivn.size()); ++k) derivn[k].clear();
    derivn.resize(3, linsize);
    for (int k = 0; k < (int)(derivtxi.size()); ++k) derivtxi[k].clear();
    derivtxi.resize(3, linsize);
    for (int k = 0; k < (int)(derivteta.size()); ++k) derivteta[k].clear();
    derivteta.resize(3, linsize);

    for (int k = 0; k < (int)(n_x_key[gid].size()); ++k)
      (cnode->Data().GetDerivN()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k = 0; k < (int)(n_y_key[gid].size()); ++k)
      (cnode->Data().GetDerivN()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k = 0; k < (int)(n_z_key[gid].size()); ++k)
      (cnode->Data().GetDerivN()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k = 0; k < (int)(txi_x_key[gid].size()); ++k)
      (cnode->Data().GetDerivTxi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k = 0; k < (int)(txi_y_key[gid].size()); ++k)
      (cnode->Data().GetDerivTxi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k = 0; k < (int)(txi_z_key[gid].size()); ++k)
      (cnode->Data().GetDerivTxi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k = 0; k < (int)(teta_x_key[gid].size()); ++k)
      (cnode->Data().GetDerivTeta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k = 0; k < (int)(teta_y_key[gid].size()); ++k)
      (cnode->Data().GetDerivTeta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k = 0; k < (int)(teta_z_key[gid].size()); ++k)
      (cnode->Data().GetDerivTeta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();

  /*std::cout << "---------- nodal normals
  ----------------------------------------------------------" << std::endl;

  // print nodal normals
  for (int p=0;p<Comm().NumProc();++p)
  {
    // one proc after the other
    if (p==Comm().MyPID())
    {
      std::cout << "\n*****\nPROC " << p << "\n*****" << std::endl;
      for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
      {
        int gid = snodecolmapbound_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("Cannot find node with gid %",gid);
        Node* cnode = dynamic_cast<Node*>(node);

        // print normal and tangents at each slave node
        std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " Normal: " << cnode->MoData().n()[0]
             << " " << cnode->MoData().n()[1] << " " << cnode->MoData().n()[2] << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " TXi: " << cnode->Data().txi()[0]
             << " " << cnode->Data().txi()[1] << " " << cnode->Data().txi()[2] << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " TEta: " << cnode->Data().teta()[0]
             << " " << cnode->Data().teta()[1] << " " << cnode->Data().teta()[2] << std::endl;

        // print linearizations at each slave node
        std::cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinN:
  "; for
  (_iter=cnode->Data().GetDerivN()[0].begin();_iter!=cnode->Data().GetDerivN()[0].end();++_iter)
          std::cout << "\n" << _iter->first << "\t" << _iter->second;
        std::cout << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinTxi:
  "; for
  (_iter=cnode->Data().GetDerivTxi()[0].begin();_iter!=cnode->Data().GetDerivTxi()[0].end();++_iter)
          std::cout << "\n" << _iter->first << "\t" << _iter->second;
        std::cout << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << "
  LinTeta: "; for
  (_iter=cnode->Data().GetDerivTeta()[0].begin();_iter!=cnode->Data().GetDerivTeta()[0].end();++_iter)
          std::cout << "\n" << _iter->first << "\t" << _iter->second;
        std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;
    }

    // barrier
    Comm().Barrier();
  }*/

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 10/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::EvaluateSearchBinarytree()
{
  // *********************************************************************
  // self contact:
  // *********************************************************************
  // We call EvaluateSearch(), which does both, the bottom-up update (on the whole interface) and
  // the search. Then the dynamic master/slave assignment routine UpdateMasterSlaveSets() is called
  // and the new slave nodes' data containers are initialized.
  // *********************************************************************
  if (SelfContact())
  {
    // evaluate search itself
    binarytreeself_->EvaluateSearch();

    // update master/slave sets of interface
    UpdateMasterSlaveSets();

    // initialize node data container
    // (include slave side boundary nodes / crosspoints)
    for (int i = 0; i < SlaveColNodesBound()->NumMyElements(); ++i)
    {
      int gid = SlaveColNodesBound()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("Cannot find node with gid %i", gid);
      MORTAR::Node* mnode = dynamic_cast<MORTAR::Node*>(node);

      // initialize container if not yet initialized before
      mnode->InitializeDataContainer();
    }

    // no initialization of element data container as this would
    // possibly destroy the information on search elements again
    // (this was already done in SetElementAreas())
  }

  else
  {
    // call mortar routine
    MORTAR::Interface::EvaluateSearchBinarytree();
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-line coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateSTL()
{
  // check
  if (Dim() == 2) dserror("LTS algorithm only for 3D simulations!");

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // guarantee uniquness
    std::set<std::pair<int, int>> donebefore;

    // loop over found meles
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      if (melement->Shape() == CORE::FE::CellType::quad4)
      {
        for (int j = 0; j < 4; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = melement->NodeIds()[3];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

          // if not then create ele
          if (iter == donebefore.end() and itertw == donebefore.end())
          {
            // add to set of processed nodes
            donebefore.insert(actIDs);
            donebefore.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
                j, melement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<DRT::Node*, 2> nodes = {
                melement->Nodes()[nodeLIds[0]], melement->Nodes()[nodeLIds[1]]};
            lineEle->BuildNodalPointers(nodes.data());

            // init data container for dual shapes
            lineEle->InitializeDataContainer();

            std::vector<Element*> seleElements;
            seleElements.push_back(selement);

            // create coupling object
            LineToSurfaceCoupling3d coup(*idiscret_, 3, InterfaceParams(), *melement, lineEle,
                seleElements, LineToSurfaceCoupling3d::stl);

            // perform evaluate!
            coup.EvaluateCoupling();
          }
        }  // end edge loop
      }
      else
        dserror("LTS only for quad4!");

    }  // end found mele loop
  }    // end slave ele loop

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type ndoe-to-segment coupl             farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateNTSMaster()
{
  // create one interpolator instance which is valid for all nodes!
  Teuchos::RCP<NTS::Interpolator> interpolator =
      Teuchos::rcp(new NTS::Interpolator(InterfaceParams(), Dim()));

  // guarantee uniquness
  std::set<int> donebefore;

  // loop over all slave col elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      for (int n = 0; n < (int)melement->NumNode(); ++n)
      {
        MORTAR::Node* cnode = dynamic_cast<MORTAR::Node*>(melement->Nodes()[n]);

        std::set<int>::iterator iter = donebefore.find(cnode->Id());

        // if not evaluated before
        if (iter == donebefore.end())
        {
          std::vector<MORTAR::Element*> dummy;
          MORTAR::Element* sele = dynamic_cast<MORTAR::Element*>(selement);
          dummy.push_back(sele);

          // call interpolation functions
          bool success = interpolator->Interpolate(*cnode, dummy);
          if (success) donebefore.insert(cnode->Id());
        }

      }  // node loop
    }    // found contact eles
  }      // sele loop

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateLTSMaster()
{
  // check
  if (Dim() == 2) dserror("LTS algorithm only for 3D simulations!");

  // guarantee uniquness
  std::set<std::pair<int, int>> donebefore;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // ele check
    if (selement->Shape() != CORE::FE::CellType::quad4 and
        selement->Shape() != CORE::FE::CellType::tri3)
      dserror("LTS algorithm only for tri3/quad4!");

    // empty vector of master element pointers
    std::vector<Teuchos::RCP<MORTAR::Element>> lineElements;
    std::vector<Element*> meleElements;

    // compute slave normal
    double slaveN[3] = {0.0, 0.0, 0.0};
    double loccenter[2] = {0.0, 0.0};

    CORE::FE::CellType dt = selement->Shape();
    if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
    {
      loccenter[0] = 1.0 / 3.0;
      loccenter[1] = 1.0 / 3.0;
    }
    else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
             dt == CORE::FE::CellType::quad9)
    {
      loccenter[0] = 0.0;
      loccenter[1] = 0.0;
    }
    else
      dserror("AuxiliaryPlane called for unknown element type");

    // we then compute the unit normal vector at the element center
    selement->ComputeUnitNormalAtXi(loccenter, slaveN);

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      // check orientation
      // compute slave normal
      double masterN[3] = {0.0, 0.0, 0.0};
      double loccenterM[2] = {0.0, 0.0};

      CORE::FE::CellType dt = melement->Shape();
      if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
      {
        loccenterM[0] = 1.0 / 3.0;
        loccenterM[1] = 1.0 / 3.0;
      }
      else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
               dt == CORE::FE::CellType::quad9)
      {
        loccenterM[0] = 0.0;
        loccenterM[1] = 0.0;
      }
      else
        dserror("AuxiliaryPlane called for unknown element type");

      // we then compute the unit normal vector at the element center
      melement->ComputeUnitNormalAtXi(loccenterM, masterN);

      double scaprod = slaveN[0] * masterN[0] + slaveN[1] * masterN[1] + slaveN[2] * masterN[2];

      // tolerance for line clipping
      const double sminedge = selement->MinEdgeSize();
      const double mminedge = melement->MinEdgeSize();
      const double tol = 0.001 * std::min(sminedge, mminedge);
      if (abs(scaprod) < tol) continue;

      // if orientation is okay
      meleElements.push_back(melement);
    }

    // no valid maste elements?
    if (meleElements.size() < 1) continue;

    for (int m = 0; m < (int)meleElements.size(); ++m)
    {
      // loop over slave edges -> match node number for tri3/quad4
      for (int j = 0; j < meleElements[m]->NumNode(); ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (meleElements[m]->Shape() == CORE::FE::CellType::quad4)
        {
          if (j == 0)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[0];
            nodeIds[1] = meleElements[m]->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[1];
            nodeIds[1] = meleElements[m]->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[2];
            nodeIds[1] = meleElements[m]->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[3];
            nodeIds[1] = meleElements[m]->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            dserror("loop counter and edge number do not match!");
        }
        else if (meleElements[m]->Shape() == CORE::FE::CellType::tri3)
        {
          if (j == 0)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[0];
            nodeIds[1] = meleElements[m]->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[1];
            nodeIds[1] = meleElements[m]->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[2];
            nodeIds[1] = meleElements[m]->NodeIds()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }
          else
            dserror("loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge =
            dynamic_cast<MORTAR::Node*>(meleElements[m]->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge =
            dynamic_cast<MORTAR::Node*>(meleElements[m]->Nodes()[nodeLIds[1]])->IsOnEdge();

        if (nonSmoothContact_ and (!node0Edge or !node1Edge)) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

        // if not then create ele
        if (iter == donebefore.end() and itertw == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(actIDs);
          donebefore.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
              j, meleElements[m]->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<DRT::Node*, 2> nodes = {
              meleElements[m]->Nodes()[nodeLIds[0]], meleElements[m]->Nodes()[nodeLIds[1]]};
          lineEle->BuildNodalPointers(nodes.data());

          // init data container for dual shapes
          lineEle->InitializeDataContainer();

          // push back into vector
          lineElements.push_back(lineEle);

          std::vector<Element*> dummy;
          dummy.push_back(selement);

          // create coupling object
          LineToSurfaceCoupling3d coup(*idiscret_, 3, InterfaceParams(), *meleElements[m], lineEle,
              dummy, LineToSurfaceCoupling3d::lts);


          // perform evaluate!
          coup.EvaluateCoupling();
        }
      }  // end edge loop
    }
  }  // slave ele loop
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateLTS()
{
  // check
  if (Dim() == 2) dserror("LTS algorithm only for 3D simulations!");

  // guarantee uniquness
  std::set<std::pair<int, int>> donebefore;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // ele check
    if (selement->Shape() != CORE::FE::CellType::quad4 and
        selement->Shape() != CORE::FE::CellType::tri3)
      dserror("LTS algorithm only for tri3/quad4!");

    // empty vector of master element pointers
    std::vector<Teuchos::RCP<MORTAR::Element>> lineElements;
    std::vector<Element*> meleElements;

    // compute slave normal
    double slaveN[3] = {0.0, 0.0, 0.0};
    double loccenter[2] = {0.0, 0.0};

    CORE::FE::CellType dt = selement->Shape();
    if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
    {
      loccenter[0] = 1.0 / 3.0;
      loccenter[1] = 1.0 / 3.0;
    }
    else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
             dt == CORE::FE::CellType::quad9)
    {
      loccenter[0] = 0.0;
      loccenter[1] = 0.0;
    }
    else
      dserror("AuxiliaryPlane called for unknown element type");

    // we then compute the unit normal vector at the element center
    selement->ComputeUnitNormalAtXi(loccenter, slaveN);

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      // check orientation
      // compute slave normal
      double masterN[3] = {0.0, 0.0, 0.0};
      double loccenterM[2] = {0.0, 0.0};

      CORE::FE::CellType dt = melement->Shape();
      if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
      {
        loccenterM[0] = 1.0 / 3.0;
        loccenterM[1] = 1.0 / 3.0;
      }
      else if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
               dt == CORE::FE::CellType::quad9)
      {
        loccenterM[0] = 0.0;
        loccenterM[1] = 0.0;
      }
      else
        dserror("AuxiliaryPlane called for unknown element type");

      // we then compute the unit normal vector at the element center
      melement->ComputeUnitNormalAtXi(loccenterM, masterN);

      double scaprod = slaveN[0] * masterN[0] + slaveN[1] * masterN[1] + slaveN[2] * masterN[2];

      // tolerance for line clipping
      const double sminedge = selement->MinEdgeSize();
      const double mminedge = melement->MinEdgeSize();
      const double tol = 0.001 * std::min(sminedge, mminedge);
      if (abs(scaprod) < tol) continue;

      // if orientation is okay
      meleElements.push_back(melement);
    }

    // no valid maste elements?
    if (meleElements.size() < 1) continue;

    // loop over slave edges -> match node number for tri3/quad4
    for (int j = 0; j < selement->NumNode(); ++j)
    {
      int nodeIds[2] = {0, 0};
      int nodeLIds[2] = {0, 0};

      if (selement->Shape() == CORE::FE::CellType::quad4)
      {
        if (j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = selement->NodeIds()[3];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }
        else
          dserror("loop counter and edge number do not match!");
      }
      else if (selement->Shape() == CORE::FE::CellType::tri3)
      {
        if (j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 2;
          nodeLIds[1] = 0;
        }
        else
          dserror("loop counter and edge number do not match!");
      }

      // check if both nodes on edge geometry
      bool node0Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
      bool node1Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

      if (nonSmoothContact_ and (!node0Edge or !node1Edge)) continue;

      // create pair
      std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
      std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

      // check if processed before
      std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
      std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

      // if not then create ele
      if (iter == donebefore.end() and itertw == donebefore.end())
      {
        // add to set of processed nodes
        donebefore.insert(actIDs);
        donebefore.insert(actIDstw);

        // create line ele:
        Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
            j, selement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

        // get nodes
        std::array<DRT::Node*, 2> nodes = {
            selement->Nodes()[nodeLIds[0]], selement->Nodes()[nodeLIds[1]]};
        lineEle->BuildNodalPointers(nodes.data());

        // init data container for dual shapes
        lineEle->InitializeDataContainer();

        // push back into vector
        lineElements.push_back(lineEle);
      }
    }  // end edge loop

    // loop over all created line elements
    for (int l = 0; l < (int)lineElements.size(); ++l)
    {
      // create coupling object
      LineToSurfaceCoupling3d coup(*idiscret_, 3, InterfaceParams(), *selement, lineElements[l],
          meleElements, LineToSurfaceCoupling3d::lts);

      // perform evaluate!
      coup.EvaluateCoupling();
    }
  }  // slave ele loop

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-line coupl                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateLTL()
{
  // check
  if (Dim() == 2) dserror("LTL algorithm only for 3D simulations!");

  // guarantee uniquness of slave edges
  std::set<std::pair<int, int>> donebeforeS;

  // loop over slave elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    Element* selement = dynamic_cast<Element*>(ele1);

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::Element>> lineElementsS;

    if (selement->Shape() == CORE::FE::CellType::quad4)
    {
      for (int j = 0; j < 4; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = selement->NodeIds()[3];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebeforeS.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebeforeS.find(actIDstw);

        // if not then create ele
        if (iter == donebeforeS.end() and itertw == donebeforeS.end())
        {
          // add to set of processed nodes
          donebeforeS.insert(actIDs);
          donebeforeS.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
              j, selement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<DRT::Node*, 2> nodes = {
              selement->Nodes()[nodeLIds[0]], selement->Nodes()[nodeLIds[1]]};
          lineEle->BuildNodalPointers(nodes.data());

          // init data container for dual shapes
          lineEle->InitializeDataContainer();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else if (selement->Shape() == CORE::FE::CellType::tri3)
    {
      for (int j = 0; j < 3; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 2;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::Node*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebeforeS.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebeforeS.find(actIDstw);

        // if not then create ele
        if (iter == donebeforeS.end() and itertw == donebeforeS.end())
        {
          // add to set of processed nodes
          donebeforeS.insert(actIDs);
          donebeforeS.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
              j, selement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<DRT::Node*, 2> nodes = {
              selement->Nodes()[nodeLIds[0]], selement->Nodes()[nodeLIds[1]]};
          lineEle->BuildNodalPointers(nodes.data());

          // init data container for dual shapes
          lineEle->InitializeDataContainer();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else
      dserror("LTL only for quad4 and tri3!");

    // guarantee uniquness of master edges
    std::set<std::pair<int, int>> donebeforeM;

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::Element>> lineElementsM;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < selement->MoData().NumSearchElements(); ++k)
    {
      int gid2 = selement->MoData().SearchElements()[k];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      Element* melement = dynamic_cast<Element*>(ele2);

      if (melement->Shape() == CORE::FE::CellType::quad4)
      {
        for (int j = 0; j < 4; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if (j == 3)
          {
            nodeIds[0] = melement->NodeIds()[3];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge = dynamic_cast<MORTAR::Node*>(melement->Nodes()[nodeLIds[0]])->IsOnEdge();
          bool node1Edge = dynamic_cast<MORTAR::Node*>(melement->Nodes()[nodeLIds[1]])->IsOnEdge();

          if (!node0Edge or !node1Edge) continue;

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebeforeM.find(actIDstw);

          // if not then create ele
          if (iter == donebeforeM.end() and itertw == donebeforeM.end())
          {
            // add to set of processed nodes
            donebeforeM.insert(actIDs);
            donebeforeM.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
                j, melement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<DRT::Node*, 2> nodes = {
                melement->Nodes()[nodeLIds[0]], melement->Nodes()[nodeLIds[1]]};
            lineEle->BuildNodalPointers(nodes.data());

            // init data container for dual shapes
            lineEle->InitializeDataContainer();

            // push back into vector
            lineElementsM.push_back(lineEle);
          }
        }  // end edge loop
      }
      else if (melement->Shape() == CORE::FE::CellType::tri3)
      {
        for (int j = 0; j < 3; ++j)
        {
          int nodeIds[2] = {0, 0};
          int nodeLIds[2] = {0, 0};

          if (j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if (j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if (j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge = dynamic_cast<MORTAR::Node*>(melement->Nodes()[nodeLIds[0]])->IsOnEdge();
          bool node1Edge = dynamic_cast<MORTAR::Node*>(melement->Nodes()[nodeLIds[1]])->IsOnEdge();

          if (!node0Edge or !node1Edge) continue;

          // create pair
          std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
          std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

          // check if processed before
          std::set<std::pair<int, int>>::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int, int>>::iterator itertw = donebeforeM.find(actIDstw);

          // if not then create ele
          if (iter == donebeforeM.end() and itertw == donebeforeM.end())
          {
            // add to set of processed nodes
            donebeforeM.insert(actIDs);
            donebeforeM.insert(actIDstw);

            // create line ele:
            Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(new MORTAR::Element(
                j, melement->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

            // get nodes
            std::array<DRT::Node*, 2> nodes = {
                melement->Nodes()[nodeLIds[0]], melement->Nodes()[nodeLIds[1]]};
            lineEle->BuildNodalPointers(nodes.data());

            // init data container for dual shapes
            lineEle->InitializeDataContainer();

            // push back into vector
            lineElementsM.push_back(lineEle);
          }
        }  // end edge loop
      }
      else
        dserror("LTL only for quad4 and tri3!");
    }  // end found mele loop

    // loop over slave edges
    for (int s = 0; s < (int)lineElementsS.size(); ++s)
    {
      // loop over master edges
      for (int m = 0; m < (int)lineElementsM.size(); ++m)
      {
        // create coupling object
        LineToLineCouplingPoint3d coup(
            *idiscret_, 3, InterfaceParams(), lineElementsS[s], lineElementsM[m]);

        // perform evaluate!
        coup.EvaluateCoupling();
      }
    }

  }  // end slave loop

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type node-to-segment coupl             farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateNTS()
{
  // create one interpolator instance which is valid for all nodes!
  Teuchos::RCP<NTS::Interpolator> interpolator =
      Teuchos::rcp(new NTS::Interpolator(InterfaceParams(), Dim()));

  // loop over slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    MORTAR::Node* mrtrnode = dynamic_cast<MORTAR::Node*>(node);

    if (!mrtrnode->IsOnCorner() and nonSmoothContact_) continue;

    if (mrtrnode->Owner() != Comm().MyPID()) dserror("Node ownership inconsistency!");

    // vector with possible contacting master eles
    std::vector<MORTAR::Element*> meles;

    // fill vector with possibly contacting meles
    FindMEles(*mrtrnode, meles);

    // skip calculation if no meles vector is empty
    if (meles.size() < 1) continue;

    // call interpolation functions
    interpolator->Interpolate(*mrtrnode, meles);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlaps     popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::MortarCoupling(MORTAR::Element* sele, std::vector<MORTAR::Element*> mele,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // do stuff before the actual coupling is going to be evaluated
  PreMortarCoupling(sele, mele, mparams_ptr);

  // increase counter of slave/master pairs
  smpairs_ += (int)mele.size();

  // check if quadratic interpolation is involved
  bool quadratic = false;
  if (sele->IsQuad()) quadratic = true;
  for (int m = 0; m < (int)mele.size(); ++m)
    if (mele[m]->IsQuad()) quadratic = true;

  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (Dim() == 2)
  {
    // *************************************************** linear 2D ***
    // ************************************************ quadratic 2D ***
    // neither quadratic interpolation nor mixed linear and quadratic
    // interpolation need any special treatment in the 2d case

    // create Coupling2dManager
    CONTACT::Coupling2dManager coup(Discret(), Dim(), quadratic, InterfaceParams(), sele, mele);
    // evaluate
    coup.EvaluateCoupling(mparams_ptr);

    // increase counter of slave/master integration pairs and intcells
    smintpairs_ += (int)mele.size();
    intcells_ += (int)mele.size();
  }
  // ************************************************************** 3D ***
  else if (Dim() == 3)
  {
    // *************************************************** linear 3D ***
    if (!quadratic)
    {
      // create Coupling3dManager
      CONTACT::Coupling3dManager coup(Discret(), Dim(), quadratic, InterfaceParams(), sele, mele);
      // evaluate
      coup.EvaluateCoupling(mparams_ptr);

      // increase counter of slave/master integration pairs and intcells
      smintpairs_ += (int)mele.size();
      intcells_ += coup.IntegrationCells();
    }

    // ************************************************** quadratic 3D ***
    else
    {
      // create Coupling3dQuadManager
      CONTACT::Coupling3dQuadManager coup(
          Discret(), Dim(), quadratic, InterfaceParams(), sele, mele);
      // evaluate
      coup.EvaluateCoupling(mparams_ptr);
    }  // quadratic
  }    // 3D
  else
    dserror("Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  // do stuff after the coupling evaluation
  PostMortarCoupling(sele, mele, mparams_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate penalty scaling factor kapp (public)            popp 11/09|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::IntegrateKappaPenalty(CONTACT::Element& sele)
{
  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape() == CORE::FE::CellType::tri3 || sele.Shape() == CORE::FE::CellType::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0;
    sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0;
    sxia[1] = -1.0;
    sxib[0] = 1.0;
    sxib[1] = 1.0;
  }

  // ************************************************** quadratic 3D ***
  if (Dim() == 3 && sele.IsQuad())
  {
    // get LM interpolation and testing type
    INPAR::MORTAR::LagMultQuad lmtype =
        CORE::UTILS::IntegralValue<INPAR::MORTAR::LagMultQuad>(InterfaceParams(), "LM_QUAD");

    // build linear integration elements from quadratic CElements
    std::vector<Teuchos::RCP<MORTAR::IntElement>> sauxelements(0);
    SplitIntElements(sele, sauxelements);

    // different options for mortar integration
    if (lmtype == INPAR::MORTAR::lagmult_quad || lmtype == INPAR::MORTAR::lagmult_lin)
    {
      // do the element integration of kappa and store into gap
      int nrow = sele.NumNode();
      Teuchos::RCP<CORE::LINALG::SerialDenseVector> gseg =
          Teuchos::rcp(new CORE::LINALG::SerialDenseVector(nrow));

      // create a CONTACT integrator instance with correct NumGP and Dim
      CONTACT::Integrator integrator(imortar_, sele.Shape(), Comm());
      integrator.IntegrateKappaPenalty(sele, sxia, sxib, gseg);

      // do the assembly into the slave nodes
      integrator.AssembleG(Comm(), sele, *gseg);
    }

    else if (lmtype == INPAR::MORTAR::lagmult_pwlin)
    {
      // integrate each int element seperately
      for (int i = 0; i < (int)sauxelements.size(); ++i)
      {
        // do the int element integration of kappa and store into gap
        int nrow = sauxelements[i]->NumNode();
        Teuchos::RCP<CORE::LINALG::SerialDenseVector> gseg =
            Teuchos::rcp(new CORE::LINALG::SerialDenseVector(nrow));

        // create a CONTACT integrator instance with correct NumGP and Dim
        CONTACT::Integrator integrator(imortar_, sauxelements[i]->Shape(), Comm());
        integrator.IntegrateKappaPenalty(sele, *(sauxelements[i]), sxia, sxib, gseg);

        // do the assembly into the slave nodes
        integrator.AssembleG(Comm(), *(sauxelements[i]), *gseg);
      }
    }

    else
    {
      dserror("IntegrateKappaPenalty: Invalid case for 3D mortar contact LM interpolation");
    }
  }

  // *************************************************** other cases ***
  else
  {
    // do the element integration of kappa and store into gap
    int nrow = sele.NumNode();
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> gseg =
        Teuchos::rcp(new CORE::LINALG::SerialDenseVector(nrow));

    // create a CONTACT integrator instance with correct NumGP and Dim
    CONTACT::Integrator integrator(imortar_, sele.Shape(), Comm());
    integrator.IntegrateKappaPenalty(sele, sxia, sxib, gseg);

    // do the assembly into the slave nodes
    integrator.AssembleG(Comm(), sele, *gseg);
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate relative movement (jump) of a slave node     gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateRelMov(const Teuchos::RCP<Epetra_Vector> xsmod,
    const Teuchos::RCP<CORE::LINALG::SparseMatrix> dmatrixmod,
    const Teuchos::RCP<CORE::LINALG::SparseMatrix> doldmod)
{
  if (friction_ == false)
    dserror("Error in Interface::EvaluateRelMov(): Only evaluated for frictional contact");

  // parameters
  double pp = InterfaceParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);
    double cn = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];

    // get some informatiom form the node
    double gap = cnode->Data().Getg();

    const int dim = cnode->NumDof();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k = 0; k < 3; ++k) nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

    std::vector<double> jump(dim);
    for (int dim = 0; dim < Dim(); dim++) jump[dim] = 0;

    double lmuzawan = 0.0;
    for (int k = 0; k < dim; ++k) lmuzawan += cnode->MoData().lmuzawa()[k] * cnode->MoData().n()[k];

    double kappa = cnode->Data().Kappa();

    // evaluate jump (relative displacement) of this node
    // only when the node is going to be active, otherwise,
    // this value isn't needed.
    bool activeinfuture = false;

    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
            InterfaceParams(), "STRATEGY") == INPAR::CONTACT::solution_penalty ||
        CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
            InterfaceParams(), "STRATEGY") == INPAR::CONTACT::solution_multiscale)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
                 InterfaceParams(), "STRATEGY") == INPAR::CONTACT::solution_lagmult and
             CORE::UTILS::IntegralValue<int>(InterfaceParams(), "SEMI_SMOOTH_NEWTON") != 1)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
                 InterfaceParams(), "STRATEGY") == INPAR::CONTACT::solution_lagmult and
             CORE::UTILS::IntegralValue<int>(InterfaceParams(), "SEMI_SMOOTH_NEWTON") == 1)
    {
      if ((nz - cn * gap > 0) or cnode->Active()) activeinfuture = true;
    }
    else if (CORE::UTILS::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
                 InterfaceParams(), "STRATEGY") == INPAR::CONTACT::solution_uzawa)
    {
      if (lmuzawan - kappa * pp * gap >= 0) activeinfuture = true;
    }
    else
      dserror("Error in Interface::EvaluateRelMov(): Solution strategy not known!");

    if (activeinfuture == true)
    {
      CORE::GEN::pairedvector<int, double>& dmap = cnode->MoData().GetD();
      CORE::GEN::pairedvector<int, double>& dmapold = cnode->FriData().GetDOld();

      std::set<int> snodes = cnode->FriData().GetSNodes();

      // check if there are entries in the old D map
      if (dmapold.size() < 1) dserror("Error in Interface::EvaluateRelMov(): No old D-Map!");

      std::map<int, double>::iterator colcurr;
      std::set<int>::iterator scurr;

      // loop over all slave nodes with an entry adjacent to this node
      for (scurr = snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("Cannot find node with gid %", gid);
        Node* csnode = dynamic_cast<Node*>(snode);

        double dik = dmap[csnode->Id()];
        double dikold = dmapold[csnode->Id()];

        std::map<int, double>::iterator mcurr;

        for (int dim = 0; dim < csnode->NumDof(); ++dim)
        {
          int locid = (xsmod->Map()).LID(csnode->Dofs()[dim]);
          jump[dim] -= (dik - dikold) * (*xsmod)[locid];
        }
      }  //  loop over adjacent slave nodes

      std::map<int, double>& mmap = cnode->MoData().GetM();
      std::map<int, double>& mmapold = cnode->FriData().GetMOld();

      const std::set<int>& mnodescurrent = cnode->FriData().GetMNodes();
      const std::set<int>& mnodesold = cnode->FriData().GetMNodesOld();

      // check if there are entries in the M map
      if (mmap.size() < 1) dserror("Error in Interface::EvaluateRelMov(): No M-Map!");

      // check if there are entries in the old M map
      if (mmapold.size() < 1) dserror("Error in Interface::EvaluateRelMov(): No old M-Map!");

      if (mnodesold.size() < 1) dserror("Error in Interface::EvaluateRelMov(): No old M-Set!");

      std::set<int> mnodes;
      std::set<int>::iterator mcurr;

      for (mcurr = mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr = mnodesold.begin(); mcurr != mnodesold.end(); mcurr++) mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("Cannot find node with gid %", gid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        double mik = mmap[cmnode->Id()];
        double mikold = mmapold[cmnode->Id()];

        std::map<int, double>::iterator mcurr;

        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          jump[dim] += (mik - mikold) * (cmnode->xspatial()[dim]);
        }
      }  //  loop over master nodes

      // write it to nodes
      for (int dim = 0; dim < Dim(); dim++) cnode->FriData().jump()[dim] = jump[dim];

      // linearization of jump vector

      // reset derivative map of jump
      for (int j = 0; j < (int)((cnode->FriData().GetDerivJump()).size()); ++j)
        (cnode->FriData().GetDerivJump())[j].clear();
      (cnode->FriData().GetDerivJump()).resize(0);

      /*** 01  **********************************************************/

      if (dmatrixmod == Teuchos::null)
      {
        // loop over according slave nodes
        for (scurr = snodes.begin(); scurr != snodes.end(); scurr++)
        {
          int gid = *scurr;
          DRT::Node* snode = idiscret_->gNode(gid);
          if (!snode) dserror("Cannot find node with gid %", gid);
          Node* csnode = dynamic_cast<Node*>(snode);

          double dik = dmap[csnode->Id()];
          double dikold = dmapold[csnode->Id()];

          for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
          {
            int col = csnode->Dofs()[dimrow];
            double val = -(dik - dikold);
            if (abs(val) > 1e-14) cnode->AddDerivJumpValue(dimrow, col, val);
          }
        }
      }
      // in the 3D quadratic case, the values are obtained from the
      // global matrices Dmod and Doldmod
      else
      {
        // loop over dimension of the node
        for (int dim = 0; dim < cnode->NumDof(); ++dim)
        {
          int NumEntries = 0;
          int NumEntriesOld = 0;
          std::vector<double> Values((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<int> Indices((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<double> ValuesOld((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<int> IndicesOld((dmatrixmod->EpetraMatrix())->MaxNumEntries());

          // row
          int row = cnode->Dofs()[dim];

          // extract entries of this row from matrix
          int err = (dmatrixmod->EpetraMatrix())
                        ->ExtractGlobalRowCopy(row, (dmatrixmod->EpetraMatrix())->MaxNumEntries(),
                            NumEntries, Values.data(), Indices.data());
          if (err) dserror("ExtractMyRowView failed: err=%d", err);

          int errold = (doldmod->EpetraMatrix())
                           ->ExtractGlobalRowCopy(row, (doldmod->EpetraMatrix())->MaxNumEntries(),
                               NumEntriesOld, ValuesOld.data(), IndicesOld.data());
          if (errold) dserror("ExtractMyRowView failed: err=%d", err);

          // loop over entries of this vector
          for (int j = 0; j < NumEntries; ++j)
          {
            double ValueOld = 0;
            bool found = false;

            // find value with the same index in vector of Dold
            for (int k = 0; k < NumEntriesOld; ++k)
            {
              if (Indices[k] == Indices[j])
              {
                ValueOld = ValuesOld[k];
                found = true;
                break;
              }
            }

            if (found == false or abs(ValueOld) < 1e-12)
              dserror("Error in EvaluareRelMov(): No old D value exists");

            // write to node
            cnode->AddDerivJumpValue(dim, Indices[j], (Values[j] - ValueOld));
          }
        }
      }

      /*** 02  **********************************************************/
      // loop over according master nodes
      for (mcurr = mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("Cannot find node with gid %", gid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        double mik = mmap[cmnode->Id()];
        double mikold = mmapold[cmnode->Id()];

        for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
        {
          int col = cmnode->Dofs()[dimrow];
          double val = (mik - mikold);
          if (abs(val) > 1e-14) cnode->AddDerivJumpValue(dimrow, col, val);
        }
      }

      /*** 03 ***********************************************************/
      // we need the Lin(D-matrix) entries of this node
      std::map<int, std::map<int, double>>& ddmap = cnode->Data().GetDerivD();
      std::map<int, std::map<int, double>>::iterator dscurr;

      // loop over all slave nodes in the DerivM-map of the stick slave node
      for (dscurr = ddmap.begin(); dscurr != ddmap.end(); ++dscurr)
      {
        int gid = dscurr->first;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("Cannot find node with gid %", gid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // compute entry of the current stick node / slave node pair
        std::map<int, double>& thisdmmap = cnode->Data().GetDerivD(gid);

        // loop over all entries of the current derivative map
        for (colcurr = thisdmmap.begin(); colcurr != thisdmmap.end(); ++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for (int dim = 0; dim < cnode->NumDof(); ++dim)
          {
            int locid = (xsmod->Map()).LID(csnode->Dofs()[dim]);
            double val = -colcurr->second * (*xsmod)[locid];
            if (abs(val) > 1e-14) cnode->AddDerivJumpValue(dim, col, val);
          }
        }
      }

      /*** 04 ***********************************************************/
      // we need the Lin(M-matrix) entries of this node
      std::map<int, std::map<int, double>>& dmmap = cnode->Data().GetDerivM();
      std::map<int, std::map<int, double>>::iterator dmcurr;

      // loop over all master nodes in the DerivM-map of the stick slave node
      for (dmcurr = dmmap.begin(); dmcurr != dmmap.end(); ++dmcurr)
      {
        int gid = dmcurr->first;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("Cannot find node with gid %", gid);
        Node* cmnode = dynamic_cast<Node*>(mnode);
        double* mxi = cmnode->xspatial();

        // compute entry of the current stick node / master node pair
        std::map<int, double>& thisdmmap = cnode->Data().GetDerivM(gid);

        // loop over all entries of the current derivative map
        for (colcurr = thisdmmap.begin(); colcurr != thisdmmap.end(); ++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for (int dimrow = 0; dimrow < cnode->NumDof(); ++dimrow)
          {
            double val = colcurr->second * mxi[dimrow];
            if (abs(val) > 1e-14) cnode->AddDerivJumpValue(dimrow, col, val);
          }
        }
      }

      if (constr_direction_ == INPAR::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); j++)
          if (cnode->DbcDofs()[j] == true)
          {
            cnode->FriData().jump()[j] = 0.;
            cnode->FriData().GetDerivJump()[j].clear();
          }

    }  // active nodes
  }    // loop over slave nodes
  return;
}

/*----------------------------------------------------------------------*
 |  calculate nodal distances (public)                     pfaller Jan15|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateDistances(const Teuchos::RCP<const Epetra_Vector>& vec,
    std::map<int, std::vector<double>>& mynormals,
    std::map<int, std::vector<CORE::GEN::pairedvector<int, double>>>& dmynormals,
    std::map<int, double>& mygap, std::map<int, std::map<int, double>>& dmygap)
{
  SetState(MORTAR::state_new_displacement, *vec);
  Initialize();

  // interface needs to be complete
  if (!Filled() && Comm().MyPID() == 0) dserror("FillComplete() not called on interface %", id_);

  // create an interpolator instance
  Teuchos::RCP<NTS::Interpolator> interpolator =
      Teuchos::rcp(new NTS::Interpolator(imortar_, Dim()));

  // create normals
  PreEvaluate(-1, -1);  // dummy values

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  Comm().Barrier();

  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("Cannot find slave element with gid %", gid1);
    CONTACT::Element* selement = dynamic_cast<CONTACT::Element*>(ele1);

    if (selement->MoData().NumSearchElements() < 1)
    {
      std::cout << "WARNING: No elements found!" << std::endl;
      continue;
    }

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized()) continue;

    // empty vector of master element pointers
    std::vector<CONTACT::Element*> melements;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("Cannot find master element with gid %", gid2);
      CONTACT::Element* melement = dynamic_cast<CONTACT::Element*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      melements.push_back(melement);
    }

    //**************************************************************
    //                loop over all Slave nodes
    //**************************************************************
    for (int snodes = 0; snodes < selement->NumNode(); ++snodes)
    {
      CONTACT::Node* mynode = dynamic_cast<CONTACT::Node*>(selement->Nodes()[snodes]);

      // skip this node if already considered
      if (mynode->HasProj()) continue;

      //                store node normals
      //**************************************************************
      //      int gid = snoderowmapbound_->GID(snodes);
      int gid = mynode->Id();

      int numdofs = mynode->NumDof();
      std::vector<double> temp(numdofs, 0.0);
      for (int kk = 0; kk < numdofs; kk++)
      {
        temp[kk] = mynode->MoData().n()[kk];
      }
      mynormals.insert(std::pair<int, std::vector<double>>(gid, temp));
      dmynormals.insert(std::pair<int, std::vector<CORE::GEN::pairedvector<int, double>>>(
          gid, mynode->Data().GetDerivN()));

      //**************************************************************
      double sxi[2] = {0.0, 0.0};

      if (selement->Shape() == CORE::FE::CellType::quad4 or
          selement->Shape() == CORE::FE::CellType::quad8 or
          selement->Shape() == CORE::FE::CellType::quad9)
      {
        // TODO (pfaller): switch case
        if (snodes == 0)
        {
          sxi[0] = -1;
          sxi[1] = -1;
        }
        else if (snodes == 1)
        {
          sxi[0] = 1;
          sxi[1] = -1;
        }
        else if (snodes == 2)
        {
          sxi[0] = 1;
          sxi[1] = 1;
        }
        else if (snodes == 3)
        {
          sxi[0] = -1;
          sxi[1] = 1;
        }
        else if (snodes == 4)
        {
          sxi[0] = 0;
          sxi[1] = -1;
        }
        else if (snodes == 5)
        {
          sxi[0] = 1;
          sxi[1] = 0;
        }
        else if (snodes == 6)
        {
          sxi[0] = 0;
          sxi[1] = 1;
        }
        else if (snodes == 7)
        {
          sxi[0] = -1;
          sxi[1] = 0;
        }
        else if (snodes == 8)
        {
          sxi[0] = 0;
          sxi[1] = 0;
        }
        else
          dserror("ERORR: wrong node LID");
      }
      else if (selement->Shape() == CORE::FE::CellType::tri3 or
               selement->Shape() == CORE::FE::CellType::tri6)
      {
        if (snodes == 0)
        {
          sxi[0] = 0;
          sxi[1] = 0;
        }
        else if (snodes == 1)
        {
          sxi[0] = 1;
          sxi[1] = 0;
        }
        else if (snodes == 2)
        {
          sxi[0] = 0;
          sxi[1] = 1;
        }
        else if (snodes == 3)
        {
          sxi[0] = 0.5;
          sxi[1] = 0;
        }
        else if (snodes == 4)
        {
          sxi[0] = 0.5;
          sxi[1] = 0.5;
        }
        else if (snodes == 5)
        {
          sxi[0] = 0;
          sxi[1] = 0.5;
        }
        else
          dserror("ERORR: wrong node LID");
      }
      else
      {
        dserror("Chosen element type not supported for NTS!");
      }

      //**************************************************************
      //                loop over all Master Elements
      //**************************************************************
      // create vectors to store projection information for several master elements in case
      // projection is not unique
      std::vector<double> gap_vec;
      std::vector<std::map<int, double>> dgap_vec;

      for (int nummaster = 0; nummaster < (int)melements.size(); ++nummaster)
      {
        // project Gauss point onto master element
        double mxi[2] = {0.0, 0.0};
        double projalpha = 0.0;
        bool is_projected =
            MORTAR::Projector::Impl(*selement, *melements[nummaster])
                ->ProjectGaussPoint3D(*selement, sxi, *melements[nummaster], mxi, projalpha);

        bool is_on_mele = true;

        // check GP projection
        CORE::FE::CellType dt = melements[nummaster]->Shape();
        const double tol = 1e-8;
        if (dt == CORE::FE::CellType::quad4 || dt == CORE::FE::CellType::quad8 ||
            dt == CORE::FE::CellType::quad9)
        {
          if (mxi[0] < -1.0 - tol || mxi[1] < -1.0 - tol || mxi[0] > 1.0 + tol ||
              mxi[1] > 1.0 + tol)
          {
            is_on_mele = false;
          }
        }
        else
        {
          if (mxi[0] < -tol || mxi[1] < -tol || mxi[0] > 1.0 + tol || mxi[1] > 1.0 + tol ||
              mxi[0] + mxi[1] > 1.0 + 2 * tol)
          {
            is_on_mele = false;
          }
        }

        // node on mele?
        if (is_on_mele && is_projected)
        {
          // store information of projection so that this node is not considered again
          mynode->HasProj() = true;

          int ndof = 3;
          int ncol = melements[nummaster]->NumNode();
          CORE::LINALG::SerialDenseVector mval(ncol);
          CORE::LINALG::SerialDenseMatrix mderiv(ncol, 2);
          melements[nummaster]->EvaluateShape(mxi, mval, mderiv, ncol, false);

          //          int linsize    = mynode->GetLinsize();
          int linsize = 100;
          double gpn[3] = {0.0, 0.0, 0.0};
          //**************************************************************

          // evalute the GP slave coordinate derivatives --> no entries
          std::vector<CORE::GEN::pairedvector<int, double>> dsxi(2, 0);
          std::vector<CORE::GEN::pairedvector<int, double>> dmxi(2, 4 * linsize + ncol * ndof);

          (*interpolator)
              .DerivXiGP3D(*selement, *melements[nummaster], sxi, mxi, dsxi, dmxi, projalpha);
          (*interpolator).nwGap3D(*mynode, *melements[nummaster], mval, mderiv, dmxi, gpn);

          // store linearization for node
          std::map<int, double> dgap = mynode->Data().GetDerivGnts();  // (dof,value)

          // store gap information
          gap_vec.push_back(mynode->Data().Getgnts());
          dgap_vec.push_back(dgap);

          // reset nodal weighted gap and derivative
          mynode->Data().Getgnts() = 1.0e12;
          (mynode->Data().GetDerivGnts()).clear();
        }  // End hit ele
      }    // End Loop over all Master Elements

      if (gap_vec.size() > 0)
      {
        // find projection with smallest absoluate value of gap
        std::vector<double>::iterator iter_min =
            std::min_element(gap_vec.begin(), gap_vec.end(), abs_compare);
        const int i_min = std::distance(gap_vec.begin(), iter_min);

        // save to map at GID
        mygap.insert(std::pair<int, double>(gid, gap_vec[i_min]));
        dmygap.insert(std::pair<int, std::map<int, double>>(gid, dgap_vec[i_min]));
      }
    }
  }

  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate L2 Norm of tangential contact conditions     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvaluateTangentNorm(double& cnormtan)
{
  // friction coefficient
  double frcoeff = InterfaceParams().get<double>("FRCOEFF");

  // loop over all slave row nodes on the current interface
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some information from node
    double* n = cnode->MoData().n();
    const int dim = cnode->NumDof();

    // tangential plane
    CORE::LINALG::SerialDenseMatrix tanplane(dim, dim);
    if (dim == 3)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);
      tanplane(0, 2) = -(n[0] * n[2]);
      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
      tanplane(1, 2) = -(n[1] * n[2]);

      tanplane(2, 0) = -(n[2] * n[0]);
      tanplane(2, 1) = -(n[2] * n[1]);
      tanplane(2, 2) = 1 - (n[2] * n[2]);
    }
    else if (dim == 2)
    {
      tanplane(0, 0) = 1 - (n[0] * n[0]);
      tanplane(0, 1) = -(n[0] * n[1]);

      tanplane(1, 0) = -(n[1] * n[0]);
      tanplane(1, 1) = 1 - (n[1] * n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // jump vector
    CORE::LINALG::SerialDenseMatrix jumpvec(dim, 1);
    for (int i = 0; i < dim; i++) jumpvec(i, 0) = cnode->FriData().jump()[i];

    // evaluate jump in tangential direction
    CORE::LINALG::SerialDenseMatrix jumptan(dim, 1);
    CORE::LINALG::multiply(jumptan, tanplane, jumpvec);

    // force vector
    CORE::LINALG::SerialDenseMatrix forcevec(dim, 1);
    for (int i = 0; i < dim; i++) forcevec(i, 0) = cnode->MoData().lm()[i];

    // evaluate force in normal direction
    double forcen = 0.0;
    for (int k = 0; k < dim; ++k) forcen += forcevec(k, 0) * n[k];

    // norm of constraint violation for stick nodes
    if (cnode->Active() == true and cnode->FriData().Slip() == false)
    {
      for (int j = 0; j < dim; j++) cnormtan += jumptan(j, 0) * jumptan(j, 0);
    }
    // norm of constraint violation for slip nodes
    else if (cnode->Active() == true and cnode->FriData().Slip() == true)
    {
      double part1 = 0.0;
      double jumpnorm = 0.0;

      for (int j = 0; j < dim; j++)
      {
        jumpnorm += jumptan(j, 0) * jumptan(j, 0);
        part1 += jumptan(j, 0) * forcevec(j, 0);
      }

      jumpnorm = sqrt(jumpnorm);
      cnormtan += (part1 - frcoeff * forcen * jumpnorm) * (part1 - frcoeff * forcen * jumpnorm);
    }
  }  // loop over slave nodes

  // get cnorm from all procs
  double sumcnormtanallprocs = 0.0;
  Comm().SumAll(&cnormtan, &sumcnormtanallprocs, 1);
  cnormtan = sumcnormtanallprocs;

  return;
}


/*----------------------------------------------------------------------*
 |  Update active set and check for convergence (public)      popp 06/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::UpdateActiveSetSemiSmooth()
{
  // get input parameter ftype
  INPAR::CONTACT::FrictionType ftype =
      CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(), "FRICTION");

  // this is the complementarity parameter we use for the decision.
  // it might be scaled with a mesh-size dependent factor
  double cn = 0.;
  double ct = 0.;

  // assume that active set has converged and check for opposite
  int localcheck = true;

  // loop over all slave nodes on the current interface
  for (int j = 0; j < SlaveRowNodes()->NumMyElements(); ++j)
  {
    int gid = SlaveRowNodes()->GID(j);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);

    Node* cnode = dynamic_cast<Node*>(node);

    cn = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];
    if (friction_) ct = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

    // get weighted gap
    double wgap = cnode->Data().Getg();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k = 0; k < 3; ++k)
    {
      nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];
    }

    // friction
    std::vector<double> tz(Dim() - 1, 0);
    std::vector<double> tjump(Dim() - 1, 0);
    double euclidean = 0.0;

    if (friction_)
    {
      // static cast
      FriNode* frinode = dynamic_cast<FriNode*>(cnode);

      // compute tangential parts and of Lagrange multiplier and incremental jumps
      for (int i = 0; i < Dim(); ++i)
      {
        tz[0] += frinode->Data().txi()[i] * frinode->MoData().lm()[i];
        if (Dim() == 3) tz[1] += frinode->Data().teta()[i] * frinode->MoData().lm()[i];

        if (CORE::UTILS::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == false)
        {
          tjump[0] += frinode->Data().txi()[i] * frinode->FriData().jump()[i];
          if (Dim() == 3) tjump[1] += frinode->Data().teta()[i] * frinode->FriData().jump()[i];
        }
      }

      if (CORE::UTILS::IntegralValue<int>(InterfaceParams(), "GP_SLIP_INCR") == true)
      {
        tjump[0] = frinode->FriData().jump_var()[0];
        if (Dim() == 3) tjump[1] = frinode->FriData().jump_var()[1];
      }

      // evaluate euclidean norm |tz+ct.tjump|
      std::vector<double> sum(Dim() - 1, 0);
      sum[0] = tz[0] + ct * tjump[0];
      if (Dim() == 3) sum[1] = tz[1] + ct * tjump[1];
      if (Dim() == 2) euclidean = abs(sum[0]);
      if (Dim() == 3) euclidean = sqrt(sum[0] * sum[0] + sum[1] * sum[1]);
    }

    // adhesion
    double adhbound = 0.0;
    if (CORE::UTILS::IntegralValue<INPAR::CONTACT::AdhesionType>(InterfaceParams(), "ADHESION") ==
        INPAR::CONTACT::adhesion_bound)
      adhbound = InterfaceParams().get<double>("ADHESION_BOUND");

    // check nodes of inactive set *************************************
    if (cnode->Active() == false)
    {
      // check for penetration and/or tensile contact forces
      if (nz - cn * wgap >
          0)  // no averaging of Lagrange multipliers
              // if ((0.5*nz+0.5*nzold) - cn*wgap > 0) // averaging of Lagrange multipliers
      {
        cnode->Active() = true;
        localcheck = false;

        // friction
        if (friction_)
        {
          // nodes coming into contact
          dynamic_cast<FriNode*>(cnode)->FriData().Slip() = true;
        }
      }
    }

    // check nodes of active set ***************************************
    else
    {
      // adhesion modification
      nz += adhbound;

      // check for tensile contact forces and/or penetration
      if (nz - cn * wgap <=
          0)  // no averaging of Lagrange multipliers
              // if ((0.5*nz+0.5*nzold) - cn*wgap <= 0) // averaging of Lagrange multipliers
      {
        cnode->Active() = false;
        localcheck = false;

        // friction
        if (friction_) dynamic_cast<FriNode*>(cnode)->FriData().Slip() = false;
      }

      // only do something for friction
      else
      {
        // friction tresca
        if (ftype == INPAR::CONTACT::friction_tresca)
        {
          FriNode* frinode = dynamic_cast<FriNode*>(cnode);

          // CAREFUL: friction bound is now interface-local (popp 08/2012)
          double frbound = InterfaceParams().get<double>("FRBOUND");

          if (frinode->FriData().Slip() == false)
          {
            // check (euclidean)-frbound <= 0
            if (euclidean - frbound <= 0)
            {
            }
            // do nothing (stick was correct)
            else
            {
              frinode->FriData().Slip() = true;
              localcheck = false;
            }
          }
          else
          {
            // check (euclidean)-frbound > 0
            if (euclidean - frbound > 0)
            {
            }
            // do nothing (slip was correct)
            else
            {
              frinode->FriData().Slip() = false;
              localcheck = false;
            }
          }
        }  // if (fytpe=="tresca")

        // friction coulomb
        if (ftype == INPAR::CONTACT::friction_coulomb)
        {
          FriNode* frinode = dynamic_cast<FriNode*>(cnode);

          // CAREFUL: friction coefficient is now interface-local (popp 08/2012)
          double frcoeff = frinode->FrCoeff(InterfaceParams().get<double>("FRCOEFF"));
          double frbound;
          static const bool regularization =
              CORE::UTILS::IntegralValue<int>(InterfaceParams(), "REGULARIZED_NORMAL_CONTACT");
          if (!regularization)
            frbound = frcoeff * (nz - cn * wgap);
          else
          {
            static const double k = 1. / InterfaceParams().get<double>("REGULARIZATION_STIFFNESS");
            static const double gmax = InterfaceParams().get<double>("REGULARIZATION_THICKNESS");
            if (cnode->MoData().GetD().size() != 1)
              dserror(
                  "we need to have a D-value for active contact nodes\nAnd exactly one due to "
                  "biorthogonality");
            double dval = cnode->MoData().GetD().at(cnode->Id());
            const double gLM = gmax * (1. - exp(-k / gmax * nz));
            frbound = frcoeff * std::max(0., nz - cn * (wgap + dval * gLM));
          }

          if (frinode->FriData().Slip() == false)
          {
            // check (euclidean)-frbound <= 0
            if (euclidean - frbound <= 1e-10)
            {
            }
            // do nothing (stick was correct)
            else
            {
              frinode->FriData().Slip() = true;
              localcheck = false;
            }
          }
          else
          {
            // check (euclidean)-frbound > 0
            if (euclidean - frbound > -1e-10)
            {
            }
            // do nothing (slip was correct)
            else
            {
              frinode->FriData().Slip() = false;
              localcheck = false;
            }
          }
        }  // if (ftype == INPAR::CONTACT::friction_coulomb)
      }    // if (nz - cn*wgap <= 0)
    }      // if (cnode->Active()==false)
  }        // loop over all slave nodes

  // broadcast convergence status among processors
  int convcheck = 0;
  Comm().MinAll(&localcheck, &convcheck, 1);

  return convcheck;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::UpdateActiveSetInitialStatus() const
{
  // List of GIDs of all my slave row nodes
  int* my_slave_row_node_gids = SlaveRowNodes()->MyGlobalElements();

  // loop over all slave nodes on the current interface
  for (int j = 0; j < SlaveRowNodes()->NumMyElements(); ++j)
  {
    // Grab the current slave node
    const int gid = my_slave_row_node_gids[j];
    Node* cnode = dynamic_cast<Node*>(Discret().gNode(gid));
    if (!cnode) dserror("Cannot find node with gid %", gid);

    SetNodeInitiallyActive(*cnode);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::BuildActiveSet(bool init)
{
  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mynodegidsInactive(0);
  std::vector<int> mydofgids(0);
  std::vector<int> mydofgidsInactive(0);
  std::vector<int> myslipnodegids(0);
  std::vector<int> myslipdofgids(0);
  std::vector<int> mymnodegids(0);
  std::vector<int> mymdofgids(0);

  // loop over all slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // *******************************************************************
    // This is given by the Node member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding Nodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (init)
    {
      // flag for initialization of init active nodes with nodal gaps
      bool initcontactbygap =
          CORE::UTILS::IntegralValue<int>(InterfaceParams(), "INITCONTACTBYGAP");
      // value
      double initcontactval = InterfaceParams().get<double>("INITCONTACTGAPVALUE");

      // Either init contact by definition or by gap
      if (cnode->IsInitActive() and initcontactbygap)
        dserror("Init contact either by definition in condition or by gap!");

      // check if node is initially active or, if initialization with nodal, gap,
      // the gap is smaller than the prescribed value
      if (cnode->IsInitActive() or (initcontactbygap and cnode->Data().Getg() < initcontactval)
          //          sqrt(cnode->X()[0]*cnode->X()[0]+cnode->X()[1]*cnode->X()[1])>65.
          //          (sqrt(cnode->X()[0]*cnode->X()[0]+cnode->X()[1]*cnode->X()[1])>12. &&
          //          abs(cnode->X()[2]-10.5)<1.e-12)
          //          || abs(cnode->X()[2]-0.5)<1.e-12
      )
      {
        /*
              // **********************************************************************************
              // ***                     CONSISTENT INITIALIZATION STATE                        ***
              // **********************************************************************************
              // hiermeier 08/2013
              //
              // The consistent nonlinear complementarity function C_tau is for the special
           initialization case
              //
              //                zn == 0   and   gap == 0
              //
              // in conjunction with the Coulomb friction model no longer unique and a special
           treatment is
              // necessary. For this purpose we identify the critical nodes here and set the
           slip-state to true.
              // In the next step the tangential part of the Lagrange multiplier vector will be set
           to zero. Hence,
              // we treat these nodes as frictionless nodes! (see AssembleLinSlip)
              INPAR::CONTACT::FrictionType ftype =
                  CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(InterfaceParams(),"FRICTION");
              if (ftype == INPAR::CONTACT::friction_coulomb)
              {
              dynamic_cast<FriNode*>(cnode)->FriData().InconInit() = true;
              myslipnodegids.push_back(cnode->Id());

              for (int j=0;j<numdof;++j)
                myslipdofgids.push_back(cnode->Dofs()[j]);
              }

        */
        cnode->Active() = true;
        mynodegids.push_back(cnode->Id());

        for (int j = 0; j < numdof; ++j) mydofgids.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is initially in slip state
      if (friction_)
      {
        // do nothing: we always assume STICK at t=0
      }
    }

    // *******************************************************************
    // RE-BUILDING OF THE ACTIVE SET
    // *******************************************************************
    else
    {
      // check if node is active
      if (cnode->Active())
      {
        mynodegids.push_back(cnode->Id());

        for (int j = 0; j < numdof; ++j) mydofgids.push_back(cnode->Dofs()[j]);
      }
      else
      {
        mynodegidsInactive.push_back(cnode->Id());

        for (int j = 0; j < numdof; ++j) mydofgidsInactive.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is in slip state
      if (friction_)
      {
        if (dynamic_cast<FriNode*>(cnode)->FriData().Slip())
        {
          myslipnodegids.push_back(cnode->Id());

          for (int j = 0; j < numdof; ++j) myslipdofgids.push_back(cnode->Dofs()[j]);
        }
      }
    }
  }

  // create active node map and active dof map
  activenodes_ = CORE::LINALG::CreateMap(mynodegids, Comm());
  activedofs_ = CORE::LINALG::CreateMap(mydofgids, Comm());
  inactivenodes_ = CORE::LINALG::CreateMap(mynodegidsInactive, Comm());
  inactivedofs_ = CORE::LINALG::CreateMap(mydofgidsInactive, Comm());

  if (friction_)
  {
    // create slip node map and slip dof map
    slipnodes_ = CORE::LINALG::CreateMap(myslipnodegids, Comm());
    slipdofs_ = CORE::LINALG::CreateMap(myslipdofgids, Comm());
  }

  // split active dofs and slip dofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  split active dofs into Ndofs, Tdofs and slipTdofs         popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Interface::SplitActiveDofs()
{
  // get out of here if active set is empty
  if (activenodes_ == Teuchos::null)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return true;
  }

  else if (activenodes_->NumGlobalElements() == 0)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return true;
  }

  // define local variables
  int countN = 0;
  int countT = 0;
  std::vector<int> myNgids(activenodes_->NumMyElements());
  std::vector<int> myTgids((Dim() - 1) * activenodes_->NumMyElements());

  // dimension check
  double dimcheck = (activedofs_->NumGlobalElements()) / (activenodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all active row nodes
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    const int numdof = cnode->NumDof();

    // add first dof to Nmap
    myNgids[countN] = cnode->Dofs()[0];
    ++countN;

    // add remaining dofs to Tmap
    for (int j = 1; j < numdof; ++j)
    {
      myTgids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNgids.resize(countN);
  myTgids.resize(countT);

  // communicate countN and countT among procs
  int gcountN, gcountT;
  Comm().SumAll(&countN, &gcountN, 1);
  Comm().SumAll(&countT, &gcountT, 1);

  // check global dimensions
  if ((gcountN + gcountT) != activedofs_->NumGlobalElements())
    dserror("SplitActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = Teuchos::rcp(new Epetra_Map(gcountN, countN, myNgids.data(), 0, Comm()));
  activet_ = Teuchos::rcp(new Epetra_Map(gcountT, countT, myTgids.data(), 0, Comm()));

  // *******************************************************************
  // FRICTION - EXTRACTING TANGENTIAL DOFS FROM SLIP DOFS
  // *******************************************************************

  // get out of here if there is no friction
  if (friction_ == false) return true;

  // get out of here if slip set is empty
  if (slipnodes_ == Teuchos::null)
  {
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return true;
  }

  if (slipnodes_->NumGlobalElements() == 0)
  {
    slipt_ = Teuchos::rcp(new Epetra_Map(0, 0, Comm()));
    return true;
  }

  // define local variables
  int countslipT = 0;
  std::vector<int> myslipTgids((Dim() - 1) * slipnodes_->NumMyElements());

  // dimension check
  dimcheck = (slipdofs_->NumGlobalElements()) / (slipnodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slip row nodes
  for (int i = 0; i < slipnodes_->NumMyElements(); ++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    const int numdof = cnode->NumDof();

    // add dofs to slipTmap
    for (int j = 1; j < numdof; ++j)
    {
      myslipTgids[countslipT] = cnode->Dofs()[j];
      ++countslipT;
    }
  }

  // resize the temporary vectors
  myslipTgids.resize(countslipT);

  // communicate countslipT among procs
  int gcountslipT;
  Comm().SumAll(&countslipT, &gcountslipT, 1);

  // create Tslipmap objects
  slipt_ = Teuchos::rcp(new Epetra_Map(gcountslipT, countslipT, myslipTgids.data(), 0, Comm()));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::GetForceOfNode(CORE::LINALG::Matrix<3, 1>& nodal_force,
    const Epetra_Vector& force, const DRT::Node& node) const
{
  std::vector<int> dofs;
  idiscret_->Dof(&node, dofs);

  // reset nodal force vector
  std::fill(nodal_force.A(), nodal_force.A() + 3, 0.0);
  const double* f_vals = force.Values();

  if (dofs.size() > 3) dserror("The interface node seems to have more than 3 DOFs!");

  for (unsigned i = 0; i < dofs.size(); ++i)
  {
    const int dof = dofs[i];
    const int f_lid = force.Map().LID(dof);
    if (f_lid == -1) dserror("Couldn't find the nodal DOF %d in the global force vector!", dof);

    nodal_force(i, 0) = f_vals[f_lid];
  }
}

/*----------------------------------------------------------------------*
 |  Calculate angular interface moments                  hiermeier 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::EvalResultantMoment(const Epetra_Vector& fs, const Epetra_Vector& fm,
    CORE::LINALG::SerialDenseMatrix* conservation_data_ptr) const
{
  double lresMoSl[3] = {0.0, 0.0, 0.0};  // local slave moment
  double gresMoSl[3] = {0.0, 0.0, 0.0};  // global slave moment
  double lresMoMa[3] = {0.0, 0.0, 0.0};  // local master momemnt
  double gresMoMa[3] = {0.0, 0.0, 0.0};  // global master moment
  double gbalMo[3] = {0.0, 0.0, 0.0};    // global moment balance

  double lresFSl[3] = {0.0, 0.0, 0.0};  // local slave force
  double gresFSl[3] = {0.0, 0.0, 0.0};  // global slave force
  double lresFMa[3] = {0.0, 0.0, 0.0};  // local master force
  double gresFMa[3] = {0.0, 0.0, 0.0};  // global master force
  double gbalF[3] = {0.0, 0.0, 0.0};    // global force balance

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  CORE::LINALG::Matrix<3, 1> nforce(true);
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Node* snode = dynamic_cast<Node*>(idiscret_->gNode(gid));

    GetForceOfNode(nforce, fs, *snode);

    if (Dim() == 2)
      lresMoSl[2] += snode->xspatial()[0] * nforce(1, 0) - snode->xspatial()[1] * nforce(0, 0);
    else
    {
      lresMoSl[0] += snode->xspatial()[1] * nforce(2, 0) - snode->xspatial()[2] * nforce(1, 0);
      lresMoSl[1] += snode->xspatial()[2] * nforce(0, 0) - snode->xspatial()[0] * nforce(2, 0);
      lresMoSl[2] += snode->xspatial()[0] * nforce(1, 0) - snode->xspatial()[1] * nforce(0, 0);
    }
    for (unsigned k = 0; k < 3; ++k) lresFSl[k] += nforce(k, 0);
  }
  Comm().SumAll(lresMoSl, gresMoSl, 3);
  Comm().SumAll(lresFSl, gresFSl, 3);

  // loop over proc's master nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
  {
    int gid = mnoderowmap_->GID(i);
    Node* mnode = dynamic_cast<Node*>(idiscret_->gNode(gid));

    GetForceOfNode(nforce, fm, *mnode);

    if (Dim() == 2)
      lresMoMa[2] += mnode->xspatial()[0] * nforce(1, 0) - mnode->xspatial()[1] * nforce(0, 0);
    else
    {
      lresMoMa[0] += mnode->xspatial()[1] * nforce(2, 0) - mnode->xspatial()[2] * nforce(1, 0);
      lresMoMa[1] += mnode->xspatial()[2] * nforce(0, 0) - mnode->xspatial()[0] * nforce(2, 0);
      lresMoMa[2] += mnode->xspatial()[0] * nforce(1, 0) - mnode->xspatial()[1] * nforce(0, 0);
    }
    for (unsigned k = 0; k < 3; ++k) lresFMa[k] += nforce(k, 0);
  }
  Comm().SumAll(lresMoMa, gresMoMa, 3);
  Comm().SumAll(lresFMa, gresFMa, 3);

  for (int d = 0; d < 3; ++d)
  {
    gbalMo[d] = gresMoSl[d] + gresMoMa[d];
    gbalF[d] = gresFSl[d] + gresFMa[d];
  }
  if (Comm().MyPID() == 0)
  {
    std::cout << "SLAVE:   "
              << " [" << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[0]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[1]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[2]
              << "]" << std::endl;
    std::cout << "Master:  "
              << " [" << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[0]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[1]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[2]
              << "]" << std::endl;
    std::cout << "Balance: "
              << " [" << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[0]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[1]
              << ", " << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[2]
              << "]" << std::endl;
  }

  if (conservation_data_ptr)
  {
    CORE::LINALG::SerialDenseMatrix& conservation_data = *conservation_data_ptr;
    conservation_data.putScalar(0.0);
    if (conservation_data.numRows() < 18) dserror("conservation_data length is too short!");

    std::copy(gresFSl, gresFSl + 3, conservation_data.values());
    std::copy(gresFMa, gresFMa + 3, conservation_data.values() + 3);
    std::copy(gbalF, gbalF + 3, conservation_data.values() + 6);

    std::copy(gresMoSl, gresMoSl + 3, conservation_data.values() + 9);
    std::copy(gresMoMa, gresMoMa + 3, conservation_data.values() + 12);
    std::copy(gbalMo, gbalMo + 3, conservation_data.values() + 15);
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::Interface& CONTACT::Interface::GetMaSharingRefInterface() const
{
  return dynamic_cast<const Interface&>(interfaceData_->GetMaSharingRefInterface());
}


/*----------------------------------------------------------------------*
 | Store nodal quant. to old ones (last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::Interface::StoreToOld(MORTAR::StrategyBase::QuantityType type)
{
  // loop over all slave row nodes on the current interface
  for (int j = 0; j < SlaveColNodes()->NumMyElements(); ++j)
  {
    int gid = SlaveColNodes()->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("Cannot find node with gid %", gid);

    switch (type)
    {
      case MORTAR::StrategyBase::dm:
      {
        // store D and M entries
        dynamic_cast<FriNode*>(node)->StoreDMOld();
        break;
      }
      case MORTAR::StrategyBase::pentrac:
      {
        // store penalty tractions to old ones
        dynamic_cast<FriNode*>(node)->StoreTracOld();
        break;
      }
      case MORTAR::StrategyBase::n_old:
      {
        dynamic_cast<Node*>(node)->StoreOldNormal();
        break;
      }
      case MORTAR::StrategyBase::activeold:
      {
        dynamic_cast<Node*>(node)->Data().ActiveOld() = dynamic_cast<Node*>(node)->Active();
        break;
      }
      default:
        dserror("StoreDMToNodes: Unknown state std::string variable!");
        break;
    }  // switch
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Interface::UpdateSelfContactLagMultSet(
    const Epetra_Map& gref_lmmap, const Epetra_Map& gref_smmap)
{
  if (gref_lmmap.NumMyElements() != gref_smmap.NumMyElements()) dserror("Size mismatch!");

  const int num_sgids = sdofrowmap_->NumMyElements();
  const int* sgids = sdofrowmap_->MyGlobalElements();
  const int* ref_lmgids = gref_lmmap.MyGlobalElements();

  std::vector<int> lmdofs;
  lmdofs.reserve(num_sgids);

  for (int i = 0; i < num_sgids; ++i)
  {
    const int sgid = sgids[i];
    const int ref_lid = gref_smmap.LID(sgid);
    if (ref_lid == -1)
      dserror(
          "Couldn't find the current slave gid #%d in the reference self "
          "contact slave master map.",
          sgid);
    lmdofs.push_back(ref_lmgids[ref_lid]);
  }

  lmdofmap_ =
      Teuchos::rcp(new Epetra_Map(-1, static_cast<int>(lmdofs.size()), lmdofs.data(), 0, Comm()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::SetNodeInitiallyActive(CONTACT::Node& cnode) const
{
  static const bool init_contact_by_gap =
      CORE::UTILS::IntegralValue<int>(InterfaceParams(), "INITCONTACTBYGAP");

  const bool node_init_active = cnode.IsInitActive();

  // Either init contact by definition or by gap
  if (node_init_active and init_contact_by_gap)
    dserror("Init contact either by definition in condition or by gap!");
  else if (node_init_active)
    cnode.Active() = true;
  else if (init_contact_by_gap)
    SetNodeInitiallyActiveByGap(cnode);

#ifdef BACI_DEBUG
  if (node_init_active)
    std::cout << "Node #" << std::setw(5) << cnode.Id()
              << " is set initially active via the condition line.\n";
  else if (init_contact_by_gap)
    std::cout << "Node #" << std::setw(5) << cnode.Id() << " is set initially active by gap.\n";
#endif

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::SetNodeInitiallyActiveByGap(CONTACT::Node& cnode) const
{
  static const double initcontactval = InterfaceParams().get<double>("INITCONTACTGAPVALUE");

  if (cnode.Data().Getg() < initcontactval) cnode.Active() = true;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::SetConditionSpecificParameters()
{
  if (InterfaceParams().isSublist("ContactS2ICoupling"))
  {
    // read interface parameters and set them to the scatra boundary parameter class
    auto& s2icouplinglist = InterfaceParams().sublist("ContactS2ICoupling", true);
    DRT::ELEMENTS::ScaTraEleParameterBoundary::Instance("scatra")->SetParameters(s2icouplinglist);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Interface::PostprocessQuantities(const Teuchos::ParameterList& outputParams)
{
  using Teuchos::RCP;

  // Check if the given parameter list contains all required data to be written to output
  {
    // Vector with names of all required parameter entries
    std::vector<std::string> requiredEntries;
    requiredEntries.push_back("step");
    requiredEntries.push_back("time");
    requiredEntries.push_back("displacement");
    requiredEntries.push_back("interface traction");
    requiredEntries.push_back("slave forces");
    requiredEntries.push_back("master forces");

    CheckOutputList(outputParams, requiredEntries);
  }

  // Get the discretization writer and get ready for writing
  RCP<IO::DiscretizationWriter> writer = idiscret_->Writer();

  // Get output for this time step started
  {
    const int step = outputParams.get<int>("step");
    const double time = outputParams.get<double>("time");

    writer->ClearMapCache();
    writer->WriteMesh(step, time);
    writer->NewStep(step, time);
  }

  /* Write interface displacement
   *
   * The interface displacement has been handed in via the parameter list outParams.
   * Grab it from there, then use CORE::LINALG::Export() to extract the interface
   * portion from the global displacement vector. Finally, write the interface
   * portion using this interfaces' discretization writer.
   */
  {
    // Get full displacement vector and extract interface displacement
    RCP<const Epetra_Vector> disp = outputParams.get<RCP<const Epetra_Vector>>("displacement");
    RCP<Epetra_Vector> iDisp = CORE::LINALG::CreateVector(*idiscret_->DofRowMap());
    CORE::LINALG::Export(*disp, *iDisp);

    // Write the interface displacement field
    writer->WriteVector("displacement", iDisp, IO::VectorType::dofvector);
  }

  // Write Lagrange multiplier field
  {
    // Get full Lagrange multiplier vector and extract values of this interface
    RCP<const Epetra_Vector> lagMult =
        outputParams.get<RCP<const Epetra_Vector>>("interface traction");
    RCP<Epetra_Vector> iLagMult = CORE::LINALG::CreateVector(*idiscret_->DofRowMap());
    CORE::LINALG::Export(*lagMult, *iLagMult);

    // Write this interface's Lagrange multiplier field
    writer->WriteVector("interfacetraction", iLagMult, IO::VectorType::dofvector);
  }

  // Write normal contact stress
  {
    // Get values from parameter list and export to interface DofRowMap
    RCP<const Epetra_Vector> normalStresses =
        outputParams.get<RCP<const Epetra_Vector>>("norcontactstress");
    RCP<Epetra_Vector> iNormalStresses = CORE::LINALG::CreateVector(*idiscret_->DofRowMap());
    CORE::LINALG::Export(*normalStresses, *iNormalStresses);

    // Write this interface's normal contact stress field
    writer->WriteVector("norcontactstress", iNormalStresses, IO::VectorType::dofvector);
  }

  // Write tangential contact stress
  {
    // Get values from parameter list and export to interface DofRowMap
    RCP<const Epetra_Vector> tangentialStresses =
        outputParams.get<RCP<const Epetra_Vector>>("tancontactstress");
    RCP<Epetra_Vector> iTangentialStresses = CORE::LINALG::CreateVector(*idiscret_->DofRowMap());
    CORE::LINALG::Export(*tangentialStresses, *iTangentialStresses);

    // Write this interface's normal contact stress field
    writer->WriteVector("tancontactstress", iTangentialStresses, IO::VectorType::dofvector);
  }

  // Write nodal forces of slave side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> slaveforces =
        outputParams.get<RCP<const Epetra_Vector>>("slave forces");
    RCP<Epetra_Vector> forces = CORE::LINALG::CreateVector(*idiscret_->DofRowMap());
    CORE::LINALG::Export(*slaveforces, *forces);

    // Write to output
    writer->WriteVector("slaveforces", forces, IO::VectorType::dofvector);
  }

  // Write nodal forces of master side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> masterforces =
        outputParams.get<RCP<const Epetra_Vector>>("master forces");
    RCP<Epetra_Vector> forces = CORE::LINALG::CreateVector(*idiscret_->DofRowMap());
    CORE::LINALG::Export(*masterforces, *forces);

    // Write to output
    writer->WriteVector("masterforces", forces, IO::VectorType::dofvector);
  }


  // Nodes: node-based vector with '0' at slave nodes and '1' at master nodes
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*mnoderowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> nodeRowMap = CORE::LINALG::MergeMap(snoderowmap_, mnoderowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = CORE::LINALG::CreateVector(*nodeRowMap, true);
    CORE::LINALG::Export(*masterVec, *masterSlaveVec);

    writer->WriteVector("slavemasternodes", masterSlaveVec, IO::VectorType::nodevector);
  }

  // Write active set
  {
    // evaluate active set and slip set
    RCP<Epetra_Vector> activeset = Teuchos::rcp(new Epetra_Vector(*activenodes_));
    activeset->PutScalar(1.0);

    if (IsFriction())
    {
      RCP<Epetra_Vector> slipset = Teuchos::rcp(new Epetra_Vector(*slipnodes_));
      slipset->PutScalar(1.0);
      RCP<Epetra_Vector> slipsetexp = Teuchos::rcp(new Epetra_Vector(*activenodes_));
      CORE::LINALG::Export(*slipset, *slipsetexp);
      activeset->Update(1.0, *slipsetexp, 1.0);
    }

    // export to interface node row map
    RCP<Epetra_Vector> activesetexp = Teuchos::rcp(new Epetra_Vector(*(idiscret_->NodeRowMap())));
    CORE::LINALG::Export(*activeset, *activesetexp);

    writer->WriteVector("activeset", activesetexp, IO::VectorType::nodevector);
  }

  // Elements: element-based vector with '0' at slave elements and '1' at master elements
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*melerowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> eleRowMap = CORE::LINALG::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = CORE::LINALG::CreateVector(*eleRowMap, true);
    CORE::LINALG::Export(*masterVec, *masterSlaveVec);

    writer->WriteVector("slavemasterelements", masterSlaveVec, IO::VectorType::elementvector);
  }

  // Write element owners
  {
    RCP<const Epetra_Map> eleRowMap = CORE::LINALG::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> owner = CORE::LINALG::CreateVector(*eleRowMap);

    for (int i = 0; i < idiscret_->ElementRowMap()->NumMyElements(); ++i)
      (*owner)[i] = idiscret_->lRowElement(i)->Owner();

    writer->WriteVector("Owner", owner, IO::VectorType::elementvector);
  }
}

FOUR_C_NAMESPACE_CLOSE
