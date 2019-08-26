/*----------------------------------------------------------------------*/
/*! \file
\brief One mortar coupling interface

\level 1

\maintainer Matthias Mayr
*/
/*-----------------------------------------------------------------------*/

#include "mortar_interface.H"
#include "mortar_node.H"
#include "mortar_element.H"
#include "mortar_integrator.H"
#include "mortar_coupling2d.H"
#include "mortar_coupling3d.H"
#include "mortar_coupling3d_classes.H"
#include "mortar_dofset.H"
#include "mortar_binarytree.H"
#include "mortar_defines.H"
#include "mortar_projector.H"
#include "mortar_utils.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_utils_rebalancing.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>
#include <Epetra_SerialComm.h>

#include "../drt_binstrategy/binning_strategy.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_poroelast/poroelast_utils.H"

#include "../drt_contact/contact_interpolator.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
MORTAR::InterfaceDataContainer::InterfaceDataContainer()
    : id_(-1),
      comm_(NULL),
      lcomm_(Teuchos::null),
      redistributed_(false),
      idiscret_(Teuchos::null),
      dim_(-1),
      imortar_(Teuchos::ParameterList()),
      shapefcn_(INPAR::MORTAR::shape_undefined),
      quadslave_(false),
      redundant_(INPAR::MORTAR::redundant_none),
      oldnodecolmap_(Teuchos::null),
      oldelecolmap_(Teuchos::null),
      snoderowmap_(Teuchos::null),
      snodecolmap_(Teuchos::null),
      mnoderowmap_(Teuchos::null),
      snoderowmapbound_(Teuchos::null),
      snodecolmapbound_(Teuchos::null),
      mnoderowmapnobound_(Teuchos::null),
      mnodecolmapnobound_(Teuchos::null),
      selerowmap_(Teuchos::null),
      selecolmap_(Teuchos::null),
      melerowmap_(Teuchos::null),
      melecolmap_(Teuchos::null),
      sdofrowmap_(Teuchos::null),
      sdofcolmap_(Teuchos::null),
      mdofrowmap_(Teuchos::null),
      mdofcolmap_(Teuchos::null),
      psdofrowmap_(Teuchos::null),
      plmdofmap_(Teuchos::null),
      lmdofmap_(Teuchos::null),
      maxdofglobal_(-1),
      searchalgo_(INPAR::MORTAR::search_binarytree),
      binarytree_(Teuchos::null),
      searchparam_(-1.0),
      searchuseauxpos_(false),
      inttime_interface_(0.0),
      nurbs_(false),
      poro_(false),
      ehl_(false),
      isinit_(false)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
MORTAR::MortarInterface::MortarInterface(
    const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData)
    : interfaceData_(interfaceData),
      id_(interfaceData_->Id()),
      comm_(interfaceData_->CommPtr()),
      lcomm_(interfaceData_->lComm()),
      procmap_(interfaceData_->ProcMap()),
      redistributed_(interfaceData_->IsRedistributed()),
      idiscret_(interfaceData_->IDiscret()),
      dim_(interfaceData_->Dim()),
      imortar_(interfaceData_->IMortar()),
      shapefcn_(interfaceData_->ShapeFcn()),
      quadslave_(interfaceData_->IsQuadSlave()),
      redundant_(interfaceData_->RedundantStorage()),
      oldnodecolmap_(interfaceData_->OldNodeColMap()),
      oldelecolmap_(interfaceData_->OldEleColMap()),
      snoderowmap_(interfaceData_->SNodeRowMap()),
      snodecolmap_(interfaceData_->SNodeColMap()),
      mnoderowmap_(interfaceData_->MNodeRowMap()),
      mnodecolmap_(interfaceData_->MNodeColMap()),
      snoderowmapbound_(interfaceData_->SNodeRowMapBound()),
      snodecolmapbound_(interfaceData_->SNodeColMapBound()),
      mnoderowmapnobound_(interfaceData_->MNodeRowMapNoBound()),
      mnodecolmapnobound_(interfaceData_->MNodeColMapNoBound()),
      selerowmap_(interfaceData_->SEleRowMap()),
      selecolmap_(interfaceData_->SEleColMap()),
      melerowmap_(interfaceData_->MEleRowMap()),
      melecolmap_(interfaceData_->MEleColMap()),
      sdofrowmap_(interfaceData_->SDofRowMap()),
      sdofcolmap_(interfaceData_->SDofColMap()),
      mdofrowmap_(interfaceData_->MDofRowMap()),
      mdofcolmap_(interfaceData_->MDofColMap()),
      psdofrowmap_(interfaceData_->PSDofRowMap()),
      plmdofmap_(interfaceData_->PLmDofRowMap()),
      lmdofmap_(interfaceData_->LmDofRowMap()),
      maxdofglobal_(interfaceData_->MaxDofGlobal()),
      searchalgo_(interfaceData_->SearchAlgorithm()),
      binarytree_(interfaceData_->BinaryTree()),
      searchparam_(interfaceData_->SearchParam()),
      searchuseauxpos_(interfaceData_->SearchUseAuxPos()),
      inttime_interface_(interfaceData_->IntTimeInterface()),
      nurbs_(interfaceData_->IsNurbs()),
      poro_(interfaceData_->IsPoro()),
      ehl_(interfaceData_->IsEhl())
{
  if (not interfaceData_->IsInit())
    dserror(
        "This constructor is only allowed for already initialized "
        "interface data containers!");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<MORTAR::MortarInterface> MORTAR::MortarInterface::Create(const int id,
    const Epetra_Comm& comm, const int spatialDim, const Teuchos::ParameterList& imortar,
    INPAR::MORTAR::RedundantStorage redundant)
{
  Teuchos::RCP<MORTAR::InterfaceDataContainer> interfaceData =
      Teuchos::rcp(new MORTAR::InterfaceDataContainer());

  return Teuchos::rcp(
      new MORTAR::MortarInterface(interfaceData, id, comm, spatialDim, imortar, redundant));
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarInterface::MortarInterface(const Teuchos::RCP<InterfaceDataContainer>& interfaceData,
    const int id, const Epetra_Comm& comm, const int spatialDim,
    const Teuchos::ParameterList& imortar, INPAR::MORTAR::RedundantStorage redundant)
    : interfaceData_(interfaceData),
      id_(interfaceData_->Id()),
      comm_(interfaceData_->CommPtr()),
      lcomm_(interfaceData_->lComm()),
      procmap_(interfaceData_->ProcMap()),
      redistributed_(interfaceData_->IsRedistributed()),
      idiscret_(interfaceData_->IDiscret()),
      dim_(interfaceData_->Dim()),
      imortar_(interfaceData_->IMortar()),
      shapefcn_(interfaceData_->ShapeFcn()),
      quadslave_(interfaceData_->IsQuadSlave()),
      redundant_(interfaceData_->RedundantStorage()),
      oldnodecolmap_(interfaceData_->OldNodeColMap()),
      oldelecolmap_(interfaceData_->OldEleColMap()),
      snoderowmap_(interfaceData_->SNodeRowMap()),
      snodecolmap_(interfaceData_->SNodeColMap()),
      mnoderowmap_(interfaceData_->MNodeRowMap()),
      mnodecolmap_(interfaceData_->MNodeColMap()),
      snoderowmapbound_(interfaceData_->SNodeRowMapBound()),
      snodecolmapbound_(interfaceData_->SNodeColMapBound()),
      mnoderowmapnobound_(interfaceData_->MNodeRowMapNoBound()),
      mnodecolmapnobound_(interfaceData_->MNodeColMapNoBound()),
      selerowmap_(interfaceData_->SEleRowMap()),
      selecolmap_(interfaceData_->SEleColMap()),
      melerowmap_(interfaceData_->MEleRowMap()),
      melecolmap_(interfaceData_->MEleColMap()),
      sdofrowmap_(interfaceData_->SDofRowMap()),
      sdofcolmap_(interfaceData_->SDofColMap()),
      mdofrowmap_(interfaceData_->MDofRowMap()),
      mdofcolmap_(interfaceData_->MDofColMap()),
      psdofrowmap_(interfaceData_->PSDofRowMap()),
      plmdofmap_(interfaceData_->PLmDofRowMap()),
      lmdofmap_(interfaceData_->LmDofRowMap()),
      maxdofglobal_(interfaceData_->MaxDofGlobal()),
      searchalgo_(interfaceData_->SearchAlgorithm()),
      binarytree_(interfaceData_->BinaryTree()),
      searchparam_(interfaceData_->SearchParam()),
      searchuseauxpos_(interfaceData_->SearchUseAuxPos()),
      inttime_interface_(interfaceData_->IntTimeInterface()),
      nurbs_(interfaceData_->IsNurbs()),
      poro_(interfaceData_->IsPoro()),
      ehl_(interfaceData_->IsEhl())
{
  interfaceData_->SetIsInit(true);
  id_ = id;
  comm_ = Teuchos::rcpFromRef(comm);
  dim_ = spatialDim;
  imortar_.setParameters(imortar);
  quadslave_ = false;
  redundant_ = redundant;
  searchalgo_ =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::SearchAlgorithm>(imortar, "SEARCH_ALGORITHM");
  searchparam_ = imortar.get<double>("SEARCH_PARAM");
  searchuseauxpos_ = DRT::INPUT::IntegralValue<int>(imortar, "SEARCH_USE_AUX_POS");
  nurbs_ = imortar.get<bool>("NURBS");

  if (Dim() != 2 && Dim() != 3) dserror("Mortar problem must be 2D or 3D.");

  procmap_.clear();

  CreateInterfaceDiscretization();
  SetShapeFunctionType();

  poro_ = false;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::CreateInterfaceDiscretization()
{
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());

  // Create name for mortar interface discretization
  std::stringstream dis_name;
  dis_name << "mortar_interface_" << id_;

  // Create the required type of discretization
  if (nurbs_)
  {
    idiscret_ = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization(dis_name.str(), comm));

    /*
    Note: The NurbsDiscretization needs a Knotvector to be able to write output. This is probably
    the place, where the Knotvector of the mortar coupling surface needs to be inserted into the
    NurbsDiscretization. However, it's not clear how to compute that Knotvector, so we can't do it
    right now. In the end, it should be sufficient to extract the portion from the underlying volume
    discretization's knot vector that corresponds to the mortar interface.

    As a NurbsDiscretization can't write output for now, we don't do it and rather use the 'old'
    output style, where interface output is written by the underlying volume discretization.
    */
  }
  else
  {
    idiscret_ = Teuchos::rcp(new DRT::Discretization(dis_name.str(), comm));
  }

  // Prepare discretization writer
  idiscret_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(idiscret_)));
  dsassert(not idiscret_->Writer().is_null(), "Setup of discretization writer failed.");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::SetShapeFunctionType()
{
  INPAR::MORTAR::ShapeFcn shapefcn =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(InterfaceParams(), "LM_SHAPEFCN");
  switch (shapefcn)
  {
    case INPAR::MORTAR::shape_dual:
    {
      shapefcn_ = INPAR::MORTAR::shape_dual;
      break;
    }
    case INPAR::MORTAR::shape_petrovgalerkin:
    {
      shapefcn_ = INPAR::MORTAR::shape_petrovgalerkin;
      break;
    }
    case INPAR::MORTAR::shape_standard:
    {
      shapefcn_ = INPAR::MORTAR::shape_standard;
      break;
    }
    default:
    {
      dserror(
          "Invalid shape function type. Interface must either have dual or standard shape "
          "functions.");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const MORTAR::MortarInterface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    os << "\nMortar Interface Id " << id_ << std::endl;
    os << "Mortar Interface Discretization:" << std::endl;
  }
  os << Discret();
  return;
}

/*----------------------------------------------------------------------*
 |  check if interface is FillComplete (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::Filled() const { return idiscret_->Filled(); }

/*----------------------------------------------------------------------*
 |  print parallel distribution (public)                      popp 06/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::PrintParallelDistribution() const
{
  // how many processors
  const int numproc = Discret().Comm().NumProc();

  // only print parallel distribution if numproc > 1
  if (numproc > 1)
  {
    const int myrank = Discret().Comm().MyPID();

    std::vector<int> my_n_nodes(numproc, 0);
    std::vector<int> n_nodes(numproc, 0);
    std::vector<int> my_n_ghostnodes(numproc, 0);
    std::vector<int> n_ghostnodes(numproc, 0);
    std::vector<int> my_n_elements(numproc, 0);
    std::vector<int> n_elements(numproc, 0);
    std::vector<int> my_n_ghostele(numproc, 0);
    std::vector<int> n_ghostele(numproc, 0);
    std::vector<int> my_s_nodes(numproc, 0);
    std::vector<int> s_nodes(numproc, 0);
    std::vector<int> my_s_ghostnodes(numproc, 0);
    std::vector<int> s_ghostnodes(numproc, 0);
    std::vector<int> my_s_elements(numproc, 0);
    std::vector<int> s_elements(numproc, 0);
    std::vector<int> my_s_ghostele(numproc, 0);
    std::vector<int> s_ghostele(numproc, 0);
    std::vector<int> my_m_nodes(numproc, 0);
    std::vector<int> m_nodes(numproc, 0);
    std::vector<int> my_m_elements(numproc, 0);
    std::vector<int> m_elements(numproc, 0);
    std::vector<int> my_m_ghostnodes(numproc, 0);
    std::vector<int> m_ghostnodes(numproc, 0);
    std::vector<int> my_m_ghostele(numproc, 0);
    std::vector<int> m_ghostele(numproc, 0);

    my_n_nodes[myrank] = Discret().NumMyRowNodes();
    my_n_ghostnodes[myrank] = Discret().NumMyColNodes() - my_n_nodes[myrank];
    my_n_elements[myrank] = Discret().NumMyRowElements();
    my_n_ghostele[myrank] = Discret().NumMyColElements() - my_n_elements[myrank];

    my_s_nodes[myrank] = snoderowmap_->NumMyElements();
    my_s_ghostnodes[myrank] = snodecolmap_->NumMyElements() - my_s_nodes[myrank];
    my_s_elements[myrank] = selerowmap_->NumMyElements();
    my_s_ghostele[myrank] = selecolmap_->NumMyElements() - my_s_elements[myrank];

    my_m_nodes[myrank] = mnoderowmap_->NumMyElements();
    my_m_ghostnodes[myrank] = mnodecolmap_->NumMyElements() - my_m_nodes[myrank];
    my_m_elements[myrank] = melerowmap_->NumMyElements();
    my_m_ghostele[myrank] = melecolmap_->NumMyElements() - my_m_elements[myrank];

    // adapt output for redundant master or all redundant case
    if (Redundant() == INPAR::MORTAR::redundant_master)
    {
      my_m_ghostnodes[myrank] = mnoderowmap_->NumGlobalElements() - my_m_nodes[myrank];
      my_m_ghostele[myrank] = melerowmap_->NumGlobalElements() - my_m_elements[myrank];
    }
    else if (Redundant() == INPAR::MORTAR::redundant_all)
    {
      my_m_ghostnodes[myrank] = mnoderowmap_->NumGlobalElements() - my_m_nodes[myrank];
      my_m_ghostele[myrank] = melerowmap_->NumGlobalElements() - my_m_elements[myrank];
      my_s_ghostnodes[myrank] = snoderowmap_->NumGlobalElements() - my_s_nodes[myrank];
      my_s_ghostele[myrank] = selerowmap_->NumGlobalElements() - my_s_elements[myrank];
    }

    Discret().Comm().SumAll(&my_n_nodes[0], &n_nodes[0], numproc);
    Discret().Comm().SumAll(&my_n_ghostnodes[0], &n_ghostnodes[0], numproc);
    Discret().Comm().SumAll(&my_n_elements[0], &n_elements[0], numproc);
    Discret().Comm().SumAll(&my_n_ghostele[0], &n_ghostele[0], numproc);

    Discret().Comm().SumAll(&my_s_nodes[0], &s_nodes[0], numproc);
    Discret().Comm().SumAll(&my_s_ghostnodes[0], &s_ghostnodes[0], numproc);
    Discret().Comm().SumAll(&my_s_elements[0], &s_elements[0], numproc);
    Discret().Comm().SumAll(&my_s_ghostele[0], &s_ghostele[0], numproc);

    Discret().Comm().SumAll(&my_m_nodes[0], &m_nodes[0], numproc);
    Discret().Comm().SumAll(&my_m_ghostnodes[0], &m_ghostnodes[0], numproc);
    Discret().Comm().SumAll(&my_m_elements[0], &m_elements[0], numproc);
    Discret().Comm().SumAll(&my_m_ghostele[0], &m_ghostele[0], numproc);

    if (myrank == 0)
    {
      std::cout << std::endl;
      std::cout << "   Discretization: " << Discret().Name() << std::endl;
      printf("   +-----+-----------------+--------------+-----------------+--------------+\n");
      printf("   | PID |   n_rownodes    | n_ghostnodes |  n_rowelements  |  n_ghostele  |\n");
      printf("   +-----+-----------------+--------------+-----------------+--------------+\n");
      for (int npid = 0; npid < numproc; ++npid)
      {
        printf("   | %3d | Total %9d | %12d | Total %9d | %12d |\n", npid, n_nodes[npid],
            n_ghostnodes[npid], n_elements[npid], n_ghostele[npid]);
        printf("   |     | Slave %9d | %12d | Slave %9d | %12d |\n", s_nodes[npid],
            s_ghostnodes[npid], s_elements[npid], s_ghostele[npid]);
        printf("   |     | Master %8d | %12d | Master %8d | %12d |\n", m_nodes[npid],
            m_ghostnodes[npid], m_elements[npid], m_ghostele[npid]);
        printf("   +-----+-----------------+--------------+-----------------+--------------+\n");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  add mortar node (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AddMortarNode(Teuchos::RCP<MORTAR::MortarNode> mrtrnode)
{
  idiscret_->AddNode(mrtrnode);
  return;
}

/*----------------------------------------------------------------------*
 |  add mortar element (public)                              mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AddMortarElement(Teuchos::RCP<MORTAR::MortarElement> mrtrele)
{
  // check for quadratic 2d slave elements to be modified
  if (mrtrele->IsSlave() &&
      (mrtrele->Shape() == DRT::Element::line3 || mrtrele->Shape() == DRT::Element::nurbs3))
    quadslave_ = true;

  // check for quadratic 3d slave elements to be modified
  if (mrtrele->IsSlave() &&
      (mrtrele->Shape() == DRT::Element::quad9 || mrtrele->Shape() == DRT::Element::quad8 ||
          mrtrele->Shape() == DRT::Element::tri6 || mrtrele->Shape() == DRT::Element::nurbs8 ||
          mrtrele->Shape() == DRT::Element::nurbs9))
    quadslave_ = true;

  idiscret_->AddElement(mrtrele);
  return;
}

/*----------------------------------------------------------------------*
 |  delete all elements/nodes on one side of the interface   ghamm 11/13|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::RemoveSingleInterfaceSide(bool slaveside)
{
  Teuchos::RCP<Epetra_Map> elecolmap;
  Teuchos::RCP<Epetra_Map> nodecolmap;
  Teuchos::RCP<Epetra_Map> pointcolmap;

  if (idiscret_->Filled())
  {
    // extract maps before deleting first element/node/point
    if (slaveside == true)
    {
      elecolmap = SlaveColElements();
      nodecolmap = SlaveColNodes();
    }
    else
    {
      elecolmap = MasterColElements();
      nodecolmap = MasterColNodes();
    }

    // delete elements on desired side
    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      idiscret_->DeleteElement(gid);
    }

    // delete nodes on desired side
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      idiscret_->DeleteNode(gid);
    }
  }
  else
  {
    idiscret_->FillComplete(false, false, false);
    // extract maps before deleting first element/node/point
    elecolmap = Teuchos::rcp(new Epetra_Map(*idiscret_->ElementColMap()));
    nodecolmap = Teuchos::rcp(new Epetra_Map(*idiscret_->NodeColMap()));

    // delete elements on desired side
    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      bool isslave = dynamic_cast<MORTAR::MortarElement*>(idiscret_->gElement(gid))->IsSlave();
      if (isslave == slaveside) idiscret_->DeleteElement(gid);
    }

    // delete nodes on desired side
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      bool isslave = dynamic_cast<MORTAR::MortarNode*>(idiscret_->gNode(gid))->IsSlave();
      if (isslave == slaveside) idiscret_->DeleteNode(gid);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  finalize construction of interface (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::FillComplete(int maxdof, bool newghosting)
{
  // store maximum global dof ID handed in
  // this ID is later needed when setting up the Lagrange multiplier
  // dof map, which of course must not overlap with existing dof ranges
  maxdofglobal_ = maxdof;

  // we'd like to call idiscret_.FillComplete(true,false,false) but this
  // will assign all nodes new degrees of freedom which we don't want.
  // We would like to use the degrees of freedom that were stored in the
  // mortar nodes. To do so, we have to create and set our own
  // version of a DofSet class before we call FillComplete on the
  // interface discretization.
  // Our special dofset class will not assign new dofs but will assign the
  // dofs stored in the nodes.
  {
    Teuchos::RCP<MORTAR::MortarDofSet> mrtrdofset = Teuchos::rcp(new MORTAR::MortarDofSet());
    Discret().ReplaceDofSet(mrtrdofset);
    // do not assign dofs yet, we'll do this below after
    // shuffling around of nodes and elements (saves time)
    Discret().FillComplete(false, false, false);
  }

  // check whether crosspoints / edge nodes shall be considered or not
  InitializeCrossPoints();

  // check for const/linear interpolation of 2D/3D quadratic Lagrange multipliers
  InitializeLagMultConst();
  InitializeLagMultLin();

  // check/init corner/edge modification
  InitializeCornerEdge();

  // later we might export node and element column map to extended or even FULL overlap,
  // thus store the standard column maps first
  // get standard nodal column map (overlap=1)
  oldnodecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().NodeColMap())));
  // get standard element column map (overlap=1)
  oldelecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().ElementColMap())));

  CreateInterfaceLocalCommunicator();

  // create interface ghosting
  // (currently, the slave is kept with the standard overlap of one,
  // but the master is made fully redundant, i.e. it is exported to
  // fully overlapping column layout, for the ease of interface search)
  // (the only exceptions are self contact and coupled problems, where
  // also the slave is still made fully redundant)
  // ghosting can be skipped if the desired parallel layout is already present
  if (newghosting) ExtendInterfaceGhosting();

  // make sure discretization is complete
  Discret().FillComplete(true, false, false);

  // ghost also parent elements according to the ghosting strategy of the interface (atm just for
  // poro)
  if (newghosting && poro_)
    POROELAST::UTILS::CreateVolumeGhosting(Discret());
  else if (imortar_.isParameter("STRATEGY"))
    if (newghosting && DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(
                           imortar_, "ALGORITHM") == INPAR::MORTAR::algorithm_gpts)
      CreateVolumeGhosting();

  // need row and column maps of slave and master nodes / elements / dofs
  // separately so we can easily address them
  UpdateMasterSlaveSets();

  // initialize node and element data container
  InitializeDataContainer();

  // communicate quadslave status among ALL processors
  // (not only those participating in interface)
  int localstatus = (int)(quadslave_);
  int globalstatus = 0;
  Comm().SumAll(&localstatus, &globalstatus, 1);
  quadslave_ = (bool)(globalstatus);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::CreateInterfaceLocalCommunicator()
{
  std::vector<int> lin(Comm().NumProc());
  std::vector<int> gin(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) lin[i] = 0;

  // check ownership or ghosting of any elements / nodes
  // const Epetra_Map* nodemap = Discret().NodeColMap();
  // const Epetra_Map* elemap  = Discret().ElementColMap();

  //********************************************************************
  // NOTE: currently we choose local=global communicator, but we have
  // all structures present in the code to change this assignment any time.
  //********************************************************************
  // if (nodemap->NumMyElements() || elemap->NumMyElements())
  lin[Comm().MyPID()] = 1;

  Comm().MaxAll(&lin[0], &gin[0], Comm().NumProc());
  lin.clear();

  // build global -> local communicator PID map
  // we need this when calling Broadcast() on lComm later
  int counter = 0;
  for (int i = 0; i < Comm().NumProc(); ++i)
  {
    if (gin[i])
      procmap_[i] = counter++;
    else
      procmap_[i] = -1;
  }

  // typecast the Epetra_Comm to Epetra_MpiComm
  Teuchos::RCP<Epetra_Comm> copycomm = Teuchos::rcp(Comm().Clone());
  Epetra_MpiComm* epetrampicomm = dynamic_cast<Epetra_MpiComm*>(copycomm.get());
  if (epetrampicomm != NULL)
  {
    // split the communicator into participating and none-participating procs
    int color;
    int key = Comm().MyPID();
    // I am taking part in the new comm if I have any ownership
    if (gin[Comm().MyPID()]) color = 0;
    // I am not taking part in the new comm
    else
      color = MPI_UNDEFINED;

    // tidy up
    gin.clear();

    // free lcomm_ first
    if (lcomm_ != Teuchos::null)
    {
      MPI_Comm oldcomm = Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(lcomm_)->GetMpiComm();
      lcomm_ = Teuchos::null;
      MPI_Comm_free(&oldcomm);
    }

    // create new local communicator
    MPI_Comm mpi_global_comm = epetrampicomm->GetMpiComm();
    MPI_Comm mpi_local_comm;
    MPI_Comm_split(mpi_global_comm, color, key, &mpi_local_comm);

    // create the new Epetra_MpiComm only for participating procs
    if (mpi_local_comm != MPI_COMM_NULL) lcomm_ = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
  }
  else
  {
    // check for serial communicator
    if (Comm().NumProc() != 1) dserror("ERROR: Epetra_SerialComm can only handle 1 processor!");
    Epetra_SerialComm* serialcomm = dynamic_cast<Epetra_SerialComm*>(copycomm.get());
    if (!serialcomm) dserror("ERROR: casting Epetra_Comm -> Epetra_SerialComm failed");
    lcomm_ = Teuchos::rcp(new Epetra_SerialComm(*serialcomm));
  }
}

/*----------------------------------------------------------------------*
 |  Check and initialize corner/edge contact                 farah 07/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::InitializeCornerEdge()
{
  // if linear LM for quad displacements return!
  // TODO: this case needs a special treatment
  bool lagmultlin = (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(
                         InterfaceParams(), "LM_QUAD") == INPAR::MORTAR::lagmult_lin);

  if (lagmultlin) return;

  for (int i = 0; i < (Discret().NodeRowMap())->NumMyElements(); ++i)
  {
    // static_cast to the corresponding mortar/contact/friction/... node
    // or element would be enough in all places
    // performance loss is negligible when using a dynamic_cast instead
    // but safety is increased enormously
    MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

    // remove bound/corner/edge flag for master nodes!
    if (!node->IsSlave() && node->IsOnCorner() == true) node->SetOnCorner() = false;
    if (!node->IsSlave() && node->IsOnEdge() == true) node->SetOnEdge() = false;

    // candidates are slave nodes with only 1 adjacent MortarElement
    if (node->IsSlave() && node->IsOnCornerEdge() == true)
    {
      node->SetSlave() = false;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Check and initialize cross points                        farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::InitializeCrossPoints()
{
  // check whether crosspoints / edge nodes shall be considered or not
  bool crosspoints = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "CROSSPOINTS");

  // modify crosspoints / edge nodes
  if (crosspoints)
  {
    // only applicable for 2D problems up to now
    if (Dim() == 3) dserror("ERROR: Crosspoint / edge node modification not yet impl. for 3D");

    // ---------------------------------------------------------------------
    // Detect relevant nodes on slave side
    // ---------------------------------------------------------------------
    // A typical application are so-called crosspoints within mortar mesh
    // tying, where this approach is necessary to avoid over-constraint.
    // Otherwise these crosspoints would be active with respect to more
    // than one interface and thus the LM cannot sufficiently represent
    // all geometrical constraints. Another typical application is mortar
    // contact, when we want to make use of symmetry boundary conditions.
    // In this case, we deliberately modify so-called edge nodes of the
    // contact boundary and thus free them from any contact constraint.
    // ---------------------------------------------------------------------
    // Basically, the status of the crosspoints / edge nodes is simply
    // changed to MASTER and consequently they will NOT carry Lagrange
    // multipliers later on. In order to sustain the partition of unity
    // property of the LM shape functions on the adjacent slave elements,
    // the LM shape functions of the adjacent nodes will be modified! This
    // way, the mortar operator entries of the crosspoints / edge nodes are
    // transferred to the neighboring slave nodes!
    // ---------------------------------------------------------------------

    for (int i = 0; i < (Discret().NodeRowMap())->NumMyElements(); ++i)
    {
      // static_cast to the corresponding mortar/contact/friction/... node
      // or element would be enough in all places
      // performance loss is negligible when using a dynamic_cast instead
      // but safety is increased enormously
      MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with only 1 adjacent MortarElement
      if (node->IsSlave() && node->NumElement() == 1)
      {
        // case1: linear shape functions, boundary nodes already found
        if ((node->Elements()[0])->NumNode() == 2)
        {
          node->SetBound() = true;
          node->SetSlave() = false;
        }
        // case2: quad. shape functions, middle nodes must be sorted out
        else if (node->Id() != (node->Elements()[0])->NodeIds()[2])
        {
          node->SetBound() = true;
          node->SetSlave() = false;
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Check and initialize for lin lagmult interpolation       farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::InitializeLagMultLin()
{
  // check for linear interpolation of 2D/3D quadratic Lagrange multipliers
  bool lagmultlin = (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(
                         InterfaceParams(), "LM_QUAD") == INPAR::MORTAR::lagmult_lin);

  // modify nodes accordingly
  if (lagmultlin)
  {
    // modified treatment of vertex nodes and edge nodes
    // detect middle nodes (quadratic nodes) on slave side
    // set status of middle nodes -> MASTER
    // set status of vertex nodes -> SLAVE

    // loop over all elements
    for (int i = 0; i < Discret().NodeRowMap()->NumMyElements(); ++i)
    {
      // get node and cast to cnode
      MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with shape line3 (2D), tri6 and quad8/9 (3D)
      if (node->IsSlave())
      {
        // search the first adjacent element
        MORTAR::MortarElement::DiscretizationType shape = (node->Elements()[0])->Shape();

        // which discretization type
        switch (shape)
        {
          // line3 contact elements (= quad8/9 or tri6 discretizations)
          case MORTAR::MortarElement::line3:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0] ||
                node->Id() == (node->Elements()[0])->NodeIds()[1])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

            // tri6 contact elements (= tet10 discretizations)
          case MORTAR::MortarElement::tri6:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0] ||
                node->Id() == (node->Elements()[0])->NodeIds()[1] ||
                node->Id() == (node->Elements()[0])->NodeIds()[2])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

            // quad8 contact elements (= hex20 discretizations)
          case MORTAR::MortarElement::quad8:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0] ||
                node->Id() == (node->Elements()[0])->NodeIds()[1] ||
                node->Id() == (node->Elements()[0])->NodeIds()[2] ||
                node->Id() == (node->Elements()[0])->NodeIds()[3])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

            // quad9 contact elements (= hex27 discretizations)
          case MORTAR::MortarElement::quad9:
          {
            // case1: vertex nodes remain SLAVE
            if (node->Id() == (node->Elements()[0])->NodeIds()[0] ||
                node->Id() == (node->Elements()[0])->NodeIds()[1] ||
                node->Id() == (node->Elements()[0])->NodeIds()[2] ||
                node->Id() == (node->Elements()[0])->NodeIds()[3])
            {
              // do nothing
            }

            // case2: middle nodes must be set to MASTER
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            break;
          }

            // other cases
          default:
          {
            dserror(
                "ERROR: Lin/Lin interpolation of LM only for line3/tri6/quad8/quad9 mortar "
                "elements");
            break;
          }
        }  // switch(Shape)
      }    // if (IsSlave())
    }      // for-loop
  }

  return;
}



/*----------------------------------------------------------------------*
 |  Check and initialize for const lagmult interpolation     seitz 09/17|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::InitializeLagMultConst()
{
  if ((DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(InterfaceParams(), "LM_QUAD") ==
          INPAR::MORTAR::lagmult_const))
  {
    // modified treatment slave side nodes:
    // only the center-node carries LM

    // loop over all elements
    for (int i = 0; i < Discret().NodeRowMap()->NumMyElements(); ++i)
    {
      // get node and cast to cnode
      MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with shape line3 (2D), tri6 and quad8/9 (3D)
      if (node->IsSlave())
      {
        // search the first adjacent element
        MORTAR::MortarElement::DiscretizationType shape = (node->Elements()[0])->Shape();

        // which discretization type
        switch (shape)
        {
          // line3 contact elements (= quad8/9 or tri6 discretizations)
          case MORTAR::MortarElement::line3:
          {
            // case1: vertex nodes must be set to MASTER
            if (node->Id() == (node->Elements()[0])->NodeIds()[0] ||
                node->Id() == (node->Elements()[0])->NodeIds()[1])
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }

            // case2: middle nodes remain SLAVE
            else
            {
              // do nothing
            }

            break;
          }
          case MORTAR::MortarElement::quad9:
            if (node->Id() == (node->Elements()[0])->NodeIds()[8])
            {
              // do nothing
            }
            else
            {
              node->SetBound() = true;
              node->SetSlave() = false;
            }
            break;

            // other cases
          default:
          {
            dserror(
                "ERROR: Lin/Lin interpolation of LM only for line3/tri6/quad8/quad9 mortar "
                "elements");
            break;
          }
        }  // switch(Shape)
      }    // if (IsSlave())
    }      // for-loop
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Initialize Data Container for nodes and elements         farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::InitializeDataContainer()
{
  // initialize node data container
  // (include slave side boundary nodes / crosspoints)
  for (int i = 0; i < SlaveColNodesBound()->NumMyElements(); ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %i", gid);
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    //********************************************************
    // NOTE: depending on which kind of node this really is,
    // i.e. mortar, contact or friction node, several derived
    // versions of the InitializeDataContainer() methods will
    // be called here, apart from the base class version.
    //********************************************************

    // initialize container if not yet initialized before
    mnode->InitializeDataContainer();
    if (poro_)  // initialize just for poro contact case!
      mnode->InitializePoroDataContainer();
    if (ehl_) mnode->InitializeEhlDataContainer();
  }
  if (poro_)  // as velocities of structure and fluid exist also on master nodes!!!
  {
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(MasterRowNodes()));
    // initialize poro node data container for master nodes!!!
    for (int i = 0; i < masternodes()->NumMyElements(); ++i)
    {
      int gid = masternodes()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %i", gid);
      MortarNode* mnode = dynamic_cast<MortarNode*>(node);

      // ATM just implemented for ContactNode ... otherwise error!!!

      // initialize container if not yet initialized before
      mnode->InitializePoroDataContainer();
    }
  }

  // initialize element data container
  for (int i = 0; i < SlaveColElements()->NumMyElements(); ++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MortarElement* mele = dynamic_cast<MortarElement*>(ele);

    // initialize container if not yet initialized before
    mele->InitializeDataContainer();
  }

  if (poro_)
  {
    // initialize master element data container
    for (int i = 0; i < MasterColElements()->NumMyElements(); ++i)
    {
      int gid = MasterColElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
      MortarElement* mele = dynamic_cast<MortarElement*>(ele);

      // initialize container if not yet initialized before
      mele->InitializeDataContainer();
    }
  }

  if (InterfaceParams().isParameter("ALGORITHM"))
    if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(InterfaceParams(), "ALGORITHM") ==
        INPAR::MORTAR::algorithm_gpts)
      for (int i = 0; i < MasterColElements()->NumMyElements(); ++i)
        dynamic_cast<MortarElement*>(Discret().gElement(MasterColElements()->GID(i)))
            ->InitializeDataContainer();

  return;
}


/*----------------------------------------------------------------------*
 |  Parallel Strategy based on bin distribution (public)     farah 11/13|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::BinningStrategy(
    Teuchos::RCP<Epetra_Map> initial_elecolmap, const double vel)
{
  // Create the binning strategy
  Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy = SetupBinningStrategy(vel);

  // fill master and slave elements into bins
  std::map<int, std::set<int>> slavebinelemap;
  binningstrategy->DistributeElesToBins(Discret(), slavebinelemap, true);
  std::map<int, std::set<int>> masterbinelemap;
  binningstrategy->DistributeElesToBins(Discret(), masterbinelemap, false);

  // ghosting is extended
  Teuchos::RCP<Epetra_Map> extendedmastercolmap = binningstrategy->ExtendGhosting(
      Discret(), *initial_elecolmap, slavebinelemap, masterbinelemap);

  // adapt layout to extended ghosting in discret
  // first export the elements according to the processor local element column maps
  Discret().ExportColumnElements(*extendedmastercolmap);

  // get the node ids of the elements that are to be ghosted and create a proper node column map for
  // their export
  std::set<int> nodes;
  for (int lid = 0; lid < extendedmastercolmap->NumMyElements(); ++lid)
  {
    DRT::Element* ele = Discret().gElement(extendedmastercolmap->GID(lid));
    const int* nodeids = ele->NodeIds();
    for (int inode = 0; inode < ele->NumNode(); ++inode) nodes.insert(nodeids[inode]);
  }

  std::vector<int> colnodes(nodes.begin(), nodes.end());
  Teuchos::RCP<Epetra_Map> nodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), &colnodes[0], 0, Comm()));

  // now ghost the nodes
  Discret().ExportColumnNodes(*nodecolmap);

  // fillcomplete interface
  FillComplete();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BINSTRATEGY::BinningStrategy> MORTAR::MortarInterface::SetupBinningStrategy(
    const double vel)
{
  // Initialize eXtendedAxisAlignedBoundingBox (XAABB)
  LINALG::Matrix<3, 2> XAABB(false);
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = +1.0e12;
    XAABB(dim, 1) = -1.0e12;
  }

  // loop all slave nodes and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int lid = 0; lid < SlaveColNodes()->NumMyElements(); ++lid)
  {
    int gid = SlaveColNodes()->GID(lid);
    DRT::Node* node = Discret().gNode(gid);
    if (!node)
      dserror(
          "Cannot find node with gid %i in discretization '%s'.", gid, Discret().Name().c_str());
    MORTAR::MortarNode* mtrnode = dynamic_cast<MORTAR::MortarNode*>(node);

    for (unsigned int dim = 0; dim < 3; ++dim)
    {
      XAABB(dim, 0) = std::min(XAABB(dim, 0), mtrnode->xspatial()[dim] - MORTARPROJLIM);
      XAABB(dim, 1) = std::max(XAABB(dim, 1), mtrnode->xspatial()[dim] + MORTARPROJLIM);
    }
  }

  // local bounding box
  double locmin[3] = {XAABB(0, 0), XAABB(1, 0), XAABB(2, 0)};
  double locmax[3] = {XAABB(0, 1), XAABB(1, 1), XAABB(2, 1)};

  // global bounding box
  double globmin[3];
  double globmax[3];

  // do the necessary communication
  Comm().MinAll(&locmin[0], &globmin[0], 3);
  Comm().MaxAll(&locmax[0], &globmax[0], 3);

  // compute cutoff radius:
  double totalsme = -1.0;
  for (int lid = 0; lid < SlaveColElements()->NumMyElements(); ++lid)
  {
    int gid = SlaveColElements()->GID(lid);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele)
      dserror(
          "Cannot find element with gid %i in discretization '%s'.", gid, Discret().Name().c_str());
    MORTAR::MortarElement* mtrele = dynamic_cast<MORTAR::MortarElement*>(ele);

    // to be thought about, whether this is enough (safety = 2??)
    double sme = mtrele->MaxEdgeSize();
    totalsme = std::max(sme, totalsme);
  }

  double cutoff = -1.0;
  Comm().MaxAll(&totalsme, &cutoff, 1);

  // extend cutoff based on problem interface velocity
  // --> only for contact problems
  if (vel >= EPS12)
  {
    const double dt = InterfaceParams().get<double>("TIMESTEP");
    cutoff = cutoff + 2 * dt * vel;
  }

  // increase XAABB by 2x cutoff radius
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = globmin[dim] - cutoff;
    XAABB(dim, 1) = globmax[dim] + cutoff;
  }

  return Teuchos::rcp(new BINSTRATEGY::BinningStrategy(Comm(), cutoff, XAABB));
}


/*----------------------------------------------------------------------*
 |  redistribute interface (public)                           popp 08/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Redistribute()
{
  std::stringstream ss;
  ss << "MORTAR::MortarInterface::Redistribute of '" << Discret().Name() << "'";
  TEUCHOS_FUNC_TIME_MONITOR(ss.str());

  const Teuchos::ParameterList& mortarParallelRedistParams =
      InterfaceParams().sublist("PARALLEL REDISTRIBUTION");

  // make sure we are supposed to be here
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
          mortarParallelRedistParams, "PARALLEL_REDIST") == INPAR::MORTAR::parredist_none)
    dserror("ERROR: You are not supposed to be here...");

  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();
  Epetra_Time time(*comm);

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  // we need an arbitrary preliminary element row map
  Teuchos::RCP<Epetra_Map> sroweles = Teuchos::rcp(new Epetra_Map(*SlaveRowElements()));
  Teuchos::RCP<Epetra_Map> mroweles = Teuchos::rcp(new Epetra_Map(*MasterRowElements()));

  //**********************************************************************
  // (1) PREPARATIONS decide how many procs are used
  //**********************************************************************
  // first we assume that all procs will be used
  int sproc = numproc;
  int mproc = numproc;

  // minimum number of elements per proc
  int minele = mortarParallelRedistParams.get<int>("MIN_ELEPROC");

  // Max. relative imbalance between subdomain sizes
  const double imbalance_tol = mortarParallelRedistParams.get<double>("IMBALANCE_TOL");

  // calculate real number of procs to be used
  if (minele > 0)
  {
    sproc = static_cast<int>((sroweles->NumGlobalElements()) / minele);
    mproc = static_cast<int>((mroweles->NumGlobalElements()) / minele);
    if (sroweles->NumGlobalElements() < 2 * minele) sproc = 1;
    if (mroweles->NumGlobalElements() < 2 * minele) mproc = 1;
    if (sproc > numproc) sproc = numproc;
    if (mproc > numproc) mproc = numproc;
  }

  // print message
  if (!myrank)
  {
    std::cout << "\nRedistributing interface '" << Discret().Name() << "' using ZOLTAN.........\n";
    std::cout << "Procs used for redistribution: " << sproc << " / " << mproc << " (S / M)\n";
  }

  //**********************************************************************
  // (2) SLAVE redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> srownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> scolnodes = Teuchos::null;

  DRT::UTILS::REBALANCING::ComputeRebalancedNodeMaps(
      idiscret_, sroweles, srownodes, scolnodes, comm, false, sproc, imbalance_tol);

  //**********************************************************************
  // (3) MASTER redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> mrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> mcolnodes = Teuchos::null;

  RedistributeMasterSide(mrownodes, mcolnodes, mroweles, comm, mproc, imbalance_tol);

  //**********************************************************************
  // (4) Merge global interface node row and column map
  //**********************************************************************
  // merge node maps from slave and master parts
  Teuchos::RCP<Epetra_Map> rownodes = LINALG::MergeMap(srownodes, mrownodes, false);
  Teuchos::RCP<Epetra_Map> colnodes = LINALG::MergeMap(scolnodes, mcolnodes, false);

  //**********************************************************************
  // (5) Get partitioning information into discretization
  //**********************************************************************
  // build reasonable element maps from the already valid and final node maps
  // (note that nothing is actually redistributed in here)
  Teuchos::RCP<Epetra_Map> roweles = Teuchos::null;
  Teuchos::RCP<Epetra_Map> coleles = Teuchos::null;
  Discret().BuildElementRowColumn(*rownodes, *colnodes, roweles, coleles);

  // export nodes and elements to the row map
  Discret().ExportRowNodes(*rownodes);
  Discret().ExportRowElements(*roweles);

  // export nodes and elements to the column map (create ghosting)
  Discret().ExportColumnNodes(*colnodes);
  Discret().ExportColumnElements(*coleles);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MortarInterface::RedistributeMasterSide(Teuchos::RCP<Epetra_Map>& rownodes,
    Teuchos::RCP<Epetra_Map>& colnodes, const Teuchos::RCP<Epetra_Map>& roweles,
    const Teuchos::RCP<Epetra_Comm>& comm, const int parts, const double imbalance) const
{
  if (not HasMaSharingRefInterface())
  {
    // call ZOLTAN for parallel redistribution
    DRT::UTILS::REBALANCING::ComputeRebalancedNodeMaps(
        idiscret_, roweles, rownodes, colnodes, comm, false, parts, imbalance);
  }
  else
  {
    DRT::UTILS::RedistributeInAccordanceWithReference(
        *GetMaSharingRefInterfacePtr()->MasterRowNodes(), *MasterRowNodes(), rownodes);

    DRT::UTILS::RedistributeInAccordanceWithReference(
        *GetMaSharingRefInterfacePtr()->MasterColNodes(), *MasterColNodes(), colnodes);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::ExtendInterfaceGhosting()
{
  //*****REDUNDANT SLAVE AND MASTER STORAGE*****
  if (Redundant() == INPAR::MORTAR::redundant_all)
  {
    // to ease our search algorithms we'll afford the luxury to ghost all nodes
    // on all processors. To do so, we'll take the node row map and export it to
    // full overlap. Then we export the discretization to full overlap column map.
    // This way, also all mortar elements will be fully ghosted on all processors.
    // Note that we'll do ghosting NOT ONLY on procs that do own or ghost any of the
    // nodes in the natural distribution of idiscret_, but really on ALL procs.
    // This makes dynamic redistribution easier!

    //**********************************************************************
    // IMPORTANT NOTE:
    // In an older code version, we only did ghosting on procs that own or ghost
    // any of the interface nodes or elements  in the natural distr. of idiscret_.
    // The corresponding code lines for creating this proc list are:
    //
    // std::vector<int> stproc(0);
    // if (oldnodecolmap_->NumMyElements() || oldelecolmap_->NumMyElements())
    //   stproc.push_back(Comm().MyPID());
    // std::vector<int> rtproc(0);
    // LINALG::Gather<int>(stproc,rtproc,Comm().NumProc(),&allproc[0],Comm());
    //
    // In this case, we use "rtproc" instead of "allproc" afterwards, i.e. when
    // the node gids and element gids are gathered among procs.
    //**********************************************************************

    // we want to do full ghosting on all procs
    std::vector<int> allproc(Comm().NumProc());
    for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

    // fill my own row node ids
    const Epetra_Map* noderowmap = Discret().NodeRowMap();
    std::vector<int> sdata(noderowmap->NumMyElements());
    for (int i = 0; i < noderowmap->NumMyElements(); ++i) sdata[i] = noderowmap->GID(i);

    // gather all gids of nodes redundantly
    std::vector<int> rdata;
    LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), &allproc[0], Comm());

    // build completely overlapping map of nodes (on ALL processors)
    Teuchos::RCP<Epetra_Map> newnodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, Comm()));
    sdata.clear();
    rdata.clear();

    // fill my own row element ids
    const Epetra_Map* elerowmap = Discret().ElementRowMap();
    sdata.resize(elerowmap->NumMyElements());
    for (int i = 0; i < elerowmap->NumMyElements(); ++i) sdata[i] = elerowmap->GID(i);

    // gather all gids of elements redundantly
    rdata.resize(0);
    LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), &allproc[0], Comm());

    // build complete overlapping map of elements (on ALL processors)
    Teuchos::RCP<Epetra_Map> newelecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, Comm()));
    sdata.clear();
    rdata.clear();
    allproc.clear();

    // redistribute the discretization of the interface according to the
    // new column layout
    Discret().ExportColumnNodes(*newnodecolmap);
    Discret().ExportColumnElements(*newelecolmap);
  }

  //*****ONLY REDUNDANT MASTER STORAGE*****
  else if (Redundant() == INPAR::MORTAR::redundant_master)
  {
    // to ease our search algorithms we'll afford the luxury to ghost all master
    // nodes on all processors. To do so, we'll take the master node row map and
    // export it to full overlap. Then we export the discretization to partially
    // full overlap column map. This way, also all master elements will be fully
    // ghosted on all processors. Note that we'll do ghosting NOT ONLY on procs
    // that do own or ghost any nodes in the natural distribution of idiscret_,
    // but really on ALL procs. This makes dynamic redistribution easier!

    //**********************************************************************
    // IMPORTANT NOTE:
    // In an older code version, we only did ghosting on procs that own or ghost
    // any of the interface nodes or elements  in the natural distr. of idiscret_.
    // The corresponding code lines for creating this proc list are:
    //
    // std::vector<int> stproc(0);
    // if (oldnodecolmap_->NumMyElements() || oldelecolmap_->NumMyElements())
    //   stproc.push_back(Comm().MyPID());
    // std::vector<int> rtproc(0);
    // LINALG::Gather<int>(stproc,rtproc,Comm().NumProc(),&allproc[0],Comm());
    //
    // In this case, we use "rtproc" instead of "allproc" afterwards, i.e. when
    // the node gids and element gids are gathered among procs.
    //**********************************************************************

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
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);
      if (!mrtrnode->IsSlave()) sdata.push_back(gid);
    }

    // gather all master row node gids redundantly
    std::vector<int> rdata;
    LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), &allproc[0], Comm());

    // add my own slave column node ids (non-redundant, standard overlap)
    const Epetra_Map* nodecolmap = Discret().NodeColMap();
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);
      if (mrtrnode->IsSlave()) rdata.push_back(gid);
    }

    // build new node column map (on ALL processors)
    Teuchos::RCP<Epetra_Map> newnodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, Comm()));
    sdata.clear();
    rdata.clear();

    // fill my own master row element ids
    const Epetra_Map* elerowmap = Discret().ElementRowMap();
    sdata.resize(0);
    for (int i = 0; i < elerowmap->NumMyElements(); ++i)
    {
      int gid = elerowmap->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find element with gid %", gid);
      MortarElement* mrtrele = dynamic_cast<MortarElement*>(ele);
      if (!mrtrele->IsSlave()) sdata.push_back(gid);
    }

    // gather all gids of elements redundantly
    rdata.resize(0);
    LINALG::Gather<int>(sdata, rdata, (int)allproc.size(), &allproc[0], Comm());

    // add my own slave column node ids (non-redundant, standard overlap)
    const Epetra_Map* elecolmap = Discret().ElementColMap();
    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find element with gid %", gid);
      MortarElement* mrtrele = dynamic_cast<MortarElement*>(ele);
      if (mrtrele->IsSlave()) rdata.push_back(gid);
    }

    // build new element column map (on ALL processors)
    Teuchos::RCP<Epetra_Map> newelecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, Comm()));
    sdata.clear();
    rdata.clear();
    allproc.clear();

    // redistribute the discretization of the interface according to the
    // new node / element column layout (i.e. master = full overlap)
    Discret().ExportColumnNodes(*newnodecolmap);
    Discret().ExportColumnElements(*newelecolmap);
  }

  //*****NON-REDUNDANT STORAGE*****
  else if (Redundant() == INPAR::MORTAR::redundant_none)
  {
    // nothing to do here, we work with the given non-redundant distribution
    // of both slave and master nodes to the individual processsors. However
    // we want ALL procs to be part of the interface discretization, not only
    // the ones that do own or ghost any nodes in the natural distribution of
    // idiscret_. This makes dynamic redistribution easier!

    //**********************************************************************
    // IMPORTANT NOTE:
    // In an older code version, we only did ghosting on procs that own or ghost
    // any of the interface nodes or elements  in the natural distr. of idiscret_.
    // The corresponding code lines for creating this proc list are:
    //
    // std::vector<int> stproc(0);
    // if (oldnodecolmap_->NumMyElements() || oldelecolmap_->NumMyElements())
    //   stproc.push_back(Comm().MyPID());
    // std::vector<int> rtproc(0);
    // LINALG::Gather<int>(stproc,rtproc,Comm().NumProc(),&allproc[0],Comm());
    //
    // In this case, we use "rtproc" instead of "allproc" afterwards, i.e. when
    // the node gids and element gids are gathered among procs.
    //**********************************************************************

    // we keep the current ghosting, but we want to (formally) include
    // all processors, even if they do not own or ghost a single node or
    // element in the natural distribution of idiscret_
    std::vector<int> rdata;

    // fill my own slave and master column node ids (non-redundant)
    const Epetra_Map* nodecolmap = Discret().NodeColMap();
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      rdata.push_back(gid);
    }

    // re-build node column map (now formally on ALL processors)
    Teuchos::RCP<Epetra_Map> newnodecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, Comm()));
    rdata.clear();

    // fill my own slave and master column element ids (non-redundant)
    const Epetra_Map* elecolmap = Discret().ElementColMap();
    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      rdata.push_back(gid);
    }

    // re-build element column map (now formally on ALL processors)
    Teuchos::RCP<Epetra_Map> newelecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, Comm()));
    rdata.clear();

    // redistribute the discretization of the interface according to the
    // new (=old) node / element column layout
    Discret().ExportColumnNodes(*newnodecolmap);
    Discret().ExportColumnElements(*newelecolmap);

    const INPAR::MORTAR::GhostingStrategy ghostingStrategy =
        DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
            DRT::Problem::Instance()->MortarCouplingParams().sublist("PARALLEL REDISTRIBUTION"),
            "GHOSTING_STRATEGY");

    if (ghostingStrategy == INPAR::MORTAR::binningstrategy)
    {
      /* We have to update the row/column maps split into master/slave. We start from the new
       * node/element column maps. Since we don't have row maps at this point, we can/have to pass
       * the column map as row map.
       */
      UpdateMasterSlaveElementMaps(*newelecolmap, *newelecolmap);
      UpdateMasterSlaveNodeMaps(*newnodecolmap, *newnodecolmap);

      /* Ask the discretization to initialize the elements. We need this, since the setup of the
       * binning strategy relies on some element information.
       * Note: This should be cheap, since we don't assign new degrees of freedom.
       */
      Discret().FillComplete(false, true, false);

      // Create the binning strategy
      const double vel = 0.0;
      Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy = SetupBinningStrategy(vel);

      // fill master and slave elements into bins
      std::map<int, std::set<int>> slavebinelemap;
      binningstrategy->DistributeElesToBins(Discret(), slavebinelemap, true);
      std::map<int, std::set<int>> masterbinelemap;
      binningstrategy->DistributeElesToBins(Discret(), masterbinelemap, false);

      // Extend ghosting of the master elements
      Teuchos::RCP<Epetra_Map> extendedmastercolmap = binningstrategy->ExtendGhosting(
          Discret(), *newelecolmap, slavebinelemap, masterbinelemap);

      // adapt layout to extended ghosting in the discretization
      // first export the elements according to the processor local element column maps
      Discret().ExportColumnElements(*extendedmastercolmap);

      // get the node ids of the elements that are to be ghosted and create a proper node column map
      // for their export
      std::set<int> nodes;
      for (int lid = 0; lid < extendedmastercolmap->NumMyElements(); ++lid)
      {
        DRT::Element* ele = Discret().gElement(extendedmastercolmap->GID(lid));
        const int* nodeids = ele->NodeIds();
        for (int inode = 0; inode < ele->NumNode(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), &colnodes[0], 0, Comm()));

      // now ghost the nodes
      Discret().ExportColumnNodes(*nodecolmap);
    }
  }

  //*****INVALID CASES*****
  else
  {
    dserror("ERROR: ExtendInterfaceGhosting: Invalid redundancy type.");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::CreateSearchTree()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

    // warning
#ifdef MORTARGMSHCTN
  if (Dim() == 3 && Comm().MyPID() == 0)
  {
    std::cout << "\n*****************************************************************\n";
    std::cout << "GMSH output of all mortar tree nodes in 3D needs a lot of memory!\n";
    std::cout << "*****************************************************************\n";
  }
#endif  // MORTARGMSHCTN
  // binary tree search
  if (SearchAlg() == INPAR::MORTAR::search_binarytree)
  {
    // create fully overlapping map of all master elements
    // for non-redundant storage (RRloop) we handle the master elements
    // like the slave elements --> melecolmap_
    INPAR::MORTAR::GhostingStrategy strat =
        DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
            InterfaceParams().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");

    // get update type of binary tree
    INPAR::MORTAR::BinaryTreeUpdateType updatetype =
        DRT::INPUT::IntegralValue<INPAR::MORTAR::BinaryTreeUpdateType>(
            InterfaceParams(), "BINARYTREE_UPDATETYPE");

    Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;
    if (strat == INPAR::MORTAR::roundrobinghost || strat == INPAR::MORTAR::binningstrategy)
      melefullmap = melecolmap_;
    else if (strat == INPAR::MORTAR::ghosting_redundant)
      melefullmap = LINALG::AllreduceEMap(*melerowmap_);
    else
      dserror("Chosen parallel strategy not supported!");

    // create binary tree object for search and setup tree
    binarytree_ = Teuchos::rcp(new MORTAR::BinaryTree(
        Discret(), selecolmap_, melefullmap, Dim(), SearchParam(), updatetype, SearchUseAuxPos()));
    // initialize the binary tree
    binarytree_->Init();
  }
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                 popp 11/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveSets()
{
  UpdateMasterSlaveNodeMaps();
  UpdateMasterSlaveElementMaps();
  UpdateMasterSlaveDofMaps();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveDofMaps()
{
  // Vectors to collect GIDs to build maps
  std::vector<int> sc;  // slave column map
  std::vector<int> sr;  // slave row map
  std::vector<int> mc;  // master column map
  std::vector<int> mr;  // master row map

  const int numMyColumnDofs = Discret().NodeColMap()->NumMyElements();
  for (int i = 0; i < numMyColumnDofs; ++i)
  {
    int gid = Discret().NodeColMap()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);
    bool isslave = mrtrnode->IsSlave();

    if (isslave)
      for (int j = 0; j < mrtrnode->NumDof(); ++j) sc.push_back(mrtrnode->Dofs()[j]);
    else
      for (int j = 0; j < mrtrnode->NumDof(); ++j) mc.push_back(mrtrnode->Dofs()[j]);

    if (Discret().NodeRowMap()->MyGID(gid))
    {
      if (isslave)
        for (int j = 0; j < mrtrnode->NumDof(); ++j) sr.push_back(mrtrnode->Dofs()[j]);
      else
        for (int j = 0; j < mrtrnode->NumDof(); ++j) mr.push_back(mrtrnode->Dofs()[j]);
    }
  }

  sdofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sr.size(), &sr[0], 0, Comm()));
  sdofcolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sc.size(), &sc[0], 0, Comm()));
  mdofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mr.size(), &mr[0], 0, Comm()));
  mdofcolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mc.size(), &mc[0], 0, Comm()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveElementMaps()
{
  UpdateMasterSlaveElementMaps(*Discret().ElementRowMap(), *Discret().ElementColMap());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveElementMaps(
    const Epetra_Map& elementRowMap, const Epetra_Map& elementColumnMap)
{
  // Vectors to collect GIDs to build maps
  std::vector<int> sc;  // slave column map
  std::vector<int> sr;  // slave row map
  std::vector<int> mc;  // master column map
  std::vector<int> mr;  // master row map

  const int numMyColumnElements = elementColumnMap.NumMyElements();
  for (int i = 0; i < numMyColumnElements; ++i)
  {
    int gid = elementColumnMap.GID(i);
    bool isslave = dynamic_cast<MORTAR::MortarElement*>(Discret().gElement(gid))->IsSlave();

    if (isslave)
      sc.push_back(gid);
    else
      mc.push_back(gid);

    if (elementRowMap.MyGID(gid))
    {
      if (isslave)
        sr.push_back(gid);
      else
        mr.push_back(gid);
    }
  }

  selerowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sr.size(), &sr[0], 0, Comm()));
  selecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sc.size(), &sc[0], 0, Comm()));
  melerowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mr.size(), &mr[0], 0, Comm()));
  melecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mc.size(), &mc[0], 0, Comm()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveNodeMaps()
{
  UpdateMasterSlaveNodeMaps(*Discret().NodeRowMap(), *Discret().NodeColMap());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::UpdateMasterSlaveNodeMaps(
    const Epetra_Map& nodeRowMap, const Epetra_Map& nodeColumnMap)
{
  // Vectors to collect GIDs to build maps
  std::vector<int> sc;   // slave column map
  std::vector<int> sr;   // slave row map
  std::vector<int> mc;   // master column map
  std::vector<int> mr;   // master row map
  std::vector<int> srb;  // slave row map + boundary nodes
  std::vector<int> scb;  // slave column map + boundary nodes
  std::vector<int> mrb;  // master row map - boundary nodes
  std::vector<int> mcb;  // master column map - boundary nodes

  const int numMyColumnNodes = nodeColumnMap.NumMyElements();
  for (int i = 0; i < numMyColumnNodes; ++i)
  {
    int gid = nodeColumnMap.GID(i);
    bool isslave = dynamic_cast<MORTAR::MortarNode*>(Discret().gNode(gid))->IsSlave();
    bool isonbound = dynamic_cast<MORTAR::MortarNode*>(Discret().gNode(gid))->IsOnBoundorCE();

    if (isslave || isonbound)
      scb.push_back(gid);
    else
      mcb.push_back(gid);
    if (isslave)
      sc.push_back(gid);
    else
      mc.push_back(gid);

    if (nodeRowMap.MyGID(gid))
    {
      if (isslave || isonbound)
        srb.push_back(gid);
      else
        mrb.push_back(gid);
      if (isslave)
        sr.push_back(gid);
      else
        mr.push_back(gid);
    }
  }

  snoderowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sr.size(), &sr[0], 0, Comm()));
  snodecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sc.size(), &sc[0], 0, Comm()));
  mnoderowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mr.size(), &mr[0], 0, Comm()));
  mnodecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mc.size(), &mc[0], 0, Comm()));

  snoderowmapbound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)srb.size(), &srb[0], 0, Comm()));
  snodecolmapbound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)scb.size(), &scb[0], 0, Comm()));
  mnoderowmapnobound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mrb.size(), &mrb[0], 0, Comm()));
  mnodecolmapnobound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mcb.size(), &mcb[0], 0, Comm()));

  // build exporter
  interfaceData_->SlExporterPtr() =
      Teuchos::rcp(new DRT::Exporter(*snoderowmapbound_, *snodecolmapbound_, Comm()));
}

/*----------------------------------------------------------------------*
 |  restrict slave sets to actual meshtying zone              popp 08/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::RestrictSlaveSets()
{
  //********************************************************************
  // NODES
  //********************************************************************
  {
    std::vector<int> sc;      // slave column map
    std::vector<int> sr;      // slave row map
    std::vector<int> scfull;  // slave full map

    for (int i = 0; i < snodecolmap_->NumMyElements(); ++i)
    {
      int gid = snodecolmap_->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);
      int istied = (int)mrtrnode->IsTiedSlave();

      if (istied && snodecolmap_->MyGID(gid)) sc.push_back(gid);
      if (istied && snoderowmap_->MyGID(gid)) sr.push_back(gid);
    }

    snoderowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sr.size(), &sr[0], 0, Comm()));
    snodecolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sc.size(), &sc[0], 0, Comm()));
  }

  //********************************************************************
  // ELEMENTS
  //********************************************************************
  // no need to do this for elements, because all mortar quantities
  // are defined with respect to node or dof maps (D,M,...). As all
  // mortar stuff has already been evaluated, it would not matter if
  // we adapted the element maps as well, but we just skip it.

  //********************************************************************
  // DOFS
  //********************************************************************
  {
    std::vector<int> sc;      // slave column map
    std::vector<int> sr;      // slave row map
    std::vector<int> scfull;  // slave full map

    for (int i = 0; i < snodecolmap_->NumMyElements(); ++i)
    {
      int gid = snodecolmap_->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);
      int istied = (int)mrtrnode->IsTiedSlave();

      if (istied && snodecolmap_->MyGID(gid))
        for (int j = 0; j < mrtrnode->NumDof(); ++j) sc.push_back(mrtrnode->Dofs()[j]);

      if (istied && snoderowmap_->MyGID(gid))
        for (int j = 0; j < mrtrnode->NumDof(); ++j) sr.push_back(mrtrnode->Dofs()[j]);
    }

    sdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sr.size(), &sr[0], 0, Comm()));
    sdofcolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sc.size(), &sc[0], 0, Comm()));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  update Lagrange multiplier set (dofs)                     popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> MORTAR::MortarInterface::UpdateLagMultSets(
    int offset_if, const bool& redistributed, const Epetra_Map& ref_map) const
{
  if (redistributed)
  {
    return RedistributeLagMultSets();
  }
  //********************************************************************
  // LAGRANGE MULTIPLIER DOFS
  //********************************************************************
  // NOTE: we want no gap between the displacement dofs and the newly
  // defined Lagrange multiplier dofs!! Thus, if the maximum displacement
  // dof is 12.345, we want the LM dofs to start with 12.346. This can
  // be readily achieved, because we know that the lmdofmap will have
  // the same parallel distribution as the slavedofrowmap. The only
  // thing we need to take care of is to avoid overlapping of the LM
  // dofs among different processors. Therefore, the total number of
  // slave nodes (and thus LM nodes) of each processor is communicated
  // to ALL other processors and an offset is then determined for each
  // processor based on this information.
  //********************************************************************
  // temporary vector of LM dofs
  std::vector<int> lmdof;

  // gather information over all procs
  std::vector<int> localnumlmdof(Comm().NumProc());
  std::vector<int> globalnumlmdof(Comm().NumProc());
  localnumlmdof[Comm().MyPID()] = ref_map.NumMyElements();
  Comm().SumAll(&localnumlmdof[0], &globalnumlmdof[0], Comm().NumProc());

  // compute offset for LM dof initialization for all procs
  int offset = 0;
  for (int k = 0; k < Comm().MyPID(); ++k) offset += globalnumlmdof[k];

  // loop over all slave dofs and initialize LM dofs
  for (int i = 0; i < ref_map.NumMyElements(); ++i)
    lmdof.push_back(MaxDofGlobal() + 1 + offset_if + offset + i);

  // create interface LM map
  // (if maxdofglobal_ == 0, we do not want / need this)
  if (MaxDofGlobal() > 0)
    return Teuchos::rcp(new Epetra_Map(-1, (int)lmdof.size(), &lmdof[0], 0, Comm()));

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::StoreUnredistributedMaps()
{
  psdofrowmap_ = Teuchos::rcp(new Epetra_Map(*sdofrowmap_));
  interfaceData_->PMDofRowMap() = Teuchos::rcp(new Epetra_Map(*mdofrowmap_));
  plmdofmap_ = Teuchos::rcp(new Epetra_Map(*lmdofmap_));

  interfaceData_->PSNodeRowMap() = Teuchos::rcp(new Epetra_Map(*snoderowmap_));
  interfaceData_->PMNodeRowMap() = Teuchos::rcp(new Epetra_Map(*mnoderowmap_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> MORTAR::MortarInterface::RedistributeLagMultSets() const
{
  if (plmdofmap_.is_null()) dserror("The plmdofmap_ is not yet initialized!");
  if (psdofrowmap_.is_null()) dserror("The psdofrowmap_ is not yet initialized!");

  // new lm dofs
  std::vector<int> lmdof(sdofrowmap_->NumMyElements());

  /* get the initial correlation between the slave dofs
   * and the lagrange multiplier dofs
   *
   * There is always a tuple of two values which correlate
   * with each other. The first one is the lm-dof, the second
   * one is the slave dof. */
  std::vector<int> my_related_gids(2 * plmdofmap_->NumMyElements(), -1);
  for (int i = 0; i < plmdofmap_->NumMyElements(); ++i)
  {
    my_related_gids[2 * i] = plmdofmap_->GID(i);
    my_related_gids[2 * i + 1] = psdofrowmap_->GID(i);
  }

  std::vector<int> g_related_gids;
  Comm().Barrier();
  for (int p = 0; p < Comm().NumProc(); ++p)
  {
    int num_mygids = plmdofmap_->NumMyElements();

    Comm().Broadcast(&num_mygids, 1, p);
    // skip processors which hold no correlation info
    if (num_mygids == 0) continue;
    g_related_gids.resize(2 * num_mygids);

    /* communicate the correlation list of proc p
     * to all procs */
    if (p == Comm().MyPID())
      for (std::size_t i = 0; i < my_related_gids.size(); ++i)
        g_related_gids[i] = my_related_gids[i];
    Comm().Broadcast(&g_related_gids[0], 2 * num_mygids, p);

    for (int i = 0; i < num_mygids; ++i)
    {
      /* check in the already redistributed sdofrowmap on
       * each processor which one holds the current gid */
      int my_sllid = sdofrowmap_->LID(g_related_gids[2 * i + 1]);
      /* on the proc holding the gid, we store the corresponding
       * lm-dof-gid as well at the same lid. */
      if (my_sllid != -1) lmdof[my_sllid] = g_related_gids[2 * i];
    }
    // wait for the arrival of all procs
    Comm().Barrier();
  }

  // create deterministic interface LM map
  return Teuchos::rcp(new Epetra_Map(-1, (int)lmdof.size(), &lmdof[0], 0, Comm()));
}

/*----------------------------------------------------------------------*
 |  initialize / reset mortar interface                       popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over all nodes to reset stuff (fully overlapping column map)
  for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
  {
    MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lColNode(i));

    // reset feasible projection and segmentation status
    node->HasProj() = false;
    node->HasSegment() = false;
  }

  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i = 0; i < SlaveColNodesBound()->NumMyElements(); ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MORTAR::MortarNode* monode = dynamic_cast<MORTAR::MortarNode*>(node);

    // reset nodal normal
    for (int j = 0; j < 3; ++j) monode->MoData().n()[j] = 0.0;

    // reset nodal Mortar maps
    monode->MoData().GetD().clear();
    monode->MoData().GetM().clear();
    monode->MoData().GetMmod().clear();
  }

  // loop over all elements to reset candidates / search lists
  // (use standard slave column map)
  for (int i = 0; i < SlaveColElements()->NumMyElements(); ++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MortarElement* mele = dynamic_cast<MortarElement*>(ele);

    mele->MoData().SearchElements().resize(0);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  set current and old deformation state                      popp 12/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::SetState(const enum StateType& statetype, const Epetra_Vector& vec)
{
  // ***WARNING:*** This is commented out here, as idiscret_->SetState()
  // needs all the procs around, not only the interface local ones!
  // if (!lComm()) return;

  switch (statetype)
  {
    case state_new_displacement:
    {
      // alternative method to get vec to full overlap
      Teuchos::RCP<Epetra_Vector> global =
          Teuchos::rcp(new Epetra_Vector(*idiscret_->DofColMap(), false));
      LINALG::Export(vec, *global);

      // set displacements in interface discretization
      idiscret_->SetState(StateType2String(statetype), global);

      // loop over all nodes to set current displacement
      // (use fully overlapping column map)
      for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
      {
        MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lColNode(i));
        const int numdof = node->NumDof();
        std::vector<double> mydisp(numdof);
        std::vector<int> lm(numdof);

        for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

        DRT::UTILS::ExtractMyValues(*global, mydisp, lm);

        // add mydisp[2]=0 for 2D problems
        if (mydisp.size() < 3) mydisp.resize(3);

        // set current configuration
        for (int j = 0; j < 3; ++j) node->xspatial()[j] = node->X()[j] + mydisp[j];
      }

      // compute element areas
      SetElementAreas();
      break;
    }
    case state_lagrange_multiplier:
    {
      // alternative method to get vec to full overlap
      Teuchos::RCP<Epetra_Vector> global =
          Teuchos::rcp(new Epetra_Vector(*idiscret_->DofColMap(), false));
      LINALG::Export(vec, *global);

      // loop over all nodes to set current displacement
      // (use fully overlapping column map)
      for (int i = 0; i < SlaveColNodes()->NumMyElements(); ++i)
      {
        MORTAR::MortarNode* node =
            dynamic_cast<MORTAR::MortarNode*>(idiscret_->gNode(SlaveColNodes()->GID(i)));
        const int numdof = node->NumDof();
        std::vector<double> mydisp(numdof);
        std::vector<int> lm(numdof);

        for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

        DRT::UTILS::ExtractMyValues(*global, mydisp, lm);

        // add mydisp[2]=0 for 2D problems
        if (mydisp.size() < 3) mydisp.resize(3);

        // set current configuration
        for (int j = 0; j < 3; ++j) node->MoData().lm()[j] = mydisp[j];
      }
      break;
    }
    case state_old_displacement:
    {
      // alternative method to get vec to full overlap
      Teuchos::RCP<Epetra_Vector> global =
          Teuchos::rcp(new Epetra_Vector(*idiscret_->DofColMap(), false));
      LINALG::Export(vec, *global);

      // set displacements in interface discretization
      idiscret_->SetState(StateType2String(statetype), global);

      // loop over all nodes to set current displacement
      // (use fully overlapping column map)
      for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
      {
        MORTAR::MortarNode* node = dynamic_cast<MORTAR::MortarNode*>(idiscret_->lColNode(i));
        const int numdof = node->NumDof();
        std::vector<double> myolddisp(numdof);
        std::vector<int> lm(numdof);

        for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

        DRT::UTILS::ExtractMyValues(*global, myolddisp, lm);

        // add mydisp[2]=0 for 2D problems
        if (myolddisp.size() < 3) myolddisp.resize(3);

        // set old displacement
        for (int j = 0; j < 3; ++j) node->uold()[j] = myolddisp[j];
      }

      break;
    }
    default:
    {
      dserror(
          "The given state type is unsupported! (type = %s)", StateType2String(statetype).c_str());
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::SetElementAreas()
{
  // loop over all elements to set current element length / area
  // (use standard slave column map)
  for (int i = 0; i < SlaveColElements()->NumMyElements(); ++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i", gid);
    MortarElement* mele = dynamic_cast<MortarElement*>(ele);

    mele->MoData().Area() = mele->ComputeArea();
  }
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate geometric setting (create integration cells)    farah 01/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateGeometry(std::vector<Teuchos::RCP<MORTAR::IntCell>>& intcells)
{
  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // check
  if (Dim() == 2) dserror("ERROR: Geometry evaluation for mortar interface only for 3D problems!");

  INPAR::MORTAR::AlgorithmType algo =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  if (algo == INPAR::MORTAR::algorithm_nts)
    dserror("ERROR: Geometry evaluation only for mortar problems!");

  // interface needs to be complete
  if (!Filled() && Comm().MyPID() == 0)
    dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm()) return;

  // clear vector
  intcells.clear();

  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("ERROR: Invalid search algorithm");

  // create normals
  EvaluateNodalNormals();

  // export nodal normals to slave node column map
  // this call is very expensive and the computation
  // time scales directly with the proc number !
  ExportNodalNormals();

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %", gid1);
    MortarElement* selement = dynamic_cast<MortarElement*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %", gid2);
      MortarElement* melement = dynamic_cast<MortarElement*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      //********************************************************************
      if (selement->IsQuad())
      {
        dserror("ERROR: Geometry evaluation only implemented for first order elements!");
      }
      // noquad!
      else
      {
        MORTAR::Coupling3d coup(*idiscret_, dim_, false, imortar_, *selement, *melement);

        // do coupling
        coup.EvaluateCoupling();

        // set sele and mele id and push into global vector
        for (size_t c = 0; c < coup.Cells().size(); ++c)
        {
          coup.Cells()[c]->SetSlaveId(selement->Id());
          coup.Cells()[c]->SetMasterId(melement->Id());
          intcells.push_back(coup.Cells()[c]);
        }
      }
    }
  }  // end sele loop

  // time measurement
  Comm().Barrier();
  const double evaltime = Teuchos::Time::wallTime() - t_start;

  // time output
  std::cout << "Required time for geometry evaluation: " << evaltime << std::endl;

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate mortar coupling (public)                         popp 11/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::Evaluate(
    int rriter, const int& step, const int& iter, Teuchos::RCP<MORTAR::ParamsInterface> mparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR("MORTAR::MortarInterface::Evaluate");

  // interface needs to be complete
  if (!Filled() && Comm().MyPID() == 0)
    dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm()) return;

  //******************************************
  // Start basic evaluation part of interface
  //******************************************
  // start time measurement
  const double t_start = Teuchos::Time::wallTime();

  // evaluate nodal normals and decide in contact case if
  // this is a nonsmooth or smooth contact
  PreEvaluate(step, iter);

  // evaluation routine for coupling
  EvaluateCoupling(*selecolmap_, snoderowmap_.get(), mparams_ptr);

  // do some post operations. nothing happens for standard cases...
  PostEvaluate(step, iter);

  // end time on this proc
  const double inttime = Teuchos::Time::wallTime() - t_start;
  //******************************************
  // End basic evaluation part of interface
  //******************************************

  // store integrationtime
  inttime_interface_ = inttime;

  // bye
  return;
}


/*----------------------------------------------------------------------*
 |  protected evaluate routine                               farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateCoupling(const Epetra_Map& selecolmap,
    const Epetra_Map* snoderowmap, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // decide which type of coupling should be evaluated
  INPAR::MORTAR::AlgorithmType algo =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

  // smooth contact
  switch (algo)
  {
    //*********************************
    // Mortar Coupling (STS)    (2D/3D)
    // Gauß-Point-To-Segment (GPTS)
    //*********************************
    case INPAR::MORTAR::algorithm_mortar:
    case INPAR::MORTAR::algorithm_gpts:
    {
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      MORTAR::MortarInterface::EvaluateSTS(selecolmap, mparams_ptr);
      break;
    }
    //*********************************
    // Segment-to-Line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_stl:
    {
      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateSTL();
      break;
    }
    //*********************************
    // Line-to-Segment Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_lts:
    {
      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateLTS();
      break;
    }
    //*********************************
    // line-to-line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_ltl:
    {
      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateLTL();
      break;
    }
    //*********************************
    // Node-to-Segment Coupling (2D/3D)
    //*********************************
    case INPAR::MORTAR::algorithm_nts:
    {
      //********************************************************************
      // 1) try to project slave nodes onto master elements
      // 2) evaluate shape functions at projected positions
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateNTS();
      break;
    }
    //*********************************
    // Node-to-Line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_ntl:
    {
      dserror("ERROR: not yet implemented!");
      break;
    }
    //*********************************
    // Default case
    //*********************************
    default:
    {
      dserror("ERROR: Unknown discr. type for constraints!");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-segment coupl          farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateSTS(
    const Epetra_Map& selecolmap, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR("MORTAR::MortarInterface::EvaluateSTS");

  // loop over all slave col elements
  for (int i = 0; i < selecolmap.NumMyElements(); ++i)
  {
    const int gid1 = selecolmap.GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %d", gid1);

    MortarElement* selement = dynamic_cast<MortarElement*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized()) continue;

    // empty vector of master element pointers
    std::vector<MortarElement*> melements;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %d", gid2);
      MortarElement* melement = dynamic_cast<MortarElement*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      melements.push_back(melement);
    }

    // concrete coupling evaluation routine
    MortarCoupling(selement, melements, mparams_ptr);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type node-to-segment coupl             farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateNTS()
{
  // loop over slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) dserror("ERROR: Node ownership inconsistency!");

    // vector with possible contacting master eles
    std::vector<MORTAR::MortarElement*> meles;

    // fill vector with possibly contacting meles
    FindMEles(*mrtrnode, meles);

    // skip calculation if no meles vector is empty
    if (meles.size() < 1) continue;

    // call interpolation functions
    NTS::MTInterpolator::Impl(meles)->Interpolate(*mrtrnode, meles);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateLTS()
{
  dserror("ERROR: LTS NOT FOR MESHTYING!!!!");
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-line coupl                farah 07/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateLTL()
{
  dserror("ERROR: LTL NOT FOR MESHTYING!!!!");
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-line coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateSTL()
{
  dserror("ERROR: LTS NOT FOR MESHTYING!!!!");
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate nodal normals (public)                           popp 10/11|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateNodalNormals() const
{
  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();
  }

  return;
}

/*----------------------------------------------------------------------*
 |  pre evaluate to calc normals                            farah 02/16 |
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::PreEvaluate(const int& step, const int& iter)
{
  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("ERROR: Invalid search algorithm");

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

  // evaluate averaged nodal normals on slave side
  EvaluateNodalNormals();

  // export nodal normals to slave node column map
  // this call is very expensive and the computation
  // time scales directly with the proc number !
  ExportNodalNormals();

  return;
}


/*----------------------------------------------------------------------*
 |  post evaluate                                           farah 02/16 |
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::PostEvaluate(const int step, const int iter)
{
  // nothing to do...

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
 |  find meles to snode                                     farah 01/16 |
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::FindMEles(
    MortarNode& mrtrnode, std::vector<MortarElement*>& meles) const
{
  // clear vector
  meles.clear();

  // get adjacent elements for this node
  DRT::Element** adjeles = mrtrnode.Elements();

  // empty vector of master element pointers
  std::set<int> donebefore;

  for (int j = 0; j < mrtrnode.NumElement(); ++j)
  {
    MortarElement* adjcele = dynamic_cast<MortarElement*>(adjeles[j]);

    // skip zero-sized nurbs elements (slave)
    if (adjcele->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < adjcele->MoData().NumSearchElements(); ++k)
    {
      int gid2 = adjcele->MoData().SearchElements()[k];
      DRT::Element* mele = idiscret_->gElement(gid2);
      if (!mele) dserror("ERROR: Cannot find master element with gid %", gid2);
      MortarElement* melement = dynamic_cast<MortarElement*>(mele);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      // check uniqueness
      std::set<int>::iterator iter = donebefore.find(melement->Id());
      if (iter != donebefore.end()) continue;

      donebefore.insert(melement->Id());

      // fill vector
      meles.push_back(melement);
    }  // found eles
  }    // loop over adjacent slave elements

  return;
}


/*----------------------------------------------------------------------*
 |  find mnodes to snode                                    farah 01/16 |
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::FindMNodes(
    MortarNode& mrtrnode, std::vector<MortarElement*>& meles, std::vector<MortarNode*>& mnodes)
{
  // clear vector
  mnodes.clear();

  // check meles
  if (meles.size() < 1) return;

  // set object to guarantee uniqueness of found mnodes
  std::set<int> donebefore;

  for (size_t j = 0; j < meles.size(); ++j)
  {
    MortarElement* mele = dynamic_cast<MortarElement*>(meles[j]);

    // skip zero-sized nurbs elements (master)
    if (mele->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < mele->NumNode(); ++k)
    {
      DRT::Node* node = mele->Nodes()[k];
      if (!node) dserror("ERROR: Cannot find master node");
      MortarNode* mnode = dynamic_cast<MortarNode*>(node);

      // check uniqueness
      std::set<int>::iterator iter = donebefore.find(mnode->Id());
      if (iter != donebefore.end()) continue;

      donebefore.insert(mnode->Id());

      // fill vector
      mnodes.push_back(mnode);
    }  // found eles
  }    // loop over adjacent slave elements

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate nodal normals and store them in map (public)      jb 07/14 |
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateNodalNormals(std::map<int, std::vector<double>>& mynormals)
{
  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();

    int numdofs = mrtrnode->NumDof();
    std::vector<double> temp(numdofs, 0.0);
    for (int i = 0; i < numdofs; i++)
    {
      temp[i] = mrtrnode->MoData().n()[i];
    }
    mynormals.insert(std::pair<int, std::vector<double>>(gid, temp));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  export nodal normals (public)                             popp 11/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::ExportNodalNormals() const
{
  // create empty data objects
  std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>> triad;

  // build info on row map
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    // fill nodal matrix
    Teuchos::RCP<Epetra_SerialDenseMatrix> loc = Teuchos::rcp(new Epetra_SerialDenseMatrix(3, 1));
    (*loc)(0, 0) = mrtrnode->MoData().n()[0];
    (*loc)(1, 0) = mrtrnode->MoData().n()[1];
    (*loc)(2, 0) = mrtrnode->MoData().n()[2];

    triad[gid] = loc;
  }

  // communicate from slave node row to column map

  interfaceData_->Exporter().Export(triad);

  // extract info on column map
  for (int i = 0; i < snodecolmapbound_->NumMyElements(); ++i)
  {
    // only do something for ghosted nodes
    int gid = snodecolmapbound_->GID(i);
    if (snoderowmapbound_->MyGID(gid)) continue;

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    // extract info
    Teuchos::RCP<Epetra_SerialDenseMatrix> loc = triad[gid];
    mrtrnode->MoData().n()[0] = (*loc)(0, 0);
    mrtrnode->MoData().n()[1] = (*loc)(1, 0);
    mrtrnode->MoData().n()[2] = (*loc)(2, 0);
  }

  // free memory
  triad.clear();

  /*// print nodal normals
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
   if (!node) dserror("ERROR: Cannot find node with gid %",gid);
   MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

   // print averaged normal at each slave node
   std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << mrtrnode->Owner()
   << " Normal: " << mrtrnode->MoData().n()[0]
   << " " << mrtrnode->MoData().n()[1] << " " << mrtrnode->MoData().n()[2] << std::endl;
   }
   std::cout << std::endl << std::endl;
   }

   // barrier
   Comm().Barrier();
   }*/

  return;
}

/*----------------------------------------------------------------------*
 |  Search element-based "brute force" (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::EvaluateSearchBruteForce(const double& eps)
{
  /**********************************************************************/
  /* SEARCH ALGORITHM:                                                  */
  /* The idea of the search is to reduce the number of master / slave   */
  /* element pairs that are checked for overlap and coupling by intro-  */
  /* ducing information about proximity and maybe history!              */
  /* This old version is already based on bounding volumes, but still   */
  /* brute force for finding the overlap of these bounding volumes,     */
  /* so it has been replaced by a more efficient approach (binary tree).*/
  /**********************************************************************/

  // calculate minimal element length
  double lmin = 1.0e12;
  double enlarge = 0.0;

  // create fully overlapping map of all master elements
  // for non-redundant storage (RRloop) we handle the master elements
  // like the slave elements --> melecolmap_
  INPAR::MORTAR::GhostingStrategy strat =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::GhostingStrategy>(
          InterfaceParams().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");
  Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;

  if (strat == INPAR::MORTAR::ghosting_redundant)
    melefullmap = LINALG::AllreduceEMap(*melerowmap_);
  else if (strat == INPAR::MORTAR::roundrobinghost)
    melefullmap = melerowmap_;
  else if (strat == INPAR::MORTAR::binningstrategy)
    melefullmap = melecolmap_;
  else
    dserror("Choosen parallel strategy not supported!");

  // loop over all slave elements on this proc.
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    DRT::Element* element = idiscret_->gElement(selecolmap_->GID(i));
    if (!element) dserror("ERROR: Cannot find element with gid %\n", selecolmap_->GID(i));
    MORTAR::MortarElement* mrtrelement = dynamic_cast<MortarElement*>(element);
    if (mrtrelement->MinEdgeSize() < lmin) lmin = mrtrelement->MinEdgeSize();
  }

  // loop over all master elements on this proc.
  for (int i = 0; i < melefullmap->NumMyElements(); ++i)
  {
    DRT::Element* element = idiscret_->gElement(melefullmap->GID(i));
    if (!element) dserror("ERROR: Cannot find element with gid %\n", melefullmap->GID(i));
    MORTAR::MortarElement* mrtrelement = dynamic_cast<MortarElement*>(element);
    if (mrtrelement->MinEdgeSize() < lmin) lmin = mrtrelement->MinEdgeSize();
  }

  // compute DOP inflation length
  enlarge = eps * lmin;

  // define dopnormals
  Epetra_SerialDenseMatrix dopnormals;
  int kdop = 0;

  if (dim_ == 2)
  {
    kdop = 8;

    // setup normals for DOP
    dopnormals.Reshape(4, 3);
    dopnormals(0, 0) = 1;
    dopnormals(0, 1) = 0;
    dopnormals(0, 2) = 0;
    dopnormals(1, 0) = 0;
    dopnormals(1, 1) = 1;
    dopnormals(1, 2) = 0;
    dopnormals(2, 0) = 1;
    dopnormals(2, 1) = 1;
    dopnormals(2, 2) = 0;
    dopnormals(3, 0) = -1;
    dopnormals(3, 1) = 1;
    dopnormals(3, 2) = 0;
  }
  else if (dim_ == 3)
  {
    kdop = 18;

    // setup normals for DOP
    dopnormals.Reshape(9, 3);
    dopnormals(0, 0) = 1;
    dopnormals(0, 1) = 0;
    dopnormals(0, 2) = 0;
    dopnormals(1, 0) = 0;
    dopnormals(1, 1) = 1;
    dopnormals(1, 2) = 0;
    dopnormals(2, 0) = 0;
    dopnormals(2, 1) = 0;
    dopnormals(2, 2) = 1;
    dopnormals(3, 0) = 1;
    dopnormals(3, 1) = 1;
    dopnormals(3, 2) = 0;
    dopnormals(4, 0) = 1;
    dopnormals(4, 1) = 0;
    dopnormals(4, 2) = 1;
    dopnormals(5, 0) = 0;
    dopnormals(5, 1) = 1;
    dopnormals(5, 2) = 1;
    dopnormals(6, 0) = 1;
    dopnormals(6, 1) = 0;
    dopnormals(6, 2) = -1;
    dopnormals(7, 0) = -1;
    dopnormals(7, 1) = 1;
    dopnormals(7, 2) = 0;
    dopnormals(8, 0) = 0;
    dopnormals(8, 1) = -1;
    dopnormals(8, 2) = 1;
  }
  else
    dserror("ERROR: Problem dimension must be 2D or 3D!");

  // decide whether auxiliary positions are used when computing dops
  const bool useauxpos = SearchUseAuxPos();

  // define slave and master slabs
  Epetra_SerialDenseMatrix sslabs(kdop / 2, 2);
  Epetra_SerialDenseMatrix mslabs(kdop / 2, 2);

  //**********************************************************************
  // perform brute-force search (element-based)
  //**********************************************************************
  // for every slave element
  for (int i = 0; i < selecolmap_->NumMyElements(); i++)
  {
    // calculate slabs
    double dcurrent = 0.0;

    // initialize slabs with first node
    int sgid = selecolmap_->GID(i);
    DRT::Element* element = idiscret_->gElement(sgid);
    if (!element) dserror("ERROR: Cannot find element with gid %\n", sgid);
    DRT::Node** node = element->Nodes();
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node[0]);
    const double* posnode = mrtrnode->xspatial();

    // calculate slabs initialization
    for (int j = 0; j < kdop / 2; j++)
    {
      //= ax+by+cz=d/sqrt(aa+bb+cc)
      sslabs(j, 0) = sslabs(j, 1) =
          (dopnormals(j, 0) * posnode[0] + dopnormals(j, 1) * posnode[1] +
              dopnormals(j, 2) * posnode[2]) /
          sqrt((dopnormals(j, 0) * dopnormals(j, 0)) + (dopnormals(j, 1) * dopnormals(j, 1)) +
               (dopnormals(j, 2) * dopnormals(j, 2)));
    }

    // for int j=1, because of initialization done before
    for (int j = 1; j < element->NumNode(); j++)
    {
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node[j]);
      posnode = mrtrnode->xspatial();

      for (int k = 0; k < kdop / 2; k++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        dcurrent =
            (dopnormals(k, 0) * posnode[0] + dopnormals(k, 1) * posnode[1] +
                dopnormals(k, 2) * posnode[2]) /
            sqrt((dopnormals(k, 0) * dopnormals(k, 0)) + (dopnormals(k, 1) * dopnormals(k, 1)) +
                 (dopnormals(k, 2) * dopnormals(k, 2)));
        if (dcurrent > sslabs(k, 1)) sslabs(k, 1) = dcurrent;
        if (dcurrent < sslabs(k, 0)) sslabs(k, 0) = dcurrent;
      }
    }

    // add auxiliary positions
    // (last converged positions for all slave nodes)
    if (useauxpos == true)
    {
      for (int j = 0; j < element->NumNode(); j++)
      {
        // get pointer to slave node
        MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node[j]);

        double auxpos[3] = {0.0, 0.0, 0.0};
        double scalar = 0.0;
        for (int k = 0; k < dim_; k++)
          scalar += (mrtrnode->X()[k] + mrtrnode->uold()[k] - mrtrnode->xspatial()[k]) *
                    mrtrnode->MoData().n()[k];
        for (int k = 0; k < dim_; k++)
          auxpos[k] = mrtrnode->xspatial()[k] + scalar * mrtrnode->MoData().n()[k];

        for (int j = 0; j < kdop / 2; j++)
        {
          //= ax+by+cz=d/sqrt(aa+bb+cc)
          dcurrent =
              (dopnormals(j, 0) * auxpos[0] + dopnormals(j, 1) * auxpos[1] +
                  dopnormals(j, 2) * auxpos[2]) /
              sqrt((dopnormals(j, 0) * dopnormals(j, 0)) + (dopnormals(j, 1) * dopnormals(j, 1)) +
                   (dopnormals(j, 2) * dopnormals(j, 2)));
          if (dcurrent > sslabs(j, 1)) sslabs(j, 1) = dcurrent;
          if (dcurrent < sslabs(j, 0)) sslabs(j, 0) = dcurrent;
        }
      }
    }

    // enlarge slabs with scalar factor
    for (int j = 0; j < kdop / 2; j++)
    {
      sslabs(j, 0) = sslabs(j, 0) - enlarge;
      sslabs(j, 1) = sslabs(j, 1) + enlarge;
    }

    // for every master element
    for (int j = 0; j < melefullmap->NumMyElements(); j++)
    {
      // calculate slabs
      double dcurrent = 0.0;

      // initialize slabs with first node
      int mgid = melefullmap->GID(j);
      DRT::Element* element = idiscret_->gElement(mgid);
      if (!element) dserror("ERROR: Cannot find element with gid %\n", mgid);
      DRT::Node** node = element->Nodes();
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node[0]);
      const double* posnode = mrtrnode->xspatial();

      // calculate slabs initialization
      for (int k = 0; k < kdop / 2; k++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        mslabs(k, 0) = mslabs(k, 1) =
            (dopnormals(k, 0) * posnode[0] + dopnormals(k, 1) * posnode[1] +
                dopnormals(k, 2) * posnode[2]) /
            sqrt((dopnormals(k, 0) * dopnormals(k, 0)) + (dopnormals(k, 1) * dopnormals(k, 1)) +
                 (dopnormals(k, 2) * dopnormals(k, 2)));
      }

      // for int k=1, because of initialization done before
      for (int k = 1; k < element->NumNode(); k++)
      {
        MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node[k]);
        posnode = mrtrnode->xspatial();

        for (int l = 0; l < kdop / 2; l++)
        {
          //= d=ax+by+cz/sqrt(aa+bb+cc)
          dcurrent =
              (dopnormals(l, 0) * posnode[0] + dopnormals(l, 1) * posnode[1] +
                  dopnormals(l, 2) * posnode[2]) /
              sqrt((dopnormals(l, 0) * dopnormals(l, 0)) + (dopnormals(l, 1) * dopnormals(l, 1)) +
                   (dopnormals(l, 2) * dopnormals(l, 2)));
          if (dcurrent > mslabs(l, 1)) mslabs(l, 1) = dcurrent;
          if (dcurrent < mslabs(l, 0)) mslabs(l, 0) = dcurrent;
        }
      }

      // enlarge slabs with scalar factor
      for (int k = 0; k < kdop / 2; k++)
      {
        mslabs(k, 0) = mslabs(k, 0) - enlarge;
        mslabs(k, 1) = mslabs(k, 1) + enlarge;
      }

      // check if slabs of current master and slave element intercept
      int nintercepts = 0;
      for (int k = 0; k < kdop / 2; k++)
      {
        if ((sslabs(k, 0) <= mslabs(k, 0) && sslabs(k, 1) >= mslabs(k, 0)) ||
            (mslabs(k, 1) >= sslabs(k, 0) && mslabs(k, 0) <= sslabs(k, 0)) ||
            (sslabs(k, 0) <= mslabs(k, 0) && sslabs(k, 1) >= mslabs(k, 1)) ||
            (sslabs(k, 0) >= mslabs(k, 0) && mslabs(k, 1) >= sslabs(k, 1)))
        {
          nintercepts++;
        }
      }

      // std::cout <<"\n"<< Comm().MyPID() << " Number of intercepts found: " << nintercepts ;

      // slabs of current master and slave element do intercept
      if (nintercepts == kdop / 2)
      {
        // std::cout << Comm().MyPID() << " Coupling found between slave element: " << sgid <<" and
        // master element: "<< mgid << std::endl;
        DRT::Element* element = idiscret_->gElement(sgid);
        MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(element);
        selement->AddSearchElements(mgid);
      }
    }  // for all master elements
  }    // for all slave elements

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially coupling sl/ma pairs (public)      popp 10/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::EvaluateSearchBinarytree()
{
  // get out of here if not participating in interface
  if (!lComm()) return true;

  binarytree_->EvaluateSearch();

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlap      popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::MortarCoupling(MORTAR::MortarElement* sele,
    std::vector<MORTAR::MortarElement*> mele,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  PreMortarCoupling(sele, mele, mparams_ptr);

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

    // create Coupling2dManager and evaluate
    MORTAR::Coupling2dManager(Discret(), Dim(), quadratic, InterfaceParams(), sele, mele)
        .EvaluateCoupling(mparams_ptr);
  }
  // ************************************************************** 3D ***
  else if (Dim() == 3)
  {
    // *************************************************** linear 3D ***
    if (!quadratic)
    {
      // create Coupling3dManager and evaluate
      MORTAR::Coupling3dManager(Discret(), Dim(), false, InterfaceParams(), sele, mele)
          .EvaluateCoupling(mparams_ptr);
    }

    // ************************************************** quadratic 3D ***
    else
    {
      // create Coupling3dQuadManager and evaluate
      MORTAR::Coupling3dQuadManager(Discret(), Dim(), false, InterfaceParams(), sele, mele)
          .EvaluateCoupling(mparams_ptr);
    }  // quadratic
  }    // 3D
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  PostMortarCoupling(sele, mele, mparams_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 | Split MortarElements->IntElements for 3D quad. coupling    popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarInterface::SplitIntElements(
    MORTAR::MortarElement& ele, std::vector<Teuchos::RCP<MORTAR::IntElement>>& auxele)
{
  // *********************************************************************
  // do splitting for given element
  // *********************************************************** quad9 ***
  if (ele.Shape() == DRT::Element::quad9)
  {
    // split into for quad4 elements
    int numnode = 4;
    DRT::Element::DiscretizationType dt = DRT::Element::quad4;

    // first integration element
    // containing parent nodes 0,4,8,7
    int nodeids[4] = {0, 0, 0, 0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[8];
    nodeids[3] = ele.NodeIds()[7];

    std::vector<DRT::Node*> nodes(4);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[4];
    nodes[2] = ele.Nodes()[8];
    nodes[3] = ele.Nodes()[7];

    auxele.push_back(Teuchos::rcp(new IntElement(
        0, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));

    // second integration element
    // containing parent nodes 4,1,5,8
    nodeids[0] = ele.NodeIds()[4];
    nodeids[1] = ele.NodeIds()[1];
    nodeids[2] = ele.NodeIds()[5];
    nodeids[3] = ele.NodeIds()[8];

    nodes[0] = ele.Nodes()[4];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[5];
    nodes[3] = ele.Nodes()[8];

    auxele.push_back(Teuchos::rcp(new IntElement(
        1, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));

    // third integration element
    // containing parent nodes 8,5,2,6
    nodeids[0] = ele.NodeIds()[8];
    nodeids[1] = ele.NodeIds()[5];
    nodeids[2] = ele.NodeIds()[2];
    nodeids[3] = ele.NodeIds()[6];

    nodes[0] = ele.Nodes()[8];
    nodes[1] = ele.Nodes()[5];
    nodes[2] = ele.Nodes()[2];
    nodes[3] = ele.Nodes()[6];

    auxele.push_back(Teuchos::rcp(new IntElement(
        2, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));

    // fourth integration element
    // containing parent nodes 7,8,6,3
    nodeids[0] = ele.NodeIds()[7];
    nodeids[1] = ele.NodeIds()[8];
    nodeids[2] = ele.NodeIds()[6];
    nodeids[3] = ele.NodeIds()[3];

    nodes[0] = ele.Nodes()[7];
    nodes[1] = ele.Nodes()[8];
    nodes[2] = ele.Nodes()[6];
    nodes[3] = ele.Nodes()[3];

    auxele.push_back(Teuchos::rcp(new IntElement(
        3, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));
  }

  // *********************************************************** quad8 ***
  else if (ele.Shape() == DRT::Element::quad8)
  {
    // split into four tri3 elements and one quad4 element
    int numnodetri = 3;
    int numnodequad = 4;
    DRT::Element::DiscretizationType dttri = DRT::Element::tri3;
    DRT::Element::DiscretizationType dtquad = DRT::Element::quad4;

    // first integration element
    // containing parent nodes 0,4,7
    int nodeids[3] = {0, 0, 0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[7];

    std::vector<DRT::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[4];
    nodes[2] = ele.Nodes()[7];

    auxele.push_back(Teuchos::rcp(new IntElement(
        0, ele.Id(), ele.Owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.IsSlave(), false)));

    // second integration element
    // containing parent nodes 1,5,4
    nodeids[0] = ele.NodeIds()[1];
    nodeids[1] = ele.NodeIds()[5];
    nodeids[2] = ele.NodeIds()[4];

    nodes[0] = ele.Nodes()[1];
    nodes[1] = ele.Nodes()[5];
    nodes[2] = ele.Nodes()[4];

    auxele.push_back(Teuchos::rcp(new IntElement(
        1, ele.Id(), ele.Owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.IsSlave(), false)));

    // third integration element
    // containing parent nodes 2,6,5
    nodeids[0] = ele.NodeIds()[2];
    nodeids[1] = ele.NodeIds()[6];
    nodeids[2] = ele.NodeIds()[5];

    nodes[0] = ele.Nodes()[2];
    nodes[1] = ele.Nodes()[6];
    nodes[2] = ele.Nodes()[5];

    auxele.push_back(Teuchos::rcp(new IntElement(
        2, ele.Id(), ele.Owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.IsSlave(), false)));

    // fourth integration element
    // containing parent nodes 3,7,6
    nodeids[0] = ele.NodeIds()[3];
    nodeids[1] = ele.NodeIds()[7];
    nodeids[2] = ele.NodeIds()[6];

    nodes[0] = ele.Nodes()[3];
    nodes[1] = ele.Nodes()[7];
    nodes[2] = ele.Nodes()[6];

    auxele.push_back(Teuchos::rcp(new IntElement(
        3, ele.Id(), ele.Owner(), &ele, dttri, numnodetri, nodeids, nodes, ele.IsSlave(), false)));

    // fifth integration element
    // containing parent nodes 4,5,6,7
    int nodeidsquad[4] = {0, 0, 0, 0};
    nodeidsquad[0] = ele.NodeIds()[4];
    nodeidsquad[1] = ele.NodeIds()[5];
    nodeidsquad[2] = ele.NodeIds()[6];
    nodeidsquad[3] = ele.NodeIds()[7];

    std::vector<DRT::Node*> nodesquad(4);
    nodesquad[0] = ele.Nodes()[4];
    nodesquad[1] = ele.Nodes()[5];
    nodesquad[2] = ele.Nodes()[6];
    nodesquad[3] = ele.Nodes()[7];

    auxele.push_back(Teuchos::rcp(new IntElement(4, ele.Id(), ele.Owner(), &ele, dtquad,
        numnodequad, nodeidsquad, nodesquad, ele.IsSlave(), false)));
  }

  // ************************************************************ tri6 ***
  else if (ele.Shape() == DRT::Element::tri6)
  {
    // split into four tri3 elements
    int numnode = 3;
    DRT::Element::DiscretizationType dt = DRT::Element::tri3;

    // first integration element
    // containing parent nodes 0,3,5
    int nodeids[3] = {0, 0, 0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[3];
    nodeids[2] = ele.NodeIds()[5];

    std::vector<DRT::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[3];
    nodes[2] = ele.Nodes()[5];

    auxele.push_back(Teuchos::rcp(new IntElement(
        0, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));

    // second integration element
    // containing parent nodes 3,1,4
    nodeids[0] = ele.NodeIds()[3];
    nodeids[1] = ele.NodeIds()[1];
    nodeids[2] = ele.NodeIds()[4];

    nodes[0] = ele.Nodes()[3];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[4];

    auxele.push_back(Teuchos::rcp(new IntElement(
        1, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));

    // third integration element
    // containing parent nodes 5,4,2
    nodeids[0] = ele.NodeIds()[5];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[2];

    nodes[0] = ele.Nodes()[5];
    nodes[1] = ele.Nodes()[4];
    nodes[2] = ele.Nodes()[2];

    auxele.push_back(Teuchos::rcp(new IntElement(
        2, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));

    // fourth integration element
    // containing parent nodes 4,5,3
    nodeids[0] = ele.NodeIds()[4];
    nodeids[1] = ele.NodeIds()[5];
    nodeids[2] = ele.NodeIds()[3];

    nodes[0] = ele.Nodes()[4];
    nodes[1] = ele.Nodes()[5];
    nodes[2] = ele.Nodes()[3];

    auxele.push_back(Teuchos::rcp(new IntElement(
        3, ele.Id(), ele.Owner(), &ele, dt, numnode, nodeids, nodes, ele.IsSlave(), false)));
  }

  // *********************************************************** quad4 ***
  else if (ele.Shape() == DRT::Element::quad4)
  {
    // 1:1 conversion to IntElement
    std::vector<DRT::Node*> nodes(4);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[2];
    nodes[3] = ele.Nodes()[3];

    auxele.push_back(Teuchos::rcp(new IntElement(0, ele.Id(), ele.Owner(), &ele, ele.Shape(),
        ele.NumNode(), ele.NodeIds(), nodes, ele.IsSlave(), false)));
  }

  // ************************************************************ tri3 ***
  else if (ele.Shape() == DRT::Element::tri3)
  {
    // 1:1 conversion to IntElement
    std::vector<DRT::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[2];

    auxele.push_back(Teuchos::rcp(new IntElement(0, ele.Id(), ele.Owner(), &ele, ele.Shape(),
        ele.NumNode(), ele.NodeIds(), nodes, ele.IsSlave(), false)));
  }

  // ********************************************************* invalid ***
  else
    dserror("ERROR: SplitIntElements called for unknown element shape!");

  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble geometry-dependent lagrange multipliers (global)      popp 05/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleLM(Epetra_Vector& zglobal)
{
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    int dim = mrtrnode->NumDof();
    double* lm = mrtrnode->MoData().lm();

    Epetra_SerialDenseVector lmnode(dim);
    std::vector<int> lmdof(dim);
    std::vector<int> lmowner(dim);

    for (int k = 0; k < dim; ++k)
    {
      lmnode(k) = lm[k];
      lmdof[k] = mrtrnode->Dofs()[k];
      lmowner[k] = mrtrnode->Owner();
    }

    // do assembly
    LINALG::Assemble(zglobal, lmnode, lmdof, lmowner);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble Mortar D matrix                                  popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleD(LINALG::SparseMatrix& dglobal)
{
  const bool nonsmooth = DRT::INPUT::IntegralValue<int>(InterfaceParams(), "NONSMOOTH_GEOMETRIES");
  const bool lagmultlin = (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(
                               InterfaceParams(), "LM_QUAD") == INPAR::MORTAR::lagmult_lin);

  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDM: Node ownership inconsistency!");

    /**************************************************** D-matrix ******/
    if ((mrtrnode->MoData().GetD()).size() > 0)
    {
      const GEN::pairedvector<int, double>& dmap = mrtrnode->MoData().GetD();
      int rowsize = mrtrnode->NumDof();

      GEN::pairedvector<int, double>::const_iterator colcurr;

      for (colcurr = dmap.begin(); colcurr != dmap.end(); ++colcurr)
      {
        double val = colcurr->second;

        DRT::Node* knode = Discret().gNode(colcurr->first);
        if (!knode) dserror("node not found");
        MortarNode* kcnode = dynamic_cast<MortarNode*>(knode);
        if (!kcnode) dserror("node not found");

        for (int j = 0; j < rowsize; ++j)
        {
          int row = mrtrnode->Dofs()[j];
          int col = kcnode->Dofs()[j];

          // do the assembly into global D matrix
          if (!nonsmooth and (shapefcn_ == INPAR::MORTAR::shape_dual or
                                 shapefcn_ == INPAR::MORTAR::shape_petrovgalerkin))
          {
            if (lagmultlin)
            {
              // do lumping of D-matrix
              // create an explicitly diagonal d matrix
              dglobal.Assemble(val, row, row);
            }
            else
            {
              // check for diagonality
              if (row != col && abs(val) > 1.0e-12)
                dserror("ERROR: AssembleDM: D-Matrix is not diagonal!");

              // create an explicitly diagonal d matrix
              if (row == col) dglobal.Assemble(val, row, col);
            }
          }
          else if (nonsmooth or shapefcn_ == INPAR::MORTAR::shape_standard)
          {
            // don't check for diagonality
            // since for standard shape functions, as in general when using
            // arbitrary shape function types, this is not the case

            // create the d matrix, do not assemble zeros
            dglobal.Assemble(val, row, col);
          }
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble Mortar M matrix                                  popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleM(LINALG::SparseMatrix& mglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDM: Node ownership inconsistency!");

    /**************************************************** M-matrix ******/
    if ((mrtrnode->MoData().GetM()).size() > 0)
    {
      const std::map<int, double>& mmap = mrtrnode->MoData().GetM();
      int rowsize = mrtrnode->NumDof();

      std::map<int, double>::const_iterator colcurr;

      for (colcurr = mmap.begin(); colcurr != mmap.end(); ++colcurr)
      {
        DRT::Node* knode = Discret().gNode(colcurr->first);
        if (!knode) dserror("node not found");
        MortarNode* kcnode = dynamic_cast<MortarNode*>(knode);
        if (!kcnode) dserror("node not found");

        double val = colcurr->second;

        for (int j = 0; j < rowsize; ++j)
        {
          int row = mrtrnode->Dofs()[j];
          int col = kcnode->Dofs()[j];

          // do not assemble zeros into m matrix
          //          if (abs(val) > 1.0e-12)
          mglobal.Assemble(val, row, col);
        }
      }
    }

    /************************************************* Mmod-matrix ******/
    if ((mrtrnode->MoData().GetMmod()).size() > 0)
    {
      std::map<int, double>& mmap = mrtrnode->MoData().GetMmod();
      int rowsize = mrtrnode->NumDof();
      int colsize = (int)mmap.size() * rowsize;

      Epetra_SerialDenseMatrix Mnode(rowsize, colsize);
      std::vector<int> lmrow(rowsize);
      std::vector<int> lmcol(colsize);
      std::vector<int> lmrowowner(rowsize);
      std::map<int, double>::const_iterator colcurr;
      int k = 0;

      for (colcurr = mmap.begin(); colcurr != mmap.end(); ++colcurr)
      {
        DRT::Node* knode = Discret().gNode(colcurr->first);
        if (!knode) dserror("node not found");
        MortarNode* kcnode = dynamic_cast<MortarNode*>(knode);
        if (!kcnode) dserror("node not found");

        for (int j = 0; j < rowsize; ++j)
        {
          int row = mrtrnode->Dofs()[j];
          lmrow[j] = row;
          lmrowowner[j] = mrtrnode->Owner();

          int col = kcnode->Dofs()[j];
          double val = colcurr->second;
          lmcol[k] = col;

          Mnode(j, k) = val;
          ++k;
        }
      }

      mglobal.Assemble(-1, Mnode, lmrow, lmrowowner, lmcol);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble Mortar matrices                                 farah 02/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleDM(
    LINALG::SparseMatrix& dglobal, LINALG::SparseMatrix& mglobal)
{
  // call subroutines:

  // assemble mortar matrix D (slave side)
  AssembleD(dglobal);

  // assemble mortar matrix M (master side)
  AssembleM(mglobal);

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrix of normals                                popp 10/11|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleNormals(LINALG::SparseMatrix& nglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDM: Node ownership inconsistency!");

    // nodal normal
    double* nodalnormal = mrtrnode->MoData().n();

    // add normal to corresponding row in global matrix
    for (int k = 0; k < mrtrnode->NumDof(); ++k)
      nglobal.Assemble(nodalnormal[k], gid, mrtrnode->Dofs()[k]);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble interface displacement trafo matrices            popp 06/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::AssembleTrafo(
    LINALG::SparseMatrix& trafo, LINALG::SparseMatrix& invtrafo, std::set<int>& donebefore)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // check for dual shape functions and quadratic slave elements
  if (shapefcn_ == INPAR::MORTAR::shape_standard || quadslave_ == false)
    dserror("ERROR: AssembleTrafo -> you should not be here...");

  // check whether locally linear LM interpolation is used
  const bool lagmultlin = (DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(
                               InterfaceParams(), "LM_QUAD") == INPAR::MORTAR::lagmult_lin);

  //********************************************************************
  //********************************************************************
  // LOOP OVER ALL SLAVE NODES
  //********************************************************************
  //********************************************************************
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleTrafo: Node ownership inconsistency!");

    // find out whether this is a corner node (no trafo), an edge
    // node (trafo of displacement DOFs) or a center node (only for
    // quad9, theoretically trafo of displacement DOFs but not used)
    // and also store transformation factor theta
    enum NodeType
    {
      corner,
      edge,
      center,
      undefined
    };
    NodeType nt = undefined;
    double theta = 0.0;

    // search within the first adjacent element
    MortarElement* mrtrele = dynamic_cast<MortarElement*>(mrtrnode->Elements()[0]);
    MORTAR::MortarElement::DiscretizationType shape = mrtrele->Shape();

    // which discretization type
    switch (shape)
    {
      // line3 contact elements (= tri6||quad8||quad9 discretizations)
      case MORTAR::MortarElement::line3:
      {
        // modification factor
        if (INPAR::MORTAR::LagMultQuad() == INPAR::MORTAR::lagmult_lin)
          theta = 1.0 / 2.0;
        else
          theta = 1.0 / 5.0;

        // corner nodes
        if (mrtrnode->Id() == mrtrele->NodeIds()[0] || mrtrnode->Id() == mrtrele->NodeIds()[1])
        {
          nt = corner;
        }

        // edge nodes
        else if (mrtrnode->Id() == mrtrele->NodeIds()[2])
        {
          nt = edge;
        }

        break;
      }

      // tri6 contact elements (= tet10 discretizations)
      case MORTAR::MortarElement::tri6:
      {
        // modification factor
        theta = 1.0 / 5.0;

        // corner nodes
        if (mrtrnode->Id() == mrtrele->NodeIds()[0] || mrtrnode->Id() == mrtrele->NodeIds()[1] ||
            mrtrnode->Id() == mrtrele->NodeIds()[2])
        {
          nt = corner;
        }

        // edge nodes
        else if (mrtrnode->Id() == mrtrele->NodeIds()[3] ||
                 mrtrnode->Id() == mrtrele->NodeIds()[4] || mrtrnode->Id() == mrtrele->NodeIds()[5])
        {
          nt = edge;
        }

        break;
      }

        // quad8 contact elements (= hex20 discretizations)
      case MORTAR::MortarElement::quad8:
      {
        // modification factor
        theta = 1.0 / 5.0;

        // corner nodes
        if (mrtrnode->Id() == mrtrele->NodeIds()[0] || mrtrnode->Id() == mrtrele->NodeIds()[1] ||
            mrtrnode->Id() == mrtrele->NodeIds()[2] || mrtrnode->Id() == mrtrele->NodeIds()[3])
        {
          nt = corner;
        }

        // edge nodes
        else if (mrtrnode->Id() == mrtrele->NodeIds()[4] ||
                 mrtrnode->Id() == mrtrele->NodeIds()[5] ||
                 mrtrnode->Id() == mrtrele->NodeIds()[6] || mrtrnode->Id() == mrtrele->NodeIds()[7])
        {
          nt = edge;
        }

        break;
      }

        // quad9 contact elements (= hex27 discretizations)
        // *************************************************
        // ** currently we only use this modification for **
        // ** tri6 and quad8 surfaces, but NOT for quad9  **
        // ** as in this case, there is no real need!     **
        // ** (positivity of shape function integrals)    **
        // ** thus, we simply want to assemble the        **
        // ** identity matrix here, which we achieve by   **
        // ** setting the trafo factor theta = 0.0!       **
        // *************************************************
      case MORTAR::MortarElement::quad9:
      {
        // modification factor
        theta = 0.0;

        // corner nodes
        if (mrtrnode->Id() == mrtrele->NodeIds()[0] || mrtrnode->Id() == mrtrele->NodeIds()[1] ||
            mrtrnode->Id() == mrtrele->NodeIds()[2] || mrtrnode->Id() == mrtrele->NodeIds()[3])
        {
          nt = corner;
        }

        // edge nodes
        else if (mrtrnode->Id() == mrtrele->NodeIds()[4] ||
                 mrtrnode->Id() == mrtrele->NodeIds()[5] ||
                 mrtrnode->Id() == mrtrele->NodeIds()[6] || mrtrnode->Id() == mrtrele->NodeIds()[7])
        {
          nt = edge;
        }

        // center node
        else if (mrtrnode->Id() == mrtrele->NodeIds()[8])
        {
          nt = center;
        }

        break;
      }

        // other cases
      default:
      {
        dserror("ERROR: Trafo matrix only for line3/tri6/quad8/quad9 contact elements");
        break;
      }
    }  // switch(Shape)

    //********************************************************************
    // CASE 1: CORNER NODES AND CENTER NODE
    //********************************************************************
    if (nt == corner || nt == center)
    {
      // check if processed before
      std::set<int>::iterator iter = donebefore.find(gid);

      // if not then assemble trafo matrix block
      if (iter == donebefore.end())
      {
        // add to set of processed nodes
        donebefore.insert(gid);

        // add transformation matrix block (unity block!)
        for (int k = 0; k < mrtrnode->NumDof(); ++k)
        {
          // assemble diagonal values
          trafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
          invtrafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
        }
      }
    }

    //********************************************************************
    // CASE 2: EDGE NODES
    //********************************************************************
    else if (nt == edge)
    {
      // check if processed before
      std::set<int>::iterator iter = donebefore.find(gid);

      // if not then assemble trafo matrix block
      if (iter == donebefore.end())
      {
        // add to set of processed nodes
        donebefore.insert(gid);

        // find adjacent corner nodes locally
        int index1 = 0;
        int index2 = 0;
        int hoindex = mrtrele->GetLocalNodeId(gid);
        DRT::UTILS::getCornerNodeIndices(index1, index2, hoindex, shape);

        // find adjacent corner nodes globally
        int gindex1 = mrtrele->NodeIds()[index1];
        int gindex2 = mrtrele->NodeIds()[index2];
        // std::cout << "-> adjacent corner nodes: " << gindex1 << " " << gindex2 << std::endl;
        DRT::Node* adjnode1 = idiscret_->gNode(gindex1);
        if (!adjnode1) dserror("ERROR: Cannot find node with gid %", gindex1);
        MortarNode* adjmrtrnode1 = dynamic_cast<MortarNode*>(adjnode1);
        DRT::Node* adjnode2 = idiscret_->gNode(gindex2);
        if (!adjnode2) dserror("ERROR: Cannot find node with gid %", gindex2);
        MortarNode* adjmrtrnode2 = dynamic_cast<MortarNode*>(adjnode2);

        // add transformation matrix block
        for (int k = 0; k < mrtrnode->NumDof(); ++k)
        {
          // assemble diagonal values
          trafo.Assemble(1.0 - 2 * theta, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
          invtrafo.Assemble(1.0 / (1.0 - 2 * theta), mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);

          // assemble off-diagonal values
          trafo.Assemble(theta, mrtrnode->Dofs()[k], adjmrtrnode1->Dofs()[k]);
          trafo.Assemble(theta, mrtrnode->Dofs()[k], adjmrtrnode2->Dofs()[k]);
          invtrafo.Assemble(
              -theta / (1.0 - 2 * theta), mrtrnode->Dofs()[k], adjmrtrnode1->Dofs()[k]);
          invtrafo.Assemble(
              -theta / (1.0 - 2 * theta), mrtrnode->Dofs()[k], adjmrtrnode2->Dofs()[k]);
        }
      }
    }

    //********************************************************************
    // CASE 3: UNDEFINED NODES
    //********************************************************************
    else
    {
      dserror("ERROR: Undefined node type (corner, edge, center)");
    }
  }

  // assembly for locally linear LM interpolation
  if (lagmultlin)
  {
    //********************************************************************
    //********************************************************************
    // LOOP OVER ALL MASTER NODES
    //********************************************************************
    //********************************************************************
    // loop over proc's master nodes of the interface for assembly
    // use standard row map to assemble each node only once
    for (int i = 0; i < mnoderowmap_->NumMyElements(); ++i)
    {
      int gid = mnoderowmap_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %", gid);
      MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

      if (mrtrnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleTrafo: Node ownership inconsistency!");

      // find out whether this is a "real" master node (no trafo), a former
      // slave edge node (trafo of displacement DOFs) or a former slave
      // center node (only for quadd9, trafo of displacement DOFs)
      enum NodeType
      {
        master,
        slaveedge,
        slavecenter,
        undefined
      };
      NodeType nt = undefined;

      // search within the first adjacent element
      MortarElement* mrtrele = dynamic_cast<MortarElement*>(mrtrnode->Elements()[0]);
      MORTAR::MortarElement::DiscretizationType shape = mrtrele->Shape();

      // real master nodes are easily identified
      if (!mrtrnode->IsOnBound()) nt = master;
      // former slave node type depends on discretization
      else
      {
        switch (shape)
        {
          case MORTAR::MortarElement::line3:
          {
            // edge node
            if (mrtrnode->Id() == mrtrele->NodeIds()[2]) nt = slaveedge;

            break;
          }

          // tri6 contact elements (= tet10 discretizations)
          case MORTAR::MortarElement::tri6:
          {
            // edge nodes
            if (mrtrnode->Id() == mrtrele->NodeIds()[3] ||
                mrtrnode->Id() == mrtrele->NodeIds()[4] || mrtrnode->Id() == mrtrele->NodeIds()[5])
              nt = slaveedge;

            break;
          }

          // quad8 contact elements (= hex20 discretizations)
          case MORTAR::MortarElement::quad8:
          {
            // edge nodes
            if (mrtrnode->Id() == mrtrele->NodeIds()[4] ||
                mrtrnode->Id() == mrtrele->NodeIds()[5] ||
                mrtrnode->Id() == mrtrele->NodeIds()[6] || mrtrnode->Id() == mrtrele->NodeIds()[7])
              nt = slaveedge;

            break;
          }

          // quad9 contact elements (= hex27 discretizations)
          case MORTAR::MortarElement::quad9:
          {
            // edge nodes
            if (mrtrnode->Id() == mrtrele->NodeIds()[4] ||
                mrtrnode->Id() == mrtrele->NodeIds()[5] ||
                mrtrnode->Id() == mrtrele->NodeIds()[6] || mrtrnode->Id() == mrtrele->NodeIds()[7])
              nt = slaveedge;
            // center node
            else if (mrtrnode->Id() == mrtrele->NodeIds()[8])
              nt = slavecenter;

            break;
          }

          // other cases
          default:
          {
            dserror("ERROR: Trafo matrix only for line3/tri6/quad8/quad9 contact elements");
            break;
          }
        }  // switch(Shape)
      }

      //********************************************************************
      // CASE 1: REAL MASTER NODE
      //********************************************************************
      if (nt == master)
      {
        // check if processed before
        std::set<int>::iterator iter = donebefore.find(gid);

        // if not then assemble trafo matrix block
        if (iter == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(gid);

          // add transformation matrix block (unity block!)
          for (int k = 0; k < mrtrnode->NumDof(); ++k)
          {
            // assemble diagonal values
            trafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
            invtrafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
          }
        }
      }

      //********************************************************************
      // CASE 2: FORMER SLAVE EDGE NODE
      // (for linear LM interpolation -> full distribution of edge nodes)
      // (nevertheless, we keep the 1.0 on the main diagonal -> no PoU!)
      //********************************************************************
      else if (nt == slaveedge)
      {
        // check if processed before
        std::set<int>::iterator iter = donebefore.find(gid);

        // if not then assemble trafo matrix block
        if (iter == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(gid);

          // find adjacent corner nodes locally
          int index1 = 0;
          int index2 = 0;
          int hoindex = mrtrele->GetLocalNodeId(gid);
          DRT::UTILS::getCornerNodeIndices(index1, index2, hoindex, shape);

          // find adjacent corner nodes globally
          int gindex1 = mrtrele->NodeIds()[index1];
          int gindex2 = mrtrele->NodeIds()[index2];
          // std::cout << "-> adjacent corner nodes: " << gindex1 << " " << gindex2 << std::endl;
          DRT::Node* adjnode1 = idiscret_->gNode(gindex1);
          if (!adjnode1) dserror("ERROR: Cannot find node with gid %", gindex1);
          MortarNode* adjmrtrnode1 = dynamic_cast<MortarNode*>(adjnode1);
          DRT::Node* adjnode2 = idiscret_->gNode(gindex2);
          if (!adjnode2) dserror("ERROR: Cannot find node with gid %", gindex2);
          MortarNode* adjmrtrnode2 = dynamic_cast<MortarNode*>(adjnode2);

          // add transformation matrix block
          for (int k = 0; k < mrtrnode->NumDof(); ++k)
          {
            // assemble diagonal values
            trafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
            invtrafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);

            // assemble off-diagonal values
            trafo.Assemble(0.5, mrtrnode->Dofs()[k], adjmrtrnode1->Dofs()[k]);
            trafo.Assemble(0.5, mrtrnode->Dofs()[k], adjmrtrnode2->Dofs()[k]);
            invtrafo.Assemble(-0.5, mrtrnode->Dofs()[k], adjmrtrnode1->Dofs()[k]);
            invtrafo.Assemble(-0.5, mrtrnode->Dofs()[k], adjmrtrnode2->Dofs()[k]);
          }
        }
      }

      //********************************************************************
      // CASE 3: FORMER SLAVE CENTER NODE (QUAD9)
      // (for linear LM interpolation -> full distribution of corner nodes)
      // (nevertheless, we keep the 1.0 on the main diagonal -> no PoU!)
      //********************************************************************
      else if (nt == slavecenter)
      {
        // check if processed before
        std::set<int>::iterator iter = donebefore.find(gid);

        // if not then assemble trafo matrix block
        if (iter == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(gid);

          // find adjacent corner nodes globally
          int gindex1 = mrtrele->NodeIds()[0];
          int gindex2 = mrtrele->NodeIds()[1];
          int gindex3 = mrtrele->NodeIds()[2];
          int gindex4 = mrtrele->NodeIds()[3];
          // std::cout << "-> adjacent corner nodes: " << gindex1 << " " << gindex2 << std::endl;
          // std::cout << "-> adjacent corner nodes: " << gindex3 << " " << gindex4 << std::endl;
          DRT::Node* adjnode1 = idiscret_->gNode(gindex1);
          if (!adjnode1) dserror("ERROR: Cannot find node with gid %", gindex1);
          MortarNode* adjmrtrnode1 = dynamic_cast<MortarNode*>(adjnode1);
          DRT::Node* adjnode2 = idiscret_->gNode(gindex2);
          if (!adjnode2) dserror("ERROR: Cannot find node with gid %", gindex2);
          MortarNode* adjmrtrnode2 = dynamic_cast<MortarNode*>(adjnode2);
          DRT::Node* adjnode3 = idiscret_->gNode(gindex3);
          if (!adjnode3) dserror("ERROR: Cannot find node with gid %", gindex3);
          MortarNode* adjmrtrnode3 = dynamic_cast<MortarNode*>(adjnode3);
          DRT::Node* adjnode4 = idiscret_->gNode(gindex4);
          if (!adjnode4) dserror("ERROR: Cannot find node with gid %", gindex4);
          MortarNode* adjmrtrnode4 = dynamic_cast<MortarNode*>(adjnode4);

          // add transformation matrix block
          for (int k = 0; k < mrtrnode->NumDof(); ++k)
          {
            // assemble diagonal values
            trafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);
            invtrafo.Assemble(1.0, mrtrnode->Dofs()[k], mrtrnode->Dofs()[k]);

            // assemble off-diagonal values
            trafo.Assemble(0.25, mrtrnode->Dofs()[k], adjmrtrnode1->Dofs()[k]);
            trafo.Assemble(0.25, mrtrnode->Dofs()[k], adjmrtrnode2->Dofs()[k]);
            trafo.Assemble(0.25, mrtrnode->Dofs()[k], adjmrtrnode3->Dofs()[k]);
            trafo.Assemble(0.25, mrtrnode->Dofs()[k], adjmrtrnode4->Dofs()[k]);
            invtrafo.Assemble(-0.25, mrtrnode->Dofs()[k], adjmrtrnode1->Dofs()[k]);
            invtrafo.Assemble(-0.25, mrtrnode->Dofs()[k], adjmrtrnode2->Dofs()[k]);
            invtrafo.Assemble(-0.25, mrtrnode->Dofs()[k], adjmrtrnode3->Dofs()[k]);
            invtrafo.Assemble(-0.25, mrtrnode->Dofs()[k], adjmrtrnode4->Dofs()[k]);
          }
        }
      }

      //********************************************************************
      // CASE 4: UNDEFINED NODES
      //********************************************************************
      else
        dserror("ERROR: Undefined node type (corner, edge, center)");
    }
  }  // end of assembly for locally linear LM interpolation

  return;
}

/*----------------------------------------------------------------------*
 |  Detect actual meshtying zone (node by node)               popp 08/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::DetectTiedSlaveNodes(int& founduntied)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  //**********************************************************************
  // STEP 1: Build tying info for slave node row map (locally+globally)
  //**********************************************************************
  // global vector for tying info
  Teuchos::RCP<Epetra_Vector> rowtied = Teuchos::rcp(new Epetra_Vector(*snoderowmap_));

  // loop over proc's slave row nodes of the interface for detection
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    // perform detection
    const GEN::pairedvector<int, double>& dmap = mrtrnode->MoData().GetD();
    const std::map<int, double>& mmap = mrtrnode->MoData().GetM();
    int sized = dmap.size();
    int sizem = mmap.size();

    // found untied node
    if (sized == 0 && sizem == 0)
    {
      // increase counter
      founduntied += 1;

      // set node status to untied slave
      mrtrnode->SetTiedSlave() = false;

      // set vector entry (tiedtoggle)
      (*rowtied)[i] = 1.0;
    }

    // found tied node
    else if (sized > 0 && sizem > 0)
    {
      // do nothing
    }

    // found inconsistency
    else
    {
      dserror("ERROR: Inconsistency in tied/untied node detection");
    }
  }

  //**********************************************************************
  // STEP 2: Export tying info to slave node column map (globally)
  //**********************************************************************
  // export tying information to standard column map
  Teuchos::RCP<Epetra_Vector> coltied = Teuchos::rcp(new Epetra_Vector(*snodecolmap_));
  LINALG::Export(*rowtied, *coltied);

  //**********************************************************************
  // STEP 3: Extract tying info for slave node column map (locally)
  //**********************************************************************
  // loop over proc's slave col nodes of the interface for storage
  for (int i = 0; i < snodecolmap_->NumMyElements(); ++i)
  {
    int gid = snodecolmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    MortarNode* mrtrnode = dynamic_cast<MortarNode*>(node);

    // check if this node is untied
    if ((*coltied)[i] == 1.0) mrtrnode->SetTiedSlave() = false;
  }

  return;
}

/*----------------------------------------------------------------------*
 | create volume ghosting (public)                            ager 06/15|
 *----------------------------------------------------------------------*/
void MORTAR::MortarInterface::CreateVolumeGhosting()
{
  INPAR::CONTACT::Problemtype prb = (INPAR::CONTACT::Problemtype)InterfaceParams().get<int>(
      "PROBTYPE", (int)INPAR::CONTACT::other);

  switch (prb)
  {
    case INPAR::CONTACT::tsi:
    {
      std::vector<std::string> tar_dis;
      tar_dis.push_back("structure");
      tar_dis.push_back("thermo");
      std::vector<std::pair<int, int>> material_map;
      material_map.push_back(std::pair<int, int>(0, 1));
      material_map.push_back(std::pair<int, int>(1, 0));

      MORTAR::UTILS::CreateVolumeGhosting(Discret(), tar_dis, material_map);
      break;
    }
    default:
    {
      std::vector<std::string> tar_dis;
      tar_dis.push_back("structure");
      MORTAR::UTILS::CreateVolumeGhosting(Discret(), tar_dis, std::vector<std::pair<int, int>>(0));

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool MORTAR::MortarInterface::HasMaSharingRefInterface() const
{
  return (interfaceData_->GetMaSharingRefInterfacePtr() != NULL);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const MORTAR::MortarInterface* MORTAR::MortarInterface::GetMaSharingRefInterfacePtr() const
{
  return interfaceData_->GetMaSharingRefInterfacePtr();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MortarInterface::AddMaSharingRefInterface(const MortarInterface* ref_interface)
{
  // avoid non-uniqueness and closed loops
  if (ref_interface->HasMaSharingRefInterface())
  {
    if (ref_interface->GetMaSharingRefInterfacePtr() == this) return;
  }

  /* The following test is valid, since this interface must be a FULL subset of
   * the reference interface, i.e. all master elements of this interface must
   * be also contained in the reference master element map. Therefore, if a new
   * reference interface candidate is supposed to replace the current one, it
   * must have more global entries in its master element map. Otherwise, it is
   * as well a sub-set of the current reference interface.
   *
   * Again: The last assumption holds only if no partial overlaps are allowed.
   *                                                          hiermeier 01/18 */
  if (HasMaSharingRefInterface())
  {
    const int size_curr_ref_interface =
        GetMaSharingRefInterfacePtr()->MasterRowElements()->NumGlobalElements();
    const int size_new_ref_interface = ref_interface->MasterRowElements()->NumGlobalElements();

    if (size_curr_ref_interface >= size_new_ref_interface) return;
  }

  interfaceData_->SetMaSharingRefInterfacePtr(ref_interface);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MORTAR::MortarInterface::PostprocessQuantities(const Teuchos::ParameterList& outputParams)
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
   * Grab it from there, then use LINALG::Export() to extract the interface
   * portion from the global displacement vector. Finally, write the interface
   * portion using this interfaces' discretization writer.
   */
  {
    // Get full displacement vector and extract interface displacement
    RCP<const Epetra_Vector> disp = outputParams.get<RCP<const Epetra_Vector>>("displacement");
    RCP<Epetra_Vector> iDisp = LINALG::CreateVector(*idiscret_->DofRowMap());
    LINALG::Export(*disp, *iDisp);

    // Write the interface displacement field
    writer->WriteVector("displacement", iDisp, IO::VectorType::dofvector);
  }

  // Write Lagrange multiplier field
  {
    // Get full Lagrange multiplier vector and extract values of this interface
    RCP<const Epetra_Vector> lagMult =
        outputParams.get<RCP<const Epetra_Vector>>("interface traction");
    RCP<Epetra_Vector> iLagMult = LINALG::CreateVector(*idiscret_->DofRowMap());
    LINALG::Export(*lagMult, *iLagMult);

    // Write this interface's Lagrange multiplier field
    writer->WriteVector("interfacetraction", iLagMult, IO::VectorType::dofvector);
  }

  // Write nodal forces of slave side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> slaveforces =
        outputParams.get<RCP<const Epetra_Vector>>("slave forces");
    RCP<Epetra_Vector> forces = LINALG::CreateVector(*idiscret_->DofRowMap());
    LINALG::Export(*slaveforces, *forces);

    // Write to output
    writer->WriteVector("slaveforces", forces, IO::VectorType::dofvector);
  }

  // Write nodal forces of master side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> masterforces =
        outputParams.get<RCP<const Epetra_Vector>>("master forces");
    RCP<Epetra_Vector> forces = LINALG::CreateVector(*idiscret_->DofRowMap());
    LINALG::Export(*masterforces, *forces);

    // Write to output
    writer->WriteVector("masterforces", forces, IO::VectorType::dofvector);
  }


  // Nodes: node-based vector with '0' at slave nodes and '1' at master nodes
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*mnoderowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> nodeRowMap = LINALG::MergeMap(snoderowmap_, mnoderowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = LINALG::CreateVector(*nodeRowMap, true);
    LINALG::Export(*masterVec, *masterSlaveVec);

    writer->WriteVector("slavemasternodes", masterSlaveVec, IO::VectorType::nodevector);
  }

  // Elements: element-based vector with '0' at slave elements and '1' at master elements
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*melerowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> eleRowMap = LINALG::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = LINALG::CreateVector(*eleRowMap, true);
    LINALG::Export(*masterVec, *masterSlaveVec);

    writer->WriteVector("slavemasterelements", masterSlaveVec, IO::VectorType::elementvector);
  }

  // Write element owners
  {
    RCP<const Epetra_Map> eleRowMap = LINALG::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> owner = LINALG::CreateVector(*eleRowMap);

    for (int i = 0; i < idiscret_->ElementRowMap()->NumMyElements(); ++i)
      (*owner)[i] = idiscret_->lRowElement(i)->Owner();

    writer->WriteVector("Owner", owner, IO::VectorType::elementvector);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool MORTAR::MortarInterface::CheckOutputList(
    const Teuchos::ParameterList& outParams, const std::vector<std::string>& requiredEntries) const
{
  // Check for each required parameter entry if it exists
  for (std::vector<std::string>::const_iterator it = requiredEntries.begin();
       it < requiredEntries.end(); ++it)
  {
    if (not outParams.isParameter(*it))
    {
      dserror("Parameter list is missing the required entry '%s'.", (*it).c_str());
      return false;
    }
  }

  // We only make it to here, if all checks passed. So it's safe to return 'true'.
  return true;
}
