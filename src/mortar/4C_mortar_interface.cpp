/*----------------------------------------------------------------------*/
/*! \file
\brief One mortar coupling interface

\level 1

*/
/*-----------------------------------------------------------------------*/

#include "4C_mortar_interface.hpp"

#include "4C_binstrategy.hpp"
#include "4C_contact_interpolator.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_binarytree.hpp"
#include "4C_mortar_coupling2d.hpp"
#include "4C_mortar_coupling3d.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_dofset.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_integrator.hpp"
#include "4C_mortar_interface_utils.hpp"
#include "4C_mortar_node.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_rebalance_graph_based.hpp"

#include <Epetra_Map.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Vector.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Mortar::InterfaceDataContainer::InterfaceDataContainer()
    : id_(-1),
      comm_(nullptr),
      redistributed_(false),
      idiscret_(Teuchos::null),
      dim_(-1),
      imortar_(Teuchos::ParameterList()),
      shapefcn_(Inpar::Mortar::shape_undefined),
      quadslave_(false),
      extendghosting_(Inpar::Mortar::ExtendGhosting::redundant_master),
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
      searchalgo_(Inpar::Mortar::search_binarytree),
      binarytree_(Teuchos::null),
      searchparam_(-1.0),
      searchuseauxpos_(false),
      inttime_interface_(0.0),
      nurbs_(false),
      poro_(false),
      porotype_(Inpar::Mortar::other),
      ehl_(false),
      isinit_(false)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Mortar::Interface::Interface(Teuchos::RCP<Mortar::InterfaceDataContainer> interfaceData)
    : interface_data_(std::move(interfaceData)),
      id_(interface_data_->Id()),
      comm_(interface_data_->CommPtr()),
      procmap_(interface_data_->ProcMap()),
      redistributed_(interface_data_->IsRedistributed()),
      idiscret_(interface_data_->i_discret()),
      dim_(interface_data_->Dim()),
      imortar_(interface_data_->IMortar()),
      shapefcn_(interface_data_->ShapeFcn()),
      quadslave_(interface_data_->IsQuadSlave()),
      oldnodecolmap_(interface_data_->OldNodeColMap()),
      oldelecolmap_(interface_data_->OldEleColMap()),
      snoderowmap_(interface_data_->SNodeRowMap()),
      snodecolmap_(interface_data_->SNodeColMap()),
      mnoderowmap_(interface_data_->MNodeRowMap()),
      mnodecolmap_(interface_data_->MNodeColMap()),
      snoderowmapbound_(interface_data_->SNodeRowMapBound()),
      snodecolmapbound_(interface_data_->SNodeColMapBound()),
      mnoderowmapnobound_(interface_data_->MNodeRowMapNoBound()),
      mnodecolmapnobound_(interface_data_->MNodeColMapNoBound()),
      selerowmap_(interface_data_->SEleRowMap()),
      selecolmap_(interface_data_->SEleColMap()),
      melerowmap_(interface_data_->MEleRowMap()),
      melecolmap_(interface_data_->MEleColMap()),
      sdofrowmap_(interface_data_->SDofRowMap()),
      sdofcolmap_(interface_data_->SDofColMap()),
      mdofrowmap_(interface_data_->MDofRowMap()),
      mdofcolmap_(interface_data_->MDofColMap()),
      psdofrowmap_(interface_data_->PSDofRowMap()),
      plmdofmap_(interface_data_->PLmDofRowMap()),
      lmdofmap_(interface_data_->LmDofRowMap()),
      maxdofglobal_(interface_data_->MaxDofGlobal()),
      searchalgo_(interface_data_->SearchAlgorithm()),
      binarytree_(interface_data_->BinaryTree()),
      searchparam_(interface_data_->SearchParam()),
      searchuseauxpos_(interface_data_->SearchUseAuxPos()),
      inttime_interface_(interface_data_->IntTimeInterface()),
      nurbs_(interface_data_->IsNurbs()),
      ehl_(interface_data_->IsEhl())
{
  if (not interface_data_->is_init())
  {
    FOUR_C_THROW(
        "This constructor is only allowed for already initialized "
        "interface data containers!");
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Mortar::Interface> Mortar::Interface::Create(const int id, const Epetra_Comm& comm,
    const int spatialDim, const Teuchos::ParameterList& imortar)
{
  Teuchos::RCP<Mortar::InterfaceDataContainer> interfaceData =
      Teuchos::rcp(new Mortar::InterfaceDataContainer());

  return Teuchos::rcp(new Mortar::Interface(interfaceData, id, comm, spatialDim, imortar));
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
Mortar::Interface::Interface(Teuchos::RCP<InterfaceDataContainer> interfaceData, const int id,
    const Epetra_Comm& comm, const int spatialDim, const Teuchos::ParameterList& imortar)
    : interface_data_(std::move(interfaceData)),
      id_(interface_data_->Id()),
      comm_(interface_data_->CommPtr()),
      procmap_(interface_data_->ProcMap()),
      redistributed_(interface_data_->IsRedistributed()),
      idiscret_(interface_data_->i_discret()),
      dim_(interface_data_->Dim()),
      imortar_(interface_data_->IMortar()),
      shapefcn_(interface_data_->ShapeFcn()),
      quadslave_(interface_data_->IsQuadSlave()),
      oldnodecolmap_(interface_data_->OldNodeColMap()),
      oldelecolmap_(interface_data_->OldEleColMap()),
      snoderowmap_(interface_data_->SNodeRowMap()),
      snodecolmap_(interface_data_->SNodeColMap()),
      mnoderowmap_(interface_data_->MNodeRowMap()),
      mnodecolmap_(interface_data_->MNodeColMap()),
      snoderowmapbound_(interface_data_->SNodeRowMapBound()),
      snodecolmapbound_(interface_data_->SNodeColMapBound()),
      mnoderowmapnobound_(interface_data_->MNodeRowMapNoBound()),
      mnodecolmapnobound_(interface_data_->MNodeColMapNoBound()),
      selerowmap_(interface_data_->SEleRowMap()),
      selecolmap_(interface_data_->SEleColMap()),
      melerowmap_(interface_data_->MEleRowMap()),
      melecolmap_(interface_data_->MEleColMap()),
      sdofrowmap_(interface_data_->SDofRowMap()),
      sdofcolmap_(interface_data_->SDofColMap()),
      mdofrowmap_(interface_data_->MDofRowMap()),
      mdofcolmap_(interface_data_->MDofColMap()),
      psdofrowmap_(interface_data_->PSDofRowMap()),
      plmdofmap_(interface_data_->PLmDofRowMap()),
      lmdofmap_(interface_data_->LmDofRowMap()),
      maxdofglobal_(interface_data_->MaxDofGlobal()),
      searchalgo_(interface_data_->SearchAlgorithm()),
      binarytree_(interface_data_->BinaryTree()),
      searchparam_(interface_data_->SearchParam()),
      searchuseauxpos_(interface_data_->SearchUseAuxPos()),
      inttime_interface_(interface_data_->IntTimeInterface()),
      nurbs_(interface_data_->IsNurbs()),
      ehl_(interface_data_->IsEhl())
{
  interface_data_->set_is_init(true);
  id_ = id;
  comm_ = Teuchos::rcpFromRef(comm);
  dim_ = spatialDim;
  imortar_.setParameters(imortar);
  quadslave_ = false;
  interface_data_->SetExtendGhosting(Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
      imortar.sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY"));
  searchalgo_ =
      Core::UTILS::IntegralValue<Inpar::Mortar::SearchAlgorithm>(imortar, "SEARCH_ALGORITHM");
  searchparam_ = imortar.get<double>("SEARCH_PARAM");
  searchuseauxpos_ = Core::UTILS::IntegralValue<int>(imortar, "SEARCH_USE_AUX_POS");
  nurbs_ = imortar.get<bool>("NURBS");

  if (Dim() != 2 && Dim() != 3) FOUR_C_THROW("Mortar problem must be 2D or 3D.");

  procmap_.clear();

  create_interface_discretization();
  set_shape_function_type();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::create_interface_discretization()
{
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());

  // Create name for mortar interface discretization
  std::stringstream dis_name;
  dis_name << "mortar_interface_" << id_;

  // Create the required type of discretization
  if (nurbs_)
  {
    idiscret_ = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
        dis_name.str(), comm, Global::Problem::Instance()->NDim()));

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
    idiscret_ = Teuchos::rcp(
        new Core::FE::Discretization(dis_name.str(), comm, Global::Problem::Instance()->NDim()));
  }

  // Prepare discretization writer
  idiscret_->SetWriter(Teuchos::rcp(new Core::IO::DiscretizationWriter(idiscret_,
      Global::Problem::Instance()->OutputControlFile(),
      Global::Problem::Instance()->spatial_approximation_type())));
  FOUR_C_ASSERT(not idiscret_->Writer().is_null(), "Setup of discretization writer failed.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::set_shape_function_type()
{
  auto shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(interface_params(), "LM_SHAPEFCN");
  switch (shapefcn)
  {
    case Inpar::Mortar::shape_dual:
    {
      shapefcn_ = Inpar::Mortar::shape_dual;
      break;
    }
    case Inpar::Mortar::shape_petrovgalerkin:
    {
      shapefcn_ = Inpar::Mortar::shape_petrovgalerkin;
      break;
    }
    case Inpar::Mortar::shape_standard:
    {
      shapefcn_ = Inpar::Mortar::shape_standard;
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Invalid shape function type. Interface must either have dual or standard shape "
          "functions.");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Mortar::Interface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    os << "\nMortar Interface Id " << id_ << std::endl;
    os << "Mortar Interface discretization:" << std::endl;
  }
  os << Discret();
}

/*----------------------------------------------------------------------*
 |  check if interface is fill_complete (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
bool Mortar::Interface::Filled() const { return idiscret_->Filled(); }

/*----------------------------------------------------------------------*
 |  print parallel distribution (public)                      popp 06/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::print_parallel_distribution() const
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
    if (interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::redundant_master)
    {
      my_m_ghostnodes[myrank] = mnoderowmap_->NumGlobalElements() - my_m_nodes[myrank];
      my_m_ghostele[myrank] = melerowmap_->NumGlobalElements() - my_m_elements[myrank];
    }
    else if (interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::redundant_all)
    {
      my_m_ghostnodes[myrank] = mnoderowmap_->NumGlobalElements() - my_m_nodes[myrank];
      my_m_ghostele[myrank] = melerowmap_->NumGlobalElements() - my_m_elements[myrank];
      my_s_ghostnodes[myrank] = snoderowmap_->NumGlobalElements() - my_s_nodes[myrank];
      my_s_ghostele[myrank] = selerowmap_->NumGlobalElements() - my_s_elements[myrank];
    }

    Discret().Comm().SumAll(my_n_nodes.data(), n_nodes.data(), numproc);
    Discret().Comm().SumAll(my_n_ghostnodes.data(), n_ghostnodes.data(), numproc);
    Discret().Comm().SumAll(my_n_elements.data(), n_elements.data(), numproc);
    Discret().Comm().SumAll(my_n_ghostele.data(), n_ghostele.data(), numproc);

    Discret().Comm().SumAll(my_s_nodes.data(), s_nodes.data(), numproc);
    Discret().Comm().SumAll(my_s_ghostnodes.data(), s_ghostnodes.data(), numproc);
    Discret().Comm().SumAll(my_s_elements.data(), s_elements.data(), numproc);
    Discret().Comm().SumAll(my_s_ghostele.data(), s_ghostele.data(), numproc);

    Discret().Comm().SumAll(my_m_nodes.data(), m_nodes.data(), numproc);
    Discret().Comm().SumAll(my_m_ghostnodes.data(), m_ghostnodes.data(), numproc);
    Discret().Comm().SumAll(my_m_elements.data(), m_elements.data(), numproc);
    Discret().Comm().SumAll(my_m_ghostele.data(), m_ghostele.data(), numproc);

    if (myrank == 0)
    {
      std::cout << std::endl;
      std::cout << "  discretization: " << Discret().Name() << std::endl;

      // Compute and print statistics
      {
        std::cout << "\n"
                  << "    Statistics of parallel distribution across " << numproc
                  << " ranks:" << std::endl;
        printf(
            "    "
            "+----------------------+----------------+----------------+-----------------+----------"
            "--------+\n");
        printf(
            "    | Type                 | min over procs | max over procs | mean over procs | "
            "max-to-min ratio |\n");
        printf(
            "    "
            "+----------------------+----------------+----------------+-----------------+----------"
            "--------+\n");
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "nodes (s+m)", n_nodes, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "ghost nodes (s+m)", n_ghostnodes, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "elements (s+m)", n_elements, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "ghost elements (s+m)", n_ghostele, myrank == 0);
        printf(
            "    "
            "+----------------------+----------------+----------------+-----------------+----------"
            "--------+\n");
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "nodes (s)", s_nodes, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "ghost nodes (s)", s_ghostnodes, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "elements (s)", s_elements, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "ghost elements (s)", s_ghostele, myrank == 0);
        printf(
            "    "
            "+----------------------+----------------+----------------+-----------------+----------"
            "--------+\n");
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "nodes (m)", m_nodes, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "ghost nodes (m)", m_ghostnodes, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "elements (m)", m_elements, myrank == 0);
        Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
            "ghost elements (m)", m_ghostele, myrank == 0);
        printf(
            "    "
            "+----------------------+----------------+----------------+-----------------+----------"
            "--------+\n");
      }

      // Print details of parallel distribution for each proc if requested by the user
      const bool printDetails = Core::UTILS::IntegralValue<bool>(
          interface_params().sublist("PARALLEL REDISTRIBUTION"), "PRINT_DISTRIBUTION");
      if (printDetails)
      {
        std::cout << std::endl;
        std::cout << "    Detailed distribution:" << std::endl;
        printf("    +-----+-----------------+--------------+-----------------+--------------+\n");
        printf("    | PID |   n_rownodes    | n_ghostnodes |  n_rowelements  |  n_ghostele  |\n");
        printf("    +-----+-----------------+--------------+-----------------+--------------+\n");
        for (int npid = 0; npid < numproc; ++npid)
        {
          printf("    | %3d | Total %9d | %12d | Total %9d | %12d |\n", npid, n_nodes[npid],
              n_ghostnodes[npid], n_elements[npid], n_ghostele[npid]);
          printf("    |     | Slave %9d | %12d | Slave %9d | %12d |\n", s_nodes[npid],
              s_ghostnodes[npid], s_elements[npid], s_ghostele[npid]);
          printf("    |     | Master %8d | %12d | Master %8d | %12d |\n", m_nodes[npid],
              m_ghostnodes[npid], m_elements[npid], m_ghostele[npid]);
          printf("    +-----+-----------------+--------------+-----------------+--------------+\n");
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  add mortar node (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AddMortarNode(Teuchos::RCP<Mortar::Node> mrtrnode)
{
  idiscret_->AddNode(mrtrnode);
}

/*----------------------------------------------------------------------*
 |  add mortar element (public)                              mwgee 10/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AddMortarElement(Teuchos::RCP<Mortar::Element> mrtrele)
{
  // check for quadratic 2d slave elements to be modified
  if (mrtrele->IsSlave() && (mrtrele->Shape() == Core::FE::CellType::line3 ||
                                mrtrele->Shape() == Core::FE::CellType::nurbs3))
    quadslave_ = true;

  // check for quadratic 3d slave elements to be modified
  if (mrtrele->IsSlave() && (mrtrele->Shape() == Core::FE::CellType::quad9 ||
                                mrtrele->Shape() == Core::FE::CellType::quad8 ||
                                mrtrele->Shape() == Core::FE::CellType::tri6 ||
                                mrtrele->Shape() == Core::FE::CellType::nurbs8 ||
                                mrtrele->Shape() == Core::FE::CellType::nurbs9))
    quadslave_ = true;

  idiscret_->add_element(mrtrele);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::fill_complete_new(const bool isFinalParallelDistribution, const int maxdof)
{
  FOUR_C_THROW("Not implemented for meshtying.");
}

/*----------------------------------------------------------------------*
 |  finalize construction of interface (public)              mwgee 10/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::fill_complete(
    const bool isFinalParallelDistribution, const int maxdof, const double meanVelocity)
{
  TEUCHOS_FUNC_TIME_MONITOR("Mortar::Interface::fill_complete");

  // store maximum global dof ID handed in
  // this ID is later needed when setting up the Lagrange multiplier
  // dof map, which of course must not overlap with existing dof ranges
  maxdofglobal_ = maxdof;

  // we'd like to call idiscret_.fill_complete(true,false,false) but this
  // will assign all nodes new degrees of freedom which we don't want.
  // We would like to use the degrees of freedom that were stored in the
  // mortar nodes. To do so, we have to create and set our own
  // version of a DofSet class before we call fill_complete on the
  // interface discretization.
  // Our special dofset class will not assign new dofs but will assign the
  // dofs stored in the nodes.
  {
    Teuchos::RCP<Mortar::DofSet> mrtrdofset = Teuchos::rcp(new Mortar::DofSet());
    Discret().ReplaceDofSet(mrtrdofset);
    // do not assign dofs yet, we'll do this below after
    // shuffling around of nodes and elements (saves time)
    Discret().fill_complete(false, false, false);
  }

  // check whether crosspoints / edge nodes shall be considered or not
  initialize_cross_points();

  // check for const/linear interpolation of 2D/3D quadratic Lagrange multipliers
  initialize_lag_mult_const();
  initialize_lag_mult_lin();

  // check/init corner/edge modification
  initialize_corner_edge();

  // later we might export node and element column map to extended or even FULL overlap,
  // thus store the standard column maps first
  // get standard nodal column map (overlap=1)
  oldnodecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().NodeColMap())));
  // get standard element column map (overlap=1)
  oldelecolmap_ = Teuchos::rcp(new Epetra_Map(*(Discret().ElementColMap())));

  extend_interface_ghosting(isFinalParallelDistribution, meanVelocity);

  // make sure discretization is complete
  Discret().fill_complete(isFinalParallelDistribution, false, false);

  // ghost also parent elements according to the ghosting strategy of the interface (atm just for
  // poro)
  if (interface_data_->IsPoro())
  {
    if (interface_data_->PoroType() == Inpar::Mortar::poroscatra)
      PoroElastScaTra::UTILS::create_volume_ghosting(Discret());
    else
      PoroElast::UTILS::create_volume_ghosting(Discret());
  }
  else if (imortar_.isParameter("STRATEGY"))
  {
    if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM") ==
        Inpar::Mortar::algorithm_gpts)
      create_volume_ghosting();
  }

  // need row and column maps of slave and master nodes / elements / dofs
  // separately so we can easily address them
  update_master_slave_sets();

  // initialize node and element data container
  initialize_data_container();

  // Communicate quadslave status among ALL processors
  communicate_quad_slave_status_among_all_procs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::communicate_quad_slave_status_among_all_procs()
{
  int localstatus = static_cast<int>(quadslave_);
  int globalstatus = 0;
  Comm().SumAll(&localstatus, &globalstatus, 1);
  quadslave_ = static_cast<bool>(globalstatus);
}

/*----------------------------------------------------------------------*
 |  Check and initialize corner/edge contact                 farah 07/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::initialize_corner_edge()
{
  // if linear LM for quad displacements return!
  // TODO: this case needs a special treatment
  bool lagmultlin = (Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(
                         interface_params(), "LM_QUAD") == Inpar::Mortar::lagmult_lin);

  if (lagmultlin) return;

  for (int i = 0; i < (Discret().NodeRowMap())->NumMyElements(); ++i)
  {
    // static_cast to the corresponding mortar/contact/friction/... node
    // or element would be enough in all places
    // performance loss is negligible when using a dynamic_cast instead
    // but safety is increased enormously
    auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lRowNode(i));

    // remove bound/corner/edge flag for master nodes!
    if (!node->IsSlave() && node->IsOnCorner()) node->SetOnCorner() = false;
    if (!node->IsSlave() && node->IsOnEdge()) node->SetOnEdge() = false;

    // candidates are slave nodes with only 1 adjacent Mortar::Element
    if (node->IsSlave() && node->IsOnCornerEdge())
    {
      node->SetSlave() = false;
    }
  }
}


/*----------------------------------------------------------------------*
 |  Check and initialize cross points                        farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::initialize_cross_points()
{
  // check whether crosspoints / edge nodes shall be considered or not
  bool crosspoints = Core::UTILS::IntegralValue<int>(interface_params(), "CROSSPOINTS");

  // modify crosspoints / edge nodes
  if (crosspoints)
  {
    // only applicable for 2D problems up to now
    if (Dim() == 3)
      FOUR_C_THROW("Crosspoint / edge node modification not yet implemented for 3D problems.");

    // Detect relevant nodes on slave side
    const int numRowNodes = Discret().NodeRowMap()->NumMyElements();
    for (int i = 0; i < numRowNodes; ++i)
    {
      // static_cast to the corresponding mortar/contact/friction/... node
      // or element would be enough in all places
      // performance loss is negligible when using a dynamic_cast instead
      // but safety is increased enormously
      auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with only 1 adjacent Mortar::Element
      if (node->IsSlave() && node->NumElement() == 1)
      {
        // case1: linear shape functions, boundary nodes already found
        if ((node->Elements()[0])->num_node() == 2)
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
}


/*----------------------------------------------------------------------*
 |  Check and initialize for lin lagmult interpolation       farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::initialize_lag_mult_lin()
{
  // check for linear interpolation of 2D/3D quadratic Lagrange multipliers
  bool lagmultlin = (Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(
                         interface_params(), "LM_QUAD") == Inpar::Mortar::lagmult_lin);

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
      auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with shape line3 (2D), tri6 and quad8/9 (3D)
      if (node->IsSlave())
      {
        // search the first adjacent element
        Core::FE::CellType shape = (node->Elements()[0])->Shape();

        // which discretization type
        switch (shape)
        {
          // line3 contact elements (= quad8/9 or tri6 discretizations)
          case Core::FE::CellType::line3:
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
          case Core::FE::CellType::tri6:
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
          case Core::FE::CellType::quad8:
            // quad9 contact elements (= hex27 discretizations)
          case Core::FE::CellType::quad9:
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
            FOUR_C_THROW(
                "Lin/Lin interpolation of LM only for line3/tri6/quad8/quad9 mortar elements");
            break;
          }
        }  // switch(Shape)
      }    // if (IsSlave())
    }      // for-loop
  }
}



/*----------------------------------------------------------------------*
 |  Check and initialize for const lagmult interpolation     seitz 09/17|
 *----------------------------------------------------------------------*/
void Mortar::Interface::initialize_lag_mult_const()
{
  if ((Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(interface_params(), "LM_QUAD") ==
          Inpar::Mortar::lagmult_const))
  {
    // modified treatment slave side nodes:
    // only the center-node carries LM

    // loop over all elements
    for (int i = 0; i < Discret().NodeRowMap()->NumMyElements(); ++i)
    {
      // get node and cast to cnode
      auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lRowNode(i));

      // candidates are slave nodes with shape line3 (2D), tri6 and quad8/9 (3D)
      if (node->IsSlave())
      {
        // search the first adjacent element
        Core::FE::CellType shape = (node->Elements()[0])->Shape();

        // which discretization type
        switch (shape)
        {
          // line3 contact elements (= quad8/9 or tri6 discretizations)
          case Core::FE::CellType::line3:
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
          case Core::FE::CellType::quad9:
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
            FOUR_C_THROW(
                "Lin/Lin interpolation of LM only for line3/tri6/quad8/quad9 mortar elements");
            break;
          }
        }  // switch(Shape)
      }    // if (IsSlave())
    }      // for-loop
  }
}


/*----------------------------------------------------------------------*
 |  Initialize Data Container for nodes and elements         farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::initialize_data_container()
{
  // initialize node data container
  // (include slave side boundary nodes / crosspoints)
  const int numMySlaveColumnNodes = SlaveColNodesBound()->NumMyElements();
  for (int i = 0; i < numMySlaveColumnNodes; ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    Core::Nodes::Node* node = Discret().gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %i", gid);
    auto* mnode = dynamic_cast<Node*>(node);

    //********************************************************
    // NOTE: depending on which kind of node this really is,
    // i.e. mortar, contact or friction node, several derived
    // versions of the initialize_data_container() methods will
    // be called here, apart from the base class version.
    //********************************************************

    // initialize container if not yet initialized before
    mnode->initialize_data_container();
    if (interface_data_->IsPoro())  // initialize just for poro contact case!
      mnode->initialize_poro_data_container();
    if (ehl_) mnode->initialize_ehl_data_container();
  }
  if (interface_data_
          ->IsPoro())  // as velocities of structure and fluid exist also on master nodes!!!
  {
    const Teuchos::RCP<Epetra_Map> masternodes = Core::LinAlg::AllreduceEMap(*(MasterRowNodes()));
    // initialize poro node data container for master nodes!!!
    for (int i = 0; i < masternodes()->NumMyElements(); ++i)
    {
      int gid = masternodes()->GID(i);
      Core::Nodes::Node* node = Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %i", gid);
      auto* mnode = dynamic_cast<Node*>(node);

      // ATM just implemented for ContactNode ... otherwise error!!!

      // initialize container if not yet initialized before
      mnode->initialize_poro_data_container();
    }
  }

  // initialize element data container
  const int numMySlaveColumnElements = SlaveColElements()->NumMyElements();
  for (int i = 0; i < numMySlaveColumnElements; ++i)
  {
    int gid = SlaveColElements()->GID(i);
    Core::Elements::Element* ele = Discret().gElement(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    auto* mele = dynamic_cast<Mortar::Element*>(ele);

    // initialize container if not yet initialized before
    mele->initialize_data_container();
  }

  if (interface_data_->IsPoro())
  {
    // initialize master element data container
    for (int i = 0; i < MasterColElements()->NumMyElements(); ++i)
    {
      int gid = MasterColElements()->GID(i);
      Core::Elements::Element* ele = Discret().gElement(gid);
      if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
      auto* mele = dynamic_cast<Mortar::Element*>(ele);

      // initialize container if not yet initialized before
      mele->initialize_data_container();
    }
  }

  if (interface_params().isParameter("ALGORITHM"))
  {
    if (Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(interface_params(), "ALGORITHM") ==
        Inpar::Mortar::algorithm_gpts)
    {
      const int numMyMasterColumnElements = MasterColElements()->NumMyElements();
      for (int i = 0; i < numMyMasterColumnElements; ++i)
        dynamic_cast<Mortar::Element*>(Discret().gElement(MasterColElements()->GID(i)))
            ->initialize_data_container();
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BINSTRATEGY::BinningStrategy> Mortar::Interface::setup_binning_strategy(
    const double meanVelocity)
{
  // Initialize eXtendedAxisAlignedBoundingBox (XAABB)
  Core::LinAlg::Matrix<3, 2> XAABB(false);
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    XAABB(dim, 0) = +1.0e12;
    XAABB(dim, 1) = -1.0e12;
  }

  // loop all slave nodes and merge XAABB with their eXtendedAxisAlignedBoundingBox
  for (int lid = 0; lid < SlaveColNodes()->NumMyElements(); ++lid)
  {
    int gid = SlaveColNodes()->GID(lid);
    Core::Nodes::Node* node = Discret().gNode(gid);
    if (!node)
      FOUR_C_THROW(
          "Cannot find node with gid %i in discretization '%s'.", gid, Discret().Name().c_str());
    auto* mtrnode = dynamic_cast<Mortar::Node*>(node);

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
  Comm().MinAll(locmin, globmin, 3);
  Comm().MaxAll(locmax, globmax, 3);

  // compute cutoff radius:
  double global_slave_max_edge_size = -1.0;
  for (int lid = 0; lid < SlaveColElements()->NumMyElements(); ++lid)
  {
    int gid = SlaveColElements()->GID(lid);
    Core::Elements::Element* ele = Discret().gElement(gid);
    if (!ele)
      FOUR_C_THROW(
          "Cannot find element with gid %i in discretization '%s'.", gid, Discret().Name().c_str());
    auto* mtrele = dynamic_cast<Mortar::Element*>(ele);

    // to be thought about, whether this is enough (safety = 2??)
    double slave_max_edge_size = mtrele->MaxEdgeSize();
    global_slave_max_edge_size = std::max(slave_max_edge_size, global_slave_max_edge_size);
  }

  double cutoff = -1.0;
  Comm().MaxAll(&global_slave_max_edge_size, &cutoff, 1);

  // extend cutoff based on problem interface velocity
  // --> only for contact problems
  if (meanVelocity >= 1e-12)
  {
    const double dt = interface_params().get<double>("TIMESTEP");
    cutoff = cutoff + 2 * dt * meanVelocity;
  }

  // increase XAABB by 2x cutoff radius
  std::stringstream domain_bounding_box_stream;
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    domain_bounding_box_stream << globmin[dim] - cutoff << " ";
  }
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    domain_bounding_box_stream << globmax[dim] + cutoff << " ";
  }

  Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();

  binning_params.set<double>("BIN_SIZE_LOWER_BOUND", cutoff);
  binning_params.set<std::string>("DOMAINBOUNDINGBOX", domain_bounding_box_stream.str());
  Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
      "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
      binning_params);
  Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
      Teuchos::rcp(new BINSTRATEGY::BinningStrategy(binning_params,
          Global::Problem::Instance()->OutputControlFile(), Comm(), Comm().MyPID()));

  return binningstrategy;
}


/*----------------------------------------------------------------------*
 |  redistribute interface (public)                           popp 08/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::Redistribute()
{
  std::stringstream ss;
  ss << "Mortar::Interface::Redistribute of '" << Discret().Name() << "'";
  TEUCHOS_FUNC_TIME_MONITOR(ss.str());

  const Teuchos::ParameterList& mortarParallelRedistParams =
      interface_params().sublist("PARALLEL REDISTRIBUTION");

  // make sure we are supposed to be here
  if (Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(mortarParallelRedistParams,
          "PARALLEL_REDIST") == Inpar::Mortar::ParallelRedist::redist_none)
    FOUR_C_THROW("You are not supposed to be here...");

  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());
  const int myrank = comm->MyPID();
  const int numproc = comm->NumProc();
  Teuchos::Time time("", true);

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
    std::cout << "\nRedistributing interface '" << Discret().Name() << "' .........\n";
    std::cout << "Procs used for redistribution: " << sproc << " / " << mproc << " (S / M)\n";
  }

  //**********************************************************************
  // (2) SLAVE redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> srownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> scolnodes = Teuchos::null;

  {
    std::stringstream ss_slave;
    ss_slave << "Mortar::Interface::Redistribute of '" << Discret().Name() << "' (slave)";
    TEUCHOS_FUNC_TIME_MONITOR(ss_slave.str());

    Teuchos::RCP<const Epetra_CrsGraph> snodegraph =
        Core::Rebalance::BuildGraph(idiscret_, sroweles);

    Teuchos::ParameterList rebalanceParams;
    rebalanceParams.set<std::string>("num parts", std::to_string(sproc));
    rebalanceParams.set<std::string>("imbalance tol", std::to_string(imbalance_tol));

    std::tie(srownodes, scolnodes) =
        Core::Rebalance::RebalanceNodeMaps(snodegraph, rebalanceParams);
  }

  //**********************************************************************
  // (3) MASTER redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> mrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> mcolnodes = Teuchos::null;

  {
    std::stringstream ss_master;
    ss_master << "Mortar::Interface::Redistribute of '" << Discret().Name() << "' (master)";
    TEUCHOS_FUNC_TIME_MONITOR(ss_master.str());

    redistribute_master_side(mrownodes, mcolnodes, mroweles, comm, mproc, imbalance_tol);
  }

  //**********************************************************************
  // (4) Merge global interface node row and column map
  //**********************************************************************
  // merge node maps from slave and master parts
  Teuchos::RCP<Epetra_Map> rownodes = Core::LinAlg::MergeMap(srownodes, mrownodes, false);
  Teuchos::RCP<Epetra_Map> colnodes = Core::LinAlg::MergeMap(scolnodes, mcolnodes, false);

  //**********************************************************************
  // (5) Get partitioning information into discretization
  //**********************************************************************
  {
    std::stringstream ss_comm;
    ss_comm << "Mortar::Interface::Redistribute of '" << Discret().Name() << "' (communicate)";
    TEUCHOS_FUNC_TIME_MONITOR(ss_comm.str());

    // build reasonable element maps from the already valid and final node maps
    // (note that nothing is actually redistributed in here)
    auto const& [roweles, coleles] = Discret().build_element_row_column(*rownodes, *colnodes);

    // export nodes and elements to the row map
    Discret().ExportRowNodes(*rownodes);
    Discret().ExportRowElements(*roweles);

    // export nodes and elements to the column map (create ghosting)
    Discret().ExportColumnNodes(*colnodes);
    Discret().export_column_elements(*coleles);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::Interface::redistribute_master_side(Teuchos::RCP<Epetra_Map>& rownodes,
    Teuchos::RCP<Epetra_Map>& colnodes, const Teuchos::RCP<Epetra_Map>& roweles,
    const Teuchos::RCP<Epetra_Comm>& comm, const int parts, const double imbalance) const
{
  if (not has_ma_sharing_ref_interface())
  {
    // call parallel redistribution
    Teuchos::RCP<const Epetra_CrsGraph> nodegraph = Core::Rebalance::BuildGraph(idiscret_, roweles);

    Teuchos::ParameterList rebalanceParams;
    rebalanceParams.set<std::string>("num parts", std::to_string(parts));
    rebalanceParams.set<std::string>("imbalance tol", std::to_string(imbalance));

    std::tie(rownodes, colnodes) = Core::Rebalance::RebalanceNodeMaps(nodegraph, rebalanceParams);
  }
  else
  {
    Core::Rebalance::RebalanceInAccordanceWithReference(
        *get_ma_sharing_ref_interface_ptr()->MasterRowNodes(), *MasterRowNodes(), rownodes);

    Core::Rebalance::RebalanceInAccordanceWithReference(
        *get_ma_sharing_ref_interface_ptr()->MasterColNodes(), *MasterColNodes(), colnodes);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_parallel_layout_and_data_structures(const bool perform_rebalancing,
    const bool enforce_ghosting_update, const int maxdof, const double meanVelocity)
{
  FOUR_C_THROW("Not implemented for meshtying.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::extend_interface_ghosting_safely(const double meanVelocity)
{
  FOUR_C_THROW("Not implemented for meshtying.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::extend_interface_ghosting(
    const bool isFinalParallelDistribution, const double meanVelocity)
{
  //*****REDUNDANT SLAVE AND MASTER STORAGE*****
  if (interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::redundant_all)
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
    // Core::LinAlg::Gather<int>(stproc,rtproc,Comm().NumProc(),allproc.data(),Comm());
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
    Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

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
    Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

    // build complete overlapping map of elements (on ALL processors)
    Teuchos::RCP<Epetra_Map> newelecolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
    sdata.clear();
    rdata.clear();
    allproc.clear();

    // redistribute the discretization of the interface according to the
    // new column layout
    Discret().ExportColumnNodes(*newnodecolmap);
    Discret().export_column_elements(*newelecolmap);
  }

  //*****ONLY REDUNDANT MASTER STORAGE*****
  else if (interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::redundant_master)
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
    // Core::LinAlg::Gather<int>(stproc,rtproc,Comm().NumProc(),allproc.data(),Comm());
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
      Core::Nodes::Node* node = Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);
      if (!mrtrnode->IsSlave()) sdata.push_back(gid);
    }

    // gather all master row node gids redundantly
    std::vector<int> rdata;
    Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

    // add my own slave column node ids (non-redundant, standard overlap)
    const Epetra_Map* nodecolmap = Discret().NodeColMap();
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      Core::Nodes::Node* node = Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);
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
      Core::Elements::Element* ele = Discret().gElement(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      auto* mrtrele = dynamic_cast<Mortar::Element*>(ele);
      if (!mrtrele->IsSlave()) sdata.push_back(gid);
    }

    // gather all gids of elements redundantly
    rdata.resize(0);
    Core::LinAlg::Gather<int>(sdata, rdata, (int)allproc.size(), allproc.data(), Comm());

    // add my own slave column node ids (non-redundant, standard overlap)
    const Epetra_Map* elecolmap = Discret().ElementColMap();
    for (int i = 0; i < elecolmap->NumMyElements(); ++i)
    {
      int gid = elecolmap->GID(i);
      Core::Elements::Element* ele = Discret().gElement(gid);
      if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
      auto* mrtrele = dynamic_cast<Mortar::Element*>(ele);
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
    Discret().export_column_elements(*newelecolmap);
  }

  //*****NON-REDUNDANT STORAGE*****
  else if (interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::roundrobin ||
           interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::binning)
  {
    // nothing to do here, we work with the given non-redundant distribution
    // of both slave and master nodes to the individual processors. However
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
    // Core::LinAlg::Gather<int>(stproc,rtproc,Comm().NumProc(),allproc.data(),Comm());
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
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
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
        Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, Comm()));
    rdata.clear();

    // redistribute the discretization of the interface according to the
    // new (=old) node / element column layout
    Discret().ExportColumnNodes(*newnodecolmap);
    Discret().export_column_elements(*newelecolmap);

    if (interface_data_->GetExtendGhosting() == Inpar::Mortar::ExtendGhosting::binning)
    {
      /* We have to update the row/column maps split into master/slave. We start from the new
       * node/element column maps. Since we don't have row maps at this point, we can/have to pass
       * the column map as row map.
       */
      update_master_slave_element_maps(*newelecolmap, *newelecolmap);
      update_master_slave_node_maps(*newnodecolmap, *newnodecolmap);

      /* Ask the discretization to initialize the elements. We need this, since the setup of the
       * binning strategy relies on some element information.
       * Note: This should be cheap, since we don't assign new degrees of freedom. Still, we skip
       * initialization of elements, if we know, that the discretization will be redistributed
       * again.
       */
      Discret().fill_complete(false, isFinalParallelDistribution, false);

      // Create the binning strategy
      Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
          setup_binning_strategy(meanVelocity);

      // fill master and slave elements into bins
      std::map<int, std::set<int>> slavebinelemap;
      binningstrategy->distribute_eles_to_bins(Discret(), slavebinelemap, true);
      std::map<int, std::set<int>> masterbinelemap;
      binningstrategy->distribute_eles_to_bins(Discret(), masterbinelemap, false);

      // Extend ghosting of the master elements
      std::map<int, std::set<int>> ext_bin_to_ele_map;
      Teuchos::RCP<const Epetra_Map> extendedmastercolmap =
          binningstrategy->ExtendElementColMap(slavebinelemap, masterbinelemap, ext_bin_to_ele_map,
              Teuchos::null, Teuchos::null, newelecolmap.get());

      // adapt layout to extended ghosting in the discretization
      // first export the elements according to the processor local element column maps
      Discret().export_column_elements(*extendedmastercolmap);

      // get the node ids of the elements that are to be ghosted and create a proper node column map
      // for their export
      std::set<int> nodes;
      for (int lid = 0; lid < extendedmastercolmap->NumMyElements(); ++lid)
      {
        Core::Elements::Element* ele = Discret().gElement(extendedmastercolmap->GID(lid));
        const int* nodeids = ele->NodeIds();
        for (int inode = 0; inode < ele->num_node(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), colnodes.data(), 0, Comm()));

      // now ghost the nodes
      Discret().ExportColumnNodes(*nodecolmap);
    }
  }

  //*****INVALID CASES*****
  else
  {
    FOUR_C_THROW("Invalid redundancy type.");
  }
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::CreateSearchTree()
{
  // warning
#ifdef MORTARGMSHCTN
  if (Dim() == 3 && Comm().MyPID() == 0)
  {
    std::cout << "\n*****************************************************************\n";
    std::cout << "GMSH output of all mortar tree nodes in 3D needs a lot of memory!\n";
    std::cout << "*****************************************************************\n";
  }
#endif
  // binary tree search
  if (SearchAlg() == Inpar::Mortar::search_binarytree)
  {
    // create fully overlapping map of all master elements
    // for non-redundant storage (RRloop) we handle the master elements
    // like the slave elements --> melecolmap_
    auto strat = Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
        interface_params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");

    // get update type of binary tree
    auto updatetype = Core::UTILS::IntegralValue<Inpar::Mortar::BinaryTreeUpdateType>(
        interface_params(), "BINARYTREE_UPDATETYPE");

    Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;
    switch (strat)
    {
      case Inpar::Mortar::ExtendGhosting::roundrobin:
      case Inpar::Mortar::ExtendGhosting::binning:
      {
        melefullmap = melecolmap_;
        break;
      }
      case Inpar::Mortar::ExtendGhosting::redundant_all:
      case Inpar::Mortar::ExtendGhosting::redundant_master:
      {
        melefullmap = Core::LinAlg::AllreduceEMap(*melerowmap_);
        break;
      }
      default:
      {
        FOUR_C_THROW("Unknown strategy to deal with interface ghosting.");
        break;
      }
    }

    // create binary tree object for search and setup tree
    binarytree_ = Teuchos::rcp(new Mortar::BinaryTree(
        Discret(), selecolmap_, melefullmap, Dim(), SearchParam(), updatetype, SearchUseAuxPos()));
    // initialize the binary tree
    binarytree_->Init();
  }
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                 popp 11/09|
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_master_slave_sets()
{
  update_master_slave_node_maps();
  update_master_slave_element_maps();
  update_master_slave_dof_maps();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_master_slave_dof_maps()
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
    Core::Nodes::Node* node = Discret().gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);
    bool isslave = mrtrnode->IsSlave();
    const int numdof = mrtrnode->NumDof();

    if (isslave)
      for (int j = 0; j < numdof; ++j) sc.push_back(mrtrnode->Dofs()[j]);
    else
      for (int j = 0; j < numdof; ++j) mc.push_back(mrtrnode->Dofs()[j]);

    if (Discret().NodeRowMap()->MyGID(gid))
    {
      if (isslave)
        for (int j = 0; j < numdof; ++j) sr.push_back(mrtrnode->Dofs()[j]);
      else
        for (int j = 0; j < numdof; ++j) mr.push_back(mrtrnode->Dofs()[j]);
    }
  }

  sdofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sr.size(), sr.data(), 0, Comm()));
  sdofcolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sc.size(), sc.data(), 0, Comm()));
  mdofrowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mr.size(), mr.data(), 0, Comm()));
  mdofcolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mc.size(), mc.data(), 0, Comm()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_master_slave_element_maps()
{
  update_master_slave_element_maps(*Discret().ElementRowMap(), *Discret().ElementColMap());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_master_slave_element_maps(
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
    bool isslave = dynamic_cast<Mortar::Element*>(Discret().gElement(gid))->IsSlave();

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

  selerowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sr.size(), sr.data(), 0, Comm()));
  selecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sc.size(), sc.data(), 0, Comm()));
  melerowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mr.size(), mr.data(), 0, Comm()));
  melecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mc.size(), mc.data(), 0, Comm()));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_master_slave_node_maps()
{
  update_master_slave_node_maps(*Discret().NodeRowMap(), *Discret().NodeColMap());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::update_master_slave_node_maps(
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
    bool isslave = dynamic_cast<Mortar::Node*>(Discret().gNode(gid))->IsSlave();
    bool isonbound = dynamic_cast<Mortar::Node*>(Discret().gNode(gid))->IsOnBoundorCE();

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

  snoderowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sr.size(), sr.data(), 0, Comm()));
  snodecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)sc.size(), sc.data(), 0, Comm()));
  mnoderowmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mr.size(), mr.data(), 0, Comm()));
  mnodecolmap_ = Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mc.size(), mc.data(), 0, Comm()));

  snoderowmapbound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)srb.size(), srb.data(), 0, Comm()));
  snodecolmapbound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)scb.size(), scb.data(), 0, Comm()));
  mnoderowmapnobound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mrb.size(), mrb.data(), 0, Comm()));
  mnodecolmapnobound_ =
      Teuchos::rcp<Epetra_Map>(new Epetra_Map(-1, (int)mcb.size(), mcb.data(), 0, Comm()));

  // build exporter
  interface_data_->sl_exporter_ptr() = Teuchos::rcp(
      new Core::Communication::Exporter(*snoderowmapbound_, *snodecolmapbound_, Comm()));
}

/*----------------------------------------------------------------------*
 |  restrict slave sets to actual meshtying zone              popp 08/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::RestrictSlaveSets()
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
      Core::Nodes::Node* node = Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);
      int istied = (int)mrtrnode->IsTiedSlave();

      if (istied && snodecolmap_->MyGID(gid)) sc.push_back(gid);
      if (istied && snoderowmap_->MyGID(gid)) sr.push_back(gid);
    }

    snoderowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sr.size(), sr.data(), 0, Comm()));
    snodecolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sc.size(), sc.data(), 0, Comm()));
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
      Core::Nodes::Node* node = Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);
      const bool istied = mrtrnode->IsTiedSlave();
      const int numdof = mrtrnode->NumDof();

      if (istied && snodecolmap_->MyGID(gid))
        for (int j = 0; j < numdof; ++j) sc.push_back(mrtrnode->Dofs()[j]);

      if (istied && snoderowmap_->MyGID(gid))
        for (int j = 0; j < numdof; ++j) sr.push_back(mrtrnode->Dofs()[j]);
    }

    sdofrowmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sr.size(), sr.data(), 0, Comm()));
    sdofcolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)sc.size(), sc.data(), 0, Comm()));
  }
}

/*----------------------------------------------------------------------*
 |  update Lagrange multiplier set (dofs)                     popp 08/10|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> Mortar::Interface::UpdateLagMultSets(
    int offset_if, const bool& redistributed, const Epetra_Map& ref_map) const
{
  if (redistributed)
  {
    return redistribute_lag_mult_sets();
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
  Comm().SumAll(localnumlmdof.data(), globalnumlmdof.data(), Comm().NumProc());

  // compute offset for LM dof initialization for all procs
  int offset = 0;
  for (int k = 0; k < Comm().MyPID(); ++k) offset += globalnumlmdof[k];

  // loop over all slave dofs and initialize LM dofs
  lmdof.reserve(ref_map.NumMyElements());
  for (int i = 0; i < ref_map.NumMyElements(); ++i)
    lmdof.push_back(MaxDofGlobal() + 1 + offset_if + offset + i);

  // create interface LM map
  // (if maxdofglobal_ == 0, we do not want / need this)
  if (MaxDofGlobal() > 0)
    return Teuchos::rcp(new Epetra_Map(-1, (int)lmdof.size(), lmdof.data(), 0, Comm()));

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mortar::Interface::store_unredistributed_maps()
{
  psdofrowmap_ = Teuchos::rcp(new Epetra_Map(*sdofrowmap_));
  interface_data_->PMDofRowMap() = Teuchos::rcp(new Epetra_Map(*mdofrowmap_));
  plmdofmap_ = Teuchos::rcp(new Epetra_Map(*lmdofmap_));

  interface_data_->PSNodeRowMap() = Teuchos::rcp(new Epetra_Map(*snoderowmap_));
  interface_data_->PMNodeRowMap() = Teuchos::rcp(new Epetra_Map(*mnoderowmap_));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> Mortar::Interface::redistribute_lag_mult_sets() const
{
  if (plmdofmap_.is_null()) FOUR_C_THROW("The plmdofmap_ is not yet initialized!");
  if (psdofrowmap_.is_null()) FOUR_C_THROW("The psdofrowmap_ is not yet initialized!");

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
    Comm().Broadcast(g_related_gids.data(), 2 * num_mygids, p);

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
  return Teuchos::rcp(new Epetra_Map(-1, (int)lmdof.size(), lmdof.data(), 0, Comm()));
}

/*----------------------------------------------------------------------*
 |  initialize / reset mortar interface                       popp 01/08|
 *----------------------------------------------------------------------*/
void Mortar::Interface::Initialize()
{
  // loop over all nodes to reset stuff (fully overlapping column map)
  for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
  {
    auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lColNode(i));

    // reset feasible projection and segmentation status
    node->HasProj() = false;
    node->HasSegment() = false;
  }

  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i = 0; i < SlaveColNodesBound()->NumMyElements(); ++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    Core::Nodes::Node* node = Discret().gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* monode = dynamic_cast<Mortar::Node*>(node);

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
    Core::Elements::Element* ele = Discret().gElement(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    auto* mele = dynamic_cast<Mortar::Element*>(ele);

    mele->MoData().SearchElements().resize(0);
  }
}


/*----------------------------------------------------------------------*
 |  set current and old deformation state                      popp 12/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::set_state(const enum StateType& statetype, const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case state_new_displacement:
    {
      // alternative method to get vec to full overlap
      Teuchos::RCP<Epetra_Vector> global =
          Teuchos::rcp(new Epetra_Vector(*idiscret_->DofColMap(), false));
      Core::LinAlg::Export(vec, *global);

      // set displacements in interface discretization
      idiscret_->set_state(StateType2String(statetype), global);

      // loop over all nodes to set current displacement
      // (use fully overlapping column map)
      for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
      {
        auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lColNode(i));
        const int numdof = node->NumDof();
        std::vector<double> mydisp(numdof);
        std::vector<int> lm(numdof);

        for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

        Core::FE::ExtractMyValues(*global, mydisp, lm);

        // add mydisp[2]=0 for 2D problems
        if (mydisp.size() < 3) mydisp.resize(3);

        // set current configuration
        for (int j = 0; j < 3; ++j) node->xspatial()[j] = node->X()[j] + mydisp[j];
      }

      // compute element areas
      set_element_areas();
      break;
    }
    case state_lagrange_multiplier:
    {
      // alternative method to get vec to full overlap
      Teuchos::RCP<Epetra_Vector> global =
          Teuchos::rcp(new Epetra_Vector(*idiscret_->DofColMap(), false));
      Core::LinAlg::Export(vec, *global);

      // loop over all nodes to set current displacement
      // (use fully overlapping column map)
      for (int i = 0; i < SlaveColNodes()->NumMyElements(); ++i)
      {
        auto* node = dynamic_cast<Mortar::Node*>(idiscret_->gNode(SlaveColNodes()->GID(i)));
        const int numdof = node->NumDof();
        std::vector<double> mydisp(numdof);
        std::vector<int> lm(numdof);

        for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

        Core::FE::ExtractMyValues(*global, mydisp, lm);

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
      Core::LinAlg::Export(vec, *global);

      // set displacements in interface discretization
      idiscret_->set_state(StateType2String(statetype), global);

      // loop over all nodes to set current displacement
      // (use fully overlapping column map)
      for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
      {
        auto* node = dynamic_cast<Mortar::Node*>(idiscret_->lColNode(i));
        const int numdof = node->NumDof();
        std::vector<double> myolddisp(numdof);
        std::vector<int> lm(numdof);

        for (int j = 0; j < numdof; ++j) lm[j] = node->Dofs()[j];

        Core::FE::ExtractMyValues(*global, myolddisp, lm);

        // add mydisp[2]=0 for 2D problems
        if (myolddisp.size() < 3) myolddisp.resize(3);

        // set old displacement
        for (int j = 0; j < 3; ++j) node->uold()[j] = myolddisp[j];
      }

      break;
    }
    default:
    {
      FOUR_C_THROW(
          "The given state type is unsupported! (type = %s)", StateType2String(statetype).c_str());
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::set_element_areas()
{
  // loop over all elements to set current element length / area
  // (use standard slave column map)
  for (int i = 0; i < SlaveColElements()->NumMyElements(); ++i)
  {
    int gid = SlaveColElements()->GID(i);
    Core::Elements::Element* ele = Discret().gElement(gid);
    if (!ele) FOUR_C_THROW("Cannot find ele with gid %i", gid);
    auto* mele = dynamic_cast<Mortar::Element*>(ele);

    mele->MoData().Area() = mele->ComputeArea();
  }
}


/*----------------------------------------------------------------------*
 |  evaluate geometric setting (create integration cells)    farah 01/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::EvaluateGeometry(std::vector<Teuchos::RCP<Mortar::IntCell>>& intcells)
{
  // time measurement
  Comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // check
  if (Dim() == 2) FOUR_C_THROW("Geometry evaluation for mortar interface only for 3D problems!");

  auto algo = Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  if (algo == Inpar::Mortar::algorithm_nts)
    FOUR_C_THROW("Geometry evaluation only for mortar problems!");

  // interface needs to be complete
  if (!Filled() && Comm().MyPID() == 0)
    FOUR_C_THROW("fill_complete() not called on interface %", id_);

  // clear vector
  intcells.clear();

  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (SearchAlg() == Inpar::Mortar::search_bfele)
    evaluate_search_brute_force(SearchParam());
  else if (SearchAlg() == Inpar::Mortar::search_binarytree)
    evaluate_search_binarytree();
  else
    FOUR_C_THROW("Invalid search algorithm");

  // create normals
  evaluate_nodal_normals();

  // export nodal normals to slave node column map
  // this call is very expensive and the computation
  // time scales directly with the proc number !
  export_nodal_normals();

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    Core::Elements::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %", gid1);
    auto* selement = dynamic_cast<Mortar::Element*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      Core::Elements::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      auto* melement = dynamic_cast<Mortar::Element*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      //********************************************************************
      if (selement->IsQuad())
      {
        FOUR_C_THROW("Geometry evaluation only implemented for first order elements!");
      }
      // noquad!
      else
      {
        Mortar::Coupling3d coup(*idiscret_, dim_, false, imortar_, *selement, *melement);

        // do coupling
        coup.evaluate_coupling();

        // set sele and mele id and push into global vector
        for (auto& coupcell : coup.Cells())
        {
          coupcell->SetSlaveId(selement->Id());
          coupcell->SetMasterId(melement->Id());
          intcells.push_back(coupcell);
        }
      }
    }
  }  // end sele loop

  // time measurement
  Comm().Barrier();
  const double evaltime = Teuchos::Time::wallTime() - t_start;

  // time output
  std::cout << "Required time for geometry evaluation: " << evaltime << std::endl;
}


/*----------------------------------------------------------------------*
 |  evaluate mortar coupling (public)                         popp 11/07|
 *----------------------------------------------------------------------*/
void Mortar::Interface::Evaluate(
    int rriter, const int& step, const int& iter, Teuchos::RCP<Mortar::ParamsInterface> mparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR("Mortar::Interface::Evaluate");

  // interface needs to be complete
  if (!Filled()) FOUR_C_THROW("fill_complete() not called on interface %", id_);

  //******************************************
  // Start basic evaluation part of interface
  //******************************************
  // start time measurement
  const double t_start = Teuchos::Time::wallTime();

  // evaluate nodal normals and decide in contact case if
  // this is a nonsmooth or smooth contact
  pre_evaluate(step, iter);

  // evaluation routine for coupling
  evaluate_coupling(*selecolmap_, snoderowmap_.get(), mparams_ptr);

  // do some post operations. nothing happens for standard cases...
  post_evaluate(step, iter);

  // end time on this proc
  const double inttime = Teuchos::Time::wallTime() - t_start;
  //******************************************
  // End basic evaluation part of interface
  //******************************************

  // store integrationtime
  inttime_interface_ = inttime;
}


/*----------------------------------------------------------------------*
 |  protected evaluate routine                               farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_coupling(const Epetra_Map& selecolmap,
    const Epetra_Map* snoderowmap, const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  // decide which type of coupling should be evaluated
  auto algo = Core::UTILS::IntegralValue<Inpar::Mortar::AlgorithmType>(imortar_, "ALGORITHM");

  // smooth contact
  switch (algo)
  {
    //*********************************
    // Mortar Coupling (STS)    (2D/3D)
    // Gauss-Point-To-Segment (GPTS)
    //*********************************
    case Inpar::Mortar::algorithm_mortar:
    case Inpar::Mortar::algorithm_gpts:
    {
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      Mortar::Interface::evaluate_sts(selecolmap, mparams_ptr);
      break;
    }
    //*********************************
    // Segment-to-Line Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_stl:
    {
      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      evaluate_stl();
      break;
    }
    //*********************************
    // Line-to-Segment Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_lts:
    {
      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      evaluate_lts();
      break;
    }
    //*********************************
    // line-to-line Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_ltl:
    {
      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //    (only for contact setting)
      //********************************************************************
      evaluate_ltl();
      break;
    }
    //*********************************
    // Node-to-Segment Coupling (2D/3D)
    //*********************************
    case Inpar::Mortar::algorithm_nts:
    {
      //********************************************************************
      // 1) try to project slave nodes onto master elements
      // 2) evaluate shape functions at projected positions
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      evaluate_nts();
      break;
    }
    //*********************************
    // Node-to-Line Coupling (3D)
    //*********************************
    case Inpar::Mortar::algorithm_ntl:
    {
      FOUR_C_THROW("not yet implemented!");
      break;
    }
    //*********************************
    // Default case
    //*********************************
    default:
    {
      FOUR_C_THROW("Unknown discretization type for constraints!");
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-segment coupl          farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_sts(
    const Epetra_Map& selecolmap, const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR("Mortar::Interface::EvaluateSTS");

  // loop over all slave col elements
  for (int i = 0; i < selecolmap.NumMyElements(); ++i)
  {
    const int gid1 = selecolmap.GID(i);
    Core::Elements::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) FOUR_C_THROW("Cannot find slave element with gid %d", gid1);

    auto* selement = dynamic_cast<Mortar::Element*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized()) continue;

    // empty vector of master element pointers
    std::vector<Mortar::Element*> melements;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      Core::Elements::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) FOUR_C_THROW("Cannot find master element with gid %d", gid2);
      auto* melement = dynamic_cast<Mortar::Element*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      melements.push_back(melement);
    }

    // concrete coupling evaluation routine
    MortarCoupling(selement, melements, mparams_ptr);
  }
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type node-to-segment coupl             farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_nts()
{
  // loop over slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Mortar::Node*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    // vector with possible contacting master eles
    std::vector<Mortar::Element*> meles;

    // fill vector with possibly contacting meles
    FindMEles(*mrtrnode, meles);

    // skip calculation if no meles vector is empty
    if (meles.size() < 1) continue;

    // call interpolation functions
    NTS::MTInterpolator::Impl(meles)->Interpolate(*mrtrnode, meles);
  }
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_lts()
{
  FOUR_C_THROW("Line -to-segment is not available for meshtying problems.");
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-line coupl                farah 07/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_ltl()
{
  FOUR_C_THROW("Line-to-line is not available for meshtying problems.");
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-line coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_stl()
{
  FOUR_C_THROW("Segment-to-line is not available for meshtying problems.");
}

/*----------------------------------------------------------------------*
 |  evaluate nodal normals (public)                           popp 10/11|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_nodal_normals() const
{
  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();
  }
}

/*----------------------------------------------------------------------*
 |  pre evaluate to calc normals                            farah 02/16 |
 *----------------------------------------------------------------------*/
void Mortar::Interface::pre_evaluate(const int& step, const int& iter)
{
  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (SearchAlg() == Inpar::Mortar::search_bfele)
    evaluate_search_brute_force(SearchParam());
  else if (SearchAlg() == Inpar::Mortar::search_binarytree)
    evaluate_search_binarytree();
  else
    FOUR_C_THROW("Invalid search algorithm");

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
  evaluate_nodal_normals();

  // export nodal normals to slave node column map
  // this call is very expensive and the computation
  // time scales directly with the proc number !
  export_nodal_normals();
}


/*----------------------------------------------------------------------*
 |  post evaluate                                           farah 02/16 |
 *----------------------------------------------------------------------*/
void Mortar::Interface::post_evaluate(const int step, const int iter)
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
    FOUR_C_THROW("Gmsh output implemented for a maximum of 99.999 time steps");
  newfilename << step;

  // second index = Newton iteration index
  newfilename << "_";
  if (iter < 10)
    newfilename << 0;
  else if (iter > 99)
    FOUR_C_THROW("Gmsh output implemented for a maximum of 99 iterations");
  newfilename << iter << "_p" << proc << ".pos";

  // rename file
  rename(filename.str().c_str(), newfilename.str().c_str());
#endif  // #ifdef MORTARGMSHCELLS
}


/*----------------------------------------------------------------------*
 |  find meles to snode                                     farah 01/16 |
 *----------------------------------------------------------------------*/
void Mortar::Interface::FindMEles(Node& mrtrnode, std::vector<Mortar::Element*>& meles) const
{
  // clear vector
  meles.clear();

  // get adjacent elements for this node
  Core::Elements::Element** adjeles = mrtrnode.Elements();

  // empty vector of master element pointers
  std::set<int> donebefore;

  for (int j = 0; j < mrtrnode.NumElement(); ++j)
  {
    auto* adjcele = dynamic_cast<Mortar::Element*>(adjeles[j]);

    // skip zero-sized nurbs elements (slave)
    if (adjcele->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < adjcele->MoData().NumSearchElements(); ++k)
    {
      int gid2 = adjcele->MoData().SearchElements()[k];
      Core::Elements::Element* mele = idiscret_->gElement(gid2);
      if (!mele) FOUR_C_THROW("Cannot find master element with gid %", gid2);
      auto* melement = dynamic_cast<Mortar::Element*>(mele);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized()) continue;

      // check uniqueness
      auto iter = donebefore.find(melement->Id());
      if (iter != donebefore.end()) continue;

      donebefore.insert(melement->Id());

      // fill vector
      meles.push_back(melement);
    }  // found eles
  }    // loop over adjacent slave elements
}


/*----------------------------------------------------------------------*
 |  find mnodes to snode                                    farah 01/16 |
 *----------------------------------------------------------------------*/
void Mortar::Interface::find_m_nodes(
    Node& mrtrnode, std::vector<Mortar::Element*>& meles, std::vector<Node*>& mnodes)
{
  // clear vector
  mnodes.clear();

  // check meles
  if (meles.size() < 1) return;

  // set object to guarantee uniqueness of found mnodes
  std::set<int> donebefore;

  for (auto& mele : meles)
  {
    // skip zero-sized nurbs elements (master)
    if (mele->ZeroSized()) continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < mele->num_node(); ++k)
    {
      Core::Nodes::Node* node = mele->Nodes()[k];
      if (!node) FOUR_C_THROW("Cannot find master node");
      auto* mnode = dynamic_cast<Node*>(node);

      // check uniqueness
      auto iter = donebefore.find(mnode->Id());
      if (iter != donebefore.end()) continue;

      donebefore.insert(mnode->Id());

      // fill vector
      mnodes.push_back(mnode);
    }  // found eles
  }    // loop over adjacent slave elements
}


/*----------------------------------------------------------------------*
 |  evaluate nodal normals and store them in map (public)      jb 07/14 |
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_nodal_normals(std::map<int, std::vector<double>>& mynormals)
{
  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();

    int numdofs = mrtrnode->NumDof();
    std::vector<double> temp(numdofs, 0.0);
    for (int j = 0; j < numdofs; j++)
    {
      temp[j] = mrtrnode->MoData().n()[j];
    }
    mynormals.insert(std::pair<int, std::vector<double>>(gid, temp));
  }
}

/*----------------------------------------------------------------------*
 |  export nodal normals (public)                             popp 11/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::export_nodal_normals() const
{
  // create empty data objects
  std::map<int, Teuchos::RCP<Core::LinAlg::SerialDenseMatrix>> triad;

  // build info on row map
  for (int i = 0; i < snoderowmapbound_->NumMyElements(); ++i)
  {
    int gid = snoderowmapbound_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    // fill nodal matrix
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> loc =
        Teuchos::rcp(new Core::LinAlg::SerialDenseMatrix(3, 1));
    (*loc)(0, 0) = mrtrnode->MoData().n()[0];
    (*loc)(1, 0) = mrtrnode->MoData().n()[1];
    (*loc)(2, 0) = mrtrnode->MoData().n()[2];

    triad[gid] = loc;
  }

  // communicate from slave node row to column map

  interface_data_->Exporter().Export(triad);

  // extract info on column map
  for (int i = 0; i < snodecolmapbound_->NumMyElements(); ++i)
  {
    // only do something for ghosted nodes
    int gid = snodecolmapbound_->GID(i);
    if (snoderowmapbound_->MyGID(gid)) continue;

    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    // extract info
    Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> loc = triad[gid];
    mrtrnode->MoData().n()[0] = (*loc)(0, 0);
    mrtrnode->MoData().n()[1] = (*loc)(1, 0);
    mrtrnode->MoData().n()[2] = (*loc)(2, 0);
  }
}

/*----------------------------------------------------------------------*
 |  Search element-based "brute force" (public)               popp 10/08|
 *----------------------------------------------------------------------*/
void Mortar::Interface::evaluate_search_brute_force(const double& eps)
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
  auto strat = Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
      interface_params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");
  Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;

  switch (strat)
  {
    case Inpar::Mortar::ExtendGhosting::redundant_all:
    case Inpar::Mortar::ExtendGhosting::redundant_master:
    {
      melefullmap = Core::LinAlg::AllreduceEMap(*melerowmap_);
      break;
    }
    case Inpar::Mortar::ExtendGhosting::roundrobin:
    {
      melefullmap = melerowmap_;
      break;
    }
    case Inpar::Mortar::ExtendGhosting::binning:
    {
      melefullmap = melecolmap_;
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown strategy to deal with interface ghosting.");
      break;
    }
  }

  // loop over all slave elements on this proc.
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    Core::Elements::Element* element = idiscret_->gElement(selecolmap_->GID(i));
    if (!element) FOUR_C_THROW("Cannot find element with gid %\n", selecolmap_->GID(i));
    auto* mrtrelement = dynamic_cast<Mortar::Element*>(element);
    if (mrtrelement->MinEdgeSize() < lmin) lmin = mrtrelement->MinEdgeSize();
  }

  // loop over all master elements on this proc.
  for (int i = 0; i < melefullmap->NumMyElements(); ++i)
  {
    Core::Elements::Element* element = idiscret_->gElement(melefullmap->GID(i));
    if (!element) FOUR_C_THROW("Cannot find element with gid %\n", melefullmap->GID(i));
    auto* mrtrelement = dynamic_cast<Mortar::Element*>(element);
    if (mrtrelement->MinEdgeSize() < lmin) lmin = mrtrelement->MinEdgeSize();
  }

  // compute DOP inflation length
  enlarge = eps * lmin;

  // define dopnormals
  Core::LinAlg::SerialDenseMatrix dopnormals;
  int kdop = 0;

  if (dim_ == 2)
  {
    kdop = 8;

    // setup normals for DOP
    dopnormals.reshape(4, 3);
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
    dopnormals.reshape(9, 3);
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
    FOUR_C_THROW("Problem dimension must be either 2D or 3D.");

  // decide whether auxiliary positions are used when computing dops
  const bool useauxpos = SearchUseAuxPos();

  // define slave and master slabs
  Core::LinAlg::SerialDenseMatrix sslabs(kdop / 2, 2);
  Core::LinAlg::SerialDenseMatrix mslabs(kdop / 2, 2);

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
    Core::Elements::Element* element = idiscret_->gElement(sgid);
    if (!element) FOUR_C_THROW("Cannot find element with gid %\n", sgid);
    Core::Nodes::Node** node = element->Nodes();
    auto* mrtrnode = dynamic_cast<Node*>(node[0]);
    const double* posnode = mrtrnode->xspatial();

    // calculate slabs initialization
    for (int j = 0; j < kdop / 2; j++)
    {
      //= ax+by+cz=d/sqrt(aa+bb+cc)
      sslabs(j, 0) = sslabs(j, 1) =
          (dopnormals(j, 0) * posnode[0] + dopnormals(j, 1) * posnode[1] +
              dopnormals(j, 2) * posnode[2]) /
          std::sqrt((dopnormals(j, 0) * dopnormals(j, 0)) + (dopnormals(j, 1) * dopnormals(j, 1)) +
                    (dopnormals(j, 2) * dopnormals(j, 2)));
    }

    // for int j=1, because of initialization done before
    for (int j = 1; j < element->num_node(); j++)
    {
      auto* mrtrnode = dynamic_cast<Node*>(node[j]);
      posnode = mrtrnode->xspatial();

      for (int k = 0; k < kdop / 2; k++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        dcurrent = (dopnormals(k, 0) * posnode[0] + dopnormals(k, 1) * posnode[1] +
                       dopnormals(k, 2) * posnode[2]) /
                   std::sqrt((dopnormals(k, 0) * dopnormals(k, 0)) +
                             (dopnormals(k, 1) * dopnormals(k, 1)) +
                             (dopnormals(k, 2) * dopnormals(k, 2)));
        if (dcurrent > sslabs(k, 1)) sslabs(k, 1) = dcurrent;
        if (dcurrent < sslabs(k, 0)) sslabs(k, 0) = dcurrent;
      }
    }

    // add auxiliary positions
    // (last converged positions for all slave nodes)
    if (useauxpos)
    {
      for (int j = 0; j < element->num_node(); j++)
      {
        // get pointer to slave node
        auto* mrtrnode = dynamic_cast<Node*>(node[j]);

        std::array<double, 3> auxpos = {0.0, 0.0, 0.0};
        double scalar = 0.0;
        for (int k = 0; k < dim_; k++)
          scalar += (mrtrnode->X()[k] + mrtrnode->uold()[k] - mrtrnode->xspatial()[k]) *
                    mrtrnode->MoData().n()[k];
        for (int k = 0; k < dim_; k++)
          auxpos[k] = mrtrnode->xspatial()[k] + scalar * mrtrnode->MoData().n()[k];

        for (int l = 0; l < kdop / 2; l++)
        {
          //= ax+by+cz=d/sqrt(aa+bb+cc)
          dcurrent = (dopnormals(l, 0) * auxpos[0] + dopnormals(l, 1) * auxpos[1] +
                         dopnormals(l, 2) * auxpos[2]) /
                     std::sqrt((dopnormals(l, 0) * dopnormals(l, 0)) +
                               (dopnormals(l, 1) * dopnormals(l, 1)) +
                               (dopnormals(l, 2) * dopnormals(l, 2)));
          if (dcurrent > sslabs(l, 1)) sslabs(l, 1) = dcurrent;
          if (dcurrent < sslabs(l, 0)) sslabs(l, 0) = dcurrent;
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
      Core::Elements::Element* element = idiscret_->gElement(mgid);
      if (!element) FOUR_C_THROW("Cannot find element with gid %\n", mgid);
      Core::Nodes::Node** node = element->Nodes();
      auto* mrtrnode = dynamic_cast<Node*>(node[0]);
      const double* posnode = mrtrnode->xspatial();

      // calculate slabs initialization
      for (int k = 0; k < kdop / 2; k++)
      {
        //= ax+by+cz=d/sqrt(aa+bb+cc)
        mslabs(k, 0) = mslabs(k, 1) =
            (dopnormals(k, 0) * posnode[0] + dopnormals(k, 1) * posnode[1] +
                dopnormals(k, 2) * posnode[2]) /
            std::sqrt((dopnormals(k, 0) * dopnormals(k, 0)) +
                      (dopnormals(k, 1) * dopnormals(k, 1)) +
                      (dopnormals(k, 2) * dopnormals(k, 2)));
      }

      // for int k=1, because of initialization done before
      for (int k = 1; k < element->num_node(); k++)
      {
        auto* mrtrnode = dynamic_cast<Node*>(node[k]);
        posnode = mrtrnode->xspatial();

        for (int l = 0; l < kdop / 2; l++)
        {
          //= d=ax+by+cz/sqrt(aa+bb+cc)
          dcurrent = (dopnormals(l, 0) * posnode[0] + dopnormals(l, 1) * posnode[1] +
                         dopnormals(l, 2) * posnode[2]) /
                     std::sqrt((dopnormals(l, 0) * dopnormals(l, 0)) +
                               (dopnormals(l, 1) * dopnormals(l, 1)) +
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
        Core::Elements::Element* element = idiscret_->gElement(sgid);
        auto* selement = dynamic_cast<Mortar::Element*>(element);
        selement->AddSearchElements(mgid);
      }
    }  // for all master elements
  }    // for all slave elements
}

/*----------------------------------------------------------------------*
 |  Search for potentially coupling sl/ma pairs (public)      popp 10/08|
 *----------------------------------------------------------------------*/
bool Mortar::Interface::evaluate_search_binarytree()
{
  binarytree_->evaluate_search();

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlap      popp 11/08|
 *----------------------------------------------------------------------*/
bool Mortar::Interface::MortarCoupling(Mortar::Element* sele, std::vector<Mortar::Element*> mele,
    const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr)
{
  pre_mortar_coupling(sele, mele, mparams_ptr);

  // check if quadratic interpolation is involved
  bool quadratic = false;
  if (sele->IsQuad()) quadratic = true;
  for (auto& m : mele)
    if (m->IsQuad()) quadratic = true;

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
    Mortar::Coupling2dManager(Discret(), Dim(), quadratic, interface_params(), sele, mele)
        .evaluate_coupling(mparams_ptr);
  }
  // ************************************************************** 3D ***
  else if (Dim() == 3)
  {
    // *************************************************** linear 3D ***
    if (!quadratic)
    {
      // create Coupling3dManager and evaluate
      Mortar::Coupling3dManager(Discret(), Dim(), false, interface_params(), sele, mele)
          .evaluate_coupling(mparams_ptr);
    }

    // ************************************************** quadratic 3D ***
    else
    {
      // create Coupling3dQuadManager and evaluate
      Mortar::Coupling3dQuadManager(Discret(), Dim(), false, interface_params(), sele, mele)
          .evaluate_coupling(mparams_ptr);
    }  // quadratic
  }    // 3D
  else
    FOUR_C_THROW("Dimension for Mortar coupling must be either 2D or 3D.");
  // *********************************************************************

  post_mortar_coupling(sele, mele, mparams_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 | Split Mortar::Elements->IntElements for 3D quad. coupling    popp 03/09|
 *----------------------------------------------------------------------*/
bool Mortar::Interface::split_int_elements(
    Mortar::Element& ele, std::vector<Teuchos::RCP<Mortar::IntElement>>& auxele)
{
  // *********************************************************************
  // do splitting for given element
  // *********************************************************** quad9 ***
  if (ele.Shape() == Core::FE::CellType::quad9)
  {
    // split into for quad4 elements
    int numnode = 4;
    Core::FE::CellType dt = Core::FE::CellType::quad4;

    // first integration element
    // containing parent nodes 0,4,8,7
    int nodeids[4] = {0, 0, 0, 0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[8];
    nodeids[3] = ele.NodeIds()[7];

    std::vector<Core::Nodes::Node*> nodes(4);
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
  else if (ele.Shape() == Core::FE::CellType::quad8)
  {
    // split into four tri3 elements and one quad4 element
    int numnodetri = 3;
    int numnodequad = 4;
    Core::FE::CellType dttri = Core::FE::CellType::tri3;
    Core::FE::CellType dtquad = Core::FE::CellType::quad4;

    // first integration element
    // containing parent nodes 0,4,7
    int nodeids[3] = {0, 0, 0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[4];
    nodeids[2] = ele.NodeIds()[7];

    std::vector<Core::Nodes::Node*> nodes(3);
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

    std::vector<Core::Nodes::Node*> nodesquad(4);
    nodesquad[0] = ele.Nodes()[4];
    nodesquad[1] = ele.Nodes()[5];
    nodesquad[2] = ele.Nodes()[6];
    nodesquad[3] = ele.Nodes()[7];

    auxele.push_back(Teuchos::rcp(new IntElement(4, ele.Id(), ele.Owner(), &ele, dtquad,
        numnodequad, nodeidsquad, nodesquad, ele.IsSlave(), false)));
  }

  // ************************************************************ tri6 ***
  else if (ele.Shape() == Core::FE::CellType::tri6)
  {
    // split into four tri3 elements
    int numnode = 3;
    Core::FE::CellType dt = Core::FE::CellType::tri3;

    // first integration element
    // containing parent nodes 0,3,5
    int nodeids[3] = {0, 0, 0};
    nodeids[0] = ele.NodeIds()[0];
    nodeids[1] = ele.NodeIds()[3];
    nodeids[2] = ele.NodeIds()[5];

    std::vector<Core::Nodes::Node*> nodes(3);
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
  else if (ele.Shape() == Core::FE::CellType::quad4)
  {
    // 1:1 conversion to IntElement
    std::vector<Core::Nodes::Node*> nodes(4);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[2];
    nodes[3] = ele.Nodes()[3];

    auxele.push_back(Teuchos::rcp(new IntElement(0, ele.Id(), ele.Owner(), &ele, ele.Shape(),
        ele.num_node(), ele.NodeIds(), nodes, ele.IsSlave(), false)));
  }

  // ************************************************************ tri3 ***
  else if (ele.Shape() == Core::FE::CellType::tri3)
  {
    // 1:1 conversion to IntElement
    std::vector<Core::Nodes::Node*> nodes(3);
    nodes[0] = ele.Nodes()[0];
    nodes[1] = ele.Nodes()[1];
    nodes[2] = ele.Nodes()[2];

    auxele.push_back(Teuchos::rcp(new IntElement(0, ele.Id(), ele.Owner(), &ele, ele.Shape(),
        ele.num_node(), ele.NodeIds(), nodes, ele.IsSlave(), false)));
  }

  // ********************************************************* invalid ***
  else
    FOUR_C_THROW("split_int_elements called for unknown element shape!");

  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  Assemble geometry-dependent lagrange multipliers (global)      popp 05/09|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AssembleLM(Epetra_Vector& zglobal)
{
  // loop over all slave nodes
  for (int j = 0; j < snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    int dim = mrtrnode->NumDof();
    double* lm = mrtrnode->MoData().lm();

    Core::LinAlg::SerialDenseVector lmnode(dim);
    std::vector<int> lmdof(dim);
    std::vector<int> lmowner(dim);

    for (int k = 0; k < dim; ++k)
    {
      lmnode(k) = lm[k];
      lmdof[k] = mrtrnode->Dofs()[k];
      lmowner[k] = mrtrnode->Owner();
    }

    // do assembly
    Core::LinAlg::Assemble(zglobal, lmnode, lmdof, lmowner);
  }
}


/*----------------------------------------------------------------------*
 |  Assemble Mortar D matrix                                  popp 01/08|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AssembleD(Core::LinAlg::SparseMatrix& dglobal)
{
  const bool nonsmooth =
      Core::UTILS::IntegralValue<int>(interface_params(), "NONSMOOTH_GEOMETRIES");
  const bool lagmultlin = (Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(
                               interface_params(), "LM_QUAD") == Inpar::Mortar::lagmult_lin);

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    /**************************************************** D-matrix ******/
    if ((mrtrnode->MoData().GetD()).size() > 0)
    {
      const Core::Gen::Pairedvector<int, double>& dmap = mrtrnode->MoData().GetD();
      int rowsize = mrtrnode->NumDof();

      Core::Gen::Pairedvector<int, double>::const_iterator colcurr;

      for (colcurr = dmap.begin(); colcurr != dmap.end(); ++colcurr)
      {
        double val = colcurr->second;

        Core::Nodes::Node* knode = Discret().gNode(colcurr->first);
        if (!knode) FOUR_C_THROW("node not found");
        auto* kcnode = dynamic_cast<Node*>(knode);
        if (!kcnode) FOUR_C_THROW("node not found");

        for (int j = 0; j < rowsize; ++j)
        {
          int row = mrtrnode->Dofs()[j];
          int col = kcnode->Dofs()[j];

          // do the assembly into global D matrix
          if (!nonsmooth and (shapefcn_ == Inpar::Mortar::shape_dual or
                                 shapefcn_ == Inpar::Mortar::shape_petrovgalerkin))
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
              if (row != col && abs(val) > 1.0e-12) FOUR_C_THROW("D-Matrix is not diagonal!");

              // create an explicitly diagonal d matrix
              if (row == col) dglobal.Assemble(val, row, col);
            }
          }
          else if (nonsmooth or shapefcn_ == Inpar::Mortar::shape_standard)
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
}


/*----------------------------------------------------------------------*
 |  Assemble Mortar M matrix                                  popp 01/08|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AssembleM(Core::LinAlg::SparseMatrix& mglobal)
{
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    /**************************************************** M-matrix ******/
    if ((mrtrnode->MoData().GetM()).size() > 0)
    {
      const std::map<int, double>& mmap = mrtrnode->MoData().GetM();
      int rowsize = mrtrnode->NumDof();

      std::map<int, double>::const_iterator colcurr;

      for (colcurr = mmap.begin(); colcurr != mmap.end(); ++colcurr)
      {
        Core::Nodes::Node* knode = Discret().gNode(colcurr->first);
        if (!knode) FOUR_C_THROW("node not found");
        auto* kcnode = dynamic_cast<Node*>(knode);
        if (!kcnode) FOUR_C_THROW("node not found");

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
      int colsize = static_cast<int>(mmap.size()) * rowsize;

      Core::LinAlg::SerialDenseMatrix Mnode(rowsize, colsize);
      std::vector<int> lmrow(rowsize);
      std::vector<int> lmcol(colsize);
      std::vector<int> lmrowowner(rowsize);
      std::map<int, double>::const_iterator colcurr;
      int k = 0;

      for (colcurr = mmap.begin(); colcurr != mmap.end(); ++colcurr)
      {
        Core::Nodes::Node* knode = Discret().gNode(colcurr->first);
        if (!knode) FOUR_C_THROW("node not found");
        auto* kcnode = dynamic_cast<Node*>(knode);
        if (!kcnode) FOUR_C_THROW("node not found");

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
}


/*----------------------------------------------------------------------*
 |  Assemble Mortar matrices                                 farah 02/16|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AssembleDM(
    Core::LinAlg::SparseMatrix& dglobal, Core::LinAlg::SparseMatrix& mglobal)
{
  // call subroutines:

  // assemble mortar matrix D (slave side)
  AssembleD(dglobal);

  // assemble mortar matrix M (master side)
  AssembleM(mglobal);
}


/*----------------------------------------------------------------------*
 |  Assemble matrix of normals                                popp 10/11|
 *----------------------------------------------------------------------*/
void Mortar::Interface::assemble_normals(Core::LinAlg::SparseMatrix& nglobal)
{
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

    // nodal normal
    double* nodalnormal = mrtrnode->MoData().n();

    // add normal to corresponding row in global matrix
    for (int k = 0; k < mrtrnode->NumDof(); ++k)
      nglobal.Assemble(nodalnormal[k], gid, mrtrnode->Dofs()[k]);
  }
}

/*----------------------------------------------------------------------*
 |  Assemble interface displacement trafo matrices            popp 06/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::AssembleTrafo(Core::LinAlg::SparseMatrix& trafo,
    Core::LinAlg::SparseMatrix& invtrafo, std::set<int>& donebefore)
{
  // check for dual shape functions and quadratic slave elements
  if (shapefcn_ == Inpar::Mortar::shape_standard || !quadslave_)
    FOUR_C_THROW("AssembleTrafo -> you should not be here...");

  // check whether locally linear LM interpolation is used
  const bool lagmultlin = (Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(
                               interface_params(), "LM_QUAD") == Inpar::Mortar::lagmult_lin);

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
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    if (mrtrnode->Owner() != Comm().MyPID()) FOUR_C_THROW("Node ownership inconsistency!");

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
    auto* mrtrele = dynamic_cast<Mortar::Element*>(mrtrnode->Elements()[0]);
    Core::FE::CellType shape = mrtrele->Shape();

    // which discretization type
    switch (shape)
    {
      // line3 contact elements (= tri6||quad8||quad9 discretizations)
      case Core::FE::CellType::line3:
      {
        // modification factor
        if (Inpar::Mortar::LagMultQuad() == Inpar::Mortar::lagmult_lin)
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
      case Core::FE::CellType::tri6:
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
      case Core::FE::CellType::quad8:
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
      case Core::FE::CellType::quad9:
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
        FOUR_C_THROW("Trafo matrix only for line3/tri6/quad8/quad9 contact elements");
        break;
      }
    }  // switch(Shape)

    //********************************************************************
    // CASE 1: CORNER NODES AND CENTER NODE
    //********************************************************************
    if (nt == corner || nt == center)
    {
      // check if processed before
      auto iter = donebefore.find(gid);

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
      auto iter = donebefore.find(gid);

      // if not then assemble trafo matrix block
      if (iter == donebefore.end())
      {
        // add to set of processed nodes
        donebefore.insert(gid);

        // find adjacent corner nodes locally
        int index1 = 0;
        int index2 = 0;
        int hoindex = mrtrele->GetLocalNodeId(gid);
        Core::FE::getCornerNodeIndices(index1, index2, hoindex, shape);

        // find adjacent corner nodes globally
        int gindex1 = mrtrele->NodeIds()[index1];
        int gindex2 = mrtrele->NodeIds()[index2];
        // std::cout << "-> adjacent corner nodes: " << gindex1 << " " << gindex2 << std::endl;
        Core::Nodes::Node* adjnode1 = idiscret_->gNode(gindex1);
        if (!adjnode1) FOUR_C_THROW("Cannot find node with gid %", gindex1);
        auto* adjmrtrnode1 = dynamic_cast<Node*>(adjnode1);
        Core::Nodes::Node* adjnode2 = idiscret_->gNode(gindex2);
        if (!adjnode2) FOUR_C_THROW("Cannot find node with gid %", gindex2);
        auto* adjmrtrnode2 = dynamic_cast<Node*>(adjnode2);

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
      FOUR_C_THROW("Undefined node type (corner, edge, center)");
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
      Core::Nodes::Node* node = idiscret_->gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      auto* mrtrnode = dynamic_cast<Node*>(node);

      if (mrtrnode->Owner() != Comm().MyPID())
        FOUR_C_THROW("AssembleTrafo: Node ownership inconsistency!");

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
      auto* mrtrele = dynamic_cast<Mortar::Element*>(mrtrnode->Elements()[0]);
      Core::FE::CellType shape = mrtrele->Shape();

      // real master nodes are easily identified
      if (!mrtrnode->IsOnBound()) nt = master;
      // former slave node type depends on discretization
      else
      {
        switch (shape)
        {
          case Core::FE::CellType::line3:
          {
            // edge node
            if (mrtrnode->Id() == mrtrele->NodeIds()[2]) nt = slaveedge;

            break;
          }

          // tri6 contact elements (= tet10 discretizations)
          case Core::FE::CellType::tri6:
          {
            // edge nodes
            if (mrtrnode->Id() == mrtrele->NodeIds()[3] ||
                mrtrnode->Id() == mrtrele->NodeIds()[4] || mrtrnode->Id() == mrtrele->NodeIds()[5])
              nt = slaveedge;

            break;
          }

          // quad8 contact elements (= hex20 discretizations)
          case Core::FE::CellType::quad8:
          {
            // edge nodes
            if (mrtrnode->Id() == mrtrele->NodeIds()[4] ||
                mrtrnode->Id() == mrtrele->NodeIds()[5] ||
                mrtrnode->Id() == mrtrele->NodeIds()[6] || mrtrnode->Id() == mrtrele->NodeIds()[7])
              nt = slaveedge;

            break;
          }

          // quad9 contact elements (= hex27 discretizations)
          case Core::FE::CellType::quad9:
          {
            // edge nodes
            if (mrtrnode->Id() == mrtrele->NodeIds()[4] ||
                mrtrnode->Id() == mrtrele->NodeIds()[5] ||
                mrtrnode->Id() == mrtrele->NodeIds()[6] || mrtrnode->Id() == mrtrele->NodeIds()[7])
            {
              nt = slaveedge;
            }

            // center node
            else if (mrtrnode->Id() == mrtrele->NodeIds()[8])
              nt = slavecenter;

            break;
          }

          // other cases
          default:
          {
            FOUR_C_THROW("Trafo matrix only for line3/tri6/quad8/quad9 contact elements");
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
        auto iter = donebefore.find(gid);

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
        auto iter = donebefore.find(gid);

        // if not then assemble trafo matrix block
        if (iter == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(gid);

          // find adjacent corner nodes locally
          int index1 = 0;
          int index2 = 0;
          int hoindex = mrtrele->GetLocalNodeId(gid);
          Core::FE::getCornerNodeIndices(index1, index2, hoindex, shape);

          // find adjacent corner nodes globally
          int gindex1 = mrtrele->NodeIds()[index1];
          int gindex2 = mrtrele->NodeIds()[index2];
          // std::cout << "-> adjacent corner nodes: " << gindex1 << " " << gindex2 << std::endl;
          Core::Nodes::Node* adjnode1 = idiscret_->gNode(gindex1);
          if (!adjnode1) FOUR_C_THROW("Cannot find node with gid %", gindex1);
          auto* adjmrtrnode1 = dynamic_cast<Node*>(adjnode1);
          Core::Nodes::Node* adjnode2 = idiscret_->gNode(gindex2);
          if (!adjnode2) FOUR_C_THROW("Cannot find node with gid %", gindex2);
          auto* adjmrtrnode2 = dynamic_cast<Node*>(adjnode2);

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
        auto iter = donebefore.find(gid);

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
          Core::Nodes::Node* adjnode1 = idiscret_->gNode(gindex1);
          if (!adjnode1) FOUR_C_THROW("Cannot find node with gid %", gindex1);
          auto* adjmrtrnode1 = dynamic_cast<Node*>(adjnode1);
          Core::Nodes::Node* adjnode2 = idiscret_->gNode(gindex2);
          if (!adjnode2) FOUR_C_THROW("Cannot find node with gid %", gindex2);
          auto* adjmrtrnode2 = dynamic_cast<Node*>(adjnode2);
          Core::Nodes::Node* adjnode3 = idiscret_->gNode(gindex3);
          if (!adjnode3) FOUR_C_THROW("Cannot find node with gid %", gindex3);
          auto* adjmrtrnode3 = dynamic_cast<Node*>(adjnode3);
          Core::Nodes::Node* adjnode4 = idiscret_->gNode(gindex4);
          if (!adjnode4) FOUR_C_THROW("Cannot find node with gid %", gindex4);
          auto* adjmrtrnode4 = dynamic_cast<Node*>(adjnode4);

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
        FOUR_C_THROW("Undefined node type (corner, edge, center)");
    }
  }  // end of assembly for locally linear LM interpolation
}

/*----------------------------------------------------------------------*
 |  Detect actual meshtying zone (node by node)               popp 08/10|
 *----------------------------------------------------------------------*/
void Mortar::Interface::detect_tied_slave_nodes(int& founduntied)
{
  //**********************************************************************
  // STEP 1: Build tying info for slave node row map (locally+globally)
  //**********************************************************************
  // global vector for tying info
  Teuchos::RCP<Epetra_Vector> rowtied = Teuchos::rcp(new Epetra_Vector(*snoderowmap_));

  // loop over proc's slave row nodes of the interface for detection
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    // perform detection
    const Core::Gen::Pairedvector<int, double>& dmap = mrtrnode->MoData().GetD();
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
      FOUR_C_THROW("Inconsistency in tied/untied node detection");
    }
  }

  //**********************************************************************
  // STEP 2: Export tying info to slave node column map (globally)
  //**********************************************************************
  // export tying information to standard column map
  Teuchos::RCP<Epetra_Vector> coltied = Teuchos::rcp(new Epetra_Vector(*snodecolmap_));
  Core::LinAlg::Export(*rowtied, *coltied);

  //**********************************************************************
  // STEP 3: Extract tying info for slave node column map (locally)
  //**********************************************************************
  // loop over proc's slave col nodes of the interface for storage
  for (int i = 0; i < snodecolmap_->NumMyElements(); ++i)
  {
    int gid = snodecolmap_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<Node*>(node);

    // check if this node is untied
    if ((*coltied)[i] == 1.0) mrtrnode->SetTiedSlave() = false;
  }
}

/*----------------------------------------------------------------------*
 | create volume ghosting (public)                            ager 06/15|
 *----------------------------------------------------------------------*/
void Mortar::Interface::create_volume_ghosting()
{
  Inpar::CONTACT::Problemtype prb = (Inpar::CONTACT::Problemtype)interface_params().get<int>(
      "PROBTYPE", (int)Inpar::CONTACT::other);

  switch (prb)
  {
    case Inpar::CONTACT::ssi:
    case Inpar::CONTACT::ssi_elch:
    {
      std::vector<std::string> tar_dis;
      tar_dis.emplace_back("structure");
      tar_dis.emplace_back("scatra");
      std::vector<std::pair<int, int>> material_map;
      material_map.emplace_back(std::pair<int, int>(0, 1));
      material_map.emplace_back(std::pair<int, int>(1, 0));

      Mortar::UTILS::create_volume_ghosting(Discret(), tar_dis, material_map);

      // we need to redistribute the scalar field since distribution has changed during setup
      auto structdis = Global::Problem::Instance()->GetDis("structure");
      structdis->RedistributeState(1, "scalarfield");

      break;
    }
    case Inpar::CONTACT::tsi:
    {
      std::vector<std::string> tar_dis;
      tar_dis.emplace_back("structure");
      tar_dis.emplace_back("thermo");
      std::vector<std::pair<int, int>> material_map;
      material_map.emplace_back(std::pair<int, int>(0, 1));
      material_map.emplace_back(std::pair<int, int>(1, 0));

      Mortar::UTILS::create_volume_ghosting(Discret(), tar_dis, material_map);
      break;
    }
    default:
    {
      std::vector<std::string> tar_dis;
      tar_dis.emplace_back("structure");
      Mortar::UTILS::create_volume_ghosting(
          Discret(), tar_dis, std::vector<std::pair<int, int>>(0));

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Mortar::Interface::has_ma_sharing_ref_interface() const
{
  return (interface_data_->get_ma_sharing_ref_interface_ptr() != nullptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Mortar::Interface* Mortar::Interface::get_ma_sharing_ref_interface_ptr() const
{
  return interface_data_->get_ma_sharing_ref_interface_ptr();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::Interface::add_ma_sharing_ref_interface(const Interface* ref_interface)
{
  // avoid non-uniqueness and closed loops
  if (ref_interface->has_ma_sharing_ref_interface())
  {
    if (ref_interface->get_ma_sharing_ref_interface_ptr() == this) return;
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
  if (has_ma_sharing_ref_interface())
  {
    const int size_curr_ref_interface =
        get_ma_sharing_ref_interface_ptr()->MasterRowElements()->NumGlobalElements();
    const int size_new_ref_interface = ref_interface->MasterRowElements()->NumGlobalElements();

    if (size_curr_ref_interface >= size_new_ref_interface) return;
  }

  interface_data_->set_ma_sharing_ref_interface_ptr(ref_interface);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Mortar::Interface::postprocess_quantities(const Teuchos::ParameterList& outputParams)
{
  using Teuchos::RCP;

  // Check if the given parameter list contains all required data to be written to output
  {
    // Vector with names of all required parameter entries
    std::vector<std::string> requiredEntries;
    requiredEntries.emplace_back("step");
    requiredEntries.emplace_back("time");
    requiredEntries.emplace_back("displacement");
    requiredEntries.emplace_back("interface traction");
    requiredEntries.emplace_back("slave forces");
    requiredEntries.emplace_back("master forces");

    check_output_list(outputParams, requiredEntries);
  }

  // Get the discretization writer and get ready for writing
  RCP<Core::IO::DiscretizationWriter> writer = idiscret_->Writer();

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
   * Grab it from there, then use Core::LinAlg::Export() to extract the interface
   * portion from the global displacement vector. Finally, write the interface
   * portion using this interfaces' discretization writer.
   */
  {
    // Get full displacement vector and extract interface displacement
    RCP<const Epetra_Vector> disp = outputParams.get<RCP<const Epetra_Vector>>("displacement");
    RCP<Epetra_Vector> iDisp = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*disp, *iDisp);

    // Write the interface displacement field
    writer->WriteVector("displacement", iDisp, Core::IO::VectorType::dofvector);
  }

  // Write Lagrange multiplier field
  {
    // Get full Lagrange multiplier vector and extract values of this interface
    RCP<const Epetra_Vector> lagMult =
        outputParams.get<RCP<const Epetra_Vector>>("interface traction");
    RCP<Epetra_Vector> iLagMult = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*lagMult, *iLagMult);

    // Write this interface's Lagrange multiplier field
    writer->WriteVector("interfacetraction", iLagMult, Core::IO::VectorType::dofvector);
  }

  // Write nodal forces of slave side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> slaveforces =
        outputParams.get<RCP<const Epetra_Vector>>("slave forces");
    RCP<Epetra_Vector> forces = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*slaveforces, *forces);

    // Write to output
    writer->WriteVector("slaveforces", forces, Core::IO::VectorType::dofvector);
  }

  // Write nodal forces of master side
  {
    // Get nodal forces
    RCP<const Epetra_Vector> masterforces =
        outputParams.get<RCP<const Epetra_Vector>>("master forces");
    RCP<Epetra_Vector> forces = Core::LinAlg::CreateVector(*idiscret_->dof_row_map());
    Core::LinAlg::Export(*masterforces, *forces);

    // Write to output
    writer->WriteVector("masterforces", forces, Core::IO::VectorType::dofvector);
  }


  // Nodes: node-based vector with '0' at slave nodes and '1' at master nodes
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*mnoderowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> nodeRowMap = Core::LinAlg::MergeMap(snoderowmap_, mnoderowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = Core::LinAlg::CreateVector(*nodeRowMap, true);
    Core::LinAlg::Export(*masterVec, *masterSlaveVec);

    writer->WriteVector("slavemasternodes", masterSlaveVec, Core::IO::VectorType::nodevector);
  }

  // Elements: element-based vector with '0' at slave elements and '1' at master elements
  {
    RCP<Epetra_Vector> masterVec = Teuchos::rcp(new Epetra_Vector(*melerowmap_));
    masterVec->PutScalar(1.0);

    RCP<const Epetra_Map> eleRowMap = Core::LinAlg::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> masterSlaveVec = Core::LinAlg::CreateVector(*eleRowMap, true);
    Core::LinAlg::Export(*masterVec, *masterSlaveVec);

    writer->WriteVector("slavemasterelements", masterSlaveVec, Core::IO::VectorType::elementvector);
  }

  // Write element owners
  {
    RCP<const Epetra_Map> eleRowMap = Core::LinAlg::MergeMap(selerowmap_, melerowmap_, false);
    RCP<Epetra_Vector> owner = Core::LinAlg::CreateVector(*eleRowMap);

    for (int i = 0; i < idiscret_->ElementRowMap()->NumMyElements(); ++i)
      (*owner)[i] = idiscret_->lRowElement(i)->Owner();

    writer->WriteVector("Owner", owner, Core::IO::VectorType::elementvector);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Mortar::Interface::check_output_list(
    const Teuchos::ParameterList& outParams, const std::vector<std::string>& requiredEntries) const
{
  // Check for each required parameter entry if it exists
  for (const auto& requiredEntry : requiredEntries)
  {
    if (not outParams.isParameter(requiredEntry))
    {
      FOUR_C_THROW("Parameter list is missing the required entry '%s'.", (requiredEntry).c_str());
      return false;
    }
  }

  // We only make it to here, if all checks passed. So it's safe to return 'true'.
  return true;
}

FOUR_C_NAMESPACE_CLOSE
