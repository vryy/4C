/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality for reading nodes

\level 0

*/
/*---------------------------------------------------------------------*/

#include "baci_lib_meshreader.H"

#include "baci_io_pstream.H"
#include "baci_lib_discret.H"
#include "baci_lib_domainreader.H"
#include "baci_lib_elementreader.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_inputreader.H"
#include "baci_rebalance.H"
#include "baci_rebalance_utils.H"

#include <string>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::INPUT::MeshReader::MeshReader(const Teuchos::RCP<Epetra_Comm> comm) : comm_(comm) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::MeshReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
    const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
    const std::set<std::string>& elementtypes, const INPAR::GeometryType geometrysource,
    const std::string* geofilepath)
{
  switch (geometrysource)
  {
    case INPAR::geometry_full:
    {
      std::string fullsectionname("--" + sectionname + " ELEMENTS");
      ElementReader er = DRT::INPUT::ElementReader(dis, reader, fullsectionname, elementtypes);
      element_readers_.emplace_back(er);
      break;
    }
    case INPAR::geometry_box:
    {
      std::string fullsectionname("--" + sectionname + " DOMAIN");
      DomainReader dr = DRT::INPUT::DomainReader(dis, reader, fullsectionname);
      domain_readers_.emplace_back(dr);
      break;
    }
    case INPAR::geometry_file:
    {
      dserror("Unfortunately not yet implemented, but feel free ...");
      break;
    }
    default:
      dserror("Unknown geometry source");
      break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::MeshReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
    const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
    const INPAR::GeometryType geometrysource, const std::string* geofilepath)
{
  std::set<std::string> dummy;
  AddAdvancedReader(dis, reader, sectionname, dummy, geometrysource, geofilepath);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::MeshReader::ReadAndPartition()
{
  // We need to track the max global node ID to offset node numbering and for sanity checks
  int max_node_id = 0;

  graph_.resize(element_readers_.size());

  ReadMeshFromDatFile(max_node_id);
  Rebalance();
  CreateInlineMesh(max_node_id);

  // last check if there are enough nodes
  node_reader_->ThrowIfNotEnoughNodes(max_node_id);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::MeshReader::ReadMeshFromDatFile(int& max_node_id)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::INPUT::MeshReader::ReadMeshFromDatFile");

  // read element information
  for (auto& element_reader : element_readers_) element_reader.ReadAndDistribute();

  // read nodes based on the element information
  node_reader_->Read(element_readers_, max_node_id);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::MeshReader::Rebalance()
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::INPUT::MeshReader::Rebalance");

  // do the real partitioning and distribute maps
  for (size_t i = 0; i < element_readers_.size(); i++)
  {
    // global node ids --- this will be a fully redundant vector!
    int numnodes = static_cast<int>(element_readers_[i].GetUniqueNodes().size());
    comm_->Broadcast(&numnodes, 1, 0);

    // We want to be able to read empty fields. If we have such a beast
    // just skip the building of the node  graph and do a proper initialization
    if (numnodes)
      graph_[i] =
          REBALANCE::BuildGraph(element_readers_[i].GetDis(), element_readers_[i].GetRowElements());
    else
      graph_[i] = Teuchos::null;

    // create partitioning parameters
    const double imbalance_tol =
        DRT::Problem::Instance()->MeshPartitioningParams().get<double>("IMBALANCE_TOL");

    Teuchos::RCP<Teuchos::ParameterList> rebalanceParams =
        Teuchos::rcp(new Teuchos::ParameterList());
    rebalanceParams->set<std::string>("imbalance tol", std::to_string(imbalance_tol));

    const auto rebalanceMethod = Teuchos::getIntegralValue<INPAR::REBALANCE::RebalanceType>(
        DRT::Problem::Instance()->MeshPartitioningParams(), "METHOD");

    Teuchos::RCP<Epetra_Map> rowmap, colmap;

    switch (rebalanceMethod)
    {
      case INPAR::REBALANCE::RebalanceType::hypergraph:
      {
        rebalanceParams->set("partitioning method", "HYPERGRAPH");
        if (!graph_[i].is_null())
        {
          // here we can reuse the graph, which was calculated before, this saves us some time
          std::tie(rowmap, colmap) = REBALANCE::RebalanceNodeMaps(graph_[i], *rebalanceParams);
        }
        else
          rowmap = colmap = Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, *comm_));
        break;
      }
      case INPAR::REBALANCE::RebalanceType::recursive_coordinate_bisection:
      {
        rebalanceParams->set("partitioning method", "RCB");
        if (!graph_[i].is_null())
        {
          // here we can reuse the graph, which was calculated before, this saves us some time and
          // in addition calculate geometric information based on the coordinates of the
          // discretization
          rowmap = Teuchos::rcp(new Epetra_Map(-1, graph_[i]->RowMap().NumMyElements(),
              graph_[i]->RowMap().MyGlobalElements(), 0, *comm_));
          colmap = Teuchos::rcp(new Epetra_Map(-1, graph_[i]->ColMap().NumMyElements(),
              graph_[i]->ColMap().MyGlobalElements(), 0, *comm_));
          element_readers_[i].GetDis()->Redistribute(*rowmap, *colmap, false, false, false);
          Teuchos::RCP<Epetra_MultiVector> coordinates =
              element_readers_[i].GetDis()->BuildNodeCoordinates();
          std::tie(rowmap, colmap) = REBALANCE::RebalanceNodeMaps(
              graph_[i], *rebalanceParams, Teuchos::null, Teuchos::null, coordinates);
        }
        else
          rowmap = colmap = Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, *comm_));
        break;
      }
      default:
      {
        dserror("Appropriate partitioning has to be set!");
        break;
      }
    }

    element_readers_[i].GetDis()->Redistribute(*rowmap, *colmap, false, false, false);

    REBALANCE::UTILS::PrintParallelDistribution(*element_readers_[i].GetDis());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::INPUT::MeshReader::CreateInlineMesh(int& max_node_id)
{
  for (const auto& domain_reader : domain_readers_)
  {
    // communicate node offset to all procs
    int local_max_node_id = max_node_id;
    comm_->MaxAll(&local_max_node_id, &max_node_id, 1);

    domain_reader.CreatePartitionedMesh(max_node_id);
    domain_reader.Complete();
    max_node_id = domain_reader.MyDis()->NodeRowMap()->MaxAllGID() + 1;
  }
}
