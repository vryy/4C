/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality for reading nodes

\level 0

*/
/*---------------------------------------------------------------------*/

#include "4C_io_meshreader.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_domainreader.hpp"
#include "4C_io_elementreader.hpp"
#include "4C_io_inputreader.hpp"
#include "4C_io_nodereader.hpp"
#include "4C_rebalance.hpp"
#include "4C_rebalance_graph_based.hpp"
#include "4C_rebalance_print.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Core::IO::MeshReader::MeshReader(
    Core::IO::DatFileReader& reader, std::string node_section_name, MeshReaderParameters parameters)
    : comm_(reader.get_comm()),
      reader_(reader),
      node_section_name_(std::move(node_section_name)),
      parameters_(std::move(parameters))
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::add_advanced_reader(Teuchos::RCP<Core::FE::Discretization> dis,
    const Core::IO::DatFileReader& reader, const std::string& sectionname,
    const Core::IO::GeometryType geometrysource, const std::string* geofilepath)
{
  std::set<std::string> elementtypes;
  switch (geometrysource)
  {
    case Core::IO::geometry_full:
    {
      std::string fullsectionname("--" + sectionname + " ELEMENTS");
      ElementReader er = ElementReader(dis, reader, fullsectionname, elementtypes);
      element_readers_.emplace_back(er);
      break;
    }
    case Core::IO::geometry_box:
    {
      std::string fullsectionname("--" + sectionname + " DOMAIN");
      DomainReader dr = DomainReader(dis, reader, fullsectionname);
      domain_readers_.emplace_back(dr);
      break;
    }
    case Core::IO::geometry_file:
    {
      FOUR_C_THROW("Unfortunately not yet implemented, but feel free ...");
      break;
    }
    default:
      FOUR_C_THROW("Unknown geometry source");
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::read_and_partition()
{
  // We need to track the max global node ID to offset node numbering and for sanity checks
  int max_node_id = 0;

  graph_.resize(element_readers_.size());

  read_mesh_from_dat_file(max_node_id);
  rebalance();
  create_inline_mesh(max_node_id);

  // last check if there are enough nodes
  {
    int local_max_node_id = max_node_id;
    comm_->MaxAll(&local_max_node_id, &max_node_id, 1);

    if (max_node_id < comm_->NumProc() && reader_.excluded_section_length(node_section_name_) != 0)
      FOUR_C_THROW("Bad idea: Simulation with %d procs for problem with %d nodes", comm_->NumProc(),
          max_node_id);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::read_mesh_from_dat_file(int& max_node_id)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::IO::MeshReader::read_mesh_from_dat_file");

  // read element information
  for (auto& element_reader : element_readers_) element_reader.read_and_distribute();

  // read nodes based on the element information
  read_nodes(reader_, node_section_name_, element_readers_, max_node_id);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::rebalance()
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::IO::MeshReader::Rebalance");

  // do the real partitioning and distribute maps
  for (size_t i = 0; i < element_readers_.size(); i++)
  {
    // global node ids --- this will be a fully redundant vector!
    int numnodes = static_cast<int>(element_readers_[i].get_unique_nodes().size());
    comm_->Broadcast(&numnodes, 1, 0);

    const auto discret = element_readers_[i].get_dis();

    // We want to be able to read empty fields. If we have such a beast
    // just skip the building of the node  graph and do a proper initialization
    if (numnodes)
      graph_[i] = Core::Rebalance::build_graph(discret, element_readers_[i].get_row_elements());
    else
      graph_[i] = Teuchos::null;

    // create partitioning parameters
    const double imbalance_tol =
        parameters_.mesh_paritioning_parameters.get<double>("IMBALANCE_TOL");

    Teuchos::RCP<Teuchos::ParameterList> rebalanceParams =
        Teuchos::rcp(new Teuchos::ParameterList());
    rebalanceParams->set<std::string>("imbalance tol", std::to_string(imbalance_tol));

    const auto rebalanceMethod = Teuchos::getIntegralValue<Core::Rebalance::RebalanceType>(
        parameters_.mesh_paritioning_parameters, "METHOD");

    Teuchos::RCP<Epetra_Map> rowmap, colmap;

    if (!graph_[i].is_null())
    {
      switch (rebalanceMethod)
      {
        case Core::Rebalance::RebalanceType::hypergraph:
        {
          rebalanceParams->set("partitioning method", "HYPERGRAPH");

          // here we can reuse the graph, which was calculated before, this saves us some time
          std::tie(rowmap, colmap) =
              Core::Rebalance::rebalance_node_maps(graph_[i], *rebalanceParams);

          break;
        }
        case Core::Rebalance::RebalanceType::recursive_coordinate_bisection:
        {
          rebalanceParams->set("partitioning method", "RCB");

          // here we can reuse the graph, which was calculated before, this saves us some time and
          // in addition calculate geometric information based on the coordinates of the
          // discretization
          rowmap = Teuchos::rcp(new Epetra_Map(-1, graph_[i]->RowMap().NumMyElements(),
              graph_[i]->RowMap().MyGlobalElements(), 0, *comm_));
          colmap = Teuchos::rcp(new Epetra_Map(-1, graph_[i]->ColMap().NumMyElements(),
              graph_[i]->ColMap().MyGlobalElements(), 0, *comm_));

          discret->redistribute(*rowmap, *colmap, false, false, false);

          Teuchos::RCP<Epetra_MultiVector> coordinates = discret->build_node_coordinates();

          std::tie(rowmap, colmap) = Core::Rebalance::rebalance_node_maps(
              graph_[i], *rebalanceParams, Teuchos::null, Teuchos::null, coordinates);

          break;
        }
        case Core::Rebalance::RebalanceType::monolithic:
        {
          rebalanceParams->set("partitioning method", "HYPERGRAPH");

          rowmap = Teuchos::rcp(new Epetra_Map(-1, graph_[i]->RowMap().NumMyElements(),
              graph_[i]->RowMap().MyGlobalElements(), 0, *comm_));
          colmap = Teuchos::rcp(new Epetra_Map(-1, graph_[i]->ColMap().NumMyElements(),
              graph_[i]->ColMap().MyGlobalElements(), 0, *comm_));

          discret->redistribute(*rowmap, *colmap, true, true, false);

          Teuchos::RCP<const Epetra_CrsGraph> enriched_graph =
              Core::Rebalance::build_monolithic_node_graph(*discret,
                  Core::GeometricSearch::GeometricSearchParams(
                      parameters_.geometric_search_parameters, parameters_.io_parameters));

          std::tie(rowmap, colmap) =
              Core::Rebalance::rebalance_node_maps(enriched_graph, *rebalanceParams);

          break;
        }
        default:
          FOUR_C_THROW("Appropriate partitioning has to be set!");
      }
    }
    else
    {
      rowmap = colmap = Teuchos::rcp(new Epetra_Map(-1, 0, nullptr, 0, *comm_));
    }

    discret->redistribute(*rowmap, *colmap, false, false, false);

    Core::Rebalance::UTILS::print_parallel_distribution(*discret);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::IO::MeshReader::create_inline_mesh(int& max_node_id)
{
  for (const auto& domain_reader : domain_readers_)
  {
    // communicate node offset to all procs
    int local_max_node_id = max_node_id;
    comm_->MaxAll(&local_max_node_id, &max_node_id, 1);

    domain_reader.create_partitioned_mesh(max_node_id);
    domain_reader.complete();
    max_node_id = domain_reader.my_dis()->node_row_map()->MaxAllGID() + 1;
  }
}

FOUR_C_NAMESPACE_CLOSE
