/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality for reading nodes

\level 0

*/
/*---------------------------------------------------------------------*/

#include "lib_discret.H"
#include "lib_meshreader.H"
#include "lib_domainreader.H"
#include "lib_elementreader.H"
#include "lib_globalproblem.H"
#include "lib_inputreader.H"
#include "nurbs_discret_control_point.H"
#include "immersed_problem_immersed_node.H"
#include "fiber_node.H"
#include "io_pstream.H"

#include <string>

namespace DRT::INPUT
{
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  MeshReader::MeshReader(const DRT::INPUT::DatFileReader& reader) : comm_(reader.Comm()) {}


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
      const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
      const std::set<std::string>& elementtypes, const INPAR::GeometryType geometrysource,
      const std::string* geofilepath)
  {
    switch (geometrysource)
    {
      case INPAR::geometry_full:
      {
        std::string fullsectionname("--" + sectionname + " ELEMENTS");
        Teuchos::RCP<ElementReader> er =
            Teuchos::rcp(new DRT::INPUT::ElementReader(dis, reader, fullsectionname, elementtypes));
        element_readers_.push_back(er);
        break;
      }
      case INPAR::geometry_box:
      {
        std::string fullsectionname("--" + sectionname + " DOMAIN");
        Teuchos::RCP<DomainReader> dr =
            Teuchos::rcp(new DRT::INPUT::DomainReader(dis, reader, fullsectionname));
        domain_readers_.push_back(dr);
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
  void MeshReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
      const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
      const INPAR::GeometryType geometrysource, const std::string* geofilepath)
  {
    std::set<std::string> dummy;
    AddAdvancedReader(dis, reader, sectionname, dummy, geometrysource, geofilepath);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::ReadAndPartition()
  {
    // We need to track the max global node ID to offset node numbering and for sanity checks
    int max_node_id = 0;

    ReadMeshFromDatFile(max_node_id);
    CreateInlineMesh(max_node_id);

    // last check if there are enough nodes
    node_reader->ThrowIfNotEnoughNodes(max_node_id);
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::ReadMeshFromDatFile(int& max_node_id)
  {
    TEUCHOS_FUNC_TIME_MONITOR("MeshReader::ReadMeshFromDatFile");

    // read and partition element information
    for (const auto& element_reader : element_readers_) element_reader->ReadAndPartition();

    // TODO: here we should place the partitioning -> extract from element reader

    // read nodes based on the element information
    node_reader->Read(element_readers_, max_node_id);

    // last thing to do here is to produce nodal ghosting/overlap
    for (const auto& element_reader : element_readers_)
    {
      element_reader->MyDis()->ExportColumnNodes(*element_reader->MyColNodes());
      element_reader->Complete();
    }
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void MeshReader::CreateInlineMesh(int& max_node_id)
  {
    for (const auto& domain_reader : domain_readers_)
    {
      // communicate node offset to all procs
      int local_max_node_id = max_node_id;
      comm_->MaxAll(&local_max_node_id, &max_node_id, 1);

      domain_reader->CreatePartitionedMesh(max_node_id);
      domain_reader->Complete();
      max_node_id = domain_reader->MyDis()->NodeRowMap()->MaxAllGID() + 1;
    }
  }
}  // namespace DRT::INPUT
