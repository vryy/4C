/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality to read a mesh from a dat file

\level 0

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_MESHREADER_H
#define BACI_LIB_MESHREADER_H

#include "baci_inpar.H"
#include "baci_io_inputreader.H"
#include "baci_lib_domainreader.H"
#include "baci_lib_elementreader.H"
#include "baci_lib_nodereader.H"

#include <Epetra_CrsGraph.h>

BACI_NAMESPACE_OPEN

namespace DRT::INPUT
{
  /*!
    \brief helper class to read a mesh

    This is an interface for handling node, element and domain readers
    and to set up a discretization from a dat file which is FillComplete().
   */
  class MeshReader
  {
   public:
    /// construct a mesh reader
    MeshReader(const Teuchos::RCP<Epetra_Comm> comm);

    /// add an element reader for each discretization
    /*!
      Each discretization needs its own ElementReader. These readers
      have to be registered at the MeshReader.

      \param er (i) a reader of one discretization that uses (a fraction of) our nodes
     */
    void AddElementReader(const ElementReader& er) { element_readers_.emplace_back(er); }

    /// add an node reader
    /*!
      The reader has to be registered at the MeshReader.

    \param nr a node reader of the overall problem
    */
    void SetNodeReader(Teuchos::RCP<NodeReader> nr) { node_reader_ = nr; }

    /*!
     * \brief Adds the selected reader to this meshreader
     *
     * \param dis            [in] This discretization will be passed on
     * \param reader         [in] This reader will be passed on
     * \param sectionname    [in] This will be passed on element/domain readers only (not used for
     *                            file reader)
     * \param elementtypes   [in] This set of types will be passed on
     * \param geometrysource [in] selects which reader will be created
     * \param geofilepath    [in] path to the file for the file reader (not used for the others)
     */
    void AddAdvancedReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
        const std::set<std::string>& elementtypes, const INPAR::GeometryType geometrysource,
        const std::string* geofilepath);

    /*!
     * \brief Adds the selected reader to this meshreader
     *
     * This is a version without specific elementtypes. It just calls the full
     * version with a dummy set.
     *
     * \param dis            [in] This discretization will be passed on
     * \param reader         [in] This reader will be passed on
     * \param sectionname    [in] This will be passed on element/domain readers only (not used for
     *                            file reader)
     * \param geometrysource [in] selects which reader will be created
     * \param geofilepath    [in] path to the file for the file reader (not used for the others)
     */
    void AddAdvancedReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
        const INPAR::GeometryType geometrysource, const std::string* geofilepath);

    /// do the actual reading
    /*!
      This method contains the whole machinery. The reading consists of
      three steps:

      - Reading and distributing elements using each registered
        ElementReader. This includes creating the connectivity graph,
        building the node row and column maps.

      - Reading and distributing all nodes. Each node gets assigned to
        its discretization.

      - Finalizing the discretizations.

      Actually most of the work gets done by the ElementReader. The
      reading of both elements and nodes happens in blocks on processor
      0. After each block read the discretizations are redistributed.

     */
    void ReadAndPartition();

   private:
    /*!
    \brief Read pre-generated mesh from dat-file and generate related FillCompleted()
    discretizations

    \param[in/out] max_node_id Maximum node id in a given discretization. To be used as global
                               offset to start node numbering (based on already existing nodes)
    */
    void ReadMeshFromDatFile(int& max_node_id);

    /*!
    \brief Rebalance discretizations built in ReadMeshFromDatFile()
    */
    void Rebalance();

    /*!
    \brief Create inline mesh

    Ask the DomainReader to process input data and create a box shaped mesh at runtime without
    having to read or process nodes and elements from the input file.

    \param[in/out] max_node_id Maximum node id in a given discretization. To be used as global
                               offset to start node numbering (based on already existing nodes)
    */
    void CreateInlineMesh(int& max_node_id);

    /// my comm
    Teuchos::RCP<Epetra_Comm> comm_;

    //! graphs of each discretization
    std::vector<Teuchos::RCP<const Epetra_CrsGraph>> graph_;

    /// my element readers
    std::vector<ElementReader> element_readers_;

    /// my domain readers
    std::vector<DomainReader> domain_readers_;

    /// my node reader
    Teuchos::RCP<NodeReader> node_reader_;
  };
}  // namespace DRT::INPUT

BACI_NAMESPACE_CLOSE

#endif
