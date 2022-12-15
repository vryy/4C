/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality to read nodes from file

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef DRT_MESHREADER_H
#define DRT_MESHREADER_H

#include "inpar.H"

namespace DRT
{
  namespace INPUT
  {
    class ElementReader;
    class DomainReader;
    class DatFileReader;

    /*----------------------------------------------------------------------*/
    /*!
      \brief helper class to read all nodes shared by many discretizations

      The interface class to the reading subsystem.
     */
    /*----------------------------------------------------------------------*/
    class MeshReader
    {
     public:
      /// construct a reader that reads a given section
      MeshReader(const DatFileReader& reader, std::string sectionname);

      /// add an element reader for each discretization
      /*!
        Each discretization needs its own ElementReader. These readers
        have to be registered at the MeshReader.

        \param er (i) a reader of one discretization that uses (a fraction of) our nodes
       */
      void AddElementReader(Teuchos::RCP<ElementReader> er) { element_readers_.push_back(er); }

      /*!
       * \brief Adds the selected element/domain/file reader to this nodereader
       *
       * \param dis            [in] This discretization will be passed on
       * \param reader         [in] This reader will be passed on
       * \param sectionname    [in] This will be passed on element/domain readers only (not used for
       * file reader) \param elementtypes   [in] This set of types will be passed on \param
       * geometrysource [in] selects which reader will be created \param geofilepath    [in] path to
       * the file for the file reader (not used for the others)
       */
      void AddAdvancedReader(Teuchos::RCP<Discretization> dis,
          const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
          const std::set<std::string>& elementtypes, const INPAR::GeometryType geometrysource,
          const std::string* geofilepath);

      /*!
       * \brief Adds the selected element/domain/file reader to this nodereader
       *
       * This is a version without specific elementtypes. It just calls the full
       * version with a dummy set.
       *
       * \param dis            [in] This discretization will be passed on
       * \param reader         [in] This reader will be passed on
       * \param sectionname    [in] This will be passed on element/domain readers only (not used for
       * file reader) \param geometrysource [in] selects which reader will be created \param
       * geofilepath    [in] path to the file for the file reader (not used for the others)
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
          partitioning, building the node row and column maps.

        - Reading and distributing all nodes. Each node gets assigned to
          its discretization.

        - Finalizing the discretizations.

        Actually most of the work gets done by the ElementReader. The
        reading of both elements and nodes happens in blocks on processor
        0. After each block read the discretizations are redistributed.

       */
      void ReadAndPartition();

     private:
      /// find the discretization that contains the given node
      std::vector<Teuchos::RCP<DRT::Discretization>> FindDisNode(int nodeid);

      /*!
      \brief Read pre-generated mesh from dat-file

      \param[in/out] node_id_offset Global offset to start node numbering (based on already existing
      nodes)
      */
      void ReadMeshFromDatFile(int& node_id_offset);

      /*!
      \brief Create inline mesh

      Ask the DomainReader to process input data and create a mesh at runtime.

      \param[in/out] node_id_offset Global offset to start node numbering (based on already existing
      nodes)
      */
      void CreateInlineMesh(int& node_id_offset);

      /*!
      \brief Check for sufficient number of nodes and throw if too few

      We need at least as many nodes as there are processors. If there are not enough nodes, throw
      an error.
      */
      void ThrowIfNotEnoughNodes(int max_node_id) const;

      /// the main dat file reader
      const DatFileReader& reader_;

      /// my comm
      Teuchos::RCP<Epetra_Comm> comm_;

      /// my section to read
      std::string sectionname_;

      /// my element readers
      std::vector<Teuchos::RCP<ElementReader>> element_readers_;

      /// my domain readers
      std::vector<Teuchos::RCP<DomainReader>> domain_readers_;
    };

  }  // namespace INPUT
}  // namespace DRT

#endif
