/*----------------------------------------------------------------------*/
/*! \file

\brief Read element sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_IO_ELEMENTREADER_HPP
#define BACI_IO_ELEMENTREADER_HPP

#include "baci_config.hpp"

#include "baci_io_inputreader.hpp"

#include <Epetra_Map.h>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace IO
{
  /*----------------------------------------------------------------------*/
  /*!
    \brief helper class to read the elements of a discretization

    Together with MeshReader this class constitutes a (almost) parallel
    and efficient reading mechanism for discretizations from dat files.

    We face the following problem:

    - There are elements and nodes. One set of elements per
      discretization. One set of nodes with the nodes from all
      discretizations.

    - Elements and nodes have ids. These are unique but otherwise
      arbitrary.

    - We cannot afford to read all elements or all nodes on one
      processor.

    - Only processor 0 can actually read the (ascii) input file

    - We do not want to setup (that is read) elements more than once.

    The idea is to read blocks of elements and nodes and distribute them
    to different processors at first. Afterwards a reasonable
    distribution can be calculated and the discretizations can be
    redistributed. How this work is done is a mere
    technicality. However, we need to be able to use the discretization
    in a partially constructed state. In particular we need to read, add
    and distribute elements even if nodes not yet are available.
   */
  /*----------------------------------------------------------------------*/
  class ElementReader
  {
   public:
    /*!
    \brief Construct element reader for a given field that reads a given section

    Create empty discretization and append it to given field.

    \param dis (i) the new discretization
    \param comm (i) our communicator
    \param sectionname (i) the section that contains the element lines
    */
    ElementReader(Teuchos::RCP<DRT::Discretization> dis, const INPUT::DatFileReader& reader,
        std::string sectionname);

    /*!
    \brief Construct element reader for a given field that reads a given section

    Create empty discretization and append it to given field.

    \param dis (i) the new discretization
    \param comm (i) our communicator
    \param sectionname (i) the section that contains the element lines
    \param elementtype (i) element type name to read in this discretization
    */
    ElementReader(Teuchos::RCP<DRT::Discretization> dis, const INPUT::DatFileReader& reader,
        std::string sectionname, std::string elementtype);

    /*!
    \brief Construct element reader for a given field that reads a given section

    Create empty discretization and append it to given field.

    \param dis (i) the new discretization
    \param comm (i) our communicator
    \param sectionname (i) the section that contains the element lines
    \param elementtypes (i) element type names to read in this discretization
    */
    ElementReader(Teuchos::RCP<DRT::Discretization> dis, const INPUT::DatFileReader& reader,
        std::string sectionname, const std::set<std::string>& elementtypes);

    //! Destructor
    virtual ~ElementReader() = default;
    // Return list of unique nodes
    std::set<int> GetUniqueNodes() const { return nodes_; }

    /// give the discretization this reader fills
    Teuchos::RCP<DRT::Discretization> GetDis() const { return dis_; }

    /// Return the list of row elements
    Teuchos::RCP<Epetra_Map> GetRowElements() const { return roweles_; }

    /*! Read elements and partition the node graph

    - read global ids of elements of this discretization
      (this is one fully redundant vector for elements)
    - determine a preliminary element distribution. The fully redundant
      vector is trashed after the construction.
    - define blocksizes for blocks of elements we read (not necessarily
      the same as it was used to construct the map --- we may have a
      smaller blocksize here).
    - read elements of this discretization and distribute according
      to a linear map. While reading, remember node gids and assemble
      them into a second fully redundant vector (mapping node id->gid).
      In addition, we keep them in a fully redundant set (required by
      node reader). Construct reverse lookup from gids to node ids.
      Again, this is a global, fully redundant map!
    - define preliminary linear distributed nodal row map
    - determine adjacency array (i.e. the infos for the node graph)
      using the nodal row distribution and a round robin communication
      of element connectivity information.
      Use adjacency array to build an initial Crsgraph on the linear map.
    - do partitioning using parmetis
      Results are distributed to other procs using two global vectors!
    - build final nodal row map, export graph to the new map
    */
    virtual void ReadAndDistribute();

    /*!
    \brief Tell whether the given node belongs to us

    \note This is based on the redundant nodes_ set and only available on processor 0.
    */
    bool HasNode(const int nodeid) const { return nodes_.find(nodeid) != nodes_.end(); }

   private:
    /// Get the overall number of elements and their corresponding global IDs
    std::pair<int, std::vector<int>> GetElementSizeAndIDs() const;

    /// Read the file and get element information, distribute them to each processor
    void GetAndDistributeElements(const int nblock, const int bsize);

    /// discretization name
    std::string name_;

    /// the main dat file reader
    const INPUT::DatFileReader& reader_;

    /// my comm
    Teuchos::RCP<Epetra_Comm> comm_;

    /// my section to read
    std::string sectionname_;

    /*!
    \brief All global node ids of a discretization on processor 0

    This is a redundant set of all node numbers. But it is only valid
    on processor 0. We need it to easily figure out to which
    discretization a node belongs.
    */
    std::set<int> nodes_;

    /// my discretization
    Teuchos::RCP<DRT::Discretization> dis_;

    /// node row map
    Teuchos::RCP<Epetra_Map> rownodes_;

    /// node col map
    Teuchos::RCP<Epetra_Map> colnodes_;

    /// element row map
    Teuchos::RCP<Epetra_Map> roweles_;

    /// element col map
    Teuchos::RCP<Epetra_Map> coleles_;

    /// element type names to read
    std::set<std::string> elementtypes_;
  };
}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif
