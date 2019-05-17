/*----------------------------------------------------------------------*/
/*!

\brief Read domain sections of dat files.

\level 0

\maintainer  Christoph Ager

*/
/*----------------------------------------------------------------------*/



#ifndef DRT_DOMAINREADER_H
#define DRT_DOMAINREADER_H

#include <set>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include "drt_inputreader.H"

namespace DRT
{
  class Discretization;

  namespace INPUT
  {
    class LineDefinition;

    /*----------------------------------------------------------------------*/
    /*!
      \brief helper class to read the elements of a discretization

      Together with NodeReader this class constitutes a (almost) parallel
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

      \author u.kue
      \date 05/07
     */
    /*----------------------------------------------------------------------*/
    class DomainReader
    {
      friend class NodeReader;

     public:
      /// construct element reader for a given field that reads a given section
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param reader (i) the dat file reader
        \param sectionname (i) the section that contains the element lines
       */
      DomainReader(
          Teuchos::RCP<Discretization> dis, const DatFileReader& reader, std::string sectionname);

      /// construct element reader for a given field that reads a given section
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param reader (i) the dat file reader
        \param sectionname (i) the section that contains the element lines
        \param elementtype (i) element type name to read in this discretization
       */
      DomainReader(Teuchos::RCP<Discretization> dis, const DatFileReader& reader,
          std::string sectionname, std::string elementtype);

      /// construct element reader for a given field that reads a given section
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param reader (i) the dat file reader
        \param sectionname (i) the section that contains the element lines
        \param elementtypes (i) element type names to read in this discretization
       */
      DomainReader(Teuchos::RCP<Discretization> dis, const DatFileReader& reader,
          std::string sectionname, const std::set<std::string>& elementtypes);

      //! Destructor
      virtual ~DomainReader() { ; }

      /// give the discretization this reader fills
      Teuchos::RCP<Discretization> MyDis() const { return dis_; }

     protected:
      /*! \brief generate elements, partition node graph, create nodes
       *
       * This method processes DOMAIN SECTIONs. It
       * 1) generates elements
       * 2) partitions the element graph
       * 3) exports the col elements
       * 4) generates the node graph
       * 5) creates the nodes
       * 6) exports the col nodes
       *
       * \param nodeoffset [in/out] Node GIDs may not overlap, so this must/will be MaxGID+1
       */
      virtual void Partition(int* nodeoffset);

      /// finialize reading. FillComplete(false,false,false), that is, do not
      /// initialize elements. This is done later after reading boundary conditions.
      void Complete();

      /// discretization name
      std::string name_;

      /// the main dat file reader
      const DatFileReader& reader_;

      /// my comm
      Teuchos::RCP<Epetra_Comm> comm_;

      /// my section to read
      std::string sectionname_;

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

      /// lower boundary of the domain to be created
      double lower_bound_[3];

      /// upper boundary of the domain to be created
      double upper_bound_[3];

      /// intervals of the domain to be created
      int interval_[3];

      /// rotation angles of the box geometry
      double rotation_angle_[3];

      /// element type of the domain to be created
      std::string elementtype_;

      /// discretization type of the domain to be created
      std::string distype_;

      /// further arguments to the elements to be created
      std::string elearguments_;

     private:
      /// create HEX type elements for the partition
      void CreateHexElement(
          int eleid, int* nodeoffset, int myrank, DRT::INPUT::LineDefinition* linedef);

      /*!
         Create WEDGE type elements for the partition.
         For even eleids -> create 1st part of HEX equivalent, odd -> 2nd part of HEX equivalent.
         Wedges aligned in z-direction.
       */
      void CreateWedgeElement(
          int eleid, int* nodeoffset, int myrank, DRT::INPUT::LineDefinition* linedef);
    };


    /*----------------------------------------------------------------------*/
    /*!
    Particle discretizations do not need elements. Hence, the functionality
    of the domain reader must be specified for particles separately

    \author ghamm
    \date 11/14
    */
    /*----------------------------------------------------------------------*/
    class ParticleDomainReader : public DomainReader
    {
     public:
      /// construct particle reader for a given field
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param reader (i) the dat file reader
       */
      ParticleDomainReader(
          Teuchos::RCP<Discretization> dis, const DatFileReader& reader, std::string sectionname);

      //! Destructor
      virtual ~ParticleDomainReader() { ; }

     protected:
      /// partition particles
      virtual void Partition(int* nodeoffset);
    };

  }  // namespace INPUT
}  // namespace DRT
#endif
