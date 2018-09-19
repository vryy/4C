/*----------------------------------------------------------------------*/
/*!
\file drt_elementreader.H

\brief Read element sections of dat files.

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef DRT_ELEMENTREADER_H
#define DRT_ELEMENTREADER_H

#include <set>
#include <Epetra_Map.h>

#include "drt_inputreader.H"

namespace DRT
{
  class Discretization;

  namespace INPUT
  {
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
    class ElementReader
    {
      friend class NodeReader;

     public:
      /// construct element reader for a given field that reads a given section
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param comm (i) our communicator
        \param sectionname (i) the section that contains the element lines
       */
      ElementReader(
          Teuchos::RCP<Discretization> dis, const DatFileReader& reader, std::string sectionname);

      /// construct element reader for a given field that reads a given section
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param comm (i) our communicator
        \param sectionname (i) the section that contains the element lines
        \param elementtype (i) element type name to read in this discretization
       */
      ElementReader(Teuchos::RCP<Discretization> dis, const DatFileReader& reader,
          std::string sectionname, std::string elementtype);

      /// construct element reader for a given field that reads a given section
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param comm (i) our communicator
        \param sectionname (i) the section that contains the element lines
        \param elementtypes (i) element type names to read in this discretization
       */
      ElementReader(Teuchos::RCP<Discretization> dis, const DatFileReader& reader,
          std::string sectionname, const std::set<std::string>& elementtypes);

      //! Destructor
      virtual ~ElementReader() { ; }

      /// give the discretization this reader fills
      Teuchos::RCP<Discretization> MyDis() const { return dis_; }

     protected:
      /// read elements, partition node graph
      virtual void Partition();

      /// tell whether the given node belongs to us
      /*!
        \note This is based on the redundant nodes_ set and only available
        on processor 0.
       */
      bool HasNode(int nodeid) { return nodes_.find(nodeid) != nodes_.end(); }

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

      /// all global node ids of a discretization on processor 0
      /*!
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



    /*----------------------------------------------------------------------*/
    /*!
    Particle discretizations do not need elements. Hence, the functionality
    of the element reader must be specified for particles separately

    \author ghamm
    \date 11/12
    */
    /*----------------------------------------------------------------------*/
    class OldParticleReader : public ElementReader
    {
     public:
      /// construct particle reader for a given field
      /*!
        Create empty discretization and append it to given field.

        \param dis (i) the new discretization
        \param comm (i) our communicator
       */
      OldParticleReader(Teuchos::RCP<Discretization> dis, const DatFileReader& reader);

      //! Destructor
      virtual ~OldParticleReader() { ; }

     protected:
      /// partition particles
      virtual void Partition();
    };

  }  // namespace INPUT
}  // namespace DRT

#endif
