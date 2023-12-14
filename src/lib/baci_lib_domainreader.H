/*----------------------------------------------------------------------*/
/*! \file

\brief Read domain sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/



#ifndef BACI_LIB_DOMAINREADER_H
#define BACI_LIB_DOMAINREADER_H

#include "baci_io_inputreader.H"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <memory>
#include <set>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace GRIDGENERATOR
  {
    class RectangularCuboidInputs;
  }

  namespace INPUT
  {
    class LineDefinition;

    /*----------------------------------------------------------------------*/
    /*!
      \brief helper class to create a rectangular cuboid discretization of (currently) hex or wedge
      elements based on a few input parameters
     */
    /*----------------------------------------------------------------------*/
    class DomainReader
    {
      friend class MeshReader;

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

      //! Destructor
      virtual ~DomainReader() = default;
      /// give the discretization this reader fills
      Teuchos::RCP<Discretization> MyDis() const { return dis_; }

     private:
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
       * \param nodeoffset [in] Node GID of first newly created node.
       * \note Node GIDs may not overlap
       */
      void CreatePartitionedMesh(int nodeGIdOfFirstNewNode) const;

      /*! \brief read input parameters from input file
         \return class holding all input parameters for rectangular cuboid domain
       */
      DRT::GRIDGENERATOR::RectangularCuboidInputs ReadRectangularCuboidInputData() const;

      /// finalize reading. FillComplete(false,false,false), that is, do not
      /// initialize elements. This is done later after reading boundary conditions.
      void Complete() const;

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
    };

  }  // namespace INPUT
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif
