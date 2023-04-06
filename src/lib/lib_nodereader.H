/*----------------------------------------------------------------------*/
/*! \file

\brief Read node sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef LIB_NODEREADER_H
#define LIB_NODEREADER_H

#include "lib_inputreader.H"
#include "lib_elementreader.H"

#include <Epetra_Map.h>

namespace DRT
{
  class Discretization;

  namespace INPUT
  {
    class DatFileReader;

    class NodeReader
    {
     public:
      /*!
      \brief Construct node reader for a given field that reads a given section

      Create empty discretization and append it to given field.

      \param comm (i) our communicator
      \param sectionname (i) the section that contains the node lines
      */
      NodeReader(const DatFileReader& reader, std::string sectionname);

      void Read(const std::vector<ElementReader>& element_readers, int& max_node_id);

      void ThrowIfNotEnoughNodes(int max_node_id) const;

      std::vector<Teuchos::RCP<DRT::Discretization>> FindDisNode(
          const std::vector<ElementReader>& element_readers, int global_node_id);

     private:
      /// the main dat file reader
      const DatFileReader& reader_;

      /// my comm
      Teuchos::RCP<Epetra_Comm> comm_;

      /// my section to read
      std::string sectionname_;
    };
  }  // namespace INPUT
}  // namespace DRT

#endif
