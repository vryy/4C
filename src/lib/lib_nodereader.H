/*----------------------------------------------------------------------*/
/*! \file

\brief Read node sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef NODEREADER_H
#define NODEREADER_H

#include "lib_inputreader.H"

#include "Epetra_Map.h"

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
      \brief Construct element reader for a given field that reads a given section

      Create empty discretization and append it to given field.

      \param comm (i) our communicator
      \param sectionname (i) the section that contains the node lines
      */
      NodeReader(const DatFileReader& reader, std::string sectionname);

      void Read();

      void ThrowIfNotEnoughNodes(int max_node_id) const;

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
