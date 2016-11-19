/*!---------------------------------------------------------------------*
\file DatFileReader_stub.H

\brief A stub for reading from fake datfiles

\level 3

\maintainer Karl-Robert Wichmann
*-----------------------------------------------------------------------*/

#ifndef DATFILEREADER_STUB_H
#define DATFILEREADER_STUB_H

#include "src/drt_lib/drt_inputreader.H"

namespace DRT
{
namespace INPUT
{

class DatFileReaderStub : public DatFileReader
{
 public:
  DatFileReaderStub();
  Teuchos::RCP<Epetra_Comm> Comm() const;
};

}
}

#endif
