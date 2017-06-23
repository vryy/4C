/*!---------------------------------------------------------------------*
\file DatFileReader_stub.cpp

\brief A stub for reading from fake datfiles

\level 3

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/

#include "DatFileReader_stub.H"
#include <Epetra_SerialComm.h>

namespace DRT
{
namespace INPUT
{

DatFileReaderStub::DatFileReaderStub()
{}

Teuchos::RCP<Epetra_Comm> DatFileReaderStub::Comm() const
{
  return Teuchos::rcp(new Epetra_SerialComm);
}

}
}
