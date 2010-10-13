#if 0
#ifdef STKADAPTIVE

#ifndef STK_GMSH_H
#define STK_GMSH_H

#include <string>

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Discretization;
}

namespace STK
{
  void DumpGmsh( DRT::Discretization& dis, std::string name, Epetra_Vector & f, int step );
}

#endif
#endif
#endif
