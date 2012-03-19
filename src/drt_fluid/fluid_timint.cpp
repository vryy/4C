/*!----------------------------------------------------------------------
\file fluid_timint.cpp
\brief Baseclass for all fluid time integrations

<pre>
Maintainers: Volker Gravemeier & Andreas Ehrl
             {vgravem,ehrl}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15245/-252
</pre>

*----------------------------------------------------------------------*/

#include "fluid_timint.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FLD::TimInt::~TimInt()
{
}

Teuchos::RCP<const Epetra_Map> FLD::TimInt::DofRowMap(unsigned nds)
{
  return Teuchos::rcp(Discretization()->DofRowMap(nds), false);
}

FLD::UTILS::FluidXFluidMapExtractor FLD::TimInt::XFluidFluidMapExtractor()
{
  dserror("Not implemented in the base class, may be overridden by a subclass.");
  FLD::UTILS::FluidXFluidMapExtractor compilerBeQuiet;
  return compilerBeQuiet;
}

