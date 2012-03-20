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
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_inpar/inpar_fluid.H"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

FLD::TimInt::TimInt(Teuchos::RCP<Teuchos::ParameterList> params)
  :time_    (0.0),
   step_    (0),
   dta_     (params->get<double> ("time step size")),
   timealgo_(DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params, "time int algo"))
{}

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

