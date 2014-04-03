/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.cpp

\brief Fluid field adapter for fpsi

Can only be used in conjunction with #FluidImplicitTimeInt

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fluid_fpsi.H"
#include "ad_fld_fluid_fsi.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

/* constructor */
ADAPTER::FluidFPSI::FluidFPSI(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output,
    bool isale,
    bool dirichletcond)
: FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond),
  fpsiinterface_(Teuchos::rcp(new FLD::UTILS::MapExtractor()))
{
  return;
} // constructor


/* initialization */
void ADAPTER::FluidFPSI::Init()
{
  // call base class init
  FluidFSI::Init();

  fpsiinterface_->Setup(*dis_,true);

  return;
}


