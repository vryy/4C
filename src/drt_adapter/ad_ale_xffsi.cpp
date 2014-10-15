/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_xffsi.cpp
 <pre>
       Maintainer: Matthias Mayr
       mayr@mhpc.mw.tum.de
       089 - 289-10362
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_xffsi.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleXFFsiWrapper::AleXFFsiWrapper(Teuchos::RCP<Ale> ale)
  : AleFsiWrapper(ale)
{
  return;
}

//! ToDo (mayr) move this to XFluidFluid adapter
void ADAPTER::AleXFFsiWrapper::SolveAleXFluidFluidFSI()
{
  dserror("XFFSI Not implemented yet!");
}

