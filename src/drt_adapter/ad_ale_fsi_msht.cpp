/*--------------------------------------------------------------------------*/
/*!
\file ad_ale_fsi_msht.cpp

\brief FSI Wrapper for the ALE time integration with internal mesh tying or mesh sliding interface

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fsi_msht.H"

#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFsiMshtWrapper::AleFsiMshtWrapper(Teuchos::RCP<Ale> ale) : AleFsiWrapper(ale)
{
  // create the FSI interface
  fsiinterface_ = Teuchos::rcp(new ALE::UTILS::FsiMapExtractor);
  fsiinterface_->Setup(*Discretization());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::FsiMapExtractor> ADAPTER::AleFsiMshtWrapper::FsiInterface() const
{
  return fsiinterface_;
}
