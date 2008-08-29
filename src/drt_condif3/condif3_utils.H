/*!
\file condif3_utils.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONDIF3_UTILS_H
#define CONDIF3_UTILS_H


#include "condif3.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_discret.H"


namespace SCATRA
{
  //! check for higher order derivatives for shape functions
  bool isHigherOrderElement
  (
   const DRT::Element::DiscretizationType& distype ///< distype of the element
  );


  //! get optimal gaussrule for discretization type
  DRT::UTILS::GaussRule3D getOptimalGaussrule
  (
   const DRT::Element::DiscretizationType& distype ///< distype of the element
  );


} // namespace SCATRA

#endif
#endif  // #ifdef CCADISCRET
