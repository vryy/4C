/*!
\file condif2_utils.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONDIF2_UTILS_H
#define CONDIF2_UTILS_H


#include "condif2.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_discret.H"


namespace SCATRA
{
  //! check for higher order derivatives for shape functions
  bool is2DHigherOrderElement
  (
   const DRT::Element::DiscretizationType& distype ///< distype of the element
  );


  //! get optimal gaussrule for discretization type
  DRT::UTILS::GaussRule2D get2DOptimalGaussrule
  (
   const DRT::Element::DiscretizationType& distype ///< distype of the element
  );


} // namespace SCATRA

#endif
#endif  // #ifdef CCADISCRET
