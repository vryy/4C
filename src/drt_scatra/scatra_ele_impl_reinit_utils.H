/*!
\file scatra_ele_impl_reinit_utils.H

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_IMPL_REINIT_UTILS_H
#define SCATRA_ELE_IMPL_REINIT_UTILS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"


namespace SCATRA
{
  /*!
  \brief Decide, whether second derivatives are needed  (template version)
   *  In characteristic Galerkin methods (for reinitialization of level-set) we use mixed 2nd order derivatives which do not vanish */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct Use2ndDerivsReinitialization {};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::hex8>    {static const bool use = true;}; // modified by Benedikt!!!!
  template<> struct Use2ndDerivsReinitialization<DRT::Element::tet4>    {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::wedge6>  {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::pyramid5>{static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::nurbs8>  {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::quad4>   {static const bool use = true;}; // modified by Benedikt!!!!
  template<> struct Use2ndDerivsReinitialization<DRT::Element::nurbs4>  {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::tri3>    {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::line2>   {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::nurbs2>  {static const bool use = false;};
  template<> struct Use2ndDerivsReinitialization<DRT::Element::hex20>   {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::hex27>   {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::nurbs27> {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::tet10>   {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::quad8>   {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::quad9>   {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::nurbs9>  {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::tri6>    {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::line3>   {static const bool use = true; };
  template<> struct Use2ndDerivsReinitialization<DRT::Element::nurbs3>  {static const bool use = true;};

} // namespace SCATRA

#endif
