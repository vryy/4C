/*----------------------------------------------------------------------*/
/*!
\file drt_meshfree_utils.H

\brief service methods for a given meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_MESHFREE_UTILS_H
#define DRT_MESHFREE_UTILS_H

#include "../drt_lib/drt_element.H"

namespace DRT
{

namespace MESHFREE
{

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Determines whether a given element is a MeshfreeElement or not.
   *
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
inline bool IsMeshfree(const enum DRT::Element::DiscretizationType& distype)
{
  dserror("discretisation type for meshfree cell still to be defined!");
  return false;
//  switch(distype)
//  {
//  case DRT::Element::nurbs8:
//  {
//    return true;
//    break;
//  }
//  default:
//    return false;
//  }
};

} // namespace MESHFREE

} // namespace DRT

#endif // ifndef DRT_MESHFREE_UTILS_H
