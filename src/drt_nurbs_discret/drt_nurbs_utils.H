/*----------------------------------------------------------------------*/
/*!
\file drt_nurbs_utils.H

\brief A service method allowing the application of initial conditions
       for nurbs discretisations.

Since nurbs shape functions are not interpolating, it is not as
straightforward to apply initial conditions to the degrees of freedom.
(dofs are always associated with control points, i.e. the location
associated with the 'node'=control point is not the physical location
and the value at the control point is not the prescribed value at this
position since dofs associated with neighbouring control points influence
the function value as well)


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_NURBS_UTILS_H
#define DRT_NURBS_UTILS_H

#include "../drt_nurbs_discret/drt_nurbs_discret.H"

namespace DRT
{

namespace NURBS
{
/*----------------------------------------------------------------------*/
/*!
\brief A service method for accessing knotvector and weights for
       an isogeometric element

</pre>

\param dis         (i) the discretisation
\param ele         (i) a pointer to the element
\param myknots     (o) knot vector (to be filled)
\param weights     (o) knot vector (to be filled)

\date 12/10
*/
template <class WG>
bool GetMyNurbsKnotsAndWeights(
    DRT::Discretization&                      dis,
    DRT::Element*                             ele,
    std::vector<Epetra_SerialDenseVector>&    myknots,
    WG&                                       weights
)
{
  // try to cast dis to a nurbs discretisation --- if possible, proceed
  // with setting initial conditions. Otherwise we have a problem.
  DRT::NURBS::NurbsDiscretization* nurbsdis
  = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));
  if(nurbsdis==NULL)
    dserror("Received discretization which is not Nurbs!");

  // get local knot vector entries and check for zero sized elements
  const bool zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());

  // if we have a zero sized element due to a interpolated
  // point --- exit here and tell the outside world about that
  if(zero_size)
  {
    return(zero_size);
  }
  // you are still here? So get the node weights for the nurbs element as well
  DRT::Node** nodes = ele->Nodes();
  const int nen = ele->NumNode();
  for (int inode=0; inode<nen; inode++)
  {
    DRT::NURBS::ControlPoint* cp
    = dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
    weights(inode) = cp->W();
  }

  // goodbye
  return zero_size;
}


//! determine whether this is an NURBS element or not
inline bool IsNurbs(DRT::Element::DiscretizationType distype)
{
  switch(distype)
  {
  case DRT::Element::nurbs8:
  case DRT::Element::nurbs27:
  case DRT::Element::nurbs4:
  case DRT::Element::nurbs9:
  case DRT::Element::nurbs2:
  case DRT::Element::nurbs3:
  {
    return true;
    break;
  }
  default:
    return false;
  }
};

} // namespace NURBS

} //namespace DRT

#endif // ifndef DRT_NURBS_UTILS_H
#endif // #ifdef CCADISCRET
