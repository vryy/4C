/*!
\file xdiff3_utils.H

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XDIFF3_UTILS_H
#define XDIFF3_UTILS_H

#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_xfem/xfem_utils.H"
#include "../drt_geometry/intersection_service.H"
#include <Epetra_SerialDenseSolver.h>

namespace XDIFF
{

template<class M1, class M2>
void fill_flux(
    const size_t         numparamtau,
    const M1&            shp_tau,
    const M2&            eflux,
    LINALG::Matrix<3,1>& flux
)
{
  flux.Clear();
  for (size_t i = 0; i < numparamtau; ++i)
  {
    const double fac = shp_tau(i);
    flux(0,0) += fac*eflux(0,i);
    flux(1,0) += fac*eflux(1,i);
    flux(2,0) += fac*eflux(2,i);
  }
}

/*!
 * \brief calculate the body force from the nodes condition
 *
 *  get the body force in the nodes of the element (private) gammi 04/07
 *  the Neumann condition associated with the nodes is stored in the
 *  array edeadng only if all nodes have a VolumeNeumann condition
 *
 */
template <DRT::Element::DiscretizationType DISTYPE>
LINALG::SerialDenseMatrix BodyForce(DRT::Node** nodes, const double pseudotime)
{
  dserror("not adapted to xfem (, yet)!!!");

  std::vector<DRT::Condition*> myneumcond;
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  // check whether all nodes have a unique VolumeNeumann condition
  size_t nodecount = 0;
  for (size_t inode=0;inode<numnode;++inode)
  {
    nodes[inode]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      nodecount++;
    }
  }

  LINALG::SerialDenseMatrix edeadng_(3,numnode);
  if (nodecount == numnode)
  {
    // find out whether we will use a (pseudo-)time curve
    const std::vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a (pseudo-)timecurve
    {
      // factor for the intermediate step
      if(pseudotime >= 0.0)
      {
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(pseudotime);
      }
      else
      {
    // do not compute an "alternative" curvefac here since a negative pseudotime value
    // indicates an error.
        dserror("Negative pseudotime value in body force calculation: time = %f",pseudotime);
        //curvefac = DRT::Problem::Instance()->Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a (pseudo-)timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for (size_t jnode=0; jnode<numnode; ++jnode)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const std::vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const std::vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(size_t isd=0;isd<3;++isd)
      {
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_.Zero();
  }

  return edeadng_;
}

//! interpolate from nodal vector array to integration point vector using the shape function
template <class M, class VS>
double interpolateScalarFieldToIntPoint(
    const M&  eleVectorField,       ///< array with nodal vector values
    const VS& shp,                  ///< array with nodal shape function
    const int numparam              ///< number of parameters
    )
{
  double v = 0.0;
  for (int iparam = 0; iparam < numparam; ++iparam)
      v += eleVectorField(iparam)*shp(iparam);
  return v;
}

} // namespace XFLUID

#endif
#endif  // #ifdef CCADISCRET
