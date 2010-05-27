/*!
\file xfluid3_utils.H

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFLUID3_UTILS_H
#define XFLUID3_UTILS_H

#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_xfem/xfem_enums.H"
#include "../drt_geometry/intersection_service.H"
#include <Epetra_SerialDenseSolver.h>

namespace XFLUID
{

template <size_t numparam>
struct EnrViscs2
{
  LINALG::Matrix<numparam,1> xx;
  LINALG::Matrix<numparam,1> xy;
  LINALG::Matrix<numparam,1> xz;
  LINALG::Matrix<numparam,1> yx;
  LINALG::Matrix<numparam,1> yy;
  LINALG::Matrix<numparam,1> yz;
  LINALG::Matrix<numparam,1> zx;
  LINALG::Matrix<numparam,1> zy;
  LINALG::Matrix<numparam,1> zz;
};

template <size_t numparam>
struct DnDx
{
  LINALG::Matrix<numparam,1> dnxdx;
  LINALG::Matrix<numparam,1> dnxdy;
  LINALG::Matrix<numparam,1> dnxdz;
  LINALG::Matrix<numparam,1> dnydx;
  LINALG::Matrix<numparam,1> dnydy;
  LINALG::Matrix<numparam,1> dnydz;
  LINALG::Matrix<numparam,1> dnzdx;
  LINALG::Matrix<numparam,1> dnzdy;
  LINALG::Matrix<numparam,1> dnzdz;
};

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

} // namespace XFLUID

#endif
#endif  // #ifdef CCADISCRET
