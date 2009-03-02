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
#include "../drt_xfem/xfem_enums.H"
#include "../drt_geometry/intersection_service.H"
#include <Epetra_SerialDenseSolver.h>

namespace XFLUID
{

template<class M1, class M2>
void fill_tau(const int numparamtau, const M1&  shp_tau, const M2& etau,
    LINALG::Matrix<3,3>& tau)
{
    tau.Clear();
    for (int i = 0; i < numparamtau; ++i)
    {
        const double fac = shp_tau(i);
        tau(0,0) += fac*etau(0,i);
        tau(1,1) += fac*etau(1,i);
        tau(2,2) += fac*etau(2,i);
        tau(0,1) += fac*etau(3,i);
        tau(0,2) += fac*etau(4,i);
        tau(1,2) += fac*etau(5,i);
    }
    tau(1,0) = tau(0,1);
    tau(2,0) = tau(0,2);
    tau(2,1) = tau(1,2);
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
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<numnode;++inode)
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
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(pseudotime);
      }
      else
      {
    // do not compute an "alternative" curvefac here since a negative pseudotime value
    // indicates an error.
        dserror("Negative pseudotime value in body force calculation: time = %f",pseudotime);
        //curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a (pseudo-)timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for (int jnode=0; jnode<numnode; ++jnode)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const std::vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const std::vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<3;++isd)
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


/*!
 \brief get optimal gaussrule for given discretization type
        template version
 \return enum with optimal gausrule for this element
 */
inline DRT::UTILS::GaussRule3D getOptimalVelocityGaussrule(DRT::Element::DiscretizationType  distype)
{
    DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    switch (distype)
    {
    case DRT::Element::hex8:
        rule = DRT::UTILS::intrule_hex_8point;
        break;
    case DRT::Element::hex20: case DRT::Element::hex27:
        rule = DRT::UTILS::intrule_hex_27point;
        break;
    case DRT::Element::tet4:
        rule = DRT::UTILS::intrule_tet_4point;
        break;
    case DRT::Element::tet10:
        rule = DRT::UTILS::intrule_tet_11point;
        break;
    default:
        dserror("unknown distype for gaussrule initialization");
  }
  return rule;
}

/*!
 \brief get optimal gaussrule for given discretization type

 \return enum with optimal gausrule for this element
 */
template<DRT::Element::DiscretizationType XELEDISTYPE, class M>
DRT::UTILS::GaussRule3D getXFEMGaussrule(
    const DRT::Element*                     xele,
    const M&                                xyze_element,
    const bool                              intersected,
    const DRT::Element::DiscretizationType  celldistype
    )
{
  DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
  if (not intersected)
  {
    rule = getOptimalVelocityGaussrule(celldistype);
  }
  else
  {
    const GEO::EleGeoType xeleGeoType = GEO::checkGeoTypeT<XELEDISTYPE>(xele,xyze_element);
    
    if (celldistype == DRT::Element::tet4 or celldistype == DRT::Element::tet10)
    {
      switch (XELEDISTYPE)
      {
        case DRT::Element::hex8:
        {
          if (xeleGeoType == GEO::CARTESIAN) // constant jacobian
            rule = DRT::UTILS::intrule_tet_15point;
          else
            rule = DRT::UTILS::intrule_tet_24point;
          break;
        }
        case DRT::Element::hex20:
        {
          if (xeleGeoType == GEO::CARTESIAN) // constant jacobian
            rule = DRT::UTILS::intrule_tet_24point;
          else
            rule = DRT::UTILS::intrule_tet_45point;
          break;
        }
        case DRT::Element::hex27:
        {
          if (xeleGeoType == GEO::CARTESIAN) // constant jacobian
            rule = DRT::UTILS::intrule_tet_45point;
          else
            rule = DRT::UTILS::intrule_tet_45point;
          break;
        }
        case DRT::Element::tet4:
        {
          rule = DRT::UTILS::intrule_tet_4point;
          break;
        }
        case DRT::Element::tet10:
        {
          rule = DRT::UTILS::intrule_tet_11point;
          break;
        }
        default:
          dserror("no rule defined for this xele");
      }
    }
    else
    {
      dserror("gaussrule only for tetrahedral integration cells defined. Please add additional rules here!");
    }
  }

  return rule;
}

/*!
  tell, whether second derivatives have to be calculated
 */
template <DRT::Element::DiscretizationType  DISTYPE>
bool secondDerivativesAvailable()
{
  bool hoel = true;
  switch (DISTYPE)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
  case DRT::Element::pyramid5:
    hoel = true;
    break;
  case DRT::Element::tet4:
    hoel = false;
    break;
  default:
    dserror("DISTYPE unknown!");
  }
  return hoel;
}

//! decide, whether second derivatives can occur
template <DRT::Element::DiscretizationType  DISTYPE>
bool isHigherOrderElement()
{
  bool hoel = true;
  switch (DISTYPE)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge15:
    hoel = true;
    break;
  case DRT::Element::tet4:
  case DRT::Element::wedge6:
  case DRT::Element::pyramid5:
    hoel = false;
    break;
  default:
    dserror("DISTYPE unknown!");
  }
  return hoel;
}


//! interpolate from nodal vector array to integration point vector using the shape function
template <class M, class VS>
LINALG::Matrix<3,1> interpolateVectorFieldToIntPoint(
    const M&  eleVectorField,       ///< array with nodal vector values
    const VS& shp,                  ///< array with nodal shape function
    const int numparam              ///< number of parameters
    )
{
  const int nsd = 3;
  LINALG::Matrix<nsd,1> v(true);
  for (int iparam = 0; iparam < numparam; ++iparam)
    for (int isd = 0; isd < nsd; ++isd)
      v(isd) += eleVectorField(isd,iparam)*shp(iparam);
  return v;
}

} // namespace XFLUID

#endif
#endif  // #ifdef CCADISCRET
