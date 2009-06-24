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
#include "../drt_xfem/xfem_enums.H"
#include "../drt_geometry/intersection_service.H"
#include <Epetra_SerialDenseSolver.h>

namespace XDIFF
{

  /// hold arrays for all - potentially enriched - approximation functions
  /// and their derivatives with respect to the physical coordinates
  ///
  /// Let P be the approximation function and N the shape function.
  /// Then, for standard FEM, P = N
  ///
  /// Two examples:
  /// 
  /// u^h = SUM_I P^I u^I
  ///  u^I are the free parameters (nodal unknowns in the standard FEM) and
  ///  normal and enriched parameters combined in the XFEM 
  ///
  /// | P1 |   | N1     |
  /// | P2 |   | N2     |
  /// | P3 | = | N1*psi |
  /// | P4 |   | N2*psi |
  ///   ...
  /// similar
  /// du^h/dx = SUM_I dP^I/dx u^I 
  ///   P is the approximation function
  ///
  /// | P1/dx |   | N1/dx                   |
  /// | P2/dx |   | N2/dx                   |
  /// | P3/dx | = | dN1/dx*psi + N1*dpsi/dx |
  /// | P4/dx |   | dN2/dx*psi + N2*dpsi/dx |
  ///   ...
  ///   Note the chain rule above!
  template <size_t numparam>
  struct ApproxFunc
  {
    LINALG::Matrix<numparam,1>  d0;
    LINALG::Matrix<numparam,1>  dx;
    LINALG::Matrix<numparam,1>  dy;
    LINALG::Matrix<numparam,1>  dz;
//    LINALG::Matrix<numparam,1>  dxdx;
//    LINALG::Matrix<numparam,1>  dxdy;
//    LINALG::Matrix<numparam,1>  dxdz;
//    LINALG::Matrix<numparam,1>  dydx;
//    LINALG::Matrix<numparam,1>  dydy;
//    LINALG::Matrix<numparam,1>  dydz;
//    LINALG::Matrix<numparam,1>  dzdx;
//    LINALG::Matrix<numparam,1>  dzdy;
//    LINALG::Matrix<numparam,1>  dzdz;
  };
  
template<class M1, class M2>
void fill_tau(
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
        rule = DRT::UTILS::intrule_tet_24point;
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
          rule = DRT::UTILS::intrule_tet_24point;
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
