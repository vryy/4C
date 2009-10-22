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


template<class M1, class M2>
void fill_tau(
    const size_t         numparamtau,
    const M1&            shp_tau,
    const M2&            etau,
    LINALG::Matrix<3,3>& tau
    )
{
    tau.Clear();
    for (size_t i = 0; i < numparamtau; ++i)
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

/*!
 * \brief calculate the body force from the nodes condition
 *
 *  get the body force in the nodes of the element (private) gammi 04/07
 *  the Neumann condition associated with the nodes is stored in the
 *  array edeadng only if all nodes have a VolumeNeumann condition
 *
 *  modified for two-phase flows
 *
 *  \author rasthofer 08/09
 *
 */
template <DRT::Element::DiscretizationType DISTYPE>
LINALG::SerialDenseMatrix BodyForceTwoPhaseFlow(
     const DRT::Element* ele,
     const double time)
{
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    LINALG::SerialDenseMatrix edeadng(3,numnode);
    edeadng.Zero();
    vector<DRT::Condition*> myneumcond;

    // check whether all nodes have a unique VolumeNeumann condition
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

    if (myneumcond.size()>1)
      dserror("more than one VolumeNeumann cond on one node");

    if (myneumcond.size()==1)
    {
      // find out whether we will use a time curve
      const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
      int curvenum = -1;

      if (curve) curvenum = (*curve)[0];

      // initialisation
      double curvefac    = 0.0;

      if (curvenum >= 0) // yes, we have a timecurve
      {
        // time factor for the intermediate step
        if(time >= 0.0)
        {
          curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
        }
        else
        {
        // do not compute an "alternative" curvefac here since a negative time value
        // indicates an error.
          dserror("Negative time value in body force calculation: time = %f",time);
          //curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
        }
      }
      else // we do not have a timecurve --- timefactors are constant equal 1
      {
        curvefac = 1.0;
      }

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );
      const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

      // factor given by spatial function
      double functionfac = 1.0;
      int functnum = -1;

      // set this condition to the edeadng array
      for(int isd=0;isd<3;isd++)
      {
        // get factor given by spatial function
        if (functions) functnum = (*functions)[isd];
        else functnum = -1;

        double num = (*onoff)[isd]*(*val)[isd]*curvefac;

        for (std::size_t jnode=0; jnode<numnode; jnode++)
        {
          if (functnum>0)
          {
            // evaluate function at the position of the current node
            functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X(),time,NULL);
          }
          else functionfac = 1.0;

          edeadng(isd,jnode) = num*functionfac;
        }
      }
    }

    return edeadng;
}

/*!
 \brief get optimal gaussrule for given discretization type
        template version
 \return enum with optimal gaussrule for this element
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
    const DRT::Element::DiscretizationType  celldistype,
    const bool                              fast_integration
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
          if (fast_integration)
          {
            if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_5point;
            else                               rule = DRT::UTILS::intrule_tet_5point;
          }
          else
          {
            if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_15point;
            else                               rule = DRT::UTILS::intrule_tet_24point;
          }
          break;
        }
        case DRT::Element::hex20:
        {
          if (fast_integration)
          {
            if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_11point;
            else                               rule = DRT::UTILS::intrule_tet_11point;
          }
          else
          {
            if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_24point;
            else                               rule = DRT::UTILS::intrule_tet_45point;
          }
          break;
        }
        case DRT::Element::hex27:
        {
          if (fast_integration)
          {
            if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_45point;
            else                               rule = DRT::UTILS::intrule_tet_45point;
          }
          else
          {
            if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_343point_peano;
            else                               rule = DRT::UTILS::intrule_tet_343point_peano;
          }
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
    // TODO: Clarify this! hex integration cells for uncut refinement cells in an element
    else if (celldistype == DRT::Element::hex8)
    {
      switch (XELEDISTYPE)
      {
        case DRT::Element::hex8:
        {
          rule = getOptimalVelocityGaussrule(celldistype);
          // Ich habe xeleGeoType == GEO::LINEAR
//            if (xeleGeoType == GEO::CARTESIAN) // constant jacobian
//              rule = DRT::UTILS::intrule_hex_8point; //27point?
//            else
//              rule = DRT::UTILS::intrule_hex_8point;
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
 \brief get optimal gaussrule for given discretization type

 \return enum with optimal gausrule for this element

 for kink enrichment

 \author rasthofer 08/09
 */
template<DRT::Element::DiscretizationType XELEDISTYPE, class M>
DRT::UTILS::GaussRule3D getXFEMGaussruleKinkEnr(
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

    if (celldistype == DRT::Element::tet4)
    {
      switch (XELEDISTYPE)
      {
        case DRT::Element::hex8:
        {
          if (xeleGeoType == GEO::CARTESIAN) // constant jacobian
            rule = DRT::UTILS::intrule_tet_45point;
            //rule = DRT::UTILS::intrule_tet_24point;
            //rule = DRT::UTILS::intrule_tet_343point_peano;
          else
            rule = DRT::UTILS::intrule_tet_45point;
            //rule = DRT::UTILS::intrule_tet_24point;
            //rule = DRT::UTILS::intrule_tet_343point_peano;
          break;
        }
        default:
          dserror("no rule defined for this xele");
      }
    }
    // TODO: Clarify this! hex integration cells for uncut refinement cells in an element
    else if (celldistype == DRT::Element::hex8)
    {
      switch (XELEDISTYPE)
      {
        case DRT::Element::hex8:
        {
          //rule = getOptimalVelocityGaussrule(celldistype);
          rule = DRT::UTILS::intrule_hex_27point;
            break;
        }
        default:
          dserror("no rule defined for this xele");
      }
    }
    else
    {
      dserror("gaussrule only for tetrahedral and hexahedral integration cells defined. Please add additional rules here!");
    }
  }

  return rule;
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
