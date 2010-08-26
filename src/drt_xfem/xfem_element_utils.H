/*----------------------------------------------------------------------*/
/*!
\file xfem_element_utils.H

\brief collection of routines for XFEM elements

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef XFEM_ELEMENT_UTILS_H
#define XFEM_ELEMENT_UTILS_H

#include "../drt_geometry/intersection_service.H"

namespace XFEM
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
  ///  psi denotes the enrichment function, e.g. a step function H
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
  template <unsigned int numdiv, size_t numparam>
  struct ApproxFunc
  {
    // intentionally left empty -> Compiler cries, if numdiv is not one of the following structs
  };

  template <size_t numparam>
  struct ApproxFunc<0, numparam>
  {
    LINALG::Matrix<numparam,1>  d0;
  };

  template <size_t numparam>
  struct ApproxFunc<1, numparam>
  {
    LINALG::Matrix<numparam,1>  d0;
    LINALG::Matrix<numparam,1>  dx;
    LINALG::Matrix<numparam,1>  dy;
    LINALG::Matrix<numparam,1>  dz;
  };

  template <size_t numparam>
  struct ApproxFunc<2, numparam>
  {
    LINALG::Matrix<numparam,1>  d0;
    LINALG::Matrix<numparam,1>  dx;
    LINALG::Matrix<numparam,1>  dy;
    LINALG::Matrix<numparam,1>  dz;
    LINALG::Matrix<numparam,1>  dxdx;
    LINALG::Matrix<numparam,1>  dxdy;
    LINALG::Matrix<numparam,1>  dxdz;
    LINALG::Matrix<numparam,1>  dydx;
    LINALG::Matrix<numparam,1>  dydy;
    LINALG::Matrix<numparam,1>  dydz;
    LINALG::Matrix<numparam,1>  dzdx;
    LINALG::Matrix<numparam,1>  dzdy;
    LINALG::Matrix<numparam,1>  dzdz;
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
          if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_24point;
          else                               rule = DRT::UTILS::intrule_tet_24point;
          break;
        }
        case DRT::Element::hex20:
        {
          if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_45point;
          else                               rule = DRT::UTILS::intrule_tet_45point;
          break;
        }
        case DRT::Element::hex27:
        {
          if (xeleGeoType == GEO::CARTESIAN) rule = DRT::UTILS::intrule_tet_343point_peano;
          else                               rule = DRT::UTILS::intrule_tet_343point_peano;
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

    if (not DRT::UTILS::IntRuleIsOpen(rule))
      dserror("You have to use an open integration rule for XFEM volume integrals!");

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
    LINALG::Matrix<3,1> v(true);
    const int nsd = 3;
    for (int iparam = 0; iparam < numparam; ++iparam)
      for (int isd = 0; isd < nsd; ++isd)
        v(isd) += eleVectorField(isd,iparam)*shp(iparam);
    return v;
  }

}

#endif
#endif
