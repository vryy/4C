/*----------------------------------------------------------------------*/
/*!
\file xfem_element_utils.H

\brief collection of routines for XFEM elements


\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_ELEMENT_UTILS_H
#define XFEM_ELEMENT_UTILS_H


#include "../drt_fem_general/drt_utils_integration.H"
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

  // fill stress/strain tensor
  template<class M1, class M2>
  void fill_tensor(
      const size_t         numparamstress,
      const M1&            shp_stress,
      const M2&            etensor,
      LINALG::Matrix<3,3>& tensor
  )
  {
    tensor.Clear();
    for (size_t i = 0; i < numparamstress; ++i)
    {
      const double fac = shp_stress(i);
      tensor(0,0) += fac*etensor(0,i);
      tensor(1,1) += fac*etensor(1,i);
      tensor(2,2) += fac*etensor(2,i);
      tensor(0,1) += fac*etensor(3,i);
      tensor(0,2) += fac*etensor(4,i);
      tensor(1,2) += fac*etensor(5,i);
    }
    tensor(1,0) = tensor(0,1);
    tensor(2,0) = tensor(0,2);
    tensor(2,1) = tensor(1,2);
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
        template version
 \return enum with optimal gaussrule for this element
   */
  inline DRT::UTILS::GaussRule2D getOptimalFaceGaussrule(DRT::Element::DiscretizationType  distype)
  {
    DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
    switch (distype)
    {
    case DRT::Element::quad4:
      rule = DRT::UTILS::intrule_quad_4point;
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
      else if (celldistype == DRT::Element::hex8)
      {
        switch (XELEDISTYPE)
        {
        case DRT::Element::hex8:
        {
          rule = DRT::UTILS::intrule_hex_64point;
          break;
        }
        case DRT::Element::hex20:
        {
          rule = DRT::UTILS::intrule_hex_125point;
          break;
        }
        case DRT::Element::hex27:
        {
          rule = DRT::UTILS::intrule_hex_125point;
          break;
        }
        default:
          dserror("no rule defined for this xele");
        }
      }
      else if (celldistype == DRT::Element::wedge6)
      {
        rule = DRT::UTILS::intrule_wedge_9point;
      }
      else if (celldistype == DRT::Element::pyramid5)
      {
        rule = DRT::UTILS::intrule_pyramid_8point;
      }
      else
      {
        dserror("Please add additional rules here!");
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
          //rule = DRT::UTILS::intrule_hex_1point;
          //rule = DRT::UTILS::intrule_hex_8point;
          rule = DRT::UTILS::intrule_hex_27point;
          //rule = DRT::UTILS::intrule_hex_64point;
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

  //! interpolate from nodal vector array to integration point vector using the shape function
  template <class M, class VS>
  double interpolateScalarFieldToIntPoint(
      const M&  eleVectorField,       ///< array with nodal vector values
      const VS& shp,                  ///< array with nodal shape function
      const int numparam              ///< number of parameters
      )
  {
    double s = 0.0;
    for (int iparam = 0; iparam < numparam; ++iparam)
        s += eleVectorField(iparam,0)*shp(iparam);
    return s;
  }

template <size_t numparam>
class ApproxFuncNormalPair
{
public:
  /// default constructor
  explicit ApproxFuncNormalPair()
  {
    // initialize approximation functions
    for(size_t irow=0;irow<numparam;++irow)
    {
      s(irow,0)=0.0;
      n(irow,0)=0.0;
    }
  };

  /// constructor
  explicit ApproxFuncNormalPair(bool setzero)
  {
    if (setzero)
    {
      // initialize approximation functions
      for(size_t irow=0;irow<numparam;++irow)
      {
        s(irow,0)=0.0;
        n(irow,0)=0.0;
      }
    }
  };

  void Clear()
  {
    for(size_t irow=0;irow<numparam;++irow)
    {
      s(irow,0)=0.0;
      n(irow,0)=0.0;
    }
  };

  void Print()
  {
    std::cout << "s ";
    for(size_t irow=0;irow<numparam;++irow)
    {
      std::cout << s(irow,0) << " ";
    }
    std::cout << std::endl;
    std::cout << "n ";
    for(size_t irow=0;irow<numparam;++irow)
    {
      std::cout << n(irow,0) << " ";
    }
    std::cout << std::endl;
  };

  /// approximation functions
  LINALG::Matrix<numparam,1> s;
  LINALG::Matrix<numparam,1> n;

private:
  /// hide copy constructor
  explicit ApproxFuncNormalPair(const ApproxFuncNormalPair&){}
};

template <unsigned int numdiv, size_t numparam>
class ApproxFuncNormalSet
{
private:
  /// hide default constructor
  explicit ApproxFuncNormalSet(){}
  /// hide default copy constructor
  explicit ApproxFuncNormalSet(const ApproxFuncNormalSet&){}
};

template <size_t numparam>
class ApproxFuncNormalSet<0,numparam>
{
public:
  /// default constructor
  explicit ApproxFuncNormalSet(){}
  /// constructor
  explicit ApproxFuncNormalSet(bool setzero) :
    d0(setzero){}

  void Clear()
  {
    d0.Clear();
  }

  void Print()
  {
    std::cout << "  d0" << std::endl; d0.Print();
  }

  /// pair of approximation functions belonging to this set
  ApproxFuncNormalPair<numparam> d0;

private:
  /// hide default copy constructor
  explicit ApproxFuncNormalSet(const ApproxFuncNormalSet&){}
};

template <size_t numparam>
class ApproxFuncNormalSet<1,numparam>
{
public:
  /// default constructor
  explicit ApproxFuncNormalSet(){}
  /// constructor
  explicit ApproxFuncNormalSet(bool setzero) :
    d0(setzero),
    dx(setzero),
    dy(setzero),
    dz(setzero){}

  void Clear()
  {
    d0.Clear();
    dx.Clear();
    dy.Clear();
    dz.Clear();
  }

  void Print()
  {
    std::cout << "  d0" << std::endl; d0.Print();
    std::cout << "  dx" << std::endl; dx.Print();
    std::cout << "  dy" << std::endl; dy.Print();
    std::cout << "  dz" << std::endl; dz.Print();
  }

  /// pairs of approximation functions belonging to this set
  ApproxFuncNormalPair<numparam>  d0;
  ApproxFuncNormalPair<numparam>  dx;
  ApproxFuncNormalPair<numparam>  dy;
  ApproxFuncNormalPair<numparam>  dz;

private:
  /// hide default copy constructor
  explicit ApproxFuncNormalSet(const ApproxFuncNormalSet&){}
};

template <size_t numparam>
class ApproxFuncNormalSet<2,numparam>
{
public:
  /// default constructor
  explicit ApproxFuncNormalSet(){}
  /// constructor
  explicit ApproxFuncNormalSet(bool setzero) :
    d0(setzero),
    dx(setzero),
    dy(setzero),
    dz(setzero),
    dxdx(setzero),
    dxdy(setzero),
    dxdz(setzero),
    dydx(setzero),
    dydy(setzero),
    dydz(setzero),
    dzdx(setzero),
    dzdy(setzero),
    dzdz(setzero){}

  void Clear()
  {
    d0.Clear();
    dx.Clear();
    dy.Clear();
    dz.Clear();
    dxdx.Clear();
    dxdy.Clear();
    dxdz.Clear();
    dydx.Clear();
    dydy.Clear();
    dydz.Clear();
    dzdx.Clear();
    dzdy.Clear();
    dzdz.Clear();
  };

  void Print()
  {
    std::cout << "  d0"   << std::endl;   d0.Print();
    std::cout << "  dx"   << std::endl;   dx.Print();
    std::cout << "  dy"   << std::endl;   dy.Print();
    std::cout << "  dz"   << std::endl;   dz.Print();
    std::cout << "  dxdx" << std::endl; dxdx.Print();
    std::cout << "  dxdy" << std::endl; dxdy.Print();
    std::cout << "  dxdz" << std::endl; dxdz.Print();
    std::cout << "  dydx" << std::endl; dydx.Print();
    std::cout << "  dydy" << std::endl; dydy.Print();
    std::cout << "  dydz" << std::endl; dydz.Print();
    std::cout << "  dzdx" << std::endl; dzdx.Print();
    std::cout << "  dzdy" << std::endl; dzdy.Print();
    std::cout << "  dzdz" << std::endl; dzdz.Print();
  };

  /// pairs of approximation functions belonging to this set
  ApproxFuncNormalPair<numparam>  d0;
  ApproxFuncNormalPair<numparam>  dx;
  ApproxFuncNormalPair<numparam>  dy;
  ApproxFuncNormalPair<numparam>  dz;
  ApproxFuncNormalPair<numparam>  dxdx;
  ApproxFuncNormalPair<numparam>  dxdy;
  ApproxFuncNormalPair<numparam>  dxdz;
  ApproxFuncNormalPair<numparam>  dydx;
  ApproxFuncNormalPair<numparam>  dydy;
  ApproxFuncNormalPair<numparam>  dydz;
  ApproxFuncNormalPair<numparam>  dzdx;
  ApproxFuncNormalPair<numparam>  dzdy;
  ApproxFuncNormalPair<numparam>  dzdz;

private:
  /// hide default copy constructor
  explicit ApproxFuncNormalSet(const ApproxFuncNormalSet&){}
};

template <unsigned int numdiv, size_t numparam>
class ApproxFuncNormalVector
{
public:
  /// default constructor
  explicit ApproxFuncNormalVector(){}
  /// constructor
  explicit ApproxFuncNormalVector(bool setzero) :
    velx(setzero),
    vely(setzero),
    velz(setzero){}

  void Clear()
  {
    velx.Clear();
    vely.Clear();
    velz.Clear();
  };

  void Print()
  {
    std::cout << "velx" << std::endl; velx.Print();
    std::cout << "vely" << std::endl; vely.Print();
    std::cout << "velz" << std::endl; velz.Print();
  };

  /// sets of approximation functions belonging to this vector
  ApproxFuncNormalSet<numdiv,numparam> velx;
  ApproxFuncNormalSet<numdiv,numparam> vely;
  ApproxFuncNormalSet<numdiv,numparam> velz;

private:
  /// hide default copy constructor
  explicit ApproxFuncNormalVector(const ApproxFuncNormalVector&){}
};

template <size_t numparam>
struct EnrViscs2N
{
  ApproxFuncNormalPair<numparam> xx;
  ApproxFuncNormalPair<numparam> xy;
  ApproxFuncNormalPair<numparam> xz;
  ApproxFuncNormalPair<numparam> yx;
  ApproxFuncNormalPair<numparam> yy;
  ApproxFuncNormalPair<numparam> yz;
  ApproxFuncNormalPair<numparam> zx;
  ApproxFuncNormalPair<numparam> zy;
  ApproxFuncNormalPair<numparam> zz;
};


}

#endif
