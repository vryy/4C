/*----------------------------------------------------------------------*/
/*!
\file combust3_interpolation.H

\brief rotuines to assemble element stiffness matrix and vector

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 */
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET


#ifndef COMBUST3_INTERPOLATION_H
#define COMBUST3_INTERPOLATION_H

#include "../drt_lib/drt_element.H"
#include "../drt_xfem/physics.H"


namespace COMBUST
{
  
  /*!
   * \brief Returns the distype for stress interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for stress interpolation
   */
  template<DRT::Element::DiscretizationType veldistype>
  struct StressInterpolation3D {};
  /// @cond
  template<> struct StressInterpolation3D<DRT::Element::tet4>  {static const DRT::Element::DiscretizationType distype = DRT::Element::tet4;};
  template<> struct StressInterpolation3D<DRT::Element::tet10> {static const DRT::Element::DiscretizationType distype = DRT::Element::tet10;};
  template<> struct StressInterpolation3D<DRT::Element::hex8>  {static const DRT::Element::DiscretizationType distype = DRT::Element::hex8;};
  template<> struct StressInterpolation3D<DRT::Element::hex20> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex20;};
  template<> struct StressInterpolation3D<DRT::Element::hex27> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex27;};
  /// @endcond
  
  
  /*!
   * \brief Returns the distype for discont. pressure interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for discont. pressure interpolation
   */
  template<DRT::Element::DiscretizationType veldistype>
  struct DiscPressureInterpolation3D {};
  /// @cond
  template<> struct DiscPressureInterpolation3D<DRT::Element::tet4>  {static const DRT::Element::DiscretizationType distype = DRT::Element::tet4;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::tet10> {static const DRT::Element::DiscretizationType distype = DRT::Element::tet10;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::hex8>  {static const DRT::Element::DiscretizationType distype = DRT::Element::hex8;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::hex20> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex20;};
  template<> struct DiscPressureInterpolation3D<DRT::Element::hex27> {static const DRT::Element::DiscretizationType distype = DRT::Element::hex27;};
  /// @endcond
  
  /*!
   * \brief Returns the distype for stress interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for stress interpolation
   */  
  inline DRT::Element::DiscretizationType getStressInterpolationType3D(DRT::Element::DiscretizationType veldistype)
  {
    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
    switch (veldistype)
    {
      case DRT::Element::tet4:
        distype = StressInterpolation3D<DRT::Element::tet4>::distype;
        break;
      case DRT::Element::tet10:
        distype = StressInterpolation3D<DRT::Element::tet10>::distype;
        break;
      case DRT::Element::hex8:
        distype = StressInterpolation3D<DRT::Element::hex8>::distype;
        break;
      case DRT::Element::hex20:
        distype = StressInterpolation3D<DRT::Element::hex20>::distype;
        break;
      case DRT::Element::hex27:
        distype = StressInterpolation3D<DRT::Element::hex27>::distype;
        break;
      default:
        dserror("no stressdistype defined for velocity distype");
    }
    return distype;
  }
  
  /*!
   * \brief Returns the distype for discont. pressure interpolation
   *        for each parent element discretization type in 3d
   * 
   * \return distype for discont. pressure interpolation
   */
  inline DRT::Element::DiscretizationType getDiscPressureInterpolationType3D(DRT::Element::DiscretizationType veldistype)
  {
    DRT::Element::DiscretizationType distype = DRT::Element::dis_none;
    switch (veldistype)
    {
      case DRT::Element::tet4:
        distype = DiscPressureInterpolation3D<DRT::Element::tet4>::distype;
        break;
      case DRT::Element::tet10:
        distype = DiscPressureInterpolation3D<DRT::Element::tet10>::distype;
        break;
      case DRT::Element::hex8:
        distype = DiscPressureInterpolation3D<DRT::Element::hex8>::distype;
        break;
      case DRT::Element::hex20:
        distype = DiscPressureInterpolation3D<DRT::Element::hex20>::distype;
        break;
      case DRT::Element::hex27:
        distype = DiscPressureInterpolation3D<DRT::Element::hex27>::distype;
        break;
      default:
        dserror("no discpressuredistype defined for velocity distype");
    }
    return distype;
  }
  
  //! for a given velocity interpolation, a map is produced which has the number of element unknowns per physical field
  inline map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> getElementAnsatz(DRT::Element::DiscretizationType eledistype)
  {
    map<XFEM::PHYSICS::Field, DRT::Element::DiscretizationType> element_ansatz;
    //          element_ansatz[XFEM::PHYSICS::Tauxx] = COMBUST::getStressInterpolationType3D(eledistype);
    //          element_ansatz[XFEM::PHYSICS::Tauyy] = COMBUST::getStressInterpolationType3D(eledistype);
    //          element_ansatz[XFEM::PHYSICS::Tauzz] = COMBUST::getStressInterpolationType3D(eledistype);
    //          element_ansatz[XFEM::PHYSICS::Tauxy] = COMBUST::getStressInterpolationType3D(eledistype);
    //          element_ansatz[XFEM::PHYSICS::Tauxz] = COMBUST::getStressInterpolationType3D(eledistype);
    //          element_ansatz[XFEM::PHYSICS::Tauyz] = COMBUST::getStressInterpolationType3D(eledistype);
    //element_ansatz[XFEM::PHYSICS::DiscPres] = COMBUST::getDiscPressureInterpolationType3D(eledistype);
    element_ansatz[XFEM::PHYSICS::Sigmaxx] = COMBUST::getStressInterpolationType3D(eledistype);
    element_ansatz[XFEM::PHYSICS::Sigmayy] = COMBUST::getStressInterpolationType3D(eledistype);
    element_ansatz[XFEM::PHYSICS::Sigmazz] = COMBUST::getStressInterpolationType3D(eledistype);
    element_ansatz[XFEM::PHYSICS::Sigmaxy] = COMBUST::getStressInterpolationType3D(eledistype);
    element_ansatz[XFEM::PHYSICS::Sigmaxz] = COMBUST::getStressInterpolationType3D(eledistype);
    element_ansatz[XFEM::PHYSICS::Sigmayz] = COMBUST::getStressInterpolationType3D(eledistype);
    return element_ansatz;
  }
}


#endif /*COMBUST3_INTERPOLATION_H*/
#endif
