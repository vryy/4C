/*----------------------------------------------------------------------*/
/*! \file
 \brief helpful methods and template definitions for the porofluidmultiphase element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROFLUIDMULTIPHASE_ELE_CALC_UTILS_H
#define BACI_POROFLUIDMULTIPHASE_ELE_CALC_UTILS_H

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_element.H"

namespace MAT
{
  class Material;
  class FluidPoroSinglePhase;
  class FluidPoroSingleVolFrac;
  class FluidPoroVolFracPressure;
  class FluidPoroMultiPhase;
  class FluidPoroMultiPhaseReactions;
  class FluidPoroSingleReaction;
}  // namespace MAT

namespace POROFLUIDMULTIPHASE
{
  namespace ELEUTILS
  {
    //! get the single phase material from the element material
    const MAT::FluidPoroSinglePhase& GetSinglePhaseMatFromMaterial(
        const MAT::Material& material, int phasenum);

    //! get the single phase material from the element multiphase material
    const MAT::FluidPoroSinglePhase& GetSinglePhaseMatFromMultiMaterial(
        const MAT::FluidPoroMultiPhase& multiphasemat, int phasenum);

    //! get the single volume fraction material from the element material
    const MAT::FluidPoroSingleVolFrac& GetSingleVolFracMatFromMaterial(
        const MAT::Material& material, int volfracnum);

    //! get the single volume fraction material from the element multiphase material
    const MAT::FluidPoroSingleVolFrac& GetSingleVolFracMatFromMultiMaterial(
        const MAT::FluidPoroMultiPhase& multiphasemat, int volfracnum);

    //! get the volume fraction pressure material from the element material
    const MAT::FluidPoroVolFracPressure& GetVolFracPressureMatFromMaterial(
        const MAT::Material& material, int volfracnum);

    //! get the volume fraction pressure material from the element multiphase material
    const MAT::FluidPoroVolFracPressure& GetVolFracPressureMatFromMultiMaterial(
        const MAT::FluidPoroMultiPhase& multiphasemat, int volfracnum);

    //! get the single phase material from the element multiphase reactions material
    MAT::FluidPoroSingleReaction& GetSingleReactionMatFromMultiReactionsMaterial(
        const MAT::FluidPoroMultiPhaseReactions& multiphasereacmat, int phasenum);

    /*!
    \brief Decide, whether second derivatives are needed  (template version)
     *  In convection-diffusion problems, ONLY N,xx , N,yy and N,zz are needed
     *  to evaluate the laplacian operator for the residual-based stabilization.
     *  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5
     *  return false although they have non-zero MIXED second derivatives.*/
    template <CORE::FE::CellType DISTYPE>
    struct Use2ndDerivs
    {
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::hex8>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::tet4>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::wedge6>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::pyramid5>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::nurbs8>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::quad4>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::nurbs4>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::tri3>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::line2>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::nurbs2>
    {
      static const bool use = false;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::hex20>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::hex27>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::nurbs27>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::tet10>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::quad8>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::quad9>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::nurbs9>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::tri6>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::line3>
    {
      static const bool use = true;
    };
    template <>
    struct Use2ndDerivs<CORE::FE::CellType::nurbs3>
    {
      static const bool use = true;
    };


    //! Template Meta Programming version of switch over discretization type
    template <CORE::FE::CellType DISTYPE>
    struct DisTypeToOptGaussRule
    {
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::hex8>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::hex20>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::hex27>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::tet4>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::tet10>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_5point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::wedge6>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::wedge_6point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::pyramid5>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule =
          CORE::DRT::UTILS::GaussRule3D::pyramid_8point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs8>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs27>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::quad4>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::quad8>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::quad9>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::tri3>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_3point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::tri6>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_6point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs4>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs9>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::line2>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::line3>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs2>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs3>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
    };


    //! Template Meta Programming version of switch over discretization type
    template <CORE::FE::CellType DISTYPE>
    struct DisTypeToGaussRuleForExactSol
    {
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex8>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex20>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex27>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tet4>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tet10>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::wedge6>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::pyramid5>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs8>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs27>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad4>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad8>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad9>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri3>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri6>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs4>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs9>
    {
      static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line2>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
    };  // not tested
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line3>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs2>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs3>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
    };

  }  // namespace ELEUTILS

}  // namespace POROFLUIDMULTIPHASE



#endif  // POROFLUIDMULTIPHASE_ELE_CALC_UTILS_H
