/*! \file

\brief Definition of routines for calculation of solid poro element with pressure based
implementation

\level 1
*/

#ifndef BACI_SOLID_ELE_PORO_CALC_PRESSURE_BASED_H
#define BACI_SOLID_ELE_PORO_CALC_PRESSURE_BASED_H


#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_element.H"
#include "baci_lib_element_integration_select.H"
#include "baci_solid_ele_poro_calc_interface.H"
#include "baci_solid_ele_utils.H"
#include "baci_utils_singleton_owner.H"

namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace DRT
{
  namespace UTILS
  {
    template <DRT::Element::DiscretizationType>
    struct DisTypeToNumNodePerEle;
    template <DRT::Element::DiscretizationType>
    struct DisTypeToDim;
  }  // namespace UTILS

  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class SolidPoroPressureBasedEleCalc : public virtual SolidPoroEleCalcInterface
    {
     public:
      SolidPoroPressureBasedEleCalc();

      /// Singleton access method
      static SolidPoroPressureBasedEleCalc<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      void EvaluateNonlinearForceStiffness(const DRT::Element& ele, MAT::StructPoro& porostructmat,
          MAT::FluidPoroMultiPhase& porofluidmat, const INPAR::STR::KinemType& kinematictype,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix) override;

      void CouplingPoroelast(const DRT::Element& ele, MAT::StructPoro& porostructmat,
          MAT::FluidPoroMultiPhase& porofluidmat, const INPAR::STR::KinemType& kinematictype,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseMatrix& stiffness_matrix) override;

      void CouplingStress(const DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params) override;

      void PoroSetup(MAT::StructPoro& porostructmat, DRT::INPUT::LineDefinition* linedef) override;

     private:
      /// static values for matrix sizes
      static constexpr int nen_ =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static constexpr int nsd_ = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;
      static constexpr int numdofperelement_ = nen_ * nsd_;
      static constexpr int numstr_ = nsd_ * (nsd_ + 1) / 2;

      CORE::DRT::UTILS::GaussIntegration gauss_integration_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SOLID_ELE_PORO_CALC_PRESSURE_BASED_H