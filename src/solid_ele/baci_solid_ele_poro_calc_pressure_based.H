/*! \file

\brief Definition of routines for calculation of solid poro element with pressure based
implementation

\level 1
*/

#ifndef BACI_SOLID_ELE_PORO_CALC_PRESSURE_BASED_H
#define BACI_SOLID_ELE_PORO_CALC_PRESSURE_BASED_H


#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_inpar_structure.H"
#include "baci_lib_element.H"

namespace MAT
{
  class FluidPoroMultiPhase;
  class StructPoro;
}  // namespace MAT
namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace DRT
{

  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class SolidPoroPressureBasedEleCalc
    {
     public:
      SolidPoroPressureBasedEleCalc();

      void EvaluateNonlinearForceStiffness(const DRT::Element& ele, MAT::StructPoro& porostructmat,
          MAT::FluidPoroMultiPhase& porofluidmat, const INPAR::STR::KinemType& kinematictype,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix);

      void CouplingPoroelast(const DRT::Element& ele, MAT::StructPoro& porostructmat,
          MAT::FluidPoroMultiPhase& porofluidmat, const INPAR::STR::KinemType& kinematictype,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& stiffness_matrix);

      void CouplingStress(const DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params);

      void PoroSetup(MAT::StructPoro& porostructmat, DRT::INPUT::LineDefinition* linedef);

     private:
      /// static values for matrix sizes
      static constexpr int num_nodes_ = CORE::FE::num_nodes<distype>;
      static constexpr int num_dim_ = CORE::FE::dim<distype>;
      static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
      static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

      CORE::DRT::UTILS::GaussIntegration gauss_integration_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // SOLID_ELE_PORO_CALC_PRESSURE_BASED_H