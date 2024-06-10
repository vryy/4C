/*! \file

\brief Definition of routines for calculation of solid poro element with pressure based
implementation

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_PRESSURE_BASED_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_PRESSURE_BASED_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class FluidPoroMultiPhase;
  class StructPoro;
}  // namespace Mat
namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace Discret
{

  namespace ELEMENTS
  {
    template <Core::FE::CellType celltype>
    class SolidPoroPressureBasedEleCalc
    {
     public:
      SolidPoroPressureBasedEleCalc();

      void evaluate_nonlinear_force_stiffness(const Core::Elements::Element& ele,
          Mat::StructPoro& porostructmat, Mat::FluidPoroMultiPhase& porofluidmat,
          const Inpar::STR::KinemType& kinematictype,
          const Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
          Core::LinAlg::SerialDenseVector* force_vector,
          Core::LinAlg::SerialDenseMatrix* stiffness_matrix);

      void coupling_poroelast(const Core::Elements::Element& ele, Mat::StructPoro& porostructmat,
          Mat::FluidPoroMultiPhase& porofluidmat, const Inpar::STR::KinemType& kinematictype,
          const Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
          Core::LinAlg::SerialDenseMatrix& stiffness_matrix);

      void CouplingStress(const Core::Elements::Element& ele,
          const Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      void PoroSetup(Mat::StructPoro& porostructmat, Input::LineDefinition* linedef);

     private:
      /// static values for matrix sizes
      static constexpr int num_nodes_ = Core::FE::num_nodes<celltype>;
      static constexpr int num_dim_ = Core::FE::dim<celltype>;
      static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
      static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

      Core::FE::GaussIntegration gauss_integration_;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif