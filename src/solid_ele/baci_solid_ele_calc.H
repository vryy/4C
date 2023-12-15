/*! \file

\brief Declaration of routines for calculation of solid element simple displacement based

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_H
#define BACI_SOLID_ELE_CALC_H

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_calc_interface.H"

#include <memory>
#include <string>
#include <unordered_map>

BACI_NAMESPACE_OPEN

namespace MAT
{
  class So3Material;
}
namespace STR::MODELEVALUATOR
{
  class GaussPointDataOutputManager;
}
namespace DRT
{

  namespace ELEMENTS
  {
    template <CORE::FE::CellType celltype>
    class SolidEleCalc
    {
     public:
      SolidEleCalc();

      void Setup(MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef);

      void MaterialPostSetup(const DRT::Element& ele, MAT::So3Material& solid_material);

      void EvaluateNonlinearForceStiffnessMass(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
          CORE::LINALG::SerialDenseMatrix* mass_matrix);

      void EvaluateNonlinearForceStiffnessMassGEMM(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
          CORE::LINALG::SerialDenseMatrix* mass_matrix);

      void Recover(const DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params);

      void CalculateStress(const DRT::Element& ele, MAT::So3Material& solid_material,
          const StressIO& stressIO, const StrainIO& strainIO,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      double CalculateInternalEnergy(const DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      void Update(const DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      void InitializeGaussPointDataOutput(const DRT::Element& ele,
          const MAT::So3Material& solid_material,
          STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

      void EvaluateGaussPointDataOutput(const DRT::Element& ele,
          const MAT::So3Material& solid_material,
          STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

      void ResetToLastConverged(const DRT::Element& ele, MAT::So3Material& solid_material);

     private:
      /// static values for matrix sizes
      static constexpr int num_nodes_ = CORE::FE::num_nodes<celltype>;
      static constexpr int num_dim_ = CORE::FE::dim<celltype>;
      static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
      static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

      CORE::DRT::UTILS::GaussIntegration stiffness_matrix_integration_;
      CORE::DRT::UTILS::GaussIntegration mass_matrix_integration_;

    };  // class SolidEleCalc
  }     // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // SOLID_ELE_CALC_H
