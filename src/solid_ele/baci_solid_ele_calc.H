/*! \file

\brief Declaration of routines for calculation of solid element simple displacement based

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_H
#define BACI_SOLID_ELE_CALC_H

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_element.H"
#include "baci_lib_element_integration_select.H"
#include "baci_solid_ele_calc_interface.H"
#include "baci_solid_ele_utils.H"

#include <memory>
#include <string>
#include <unordered_map>


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
    class SolidEleCalc : public virtual SolidEleCalcInterface
    {
     public:
      SolidEleCalc();

      void Setup(MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef) override;

      void MaterialPostSetup(const DRT::Element& ele, MAT::So3Material& solid_material) override;

      void EvaluateNonlinearForceStiffnessMass(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
          CORE::LINALG::SerialDenseMatrix* mass_matrix) override;

      void EvaluateNonlinearForceStiffnessMassGEMM(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
          CORE::LINALG::SerialDenseMatrix* mass_matrix) override;

      void Recover(const DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params) override;

      void CalculateStress(const DRT::Element& ele, MAT::So3Material& solid_material,
          const StressIO& stressIO, const StrainIO& strainIO,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) override;

      double CalculateInternalEnergy(const DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) override;

      void Update(const DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) override;

      void InitializeGaussPointDataOutput(const DRT::Element& ele,
          const MAT::So3Material& solid_material,
          STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const override;

      void EvaluateGaussPointDataOutput(const DRT::Element& ele,
          const MAT::So3Material& solid_material,
          STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const override;

      void ResetAll(const DRT::Element& ele, MAT::So3Material& solid_material) override;

      void ResetToLastConverged(const DRT::Element& ele, MAT::So3Material& solid_material) override;

     private:
      /// static values for matrix sizes
      static constexpr int num_nodes_ =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static constexpr int num_dim_ = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;
      static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
      static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

      CORE::DRT::UTILS::GaussIntegration stiffness_matrix_integration_;
      CORE::DRT::UTILS::GaussIntegration mass_matrix_integration_;

    };  // class SolidEleCalc
  }     // namespace ELEMENTS
}  // namespace DRT

#endif  // SOLID_ELE_CALC_H
