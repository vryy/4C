/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element with EAS element technology
\level 1

*----------------------------------------------------------------------*/

#ifndef BACI_SOLID_ELE_CALC_EAS_H
#define BACI_SOLID_ELE_CALC_EAS_H

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_element.H"
#include "baci_lib_element_integration_select.H"
#include "baci_solid_ele_calc_interface.H"
#include "baci_solid_ele_interface_serializable.H"
#include "baci_solid_ele_utils.H"

#include <memory>
#include <string>
#include <unordered_map>

namespace STR::ELEMENTS
{
  class ParamsInterface;

  enum class EasType
  {
    soh8_easnone,
    eastype_h8_9,
    eastype_h8_21,
    eastype_sh8_7,
    eastype_sh18_9,
    eastype_undefined
  };

  template <STR::ELEMENTS::EasType eastype>
  struct EasTypeToNumEas
  {
  };
  template <>
  struct EasTypeToNumEas<STR::ELEMENTS::EasType::eastype_h8_9>
  {
    static const int neas = 9;
  };
  template <>
  struct EasTypeToNumEas<STR::ELEMENTS::EasType::eastype_h8_21>
  {
    static const int neas = 21;
  };
  template <>
  struct EasTypeToNumEas<STR::ELEMENTS::EasType::eastype_sh8_7>
  {
    static const int neas = 7;
  };
  template <>
  struct EasTypeToNumEas<STR::ELEMENTS::EasType::eastype_sh18_9>
  {
    static const int neas = 9;
  };
  template <>
  struct EasTypeToNumEas<STR::ELEMENTS::EasType::eastype_undefined>
  {
  };

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
    /// struct for EAS matrices and vectors to be stored between iterations
    template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
    struct EasIterationData
    {
      constexpr static int neas = STR::ELEMENTS::EasTypeToNumEas<eastype>::neas;

      /// inverse EAS matrix K_{alpha alpha}, sometimes called Dtilde
      CORE::LINALG::Matrix<neas, neas> invKaa_{true};

      /// EAS matrix K_{d alpha}, sometimes called Ltilde
      CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement *
                               CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
          neas>
          Kda_{true};

      /// EAS portion of internal forces f_{eas}, sometimes called enhacement vector s or Rtilde
      CORE::LINALG::Matrix<neas, 1> feas_{true};

      /// discrete enhanced strain scalars alpha
      CORE::LINALG::Matrix<neas, 1> alpha_{true};
    };

    template <DRT::Element::DiscretizationType distype, STR::ELEMENTS::EasType eastype>
    class SolidEleCalcEas : public SolidEleCalcInterface, public Serializable
    {
     public:
      SolidEleCalcEas();

      void Setup(MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef) override;

      void Pack(DRT::PackBuffer& data) const override;

      void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data) override;

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
      /// EAS matrices and vectors to be stored between iterations
      DRT::ELEMENTS::EasIterationData<distype, eastype> eas_iteration_data_ = {};

      /// static values for matrix sizes
      static constexpr int nen_ =
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      static constexpr int nsd_ = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;
      static constexpr int numdofperelement_ = nen_ * nsd_;
      static constexpr int numstr_ = nsd_ * (nsd_ + 1) / 2;

      CORE::DRT::UTILS::GaussIntegration stiffness_matrix_integration_;
      CORE::DRT::UTILS::GaussIntegration mass_matrix_integration_;
    };  // class SolidEleCalcEas
  }     // namespace ELEMENTS
}  // namespace DRT

#endif  // SOLID_ELE_CALC_EAS_H
