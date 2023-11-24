/*! \file

\brief Declaration of routines for calculation of solid element displacement
based with MULF prestressing

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_MULF_H
#define BACI_SOLID_ELE_CALC_MULF_H

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_discret.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_calc_interface.H"
#include "baci_solid_ele_interface_serializable.H"
#include "baci_utils_demangle.H"

#include <memory>
#include <string>
#include <unordered_map>
namespace MAT
{
  class So3Material;
}
namespace STR::MODELEVALUATOR
{
  class GaussPointDataOutputManager;
}
namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType distype>
  struct MulfHistoryData
  {
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
        inverse_jacobian;
    CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
        CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
        deformation_gradient;
    bool is_setup = false;

    MulfHistoryData()
    {
      for (int i = 0; i < CORE::DRT::UTILS::DisTypeToDim<distype>::dim; ++i)
      {
        for (int j = 0; j < CORE::DRT::UTILS::DisTypeToDim<distype>::dim; ++j)
        {
          inverse_jacobian(i, j) = static_cast<double>(i == j);
          deformation_gradient(i, j) = static_cast<double>(i == j);
        }
      }
    }
  };

  template <CORE::FE::CellType distype>
  class SolidEleCalcMulf
  {
   public:
    SolidEleCalcMulf();

    void Pack(DRT::PackBuffer& data) const;

    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

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

    void UpdatePrestress(const DRT::Element& ele, MAT::So3Material& solid_material,
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
    static constexpr int num_nodes_ = CORE::FE::num_nodes<distype>;
    static constexpr int num_dim_ = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;
    static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
    static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

    std::vector<MulfHistoryData<distype>> history_data_;

    CORE::DRT::UTILS::GaussIntegration stiffness_matrix_integration_;
    CORE::DRT::UTILS::GaussIntegration mass_matrix_integration_;

  };  // class SolidEleCalc

  template <typename T, typename AlwaysVoid = void>
  constexpr bool IsPrestressUpdateable = false;

  template <typename T>
  constexpr bool IsPrestressUpdateable<T,
      std::void_t<decltype(std::declval<T>()->UpdatePrestress(std::declval<const DRT::Element&>(),
          std::declval<MAT::So3Material&>(), std::declval<const DRT::Discretization&>(),
          std::declval<const std::vector<int>&>(), std::declval<Teuchos::ParameterList&>()))>> =
      true;

  namespace DETAILS
  {
    struct UpdatePrestressAction
    {
      UpdatePrestressAction(const DRT::Element& e, MAT::So3Material& m,
          const DRT::Discretization& d, const std::vector<int>& lmvec, Teuchos::ParameterList& p)
          : element(e), mat(m), discretization(d), lm(lmvec), params(p)
      {
      }

      template <typename T, std::enable_if_t<IsPrestressUpdateable<T&>, bool> = true>
      void operator()(T& updateable)
      {
        updateable->UpdatePrestress(element, mat, discretization, lm, params);
      }

      template <typename T, std::enable_if_t<!IsPrestressUpdateable<T&>, bool> = true>
      void operator()(T& other)
      {
        dserror(
            "Your element evaluation %s does not allow to update prestress. You may need to add "
            "MULF to your element line definitions.",
            CORE::UTILS::TryDemangle(typeid(T).name()).c_str());
      }

      const DRT::Element& element;
      MAT::So3Material& mat;
      const DRT::Discretization& discretization;
      const std::vector<int>& lm;
      Teuchos::ParameterList& params;
    };
  }  // namespace DETAILS

  template <typename VariantType>
  void UpdatePrestress(VariantType& variant, const DRT::Element& element, MAT::So3Material& mat,
      const DRT::Discretization& discretization, const std::vector<int>& lm,
      Teuchos::ParameterList& params)
  {
    std::visit(DETAILS::UpdatePrestressAction(element, mat, discretization, lm, params), variant);
  }
}  // namespace DRT::ELEMENTS

#endif  // BACI_SOLID_ELE_CALC_MULF_H
