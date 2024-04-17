/*! \file

\brief Declaration of routines for calculation of solid element simple displacement based

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_gausspoints.hpp"
#include "baci_lib_element.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_solid_3D_ele_calc_interface.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"
#include "baci_solid_3D_ele_calc_lib_integration.hpp"
#include "baci_solid_3D_ele_calc_lib_io.hpp"
#include "baci_solid_3D_ele_calc_lib_nitsche.hpp"
#include "baci_solid_3D_ele_formulation.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

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

  /*!
   * @brief A solid evaluator for different element formulations
   *
   * This is a class that takes the element formulation @p ElementFormulation as a template
   * parameter and generates all necessary code that is needed such that a solid element can be
   * evaluated. The two additional template parameter @p PreparationData and @p HistoryData store
   * preparation data and history data.
   *
   * The type @p ElementFormulation defines static member function @p Evaluate(...) for the
   * evaluation of the deformation gradient at a specific position @p xi. Additionally, it has to
   * provide @p Prepare(...) for eventual preparation of data that is constant. It should return the
   * type @p PreparationData that is passed to every subsequent call. The static functions
   * @p AddInternalForceVector(...) and @p AddStiffnessMatrix(...) are adding the gauss point
   * contribution to the force vector or stiffness matrix, respectively. Some element formulations
   * require to store data @p HistoryData. This type needs to be default constructable. The static
   * member functions @p Pack(...) and @p Unpack(...) need to be defined to allow parallel
   * distribution and restarts.
   *
   * @tparam celltype : celltype of the evaluator
   * @tparam ElementFormulation : A type that defines the element formulation in static member
   * functions
   */
  template <CORE::FE::CellType celltype, typename ElementFormulation>
  class SolidEleCalc
  {
   public:
    SolidEleCalc();

    void Pack(CORE::COMM::PackBuffer& data) const;

    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

    void Setup(MAT::So3Material& solid_material, INPUT::LineDefinition* linedef);

    void MaterialPostSetup(const DRT::Element& ele, MAT::So3Material& solid_material);

    void EvaluateNonlinearForceStiffnessMass(const DRT::Element& ele,
        MAT::So3Material& solid_material, const DRT::Discretization& discretization,
        const std::vector<int>& lm, Teuchos::ParameterList& params,
        CORE::LINALG::SerialDenseVector* force_vector,
        CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
        CORE::LINALG::SerialDenseMatrix* mass_matrix);

    void Recover(const DRT::Element& ele, const DRT::Discretization& discretization,
        const std::vector<int>& lm, Teuchos::ParameterList& params);

    void Update(const DRT::Element& ele, MAT::So3Material& solid_material,
        const DRT::Discretization& discretization, const std::vector<int>& lm,
        Teuchos::ParameterList& params);

    void CalculateStress(const DRT::Element& ele, MAT::So3Material& solid_material,
        const StressIO& stressIO, const StrainIO& strainIO,
        const DRT::Discretization& discretization, const std::vector<int>& lm,
        Teuchos::ParameterList& params);

    double CalculateInternalEnergy(const DRT::Element& ele, MAT::So3Material& solid_material,
        const DRT::Discretization& discretization, const std::vector<int>& lm,
        Teuchos::ParameterList& params);

    void InitializeGaussPointDataOutput(const DRT::Element& ele,
        const MAT::So3Material& solid_material,
        STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

    void EvaluateGaussPointDataOutput(const DRT::Element& ele,
        const MAT::So3Material& solid_material,
        STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

    void UpdatePrestress(const DRT::Element& ele, MAT::So3Material& solid_material,
        const DRT::Discretization& discretization, const std::vector<int>& lm,
        Teuchos::ParameterList& params);

    void ResetToLastConverged(const DRT::Element& ele, MAT::So3Material& solid_material);

    CauchyNDirAndLinearization<3> GetCauchyNDirAndDerivativesAtXi(const DRT::Element& ele,
        MAT::So3Material& solid_material, const std::vector<double>& disp,
        const CORE::LINALG::Matrix<3, 1>& xi, const CORE::LINALG::Matrix<3, 1>& n,
        const CORE::LINALG::Matrix<3, 1>& dir);

   private:
    /// static values for matrix sizes
    static constexpr int num_nodes_ = CORE::FE::num_nodes<celltype>;
    static constexpr int num_dim_ = CORE::FE::dim<celltype>;
    static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
    static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

    CORE::FE::GaussIntegration stiffness_matrix_integration_;
    CORE::FE::GaussIntegration mass_matrix_integration_;

    SolidFormulationHistory<ElementFormulation> history_data_{};

  };  // class SolidEleCalc
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif