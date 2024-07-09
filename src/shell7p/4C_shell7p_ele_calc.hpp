/*! \file

\brief Declaration of routines for calculation of shell element
       simple displacement based

\level 3
*/

#ifndef FOUR_C_SHELL7P_ELE_CALC_HPP
#define FOUR_C_SHELL7P_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_shell7p_ele_calc_interface.hpp"
#include "4C_shell7p_ele_calc_lib.hpp"
#include "4C_shell7p_ele_interface_serializable.hpp"

#include <memory>
#include <string>
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace Solid::ELEMENTS
{
  class ParamsInterface;
}  // namespace Solid::ELEMENTS

namespace Discret
{

  namespace ELEMENTS
  {
    template <Core::FE::CellType distype>
    class Shell7pEleCalc : public Shell7pEleCalcInterface, public Shell::Serializable
    {
     public:
      Shell7pEleCalc();

      void setup(Core::Elements::Element& ele, Mat::So3Material& solid_material,
          Input::LineDefinition* linedef, const Solid::ELEMENTS::ShellLockingTypes& locking_types,
          const Solid::ELEMENTS::ShellData& shell_data) override;

      void pack(Core::Communication::PackBuffer& data) const override;

      void unpack(std::vector<char>::size_type& position, const std::vector<char>& data) override;

      void material_post_setup(
          Core::Elements::Element& ele, Mat::So3Material& solid_material) override;

      void evaluate_nonlinear_force_stiffness_mass(Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
          const Core::LinAlg::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params,
          Core::LinAlg::SerialDenseVector* force_vector,
          Core::LinAlg::SerialDenseMatrix* stiffness_matrix,
          Core::LinAlg::SerialDenseMatrix* mass_matrix) override;

      void recover(Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params,
          Solid::ELEMENTS::ParamsInterface& str_interface) override;

      void calculate_stresses_strains(Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const ShellStressIO& stressIO,
          const ShellStrainIO& strainIO, const Core::FE::Discretization& discretization,
          const Core::LinAlg::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) override;

      double calculate_internal_energy(Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
          const Core::LinAlg::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) override;

      void update(Core::Elements::Element& ele, Mat::So3Material& solid_material,
          const Core::FE::Discretization& discretization,
          const Core::LinAlg::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) override;

      void reset_to_last_converged(
          Core::Elements::Element& ele, Mat::So3Material& solid_material) override;

      void vis_data(const std::string& name, std::vector<double>& data) override;

     private:
      //! number of integration points in thickness direction (note: currently they are fixed to 2,
      //! otherwise the element would suffer from nonlinear poisson stiffening)
      const Core::FE::IntegrationPoints1D intpoints_thickness_ =
          Core::FE::IntegrationPoints1D(Core::FE::GaussRule1D::line_2point);

      //! integration points in mid-surface
      Core::FE::IntegrationPoints2D intpoints_midsurface_;

      //! shell data (thickness, SDC, number of ANS parameter)
      Solid::ELEMENTS::ShellData shell_data_ = {};

      //! shell thickness at gauss point in spatial frame
      std::vector<double> cur_thickness_;

    };  // class Shell7pEleCalc
  }     // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
