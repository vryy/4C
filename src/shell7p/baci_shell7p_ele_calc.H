/*! \file

\brief Declaration of routines for calculation of shell element
       simple displacement based

\level 3
*/

#ifndef BACI_SHELL7P_ELE_CALC_H
#define BACI_SHELL7P_ELE_CALC_H

#include "baci_config.H"

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_element_integration_select.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_shell7p_ele_calc_interface.H"
#include "baci_shell7p_ele_calc_lib.H"
#include "baci_shell7p_ele_interface_serializable.H"

#include <memory>
#include <string>
#include <unordered_map>

BACI_NAMESPACE_OPEN

namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace DRT
{

  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class Shell7pEleCalc : public Shell7pEleCalcInterface, public SHELL::Serializable
    {
     public:
      Shell7pEleCalc();

      void Setup(DRT::Element& ele, MAT::So3Material& solid_material,
          INPUT::LineDefinition* linedef, const STR::ELEMENTS::ShellLockingTypes& locking_types,
          const STR::ELEMENTS::ShellData& shell_data) override;

      void Pack(CORE::COMM::PackBuffer& data) const override;

      void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data) override;

      void MaterialPostSetup(DRT::Element& ele, MAT::So3Material& solid_material) override;

      void EvaluateNonlinearForceStiffnessMass(DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
          CORE::LINALG::SerialDenseMatrix* mass_matrix) override;

      void Recover(DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params,
          STR::ELEMENTS::ParamsInterface& str_interface) override;

      void CalculateStressesStrains(DRT::Element& ele, MAT::So3Material& solid_material,
          const ShellStressIO& stressIO, const ShellStrainIO& strainIO,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) override;

      double CalculateInternalEnergy(DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) override;

      void Update(DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) override;

      void ResetToLastConverged(DRT::Element& ele, MAT::So3Material& solid_material) override;

      void VisData(const std::string& name, std::vector<double>& data) override;

     private:
      //! number of integration points in thickness direction (note: currently they are fixed to 2,
      //! otherwise the element would suffer from nonlinear poisson stiffening)
      const CORE::FE::IntegrationPoints1D intpoints_thickness_ =
          CORE::FE::IntegrationPoints1D(CORE::FE::GaussRule1D::line_2point);

      //! integration points in mid-surface
      CORE::FE::IntegrationPoints2D intpoints_midsurface_;

      //! shell data (thickness, SDC, number of ANS parameter)
      STR::ELEMENTS::ShellData shell_data_ = {};

      //! shell thickness at gauss point in spatial frame
      std::vector<double> cur_thickness_;

    };  // class Shell7pEleCalc
  }     // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
