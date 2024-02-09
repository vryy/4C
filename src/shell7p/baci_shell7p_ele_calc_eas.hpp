/*! \file

\brief Declaration of routines for calculation of shell element with EAS element technology

\level 3
*/
#ifndef BACI_SHELL7P_ELE_CALC_EAS_HPP
#define BACI_SHELL7P_ELE_CALC_EAS_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_gausspoints.hpp"
#include "baci_lib_element_integration_select.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_shell7p_ele_calc_interface.hpp"
#include "baci_shell7p_ele_calc_lib.hpp"
#include "baci_shell7p_ele_interface_serializable.hpp"

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
    /*!
    \brief A struct holding the iteration data used within the Enhanced Assumed Strain element
    technology
    */
    struct ShellEASIterationData
    {
      // EAS portion of internal forces Rtilde (dimension neas x 1)
      CORE::LINALG::SerialDenseMatrix RTilde_;
      // inverse EAS matrix Dtilde (dimension neas x neas)
      CORE::LINALG::SerialDenseMatrix invDTilde_;
      // transposed EAS matrix  L^T (dimension neas x numdofperelement)
      CORE::LINALG::SerialDenseMatrix transL_;
      // alpha (dimension neas x 1)
      CORE::LINALG::SerialDenseMatrix alpha_;
    };

    template <CORE::FE::CellType distype>
    class Shell7pEleCalcEas : public Shell7pEleCalcInterface, public SHELL::Serializable
    {
     public:
      Shell7pEleCalcEas();

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
          STR::ELEMENTS::ParamsInterface& interface_ptr) override;

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
      //! EAS matrices and vectors to be stored between iterations
      DRT::ELEMENTS::ShellEASIterationData eas_iteration_data_ = {};

      //! number of integration points in thickness direction (note: currently they are fixed to 2,
      //! otherwise the element would suffer from nonlinear poisson stiffening)
      const CORE::FE::IntegrationPoints1D intpoints_thickness_ =
          CORE::FE::IntegrationPoints1D(CORE::FE::GaussRule1D::line_2point);

      //! integration points on shell mid-surface
      CORE::FE::IntegrationPoints2D intpoints_midsurface_;

      //! shell data (thickness, SDC, number of ANS parameter)
      STR::ELEMENTS::ShellData shell_data_ = {};

      //! different Locking types that we want to alleviate
      STR::ELEMENTS::ShellLockingTypes locking_types_ = {};

      //! old step length
      double old_step_length_;

      //! shell thickness at gauss point in spatial frame
      std::vector<double> cur_thickness_;

    };  // class Shell7pEleCalcEas
  }     // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
