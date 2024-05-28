/*! \file

\brief Interface of Shell 7-Parameter Model Finite Elements

\level 1
*/

#ifndef FOUR_C_SHELL7P_ELE_CALC_INTERFACE_HPP
#define FOUR_C_SHELL7P_ELE_CALC_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_lib_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_shell7p_utils.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace STR::MODELEVALUATOR
{
  class GaussPointDataOutputManager;
}

namespace INPUT
{
  class LineDefinition;
}
namespace DRT
{
  // forward declaration
  class Discretization;

  namespace ELEMENTS
  {
    // forward declaration
    class Shell7p;
    struct ShellStressIO
    {
      INPAR::STR::StressType type;
      std::vector<char>& mutable_data;
    };

    struct ShellStrainIO
    {
      INPAR::STR::StrainType type;
      std::vector<char>& mutable_data;
    };

    class Shell7pEleCalcInterface
    {
     public:
      //! destructor
      virtual ~Shell7pEleCalcInterface() = default;

      /*!
       * @brief Setup routine for the element
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param locking_types (in) : Shell locking types
       * @param shell_data (in) : Shell data (thickness, number of ANS, SDC)
       * @param linedef (in) : Input line of the corresponding element
       */
      virtual void Setup(DRT::Element& ele, MAT::So3Material& solid_material,
          INPUT::LineDefinition* linedef, const STR::ELEMENTS::ShellLockingTypes& locking_types,
          const STR::ELEMENTS::ShellData& shell_data) = 0;

      /*!
       * @brief A setup routine for the materials after the whole input is read before evaluation.
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       *
       */
      virtual void material_post_setup(DRT::Element& ele, MAT::So3Material& solid_material) = 0;

      /*!
       * @brief Evaluate the force vector, stiffness matrix and mass matrix of the element
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param discretization  (in) : Reference to the discretization
       * @param nodal_directors  (in) : Nodal directors in reference frame
       * @param dof_index_array (in) : Location vector of the owned dofs
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       * @param force_vector (out) : Pointer to force vector or nullptr
       * @param stiffness_matrix (out) : Pointer to stiffness matrix or nullptr
       * @param mass_matrix (out) : Pointer to mass matrix or nullptr
       */
      virtual void evaluate_nonlinear_force_stiffness_mass(DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
          CORE::LINALG::SerialDenseMatrix* mass_matrix) = 0;

      /*!
       * @brief Recover condensed EAS variables
       *
       * @param ele  (in) : Reference to the element
       * @param discretization  (in) : Reference to the discretization
       * @param dof_index_array (in) : Location vector of the owned dofs
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       */
      virtual void Recover(DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params,
          STR::ELEMENTS::ParamsInterface& str_interface) = 0;

      /*!
       * @brief Evaluates the stresses and strains
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param strainIO  (in) : Strain input/output type
       * @param stressIO (in) : Stress input/output type
       * @param discretization  (in) : Reference to the discretization
       * @param nodal_directors  (in) : Nodal directors in reference frame
       * @param dof_index_array (in) : Location vector of the owned dofs
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       */
      virtual void calculate_stresses_strains(DRT::Element& ele, MAT::So3Material& solid_material,
          const ShellStressIO& stressIO, const ShellStrainIO& strainIO,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) = 0;

      /*!
       * @brief Calculates the internal energy
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param discretization (in) : Reference to the discretization
       * @param nodal_directors  (in) : Nodal directors in reference frame
       * @param dof_index_array (in) : Location vector of the owned dofs
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       */
      virtual double calculate_internal_energy(DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) = 0;

      /*!
       * @brief The update routine of the element
       *
       * @param ele  (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param discretization (in) : Reference to the discretization
       * @param nodal_directors  (in) : Nodal directors in reference frame
       * @param dof_index_array (in) : The location array of the owned dofs
       * @param params (in/[out]) : A ParameterList to pass values from the time integrator to
       * the elements/materials
       */
      virtual void Update(DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization,
          const CORE::LINALG::SerialDenseMatrix& nodal_directors,
          const std::vector<int>& dof_index_array, Teuchos::ParameterList& params) = 0;

      /*!
       * @brief Reset time step of material (for time adaptivity)
       *
       * @param ele  (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       */
      virtual void reset_to_last_converged(DRT::Element& ele, MAT::So3Material& solid_material) = 0;

      /*!
       * \brief Query data to be visualized using BINIO of a given name
       *
       * @param name (in):   Name of data that is currently processed for visualization
       * @param data (out):  data to be filled by element if it recognizes the name
       */
      virtual void VisData(const std::string& name, std::vector<double>& data) = 0;
    };  // class Shell7pEleInterface
  }     // namespace ELEMENTS

}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
