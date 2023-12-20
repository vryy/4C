/*! \file

\brief Interface of Shell 7-Parameter Model Finite Elements

\level 1
*/

#ifndef BACI_SHELL7P_ELE_CALC_INTERFACE_H
#define BACI_SHELL7P_ELE_CALC_INTERFACE_H

#include "baci_config.H"

#include "baci_lib_element.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_mat_so3_material.H"
#include "baci_shell7p_utils.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace STR::MODELEVALUATOR
{
  class GaussPointDataOutputManager;
}

namespace DRT
{
  // forward declaration
  class Discretization;
  class Condition;

  namespace INPUT
  {
    class LineDefinition;
  }

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
          DRT::INPUT::LineDefinition* linedef,
          const STR::ELEMENTS::ShellLockingTypes& locking_types,
          const STR::ELEMENTS::ShellData& shell_data) = 0;

      /*!
       * @brief A setup routine for the materials after the whole input is read before evaluation.
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       *
       */
      virtual void MaterialPostSetup(DRT::Element& ele, MAT::So3Material& solid_material) = 0;

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
      virtual void EvaluateNonlinearForceStiffnessMass(DRT::Element& ele,
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
      virtual void CalculateStressesStrains(DRT::Element& ele, MAT::So3Material& solid_material,
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
      virtual double CalculateInternalEnergy(DRT::Element& ele, MAT::So3Material& solid_material,
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
      virtual void ResetToLastConverged(DRT::Element& ele, MAT::So3Material& solid_material) = 0;

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

BACI_NAMESPACE_CLOSE

#endif
