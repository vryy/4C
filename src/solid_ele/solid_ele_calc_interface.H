/*----------------------------------------------------------------------*/
/*! \file

\brief Interface of solid elements

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_CALC_INTERFACE_H
#define SOLID_ELE_CALC_INTERFACE_H

#include <Teuchos_ParameterList.hpp>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "utils_exceptions.H"
#include "lib_element.H"
#include "mat_so3_material.H"

namespace STR::MODELEVALUATOR
{
  class GaussPointDataOutputManager;
}

namespace DRT
{
  // forward declaration
  class Discretization;
  class Condition;
  namespace UTILS
  {
    class GaussIntegration;
  }
  namespace INPUT
  {
    class LineDefinition;
  }

  namespace ELEMENTS
  {
    // forward declaration
    class Solid;

    struct StressIO
    {
      INPAR::STR::StressType type;
      std::vector<char>& mutable_data;
    };

    struct StrainIO
    {
      INPAR::STR::StrainType type;
      std::vector<char>& mutable_data;
    };

    class SolidEleCalcInterface
    {
     public:
      virtual ~SolidEleCalcInterface() = default;

      /*!
       * @brief Setup routine for the element
       *
       * @param solid_material (in) : Solid material of the element
       * @param linedef (in) : 778
       */
      virtual void Setup(MAT::So3Material& solid_material, DRT::INPUT::LineDefinition* linedef) = 0;

      /*!
       * @brief A setup routine for the materials after the whole input is read before evaluation.
       *
       * @param ele (in) : reference to the element
       * @param solid_material (in) : Solid material of the element
       *
       */
      virtual void MaterialPostSetup(const DRT::Element& ele, MAT::So3Material& solid_material) = 0;

      /*!
       * @brief Evaluate the force vector, stiffness matrix and mass matrix of the element
       *
       * @param ele (in) : reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param discretization  (in) : Reference to the discretization
       * @param lm (in) : Location vector
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       * @param force_vector (out) : pointer to force vector or nullptr
       * @param stiffness_matrix (out) : pointer to stiffness matrix or nullptr
       * @param mass_matrix (out) : pointer to mass matrix or nullptr
       */
      virtual void EvaluateNonlinearForceStiffnessMass(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseVector* force_vector, Epetra_SerialDenseMatrix* stiffness_matrix,
          Epetra_SerialDenseMatrix* mass_matrix) = 0;

      /*!
       * @brief Evaluate the force vector, stiffness matrix and mass matrix of the element with GEMM
       *
       * @param ele (in) : reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param discretization  (in) : Reference to the discretization
       * @param lm (in) : Location vector
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       * @param force_vector (out) : pointer to force vector or nullptr
       * @param stiffness_matrix (out) : pointer to stiffness matrix or nullptr
       * @param mass_matrix (out) : pointer to mass matrix or nullptr
       */
      virtual void EvaluateNonlinearForceStiffnessMassGEMM(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseVector* force_vector, Epetra_SerialDenseMatrix* stiffness_matrix,
          Epetra_SerialDenseMatrix* mass_matrix) = 0;

      virtual void Recover(const DRT::Element& ele, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params) = 0;

      virtual void CalculateStress(const DRT::Element& ele, MAT::So3Material& solid_material,
          const StressIO& stressIO, const StrainIO& strainIO,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;

      virtual double CalculateInternalEnergy(const DRT::Element& ele,
          MAT::So3Material& solid_material, const DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params) = 0;

      /*!
       * @brief The update routine of the element
       *
       * @param ele  (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param discretization (in) : Reference to the discretization
       * @param lm (in) : The location array of the owned dofs
       * @param params (in/[out]) : A ParameterList to pass values from the time integrator to
       * the elements/materials
       */
      virtual void Update(const DRT::Element& ele, MAT::So3Material& solid_material,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;

      /*!
       * @brief Initializes the Gauss point data output
       *
       * The element must register the number of used Gauss points, the names and the the rank of
       * each quantity. It is also supposed to ask the material about their quantities.
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param gp_data_output_manager (in) : Gauss point data output manager
       *                                    (only for new structure time integration)
       */
      virtual void InitializeGaussPointDataOutput(const DRT::Element& ele,
          const MAT::So3Material& solid_material,
          STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const = 0;

      /*!
       * @brief Evaluation of the Gauss point data output
       *
       * The element must go over all quantities and check whether it has some of the values. In
       * case there are some values, they must be provided at the specified location (element
       * center, nodes or gauss points)
       *
       * @param ele (in) : Reference to the element
       * @param solid_material (in) : Solid material of the element
       * @param gp_data_output_manager (in) : Gauss point data output manager
       *                                    (only for new structure time integration)
       */
      virtual void EvaluateGaussPointDataOutput(const DRT::Element& ele,
          const MAT::So3Material& solid_material,
          STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const = 0;

      virtual void ResetAll(const DRT::Element& ele, MAT::So3Material& solid_material) = 0;

      virtual void ResetToLastConverged(
          const DRT::Element& ele, MAT::So3Material& solid_material) = 0;

    };  // class SolidEleInterface
  }     // namespace ELEMENTS

}  // namespace DRT

#endif  // SOLID_ELE_CALC_INTERFACE_H
