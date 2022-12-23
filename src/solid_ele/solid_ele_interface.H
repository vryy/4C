/*----------------------------------------------------------------------*/
/*! \file

\brief Interface of solid elements

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_INTERFACE_H
#define SOLID_ELE_INTERFACE_H

#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "dserror.H"
#include "element.H"

namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

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

    class SolidEleInterface
    {
     public:
      virtual ~SolidEleInterface() = default;

      /*!
       * @brief Setup routine for the element
       *
       * @param ele (in) : reference to the element
       * @param linedef (in) : 778
       */
      virtual void Setup(const DRT::ELEMENTS::Solid& ele, DRT::INPUT::LineDefinition* linedef) = 0;

      /*!
       * @brief A setup routine for the materials after the whole input is read before evaluation.
       *
       * @param ele (in) : reference to the element
       */
      virtual void MaterialPostSetup(const DRT::ELEMENTS::Solid& ele) = 0;

      /*!
       * @brief Evaluate the force vector, stiffness matrix and mass matrix of the element
       *
       * @param ele (in) : reference to the element
       * @param discretization  (in) : Reference to the discretization
       * @param lm (in) : Location vector
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       * @param force_vector (out) : pointer to force vector or nullptr
       * @param stiffness_matrix (out) : pointer to stiffness matrix or nullptr
       * @param mass_matrix (out) : pointer to mass matrix or nullptr
       */
      virtual void EvaluateNonlinearForceStiffnessMass(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseVector* force_vector,
          Epetra_SerialDenseMatrix* stiffness_matrix, Epetra_SerialDenseMatrix* mass_matrix) = 0;

      virtual void Recover(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;

      virtual void CalculateStress(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;

      virtual void PostProcessStressStrain(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;

      /*!
       * @brief The update routine of the element
       *
       * @param ele  (in) : Reference to the solid element
       * @param discretization (in) : Reference to the discretization
       * @param lm (in) : The location array of the owned dofs
       * @param params (in/[out]) : A ParameterList to pass values from the time integrator to the
       * elements/materials
       */
      virtual void Update(const DRT::ELEMENTS::Solid& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;

      /*!
       * @brief Initializes the Gauss point data output
       *
       * The element must register the number of used Gauss points, the names and the the rank of
       * each quantity. It is also supposed to ask the material about their quantities.
       *
       * @param ele (in) : Reference to the solid element
       */
      virtual void InitializeGaussPointDataOutput(const DRT::ELEMENTS::Solid& ele) const = 0;

      /*!
       * @brief Evaluation of the Gauss point data output
       *
       * The element must go over all quantities and check whether it has some of the values. In
       * case there are some values, they must be provided at the specified location (element
       * center, nodes or gauss points)
       *
       * @param ele (in) : Reference to the solid element
       */
      virtual void EvaluateGaussPointDataOutput(const DRT::ELEMENTS::Solid& ele) const = 0;


    };  // class SolidEleInterface
  }     // namespace ELEMENTS

}  // namespace DRT

#endif  // SOLID_ELE_INTERFACE_H
