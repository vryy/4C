/*! \file

\brief Interface of solid-poro elements

\level 1
*/

#ifndef BACI_SOLID_ELE_PORO_CALC_INTERFACE_H
#define BACI_SOLID_ELE_PORO_CALC_INTERFACE_H

#include "baci_lib_element.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_mat_fluidporo_multiphase.H"
#include "baci_mat_structporo.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_ParameterList.hpp>

namespace DRT
{
  namespace UTILS
  {
    class GaussIntegration;
  }

  namespace ELEMENTS
  {
    class SolidPoroEleCalcInterface
    {
     public:
      virtual ~SolidPoroEleCalcInterface() = default;

      /*!
       * @brief Setup routine for the element
       *
       * @param ele (in) : reference to the element
       * @param linedef (in) : definition of one (possible) line in a dat file
       */
      virtual void PoroSetup(
          MAT::StructPoro& porostructmat, DRT::INPUT::LineDefinition* linedef) = 0;

      /*!
       * @brief Evaluate the force vector, stiffness matrix and mass matrix of the element for
       * poroelasticity problems
       *
       * @param ele (in) : reference to the element
       * @param porostructmat (in) : Solid material of the element
       * @param porofluidmat (in): fluid multiphase material
       * @param kinematictype (in): kinematic type of element
       * @param discretization  (in) : Reference to the discretization
       * @param la (in) : Location vector
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       * @param force_vector (out) : pointer to force vector or nullptr
       * @param stiffness_matrix (out) : pointer to stiffness matrix or nullptr
       */
      virtual void EvaluateNonlinearForceStiffness(const DRT::Element& ele,
          MAT::StructPoro& porostructmat, MAT::FluidPoroMultiPhase& porofluidmat,
          const INPAR::STR::KinemType& kinematictype, const DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Teuchos::ParameterList& params,
          CORE::LINALG::SerialDenseVector* force_vector,
          CORE::LINALG::SerialDenseMatrix* stiffness_matrix) = 0;

      /*!
       * @brief Calculate coupling terms in nonlinear stiffness and internal force for
       * poroelasticity problems
       *
       * @param ele (in) : reference to the element
       * @param porostructmat (in) : Solid material of the element
       * @param porofluidmat (in): fluid multiphase material
       * @param kinematictype (in): kinematic type of element
       * @param discretization  (in) : Reference to the discretization
       * @param la (in) : Location vector
       * @param params (in) : ParameterList for communication between control routine, elements and
       * materials
       * @param stiffness_matrix (out) : pointer to stiffness matrix or nullptr
       */
      virtual void CouplingPoroelast(const DRT::Element& ele, MAT::StructPoro& porostructmat,
          MAT::FluidPoroMultiPhase& porofluidmat, const INPAR::STR::KinemType& kinematictype,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& stiffness_matrix) = 0;

      //! Calculate coupling stress for poroelasticity problems
      virtual void CouplingStress(const DRT::Element& ele,
          const DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif  // SOLID_ELE_PORO_CALC_INTERFACE_H
