/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class defining the fluid boundary element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_INTERFACE_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class FluidBoundary;

    /// Interface base class for FluidEleBoundaryCalc
    /*!
      This class exists to provide a common interface for all template
      versions of FluidBoundaryCalc.
     */
    //  class FluidBoundaryImplInterface
    class FluidBoundaryInterface
    {
     public:
      /// Empty constructor
      FluidBoundaryInterface(){};

      /// Empty destructor
      virtual ~FluidBoundaryInterface() = default;

      virtual void evaluate_action(DRT::ELEMENTS::FluidBoundary* ele1,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
          CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) = 0;

      virtual int evaluate_neumann(DRT::ELEMENTS::FluidBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          CORE::Conditions::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseMatrix* elemat1) = 0;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
