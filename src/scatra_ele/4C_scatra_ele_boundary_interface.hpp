/*----------------------------------------------------------------------*/
/*! \file

\brief Interface of scatra boundary elements

\level 1

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_INTERFACE_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::Elements
{
  class FaceElement;
}

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class TransportBoundary;

    /// Interface base class for ScaTraEleBoundaryCalc
    /*!
      This class exists to provide a common interface for all template
      versions of ScaTraEleBoundaryCalc.
     */
    class ScaTraBoundaryInterface
    {
     public:
      /// Virtual destructor.
      virtual ~ScaTraBoundaryInterface() = default;

      virtual int Evaluate(CORE::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
          CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) = 0;

      virtual int evaluate_neumann(CORE::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
          CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseVector& elevec1,
          const double scalar) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
