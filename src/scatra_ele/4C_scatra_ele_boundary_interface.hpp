/*----------------------------------------------------------------------*/
/*! \file

\brief Interface of scatra boundary elements

\level 1

 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_ELE_BOUNDARY_INTERFACE_HPP
#define FOUR_C_SCATRA_ELE_BOUNDARY_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class FaceElement;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
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

      virtual int evaluate(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
          Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) = 0;

      virtual int evaluate_neumann(Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
          Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
          const double scalar) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
