/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class defining the fluid boundary element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_BOUNDARY_INTERFACE_HPP
#define FOUR_C_FLUID_ELE_BOUNDARY_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
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

      virtual void evaluate_action(Discret::ELEMENTS::FluidBoundary* ele1,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3) = 0;

      virtual int evaluate_neumann(Discret::ELEMENTS::FluidBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Conditions::Condition& condition, std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseMatrix* elemat1) = 0;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
