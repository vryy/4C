/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class defining the fluid boundary element

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_BOUNDARY_INTERFACE_HPP
#define BACI_FLUID_ELE_BOUNDARY_INTERFACE_HPP

#include "baci_config.hpp"

#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Condition;

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

      virtual void EvaluateAction(DRT::ELEMENTS::FluidBoundary* ele1,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
          CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) = 0;

      virtual int EvaluateNeumann(DRT::ELEMENTS::FluidBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseMatrix* elemat1) = 0;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
