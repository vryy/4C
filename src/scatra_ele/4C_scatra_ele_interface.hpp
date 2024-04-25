/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of scatra elements

\level 1

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_INTERFACE_HPP
#define FOUR_C_SCATRA_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_lib_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
  class Element;

  namespace ELEMENTS
  {
    /// Interface base class for ScaTraEleCalc
    /*!
      This class exists to provide a common interface for all template
      versions of ScaTraEleCalc.
     */
    class ScaTraEleInterface
    {
     public:
      /// Virtual destructor.
      virtual ~ScaTraEleInterface() = default;

      /// Setup element evaluation
      virtual int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization) = 0;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in ScatraEleCalc.
       */
      virtual int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;

      virtual int EvaluateService(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;

      virtual int EvaluateOD(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
