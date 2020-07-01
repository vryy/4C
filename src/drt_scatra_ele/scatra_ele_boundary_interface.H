/*----------------------------------------------------------------------*/
/*! \file

\brief Interface of scatra boundary elements

\level 1

 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_INTERFACE_H
#define SCATRA_ELE_BOUNDARY_INTERFACE_H

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_lib/drt_element.H"

namespace DRT
{
  class Discretization;
  class Condition;
  class FaceElement;

  namespace ELEMENTS
  {
    class TransportBoundary;

    /// Interface base class for ScaTraEleBoundaryCalc
    /*!
      This class exists to provide a common interface for all template
      versions of ScaTraEleBoundaryCalc.
     */
    class ScaTraBoundaryInterface : DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      ScaTraBoundaryInterface() { return; };

      /// Empty destructor
      virtual ~ScaTraBoundaryInterface() { return; };

      virtual int Evaluate(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3) = 0;

      virtual int EvaluateNeumann(DRT::FaceElement* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition,
          DRT::Element::LocationArray& la, Epetra_SerialDenseVector& elevec1,
          const double scalar) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
