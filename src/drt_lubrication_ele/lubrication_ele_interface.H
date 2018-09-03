/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_interface.H

\brief Interface of Lubrication element evaluation routine

\level 3

\maintainer Alexander Seitz
*/
/*--------------------------------------------------------------------------*/

#ifndef SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_INTERFACE_H_
#define SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_INTERFACE_H_

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_lib/drt_element.H"

namespace DRT
{
  class Discretization;
  class Element;

  namespace ELEMENTS
  {
    /// Interface base class for LubricationEleCalc
    /*!
      This class exists to provide a common interface for all template
      versions of LubricationEleCalc.
     */
    class LubricationEleInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      LubricationEleInterface() { return; };

      /// Setup element evaluation
      virtual int SetupCalc(DRT::Element* ele, DRT::Discretization& discretization) = 0;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in LubricationEleCalc.
        The Evaluate() method is meant only for the assembling of the
        linearized matrix and the right hand side
       */
      virtual int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra) = 0;

      virtual int EvaluateEHLMon(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra) = 0;

      /*!
        This class does not provide a definition for this function; it
        must be defined in LubricationEleCalc.
        The EvaluateService() method is meant for everything not related to
        the assembling of the linearized matrix and the right hand side
      */
      virtual int EvaluateService(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_INTERFACE_H_ */
