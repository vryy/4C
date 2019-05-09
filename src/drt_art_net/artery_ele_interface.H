/*----------------------------------------------------------------------*/
/*!

 \brief generic interface for implementations of the artery element
        evaluation routines

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_ART_NET_ARTERY_ELE_INTERFACE_H_
#define SRC_DRT_ART_NET_ARTERY_ELE_INTERFACE_H_

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_lib/drt_element.H"
#include "artery.H"
#include "artery_ele_action.H"

namespace DRT
{
  class Discretization;
  class Element;

  namespace ELEMENTS
  {
    /// Interface base class for ArteryEleCalc
    /*!
      This class exists to provide a common interface for all template
      versions of ArteryEleCalc.
     */
    class ArteryEleInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      ArteryEleInterface() { return; };

      /// Empty destructor
      virtual ~ArteryEleInterface() { return; };

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in ScatraEleCalc.
       */
      virtual int Evaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual int EvaluateService(Artery* ele, const ARTERY::Action action,
          Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      virtual int ScatraEvaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT



#endif /* SRC_DRT_ART_NET_ARTERY_ELE_INTERFACE_H_ */
