/*----------------------------------------------------------------------*/
/*! \file

 \brief generic interface for implementations of the artery element
        evaluation routines

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_ART_NET_ARTERY_ELE_INTERFACE_H
#define BACI_ART_NET_ARTERY_ELE_INTERFACE_H

#include "baci_art_net_artery.H"
#include "baci_art_net_artery_ele_action.H"
#include "baci_lib_element.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"

#include <Teuchos_ParameterList.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

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
    class ArteryEleInterface
    {
     public:
      /// Virtual destructor
      virtual ~ArteryEleInterface() = default;

      //! evaluate the element
      virtual int Evaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      //! evaluate service (other quantities apart from rhs and matrix)
      virtual int EvaluateService(Artery* ele, const ARTERY::Action action,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;

      //! evaluate scalar transport (only lin-exp formulation uses this)
      virtual int ScatraEvaluate(Artery* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT



BACI_NAMESPACE_CLOSE

#endif  // ART_NET_ARTERY_ELE_INTERFACE_H
