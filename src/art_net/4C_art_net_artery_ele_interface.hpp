/*----------------------------------------------------------------------*/
/*! \file

 \brief generic interface for implementations of the artery element
        evaluation routines

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_ART_NET_ARTERY_ELE_INTERFACE_HPP
#define FOUR_C_ART_NET_ARTERY_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_art_net_artery.hpp"
#include "4C_art_net_artery_ele_action.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Discret
{
  class Discretization;

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
          Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      //! evaluate service (other quantities apart from rhs and matrix)
      virtual int EvaluateService(Artery* ele, const Arteries::Action action,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      //! evaluate scalar transport (only lin-exp formulation uses this)
      virtual int ScatraEvaluate(Artery* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, std::vector<int>& lm,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
