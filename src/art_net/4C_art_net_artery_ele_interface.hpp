// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ART_NET_ARTERY_ELE_INTERFACE_HPP
#define FOUR_C_ART_NET_ARTERY_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_art_net_artery.hpp"
#include "4C_art_net_artery_ele_action.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{

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
      virtual int evaluate(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      //! evaluate service (other quantities apart from rhs and matrix)
      virtual int evaluate_service(Artery* ele, const Arteries::Action action,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra,
          Teuchos::RCP<Core::Mat::Material> mat) = 0;

      //! evaluate scalar transport (only lin-exp formulation uses this)
      virtual int scatra_evaluate(Artery* ele, Teuchos::ParameterList& params,
          Core::FE::Discretization& discretization, std::vector<int>& lm,
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
