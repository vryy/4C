// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ELEMAG_ELE_INTERFACE_HPP
#define FOUR_C_ELEMAG_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_cut_utils.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class Material;
}

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Elemag;

    class ElemagEleInterface
    {
     public:
      /// Virtual destructor
      virtual ~ElemagEleInterface() = default;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in ElemagEleCalc.
       */
      virtual int evaluate(Discret::ELEMENTS::Elemag* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<Core::Mat::Material>& mat, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) = 0;

      /// Integrate shape function
      /*!
        This class does not provide a definition for this function; it
        must be defined in ElemagEleCalc.
       */
      virtual int integrate_shape_function(Discret::ELEMENTS::Elemag* ele,
          Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Core::LinAlg::SerialDenseVector& elevec1) = 0;
    };


  }  // namespace ELEMENTS

}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
