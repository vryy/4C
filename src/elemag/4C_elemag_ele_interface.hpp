/*--------------------------------------------------------------------------*/
/*! \file

\brief Interface of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_ELE_INTERFACE_HPP
#define FOUR_C_ELEMAG_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_cut_utils.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
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
      virtual int Evaluate(Discret::ELEMENTS::Elemag* ele, Core::FE::Discretization& discretization,
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
