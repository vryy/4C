/*--------------------------------------------------------------------------*/
/*! \file

\brief Interface of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_ELEMAG_ELE_INTERFACE_HPP
#define FOUR_C_ELEMAG_ELE_INTERFACE_HPP

#include "baci_config.hpp"

#include "baci_cut_utils.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace MAT
{
  class Material;
}

namespace DRT
{
  class Discretization;

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
      virtual int Evaluate(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag = false) = 0;

      /// Integrate shape function
      /*!
        This class does not provide a definition for this function; it
        must be defined in ElemagEleCalc.
       */
      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Elemag* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1) = 0;
    };


  }  // namespace ELEMENTS

}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_ELE_INTERFACE_H
