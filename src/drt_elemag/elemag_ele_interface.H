/*--------------------------------------------------------------------------*/
/*! \file

\brief Interface of electromagnetic elements

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef ELEMAG_ELE_INTERFACE_H
#define ELEMAG_ELE_INTERFACE_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"

#include "../drt_cut/cut_utils.H"

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

    class ElemagEleInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty Constructor
      ElemagEleInterface() {}

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in ElemagEleCalc.
       */
      virtual int Evaluate(DRT::ELEMENTS::Elemag* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false) = 0;

      /// Integrate shape function
      /*!
        This class does not provide a definition for this function; it
        must be defined in ElemagEleCalc.
       */
      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Elemag* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1) = 0;
    };


  }  // namespace ELEMENTS

}  // namespace DRT


#endif /* ELEMAG_ELE_INTERFACE_H */
