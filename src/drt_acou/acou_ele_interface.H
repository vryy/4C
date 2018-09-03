/*--------------------------------------------------------------------------*/
/*!
\file acou_ele_interface.H

\brief Interface of acoustic elements

<pre>
\level 2

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef ACOU_ELE_INTERFACE_H
#define ACOU_ELE_INTERFACE_H

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
    class Acou;

    class AcouEleInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty Constructor
      AcouEleInterface() {}

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in AcouEleCalc.
       */
      virtual int Evaluate(DRT::ELEMENTS::Acou* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
          bool offdiag = false) = 0;

      /// Integrate shape function
      /*!
        This class does not provide a definition for this function; it
        must be defined in AcouEleCalc.
       */
      virtual int IntegrateShapeFunction(DRT::ELEMENTS::Acou* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1) = 0;
    };


  }  // namespace ELEMENTS

}  // namespace DRT


#endif /* ACOU_ELE_INTERFACE_H */
