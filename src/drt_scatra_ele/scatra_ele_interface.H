/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_interface.H

\brief Factory of scatra elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_INTERFACE_H
#define SCATRA_ELE_INTERFACE_H

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
  class Discretization;

namespace ELEMENTS
{
  class Transport;

  /// Interface base class for ScaTraEleCalc
  /*!
    This class exists to provide a common interface for all template
    versions of FluidEleCalc. The only function this class actually defines
    is Ele, which returns a pointer to the appropriate version of FluidEleCalc.
   */
  class ScaTraEleInterface : public DRT::SingletonDestruction
  {
  public:

    /// Empty constructor
    ScaTraEleInterface() {}

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function; it
      must be defined in ScatraEleCalc.
     */
    virtual int Evaluate(DRT::ELEMENTS::Transport*     ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization &         discretization,
                         const std::vector<int> &      lm,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra) = 0;

    virtual int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization &         discretization,
                         const std::vector<int> &      lm,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra) = 0;

  };

} //namespace Elements
} //namespace DRT

#endif

