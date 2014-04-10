/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_interface.H

\brief Interface of scatra elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_INTERFACE_H
#define SCATRA_ELE_BOUNDARY_INTERFACE_H

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
  class Discretization;
  class Condition;

namespace ELEMENTS
{
class TransportBoundary;

  /// Interface base class for ScaTraEleBoundaryCalc
  /*!
    This class exists to provide a common interface for all template
    versions of FluidBoundaryCalc.
   */
  //  class FluidBoundaryImplInterface : DRT::SingletonDestruction
  class ScaTraBoundaryInterface : DRT::SingletonDestruction
  {

  public:

    /// Empty constructor
    ScaTraBoundaryInterface() {};

    /// Empty destructor
    virtual ~ScaTraBoundaryInterface() {};

    virtual int EvaluateAction(
      DRT::ELEMENTS::TransportBoundary* ele1,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      std::vector<int>&             lm,
      Epetra_SerialDenseMatrix&     elemat1,
      Epetra_SerialDenseMatrix&     elemat2,
      Epetra_SerialDenseVector&     elevec1,
      Epetra_SerialDenseVector&     elevec2,
      Epetra_SerialDenseVector&     elevec3) = 0;

    virtual int EvaluateNeumann(
      DRT::ELEMENTS::TransportBoundary* ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Condition&               condition,
      std::vector<int>&             lm,
      Epetra_SerialDenseVector&     elevec1_epetra
    ) = 0;

  };

} //namespace Elements
} //namespace DRT

#endif
