/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_boundary_interface.H

 \brief generic interface for implementations of the porofluidmultiphase
        boundary element evaluation routines

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_INTERFACE_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_INTERFACE_H_


#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_lib/drt_element.H"

namespace DRT
{
  class Discretization;
  class Condition;
  class FaceElement;

namespace ELEMENTS
{

  /// Interface base class for PoroFluidMultiPhaseEleBoundaryCalc
  /*!
    This class exists to provide a common interface for all template
    versions of PoroFluidMultiPhaseEleBoundaryCalc.
   */
  class PoroFluidMultiPhaseEleBoundaryInterface : DRT::SingletonDestruction
  {

  public:

    /// Empty constructor
    PoroFluidMultiPhaseEleBoundaryInterface(){return;};

    /// Empty destructor
    virtual ~PoroFluidMultiPhaseEleBoundaryInterface(){return;};

    /// evaluate the element
    virtual int Evaluate(
      DRT::FaceElement*                   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      DRT::Element::LocationArray&        la,
      Epetra_SerialDenseMatrix&           elemat1,
      Epetra_SerialDenseMatrix&           elemat2,
      Epetra_SerialDenseVector&           elevec1,
      Epetra_SerialDenseVector&           elevec2,
      Epetra_SerialDenseVector&           elevec3
    ) = 0;

  };
} // namespace Elements
} // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_BOUNDARY_INTERFACE_H_ */
