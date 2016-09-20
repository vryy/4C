/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_interface.H

 \brief generic interface for implementations of the porofluidmultiphase
        boundary element evaluation routines

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_INTERFACE_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_INTERFACE_H_



#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_lib/drt_element.H"

namespace DRT
{
  class Discretization;
  class Element;

namespace ELEMENTS
{

  /// Interface base class for PoroFluidMultiPhaseEleCalc
  /*!
    This class exists to provide a common interface for all template
    versions of PoroFluidMultiPhaseEleCalc.
   */
  class PoroFluidMultiPhaseEleInterface : public DRT::SingletonDestruction
  {
  public:

    /// Empty constructor
    PoroFluidMultiPhaseEleInterface(){return;};

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function; it
      must be defined in PoroFluidMultiPhaseEleCalc.
      The Evaluate() method is meant only for the assembling of the
      linearized matrix and the right hand side
     */
    virtual int Evaluate(DRT::Element*                 ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization&          discretization,
                         DRT::Element::LocationArray&  la,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra) = 0;

    /*!
      This class does not provide a definition for this function; it
      must be defined in PoroFluidMultiPhaseEleCalc.
      The EvaluateService() method is meant for everything not related to
      the assembling of the linearized matrix and the right hand side
    */
    virtual int EvaluateService(DRT::Element*                 ele,
                                Teuchos::ParameterList&       params,
                                DRT::Discretization&          discretization,
                                DRT::Element::LocationArray&  la,
                                Epetra_SerialDenseMatrix&     elemat1_epetra,
                                Epetra_SerialDenseMatrix&     elemat2_epetra,
                                Epetra_SerialDenseVector&     elevec1_epetra,
                                Epetra_SerialDenseVector&     elevec2_epetra,
                                Epetra_SerialDenseVector&     elevec3_epetra) = 0;

  };
} // namespace Elements
} // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_INTERFACE_H_ */
