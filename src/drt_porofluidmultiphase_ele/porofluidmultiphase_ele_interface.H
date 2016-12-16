/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_interface.H

 \brief generic interface for implementations of the porofluidmultiphase
        boundary element evaluation routines

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
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

    /// Empty destructor
    virtual ~PoroFluidMultiPhaseEleInterface(){return;};

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function; it
      must be defined in PoroFluidMultiPhaseEleCalc.
     */
    virtual int Evaluate(DRT::Element*                 ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization&          discretization,
                         DRT::Element::LocationArray&  la,
                         std::vector<Epetra_SerialDenseMatrix*>& elemat,
                         std::vector<Epetra_SerialDenseVector*>& elevec
                         ) = 0;

  };
} // namespace Elements
} // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_INTERFACE_H_ */
