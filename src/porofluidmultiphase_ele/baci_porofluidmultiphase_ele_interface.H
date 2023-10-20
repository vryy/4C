/*----------------------------------------------------------------------*/
/*! \file
 \brief generic interface for implementations of the porofluidmultiphase
        boundary element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROFLUIDMULTIPHASE_ELE_INTERFACE_H
#define BACI_POROFLUIDMULTIPHASE_ELE_INTERFACE_H



#include "baci_lib_element.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"

#include <Teuchos_ParameterList.hpp>

#include <vector>

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
    class PoroFluidMultiPhaseEleInterface
    {
     public:
      /// Empty constructor
      PoroFluidMultiPhaseEleInterface() { return; };

      /// Empty destructor
      virtual ~PoroFluidMultiPhaseEleInterface() = default;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in PoroFluidMultiPhaseEleCalc.
       */
      virtual int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          std::vector<CORE::LINALG::SerialDenseMatrix*>& elemat,
          std::vector<CORE::LINALG::SerialDenseVector*>& elevec) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // POROFLUIDMULTIPHASE_ELE_INTERFACE_H
