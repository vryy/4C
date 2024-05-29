/*----------------------------------------------------------------------*/
/*! \file
 \brief factory class providing the implementations of the porofluidmultiphase
        element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_FACTORY_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_FACTORY_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class PoroFluidMultiPhaseEleInterface;

    class PoroFluidMultiPhaseFactory
    {
     public:
      //! ctor
      PoroFluidMultiPhaseFactory() { return; }

      //! dtor
      virtual ~PoroFluidMultiPhaseFactory() = default;
      //! ProvideImpl
      static PoroFluidMultiPhaseEleInterface* ProvideImpl(
          CORE::FE::CellType distype, const int numdofpernode, const std::string& disname);

     private:
      //! define PoroFluidMultiPhaseEle instances dependent on problem
      template <CORE::FE::CellType distype>
      static PoroFluidMultiPhaseEleInterface* define_problem_type(
          const int numdofpernode, const std::string& disname);


    };  // end class PoroFluidMultiPhaseFactory

  }  // namespace ELEMENTS

}  // namespace DRT



FOUR_C_NAMESPACE_CLOSE

#endif
