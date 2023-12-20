/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class going from the generic evaluation routines to the ones
  templated by the element shape and specialization

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_FACTORY_H
#define BACI_FLUID_ELE_FACTORY_H

#include "baci_config.H"

#include "baci_lib_element.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidEleInterface;

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief to do
     *
     *
     * \date March, 2012
     */
    /*--------------------------------------------------------------------------*/
    class FluidFactory
    {
     public:
      //! ctor
      FluidFactory() { return; }

      //! dtor
      virtual ~FluidFactory() = default;
      //! ProvideImpl
      static FluidEleInterface* ProvideImpl(CORE::FE::CellType distype, std::string problem);

      //! special ProvideImpl for XFEM problems to reduce created template combinations
      static FluidEleInterface* ProvideImplXFEM(CORE::FE::CellType distype, std::string problem);

     private:
      //! define FluidEle instances dependent on problem
      template <CORE::FE::CellType distype>
      static FluidEleInterface* DefineProblemType(std::string problem);

      //! special DefineProblemTypeXFEM for XFEM problems
      template <CORE::FE::CellType distype>
      static FluidEleInterface* DefineProblemTypeXFEM(std::string problem);

    };  // end class FluidFactory

  }  // namespace ELEMENTS

}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
