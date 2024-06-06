/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of Lubrication elements

\level 3


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_LUBRICATION_ELE_FACTORY_HPP
#define FOUR_C_LUBRICATION_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    // forward declaration
    class LubricationEleInterface;

    class LubricationFactory
    {
     public:
      //! ctor
      LubricationFactory() { return; }

      //! dtor
      virtual ~LubricationFactory() = default;
      //! ProvideImpl
      static LubricationEleInterface* ProvideImpl(
          Core::FE::CellType distype, const std::string& disname);

     private:
      //! define LubricationEle instances dependent on problem
      template <Core::FE::CellType distype, int probdim>
      static LubricationEleInterface* define_problem_type(const std::string& disname);


    };  // end class LubricationFactory

  }  // namespace ELEMENTS

}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
