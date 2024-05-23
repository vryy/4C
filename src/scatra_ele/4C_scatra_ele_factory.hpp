/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of scatra elements

\level 1

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_FACTORY_HPP
#define FOUR_C_SCATRA_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_lib_element.hpp"
#include "4C_scatra_ele_interface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraFactory
    {
     public:
      //! ctor
      ScaTraFactory() = default;

      //! dtor
      virtual ~ScaTraFactory() = default;

      //! ProvideImpl
      static ScaTraEleInterface* ProvideImpl(CORE::FE::CellType distype,
          INPAR::SCATRA::ImplType problem, const int numdofpernode, const int numscal,
          const std::string& disname);

      //! ProvideImplHDG
      static ScaTraEleInterface* ProvideImplHDG(CORE::FE::CellType distype,
          INPAR::SCATRA::ImplType problem, const int numdofpernode, const int numscal,
          const std::string& disname);

     private:
      //! define ScatraEle instances dependent on problem
      template <CORE::FE::CellType distype, int probdim>
      static ScaTraEleInterface* DefineProblemType(INPAR::SCATRA::ImplType problem,
          const int numdofpernode, const int numscal, const std::string& disname);

      //! define ScatraEle instances dependent on problem
      template <CORE::FE::CellType distype, int probdim>
      static ScaTraEleInterface* define_problem_type_hdg(INPAR::SCATRA::ImplType problem,
          const int numdofpernode, const int numscal, const std::string& disname);

    };  // end class ScaTraFactory

  }  // namespace ELEMENTS

}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
