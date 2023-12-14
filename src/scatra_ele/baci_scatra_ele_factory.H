/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of scatra elements

\level 1

*/
/*--------------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_FACTORY_H
#define BACI_SCATRA_ELE_FACTORY_H

#include "baci_inpar_scatra.H"
#include "baci_lib_element.H"
#include "baci_scatra_ele_interface.H"

BACI_NAMESPACE_OPEN

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
      static ScaTraEleInterface* DefineProblemTypeHDG(INPAR::SCATRA::ImplType problem,
          const int numdofpernode, const int numscal, const std::string& disname);

    };  // end class ScaTraFactory

  }  // namespace ELEMENTS

}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
