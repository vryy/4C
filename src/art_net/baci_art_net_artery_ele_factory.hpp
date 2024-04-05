/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of artery elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef FOUR_C_ART_NET_ARTERY_ELE_FACTORY_HPP
#define FOUR_C_ART_NET_ARTERY_ELE_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_art_net_artery_ele_interface.hpp"
#include "baci_inpar_bio.hpp"
#include "baci_lib_element.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ArteryEleInterface;

    class ArtNetFactory
    {
     public:
      //! ctor
      ArtNetFactory() { return; }

      //! dtor
      virtual ~ArtNetFactory() = default;
      //! ProvideImpl
      static ArteryEleInterface* ProvideImpl(
          CORE::FE::CellType distype, INPAR::ARTDYN::ImplType problem, const std::string& disname);

     private:
      //! define ArteryEle instances dependent on problem
      template <CORE::FE::CellType distype>
      static ArteryEleInterface* DefineProblemType(
          INPAR::ARTDYN::ImplType problem, const std::string& disname);


    };  // end class ArtNetFactory

  }  // namespace ELEMENTS

}  // namespace DRT



BACI_NAMESPACE_CLOSE

#endif  // ART_NET_ARTERY_ELE_FACTORY_H
