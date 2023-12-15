/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of artery elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef BACI_ART_NET_ARTERY_ELE_FACTORY_H
#define BACI_ART_NET_ARTERY_ELE_FACTORY_H

#include "baci_art_net_artery_ele_interface.H"
#include "baci_inpar_bio.H"
#include "baci_lib_element.H"

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
