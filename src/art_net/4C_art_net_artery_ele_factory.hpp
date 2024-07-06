/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of artery elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef FOUR_C_ART_NET_ARTERY_ELE_FACTORY_HPP
#define FOUR_C_ART_NET_ARTERY_ELE_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_art_net_artery_ele_interface.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_inpar_bio.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
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
      static ArteryEleInterface* provide_impl(
          Core::FE::CellType distype, Inpar::ArtDyn::ImplType problem, const std::string& disname);

     private:
      //! define ArteryEle instances dependent on problem
      template <Core::FE::CellType distype>
      static ArteryEleInterface* define_problem_type(
          Inpar::ArtDyn::ImplType problem, const std::string& disname);


    };  // end class ArtNetFactory

  }  // namespace ELEMENTS

}  // namespace Discret



FOUR_C_NAMESPACE_CLOSE

#endif
