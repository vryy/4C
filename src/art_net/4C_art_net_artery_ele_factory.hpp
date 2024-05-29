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
#include "4C_discretization_fem_general_element.hpp"
#include "4C_inpar_bio.hpp"

FOUR_C_NAMESPACE_OPEN

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
      static ArteryEleInterface* define_problem_type(
          INPAR::ARTDYN::ImplType problem, const std::string& disname);


    };  // end class ArtNetFactory

  }  // namespace ELEMENTS

}  // namespace DRT



FOUR_C_NAMESPACE_CLOSE

#endif
