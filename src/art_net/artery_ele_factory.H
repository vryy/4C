/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of artery elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef SRC_ART_NET_ARTERY_ELE_FACTORY_H_
#define SRC_ART_NET_ARTERY_ELE_FACTORY_H_

#include "inpar_bio.H"

#include "element.H"

#include "artery_ele_interface.H"

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
      virtual ~ArtNetFactory() { return; }

      //! ProvideImpl
      static ArteryEleInterface* ProvideImpl(DRT::Element::DiscretizationType distype,
          INPAR::ARTDYN::ImplType problem, const std::string& disname);

     private:
      //! define ArteryEle instances dependent on problem
      template <DRT::Element::DiscretizationType distype>
      static ArteryEleInterface* DefineProblemType(
          INPAR::ARTDYN::ImplType problem, const std::string& disname);


    };  // end class ArtNetFactory

  }  // namespace ELEMENTS

}  // namespace DRT



#endif /* SRC_ART_NET_ARTERY_ELE_FACTORY_H_ */
