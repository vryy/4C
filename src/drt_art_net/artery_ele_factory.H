/*--------------------------------------------------------------------------*/
/*!

\brief Factory of artery elements

\level 3

\maintainer Johannes Kremheller
*/
/*--------------------------------------------------------------------------*/


#ifndef SRC_DRT_ART_NET_ARTERY_ELE_FACTORY_H_
#define SRC_DRT_ART_NET_ARTERY_ELE_FACTORY_H_

#include "../drt_inpar/inpar_bio.H"

#include "../drt_lib/drt_element.H"

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



#endif /* SRC_DRT_ART_NET_ARTERY_ELE_FACTORY_H_ */
