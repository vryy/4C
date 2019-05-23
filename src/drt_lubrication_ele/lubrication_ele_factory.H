/*--------------------------------------------------------------------------*/
/*!

\brief Factory of Lubrication elements

\level 3

\maintainer Mostafa Faraji

*/
/*--------------------------------------------------------------------------*/

#ifndef SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_FACTORY_H_
#define SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_FACTORY_H_

#include "../drt_lib/drt_element.H"

namespace DRT
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
      virtual ~LubricationFactory() { return; }

      //! ProvideImpl
      static LubricationEleInterface* ProvideImpl(
          DRT::Element::DiscretizationType distype, const std::string& disname);

     private:
      //! define LubricationEle instances dependent on problem
      template <DRT::Element::DiscretizationType distype, int probdim>
      static LubricationEleInterface* DefineProblemType(const std::string& disname);


    };  // end class LubricationFactory

  }  // namespace ELEMENTS

}  // namespace DRT


#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_FACTORY_H_ */
