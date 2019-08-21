/*--------------------------------------------------------------------------*/
/*! \file

\brief Factory of acoustic elements

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*--------------------------------------------------------------------------*/

#ifndef ACOU_ELE_FACTORY_H
#define ACOU_ELE_FACTORY_H

#include "../drt_lib/drt_element.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class AcouEleInterface;

    class AcouFactory
    {
     public:
      //! ctor
      AcouFactory() { return; }

      //! dtor
      virtual ~AcouFactory() { return; }

      //! ProvideImpl
      static AcouEleInterface* ProvideImpl(
          DRT::Element::DiscretizationType distype, std::string problem);

     private:
      //! define AcouEle instances dependent on problem
      template <DRT::Element::DiscretizationType distype>
      static AcouEleInterface* DefineProblemType(std::string problem);
    };

  }  // namespace ELEMENTS

}  // namespace DRT

#endif /* ACOU_ELE_FACTORY_H_ */
