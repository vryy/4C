/*----------------------------------------------------------------------*/
/*!
 \file volmortar_utilis.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef VOLMORTAR_UTILIS_H_
#define VOLMORTAR_UTILIS_H_

namespace DRT
{
  class Element;
  class Discretization;
}

namespace VOLMORTAR
{
  class VolMortarCoupl;

  namespace UTILS
  {
    class DefaultMaterialStrategy
    {
    public:
      DefaultMaterialStrategy(){};

      virtual ~DefaultMaterialStrategy(){};

      virtual void AssignMaterialBToA(const VOLMORTAR::VolMortarCoupl* volmortar,
                                  DRT::Element* Aele,
                                  const std::vector<int>& Bids,
                                  Teuchos::RCP<DRT::Discretization> disA,
                                  Teuchos::RCP<DRT::Discretization> disB);

      virtual void AssignMaterialAToB(const VOLMORTAR::VolMortarCoupl* volmortar,
                                  DRT::Element* Bele,
                                  const std::vector<int>& Aids,
                                  Teuchos::RCP<DRT::Discretization> disA,
                                  Teuchos::RCP<DRT::Discretization> disB);
    };
  }
}


#endif /* VOLMORTAR_UTILIS_H_ */
