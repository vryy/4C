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

/*---------------------------------------------------------------------*
 | forward declarations                                    vuong 09/14 |
 *---------------------------------------------------------------------*/
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
    /// Helper class for assigning materials for volumetric coupling of non conforming meshes
    /*!
     When coupling two overlapping discretizations, most often one discretization needs access
     to the corresponding element/material on the other side. For conforming meshes this is straight forward
     as there is one unique element on the other side and therefore one unique material, which can
     be accessed. However, for non conforming meshes there are potentially several elements overlapping.
     Therefore, some rule for assigning materials is needed. This class is meant to do that.
     It gets the element to which it shall assign a material and a vector of IDs of the overlapping elements
     of the other discretization.

     The default strategy will just assign the material of the first element in the vector to the other element.
     This is fine for constant material properties, for instance. If there is furhter work to be done it is meant
     to derive from this class (see e.g. TSIMaterialStrategy).

     \author vuong 10/14
     */

    class DefaultMaterialStrategy
    {
    public:
      //! standard constructor
      DefaultMaterialStrategy(){};

      //! virtual destructor
      virtual ~DefaultMaterialStrategy(){};

      //! assign material of discretization B
      virtual void AssignMaterialBToA(const VOLMORTAR::VolMortarCoupl* volmortar,
                                  DRT::Element* Aele,
                                  const std::vector<int>& Bids,
                                  Teuchos::RCP<DRT::Discretization> disA,
                                  Teuchos::RCP<DRT::Discretization> disB);

      //! assign material of discretization B
      virtual void AssignMaterialAToB(const VOLMORTAR::VolMortarCoupl* volmortar,
                                  DRT::Element* Bele,
                                  const std::vector<int>& Aids,
                                  Teuchos::RCP<DRT::Discretization> disA,
                                  Teuchos::RCP<DRT::Discretization> disB);
    };
  }
}


#endif /* VOLMORTAR_UTILIS_H_ */
