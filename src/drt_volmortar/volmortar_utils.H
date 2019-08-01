/*----------------------------------------------------------------------*/
/*!

 \brief utils functions related to volmortar

\level 2

\maintainer Matthias Mayr
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
}  // namespace DRT

namespace VOLMORTAR
{
  class VolMortarCoupl;

  namespace UTILS
  {
    /// Helper class for assigning materials for volumetric coupling of non conforming meshes
    /*!
     When coupling two overlapping discretizations, most often one discretization needs access
     to the corresponding element/material on the other side. For conforming meshes this is straight
     forward as there is one unique element on the other side and therefore one unique material,
     which can be accessed. However, for non conforming meshes there are potentially several
     elements overlapping. Therefore, some rule for assigning materials is needed. This class is
     meant to do that. It gets the element to which it shall assign a material and a vector of IDs
     of the overlapping elements of the other discretization.

     The default strategy will just assign the material of the first element in the vector to the
     other element. This is fine for constant material properties, for instance. If there is furhter
     work to be done it is meant to derive from this class (see e.g. TSIMaterialStrategy).

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
      virtual void AssignMaterial2To1(const VOLMORTAR::VolMortarCoupl* volmortar,
          DRT::Element* ele1, const std::vector<int>& ids_2, Teuchos::RCP<DRT::Discretization> dis1,
          Teuchos::RCP<DRT::Discretization> dis2);

      //! assign material of discretization B
      virtual void AssignMaterial1To2(const VOLMORTAR::VolMortarCoupl* volmortar,
          DRT::Element* ele2, const std::vector<int>& ids_1, Teuchos::RCP<DRT::Discretization> dis1,
          Teuchos::RCP<DRT::Discretization> dis2);
    };
  }  // namespace UTILS
}  // namespace VOLMORTAR


#endif /* VOLMORTAR_UTILIS_H_ */
