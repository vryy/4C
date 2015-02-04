/*!----------------------------------------------------------------------
\file xfluid_utils.H
\brief Basic functions used for xfluid applications

<pre>
Maintainer:  Magnus Winter
             winter@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

#ifndef XFLUID_UTILS_H_
#define XFLUID_UTILS_H_

#include "../drt_cut/cut_point.H"


namespace XFLUID
{
  namespace UTILS
  {
  //! @name Get material properties for the Volume Cell
  /*!

  \brief Element material for the volume cell, depending on element and position.
         If an element which is not a material list is given, the provided material is chosen.
         If however a material list is given the material chosen for the volume cell is depending on the point position.

   */
  void GetVolumeCellMaterial(DRT::Element* actele,                        // element for volume cell INPUT
                             Teuchos::RCP<MAT::Material> & mat,           // material of volume cell OUTPUT
                             const GEO::CUT::Point::PointPosition position = GEO::CUT::Point::outside //position of volume cell INPUT
                             );


  //! @name Check whether materials are identical
  /*!

  \brief A Safety check is done for XFEM-type problems. Is utilized in the edgebased framework.

   */
  void SafetyCheckMaterials(Teuchos::RCP<MAT::Material> &          pmat,
                            Teuchos::RCP<MAT::Material> &          nmat
                            );

  }

}


#endif /* XFEM_UTILS_H_ */
