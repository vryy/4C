/*----------------------------------------------------------------------*/
/*!
\file xfem_utils.H
\brief Basic tools used in XFEM routines

\level 3

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>

\warning this file should be cleaned up
*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_UTILS_H_
#define XFEM_UTILS_H_

#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_point.H"

namespace XFEM
{
  namespace UTILS
  {
    //! extract the nodal vectors and store them in node-vector-map
    //! \author schott \date 01/13
    void ExtractNodeVectors(Teuchos::RCP<DRT::Discretization> dis,
        std::map<int, LINALG::Matrix<3, 1>>& nodevecmap, Teuchos::RCP<Epetra_Vector> idispnp);

    //! @name Get material properties for the Volume Cell

    /*!

    \brief Element material for the volume cell, depending on element and position.
           If an element which is not a material list is given, the provided material is chosen.
           If however a material list is given the material chosen for the volume cell is depending
    on the point position.

     */
    void GetVolumeCellMaterial(DRT::Element* actele,  // element for volume cell INPUT
        Teuchos::RCP<MAT::Material>& mat,             // material of volume cell OUTPUT
        GEO::CUT::Point::PointPosition position =
            GEO::CUT::Point::outside  // position of volume cell INPUT to determine position
    );


    //! @name Check whether materials are identical
    /*!

    \brief A Safety check is done for XFEM-type problems. Is utilized in the edgebased framework.

     */
    void SafetyCheckMaterials(Teuchos::RCP<MAT::Material>& pmat, Teuchos::RCP<MAT::Material>& nmat);

    //! @name Extract quantities on a element
    /*!
    \brief Needs a column-vector to extract correctly in parallel
     */
    void ExtractQuantityAtElement(Epetra_SerialDenseMatrix& element_vector,
        const DRT::Element* element,
        const Teuchos::RCP<const Epetra_MultiVector>& global_col_vector,
        Teuchos::RCP<DRT::Discretization>& dis, const int nds_vector, const int nsd);

    //! @name Extract quantities on a node
    /*!
    \brief Needs a column-vector to extract correctly in parallel
     */
    void ExtractQuantityAtNode(Epetra_SerialDenseMatrix& element_vector, const DRT::Node* node,
        const Teuchos::RCP<const Epetra_MultiVector>& global_col_vector,
        Teuchos::RCP<DRT::Discretization>& dis, const int nds_vector, const unsigned int nsd);

  }  // namespace UTILS
}  // namespace XFEM


#endif /* XFEM_UTILS_H_ */
