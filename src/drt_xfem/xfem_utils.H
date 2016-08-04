/*!----------------------------------------------------------------------
\file xfem_utils.H
\brief Basic tools used in XFEM routines

\level 2

<pre>
\maintainer Magnus Winter
             winter@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

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
    void ExtractNodeVectors(
      Teuchos::RCP<DRT::Discretization>    dis,
      std::map<int, LINALG::Matrix<3,1> >& nodevecmap,
      Teuchos::RCP<Epetra_Vector>          idispnp
    );

    //! @name Get material properties for the Volume Cell

    /*!

    \brief Element material for the volume cell, depending on element and position.
           If an element which is not a material list is given, the provided material is chosen.
           If however a material list is given the material chosen for the volume cell is depending on the point position.

     */
    void GetVolumeCellMaterial( DRT::Element* actele,                        // element for volume cell INPUT
                                Teuchos::RCP<MAT::Material> & mat,           // material of volume cell OUTPUT
                                GEO::CUT::Point::PointPosition position = GEO::CUT::Point::outside // position of volume cell INPUT to determine position
                               );


    //! @name Check whether materials are identical
    /*!

    \brief A Safety check is done for XFEM-type problems. Is utilized in the edgebased framework.

     */
    void SafetyCheckMaterials(Teuchos::RCP<MAT::Material> &          pmat,
                              Teuchos::RCP<MAT::Material> &          nmat
                              );

    //! @name GetNavierSlipStabilizationParameters
    /*!

    \brief Get NavierSlip Stabilization Parameters for tangential direction

     */
    void GetNavierSlipStabilizationParameters(
      const double &                                NIT_full_stab_fac,            ///< full Nitsche stab fac
      const double &                                NIT_visc_stab_fac,            ///< viscous Nitsche stab fac
      double &                                      dynvisc,                      ///< average dynamic viscosity
      double &                                      sliplength,                   ///< sliplength
      double &                                      stabnit,                      ///< stabilization factor NIT_Penalty
      double &                                      stabadj                       ///< stabilization factor Adjoint
    );
  }
}


#endif /* XFEM_UTILS_H_ */
