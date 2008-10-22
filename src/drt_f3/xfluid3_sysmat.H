/*----------------------------------------------------------------------*/
/*!
\file xfluid3_sysmat.H

\brief element formulations for 3d XFEM fluid element

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef XFLUID3_SYSMAT_H
#define XFLUID3_SYSMAT_H

#include "xfluid3.H"
#include "../drt_xfem/xfem_enums.H"

namespace XFLUID
{
  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmatTP1(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
          const Teuchos::RCP<const Epetra_Vector> ivelcol,
          const Teuchos::RCP<Epetra_Vector> iforcecol,     ///< reaction force due to given interface velocity
          Epetra_SerialDenseMatrix&         estif,
          Epetra_SerialDenseVector&         eforce,
          const struct _MATERIAL*           material,
          const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
          const double                      dt,            ///< delta t (time step size)
          const double                      theta,         ///< factor for one step theta scheme
          const bool                        newton ,
          const bool                        pstab  ,
          const bool                        supg   ,
          const bool                        cstab  ,
          const bool                        instationary,
          const bool                        ifaceForceContribution
          );

  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmat3(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const std::vector<double>&        locval,
          const std::vector<double>&        locval_hist,
          const Teuchos::RCP<const Epetra_Vector> ivelcol,
          const Teuchos::RCP<Epetra_Vector> iforcecol,     ///< reaction force due to given interface velocity
          Epetra_SerialDenseMatrix&         estif,
          Epetra_SerialDenseVector&         eforce,
          const struct _MATERIAL*           material,
          const double                      time,          ///< current time (pseudotime for stationary formulation)
          const double                      timefac,       ///< One-step-Theta: theta*dt, BDF2: 2/3 * dt
          const bool                        newton ,
          const bool                        pstab  ,
          const bool                        supg   ,
          const bool                        cstab  ,
          const bool                        instationary
          );
  
  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmat4(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
          const Teuchos::RCP<const Epetra_Vector> ivelcol,
          const Teuchos::RCP<Epetra_Vector> iforcecol,     ///< reaction force due to given interface velocity
          Epetra_SerialDenseMatrix&         estif,
          Epetra_SerialDenseVector&         eforce,
          const struct _MATERIAL*           material,
          const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
          const double                      dt,            ///< delta t (time step size)
          const double                      theta,         ///< factor for one step theta scheme
          const bool                        newton ,
          const bool                        pstab  ,
          const bool                        supg   ,
          const bool                        cstab  ,
          const bool                        instationary,
          const bool                        ifaceForceContribution
          );
}

#endif
#endif
#endif
