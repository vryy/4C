/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat.H

\brief element formulations for 3d Combust fluid element

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef COMBUST3_SYSMAT_H
#define COMBUST3_SYSMAT_H

#include "combust3.H"
#include "../drt_xfem/interface.H"

namespace COMBUST
{
  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template intantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmat3(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::Combust3*    ele,
          const Teuchos::RCP<XFEM::InterfaceHandle>  ih,
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
   * at one point, one has to call specific template intantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmat4(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::Combust3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandle>  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::Combust3::MyState  mystate,   ///< element state variables
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
          const bool                        instationary
          );
}

#endif
#endif
#endif
