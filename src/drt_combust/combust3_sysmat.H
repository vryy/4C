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
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_xfem/xfem_enums.H"

namespace COMBUST
{
  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmat(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::Combust3*     ele,
          const Teuchos::RCP<COMBUST::InterfaceHandleCombust>&  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::Combust3::MyState&  mystate,   ///< element state variables
          Epetra_SerialDenseMatrix&         estif,
          Epetra_SerialDenseVector&         eforce,
          Teuchos::RCP<const MAT::Material> material,
          const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
          const double                      dt,            ///< delta t (time step size)
          const double                      theta,         ///< factor for one step theta scheme
          const bool                        newton,
          const bool                        pstab,
          const bool                        supg,
          const bool                        cstab,
          const INPAR::FLUID::TauType       tautype,       ///< stabilization parameter definition
          const bool                        instationary,
          const INPAR::COMBUST::CombustionType combusttype,
          const double                       flamespeed,
          const double                       nitschevel,
          const double                       nitschepres
          );
}

#endif
#endif
#endif
