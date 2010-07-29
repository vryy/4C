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
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

namespace XFLUID
{
  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmatTauPressure(
          ParameterList&                    params,
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
          const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
          Epetra_SerialDenseMatrix&         estif,
          Epetra_SerialDenseVector&         eforce,
          Teuchos::RCP<const MAT::Material> material,
          const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
          const double                      dt,            ///< delta t (time step size)
          const double                      theta,         ///< factor for one step theta scheme
          const bool                        newton ,
          const bool                        pstab  ,
          const bool                        supg   ,
          const bool                        cstab  ,
          const bool                        ifaceForceContribution, ///< is xfem element a row element?
          const bool                        monolithic_FSI,
          double&                           L2,
          DRT::ELEMENTS::XFluid3::FluidFluidCouplingMatrices fluidfluidmatrices
          );

  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmatSigma(
          ParameterList&                    params,
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
          const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
          Epetra_SerialDenseMatrix&         estif,         ///< element stiffness matrix
          Epetra_SerialDenseVector&         eforce,        ///< element rhs vector
          Teuchos::RCP<const MAT::Material> material,
          const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
          const double                      dt,            ///< delta t (time step size)
          const double                      theta,         ///< factor for one step theta scheme
          const bool                        newton ,
          const bool                        pstab  ,
          const bool                        supg   ,
          const bool                        cstab  ,
          const bool                        ifaceForceContribution, ///< is xfem element a row element?
          const bool                        monolithic_FSI,
          double&                           L2,
          DRT::ELEMENTS::XFluid3::FluidFluidCouplingMatrices fluidfluidmatrices
          );


  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  inline void callSysmat(
          INPAR::XFEM::BoundaryIntegralType btype,
          ParameterList&                    params,
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
          const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
          Epetra_SerialDenseMatrix&         estif,
          Epetra_SerialDenseVector&         eforce,
          Teuchos::RCP<const MAT::Material> material,
          const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
          const double                      dt,            ///< delta t (time step size)
          const double                      theta,         ///< factor for one step theta scheme
          const bool                        newton ,
          const bool                        pstab  ,
          const bool                        supg   ,
          const bool                        cstab  ,
          const bool                        ifaceForceContribution, ///< is xfem element a row element?
          const bool                        monolithic_FSI,
          double&                           L2,
          DRT::ELEMENTS::XFluid3::FluidFluidCouplingMatrices fluidfluidmatrices
          )
  {
    switch (btype)
    {
    case INPAR::XFEM::BoundaryTypeSigma:
      XFLUID::callSysmatSigma(params, assembly_type,
          ele, ih, eleDofManager, mystate, iforcecol, estif, eforce,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, ifaceForceContribution, monolithic_FSI, L2, fluidfluidmatrices);
      break;
    case INPAR::XFEM::BoundaryTypeTauPressure:
      XFLUID::callSysmatTauPressure(params, assembly_type,
          ele, ih, eleDofManager, mystate, iforcecol, estif, eforce,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, ifaceForceContribution, monolithic_FSI, L2, fluidfluidmatrices);
      break;
    default:
      dserror("unknown boundary type");
    }
  }

  /*!
   * \brief entry point for Sysmat call
   * at one point, one has to call specific template instantiations of Sysmat using the current Shape() of the element.
   * This is the point.
   */
  void callSysmatProjection(
          const XFEM::AssemblyType          assembly_type,
          const DRT::ELEMENTS::XFluid3*     ele,
          const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
          Epetra_SerialDenseMatrix&         estif_veln,         ///< element matrix to calculate
          Epetra_SerialDenseMatrix&         estif_accn,         ///< element matrix to calculate
          Epetra_SerialDenseVector&         eforce_veln,        ///< element rhs to calculate
          Epetra_SerialDenseVector&         eforce_accn,        ///< element rhs to calculate
          const bool                        pstab  ,
          const bool                        ifaceForceContribution ///< is xfem element a row element?
          );

}

#endif
#endif
#endif
