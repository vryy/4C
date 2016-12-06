/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_premixed_nitsche.H

\brief system matrix formulation with Nitsche's method
       for premixed combustion problem / 3D combustion element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_PREMIXED_NITSCHE_H
#define COMBUST3_SYSMAT_PREMIXED_NITSCHE_H


#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat_standard.H"
#include "combust3_sysmat_surfacetension.H"
#include "combust3_sysmat_turbulence.H"
#include "combust_defines.H"
#include "combust3_stabilization.H"
#ifdef DEBUG
#include "../drt_mat/material.H"
#include "../drt_inpar/inpar_material.H"
#endif


namespace COMBUST
{


/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class V2>
void SysmatDomainNitsche(
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const M1&                      evelnp,       ///< nodal velocity n+1 values for this element
    const M1&                      eveln,        ///< nodal velocity n   values for this element
    const M1&                      evelnm,       ///< nodal velocity n-1 values for this element
    const M1&                      eaccn,        ///< nodal acceleration n values for this element
    const M1&                      eaccam,       ///< nodal acceleration n values for this element
    const M1&                      efsvel,       ///< nodal fine-scale velocity n+1 values for this element
    const V1&                      eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                      ephi,         ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material> material,  ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,     ///< time discretization scheme
    const double                   time,         ///< current time
    const double                   dt,           ///< delta t (time step size)
    const double                   theta,        ///< parameter for one step theta scheme
    const double                   ga_alphaF,
    const double                   ga_alphaM,
    const double                   ga_gamma,
    const bool                     newton,       ///< full Newton or fixed-point-like
    const bool                     pstab,        ///< turn on/off PSPG stabilization
    const bool                     supg,         ///< turn on/off SUPG stabilization
    const bool                     graddiv,        ///< turn on/off BV   stabilization
    const INPAR::FLUID::TauType    tautype,      //< stabilization parameter definition
    const bool                     instationary, ///< switch between stationary and instationary formulation
    const bool                     genalpha,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    double&                        ele_meas_plus,      // set measure of element in plus domain
    double&                        ele_meas_minus,     // set measure of element in minus domain
    const INPAR::FLUID::TurbModelAction    turbmodel,
    const INPAR::FLUID::FineSubgridVisc    fssgv,
    const double                           Cs,
    const double                           Csgs,
    const double                           alpha,
    const bool                             CalcN,
    const double                           N,
    const INPAR::FLUID::RefVelocity        refvel,
    const INPAR::FLUID::RefLength          reflength,
    const double                           c_nu,
    const bool                             near_wall_limit,
    const bool                             B_gp,
    const bool                             mfs_is_conservative
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - domain");

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration factor
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of this element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // read body force (dead load) from input file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> edead(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, "DeadOrLive", edead);

  // read the pressure gradient (body force w/out scaling w/ density)
  // from inut file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> epgrad(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, "neum_pgrad", epgrad);

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // stabilization parameter
  const double hk = COMBUST::UTILS::HK<DISTYPE>(xyze);
  const double mk = COMBUST::UTILS::MK<DISTYPE>();

  // element volume
  const double elevol = COMBUST::UTILS::VolumeViaNumIntegration<DISTYPE>(xyze);

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  //-----------------------------------
  // loop over domain integration cells
  //-----------------------------------
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    //--------------------------------------------------
    // get material parameters for this integration cell
    //--------------------------------------------------
    // density \rho^{n+1}
    double densaf  = -1.0;
    // dynamic viscosity \mu^{n+1}
    double dynvisc = -1.0;
    // set density and dynamic viscosity
    {
      double densplus = -1.0;
      double densminus = -1.0;
      double dynviscplus = -1.0;
      double dynviscminus = -1.0;
      COMBUST::GetMaterialParams(material, true, densplus, dynviscplus);
      COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);

      // blend smoothly to minus (gas bubble) domain
      COMBUST::BlendMaterial(ele, cell->getDomainPlus(), time, densaf, dynvisc, densplus, dynviscplus, densminus, dynviscminus);
    }
    // density \rho^{n}
    // until the XFEM time integration is done properly, densities from previous and
    // current time steps are identical
    const double densn = densaf;

    //-----------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //-----------------------------------------------------------
#ifndef COMBUST_SXFEM
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);
#endif
    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Splited(), cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of the current integration point in cell coordinates \eta^domain
      const LINALG::Matrix<nsd,1> pos_eta_domain(intpoints.qxg[iquad]);

      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
      const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
#ifdef DEBUG
      if (detcell < 0.0)
      {
        std::cout << "detcell :  " << detcell << std::endl;
        dserror("negative detcell!");
      }
#endif
      //----------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //----------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //-------------------------------------
      // procedures involving Jacobian matrix
      //-------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      //-------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //-------------------------------------------------------------------
      static LINALG::Matrix<3,numnode> derxy;
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

#ifdef COMBUST_SXFEM
      const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,ephi,*cell,funct,derxy,derxy2);
      //const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);
#endif

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSize,1> enr_funct;
        enr_funct.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy;
        enr_derxy.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2;
        enr_derxy2.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,
            funct,
            derxy,
            derxy2,
            enr_funct,
            enr_derxy,
            enr_derxy2);

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shp.d0(iparam) = enr_funct(iparam);
          shp.dx(iparam) = enr_derxy(0,iparam);
          shp.dy(iparam) = enr_derxy(1,iparam);
          shp.dz(iparam) = enr_derxy(2,iparam);
          shp.dxdx(iparam) = enr_derxy2(0,iparam);
          shp.dxdy(iparam) = enr_derxy2(3,iparam);
          shp.dxdz(iparam) = enr_derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = enr_derxy2(1,iparam);
          shp.dydz(iparam) = enr_derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = enr_derxy2(2,iparam);
        }
      }
      else // not xfem_assembly i.e. standard assembly
      {
        // fill approximation functions for standard FEM
        // remark: numparamvelx == numnode, for standard FEM
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shp.d0(iparam) = funct(iparam);
          shp.dx(iparam) = derxy(0,iparam);
          shp.dy(iparam) = derxy(1,iparam);
          shp.dz(iparam) = derxy(2,iparam);
          shp.dxdx(iparam) = derxy2(0,iparam);
          shp.dxdy(iparam) = derxy2(3,iparam);
          shp.dxdz(iparam) = derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = derxy2(1,iparam);
          shp.dydz(iparam) = derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = derxy2(2,iparam);
        }
      }
      // set measure of element in plus and minus domain
      if(cell->getDomainPlus() == true) {ele_meas_plus  += fac;}
      else  {ele_meas_minus += fac;}
      //------------------------------------------------------
      // get velocities and accelerations at integration point
      //------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);
      // get acceleration at n+alpha_m for generalized alpha time integration scheme
      LINALG::Matrix<3,1> gpaccam(true);
      if (genalpha)
        gpaccam  = XFEM::interpolateVectorFieldToIntPoint(eaccam , shp.d0, numparamvelx);

      //----------------------------------------------
      // get velocity derivatives at integration point
      //----------------------------------------------
      static LINALG::Matrix<3,nsd> vderxy;
      vderxy.Clear();
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy(isd,0) += evelnp(isd,iparam) * shp.dx(iparam);
          vderxy(isd,1) += evelnp(isd,iparam) * shp.dy(iparam);
          vderxy(isd,2) += evelnp(isd,iparam) * shp.dz(iparam);
        }
      }

      //-----------------------------------------------------
      // get second velocity derivatives at integration point
      //-----------------------------------------------------
      static LINALG::Matrix<3,6> vderxy2;
      vderxy2.Clear();
      //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
      if (higher_order_ele)
      {
        for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          for (size_t isd = 0; isd < nsd; ++isd)
          {
            vderxy2(isd,0) += evelnp(isd,iparam) * shp.dxdx(iparam);
            vderxy2(isd,1) += evelnp(isd,iparam) * shp.dydy(iparam);
            vderxy2(isd,2) += evelnp(isd,iparam) * shp.dzdz(iparam);
            vderxy2(isd,3) += evelnp(isd,iparam) * shp.dxdy(iparam);
            vderxy2(isd,4) += evelnp(isd,iparam) * shp.dxdz(iparam);
            vderxy2(isd,5) += evelnp(isd,iparam) * shp.dydz(iparam);
          }
        }
      }

      //----------------------------------
      // get pressure at integration point
      //----------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp.d0(iparam);

      //--------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp.dx(iparam);
        gradp(1) += eprenp(iparam) * shp.dy(iparam);
        gradp(2) += eprenp(iparam) * shp.dz(iparam);
      }

      //------------------------------------------
      // compute stabilization parameters (3 taus)
      //------------------------------------------
      const double vel_norm = gpvelnp.Norm2();
      //const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_Mu = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;
      // remark: stabilization parameters 'tau_stab_Mu' and 'tau_stab_Mp' are both computed based on
      //         'hk' instead of using the streamlength for 'tau_stab_Mu'
      COMBUST::UTILS::computeStabilizationParams(gpvelnp, xji,
          instationary, dynvisc, densaf, vel_norm, hk, hk, mk, timefac, dt,
          tautype, tau_stab_Mu, tau_stab_Mp, tau_stab_C);

      //------------------------------------------
      // prepare turbulence models
      //------------------------------------------

      static LINALG::Matrix<3,nsd> fsvderxy(true);
      fsvderxy.Clear();
      if(fssgv != INPAR::FLUID::no_fssgv or turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
      {
        // get fine scale velocity (np,i) derivatives at integration point
         for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
         {
           for (size_t isd = 0; isd < nsd; ++isd)
           {
             fsvderxy(isd,0) += efsvel(isd,iparam) * shp.dx(iparam);
             fsvderxy(isd,1) += efsvel(isd,iparam) * shp.dy(iparam);
             fsvderxy(isd,2) += efsvel(isd,iparam) * shp.dz(iparam);
           }
         }
      }

      //----------------
      // calculate finscale subgrid viscosity
      // for turbulent two-phase flows with AVM3 approach
      //----------------
      // kinematic fine-scale viscosity
      // (multiplication by density below in element term)
      double fssgvisc=0.0;
      if (fssgv != INPAR::FLUID::no_fssgv)
        fssgvisc = COMBUST::CalcFineScaleSubgrVisc<DISTYPE>(vderxy,fsvderxy,elevol, Cs,fssgv);

      //----------------
      // multifractal subgrid-scales
      //----------------

      LINALG::Matrix<nsd,1> gpmfsvelnp(true);
      LINALG::Matrix<3,nsd> mfsvderxy(true);
      double mfsvdiv = 0.0;
      double vdiv = 0.0;

      if (turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
      {
        // get fine-scale velocity at integration point
        LINALG::Matrix<nsd,1> gpfsvelnp = XFEM::interpolateVectorFieldToIntPoint(efsvel, shp.d0, numparamvelx);

        // calculate parameters of multifractal subgrid-scales
        double B_mfs = 0.0;
        PrepareMultifractalSubgrScales<DISTYPE>(B_mfs, gpvelnp, gpfsvelnp, vderxy, elevol,
          Csgs, alpha, CalcN, N, reflength, refvel, c_nu, near_wall_limit, densaf, dynvisc);

        // calculate fine-scale velocity, its derivative and divergence for multifractal subgrid-scale modeling
        for (std::size_t idim=0; idim<nsd; idim++)
            gpmfsvelnp(idim,0) = gpfsvelnp(idim,0) * B_mfs;

        for (std::size_t idim=0; idim<nsd; idim++)
        {
          for (std::size_t jdim=0; jdim<nsd; jdim++)
              mfsvderxy(idim,jdim) = fsvderxy(idim,jdim) * B_mfs;
        }

        mfsvdiv = mfsvderxy(0,0) + mfsvderxy(1,1) + mfsvderxy(2,2);
        vdiv = vderxy(0,0) + vderxy(1,1) + vderxy(2,2);
      }

      //------------------------------------------------------------------------
      // build RHS vector at integration point (scaled with timefac and density)
      //------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> rhsint;
      rhsint.Clear();

      // time RHS vector from previous time step
      // remark: histvec = \rho^{n}*(u^n + delta t (1-theta) a^n)
      LINALG::Matrix<nsd,1> histvec = COMBUST::TIMEINT::GetOldPartOfRighthandside(
          gpveln, gpvelnm, gpaccn, timealgo, dt, theta);
      histvec.Scale(densn);

      // body force vector
      // remark: unit is acceleration (not density weighted)
      static LINALG::Matrix<nsd,1> bodyforce;
      bodyforce.Clear();
      static LINALG::Matrix<nsd,1> noscalebodyforce;
      noscalebodyforce.Clear();
      for (std::size_t isd = 0; isd < nsd; isd++)
      {
        for (std::size_t inode = 0; inode < numnode; inode++)
        {
          bodyforce(isd) += edead(isd,inode) * funct(inode);
          noscalebodyforce(isd) += epgrad(isd,inode) * funct(inode);
        }
      }
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densaf);

      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = bodyforce(isd) + noscalebodyforce(isd);
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = histvec(isd) + timefac*(bodyforce(isd) + noscalebodyforce(isd));
      }

      //-----------------------------------------------------------------------
      // compute numerical representation of some single operators (not scaled)
      //-----------------------------------------------------------------------

      // convective term:  u_old * grad u_old:
      // remark: conv_old = vderxy(i,j)*gpvelnp(j);
      LINALG::Matrix<nsd,1> conv_old(true);
      conv_old.Multiply(vderxy,gpvelnp);

      // viscous term:  div epsilon(u_old)
      LINALG::Matrix<nsd,1> visc_old;
      visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
      visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
      visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

      // convective part of linearized convective term: u_old * grad (funct)
      // remark: u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      static LINALG::Matrix<shpVecSize,1> enr_conv_c;
      enr_conv_c.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_conv_c(iparam) += shp.dx(iparam)*gpvelnp(0);
        enr_conv_c(iparam) += shp.dy(iparam)*gpvelnp(1);
        enr_conv_c(iparam) += shp.dz(iparam)*gpvelnp(2);
      }

      // viscous term  - grad * epsilon(u)
      /*    /                                              \
                 |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
               1 |                                                |
               - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
               2 |                                                |
                 |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
                  \                                              /
       */
      static COMBUST::EnrViscs2<shpVecSize> enr_viscs2;
      //            enr_viscs2.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_viscs2.xx(iparam) = 0.5 * (2.0 * shp.dxdx(iparam) + shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.xy(iparam) = 0.5 *  shp.dydx(iparam);
        enr_viscs2.xz(iparam) = 0.5 *  shp.dzdx(iparam);
        enr_viscs2.yx(iparam) = 0.5 *  shp.dxdy(iparam);
        enr_viscs2.yy(iparam) = 0.5 * (shp.dxdx(iparam) + 2.0 * shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.yz(iparam) = 0.5 *  shp.dzdy(iparam);
        enr_viscs2.zx(iparam) = 0.5 *  shp.dxdz(iparam);
        enr_viscs2.zy(iparam) = 0.5 *  shp.dydz(iparam);
        enr_viscs2.zz(iparam) = 0.5 * (shp.dxdx(iparam) + shp.dydy(iparam) + 2.0 * shp.dzdz(iparam));
      }

      //----------------------------------------------------
      // evaluate residual (scaled with timefac and density)
      // remark: used for stabilization right hand sides
      //----------------------------------------------------
      //            LINALG::Matrix<nsd,1> res_old(true);
      //            if (instationary)
      //            {
      //              // timefac * (\rho*du/dt + (\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              for (size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = densaf*gpvelnp(isd) - rhsint(isd) + timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));
      //            }
      //            else // stationary
      //            {
      //              // timefac * ((\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              // remark: Although calculation is stationary, scaling with timefac is kept to be
      //              //         consistent with formulation of matrix entries, which are also scaled.
      //              for(size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd)-bodyforce(isd));
      //            }

      // evaluate residual once for all stabilisation right hand sides
      LINALG::Matrix<nsd,1> res_old(true);
      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = densaf*conv_old(isd) + gradp(isd) - 2.0*dynvisc*visc_old(isd) - bodyforce(isd);
        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpaccam(isd);
        }
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpvelnp(isd);
        }
      }
      // integration factors and coefficients of single terms
      double rhsfac     = fac;
      if (genalpha)
        rhsfac *= timefac / ga_alphaF; // (fac * gamma * dt ) / alpha_M

      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize>(
          assembler, shp, shp, fac, timefac, rhsfac, densaf, dynvisc,
          gpvelnp, pres, gpaccam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_conv_c, enr_viscs2,
          instationary, genalpha, newton, pstab, supg, graddiv,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);

      //-------------------------------------------
      // build stiffness matrix turbulence entries
      //-------------------------------------------
      if (fssgv != INPAR::FLUID::no_fssgv or turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
        COMBUST::BuildTurbulenceDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize>(
               fssgv, turbmodel,
               assembler, shp, fsvderxy, gpvelnp, gpmfsvelnp, vderxy, mfsvderxy, vdiv, mfsvdiv,
               genalpha, fac, timefac, rhsfac, densaf, fssgvisc, mfs_is_conservative);
    } // end loop over gauss points
  } // end loop over integration cells

  return;
}


/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class V2>
void SysmatDomainNitscheGalerkin(
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const M1&                      evelnp,       ///< nodal velocity n+1 values for this element
    const M1&                      eveln,        ///< nodal velocity n   values for this element
    const M1&                      evelnm,       ///< nodal velocity n-1 values for this element
    const M1&                      eaccn,        ///< nodal acceleration n values for this element
    const M1&                      eaccam,       ///< nodal acceleration n values for this element
    const V1&                      eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                      ephi,         ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material> material,  ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,     ///< time discretization scheme
    const double                   time,         ///< current time
    const double                   dt,           ///< delta t (time step size)
    const double                   theta,        ///< parameter for one step theta scheme
    const double                   ga_alphaF,
    const double                   ga_alphaM,
    const double                   ga_gamma,
    const bool                     newton,       ///< full Newton or fixed-point-like
    const bool                     pstab,        ///< turn on/off PSPG stabilization
    const bool                     supg,         ///< turn on/off SUPG stabilization
    const bool                     graddiv,        ///< turn on/off BV   stabilization
    const INPAR::FLUID::TauType    tautype,      //< stabilization parameter definition
    const bool                     instationary, ///< switch between stationary and instationary formulation
    const bool                     genalpha,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    double&                        ele_meas_plus,      // set measure of element in plus domain
    double&                        ele_meas_minus      // set measure of element in minus domain
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - domain");

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration factor
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of this element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // read body force (dead load) from input file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> edead(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, "DeadOrLive", edead);

  // read the pressure gradient (body force w/out scaling w/ density)
  // from inut file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> epgrad(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, "neum_pgrad", epgrad);

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // stabilization parameter
  const double hk = COMBUST::UTILS::HK<DISTYPE>(evelnp,xyze);
  const double mk = COMBUST::UTILS::MK<DISTYPE>();

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  //-----------------------------------
  // loop over domain integration cells
  //-----------------------------------
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    //--------------------------------------------------
    // get material parameters for this integration cell
    //--------------------------------------------------
    // density \rho^{n+1}
    double densaf  = -1.0;
    // dynamic viscosity \mu^{n+1}
    double dynvisc = -1.0;
    // set density and dynamic viscosity
    {
      double densplus = -1.0;
      double densminus = -1.0;
      double dynviscplus = -1.0;
      double dynviscminus = -1.0;
      COMBUST::GetMaterialParams(material, true, densplus, dynviscplus);
      COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);

      // blend smoothly to minus (gas bubble) domain
      COMBUST::BlendMaterial(ele, cell->getDomainPlus(), time, densaf, dynvisc, densplus, dynviscplus, densminus, dynviscminus);
    }
    // density \rho^{n}
    // until the XFEM time integration is done properly, densities from previous and
    // current time steps are identical
    const double densn = densaf;

    //-----------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //-----------------------------------------------------------
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);

    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Splited(), cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of the current integration point in cell coordinates \eta^domain
      const LINALG::Matrix<nsd,1> pos_eta_domain(intpoints.qxg[iquad]);

      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
      const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
#ifdef DEBUG
      if (detcell < 0.0)
      {
        std::cout << "detcell :  " << detcell << std::endl;
        dserror("negative detcell!");
      }
#endif
      //----------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //----------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //-------------------------------------
      // procedures involving Jacobian matrix
      //-------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      //-------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //-------------------------------------------------------------------
      static LINALG::Matrix<3,numnode> derxy;
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSize,1> enr_funct;
        enr_funct.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy;
        enr_derxy.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2;
        enr_derxy2.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,
            funct,
            derxy,
            derxy2,
            enr_funct,
            enr_derxy,
            enr_derxy2);

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shp.d0(iparam) = enr_funct(iparam);
          shp.dx(iparam) = enr_derxy(0,iparam);
          shp.dy(iparam) = enr_derxy(1,iparam);
          shp.dz(iparam) = enr_derxy(2,iparam);
          shp.dxdx(iparam) = enr_derxy2(0,iparam);
          shp.dxdy(iparam) = enr_derxy2(3,iparam);
          shp.dxdz(iparam) = enr_derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = enr_derxy2(1,iparam);
          shp.dydz(iparam) = enr_derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = enr_derxy2(2,iparam);
        }
      }
      else // not xfem_assembly i.e. standard assembly
      {
        // fill approximation functions for standard FEM
        // remark: numparamvelx == numnode, for standard FEM
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shp.d0(iparam) = funct(iparam);
          shp.dx(iparam) = derxy(0,iparam);
          shp.dy(iparam) = derxy(1,iparam);
          shp.dz(iparam) = derxy(2,iparam);
          shp.dxdx(iparam) = derxy2(0,iparam);
          shp.dxdy(iparam) = derxy2(3,iparam);
          shp.dxdz(iparam) = derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = derxy2(1,iparam);
          shp.dydz(iparam) = derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = derxy2(2,iparam);
        }
      }
      // set measure of element in plus and minus domain
      if(cell->getDomainPlus() == true) {ele_meas_plus  += fac;}
      else  {ele_meas_minus += fac;}
      //------------------------------------------------------
      // get velocities and accelerations at integration point
      //------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);
      // get acceleration at n+alpha_m for generalized alpha time integration scheme
      LINALG::Matrix<3,1> gpaccam(true);
      if (genalpha)
        gpaccam  = XFEM::interpolateVectorFieldToIntPoint(eaccam , shp.d0, numparamvelx);

      //----------------------------------------------
      // get velocity derivatives at integration point
      //----------------------------------------------
      static LINALG::Matrix<3,nsd> vderxy;
      vderxy.Clear();
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy(isd,0) += evelnp(isd,iparam) * shp.dx(iparam);
          vderxy(isd,1) += evelnp(isd,iparam) * shp.dy(iparam);
          vderxy(isd,2) += evelnp(isd,iparam) * shp.dz(iparam);
        }
      }

      //-----------------------------------------------------
      // get second velocity derivatives at integration point
      //-----------------------------------------------------
      static LINALG::Matrix<3,6> vderxy2;
      vderxy2.Clear();
      //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
      if (higher_order_ele)
      {
        for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          for (size_t isd = 0; isd < nsd; ++isd)
          {
            vderxy2(isd,0) += evelnp(isd,iparam) * shp.dxdx(iparam);
            vderxy2(isd,1) += evelnp(isd,iparam) * shp.dydy(iparam);
            vderxy2(isd,2) += evelnp(isd,iparam) * shp.dzdz(iparam);
            vderxy2(isd,3) += evelnp(isd,iparam) * shp.dxdy(iparam);
            vderxy2(isd,4) += evelnp(isd,iparam) * shp.dxdz(iparam);
            vderxy2(isd,5) += evelnp(isd,iparam) * shp.dydz(iparam);
          }
        }
      }

      //----------------------------------
      // get pressure at integration point
      //----------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp.d0(iparam);

      //--------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp.dx(iparam);
        gradp(1) += eprenp(iparam) * shp.dy(iparam);
        gradp(2) += eprenp(iparam) * shp.dz(iparam);
      }

      //------------------------------------------
      // compute stabilization parameters (3 taus)
      //------------------------------------------
      const double vel_norm = gpvelnp.Norm2();
      //const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_Mu = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;
      // remark: stabilization parameters 'tau_stab_Mu' and 'tau_stab_Mp' are both computed based on
      //         'hk' instead of using the streamlength for 'tau_stab_Mu'
      COMBUST::UTILS::computeStabilizationParams(gpvelnp, xji,
          instationary, dynvisc, densaf, vel_norm, hk, hk, mk, timefac, dt,
          tautype, tau_stab_Mu, tau_stab_Mp, tau_stab_C);

      //------------------------------------------------------------------------
      // build RHS vector at integration point (scaled with timefac and density)
      //------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> rhsint;
      rhsint.Clear();

      // time RHS vector from previous time step
      // remark: histvec = \rho^{n}*(u^n + delta t (1-theta) a^n)
      LINALG::Matrix<nsd,1> histvec = COMBUST::TIMEINT::GetOldPartOfRighthandside(
          gpveln, gpvelnm, gpaccn, timealgo, dt, theta);
      histvec.Scale(densn);

      // body force vector
      // remark: unit is acceleration (not density weighted)
      static LINALG::Matrix<nsd,1> bodyforce;
      bodyforce.Clear();
      static LINALG::Matrix<nsd,1> noscalebodyforce;
      noscalebodyforce.Clear();
      for (std::size_t isd = 0; isd < nsd; isd++)
      {
        for (std::size_t inode = 0; inode < numnode; inode++)
        {
          bodyforce(isd) += edead(isd,inode) * funct(inode);
          noscalebodyforce(isd) += epgrad(isd,inode) * funct(inode);
        }
      }
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densaf);

      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = bodyforce(isd) + noscalebodyforce(isd);
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = histvec(isd) + timefac*(bodyforce(isd) + noscalebodyforce(isd));
      }

      //-----------------------------------------------------------------------
      // compute numerical representation of some single operators (not scaled)
      //-----------------------------------------------------------------------

      // convective term:  u_old * grad u_old:
      // remark: conv_old = vderxy(i,j)*gpvelnp(j);
      LINALG::Matrix<nsd,1> conv_old(true);
      conv_old.Multiply(vderxy,gpvelnp);

      // viscous term:  div epsilon(u_old)
      LINALG::Matrix<nsd,1> visc_old;
      visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
      visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
      visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

      // convective part of linearized convective term: u_old * grad (funct)
      // remark: u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      static LINALG::Matrix<shpVecSize,1> enr_conv_c;
      enr_conv_c.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_conv_c(iparam) += shp.dx(iparam)*gpvelnp(0);
        enr_conv_c(iparam) += shp.dy(iparam)*gpvelnp(1);
        enr_conv_c(iparam) += shp.dz(iparam)*gpvelnp(2);
      }

      // viscous term  - grad * epsilon(u)
      /*    /                                              \
                 |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
               1 |                                                |
               - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
               2 |                                                |
                 |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
                  \                                              /
       */
      static COMBUST::EnrViscs2<shpVecSize> enr_viscs2;
      //            enr_viscs2.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_viscs2.xx(iparam) = 0.5 * (2.0 * shp.dxdx(iparam) + shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.xy(iparam) = 0.5 *  shp.dydx(iparam);
        enr_viscs2.xz(iparam) = 0.5 *  shp.dzdx(iparam);
        enr_viscs2.yx(iparam) = 0.5 *  shp.dxdy(iparam);
        enr_viscs2.yy(iparam) = 0.5 * (shp.dxdx(iparam) + 2.0 * shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.yz(iparam) = 0.5 *  shp.dzdy(iparam);
        enr_viscs2.zx(iparam) = 0.5 *  shp.dxdz(iparam);
        enr_viscs2.zy(iparam) = 0.5 *  shp.dydz(iparam);
        enr_viscs2.zz(iparam) = 0.5 * (shp.dxdx(iparam) + shp.dydy(iparam) + 2.0 * shp.dzdz(iparam));
      }

      //----------------------------------------------------
      // evaluate residual (scaled with timefac and density)
      // remark: used for stabilization right hand sides
      //----------------------------------------------------
      //            LINALG::Matrix<nsd,1> res_old(true);
      //            if (instationary)
      //            {
      //              // timefac * (\rho*du/dt + (\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              for (size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = densaf*gpvelnp(isd) - rhsint(isd) + timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));
      //            }
      //            else // stationary
      //            {
      //              // timefac * ((\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              // remark: Although calculation is stationary, scaling with timefac is kept to be
      //              //         consistent with formulation of matrix entries, which are also scaled.
      //              for(size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd)-bodyforce(isd));
      //            }

      // evaluate residual once for all stabilisation right hand sides
      LINALG::Matrix<nsd,1> res_old(true);
      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = densaf*conv_old(isd) + gradp(isd) - 2.0*dynvisc*visc_old(isd) - bodyforce(isd);
        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpaccam(isd);
        }
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpvelnp(isd);
        }
      }
      // integration factors and coefficients of single terms
      double rhsfac = fac;
      if (genalpha)
        rhsfac *= timefac / ga_alphaF; // (fac * gamma * dt ) / alpha_M

      if (ASSTYPE== XFEM::xfem_assembly)
      {
        if(ele->Splited())
        {
//          std::cout << "funct" << funct << std::endl;
//          std::cout << "evelnp " << evelnp << std::endl;
//          std::cout << "eprenp " << eprenp << std::endl;
//          std::cout << "shp" << std::endl;
//          for(size_t iparam=0;iparam<shpVecSize;++iparam)
//          {
//            std::cout << " " << shp.d0(iparam);
//          } std::cout << std::endl;
//          std::cout << "gpvelnp" << gpvelnp << std::endl;
//          std::cout << "vderxy " << vderxy << std::endl;
//          std::cout << "vderxy2" << vderxy2 << std::endl;
//          std::cout << "pres " << pres << std::endl;
//          std::cout << "gradp " << gradp << std::endl;
//          std::cout << "visc_old " << visc_old << std::endl;
//          std::cout << "enr_conv_c" << std::endl;
          //enr_conv_c.Print();

//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.xx " << enr_viscs2.xx.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.xx " << enr_viscs2.xx.n(iparam);
//          }
//          std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.xy " << enr_viscs2.xy.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.xy " << enr_viscs2.xy.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.xz " << enr_viscs2.xz.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.xz " << enr_viscs2.xz.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.yx " << enr_viscs2.yx.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.yx " << enr_viscs2.yx.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.yy " << enr_viscs2.yy.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.yy " << enr_viscs2.yy.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.yz " << enr_viscs2.yz.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.yz " << enr_viscs2.yz.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.zx " << enr_viscs2.zx.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.zx " << enr_viscs2.zx.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.zy " << enr_viscs2.zy.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.zy " << enr_viscs2.zy.n(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.zz " << enr_viscs2.zz.s(iparam);
//          } std::cout << std::endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            std::cout << "enr_viscs2.zz " << enr_viscs2.zz.n(iparam);
//          } std::cout << std::endl;
        }
      }
      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize>(
          assembler, shp, shp, fac, timefac, rhsfac, densaf, dynvisc,
          gpvelnp, pres, gpaccam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_conv_c, enr_viscs2,
          instationary, genalpha, newton, pstab, supg, graddiv,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
    } // end loop over gauss points
  } // end loop over integration cells

  return;
}


/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class V2>
void SysmatDomainNitscheStabHexRule(
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const M1&                      evelnp,       ///< nodal velocity n+1 values for this element
    const M1&                      eveln,        ///< nodal velocity n   values for this element
    const M1&                      evelnm,       ///< nodal velocity n-1 values for this element
    const M1&                      eaccn,        ///< nodal acceleration n values for this element
    const M1&                      eaccam,       ///< nodal acceleration n values for this element
    const V1&                      eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                      ephi,         ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material> material,  ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,     ///< time discretization scheme
    const double                   time,         ///< current time
    const double                   dt,           ///< delta t (time step size)
    const double                   theta,        ///< parameter for one step theta scheme
    const double                   ga_alphaF,
    const double                   ga_alphaM,
    const double                   ga_gamma,
    const bool                     newton,       ///< full Newton or fixed-point-like
    const bool                     pstab,        ///< turn on/off PSPG stabilization
    const bool                     supg,         ///< turn on/off SUPG stabilization
    const bool                     graddiv,        ///< turn on/off BV   stabilization
    const INPAR::FLUID::TauType    tautype,      //< stabilization parameter definition
    const bool                     instationary, ///< switch between stationary and instationary formulation
    const bool                     genalpha,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    bool                           domainplus
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - domain");

  if (ASSTYPE != XFEM::xfem_assembly)
    dserror("you should not be here");

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration factor
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of this element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // read body force (dead load) from input file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> edead(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, "DeadOrLive", edead);

  // read the pressure gradient (body force w/out scaling w/ density)
  // from inut file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> epgrad(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, "neum_pgrad", epgrad);

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // stabilization parameter
  const double hk = COMBUST::UTILS::HK<DISTYPE>(evelnp,xyze);
  const double mk = COMBUST::UTILS::MK<DISTYPE>();

  // get domain integration cells for this element
  //const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));

  //-----------------------------------
  // loop over domain integration cells
  //-----------------------------------
  //for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  //{
    //--------------------------------------------------
    // get material parameters for this integration cell
    //--------------------------------------------------
    // density \rho^{n+1}
    double densaf  = -1.0;
    // dynamic viscosity \mu^{n+1}
    double dynvisc = -1.0;
    // set density and dynamic viscosity
    COMBUST::GetMaterialParams(material, domainplus, densaf, dynvisc);

//    // blend smoothly to plus (burnt) domain
//    if (cell->getDomainPlus())
//    {
//      double densminus = -1.0;
//      double dynviscminus = -1.0;
//      COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);
//
//      COMBUST::BlendMaterial(ele, time, densaf, dynvisc, densminus, dynviscminus);
//    }
    // density \rho^{n}
    // until the XFEM time integration is done properly, densities from previous and
    // current time steps are identical
    const double densn = densaf;

    //-----------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //-----------------------------------------------------------
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,ephi,domainplus);

    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, false,  DRT::Element::hex8);
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
if (intpoints.nquad != 8)
  dserror("8 Gausspoints expected");
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of the current integration point in cell coordinates \eta^domain
      const LINALG::Matrix<nsd,1> posXiDomain(intpoints.qxg[iquad]);
//std::cout << "position of Gauss point " << iquad << " is "<< posXiDomain << std::endl;
      // coordinates of the current integration point in element coordinates \xi
      //static LINALG::Matrix<nsd,1> posXiDomain;
      //GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
      const double detcell = 1.0; // GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
#ifdef DEBUG
      if (detcell < 0.0)
      {
        std::cout << "detcell :  " << detcell << std::endl;
        dserror("negative detcell!");
      }
#endif
      //----------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //----------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //-------------------------------------
      // procedures involving Jacobian matrix
      //-------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      //-------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //-------------------------------------------------------------------
      static LINALG::Matrix<3,numnode> derxy;
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSize,1> enr_funct;
        enr_funct.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy;
        enr_derxy.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2;
        enr_derxy2.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,
            funct,
            derxy,
            derxy2,
            enr_funct,
            enr_derxy,
            enr_derxy2);

//std::cout << "funct " << funct << std::endl;
//std::cout << "enr_funct " << enr_funct << std::endl;

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shp.d0(iparam) = enr_funct(iparam);
          shp.dx(iparam) = enr_derxy(0,iparam);
          shp.dy(iparam) = enr_derxy(1,iparam);
          shp.dz(iparam) = enr_derxy(2,iparam);
          shp.dxdx(iparam) = enr_derxy2(0,iparam);
          shp.dxdy(iparam) = enr_derxy2(3,iparam);
          shp.dxdz(iparam) = enr_derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = enr_derxy2(1,iparam);
          shp.dydz(iparam) = enr_derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = enr_derxy2(2,iparam);
        }
      }
      else // not xfem_assembly i.e. standard assembly
      {
        dserror("das gibts gar nicht");
        // fill approximation functions for standard FEM
        // remark: numparamvelx == numnode, for standard FEM
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shp.d0(iparam) = funct(iparam);
          shp.dx(iparam) = derxy(0,iparam);
          shp.dy(iparam) = derxy(1,iparam);
          shp.dz(iparam) = derxy(2,iparam);
          shp.dxdx(iparam) = derxy2(0,iparam);
          shp.dxdy(iparam) = derxy2(3,iparam);
          shp.dxdz(iparam) = derxy2(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = derxy2(1,iparam);
          shp.dydz(iparam) = derxy2(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = derxy2(2,iparam);
        }
      }
      //------------------------------------------------------
      // get velocities and accelerations at integration point
      //------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);
      // get acceleration at n+alpha_m for generalized alpha time integration scheme
      LINALG::Matrix<3,1> gpaccam(true);
      if (genalpha)
        gpaccam  = XFEM::interpolateVectorFieldToIntPoint(eaccam , shp.d0, numparamvelx);

      //----------------------------------------------
      // get velocity derivatives at integration point
      //----------------------------------------------
      static LINALG::Matrix<3,nsd> vderxy;
      vderxy.Clear();
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy(isd,0) += evelnp(isd,iparam) * shp.dx(iparam);
          vderxy(isd,1) += evelnp(isd,iparam) * shp.dy(iparam);
          vderxy(isd,2) += evelnp(isd,iparam) * shp.dz(iparam);
        }
      }

      //-----------------------------------------------------
      // get second velocity derivatives at integration point
      //-----------------------------------------------------
      static LINALG::Matrix<3,6> vderxy2;
      vderxy2.Clear();
      //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
      if (higher_order_ele)
      {
        for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          for (size_t isd = 0; isd < nsd; ++isd)
          {
            vderxy2(isd,0) += evelnp(isd,iparam) * shp.dxdx(iparam);
            vderxy2(isd,1) += evelnp(isd,iparam) * shp.dydy(iparam);
            vderxy2(isd,2) += evelnp(isd,iparam) * shp.dzdz(iparam);
            vderxy2(isd,3) += evelnp(isd,iparam) * shp.dxdy(iparam);
            vderxy2(isd,4) += evelnp(isd,iparam) * shp.dxdz(iparam);
            vderxy2(isd,5) += evelnp(isd,iparam) * shp.dydz(iparam);
          }
        }
      }

      //----------------------------------
      // get pressure at integration point
      //----------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp.d0(iparam);

      //--------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp.dx(iparam);
        gradp(1) += eprenp(iparam) * shp.dy(iparam);
        gradp(2) += eprenp(iparam) * shp.dz(iparam);
      }

      //------------------------------------------
      // compute stabilization parameters (3 taus)
      //------------------------------------------
      const double vel_norm = gpvelnp.Norm2();
      //const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_Mu = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;
      // remark: stabilization parameters 'tau_stab_Mu' and 'tau_stab_Mp' are both computed based on
      //         'hk' instead of using the streamlength for 'tau_stab_Mu'
      COMBUST::UTILS::computeStabilizationParams(gpvelnp, xji,
          instationary, dynvisc, densaf, vel_norm, hk, hk, mk, timefac, dt,
          tautype, tau_stab_Mu, tau_stab_Mp, tau_stab_C);

//std::cout << "tau_stab_Mu " << tau_stab_Mu << std::endl;
//std::cout << "tau_stab_Mp " << tau_stab_Mp << std::endl;
//std::cout << "tau_stab_C " << tau_stab_C << std::endl << std::endl;

      //------------------------------------------------------------------------
      // build RHS vector at integration point (scaled with timefac and density)
      //------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> rhsint;
      rhsint.Clear();

      // time RHS vector from previous time step
      // remark: histvec = \rho^{n}*(u^n + delta t (1-theta) a^n)
      LINALG::Matrix<nsd,1> histvec = COMBUST::TIMEINT::GetOldPartOfRighthandside(
          gpveln, gpvelnm, gpaccn, timealgo, dt, theta);
      histvec.Scale(densn);

      // body force vector
      // remark: unit is acceleration (not density weighted)
      static LINALG::Matrix<nsd,1> bodyforce;
      bodyforce.Clear();
      static LINALG::Matrix<nsd,1> noscalebodyforce;
      noscalebodyforce.Clear();
      for (std::size_t isd = 0; isd < nsd; isd++)
      {
        for (std::size_t inode = 0; inode < numnode; inode++)
        {
          bodyforce(isd) += edead(isd,inode) * funct(inode);
          noscalebodyforce(isd) += epgrad(isd,inode) * funct(inode);
        }
      }
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densaf);

      for (size_t isd = 0; isd < nsd; ++isd)
        rhsint(isd) = histvec(isd) + timefac*(bodyforce(isd) + noscalebodyforce(isd));

      //-----------------------------------------------------------------------
      // compute numerical representation of some single operators (not scaled)
      //-----------------------------------------------------------------------

      // convective term:  u_old * grad u_old:
      // remark: conv_old = vderxy(i,j)*gpvelnp(j);
      LINALG::Matrix<nsd,1> conv_old(true);
      conv_old.Multiply(vderxy,gpvelnp);

      // viscous term:  div epsilon(u_old)
      LINALG::Matrix<nsd,1> visc_old;
      visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
      visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
      visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

      // convective part of linearized convective term: u_old * grad (funct)
      // remark: u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      static LINALG::Matrix<shpVecSize,1> enr_conv_c;
      enr_conv_c.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_conv_c(iparam) += shp.dx(iparam)*gpvelnp(0);
        enr_conv_c(iparam) += shp.dy(iparam)*gpvelnp(1);
        enr_conv_c(iparam) += shp.dz(iparam)*gpvelnp(2);
      }

      // viscous term  - grad * epsilon(u)
      /*    /                                              \
                 |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
               1 |                                                |
               - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
               2 |                                                |
                 |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
                  \                                              /
       */
      static COMBUST::EnrViscs2<shpVecSize> enr_viscs2;
      //            enr_viscs2.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_viscs2.xx(iparam) = 0.5 * (2.0 * shp.dxdx(iparam) + shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.xy(iparam) = 0.5 *  shp.dydx(iparam);
        enr_viscs2.xz(iparam) = 0.5 *  shp.dzdx(iparam);
        enr_viscs2.yx(iparam) = 0.5 *  shp.dxdy(iparam);
        enr_viscs2.yy(iparam) = 0.5 * (shp.dxdx(iparam) + 2.0 * shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.yz(iparam) = 0.5 *  shp.dzdy(iparam);
        enr_viscs2.zx(iparam) = 0.5 *  shp.dxdz(iparam);
        enr_viscs2.zy(iparam) = 0.5 *  shp.dydz(iparam);
        enr_viscs2.zz(iparam) = 0.5 * (shp.dxdx(iparam) + shp.dydy(iparam) + 2.0 * shp.dzdz(iparam));
      }

      //----------------------------------------------------
      // evaluate residual (scaled with timefac and density)
      // remark: used for stabilization right hand sides
      //----------------------------------------------------
      //            LINALG::Matrix<nsd,1> res_old(true);
      //            if (instationary)
      //            {
      //              // timefac * (\rho*du/dt + (\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              for (size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = densaf*gpvelnp(isd) - rhsint(isd) + timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));
      //            }
      //            else // stationary
      //            {
      //              // timefac * ((\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              // remark: Although calculation is stationary, scaling with timefac is kept to be
      //              //         consistent with formulation of matrix entries, which are also scaled.
      //              for(size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd)-bodyforce(isd));
      //            }

      // evaluate residual once for all stabilisation right hand sides
      LINALG::Matrix<nsd,1> res_old(true);
      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = dt*(densn*gpaccam(isd)+densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd)-densaf*bodyforce(isd));
        // for rho(t): dt*(densam*gpaccam(isd)+densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd)-densaf*bodyforce(isd));
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpvelnp(isd);
        }
      }
      // integration factors and coefficients of single terms
      const double timefacfac = timefac * fac;
      const double rhsfac     = timefacfac / ga_alphaF;

      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegralsStab<DISTYPE,ASSTYPE,NUMDOF,shpVecSize>(
          assembler, shp, shp, fac, timefac, timefacfac, rhsfac, ga_alphaM, densaf, dynvisc,
          gpvelnp, pres, gpaccam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_conv_c, enr_viscs2,
          instationary, genalpha, newton, pstab, supg, graddiv,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
    } // end loop over gauss points
  //} // end loop over integration cells

  return;
}

/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType CELLDISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          class M1, class M2, class M3, class M4, class V1, class V2>
void Nitsche_SysmatBoundaryCell(
    GEO::BoundaryIntCell           cell,
    const DRT::UTILS::GaussRule2D& intrule2D,
    M2&                            xyze,
    const DRT::ELEMENTS::Combust3* ele,                ///< the element those matrix is calculated
    const bool                     higher_order_ele,
    const M1&                      evelnp,
    const V1&                      eprenp,
    const M3&                      grad_ephi_smoothed, ///< smoothed nodal G-function gradient values for this element
    const M4&                      grad_ephi2_smoothed, ///< smoothed nodal G-function second derivative values for this element
    const V2&                      ecurv,
    const double                   flamespeed,
    const double                   marksteinlength,
    const double                   timefac,
    const bool                     instationary,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
#ifndef COMBUST_SXFEM
    const XFEM::ElementEnrichmentValues&       enrvals_plus,
    const XFEM::ElementEnrichmentValues&       enrvals_minus,
#else
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const V2&                      ephi,         ///< nodal G-function values for this element
#endif
    const double                               kappa_plus,          ///< weights of average in plus domain
    const double                               kappa_minus,         ///< weights of average in minus domain
    const double                               visc_plus,           ///< viscosity in plus domain
    const double                               visc_minus,          ///< viscosity in minus domain
    const double                               dens_plus,           ///< density in plus domain
    const double                               dens_minus,          ///< density in minus domain
    const double                               alphau,              ///< Nitsche's stabilization parameter
    const double                               wallfac,             ///< wall factor to fade jump at walls (for ORACLES)
    const double                               hk_ele,              ///< element length
    const double                               nitschevel,          ///< Nitsche's velocity
    const size_t                               numparamvelx,
    const size_t                               numparampres,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox,      ///<
    const double                               variablesurftens,
    const bool                                 second_deriv,
    const double                               surftenscoeff,       ///<
    const bool                                 connected_interface, ///<
    const INPAR::COMBUST::VelocityJumpType     veljumptype,
    const INPAR::COMBUST::FluxJumpType         fluxjumptype,
    const bool                                 smoothed_boundary_integration,
    const bool                                 nitsche_convflux,
    const bool                                 nitsche_convstab,
    const bool                                 nitsche_convpenalty,
    const bool                                 nitsche_mass
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  //----------------------------------------------
  // get Gaussian points for this integration cell
  //----------------------------------------------
  const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

  // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
  LINALG::SerialDenseMatrix cellXiDomaintmp = cell.CellNodalPosXiDomain();
  // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
  LINALG::SerialDenseMatrix cellXYZDomaintmp = cell.CellNodalPosXYZ();
  // transform to fixed size format
  const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
  const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

  // assembly line integrals for surface tension if interface is not connected
  // and has intersection with domain boundary
  // no line integrals for not laplace-beltrami surface tension approximation types
  if(!connected_interface
      && (   fluxjumptype == INPAR::COMBUST::flux_jump_surface_tension)
      && ( surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami ||
           surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed )
  ) // only for laplace_beltrami approximation
  {
#ifndef COMBUST_SXFEM
    COMBUST::RHS_LaplaceBeltrami_1DLineIntegrals_Nitsche<DISTYPE,CELLDISTYPE,ASSTYPE,NUMDOF>(
        cellXiDomain, xyze, grad_ephi_smoothed,
        enrvals_plus,enrvals_minus,
        kappa_plus,kappa_minus,numparamvelx,timefac,surftensapprox,surftenscoeff,assembler);
#else
    std::cout << "RHS_LaplaceBeltrami_1DLineIntegrals_Nitsche deactivated - dont do two-phase flow intersecting the boundary!" << std::endl;
#endif
  }

  //--------------------------------------
  // integration loop over Gaussian points
  //--------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // new transformation for boundary integrals
    // 1. define a coupled transformation x_3D(xi_3D(eta_2D)): transformation from 2D->3D
    // 2. compute the corresponding Jacobian J_eta2D->x_3D and
    // 3. the corresponding surface integral factor sqrt(det(J_eta2D->x_3D^T * J_eta2D->x_3D))
    // 4. approximate integral with Gauss rule in eta coordinate system
    // 5. evaluate the transformed integrand f(x(xi(eta)))

    const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

    // jacobian for coupled transformation
    // get derivatives dxi_3D/deta_2D
    static LINALG::Matrix<2,numvertices> deriv_eta2D;
    DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),CELLDISTYPE);

    // calculate dxi3Ddeta2D
    static LINALG::Matrix<3,2> dXi3Ddeta2D;
    dXi3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < (int)numvertices; k++)
          dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);

    // transform Gauss point to xi3D space (element parameter space)
    static LINALG::Matrix<nsd,1> gpinXi3D;
    gpinXi3D.Clear();
    // coordinates of this integration point in element coordinates \xi^domain
    GEO::mapEtaBToXiD(cell, gpinEta2D, gpinXi3D);

    static LINALG::Matrix<3,numnode> deriv_xi3D;
    DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DISTYPE);

    // calculate dx3Ddxi3D
    static LINALG::Matrix<3,3> dX3DdXi3D;
    dX3DdXi3D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 3; j++) // derivatives
        for (int k = 0; k < (int)numnode; k++)
          dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

    // get the coupled Jacobian dx3Ddeta2D
    static LINALG::Matrix<3,2> dx3Ddeta2D;
    dx3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < 3; k++)
          dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);

    // get deformation factor
    static LINALG::Matrix<2,2> Jac_tmp; // J^T*J
    Jac_tmp.Clear();
    Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

    if(Jac_tmp.Determinant() <= 0.0) dserror("deformation factor for boundary integration is zero");

    const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

    const double fac = intpoints.qwgt[iquad]*deform_factor;

    LINALG::Matrix<2,1> posEtaBoundary;
    posEtaBoundary.Clear();
    for (int i= 0; i< 2; i++){
      posEtaBoundary(i,0) = gpinEta2D(i,0);
    }

    LINALG::Matrix<nsd,1> posXiDomain;
    posXiDomain.Clear();
    for (int i= 0; i< 3; i++){
      posXiDomain(i,0) = gpinXi3D(i,0);
    }

    // compute total (time and spatial) integration factor (and coefficients of single terms)?
    const double timefacfac = timefac * fac;

    //--------------------------------------------------------------------------------------------
    // compute normal vector (normed, in physical coordinates)
    // remark: for linear boundary integrastion cells this could be done before the loop over all
    //         Gaussian points
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    static LINALG::Matrix<nsd,1> tang1(true);
    static LINALG::Matrix<nsd,1> tang2(true);
    normal.Clear();
    tang1.Clear();
    tang2.Clear();
    // normal vector points from Omega^+ to Omega^- (n = n^+ = -n^-)
    GEO::computeNormalAndTangentToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal, tang1, tang2);

    static LINALG::Matrix<nsd,1> normal_smooth(true);
    static LINALG::Matrix<nsd,1> tang1_smooth(true);
    static LINALG::Matrix<nsd,1> tang2_smooth(true);
    // arbitrary auxiliary vector
    static LINALG::Matrix<nsd,1> arbvec(true);

    normal_smooth.Clear();
    tang1_smooth.Clear();
    tang2_smooth.Clear();
    arbvec.Clear();

    static double curv;
    curv = 0.0;

    // call smoothed normals for an accurate evaluation of boundary integration terms
    if ((smoothed_boundary_integration == true) and
        (surftensapprox != INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed
         and surftensapprox != INPAR::COMBUST::surface_tension_approx_divgrad_normal
         and surftensapprox != INPAR::COMBUST::surface_tension_approx_nodal_curvature))
    {
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      // loop over nodes
      for(size_t i = 0; i< numnode; i++)
      {
        LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
        nodal_grad_xyz_tmp.Clear();

        // compute smoothed normal vector
        normal_smooth(0)+=grad_ephi_smoothed(0,i)*funct_gp(i);
        normal_smooth(1)+=grad_ephi_smoothed(1,i)*funct_gp(i);
        normal_smooth(2)+=grad_ephi_smoothed(2,i)*funct_gp(i);

        // interpolate curvature to Gaussian point
        curv += ecurv(i)*funct_gp(i);
      }
#ifdef COMBUST_2D
      normal_smooth(2) = 0.0;
#endif
      // normalize final gradient
      // remark: n = -grad(\phi)
      const double norm = normal_smooth.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal_smooth.Scale(-1.0/norm);

      arbvec(0) = 1.0;
      arbvec(1) = 0.0;
      arbvec(2) = 0.0;
      // check if normal and tang1 are almost parallel
      if (fabs(arbvec(0)*normal_smooth(0)+arbvec(1)*normal_smooth(1)+arbvec(2)*normal_smooth(2))>0.5)
      {
        arbvec(0) = 0.0;
        arbvec(1) = 1.0;
        arbvec(2) = 0.0;
      }
      tang1_smooth(0) = arbvec(1)*normal_smooth(2) - arbvec(2)*normal_smooth(1);
      tang1_smooth(1) = arbvec(2)*normal_smooth(0) - arbvec(0)*normal_smooth(2);
      tang1_smooth(2) = arbvec(0)*normal_smooth(1) - arbvec(1)*normal_smooth(0);

      tang2_smooth(0) = normal_smooth(1)*tang1_smooth(2) - normal_smooth(2)*tang1_smooth(1);
      tang2_smooth(1) = normal_smooth(2)*tang1_smooth(0) - normal_smooth(0)*tang1_smooth(2);
      tang2_smooth(2) = normal_smooth(0)*tang1_smooth(1) - normal_smooth(1)*tang1_smooth(0);

      // scale
      // remark: tang1_smooth is already normed
      const double norm_tang1 = tang1_smooth.Norm2();
      const double norm_tang2 = tang2_smooth.Norm2();
      if (norm_tang1 == 0.0) dserror("norm of tangential vector 1 is zero!");
      if (norm_tang2 == 0.0) dserror("norm of tangential vector 2 is zero!");
      tang1_smooth.Scale(1.0 / norm_tang1);
      tang2_smooth.Scale(1.0 / norm_tang2);
    }
    else
    {
      normal_smooth = normal;
      tang1_smooth = tang1;
      tang2_smooth = tang2;
    }

#ifdef COLLAPSE_FLAME
      normal_smooth.Clear();
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      for (unsigned i=0;i<numnode;i++)
      {
        normal_smooth(0) += funct_gp(i)*xyze(0,i);
        normal_smooth(1) += funct_gp(i)*xyze(1,i);
        normal_smooth(2) += funct_gp(i)*xyze(2,i);
      }
#ifdef COMBUST_2D
      normal_smooth(2) = 0.0;
#endif
      const double norm = normal_smooth.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal_smooth.Scale(1.0/norm);
#endif

      //-----------------------------
      // prepare interface conditions
      //-----------------------------
      // mass flux through flame front
      const double mflux = flamespeed*(1.0 - marksteinlength*curv)*dens_minus*wallfac;

      // velocity jump in normal direction
#if 0
      // jump velocity vector ju_{n,t1,t2}
      LINALG::Matrix<nsd,1> ju_nt(true);
      ju_nt(0,0) = -mflux*(1.0/dens_minus - 1.0/dens_plus);
      ju_nt(1,0) = 0.0;
      ju_nt(2,0) = 0.0;
#else
      double ju_n = 0.0;
      if (veljumptype == INPAR::COMBUST::vel_jump_premixed_combustion)
      {
        ju_n = -mflux*(1.0/dens_minus - 1.0/dens_plus);
      }
      else if(veljumptype == INPAR::COMBUST::vel_jump_const)
      {
        ju_n = 5.0;
      }
      else
      {
        // velocity jump is zero
      }
#endif

      // jump in flux (normal stress)
      double jflux_n = 0.0;//mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);
      if(fluxjumptype == INPAR::COMBUST::flux_jump_premixed_combustion){
        jflux_n = mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);
      }
      else if(fluxjumptype == INPAR::COMBUST::flux_jump_const){
        jflux_n = 3.0;
      }
      else{
        // flux jump is zero
      }

      // jump in velocity
      LINALG::Matrix<nsd,1> ju(true);
#if 0
    //----------------------------------------------------------------------------------------
    // remark: does not make much sense; multiplication of 'ju' by 'normal' yields same result
    //----------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> tang1(true);
    tang1.Clear();
    tang1(0) = cellXYZDomain(0,1) - cellXYZDomain(0,0);
    tang1(1) = cellXYZDomain(1,1) - cellXYZDomain(1,0);
    tang1(2) = cellXYZDomain(2,1) - cellXYZDomain(2,0);
    // compute unit (normed) tang1 vector
    double norm = sqrt(tang1(0)*tang1(0) + tang1(1)*tang1(1) + tang1(2)*tang1(2));
    if (norm == 0.0) dserror("norm of tang1 vector is zero!");
    tang1.Scale(1.0/norm);

    static LINALG::Matrix<nsd,1> tang2(true);
    tang2.Clear();
    tang2(0) = cellXYZDomain(0,2) - cellXYZDomain(0,0);
    tang2(1) = cellXYZDomain(1,2) - cellXYZDomain(1,0);
    tang2(2) = cellXYZDomain(2,2) - cellXYZDomain(2,0);
    // compute unit (normed) tang1 vector
    norm = sqrt(tang2(0)*tang2(0) + tang2(1)*tang2(1) + tang2(2)*tang2(2));
    if (norm == 0.0) dserror("norm of tang1 vector is zero!");
    tang2.Scale(1.0/norm);

    // coordinate transform of ju vector
    //              | n_1, t1_1, t2_1 |
    // ju_{x,y,z} = | n_2, t1_2, t2_2 | * ju_{n,t1,t2}
    //              | n_3, t1_3, t2_3 |

    LINALG::Matrix<nsd,nsd> trafo(true);
    for (size_t idim=0;idim<nsd;idim++)
    {
      trafo(idim,0) = normal(idim,0);
      trafo(idim,1) = tang1(idim,0);
      trafo(idim,2) = tang2(idim,0);
    }
    //ju = trafo*ju_nt;
    ju.Multiply(trafo,ju_nt);
#else
    if (veljumptype == INPAR::COMBUST::vel_jump_premixed_combustion or
        veljumptype == INPAR::COMBUST::vel_jump_const)
    {
      // TODO: macht das wirklich Sinn mit dem geglaetten Normalenvektor? Ursula
      // use smoothes normal for velocity jump
      ju(0) = ju_n*normal_smooth(0);
      ju(1) = ju_n*normal_smooth(1);
      ju(2) = ju_n*normal_smooth(2);
    //std::cout << "ju " << ju << std::endl;
    }
    //else: two-phase flow with zero jump
    //std::cout << "ju " << ju << std::endl;
#ifdef COMBUST_TESTCOUETTEFLOW
    // overwrite velocity jump condition for shear flow test case
    ju(0) = -2.0;
    ju(1) = 0.0;
    ju(2) = 0.0;
#endif
#endif

    // jump of traction vector (stress tensor * normal)
    LINALG::Matrix<nsd,1> jflux(true);
    if(fluxjumptype == INPAR::COMBUST::flux_jump_const ||
       fluxjumptype == INPAR::COMBUST::flux_jump_premixed_combustion)
    {
      jflux(0) = jflux_n*normal(0);
      jflux(1) = jflux_n*normal(1);
      jflux(2) = jflux_n*normal(2);
      //std::cout << "jflux " << jflux << std::endl;
      // overwrite velocity jump condition (e.g. for specific test case)
      //jflux(0) = 0.0;
      //jflux(1) = 0.0;
      //jflux(2) = 0.0;
    }

    //----------------------------------------------------------------------------
    // evaluate shape functions and their first derivatives at this Gaussian point
    //----------------------------------------------------------------------------
    // shape functions for fields with dofs at nodes
    static LINALG::Matrix<numnode,1> funct;
    funct.Clear();
    DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
    // derivatives of shape functions for fields with dofs at nodes
    static LINALG::Matrix<nsd,numnode> deriv;
    deriv.Clear();
    DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    //--------------------------------------------------------
    // procedures involving Jacobian matrix for domain mapping
    //--------------------------------------------------------
    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<nsd,nsd> xjm;
    xjm.Clear();
    xjm.MultiplyNT(deriv,xyze);
    // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
    const double detXtoXi = xjm.Determinant();

    if (detXtoXi < 0.0)
    {
      dserror("negative Jacobian determinant detXtoXi: %f for element %d!", detXtoXi, ele->Id());
    }

    // inverse of jacobian
    static LINALG::Matrix<nsd,nsd> xji;
    xji.Clear();
    xji.Invert(xjm);

    //-------------------------------------------------------------------
    // compute global derivates of shape functions at this Gaussian point
    //-------------------------------------------------------------------
    // compute first global derivative
    static LINALG::Matrix<3,numnode> derxy;
    derxy.Clear();
    // derxy(i,j) = xji(i,k) * deriv(k,j)
    derxy.Multiply(xji,deriv);

    // compute second global derivative
    static LINALG::Matrix<6,numnode> derxy2;
    if (higher_order_ele)
    {
      static LINALG::Matrix<6,numnode> deriv2;
      DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
    }
    else
    {
      derxy2.Clear();
    }

#ifdef COMBUST_SXFEM
    const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi,funct,derxy,derxy2);
    const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi,funct,derxy,derxy2);
    //const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
    //const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);
#endif
    //--------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------
    const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
    static LINALG::Matrix<shpVecSize,1> enrfunct_plus;
    enrfunct_plus.Clear();
    static LINALG::Matrix<shpVecSize,1> enrfunct_minus;
    enrfunct_minus.Clear();

    static LINALG::Matrix<3,shpVecSize> enrderxy_plus;
    enrderxy_plus.Clear();
    static LINALG::Matrix<3,shpVecSize> enrderxy_minus;
    enrderxy_minus.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_plus;
    enrderxy2_plus.Clear();
    static LINALG::Matrix<6,shpVecSize> enrderxy2_minus;
    enrderxy2_minus.Clear();

    // shape functions for nodal parameters (dofs) on plus and minus side
    enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
        enrfunct_plus, enrderxy_plus, enrderxy2_plus);
    enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
        enrfunct_minus, enrderxy_minus, enrderxy2_minus);

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<0,shpVecSize> shp_jump;      // [[ ]] notation
    static XFEM::ApproxFunc<0,shpVecSize> shp_mean;      //  { }  notation
    static XFEM::ApproxFunc<1,shpVecSize> shp_mean_inv;  //  < >  notation
    static XFEM::ApproxFunc<2,shpVecSize> shp_mean_visc;
    static XFEM::ApproxFunc<2,shpVecSize> shp_jump_visc;

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_jump.d0(iparam) = enrfunct_minus(iparam) - enrfunct_plus(iparam);

      shp_mean.d0(iparam) = kappa_plus * enrfunct_plus(iparam) + kappa_minus * enrfunct_minus(iparam);

      shp_mean_inv.d0(iparam) = kappa_minus * enrfunct_plus(iparam) + kappa_plus * enrfunct_minus(iparam);

      shp_mean_inv.dx(iparam) = kappa_minus * enrderxy_plus(0,iparam) + kappa_plus * enrderxy_minus(0,iparam);
      shp_mean_inv.dy(iparam) = kappa_minus * enrderxy_plus(1,iparam) + kappa_plus * enrderxy_minus(1,iparam);
      shp_mean_inv.dz(iparam) = kappa_minus * enrderxy_plus(2,iparam) + kappa_plus * enrderxy_minus(2,iparam);

      shp_mean_visc.d0(iparam) = kappa_plus * visc_plus*enrfunct_plus(iparam) + kappa_minus * visc_minus*enrfunct_minus(iparam);
      shp_mean_visc.dx(iparam) = kappa_plus * visc_plus*enrderxy_plus(0,iparam) + kappa_minus * visc_minus*enrderxy_minus(0,iparam);
      shp_mean_visc.dy(iparam) = kappa_plus * visc_plus*enrderxy_plus(1,iparam) + kappa_minus * visc_minus*enrderxy_minus(1,iparam);
      shp_mean_visc.dz(iparam) = kappa_plus * visc_plus*enrderxy_plus(2,iparam) + kappa_minus * visc_minus*enrderxy_minus(2,iparam);
      shp_mean_visc.dxdx(iparam) = kappa_plus * visc_plus*enrderxy2_plus(0,iparam) + kappa_minus * visc_minus*enrderxy2_minus(0,iparam);
      shp_mean_visc.dxdy(iparam) = kappa_plus * visc_plus*enrderxy2_plus(3,iparam) + kappa_minus * visc_minus*enrderxy2_minus(3,iparam);
      shp_mean_visc.dxdz(iparam) = kappa_plus * visc_plus*enrderxy2_plus(4,iparam) + kappa_minus * visc_minus*enrderxy2_minus(4,iparam);
      shp_mean_visc.dydx(iparam) = shp_mean_visc.dxdy(iparam);
      shp_mean_visc.dydy(iparam) = kappa_plus * visc_plus*enrderxy2_plus(1,iparam) + kappa_minus * visc_minus*enrderxy2_minus(1,iparam);
      shp_mean_visc.dydz(iparam) = kappa_plus * visc_plus*enrderxy2_plus(5,iparam) + kappa_minus * visc_minus*enrderxy2_minus(5,iparam);
      shp_mean_visc.dzdx(iparam) = shp_mean_visc.dxdz(iparam);
      shp_mean_visc.dzdy(iparam) = shp_mean_visc.dydz(iparam);
      shp_mean_visc.dzdz(iparam) = kappa_plus * visc_plus*enrderxy2_plus(2,iparam) + kappa_minus * visc_minus*enrderxy2_minus(2,iparam);

      shp_jump_visc.d0(iparam) = visc_minus*enrfunct_minus(iparam) - visc_plus*enrfunct_plus(iparam);
      shp_jump_visc.dx(iparam) = visc_minus*enrderxy_minus(0,iparam) - visc_plus*enrderxy_plus(0,iparam);
      shp_jump_visc.dy(iparam) = visc_minus*enrderxy_minus(1,iparam) - visc_plus*enrderxy_plus(1,iparam);
      shp_jump_visc.dz(iparam) = visc_minus*enrderxy_minus(2,iparam) - visc_plus*enrderxy_plus(2,iparam);
      shp_jump_visc.dxdx(iparam) = visc_minus*enrderxy2_minus(0,iparam) - visc_plus*enrderxy2_plus(0,iparam);
      shp_jump_visc.dxdy(iparam) = visc_minus*enrderxy2_minus(3,iparam) - visc_plus*enrderxy2_plus(3,iparam);
      shp_jump_visc.dxdz(iparam) = visc_minus*enrderxy2_minus(4,iparam) - visc_plus*enrderxy2_plus(4,iparam);
      shp_jump_visc.dydx(iparam) = shp_jump_visc.dxdy(iparam);
      shp_jump_visc.dydy(iparam) = visc_minus*enrderxy2_minus(1,iparam) - visc_plus*enrderxy2_plus(1,iparam);
      shp_jump_visc.dydz(iparam) = visc_minus*enrderxy2_minus(5,iparam) - visc_plus*enrderxy2_plus(5,iparam);
      shp_jump_visc.dzdx(iparam) = shp_jump_visc.dxdz(iparam);
      shp_jump_visc.dzdy(iparam) = shp_jump_visc.dydz(iparam);
      shp_jump_visc.dzdz(iparam) = visc_minus*enrderxy2_minus(2,iparam) - visc_plus*enrderxy2_plus(2,iparam);
    }
    //---------------------------------------
    // compute data at Gaussian point for rhs
    //---------------------------------------

    // velocity jump
    static LINALG::Matrix<nsd,1> veljump(true);
    veljump.Clear();
    veljump = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);

    // mean velocity
    static LINALG::Matrix<nsd,1> velmean(true);
    velmean.Clear();
    velmean = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_mean.d0, numparamvelx);

    // get velocity (np,i) derivatives at integration point
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    static LINALG::Matrix<nsd,nsd> vderxy_mean_visc(true);
    vderxy_mean_visc.Clear();
    for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      for (size_t isd = 0; isd < nsd; ++isd)
      {
        vderxy_mean_visc(isd,0) += evelnp(isd,iparam) * shp_mean_visc.dx(iparam);
        vderxy_mean_visc(isd,1) += evelnp(isd,iparam) * shp_mean_visc.dy(iparam);
        vderxy_mean_visc(isd,2) += evelnp(isd,iparam) * shp_mean_visc.dz(iparam);
      }

    // get velocity (np,i) derivatives at integration point
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    static LINALG::Matrix<nsd,nsd> vderxy_jump_visc(true);
    vderxy_jump_visc.Clear();
    for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      for (size_t isd = 0; isd < nsd; ++isd)
      {
        vderxy_jump_visc(isd,0) += evelnp(isd,iparam) * shp_jump_visc.dx(iparam);
        vderxy_jump_visc(isd,1) += evelnp(isd,iparam) * shp_jump_visc.dy(iparam);
        vderxy_jump_visc(isd,2) += evelnp(isd,iparam) * shp_jump_visc.dz(iparam);
      }

    // get pressure jump
    static double pjump;
    pjump = 0.0;
    for (size_t iparam = 0; iparam != numparampres; ++iparam)
      pjump += shp_jump.d0(iparam)*eprenp(iparam);

    static double pmean;
    pmean = 0.0;
    for (size_t iparam = 0; iparam != numparampres; ++iparam)
      pmean += shp_mean.d0(iparam)*eprenp(iparam);

    static XFEM::ApproxFunc<0,shpVecSize> shp_plus;
    static XFEM::ApproxFunc<0,shpVecSize> shp_minus;

    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_plus.d0(iparam) = enrfunct_plus(iparam);
      shp_minus.d0(iparam) = enrfunct_minus(iparam);
    }

    // inflow velocity
    LINALG::Matrix<3,1> vel_minus(true);
    for (unsigned iparam = 0; iparam < numparamvelx; ++iparam)
      for (int isd = 0; isd < 3; ++isd)
        vel_minus(isd) += evelnp(isd,iparam)*enrfunct_minus(iparam);

    LINALG::Matrix<3,1> vel_plus(true);
    for (unsigned iparam = 0; iparam < numparamvelx; ++iparam)
      for (int isd = 0; isd < 3; ++isd)
        vel_plus(isd) += evelnp(isd,iparam)*enrfunct_plus(iparam);

    double velproj_plus  = vel_plus(0)*normal(0) + vel_plus(1)*normal(1) + vel_plus(2)*normal(2);
    double velproj_minus = vel_minus(0)*normal(0) + vel_minus(1)*normal(1) + vel_minus(2)*normal(2);
    const double absvelproj = fabs(velproj_plus);

    if (velproj_plus > 0.0) // outflow plus domain
      velproj_plus = 0.0;
    if (velproj_minus < 0.0) // outflow minus domain
      velproj_minus = 0.0;

    //--------------------------------------------------------------------------------------------
    // build single boundary integral stiffness terms
    //
    // remarks: || x || stands for the jump operator ( || x || = x^- - x^+ )
    //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
    //             n    stands for the normal vector on the interface pointing from the positive (+)
    //                  domain into the negative (-) domain (n = n^+ = -n^-)
    //--------------------------------------------------------------------------------------------

    const double halftimefacfacplus   = dens_plus*0.5*timefacfac;
    const double halftimefacfacminus  = dens_minus*0.5*timefacfac;

    double alphatimefacfac = 0.0;
    if (not nitsche_mass)
      alphatimefacfac = alphau*timefacfac;
    else
    {
      // max norm of velocity for plus and minus domain at Gauss point
      double vel_max_norm_plus = 0.0;
      for (std::size_t idim=0; idim<3; idim++)
      {
        if (std::abs(vel_plus(idim,0)) > vel_max_norm_plus)
          vel_max_norm_plus = std::abs(vel_plus(idim,0));
      }
      double vel_max_norm_minus = 0.0;
      for (std::size_t idim=0; idim<3; idim++)
      {
        if (std::abs(vel_minus(idim,0)) > vel_max_norm_minus)
          vel_max_norm_minus = std::abs(vel_minus(idim,0));
      }

      if (instationary)
        alphatimefacfac = nitschevel * (kappa_minus* (visc_minus / hk_ele + dens_minus*vel_max_norm_minus / 6.0 + dens_minus * hk_ele / (timefac * 12.0))
                                       +kappa_plus* (visc_plus / hk_ele + dens_plus*vel_max_norm_plus / 6.0 + dens_plus * hk_ele / (timefac * 12.0))) *timefacfac;
      else
        alphatimefacfac = nitschevel * (kappa_minus* (visc_minus / hk_ele + dens_minus*vel_max_norm_minus / 6.0)
                                       +kappa_plus* (visc_plus / hk_ele + dens_plus*vel_max_norm_plus / 6.0)) *timefacfac;
    }

    // TODO: convstabtimefacfac = alpha_c * timefacfac
    //       alpha_c = 0.5 * rho_max * absvelproj or rho_inflow_domain?
    // due to the convection-dependence of the penalty parameter, this term should not be required for the time being
    double rho_max = std::max(dens_plus,dens_minus);
    const double convstabtimefacfac = 0.5*absvelproj*rho_max*timefacfac;
    // TODO: currently, a mass penalty term acting merely in normal direction is not used
#if 0
    const double B1 = 0.05;
    const double B2 = 0.05;
    const double massstabtimefacfac = (B1/hk_ele * 0.5 * (visc_plus+visc_minus)
                                    + B2 * 0.5 * (dens_plus*vel_plus.Norm2()+dens_minus*vel_minus.Norm2())) * timefacfac;
#endif


#ifdef COMBUST_NITSCHE_VISCCONSISTMEAN
    //-------------------------    |                                      |
    // viscous consistency term  + |  || v || , { 2\mu epsilon( Du ) } n  |
    //-------------------------    |                                      |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, 2.0*timefacfac*normal(0), shp_mean_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dx);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, 2.0*timefacfac*normal(1), shp_mean_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dy);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, 2.0*timefacfac*normal(2), shp_mean_visc.dz);

    //   |                                       |
    // - |  || v || , { 2\mu epsilon( u_i ) } n  |
    //   |                                       |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(0,0) + vderxy_mean_visc(0,0))*normal(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(0,1) + vderxy_mean_visc(1,0))*normal(1));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(0,2) + vderxy_mean_visc(2,0))*normal(2));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(1,0) + vderxy_mean_visc(0,1))*normal(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(1,1) + vderxy_mean_visc(1,1))*normal(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(1,2) + vderxy_mean_visc(2,1))*normal(2));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(2,0) + vderxy_mean_visc(0,2))*normal(0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(2,1) + vderxy_mean_visc(1,2))*normal(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(2,2) + vderxy_mean_visc(2,2))*normal(2));
#endif

#ifdef COMBUST_NITSCHE_VISCADJCONSIST
    //---------------------------------    |                                     |
    // viscous adjoint consistency term  + |  { 2\mu epsilon( v ) } n, || Du ||  |
    //---------------------------------    |                                     |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     timefacfac*normal(2), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_mean_visc.dx,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_mean_visc.dx,     timefacfac*normal(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_mean_visc.dy,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     timefacfac*normal(2), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_mean_visc.dy,     timefacfac*normal(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_mean_visc.dz,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_mean_visc.dz,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_visc.dz, 2.0*timefacfac*normal(2), shp_jump.d0);

    //   |                                      |
    // - |  { 2\mu epsilon( v ) } n, || u_i ||  |
    //   |                                      |

    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, -2.0* normal(0)           *timefacfac*veljump(0,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy, -    (normal(0)+normal(1))*timefacfac*veljump(0,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz, -    (normal(0)+normal(2))*timefacfac*veljump(0,0));

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, -2.0*timefacfac*normal(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy, -    timefacfac*normal(1)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz, -    timefacfac*normal(2)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx, -    timefacfac*normal(1)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx, -    timefacfac*normal(2)*veljump(0,0));

    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx, -    (normal(1)+normal(0))*timefacfac*veljump(1,0));
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, -2.0* normal(1)           *timefacfac*veljump(1,0));
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz, -    (normal(1)+normal(2))*timefacfac*veljump(1,0));

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy, -    timefacfac*normal(0)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx, -    timefacfac*normal(0)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, -2.0*timefacfac*normal(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz, -    timefacfac*normal(2)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy, -    timefacfac*normal(2)*veljump(1,0));

    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx, -    (normal(2)+normal(0))*timefacfac*veljump(2,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy, -    (normal(2)+normal(1))*timefacfac*veljump(2,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, -2.0* normal(2)           *timefacfac*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz, -    timefacfac*normal(0)*veljump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz, -    timefacfac*normal(1)*veljump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx, -    timefacfac*normal(0)*veljump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy, -    timefacfac*normal(1)*veljump(2,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, -2.0*timefacfac*normal(2)*veljump(2,0));

    //-------------------------------------    |                                    |
    // viscous adjoint consistency term RHS  + |  { 2\mu epsilon( v ) } n, j_u  |
    //-------------------------------------    |                                    |
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, 2.0* normal(0)           *timefacfac*normal(0)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     (normal(0)+normal(1))*timefacfac*normal(0)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     (normal(0)+normal(2))*timefacfac*normal(0)*ju);

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(1)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     timefacfac*normal(2)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(1)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     timefacfac*normal(2)*ju(0));
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(1)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     timefacfac*normal(2)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(1)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     timefacfac*normal(2)*normal(0)*ju_n);

    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     (normal(1)+normal(0))*timefacfac*normal(1)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, 2.0* normal(1)           *timefacfac*normal(1)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     (normal(1)+normal(2))*timefacfac*normal(1)*ju);

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(0)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(0)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     timefacfac*normal(2)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     timefacfac*normal(2)*ju(1));
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(0)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(0)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     timefacfac*normal(2)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     timefacfac*normal(2)*normal(1)*ju_n);

    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     (normal(2)+normal(0))*timefacfac*normal(2)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     (normal(2)+normal(1))*timefacfac*normal(2)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, 2.0* normal(2)           *timefacfac*normal(2)*ju);

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     timefacfac*normal(0)*ju(2));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     timefacfac*normal(1)*ju(2));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     timefacfac*normal(0)*ju(2));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     timefacfac*normal(1)*ju(2));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, 2.0*timefacfac*normal(2)*ju(2));
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     timefacfac*normal(0)*normal(2)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     timefacfac*normal(1)*normal(2)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     timefacfac*normal(0)*normal(2)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     timefacfac*normal(1)*normal(2)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, 2.0*timefacfac*normal(2)*normal(2)*ju_n);
#endif

#ifdef COMBUST_NITSCHE_PRESCONSISTMEAN
    //--------------------------    |                     |
    // pressure consistency term  - |  || v ||, { Dp } n  |
    //--------------------------    |                     |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Pres>(shp_jump.d0, -timefacfac*normal(0), shp_mean.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Pres>(shp_jump.d0, -timefacfac*normal(1), shp_mean.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Pres>(shp_jump.d0, -timefacfac*normal(2), shp_mean.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Pres>(shp_jump.d0, -timefacfac*normal(0), shp_presproj_mean);
    //      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Pres>(shp_jump.d0, -timefacfac*normal(1), shp_presproj_mean);
    //      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Pres>(shp_jump.d0, -timefacfac*normal(2), shp_presproj_mean);

    //   |                      |
    // + |  || v ||, { p_i } n  |
    //   |                      |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, timefacfac*normal(0)*pmean);
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, timefacfac*normal(1)*pmean);
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, timefacfac*normal(2)*pmean);
#endif

#ifdef COMBUST_NITSCHE_PRESADJCONSIST
    //----------------------------------    |                     |
    // pressure adjoint consistency term  + |  { q }, || Du || n  |
    //----------------------------------    |                     |

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(2), shp_jump.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velx>(shp_presproj_mean, +timefacfac*normal(0), shp_jump.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Vely>(shp_presproj_mean, +timefacfac*normal(1), shp_jump.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velz>(shp_presproj_mean, +timefacfac*normal(2), shp_jump.d0);

    //   |                      |
    // - |  { q }, || u_i || n  |
    //   |                      |

    assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, -timefacfac*(veljump(0,0)*normal(0)
                                                                            +veljump(1,0)*normal(1)
                                                                            +veljump(2,0)*normal(2)));
    //      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_presproj_mean, -timefacfac*(vjump(0,0)*normal(0)
    //          +vjump(1,0)*normal(1)
    //          +vjump(2,0)*normal(2)));

    //--------------------------------------    |                |
    // pressure adjoint consistency term RHS  + |  { q }, j_u*n  |
    //--------------------------------------    |                |

    assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, timefacfac*(ju(0)*normal(0)
                                                                           +ju(1)*normal(1)
                                                                           +ju(2)*normal(2)));
    //assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, timefacfac*ju_n);

    //      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_presproj_mean, +timefacfac*(jump_vel(0,0)*normal(0)
    //          +jump_vel(1,0)*normal(1)
    //          +jump_vel(2,0)*normal(2)));
#endif

#ifdef COMBUST_NITSCHE_FLUXJUMP

    if(fluxjumptype == INPAR::COMBUST::flux_jump_surface_tension)
    {
      //------------------------------------------------
      // assemble boundary integrals for surface tension
      //------------------------------------------------
      if ((smoothed_boundary_integration == true) and
          (surftensapprox != INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed
           and surftensapprox != INPAR::COMBUST::surface_tension_approx_divgrad_normal
           and surftensapprox != INPAR::COMBUST::surface_tension_approx_nodal_curvature))
     {
        COMBUST::RHS_SurfaceTension2DIntegrals<DISTYPE,ASSTYPE,NUMDOF>(posXiDomain, shp_mean_inv, ele,
                    xyze, grad_ephi_smoothed, grad_ephi2_smoothed, ecurv, normal_smooth, timefacfac,
                    surftensapprox, variablesurftens, second_deriv, surftenscoeff, assembler);
      }
      else
        COMBUST::RHS_SurfaceTension2DIntegrals<DISTYPE,ASSTYPE,NUMDOF>(posXiDomain, shp_mean_inv, ele,
            xyze, grad_ephi_smoothed, grad_ephi2_smoothed, ecurv, normal, timefacfac,
            surftensapprox, variablesurftens, second_deriv, surftenscoeff, assembler);
    }
    else
    {
      //---------------------------------------------------    |                 |
      // flux (traction/stress) jump (consistency) term RHS  - |  < v >, j_flux  |
      //---------------------------------------------------    |                 |

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, -timefacfac*jflux(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, -timefacfac*jflux(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, -timefacfac*jflux(2));
    }
#endif

#ifdef COMBUST_NITSCHE_NITSCHEVEL
    //------------------------    |                              |
    // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
    //------------------------    |                              |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, alphatimefacfac, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, alphatimefacfac, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, alphatimefacfac, shp_jump.d0);

    //    |                               |
    //  - |  \alpha_u || v ||, || u_i ||  |
    //    |                               |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -alphatimefacfac*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -alphatimefacfac*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -alphatimefacfac*veljump(2,0));

    //----------------------------    |                         |
    // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, j_u  |
    //----------------------------    |                         |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, alphatimefacfac*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, alphatimefacfac*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, alphatimefacfac*ju(2));
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, alphatimefacfac*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, alphatimefacfac*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, alphatimefacfac*normal(2)*ju_n);
#endif

//#ifdef COMBUST_NITSCHE_CONVSTAB
if(nitsche_convstab)
{
    //--------------------     |                                     |
    // normal penalty term   + |  \alpha_c *|u*n| || v ||, || Du ||  |
    //--------------------     |                                     |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, convstabtimefacfac, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, convstabtimefacfac, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, convstabtimefacfac, shp_jump.d0);

    //    |                                   |
    //  - |  1/(2h)*|u*n| || v ||, || u_i ||  |
    //    |                                   |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -convstabtimefacfac*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -convstabtimefacfac*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -convstabtimefacfac*veljump(2,0));

    //------------------------    |                             |
    // normal penalty term RHS  + |  1/(2h)*|u*n| || v ||, j_u  |
    //------------------------    |                             |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, convstabtimefacfac*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, convstabtimefacfac*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, convstabtimefacfac*ju(2));
}
//#endif

//#ifdef COMBUST_NITSCHE_CONVPENALTY
if(nitsche_convpenalty)
{
    dserror("Either you know what this term does or you should not use it!");

    const double pentimefacfac   = 1.0/hk_ele * timefacfac; // note: in Florian's Version, pen_fac was set to 1/h
    //--------------------     |                             |
    // normal penalty term   + |  1/h || v ||*n, || Du ||*n  |
    //--------------------     |                             |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, normal_smooth(0)*pentimefacfac*normal_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_jump.d0, normal_smooth(0)*pentimefacfac*normal_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_jump.d0, normal_smooth(0)*pentimefacfac*normal_smooth(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_jump.d0, normal_smooth(1)*pentimefacfac*normal_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, normal_smooth(1)*pentimefacfac*normal_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_jump.d0, normal_smooth(1)*pentimefacfac*normal_smooth(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_jump.d0, normal_smooth(2)*pentimefacfac*normal_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_jump.d0, normal_smooth(2)*pentimefacfac*normal_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, normal_smooth(2)*pentimefacfac*normal_smooth(2), shp_jump.d0);

    //    |                              |
    //  - |  1/h || v ||*n, || u_i ||*n  |
    //    |                              |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -normal_smooth(0)*pentimefacfac*normal_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -normal_smooth(0)*pentimefacfac*normal_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -normal_smooth(0)*pentimefacfac*normal_smooth(2)*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -normal_smooth(1)*pentimefacfac*normal_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -normal_smooth(1)*pentimefacfac*normal_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -normal_smooth(1)*pentimefacfac*normal_smooth(2)*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -normal_smooth(2)*pentimefacfac*normal_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -normal_smooth(2)*pentimefacfac*normal_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -normal_smooth(2)*pentimefacfac*normal_smooth(2)*veljump(2,0));

    //------------------------    |                        |
    // normal penalty term RHS  + |  1/h || v ||*n, j_u*n  |
    //------------------------    |                        |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, normal_smooth(0)*pentimefacfac*normal_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, normal_smooth(0)*pentimefacfac*normal_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, normal_smooth(0)*pentimefacfac*normal_smooth(2)*ju(2));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, normal_smooth(1)*pentimefacfac*normal_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, normal_smooth(1)*pentimefacfac*normal_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, normal_smooth(1)*pentimefacfac*normal_smooth(2)*ju(2));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, normal_smooth(2)*pentimefacfac*normal_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, normal_smooth(2)*pentimefacfac*normal_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, normal_smooth(2)*pentimefacfac*normal_smooth(2)*ju(2));

    //--------------------------    |                                 |
    // tangential penalty term 1  + |  1/h || v ||*t_1, || Du ||*t_1  |
    //--------------------------    |                                 |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, tang1_smooth(0)*pentimefacfac*tang1_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_jump.d0, tang1_smooth(0)*pentimefacfac*tang1_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_jump.d0, tang1_smooth(0)*pentimefacfac*tang1_smooth(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_jump.d0, tang1_smooth(1)*pentimefacfac*tang1_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, tang1_smooth(1)*pentimefacfac*tang1_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_jump.d0, tang1_smooth(1)*pentimefacfac*tang1_smooth(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_jump.d0, tang1_smooth(2)*pentimefacfac*tang1_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_jump.d0, tang1_smooth(2)*pentimefacfac*tang1_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, tang1_smooth(2)*pentimefacfac*tang1_smooth(2), shp_jump.d0);


    //    |                                  |
    //  - |  1/h || v ||*t_1, || u_i ||*t_1  |
    //    |                                  |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -tang1_smooth(0)*pentimefacfac*tang1_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -tang1_smooth(0)*pentimefacfac*tang1_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -tang1_smooth(0)*pentimefacfac*tang1_smooth(2)*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -tang1_smooth(1)*pentimefacfac*tang1_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -tang1_smooth(1)*pentimefacfac*tang1_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -tang1_smooth(1)*pentimefacfac*tang1_smooth(2)*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -tang1_smooth(2)*pentimefacfac*tang1_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -tang1_smooth(2)*pentimefacfac*tang1_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -tang1_smooth(2)*pentimefacfac*tang1_smooth(2)*veljump(2,0));

    //------------------------------    |                            |
    // tangential penalty term 1 RHS  + |  1/h || v ||*t_1, j_u*t_1  |
    //------------------------------    |                            |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, tang1_smooth(0)*pentimefacfac*tang1_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, tang1_smooth(0)*pentimefacfac*tang1_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, tang1_smooth(0)*pentimefacfac*tang1_smooth(2)*ju(2));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, tang1_smooth(1)*pentimefacfac*tang1_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, tang1_smooth(1)*pentimefacfac*tang1_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, tang1_smooth(1)*pentimefacfac*tang1_smooth(2)*ju(2));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, tang1_smooth(2)*pentimefacfac*tang1_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, tang1_smooth(2)*pentimefacfac*tang1_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, tang1_smooth(2)*pentimefacfac*tang1_smooth(2)*ju(2));

    //--------------------------    |                                 |
    // tangential penalty term 2  + |  1/h || v ||*t_2, || Du ||*t_2  |
    //--------------------------    |                                 |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, tang2_smooth(0)*pentimefacfac*tang2_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_jump.d0, tang2_smooth(0)*pentimefacfac*tang2_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_jump.d0, tang2_smooth(0)*pentimefacfac*tang2_smooth(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_jump.d0, tang2_smooth(1)*pentimefacfac*tang2_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, tang2_smooth(1)*pentimefacfac*tang2_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_jump.d0, tang2_smooth(1)*pentimefacfac*tang2_smooth(2), shp_jump.d0);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_jump.d0, tang2_smooth(2)*pentimefacfac*tang2_smooth(0), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_jump.d0, tang2_smooth(2)*pentimefacfac*tang2_smooth(1), shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, tang2_smooth(2)*pentimefacfac*tang2_smooth(2), shp_jump.d0);


    //    |                                  |
    //  - |  1/h || v ||*t_2, || u_i ||*t_2  |
    //    |                                  |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -tang2_smooth(0)*pentimefacfac*tang2_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -tang2_smooth(0)*pentimefacfac*tang2_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -tang2_smooth(0)*pentimefacfac*tang2_smooth(2)*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -tang2_smooth(1)*pentimefacfac*tang2_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -tang2_smooth(1)*pentimefacfac*tang2_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -tang2_smooth(1)*pentimefacfac*tang2_smooth(2)*veljump(2,0));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -tang2_smooth(2)*pentimefacfac*tang2_smooth(0)*veljump(0,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -tang2_smooth(2)*pentimefacfac*tang2_smooth(1)*veljump(1,0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -tang2_smooth(2)*pentimefacfac*tang2_smooth(2)*veljump(2,0));

    //------------------------------    |                            |
    // tangential penalty term 2 RHS  + |  1/h || v ||*t_2, j_u*t_2  |
    //------------------------------    |                            |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, tang2_smooth(0)*pentimefacfac*tang2_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, tang2_smooth(0)*pentimefacfac*tang2_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, tang2_smooth(0)*pentimefacfac*tang2_smooth(2)*ju(2));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, tang2_smooth(1)*pentimefacfac*tang2_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, tang2_smooth(1)*pentimefacfac*tang2_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, tang2_smooth(1)*pentimefacfac*tang2_smooth(2)*ju(2));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, tang2_smooth(2)*pentimefacfac*tang2_smooth(0)*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, tang2_smooth(2)*pentimefacfac*tang2_smooth(1)*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, tang2_smooth(2)*pentimefacfac*tang2_smooth(2)*ju(2));
}
//#endif


//#ifdef COMBUST_NITSCHE_CONVFLUX
if(nitsche_convflux)
{
#if 1 // this is Florian's version
    //---------------------------------    |                                  |
    // convective flux inflow term plus  + |  v^{+} (u^{+}*n)_{in}, || Du ||  |
    //---------------------------------    |                                  |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfacplus*velproj_plus, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfacplus*velproj_plus, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfacplus*velproj_plus, shp_jump.d0);

    //   |                                   |
    // - |  v^{+} (u^{+}*n)_{in}, || u_i ||  |
    //   |                                   |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, -halftimefacfacplus*velproj_plus*veljump(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, -halftimefacfacplus*velproj_plus*veljump(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, -halftimefacfacplus*velproj_plus*veljump(2));

    //   |                             |
    // + |  v^{+} (u^{+}*n)_{in}, j_u  |
    //   |                             |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfacplus*velproj_plus*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfacplus*velproj_plus*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfacplus*velproj_plus*ju(2));

    //----------------------------------    |                                  |
    // convective flux inflow term minus  + |  v^{-} (u^{-}*n)_{in}, || Du ||  |
    //----------------------------------    |                                  |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfacminus*velproj_minus, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfacminus*velproj_minus, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfacminus*velproj_minus, shp_jump.d0);

    //   |                                   |
    // - |  v^{-} (u^{-}*n)_{in}, || u_i ||  |
    //   |                                   |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, -halftimefacfacminus*velproj_minus*veljump(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, -halftimefacfacminus*velproj_minus*veljump(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, -halftimefacfacminus*velproj_minus*veljump(2));

    //   |                             |
    // + |  v^{-} (u^{-}*n)_{in}, j_u  |
    //   |                             |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfacminus*velproj_minus*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfacminus*velproj_minus*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfacminus*velproj_minus*ju(2));


//    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfac*velproj_plus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfac*velproj_plus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfac*velproj_plus, shp_jump.d0);
//
//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, -halftimefacfac*velproj_plus*veljump(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, -halftimefacfac*velproj_plus*veljump(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, -halftimefacfac*velproj_plus*veljump(2));
//
//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfac*velproj_plus*ju(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfac*velproj_plus*ju(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfac*velproj_plus*ju(2));
//
//    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfac*velproj_minus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfac*velproj_minus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfac*velproj_minus, shp_jump.d0);
//
//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, -halftimefacfac*velproj_minus*veljump(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, -halftimefacfac*velproj_minus*veljump(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, -halftimefacfac*velproj_minus*veljump(2));
//
//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_inv.d0, halftimefacfac*velproj_minus*ju(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_inv.d0, halftimefacfac*velproj_minus*ju(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_inv.d0, halftimefacfac*velproj_minus*ju(2));


//    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, -halftimefacfac*velproj_plus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, -halftimefacfac*velproj_plus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, -halftimefacfac*velproj_plus, shp_jump.d0);

//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, halftimefacfac*velproj_plus*veljump(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, halftimefacfac*velproj_plus*veljump(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, halftimefacfac*velproj_plus*veljump(2));

//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -halftimefacfac*velproj_plus*ju(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -halftimefacfac*velproj_plus*ju(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -halftimefacfac*velproj_plus*ju(2));

//    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, halftimefacfac*velproj_minus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, halftimefacfac*velproj_minus, shp_jump.d0);
//    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, halftimefacfac*velproj_minus, shp_jump.d0);

//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -halftimefacfac*velproj_minus*veljump(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -halftimefacfac*velproj_minus*veljump(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -halftimefacfac*velproj_minus*veljump(2));

//    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, halftimefacfac*velproj_minus*ju(0));
//    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, halftimefacfac*velproj_minus*ju(1));
//    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, halftimefacfac*velproj_minus*ju(2));
#endif
}

#if 0
  //COMBUST_NITSCHE_CONTINUITY
if(nitsche_mass)
{
  //--------------------     |                                                                           |
  // normal penalty term   + | (B_1/h*{\mu}_m + B_2{\rho* norm2(u)}_m) || v || * n_int, || Du || * n_int |
  //--------------------     |                                                                           |

  assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, massstabtimefacfac*normal(0)*normal(0), shp_jump.d0);
  assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_jump.d0, massstabtimefacfac*normal(0)*normal(1), shp_jump.d0);
  assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_jump.d0, massstabtimefacfac*normal(0)*normal(2), shp_jump.d0);

  assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_jump.d0, massstabtimefacfac*normal(1)*normal(0), shp_jump.d0);
  assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, massstabtimefacfac*normal(1)*normal(1), shp_jump.d0);
  assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_jump.d0, massstabtimefacfac*normal(1)*normal(2), shp_jump.d0);

  assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_jump.d0, massstabtimefacfac*normal(2)*normal(0), shp_jump.d0);
  assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_jump.d0, massstabtimefacfac*normal(2)*normal(1), shp_jump.d0);
  assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, massstabtimefacfac*normal(2)*normal(2), shp_jump.d0);

  //   |                                                                           |
  // - |  (B_1/h*{\mu}_m + B_2{\rho* norm2(u)}_m) || v || * n_int, || u || * n_int |
  //   |                                                                           |

  assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, - massstabtimefacfac*normal(0)*(veljump(0,0)*normal(0)
                                                                                             +veljump(1,0)*normal(1)
                                                                                             +veljump(2,0)*normal(2)));
  assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, - massstabtimefacfac*normal(1)*(veljump(0,0)*normal(0)
                                                                                             +veljump(1,0)*normal(1)
                                                                                             +veljump(2,0)*normal(2)));
  assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, - massstabtimefacfac*normal(2)*(veljump(0,0)*normal(0)
                                                                                             +veljump(1,0)*normal(1)
                                                                                             +veljump(2,0)*normal(2)));

  //   |                                                                       |
  // + |  (B_1/h*{\mu}_m + B_2{\rho* norm2(u)}_m) || v || * n_int, j_u * n_int |
  //   |                                                                       |

  assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0,   massstabtimefacfac*normal(0)*(ju(0,0)*normal(0)
                                                                                             +ju(1,0)*normal(1)
                                                                                             +ju(2,0)*normal(2)));
  assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0,   massstabtimefacfac*normal(1)*(ju(0,0)*normal(0)
                                                                                             +ju(1,0)*normal(1)
                                                                                             +ju(2,0)*normal(2)));
  assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0,   massstabtimefacfac*normal(2)*(ju(0,0)*normal(0)
                                                                                             +ju(1,0)*normal(1)
                                                                                             +ju(2,0)*normal(2)));
}
#endif

//#endif
} // loop Gaussian points

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class M3, class V1, class V2>
void SysmatBoundaryNitsche(
    const DRT::ELEMENTS::Combust3*             ele,            ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*&  ih,  ///< connection to the interface handler
    const XFEM::ElementDofManager&             dofman,         ///< dofmanager of the current element
    const M1&                                  evelnp,         ///<
    const V1&                                  eprenp,         ///<
    const V2&                                  ephi,           ///<
    const M2&                                  egradphi,
    const M3&                                  egradphi2,
    const V2&                                  ecurv,
    Teuchos::RCP<const MAT::Material>          material,       ///< fluid material
    const double                               time,         ///< current time
    const INPAR::FLUID::TimeIntegrationScheme  timealgo,       ///< time discretization type
    const double&                              dt,             ///< delta t (time step size)
    const double&                              theta,          ///< factor for one step theta scheme
    const double                               ga_alphaF,
    const double                               ga_alphaM,
    const double                               ga_gamma,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler,      ///<
    const double                               flamespeed,     ///<
    const double                               marksteinlength,///<
    const double                               nitschevel,     ///<
    const double                               nitschepres,    ///<
    const double                               ele_meas_plus,  ///< measure of element in plus domain
    const double                               ele_meas_minus, ///< measure of element in minus domain
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               variablesurftens,
    const bool                                 second_deriv,
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType     veljumptype,
    const INPAR::COMBUST::FluxJumpType         fluxjumptype,
    const bool                                 smoothed_boundary_integration,
    const INPAR::COMBUST::WeightType           weighttype,
    const bool                                 nitsche_convflux,
    const bool                                 nitsche_convstab,
    const bool                                 nitsche_convpenalty,
    const bool                                 nitsche_mass
    )
{
  if (timealgo == INPAR::FLUID::timeint_afgenalpha)
    dserror("SysmatBoundaryNitsche: Genalpha does not function properly in this method.");

  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - boundary");

#ifdef DEBUG
  if (ASSTYPE != XFEM::xfem_assembly)
    dserror("boundary integration integrals must only be added to intersected XFEM elements");
#endif

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // flag for instationary problems
  bool instationary = false;
  if (timealgo != INPAR::FLUID::timeint_stationary)
    instationary = true;

  // time integration constant
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  //---------------------------------------------------------------------------
  // get material parameters for all boundary integration cells of this element
  //---------------------------------------------------------------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif

  //---------------------------------------------------------------------------
  // get material parameters for all boundary integration cells of this element
  //---------------------------------------------------------------------------
  // density \rho
  double dens_plus = -1.0;
  // dynamic viscosity \mu
  double dynvisc_plus = -1.0;
  // density \rho
  double dens_minus = -1.0;
  // dynamic viscosity \mu
  double dynvisc_minus = -1.0;
  // set density and dynamic viscosity
  COMBUST::GetMaterialParams(material, dens_plus, dynvisc_plus, dens_minus, dynvisc_minus);

  // to consistently take into account the blend material option here, a small hack (due to assumption) is best added
  // I commented this part since it assumes that the minus domain is the blending domain.
  // Moreover, nobody will use this option in the old framework anymore.
  // Hence, this should not cause any problems.
#if 0
  double dynvisc = -1.0;
  double dens = -1.0;
  // blend smoothly to minus (gas bubble) domain (This is assumed here!!!!!)
  COMBUST::BlendMaterial(ele, false, time, dens, dynvisc, dens_plus, dynvisc_plus, dens_minus, dynvisc_minus);
  dens_minus = dens;
  dynvisc_minus = dynvisc;
#endif

  double wallfac = 1.0;
#ifdef ORACLES
  //--------------------------------------------------------
  // get physical coordinates of element center
  // remark: we want to blend the flame speed close to walls
  //--------------------------------------------------------
  LINALG::Matrix<3,1> center(true);
  // the element center is at 0,0,0 (we assume hex elements, no tets!)
  LINALG::Matrix<numnode,1> funct(true);
  DRT::UTILS::shape_function_3D(funct,0.0,0.0,0.0,DISTYPE);
  // compute coordinates of element center
  center.MultiplyNN(xyze,funct);

  //    wall
  // 1.0 |     ______
  //     |    /
  // 0.0 |___/ ,
  //     |     H/6
  const double wallzone = 0.0299/6.0;
  if (center(0) > 0.0) // inside combustion chamber
  {
    if ( (0.0653-abs(center(1))) < wallzone or // close to top or bottom wall
                       center(0) < wallzone )  // close to step
    {
      // wall factor is 0 at the wall and 1 at H/6 or further away from the wall
      wallfac = 6.0/0.0299 * std::min(0.0653-abs(center(1)),center(0));
      if (wallfac < 0.1) // cut off the last 10% to guarantee a zero jump at the wall
        wallfac = 0.0;
    }
  }
#endif

  //-----------------------------------------------------------------------------------
  // get surface tension coefficient for all boundary integration cells of this element
  //-----------------------------------------------------------------------------------
  double surftenscoeff = 0.0;
  // set surface tension
  COMBUST::GetMaterialParams(material, surftenscoeff);

  //-----------------------------------------------------------------------------------------
  // get weights for average operators { () } = kappa_plus * ()_plus + kappa_minus * ()_minus
  //                               and < () > = kappa_minus* ()_plus + kappa_plus  * ()_minus
  //-----------------------------------------------------------------------------------------

  const double ele_meas = ele_meas_plus + ele_meas_minus;
  if(ele_meas == 0.0) dserror("element measure is zero!!! Weighted average operators can not get defined!");
  //double kappa_plus  = ele_meas_plus/ele_meas;
  //double kappa_minus = ele_meas_minus/ele_meas;
  const double a_plus = ele_meas_plus/ele_meas;
  const double a_minus = ele_meas_minus/ele_meas;

  double kappa_plus = 0.0;
  double kappa_minus = 0.0;

  if (weighttype == INPAR::COMBUST::weight_intersection_visc_based_harmonic)
  {
    kappa_plus = dynvisc_minus*a_plus / (dynvisc_plus*a_minus + dynvisc_minus*a_plus);
    kappa_minus = dynvisc_plus*a_minus / (dynvisc_plus*a_minus + dynvisc_minus*a_plus);
  }
  else if (weighttype == INPAR::COMBUST::weight_visc_based_harmonic)
  {
    kappa_plus = dynvisc_minus / (dynvisc_plus + dynvisc_minus);
    kappa_minus = dynvisc_plus / (dynvisc_plus + dynvisc_minus);
  }
  else dserror("Weighting type not supported!");

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  //----------------------------------------
  // compute Nitsche stabilization patameter
  //----------------------------------------
  // calculate element diameter
  const double hk_eleDiam = COMBUST::getEleDiameter<DISTYPE>(xyze);

  double alphau = 0.0;

  if (weighttype == INPAR::COMBUST::weight_visc_based_harmonic)
  {
    alphau = nitschevel/hk_eleDiam * dynvisc_minus * dynvisc_plus / (dynvisc_plus + dynvisc_minus);
  }
  else if (weighttype == INPAR::COMBUST::weight_intersection_visc_based_harmonic)
  {
    double area = 0.0;
    // get domain integration cells for this element
    const GEO::BoundaryIntCells& elementBoundaryIntCells = ih->ElementBoundaryIntCells(ele->Id());
    GEO::BoundaryIntCells::const_iterator itercell;
    for(itercell = elementBoundaryIntCells.begin(); itercell != elementBoundaryIntCells.end(); ++itercell )
    {
      if (!(itercell->Shape() == DRT::Element::tri3 or itercell->Shape() == DRT::Element::quad4))
        dserror("invalid type of boundary integration cell for surface area calculation");

      // get coordinates of vertices defining flame front patch
      const LINALG::SerialDenseMatrix& coords = itercell->CellNodalPosXYZ();

      // first point of flame front patch
      LINALG::Matrix<3,1> pointA;
      pointA(0) = coords(0,0);
      pointA(1) = coords(1,0);
      pointA(2) = coords(2,0);

      // second point of flame front patch
      LINALG::Matrix<3,1> pointB;
      pointB(0) = coords(0,1);
      pointB(1) = coords(1,1);
      pointB(2) = coords(2,1);

      // first edge of flame front patch
      LINALG::Matrix<3,1> edgeBA;
      edgeBA.Update(1.0, pointA, -1.0, pointB);

      // third point of flame front patch
      LINALG::Matrix<3,1> pointC;
      pointC(0) = coords(0,2);
      pointC(1) = coords(1,2);
      pointC(2) = coords(2,2);

      // second edge of flame front patch
      LINALG::Matrix<3,1> edgeBC;
      edgeBC.Update(1.0, pointC, -1.0, pointB);

      LINALG::Matrix<3,1> crossP;
      crossP(0) = edgeBA(1)*edgeBC(2) - edgeBA(2)*edgeBC(1);
      crossP(1) = edgeBA(2)*edgeBC(0) - edgeBA(0)*edgeBC(2);
      crossP(2) = edgeBA(0)*edgeBC(1) - edgeBA(1)*edgeBC(0);

      if (itercell->Shape() == DRT::Element::tri3)
        area += crossP.Norm2() / 2.0;
      else
        area += crossP.Norm2();
     }

    const double gamma = hk_eleDiam * area / ele_meas;
    alphau = nitschevel/hk_eleDiam * ((gamma * dynvisc_minus * dynvisc_plus)/(dynvisc_plus*a_minus + dynvisc_minus*a_plus));
  }

#if 0 //Florian
  // get maximum of the viscosities of all fluid materials
  const double dynvisc_max = std::max(dynvisc_plus,dynvisc_minus);

  // Nitsche parameter velocity ( alph_u := 1/h_k * mu_max * gamma_u ) see Becker, Hansbo 2009
  alphau = 1.0/hk_eleDiam * dynvisc_max * nitschevel;
  //const double alphau = 1.0/hk_eleDiam * 4.0*nitschevel*(dynvisc_plus*dynvisc_minus);
#endif

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->ElementBoundaryIntCells(ele->Id());

  //-------------------------------------------------------------------------------
  // evaluate the enrichment function at the interface (boundary integration cells)
  //-------------------------------------------------------------------------------
#ifndef COMBUST_SXFEM
  const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
  const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);
#endif

  //-------------------------------------
  // loop over boundary integration cells
  //-------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    switch (cell->Shape())
    {

    case DRT::Element::tri3:
      COMBUST::Nitsche_SysmatBoundaryCell<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
          *cell, DRT::UTILS::intrule_tri_37point, xyze, ele, higher_order_ele,
          evelnp, eprenp, egradphi, egradphi2, ecurv,
          flamespeed, marksteinlength,
          timefac, instationary, assembler,
#ifndef COMBUST_SXFEM
          enrvals_plus, enrvals_minus,
#else
          dofman, ephi,
#endif
          kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
          alphau, wallfac, hk_eleDiam, nitschevel,
          numparamvelx, numparampres, surftensapprox, variablesurftens, second_deriv, surftenscoeff, connected_interface,
          veljumptype, fluxjumptype,
          smoothed_boundary_integration,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass
      );
      break;
    case DRT::Element::quad4:
      COMBUST::Nitsche_SysmatBoundaryCell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          *cell, DRT::UTILS::intrule_quad_25point , xyze, ele, higher_order_ele,
          evelnp, eprenp, egradphi, egradphi2, ecurv,
          flamespeed, marksteinlength,
          timefac, instationary, assembler,
#ifndef COMBUST_SXFEM
          enrvals_plus, enrvals_minus,
#else
          dofman, ephi,
#endif
          kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus, dens_plus, dens_minus,
          alphau, wallfac, hk_eleDiam, nitschevel,
          numparamvelx, numparampres, surftensapprox, variablesurftens, second_deriv, surftenscoeff, connected_interface,
          veljumptype, fluxjumptype,
          smoothed_boundary_integration,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass
      );
      break;
    default:
      dserror("invalid type of boundary integration cell");
    }
  } // loop boundary integration cells

  return;
}

/*!
  Calculate boundary integrals for Neumann inflow term
 */
template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType CELLDISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          class M1, class M2, class V1>
void Nitsche_SysmatNeumannInflow(
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const DRT::ELEMENTS::Combust3*         ele,           ///< the element those matrix is calculated
    const DRT::ELEMENTS::Combust3Surface*  elesurf,
    const XFEM::ElementDofManager&         dofman,
    GEO::BoundaryIntCell                   cell,
    const DRT::UTILS::GaussRule2D&         intrule2D,
    const M2&                              xyze,
    const M1&                              evelnp,
    const V1&                              ephi,
    const size_t&                          numparamvelx,
    const bool&                            newton,
    Teuchos::RCP<const MAT::Material>      material,
    const double&                          time,
    const double&                          timefac
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  //--------------------------------------------------
  // get material parameters for this integration cell
  //--------------------------------------------------
  // density \rho^{n+1}
  double densaf  = -1.0;
  // dynamic viscosity \mu^{n+1}
  double dynvisc = -1.0;
  // set density and dynamic viscosity
  {
    double densplus = -1.0;
    double densminus = -1.0;
    double dynviscplus = -1.0;
    double dynviscminus = -1.0;
    COMBUST::GetMaterialParams(material, true, densplus, dynviscplus);
    COMBUST::GetMaterialParams(material, false, densminus, dynviscminus);

    // blend smoothly to minus (gas bubble) domain
    COMBUST::BlendMaterial(ele, cell.getDomainPlus(), time, densaf, dynvisc, densplus, dynviscplus, densminus, dynviscminus);
  }

  //-----------------------------------------------------------
  // evaluate the enrichment function for this integration cell
  //-----------------------------------------------------------
  const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,cell,ephi);

  //----------------------------------------------
  // get Gaussian points for this integration cell
  //----------------------------------------------
  const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

  // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
  LINALG::SerialDenseMatrix cellXiDomaintmp = cell.CellNodalPosXiDomain();
  // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
  LINALG::SerialDenseMatrix cellXYZDomaintmp = cell.CellNodalPosXYZ();
  // transform to fixed size format
  const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
  const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

  //--------------------------------------
  // integration loop over Gaussian points
  //--------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // new transformation for boundary integrals
    // 1. define a coupled transformation x_3D(xi_3D(eta_2D)): transformation from 2D->3D
    // 2. compute the corresponding Jacobian J_eta2D->x_3D and
    // 3. the corresponding surface integral factor sqrt(det(J_eta2D->x_3D^T * J_eta2D->x_3D))
    // 4. approximate integral with Gauss rule in eta coordinate system
    // 5. evaluate the transformed integrand f(x(xi(eta)))

    const LINALG::Matrix<2,1> gpinEta2D(intpoints.qxg[iquad]);

    // jacobian for coupled transformation
    // get derivatives dxi_3D/deta_2D
    static LINALG::Matrix<2,numvertices> deriv_eta2D;
    DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),CELLDISTYPE);

    // calculate dxi3Ddeta2D
    static LINALG::Matrix<3,2> dXi3Ddeta2D;
    dXi3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < (int)numvertices; k++)
          dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);

    // transform Gauss point to xi3D space (element parameter space)
    static LINALG::Matrix<nsd,1> gpinXi3D;
    gpinXi3D.Clear();
    // coordinates of this integration point in element coordinates \xi^domain
    GEO::mapEtaBToXiD(cell, gpinEta2D, gpinXi3D);

    static LINALG::Matrix<3,numnode> deriv_xi3D;
    DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), DISTYPE);

    // calculate dx3Ddxi3D
    static LINALG::Matrix<3,3> dX3DdXi3D;
    dX3DdXi3D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 3; j++) // derivatives
        for (int k = 0; k < (int)numnode; k++)
          dX3DdXi3D(i,j) += xyze(i,k)*deriv_xi3D(j,k);

    // get the coupled Jacobian dx3Ddeta2D
    static LINALG::Matrix<3,2> dx3Ddeta2D;
    dx3Ddeta2D.Clear();
    for (int i = 0; i < 3; i++)   // dimensions
      for (int j = 0; j < 2; j++) // derivatives
        for (int k = 0; k < 3; k++)
          dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);

    // get deformation factor
    static LINALG::Matrix<2,2> Jac_tmp; // J^T*J
    Jac_tmp.Clear();
    Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

    if(Jac_tmp.Determinant() <= 0.0) dserror("deformation factor for boundary integration is zero");

    const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

    const double fac = intpoints.qwgt[iquad]*deform_factor;

    LINALG::Matrix<2,1> posEtaBoundary;
    posEtaBoundary.Clear();
    for (int i= 0; i< 2; i++){
      posEtaBoundary(i,0) = gpinEta2D(i,0);
    }

    LINALG::Matrix<nsd,1> posXiDomain;
    posXiDomain.Clear();
    for (int i= 0; i< 3; i++){
      posXiDomain(i,0) = gpinXi3D(i,0);
    }

    //--------------------------------------------------------------------------------------------
    // compute normal vector (normed, in physical coordinates)
    // remark: for linear boundary integration cells this could be done before the loop over all
    //         Gaussian points
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    normal.Clear();
    // outside pointing normal vector
    GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);

    //----------------------------------------------------------------------------
    // evaluate shape functions and their first derivatives at this Gaussian point
    //----------------------------------------------------------------------------
    // shape functions for fields with dofs at nodes
    static LINALG::Matrix<numnode,1> funct;
    funct.Clear();
    DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
    // derivatives of shape functions for fields with dofs at nodes
    static LINALG::Matrix<nsd,numnode> deriv;
    deriv.Clear();
    DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    //--------------------------------------------------------
    // procedures involving Jacobian matrix for domain mapping
    //--------------------------------------------------------
    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<nsd,nsd> xjm;
    xjm.Clear();
    xjm.MultiplyNT(deriv,xyze);
    // determinant for mapping from physical space (X^3D/domain) to element (Xi^3D/domain)
    const double detXtoXi = xjm.Determinant();
    //#ifdef DEBUG
    if (detXtoXi < 0.0)
    {
      dserror("negative Jacobian determinant detXtoXi: %f for element %d!", detXtoXi, ele->Id());
    }
    //#endif
    // inverse of jacobian
    static LINALG::Matrix<nsd,nsd> xji;
    xji.Clear();
    xji.Invert(xjm);

    //----------------------------------------------------------------
    // rearrange (enriched) shape functions as approximation functions
    //----------------------------------------------------------------
    const size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<0,shpVecSize> shp;

    if (ASSTYPE == XFEM::xfem_assembly)
    {
      // temporary arrays holding enriched shape functions (N * \Psi)
      static LINALG::Matrix<shpVecSize,1> enr_funct;
      enr_funct.Clear();

      // shape functions and derivatives for nodal parameters (dofs)
      enrvals.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, enr_funct);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        shp.d0(iparam) = enr_funct(iparam);
      }
    }
    else // not xfem_assembly i.e. standard assembly
    {
      // fill approximation functions for standard FEM
      // remark: numparamvelx == numnode, for standard FEM
      for (size_t iparam = 0; iparam < numnode; ++iparam)
      {
        shp.d0(iparam) = funct(iparam);
      }
    }
    //------------------------------------
    // get velocities at integration point
    //------------------------------------
    const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
    //const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
    //const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);

    // compute velocity vector and normal velocity at integration point
    const double inflowvelnorm = gpvelnp.Dot(normal);

    // check normal velocity -> further computation only required for
    // negative normal velocity, that is, inflow at this Neumann boundary
    if (inflowvelnorm < -1.0e-8)
    {
#if 1 //def COMBUST_NEUMANN_INFLOW

      //---------------------------------------------------    |                 |
      // Neumann inflow term
      //---------------------------------------------------    |                 |
      // compute total (time and spatial) integration factor
      const double denstimefacfac = densaf * timefac * fac;

      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp.d0, -denstimefacfac*inflowvelnorm, shp.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp.d0, -denstimefacfac*inflowvelnorm, shp.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp.d0, -denstimefacfac*inflowvelnorm, shp.d0);

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp.d0, denstimefacfac*inflowvelnorm*gpvelnp(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp.d0, denstimefacfac*inflowvelnorm*gpvelnp(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp.d0, denstimefacfac*inflowvelnorm*gpvelnp(2));

      if(newton)
      {
        assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp.d0, -denstimefacfac*gpvelnp(0)*normal(0), shp.d0);
        assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp.d0, -denstimefacfac*gpvelnp(0)*normal(1), shp.d0);
        assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp.d0, -denstimefacfac*gpvelnp(0)*normal(2), shp.d0);

        assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp.d0, -denstimefacfac*gpvelnp(1)*normal(0), shp.d0);
        assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp.d0, -denstimefacfac*gpvelnp(1)*normal(1), shp.d0);
        assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp.d0, -denstimefacfac*gpvelnp(1)*normal(2), shp.d0);

        assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp.d0, -denstimefacfac*gpvelnp(2)*normal(0), shp.d0);
        assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp.d0, -denstimefacfac*gpvelnp(2)*normal(1), shp.d0);
        assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp.d0, -denstimefacfac*gpvelnp(2)*normal(2), shp.d0);
      }
#endif
    }

  } // loop Gaussian points
  return;
}

} // namespace COMBUST


#endif

