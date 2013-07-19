/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_premixed_nitsche.H

\brief system matrix formulation with Nitsche's method
       for premixed combustion problem / 3D combustion element

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_PREMIXED_NITSCHE_NORMAL_H
#define COMBUST3_SYSMAT_PREMIXED_NITSCHE_NORMAL_H


#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat_standard.H"
#include "combust3_sysmat_surfacetension.H"
#include "combust_defines.H"
#include "combust3_stabilization.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_normals.H"


namespace COMBUST
{


/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2>
void SysmatDomainNitscheNormal(
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
    const M2&                      egradphi,     ///<
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

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  const size_t numparamveln = XFEM::NumParam<0,ASSTYPE>::get(dofman, XFEM::PHYSICS::Veln);

  // stabilization parameter
  const double hk = COMBUST::UTILS::HK<DISTYPE>(evelnp,xyze);
  const double mk = COMBUST::UTILS::MK<DISTYPE>();

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  //-----------------------------------
  // loop over domain integration cells
  //-----------------------------------
  int icell = 0;
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
    COMBUST::GetMaterialParams(material, cell->getDomainPlus(), densaf, dynvisc);
    // density \rho^{n}
    // TODO until the XFEM time integration is done properly, densities from previous and
    // current time steps are identical
    const double densn = densaf;

    //-----------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //-----------------------------------------------------------
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);

    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Bisected(), cell->Shape());
    //const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_4point;
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
        cout << "detcell :  " << detcell << endl;
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

      if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

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
        DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const size_t shpVecSizeVel    = COMBUST::SizeFacVel<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const size_t shpVecSizePres   = COMBUST::SizeFacPres<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      // second derivatives needed for velocity
      static XFEM::ApproxFuncNormalVector<2,shpVecSizeVel> shp(false);
      shp.Clear();

      // fill approximation functions for standard FEM and XFEM
      for (size_t iparam = 0; iparam < numnode; ++iparam)
      {
        shp.velx.d0.s(iparam) = funct(iparam);
        shp.vely.d0.s(iparam) = funct(iparam);
        shp.velz.d0.s(iparam) = funct(iparam);

        shp.velx.dx.s(iparam) = derxy(0,iparam);
        shp.vely.dx.s(iparam) = derxy(0,iparam);
        shp.velz.dx.s(iparam) = derxy(0,iparam);

        shp.velx.dy.s(iparam) = derxy(1,iparam);
        shp.vely.dy.s(iparam) = derxy(1,iparam);
        shp.velz.dy.s(iparam) = derxy(1,iparam);

        shp.velx.dz.s(iparam) = derxy(2,iparam);
        shp.vely.dz.s(iparam) = derxy(2,iparam);
        shp.velz.dz.s(iparam) = derxy(2,iparam);

        shp.velx.dxdx.s(iparam) = derxy2(0,iparam);
        shp.vely.dxdx.s(iparam) = derxy2(0,iparam);
        shp.velz.dxdx.s(iparam) = derxy2(0,iparam);

        shp.velx.dxdy.s(iparam) = derxy2(3,iparam);
        shp.vely.dxdy.s(iparam) = derxy2(3,iparam);
        shp.velz.dxdy.s(iparam) = derxy2(3,iparam);

        shp.velx.dxdz.s(iparam) = derxy2(4,iparam);
        shp.vely.dxdz.s(iparam) = derxy2(4,iparam);
        shp.velz.dxdz.s(iparam) = derxy2(4,iparam);

        shp.velx.dydx.s(iparam) = shp.velx.dxdy.s(iparam);
        shp.vely.dydx.s(iparam) = shp.vely.dxdy.s(iparam);
        shp.velz.dydx.s(iparam) = shp.velz.dxdy.s(iparam);

        shp.velx.dydy.s(iparam) = derxy2(1,iparam);
        shp.vely.dydy.s(iparam) = derxy2(1,iparam);
        shp.velz.dydy.s(iparam) = derxy2(1,iparam);

        shp.velx.dydz.s(iparam) = derxy2(5,iparam);
        shp.vely.dydz.s(iparam) = derxy2(5,iparam);
        shp.velz.dydz.s(iparam) = derxy2(5,iparam);

        shp.velx.dzdx.s(iparam) = shp.velx.dxdz.s(iparam);
        shp.vely.dzdx.s(iparam) = shp.vely.dxdz.s(iparam);
        shp.velz.dzdx.s(iparam) = shp.velz.dxdz.s(iparam);

        shp.velx.dzdy.s(iparam) = shp.velx.dydz.s(iparam);
        shp.vely.dzdy.s(iparam) = shp.vely.dydz.s(iparam);
        shp.velz.dzdy.s(iparam) = shp.velz.dydz.s(iparam);

        shp.velx.dzdz.s(iparam) = derxy2(2,iparam);
        shp.vely.dzdz.s(iparam) = derxy2(2,iparam);
        shp.velz.dzdz.s(iparam) = derxy2(2,iparam);
      }

      static XFEM::ApproxFunc<1,shpVecSizePres> shp_pres;

      //--------------
      // XFEM assembly
      //--------------
      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSizePres,1> enr_funct_pres;
        enr_funct_pres.Clear();
        static LINALG::Matrix<3,shpVecSizePres> enr_derxy_pres;
        enr_derxy_pres.Clear();
        static LINALG::Matrix<6,shpVecSizePres> enr_derxy2_pres;
        enr_derxy2_pres.Clear();

#ifdef COLLAPSE_FLAME
        LINALG::Matrix<3,1> normal(true);
        for (unsigned inode=0;inode<numnode;inode++)
        {
          normal(0) += funct(inode)*xyze(0,inode);
          normal(1) += funct(inode)*xyze(1,inode);
        }

        LINALG::Matrix<3,1> point(true);
        point(0) = normal(0);
        point(1) = normal(1);

        const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
        if (norm == 0.0) dserror("norm of normal vector is zero!");
        normal.Scale(1.0/norm);
#endif
//      cout << xyze << endl;
//      cout << "normal "<< normal << endl;
//if (ele->Id()==116)
////  or
////    ele->Id()==226 or
////    ele->Id()==326 or
////    ele->Id()==216)
//{
//  cout << "Normalenvektor Domain Sysmat Viertelposition" << normal << endl;
//}
//if (ele->Id()==141)
////  or
////    ele->Id()==308 or
////    ele->Id()==201 or
////    ele->Id()==133)
//{
//  cout << "Normalenvektor Domain Sysmat Achtelposition" << normal << endl;
//}
        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeNormalShapeFunction(
            funct,
            derxy,
            derxy2,
            egradphi,
#ifdef COLLAPSE_FLAME
        normal,
        point,
#endif
            shp);

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Pres,
            funct,
            derxy,
            derxy2,
            enr_funct_pres,
            enr_derxy_pres,
            enr_derxy2_pres);

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparampres; ++iparam)
        {
          shp_pres.d0(iparam) = enr_funct_pres(iparam);
          shp_pres.dx(iparam) = enr_derxy_pres(0,iparam);
          shp_pres.dy(iparam) = enr_derxy_pres(1,iparam);
          shp_pres.dz(iparam) = enr_derxy_pres(2,iparam);
        }
      }
      //------------------
      // standard assembly
      //------------------
      else // not xfem_assembly i.e. standard assembly
      {
#ifdef DEBUG
        // for standard FEM numparamvelx == numnode must hold
        dsassert(numparamvelx == numnode,"numparamvelx = numnode, for standard FEM assembly");
#endif
        for (size_t inode = 0; inode < numnode; ++inode)
        {
          shp_pres.d0(inode) = funct(inode);
          shp_pres.dx(inode) = derxy(0,inode);
          shp_pres.dy(inode) = derxy(1,inode);
          shp_pres.dz(inode) = derxy(2,inode);
        }
      }

      // set measure of element in plus and minus domain
      if(cell->getDomainPlus() == true)
      {
        ele_meas_plus  = fac;
      }
      else
      {
        ele_meas_minus = fac;
      }

//if(ele->Id()==2 and icell == 0 and iquad==0)
//{
//  cout << "ele id " << ele->Id() << endl;
//  cout << posXiDomain << endl;
//if (cell->getDomainPlus()) cout << "plus" << endl;
//else cout << "minus" << endl;
//if (ASSTYPE == XFEM::xfem_assembly) cout << "XFEM assembly" << endl;
//if (ASSTYPE == XFEM::standard_assembly) cout << "standard assembly" << endl;
//shp.Print();
//  cout << "pressure d0 " << shp_pres.d0 << endl;
//  cout << "pressure dx " << shp_pres.dx << endl;
//  cout << "pressure dy " << shp_pres.dy << endl;
//  cout << "pressure dz " << shp_pres.dz << endl;
//
//  const std::vector<int>& velxdof (dofman.LocalDofPosPerField<SymmetryMangle<XFEM::PHYSICS::Velx>::name>());
//  const std::vector<int>& velydof (dofman.LocalDofPosPerField<SymmetryMangle<XFEM::PHYSICS::Vely>::name>());
//  const std::vector<int>& velzdof (dofman.LocalDofPosPerField<SymmetryMangle<XFEM::PHYSICS::Velz>::name>());
//  const std::vector<int>& velndof (dofman.LocalDofPosPerField<SymmetryMangle<XFEM::PHYSICS::Veln>::name>());
//  const std::vector<int>& presdof (dofman.LocalDofPosPerField<SymmetryMangle<XFEM::PHYSICS::Pres>::name>());
//
//  const std::size_t numparamvelxloc = velxdof.size();
//  const std::size_t numparamvelyloc = velydof.size();
//  const std::size_t numparamvelzloc = velzdof.size();
//  const std::size_t numparamvelnloc = velndof.size();
//  const std::size_t numparampresloc = presdof.size();
//
//  cout << "numparamvelxloc " << numparamvelxloc << endl;
//  cout << "numparamvelyloc " << numparamvelyloc << endl;
//  cout << "numparamvelzloc " << numparamvelzloc << endl;
//  cout << "numparamvelnloc " << numparamvelnloc << "numparamveln "<< numparamveln << endl;
//  cout << "numparampresloc " << numparampresloc << endl;
//
//}

      //------------------------------------------------------
      // get velocities and accelerations at integration point
      //------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, shp, numparamvelx, ele, dofman);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPointNormal(eveln , shp, numparamvelx, ele, dofman);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPointNormal(evelnm, shp, numparamvelx, ele, dofman);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPointNormal(eaccn , shp, numparamvelx, ele, dofman);
      // get acceleration at n+alpha_m for generalized alpha time integration scheme
      LINALG::Matrix<3,1> gpaccam(true);
      if (genalpha)
        gpaccam  = XFEM::interpolateVectorFieldToIntPoint(eaccam , shp, numparamvelx);

      //----------------------------------------------
      // get velocity derivatives at integration point
      //----------------------------------------------
      static LINALG::Matrix<3,3> vderxy;
      vderxy.Clear();
      {
        const int* nodeids = ele->NodeIds();
        std::size_t velncounter = 0;

        // vderxy = enr_derxy(j,k)*evelnp(i,k);
        for (size_t inode = 0; inode < numparamvelx; ++inode)
        {
          // standard shape functions are identical for all vector components
          // e.g. shp.velx.dx.s == shp.vely.dx.s == shp.velz.dx.s
          vderxy(0,0) += evelnp(0,inode)*shp.velx.dx.s(inode);
          vderxy(0,1) += evelnp(0,inode)*shp.velx.dy.s(inode);
          vderxy(0,2) += evelnp(0,inode)*shp.velx.dz.s(inode);

          vderxy(1,0) += evelnp(1,inode)*shp.vely.dx.s(inode);
          vderxy(1,1) += evelnp(1,inode)*shp.vely.dy.s(inode);
          vderxy(1,2) += evelnp(1,inode)*shp.vely.dz.s(inode);

          vderxy(2,0) += evelnp(2,inode)*shp.velz.dx.s(inode);
          vderxy(2,1) += evelnp(2,inode)*shp.velz.dy.s(inode);
          vderxy(2,2) += evelnp(2,inode)*shp.velz.dz.s(inode);

          const int gid = nodeids[inode];
          const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);

          for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
              enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
          {
            if (enrfield->getField() == XFEM::PHYSICS::Veln)
            {
              vderxy(0,0) += evelnp(3,velncounter)*shp.velx.dx.n(velncounter);
              vderxy(0,1) += evelnp(3,velncounter)*shp.velx.dy.n(velncounter);
              vderxy(0,2) += evelnp(3,velncounter)*shp.velx.dz.n(velncounter);

              vderxy(1,0) += evelnp(3,velncounter)*shp.vely.dx.n(velncounter);
              vderxy(1,1) += evelnp(3,velncounter)*shp.vely.dy.n(velncounter);
              vderxy(1,2) += evelnp(3,velncounter)*shp.vely.dz.n(velncounter);

              vderxy(2,0) += evelnp(3,velncounter)*shp.velz.dx.n(velncounter);
              vderxy(2,1) += evelnp(3,velncounter)*shp.velz.dy.n(velncounter);
              vderxy(2,2) += evelnp(3,velncounter)*shp.velz.dz.n(velncounter);

              velncounter += 1;
            }
          }
        }
        if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
      }

      //-----------------------------------------------------
      // get second velocity derivatives at integration point
      //-----------------------------------------------------
      static LINALG::Matrix<3,6> vderxy2;
      vderxy2.Clear();
      //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
      if (higher_order_ele)
      {
        const int* nodeids = ele->NodeIds();
        std::size_t velncounter = 0;

        // vderxy = enr_derxy(j,k)*evelnp(i,k);
        for (size_t inode = 0; inode < numparamvelx; ++inode)
        {
          // standard shape functions are identical for all vector components
          // e.g. shp.velx.dxdx.s == shp.vely.dxdx.s == shp.velz.dxdx.s
          vderxy2(0,0) += evelnp(0,inode)*shp.velx.dxdx.s(inode);
          vderxy2(0,1) += evelnp(0,inode)*shp.velx.dydy.s(inode);
          vderxy2(0,2) += evelnp(0,inode)*shp.velx.dzdz.s(inode);
          vderxy2(0,3) += evelnp(0,inode)*shp.velx.dxdy.s(inode);
          vderxy2(0,4) += evelnp(0,inode)*shp.velx.dxdz.s(inode);
          vderxy2(0,5) += evelnp(0,inode)*shp.velx.dydz.s(inode);

          vderxy2(1,0) += evelnp(1,inode)*shp.vely.dxdx.s(inode);
          vderxy2(1,1) += evelnp(1,inode)*shp.vely.dydy.s(inode);
          vderxy2(1,2) += evelnp(1,inode)*shp.vely.dzdz.s(inode);
          vderxy2(1,3) += evelnp(1,inode)*shp.vely.dxdy.s(inode);
          vderxy2(1,4) += evelnp(1,inode)*shp.vely.dxdz.s(inode);
          vderxy2(1,5) += evelnp(1,inode)*shp.vely.dydz.s(inode);

          vderxy2(2,0) += evelnp(2,inode)*shp.velz.dxdx.s(inode);
          vderxy2(2,1) += evelnp(2,inode)*shp.velz.dydy.s(inode);
          vderxy2(2,2) += evelnp(2,inode)*shp.velz.dzdz.s(inode);
          vderxy2(2,3) += evelnp(2,inode)*shp.velz.dxdy.s(inode);
          vderxy2(2,4) += evelnp(2,inode)*shp.velz.dxdz.s(inode);
          vderxy2(2,5) += evelnp(2,inode)*shp.velz.dydz.s(inode);

          const int gid = nodeids[inode];
          const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);

          for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
              enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
          {
            if (enrfield->getField() == XFEM::PHYSICS::Veln)
            {
              vderxy2(0,0) += evelnp(3,velncounter)*shp.velx.dxdx.n(velncounter);
              vderxy2(0,1) += evelnp(3,velncounter)*shp.velx.dydy.n(velncounter);
              vderxy2(0,2) += evelnp(3,velncounter)*shp.velx.dzdz.n(velncounter);
              vderxy2(0,3) += evelnp(3,velncounter)*shp.velx.dxdy.n(velncounter);
              vderxy2(0,4) += evelnp(3,velncounter)*shp.velx.dxdz.n(velncounter);
              vderxy2(0,5) += evelnp(3,velncounter)*shp.velx.dydz.n(velncounter);

              vderxy2(1,0) += evelnp(3,velncounter)*shp.vely.dxdx.n(velncounter);
              vderxy2(1,1) += evelnp(3,velncounter)*shp.vely.dydy.n(velncounter);
              vderxy2(1,2) += evelnp(3,velncounter)*shp.vely.dzdz.n(velncounter);
              vderxy2(1,3) += evelnp(3,velncounter)*shp.vely.dxdy.n(velncounter);
              vderxy2(1,4) += evelnp(3,velncounter)*shp.vely.dxdz.n(velncounter);
              vderxy2(1,5) += evelnp(3,velncounter)*shp.vely.dydz.n(velncounter);

              vderxy2(2,0) += evelnp(3,velncounter)*shp.velz.dxdx.n(velncounter);
              vderxy2(2,1) += evelnp(3,velncounter)*shp.velz.dydy.n(velncounter);
              vderxy2(2,2) += evelnp(3,velncounter)*shp.velz.dzdz.n(velncounter);
              vderxy2(2,3) += evelnp(3,velncounter)*shp.velz.dxdy.n(velncounter);
              vderxy2(2,4) += evelnp(3,velncounter)*shp.velz.dxdz.n(velncounter);
              vderxy2(2,5) += evelnp(3,velncounter)*shp.velz.dydz.n(velncounter);

              velncounter += 1;
            }
          }
        }
        if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
      }

      //----------------------------------
      // get pressure at integration point
      //----------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp_pres.d0(iparam);

      //--------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp_pres.dx(iparam);
        gradp(1) += eprenp(iparam) * shp_pres.dy(iparam);
        gradp(2) += eprenp(iparam) * shp_pres.dz(iparam);
      }

      //------------------------------------------
      // compute stabilization parameters (3 taus)
      //------------------------------------------
      //const double vel_norm = gpvelnp.Norm2();
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
      for (std::size_t isd = 0; isd < nsd; isd++)
      {
        for (std::size_t inode = 0; inode < numnode; inode++)
        {
          bodyforce(isd) += edead(isd,inode) * funct(inode);
        }
      }
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densaf);

      for (size_t isd = 0; isd < 3; ++isd)
        rhsint(isd) = histvec(isd) + timefac*bodyforce(isd);

      //-----------------------------------------------------------------------
      // compute numerical representation of some single operators (not scaled)
      //-----------------------------------------------------------------------

      // convective term:  u_old * grad u_old:
      // remark: conv_old = vderxy(i,j)*gpvelnp(j);
      static LINALG::Matrix<3,1> conv_old;
      conv_old.Clear();
      conv_old.Multiply(vderxy,gpvelnp);

      // viscous term:  div epsilon(u_old)
      LINALG::Matrix<3,1> visc_old;
      visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
      visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
      visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

      // convective part of linearized convective term: u_old * grad (funct)
      // remark: u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> enr_conv_c(false);
      // remark: initialization is essentiual here, since contributions are added up!
      enr_conv_c.Clear();
      enr_conv_c.velx.d0.s.Clear();
      enr_conv_c.vely.d0.s.Clear();
      enr_conv_c.velz.d0.s.Clear();
      enr_conv_c.velx.d0.n.Clear();
      enr_conv_c.vely.d0.n.Clear();
      enr_conv_c.velz.d0.n.Clear();

      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        enr_conv_c.velx.d0.s(iparam) += shp.velx.dx.s(iparam)*gpvelnp(0);
        enr_conv_c.velx.d0.s(iparam) += shp.velx.dy.s(iparam)*gpvelnp(1);
        enr_conv_c.velx.d0.s(iparam) += shp.velx.dz.s(iparam)*gpvelnp(2);
      }
      enr_conv_c.vely.d0.s = enr_conv_c.velx.d0.s;
      enr_conv_c.velz.d0.s = enr_conv_c.velx.d0.s;
      for (size_t iparam = 0; iparam < numparamveln; ++iparam)
      {
        enr_conv_c.velx.d0.n(iparam) += shp.velx.dx.n(iparam)*gpvelnp(0);
        enr_conv_c.velx.d0.n(iparam) += shp.velx.dy.n(iparam)*gpvelnp(1);
        enr_conv_c.velx.d0.n(iparam) += shp.velx.dz.n(iparam)*gpvelnp(2);

        enr_conv_c.vely.d0.n(iparam) += shp.vely.dx.n(iparam)*gpvelnp(0);
        enr_conv_c.vely.d0.n(iparam) += shp.vely.dy.n(iparam)*gpvelnp(1);
        enr_conv_c.vely.d0.n(iparam) += shp.vely.dz.n(iparam)*gpvelnp(2);

        enr_conv_c.velz.d0.n(iparam) += shp.velz.dx.n(iparam)*gpvelnp(0);
        enr_conv_c.velz.d0.n(iparam) += shp.velz.dy.n(iparam)*gpvelnp(1);
        enr_conv_c.velz.d0.n(iparam) += shp.velz.dz.n(iparam)*gpvelnp(2);
      }

//if(ele->Id()==2 and icell == 0 and iquad==0)
//{
//  cout << "ele id " << ele->Id() << endl;
//  cout << posXiDomain << endl;
//if (cell->getDomainPlus()) cout << "plus" << endl;
//else cout << "minus" << endl;
//if (ASSTYPE == XFEM::xfem_assembly) cout << "XFEM assembly" << endl;
//if (ASSTYPE == XFEM::standard_assembly) cout << "standard assembly" << endl;
//enr_conv_c.Print();
//cout << "gpvelnp " << gpvelnp << endl;
//cout << "evelnp " << evelnp << endl;
//}
      // viscous term  - grad * epsilon(u)
      /*    /                                              \
                 |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
               1 |                                                |
               - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
               2 |                                                |
                 |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
                  \                                              /
       */
      static XFEM::EnrViscs2N<shpVecSizeVel> enr_viscs2;
      //            enr_viscs2.Clear();
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        enr_viscs2.xx.s(iparam) = 0.5 * (2.0 * shp.velx.dxdx.s(iparam) + shp.velx.dydy.s(iparam) + shp.velx.dzdz.s(iparam));
        enr_viscs2.xy.s(iparam) = 0.5 *  shp.vely.dydx.s(iparam);
        enr_viscs2.xz.s(iparam) = 0.5 *  shp.velz.dzdx.s(iparam);
        enr_viscs2.yx.s(iparam) = 0.5 *  shp.velx.dxdy.s(iparam);
        enr_viscs2.yy.s(iparam) = 0.5 * (shp.vely.dxdx.s(iparam) + 2.0 * shp.vely.dydy.s(iparam) + shp.vely.dzdz.s(iparam));
        enr_viscs2.yz.s(iparam) = 0.5 *  shp.velz.dzdy.s(iparam);
        enr_viscs2.zx.s(iparam) = 0.5 *  shp.velx.dxdz.s(iparam);
        enr_viscs2.zy.s(iparam) = 0.5 *  shp.vely.dydz.s(iparam);
        enr_viscs2.zz.s(iparam) = 0.5 * (shp.velz.dxdx.s(iparam) + shp.velz.dydy.s(iparam) + 2.0 * shp.velz.dzdz.s(iparam));
      }
      for (size_t iparam = 0; iparam < numparamveln; ++iparam)
      {
        enr_viscs2.xx.n(iparam) = 0.5 * (2.0 * shp.velx.dxdx.n(iparam) + shp.velx.dydy.n(iparam) + shp.velx.dzdz.n(iparam));
        enr_viscs2.xy.n(iparam) = 0.5 *  shp.vely.dydx.n(iparam);
        enr_viscs2.xz.n(iparam) = 0.5 *  shp.velz.dzdx.n(iparam);
        enr_viscs2.yx.n(iparam) = 0.5 *  shp.velx.dxdy.n(iparam);
        enr_viscs2.yy.n(iparam) = 0.5 * (shp.vely.dxdx.n(iparam) + 2.0 * shp.vely.dydy.n(iparam) + shp.vely.dzdz.n(iparam));
        enr_viscs2.yz.n(iparam) = 0.5 *  shp.velz.dzdy.n(iparam);
        enr_viscs2.zx.n(iparam) = 0.5 *  shp.velx.dxdz.n(iparam);
        enr_viscs2.zy.n(iparam) = 0.5 *  shp.vely.dydz.n(iparam);
        enr_viscs2.zz.n(iparam) = 0.5 * (shp.velz.dxdx.n(iparam) + shp.velz.dydy.n(iparam) + 2.0 * shp.velz.dzdz.n(iparam));
      }

//      if (ASSTYPE== XFEM::xfem_assembly)
//      {
//        if(ele->Bisected())
//        {
//          cout << "eprenp " << eprenp << endl;
//          cout << "evelnp " << evelnp << endl;
//          cout << "funct" << funct << endl;
//          cout << "derxy" << derxy << endl;
//          cout << "shp" << endl;
//          shp.Print();
//          cout << "gpvelnp" << gpvelnp << endl;
//          cout << "vderxy " << vderxy << endl;
//          cout << "vderxy2" << vderxy2 << endl;
//          cout << "pres " << pres << endl;
//          cout << "gradp " << gradp << endl;
//          cout << "visc_old " << visc_old << endl;
//          cout << "enr_conv_c" << endl;
//          enr_conv_c.Print();
//
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.xx " << enr_viscs2.xx.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.xx " << enr_viscs2.xx.n(iparam);
//          }
//          cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.xy " << enr_viscs2.xy.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.xy " << enr_viscs2.xy.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.xz " << enr_viscs2.xz.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.xz " << enr_viscs2.xz.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.yx " << enr_viscs2.yx.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.yx " << enr_viscs2.yx.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.yy " << enr_viscs2.yy.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.yy " << enr_viscs2.yy.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.yz " << enr_viscs2.yz.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.yz " << enr_viscs2.yz.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.zx " << enr_viscs2.zx.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.zx " << enr_viscs2.zx.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.zy " << enr_viscs2.zy.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.zy " << enr_viscs2.zy.n(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.zz " << enr_viscs2.zz.s(iparam);
//          } cout << endl;
//          for(size_t iparam=0;iparam<shpVecSizeVel;++iparam)
//          {
//            cout << "enr_viscs2.zz " << enr_viscs2.zz.n(iparam);
//          } cout << endl;
//        }
//      }
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
      LINALG::Matrix<3,1> res_old(true);
      for (size_t isd = 0; isd < nsd; ++isd)
        res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

      if (instationary)
      {
        for(size_t isd=0; isd < nsd; ++isd)
          res_old(isd) += densaf * gpvelnp(isd);
      }

      // integration factors and coefficients of single terms
      const double timefacfac = timefac * fac;

      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegralsNormal<DISTYPE,ASSTYPE,NUMDOF,shpVecSizeVel,shpVecSizePres>(
          assembler, shp, shp_pres, fac, timefac, timefacfac, densaf, dynvisc,
          gpvelnp, pres, gpaccam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_viscs2,
          instationary, genalpha, newton, pstab, supg, graddiv,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
    } // end loop over gauss points
  ++icell;
  } // end loop over integration cells

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          DRT::Element::DiscretizationType CELLDISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          class M1, class M2, class M3, class V1, class V2>
void Nitsche_SysmatBoundaryCellNormal(
    GEO::BoundaryIntCell              cell,
    const DRT::UTILS::GaussRule2D&    intrule2D,
    M2&                               xyze,
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const bool&                       higher_order_ele,
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M3&                         egradphi, ///<smoothed nodal G-function gradient values for this element
    const double&                     ju_n,
    const double&                     jflux_n,
    const double&                       timefac,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const XFEM::ElementEnrichmentValues&       enrvals_plus,
    const XFEM::ElementEnrichmentValues&       enrvals_minus,
    const double&                        kappa_plus,    /// weights of average in plus domain
    const double&                        kappa_minus,    /// weights of average in minus domain
    const double&                        visc_plus,    /// viscosity in plus domain
    const double&                        visc_minus,   /// viscosity in minus domain
    const double&                        alphau,       /// Nitsche's stabilization parameter
    const double&                        alpha_conv,   /// stabilization parameter for convection dominated flows
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType     veljumptype,
    const INPAR::COMBUST::FluxJumpType         fluxjumptype,
    const bool                                 smoothed_boundary_integration
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparamveln = XFEM::NumParam<0,ASSTYPE>::get(dofman, XFEM::PHYSICS::Veln);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  //----------------------------------------------
  // get Gaussian points for this integration cell
  //----------------------------------------------
  // TODO: are 6 Gaussian points enough?
  //const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
  //const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_37point);
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
//    COMBUST::RHS_LaplaceBeltrami_1DLineIntegrals_Nitsche<DISTYPE,CELLDISTYPE,ASSTYPE,NUMDOF>(
//        cellXiDomain, xyze, egradphi,
//        enrvals_plus,enrvals_minus,
//        kappa_plus,kappa_minus,numparamvelx,timefac,surftensapprox,surftenscoeff,assembler);
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

    if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");

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
    normal.Clear();

    // call smoothed normals for an accurate evaluation of boundary integration terms
    if ((smoothed_boundary_integration == true)
        and (surftensapprox != INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed
        or surftensapprox != INPAR::COMBUST::surface_tension_approx_divgrad_normal))
    {
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      normal.Clear();

      // loop over nodes
      for(size_t i = 0; i< numnode; i++)
      {
        LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
        nodal_grad_xyz_tmp.Clear();

        normal(0)+=egradphi(0,i)*funct_gp(i);
        normal(1)+=egradphi(1,i)*funct_gp(i);
        normal(2)+=egradphi(2,i)*funct_gp(i);
      }

      // 4. normalize final gradient
      // remark: n = -grad(\phi)
#ifdef COMBUST_2D
      normal(2) = 0.0;
#endif
      const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal.Scale(-1.0/norm);
    }
    else
    {
      // normal vector points from Omega^+ to Omega^- (n = n^+ = -n^-)
      GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);
    }

#ifdef COLLAPSE_FLAME
      normal.Clear();
      // get shape functions at gaussian point
      static LINALG::Matrix<numnode,1> funct_gp;
      funct_gp.Clear();
      DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      for (unsigned i=0;i<numnode;i++)
      {
        normal(0) += funct_gp(i)*xyze(0,i);
        normal(1) += funct_gp(i)*xyze(1,i);
      }
      const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal.Scale(1.0/norm);

      LINALG::Matrix<3,1> point(true);
      for (unsigned inode=0;inode<numnode;++inode)
      {
        point(0) += funct_gp(inode)*xyze(0,inode);
        point(1) += funct_gp(inode)*xyze(1,inode);
      }
#endif
//if (ele->Id()==116)
////  or
////    ele->Id()==226 or
////    ele->Id()==326 or
////    ele->Id()==216)
//{
//  cout << "Normalenvektor Boundary Sysmat Viertelposition" << normal << endl;
//}
//if (ele->Id()==141)
////    or
////    ele->Id()==308 or
////    ele->Id()==201 or
////    ele->Id()==133)
//{
//  cout << "Normalenvektor Boundary Sysmat Achtelposition" << normal << endl;
//}
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
    ju(0) = ju_n*normal(0);
    ju(1) = ju_n*normal(1);
    ju(2) = ju_n*normal(2);
    //cout << "ju " << ju << endl;
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
      //cout << "jflux " << jflux << endl;
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
      DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
    }
    else
    {
      derxy2.Clear();
    }

    //--------------------------------------------------------------------------------
    // rearrange (enriched) shape functions and derivatives as approximation functions
    //--------------------------------------------------------------------------------
    const size_t shpVecSizeVel = COMBUST::SizeFacVel<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    const size_t shpVecSizePres = COMBUST::SizeFacPres<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
    static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> enrfunct_vel_plus(false);
    enrfunct_vel_plus.Clear();
    static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> enrfunct_vel_minus(false);
    enrfunct_vel_minus.Clear();

    static XFEM::ApproxFuncNormalVector<2,shpVecSizeVel> shp_plus(false);
    shp_plus.Clear();
    static XFEM::ApproxFuncNormalVector<2,shpVecSizeVel> shp_minus(false);
    shp_minus.Clear();

    static LINALG::Matrix<shpVecSizePres,1> enrfunct_pres_plus;
    enrfunct_pres_plus.Clear();
    static LINALG::Matrix<shpVecSizePres,1> enrfunct_pres_minus;
    enrfunct_pres_minus.Clear();

    // shape functions for nodal parameters (dofs) on plus and minus side
    enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Pres, funct, enrfunct_pres_plus);
    enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Pres, funct, enrfunct_pres_minus);

    // fill approximation functions for XFEM
    for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
    {
      shp_plus.velx.d0.s(iparam) = funct(iparam);
      shp_plus.vely.d0.s(iparam) = funct(iparam);
      shp_plus.velz.d0.s(iparam) = funct(iparam);

      shp_plus.velx.dx.s(iparam) = derxy(0,iparam);
      shp_plus.vely.dx.s(iparam) = derxy(0,iparam);
      shp_plus.velz.dx.s(iparam) = derxy(0,iparam);

      shp_plus.velx.dy.s(iparam) = derxy(1,iparam);
      shp_plus.vely.dy.s(iparam) = derxy(1,iparam);
      shp_plus.velz.dy.s(iparam) = derxy(1,iparam);

      shp_plus.velx.dz.s(iparam) = derxy(2,iparam);
      shp_plus.vely.dz.s(iparam) = derxy(2,iparam);
      shp_plus.velz.dz.s(iparam) = derxy(2,iparam);

      shp_plus.velx.dxdx.s(iparam) = derxy2(0,iparam);
      shp_plus.vely.dxdx.s(iparam) = derxy2(0,iparam);
      shp_plus.velz.dxdx.s(iparam) = derxy2(0,iparam);

      shp_plus.velx.dxdy.s(iparam) = derxy2(3,iparam);
      shp_plus.vely.dxdy.s(iparam) = derxy2(3,iparam);
      shp_plus.velz.dxdy.s(iparam) = derxy2(3,iparam);

      shp_plus.velx.dxdz.s(iparam) = derxy2(4,iparam);
      shp_plus.vely.dxdz.s(iparam) = derxy2(4,iparam);
      shp_plus.velz.dxdz.s(iparam) = derxy2(4,iparam);

      shp_plus.velx.dydx.s(iparam) = shp_plus.velx.dxdy.s(iparam);
      shp_plus.vely.dydx.s(iparam) = shp_plus.vely.dxdy.s(iparam);
      shp_plus.velz.dydx.s(iparam) = shp_plus.velz.dxdy.s(iparam);

      shp_plus.velx.dydy.s(iparam) = derxy2(1,iparam);
      shp_plus.vely.dydy.s(iparam) = derxy2(1,iparam);
      shp_plus.velz.dydy.s(iparam) = derxy2(1,iparam);

      shp_plus.velx.dydz.s(iparam) = derxy2(5,iparam);
      shp_plus.vely.dydz.s(iparam) = derxy2(5,iparam);
      shp_plus.velz.dydz.s(iparam) = derxy2(5,iparam);

      shp_plus.velx.dzdx.s(iparam) = shp_plus.velx.dxdz.s(iparam);
      shp_plus.vely.dzdx.s(iparam) = shp_plus.vely.dxdz.s(iparam);
      shp_plus.velz.dzdx.s(iparam) = shp_plus.velz.dxdz.s(iparam);

      shp_plus.velx.dzdy.s(iparam) = shp_plus.velx.dydz.s(iparam);
      shp_plus.vely.dzdy.s(iparam) = shp_plus.vely.dydz.s(iparam);
      shp_plus.velz.dzdy.s(iparam) = shp_plus.velz.dydz.s(iparam);

      shp_plus.velx.dzdz.s(iparam) = derxy2(2,iparam);
      shp_plus.vely.dzdz.s(iparam) = derxy2(2,iparam);
      shp_plus.velz.dzdz.s(iparam) = derxy2(2,iparam);
    }
    // fill approximation functions for XFEM
    for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
    {
      shp_minus.velx.d0.s(iparam) = funct(iparam);
      shp_minus.vely.d0.s(iparam) = funct(iparam);
      shp_minus.velz.d0.s(iparam) = funct(iparam);

      shp_minus.velx.dx.s(iparam) = derxy(0,iparam);
      shp_minus.vely.dx.s(iparam) = derxy(0,iparam);
      shp_minus.velz.dx.s(iparam) = derxy(0,iparam);

      shp_minus.velx.dy.s(iparam) = derxy(1,iparam);
      shp_minus.vely.dy.s(iparam) = derxy(1,iparam);
      shp_minus.velz.dy.s(iparam) = derxy(1,iparam);

      shp_minus.velx.dz.s(iparam) = derxy(2,iparam);
      shp_minus.vely.dz.s(iparam) = derxy(2,iparam);
      shp_minus.velz.dz.s(iparam) = derxy(2,iparam);

      shp_minus.velx.dxdx.s(iparam) = derxy2(0,iparam);
      shp_minus.vely.dxdx.s(iparam) = derxy2(0,iparam);
      shp_minus.velz.dxdx.s(iparam) = derxy2(0,iparam);

      shp_minus.velx.dxdy.s(iparam) = derxy2(3,iparam);
      shp_minus.vely.dxdy.s(iparam) = derxy2(3,iparam);
      shp_minus.velz.dxdy.s(iparam) = derxy2(3,iparam);

      shp_minus.velx.dxdz.s(iparam) = derxy2(4,iparam);
      shp_minus.vely.dxdz.s(iparam) = derxy2(4,iparam);
      shp_minus.velz.dxdz.s(iparam) = derxy2(4,iparam);

      shp_minus.velx.dydx.s(iparam) = shp_minus.velx.dxdy.s(iparam);
      shp_minus.vely.dydx.s(iparam) = shp_minus.vely.dxdy.s(iparam);
      shp_minus.velz.dydx.s(iparam) = shp_minus.velz.dxdy.s(iparam);

      shp_minus.velx.dydy.s(iparam) = derxy2(1,iparam);
      shp_minus.vely.dydy.s(iparam) = derxy2(1,iparam);
      shp_minus.velz.dydy.s(iparam) = derxy2(1,iparam);

      shp_minus.velx.dydz.s(iparam) = derxy2(5,iparam);
      shp_minus.vely.dydz.s(iparam) = derxy2(5,iparam);
      shp_minus.velz.dydz.s(iparam) = derxy2(5,iparam);

      shp_minus.velx.dzdx.s(iparam) = shp_minus.velx.dxdz.s(iparam);
      shp_minus.vely.dzdx.s(iparam) = shp_minus.vely.dxdz.s(iparam);
      shp_minus.velz.dzdx.s(iparam) = shp_minus.velz.dxdz.s(iparam);

      shp_minus.velx.dzdy.s(iparam) = shp_minus.velx.dydz.s(iparam);
      shp_minus.vely.dzdy.s(iparam) = shp_minus.vely.dydz.s(iparam);
      shp_minus.velz.dzdy.s(iparam) = shp_minus.velz.dydz.s(iparam);

      shp_minus.velx.dzdz.s(iparam) = derxy2(2,iparam);
      shp_minus.vely.dzdz.s(iparam) = derxy2(2,iparam);
      shp_minus.velz.dzdz.s(iparam) = derxy2(2,iparam);
    }
    enrvals_plus.ComputeNormalShapeFunction(funct, derxy, derxy2, egradphi,
#ifdef COLLAPSE_FLAME
          normal,
          point,
#endif
        shp_plus);
    enrvals_minus.ComputeNormalShapeFunction(funct, derxy, derxy2, egradphi,
#ifdef COLLAPSE_FLAME
        normal,
        point,
#endif
        shp_minus);

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> shp_vel_jump(false);     // [[ ]] notation
    shp_vel_jump.Clear();
    static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> shp_vel_mean(false);     //  { }  notation
    shp_vel_mean.Clear();
    static XFEM::ApproxFuncNormalVector<1,shpVecSizeVel> shp_vel_mean_inv(false); //  < >  notation
    shp_vel_mean_inv.Clear();
    static XFEM::ApproxFunc<0,shpVecSizePres> shp_pres_jump;           // [[ ]] notation
    static XFEM::ApproxFunc<0,shpVecSizePres> shp_pres_mean;           //  { }  notation

    static XFEM::ApproxFuncNormalVector<2,shpVecSizeVel> shp_mean_visc(false); //  { }  notation
    shp_mean_visc.Clear();
    static XFEM::ApproxFuncNormalVector<2,shpVecSizeVel> shp_jump_visc(false); // [[ ]] notation
    shp_jump_visc.Clear();

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparampres; ++iparam)
    {
      shp_pres_jump.d0(iparam) = enrfunct_pres_minus(iparam) - enrfunct_pres_plus(iparam);
      shp_pres_mean.d0(iparam) = kappa_plus*enrfunct_pres_plus(iparam) + kappa_minus*enrfunct_pres_minus(iparam);
    }

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      // standard jump shape functions are zero by definition
      //shp_vel_jump.velx.d0.s(iparam) = 0.0;
      //shp_vel_jump.vely.d0.s(iparam) = 0.0;
      //shp_vel_jump.velz.d0.s(iparam) = 0.0;

      shp_vel_mean.velx.d0.s(iparam) = funct(iparam);
      shp_vel_mean.vely.d0.s(iparam) = funct(iparam);
      shp_vel_mean.velz.d0.s(iparam) = funct(iparam);

//      shp_vel_mean_inv.velx.d0.s(iparam) = funct(iparam);
//      shp_vel_mean_inv.vely.d0.s(iparam) = funct(iparam);
//      shp_vel_mean_inv.velz.d0.s(iparam) = funct(iparam);
//
//      shp_vel_mean_inv.velx.dx.s(iparam) = derxy(0,iparam);
//      shp_vel_mean_inv.vely.dx.s(iparam) = derxy(0,iparam);
//      shp_vel_mean_inv.velz.dx.s(iparam) = derxy(0,iparam);
//
//      shp_vel_mean_inv.velx.dy.s(iparam) = derxy(1,iparam);
//      shp_vel_mean_inv.vely.dy.s(iparam) = derxy(1,iparam);
//      shp_vel_mean_inv.velz.dy.s(iparam) = derxy(1,iparam);
//
//      shp_vel_mean_inv.velx.dz.s(iparam) = derxy(2,iparam);
//      shp_vel_mean_inv.vely.dz.s(iparam) = derxy(2,iparam);
//      shp_vel_mean_inv.velz.dz.s(iparam) = derxy(2,iparam);

      shp_vel_mean_inv.velx.d0.s(iparam) = kappa_minus*shp_plus.velx.d0.s(iparam) + kappa_plus*shp_minus.velx.d0.s(iparam);
      shp_vel_mean_inv.vely.d0.s(iparam) = kappa_minus*shp_plus.vely.d0.s(iparam) + kappa_plus*shp_minus.vely.d0.s(iparam);
      shp_vel_mean_inv.velz.d0.s(iparam) = kappa_minus*shp_plus.velz.d0.s(iparam) + kappa_plus*shp_minus.velz.d0.s(iparam);

      shp_vel_mean_inv.velx.dx.s(iparam) = kappa_minus*shp_plus.velx.dx.s(iparam) + kappa_plus*shp_minus.velx.dx.s(iparam);
      shp_vel_mean_inv.vely.dx.s(iparam) = kappa_minus*shp_plus.vely.dx.s(iparam) + kappa_plus*shp_minus.vely.dx.s(iparam);
      shp_vel_mean_inv.velz.dx.s(iparam) = kappa_minus*shp_plus.velz.dx.s(iparam) + kappa_plus*shp_minus.velz.dx.s(iparam);

      shp_vel_mean_inv.velx.dy.s(iparam) = kappa_minus*shp_plus.velx.dy.s(iparam) + kappa_plus*shp_minus.velx.dy.s(iparam);
      shp_vel_mean_inv.vely.dy.s(iparam) = kappa_minus*shp_plus.vely.dy.s(iparam) + kappa_plus*shp_minus.vely.dy.s(iparam);
      shp_vel_mean_inv.velz.dy.s(iparam) = kappa_minus*shp_plus.velz.dy.s(iparam) + kappa_plus*shp_minus.velz.dy.s(iparam);

      shp_vel_mean_inv.velx.dz.s(iparam) = kappa_minus*shp_plus.velx.dz.s(iparam) + kappa_plus*shp_minus.velx.dz.s(iparam);
      shp_vel_mean_inv.vely.dz.s(iparam) = kappa_minus*shp_plus.vely.dz.s(iparam) + kappa_plus*shp_minus.vely.dz.s(iparam);
      shp_vel_mean_inv.velz.dz.s(iparam) = kappa_minus*shp_plus.velz.dz.s(iparam) + kappa_plus*shp_minus.velz.dz.s(iparam);


      shp_mean_visc.velx.d0.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.d0.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.d0.s(iparam);
      shp_mean_visc.vely.d0.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.d0.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.d0.s(iparam);
      shp_mean_visc.velz.d0.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.d0.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.d0.s(iparam);

      shp_mean_visc.velx.dx.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dx.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dx.s(iparam);
      shp_mean_visc.vely.dx.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dx.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dx.s(iparam);
      shp_mean_visc.velz.dx.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dx.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dx.s(iparam);

      shp_mean_visc.velx.dy.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dy.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dy.s(iparam);
      shp_mean_visc.vely.dy.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dy.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dy.s(iparam);
      shp_mean_visc.velz.dy.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dy.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dy.s(iparam);

      shp_mean_visc.velx.dz.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dz.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dz.s(iparam);
      shp_mean_visc.vely.dz.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dz.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dz.s(iparam);
      shp_mean_visc.velz.dz.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dz.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dz.s(iparam);

      shp_mean_visc.velx.dxdx.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dxdx.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dxdx.s(iparam);
      shp_mean_visc.vely.dxdx.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dxdx.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dxdx.s(iparam);
      shp_mean_visc.velz.dxdx.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dxdx.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dxdx.s(iparam);

      shp_mean_visc.velx.dxdy.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dxdy.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dxdy.s(iparam);
      shp_mean_visc.vely.dxdy.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dxdy.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dxdy.s(iparam);
      shp_mean_visc.velz.dxdy.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dxdy.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dxdy.s(iparam);

      shp_mean_visc.velx.dxdz.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dxdz.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dxdz.s(iparam);
      shp_mean_visc.vely.dxdz.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dxdz.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dxdz.s(iparam);
      shp_mean_visc.velz.dxdz.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dxdz.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dxdz.s(iparam);

      shp_mean_visc.velx.dydx.s(iparam) = shp_mean_visc.velx.dxdy.s(iparam);
      shp_mean_visc.vely.dydx.s(iparam) = shp_mean_visc.vely.dxdy.s(iparam);
      shp_mean_visc.velz.dydx.s(iparam) = shp_mean_visc.velz.dxdy.s(iparam);

      shp_mean_visc.velx.dydy.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dydy.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dydy.s(iparam);
      shp_mean_visc.vely.dydy.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dydy.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dydy.s(iparam);
      shp_mean_visc.velz.dydy.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dydy.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dydy.s(iparam);

      shp_mean_visc.velx.dydz.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dydz.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dydz.s(iparam);
      shp_mean_visc.vely.dydz.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dydz.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dydz.s(iparam);
      shp_mean_visc.velz.dydz.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dydz.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dydz.s(iparam);

      shp_mean_visc.velx.dzdx.s(iparam) = shp_mean_visc.velx.dxdz.s(iparam);
      shp_mean_visc.vely.dzdx.s(iparam) = shp_mean_visc.vely.dxdz.s(iparam);
      shp_mean_visc.velz.dzdx.s(iparam) = shp_mean_visc.velz.dxdz.s(iparam);

      shp_mean_visc.velx.dzdy.s(iparam) = shp_mean_visc.velx.dydz.s(iparam);
      shp_mean_visc.vely.dzdy.s(iparam) = shp_mean_visc.vely.dydz.s(iparam);
      shp_mean_visc.velz.dzdy.s(iparam) = shp_mean_visc.velz.dydz.s(iparam);

      shp_mean_visc.velx.dzdz.s(iparam) = kappa_plus*visc_plus*shp_plus.velx.dzdz.s(iparam) + kappa_minus*visc_minus*shp_minus.velx.dzdz.s(iparam);
      shp_mean_visc.vely.dzdz.s(iparam) = kappa_plus*visc_plus*shp_plus.vely.dzdz.s(iparam) + kappa_minus*visc_minus*shp_minus.vely.dzdz.s(iparam);
      shp_mean_visc.velz.dzdz.s(iparam) = kappa_plus*visc_plus*shp_plus.velz.dzdz.s(iparam) + kappa_minus*visc_minus*shp_minus.velz.dzdz.s(iparam);

      shp_jump_visc.velx.d0.s(iparam) = visc_minus*shp_minus.velx.d0.s(iparam) - visc_plus*shp_plus.velx.d0.s(iparam);
      shp_jump_visc.vely.d0.s(iparam) = visc_minus*shp_minus.vely.d0.s(iparam) - visc_plus*shp_plus.vely.d0.s(iparam);
      shp_jump_visc.velz.d0.s(iparam) = visc_minus*shp_minus.velz.d0.s(iparam) - visc_plus*shp_plus.velz.d0.s(iparam);

      shp_jump_visc.velx.dx.s(iparam) = visc_minus*shp_minus.velx.dx.s(iparam) - visc_plus*shp_plus.velx.dx.s(iparam);
      shp_jump_visc.vely.dx.s(iparam) = visc_minus*shp_minus.vely.dx.s(iparam) - visc_plus*shp_plus.vely.dx.s(iparam);
      shp_jump_visc.velz.dx.s(iparam) = visc_minus*shp_minus.velz.dx.s(iparam) - visc_plus*shp_plus.velz.dx.s(iparam);

      shp_jump_visc.velx.dy.s(iparam) = visc_minus*shp_minus.velx.dy.s(iparam) - visc_plus*shp_plus.velx.dy.s(iparam);
      shp_jump_visc.vely.dy.s(iparam) = visc_minus*shp_minus.vely.dy.s(iparam) - visc_plus*shp_plus.vely.dy.s(iparam);
      shp_jump_visc.velz.dy.s(iparam) = visc_minus*shp_minus.velz.dy.s(iparam) - visc_plus*shp_plus.velz.dy.s(iparam);

      shp_jump_visc.velx.dz.s(iparam) = visc_minus*shp_minus.velx.dz.s(iparam) - visc_plus*shp_plus.velx.dz.s(iparam);
      shp_jump_visc.vely.dz.s(iparam) = visc_minus*shp_minus.vely.dz.s(iparam) - visc_plus*shp_plus.vely.dz.s(iparam);
      shp_jump_visc.velz.dz.s(iparam) = visc_minus*shp_minus.velz.dz.s(iparam) - visc_plus*shp_plus.velz.dz.s(iparam);

      shp_jump_visc.velx.dxdx.s(iparam) = visc_minus*shp_minus.velx.dxdx.s(iparam) - visc_plus*shp_plus.velx.dxdx.s(iparam);
      shp_jump_visc.vely.dxdx.s(iparam) = visc_minus*shp_minus.vely.dxdx.s(iparam) - visc_plus*shp_plus.vely.dxdx.s(iparam);
      shp_jump_visc.velz.dxdx.s(iparam) = visc_minus*shp_minus.velz.dxdx.s(iparam) - visc_plus*shp_plus.velz.dxdx.s(iparam);

      shp_jump_visc.velx.dxdy.s(iparam) = visc_minus*shp_minus.velx.dxdy.s(iparam) - visc_plus*shp_plus.velx.dxdy.s(iparam);
      shp_jump_visc.vely.dxdy.s(iparam) = visc_minus*shp_minus.vely.dxdy.s(iparam) - visc_plus*shp_plus.vely.dxdy.s(iparam);
      shp_jump_visc.velz.dxdy.s(iparam) = visc_minus*shp_minus.velz.dxdy.s(iparam) - visc_plus*shp_plus.velz.dxdy.s(iparam);

      shp_jump_visc.velx.dxdz.s(iparam) = visc_minus*shp_minus.velx.dxdz.s(iparam) - visc_plus*shp_plus.velx.dxdz.s(iparam);
      shp_jump_visc.vely.dxdz.s(iparam) = visc_minus*shp_minus.vely.dxdz.s(iparam) - visc_plus*shp_plus.vely.dxdz.s(iparam);
      shp_jump_visc.velz.dxdz.s(iparam) = visc_minus*shp_minus.velz.dxdz.s(iparam) - visc_plus*shp_plus.velz.dxdz.s(iparam);

      shp_jump_visc.velx.dydx.s(iparam) = shp_jump_visc.velx.dxdy.s(iparam);
      shp_jump_visc.vely.dydx.s(iparam) = shp_jump_visc.vely.dxdy.s(iparam);
      shp_jump_visc.velz.dydx.s(iparam) = shp_jump_visc.velz.dxdy.s(iparam);

      shp_jump_visc.velx.dydy.s(iparam) = visc_minus*shp_minus.velx.dydy.s(iparam) - visc_plus*shp_plus.velx.dydy.s(iparam);
      shp_jump_visc.vely.dydy.s(iparam) = visc_minus*shp_minus.vely.dydy.s(iparam) - visc_plus*shp_plus.vely.dydy.s(iparam);
      shp_jump_visc.velz.dydy.s(iparam) = visc_minus*shp_minus.velz.dydy.s(iparam) - visc_plus*shp_plus.velz.dydy.s(iparam);

      shp_jump_visc.velx.dydz.s(iparam) = visc_minus*shp_minus.velx.dydz.s(iparam) - visc_plus*shp_plus.velx.dydz.s(iparam);
      shp_jump_visc.vely.dydz.s(iparam) = visc_minus*shp_minus.vely.dydz.s(iparam) - visc_plus*shp_plus.vely.dydz.s(iparam);
      shp_jump_visc.velz.dydz.s(iparam) = visc_minus*shp_minus.velz.dydz.s(iparam) - visc_plus*shp_plus.velz.dydz.s(iparam);

      shp_jump_visc.velx.dzdx.s(iparam) = shp_jump_visc.velx.dxdz.s(iparam);
      shp_jump_visc.vely.dzdx.s(iparam) = shp_jump_visc.vely.dxdz.s(iparam);
      shp_jump_visc.velz.dzdx.s(iparam) = shp_jump_visc.velz.dxdz.s(iparam);

      shp_jump_visc.velx.dzdy.s(iparam) = shp_jump_visc.velx.dydz.s(iparam);
      shp_jump_visc.vely.dzdy.s(iparam) = shp_jump_visc.vely.dydz.s(iparam);
      shp_jump_visc.velz.dzdy.s(iparam) = shp_jump_visc.velz.dydz.s(iparam);

      shp_jump_visc.velx.dzdz.s(iparam) = visc_minus*shp_minus.velx.dzdz.s(iparam) - visc_plus*shp_plus.velx.dzdz.s(iparam);
      shp_jump_visc.vely.dzdz.s(iparam) = visc_minus*shp_minus.vely.dzdz.s(iparam) - visc_plus*shp_plus.vely.dzdz.s(iparam);
      shp_jump_visc.velz.dzdz.s(iparam) = visc_minus*shp_minus.velz.dzdz.s(iparam) - visc_plus*shp_plus.velz.dzdz.s(iparam);
    }
    for (std::size_t iparam = 0; iparam < numparamveln; ++iparam)
    {
      shp_vel_jump.velx.d0.n(iparam) = shp_minus.velx.d0.n(iparam) - shp_plus.velx.d0.n(iparam);
      shp_vel_jump.vely.d0.n(iparam) = shp_minus.vely.d0.n(iparam) - shp_plus.vely.d0.n(iparam);
      shp_vel_jump.velz.d0.n(iparam) = shp_minus.velz.d0.n(iparam) - shp_plus.velz.d0.n(iparam);

      shp_vel_mean.velx.d0.n(iparam) = kappa_plus*shp_plus.velx.d0.n(iparam) + kappa_minus*shp_minus.velx.d0.n(iparam);
      shp_vel_mean.vely.d0.n(iparam) = kappa_plus*shp_plus.vely.d0.n(iparam) + kappa_minus*shp_minus.vely.d0.n(iparam);
      shp_vel_mean.velz.d0.n(iparam) = kappa_plus*shp_plus.velz.d0.n(iparam) + kappa_minus*shp_minus.velz.d0.n(iparam);

      shp_vel_mean_inv.velx.d0.n(iparam) = kappa_minus*shp_plus.velx.d0.n(iparam) + kappa_plus*shp_minus.velx.d0.n(iparam);
      shp_vel_mean_inv.vely.d0.n(iparam) = kappa_minus*shp_plus.vely.d0.n(iparam) + kappa_plus*shp_minus.vely.d0.n(iparam);
      shp_vel_mean_inv.velz.d0.n(iparam) = kappa_minus*shp_plus.velz.d0.n(iparam) + kappa_plus*shp_minus.velz.d0.n(iparam);

      shp_vel_mean_inv.velx.dx.n(iparam) = kappa_minus*shp_plus.velx.dx.n(iparam) + kappa_plus*shp_minus.velx.dx.n(iparam);
      shp_vel_mean_inv.vely.dx.n(iparam) = kappa_minus*shp_plus.vely.dx.n(iparam) + kappa_plus*shp_minus.vely.dx.n(iparam);
      shp_vel_mean_inv.velz.dx.n(iparam) = kappa_minus*shp_plus.velz.dx.n(iparam) + kappa_plus*shp_minus.velz.dx.n(iparam);

      shp_vel_mean_inv.velx.dy.n(iparam) = kappa_minus*shp_plus.velx.dy.n(iparam) + kappa_plus*shp_minus.velx.dy.n(iparam);
      shp_vel_mean_inv.vely.dy.n(iparam) = kappa_minus*shp_plus.vely.dy.n(iparam) + kappa_plus*shp_minus.vely.dy.n(iparam);
      shp_vel_mean_inv.velz.dy.n(iparam) = kappa_minus*shp_plus.velz.dy.n(iparam) + kappa_plus*shp_minus.velz.dy.n(iparam);

      shp_vel_mean_inv.velx.dz.n(iparam) = kappa_minus*shp_plus.velx.dz.n(iparam) + kappa_plus*shp_minus.velx.dz.n(iparam);
      shp_vel_mean_inv.vely.dz.n(iparam) = kappa_minus*shp_plus.vely.dz.n(iparam) + kappa_plus*shp_minus.vely.dz.n(iparam);
      shp_vel_mean_inv.velz.dz.n(iparam) = kappa_minus*shp_plus.velz.dz.n(iparam) + kappa_plus*shp_minus.velz.dz.n(iparam);

      shp_mean_visc.velx.d0.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.d0.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.d0.n(iparam);
      shp_mean_visc.vely.d0.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.d0.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.d0.n(iparam);
      shp_mean_visc.velz.d0.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.d0.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.d0.n(iparam);

      shp_mean_visc.velx.dx.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dx.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dx.n(iparam);
      shp_mean_visc.vely.dx.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dx.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dx.n(iparam);
      shp_mean_visc.velz.dx.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dx.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dx.n(iparam);

      shp_mean_visc.velx.dy.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dy.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dy.n(iparam);
      shp_mean_visc.vely.dy.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dy.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dy.n(iparam);
      shp_mean_visc.velz.dy.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dy.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dy.n(iparam);

      shp_mean_visc.velx.dz.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dz.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dz.n(iparam);
      shp_mean_visc.vely.dz.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dz.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dz.n(iparam);
      shp_mean_visc.velz.dz.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dz.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dz.n(iparam);

      shp_mean_visc.velx.dxdx.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dxdx.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dxdx.n(iparam);
      shp_mean_visc.vely.dxdx.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dxdx.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dxdx.n(iparam);
      shp_mean_visc.velz.dxdx.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dxdx.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dxdx.n(iparam);

      shp_mean_visc.velx.dxdy.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dxdy.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dxdy.n(iparam);
      shp_mean_visc.vely.dxdy.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dxdy.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dxdy.n(iparam);
      shp_mean_visc.velz.dxdy.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dxdy.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dxdy.n(iparam);

      shp_mean_visc.velx.dxdz.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dxdz.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dxdz.n(iparam);
      shp_mean_visc.vely.dxdz.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dxdz.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dxdz.n(iparam);
      shp_mean_visc.velz.dxdz.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dxdz.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dxdz.n(iparam);

      shp_mean_visc.velx.dydx.n(iparam) = shp_mean_visc.velx.dxdy.n(iparam);
      shp_mean_visc.vely.dydx.n(iparam) = shp_mean_visc.vely.dxdy.n(iparam);
      shp_mean_visc.velz.dydx.n(iparam) = shp_mean_visc.velz.dxdy.n(iparam);

      shp_mean_visc.velx.dydy.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dydy.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dydy.n(iparam);
      shp_mean_visc.vely.dydy.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dydy.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dydy.n(iparam);
      shp_mean_visc.velz.dydy.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dydy.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dydy.n(iparam);

      shp_mean_visc.velx.dydz.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dydz.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dydz.n(iparam);
      shp_mean_visc.vely.dydz.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dydz.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dydz.n(iparam);
      shp_mean_visc.velz.dydz.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dydz.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dydz.n(iparam);

      shp_mean_visc.velx.dzdx.n(iparam) = shp_mean_visc.velx.dxdz.n(iparam);
      shp_mean_visc.vely.dzdx.n(iparam) = shp_mean_visc.vely.dxdz.n(iparam);
      shp_mean_visc.velz.dzdx.n(iparam) = shp_mean_visc.velz.dxdz.n(iparam);

      shp_mean_visc.velx.dzdy.n(iparam) = shp_mean_visc.velx.dydz.n(iparam);
      shp_mean_visc.vely.dzdy.n(iparam) = shp_mean_visc.vely.dydz.n(iparam);
      shp_mean_visc.velz.dzdy.n(iparam) = shp_mean_visc.velz.dydz.n(iparam);

      shp_mean_visc.velx.dzdz.n(iparam) = kappa_plus*visc_plus*shp_plus.velx.dzdz.n(iparam) + kappa_minus*visc_minus*shp_minus.velx.dzdz.n(iparam);
      shp_mean_visc.vely.dzdz.n(iparam) = kappa_plus*visc_plus*shp_plus.vely.dzdz.n(iparam) + kappa_minus*visc_minus*shp_minus.vely.dzdz.n(iparam);
      shp_mean_visc.velz.dzdz.n(iparam) = kappa_plus*visc_plus*shp_plus.velz.dzdz.n(iparam) + kappa_minus*visc_minus*shp_minus.velz.dzdz.n(iparam);

      shp_jump_visc.velx.d0.n(iparam) = visc_minus*shp_minus.velx.d0.n(iparam) - visc_plus*shp_plus.velx.d0.n(iparam);
      shp_jump_visc.vely.d0.n(iparam) = visc_minus*shp_minus.vely.d0.n(iparam) - visc_plus*shp_plus.vely.d0.n(iparam);
      shp_jump_visc.velz.d0.n(iparam) = visc_minus*shp_minus.velz.d0.n(iparam) - visc_plus*shp_plus.velz.d0.n(iparam);

      shp_jump_visc.velx.dx.n(iparam) = visc_minus*shp_minus.velx.dx.n(iparam) - visc_plus*shp_plus.velx.dx.n(iparam);
      shp_jump_visc.vely.dx.n(iparam) = visc_minus*shp_minus.vely.dx.n(iparam) - visc_plus*shp_plus.vely.dx.n(iparam);
      shp_jump_visc.velz.dx.n(iparam) = visc_minus*shp_minus.velz.dx.n(iparam) - visc_plus*shp_plus.velz.dx.n(iparam);

      shp_jump_visc.velx.dy.n(iparam) = visc_minus*shp_minus.velx.dy.n(iparam) - visc_plus*shp_plus.velx.dy.n(iparam);
      shp_jump_visc.vely.dy.n(iparam) = visc_minus*shp_minus.vely.dy.n(iparam) - visc_plus*shp_plus.vely.dy.n(iparam);
      shp_jump_visc.velz.dy.n(iparam) = visc_minus*shp_minus.velz.dy.n(iparam) - visc_plus*shp_plus.velz.dy.n(iparam);

      shp_jump_visc.velx.dz.n(iparam) = visc_minus*shp_minus.velx.dz.n(iparam) - visc_plus*shp_plus.velx.dz.n(iparam);
      shp_jump_visc.vely.dz.n(iparam) = visc_minus*shp_minus.vely.dz.n(iparam) - visc_plus*shp_plus.vely.dz.n(iparam);
      shp_jump_visc.velz.dz.n(iparam) = visc_minus*shp_minus.velz.dz.n(iparam) - visc_plus*shp_plus.velz.dz.n(iparam);

      shp_jump_visc.velx.dxdx.n(iparam) = visc_minus*shp_minus.velx.dxdx.n(iparam) - visc_plus*shp_plus.velx.dxdx.n(iparam);
      shp_jump_visc.vely.dxdx.n(iparam) = visc_minus*shp_minus.vely.dxdx.n(iparam) - visc_plus*shp_plus.vely.dxdx.n(iparam);
      shp_jump_visc.velz.dxdx.n(iparam) = visc_minus*shp_minus.velz.dxdx.n(iparam) - visc_plus*shp_plus.velz.dxdx.n(iparam);

      shp_jump_visc.velx.dxdy.n(iparam) = visc_minus*shp_minus.velx.dxdy.n(iparam) - visc_plus*shp_plus.velx.dxdy.n(iparam);
      shp_jump_visc.vely.dxdy.n(iparam) = visc_minus*shp_minus.vely.dxdy.n(iparam) - visc_plus*shp_plus.vely.dxdy.n(iparam);
      shp_jump_visc.velz.dxdy.n(iparam) = visc_minus*shp_minus.velz.dxdy.n(iparam) - visc_plus*shp_plus.velz.dxdy.n(iparam);

      shp_jump_visc.velx.dxdz.n(iparam) = visc_minus*shp_minus.velx.dxdz.n(iparam) - visc_plus*shp_plus.velx.dxdz.n(iparam);
      shp_jump_visc.vely.dxdz.n(iparam) = visc_minus*shp_minus.vely.dxdz.n(iparam) - visc_plus*shp_plus.vely.dxdz.n(iparam);
      shp_jump_visc.velz.dxdz.n(iparam) = visc_minus*shp_minus.velz.dxdz.n(iparam) - visc_plus*shp_plus.velz.dxdz.n(iparam);

      shp_jump_visc.velx.dydx.n(iparam) = shp_jump_visc.velx.dxdy.n(iparam);
      shp_jump_visc.vely.dydx.n(iparam) = shp_jump_visc.vely.dxdy.n(iparam);
      shp_jump_visc.velz.dydx.n(iparam) = shp_jump_visc.velz.dxdy.n(iparam);

      shp_jump_visc.velx.dydy.n(iparam) = visc_minus*shp_minus.velx.dydy.n(iparam) - visc_plus*shp_plus.velx.dydy.n(iparam);
      shp_jump_visc.vely.dydy.n(iparam) = visc_minus*shp_minus.vely.dydy.n(iparam) - visc_plus*shp_plus.vely.dydy.n(iparam);
      shp_jump_visc.velz.dydy.n(iparam) = visc_minus*shp_minus.velz.dydy.n(iparam) - visc_plus*shp_plus.velz.dydy.n(iparam);

      shp_jump_visc.velx.dydz.n(iparam) = visc_minus*shp_minus.velx.dydz.n(iparam) - visc_plus*shp_plus.velx.dydz.n(iparam);
      shp_jump_visc.vely.dydz.n(iparam) = visc_minus*shp_minus.vely.dydz.n(iparam) - visc_plus*shp_plus.vely.dydz.n(iparam);
      shp_jump_visc.velz.dydz.n(iparam) = visc_minus*shp_minus.velz.dydz.n(iparam) - visc_plus*shp_plus.velz.dydz.n(iparam);

      shp_jump_visc.velx.dzdx.n(iparam) = shp_jump_visc.velx.dxdz.n(iparam);
      shp_jump_visc.vely.dzdx.n(iparam) = shp_jump_visc.vely.dxdz.n(iparam);
      shp_jump_visc.velz.dzdx.n(iparam) = shp_jump_visc.velz.dxdz.n(iparam);

      shp_jump_visc.velx.dzdy.n(iparam) = shp_jump_visc.velx.dydz.n(iparam);
      shp_jump_visc.vely.dzdy.n(iparam) = shp_jump_visc.vely.dydz.n(iparam);
      shp_jump_visc.velz.dzdy.n(iparam) = shp_jump_visc.velz.dydz.n(iparam);

      shp_jump_visc.velx.dzdz.n(iparam) = visc_minus*shp_minus.velx.dzdz.n(iparam) - visc_plus*shp_plus.velx.dzdz.n(iparam);
      shp_jump_visc.vely.dzdz.n(iparam) = visc_minus*shp_minus.vely.dzdz.n(iparam) - visc_plus*shp_plus.vely.dzdz.n(iparam);
      shp_jump_visc.velz.dzdz.n(iparam) = visc_minus*shp_minus.velz.dzdz.n(iparam) - visc_plus*shp_plus.velz.dzdz.n(iparam);
    }
    //---------------------------------------
    // compute data at Gaussian point for rhs
    //---------------------------------------

    // velocity jump
    static LINALG::Matrix<3,1> veljump(true);
    veljump.Clear();
    veljump = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, shp_vel_jump, numparamvelx, ele, dofman);

    // mean velocity
    static LINALG::Matrix<3,1> velmean(true);
    velmean.Clear();
    velmean = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, shp_vel_mean, numparamvelx, ele, dofman);

    // get velocity (np,i) derivatives at integration point
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    static LINALG::Matrix<3,3> vderxy_mean_visc(true);
    vderxy_mean_visc.Clear();
    {
      const int* nodeids = ele->NodeIds();
      std::size_t velncounter = 0;

      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (size_t inode=0; inode<numparamvelx; ++inode)
      {
        vderxy_mean_visc(0,0) += evelnp(0,inode) * shp_mean_visc.velx.dx.s(inode);
        vderxy_mean_visc(0,1) += evelnp(0,inode) * shp_mean_visc.velx.dy.s(inode);
        vderxy_mean_visc(0,2) += evelnp(0,inode) * shp_mean_visc.velx.dz.s(inode);

        vderxy_mean_visc(1,0) += evelnp(1,inode) * shp_mean_visc.vely.dx.s(inode);
        vderxy_mean_visc(1,1) += evelnp(1,inode) * shp_mean_visc.vely.dy.s(inode);
        vderxy_mean_visc(1,2) += evelnp(1,inode) * shp_mean_visc.vely.dz.s(inode);

        vderxy_mean_visc(2,0) += evelnp(2,inode) * shp_mean_visc.velz.dx.s(inode);
        vderxy_mean_visc(2,1) += evelnp(2,inode) * shp_mean_visc.velz.dy.s(inode);
        vderxy_mean_visc(2,2) += evelnp(2,inode) * shp_mean_visc.velz.dz.s(inode);

        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            vderxy_mean_visc(0,0) += evelnp(3,velncounter) * shp_mean_visc.velx.dx.n(velncounter);
            vderxy_mean_visc(0,1) += evelnp(3,velncounter) * shp_mean_visc.velx.dy.n(velncounter);
            vderxy_mean_visc(0,2) += evelnp(3,velncounter) * shp_mean_visc.velx.dz.n(velncounter);

            vderxy_mean_visc(1,0) += evelnp(3,velncounter) * shp_mean_visc.vely.dx.n(velncounter);
            vderxy_mean_visc(1,1) += evelnp(3,velncounter) * shp_mean_visc.vely.dy.n(velncounter);
            vderxy_mean_visc(1,2) += evelnp(3,velncounter) * shp_mean_visc.vely.dz.n(velncounter);

            vderxy_mean_visc(2,0) += evelnp(3,velncounter) * shp_mean_visc.velz.dx.n(velncounter);
            vderxy_mean_visc(2,1) += evelnp(3,velncounter) * shp_mean_visc.velz.dy.n(velncounter);
            vderxy_mean_visc(2,2) += evelnp(3,velncounter) * shp_mean_visc.velz.dz.n(velncounter);

            velncounter += 1;
          }
        }
      }
      if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
    }
    //cout << "veljump " << veljump << endl;
    //cout << "vderxy_mean_visc " << vderxy_mean_visc << endl;

    // get velocity (np,i) derivatives at integration point
    // vderxy = enr_derxy(j,k)*evelnp(i,k);
    static LINALG::Matrix<3,3> vderxy_jump_visc(true);
    vderxy_jump_visc.Clear();
    {
      const int* nodeids = ele->NodeIds();
      std::size_t velncounter = 0;

      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      for (size_t inode=0; inode<numparamvelx; ++inode)
      {
        vderxy_jump_visc(0,0) += evelnp(0,inode) * shp_jump_visc.velx.dx.s(inode);
        vderxy_jump_visc(0,1) += evelnp(0,inode) * shp_jump_visc.velx.dy.s(inode);
        vderxy_jump_visc(0,2) += evelnp(0,inode) * shp_jump_visc.velx.dz.s(inode);

        vderxy_jump_visc(1,0) += evelnp(1,inode) * shp_jump_visc.vely.dx.s(inode);
        vderxy_jump_visc(1,1) += evelnp(1,inode) * shp_jump_visc.vely.dy.s(inode);
        vderxy_jump_visc(1,2) += evelnp(1,inode) * shp_jump_visc.vely.dz.s(inode);

        vderxy_jump_visc(2,0) += evelnp(2,inode) * shp_jump_visc.velz.dx.s(inode);
        vderxy_jump_visc(2,1) += evelnp(2,inode) * shp_jump_visc.velz.dy.s(inode);
        vderxy_jump_visc(2,2) += evelnp(2,inode) * shp_jump_visc.velz.dz.s(inode);

        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            vderxy_jump_visc(0,0) += evelnp(3,velncounter) * shp_jump_visc.velx.dx.n(velncounter);
            vderxy_jump_visc(0,1) += evelnp(3,velncounter) * shp_jump_visc.velx.dy.n(velncounter);
            vderxy_jump_visc(0,2) += evelnp(3,velncounter) * shp_jump_visc.velx.dz.n(velncounter);

            vderxy_jump_visc(1,0) += evelnp(3,velncounter) * shp_jump_visc.vely.dx.n(velncounter);
            vderxy_jump_visc(1,1) += evelnp(3,velncounter) * shp_jump_visc.vely.dy.n(velncounter);
            vderxy_jump_visc(1,2) += evelnp(3,velncounter) * shp_jump_visc.vely.dz.n(velncounter);

            vderxy_jump_visc(2,0) += evelnp(3,velncounter) * shp_jump_visc.velz.dx.n(velncounter);
            vderxy_jump_visc(2,1) += evelnp(3,velncounter) * shp_jump_visc.velz.dy.n(velncounter);
            vderxy_jump_visc(2,2) += evelnp(3,velncounter) * shp_jump_visc.velz.dz.n(velncounter);

            velncounter += 1;
          }
        }
      }
      if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
    }

    // get pressure jump
    static double pjump;
    pjump = 0.0;
    for (size_t iparam = 0; iparam != numparampres; ++iparam)
      pjump += shp_pres_jump.d0(iparam)*eprenp(iparam);

    static double pmean;
    pmean = 0.0;
    for (size_t iparam = 0; iparam != numparampres; ++iparam)
      pmean += shp_pres_mean.d0(iparam)*eprenp(iparam);

    //--------------------------------------------------------------------------------------------
    // build single boundary integral stiffness terms
    //
    // remarks: || x || stands for the jump operator ( || x || = x^- - x^+ )
    //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
    //             n    stands for the normal vector on the interface pointing from the positive (+)
    //                  domain into the negative (-) domain (n = n^+ = -n^-)
    //--------------------------------------------------------------------------------------------
#ifdef COMBUST_NITSCHE_VISCCONSISTMEAN
    //-------------------------    |                                      |
    // viscous consistency term  + |  || v || , { 2\mu epsilon( Du ) } n  |
    //-------------------------    |                                      |

    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, 2.0*timefacfac*normal(0), shp_mean_visc.velx.dx);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0,     timefacfac*normal(1), shp_mean_visc.velx.dy);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0,     timefacfac*normal(2), shp_mean_visc.velx.dz);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_vel_jump.velx.d0,     timefacfac*normal(1), shp_mean_visc.vely.dx);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_vel_jump.velx.d0,     timefacfac*normal(2), shp_mean_visc.velz.dx);

    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_vel_jump.vely.d0,     timefacfac*normal(0), shp_mean_visc.velx.dy);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0,     timefacfac*normal(0), shp_mean_visc.vely.dx);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, 2.0*timefacfac*normal(1), shp_mean_visc.vely.dy);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0,     timefacfac*normal(2), shp_mean_visc.vely.dz);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_vel_jump.vely.d0,     timefacfac*normal(2), shp_mean_visc.velz.dy);

    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_vel_jump.velz.d0,     timefacfac*normal(0), shp_mean_visc.velx.dz);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_vel_jump.velz.d0,     timefacfac*normal(1), shp_mean_visc.vely.dz);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0,     timefacfac*normal(0), shp_mean_visc.velz.dx);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0,     timefacfac*normal(1), shp_mean_visc.velz.dy);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, 2.0*timefacfac*normal(2), shp_mean_visc.velz.dz);

    //   |                                       |
    // - |  || v || , { 2\mu epsilon( u_i ) } n  |
    //   |                                       |

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -timefacfac* (vderxy_mean_visc(0,0) + vderxy_mean_visc(0,0))*normal(0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -timefacfac* (vderxy_mean_visc(0,1) + vderxy_mean_visc(1,0))*normal(1));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -timefacfac* (vderxy_mean_visc(0,2) + vderxy_mean_visc(2,0))*normal(2));

    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -timefacfac* (vderxy_mean_visc(1,0) + vderxy_mean_visc(0,1))*normal(0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -timefacfac* (vderxy_mean_visc(1,1) + vderxy_mean_visc(1,1))*normal(1));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -timefacfac* (vderxy_mean_visc(1,2) + vderxy_mean_visc(2,1))*normal(2));

    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -timefacfac* (vderxy_mean_visc(2,0) + vderxy_mean_visc(0,2))*normal(0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -timefacfac* (vderxy_mean_visc(2,1) + vderxy_mean_visc(1,2))*normal(1));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -timefacfac* (vderxy_mean_visc(2,2) + vderxy_mean_visc(2,2))*normal(2));
#endif

#ifdef COMBUST_NITSCHE_VISCCONSISTJUMP
    //-------------------------    |                                      |
    // viscous consistency term  + |  < v > , || 2\mu epsilon( Du ) || n  |
    //-------------------------    |                                      |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_mean.d0, 2.0*timefacfac*normal(0), shp_jump_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_mean.d0,     timefacfac*normal(1), shp_jump_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_mean.d0,     timefacfac*normal(2), shp_jump_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_vel_mean.d0,     timefacfac*normal(1), shp_jump_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_vel_mean.d0,     timefacfac*normal(2), shp_jump_visc.dx);

    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_vel_mean.d0,     timefacfac*normal(0), shp_jump_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_mean.d0,     timefacfac*normal(0), shp_jump_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_mean.d0, 2.0*timefacfac*normal(1), shp_jump_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_mean.d0,     timefacfac*normal(2), shp_jump_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_vel_mean.d0,     timefacfac*normal(2), shp_jump_visc.dy);

    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_vel_mean.d0,     timefacfac*normal(0), shp_jump_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_vel_mean.d0,     timefacfac*normal(1), shp_jump_visc.dz);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_mean.d0,     timefacfac*normal(0), shp_jump_visc.dx);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_mean.d0,     timefacfac*normal(1), shp_jump_visc.dy);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_mean.d0, 2.0*timefacfac*normal(2), shp_jump_visc.dz);

    //   |                                       |
    // - |  < v > , || 2\mu epsilon( u_i ) || n  |
    //   |                                       |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(0,0) + vderxy_jump_visc(0,0))*normal(0));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(0,1) + vderxy_jump_visc(1,0))*normal(1));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(0,2) + vderxy_jump_visc(2,0))*normal(2));

    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(1,0) + vderxy_jump_visc(0,1))*normal(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(1,1) + vderxy_jump_visc(1,1))*normal(1));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(1,2) + vderxy_jump_visc(2,1))*normal(2));

    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(2,0) + vderxy_jump_visc(0,2))*normal(0));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(2,1) + vderxy_jump_visc(1,2))*normal(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_vel_mean.d0, -timefacfac* (vderxy_jump_visc(2,2) + vderxy_jump_visc(2,2))*normal(2));
#endif

#ifdef COMBUST_NITSCHE_VISCADJCONSIST
    //---------------------------------    |                                     |
    // viscous adjoint consistency term  + |  { 2\mu epsilon( v ) } n, || Du ||  |
    //---------------------------------    |                                     |

    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dx, 2.0*timefacfac*normal(0), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dy,     timefacfac*normal(1), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dz,     timefacfac*normal(2), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_mean_visc.vely.dx,     timefacfac*normal(1), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_mean_visc.velz.dx,     timefacfac*normal(2), shp_vel_jump.velx.d0);

    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_mean_visc.velx.dy,     timefacfac*normal(0), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dx,     timefacfac*normal(0), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dy, 2.0*timefacfac*normal(1), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dz,     timefacfac*normal(2), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_mean_visc.velz.dy,     timefacfac*normal(2), shp_vel_jump.vely.d0);

    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_mean_visc.velx.dz,     timefacfac*normal(0), shp_vel_jump.velz.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_mean_visc.vely.dz,     timefacfac*normal(1), shp_vel_jump.velz.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dx,     timefacfac*normal(0), shp_vel_jump.velz.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dy,     timefacfac*normal(1), shp_vel_jump.velz.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dz, 2.0*timefacfac*normal(2), shp_vel_jump.velz.d0);

    //   |                                      |
    // - |  { 2\mu epsilon( v ) } n, || u_i ||  |
    //   |                                      |

    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, -2.0* normal(0)           *timefacfac*veljump(0,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy, -    (normal(0)+normal(1))*timefacfac*veljump(0,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz, -    (normal(0)+normal(2))*timefacfac*veljump(0,0));

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dx, -2.0*timefacfac*normal(0)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dy, -    timefacfac*normal(1)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dz, -    timefacfac*normal(2)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dx, -    timefacfac*normal(1)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dx, -    timefacfac*normal(2)*veljump(0,0));

    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx, -    (normal(1)+normal(0))*timefacfac*veljump(1,0));
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, -2.0* normal(1)           *timefacfac*veljump(1,0));
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz, -    (normal(1)+normal(2))*timefacfac*veljump(1,0));

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dy, -    timefacfac*normal(0)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dx, -    timefacfac*normal(0)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dy, -2.0*timefacfac*normal(1)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dz, -    timefacfac*normal(2)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dy, -    timefacfac*normal(2)*veljump(1,0));

    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx, -    (normal(2)+normal(0))*timefacfac*veljump(2,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy, -    (normal(2)+normal(1))*timefacfac*veljump(2,0));
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, -2.0* normal(2)           *timefacfac*veljump(2,0));

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dz, -    timefacfac*normal(0)*veljump(2,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dz, -    timefacfac*normal(1)*veljump(2,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dx, -    timefacfac*normal(0)*veljump(2,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dy, -    timefacfac*normal(1)*veljump(2,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dz, -2.0*timefacfac*normal(2)*veljump(2,0));

    //-------------------------------------    |                                    |
    // viscous adjoint consistency term RHS  + |  { 2\mu epsilon( v ) } n, j_u  |
    //-------------------------------------    |                                    |
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, 2.0* normal(0)           *timefacfac*normal(0)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     (normal(0)+normal(1))*timefacfac*normal(0)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     (normal(0)+normal(2))*timefacfac*normal(0)*ju);

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dx, 2.0*timefacfac*normal(0)*ju(0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dy,     timefacfac*normal(1)*ju(0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dz,     timefacfac*normal(2)*ju(0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dx,     timefacfac*normal(1)*ju(0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dx,     timefacfac*normal(2)*ju(0));
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(1)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dz,     timefacfac*normal(2)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(1)*normal(0)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     timefacfac*normal(2)*normal(0)*ju_n);

    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     (normal(1)+normal(0))*timefacfac*normal(1)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, 2.0* normal(1)           *timefacfac*normal(1)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     (normal(1)+normal(2))*timefacfac*normal(1)*ju);

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dy,     timefacfac*normal(0)*ju(1));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dx,     timefacfac*normal(0)*ju(1));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dy, 2.0*timefacfac*normal(1)*ju(1));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dz,     timefacfac*normal(2)*ju(1));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dy,     timefacfac*normal(2)*ju(1));
    //assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean_visc.dy,     timefacfac*normal(0)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dx,     timefacfac*normal(0)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean_visc.dz,     timefacfac*normal(2)*normal(1)*ju_n);
    //assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     timefacfac*normal(2)*normal(1)*ju_n);

    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dx,     (normal(2)+normal(0))*timefacfac*normal(2)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dy,     (normal(2)+normal(1))*timefacfac*normal(2)*ju);
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean_visc.dz, 2.0* normal(2)           *timefacfac*normal(2)*ju);

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_mean_visc.velx.dz,     timefacfac*normal(0)*ju(2));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_mean_visc.vely.dz,     timefacfac*normal(1)*ju(2));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dx,     timefacfac*normal(0)*ju(2));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dy,     timefacfac*normal(1)*ju(2));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_mean_visc.velz.dz, 2.0*timefacfac*normal(2)*ju(2));
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

    assembler.template MatrixNS<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Pres>(shp_vel_jump.velx.d0, -timefacfac*normal(0), shp_pres_mean.d0);
    assembler.template MatrixNS<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Pres>(shp_vel_jump.vely.d0, -timefacfac*normal(1), shp_pres_mean.d0);
    assembler.template MatrixNS<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Pres>(shp_vel_jump.velz.d0, -timefacfac*normal(2), shp_pres_mean.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Pres>(shp_vel_jump.d0, -timefacfac*normal(0), shp_presproj_mean);
    //      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Pres>(shp_vel_jump.d0, -timefacfac*normal(1), shp_presproj_mean);
    //      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Pres>(shp_vel_jump.d0, -timefacfac*normal(2), shp_presproj_mean);

    //   |                      |
    // + |  || v ||, { p_i } n  |
    //   |                      |

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, timefacfac*normal(0)*pmean);
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, timefacfac*normal(1)*pmean);
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, timefacfac*normal(2)*pmean);
#endif

#ifdef COMBUST_NITSCHE_PRESCONSISTJUMP
    //--------------------------    |                     |
    // pressure consistency term  - |  < v >, || Dp || n  |
    //--------------------------    |                     |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Pres>(shp_vel_mean.d0, -timefacfac*normal(0), shp_pres_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Pres>(shp_vel_mean.d0, -timefacfac*normal(1), shp_pres_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Pres>(shp_vel_mean.d0, -timefacfac*normal(2), shp_pres_jump.d0);

    //   |                      |
    // + |  < v >, || p_i || n  |
    //   |                      |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_vel_mean.d0, timefacfac*normal(0)*pjump);
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_vel_mean.d0, timefacfac*normal(1)*pjump);
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_vel_mean.d0, timefacfac*normal(2)*pjump);
#endif

#ifdef COMBUST_NITSCHE_PRESADJCONSIST
    //----------------------------------    |                     |
    // pressure adjoint consistency term  + |  { q }, || Du || n  |
    //----------------------------------    |                     |

    assembler.template MatrixSN<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velx>(shp_pres_mean.d0, timefacfac*normal(0), shp_vel_jump.velx.d0);
    assembler.template MatrixSN<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Vely>(shp_pres_mean.d0, timefacfac*normal(1), shp_vel_jump.vely.d0);
    assembler.template MatrixSN<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velz>(shp_pres_mean.d0, timefacfac*normal(2), shp_vel_jump.velz.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velx>(shp_presproj_mean, +timefacfac*normal(0), shp_vel_jump.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Vely>(shp_presproj_mean, +timefacfac*normal(1), shp_vel_jump.d0);
    //      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velz>(shp_presproj_mean, +timefacfac*normal(2), shp_vel_jump.d0);

    //   |                      |
    // - |  { q }, || u_i || n  |
    //   |                      |

    assembler.template Vector<XFEM::PHYSICS::Pres>(shp_pres_mean.d0, -timefacfac*(veljump(0,0)*normal(0)
        +veljump(1,0)*normal(1)
        +veljump(2,0)*normal(2)));
    //      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_presproj_mean, -timefacfac*(vjump(0,0)*normal(0)
    //          +vjump(1,0)*normal(1)
    //          +vjump(2,0)*normal(2)));

    //--------------------------------------    |                |
    // pressure adjoint consistency term RHS  + |  { q }, j_u n  |
    //--------------------------------------    |                |

    assembler.template Vector<XFEM::PHYSICS::Pres>(shp_pres_mean.d0, timefacfac*(ju(0)*normal(0)
                                                                 +ju(1)*normal(1)
                                                                 +ju(2)*normal(2)));
    //assembler.template Vector<XFEM::PHYSICS::Pres>(shp_pres_mean.d0, timefacfac*ju_n);

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
dserror("uncomment this!");
//      COMBUST::RHS_SurfaceTension2DIntegrals<DISTYPE,ASSTYPE,NUMDOF>(posXiDomain, shp_vel_mean_inv, ele,
//          xyze, egradphi, normal, timefacfac, surftensapprox, surftenscoeff, assembler);
    }
    else
    {
      //---------------------------------------------------    |                 |
      // flux (traction/stress) jump (consistency) term RHS  - |  < v >, j_flux  |
      //---------------------------------------------------    |                 |

      assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_mean_inv.velx.d0, -timefacfac*jflux(0));
      assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_mean_inv.vely.d0, -timefacfac*jflux(1));
      assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_mean_inv.velz.d0, -timefacfac*jflux(2));

#if 0
      assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_mean.velx.d0, -timefacfac* (vderxy_jump_visc(0,0) + vderxy_jump_visc(0,0))*normal(0));
      assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_mean.velx.d0, -timefacfac* (vderxy_jump_visc(0,1) + vderxy_jump_visc(1,0))*normal(1));
      assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_mean.velx.d0, -timefacfac* (vderxy_jump_visc(0,2) + vderxy_jump_visc(2,0))*normal(2));

      assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_mean.vely.d0, -timefacfac* (vderxy_jump_visc(1,0) + vderxy_jump_visc(0,1))*normal(0));
      assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_mean.vely.d0, -timefacfac* (vderxy_jump_visc(1,1) + vderxy_jump_visc(1,1))*normal(1));
      assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_mean.vely.d0, -timefacfac* (vderxy_jump_visc(1,2) + vderxy_jump_visc(2,1))*normal(2));

      assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_mean.velz.d0, -timefacfac* (vderxy_jump_visc(2,0) + vderxy_jump_visc(0,2))*normal(0));
      assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_mean.velz.d0, -timefacfac* (vderxy_jump_visc(2,1) + vderxy_jump_visc(1,2))*normal(1));
      assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_mean.velz.d0, -timefacfac* (vderxy_jump_visc(2,2) + vderxy_jump_visc(2,2))*normal(2));
#endif
    }
#endif

#ifdef COMBUST_NITSCHE_NITSCHEVEL
    //------------------------    |                              |
    // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
    //------------------------    |                              |

    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, alphau*timefacfac, shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, alphau*timefacfac, shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, alphau*timefacfac, shp_vel_jump.velz.d0);

    //    |                               |
    //  - |  \alpha_u || v ||, || u_i ||  |
    //    |                               |

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -alphau*timefacfac*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -alphau*timefacfac*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -alphau*timefacfac*veljump(2,0));

    //----------------------------    |                         |
    // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, j_u  |
    //----------------------------    |                         |

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, alphau*timefacfac*ju(0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, alphau*timefacfac*ju(1));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, alphau*timefacfac*ju(2));
    //assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, alphau*timefacfac*normal(0)*ju_n);
    //assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, alphau*timefacfac*normal(1)*ju_n);
    //assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, alphau*timefacfac*normal(2)*ju_n);
#endif

#ifdef COMBUST_NITSCHE_CONVSTAB
    //------------------------    |                                       |
    // Nitsche term (velocity)  + |  \alpha_u_conv || v ||*n, || Du ||*n  |
    //------------------------    |                                       |

    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, alpha_conv*normal(0)*timefacfac*normal(0), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Vely>(shp_vel_jump.velx.d0, alpha_conv*normal(0)*timefacfac*normal(1), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velz>(shp_vel_jump.velx.d0, alpha_conv*normal(0)*timefacfac*normal(2), shp_vel_jump.velz.d0);

    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velx>(shp_vel_jump.vely.d0, alpha_conv*normal(1)*timefacfac*normal(0), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, alpha_conv*normal(1)*timefacfac*normal(1), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Velz>(shp_vel_jump.vely.d0, alpha_conv*normal(1)*timefacfac*normal(2), shp_vel_jump.velz.d0);

    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velx>(shp_vel_jump.velz.d0, alpha_conv*normal(2)*timefacfac*normal(0), shp_vel_jump.velx.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Vely>(shp_vel_jump.velz.d0, alpha_conv*normal(2)*timefacfac*normal(1), shp_vel_jump.vely.d0);
    assembler.template MatrixNN<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, alpha_conv*normal(2)*timefacfac*normal(2), shp_vel_jump.velz.d0);


    //    |                                        |
    //  - |  \alpha_u_conv || v ||*n, || u_i ||*n  |
    //    |                                        |

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -alpha_conv*normal(0)*timefacfac*normal(0)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -alpha_conv*normal(0)*timefacfac*normal(1)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, -alpha_conv*normal(0)*timefacfac*normal(2)*veljump(2,0));

    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -alpha_conv*normal(1)*timefacfac*normal(0)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -alpha_conv*normal(1)*timefacfac*normal(1)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, -alpha_conv*normal(1)*timefacfac*normal(2)*veljump(2,0));

    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -alpha_conv*normal(2)*timefacfac*normal(0)*veljump(0,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -alpha_conv*normal(2)*timefacfac*normal(1)*veljump(1,0));
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, -alpha_conv*normal(2)*timefacfac*normal(2)*veljump(2,0));

    //----------------------------    |                                  |
    // Nitsche term (velocity) RHS  + |  \alpha_u_conv || v ||*n, J_u*n  |
    //----------------------------    |                                  |

    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, +alpha_conv*normal(0)*timefacfac*normal(0)*normal(0)*ju_n);
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, +alpha_conv*normal(0)*timefacfac*normal(1)*normal(1)*ju_n);
    assembler.template VectorN<XFEM::PHYSICS::Velx>(shp_vel_jump.velx.d0, +alpha_conv*normal(0)*timefacfac*normal(2)*normal(2)*ju_n);

    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, +alpha_conv*normal(1)*timefacfac*normal(0)*normal(0)*ju_n);
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, +alpha_conv*normal(1)*timefacfac*normal(1)*normal(1)*ju_n);
    assembler.template VectorN<XFEM::PHYSICS::Vely>(shp_vel_jump.vely.d0, +alpha_conv*normal(1)*timefacfac*normal(2)*normal(2)*ju_n);

    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, +alpha_conv*normal(2)*timefacfac*normal(0)*normal(0)*ju_n);
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, +alpha_conv*normal(2)*timefacfac*normal(1)*normal(1)*ju_n);
    assembler.template VectorN<XFEM::PHYSICS::Velz>(shp_vel_jump.velz.d0, +alpha_conv*normal(2)*timefacfac*normal(2)*normal(2)*ju_n);

#endif

#ifdef COMBUST_NITSCHE_NITSCHEPRES
    //------------------------    |                              |
    // Nitsche term (pressure)  + |  \alpha_p || q ||, || Dp ||  |
    //------------------------    |                              |

    assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Pres>(shp_pres_jump.d0, alphap*timefacfac, shp_pres_jump.d0);

    //    |                               |
    //  - |  \alpha_p || q ||, || p_i ||  |
    //    |                               |

    assembler.template Vector<XFEM::PHYSICS::Pres>(shp_pres_jump.d0, -alphap*timefacfac*pjump);

    //----------------------------    |                         |
    // Nitsche term (pressure) RHS  + |  \alpha_p || q ||, J_p  |
    //----------------------------    |                         |

    assembler.template Vector<XFEM::PHYSICS::Pres>(shp_pres_jump.d0, alphap*timefacfac*jp);
#endif
  } // loop Gaussian points

  return;
}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2>
void SysmatBoundaryNitscheNormal(
    const DRT::ELEMENTS::Combust3*             ele,            ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*&    ih,  ///< connection to the interface handler
    const XFEM::ElementDofManager&             dofman,         ///< dofmanager of the current element
    const M1&                                  evelnp,         ///<
    const V1&                                  eprenp,         ///<
    const V2&                                  ephi,           ///<
    const M2&                                  egradphi,
    Teuchos::RCP<const MAT::Material>          material,       ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme  timealgo,       ///< time discretization type
    const double&                              dt,             ///< delta t (time step size)
    const double&                              theta,          ///< factor for one step theta scheme
    const double                               ga_alphaF,
    const double                               ga_alphaM,
    const double                               ga_gamma,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler,      ///<
    const double                               flamespeed,     ///<
    const double                               nitschevel,     ///<
    const double                               nitschepres,    ///<
    const double                               ele_meas_plus,  ///< measure of element in plus domain
    const double                               ele_meas_minus, ///< measure of element in minus domain
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType     veljumptype,
    const INPAR::COMBUST::FluxJumpType         fluxjumptype,
    const bool                                 smoothed_boundary_integration
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - boundary");

#ifdef DEBUG
  if (ASSTYPE != XFEM::xfem_assembly)
    dserror("boundary integration integrals must only be added to intersected XFEM elements");
#endif

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration constant
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  //---------------------------------------------------------------------------
  // get material parameters for all boundary integration cells of this element
  //---------------------------------------------------------------------------
  // density \rho
  double dens_plus = -1.0;
  //cout << "plus density: " << dens_plus << endl;
  // dynamic viscosity \mu
  double dynvisc_plus = -1.0;
  // density \rho
  double dens_minus = -1.0;
  //cout << "minus density: " << dens_minus << endl;
  // dynamic viscosity \mu
  double dynvisc_minus = -1.0;
  // set density and dynamic viscosity
  COMBUST::GetMaterialParams(material, dens_plus, dynvisc_plus, dens_minus, dynvisc_minus);

  //-----------------------------------------------------------------------------------
  // get surface tension coefficient for all boundary integration cells of this element
  //-----------------------------------------------------------------------------------
  double surftenscoeff = 0.0;
  // set surface tension
  COMBUST::GetMaterialParams(material, surftenscoeff);

  //-----------------------------
  // prepare interface conditions
  //-----------------------------
  // mass flux through flame front
  const double mflux = flamespeed*dens_minus;

  // velocity jump in normal direction
#if 0
  // jump velocity vector ju_{n,t1,t2}
  LINALG::Matrix<nsd,1> ju_nt(true);
  ju_nt(0,0) = -mflux*(1.0/dens_minus - 1.0/dens_plus);
  ju_nt(1,0) = 0.0;
  ju_nt(2,0) = 0.0;
#else
  double ju_n = 0.0;
  if (veljumptype == INPAR::COMBUST::vel_jump_premixed_combustion){
    ju_n = -mflux*(1.0/dens_minus - 1.0/dens_plus);
  }
  else if(veljumptype == INPAR::COMBUST::vel_jump_const){
    ju_n = 5.0;
  }
  else{
    // velocity jump is zero
  }
#endif

  // jump in flux (normal stress)
  double jflux_n = -mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);
  if(fluxjumptype == INPAR::COMBUST::flux_jump_premixed_combustion){
    jflux_n = mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);
  }
  else if(fluxjumptype == INPAR::COMBUST::flux_jump_const){
    jflux_n = 5.0;
  }
  else{
    // flux jump is zero
  }

  //----------------------------------------
  // compute Nitsche stabilization patameter
  //----------------------------------------
  // calculate element diameter
  const double hk_eleDiam = COMBUST::getEleDiameter<DISTYPE>(xyze);

  // get maximum of the viscosities of all fluid materials
  const double dynvisc_max = std::max(dynvisc_plus,dynvisc_minus);

  // Nitsche parameter velocity ( alph_u := 1/h_k * mu_max * gamma_u ) see Becker, Hansbo 2009
  const double alphau = 1.0/hk_eleDiam * dynvisc_max * nitschevel;
  //cout << "alpha_u: " << alphau << endl;

  // Benedikt:
  // the input parameter alphap is used for the stabilization term for convection dominated flows
  // new convective stabilization parameter
  const double alpha_conv = 1.0/hk_eleDiam * nitschepres;

  //-----------------------------------------------------------------------------------------
  // get weights for average operators { () } = kappa_plus * ()_plus + kappa_minus * ()_minus
  //                               and < () > = kappa_minus* ()_plus + kappa_plus  * ()_minus
  //-----------------------------------------------------------------------------------------

  double ele_meas = ele_meas_plus + ele_meas_minus;
  if(ele_meas == 0.0) dserror("element measure is zero!!! Weighted average operators can not get defined!");
  double kappa_plus  = ele_meas_plus/ele_meas;
  double kappa_minus = ele_meas_minus/ele_meas;

  if(ele->Touched() == true)
  {
    // weights are chosen inversly -> touched element as a limit case of an intersected element
    // with ele_meas_plus -> 0
    // remark: interface belongs to Omega(+)
    kappa_plus = 0.0;
    kappa_minus= 1.0;
  }

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->ElementBoundaryIntCells(ele->Id());

  //-------------------------------------------------------------------------------
  // evaluate the enrichment function at the interface (boundary integration cells)
  //-------------------------------------------------------------------------------
  const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
  const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);

  //-------------------------------------
  // loop over boundary integration cells
  //-------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    switch (cell->Shape())
    {

    // TODO GaussRule modified!
    case DRT::Element::tri3:
      COMBUST::Nitsche_SysmatBoundaryCellNormal<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
          *cell, DRT::UTILS::intrule_tri_6point, xyze, ele, dofman, higher_order_ele,
          evelnp, eprenp, ephi, egradphi,
          ju_n, jflux_n,
          timefac, assembler,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus,
          alphau, alpha_conv,
          surftensapprox, surftenscoeff, connected_interface, veljumptype, fluxjumptype,
          smoothed_boundary_integration
      );
      break;
    case DRT::Element::quad4:
      // TODO check this intrule !!!
      COMBUST::Nitsche_SysmatBoundaryCellNormal<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          *cell, DRT::UTILS::intrule_quad_9point , xyze, ele, dofman, higher_order_ele,
          evelnp, eprenp, ephi, egradphi,
          ju_n, jflux_n,
          timefac, assembler,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, dynvisc_plus, dynvisc_minus,
          alphau, alpha_conv,
          surftensapprox, surftenscoeff, connected_interface, veljumptype, fluxjumptype,
          smoothed_boundary_integration
      );
      break;
    default:
      dserror("invalid type of boundary integration cell");
    }
  } // loop boundary integration cells

  return;
}

} // namespace COMBUST

#endif

