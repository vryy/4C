/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_twophaseflow.H

\brief system matrix formulation for two-phase flow problem / 3D combustion element

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_SYSMAT_TWOPHASEFLOW_H
#define COMBUST3_SYSMAT_TWOPHASEFLOW_H


#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat.H"
#include "combust3_local_assembler.H"
#include "combust3_utils.H"
#include "combust3_interpolation.H"
#include "combust3_stabilization.H"
#include "two_phase_defines.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "combust_time_integration_element.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_xfem/enrichment.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_geometry/position_array.H"
#include "../drt_io/io_gmsh.H"
//#include "two_phase_defines.H"

namespace COMBUST
{

/*!
  Calculate integrals in matrix and rhs for stationary two-phase flow problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class M2, class V2>
void SysmatTwoPhaseFlow(
    const DRT::ELEMENTS::Combust3*      ele,           ///< the element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
    const M1&                           eveln,
    const M1&                           evelnm,
    const M1&                           eaccn,
    const M1&                           eaccam,       ///< nodal acceleration n values for this element
    const V1&                           eprenp,
    const V2&                           ephi,
    const M2&                           etensor,
    Teuchos::RCP<const MAT::Material>   material,      ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,      ///< time discretization type
    const double                        time,          ///< current time
    const double                        dt,            ///< delta t (time step size)
    const double                        theta,         ///< factor for one step theta scheme
    const double                        ga_alphaF,
    const double                        ga_alphaM,
    const double                        ga_gamma,
    //const double                        step,          ///< for gmsh output of subgrid velocity
    const bool                          newton,        ///< full Newton or fixed-point-like
    const bool                          pstab,         ///< flag for stabilization
    const bool                          supg,          ///< flag for stabilization
    const bool                          graddiv,         ///< flag for stabilization
    const INPAR::FLUID::TauType         tautype,       ///< stabilization parameter definition
    const bool                          instationary,  ///< switch between stationary and instationary formulation
    const bool                          genalpha,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluate - Sysmat - domain");

  // number of nodes for element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // space dimension for 3d fluid element
  const size_t nsd = 3;

  // time integration constant
  const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

  // get node coordinates of the current element
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

  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, Velx, numnode);
  // different enrichments for pressure and velocity possible
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  //    std::cout << "element " << ele->Bisected() << std::endl;
  //    std::cout << "velocity dofs " << numparamvelx << std::endl;
  //    std::cout << "pressure dofs " << numparampres << std::endl;

  // stabilization parameter
  const double hk = COMBUST::UTILS::HK<DISTYPE>(xyze);
  const double mk = COMBUST::UTILS::MK<DISTYPE>();

  // information about domain integration cells
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  //std::cout << "number of IntCells " <<  domainIntCells.size() << std::endl;
  // loop over integration cells
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    //--------------------------------------------------
    // get material parameters for this integration cell
    //--------------------------------------------------
    // density \rho^{n+1}
    double densaf  = -1.0;
    // dynamic viscosity \mu^{n+1}
    // remark: der Wert im input file wird hier als dynamische Viskosität interpretiert
    //         Alternative:
    //         const double kinvisc = mat->Viscosity();
    //         const double dynvisc = kinvisc * dens;
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

    // special getXFEMGaussruleKinkEnr for kink enrichment is called as parabolic shape functions are obtained
    // after multipying N and Psi
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussruleKinkEnr<DISTYPE>(ele, xyze, ele->Splited(),cell->Shape());

    // gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

    // integration loop over Gauss points
    for(int iquad=0; iquad<intpoints.nquad; ++iquad)
    {

      //          std::cout << "jetzt kommt Gausspunktschleife" << std::endl;
      //            // coordinates of the current integration point in cell coordinates \eta^domain
      //            const LINALG::Matrix<nsd,1> pos_eta_domain(intpoints.qxg[iquad]);
      // coordinates of the current integration point in cell coordinates \eta
      LINALG::Matrix<nsd,1> pos_eta_domain;
      pos_eta_domain(0) = intpoints.qxg[iquad][0];
      pos_eta_domain(1) = intpoints.qxg[iquad][1];
      pos_eta_domain(2) = intpoints.qxg[iquad][2];

      // coordinates of the current integration point in element coordinates \xi
      LINALG::Matrix<nsd,1> posXiDomain;
//      GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
//      const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
      // this is a workaround for TwoPhaseFlow with standard FEM
      // sets AssemblyType to XFEM Assembly, otherwise detEtaToXi3D is always 1.0
      // and coordinate mapping is incorrect!
      double detcell_temp=0.0;
      if(ASSTYPE==XFEM::standard_assembly and ele->Splited())
      {
        GEO::mapEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain, posXiDomain);
        detcell_temp = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain);
      }
      else
      {
        GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
        detcell_temp = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
      }
      const double detcell=detcell_temp;

      // shape functions and their first derivatives
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

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

      // compute global derivates
      static LINALG::Matrix<3,numnode> derxy;

      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);
      //TEST
      //            for (std::size_t no=0; no<numnode; no++)
      //            {
      //             std::cout << derxy(0,no) << "    " << derxy(1,no) << "    " << derxy(2,no) << std::endl;
      ////               std::cout << derxy(1,no) << std::endl;
      ////               std::cout << derxy(2,no) << std::endl;
      //            }
      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE,numnode>(xjm, derxy, deriv2, xyze, derxy2);
        //TEST
        //                for (std::size_t no=0; no<numnode; no++)
        //                {
        //                 std::cout << derxy2(0,no) << "   " << derxy2(1,no) << "   " << derxy2(2,no) << "   " << derxy2(3,no) << "   " << derxy2(4,no) << "   " << derxy2(5,no) << std::endl;
        ////                   std::cout << derxy2(1,no) << std::endl;
        ////                   std::cout << derxy2(2,no) << std::endl;
        ////                   std::cout << derxy2(3,no) << std::endl;
        ////                   std::cout << derxy2(5,no) << std::endl;
        ////                   std::cout << derxy2(5,no) << std::endl;
        //                }
      }
      else
      {
        derxy2.Clear();
      }

      // the enrichment functions depend on the gauss point
      // therefore the comutation of the enrichment functions is called here
      // the gauss point is contained in funct!


#if(0)
      // Original version of Ursula
      // here jump enrichments are called with level-set values
      const XFEM::ElementEnrichmentValues enrvals(
          *ele,
          dofman,
          ephi,
          funct,
          derxy,
          derxy2);
#else
      // schott Sep 20, 2010
      // new version of Benedikt
      // kink enrichments are called with level-set values
      // jump enrichments are called with domain cells!
      const XFEM::ElementEnrichmentValues enrvals(
          *ele,
          dofman,
          ephi,
          *cell,
          funct,
          derxy,
          derxy2);
#endif

      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;


      //static Shp<shpVecSize> shp;

      static XFEM::ApproxFunc<2,shpVecSize> shpvel;
      static XFEM::ApproxFunc<2,shpVecSize> shppres;

      //            typedef LINALG::Matrix<shpVecSize,1> ShpVec;
      //            static ShpVec shp;
      //            static ShpVec shp_dx;
      //            static ShpVec shp_dy;
      //            static ShpVec shp_dz;
      //            static ShpVec shp_dxdx;
      //            static ShpVec shp_dxdy;
      //            static ShpVec shp_dxdz;
      //            static ShpVec shp_dydx;
      //            static ShpVec shp_dydy;
      //            static ShpVec shp_dydz;
      //            static ShpVec shp_dzdx;
      //            static ShpVec shp_dzdy;
      //            static ShpVec shp_dzdz;

      static LINALG::Matrix<shpVecSize,1> enr_funct_vel;
      static LINALG::Matrix<3,shpVecSize> enr_derxy_vel;
      static LINALG::Matrix<6,shpVecSize> enr_derxy2_vel;

      static LINALG::Matrix<shpVecSize,1> enr_funct_pres;
      static LINALG::Matrix<3,shpVecSize> enr_derxy_pres;
      static LINALG::Matrix<6,shpVecSize> enr_derxy2_pres;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays
        //                static LINALG::Matrix<shpVecSize,1> enr_funct_vel;
        //                static LINALG::Matrix<3,shpVecSize> enr_derxy_vel;
        //                static LINALG::Matrix<6,shpVecSize> enr_derxy2_vel;


        // shape function for nodal dofs
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Velx,
            funct,
            derxy,
            derxy2,
            enr_funct_vel,
            enr_derxy_vel,
            enr_derxy2_vel);

        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shpvel.d0(iparam) = enr_funct_vel(iparam);
          shpvel.dx(iparam) = enr_derxy_vel(0,iparam);
          shpvel.dy(iparam) = enr_derxy_vel(1,iparam);
          shpvel.dz(iparam) = enr_derxy_vel(2,iparam);
          shpvel.dxdx(iparam) = enr_derxy2_vel(0,iparam);
          shpvel.dxdy(iparam) = enr_derxy2_vel(3,iparam);
          shpvel.dxdz(iparam) = enr_derxy2_vel(4,iparam);
          shpvel.dydx(iparam) = shpvel.dxdy(iparam);
          shpvel.dydy(iparam) = enr_derxy2_vel(1,iparam);
          shpvel.dydz(iparam) = enr_derxy2_vel(5,iparam);
          shpvel.dzdx(iparam) = shpvel.dxdz(iparam);
          shpvel.dzdy(iparam) = shpvel.dydz(iparam);
          shpvel.dzdz(iparam) = enr_derxy2_vel(2,iparam);
          //TEST
          //                  std::cout << "Geschwindigkeit " << std::endl;
          //                  std::cout << "Enr " << enr_funct_vel(iparam) << std::endl;
          //                  std::cout << "Enrderxx " << enr_derxy_vel(0,iparam) << std::endl;
          //                  std::cout << "Enrderyy " << enr_derxy_vel(1,iparam) << std::endl;
          //                  std::cout << "Enrderzz " << enr_derxy_vel(2,iparam) << std::endl;
          //                  std::cout << "Enrderxx2 " << enr_derxy2_vel(0,iparam) << std::endl;
          //                  std::cout << "Enrderyy2 " << enr_derxy2_vel(1,iparam) << std::endl;
          //                  std::cout << "Enrderzz2 " << enr_derxy2_vel(2,iparam) << std::endl;
          //                  std::cout << "Enrderxy2 " << enr_derxy2_vel(3,iparam) << std::endl;
          //                  std::cout << "Enrderxz2 " << enr_derxy2_vel(4,iparam) << std::endl;
          //                  std::cout << "Enrderyz2 " << enr_derxy2_vel(5,iparam) << std::endl;
        }

        //                static LINALG::Matrix<shpVecSize,1> enr_funct_pres;
        //                static LINALG::Matrix<3,shpVecSize> enr_derxy_pres;
        //                static LINALG::Matrix<6,shpVecSize> enr_derxy2_pres;
        // shape function for nodal dofs pressure
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Pres,
            funct,
            derxy,
            derxy2,
            enr_funct_pres,
            enr_derxy_pres,
            enr_derxy2_pres);
        for (size_t iparam = 0; iparam < numparampres; ++iparam)
        {
          shppres.d0(iparam) = enr_funct_pres(iparam);
          shppres.dx(iparam) = enr_derxy_pres(0,iparam);
          shppres.dy(iparam) = enr_derxy_pres(1,iparam);
          shppres.dz(iparam) = enr_derxy_pres(2,iparam);
          shppres.dxdx(iparam) = enr_derxy2_pres(0,iparam);
          shppres.dxdy(iparam) = enr_derxy2_pres(3,iparam);
          shppres.dxdz(iparam) = enr_derxy2_pres(4,iparam);
          shppres.dydx(iparam) = shppres.dxdy(iparam);
          shppres.dydy(iparam) = enr_derxy2_pres(1,iparam);
          shppres.dydz(iparam) = enr_derxy2_pres(5,iparam);
          shppres.dzdx(iparam) = shppres.dxdz(iparam);
          shppres.dzdy(iparam) = shppres.dydz(iparam);
          shppres.dzdz(iparam) = enr_derxy2_pres(2,iparam);

          //TEST
          //                  std::cout << "Druck " << std::endl;
          //                  std::cout << "Enr " << enr_funct_pres(iparam) << std::endl;
          //                  std::cout << "Enrderxx " << enr_derxy_pres(0,iparam) << std::endl;
          //                  std::cout << "Enrderyy " << enr_derxy_pres(1,iparam) << std::endl;
          //                  std::cout << "Enrderzz " << enr_derxy_pres(2,iparam) << std::endl;
          //                  std::cout << "Enrderxx2 " << enr_derxy2_pres(0,iparam) << std::endl;
          //                  std::cout << "Enrderyy2 " << enr_derxy2_pres(1,iparam) << std::endl;
          //                  std::cout << "Enrderzz2 " << enr_derxy2_pres(2,iparam) << std::endl;
          //                  std::cout << "Enrderxy2 " << enr_derxy2_pres(3,iparam) << std::endl;
          //                  std::cout << "Enrderxz2 " << enr_derxy2_pres(4,iparam) << std::endl;
          //                  std::cout << "Enrderyz2 " << enr_derxy2_pres(5,iparam) << std::endl;
        }
      }
      else // not xfem_assembly
      {
        //std::cout << "Shapefunction FEM" << std::endl;
        for (size_t iparam = 0; iparam < numnode; ++iparam)
        {
          shpvel.d0(iparam) = funct(iparam);
          shpvel.dx(iparam) = derxy(0,iparam);
          shpvel.dy(iparam) = derxy(1,iparam);
          shpvel.dz(iparam) = derxy(2,iparam);
          shpvel.dxdx(iparam) = derxy2(0,iparam);
          shpvel.dxdy(iparam) = derxy2(3,iparam);
          shpvel.dxdz(iparam) = derxy2(4,iparam);
          shpvel.dydx(iparam) = shpvel.dxdy(iparam);
          shpvel.dydy(iparam) = derxy2(1,iparam);
          shpvel.dydz(iparam) = derxy2(5,iparam);
          shpvel.dzdx(iparam) = shpvel.dxdz(iparam);
          shpvel.dzdy(iparam) = shpvel.dydz(iparam);
          shpvel.dzdz(iparam) = derxy2(2,iparam);

          shppres.d0(iparam) = funct(iparam);
          shppres.dx(iparam) = derxy(0,iparam);
          shppres.dy(iparam) = derxy(1,iparam);
          shppres.dz(iparam) = derxy(2,iparam);
          shppres.dxdx(iparam) = derxy2(0,iparam);
          shppres.dxdy(iparam) = derxy2(3,iparam);
          shppres.dxdz(iparam) = derxy2(4,iparam);
          shppres.dydx(iparam) = shppres.dxdy(iparam);
          shppres.dydy(iparam) = derxy2(1,iparam);
          shppres.dydz(iparam) = derxy2(5,iparam);
          shppres.dzdx(iparam) = shppres.dxdz(iparam);
          shppres.dzdy(iparam) = shppres.dydz(iparam);
          shppres.dzdz(iparam) = derxy2(2,iparam);
        }

      }

      // get velocities and accelerations at integration point
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shpvel.d0, numparamvelx);
      LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shpvel.d0, numparamvelx);
      LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shpvel.d0, numparamvelx);
      LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPoint(eaccn , shpvel.d0, numparamvelx);
      // get acceleration at n+alpha_m for generalized alpha time integration scheme
      LINALG::Matrix<3,1> gpaccam(true);
      if (genalpha)
        gpaccam  = XFEM::interpolateVectorFieldToIntPoint(eaccam , shpvel.d0, numparamvelx);

      // commenting this section out leads to problems with instationary calculations   henke 01/09
      //            if (ASSTYPE == XFEM::xfem_assembly and timealgo != timeint_stationary)
      //            {
      //              const bool valid_spacetime_cell_found = COMBUST::modifyOldTimeStepsValues<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, ivelcoln, ivelcolnm, iacccoln, gpveln, gpvelnm, gpaccn);
      //              if (not valid_spacetime_cell_found)
      //                continue;
      //            }
      //            std::cout << gpvelnp << std::endl;
      //            std::cout << evelnp << std::endl;
      //            std::cout << shp << std::endl;

      // get history data (n) at integration point
      //            LINALG::Matrix<3,1> histvec;
      //            //histvec = enr_funct(j)*evelnp_hist(i,j);
      //            for (int isd = 0; isd < nsd; ++isd)
      //            {
      //                histvec(isd) = 0.0;
      //                for (int iparam = 0; iparam < numparamvelx; ++iparam)
      //                    histvec(isd) += evelnp_hist(isd,iparam)*shp.d0(iparam);
      //            }
      LINALG::Matrix<nsd,1> histvec = COMBUST::TIMEINT::GetOldPartOfRighthandside(
          gpveln, gpvelnm, gpaccn, timealgo, dt, theta);
      histvec.Scale(densaf);

      // get velocity (np,i) derivatives at integration point
      // vderxy = enr_derxy(j,k)*evelnp(i,k);
      static LINALG::Matrix<3,nsd> vderxy;
      vderxy.Clear();
      for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
        {
          vderxy(isd,0) += evelnp(isd,iparam) * shpvel.dx(iparam);
          vderxy(isd,1) += evelnp(isd,iparam) * shpvel.dy(iparam);
          vderxy(isd,2) += evelnp(isd,iparam) * shpvel.dz(iparam);
        }
      }
      //compute divergence; test
      //if (ele->Id()==1216)
      //{
      //  const double vdiv = vderxy(0,0) + vderxy(1,1) + vderxy(2,2);
      //  std::cout << "element " << ele->Id() << "GP " << iquad << " " << vdiv << std::endl;
      //}
      //std::cout << "eps_xy" << (0.5*(vderxy(0,1)+vderxy(1,0))) << ", "<< std::endl;

      // calculate 2nd velocity derivatives at integration point
      static LINALG::Matrix<3,6> vderxy2;
      if (higher_order_ele)
      {
        //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
        vderxy2.Clear();
        for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          for (size_t isd = 0; isd < nsd; ++isd)
          {
            vderxy2(isd,0) += evelnp(isd,iparam)*shpvel.dxdx(iparam);
            vderxy2(isd,1) += evelnp(isd,iparam)*shpvel.dydy(iparam);
            vderxy2(isd,2) += evelnp(isd,iparam)*shpvel.dzdz(iparam);
            vderxy2(isd,3) += evelnp(isd,iparam)*shpvel.dxdy(iparam);
            vderxy2(isd,4) += evelnp(isd,iparam)*shpvel.dxdz(iparam);
            vderxy2(isd,5) += evelnp(isd,iparam)*shpvel.dydz(iparam);
          }
        }
      }
      else
      {
        vderxy2.Clear();
      }

      // get pressure gradients
      // gradp = enr_derxy(i,j)*eprenp(j);
      LINALG::Matrix<nsd,1> gradp(true);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += shppres.dx(iparam)*eprenp(iparam);
        gradp(1) += shppres.dy(iparam)*eprenp(iparam);
        gradp(2) += shppres.dz(iparam)*eprenp(iparam);
      }

      //            // get discont. pressure gradients
      //            LINALG::Matrix<3,1> graddiscp;
      //            //gradp = enr_derxy(i,j)*eprenp(j);
      //            for (int isd = 0; isd < nsd; ++isd)
      //            {
      //                graddiscp(isd) = 0.0;
      //                for (int iparam = 0; iparam < numparamdiscpres; ++iparam)
      //                    graddiscp(isd) += enr_derxy_discpres(isd,iparam)*ediscprenp(iparam);
      //            }

      // get pressure
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += shppres.d0(iparam)*eprenp(iparam);
      // get bodyforce in gausspoint
      //            LINALG::Matrix<3,1> bodyforce;
      //            bodyforce = 0.0;
      //            std::cout << bodyforce << std::endl;
      ///////////////LINALG::SerialDenseVector bodyforce_(enr_edeadng_(i,j)*enr_funct_(j));

      // compute stabilization parameters (3 taus)
      // compute velocity norm
      const double vel_norm = gpvelnp.Norm2();
      //const double strle = FLD::UTILS::Streamlength(shpvel.dx, shpvel.dy, shpvel.dz, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_M = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;
      // remark: stabilization parameters 'tau_stab_Mu' and 'tau_stab_Mp' are both computes based on
      //         'hk' instead of using the streamlength for 'tau_stab_Mu'
      COMBUST::UTILS::computeStabilizationParams(gpvelnp, xji,
          instationary, dynvisc, densaf, vel_norm, hk, hk, mk, timefac, dt, tautype,
          tau_stab_M, tau_stab_Mp, tau_stab_C); //tautype_franca_barrenechea_valentin_wall

      // integration factors and coefficients of single terms
      double rhsfac = fac;
      if (genalpha)
        rhsfac *= timefac / ga_alphaF; // (fac * gamma * dt ) / alpha_M

      /*------------------------- evaluate rhs vector at integration point ---*/
      LINALG::Matrix<nsd,1> rhsint;

      // body force vector
      // remark: unit is acceleration (not density weighted)
      static LINALG::Matrix<nsd,1> bodyforce;
      static LINALG::Matrix<nsd,1> noscalebodyforce;
      bodyforce.Clear();
      noscalebodyforce.Clear();
      //bodyforce(1) = -9.81;
      for (std::size_t isd = 0; isd < nsd; isd++)
      {
        for (std::size_t inode = 0; inode < numnode; inode++)
        {
          bodyforce(isd) += edead(isd,inode) * funct(inode);
          noscalebodyforce(isd) += epgrad(isd,inode) * funct(inode);
        }
      }
      // apply body force manually on Gauss point
      // remark: unit has to be acceleration (not density weighted)
      //std::cout << "BodyForce" << std::endl;
      //bodyforce(0) = 0.0;
      //bodyforce(1) = -9.81;
      //bodyforce(0) = 0.5/densaf;
      // scale body force: body force * density^{n+1}
      //std::cout << "BodyForce " << std::endl;
      //std::cout << bodyforce(0) << std::endl;
      //std::cout << bodyforce(1) << std::endl;
      //std::cout << bodyforce(2) << std::endl;
      bodyforce.Scale(densaf);

      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = bodyforce(isd) + noscalebodyforce(isd);
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = histvec(isd) + (bodyforce(isd) + noscalebodyforce(isd)) * timefac;
      }

      /*----------------- get numerical representation of single operators ---*/
      /* Convective term  u_old * grad u_old: */
      LINALG::Matrix<nsd,1> conv_old;
      //conv_old = vderxy(i, j)*gpvelnp(j);
      conv_old.Multiply(vderxy,gpvelnp);

      /* Viscous term  div epsilon(u_old) */
      LINALG::Matrix<nsd,1> visc_old;
      visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
      visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
      visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

      // evaluate residual once for all stabilisation right hand sides
      // pres?
      LINALG::Matrix<nsd,1> res_old;
      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = densaf*conv_old(isd) + gradp(isd) - 2.0*dynvisc*visc_old(isd) - bodyforce(isd) - noscalebodyforce(isd);
        if (instationary)
        {
          for (size_t isd = 0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpaccam(isd);
        }
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

        if (instationary)
        {
          for(size_t isd=0; isd<nsd; ++isd)
            res_old(isd) += densaf * gpvelnp(isd);
        }
        //res_old += gpvelnp;
      }


      /* Reactive term  u:  funct */
      /* linearise convective term */

      /*--- convective part u_old * grad (funct) --------------------------*/
      /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
             with  N .. form function matrix                                   */
      //const LINALG::SerialDenseVector enr_conv_c_(enr_derxy(j,i)*gpvelnp(j));
      static LINALG::Matrix<shpVecSize,1> enr_conv_c_;
      //static ShpVec enr_conv_c_;
      enr_conv_c_.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_conv_c_(iparam) += shpvel.dx(iparam)*gpvelnp(0);
        enr_conv_c_(iparam) += shpvel.dy(iparam)*gpvelnp(1);
        enr_conv_c_(iparam) += shpvel.dz(iparam)*gpvelnp(2);
      }

      //              /*--- convective grid part u_G * grad (funct) -----------------------*/
      //              /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
      //              enr_conv_g_ = 0.0;


      /*--- viscous term  - grad * epsilon(u): ----------------------------*/
      /*   /                                                \
               |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
             1 |                                                |
             - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
             2 |                                                |
               |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
               \                                                /

               with N_x .. x-line of N
               N_y .. y-line of N                                             */
      static COMBUST::EnrViscs2<shpVecSize> enr_viscs2;

      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_viscs2.xx(iparam) = 0.5 * (2.0 * shpvel.dxdx(iparam) + shpvel.dydy(iparam) + shpvel.dzdz(iparam));
        enr_viscs2.xy(iparam) = 0.5 *  shpvel.dxdy(iparam);
        enr_viscs2.xz(iparam) = 0.5 *  shpvel.dxdz(iparam);
        enr_viscs2.yx(iparam) = 0.5 *  shpvel.dydx(iparam);
        enr_viscs2.yy(iparam) = 0.5 * (shpvel.dxdx(iparam) + 2.0 * shpvel.dydy(iparam) + shpvel.dzdz(iparam));
        enr_viscs2.yz(iparam) = 0.5 *  shpvel.dydz(iparam);
        enr_viscs2.zx(iparam) = 0.5 *  shpvel.dzdx(iparam);
        enr_viscs2.zy(iparam) = 0.5 *  shpvel.dzdy(iparam);
        enr_viscs2.zz(iparam) = 0.5 * (shpvel.dxdx(iparam) + shpvel.dydy(iparam) + 2.0 * shpvel.dzdz(iparam));
      }

      //-------------------------------------------------------------------
      //           CROSS-/REYNODSSTRESS
      //-------------------------------------------------------------------
      bool cross_reynolds = false;
      // subgrid-scale velocity vector in gauss point
      LINALG::Matrix<nsd,1> sgvel;
      sgvel.Clear();
      // linearisation of subgrid-scale convection, convective part
      LINALG::Matrix<shpVecSize,1> sgconv_c;
      sgconv_c.Clear();
      if (cross_reynolds)
      {
        // compute subgrid-scale velocity
        sgvel.Update(-tau_stab_Mp*theta,res_old,0.0);
        // compute subgrid-scale convective operator
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          sgconv_c(iparam) += shpvel.dx(iparam)*sgvel(0);
          sgconv_c(iparam) += shpvel.dy(iparam)*sgvel(1);
          sgconv_c(iparam) += shpvel.dz(iparam)*sgvel(2);
        }
      }

// for gmsh output of subgrid velocity
#ifdef SUGRVEL_OUTPUT
      const bool screen_out = false;

      //const std::string filename = IO::GMSH::GetFileName("SubgridVelocityFluid", step, screen_out, 0);
      const std::string filename = IO::GMSH::GetFileName("SubgridVelocityFluid", 0, screen_out, 0);
      std::ofstream gmshfilecontent(filename.c_str(), ios_base::out | ios_base::app);

      LINALG::Matrix<3,1> xsi;
      for (std::size_t idim=0;idim<nsd;idim++)
          {xsi(idim) = posXiDomain(idim);}
      LINALG::SerialDenseMatrix xyze(3,numnode);
      for (std::size_t inode=0;inode<numnode;inode++)
      {
          xyze(0,inode) = ele->Nodes()[inode]->X()[0];
          xyze(1,inode) = ele->Nodes()[inode]->X()[1];
          xyze(2,inode) = ele->Nodes()[inode]->X()[2];
      }
      // transform gp from local (element) coordinates to global (physical) coordinates
      GEO::elementToCurrentCoordinatesInPlace(ele->Shape(), xyze, xsi);

      sgvel.Update(-tau_stab_Mp*theta,res_old,0.0);
      IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, sgvel, xsi, gmshfilecontent);

      gmshfilecontent.close();
//    sgvel.Clear();
//    sgconv_c.Clear();

      //const std::string filename2 = IO::GMSH::GetFileName("Residual", step, screen_out, 0);
      const std::string filename2 = IO::GMSH::GetFileName("Residual", 0, screen_out, 0);
      std::ofstream gmshfilecontent2(filename2.c_str(), ios_base::out | ios_base::app);

      LINALG::Matrix<nsd,1> res_old_sugrvel;
      for (std::size_t i=0; i<nsd; i++)
          res_old_sugrvel(i) = res_old(i) / dt;
      IO::GMSH::cellWithVectorFieldToStream(DRT::Element::point1, res_old_sugrvel, xsi, gmshfilecontent2);

      gmshfilecontent2.close();

      //const std::string filename3 = IO::GMSH::GetFileName("Tau", step, screen_out, 0);
      const std::string filename3 = IO::GMSH::GetFileName("Tau", 0, screen_out, 0);
      std::ofstream gmshfilecontent3(filename3.c_str(), ios_base::out | ios_base::app);
      IO::GMSH::cellWithScalarToStream(DRT::Element::point1, tau_stab_Mp*timefac, xsi, gmshfilecontent3);
      gmshfilecontent3.close();
#endif

//test
//      if (ele->Id()==2076) //ele->Id()==1857
//      {
////         std::cout << ele->Id() << std::endl;
////         std::cout << ele->Bisected() << std::endl;
////         std::cout << cell->getDomainPlus() << std::endl;
////         std::cout << vdiv << std::endl;
//         for (std::size_t i = 0; i<nsd; i++)
//         {
////             std::cout <<"subgrid vel " << sgvel(i) << std::endl;
////           std::cout<< "Residuum" << res_old(i) << std::endl;
////           std::cout<< "tau " << tau_stab_M*dt*theta << std::endl;
////           std::cout<< "Histvektor" << histvec(i) << std::endl;
////           std::cout<< "Geschw" << gpvelnp(i) << std::endl;
////           std::cout<< "Druck" << gradp(i) << std::endl;
////           std::cout<< "Konvek" << conv_old(i) << std::endl;
////           std::cout<< "Epsilon" << visc_old(i) << std::endl;
////           std::cout<< "Bodyforce" << bodyforce(i) << std::endl;
////           std::cout<< "dens" << dens << std::endl;
////           std::cout<< "dynvisc" << dynvisc << std::endl;
////           std::cout<< "timefac" << timefac << std::endl;
//         }
//      }

      // convective part of convection and sub-grid scale
      LINALG::Matrix<shpVecSize,1> enr_conv_c_sg;
      enr_conv_c_sg.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        enr_conv_c_sg(iparam) = enr_conv_c_(iparam) + sgconv_c(iparam);
      // Note: linearisation of sub-grid scale is taken into acstd::cout only in convective part
      //----------------------------------------------------------------------

#if 1
      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize>(
         assembler, shpvel, shppres, fac, timefac, rhsfac, densaf, dynvisc,
         gpvelnp, pres, gpaccam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c_, enr_conv_c_sg, enr_viscs2,
         instationary, genalpha, newton, pstab, supg, graddiv,
         tau_stab_M, tau_stab_Mp, tau_stab_C);
#else

      // modify stabilization
      // TODO: does this have to be modified with respect to the density, as it was before?
      tau_stab_M  /= timefac;
      tau_stab_Mp /= timefac;
      tau_stab_C /= timefac;

      //----------------------------------------------------------------------
      //                            GALERKIN PART
      //----------------------------------------------------------------------
      if (instationary)
      {
        /* inertia (contribution to mass matrix) */
        /*
                                     /           \
                                    |             |
                                    | roh Du , v  |
                                    |             |
                                     \           /
         */
        assembler.template Matrix<Velx,Velx>(shpvel.d0, fac*densaf, shpvel.d0);
        assembler.template Matrix<Vely,Vely>(shpvel.d0, fac*densaf, shpvel.d0);
        assembler.template Matrix<Velz,Velz>(shpvel.d0, fac*densaf, shpvel.d0);

        assembler.template Vector<Velx>(shpvel.d0, -fac*densaf*gpvelnp(0));
        assembler.template Vector<Vely>(shpvel.d0, -fac*densaf*gpvelnp(1));
        assembler.template Vector<Velz>(shpvel.d0, -fac*densaf*gpvelnp(2));
      }

      /* convection, convective part */
      /*
                         /                                                      \
                        |         / n+1       \                                  |
                        | v , roh| u   o nabla | Du - rho*(tau_M*R_M o nabla) Du |
                        |         \ (i)       /          subgrid velocity        |
                         \                            cross-stress              /
       */
      assembler.template Matrix<Velx,Velx>(shpvel.d0, timefacfac*densaf, enr_conv_c_sg);//enr_conv_c_+sgconv_c);
      assembler.template Matrix<Vely,Vely>(shpvel.d0, timefacfac*densaf, enr_conv_c_sg);//enr_conv_c_+sgconv_c);
      assembler.template Matrix<Velz,Velz>(shpvel.d0, timefacfac*densaf, enr_conv_c_sg);//enr_conv_c_+sgconv_c);

      assembler.template Vector<Velx>(shpvel.d0, -timefacfac*densaf*(gpvelnp(0)*vderxy(0,0) // check order
          +gpvelnp(1)*vderxy(0,1)
          +gpvelnp(2)*vderxy(0,2)+sgvel(0)*vderxy(0,0)
          +sgvel(1)*vderxy(0,1)
          +sgvel(2)*vderxy(0,2)));
      assembler.template Vector<Vely>(shpvel.d0, -timefacfac*densaf*(gpvelnp(0)*vderxy(1,0)
          +gpvelnp(1)*vderxy(1,1)
          +gpvelnp(2)*vderxy(1,2)+sgvel(0)*vderxy(1,0)
          +sgvel(1)*vderxy(1,1)
          +sgvel(2)*vderxy(1,2)));
      assembler.template Vector<Velz>(shpvel.d0, -timefacfac*densaf*(gpvelnp(0)*vderxy(2,0)
          +gpvelnp(1)*vderxy(2,1)
          +gpvelnp(2)*vderxy(2,2)+sgvel(0)*vderxy(2,0)
          +sgvel(1)*vderxy(2,1)
          +sgvel(2)*vderxy(2,2)));

      if (newton)
      {
        /*  convection, reactive part */
        /*
                       /                             \
                      |         /          \   n+1   |
                      | v ,roh | Du o nabla | u      |
                      |         \          /   (i)   |
                       \                             /
         */
         assembler.template Matrix<Velx,Velx>(shpvel.d0, timefacfac*densaf*vderxy(0,0), shpvel.d0);
         assembler.template Matrix<Velx,Vely>(shpvel.d0, timefacfac*densaf*vderxy(0,1), shpvel.d0);
         assembler.template Matrix<Velx,Velz>(shpvel.d0, timefacfac*densaf*vderxy(0,2), shpvel.d0);
         assembler.template Matrix<Vely,Velx>(shpvel.d0, timefacfac*densaf*vderxy(1,0), shpvel.d0);
         assembler.template Matrix<Vely,Vely>(shpvel.d0, timefacfac*densaf*vderxy(1,1), shpvel.d0);
         assembler.template Matrix<Vely,Velz>(shpvel.d0, timefacfac*densaf*vderxy(1,2), shpvel.d0);
         assembler.template Matrix<Velz,Velx>(shpvel.d0, timefacfac*densaf*vderxy(2,0), shpvel.d0);
         assembler.template Matrix<Velz,Vely>(shpvel.d0, timefacfac*densaf*vderxy(2,1), shpvel.d0);
         assembler.template Matrix<Velz,Velz>(shpvel.d0, timefacfac*densaf*vderxy(2,2), shpvel.d0);
      }

      /* Viskositaetsterm */
      /*
                          /                                 \
                         |       / \                  /  \   |
                         |  eps | v | , 2\mu epsilon | Du |  |
                         |       \ /                  \  /   |
                          \                                 /
       */
      // dynvisc is dynamic viscosity
      assembler.template Matrix<Velx,Velx>(shpvel.dx,   2.0*dynvisc*timefacfac, shpvel.dx);
      assembler.template Matrix<Velx,Velx>(shpvel.dy,     dynvisc*timefacfac, shpvel.dy);
      assembler.template Matrix<Velx,Vely>(shpvel.dy,     dynvisc*timefacfac, shpvel.dx);
      assembler.template Matrix<Velx,Velx>(shpvel.dz,     dynvisc*timefacfac, shpvel.dz);
      assembler.template Matrix<Velx,Velz>(shpvel.dz,     dynvisc*timefacfac, shpvel.dx);

      assembler.template Matrix<Vely,Vely>(shpvel.dx,     dynvisc*timefacfac, shpvel.dx);
      assembler.template Matrix<Vely,Velx>(shpvel.dx,     dynvisc*timefacfac, shpvel.dy);
      assembler.template Matrix<Vely,Vely>(shpvel.dy,   2.0*dynvisc*timefacfac, shpvel.dy);
      assembler.template Matrix<Vely,Vely>(shpvel.dz,     dynvisc*timefacfac, shpvel.dz);
      assembler.template Matrix<Vely,Velz>(shpvel.dz,     dynvisc*timefacfac, shpvel.dy);

      assembler.template Matrix<Velz,Velz>(shpvel.dx,     dynvisc*timefacfac, shpvel.dx);
      assembler.template Matrix<Velz,Velx>(shpvel.dx,     dynvisc*timefacfac, shpvel.dz);
      assembler.template Matrix<Velz,Velz>(shpvel.dy,     dynvisc*timefacfac, shpvel.dy);
      assembler.template Matrix<Velz,Vely>(shpvel.dy,     dynvisc*timefacfac, shpvel.dz);
      assembler.template Matrix<Velz,Velz>(shpvel.dz,   2.0*dynvisc*timefacfac, shpvel.dz);

      assembler.template Vector<Velx>(shpvel.dx,     -dynvisc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
      assembler.template Vector<Velx>(shpvel.dy,     -dynvisc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
      assembler.template Vector<Velx>(shpvel.dz,     -dynvisc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));

      assembler.template Vector<Vely>(shpvel.dx,     -dynvisc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
      assembler.template Vector<Vely>(shpvel.dy,     -dynvisc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
      assembler.template Vector<Vely>(shpvel.dz,     -dynvisc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));

      assembler.template Vector<Velz>(shpvel.dx,     -dynvisc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
      assembler.template Vector<Velz>(shpvel.dy,     -dynvisc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
      assembler.template Vector<Velz>(shpvel.dz,     -dynvisc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));

      /* Druckterm */
      /*
                            /                \
                           |                  |
                         - |  nabla o v , Dp  |
                           |                  |
                            \                /
       */
      // dynamic pressure
      assembler.template Matrix<Velx,Pres>(shpvel.dx, -timefacfac, shppres.d0);
      assembler.template Matrix<Vely,Pres>(shpvel.dy, -timefacfac, shppres.d0);
      assembler.template Matrix<Velz,Pres>(shpvel.dz, -timefacfac, shppres.d0);

      assembler.template Vector<Velx>(shpvel.dx, timefacfac*pres);
      assembler.template Vector<Vely>(shpvel.dy, timefacfac*pres);
      assembler.template Vector<Velz>(shpvel.dz, timefacfac*pres);

      /* Divergenzfreiheit - continuity equation*/
      /*
                           /              \
                          |                |
                          | q , nabla o Du |
                          |                |
                           \              /
       */
      assembler.template Matrix<Pres,Velx>(shppres.d0, timefacfac, shpvel.dx);
      assembler.template Matrix<Pres,Vely>(shppres.d0, timefacfac, shpvel.dy);
      assembler.template Matrix<Pres,Velz>(shppres.d0, timefacfac, shpvel.dz);

      //Residuum of continuity equation
      const double trace_gamma = (vderxy(0, 0) + vderxy(1, 1) + vderxy(2, 2));
      assembler.template Vector<Pres>(shppres.d0, -timefacfac*trace_gamma);

      // source term of the right hand side
      assembler.template Vector<Velx>(shpvel.d0, fac*rhsint(0));
      assembler.template Vector<Vely>(shpvel.d0, fac*rhsint(1));
      assembler.template Vector<Velz>(shpvel.d0, fac*rhsint(2));


      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART
      if(pstab)
      {
        const double timetauMp  = timefac * tau_stab_Mp * fac;
        if (instationary)
        {
          /* pressure stabilisation: inertia */
          /*
                                /                 \
                               |                   |
                               | roh Du , nabla q  |
                               |                   |
                                \                 /
           */
          //fac*densaf
          assembler.template Matrix<Pres,Velx>(shppres.dx, timetauMp*densaf, shpvel.d0);
          assembler.template Matrix<Pres,Vely>(shppres.dy, timetauMp*densaf, shpvel.d0);
          assembler.template Matrix<Pres,Velz>(shppres.dz, timetauMp*densaf, shpvel.d0);
        }
        const double ttimetauMp = timefac * timefac * tau_stab_Mp * fac;
        /* pressure stabilisation: convection, convective part */
        /*
                          /                                \
                         |                / n+1       \     |
                         | nabla q ,roh  | u   o nabla | Du |
                         |                \ i         /     |
                          \                                /
         */
        assembler.template Matrix<Pres,Velx>(shppres.dx, ttimetauMp*densaf, enr_conv_c_);
        assembler.template Matrix<Pres,Vely>(shppres.dy, ttimetauMp*densaf, enr_conv_c_);
        assembler.template Matrix<Pres,Velz>(shppres.dz, ttimetauMp*densaf, enr_conv_c_);

        if (newton)
        {
          /*  pressure stabilisation: convection, reactive part
                          /                                \
                         |              /          \   n+1  |
                         | grad q , roh| Du o nabla | u     |
                         |              \          /   (i)  |
                          \                                /
           */
          assembler.template Matrix<Pres,Velx>(shppres.dx, ttimetauMp*densaf*vderxy(0,0), shpvel.d0);
          assembler.template Matrix<Pres,Velx>(shppres.dy, ttimetauMp*densaf*vderxy(1,0), shpvel.d0);
          assembler.template Matrix<Pres,Velx>(shppres.dz, ttimetauMp*densaf*vderxy(2,0), shpvel.d0);

          assembler.template Matrix<Pres,Vely>(shppres.dx, ttimetauMp*densaf*vderxy(0,1), shpvel.d0);
          assembler.template Matrix<Pres,Vely>(shppres.dy, ttimetauMp*densaf*vderxy(1,1), shpvel.d0);
          assembler.template Matrix<Pres,Vely>(shppres.dz, ttimetauMp*densaf*vderxy(2,1), shpvel.d0);

          assembler.template Matrix<Pres,Velz>(shppres.dx, ttimetauMp*densaf*vderxy(0,2), shpvel.d0);
          assembler.template Matrix<Pres,Velz>(shppres.dy, ttimetauMp*densaf*vderxy(1,2), shpvel.d0);
          assembler.template Matrix<Pres,Velz>(shppres.dz, ttimetauMp*densaf*vderxy(2,2), shpvel.d0);
        }

        /* pressure stabilisation: viscosity (-L_visc_u) */
        /*
                           /                                      \
                          |                                  /  \  |
                        - |  nabla q , nabla o 2/mu epsilon | Du | |
                          |                                  \  /  |
                           \                                      /
         */
        assembler.template Matrix<Pres,Velx>(shppres.dx, -2.0*dynvisc*ttimetauMp, enr_viscs2.xx);
        assembler.template Matrix<Pres,Vely>(shppres.dx, -2.0*dynvisc*ttimetauMp, enr_viscs2.xy);
        assembler.template Matrix<Pres,Velz>(shppres.dx, -2.0*dynvisc*ttimetauMp, enr_viscs2.xz);

        assembler.template Matrix<Pres,Velx>(shppres.dy, -2.0*dynvisc*ttimetauMp, enr_viscs2.yx);
        assembler.template Matrix<Pres,Vely>(shppres.dy, -2.0*dynvisc*ttimetauMp, enr_viscs2.yy);
        assembler.template Matrix<Pres,Velz>(shppres.dy, -2.0*dynvisc*ttimetauMp, enr_viscs2.yz);

        assembler.template Matrix<Pres,Velx>(shppres.dz, -2.0*dynvisc*ttimetauMp, enr_viscs2.zx);
        assembler.template Matrix<Pres,Vely>(shppres.dz, -2.0*dynvisc*ttimetauMp, enr_viscs2.zy);
        assembler.template Matrix<Pres,Velz>(shppres.dz, -2.0*dynvisc*ttimetauMp, enr_viscs2.zz);

        /* pressure stabilisation: pressure( L_pres_p) */
        /*
                          /                    \
                         |                      |
                         |  nabla q , nabla Dp  |
                         |                      |
                          \                    /
         */
        assembler.template Matrix<Pres,Pres>(shppres.dx, ttimetauMp, shppres.dx);
        assembler.template Matrix<Pres,Pres>(shppres.dy, ttimetauMp, shppres.dy);
        assembler.template Matrix<Pres,Pres>(shppres.dz, ttimetauMp, shppres.dz);

        // pressure stabilization
        assembler.template Vector<Pres>(shppres.dx, -timetauMp*res_old(0));
        assembler.template Vector<Pres>(shppres.dy, -timetauMp*res_old(1));
        assembler.template Vector<Pres>(shppres.dz, -timetauMp*res_old(2));

      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART
      if(supg)
      {
        const double timetauM   = timefac * tau_stab_M * fac;
        if (instationary)
        {
          /* supg stabilisation: inertia  */
          /*
                              /                                                        \
                             |              / n+1       \                               |
                             | roh Du , roh| u   o nabla | v - rho*(tau_M*R_M o nabla)v |
                             |              \ (i)       /           subgrid velocity    |
                              \                                    reynolds-stress     /
           */
          // timetauM * densaf * densaf, da einmal Dichte von Wichtungsfunktion und einmal von Testfunktion
          // vgl. Paper bzw FEFluid-Skript bzw Blatt
          assembler.template Matrix<Velx,Velx>(enr_conv_c_sg, timetauM*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Vely,Vely>(enr_conv_c_sg, timetauM*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Velz,Velz>(enr_conv_c_sg, timetauM*densaf*densaf, shpvel.d0);

          if (newton)
          {
            /* supg stabilisation: inertia, linearisation of testfunction  */
            /*
                                   /                                  \
                                  |      n+1          /          \     |
                                  | roh u      , roh | Du o nabla | v  |
                                  |      (i)          \          /     |
                                   \                                   /

             */
            assembler.template Matrix<Velx,Velx>(shpvel.dx, timetauM*densaf*densaf*gpvelnp(0), shpvel.d0);
            assembler.template Matrix<Velx,Vely>(shpvel.dy, timetauM*densaf*densaf*gpvelnp(0), shpvel.d0);
            assembler.template Matrix<Velx,Velz>(shpvel.dz, timetauM*densaf*densaf*gpvelnp(0), shpvel.d0);
            assembler.template Matrix<Vely,Velx>(shpvel.dx, timetauM*densaf*densaf*gpvelnp(1), shpvel.d0);
            assembler.template Matrix<Vely,Vely>(shpvel.dy, timetauM*densaf*densaf*gpvelnp(1), shpvel.d0);
            assembler.template Matrix<Vely,Velz>(shpvel.dz, timetauM*densaf*densaf*gpvelnp(1), shpvel.d0);
            assembler.template Matrix<Velz,Velx>(shpvel.dx, timetauM*densaf*densaf*gpvelnp(2), shpvel.d0);
            assembler.template Matrix<Velz,Vely>(shpvel.dy, timetauM*densaf*densaf*gpvelnp(2), shpvel.d0);
            assembler.template Matrix<Velz,Velz>(shpvel.dz, timetauM*densaf*densaf*gpvelnp(2), shpvel.d0);
          }
        }

        const double ttimetauM  = timefac * timefac * tau_stab_M * fac;
        /* supg stabilisation: convective part ( L_conv_u) */
        /*
                     /                                                                           \
                    |     / n+1        \                                      / n+1        \      |
                    | roh| u    o nabla | v - rho*(tau_M*R_M o nabla)v , roh | u    o nabla | Du  |
                    |     \ (i)        /         subgrid velocity             \ (i)        /      |
                     \                           reynolds-stress                                 /
         */
        assembler.template Matrix<Velx,Velx>(enr_conv_c_sg, ttimetauM*densaf*densaf, enr_conv_c_);
        assembler.template Matrix<Vely,Vely>(enr_conv_c_sg, ttimetauM*densaf*densaf, enr_conv_c_);
        assembler.template Matrix<Velz,Velz>(enr_conv_c_sg, ttimetauM*densaf*densaf, enr_conv_c_);
        /* supg stabilisation: pressure part  ( L_pres_p) */
        /*
                          /                                                           \
                         |      / n+1       \                                          |
                         | roh | u   o nabla | v - rho*(tau_M*R_M o nabla)v , nabla Dp |
                         |      \ (i)       /         subgrid velocity                 |
                          \                            reynolds-stress                /
         */
        assembler.template Matrix<Velx,Pres>(enr_conv_c_sg, ttimetauM*densaf, shppres.dx);
        assembler.template Matrix<Vely,Pres>(enr_conv_c_sg, ttimetauM*densaf, shppres.dy);
        assembler.template Matrix<Velz,Pres>(enr_conv_c_sg, ttimetauM*densaf, shppres.dz);

        /* supg stabilisation: viscous part  (-L_visc_u) */
        /*
                      /                                                                     \
                     |               /  \       / n+1        \                               |
                   - |  nabla o eps | Du |, roh| u    o nabla | v - rho*(tau_M*R_M o nabla)v |
                     |               \  /       \ (i)        /      subgrid velocity         |
                      \                                             reynolds-stress         /
         */
        assembler.template Matrix<Velx,Velx>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.xx);
        assembler.template Matrix<Velx,Vely>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.xy);
        assembler.template Matrix<Velx,Velz>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.xz);

        assembler.template Matrix<Vely,Velx>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.yx);
        assembler.template Matrix<Vely,Vely>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.yy);
        assembler.template Matrix<Vely,Velz>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.yz);

        assembler.template Matrix<Velz,Velx>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.zx);
        assembler.template Matrix<Velz,Vely>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.zy);
        assembler.template Matrix<Velz,Velz>(enr_conv_c_sg, -2.0*dynvisc*ttimetauM*densaf, enr_viscs2.zz);

        if (newton)
        {
          /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
          /*
                               /                                                  |
                              |      |            |   n+1     | n+1          |    |
                              | roh  | Du o nabla | u    , roh| u    o nabla | v  |
                              |      |            |   (i)     | (i)          |    |
                               \                                                  |
           */
          assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(0,0), shpvel.d0);
          assembler.template Matrix<Velx,Vely>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(0,1), shpvel.d0);
          assembler.template Matrix<Velx,Velz>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(0,2), shpvel.d0);

          assembler.template Matrix<Vely,Velx>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(1,0), shpvel.d0);
          assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(1,1), shpvel.d0);
          assembler.template Matrix<Vely,Velz>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(1,2), shpvel.d0);

          assembler.template Matrix<Velz,Velx>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(2,0), shpvel.d0);
          assembler.template Matrix<Velz,Vely>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(2,1), shpvel.d0);
          assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM*densaf*densaf*vderxy(2,2), shpvel.d0);

          /*
                             /                                                  |
                            |       / n+1         |   n+1      /           |    |
                            |  roh | u    o nabla | u    ,roh | Du o nabla | v  |
                            |      | (i)          |   (i)     |            |    |
                             \                                                  |
           */
          const double con0 = ttimetauM*(gpvelnp(0)*vderxy(0,0) + gpvelnp(1)*vderxy(0,1) + gpvelnp(2)*vderxy(0,2));
          assembler.template Matrix<Velx,Velx>(shpvel.dx, con0*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Velx,Vely>(shpvel.dy, con0*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Velx,Velz>(shpvel.dz, con0*densaf*densaf, shpvel.d0);

          const double con1 = ttimetauM*(gpvelnp(0)*vderxy(1,0) + gpvelnp(1)*vderxy(1,1) + gpvelnp(2)*vderxy(1,2));
          assembler.template Matrix<Vely,Velx>(shpvel.dx, con1*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Vely,Vely>(shpvel.dy, con1*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Vely,Velz>(shpvel.dz, con1*densaf*densaf, shpvel.d0);

          const double con2 = ttimetauM*(gpvelnp(0)*vderxy(2,0) + gpvelnp(1)*vderxy(2,1) + gpvelnp(2)*vderxy(2,2));
          assembler.template Matrix<Velz,Velx>(shpvel.dx, con2*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Velz,Vely>(shpvel.dy, con2*densaf*densaf, shpvel.d0);
          assembler.template Matrix<Velz,Velz>(shpvel.dz, con2*densaf*densaf, shpvel.d0);

          /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
          /*
                                    /                                  \
                                   |         n+1       /          \     |
                                   |  nabla p    ,roh | Du o nabla | v  |
                                   |         (i)       \          /     |
                                    \                                  /
           */
          assembler.template Matrix<Velx,Velx>(shpvel.dx, ttimetauM*densaf*gradp(0), shpvel.d0);
          assembler.template Matrix<Velx,Vely>(shpvel.dy, ttimetauM*densaf*gradp(0), shpvel.d0);
          assembler.template Matrix<Velx,Velz>(shpvel.dz, ttimetauM*densaf*gradp(0), shpvel.d0);

          assembler.template Matrix<Vely,Velx>(shpvel.dx, ttimetauM*densaf*gradp(1), shpvel.d0);
          assembler.template Matrix<Vely,Vely>(shpvel.dy, ttimetauM*densaf*gradp(1), shpvel.d0);
          assembler.template Matrix<Vely,Velz>(shpvel.dz, ttimetauM*densaf*gradp(1), shpvel.d0);

          assembler.template Matrix<Velz,Velx>(shpvel.dx, ttimetauM*densaf*gradp(2), shpvel.d0);
          assembler.template Matrix<Velz,Vely>(shpvel.dy, ttimetauM*densaf*gradp(2), shpvel.d0);
          assembler.template Matrix<Velz,Velz>(shpvel.dz, ttimetauM*densaf*gradp(2), shpvel.d0);

          /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
          /*
                              /                                            \
                             |               / n+1 \       /          \     |
                           - |  nabla o eps | u     |,roh | Du o nabla | v  |
                             |               \ (i) /       \          /     |
                              \                                            /
           */
          assembler.template Matrix<Velx,Velx>(shpvel.dx, -2.0*dynvisc*ttimetauM*densaf*visc_old(0), shpvel.d0);
          assembler.template Matrix<Velx,Vely>(shpvel.dy, -2.0*dynvisc*ttimetauM*densaf*visc_old(0), shpvel.d0);
          assembler.template Matrix<Velx,Velz>(shpvel.dz, -2.0*dynvisc*ttimetauM*densaf*visc_old(0), shpvel.d0);

          assembler.template Matrix<Vely,Velx>(shpvel.dx, -2.0*dynvisc*ttimetauM*densaf*visc_old(1), shpvel.d0);
          assembler.template Matrix<Vely,Vely>(shpvel.dy, -2.0*dynvisc*ttimetauM*densaf*visc_old(1), shpvel.d0);
          assembler.template Matrix<Vely,Velz>(shpvel.dz, -2.0*dynvisc*ttimetauM*densaf*visc_old(1), shpvel.d0);

          assembler.template Matrix<Velz,Velx>(shpvel.dx, -2.0*dynvisc*ttimetauM*densaf*visc_old(2), shpvel.d0);
          assembler.template Matrix<Velz,Vely>(shpvel.dy, -2.0*dynvisc*ttimetauM*densaf*visc_old(2), shpvel.d0);
          assembler.template Matrix<Velz,Velz>(shpvel.dz, -2.0*dynvisc*ttimetauM*densaf*visc_old(2), shpvel.d0);

          /* supg stabilisation: bodyforce part, linearisation of test function */

          /*
                                  /                                \
                                 |                 /          \     |
                               - |  rhsint   , roh| Du o nabla | v  |
                                 |                 \          /     |
                                  \                                /

           */

          assembler.template Matrix<Velx,Velx>(shpvel.dx, -timetauM*densaf*rhsint(0), shpvel.d0);
          assembler.template Matrix<Velx,Vely>(shpvel.dy, -timetauM*densaf*rhsint(0), shpvel.d0);
          assembler.template Matrix<Velx,Velz>(shpvel.dz, -timetauM*densaf*rhsint(0), shpvel.d0);

          assembler.template Matrix<Vely,Velx>(shpvel.dx, -timetauM*densaf*rhsint(1), shpvel.d0);
          assembler.template Matrix<Vely,Vely>(shpvel.dy, -timetauM*densaf*rhsint(1), shpvel.d0);
          assembler.template Matrix<Vely,Velz>(shpvel.dz, -timetauM*densaf*rhsint(1), shpvel.d0);

          assembler.template Matrix<Velz,Velx>(shpvel.dx, -timetauM*densaf*rhsint(2), shpvel.d0);
          assembler.template Matrix<Velz,Vely>(shpvel.dy, -timetauM*densaf*rhsint(2), shpvel.d0);
          assembler.template Matrix<Velz,Velz>(shpvel.dz, -timetauM*densaf*rhsint(2), shpvel.d0);
        } // if newton

        // supg stabilisation
        assembler.template Vector<Velx>(enr_conv_c_sg, -timetauM*densaf*res_old(0));
        assembler.template Vector<Vely>(enr_conv_c_sg, -timetauM*densaf*res_old(1));
        assembler.template Vector<Velz>(enr_conv_c_sg, -timetauM*densaf*res_old(2));
      }


      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART
      if(graddiv)
      {
        const double timefac_timefac_tau_C=timefac*timefac*tau_stab_C * fac;
        const double timefac_timefac_tau_C_divunp=timefac_timefac_tau_C*(vderxy(0, 0)+vderxy(1, 1)+vderxy(2, 2));
        /* continuity stabilisation on left hand side */
        /*
                         /                        \
                        |                          |
                        | nabla o Du  , nabla o v  |
                        |                          |
                         \                        /
         */
        assembler.template Matrix<Velx,Velx>(shpvel.dx, timefac_timefac_tau_C, shpvel.dx);
        assembler.template Matrix<Velx,Vely>(shpvel.dx, timefac_timefac_tau_C, shpvel.dy);
        assembler.template Matrix<Velx,Velz>(shpvel.dx, timefac_timefac_tau_C, shpvel.dz);

        assembler.template Matrix<Vely,Velx>(shpvel.dy, timefac_timefac_tau_C, shpvel.dx);
        assembler.template Matrix<Vely,Vely>(shpvel.dy, timefac_timefac_tau_C, shpvel.dy);
        assembler.template Matrix<Vely,Velz>(shpvel.dy, timefac_timefac_tau_C, shpvel.dz);

        assembler.template Matrix<Velz,Velx>(shpvel.dz, timefac_timefac_tau_C, shpvel.dx);
        assembler.template Matrix<Velz,Vely>(shpvel.dz, timefac_timefac_tau_C, shpvel.dy);
        assembler.template Matrix<Velz,Velz>(shpvel.dz, timefac_timefac_tau_C, shpvel.dz);

        assembler.template Vector<Velx>(shpvel.dx, -timefac_timefac_tau_C_divunp);
        assembler.template Vector<Vely>(shpvel.dy, -timefac_timefac_tau_C_divunp);
        assembler.template Vector<Velz>(shpvel.dz, -timefac_timefac_tau_C_divunp);
      } // endif graddiv
#endif

    } // end loop over gauss points
  } // end loop over integration cells

  return;
}
}


#endif

