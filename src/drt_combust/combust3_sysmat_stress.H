/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_stress.H

\brief system matrix stress-based formulation
       for premixed combustion problem / 3D combustion element

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

#include <Teuchos_TimeMonitor.hpp>

#include "combust3_sysmat.H"
#include "combust3_local_assembler.H"
#include "combust3_utils.H"
#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_f3/xfluid3_utils.H"
#include "../drt_f3/fluid3_stabilization.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_xfem/enrichment.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"


using namespace XFEM::PHYSICS;

namespace COMBUST
{

/*!
  \brief build domain integral entries for combustion problem

  some remarks on the following sysmat entries:
  - this is an incremental formulation (details see e.g. Diss. Wall and BACI guide)
  - in comparison to standard fluid sysmat entries, here entries are scaled by the density
  - scaling has to be performed with quantities (density, dynvisc) from corresponding time step (n, n+1)
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          size_t shpVecSize,
          size_t shpVecSizeStress>
void BuildDomainIntegralsStress(
    LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFunc<shpVecSize>& shp,           // shape functions and derivatives
    const LINALG::Matrix<shpVecSizeStress,1>& shp_tau, //
    const double& fac,                                 // factor for integration in space
    const double& timefac,                             // factor for integration in time
    const double& timefacfac,                          // overall integration factor
    const double& densnp,                              // density at n+1
    const double& dynvisc,                             // dynamic viscosity at n+1
    const LINALG::Matrix<3,1>& velnp,                  // velocity at Gaussian point (GP) at n+1
    const double&              pres,                   // pressure at GP at n+1
    const LINALG::Matrix<3,1>& gradp,                  // pressure gradient at GP at n+1
    const LINALG::Matrix<3,3>& vderxy,                 // velocity derivative at GP at n+1
    const LINALG::Matrix<3,1>& rhsint,                 // time RHS (at n) and body force (at n+1) at GP
    const LINALG::Matrix<3,1>& res_old,                // residual at GP at n
    const LINALG::Matrix<3,1>& visc_old,               // viscous term at GP at n
    const LINALG::Matrix<3,3>& tau,                    // element stress unknowns
    const LINALG::Matrix<shpVecSize,1>& enr_conv_c,    // convective term at GP at n
    const COMBUST::EnrViscs2<shpVecSize>& enr_viscs2,  // viscous term (second derivatves) at GP at n
    const bool newton                                  // turn on/off linearization of convective term
)
{
//    /*
//                       /                      \
//                    - |  virt tau , eps(Dtau)  |
//                       \                      /
//     */
//
//    const double reciproke_viscfac = 1.0/(2.0*visc);
//    assembler.template Matrix<Sigmaxx,Sigmaxx>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
//    assembler.template Matrix<Sigmaxy,Sigmaxy>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
//    assembler.template Matrix<Sigmaxz,Sigmaxz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
//    assembler.template Matrix<Sigmayy,Sigmayy>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
//    assembler.template Matrix<Sigmayz,Sigmayz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
//    assembler.template Matrix<Sigmazz,Sigmazz>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
//
//    assembler.template Vector<Sigmaxx>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,0));
//    assembler.template Vector<Sigmaxy>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,1)*2.0);
//    assembler.template Vector<Sigmaxz>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,2)*2.0);
//    assembler.template Vector<Sigmayy>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,1));
//    assembler.template Vector<Sigmayz>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,2)*2.0);
//    assembler.template Vector<Sigmazz>(shp_tau,  reciproke_viscfac*timefacfac*tau(2,2));
//
//    /*             /                  \
//                    | virt tau , eps(Du) |
//                     \                  /
//     */
//    assembler.template Matrix<Sigmaxx,Velx>(shp_tau,     timefacfac    , shp.dx);
//    assembler.template Matrix<Sigmaxy,Velx>(shp_tau,     timefacfac    , shp.dy);
//    assembler.template Matrix<Sigmaxy,Vely>(shp_tau,     timefacfac    , shp.dx);
//    assembler.template Matrix<Sigmaxz,Velx>(shp_tau,     timefacfac    , shp.dz);
//    assembler.template Matrix<Sigmaxz,Velz>(shp_tau,     timefacfac    , shp.dx);
//    assembler.template Matrix<Sigmayy,Vely>(shp_tau,     timefacfac    , shp.dy);
//    assembler.template Matrix<Sigmayz,Vely>(shp_tau,     timefacfac    , shp.dz);
//    assembler.template Matrix<Sigmayz,Velz>(shp_tau,     timefacfac    , shp.dy);
//    assembler.template Matrix<Sigmazz,Velz>(shp_tau,     timefacfac    , shp.dz);
//
//    assembler.template Vector<Sigmaxx>(shp_tau,    - timefacfac*vderxy(0, 0));
//    assembler.template Vector<Sigmaxy>(shp_tau,    - timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
//    assembler.template Vector<Sigmaxz>(shp_tau,    - timefacfac*(vderxy(0, 2) + vderxy(2, 0)));
//    assembler.template Vector<Sigmayy>(shp_tau,    - timefacfac*vderxy(1, 1));
//    assembler.template Vector<Sigmayz>(shp_tau,    - timefacfac*(vderxy(1, 2) + vderxy(2, 1)));
//    assembler.template Vector<Sigmazz>(shp_tau,    - timefacfac*vderxy(2, 2));
//
//
//    /* pressure-pressure coupling, rectangular part */
//    /*
//                     /                    \
//                    |                      |
//                  - | tr(virt tau^e) , p I |
//                    |                      |
//                     \                    /
//     */
//    assembler.template Matrix<Sigmaxx,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
//    assembler.template Matrix<Sigmayy,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
//    assembler.template Matrix<Sigmazz,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
//
//    assembler.template Vector<Sigmaxx>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
//    assembler.template Vector<Sigmayy>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
//    assembler.template Vector<Sigmazz>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
}

/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class M2, class V2>
void SysmatDomainStress(
    const DRT::ELEMENTS::Combust3*      ele,          ///< element whose matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih, ///< information about the interface
    const XFEM::ElementDofManager&      dofman,       ///< dofmanager of this element
    const M1&                           evelnp,       ///< nodal velocity n+1 values for this element
    const M1&                           eveln,        ///< nodal velocity n   values for this element
    const M1&                           evelnm,       ///< nodal velocity n-1 values for this element
    const M1&                           eaccn,        ///< nodal acceleration n values for this element
    const V1&                           eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                           ephi,         ///< nodal G-function values for this element
    const M2&                           etau,         ///< element stresses (?) not used in combustion
    Teuchos::RCP<const MAT::Material>   material,     ///< fluid material
    const FLUID_TIMEINTTYPE             timealgo,     ///< time discretization scheme
    const double                        dt,           ///< delta t (time step size)
    const double                        theta,        ///< parameter for one step theta scheme
    const bool                          newton,       ///< full Newton or fixed-point-like
    const bool                          pstab,        ///< flag for stabilization
    const bool                          supg,         ///< flag for stabilization
    const bool                          cstab,        ///< flag for stabilization
    const INPAR::FLUID::TauType         tautype,      ///< stabilization parameter definition
    const bool                          instationary, ///< switch between stationary and instationary formulation
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - domain");

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration factor
  const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

  // get node coordinates of this element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // dead load in element nodes
  //LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;

  // figure out whether we have stress unknowns at all
  const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
  dsassert(tauele_unknowns_present == true, "there are no element stress unknowns!");
  //    const bool velocity_unknowns_present = (getNumParam<ASSTYPE>(dofman, Velx, 1) > 0);
  const bool pressure_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Pres, 1) > 0);
  dsassert(pressure_unknowns_present == true, "there is no element pressure unknown!");

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

  // stabilization parameter
  const double hk = FLD::UTILS::HK<DISTYPE>(evelnp,xyze);
  const double mk = FLD::UTILS::MK<DISTYPE>();

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));

  //----------------------------------------------------------------------------------------------
  // loop over domain integration cells
  //----------------------------------------------------------------------------------------------
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    //------------------------------------------------------------------------------------------
    // get material parameters for this integration cell
    //------------------------------------------------------------------------------------------
    int matid = 777;
    // check on which side of the interface the cell is located
    if(cell->getDomainPlus()) // cell belongs to burnt domain
    {
      matid = 3; // burnt material
    }
    else // cell belongs to unburnt domain
    {
      matid = 4; // unburnt material
    }
    // get material from list of materials
    Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
    // check if we really have a fluid material
    dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
    const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
    // get the kinematic viscosity \nu
    const double kinvisc = mat->Viscosity();
    // get the density \rho^{n+1}
    const double densnp = mat->Density();
    // density \rho^{n}
    // remark: until the XFEM time integration is done properly, densities from previous and
    // current time steps are identical
    const double densn = densnp;
    // compute dynamic viscosity \mu
    const double dynvisc = kinvisc * densnp;

    //------------------------------------------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //------------------------------------------------------------------------------------------
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);

    //------------------------------------------------------------------------------------------
    // get Gaussian points for this integration cell
    //------------------------------------------------------------------------------------------
    const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Intersected(), cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //------------------------------------------------------------------------------------------
    // integration loop over Gaussian points
    //------------------------------------------------------------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //--------------------------------------------------------------------------------------
      // transform coordinates of this Gaussian point
      //--------------------------------------------------------------------------------------
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
      //--------------------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      // discontinuous stress shape functions
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
      if (ASSTYPE == XFEM::xfem_assembly)
      {
        if (tauele_unknowns_present)
        {
          DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
        }
        else
        {
          funct_stress.Clear();
        }
      }

//--------------------------------------------------------------------------------------
      // procedures involving Jacobian matrix
      //--------------------------------------------------------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
      xjm.MultiplyNT(deriv,xyze);

      const double det = xjm.Determinant();
      const double fac = intpoints.qwgt[iquad]*det*detcell;
      //cout << "domain spatial factor " << fac << endl;

      if (det < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
      }

      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Invert(xjm);

      //--------------------------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //--------------------------------------------------------------------------------------
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

      //--------------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<shpVecSize> shp;
      // comment missing
      static LINALG::Matrix<shpVecSizeStress,1> shp_tau;
      if (ASSTYPE == XFEM::xfem_assembly)
      {
        //cout << "this is a fully or partially enriched element" << endl;
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSize,1> enr_funct;
        enr_funct.Clear();
        static LINALG::Matrix<3,shpVecSize> enr_derxy;
        enr_derxy.Clear();
        static LINALG::Matrix<6,shpVecSize> enr_derxy2;
        enr_derxy2.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            Velx,
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

        if (tauele_unknowns_present)
        {
          LINALG::Matrix<shpVecSizeStress,1> enr_funct_stress;

          // shape functions for element dofs
          enrvals.ComputeEnrichedElementShapefunction(
              Sigmaxx,
              funct_stress,
              enr_funct_stress);

          for (size_t iparam = 0; iparam < numparamtauxx; ++iparam)
          {
            shp_tau(iparam) = enr_funct_stress(iparam);
          }
        }
        else
        {
          shp_tau.Clear();
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
        if (tauele_unknowns_present)
        {
          dserror("no stress enrichments without xfem assembly");
        }
      }

      //--------------------------------------------------------------------------------------
      // get velocities and accelerations at integration point
      //--------------------------------------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpveln  = COMBUST::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpvelnm = COMBUST::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpaccn  = COMBUST::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);

      //--------------------------------------------------------------------------------------
      // get velocity derivatives at integration point
      //--------------------------------------------------------------------------------------
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

      //--------------------------------------------------------------------------------------
      // get second velocity derivatives at integration point
      //--------------------------------------------------------------------------------------
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

      //--------------------------------------------------------------------------------------
      // get pressure at integration point
      //--------------------------------------------------------------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp.d0(iparam);

      //--------------------------------------------------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp.dx(iparam);
        gradp(1) += eprenp(iparam) * shp.dy(iparam);
        gradp(2) += eprenp(iparam) * shp.dz(iparam);
      }

      // get viscous stress unknowns
      static LINALG::Matrix<nsd,nsd> tau;
      if (tauele_unknowns_present)
      {
        XFLUID::fill_tau(numparamtauxx, shp_tau, etau, tau);
      }
      else
      {
        tau.Clear();
      }

      //--------------------------------------------------------------------------------------
      // compute stabilization parameters (3 taus)
      //--------------------------------------------------------------------------------------
      const double vel_norm = gpvelnp.Norm2();
      const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_Mu = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;
      FLD::UTILS::computeStabilizationParams(gpvelnp, xji,
          instationary, dynvisc, densnp, vel_norm, strle, hk, mk, timefac, dt,
          INPAR::FLUID::tautype_franca_barrenechea_valentin_wall,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);

      //--------------------------------------------------------------------------------------
      // build RHS vector at integration point (scaled with timefac and density)
      //--------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> rhsint;
      rhsint.Clear();

      // time RHS vector from previous time step
      // remark: histvec = \rho^{n}*(u^n + delta t (1-theta) a^n)
      LINALG::Matrix<nsd,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
          gpveln, gpvelnm, gpaccn, timealgo, dt, theta);
      histvec.Scale(densn);

      // body force vector
      // remark: unit has to be acceleration; that is, not density weighted
      static LINALG::Matrix<nsd,1> bodyforce;
      bodyforce.Clear();
      bodyforce(0) = 0.0;
      bodyforce(1) = 0.0;
      bodyforce(2) = 0.0;
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densnp);
      //            // read body force from input file
      //            LINALG::SerialDenseMatrix edeadng = XFLUID::BodyForceTwoPhaseFlow<DISTYPE>(ele, 0.0);
      //            // remark: use BodyForce() instead?
      //            for (std::size_t isd = 0; isd < nsd; isd++)
      //            {
      //              for (std::size_t inode = 0; inode < numnode; inode++)
      //                  bodyforce(isd) += edeadng(isd,inode) * funct(inode);
      //            }

      for (size_t isd = 0; isd < nsd; ++isd)
        rhsint(isd) = histvec(isd) + timefac*bodyforce(isd);

      //--------------------------------------------------------------------------------------
      // compute numerical representation of some single operators (not scaled)
      //--------------------------------------------------------------------------------------

      // convective term:  u_old * grad u_old:
      // remark: conv_old = vderxy(i,j)*gpvelnp(j);
      LINALG::Matrix<nsd,1> conv_old;
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
        enr_viscs2.xy(iparam) = 0.5 *  shp.dxdy(iparam);
        enr_viscs2.xz(iparam) = 0.5 *  shp.dxdz(iparam);
        enr_viscs2.yx(iparam) = 0.5 *  shp.dydx(iparam);
        enr_viscs2.yy(iparam) = 0.5 * (shp.dxdx(iparam) + 2.0 * shp.dydy(iparam) + shp.dzdz(iparam));
        enr_viscs2.yz(iparam) = 0.5 *  shp.dydz(iparam);
        enr_viscs2.zx(iparam) = 0.5 *  shp.dzdx(iparam);
        enr_viscs2.zy(iparam) = 0.5 *  shp.dzdy(iparam);
        enr_viscs2.zz(iparam) = 0.5 * (shp.dxdx(iparam) + shp.dydy(iparam) + 2.0 * shp.dzdz(iparam));
      }

      //--------------------------------------------------------------------------------------
      // evaluate residual (scaled with timefac and density)
      // remark: used for stabilization right hand sides
      //--------------------------------------------------------------------------------------
      //            LINALG::Matrix<nsd,1> res_old(true);
      //            if (instationary)
      //            {
      //              // timefac * (\rho*du/dt + (\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              for (size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = densnp*gpvelnp(isd) - rhsint(isd) + timefac*(densnp*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));
      //            }
      //            else // stationary
      //            {
      //              // timefac * ((\rho u*grad) u + grad p - 2\mu\epsilon(u) - \rho g)
      //              // remark: Although calculation is stationary, scaling with timefac is kept to be
      //              //         consistent with formulation of matrix entries, which are also scaled.
      //              for(size_t isd = 0; isd < nsd; ++isd)
      //                res_old(isd) = timefac*(densnp*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd)-bodyforce(isd));
      //            }

      // evaluate residual once for all stabilisation right hand sides
      LINALG::Matrix<nsd,1> res_old(true);
      for (size_t isd = 0; isd < nsd; ++isd)
        res_old(isd) = -rhsint(isd)+timefac*(conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

      if (instationary)
      {
        res_old += gpvelnp;
      }

      // integration factors and coefficients of single terms
      const double timefacfac = timefac * fac;

      //--------------------------------------------------------------------------------------
      // build single stiffness matrix entries
      //--------------------------------------------------------------------------------------
      COMBUST::BuildDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
          assembler, shp, fac, timefac, timefacfac, densnp, dynvisc,
          gpvelnp, pres, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_viscs2,
          instationary, newton, pstab, supg, cstab,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);

      COMBUST::BuildDomainIntegralsStress<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
          assembler, shp, shp_tau, fac, timefac, timefacfac, densnp, dynvisc,
          gpvelnp, pres, gradp, vderxy, rhsint, res_old, visc_old, tau,
          enr_conv_c, enr_viscs2, newton);

    } // end loop over gauss points
  } // end loop over integration cells

  return;

}


/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2>
void SysmatBoundaryStress(
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M2&                         etau,
    Teuchos::RCP<const MAT::Material> material,      ///< fluid material
    const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
    const double&                     dt,            ///< delta t (time step size)
    const double&                     theta,         ///< factor for one step theta scheme
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const double                      flamespeed,
    const double                      nitschevel,
    const double                      nitschepres
)
{
#if 0
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
  const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  //------------------------------------------------------------------------------------------------
  // get material parameters for all boundary integration cell of this element
  //------------------------------------------------------------------------------------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif

  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // index plus == burnt material
  const int matid_plus = 3;
  // get material from list of materials
  Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid_plus);
  // check if we really have a fluid material
  dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
  // get the kinematic viscosity \nu
  const double kinvisc_plus = mat->Viscosity();
  // get the density \rho
  const double dens_plus = mat->Density();
  //cout << "plus density: " << dens_plus << endl;
  // compute dynamic viscosity \mu
  const double visc_plus = kinvisc_plus * dens_plus;

  // index minus = unburnt material
  const int matid_minus = 4;
  // get material from list of materials
  matptr = matlist->MaterialById(matid_minus);
  // check if we really have a fluid material
  dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
  // get the kinematic viscosity \nu
  const double kinvisc_minus = mat->Viscosity();
  // get the density \rho
  const double dens_minus = mat->Density();
  //cout << "minus density: " << dens_minus << endl;
  // compute dynamic viscosity \mu
  const double visc_minus = kinvisc_minus * dens_minus;

  //------------------------------------------------------------------------------------------------
  // compute jump values
  //------------------------------------------------------------------------------------------------
  // velocity jump value
  const double ju = -2.0; // -flamespeed*dens_minus*(1/dens_plus - 1/dens_minus);
  //cout << "J_u: " << ju << endl;
  // pressure jump value
  const double jp = 20.0; // (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);
  //cout << "J_p: " << jp << endl;

  // Nitsche parameter velocity
  const double alphau = nitschevel;
  //cout << "alpha_u: " << alphau << endl;
  // Nitsche parameter pressure
  const double alphap = nitschepres;
  //cout << "alpha_p: " << alphap << endl;

#if 0
  const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
  // for now, I don't try to compare to elements without stress unknowns, since they lock anyway
#endif

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
#if 0
  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, the numparam is determined by the dofmanager
  //const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
  const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
#endif

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

  //------------------------------------------------------------------------------------------
  // evaluate the enrichment function at the interface (boundary integration cells)
  //------------------------------------------------------------------------------------------
  const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
  const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);

  //------------------------------------------------------------------------------------------------
  // loop over boundary integration cells
  //------------------------------------------------------------------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    //----------------------------------------------------------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------------------------------------------------------
    // TODO: finish -> template? Only tri3 boundary integration cells allowed?
    if (cell->Shape() == DRT::Element::tri3)
    {
      // there should only be triangular boundary cells
    }
    else if (cell->Shape() == DRT::Element::quad4)
    {
      dserror("triangular boundary integration cell expected");
    }
    else
    {
      dserror("invalid type of boundary integration cell");
    }

    // here, a triangular boundary integration cell is assumed (numvertices = 3)
    const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

    // TODO: are 6 Gaussian points enough?
    //const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
    const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_37point);

    // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
    LINALG::SerialDenseMatrix cellXiDomaintmp = cell->CellNodalPosXiDomain();
    // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
    LINALG::SerialDenseMatrix cellXYZDomaintmp = cell->CellNodalPosXYZ();
    // transform to fixed size format
    const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
    const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

    //----------------------------------------------------------------------------------------------
    // integration loop over Gaussian points
    //----------------------------------------------------------------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //--------------------------------------------------------------------------------------------
      // transform coordinates of this Gaussian point
      //--------------------------------------------------------------------------------------------
      // coordinates of this integration point in boundary cell coordinates \eta^boundary
      const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);

      // coordinates of this integration point in element coordinates \xi^domain
      LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

      //--------------------------------------------------------------------------------------------
      // compute normal vector (in physical coordinates)
      // remark: for linear boundary integrastion cells this could be done before the loop over all
      //         Gaussian points
      // TODO: is this normal normed?
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> normal(true);
      normal.Clear();
      GEO::computeNormalToSurfaceElement(cell->Shape(), cellXYZDomain, posEtaBoundary, normal);

      //--------------------------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<numnode,1> funct;
      funct.Clear();
      static LINALG::Matrix<nsd,numnode> deriv;
      deriv.Clear();
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //--------------------------------------------------------------------------------------------
      // procedures involving Jacobian matrix for domain mapping
      //--------------------------------------------------------------------------------------------
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
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), detXtoXi);
      }
      //#endif
      // inverse of jacobian
      static LINALG::Matrix<nsd,nsd> xji;
      xji.Clear();
      xji.Invert(xjm);

      //--------------------------------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //--------------------------------------------------------------------------------------------
      // compute first global derivative
      static LINALG::Matrix<3,numnode> derxy;
      derxy.Clear();
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      // compute second global derivative
      static LINALG::Matrix<6,numnode> derxy2;
      if (false) // if (higher_order_ele)
      {
        static LINALG::Matrix<6,numnode> deriv2;
        DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
        DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

      //--------------------------------------------------------------------------------------------
      // compute Jacobian matrix for mapping to boundary
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,numvertices> deriv_boundary;
      deriv_boundary.Clear();
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posEtaBoundary(0),posEtaBoundary(1),cell->Shape());

      // get Jacobian matrix d \xi^domain / d \eta^boundary  (3x2)
      static LINALG::Matrix<nsd,2> dxideta;
      dxideta.Clear();
      for (std::size_t k=0; k!=numvertices; ++k)
        for (std::size_t i=0; i!=nsd; ++i)
          for (std::size_t j=0; j!=2; ++j)
            dxideta(i,j) += cellXiDomain(i,k)*deriv_boundary(j,k);

      // compute covariant metric tensor G (2x2)
      // metric = dxideta(k,i)*dxideta(k,j);
      static LINALG::Matrix<2,2> metric;
      metric.Clear();
      metric.MultiplyTN(dxideta,dxideta);
      // determinant for mapping from element (Xi^3D/domain) to integration cell (Eta^2D/boundary)
      const double detXitoEta = sqrt(metric.Determinant());
      //#ifdef DEBUG
      // actually ths makes no sence, because the square root was already taken before
      if (detXitoEta  < 0.0)
      {
        dserror("negative Jacobian determinant detXitoEta: %f!", detXitoEta );
      }
      //#endif

      //--------------------------------------------------------------------------------------------
      // compute integration factors
      //--------------------------------------------------------------------------------------------
      // compute spatial integration factor
      const double fac = intpoints.qwgt[iquad]*detXtoXi*detXitoEta;
      //cout << "boundary spatial factor " << fac << endl;
      // compute total (time and spatial) integration factor (and coefficients of single terms)?
      const double timefacfac = timefac * fac;
      //cout << "boundary time and spatial factor " << timefacfac << endl;

      //--------------------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------------------
      const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
      static LINALG::Matrix<shpVecSize,1>       enrfunct_plus;
      enrfunct_plus.Clear();
      static LINALG::Matrix<shpVecSize,1>       enrfunct_minus;
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
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
          enrfunct_plus, enrderxy_plus, enrderxy2_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
          enrfunct_minus, enrderxy_minus, enrderxy2_minus);

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<shpVecSize> shp_jump;
      static XFEM::ApproxFunc<shpVecSize> shp_mean;
      static XFEM::ApproxFunc<shpVecSize> shp_mean_visc;
      static XFEM::ApproxFunc<shpVecSize> shp_jump_visc;

      // fill approximation functions
      for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        shp_jump.d0(iparam) = enrfunct_plus(iparam) - enrfunct_minus(iparam);
        //cout << "shp_jump.d0: " << shp_jump.d0(iparam) << endl;
        shp_mean.d0(iparam) = 0.5*(enrfunct_plus(iparam) + enrfunct_minus(iparam));
        //cout << "shp_mean.d0: " << shp_mean.d0(iparam) << endl;
        shp_mean_visc.d0(iparam) = 0.5*(visc_plus*enrfunct_plus(iparam) + visc_minus*enrfunct_minus(iparam));
        shp_mean_visc.dx(iparam) = 0.5*(visc_plus*enrderxy_plus(0,iparam) + visc_minus*enrderxy_minus(0,iparam));
        shp_mean_visc.dy(iparam) = 0.5*(visc_plus*enrderxy_plus(1,iparam) + visc_minus*enrderxy_minus(1,iparam));
        shp_mean_visc.dz(iparam) = 0.5*(visc_plus*enrderxy_plus(2,iparam) + visc_minus*enrderxy_minus(2,iparam));
        shp_mean_visc.dxdx(iparam) = 0.5*(visc_plus*enrderxy2_plus(0,iparam) + visc_minus*enrderxy2_minus(0,iparam));
        shp_mean_visc.dxdy(iparam) = 0.5*(visc_plus*enrderxy2_plus(3,iparam) + visc_minus*enrderxy2_minus(3,iparam));
        shp_mean_visc.dxdz(iparam) = 0.5*(visc_plus*enrderxy2_plus(4,iparam) + visc_minus*enrderxy2_minus(4,iparam));
        shp_mean_visc.dydx(iparam) = shp_mean_visc.dxdy(iparam);
        shp_mean_visc.dydy(iparam) = 0.5*(visc_plus*enrderxy2_plus(1,iparam) + visc_minus*enrderxy2_minus(1,iparam));
        shp_mean_visc.dydz(iparam) = 0.5*(visc_plus*enrderxy2_plus(5,iparam) + visc_minus*enrderxy2_minus(5,iparam));
        shp_mean_visc.dzdx(iparam) = shp_mean_visc.dxdz(iparam);
        shp_mean_visc.dzdy(iparam) = shp_mean_visc.dydz(iparam);
        shp_mean_visc.dzdz(iparam) = 0.5*(visc_plus*enrderxy2_plus(2,iparam) + visc_minus*enrderxy2_minus(2,iparam));

        shp_jump_visc.d0(iparam) = visc_plus*enrfunct_plus(iparam) - visc_minus*enrfunct_minus(iparam);
        shp_jump_visc.dx(iparam) = visc_plus*enrderxy_plus(0,iparam) - visc_minus*enrderxy_minus(0,iparam);
        shp_jump_visc.dy(iparam) = visc_plus*enrderxy_plus(1,iparam) - visc_minus*enrderxy_minus(1,iparam);
        shp_jump_visc.dz(iparam) = visc_plus*enrderxy_plus(2,iparam) - visc_minus*enrderxy_minus(2,iparam);
        shp_jump_visc.dxdx(iparam) = visc_plus*enrderxy2_plus(0,iparam) - visc_minus*enrderxy2_minus(0,iparam);
        shp_jump_visc.dxdy(iparam) = visc_plus*enrderxy2_plus(3,iparam) - visc_minus*enrderxy2_minus(3,iparam);
        shp_jump_visc.dxdz(iparam) = visc_plus*enrderxy2_plus(4,iparam) - visc_minus*enrderxy2_minus(4,iparam);
        shp_jump_visc.dydx(iparam) = shp_jump_visc.dxdy(iparam);
        shp_jump_visc.dydy(iparam) = visc_plus*enrderxy2_plus(1,iparam) - visc_minus*enrderxy2_minus(1,iparam);
        shp_jump_visc.dydz(iparam) = visc_plus*enrderxy2_plus(5,iparam) - visc_minus*enrderxy2_minus(5,iparam);
        shp_jump_visc.dzdx(iparam) = shp_jump_visc.dxdz(iparam);
        shp_jump_visc.dzdy(iparam) = shp_jump_visc.dydz(iparam);
        shp_jump_visc.dzdz(iparam) = visc_plus*enrderxy2_plus(2,iparam) - visc_minus*enrderxy2_minus(2,iparam);
      }
      //cout << "mean enrichment function: " << shp_mean_visc.d0 << endl;
      //--------------------------------------------------------------------------------------------
      // compute data at Gaussian point for rhs
      //--------------------------------------------------------------------------------------------

      // velocity jump
      static LINALG::Matrix<nsd,1> vjump(true);
      vjump.Clear();
      vjump = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);
      //cout << "evelnp: " << evelnp << endl;
      //cout << "shp_jump.d0: " << shp_jump.d0 << endl;
      //cout << "vjump: " << vjump << endl;

      // mean velocity
      static LINALG::Matrix<nsd,1> vmean(true);
      vmean.Clear();
      vmean = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shp_mean.d0, numparamvelx);
      //cout << "vmean: " << vmean << endl;

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
      //cout << "vderxy_mean_visc: " << vderxy_mean_visc << endl;

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
      //cout << "vderxy_mean_visc: " << vderxy_mean_visc << endl;

      // get pressure jump
      static double pjump;
      pjump = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pjump += shp_jump.d0(iparam)*eprenp(iparam);
      //cout << "pjump: " << pjump << endl;

      static double pmean;
      pmean = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pmean += shp_mean.d0(iparam)*eprenp(iparam);
      //cout << "pmean: " << pmean << endl;

      //            Epetra_SerialDenseVector shp_iface(numnode_boundary*begids.size());
      //            int pos = 0;
      //            for (std::set<int>::const_iterator begid = begids.begin(); begid != begids.end();++begid)
      //            {
      //              if (*begid == boundaryele->Id())
      //              {
      //                for (std::size_t inode=0; inode < numnode_boundary; ++inode)
      //                {
      //                  shp_iface(pos+inode) = funct_boundary(inode);
      //                }
      //                break;
      //              }
      //              pos += numnode_boundary;
      //            }
      //
      //            // get velocities (n+g,i) at integration point
      //            // gpvelnp = evelnp(i,j)*shp(j);
      //            LINALG::Matrix<nsd,1> gpvelnp;
      //            gpvelnp.Clear();
      //            for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      //                for (std::size_t isd = 0; isd < nsd; ++isd)
      //                    gpvelnp(isd) += evelnp(isd,iparam)*shp(iparam);
      //
      //            // get interface velocity
      //            LINALG::Matrix<nsd,1> interface_gpvelnp;
      //            interface_gpvelnp.Clear();
      //            for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
      //                for (std::size_t isd = 0; isd < nsd; ++isd)
      //                    interface_gpvelnp(isd) += vel_boundary(isd,inode)*funct_boundary(inode);
      //
      //            // get viscous stress unknowns
      //            static LINALG::Matrix<nsd,nsd> tau;
      //            XFLUID::fill_tau(numparamtauxx, shp_tau, etau, tau);



      //--------------------------------------------------------------------------------------------
      // build single boundary integral stiffness terms
      //
      // remarks: || x || stands for the jump operator ( || x || = x^+ - x^- )
      //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
      //             n    stands for the normal vector on the interface pointing into the
      //                  burnt domain (n = n^- = -n^+)
      //--------------------------------------------------------------------------------------------

#ifdef COMBUST_VISCCONSISTMEAN
      //-------------------------    |                                      |
      // viscous consistency term  + |  || v || , < 2\mu epsilon( Du ) > n  |
      //-------------------------    |                                      |

      assembler.template Matrix<Velx,Velx>(shp_jump.d0, 2.0*timefacfac*normal(0), shp_mean_visc.dx);
      assembler.template Matrix<Velx,Velx>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dy);
      assembler.template Matrix<Velx,Velx>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dz);
      assembler.template Matrix<Velx,Vely>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dx);
      assembler.template Matrix<Velx,Velz>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dx);

      assembler.template Matrix<Vely,Velx>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dy);
      assembler.template Matrix<Vely,Vely>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dx);
      assembler.template Matrix<Vely,Vely>(shp_jump.d0, 2.0*timefacfac*normal(1), shp_mean_visc.dy);
      assembler.template Matrix<Vely,Vely>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dz);
      assembler.template Matrix<Vely,Velz>(shp_jump.d0,     timefacfac*normal(2), shp_mean_visc.dy);

      assembler.template Matrix<Velz,Velx>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dz);
      assembler.template Matrix<Velz,Vely>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dz);
      assembler.template Matrix<Velz,Velz>(shp_jump.d0,     timefacfac*normal(0), shp_mean_visc.dx);
      assembler.template Matrix<Velz,Velz>(shp_jump.d0,     timefacfac*normal(1), shp_mean_visc.dy);
      assembler.template Matrix<Velz,Velz>(shp_jump.d0, 2.0*timefacfac*normal(2), shp_mean_visc.dz);

      //   |                                       |
      // - |  || v || , < 2\mu epsilon( u_i ) > n  |
      //   |                                       |

      assembler.template Vector<Velx>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(0,0) + vderxy_mean_visc(0,0))*normal(0));
      assembler.template Vector<Velx>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(0,1) + vderxy_mean_visc(1,0))*normal(1));
      assembler.template Vector<Velx>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(0,2) + vderxy_mean_visc(2,0))*normal(2));

      assembler.template Vector<Vely>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(1,0) + vderxy_mean_visc(0,1))*normal(0));
      assembler.template Vector<Vely>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(1,1) + vderxy_mean_visc(1,1))*normal(1));
      assembler.template Vector<Vely>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(1,2) + vderxy_mean_visc(2,1))*normal(2));

      assembler.template Vector<Velz>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(2,0) + vderxy_mean_visc(0,2))*normal(0));
      assembler.template Vector<Velz>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(2,1) + vderxy_mean_visc(1,2))*normal(1));
      assembler.template Vector<Velz>(shp_jump.d0, -timefacfac* (vderxy_mean_visc(2,2) + vderxy_mean_visc(2,2))*normal(2));
#endif

#ifdef COMBUST_VISCCONSISTJUMP
      //-------------------------    |                                      |
      // viscous consistency term  + |  < v > , || 2\mu epsilon( Du ) || n  |
      //-------------------------    |                                      |

      assembler.template Matrix<Velx,Velx>(shp_mean.d0, 2.0*timefacfac*normal(0), shp_jump_visc.dx);
      assembler.template Matrix<Velx,Velx>(shp_mean.d0,     timefacfac*normal(1), shp_jump_visc.dy);
      assembler.template Matrix<Velx,Velx>(shp_mean.d0,     timefacfac*normal(2), shp_jump_visc.dz);
      assembler.template Matrix<Velx,Vely>(shp_mean.d0,     timefacfac*normal(1), shp_jump_visc.dx);
      assembler.template Matrix<Velx,Velz>(shp_mean.d0,     timefacfac*normal(2), shp_jump_visc.dx);

      assembler.template Matrix<Vely,Velx>(shp_mean.d0,     timefacfac*normal(0), shp_jump_visc.dy);
      assembler.template Matrix<Vely,Vely>(shp_mean.d0,     timefacfac*normal(0), shp_jump_visc.dx);
      assembler.template Matrix<Vely,Vely>(shp_mean.d0, 2.0*timefacfac*normal(1), shp_jump_visc.dy);
      assembler.template Matrix<Vely,Vely>(shp_mean.d0,     timefacfac*normal(2), shp_jump_visc.dz);
      assembler.template Matrix<Vely,Velz>(shp_mean.d0,     timefacfac*normal(2), shp_jump_visc.dy);

      assembler.template Matrix<Velz,Velx>(shp_mean.d0,     timefacfac*normal(0), shp_jump_visc.dz);
      assembler.template Matrix<Velz,Vely>(shp_mean.d0,     timefacfac*normal(1), shp_jump_visc.dz);
      assembler.template Matrix<Velz,Velz>(shp_mean.d0,     timefacfac*normal(0), shp_jump_visc.dx);
      assembler.template Matrix<Velz,Velz>(shp_mean.d0,     timefacfac*normal(1), shp_jump_visc.dy);
      assembler.template Matrix<Velz,Velz>(shp_mean.d0, 2.0*timefacfac*normal(2), shp_jump_visc.dz);

      //   |                                       |
      // - |  < v > , || 2\mu epsilon( u_i ) || n  |
      //   |                                       |

      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(0,0) + vderxy_jump_visc(0,0))*normal(0));
      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(0,1) + vderxy_jump_visc(1,0))*normal(1));
      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(0,2) + vderxy_jump_visc(2,0))*normal(2));

      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(1,0) + vderxy_jump_visc(0,1))*normal(0));
      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(1,1) + vderxy_jump_visc(1,1))*normal(1));
      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(1,2) + vderxy_jump_visc(2,1))*normal(2));

      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(2,0) + vderxy_jump_visc(0,2))*normal(0));
      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(2,1) + vderxy_jump_visc(1,2))*normal(1));
      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac* (vderxy_jump_visc(2,2) + vderxy_jump_visc(2,2))*normal(2));
#endif

#ifdef COMBUST_VISCADJCONSIST
      //---------------------------------    |                                     |
      // viscous adjoint consistency term  + |  < 2\mu epsilon( v ) > n, || Du ||  |
      //---------------------------------    |                                     |

      assembler.template Matrix<Velx,Velx>(shp_mean_visc.dx, -2.0*timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Velx,Velx>(shp_mean_visc.dy, -    timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Velx,Velx>(shp_mean_visc.dz, -    timefacfac*normal(2), shp_jump.d0);
      assembler.template Matrix<Vely,Velx>(shp_mean_visc.dx, -    timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Velz,Velx>(shp_mean_visc.dx, -    timefacfac*normal(2), shp_jump.d0);

      assembler.template Matrix<Velx,Vely>(shp_mean_visc.dy, -    timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_mean_visc.dx, -    timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_mean_visc.dy, -2.0*timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_mean_visc.dz, -    timefacfac*normal(2), shp_jump.d0);
      assembler.template Matrix<Velz,Vely>(shp_mean_visc.dy, -    timefacfac*normal(2), shp_jump.d0);

      assembler.template Matrix<Velx,Velz>(shp_mean_visc.dz, -    timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Vely,Velz>(shp_mean_visc.dz, -    timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_mean_visc.dx, -    timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_mean_visc.dy, -    timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_mean_visc.dz, -2.0*timefacfac*normal(2), shp_jump.d0);

      //   |                                      |
      // - |  < 2\mu epsilon( v ) > n, || u_i ||  |
      //   |                                      |

      assembler.template Vector<Velx>(shp_mean_visc.dx, +2.0* normal(0)           *timefacfac*vjump(0,0));
      assembler.template Vector<Velx>(shp_mean_visc.dy, +    (normal(0)+normal(1))*timefacfac*vjump(0,0));
      assembler.template Vector<Velx>(shp_mean_visc.dz, +    (normal(0)+normal(2))*timefacfac*vjump(0,0));

      assembler.template Vector<Vely>(shp_mean_visc.dx, +    (normal(1)+normal(0))*timefacfac*vjump(1,0));
      assembler.template Vector<Vely>(shp_mean_visc.dy, +2.0* normal(1)           *timefacfac*vjump(1,0));
      assembler.template Vector<Vely>(shp_mean_visc.dz, +    (normal(1)+normal(2))*timefacfac*vjump(1,0));

      assembler.template Vector<Velz>(shp_mean_visc.dx, +    (normal(2)+normal(0))*timefacfac*vjump(2,0));
      assembler.template Vector<Velz>(shp_mean_visc.dy, +    (normal(2)+normal(1))*timefacfac*vjump(2,0));
      assembler.template Vector<Velz>(shp_mean_visc.dz, +2.0* normal(2)           *timefacfac*vjump(2,0));

      //-------------------------------------    |                                    |
      // viscous adjoint consistency term RHS  + |  < 2\mu\rho epsilon( v ) > n, J_u  |
      //-------------------------------------    |                                    |

      assembler.template Vector<Velx>(shp_mean_visc.dx, +2.0* normal(0)           *timefacfac*normal(0)*ju);
      assembler.template Vector<Velx>(shp_mean_visc.dy, +    (normal(0)+normal(1))*timefacfac*normal(0)*ju);
      assembler.template Vector<Velx>(shp_mean_visc.dz, +    (normal(0)+normal(2))*timefacfac*normal(0)*ju);

      assembler.template Vector<Vely>(shp_mean_visc.dx, +    (normal(1)+normal(0))*timefacfac*normal(1)*ju);
      assembler.template Vector<Vely>(shp_mean_visc.dy, +2.0* normal(1)           *timefacfac*normal(1)*ju);
      assembler.template Vector<Vely>(shp_mean_visc.dz, +    (normal(1)+normal(2))*timefacfac*normal(1)*ju);

      assembler.template Vector<Velz>(shp_mean_visc.dx, +    (normal(2)+normal(0))*timefacfac*normal(2)*ju);
      assembler.template Vector<Velz>(shp_mean_visc.dy, +    (normal(2)+normal(1))*timefacfac*normal(2)*ju);
      assembler.template Vector<Velz>(shp_mean_visc.dz, +2.0* normal(2)           *timefacfac*normal(2)*ju);
#endif

#ifdef COMBUST_PRESCONSISTMEAN
      //--------------------------    |                     |
      // pressure consistency term  - |  || v ||, < Dp > n  |
      //--------------------------    |                     |

      assembler.template Matrix<Velx,Pres>(shp_jump.d0, +timefacfac*normal(0), shp_mean.d0);
      assembler.template Matrix<Vely,Pres>(shp_jump.d0, +timefacfac*normal(1), shp_mean.d0);
      assembler.template Matrix<Velz,Pres>(shp_jump.d0, +timefacfac*normal(2), shp_mean.d0);

      //   |                      |
      // + |  || v ||, < p_i > n  |
      //   |                      |

      assembler.template Vector<Velx>(shp_jump.d0, -timefacfac*normal(0)*pmean);
      assembler.template Vector<Vely>(shp_jump.d0, -timefacfac*normal(1)*pmean);
      assembler.template Vector<Velz>(shp_jump.d0, -timefacfac*normal(2)*pmean);
#endif

#ifdef COMBUST_PRESCONSISTJUMP
      //--------------------------    |                     |
      // pressure consistency term  - |  < v >, || Dp || n  |
      //--------------------------    |                     |

      assembler.template Matrix<Velx,Pres>(shp_mean.d0, -timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Vely,Pres>(shp_mean.d0, -timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Velz,Pres>(shp_mean.d0, -timefacfac*normal(2), shp_jump.d0);

      //   |                      |
      // + |  < v >, || p_i || n  |
      //   |                      |

      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(0)*pjump);
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(1)*pjump);
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(2)*pjump);
#endif

#ifdef COMBUST_PRESADJCONSIST
      //----------------------------------    |                     |
      // pressure adjoint consistency term  + |  < q > n, || Du ||  |
      //----------------------------------    |                     |

      assembler.template Matrix<Pres,Velx>(shp_mean.d0, +timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Pres,Vely>(shp_mean.d0, +timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Pres,Velz>(shp_mean.d0, +timefacfac*normal(2), shp_jump.d0);

      //   |                      |
      // - |  < q > n, || u_i ||  |
      //   |                      |

      assembler.template Vector<Pres>(shp_mean.d0, -timefacfac*(vmean(0,0)*normal(0)
          +vmean(1,0)*normal(1)
          +vmean(2,0)*normal(2)));

      //--------------------------------------    |                |
      // pressure adjoint consistency term RHS  + |  < q > n, J_u  |
      //--------------------------------------    |                |

      assembler.template Vector<Pres>(shp_mean.d0, -timefacfac*ju);
#endif

#ifdef COMBUST_PRESJUMPRHS
      //------------------------------------------    |                |
      // pressure/flux jump (consistency) term RHS  + |  < v >, J_p*n  |
      //------------------------------------------    |                |

      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(0)*jp);
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(1)*jp);
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(2)*jp);
#endif

#ifdef COMBUST_NITSCHEVEL
      //------------------------    |                              |
      // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
      //------------------------    |                              |

      assembler.template Matrix<Velx,Velx>(shp_jump.d0, +alphau*timefacfac, shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_jump.d0, +alphau*timefacfac, shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_jump.d0, +alphau*timefacfac, shp_jump.d0);

      //    |                               |
      //  - |  \alpha_u || v ||, || u_i ||  |
      //    |                               |

      assembler.template Vector<Velx>(shp_jump.d0, -alphau*timefacfac*vjump(0,0));
      assembler.template Vector<Vely>(shp_jump.d0, -alphau*timefacfac*vjump(1,0));
      assembler.template Vector<Velz>(shp_jump.d0, -alphau*timefacfac*vjump(2,0));

      //----------------------------    |                         |
      // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, J_u  |
      //----------------------------    |                         |

      assembler.template Vector<Velx>(shp_jump.d0, +alphau*timefacfac*normal(0)*ju);
      assembler.template Vector<Vely>(shp_jump.d0, +alphau*timefacfac*normal(1)*ju);
      assembler.template Vector<Velz>(shp_jump.d0, +alphau*timefacfac*normal(2)*ju);
#endif

#ifdef COMBUST_NITSCHEPRES
      //------------------------    |                              |
      // Nitsche term (pressure)  + |  \alpha_p || q ||, || Dp ||  |
      //------------------------    |                              |

      assembler.template Matrix<Pres,Pres>(shp_jump.d0, alphap*timefacfac, shp_jump.d0);

      //    |                               |
      //  - |  \alpha_p || q ||, || p_i ||  |
      //    |                               |

      assembler.template Vector<Pres>(shp_jump.d0, -alphap*timefacfac*pjump);

      //----------------------------    |                         |
      // Nitsche term (pressure) RHS  + |  \alpha_p || q ||, J_p  |
      //----------------------------    |                         |

      assembler.template Vector<Pres>(shp_jump.d0, alphap*timefacfac*jp);
#endif

    } // loop Gaussian points
  } // loop boundary integration cells
#endif
  return;
}


} // namespace COMBUST

#endif
#endif
