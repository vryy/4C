/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat_premixed_nitsche.H

\brief system matrix formulation with Nitsche's method
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

#include "combust3_sysmat_standard.H"
#include "combust3_sysmat_surfacetension.H"
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
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class V2>
void SysmatDomainNitsche(
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const M1&                      evelnp,       ///< nodal velocity n+1 values for this element
    const M1&                      eveln,        ///< nodal velocity n   values for this element
    const M1&                      evelnm,       ///< nodal velocity n-1 values for this element
    const M1&                      eaccn,        ///< nodal acceleration n values for this element
    const V1&                      eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                      ephi,         ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material> material,  ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,     ///< time discretization scheme
    const double                   dt,           ///< delta t (time step size)
    const double                   theta,        ///< parameter for one step theta scheme
    const bool                     newton,       ///< full Newton or fixed-point-like
    const bool                     pstab,        ///< turn on/off PSPG stabilization
    const bool                     supg,         ///< turn on/off SUPG stabilization
    const bool                     cstab,        ///< turn on/off BV   stabilization
    const INPAR::FLUID::TauType    tautype,      //< stabilization parameter definition
    const bool                     instationary, ///< switch between stationary and instationary formulation
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
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
  const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

  // get node coordinates of this element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // read body force (dead load) from input file and apply it to element nodes
  LINALG::Matrix<nsd,numnode> edead(true);
  COMBUST::BodyForce<DISTYPE>(ele, 0.0, edead);

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // stabilization parameter
  const double hk = FLD::UTILS::HK<DISTYPE>(evelnp,xyze);
  const double mk = FLD::UTILS::MK<DISTYPE>();

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));

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
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Intersected(), cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
#ifdef BENEDIKT
      const LINALG::Matrix<3,1> pos_eta_domain(intpoints.qxg[iquad]);

      double det_cell_EtaToXYZ = 0.0;
      static LINALG::Matrix<numnode,1> funct3Dele;
      static LINALG::Matrix<3,numnode> deriv3Dele_xyz;
      static LINALG::Matrix<6,numnode> deriv3Dele_xyz_2nd;

      static LINALG::Matrix<3,3> xji_ele_XiToXYZ;

      //        static LINALG::Matrix<3,1> posXiDomain;

      switch (cell->Shape())
      {
      case DRT::Element::tet4:
        COMBUST::GetDomainIntegrationFeatures<DISTYPE,DRT::Element::tet4,ASSTYPE>(
            ele, pos_eta_domain,xyze, *cell, det_cell_EtaToXYZ, funct3Dele, deriv3Dele_xyz, deriv3Dele_xyz_2nd,xji_ele_XiToXYZ
        );
        break;
      case DRT::Element::hex8:
        COMBUST::GetDomainIntegrationFeatures<DISTYPE,DRT::Element::hex8,ASSTYPE>(
            ele, pos_eta_domain,xyze, *cell, det_cell_EtaToXYZ, funct3Dele, deriv3Dele_xyz, deriv3Dele_xyz_2nd,xji_ele_XiToXYZ
        );
        break;
      default:
        dserror("cell distype not implemented yet");
      }

      double fac = 0.0;

      fac = intpoints.qwgt[iquad] * det_cell_EtaToXYZ;


      if (det_cell_EtaToXYZ < 0.0)
      {
        dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det_cell_EtaToXYZ);
      }


      //=================================================================================
      //      8. Get enriched shape functions and derivatives as approximation functions
      //=================================================================================
      //--------------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;
#if 0
      // comment missing
      static LINALG::Matrix<shpVecSizeStress,1> shp_epsilon;
#endif
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
            funct3Dele,
            deriv3Dele_xyz,
            deriv3Dele_xyz_2nd,
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
          shp.d0(iparam) = funct3Dele(iparam);
          shp.dx(iparam) = deriv3Dele_xyz(0,iparam);
          shp.dy(iparam) = deriv3Dele_xyz(1,iparam);
          shp.dz(iparam) = deriv3Dele_xyz(2,iparam);
          shp.dxdx(iparam) = deriv3Dele_xyz_2nd(0,iparam);
          shp.dxdy(iparam) = deriv3Dele_xyz_2nd(3,iparam);
          shp.dxdz(iparam) = deriv3Dele_xyz_2nd(4,iparam);
          shp.dydx(iparam) = shp.dxdy(iparam);
          shp.dydy(iparam) = deriv3Dele_xyz_2nd(1,iparam);
          shp.dydz(iparam) = deriv3Dele_xyz_2nd(5,iparam);
          shp.dzdx(iparam) = shp.dxdz(iparam);
          shp.dzdy(iparam) = shp.dydz(iparam);
          shp.dzdz(iparam) = deriv3Dele_xyz_2nd(2,iparam);
        }
      }

      //cout << shp.d0 << endl << endl;
#else
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
      //cout << "domain spatial factor " << fac << endl;

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
        DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
      }
      else
      {
        derxy2.Clear();
      }

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;

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
#endif //BENEDIKT

#ifdef TPF_NITSCHE_LPS_DOMAIN        
      //=======================================additional pressure stabilization======================

      // calculate
      double meas_ele_plus = 0.0;   // measure of element ele at plus side
      double meas_ele_minus = 0.0;  // measure of element ele at minus side

      //--------------------------------------------------------------------------------------------
      // shape functions for pressure projection as approximation functions
      // (constant values for each element! -> before cell loop)
      //--------------------------------------------------------------------------------------------
      const std::size_t shpPresProjVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // temporary arrays holding enriched shape functions (int^+-(N * \Psi)) on either side of the interface
      static LINALG::Matrix<shpPresProjVecSize,1>       shp_presproj_plus;
      shp_presproj_plus.Clear();
      static LINALG::Matrix<shpPresProjVecSize,1>       shp_presproj_minus;
      shp_presproj_minus.Clear();


      //  TWOPHASEFLOWJUMP::BuildPresProjShapefunction<DISTYPE,ASSTYPE,shpPresProjVecSize,NUMDOF>(meas_ele_plus, meas_ele_minus,
      //        shp_presproj_plus, shp_presproj_minus,
      //          ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, eepsilon,
      //          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);
      COMBUST::BuildPresProjShapefunction<DISTYPE,ASSTYPE,shpPresProjVecSize,NUMDOF>(
          meas_ele_plus, meas_ele_minus,
          shp_presproj_plus, shp_presproj_minus,
          ele,
          ih,
          dofman,
          ephi,
          material);
      //  cout << "meas_ele_plus:\t" << meas_ele_plus << endl;
      //  cout << "meas_ele_minus:\t" << meas_ele_minus << endl << endl;
      //
      //  for  (int iparam = 0; iparam< (int)shpPresProjVecSize; iparam++){
      //    cout << iparam << "\t" << shp_presproj_plus(iparam) << endl;
      //  }
      //
      //  cout << endl;
      //
      //  for  (int iparam = 0; iparam< (int)shpPresProjVecSize; iparam++){
      //    cout << iparam << "\t" << shp_presproj_minus(iparam) << endl;
      //  }
      //
      //  cout << endl;

      static LINALG::Matrix<shpPresProjVecSize,1>   shp_presproj_jump;
      shp_presproj_jump.Clear();

      //  static LINALG::Matrix<shpPresProjVecSize,1>  shp_presproj_mean;
      //  shp_presproj_mean.Clear();

      // TODO this is not as secure as the material access is treated elsewhere
      // get material list for this element
      const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
      // set default id in list of materials
      int matid = -1;
      // check on which side of the interface the cell is located
      matid = matlist->MatID(0); // burnt material (first material in material list)

      for (int iparam = 0; iparam < (int)shpPresProjVecSize; iparam++)
      {
        if(matid==3) shp_presproj_jump(iparam) = shp.d0(iparam) - shp_presproj_plus(iparam);
        else
          shp_presproj_jump(iparam) = shp.d0(iparam) - shp_presproj_minus(iparam);

        //shp_presproj_mean(iparam) = 0.5/meas_ele_plus * shp_presproj_plus(iparam)+ 0.5/meas_ele_minus * shp_presproj_minus(iparam);
      }
      // for each element get shp_presstab = shp_jump.d0 - shp_presproj_jump
#endif

      // TODO aufpassen!
      ele_meas_plus = 0.5;
      ele_meas_minus = 0.5;

      //------------------------------------------------------
      // get velocities and accelerations at integration point
      //------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);

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
      const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_Mu = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;

#ifdef BENEDIKT
      /*
       * author schott
       * May 7, 2010
       *
       * take BAZILEVS stabilization!!!!
       */
      FLD::UTILS::computeStabilizationParams(gpvelnp, xji_ele_XiToXYZ,
          instationary, dynvisc, densaf, vel_norm, strle, hk, mk, timefac, dt,
          INPAR::FLUID::tautype_bazilevs,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
      //     FLD::UTILS::computeStabilizationParams(gpvelnp, xji_ele_XiToXYZ,
      //               instationary, dynvisc, densnp, vel_norm, strle, hk, mk, timefac, dt,
      //               INPAR::FLUID::tautype_franca_barrenechea_valentin_wall,
      //               tau_stab_Mu, tau_stab_Mp, tau_stab_C);
// schott Jun 18, 2010
      // stabilization parameters checked
      // dimensions:
      // [tau_stab_Mp] = s/[rho]
      // [tau_stab_Mu] = s/[rho]
      // [tau_stab_C]  = [rho]*m*m/s

#else
      FLD::UTILS::computeStabilizationParams(gpvelnp, xji,
          instationary, dynvisc, densaf, vel_norm, strle, hk, mk, timefac, dt,
          INPAR::FLUID::tautype_bazilevs,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
#endif

      //------------------------------------------------------------------------
      // build RHS vector at integration point (scaled with timefac and density)
      //------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> rhsint;
      rhsint.Clear();

      // time RHS vector from previous time step
      // remark: histvec = \rho^{n}*(u^n + delta t (1-theta) a^n)
      LINALG::Matrix<nsd,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
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
#ifdef BENEDIKT
          bodyforce(isd) += edead(isd,inode) * funct3Dele(inode);
#else
          bodyforce(isd) += edead(isd,inode) * funct(inode);
#endif
        }
      }
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densaf);

      for (size_t isd = 0; isd < nsd; ++isd)
        rhsint(isd) = histvec(isd) + timefac*bodyforce(isd);

      //-----------------------------------------------------------------------
      // compute numerical representation of some single operators (not scaled)
      //-----------------------------------------------------------------------

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
      for (size_t isd = 0; isd < nsd; ++isd)
        res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

      if (instationary)
      {
        for(size_t isd=0; isd < nsd; ++isd)
          res_old(isd) += densaf * gpvelnp(isd);
      }

      // integration factors and coefficients of single terms
      const double timefacfac = timefac * fac;

#ifdef TPF_NITSCHE_LPS_DOMAIN
      //--------------------------------------------------------------------------------------
      // build single stiffness matrix entries
      //--------------------------------------------------------------------------------------
      COMBUST::BuildDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress,shpPresProjVecSize>(
          assembler, shp,shp_presproj_jump, fac, timefac, timefacfac, densaf, dynvisc,
          gpvelnp, pres, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_viscs2,
          instationary, newton, pstab, supg, cstab,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
#else

      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
          assembler, shp, fac, timefac, timefacfac, densaf, dynvisc,
          gpvelnp, pres, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_viscs2,
          instationary, newton, pstab, supg, cstab,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
#endif
    } // end loop over gauss points
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
void Nitsche_SysmatBoundaryCell(
    GEO::BoundaryIntCell              cell,
    const DRT::UTILS::GaussRule2D&    intrule2D,
    M2&                               xyze,
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const bool&                       higher_order_ele,
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M3&                         grad_ephi_smoothed, ///<smoothed nodal G-function gradient values for this element
    const double&                     ju_n,
    const double&                     jflux_n,
    const double&                       timefac,
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const XFEM::ElementEnrichmentValues&       enrvals_plus,
    const XFEM::ElementEnrichmentValues&       enrvals_minus,
    const double&                        kappa_plus,    /// weights of average in plus domain
    const double&                        kappa_minus,    /// weights of average in minus domain
    const double&                        visc_plus,    /// viscosity in plus domain
    const double&                        visc_minus,   /// viscosity in minus domain
    const double&                        alphau,       /// Nitsche's stabilization parameter
    const size_t&                              numparamvelx,
    const size_t&                              numparampres,
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType     veljumptype,
    const INPAR::COMBUST::NormalTensionJumpType normaltensionjumptype
)
{
  const size_t nsd =3;

  // number of nodes of element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<CELLDISTYPE>::numNodePerElement;

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
      && ( normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_SurfaceTension)
      && (        surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami
          || surftensapprox == INPAR::COMBUST::surface_tension_approx_laplacebeltrami_smoothed )
  ) // only for laplace_beltrami approximation
  {
    COMBUST::RHS_LaplaceBeltrami_1DLineIntegrals_Nitsche<DISTYPE,CELLDISTYPE,ASSTYPE,NUMDOF>(
        cellXYZDomain, xyze, grad_ephi_smoothed,
        enrvals_plus,enrvals_minus,
        kappa_plus,kappa_minus,numparamvelx,timefac,surftensapprox,surftenscoeff,assembler);
  }

  //--------------------------------------
  // integration loop over Gaussian points
  //--------------------------------------
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
#ifdef BENEDIKT
    // 1. cell coordinates are stored in cellXYZDomain <nsd,numvertices>
    // 2. 2D-shape functions -> evaluate at each integration point

    // coordinates of this integration point in boundary cell coordinates \eta^boundary
    const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);
    double detXYZtoEta1Eta2 = 0.0;

    static LINALG::Matrix<2,numvertices> deriv_boundary;
    DRT::UTILS::shape_function_2D_deriv1(deriv_boundary,posEtaBoundary(0,0),posEtaBoundary(1,0),CELLDISTYPE);


    // det for surface integrals  ||  dxyz/deta1 x dxyz/deta2  ||
    static LINALG::Matrix<3,2> dXYZdEta1Eta2;
    dXYZdEta1Eta2.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < (int)numvertices; k++)
          dXYZdEta1Eta2(i,j) += cellXYZDomain(i,k)*deriv_boundary(j,k);

    static LINALG::Matrix<3,1> tmp_vec;     // vector for dxyz/deta1 x dxyz/deta2
    tmp_vec.Clear();
    tmp_vec(0) = dXYZdEta1Eta2(1,0) * dXYZdEta1Eta2(2,1) - dXYZdEta1Eta2(2,0) * dXYZdEta1Eta2(1,1);
    tmp_vec(1) = dXYZdEta1Eta2(2,0) * dXYZdEta1Eta2(0,1) - dXYZdEta1Eta2(0,0) * dXYZdEta1Eta2(2,1);
    tmp_vec(2) = dXYZdEta1Eta2(0,0) * dXYZdEta1Eta2(1,1) - dXYZdEta1Eta2(1,0) * dXYZdEta1Eta2(0,1);

    detXYZtoEta1Eta2 = tmp_vec.Norm2();


    // now map posEtaBoundary in 3D-Xi-coordinates
    LINALG::Matrix<nsd,1> posXiDomain;
    posXiDomain.Clear();
    // coordinates of this integration point in element coordinates \xi^domain

    // the same results like GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain

    //TODO: check this transformation
    // schott Apr 30, 2010
    //GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

    // schott Apr 30, 2010
    static LINALG::Matrix<nsd,1> posEtaBoundaryinXYZ;
    posEtaBoundaryinXYZ.Clear();
    static LINALG::Matrix<numvertices,1> shape_boundary;
    DRT::UTILS::shape_function_2D(shape_boundary,posEtaBoundary(0,0),posEtaBoundary(1,0),cell.Shape());

    posEtaBoundaryinXYZ.Multiply(cellXYZDomain,shape_boundary);
    //get Xi-coordinates of posEtaBoundary 3D -> 3D transformation in element, not in cell
    COMBUST::MapXYZToXiDomain<DISTYPE>(xyze,posEtaBoundaryinXYZ,posXiDomain);

    //        GEO::currentToVolumeElementCoordinates(DISTYPE, xyze, posEtaBoundaryinXYZ, posXiDomain);


    //----------------------------
    // compute integration factors
    //----------------------------
    // compute spatial integration factor
    const double fac = intpoints.qwgt[iquad]*detXYZtoEta1Eta2;
    // compute total (time and spatial) integration factor (and coefficients of single terms)?
    const double timefacfac = timefac * fac;
#else
    //---------------------------------------------
    // transform coordinates of this Gaussian point
    //---------------------------------------------
    // coordinates of this integration point in boundary cell coordinates \eta^boundary
    const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);

    // coordinates of this integration point in element coordinates \xi^domain
    LINALG::Matrix<nsd,1> posXiDomain;
    GEO::mapEtaBToXiD(cell, posEtaBoundary, posXiDomain);

    //------------------------------------------------
    // compute Jacobian matrix for mapping to boundary
    //------------------------------------------------
    static LINALG::Matrix<nsd,numvertices> deriv_boundary;
    deriv_boundary.Clear();
    DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posEtaBoundary(0),posEtaBoundary(1),CELLDISTYPE);

    // get Jacobian matrix d \xi^domain / d \eta^boundary  (3x2)
    static LINALG::Matrix<nsd,2> dxideta;
    dxideta.Clear();
    for (std::size_t k=0; k!=numvertices; ++k)
      for (std::size_t i=0; i!=nsd; ++i)
        for (std::size_t j=0; j!=2; ++j)
          dxideta(i,j) += cellXYZDomain(i,k)*deriv_boundary(j,k);
    //TODO clarify
    //dxideta(i,j) += cellXiDomain(i,k)*deriv_boundary(j,k);

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
      dserror("negative Jacobian determinant detXitoEta: %f for element %d!", detXitoEta, ele->Id());
    }
    //#endif

    //----------------------------
    // compute integration factors
    //----------------------------
    // compute spatial integration factor
    const double fac = intpoints.qwgt[iquad]*detXitoEta;
    //TODO clarify
    //const double fac = intpoints.qwgt[iquad]*detXtoXi*detXitoEta;
    // compute total (time and spatial) integration factor (and coefficients of single terms)?
    const double timefacfac = timefac * fac;
#endif

    //--------------------------------------------------------------------------------------------
    // compute normal vector (normed, in physical coordinates)
    // remark: for linear boundary integrastion cells this could be done before the loop over all
    //         Gaussian points
    //--------------------------------------------------------------------------------------------
    static LINALG::Matrix<nsd,1> normal(true);
    normal.Clear();
    /*
     * author schott
     * Apr 15, 2010
     *
     * cells are triangulars or quadrilaterals (if a cell lies on a face)
     * since the normal is equal at all gausspoints
     */
    // TODO: implement also normals based on nabla(phi)/|nabla(phi)| (interpolated), necessary for surface tension //Benedikt 04/26/2010
    // here normal points from Omega^+ to Omega^-
    GEO::computeNormalToSurfaceElement(CELLDISTYPE, cellXYZDomain, posEtaBoundary, normal);

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
    jflux(0) = jflux_n*normal(0);
    jflux(1) = jflux_n*normal(1);
    jflux(2) = jflux_n*normal(2);
    //cout << "jflux " << jflux_n << endl;
    // overwrite velocity jump condition (e.g. for specific test case)
    //jflux(0) = 0.0;
    //jflux(1) = 0.0;
    //jflux(2) = 0.0;

#ifdef BENEDIKT
    //--------------------------------------------------------------------------------------------
    // evaluate jump conditions at this Gaussian point
    //--------------------------------------------------------------------------------------------
    // Benedikt 04/26/10
    // velocity jump vector, in general there is a jump in normal direction
    // [|v|] = v^+ - v^- = ju
    static LINALG::Matrix<nsd,1> jump_vel;
    jump_vel.Clear();

    if(veljumptype == INPAR::COMBUST::Vel_Jump_Const)
    {
      //          jump_vel(0,0) = normal(0)*ju;
      //          jump_vel(1,0) = normal(1)*ju;
      //          jump_vel(2,0) = normal(2)*ju;
      //
// TODO @Benedikt hier ist das 'ju' wohl anders definiert als bei mir  Flo 31/08/10
cout << "/!\\ hier ist was auskommentiert" << endl;
//      jump_vel(0,0) = ju*1.0;
//      jump_vel(1,0) = 0.0;
//      jump_vel(2,0) = 0.0;
    }

    if(veljumptype == INPAR::COMBUST::Vel_Jump_PremixedCombustion)
    {
      jump_vel(0,0) = normal(0)*ju_n;
      jump_vel(1,0) = normal(1)*ju_n;
      jump_vel(2,0) = normal(2)*ju_n;
    }


    // Benedikt 04/26/10
    // normal tension jump vector, in general there is a jump in normal direction (i.e. surface tension)
    // [|sigma*n|] = [(sigma)^+ - (sigma)^- ]*n= jump_tension

    static LINALG::Matrix<nsd,1> jump_tension;
    jump_tension.Clear();

    if(normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_Const)
    {
      // normal tension jump is set to jt
//TODO @Benedikt compiler problem, weil ich 'jt' rausgeschmissen habe    Flo 31/08/10
cout << "/!\\ hier ist was auskommentiert" << endl;
//      jump_tension(0,0) = jt * normal(0);
//      jump_tension(1,0) = jt * normal(1);
//      jump_tension(2,0) = jt * normal(2);
    }

    if(normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_PremixedCombustion)
    {
      // jump of traction vector (stress tensor * normal)
      //LINALG::Matrix<nsd,1> jtension(true);
      jump_tension(0) = jflux_n*normal(0);
      jump_tension(1) = jflux_n*normal(1);
      jump_tension(2) = jflux_n*normal(2);
    }
    // no jump value necessary for surface tension  => another implementation (jump value implicit given)

#endif
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
    enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
        enrfunct_plus, enrderxy_plus, enrderxy2_plus);
    enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Velx, funct, derxy, derxy2,
        enrfunct_minus, enrderxy_minus, enrderxy2_minus);

    // (enriched) shape functions = approximation functions (P = N * \Psi)
    static XFEM::ApproxFunc<0,shpVecSize> shp_jump;
    static XFEM::ApproxFunc<1,shpVecSize> shp_mean_inv;
    static XFEM::ApproxFunc<1,shpVecSize> shp_mean;
    static XFEM::ApproxFunc<2,shpVecSize> shp_mean_visc;
    static XFEM::ApproxFunc<2,shpVecSize> shp_jump_visc;

#ifdef TPF_NITSCHE_LPS
    static XFEM::ApproxFunc<0,shpVecSize> shp_presstab;
#endif

    // fill approximation functions
    for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
    {
      shp_jump.d0(iparam) = enrfunct_minus(iparam) - enrfunct_plus(iparam);
      // {}-average
      shp_mean.d0(iparam) = kappa_plus * enrfunct_plus(iparam) + kappa_minus * enrfunct_minus(iparam);
      // <>-average
      shp_mean_inv.d0(iparam) = kappa_minus * enrfunct_plus(iparam) + kappa_plus * enrfunct_minus(iparam);

      shp_mean.dx(iparam) = kappa_plus * enrderxy_plus(0,iparam) + kappa_minus * enrderxy_minus(0,iparam);
      shp_mean.dy(iparam) = kappa_plus * enrderxy_plus(1,iparam) + kappa_minus * enrderxy_minus(1,iparam);
      shp_mean.dz(iparam) = kappa_plus * enrderxy_plus(2,iparam) + kappa_minus * enrderxy_minus(2,iparam);

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
#ifdef TPF_NITSCHE_LPS
      shp_presstab.d0(iparam) = shp_jump.d0(iparam) - shp_presproj_jump(iparam);
      //shp_presstab.d0(iparam) = shp_presproj_jump(iparam);
      //cout << shp_jump.d0(iparam) - shp_presproj_jump(iparam) << endl;
#endif
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

#ifdef TPF_NITSCHE_LPS
    static double pstab;
    pstab = 0.0;
    for (size_t iparam = 0; iparam !=numparampres; ++iparam)
      pstab += shp_presstab.d0(iparam)*eprenp(iparam);
#endif

    static double pmean;
    pmean = 0.0;
    for (size_t iparam = 0; iparam != numparampres; ++iparam)
      pmean += shp_mean.d0(iparam)*eprenp(iparam);


    //--------------------------------------------------------------------------------------------
    // build single boundary integral stiffness terms
    //
    // remarks: || x || stands for the jump operator ( || x || = x^+ - x^- )
    //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
    //             n    stands for the normal vector on the interface pointing into the
    //                  burnt domain (n = n^- = -n^+)
    //--------------------------------------------------------------------------------------------

#if 0
    if(normaltensionjumptype == INPAR::COMBUST::NormalTension_Jump_SurfaceTension)
    {
      //================================ SURFACE TENSION ASSEMBLY =================================
      //---------------------assemble 2D surface integrals -------------------------
      COMBUST::RHS_SurfaceTension2DIntegrals<DISTYPE,ASSTYPE,NUMDOF>(posXiDomain,shp_mean_inv,ele,xyze,grad_ephi_smoothed,normal,timefacfac,surftensapprox,surftenscoeff,assembler);
      //===========================================================================================
    }
    else // all other normal tension jumps
    {
      //----------------------------------------     |                     |
      // surface tension (constant jump) term      - |  < v > , J_t >      |
      //----------------------------------------     |                     |

      assembler.template Vector<Velx>(shp_mean_inv.d0, - timefacfac*jump_tension(0,0) );
      assembler.template Vector<Vely>(shp_mean_inv.d0, - timefacfac*jump_tension(1,0) );
      assembler.template Vector<Velz>(shp_mean_inv.d0, - timefacfac*jump_tension(2,0) );
    }
#endif

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
    // - |  || v || , { 2\mu epsilon( u_i ) } n  |
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

#ifdef COMBUST_NITSCHE_VISCCONSISTJUMP
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

#ifdef COMBUST_NITSCHE_VISCADJCONSIST
    //---------------------------------    |                                     |
    // viscous adjoint consistency term  + |  { 2\mu epsilon( v ) } n, || Du ||  |
    //---------------------------------    |                                     |

    assembler.template Matrix<Velx,Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<Velx,Velx>(shp_mean_visc.dy,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<Velx,Velx>(shp_mean_visc.dz,     timefacfac*normal(2), shp_jump.d0);
    assembler.template Matrix<Vely,Velx>(shp_mean_visc.dx,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<Velz,Velx>(shp_mean_visc.dx,     timefacfac*normal(2), shp_jump.d0);

    assembler.template Matrix<Velx,Vely>(shp_mean_visc.dy,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<Vely,Vely>(shp_mean_visc.dx,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<Vely,Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<Vely,Vely>(shp_mean_visc.dz,     timefacfac*normal(2), shp_jump.d0);
    assembler.template Matrix<Velz,Vely>(shp_mean_visc.dy,     timefacfac*normal(2), shp_jump.d0);

    assembler.template Matrix<Velx,Velz>(shp_mean_visc.dz,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<Vely,Velz>(shp_mean_visc.dz,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<Velz,Velz>(shp_mean_visc.dx,     timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<Velz,Velz>(shp_mean_visc.dy,     timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<Velz,Velz>(shp_mean_visc.dz, 2.0*timefacfac*normal(2), shp_jump.d0);

    //   |                                      |
    // - |  { 2\mu epsilon( v ) } n, || u_i ||  |
    //   |                                      |

    //      assembler.template Vector<Velx>(shp_mean_visc.dx, -2.0* normal(0)           *timefacfac*veljump(0,0));
    //      assembler.template Vector<Velx>(shp_mean_visc.dy, -    (normal(0)+normal(1))*timefacfac*veljump(0,0));
    //      assembler.template Vector<Velx>(shp_mean_visc.dz, -    (normal(0)+normal(2))*timefacfac*veljump(0,0));

    assembler.template Vector<Velx>(shp_mean_visc.dx, -2.0*timefacfac*normal(0)*veljump(0,0));
    assembler.template Vector<Velx>(shp_mean_visc.dy, -    timefacfac*normal(1)*veljump(0,0));
    assembler.template Vector<Velx>(shp_mean_visc.dz, -    timefacfac*normal(2)*veljump(0,0));
    assembler.template Vector<Vely>(shp_mean_visc.dx, -    timefacfac*normal(1)*veljump(0,0));
    assembler.template Vector<Velz>(shp_mean_visc.dx, -    timefacfac*normal(2)*veljump(0,0));

    //      assembler.template Vector<Vely>(shp_mean_visc.dx, -    (normal(1)+normal(0))*timefacfac*veljump(1,0));
    //      assembler.template Vector<Vely>(shp_mean_visc.dy, -2.0* normal(1)           *timefacfac*veljump(1,0));
    //      assembler.template Vector<Vely>(shp_mean_visc.dz, -    (normal(1)+normal(2))*timefacfac*veljump(1,0));

    assembler.template Vector<Velx>(shp_mean_visc.dy, -    timefacfac*normal(0)*veljump(1,0));
    assembler.template Vector<Vely>(shp_mean_visc.dx, -    timefacfac*normal(0)*veljump(1,0));
    assembler.template Vector<Vely>(shp_mean_visc.dy, -2.0*timefacfac*normal(1)*veljump(1,0));
    assembler.template Vector<Vely>(shp_mean_visc.dz, -    timefacfac*normal(2)*veljump(1,0));
    assembler.template Vector<Velz>(shp_mean_visc.dy, -    timefacfac*normal(2)*veljump(1,0));

    //      assembler.template Vector<Velz>(shp_mean_visc.dx, -    (normal(2)+normal(0))*timefacfac*veljump(2,0));
    //      assembler.template Vector<Velz>(shp_mean_visc.dy, -    (normal(2)+normal(1))*timefacfac*veljump(2,0));
    //      assembler.template Vector<Velz>(shp_mean_visc.dz, -2.0* normal(2)           *timefacfac*veljump(2,0));

    assembler.template Vector<Velx>(shp_mean_visc.dz, -    timefacfac*normal(0)*veljump(2,0));
    assembler.template Vector<Vely>(shp_mean_visc.dz, -    timefacfac*normal(1)*veljump(2,0));
    assembler.template Vector<Velz>(shp_mean_visc.dx, -    timefacfac*normal(0)*veljump(2,0));
    assembler.template Vector<Velz>(shp_mean_visc.dy, -    timefacfac*normal(1)*veljump(2,0));
    assembler.template Vector<Velz>(shp_mean_visc.dz, -2.0*timefacfac*normal(2)*veljump(2,0));

    //-------------------------------------    |                                    |
    // viscous adjoint consistency term RHS  + |  { 2\mu epsilon( v ) } n, j_u  |
    //-------------------------------------    |                                    |
    //      assembler.template Vector<Velx>(shp_mean_visc.dx, 2.0* normal(0)           *timefacfac*normal(0)*ju);
    //      assembler.template Vector<Velx>(shp_mean_visc.dy,     (normal(0)+normal(1))*timefacfac*normal(0)*ju);
    //      assembler.template Vector<Velx>(shp_mean_visc.dz,     (normal(0)+normal(2))*timefacfac*normal(0)*ju);

    assembler.template Vector<Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0)*ju(0));
    assembler.template Vector<Velx>(shp_mean_visc.dy,     timefacfac*normal(1)*ju(0));
    assembler.template Vector<Velx>(shp_mean_visc.dz,     timefacfac*normal(2)*ju(0));
    assembler.template Vector<Vely>(shp_mean_visc.dx,     timefacfac*normal(1)*ju(0));
    assembler.template Vector<Velz>(shp_mean_visc.dx,     timefacfac*normal(2)*ju(0));
    //assembler.template Vector<Velx>(shp_mean_visc.dx, 2.0*timefacfac*normal(0)*normal(0)*ju_n);
    //assembler.template Vector<Velx>(shp_mean_visc.dy,     timefacfac*normal(1)*normal(0)*ju_n);
    //assembler.template Vector<Velx>(shp_mean_visc.dz,     timefacfac*normal(2)*normal(0)*ju_n);
    //assembler.template Vector<Vely>(shp_mean_visc.dx,     timefacfac*normal(1)*normal(0)*ju_n);
    //assembler.template Vector<Velz>(shp_mean_visc.dx,     timefacfac*normal(2)*normal(0)*ju_n);

    //      assembler.template Vector<Vely>(shp_mean_visc.dx,     (normal(1)+normal(0))*timefacfac*normal(1)*ju);
    //      assembler.template Vector<Vely>(shp_mean_visc.dy, 2.0* normal(1)           *timefacfac*normal(1)*ju);
    //      assembler.template Vector<Vely>(shp_mean_visc.dz,     (normal(1)+normal(2))*timefacfac*normal(1)*ju);

    assembler.template Vector<Velx>(shp_mean_visc.dy,     timefacfac*normal(0)*ju(1));
    assembler.template Vector<Vely>(shp_mean_visc.dx,     timefacfac*normal(0)*ju(1));
    assembler.template Vector<Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1)*ju(1));
    assembler.template Vector<Vely>(shp_mean_visc.dz,     timefacfac*normal(2)*ju(1));
    assembler.template Vector<Velz>(shp_mean_visc.dy,     timefacfac*normal(2)*ju(1));
    //assembler.template Vector<Velx>(shp_mean_visc.dy,     timefacfac*normal(0)*normal(1)*ju_n);
    //assembler.template Vector<Vely>(shp_mean_visc.dx,     timefacfac*normal(0)*normal(1)*ju_n);
    //assembler.template Vector<Vely>(shp_mean_visc.dy, 2.0*timefacfac*normal(1)*normal(1)*ju_n);
    //assembler.template Vector<Vely>(shp_mean_visc.dz,     timefacfac*normal(2)*normal(1)*ju_n);
    //assembler.template Vector<Velz>(shp_mean_visc.dy,     timefacfac*normal(2)*normal(1)*ju_n);

    //      assembler.template Vector<Velz>(shp_mean_visc.dx,     (normal(2)+normal(0))*timefacfac*normal(2)*ju);
    //      assembler.template Vector<Velz>(shp_mean_visc.dy,     (normal(2)+normal(1))*timefacfac*normal(2)*ju);
    //      assembler.template Vector<Velz>(shp_mean_visc.dz, 2.0* normal(2)           *timefacfac*normal(2)*ju);

    assembler.template Vector<Velx>(shp_mean_visc.dz,     timefacfac*normal(0)*ju(2));
    assembler.template Vector<Vely>(shp_mean_visc.dz,     timefacfac*normal(1)*ju(2));
    assembler.template Vector<Velz>(shp_mean_visc.dx,     timefacfac*normal(0)*ju(2));
    assembler.template Vector<Velz>(shp_mean_visc.dy,     timefacfac*normal(1)*ju(2));
    assembler.template Vector<Velz>(shp_mean_visc.dz, 2.0*timefacfac*normal(2)*ju(2));
    //assembler.template Vector<Velx>(shp_mean_visc.dz,     timefacfac*normal(0)*normal(2)*ju_n);
    //assembler.template Vector<Vely>(shp_mean_visc.dz,     timefacfac*normal(1)*normal(2)*ju_n);
    //assembler.template Vector<Velz>(shp_mean_visc.dx,     timefacfac*normal(0)*normal(2)*ju_n);
    //assembler.template Vector<Velz>(shp_mean_visc.dy,     timefacfac*normal(1)*normal(2)*ju_n);
    //assembler.template Vector<Velz>(shp_mean_visc.dz, 2.0*timefacfac*normal(2)*normal(2)*ju_n);
#endif

#ifdef COMBUST_NITSCHE_PRESCONSISTMEAN
    //--------------------------    |                     |
    // pressure consistency term  - |  || v ||, { Dp } n  |
    //--------------------------    |                     |

    assembler.template Matrix<Velx,Pres>(shp_jump.d0, -timefacfac*normal(0), shp_mean.d0);
    assembler.template Matrix<Vely,Pres>(shp_jump.d0, -timefacfac*normal(1), shp_mean.d0);
    assembler.template Matrix<Velz,Pres>(shp_jump.d0, -timefacfac*normal(2), shp_mean.d0);
    //      assembler.template Matrix<Velx,Pres>(shp_jump.d0, -timefacfac*normal(0), shp_presproj_mean);
    //      assembler.template Matrix<Vely,Pres>(shp_jump.d0, -timefacfac*normal(1), shp_presproj_mean);
    //      assembler.template Matrix<Velz,Pres>(shp_jump.d0, -timefacfac*normal(2), shp_presproj_mean);

    //   |                      |
    // + |  || v ||, { p_i } n  |
    //   |                      |

    assembler.template Vector<Velx>(shp_jump.d0, timefacfac*normal(0)*pmean);
    assembler.template Vector<Vely>(shp_jump.d0, timefacfac*normal(1)*pmean);
    assembler.template Vector<Velz>(shp_jump.d0, timefacfac*normal(2)*pmean);
#endif

#ifdef COMBUST_NITSCHE_PRESCONSISTJUMP
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

#ifdef COMBUST_NITSCHE_PRESADJCONSIST
    //----------------------------------    |                     |
    // pressure adjoint consistency term  + |  { q }, || Du || n  |
    //----------------------------------    |                     |

    assembler.template Matrix<Pres,Velx>(shp_mean.d0, timefacfac*normal(0), shp_jump.d0);
    assembler.template Matrix<Pres,Vely>(shp_mean.d0, timefacfac*normal(1), shp_jump.d0);
    assembler.template Matrix<Pres,Velz>(shp_mean.d0, timefacfac*normal(2), shp_jump.d0);
    //      assembler.template Matrix<Pres,Velx>(shp_presproj_mean, +timefacfac*normal(0), shp_jump.d0);
    //      assembler.template Matrix<Pres,Vely>(shp_presproj_mean, +timefacfac*normal(1), shp_jump.d0);
    //      assembler.template Matrix<Pres,Velz>(shp_presproj_mean, +timefacfac*normal(2), shp_jump.d0);

    //   |                      |
    // - |  { q }, || u_i || n  |
    //   |                      |

    assembler.template Vector<Pres>(shp_mean.d0, -timefacfac*(veljump(0,0)*normal(0)
        +veljump(1,0)*normal(1)
        +veljump(2,0)*normal(2)));
    //      assembler.template Vector<Pres>(shp_presproj_mean, -timefacfac*(vjump(0,0)*normal(0)
    //          +vjump(1,0)*normal(1)
    //          +vjump(2,0)*normal(2)));

    //--------------------------------------    |                |
    // pressure adjoint consistency term RHS  + |  { q }, j_u n  |
    //--------------------------------------    |                |

    assembler.template Vector<Pres>(shp_mean.d0, timefacfac*(ju(0)*normal(0)
        +ju(1)*normal(1)
        +ju(2)*normal(2)));
    //assembler.template Vector<Pres>(shp_mean.d0, timefacfac*ju_n);
    
    //      assembler.template Vector<Pres>(shp_presproj_mean, +timefacfac*(jump_vel(0,0)*normal(0)
    //          +jump_vel(1,0)*normal(1)
    //          +jump_vel(2,0)*normal(2)));
#endif

#ifdef COMBUST_NITSCHE_FLUXJUMP
    //---------------------------------------------------    |                 |
    // flux (traction/stress) jump (consistency) term RHS  - |  < v >, j_flux  |
    //---------------------------------------------------    |                 |

    assembler.template Vector<Velx>(shp_mean.d0, -timefacfac*jflux(0));
    assembler.template Vector<Vely>(shp_mean.d0, -timefacfac*jflux(1));
    assembler.template Vector<Velz>(shp_mean.d0, -timefacfac*jflux(2));
#if 0
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
#endif
#ifdef COMBUST_NITSCHE_NITSCHEVEL
    //------------------------    |                              |
    // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
    //------------------------    |                              |

    assembler.template Matrix<Velx,Velx>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
    assembler.template Matrix<Vely,Vely>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
    assembler.template Matrix<Velz,Velz>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);

    //    |                               |
    //  - |  \alpha_u || v ||, || u_i ||  |
    //    |                               |

    assembler.template Vector<Velx>(shp_jump.d0, -alphau*timefacfac*veljump(0,0));
    assembler.template Vector<Vely>(shp_jump.d0, -alphau*timefacfac*veljump(1,0));
    assembler.template Vector<Velz>(shp_jump.d0, -alphau*timefacfac*veljump(2,0));

    //----------------------------    |                         |
    // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, j_u  |
    //----------------------------    |                         |

    //      assembler.template Vector<Velx>(shp_jump.d0, alphau*timefacfac*ju(0));
    //      assembler.template Vector<Vely>(shp_jump.d0, alphau*timefacfac*ju(1));
    //      assembler.template Vector<Velz>(shp_jump.d0, alphau*timefacfac*ju(2));
    assembler.template Vector<Velx>(shp_jump.d0, alphau*timefacfac*normal(0)*ju_n);
    assembler.template Vector<Vely>(shp_jump.d0, alphau*timefacfac*normal(1)*ju_n);
    assembler.template Vector<Velz>(shp_jump.d0, alphau*timefacfac*normal(2)*ju_n);
#endif

#ifdef COMBUST_NITSCHE_NITSCHEPRES
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

#ifdef TPF_NITSCHE_LPS
    //------------------------------------------    |                                                           |
    // LPS-stabilization at interface (pressure)  + |  \alpha_p (|| q ||-||Proj(q)||), (|| Dp ||-||Proj(Dp)||)  |
    //------------------------------------------    |                                                           |

    assembler.template Matrix<Pres,Pres>(shp_presstab.d0, alphap*timefacfac  , shp_presstab.d0);

    //    |                                                             |
    //  - |  \alpha_p (|| q ||-||Proj(q)||), (|| p_i ||-||Proj(p_i)||)  |
    //    |                                                             |

    //cout << pstab << endl;

    assembler.template Vector<Pres>(shp_presstab.d0, - alphap*timefacfac*pstab);

    // TODO: implement RHS for old time step
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
void SysmatBoundaryNitsche(
    const DRT::ELEMENTS::Combust3*             ele,            ///< the element those matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>&  ih,  ///< connection to the interface handler
    const XFEM::ElementDofManager&             dofman,         ///< dofmanager of the current element
    const M1&                                  evelnp,         ///<
    const V1&                                  eprenp,         ///<
    const V2&                                  ephi,           ///<
    const M2&                                  egradphi,
    Teuchos::RCP<const MAT::Material>          material,       ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme  timealgo,       ///< time discretization type
    const double&                              dt,             ///< delta t (time step size)
    const double&                              theta,          ///< factor for one step theta scheme
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&  assembler,      ///<
    const double                               flamespeed,     ///<
    const double                               nitschevel,     ///<
    const double                               nitschepres,    ///<
    const double                      ele_meas_plus,    /// measure of element in plus domain
    const double                      ele_meas_minus,    /// measure of element in minus domain
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const double                               surftenscoeff,   ///<
    const bool                                 connected_interface,
    const INPAR::COMBUST::VelocityJumpType veljumptype,
    const INPAR::COMBUST::NormalTensionJumpType normaltensionjumptype
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
  const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

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

  //-----------------------------
  // prepare interface conditions
  //-----------------------------
  // mass flux through flame front
  const double mflux = flamespeed*dens_minus;

  // velocity jump in normal direction
#if 0
  // jump velocity vector ju_{n,t1,t2}
  LINALG::Matrix<nsd,1> ju_nt(true);
  ju_nt(0,0) = -flamespeed*dens_minus*(1.0/dens_minus - 1.0/dens_plus);
  ju_nt(1,0) = 0.0;
  ju_nt(2,0) = 0.0;
#else
  const double ju_n = -mflux*(1.0/dens_minus - 1.0/dens_plus);
#endif
  // jump in normal stress
  const double jflux_n = -mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);

  //----------------------------------------
  // compute Nitsche stabilization patameter
  //----------------------------------------
  // calculate element diameter
  const double hk_eleDiam = COMBUST::getEleDiameter<DISTYPE>(xyze);

  // get maximum of the viscosities of all fluid materials
  const double dynvisc_max = max(visc_plus,visc_minus);

  // Nitsche parameter velocity ( alph_u := 1/h_k * mu_max * gamma_u ) see Becker, Hansbo 2009
  const double alphau = 1.0/hk_eleDiam * dynvisc_max * nitschevel;

#ifdef COMBUST_NITSCHE_NITSCHEPRES
  // Nitsche parameter pressure
  const double alphap = nitschepres;
#endif
#ifdef TPF_NITSCHE_LPS
  const double alphap = nitschepres;
#endif
  //-----------------------------------------------------------------------------------------
  // get weights for average operators { () } = kappa_plus * ()_plus + kappa_minus * ()_minus
  //                               and < () > = kappa_minus* ()_plus + kappa_plus  * ()_minus
  //-----------------------------------------------------------------------------------------

  double ele_meas = ele_meas_plus + ele_meas_minus;
  double kappa_plus  = ele_meas_plus/ele_meas;
  double kappa_minus = ele_meas_minus/ele_meas;

  // flag for higher order elements
  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

  //-------------------------------------------------------------------------------
  // evaluate the enrichment function at the interface (boundary integration cells)
  //-------------------------------------------------------------------------------
  const XFEM::ElementEnrichmentValues enrvals_plus(*ele,dofman,XFEM::Enrichment::approachFromPlus,ephi);
  const XFEM::ElementEnrichmentValues enrvals_minus(*ele,dofman,XFEM::Enrichment::approachFromMinus,ephi);

#ifdef TPF_NITSCHE_LPS
  /*
   * author schott
   * May 6, 2010
   *
   * implements a LPS (local projection stabilization) for pressure oszillations at the interface
   * for trilinear shape functions
   */
  // calculate
  double meas_ele_plus = 0.0;   // measure of element ele at plus side
  double meas_ele_minus = 0.0;  // measure of element ele at minus side

  //--------------------------------------------------------------------------------------------
  // shape functions for pressure projection as approximation functions
  // (constant values for each element! -> before cell loop)
  //--------------------------------------------------------------------------------------------
  const std::size_t shpPresProjVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // temporary arrays holding enriched shape functions (int^+-(N * \Psi)) on either side of the interface
  static LINALG::Matrix<shpPresProjVecSize,1>       shp_presproj_plus;
  shp_presproj_plus.Clear();
  static LINALG::Matrix<shpPresProjVecSize,1>       shp_presproj_minus;
  shp_presproj_minus.Clear();


  //  TWOPHASEFLOWJUMP::BuildPresProjShapefunction<DISTYPE,ASSTYPE,shpPresProjVecSize,NUMDOF>(meas_ele_plus, meas_ele_minus,
  //        shp_presproj_plus, shp_presproj_minus,
  //          ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, eepsilon,
  //          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);

  COMBUST::BuildPresProjShapefunction<DISTYPE,ASSTYPE,shpPresProjVecSize,NUMDOF>(meas_ele_plus, meas_ele_minus,
      shp_presproj_plus, shp_presproj_minus,
      ele, ih, dofman,ephi,
      material);
  //  cout << "meas_ele_plus:\t" << meas_ele_plus << endl;
  //  cout << "meas_ele_minus:\t" << meas_ele_minus << endl << endl;
  //
  //  for  (int iparam = 0; iparam< (int)shpPresProjVecSize; iparam++){
  //    cout << iparam << "\t" << shp_presproj_plus(iparam) << endl;
  //  }
  //
  //  cout << endl;
  //
  //  for  (int iparam = 0; iparam< (int)shpPresProjVecSize; iparam++){
  //    cout << iparam << "\t" << shp_presproj_minus(iparam) << endl;
  //  }
  //
  //  cout << endl;

  static LINALG::Matrix<shpPresProjVecSize,1>   shp_presproj_jump;
  shp_presproj_jump.Clear();

  //  static LINALG::Matrix<shpPresProjVecSize,1>  shp_presproj_mean;
  //  shp_presproj_mean.Clear();

  for (int iparam = 0; iparam < (int)shpPresProjVecSize; iparam++)
  {
    shp_presproj_jump(iparam) = (1.0/meas_ele_plus * shp_presproj_plus(iparam)
        - 1.0/meas_ele_minus * shp_presproj_minus(iparam));
    //shp_presproj_mean(iparam) = 0.5/meas_ele_plus * shp_presproj_plus(iparam)+ 0.5/meas_ele_minus * shp_presproj_minus(iparam);
  }

  // for each element get shp_presstab = shp_jump.d0 - shp_presproj_jump

#endif


#if(0)
  //==================================================================================
  // schott Jun 8, 2010
  // calculate smoothed phi and grad(phi) via Least Squares Method, see Merchandise 2007
  LINALG::Matrix<numnode,1> ephi_smoothed;
  ephi_smoothed.Clear();

  LINALG::Matrix<nsd,numnode> grad_ephi_smoothed;
  grad_ephi_smoothed.Clear();

  smooth_ephi<DISTYPE>(discretization,ih,ele,ephi,ephi_smoothed,grad_ephi_smoothed);

  //  cout << "1. ephi\n"   << ephi << endl;
  //  cout << "2. ephi_smoothed\n" << ephi_smoothed << endl;
  //  cout << "3. grad_ephi_smoothed\n" << grad_ephi_smoothed << endl;
  //  cout << " geht noch eins weiter" << endl;
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
          evelnp, eprenp, ephi, egradphi,
          ju_n, jflux_n,
          timefac, assembler,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, visc_plus, visc_minus,
          alphau,
          numparamvelx, numparampres, surftensapprox, surftenscoeff, connected_interface, veljumptype, normaltensionjumptype
      );
      break;
    case DRT::Element::quad4:
      // TODO check this intrule !!!
      COMBUST::Nitsche_SysmatBoundaryCell<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
          *cell, DRT::UTILS::intrule_quad_9point , xyze, ele, higher_order_ele,
          evelnp, eprenp, ephi, egradphi,
          ju_n, jflux_n,
          timefac, assembler,
          enrvals_plus, enrvals_minus, kappa_plus, kappa_minus, visc_plus, visc_minus,
          alphau,
          numparamvelx, numparampres, surftensapprox, surftenscoeff, connected_interface, veljumptype, normaltensionjumptype
      );
      break;
    default:
      dserror("invalid type of boundary integration cell");
    }
  } // loop boundary integration cells

  return;
}



/*!
  Calculate smoothed level set field for element via least squares method
 */
template <DRT::Element::DiscretizationType DISTYPE,
          class M1, class V1>
void Nitsche_smooth_ephi(
    DRT::Discretization&              discretization,
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>&  ih, /// interface handle to get phinp
    const DRT::ELEMENTS::Combust3*    ele,                    /// < the element those matrix is calculated
    const V1&                         ephi,                   /// phi-values for each node
    V1&                               ephi_smoothed,          /// smoothed phi
    M1&                               grad_ephi_smoothed      /// smoothed grad of phi
)
{
  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element for interpolation!!!
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // number of nodes of this element (number of vertices)
  // e.g. hex20 elements has numnode=20 but numberOfNodes=8
  const int numberOfNodes = ele->NumNode(); //
  //  cout << "numberOfNodes\t" << numberOfNodes << endl;


  // get pointer to nodes of current element ele
  const int* ptToNodeIds = ele->NodeIds();

  //===========================loop over nodes of element===========================
  for (int nodeId=0; nodeId< numberOfNodes; nodeId++)
  {
    // current node
    DRT::Node* actnode = discretization.lRowNode(ptToNodeIds[nodeId]);

    //get adjacent elements to current node
    int numberOfElements = actnode->NumElement();

    //    cout << "numberOfElements\t" << numberOfElements << endl;
    // check if at least nsd+1 adjacent elements are available
    // if nsd==3 we need more than 1+nsd adjacent elements to solve LEAST SQUARES system
    //    if(numberOfElements < (int)(nsd+1))
    //TODO: !!!
    if(numberOfElements != 8)
    {
      cout << numberOfElements;
      dserror("we need at least nsd+1 adjacent elements!!!");
    }

    // get adjacent elements to current node actnode
    DRT::Element** elements = actnode->Elements();


    // we need Epetra-Matrix (numberOfElements is not a valid template parameter)
    Epetra_SerialDenseMatrix RHS_LS(numberOfElements,1);
    Epetra_SerialDenseMatrix MAT_LS(numberOfElements, nsd+1);
    //=========================loop over all adjacent elements to current node
    for (int ele_current=0; ele_current<numberOfElements; ele_current++)
    {
      // get current element
      const DRT::Element* ele_adj = elements[ele_current];

      const int* ptToNodeIds_adj = ele_adj->NodeIds();
      //const int numberOfNodes_adj = ele_adj->NumNode();


      // get phi-values of current adjacent element ele_adj
      // create vector "ephinp" holding scalar phi values for this element
      Epetra_SerialDenseVector ephinp(numnode); //local vector phi-values of adjacent element

      // which node in param space of element ele_adj has actnode
      int ID_param_space = -1;

      // get vector of node GIDs of this adjacent element -> needed for ExtractMyValues
      vector<int> nodeID_adj(numnode);
      for (size_t inode=0; inode < numnode; inode++){
        nodeID_adj[inode] = ptToNodeIds_adj[inode];
        // get local number of node actnode in ele_adj
        if(actnode->Id() == ptToNodeIds_adj[inode]) ID_param_space = inode;
      }

      //          cout << "nodeID_adj\n";
      //          for (int i=0; i< (int)numnode; i++) cout << nodeID_adj[i] << endl;
      //          cout << "ID_param_space\t" << ID_param_space << endl;

      // extract the phi-values of adjacent element with local ids from global vector *phinp
      // get pointer to vector holding G-function values at the fluid nodes
      const Teuchos::RCP<Epetra_Vector> phinp = ih->FlameFront()->Phinp();
      DRT::UTILS::ExtractMyValues(*phinp, ephinp, nodeID_adj);
      LINALG::Matrix<numnode,1> ephi_adj(ephinp);
      //          // copy ephinp into a LINALG::Matrix
      //          for(size_t iparam=0; iparam<numnode; iparam++){
      //              ephi_adj(iparam,0) = ephinp(iparam);
      //          }

      //          cout << "ephinp" << ephinp << endl;
      //          cout << "ephi_adj" << ephi_adj << endl;

      // calculate center of gravity
      static LINALG::Matrix<numnode,1> funct;
      funct.Clear();
      static LINALG::Matrix<nsd,1> centerOfGravXi;
      centerOfGravXi.Clear();

      // xi-coord = 0.0.0 for hex8
      if (DISTYPE != DRT::Element::hex8) dserror("center of gravity implemented only for hex8 elements");
      centerOfGravXi(0) = 0.0;
      centerOfGravXi(1) = 0.0;
      centerOfGravXi(2) = 0.0;

      DRT::UTILS::shape_function_3D(funct,centerOfGravXi(0),centerOfGravXi(1),centerOfGravXi(2),DISTYPE);

      // get node coordinates of this element
      static LINALG::Matrix<nsd,numnode> xyze_adj;
      GEO::fillInitialPositionArray<DISTYPE>(ele_adj, xyze_adj);

      // interpolate center of gravity
      static LINALG::Matrix<nsd,1> centerOfGravXYZ;
      centerOfGravXYZ.Clear();

      centerOfGravXYZ.Multiply(xyze_adj,funct);
      //          cout << "centerOfGravXYZ" << centerOfGravXYZ << endl;

      // calculate vector from current node to point of gravity (direction for taylor)
      static LINALG::Matrix<nsd,1> direction;
      direction.Clear();

      direction(0) = centerOfGravXYZ(0) - xyze_adj(0, ID_param_space);
      direction(1) = centerOfGravXYZ(1) - xyze_adj(1, ID_param_space);
      direction(2) = centerOfGravXYZ(2) - xyze_adj(2, ID_param_space);
      //          cout << "direction\n" << direction << endl;

      // assemble direction into matrix A

      // calculate ephi at point of gravity via interpolation
      static LINALG::Matrix<1,1> phi_adj;
      phi_adj.Clear();
      phi_adj.MultiplyTN(ephi_adj,funct);
      //          cout << "phi_adj\n" << phi_adj << endl;

      // TODO: assemble interpolated value in RHS of LGS which is to solve for each element
      // set RHS and MAT for least squares method


      RHS_LS(ele_current,0) = phi_adj(0,0);
      //set row for phi values
      MAT_LS(ele_current,0) = 1.0; // always set to 1.0
      //set rows for grad_phi values
      MAT_LS(ele_current,1) = direction(0);
      MAT_LS(ele_current,2) = direction(1);
      MAT_LS(ele_current,3) = direction(2);

      //          cout << "RHS_LS" << RHS_LS << endl;
      //          cout << "MAT_LS" << MAT_LS << endl;
      // TODO: special case for boundary elements
    } // end loop over all adjacent elements to current node

    // the system MAT_LS * phi_smoothed = RHS_LS is only solvable in a least squares manner
    // MAT_LS is not square
    // -> solve MAT_LS^T * MAT_LS * phi_smoothed = MAT_LS^T * RHS_LS
    Epetra_SerialDenseMatrix MAT_tmp(nsd+1,nsd+1);
    Epetra_SerialDenseMatrix RHS_tmp(nsd+1,1);

    MAT_tmp.Multiply('T','N', 1.0, MAT_LS,MAT_LS, 0.0);
    RHS_tmp.Multiply('T','N', 1.0, MAT_LS,RHS_LS, 0.0);

    //    cout << "MAT_tmp" << MAT_tmp << endl;
    //    cout << "RHS_tmp" << RHS_tmp << endl;
    // set LINALG-Matrix with fixed size
    // initialize element matrix for A^T*A and RHS-vector to solve Normalengleichung

    //TODO: check constructor!!!
    LINALG::Matrix<nsd+1,nsd+1> MAT(MAT_tmp);
    LINALG::Matrix<nsd+1,1> RHS(RHS_tmp);

    //    cout << "MAT" << MAT << endl;
    //    cout << "RHS" << RHS << endl;
    static LINALG::Matrix<nsd+1,1> PHI_SMOOTHED;
    PHI_SMOOTHED.Clear();


    //==============================solve A^T*A* phi_smoothed = A^T*RHS (LEAST SQUARES)

    // solve the system for current node
    LINALG::FixedSizeSerialDenseSolver<nsd+1,nsd+1,1> solver; //1 is dimension of RHS
    solver.SetMatrix(MAT);
    solver.SetVectors(PHI_SMOOTHED,RHS);
    solver.Solve();

    //    cout << "PHI_SMOOTHED" << PHI_SMOOTHED << endl;
    //    cout << "RHS" << RHS << endl;
    // TODO is nodeID the right index in parameter space of element ele
    // set smoothed node values
    ephi_smoothed(nodeId) = PHI_SMOOTHED(0,0);
    grad_ephi_smoothed(0,nodeId) = PHI_SMOOTHED(1,0);
    grad_ephi_smoothed(1,nodeId) = PHI_SMOOTHED(2,0);
    grad_ephi_smoothed(2,nodeId) = PHI_SMOOTHED(3,0);

  } // end loop over nodes of element

  // TODO: algo von GOLUB

  return;
}


/*
 * author schott
 * May 6, 2010
 *
 *
 */
/*!
  Calculate measure of element in plus and minus domain
  Calculate shape functions for pressure projections -> LPS (local projection stabilization)
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t SHPVECSIZE,
          int NUMDOF,
          class V >
void Nitsche_BuildPresProjShapefunction(
    double &                meas_ele_plus,  // reference to measure of element in plus domain
    double &                meas_ele_minus, // reference to measure of element in minus domain
    LINALG::Matrix<SHPVECSIZE,1> & shp_presproj_plus,   // reference to constant values of pressure projected shapes (plus domain)
    LINALG::Matrix<SHPVECSIZE,1> & shp_presproj_minus,  // reference to constant values of pressure projected shapes (minus domain)
    const DRT::ELEMENTS::Combust3* ele,          ///< element whose matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih, ///< information about the interface
    const XFEM::ElementDofManager& dofman,       ///< dofmanager of this element
    const V&                      ephi,         ///< nodal G-function values for this element
    Teuchos::RCP<const MAT::Material> material  ///< fluid material
)
{
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - pressure projected shape functions ");

  // =========================================================================
  // we need integrals over domain (plus and minus domain)
  // to get an average value of pressure on each side of the interface
  // =========================================================================

  // number space dimensions for 3d combustion element
  const size_t nsd = 3;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // time integration factor
  //  const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

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
  //  const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

  //  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
#if 0
  // figure out whether we have stress unknowns at all
  const bool epsilonele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
  //    const bool velocity_unknowns_present = (getNumParam<ASSTYPE>(dofman, Velx, 1) > 0);
  //    const bool pressure_unknowns_present = (getNumParam<ASSTYPE>(dofman, Pres, 1) > 0);
#endif

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  //  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  //  dsassert(numparamvelx == numparampres, "assumption violation");
#if 0
  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
#endif

  // stabilization parameter
  //  const double hk = FLD::UTILS::HK<DISTYPE>(evelnp,xyze);
  //  const double mk = FLD::UTILS::MK<DISTYPE>();

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
    //    const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
    // get the kinematic viscosity \nu
    //    const double kinvisc = mat->Viscosity();
    // get the density \rho^{n+1}
    //    const double densnp = mat->Density();
    // density \rho^{n}
    // remark: until the XFEM time integration is done properly, densities from previous and
    // current time steps are identical
    //    const double densn = densnp;
    // compute dynamic viscosity \mu
    //    const double dynvisc = kinvisc * densnp;

    //------------------------------------------------------------------------------------------
    // evaluate the enrichment function for this integration cell
    //------------------------------------------------------------------------------------------
    const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);

    //------------------------------------------------------------------------------------------
    // get Gaussian points for this integration cell
    //------------------------------------------------------------------------------------------
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Intersected(), cell->Shape());
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
    //------------------------------------------------------------------------------------------
    // integration loop over Gaussian points
    //------------------------------------------------------------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {


      const LINALG::Matrix<3,1> pos_eta_domain(intpoints.qxg[iquad]);
      /*
       * author schott
       * May 7, 2010
       *
       *
       */
      double det_cell_EtaToXYZ = 0.0;
      static LINALG::Matrix<numnode,1> funct3Dele;
      static LINALG::Matrix<3,numnode> deriv3Dele_xyz;
      static LINALG::Matrix<6,numnode> deriv3Dele_xyz_2nd;
      //const size_t NUMNODE = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      static LINALG::Matrix<3,3> xji_ele_XiToXYZ;


      switch (cell->Shape())
      {
      case DRT::Element::tet4:
        COMBUST::GetDomainIntegrationFeatures<DISTYPE,DRT::Element::tet4,ASSTYPE>(
            ele, pos_eta_domain,xyze, *cell, det_cell_EtaToXYZ, funct3Dele, deriv3Dele_xyz, deriv3Dele_xyz_2nd,xji_ele_XiToXYZ
        );
        //cout << "call tet 4";
        break;
      case DRT::Element::hex8:
        COMBUST::GetDomainIntegrationFeatures<DISTYPE,DRT::Element::hex8,ASSTYPE>(
            ele, pos_eta_domain,xyze, *cell, det_cell_EtaToXYZ, funct3Dele, deriv3Dele_xyz, deriv3Dele_xyz_2nd,xji_ele_XiToXYZ
        );
        //cout << "call hex 8";
        break;
      default:
        dserror("cell distype not implemented yet");
      }
      //cout << xyze << endl;
      double fac = 0.0;
      fac = intpoints.qwgt[iquad] * det_cell_EtaToXYZ;


      //=================================================================================
      //      8. Get enriched shape functions and derivatives as approximation functions
      //=================================================================================
      //--------------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------------
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      //const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;
#if 0
      // comment missing
      static LINALG::Matrix<shpVecSizeStress,1> shp_epsilon;
#endif
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
            funct3Dele,
            deriv3Dele_xyz,
            deriv3Dele_xyz_2nd,
            enr_funct,
            enr_derxy,
            enr_derxy2);

        // fill approximation functions for XFEM
        if (matid == 3) // burnt material
        {
          for (size_t iparam = 0; iparam != numparampres; ++iparam)
          {
            shp_presproj_plus(iparam) += enr_funct(iparam)*fac;
          }
          meas_ele_plus += fac;           // meas is int(1)= sum_gp(1*fac)
        }
        else if(matid == 4) //unburnt material
        {
          for (size_t iparam = 0; iparam != numparampres; ++iparam)
          {
            shp_presproj_minus(iparam) += enr_funct(iparam)*fac;
          }
          meas_ele_minus += fac;
        }
        else dserror("no valid matid in BuildPresProjShapefunction");

#if 0
        if (epsilonele_unknowns_present)
        {
          LINALG::Matrix<shpVecSizeStress,1> enr_funct_stress;

          // shape functions for element dofs
          enrvals.ComputeEnrichedElementShapefunction(
              Sigmaxx,
              funct_stress,
              enr_funct_stress);

          for (size_t iparam = 0; iparam < numparamepsilonxx; ++iparam)
          {
            shp_epsilon(iparam) = enr_funct_stress(iparam);
          }
        }
        else
        {
          shp_epsilon.Clear();
        }
#endif

      }
      else // not xfem_assembly i.e. standard assembly
      {
        // fill approximation functions for standard FEM
        // remark: numparamvelx == numnode, for standard FEM


#ifdef TPF_NITSCHE_LPS_DOMAIN
        if (matid == 3) // burnt material
        {
          for (size_t iparam = 0; iparam != numnode; ++iparam)
          {
            shp_presproj_plus(iparam) += funct3Dele(iparam)*fac;
          }
          meas_ele_plus += fac;           // meas is int(1)= sum_gp(1*fac)
        }
        else if(matid == 4) //unburnt material
        {
          for (size_t iparam = 0; iparam != numnode; ++iparam)
          {
            shp_presproj_minus(iparam) += funct3Dele(iparam)*fac;
          }
          meas_ele_minus += fac;
        }

#else
        dserror("should not be called");
#endif
      }
    } // loop Gaussian points
  }// loop cells
} // end BuildPresProjShapefunction


} // namespace COMBUST


#endif
#endif
