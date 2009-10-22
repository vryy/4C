/*----------------------------------------------------------------------*/
/*!
\file xfluid3_sysmat4.cpp

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

#include <Teuchos_TimeMonitor.hpp>

#include "xdiff3_sysmat.H"
#include "xdiff3_utils.H"
#include "xdiff3_local_assembler.H"
#include "xdiff3_interpolation.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_xfem/spacetime_boundary.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"

  using namespace XFEM::PHYSICS;

  //! size factor to allow fixed size arrays
  ///
  /// to allow fixed size arrays for a unknown number of unknowns, we make them bigger than necessary
  /// this factor is multiplied times numnode(distype) to get the size of many arrays
  template<XFEM::AssemblyType ASSTYPE>
  struct SizeFac {};
  /// specialization of SizeFac for XFEM::standard_assembly
  template<> struct SizeFac<XFEM::standard_assembly> {static const std::size_t fac = 1;};
  /// specialization of SizeFac for XFEM::xfem_assembly
  template<> struct SizeFac<XFEM::xfem_assembly>     {static const std::size_t fac = 3;};


  //! fill a number of local (element) arrays with unknown values
  //! from the global unknown vector given by the discretization
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            class M1, class M2>
  void fillElementUnknownsArrays4(
          const XFEM::ElementDofManager& dofman,
          const DRT::ELEMENTS::XDiff3::MyState& mystate,
          M1& evelnp,
//          M1& eveln,
//          M1& evelnm,
//          M1& eaccn,
          M2& eflux
          )
  {

      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
      const size_t numparamT = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Temp);
      // put one here to create arrays of size 1, since they are not needed anyway
      // in the xfem assembly, the numparam is determined by the dofmanager
      const size_t numparamHeatFlux_x = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::HeatFlux_x);

      const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
      const DRT::Element::DiscretizationType stressdistype = XDIFF::StressInterpolation3D<DISTYPE>::distype;
      const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      if (numparamT > shpVecSize)
      {
        dserror("increase SizeFac for nodal unknowns");
      }
      if (numparamHeatFlux_x > shpVecSizeStress)
      {
        dserror("increase SizeFac for stress unknowns");
      }

      const std::vector<int>& tempdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Temp>());

      for (size_t iparam=0; iparam<numparamT; ++iparam)
      {
          evelnp(iparam) = mystate.velnp[tempdof[iparam]];
//          if (mystate.instationary)
//          {
//              eveln( iparam) = mystate.veln[ tempdof[iparam]];
//              evelnm(iparam) = mystate.velnm[tempdof[iparam]];
//              eaccn( iparam) = mystate.accn[ tempdof[iparam]];
//          }
      }
      const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::HeatFlux_x, 0) > 0);
      if (tauele_unknowns_present)
      {
          const size_t numparamHeatFlux_y = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::HeatFlux_y, 1);
          const size_t numparamHeatFlux_z = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::HeatFlux_z, 1);
          const std::vector<int>& tauxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::HeatFlux_x>());
          const std::vector<int>& tauyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::HeatFlux_y>());
          const std::vector<int>& tauzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::HeatFlux_z>());
          for (size_t iparam=0; iparam<numparamHeatFlux_x; ++iparam)   eflux(0,iparam) = mystate.velnp[tauxxdof[iparam]];
          for (size_t iparam=0; iparam<numparamHeatFlux_y; ++iparam)   eflux(1,iparam) = mystate.velnp[tauyydof[iparam]];
          for (size_t iparam=0; iparam<numparamHeatFlux_z; ++iparam)   eflux(2,iparam) = mystate.velnp[tauzzdof[iparam]];
      }
  }

  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            size_t NUMDOF,
            size_t shpVecSize,
            size_t shpVecSizeStress>
  void BuildStiffnessMatrixEntries(
      LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>&           assembler,
      const XFEM::ApproxFunc<shpVecSize>&                     shp,
      const LINALG::Matrix<shpVecSizeStress,1>&  shp_tau,
      const double& fac,
      const double& timefac,
      const double& timefacfac,
      const double& visc,
      const LINALG::Matrix<3,1>& Tderxy,
      const double&              rhsint,
      const LINALG::Matrix<3,1>& heatflux,
      const bool tauele_unknowns_present,
      const bool instationary
        )
  {

  /*
                /                        \
               |       / \         /  \   |
               |  eps | v | , tau | Du |  |
               |       \ /         \  /   |
                \                        /
  */
  assembler.template Matrix<Temp,Temp>(shp.dx,     visc*timefacfac, shp.dx);
  assembler.template Matrix<Temp,Temp>(shp.dy,     visc*timefacfac, shp.dy);
  assembler.template Matrix<Temp,Temp>(shp.dz,     visc*timefacfac, shp.dz);

  assembler.template Vector<Temp>(shp.dx,     -visc*timefacfac*Tderxy(0));
  assembler.template Vector<Temp>(shp.dy,     -visc*timefacfac*Tderxy(1));
  assembler.template Vector<Temp>(shp.dz,     -visc*timefacfac*Tderxy(2));


  // source term of the right hand side
  assembler.template Vector<Temp>(shp.d0,     -rhsint*timefacfac);


  // Hellinger-Reissner terms
  if (tauele_unknowns_present)
  {

                       /*                     \
                    - |  virt tau , eps(Dtau)  |
                       \                     */

      const double reciproke_viscfac = 1.0/(visc);
      assembler.template Matrix<HeatFlux_x,HeatFlux_x>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
      assembler.template Matrix<HeatFlux_y,HeatFlux_y>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
      assembler.template Matrix<HeatFlux_z,HeatFlux_z>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);

      assembler.template Vector<HeatFlux_x>(shp_tau,  reciproke_viscfac*timefacfac*heatflux(0));
      assembler.template Vector<HeatFlux_y>(shp_tau,  reciproke_viscfac*timefacfac*heatflux(1));
      assembler.template Vector<HeatFlux_z>(shp_tau,  reciproke_viscfac*timefacfac*heatflux(2));

                   /*                 \
                  | virt tau , eps(Du) |
                   \                 */

      assembler.template Matrix<HeatFlux_x,Temp>(shp_tau,     -timefacfac    , shp.dx);
      assembler.template Matrix<HeatFlux_y,Temp>(shp_tau,     -timefacfac    , shp.dy);
      assembler.template Matrix<HeatFlux_z,Temp>(shp_tau,     -timefacfac    , shp.dz);

      assembler.template Vector<HeatFlux_x>(shp_tau,     timefacfac*Tderxy(0));
      assembler.template Vector<HeatFlux_y>(shp_tau,     timefacfac*Tderxy(1));
      assembler.template Vector<HeatFlux_z>(shp_tau,     timefacfac*Tderxy(2));

  }

}



/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2>
void SysmatDomain4(
    const DRT::Element*                 ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
//    const M1&                           eveln,
//    const M1&                           evelnm,
//    const M1&                           eaccn,
    const M2&                           eflux,
    Teuchos::RCP<const MAT::Material>   material,      ///< fluid material
    const FLUID_TIMEINTTYPE             timealgo,      ///< time discretization type
    const double                        dt,            ///< delta t (time step size)
    const double                        theta,         ///< factor for one step theta scheme
    const bool                          newton,        ///< full Newton or fixed-point-like
    const bool                          pstab,         ///< flag for stabilization
    const bool                          supg,          ///< flag for stabilization
    const bool                          cstab,         ///< flag for stabilization
    const bool                          instationary,  ///< switch between stationary and instationary formulation
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler,
    double&                             L2
)
{
    // number of nodes for element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // space dimension for 3d fluid element
    const size_t nsd = 3;

    // time integration constant
    const double timefac = 1.0;//FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // get node coordinates of the current element
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    // get older interface velocities and accelerations
//    const Epetra_Vector& ivelcoln  = *ih->cutterdis()->GetState("ivelcoln");
//    const Epetra_Vector& ivelcolnm = *ih->cutterdis()->GetState("ivelcolnm");
//    const Epetra_Vector& iacccoln  = *ih->cutterdis()->GetState("iacccoln");

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    // get viscosity
    // check here, if we really have a fluid !!
    dsassert(material->MaterialType() == INPAR::MAT::m_fluid, "Material law is not of type m_fluid.");
    const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(material.get());
    const double visc = actmat->Viscosity();

    // flag for higher order elements
//    const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

    const DRT::Element::DiscretizationType stressdistype = XDIFF::StressInterpolation3D<DISTYPE>::distype;

    // figure out whether we have stress unknowns at all
    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, HeatFlux_x, 0) > 0);
//    const bool velocity_unknowns_present = (getNumParam<ASSTYPE>(dofman, Velx, 1) > 0);
//    const bool pressure_unknowns_present = (getNumParam<ASSTYPE>(dofman, Pres, 1) > 0);
//    cout << endl;
//    if (ASSTYPE == XFEM::standard_assembly)
//        cout << "standard assembly" << endl;
//    else
//        cout << "xfem assembly" << endl;
//
//    cout << "stress unknowns present  : " << stress_unknowns_present << endl;
//    cout << "velocity unknowns present: " << velocity_unknowns_present << endl;
//    cout << "pressure unknowns present: " << pressure_unknowns_present << endl;


    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    //const int numparamvelx = getNumParam<ASSTYPE>(dofman, Velx, numnode);
    const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Temp);
//    const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::HeatFlux_x);

    // information about domain integration cells
    const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));
    //cout << "Element "<< ele->Id() << ": ";
    // loop over integration cells
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {
        const LINALG::Matrix<nsd,1> cellcenter_xyz(cell->GetPhysicalCenterPosition());

        int labelnp = 0;

        if (ASSTYPE == XFEM::xfem_assembly)
        {
          // integrate only in fluid integration cells (works well only with void enrichments!!!)
          labelnp = ih->PositionWithinConditionNP(cellcenter_xyz);
          const std::set<int> xlabelset(dofman.getUniqueEnrichmentLabels());
          bool compute = false;
          if (labelnp == 0) // fluid
          {
            compute = true;
          }
          if (not compute)
          {
            continue; // next integration cell
          }
        }

        const XFEM::ElementEnrichmentValues enrvals(
              *ele,
              ih,
              dofman,
              cellcenter_xyz,
              XFEM::Enrichment::approachUnknown);

        const DRT::UTILS::GaussRule3D gaussrule = XDIFF::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape());

        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            LINALG::Matrix<nsd,1> pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            LINALG::Matrix<nsd,1> posXiDomain;
            GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);

            // shape functions and their first derivatives
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

            const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
            const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

            static XFEM::ApproxFunc<shpVecSize> shp;

            static LINALG::Matrix<shpVecSizeStress,1>   shp_tau;

            if (ASSTYPE == XFEM::xfem_assembly)
            {
                // temporary arrays
                static LINALG::Matrix<shpVecSize,1> enr_funct;
                static LINALG::Matrix<3,shpVecSize> enr_derxy;

//                const std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
//                      ih,
//                      dofman.getUniqueEnrichments(),
//                      gauss_pos_xyz,
//                      XFEM::Enrichment::approachUnknown));

                // shape function for nodal dofs
                enrvals.ComputeEnrichedNodalShapefunction(
                        Temp,
                        funct,
                        derxy,
                        enr_funct,
                        enr_derxy);

                for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
                {
                  shp.d0(iparam) = enr_funct(iparam);
                  shp.dx(iparam) = enr_derxy(0,iparam);
                  shp.dy(iparam) = enr_derxy(1,iparam);
                  shp.dz(iparam) = enr_derxy(2,iparam);
                }



                if (tauele_unknowns_present)
                {
                    LINALG::Matrix<shpVecSizeStress,1> enr_funct_stress;

                    // shape functions for element dofs
                    enrvals.ComputeEnrichedElementShapefunction(
                            HeatFlux_x,
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
            else // standard assembly
            {
              // -> numparamvelx == numnode
              for (size_t iparam = 0; iparam < numnode; ++iparam)
              {
                shp.d0(iparam) = funct(iparam);
                shp.dx(iparam) = derxy(0,iparam);
                shp.dy(iparam) = derxy(1,iparam);
                shp.dz(iparam) = derxy(2,iparam);
              }

              if (tauele_unknowns_present)
              {
                dserror("no stress enrichments without xfem assembly");
              }
            }

            // get velocities and accelerations at integration point
            const double gpvelnp = XDIFF::interpolateScalarFieldToIntPoint(evelnp, shp.d0, numparamvelx);
//            LINALG::Matrix<nsd,1> gpveln  = XDIFF::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
//            LINALG::Matrix<nsd,1> gpvelnm = XDIFF::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
//            LINALG::Matrix<nsd,1> gpaccn  = XDIFF::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);


//            if (ASSTYPE == XFEM::xfem_assembly and timealgo != timeint_stationary)
//            {
//              const bool valid_spacetime_cell_found = XFLUID::modifyOldTimeStepsValues<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, ivelcoln, ivelcolnm, iacccoln, gpveln, gpvelnm, gpaccn);
//              if (not valid_spacetime_cell_found)
//                continue;
//            }
//            cout << gpvelnp << endl;
//            cout << evelnp << endl;
//            cout << shp << endl;

            // get history data (n) at integration point
//            LINALG::Matrix<3,1> histvec;
//            //histvec = enr_funct(j)*evelnp_hist(i,j);
//            for (int isd = 0; isd < nsd; ++isd)
//            {
//                histvec(isd) = 0.0;
//                for (int iparam = 0; iparam < numparamvelx; ++iparam)
//                    histvec(isd) += evelnp_hist(isd,iparam)*shp.d0(iparam);
//            }
//            const LINALG::Matrix<nsd,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
//                gpveln, gpvelnm, gpaccn, timealgo, dt, theta);

            // get velocity (np,i) derivatives at integration point
            // vderxy = enr_derxy(j,k)*evelnp(i,k);


            //cout << "eps_xy" << (0.5*(vderxy(0,1)+vderxy(1,0))) << ", "<< endl;



            // get Temperature gradients
            // gradp = enr_derxy(i,j)*eprenp(j);
            LINALG::Matrix<nsd,1> Tderxy(true);
            for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
            {
              Tderxy(0) += shp.dx(iparam)*evelnp(iparam);
              Tderxy(1) += shp.dy(iparam)*evelnp(iparam);
              Tderxy(2) += shp.dz(iparam)*evelnp(iparam);
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
//            double pres = 0.0;
//            for (size_t iparam = 0; iparam != numparampres; ++iparam)
//              pres += shp.d0(iparam)*eprenp(iparam);

            // get viscous stress unknowns
            static LINALG::Matrix<nsd,1> heatflux;
            if (tauele_unknowns_present)
            {
              XDIFF::fill_tau(numparamtauxx, shp_tau, eflux, heatflux);
            }
            else
            {
              heatflux.Clear();
            }


            // get bodyforce in gausspoint
//            LINALG::Matrix<3,1> bodyforce;
//            bodyforce = 0.0;
//            cout << bodyforce << endl;
            ///////////////LINALG::SerialDenseVector bodyforce_(enr_edeadng_(i,j)*enr_funct_(j));

            // compute stabilization parameters (3 taus)
//            double tau_stab_M  = 0.0;
//            double tau_stab_Mp = 0.0;
//            double tau_stab_C  = 0.0;
//              // function call does not work like this!
//            FLD::UTILS::computeStabilizationParams(derxy, gpvelnp, xji, numparamvelx,
//                instationary, visc, 1.0, hk, mk, timefac, dt, INPAR::FLUID::tautype_franca_barrenechea_valentin_wall,
//                tau_stab_M, tau_stab_Mp, tau_stab_C);



            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            double rhsint = -0.000;
//            LINALG::Matrix<nsd,1> bodyforce;
//            bodyforce.Clear();
//            //bodyforce(0) = 1.0;
//            for (size_t isd = 0; isd < nsd; ++isd)
//                rhsint(isd) = histvec(isd) + bodyforce(isd)*timefac;

            /*----------------- get numerical representation of single operators ---*/

            /* Convective term  u_old * grad u_old: */
//            LINALG::Matrix<nsd,1> conv_old;
            //conv_old = vderxy(i, j)*gpvelnp(j);
//            conv_old.Multiply(vderxy,gpvelnp);

//            /* Viscous term  div epsilon(u_old) */
//            LINALG::Matrix<nsd,1> visc_old;
//            visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
//            visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
//            visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

            // evaluate residual once for all stabilisation right hand sides
            double res_old = 0.0;
//            for (size_t isd = 0; isd < nsd; ++isd)
//                res_old(isd) = -rhsint(isd)+timefac*(conv_old(isd)+gradp(isd)-2.0*visc*visc_old(isd));

            if (instationary)
                res_old += gpvelnp;

            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

            LINALG::Matrix<3,1> physpos(true);
            GEO::elementToCurrentCoordinates(DISTYPE, xyze, posXiDomain, physpos);

#if 1
            // assume cylinder along z-axis
            const double radius = sqrt(physpos(0)*physpos(0) + physpos(1)*physpos(1));
            const double ri = 0.2;
            const double ra = 1.0;
            const double Ti = 5.0;
            const double Ta = 1.0;
            if (ri < radius and radius < ra)
            {
              const double T_exact = Ti -(log(radius) - log(ri))*(-(Ta-Ti)/(log(ra)-log(ri)));
              const double epsilon = (gpvelnp - T_exact);

              L2 += epsilon*epsilon*fac;
//              cout << T_exact << "   " << gpvelnp << endl;
            }
#endif
#if 0
            // block problem 3d
            const double x = physpos(0);
            const double y = physpos(1);
            const double z = physpos(2);
            const double a = 0.5;
            const double b = 1.0;
            const double c = 1.0;
            const double pi = PI;

//            rhsint = -2*sin(pi*y/b)*sin(2*pi*z/c) - x*pi*pi*(a - x)*sin(pi*y/b)*sin(2*pi*z/c)/(b*b) - 4*x*pi*pi*(a - x)*sin(pi*y/b)*sin(2*pi*z/c)/(c*c);
//            rhsint = -2*sin(pi*y/b)*sin(pi*z/c) - x*pi*pi*(a - x)*sin(pi*y/b)*sin(pi*z/c)/(b*b) - x*pi*pi*(a - x)*sin(pi*y/b)*sin(pi*z/c)/(c*c);
            rhsint = -pi*pi*sin(pi*x/a)*sin(pi*y/b)*sin(pi*z/c)/(a*a) - pi*pi*sin(pi*x/a)*sin(pi*y/b)*sin(pi*z/c)/(b*b) - pi*pi*sin(pi*x/a)*sin(pi*y/b)*sin(pi*z/c)/(c*c);
            if (0.0<x and x<a and 0.0<y and y<b and 0.0<z and z<c)
            {
              const double T_exact = sin(x/a*pi)*sin(y/b*pi)*sin(z/c*pi);
              const double epsilon = (gpvelnp - T_exact);
              L2 += epsilon*epsilon*fac;
            }
#endif
#if 0
            // block problem 2d
            const double x = physpos(0);
            const double y = physpos(1);
            const double a = 0.5;
            const double b = 1.0;
            const double pi = PI;

            rhsint = - (pi*pi*sin(pi*x/a)*sin(pi*y/b)/(a*a) + pi*pi*sin(pi*x/a)*sin(pi*y/b)/(b*b));


            if (0.0<x and x<a and 0.0<y and y<b)
            {
              const double T_exact = sin(pi*x/a)*sin(pi*y/b);
              const double epsilon = (gpvelnp - T_exact);
              L2 += epsilon*epsilon*fac;
            }
#endif


            BuildStiffnessMatrixEntries<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
                assembler, shp, shp_tau, fac, timefac, timefacfac, visc,
                Tderxy, rhsint, heatflux,
                tauele_unknowns_present, instationary);

        } // end loop over gauss points
    } // end loop over integration cells

    return;
}

/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2>
void SysmatBoundary4(
    const DRT::Element*               ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const M2&                         eflux,
    const Teuchos::RCP<const Epetra_Vector>& ivelcol,       ///< velocity for interface nodes
    const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
    const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
    const double&                     dt,            ///< delta t (time step size)
    const double&                     theta,         ///< factor for one step theta scheme
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const bool                        ifaceForceContribution
)
{
    if (ASSTYPE != XFEM::xfem_assembly) dserror("works only with xfem assembly");

    const size_t nsd = 3;
    const size_t numnodefix_boundary = 9;

    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // number of nodes for element
    const size_t numnode_xele = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const size_t numparamvelx = XFEM::NumParam<numnode_xele,ASSTYPE>::get(dofman, XFEM::PHYSICS::Temp);
    //const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    //const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::HeatFlux_x);


    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, HeatFlux_x, 0) > 0);
    // for now, I don't try to compare to elements without stress unknowns, since they lock anyway
    if (tauele_unknowns_present)
    {

    // information about boundary integration cells
    const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());


    // loop over boundary integration cells
    for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
    {

        // gaussian points
        const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_37point);

        // get the right boundary element
        const DRT::Element* boundaryele = ih->GetBoundaryEle(cell->GetSurfaceEleGid());
        const std::size_t numnode_boundary = boundaryele->NumNode();

        // get current node coordinates
//        LINALG::SerialDenseMatrix xyze_boundary(nsd,numnode_boundary);
        static LINALG::Matrix<nsd,numnodefix_boundary> xyze_boundary;
        ih->fillBoundaryNodalPositionsNP(boundaryele, xyze_boundary);

        // get interface velocities at the boundary element nodes
//        LINALG::SerialDenseMatrix vel_boundary(nsd,numnode_boundary);
        LINALG::Matrix<nsd,numnodefix_boundary> vel_boundary;
        const DRT::Node*const* nodes = boundaryele->Nodes();
        {
          std::vector<double> myvel(nsd);
          std::vector<int> gdofs(nsd);
          for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
          {
            ih->cutterdis()->Dof(nodes[inode],0,gdofs);
            DRT::UTILS::ExtractMyValues(*ivelcol,myvel,gdofs);
            vel_boundary(0,inode) = myvel[0];
            vel_boundary(1,inode) = myvel[1];
            vel_boundary(2,inode) = myvel[2];
          }
        }

//        LINALG::SerialDenseMatrix force_boundary(3,numnode_boundary,true);
        LINALG::Matrix<nsd,numnodefix_boundary> force_boundary;
        force_boundary = 0.0;

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta^\boundary
            LINALG::Matrix<2,1> pos_eta_boundary;
            pos_eta_boundary(0) = intpoints.qxg[iquad][0];
            pos_eta_boundary(1) = intpoints.qxg[iquad][1];

            // coordinates of the current integration point in element coordinates \xi^\boundary
            LINALG::Matrix<2,1> posXiBoundary;
            mapEtaBToXiB(*cell, pos_eta_boundary, posXiBoundary);

            // coordinates of the current integration point in element coordinates \xi^\domain
            LINALG::Matrix<nsd,1> posXiDomain;
            mapEtaBToXiD(*cell, pos_eta_boundary, posXiDomain);

            const double detcell = fabs(detEtaBToXiB(*cell, pos_eta_boundary)); //TODO: check normals
            if (detcell < 0.0)
            {
              cout << "detcel :  " << detcell << endl;
              dserror("negative detcell! should be a bug!");
            }

            // shape functions and their first derivatives
            dsassert((int)numnodefix_boundary >= DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()),"More than 9 nodes for boundary element - change size of fixed size array!");

            //LINALG::SerialDenseVector funct_boundary(DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()));
            static LINALG::Matrix<numnodefix_boundary,1> funct_boundary;
            DRT::UTILS::shape_function_2D(funct_boundary, posXiBoundary(0),posXiBoundary(1),boundaryele->Shape());

            //LINALG::SerialDenseMatrix deriv_boundary(nsd, DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()));
            static LINALG::Matrix<nsd,numnodefix_boundary> deriv_boundary;
            DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posXiBoundary(0),posXiBoundary(1),boundaryele->Shape());

            // shape functions and their first derivatives
            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1> funct;
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

            // stress shape function
            const DRT::Element::DiscretizationType stressdistype = XDIFF::StressInterpolation3D<DISTYPE>::distype;
            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
            DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);

            // position of the gausspoint in physical coordinates
            // gauss_pos_xyz = funct_boundary(j)*xyze_boundary(i,j);
            const LINALG::Matrix<nsd,1> gauss_pos_xyz = XDIFF::interpolateVectorFieldToIntPoint(xyze_boundary,funct_boundary,numnode_boundary);

            // get jacobian matrix d x / d \xi  (3x2)
            // dxyzdrs(i,j) = xyze_boundary(i,k)*deriv_boundary(j,k);
            static LINALG::Matrix<nsd,2> dxyzdrs;
            dxyzdrs.Clear();
            //blas.GEMM('N','T',3,2,numnode_boundary,1.0,xyze_boundary.A(),xyze_boundary.LDA(),deriv_boundary.A(),0,0.0,dxyzdrs.A(),dxyzdrs.M());
            for (std::size_t k=0; k!=numnode_boundary; ++k)
              for (std::size_t i=0; i!=nsd; ++i)
                for (std::size_t j=0; j!=2; ++j)
                  dxyzdrs(i,j) += xyze_boundary(i,k)*deriv_boundary(j,k);

            // compute covariant metric tensor G for surface element (2x2)
            // metric = dxyzdrs(k,i)*dxyzdrs(k,j);
            static LINALG::Matrix<2,2> metric;
            metric.Clear();
            metric.MultiplyTN(dxyzdrs,dxyzdrs);

            // compute global derivates
            //const LINALG::SerialDenseMatrix derxy(dxyzdrs(i,k)*deriv_boundary(k,j));
            //const LINALG::SerialDenseMatrix derxy_stress(xji(i,k)*deriv_stress(k,j));

            const double detmetric = sqrt(metric.Determinant());

            const double fac = intpoints.qwgt[iquad]*detmetric*detcell;
            if (fac < 0.0)
            {
              dserror("negative fac! should be a bug!");
            }

            const std::size_t shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
            const std::size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

            // temporary arrays
            static LINALG::Matrix<shpVecSize,1>       enr_funct;
            static LINALG::Matrix<shpVecSizeStress,1> enr_funct_stress;

//            if (dofman.getUniqueEnrichments().size() > 1)
//              dserror("for an intersected element, we assume only 1 enrichment for now!");
            const XFEM::ElementEnrichmentValues enrvals(
                  *ele,
                  ih,
                  dofman,
                  gauss_pos_xyz,
                  XFEM::Enrichment::approachFromPlus);

            // shape function for nodal dofs
            enrvals.ComputeEnrichedNodalShapefunction(
                    Temp,
                    funct,
                    enr_funct);

            // shape functions for element dofs
            enrvals.ComputeEnrichedElementShapefunction(
                    HeatFlux_x,
                    funct_stress,
                    enr_funct_stress);

            // perform integration for entire matrix and rhs
            static LINALG::Matrix<shpVecSize,1> shp;
            for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
            {
              shp(iparam) = enr_funct(iparam);
            }
            static LINALG::Matrix<shpVecSizeStress,1> shp_tau;
            for (std::size_t iparam = 0; iparam < numparamtauxx; ++iparam)
            {
              shp_tau(iparam) = enr_funct_stress(iparam);
            }

            // get normal vector (in physical coordinates) to surface element at integration point
            LINALG::Matrix<nsd,1> normalvec_solid;
            GEO::computeNormalToSurfaceElement(boundaryele->Shape(), xyze_boundary, posXiBoundary, normalvec_solid);
            LINALG::Matrix<nsd,1> normalvec_fluid(true);
            normalvec_fluid.Update(-1.0,normalvec_solid,0.0);

            // get velocities (n+g,i) at integration point
            // gpvelnp = evelnp(i,j)*shp(j);
            double gpvelnp = 0.0;
            for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
                gpvelnp += evelnp(iparam)*shp(iparam);

            // get interface velocity
            const int belegid = cell->GetSurfaceEleGid();
            const int label = ih->GetLabelPerBoundaryElementId(belegid);
            // block 2d/3d BC
            double interface_Temp = 0.0;
            // cylinder BC
            if (label == 1)
              interface_Temp = 1.0;
            else if (label == 2)
              interface_Temp = 5.0;


            // get viscous stress unknowns
            static LINALG::Matrix<nsd,1> heatflux;
            XDIFF::fill_tau(numparamtauxx, shp_tau, eflux, heatflux);

            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

               /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */

            assembler.template Matrix<HeatFlux_x,Temp>(shp_tau, timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<HeatFlux_y,Temp>(shp_tau, timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<HeatFlux_z,Temp>(shp_tau, timefacfac*normalvec_fluid(2), shp);


            assembler.template Vector<HeatFlux_x>(shp_tau, -timefacfac*normalvec_fluid(0)*gpvelnp);
            assembler.template Vector<HeatFlux_y>(shp_tau, -timefacfac*normalvec_fluid(1)*gpvelnp);
            assembler.template Vector<HeatFlux_z>(shp_tau, -timefacfac*normalvec_fluid(2)*gpvelnp);


               /*                            \
              |  (virt tau) * n^f , u^\iface  |
               \                            */

            assembler.template Vector<HeatFlux_x>(shp_tau, timefacfac*normalvec_fluid(0)*interface_Temp);
            assembler.template Vector<HeatFlux_y>(shp_tau, timefacfac*normalvec_fluid(1)*interface_Temp);
            assembler.template Vector<HeatFlux_z>(shp_tau, timefacfac*normalvec_fluid(2)*interface_Temp);


               /*               \
            - |  v , Dtau * n^f  |
               \               */

            assembler.template Matrix<Temp,HeatFlux_x>(shp, timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Temp,HeatFlux_y>(shp, timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Temp,HeatFlux_z>(shp, timefacfac*normalvec_fluid(2), shp_tau);

            double q_times_n = 0.0;
            for (std::size_t isd = 0; isd < nsd; ++isd)
              q_times_n += heatflux(isd)*normalvec_fluid(isd);
            //cout << "sigmaijnj : " << disctau_times_n << endl;
            assembler.template Vector<Temp>(shp, -timefacfac*q_times_n);
//
//            // here the interface force is integrated
//            // this is done using test shape functions of the boundary mesh
//            // hence, we can't use the local assembler here
//            for (size_t inode = 0; inode < numnode_boundary; ++inode)
//            {
//              force_boundary(0,inode) += funct_boundary(inode) * (disctau_times_n(0) * timefacfac);
//              force_boundary(1,inode) += funct_boundary(inode) * (disctau_times_n(1) * timefacfac);
//              force_boundary(2,inode) += funct_boundary(inode) * (disctau_times_n(2) * timefacfac);
//            }

        } // end loop over gauss points

        // here we need to assemble into the global force vector of the boundary discretization
        // note that we assemble into a overlapping vector, hence we add only, if we are a xfem row element
        // this way, we can later add all contributions together when exporting to interface row elements
//        if (ifaceForceContribution)
//        {
//          const Epetra_Map* dofcolmap = ih->cutterdis()->DofColMap();
//          std::vector<int> gdofs(3);
//          for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
//          {
//            ih->cutterdis()->Dof(nodes[inode],0,gdofs);
//            (*iforcecol)[dofcolmap->LID(gdofs[0])] += force_boundary(0,inode);
//            (*iforcecol)[dofcolmap->LID(gdofs[1])] += force_boundary(1,inode);
//            (*iforcecol)[dofcolmap->LID(gdofs[2])] += force_boundary(2,inode);
//          }
//        }


      } // end loop over boundary integration cells
    }
    return;
}

/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void Sysmat4(
        const DRT::Element*               ele,           ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const DRT::ELEMENTS::XDiff3::MyState&  mystate,  ///< element state variables
        const Teuchos::RCP<const Epetra_Vector>& ivelcol,       ///< velocity for interface nodes
        const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
        Epetra_SerialDenseMatrix&         estif,         ///< element matrix to calculate
        Epetra_SerialDenseVector&         eforce,        ///< element rhs to calculate
        Teuchos::RCP<const MAT::Material> material,      ///< fluid material
        const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
        const double                      dt,            ///< delta t (time step size)
        const double                      theta,         ///< factor for one step theta scheme
        const bool                        newton,        ///< full Newton or fixed-point-like
        const bool                        pstab,         ///< flag for stabilisation
        const bool                        supg,          ///< flag for stabilisation
        const bool                        cstab,         ///< flag for stabilisation
        const bool                        instationary,  ///< switch between stationary and instationary formulation
        const bool                        ifaceForceContribution,
        double&                           L2
        )
{
    // initialize arrays
    estif.Scale(0.0);
    eforce.Scale(0.0);

    const int NUMDOF = 1;

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

    // split velocity and pressure (and stress)
    const int shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    const DRT::Element::DiscretizationType stressdistype = XDIFF::StressInterpolation3D<DISTYPE>::distype;
    const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    LINALG::Matrix<shpVecSize,1> evelnp;
//    LINALG::Matrix<shpVecSize,1> eveln;
//    LINALG::Matrix<shpVecSize,1> evelnm;
//    LINALG::Matrix<shpVecSize,1> eaccn;
    LINALG::Matrix<3,shpVecSizeStress> etau;

    fillElementUnknownsArrays4<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, etau);

    SysmatDomain4<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, etau,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler, L2);

    if (ASSTYPE == XFEM::xfem_assembly)
    {
      SysmatBoundary4<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, etau, ivelcol, iforcecol,
          timealgo, dt, theta, assembler, ifaceForceContribution);
    }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XDIFF::callSysmat4(
        const XFEM::AssemblyType          assembly_type,
        const DRT::ELEMENTS::XDiff3*     ele,
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
        const XFEM::ElementDofManager&    eleDofManager,
        const DRT::ELEMENTS::XDiff3::MyState&  mystate,   ///< element state variables
        const Teuchos::RCP<const Epetra_Vector>& ivelcol,
        const Teuchos::RCP<Epetra_Vector>&  iforcecol,     ///< reaction force due to given interface velocity
        Epetra_SerialDenseMatrix&         estif,
        Epetra_SerialDenseVector&         eforce,
        Teuchos::RCP<const MAT::Material> material,
        const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
        const double                      dt,            ///< delta t (time step size)
        const double                      theta,         ///< factor for one step theta scheme
        const bool                        newton ,
        const bool                        pstab  ,
        const bool                        supg   ,
        const bool                        cstab  ,
        const bool                        instationary,
        const bool                        ifaceForceContribution,
        double&                           L2
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat4<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::hex20:
                Sysmat4<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::hex27:
                Sysmat4<DRT::Element::hex27,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::tet4:
                Sysmat4<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::tet10:
                Sysmat4<DRT::Element::tet10,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            default:
                dserror("standard_assembly Sysmat not templated yet");
        };
    }
    else
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat4<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::hex20:
                Sysmat4<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::hex27:
                Sysmat4<DRT::Element::hex27,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::tet4:
                Sysmat4<DRT::Element::tet4,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            case DRT::Element::tet10:
                Sysmat4<DRT::Element::tet10,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, L2);
                break;
            default:
                dserror("xfem_assembly Sysmat not templated yet");
        };
    }
}

#endif
#endif
