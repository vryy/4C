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

#include "xfluid3_sysmat.H"
#include "xfluid3_utils.H"
#include "xfluid3_spacetime_utils.H"
#include "fluid3_stabilization.H"
#include "xfluid3_local_assembler.H"
#include "xfluid3_local_assembler_ifacepatch.H"
#include "xfluid3_interpolation.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_xfem/spacetime_boundary.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_function.H"
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
            class M1, class V1, class M2>
  void fillElementUnknownsArrays4(
          const XFEM::ElementDofManager& dofman,
          const DRT::ELEMENTS::XFluid3::MyState& mystate,
          M1& evelnp,
          M1& eveln,
          M1& evelnm,
          M1& eaccn,
          V1& eprenp,
          M2& etau
          )
  {

      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
      //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
      const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
      const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Vely);
      const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velz);
      const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
      dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz), "assumption violation");
      // put one here to create arrays of size 1, since they are not needed anyway
      // in the xfem assembly, the numparam is determined by the dofmanager
      const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

      const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
      const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
      const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      if (numparamvelx > shpVecSize)
      {
        cout << "increase SizeFac for nodal unknowns" << endl;
      }
      if (numparamtauxx > shpVecSizeStress)
      {
        cout << "increase SizeFac for stress unknowns" << endl;
      }

      const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
      const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
      const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
      const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());

      for (size_t iparam=0; iparam<numparamvelx; ++iparam)
      {
          evelnp(0,iparam) = mystate.velnp[velxdof[iparam]];
          if (mystate.instationary)
          {
              eveln( 0,iparam) = mystate.veln[ velxdof[iparam]];
              evelnm(0,iparam) = mystate.velnm[velxdof[iparam]];
              eaccn( 0,iparam) = mystate.accn[ velxdof[iparam]];
          }
      }
      for (size_t iparam=0; iparam<numparamvely; ++iparam)
      {
          evelnp(1,iparam) = mystate.velnp[velydof[iparam]];
          if (mystate.instationary)
          {
              eveln( 1,iparam) = mystate.veln[ velydof[iparam]];
              evelnm(1,iparam) = mystate.velnm[velydof[iparam]];
              eaccn( 1,iparam) = mystate.accn[ velydof[iparam]];
          }
      }
      for (size_t iparam=0; iparam<numparamvelz; ++iparam)
      {
          evelnp(2,iparam) = mystate.velnp[velzdof[iparam]];
          if (mystate.instationary)
          {
              eveln( 2,iparam) = mystate.veln[ velzdof[iparam]];
              evelnm(2,iparam) = mystate.velnm[velzdof[iparam]];
              eaccn( 2,iparam) = mystate.accn[ velzdof[iparam]];
          }
      }
      for (size_t iparam=0; iparam<numparampres; ++iparam)
          eprenp(iparam) = mystate.velnp[presdof[iparam]];
      const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxx, 0) > 0);
      if (tauele_unknowns_present)
      {
          const size_t numparamtauyy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayy, 1);
          const size_t numparamtauzz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmazz, 1);
          const size_t numparamtauxy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxy, 1);
          const size_t numparamtauxz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxz, 1);
          const size_t numparamtauyz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayz, 1);
          const std::vector<int>& tauxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxx>());
          const std::vector<int>& tauyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayy>());
          const std::vector<int>& tauzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmazz>());
          const std::vector<int>& tauxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxy>());
          const std::vector<int>& tauxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxz>());
          const std::vector<int>& tauyzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayz>());
          for (size_t iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = mystate.velnp[tauxxdof[iparam]];
          for (size_t iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = mystate.velnp[tauyydof[iparam]];
          for (size_t iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = mystate.velnp[tauzzdof[iparam]];
          for (size_t iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = mystate.velnp[tauxydof[iparam]];
          for (size_t iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = mystate.velnp[tauxzdof[iparam]];
          for (size_t iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = mystate.velnp[tauyzdof[iparam]];
      }
  }


/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class M2>
void SysmatDomain4(
    const DRT::Element*                 ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
    const M1&                           eveln,
    const M1&                           evelnm,
    const M1&                           eaccn,
    const V1&                           eprenp,
    const M2&                           etau,
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

    // get node coordinates of the current element
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    // get interface velocities and accelerations
    const Epetra_Vector& ivelcolnp = *ih->cutterdis()->GetState("ivelcolnp");
    const Epetra_Vector& ivelcoln  = *ih->cutterdis()->GetState("ivelcoln");
    const Epetra_Vector& ivelcolnm = *ih->cutterdis()->GetState("ivelcolnm");
    const Epetra_Vector& iacccoln  = *ih->cutterdis()->GetState("iacccoln");

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    // get viscosity
    // check here, if we really have a fluid !!
    dsassert(material->MaterialType() == INPAR::MAT::m_fluid, "Material law is not of type m_fluid.");
    const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(material.get());
    const double visc = actmat->Viscosity();

    // flag for higher order elements
    const bool higher_order_ele = XFLUID::secondDerivativesAvailable<DISTYPE>();

    const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;

    // figure out whether we have stress unknowns at all
    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
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
    const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

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

        const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape());

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

            // coordinates of the current integration point in physical coordinates \xi
            LINALG::Matrix<nsd,1> posx_gp;
            GEO::elementToCurrentCoordinatesT<DISTYPE>(xyze, posXiDomain, posx_gp);

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

            const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
            const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

            static XFLUID::ApproxFunc<shpVecSize> shp;

            static LINALG::Matrix<shpVecSizeStress,1>   shp_tau;

            if (ASSTYPE == XFEM::xfem_assembly)
            {
                // temporary arrays
                static LINALG::Matrix<shpVecSize,1> enr_funct;
                static LINALG::Matrix<3,shpVecSize> enr_derxy;
                static LINALG::Matrix<6,shpVecSize> enr_derxy2;

//                const std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
//                      ih,
//                      dofman.getUniqueEnrichments(),
//                      gauss_pos_xyz,
//                      XFEM::Enrichment::approachUnknown));

                // shape function for nodal dofs
                enrvals.ComputeEnrichedNodalShapefunction(
                        Velx,
                        funct,
                        derxy,
                        derxy2,
                        enr_funct,
                        enr_derxy,
                        enr_derxy2);

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
            else // standard assembly
            {
              // -> numparamvelx == numnode
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

            // get velocities and accelerations at integration point
            const LINALG::Matrix<nsd,1> gpvelnp = XFLUID::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpveln  = XFLUID::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpvelnm = XFLUID::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpaccn  = XFLUID::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);


            const bool was_in_fluid = (ih->PositionWithinConditionN(posx_gp) == 0);

            XFLUID::TimeFormulation timeformulation = XFLUID::Eulerian;
            double dtstar = -10000.0;
//            double dtstar = dt;
            if (timealgo != timeint_stationary)
            {
              if (not was_in_fluid)
              {
                timeformulation = XFLUID::ReducedTimeStepSize;
                const bool valid_spacetime_cell_found = XFLUID::modifyOldTimeStepsValues<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, dt, ivelcolnp, ivelcoln, ivelcolnm, iacccoln, gpveln, gpvelnm, gpaccn, dtstar);
                if (not valid_spacetime_cell_found)
                {
                  cout << "not valid_spacetime_cell_found" << endl;
                  continue;
                }
              }
              else
              {
                timeformulation = XFLUID::Eulerian;
                dtstar = dt;
              }
            }
            else
            {
              dtstar = dt;
            }

            if (dtstar == -10000.0)
              dserror("something went wrong!");

            // time integration constant
            const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dtstar, theta);


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
            const LINALG::Matrix<nsd,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
                gpveln, gpvelnm, gpaccn, timealgo, dtstar, theta);

            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

            //----------------------------------------------------------------------
            //                            GALERKIN PART

            if (instationary)
            {
                /* inertia (contribution to mass matrix) */
                /*
                                     /          \
                                    |            |
                                    |  Du^n , v  |
                                    |            |
                                     \          /
                */

                assembler.template Vector<Velx>(shp.d0, -fac*gpveln(0));
                assembler.template Vector<Vely>(shp.d0, -fac*gpveln(1));
                assembler.template Vector<Velz>(shp.d0, -fac*gpveln(2));

                /* inertia (contribution to mass matrix) */
                /*
                                     /          \
                                    |            |
                                    |  Da^n , v  |
                                    |            |
                                     \          /
                */

                assembler.template Vector<Velx>(shp.d0, -fac*dt*(1-theta)*gpaccn(0));
                assembler.template Vector<Vely>(shp.d0, -fac*dt*(1-theta)*gpaccn(1));
                assembler.template Vector<Velz>(shp.d0, -fac*dt*(1-theta)*gpaccn(2));
            }
        } // end loop over gauss points
    } // end loop over integration cells

    return;
}

/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void TimeRHS(
        const DRT::Element*               ele,           ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const DRT::ELEMENTS::XFluid3::MyState&  mystate,  ///< element state variables
        const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
        Epetra_SerialDenseMatrix&         estif,         ///< element matrix to calculate
        Epetra_SerialDenseVector&         eforce,        ///< element rhs to calculate
        Epetra_SerialDenseMatrix&         Gds,
        Epetra_SerialDenseVector&         rhsd,
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
        const bool                        monolithic_FSI,
        double&                           L2
        )
{
    // initialize arrays
    estif.Scale(0.0);
    eforce.Scale(0.0);

    const int NUMDOF = 4;

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

    // split velocity and pressure (and stress)
    const int shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
    const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    LINALG::Matrix<shpVecSize,1> eprenp;
    LINALG::Matrix<3,shpVecSize> evelnp;
    LINALG::Matrix<3,shpVecSize> eveln;
    LINALG::Matrix<3,shpVecSize> evelnm;
    LINALG::Matrix<3,shpVecSize> eaccn;
    LINALG::Matrix<6,shpVecSizeStress> etau;

    fillElementUnknownsArrays4<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp, etau);

    SysmatDomain4<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, etau,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler, L2);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUID::callTimeRHS(
        const XFEM::AssemblyType          assembly_type,
        const DRT::ELEMENTS::XFluid3*     ele,
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
        const XFEM::ElementDofManager&    eleDofManager,
        const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
        const Teuchos::RCP<Epetra_Vector>&  iforcecol,     ///< reaction force due to given interface velocity
        Epetra_SerialDenseMatrix&         estif,
        Epetra_SerialDenseVector&         eforce,
        Epetra_SerialDenseMatrix&         Gds,
        Epetra_SerialDenseVector&         rhsd,
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
        const bool                        monolithic_FSI,
        double&                           L2
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                TimeRHS<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex20:
                TimeRHS<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex27:
                TimeRHS<DRT::Element::hex27,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::tet4:
                TimeRHS<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::tet10:
                TimeRHS<DRT::Element::tet10,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
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
                TimeRHS<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex20:
                TimeRHS<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex27:
                TimeRHS<DRT::Element::hex27,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::tet4:
                TimeRHS<DRT::Element::tet4,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
            case DRT::Element::tet10:
                TimeRHS<DRT::Element::tet10,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            default:
                dserror("xfem_assembly Sysmat not templated yet");
        };
    }
}

#endif
#endif
