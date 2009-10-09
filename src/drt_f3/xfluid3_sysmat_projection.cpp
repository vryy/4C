/*----------------------------------------------------------------------*/
/*!
\file xfluid3_sysmat_projection.cpp

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
            class M1, class V1>
  void fillElementUnknownsArrays4(
          const XFEM::ElementDofManager& dofman,
          const DRT::ELEMENTS::XFluid3::MyState& mystate,
          M1& evelnp,
          M1& eveln,
          M1& evelnm,
          M1& eaccn,
          V1& eprenp
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
//      const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

      const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
//      const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
//      const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      if (numparamvelx > shpVecSize)
      {
        cout << "increase SizeFac for nodal unknowns" << endl;
      }
//      if (numparamtauxx > shpVecSizeStress)
//      {
//        cout << "increase SizeFac for stress unknowns" << endl;
//      }

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
//      const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxx, 0) > 0);
//      if (tauele_unknowns_present)
//      {
//          const size_t numparamtauyy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayy, 1);
//          const size_t numparamtauzz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmazz, 1);
//          const size_t numparamtauxy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxy, 1);
//          const size_t numparamtauxz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxz, 1);
//          const size_t numparamtauyz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayz, 1);
//          const std::vector<int>& tauxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxx>());
//          const std::vector<int>& tauyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayy>());
//          const std::vector<int>& tauzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmazz>());
//          const std::vector<int>& tauxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxy>());
//          const std::vector<int>& tauxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxz>());
//          const std::vector<int>& tauyzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayz>());
//          for (size_t iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = mystate.velnp[tauxxdof[iparam]];
//          for (size_t iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = mystate.velnp[tauyydof[iparam]];
//          for (size_t iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = mystate.velnp[tauzzdof[iparam]];
//          for (size_t iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = mystate.velnp[tauxydof[iparam]];
//          for (size_t iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = mystate.velnp[tauxzdof[iparam]];
//          for (size_t iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = mystate.velnp[tauyzdof[iparam]];
//      }
  }

  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            size_t NUMDOF,
            size_t shpVecSize>
  void BuildStiffnessMatrixEntries(
      LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>&   assembler,
      const XFLUID::ApproxFunc<shpVecSize>&     shp,
      const double&                             fac,
      const double&                             visc,
      const LINALG::Matrix<3,1>&                u2_proj,
      const bool                                pstab,
      const double&                             tau_stab_Mp
        )
  {

  //----------------------------------------------------------------------
  //                            GALERKIN PART

  // linear problem -> no incremental formulation

  /* inertia (contribution to mass matrix) */
  /*
                       /          \
                      |            |
                      |  du2 , u2  |
                      |            |
                       \          /
  */
  assembler.template Matrix<Velx,Velx>(shp.d0, fac, shp.d0);
  assembler.template Matrix<Vely,Vely>(shp.d0, fac, shp.d0);
  assembler.template Matrix<Velz,Velz>(shp.d0, fac, shp.d0);


  /* inertia (contribution to mass matrix) */
  /*
                       /               \
                      |                 |
                      |  du2 , u2_proj  |
                      |                 |
                       \               /
  */

  assembler.template Vector<Velx>(shp.d0, fac*u2_proj(0));
  assembler.template Vector<Vely>(shp.d0, fac*u2_proj(1));
  assembler.template Vector<Velz>(shp.d0, fac*u2_proj(2));


   /* Druckterm */
  /*
                  /                \
                 |                  |
               - |  nabla o v , Dp  |
                 |                  |
                  \                /
  */
  assembler.template Matrix<Velx,Pres>(shp.dx, -fac, shp.d0);
  assembler.template Matrix<Vely,Pres>(shp.dy, -fac, shp.d0);
  assembler.template Matrix<Velz,Pres>(shp.dz, -fac, shp.d0);

  /* Divergenzfreiheit - continuity equation*/
  /*
                 /              \
                |                |
                | q , nabla o Du |
                |                |
                 \              /
  */
  assembler.template Matrix<Pres,Velx>(shp.d0, fac, shp.dx);
  assembler.template Matrix<Pres,Vely>(shp.d0, fac, shp.dy);
  assembler.template Matrix<Pres,Velz>(shp.d0, fac, shp.dz);

  //----------------------------------------------------------------------
  //                 PRESSURE STABILISATION PART
  if(pstab)
  {
      const double ttimetauMp = tau_stab_Mp * fac;


      /* pressure stabilisation: pressure( L_pres_p) */
      /*
                /                    \
               |                      |
               |  nabla q , nabla Dp  |
               |                      |
                \                    /
      */
      assembler.template Matrix<Pres,Pres>(shp.dx, ttimetauMp, shp.dx);
      assembler.template Matrix<Pres,Pres>(shp.dy, ttimetauMp, shp.dy);
      assembler.template Matrix<Pres,Pres>(shp.dz, ttimetauMp, shp.dz);
  }
}



/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1>
void SysmatDomainProjection(
    const DRT::Element*                 ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           eveln,
    const M1&                           eaccn,
    const bool                          pstab,         ///< flag for stabilization
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler_veln,
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler_accn
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
//    const Epetra_Vector& ivelcolnp = *ih->cutterdis()->GetState("ivelcolnp");
    const Epetra_Vector& ivelcoln  = *ih->cutterdis()->GetState("ivelcoln");
//    const Epetra_Vector& ivelcolnm = *ih->cutterdis()->GetState("ivelcolnm");
    const Epetra_Vector& iacccoln  = *ih->cutterdis()->GetState("iacccoln");

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    // get viscosity
    // check here, if we really have a fluid !!
//    dsassert(material->MaterialType() == INPAR::MAT::m_fluid, "Material law is not of type m_fluid.");
//    const MAT::NewtonianFluid* actmat = dynamic_cast<const MAT::NewtonianFluid*>(material.get());
//    const double visc = actmat->Viscosity();

//    const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;

    // figure out whether we have stress unknowns at all
    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
    if (tauele_unknowns_present)
    {
      dserror("no stress enrichments without xfem assembly");
    }

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);

    // stabilization parameter
    const double hk = FLD::UTILS::HK<DISTYPE>(eveln,xyze);
    const double mk = FLD::UTILS::MK<DISTYPE>();

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

        const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape(),true);

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

//            // discontinuous stress shape functions
//            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
//            if (ASSTYPE == XFEM::xfem_assembly)
//            {
//              if (tauele_unknowns_present)
//              {
//                DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
//              }
//              else
//              {
//                funct_stress.Clear();
//              }
//            }
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

            static XFLUID::ApproxFunc<shpVecSize> shp;

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
                        Velx,
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
            }

            // get velocities and accelerations at integration point
//            const LINALG::Matrix<nsd,1> gpvelnp = XFLUID::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> u2_proj  = XFLUID::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);\
//            cout << u2_proj << endl;
//            LINALG::Matrix<nsd,1> gpvelnm = XFLUID::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> a2_proj  = XFLUID::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);


            const bool was_in_fluid = (ih->PositionWithinConditionN(posx_gp) == 0);

            XFLUID::TimeFormulation timeformulation = XFLUID::Eulerian;
//            double dtstar = dt;

            if (not was_in_fluid)
            {
              timeformulation = XFLUID::ReducedTimeStepSize;
              const bool valid_spacetime_cell_found = XFLUID::ProjectSpaceTimeValuesToNewMesh<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, ivelcoln, iacccoln, u2_proj, a2_proj);
              if (not valid_spacetime_cell_found)
              {
                cout << "not valid_spacetime_cell_found" << endl;
                continue;
              }
            }
            else
            {
              timeformulation = XFLUID::Eulerian;
            }





            // compute stabilization parameters (3 taus)
            const double vel_norm = u2_proj.Norm2();
            const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, u2_proj, vel_norm, numparamvelx);
            double tau_stab_M  = 0.0;
            double tau_stab_Mp = 0.0;
            double tau_stab_C  = 0.0;
            FLD::UTILS::computeStabilizationParams(u2_proj, xji,
                false, 1.0, 1.0, vel_norm, strle, hk, mk, 1.0, 0.0, INPAR::FLUID::tautype_franca_barrenechea_valentin_wall,
                tau_stab_M, tau_stab_Mp, tau_stab_C);



            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

            /* inertia (contribution to mass matrix) */
            /*
                                 /          \
                                |            |
                                |  du2 , u2  |
                                |            |
                                 \          /
            */
            assembler_veln.template Matrix<Velx,Velx>(shp.d0, fac, shp.d0);
            assembler_veln.template Matrix<Vely,Vely>(shp.d0, fac, shp.d0);
            assembler_veln.template Matrix<Velz,Velz>(shp.d0, fac, shp.d0);


            /* inertia (contribution to mass matrix) */
            /*
                                 /               \
                                |                 |
                                |  du2 , u2_proj  |
                                |                 |
                                 \               /
            */

            assembler_veln.template Vector<Velx>(shp.d0, fac*u2_proj(0));
            assembler_veln.template Vector<Vely>(shp.d0, fac*u2_proj(1));
            assembler_veln.template Vector<Velz>(shp.d0, fac*u2_proj(2));


             /* Druckterm */
            /*
                            /                \
                           |                  |
                         - |  nabla o v , Dp  |
                           |                  |
                            \                /
            */
            assembler_veln.template Matrix<Velx,Pres>(shp.dx, -fac, shp.d0);
            assembler_veln.template Matrix<Vely,Pres>(shp.dy, -fac, shp.d0);
            assembler_veln.template Matrix<Velz,Pres>(shp.dz, -fac, shp.d0);

            /* Divergenzfreiheit - continuity equation*/
            /*
                           /              \
                          |                |
                          | q , nabla o Du |
                          |                |
                           \              /
            */
            assembler_veln.template Matrix<Pres,Velx>(shp.d0, fac, shp.dx);
            assembler_veln.template Matrix<Pres,Vely>(shp.d0, fac, shp.dy);
            assembler_veln.template Matrix<Pres,Velz>(shp.d0, fac, shp.dz);

            //----------------------------------------------------------------------
            //                 PRESSURE STABILISATION PART
            if(pstab)
            {
                const double ttimetauMp = tau_stab_Mp * fac;


                /* pressure stabilisation: pressure( L_pres_p) */
                /*
                          /                    \
                         |                      |
                         |  nabla q , nabla Dp  |
                         |                      |
                          \                    /
                */
                assembler_veln.template Matrix<Pres,Pres>(shp.dx, ttimetauMp, shp.dx);
                assembler_veln.template Matrix<Pres,Pres>(shp.dy, ttimetauMp, shp.dy);
                assembler_veln.template Matrix<Pres,Pres>(shp.dz, ttimetauMp, shp.dz);
            }

#if 0
            // acceleration

            /* inertia (contribution to mass matrix) */
            /*
                                 /          \
                                |            |
                                |  du2 , u2  |
                                |            |
                                 \          /
            */
            assembler_accn.template Matrix<Velx,Velx>(shp.d0, fac, shp.d0);
            assembler_accn.template Matrix<Vely,Vely>(shp.d0, fac, shp.d0);
            assembler_accn.template Matrix<Velz,Velz>(shp.d0, fac, shp.d0);


            /* inertia (contribution to mass matrix) */
            /*
                                 /               \
                                |                 |
                                |  du2 , u2_proj  |
                                |                 |
                                 \               /
            */

            assembler_accn.template Vector<Velx>(shp.d0, fac*a2_proj(0));
            assembler_accn.template Vector<Vely>(shp.d0, fac*a2_proj(1));
            assembler_accn.template Vector<Velz>(shp.d0, fac*a2_proj(2));


             /* Druckterm */
            /*
                            /                \
                           |                  |
                         - |  nabla o v , Dp  |
                           |                  |
                            \                /
            */
            assembler_accn.template Matrix<Velx,Pres>(shp.dx, -fac, shp.d0);
            assembler_accn.template Matrix<Vely,Pres>(shp.dy, -fac, shp.d0);
            assembler_accn.template Matrix<Velz,Pres>(shp.dz, -fac, shp.d0);

            /* Divergenzfreiheit - continuity equation*/
            /*
                           /              \
                          |                |
                          | q , nabla o Du |
                          |                |
                           \              /
            */
            assembler_accn.template Matrix<Pres,Velx>(shp.d0, fac, shp.dx);
            assembler_accn.template Matrix<Pres,Vely>(shp.d0, fac, shp.dy);
            assembler_accn.template Matrix<Pres,Velz>(shp.d0, fac, shp.dz);

            //----------------------------------------------------------------------
            //                 PRESSURE STABILISATION PART
            if(pstab)
            {
                const double ttimetauMp = tau_stab_Mp * fac;


                /* pressure stabilisation: pressure( L_pres_p) */
                /*
                          /                    \
                         |                      |
                         |  nabla q , nabla Dp  |
                         |                      |
                          \                    /
                */
                assembler_accn.template Matrix<Pres,Pres>(shp.dx, ttimetauMp, shp.dx);
                assembler_accn.template Matrix<Pres,Pres>(shp.dy, ttimetauMp, shp.dy);
                assembler_accn.template Matrix<Pres,Pres>(shp.dz, ttimetauMp, shp.dz);
            }
#endif
        } // end loop over gauss points
    } // end loop over integration cells

    return;
}

/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void SysmatProject(
        const DRT::Element*               ele,           ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const DRT::ELEMENTS::XFluid3::MyState&  mystate,  ///< element state variables
        Epetra_SerialDenseMatrix&         estif_veln,         ///< element matrix to calculate
        Epetra_SerialDenseMatrix&         estif_accn,         ///< element matrix to calculate
        Epetra_SerialDenseVector&         eforce_veln,        ///< element rhs to calculate
        Epetra_SerialDenseVector&         eforce_accn,        ///< element rhs to calculate
        const bool                        pstab,         ///< flag for stabilisation
        const bool                        ifaceForceContribution
        )
{
    // initialize arrays
    estif_veln.Scale(0.0);
    estif_accn.Scale(0.0);
    eforce_veln.Scale(0.0);
    eforce_accn.Scale(0.0);

    const int NUMDOF = 4;

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler_veln(dofman, estif_veln, eforce_veln);
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler_accn(dofman, estif_accn, eforce_accn);

    // split velocity and pressure (and stress)
    const int shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<shpVecSize,1> eprenp;
    static LINALG::Matrix<3,shpVecSize> evelnp;
    static LINALG::Matrix<3,shpVecSize> eveln;
    static LINALG::Matrix<3,shpVecSize> evelnm;
    static LINALG::Matrix<3,shpVecSize> eaccn;

    fillElementUnknownsArrays4<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp);

    SysmatDomainProjection<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, eveln, eaccn,
        pstab, assembler_veln, assembler_accn);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUID::callSysmatProjection(
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
        const bool                        ifaceForceContribution
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                SysmatProject<DRT::Element::hex8,XFEM::standard_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::hex20:
                SysmatProject<DRT::Element::hex20,XFEM::standard_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::hex27:
                SysmatProject<DRT::Element::hex27,XFEM::standard_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::tet4:
                SysmatProject<DRT::Element::tet4,XFEM::standard_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::tet10:
                SysmatProject<DRT::Element::tet10,XFEM::standard_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
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
                SysmatProject<DRT::Element::hex8,XFEM::xfem_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::hex20:
                SysmatProject<DRT::Element::hex20,XFEM::xfem_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::hex27:
                SysmatProject<DRT::Element::hex27,XFEM::xfem_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            case DRT::Element::tet4:
                SysmatProject<DRT::Element::tet4,XFEM::xfem_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
            case DRT::Element::tet10:
                SysmatProject<DRT::Element::tet10,XFEM::xfem_assembly>(
                    ele, ih, eleDofManager, mystate, estif_veln, estif_accn, eforce_veln, eforce_accn,
                    pstab, ifaceForceContribution);
                break;
            default:
                dserror("xfem_assembly Sysmat not templated yet");
        };
    }
}

#endif
#endif
