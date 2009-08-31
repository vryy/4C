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

  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            size_t NUMDOF,
            size_t shpVecSize,
            size_t shpVecSizeStress>
  void BuildStiffnessMatrixEntries(
      LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>&           assembler,
      const XFLUID::ApproxFunc<shpVecSize>&                     shp,
      const LINALG::Matrix<shpVecSizeStress,1>&  shp_tau,
      const double& fac,
      const double& timefac,
      const double& timefacfac,
      const double& visc,
      const LINALG::Matrix<3,1>& gpvelnp,
      const double&             pres,
      const LINALG::Matrix<3,1>& gradp,
      const LINALG::Matrix<3,3>& vderxy,
      const LINALG::Matrix<3,1>& rhsint,
      const LINALG::Matrix<3,1>& res_old,
      const LINALG::Matrix<3,1>& visc_old,
      const LINALG::Matrix<3,3>& tau,
      const LINALG::Matrix<shpVecSize,1>& enr_conv_c_,
      const XFLUID::EnrViscs2<shpVecSize>& enr_viscs2,
      const bool tauele_unknowns_present,
      const bool instationary,
      const bool newton,
      const bool pstab,
      const bool supg,
      const bool cstab,
      const double& tau_stab_Mp,
      const double& tau_stab_M,
      const double& tau_stab_C
        )
  {

  //----------------------------------------------------------------------
  //                            GALERKIN PART

  if (instationary)
  {
      /* inertia (contribution to mass matrix) */
      /*
                           /        \
                          |          |
                          |  Du , v  |
                          |          |
                           \        /
      */
      assembler.template Matrix<Velx,Velx>(shp.d0, fac, shp.d0);
      assembler.template Matrix<Vely,Vely>(shp.d0, fac, shp.d0);
      assembler.template Matrix<Velz,Velz>(shp.d0, fac, shp.d0);

      assembler.template Vector<Velx>(shp.d0, -fac*gpvelnp(0));
      assembler.template Vector<Vely>(shp.d0, -fac*gpvelnp(1));
      assembler.template Vector<Velz>(shp.d0, -fac*gpvelnp(2));
  }

  /* convection, convective part */
  /*
               /                       \
              |      / n+1       \      |
              | v , | u   o nabla | Du  |
              |      \ (i)       /      |
               \                       /
  */
  assembler.template Matrix<Velx,Velx>(shp.d0, timefacfac, enr_conv_c_);
  assembler.template Matrix<Vely,Vely>(shp.d0, timefacfac, enr_conv_c_);
  assembler.template Matrix<Velz,Velz>(shp.d0, timefacfac, enr_conv_c_);

  assembler.template Vector<Velx>(shp.d0, -timefacfac*(gpvelnp(0)*vderxy(0,0) // check order
                                                   +gpvelnp(1)*vderxy(0,1)
                                                   +gpvelnp(2)*vderxy(0,2)));
  assembler.template Vector<Vely>(shp.d0, -timefacfac*(gpvelnp(0)*vderxy(1,0)
                                                   +gpvelnp(1)*vderxy(1,1)
                                                   +gpvelnp(2)*vderxy(1,2)));
  assembler.template Vector<Velz>(shp.d0, -timefacfac*(gpvelnp(0)*vderxy(2,0)
                                                   +gpvelnp(1)*vderxy(2,1)
                                                   +gpvelnp(2)*vderxy(2,2)));

  if (newton)
  {
      /*  convection, reactive part */
      /*
             /                         \
            |      /          \   n+1   |
            | v , | Du o nabla | u      |
            |      \          /   (i)   |
             \                         /
      */
      assembler.template Matrix<Velx,Velx>(shp.d0, timefacfac*vderxy(0,0), shp.d0);
      assembler.template Matrix<Velx,Vely>(shp.d0, timefacfac*vderxy(0,1), shp.d0);
      assembler.template Matrix<Velx,Velz>(shp.d0, timefacfac*vderxy(0,2), shp.d0);
      assembler.template Matrix<Vely,Velx>(shp.d0, timefacfac*vderxy(1,0), shp.d0);
      assembler.template Matrix<Vely,Vely>(shp.d0, timefacfac*vderxy(1,1), shp.d0);
      assembler.template Matrix<Vely,Velz>(shp.d0, timefacfac*vderxy(1,2), shp.d0);
      assembler.template Matrix<Velz,Velx>(shp.d0, timefacfac*vderxy(2,0), shp.d0);
      assembler.template Matrix<Velz,Vely>(shp.d0, timefacfac*vderxy(2,1), shp.d0);
      assembler.template Matrix<Velz,Velz>(shp.d0, timefacfac*vderxy(2,2), shp.d0);
  }

  /* Viskositaetsterm */
  /*
                /                        \
               |       / \         /  \   |
               |  eps | v | , tau | Du |  |
               |       \ /         \  /   |
                \                        /
  */
  assembler.template Matrix<Velx,Velx>(shp.dx, 2.0*visc*timefacfac, shp.dx);
  assembler.template Matrix<Velx,Velx>(shp.dy,     visc*timefacfac, shp.dy);
  assembler.template Matrix<Velx,Vely>(shp.dy,     visc*timefacfac, shp.dx);
  assembler.template Matrix<Velx,Velx>(shp.dz,     visc*timefacfac, shp.dz);
  assembler.template Matrix<Velx,Velz>(shp.dz,     visc*timefacfac, shp.dx);

  assembler.template Matrix<Vely,Vely>(shp.dx,     visc*timefacfac, shp.dx);
  assembler.template Matrix<Vely,Velx>(shp.dx,     visc*timefacfac, shp.dy);
  assembler.template Matrix<Vely,Vely>(shp.dy, 2.0*visc*timefacfac, shp.dy);
  assembler.template Matrix<Vely,Vely>(shp.dz,     visc*timefacfac, shp.dz);
  assembler.template Matrix<Vely,Velz>(shp.dz,     visc*timefacfac, shp.dy);

  assembler.template Matrix<Velz,Velz>(shp.dx,     visc*timefacfac, shp.dx);
  assembler.template Matrix<Velz,Velx>(shp.dx,     visc*timefacfac, shp.dz);
  assembler.template Matrix<Velz,Velz>(shp.dy,     visc*timefacfac, shp.dy);
  assembler.template Matrix<Velz,Vely>(shp.dy,     visc*timefacfac, shp.dz);
  assembler.template Matrix<Velz,Velz>(shp.dz, 2.0*visc*timefacfac, shp.dz);

  assembler.template Vector<Velx>(shp.dx,     -visc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
  assembler.template Vector<Velx>(shp.dy,     -visc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
  assembler.template Vector<Velx>(shp.dz,     -visc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));

  assembler.template Vector<Vely>(shp.dx,     -visc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
  assembler.template Vector<Vely>(shp.dy,     -visc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
  assembler.template Vector<Vely>(shp.dz,     -visc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));

  assembler.template Vector<Velz>(shp.dx,     -visc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
  assembler.template Vector<Velz>(shp.dy,     -visc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
  assembler.template Vector<Velz>(shp.dz,     -visc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));

  /* Druckterm */
  /*
                  /                \
                 |                  |
               - |  nabla o v , Dp  |
                 |                  |
                  \                /
  */
  assembler.template Matrix<Velx,Pres>(shp.dx, -timefacfac, shp.d0);
  assembler.template Matrix<Vely,Pres>(shp.dy, -timefacfac, shp.d0);
  assembler.template Matrix<Velz,Pres>(shp.dz, -timefacfac, shp.d0);

  assembler.template Vector<Velx>(shp.dx, timefacfac*pres);
  assembler.template Vector<Vely>(shp.dy, timefacfac*pres);
  assembler.template Vector<Velz>(shp.dz, timefacfac*pres);

  /* Divergenzfreiheit - continuity equation*/
  /*
                 /              \
                |                |
                | q , nabla o Du |
                |                |
                 \              /
  */
  assembler.template Matrix<Pres,Velx>(shp.d0, timefacfac, shp.dx);
  assembler.template Matrix<Pres,Vely>(shp.d0, timefacfac, shp.dy);
  assembler.template Matrix<Pres,Velz>(shp.d0, timefacfac, shp.dz);

  const double trace_gamma = (vderxy(0, 0) + vderxy(1, 1) + vderxy(2, 2));
  assembler.template Vector<Pres>(shp.d0, -timefacfac*trace_gamma);

  // source term of the right hand side
  assembler.template Vector<Velx>(shp.d0, fac*rhsint(0));
  assembler.template Vector<Vely>(shp.d0, fac*rhsint(1));
  assembler.template Vector<Velz>(shp.d0, fac*rhsint(2));


  // Hellinger-Reissner terms
  if (tauele_unknowns_present)
  {

                       /*                     \
                    - |  virt tau , eps(Dtau)  |
                       \                     */

      const double reciproke_viscfac = 1.0/(2.0*visc);
      assembler.template Matrix<Sigmaxx,Sigmaxx>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
      assembler.template Matrix<Sigmaxy,Sigmaxy>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
      assembler.template Matrix<Sigmaxz,Sigmaxz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
      assembler.template Matrix<Sigmayy,Sigmayy>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
      assembler.template Matrix<Sigmayz,Sigmayz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
      assembler.template Matrix<Sigmazz,Sigmazz>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);

      assembler.template Vector<Sigmaxx>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,0));
      assembler.template Vector<Sigmaxy>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,1)*2.0);
      assembler.template Vector<Sigmaxz>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,2)*2.0);
      assembler.template Vector<Sigmayy>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,1));
      assembler.template Vector<Sigmayz>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,2)*2.0);
      assembler.template Vector<Sigmazz>(shp_tau,  reciproke_viscfac*timefacfac*tau(2,2));

                   /*                 \
                  | virt tau , eps(Du) |
                   \                 */

      assembler.template Matrix<Sigmaxx,Velx>(shp_tau,     timefacfac    , shp.dx);
      assembler.template Matrix<Sigmaxy,Velx>(shp_tau,     timefacfac    , shp.dy);
      assembler.template Matrix<Sigmaxy,Vely>(shp_tau,     timefacfac    , shp.dx);
      assembler.template Matrix<Sigmaxz,Velx>(shp_tau,     timefacfac    , shp.dz);
      assembler.template Matrix<Sigmaxz,Velz>(shp_tau,     timefacfac    , shp.dx);
      assembler.template Matrix<Sigmayy,Vely>(shp_tau,     timefacfac    , shp.dy);
      assembler.template Matrix<Sigmayz,Vely>(shp_tau,     timefacfac    , shp.dz);
      assembler.template Matrix<Sigmayz,Velz>(shp_tau,     timefacfac    , shp.dy);
      assembler.template Matrix<Sigmazz,Velz>(shp_tau,     timefacfac    , shp.dz);

      assembler.template Vector<Sigmaxx>(shp_tau,    - timefacfac*vderxy(0, 0));
      assembler.template Vector<Sigmaxy>(shp_tau,    - timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
      assembler.template Vector<Sigmaxz>(shp_tau,    - timefacfac*(vderxy(0, 2) + vderxy(2, 0)));
      assembler.template Vector<Sigmayy>(shp_tau,    - timefacfac*vderxy(1, 1));
      assembler.template Vector<Sigmayz>(shp_tau,    - timefacfac*(vderxy(1, 2) + vderxy(1, 2)));
      assembler.template Vector<Sigmazz>(shp_tau,    - timefacfac*vderxy(2, 2));


      /* pressure-pressure coupling, rectangular part */
      /*
                     /                    \
                    |                      |
                  - | tr(virt tau^e) , p I |
                    |                      |
                     \                    /
      */
      assembler.template Matrix<Sigmaxx,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
      assembler.template Matrix<Sigmayy,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
      assembler.template Matrix<Sigmazz,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);

      assembler.template Vector<Sigmaxx>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
      assembler.template Vector<Sigmayy>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
      assembler.template Vector<Sigmazz>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);

  }

  //----------------------------------------------------------------------
  //                 PRESSURE STABILISATION PART
  if(pstab)
  {
      const double timetauMp  = timefac * tau_stab_Mp * fac;
      if (instationary)
      {
          /* pressure stabilisation: inertia */
          /*
                      /              \
                     |                |
                     |  Du , nabla q  |
                     |                |
                      \              /
          */
          assembler.template Matrix<Pres,Velx>(shp.dx, timetauMp, shp.d0);
          assembler.template Matrix<Pres,Vely>(shp.dy, timetauMp, shp.d0);
          assembler.template Matrix<Pres,Velz>(shp.dz, timetauMp, shp.d0);
      }
      const double ttimetauMp = timefac * timefac * tau_stab_Mp * fac;
      /* pressure stabilisation: convection, convective part */
      /*
                /                             \
               |             / n+1       \     |
               | nabla q ,  | u   o nabla | Du |
               |             \ i         /     |
                \                             /
      */
      assembler.template Matrix<Pres,Velx>(shp.dx, ttimetauMp, enr_conv_c_);
      assembler.template Matrix<Pres,Vely>(shp.dy, ttimetauMp, enr_conv_c_);
      assembler.template Matrix<Pres,Velz>(shp.dz, ttimetauMp, enr_conv_c_);

      if (newton)
      {
          /*  pressure stabilisation: convection, reactive part
                /                             \
               |           /          \   n+1  |
               | grad q , | Du o nabla | u     |
               |           \          /   (i)  |
                \                             /
          */
          assembler.template Matrix<Pres,Velx>(shp.dx, ttimetauMp*vderxy(0,0), shp.d0);
          assembler.template Matrix<Pres,Velx>(shp.dy, ttimetauMp*vderxy(1,0), shp.d0);
          assembler.template Matrix<Pres,Velx>(shp.dz, ttimetauMp*vderxy(2,0), shp.d0);

          assembler.template Matrix<Pres,Vely>(shp.dx, ttimetauMp*vderxy(0,1), shp.d0);
          assembler.template Matrix<Pres,Vely>(shp.dy, ttimetauMp*vderxy(1,1), shp.d0);
          assembler.template Matrix<Pres,Vely>(shp.dz, ttimetauMp*vderxy(2,1), shp.d0);

          assembler.template Matrix<Pres,Velz>(shp.dx, ttimetauMp*vderxy(0,2), shp.d0);
          assembler.template Matrix<Pres,Velz>(shp.dy, ttimetauMp*vderxy(1,2), shp.d0);
          assembler.template Matrix<Pres,Velz>(shp.dz, ttimetauMp*vderxy(2,2), shp.d0);
      }

      /* pressure stabilisation: viscosity (-L_visc_u) */
      /*
                 /                             \
                |                         /  \  |
              - |  nabla q , nabla o tau | Du | |
                |                         \  /  |
                 \                             /
      */
      assembler.template Matrix<Pres,Velx>(shp.dx, -2.0*visc*ttimetauMp, enr_viscs2.xx);
      assembler.template Matrix<Pres,Vely>(shp.dx, -2.0*visc*ttimetauMp, enr_viscs2.xy);
      assembler.template Matrix<Pres,Velz>(shp.dx, -2.0*visc*ttimetauMp, enr_viscs2.xz);

      assembler.template Matrix<Pres,Velx>(shp.dy, -2.0*visc*ttimetauMp, enr_viscs2.xy);
      assembler.template Matrix<Pres,Vely>(shp.dy, -2.0*visc*ttimetauMp, enr_viscs2.yy);
      assembler.template Matrix<Pres,Velz>(shp.dy, -2.0*visc*ttimetauMp, enr_viscs2.yz);

      assembler.template Matrix<Pres,Velx>(shp.dz, -2.0*visc*ttimetauMp, enr_viscs2.xz);
      assembler.template Matrix<Pres,Vely>(shp.dz, -2.0*visc*ttimetauMp, enr_viscs2.yz);
      assembler.template Matrix<Pres,Velz>(shp.dz, -2.0*visc*ttimetauMp, enr_viscs2.zz);

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

      // pressure stabilization
      assembler.template Vector<Pres>(shp.dx, -timetauMp*res_old(0));
      assembler.template Vector<Pres>(shp.dy, -timetauMp*res_old(1));
      assembler.template Vector<Pres>(shp.dz, -timetauMp*res_old(2));

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
                    /                        \
                   |        / n+1       \     |
                   |  Du , | u   o nabla | v  |
                   |        \ (i)       /     |
                    \                        /
          */
          assembler.template Matrix<Velx,Velx>(enr_conv_c_, timetauM, shp.d0);
          assembler.template Matrix<Vely,Vely>(enr_conv_c_, timetauM, shp.d0);
          assembler.template Matrix<Velz,Velz>(enr_conv_c_, timetauM, shp.d0);

          if (newton)
          {
              /* supg stabilisation: inertia, linearisation of testfunction  */
              /*
                         /                           \
                        |   n+1      /          \     |
                        |  u      , | Du o nabla | v  |
                        |   (i)      \          /     |
                         \                           /

              */
              assembler.template Matrix<Velx,Velx>(shp.dx, timetauM*gpvelnp(0), shp.d0);
              assembler.template Matrix<Velx,Vely>(shp.dy, timetauM*gpvelnp(0), shp.d0);
              assembler.template Matrix<Velx,Velz>(shp.dz, timetauM*gpvelnp(0), shp.d0);
              assembler.template Matrix<Vely,Velx>(shp.dx, timetauM*gpvelnp(1), shp.d0);
              assembler.template Matrix<Vely,Vely>(shp.dy, timetauM*gpvelnp(1), shp.d0);
              assembler.template Matrix<Vely,Velz>(shp.dz, timetauM*gpvelnp(1), shp.d0);
              assembler.template Matrix<Velz,Velx>(shp.dx, timetauM*gpvelnp(2), shp.d0);
              assembler.template Matrix<Velz,Vely>(shp.dy, timetauM*gpvelnp(2), shp.d0);
              assembler.template Matrix<Velz,Velz>(shp.dz, timetauM*gpvelnp(2), shp.d0);
          }
      }
      const double ttimetauM  = timefac * timefac * tau_stab_M * fac;
      /* supg stabilisation: convective part ( L_conv_u) */
      /*
           /                                          \
          |  / n+1        \        / n+1        \      |
          | | u    o nabla | v ,  | u    o nabla | Du  |
          |  \ (i)        /        \ (i)        /      |
           \                                          /
      */
      assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM, enr_conv_c_);
      assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM, enr_conv_c_);
      assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM, enr_conv_c_);
      /* supg stabilisation: pressure part  ( L_pres_p) */
      /*
                /                             \
               |   / n+1       \               |
               |  | u   o nabla | v , nabla Dp |
               |   \ (i)       /               |
                \                             /
      */
      assembler.template Matrix<Velx,Pres>(enr_conv_c_, ttimetauM, shp.dx);
      assembler.template Matrix<Vely,Pres>(enr_conv_c_, ttimetauM, shp.dy);
      assembler.template Matrix<Velz,Pres>(enr_conv_c_, ttimetauM, shp.dz);

      /* supg stabilisation: viscous part  (-L_visc_u) */
      /*
            /                                        \
           |               /  \    / n+1        \     |
         - |  nabla o eps | Du |, | u    o nabla | v  |
           |               \  /    \ (i)        /     |
            \                                        /
      */
      assembler.template Matrix<Velx,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.xx);
      assembler.template Matrix<Velx,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.xy);
      assembler.template Matrix<Velx,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.xz);

      assembler.template Matrix<Vely,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.yx);
      assembler.template Matrix<Vely,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.yy);
      assembler.template Matrix<Vely,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.yz);

      assembler.template Matrix<Velz,Velx>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.zx);
      assembler.template Matrix<Velz,Vely>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.zy);
      assembler.template Matrix<Velz,Velz>(enr_conv_c_, -2.0*visc*ttimetauM, enr_viscs2.zz);

      if (newton)
      {
          /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
          /*
                     /                                           \
                    |    /          \   n+1    / n+1        \     |
                    |   | Du o nabla | u    , | u    o nabla | v  |
                    |    \          /   (i)    \ (i)        /     |
                     \                                           /
          */
          assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM*vderxy(0,0), shp.d0);
          assembler.template Matrix<Velx,Vely>(enr_conv_c_, ttimetauM*vderxy(0,1), shp.d0);
          assembler.template Matrix<Velx,Velz>(enr_conv_c_, ttimetauM*vderxy(0,2), shp.d0);

          assembler.template Matrix<Vely,Velx>(enr_conv_c_, ttimetauM*vderxy(1,0), shp.d0);
          assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM*vderxy(1,1), shp.d0);
          assembler.template Matrix<Vely,Velz>(enr_conv_c_, ttimetauM*vderxy(1,2), shp.d0);

          assembler.template Matrix<Velz,Velx>(enr_conv_c_, ttimetauM*vderxy(2,0), shp.d0);
          assembler.template Matrix<Velz,Vely>(enr_conv_c_, ttimetauM*vderxy(2,1), shp.d0);
          assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM*vderxy(2,2), shp.d0);

          /*
                   /                                           \
                  |    / n+1        \   n+1    /          \     |
                  |   | u    o nabla | u    , | Du o nabla | v  |
                  |    \ (i)        /   (i)    \          /     |
                   \                                           /
          */
          const double con0 = ttimetauM*(gpvelnp(0)*vderxy(0,0) + gpvelnp(1)*vderxy(0,1) + gpvelnp(2)*vderxy(0,2));
          assembler.template Matrix<Velx,Velx>(shp.dx, con0, shp.d0);
          assembler.template Matrix<Velx,Vely>(shp.dy, con0, shp.d0);
          assembler.template Matrix<Velx,Velz>(shp.dz, con0, shp.d0);

          const double con1 = ttimetauM*(gpvelnp(0)*vderxy(1,0) + gpvelnp(1)*vderxy(1,1) + gpvelnp(2)*vderxy(1,2));
          assembler.template Matrix<Vely,Velx>(shp.dx, con1, shp.d0);
          assembler.template Matrix<Vely,Vely>(shp.dy, con1, shp.d0);
          assembler.template Matrix<Vely,Velz>(shp.dz, con1, shp.d0);

          const double con2 = ttimetauM*(gpvelnp(0)*vderxy(2,0) + gpvelnp(1)*vderxy(2,1) + gpvelnp(2)*vderxy(2,2));
          assembler.template Matrix<Velz,Velx>(shp.dx, con2, shp.d0);
          assembler.template Matrix<Velz,Vely>(shp.dy, con2, shp.d0);
          assembler.template Matrix<Velz,Velz>(shp.dz, con2, shp.d0);

          /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
          /*
                          /                               \
                         |         n+1    /          \     |
                         |  nabla p    , | Du o nabla | v  |
                         |         (i)    \          /     |
                          \                               /
          */
          assembler.template Matrix<Velx,Velx>(shp.dx, ttimetauM*gradp(0), shp.d0);
          assembler.template Matrix<Velx,Vely>(shp.dy, ttimetauM*gradp(0), shp.d0);
          assembler.template Matrix<Velx,Velz>(shp.dz, ttimetauM*gradp(0), shp.d0);

          assembler.template Matrix<Vely,Velx>(shp.dx, ttimetauM*gradp(1), shp.d0);
          assembler.template Matrix<Vely,Vely>(shp.dy, ttimetauM*gradp(1), shp.d0);
          assembler.template Matrix<Vely,Velz>(shp.dz, ttimetauM*gradp(1), shp.d0);

          assembler.template Matrix<Velz,Velx>(shp.dx, ttimetauM*gradp(2), shp.d0);
          assembler.template Matrix<Velz,Vely>(shp.dy, ttimetauM*gradp(2), shp.d0);
          assembler.template Matrix<Velz,Velz>(shp.dz, ttimetauM*gradp(2), shp.d0);

            /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
            /*
                    /                                         \
                   |               / n+1 \    /          \     |
                 - |  nabla o eps | u     |, | Du o nabla | v  |
                   |               \ (i) /    \          /     |
                    \                                         /
            */
          assembler.template Matrix<Velx,Velx>(shp.dx, -2.0*visc*ttimetauM*visc_old(0), shp.d0);
          assembler.template Matrix<Velx,Vely>(shp.dy, -2.0*visc*ttimetauM*visc_old(0), shp.d0);
          assembler.template Matrix<Velx,Velz>(shp.dz, -2.0*visc*ttimetauM*visc_old(0), shp.d0);

          assembler.template Matrix<Vely,Velx>(shp.dx, -2.0*visc*ttimetauM*visc_old(1), shp.d0);
          assembler.template Matrix<Vely,Vely>(shp.dy, -2.0*visc*ttimetauM*visc_old(1), shp.d0);
          assembler.template Matrix<Vely,Velz>(shp.dz, -2.0*visc*ttimetauM*visc_old(1), shp.d0);

          assembler.template Matrix<Velz,Velx>(shp.dx, -2.0*visc*ttimetauM*visc_old(2), shp.d0);
          assembler.template Matrix<Velz,Vely>(shp.dy, -2.0*visc*ttimetauM*visc_old(2), shp.d0);
          assembler.template Matrix<Velz,Velz>(shp.dz, -2.0*visc*ttimetauM*visc_old(2), shp.d0);

          /* supg stabilisation: bodyforce part, linearisation of test function */

          /*
                        /                             \
                       |              /          \     |
                     - |  rhsint   , | Du o nabla | v  |
                       |              \          /     |
                        \                             /

          */
          assembler.template Matrix<Velx,Velx>(shp.dx, -timetauM*rhsint(0), shp.d0);
          assembler.template Matrix<Velx,Vely>(shp.dy, -timetauM*rhsint(0), shp.d0);
          assembler.template Matrix<Velx,Velz>(shp.dz, -timetauM*rhsint(0), shp.d0);

          assembler.template Matrix<Vely,Velx>(shp.dx, -timetauM*rhsint(1), shp.d0);
          assembler.template Matrix<Vely,Vely>(shp.dy, -timetauM*rhsint(1), shp.d0);
          assembler.template Matrix<Vely,Velz>(shp.dz, -timetauM*rhsint(1), shp.d0);

          assembler.template Matrix<Velz,Velx>(shp.dx, -timetauM*rhsint(2), shp.d0);
          assembler.template Matrix<Velz,Vely>(shp.dy, -timetauM*rhsint(2), shp.d0);
          assembler.template Matrix<Velz,Velz>(shp.dz, -timetauM*rhsint(2), shp.d0);
      } // if newton

      // supg stabilisation
      assembler.template Vector<Velx>(enr_conv_c_, -timetauM*res_old(0));
      assembler.template Vector<Vely>(enr_conv_c_, -timetauM*res_old(1));
      assembler.template Vector<Velz>(enr_conv_c_, -timetauM*res_old(2));
  }


  //----------------------------------------------------------------------
  //                     STABILISATION, CONTINUITY PART
  if(cstab)
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
      assembler.template Matrix<Velx,Velx>(shp.dx, timefac_timefac_tau_C, shp.dx);
      assembler.template Matrix<Velx,Vely>(shp.dx, timefac_timefac_tau_C, shp.dy);
      assembler.template Matrix<Velx,Velz>(shp.dx, timefac_timefac_tau_C, shp.dz);

      assembler.template Matrix<Vely,Velx>(shp.dy, timefac_timefac_tau_C, shp.dx);
      assembler.template Matrix<Vely,Vely>(shp.dy, timefac_timefac_tau_C, shp.dy);
      assembler.template Matrix<Vely,Velz>(shp.dy, timefac_timefac_tau_C, shp.dz);

      assembler.template Matrix<Velz,Velx>(shp.dz, timefac_timefac_tau_C, shp.dx);
      assembler.template Matrix<Velz,Vely>(shp.dz, timefac_timefac_tau_C, shp.dy);
      assembler.template Matrix<Velz,Velz>(shp.dz, timefac_timefac_tau_C, shp.dz);

      assembler.template Vector<Velx>(shp.dx, -timefac_timefac_tau_C_divunp);
      assembler.template Vector<Vely>(shp.dy, -timefac_timefac_tau_C_divunp);
      assembler.template Vector<Velz>(shp.dz, -timefac_timefac_tau_C_divunp);
  } // endif cstab
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
    const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

    // stabilization parameter
    const double hk = XFLUID::HK<DISTYPE>(evelnp,xyze);
    const double mk = XFLUID::MK<DISTYPE>();

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

            // get velocity (np,i) derivatives at integration point
            // vderxy = enr_derxy(j,k)*evelnp(i,k);
            static LINALG::Matrix<3,nsd> vderxy;
            vderxy.Clear();
            for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
            {
              for (size_t isd = 0; isd < nsd; ++isd)
              {
                vderxy(isd,0) += evelnp(isd,iparam) * shp.dx(iparam);
                vderxy(isd,1) += evelnp(isd,iparam) * shp.dy(iparam);
                vderxy(isd,2) += evelnp(isd,iparam) * shp.dz(iparam);
              }
            }

            //cout << "eps_xy" << (0.5*(vderxy(0,1)+vderxy(1,0))) << ", "<< endl;

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
                  vderxy2(isd,0) += evelnp(isd,iparam)*shp.dxdx(iparam);
                  vderxy2(isd,1) += evelnp(isd,iparam)*shp.dydy(iparam);
                  vderxy2(isd,2) += evelnp(isd,iparam)*shp.dzdz(iparam);
                  vderxy2(isd,3) += evelnp(isd,iparam)*shp.dxdy(iparam);
                  vderxy2(isd,4) += evelnp(isd,iparam)*shp.dxdz(iparam);
                  vderxy2(isd,5) += evelnp(isd,iparam)*shp.dydz(iparam);
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
              gradp(0) += shp.dx(iparam)*eprenp(iparam);
              gradp(1) += shp.dy(iparam)*eprenp(iparam);
              gradp(2) += shp.dz(iparam)*eprenp(iparam);
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
              pres += shp.d0(iparam)*eprenp(iparam);

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


            // get bodyforce in gausspoint
//            LINALG::Matrix<3,1> bodyforce;
//            bodyforce = 0.0;
//            cout << bodyforce << endl;
            ///////////////LINALG::SerialDenseVector bodyforce_(enr_edeadng_(i,j)*enr_funct_(j));

            // compute stabilization parameters (3 taus)
            double tau_stab_M  = 0.0;
            double tau_stab_Mp = 0.0;
            double tau_stab_C  = 0.0;
            XFLUID::computeStabilization(derxy, gpvelnp, numparamvelx, instationary, visc, hk, mk, FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dtstar, theta),
                tau_stab_M, tau_stab_Mp, tau_stab_C);



            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            LINALG::Matrix<nsd,1> rhsint;
            LINALG::Matrix<nsd,1> bodyforce;
            bodyforce.Clear();
//            bodyforce(0) = 1.0;
            for (size_t isd = 0; isd < nsd; ++isd)
                rhsint(isd) = histvec(isd) + bodyforce(isd)*timefac;

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
            LINALG::Matrix<nsd,1> res_old;
            for (size_t isd = 0; isd < nsd; ++isd)
                res_old(isd) = -rhsint(isd)+timefac*(conv_old(isd)+gradp(isd)-2.0*visc*visc_old(isd));

            if (instationary)
                res_old += gpvelnp;

            /* Reactive term  u:  funct */
            /* linearise convective term */

            /*--- convective part u_old * grad (funct) --------------------------*/
            /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
             with  N .. form function matrix                                   */
            //const LINALG::SerialDenseVector enr_conv_c_(enr_derxy(j,i)*gpvelnp(j));
            static LINALG::Matrix<shpVecSize,1> enr_conv_c_;
            enr_conv_c_.Clear();
            for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
            {
                enr_conv_c_(iparam) += shp.dx(iparam)*gpvelnp(0);
                enr_conv_c_(iparam) += shp.dy(iparam)*gpvelnp(1);
                enr_conv_c_(iparam) += shp.dz(iparam)*gpvelnp(2);
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
            static XFLUID::EnrViscs2<shpVecSize> enr_viscs2;

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

#if 0
            // for Jeffery-Hamel Flow
            LINALG::Matrix<3,1> physpos(true);
            GEO::elementToCurrentCoordinates(DISTYPE, xyze, posXiDomain, physpos);
            
            const double x = physpos(0);
            const double y = physpos(1);
            
            const double alpha = atan(y/x);
            const double u_alpha = DRT::UTILS::JefferyHamelFlowFunction::RadialVelocity(alpha);
            
            const double nu = 1;
            const double u_exact_x = nu * (u_alpha/(x*x+y*y))*x;
            const double u_exact_y = nu * (u_alpha/(x*x+y*y))*y;

            if (1.0 < x and x < 2.0 and 0.0 < y and y < x)
            {
              const double epsilon_x = (gpvelnp(0) - u_exact_x);
              const double epsilon_y = (gpvelnp(1) - u_exact_y);

              L2 += (epsilon_x*epsilon_x + epsilon_y*epsilon_y)*fac;
            }
#endif

            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

            BuildStiffnessMatrixEntries<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
                assembler, shp, shp_tau, fac, timefac, timefacfac, visc,
                gpvelnp, pres, gradp, vderxy, rhsint, res_old, visc_old, tau,
                enr_conv_c_, enr_viscs2,
                tauele_unknowns_present, instationary, newton, pstab, supg, cstab,
                tau_stab_Mp, tau_stab_M, tau_stab_C);

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
    const M2&                         etau,
    const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
    Epetra_SerialDenseMatrix&         Gds,
    Epetra_SerialDenseVector&         rhsd,
    const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
    const double&                     dt,            ///< delta t (time step size)
    const double&                     theta,         ///< factor for one step theta scheme
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const bool                        ifaceForceContribution,
    const bool                        monolithic_FSI
)
{
    if (ASSTYPE != XFEM::xfem_assembly) dserror("works only with xfem assembly");

    const size_t nsd = 3;
    const size_t numnodefix_boundary = 9;

#if 0 
    // get node coordinates of the current element
    // number of nodes for element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);
#endif
    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // get interface velocities
    Teuchos::RCP<const Epetra_Vector> ivelcolnp = ih->cutterdis()->GetState("ivelcolnp");

    // number of nodes for element
    const size_t numnode_xele = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const size_t numparamvelx = XFEM::NumParam<numnode_xele,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    //const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    //const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Sigmaxx, 1);
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

    IFacePatchLocalAssembler<DISTYPE, NUMDOF> patchassembler(dofman);
    const std::set<int> begids = ih->GetIntersectingBoundaryElementsGID(ele->Id());

    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
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
            DRT::UTILS::ExtractMyValues(*ivelcolnp,myvel,gdofs);
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
            const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
            DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);

            // position of the gausspoint in physical coordinates
            // gauss_pos_xyz = funct_boundary(j)*xyze_boundary(i,j);
            const LINALG::Matrix<nsd,1> gauss_pos_xyz = XFLUID::interpolateVectorFieldToIntPoint(xyze_boundary,funct_boundary,numnode_boundary);

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
                    Velx,
                    funct,
                    enr_funct);

            // shape functions for element dofs
            enrvals.ComputeEnrichedElementShapefunction(
                    Sigmaxx,
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


            Epetra_SerialDenseVector shp_iface(numnode_boundary*begids.size());
            int pos = 0;
            for (std::set<int>::const_iterator begid = begids.begin(); begid != begids.end();++begid)
            {
              if (*begid == boundaryele->Id())
              {
                for (std::size_t inode=0; inode < numnode_boundary; ++inode)
                {
                  shp_iface(pos+inode) = funct_boundary(inode);
                }
                break;
              }
              pos += numnode_boundary;
            }

            // get normal vector (in physical coordinates) to surface element at integration point
            LINALG::Matrix<nsd,1> normalvec_solid;
            GEO::computeNormalToSurfaceElement(boundaryele->Shape(), xyze_boundary, posXiBoundary, normalvec_solid);
            LINALG::Matrix<nsd,1> normalvec_fluid(true);
            normalvec_fluid.Update(-1.0,normalvec_solid,0.0);

            // get velocities (n+g,i) at integration point
            // gpvelnp = evelnp(i,j)*shp(j);
            LINALG::Matrix<nsd,1> gpvelnp;
            gpvelnp.Clear();
            for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
                for (std::size_t isd = 0; isd < nsd; ++isd)
                    gpvelnp(isd) += evelnp(isd,iparam)*shp(iparam);

            // get interface velocity
            LINALG::Matrix<nsd,1> interface_gpvelnp;
            interface_gpvelnp.Clear();
            if (timealgo != timeint_stationary)
                for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
                    for (std::size_t isd = 0; isd < nsd; ++isd)
                        interface_gpvelnp(isd) += vel_boundary(isd,inode)*funct_boundary(inode);
            
#if 0
            // for Jeffery-Hamel Flow
            LINALG::Matrix<3,1> physpos(true);
            GEO::elementToCurrentCoordinates(DISTYPE, xyze, posXiDomain, physpos);
            
            const double x = physpos(0);
            const double y = physpos(1);
            
            const double alpha = atan(y/x);
            const double u_alpha = DRT::UTILS::JefferyHamelFlowFunction::RadialVelocity(alpha);
            
            const double nu = 1;
            interface_gpvelnp(0) = nu * (u_alpha/(x*x+y*y))*x;
            interface_gpvelnp(1) = nu * (u_alpha/(x*x+y*y))*y;
            interface_gpvelnp(2) = 0.0;
#endif
            

            // get viscous stress unknowns
            static LINALG::Matrix<nsd,nsd> tau;
            XFLUID::fill_tau(numparamtauxx, shp_tau, etau, tau);

            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

               /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */

            assembler.template Matrix<Sigmaxx,Velx>(shp_tau, -timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Sigmaxy,Velx>(shp_tau, -timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Sigmaxz,Velx>(shp_tau, -timefacfac*normalvec_fluid(2), shp);
            assembler.template Matrix<Sigmayx,Vely>(shp_tau, -timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Sigmayy,Vely>(shp_tau, -timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Sigmayz,Vely>(shp_tau, -timefacfac*normalvec_fluid(2), shp);
            assembler.template Matrix<Sigmazx,Velz>(shp_tau, -timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Sigmazy,Velz>(shp_tau, -timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Sigmazz,Velz>(shp_tau, -timefacfac*normalvec_fluid(2), shp);

            assembler.template Vector<Sigmaxx>(shp_tau, timefacfac*normalvec_fluid(0)*gpvelnp(0));
            assembler.template Vector<Sigmaxy>(shp_tau, timefacfac*normalvec_fluid(1)*gpvelnp(0));
            assembler.template Vector<Sigmaxz>(shp_tau, timefacfac*normalvec_fluid(2)*gpvelnp(0));
            assembler.template Vector<Sigmayx>(shp_tau, timefacfac*normalvec_fluid(0)*gpvelnp(1));
            assembler.template Vector<Sigmayy>(shp_tau, timefacfac*normalvec_fluid(1)*gpvelnp(1));
            assembler.template Vector<Sigmayz>(shp_tau, timefacfac*normalvec_fluid(2)*gpvelnp(1));
            assembler.template Vector<Sigmazx>(shp_tau, timefacfac*normalvec_fluid(0)*gpvelnp(2));
            assembler.template Vector<Sigmazy>(shp_tau, timefacfac*normalvec_fluid(1)*gpvelnp(2));
            assembler.template Vector<Sigmazz>(shp_tau, timefacfac*normalvec_fluid(2)*gpvelnp(2));


               /*                            \
              |  (virt tau) * n^f , u^\iface  |
               \                            */

            assembler.template Vector<Sigmaxx>(shp_tau, -timefacfac*normalvec_fluid(0)*interface_gpvelnp(0));
            assembler.template Vector<Sigmaxy>(shp_tau, -timefacfac*normalvec_fluid(1)*interface_gpvelnp(0));
            assembler.template Vector<Sigmaxz>(shp_tau, -timefacfac*normalvec_fluid(2)*interface_gpvelnp(0));
            assembler.template Vector<Sigmayx>(shp_tau, -timefacfac*normalvec_fluid(0)*interface_gpvelnp(1));
            assembler.template Vector<Sigmayy>(shp_tau, -timefacfac*normalvec_fluid(1)*interface_gpvelnp(1));
            assembler.template Vector<Sigmayz>(shp_tau, -timefacfac*normalvec_fluid(2)*interface_gpvelnp(1));
            assembler.template Vector<Sigmazx>(shp_tau, -timefacfac*normalvec_fluid(0)*interface_gpvelnp(2));
            assembler.template Vector<Sigmazy>(shp_tau, -timefacfac*normalvec_fluid(1)*interface_gpvelnp(2));
            assembler.template Vector<Sigmazz>(shp_tau, -timefacfac*normalvec_fluid(2)*interface_gpvelnp(2));

               /*               \
            - |  v , Dtau * n^f  |
               \               */

            assembler.template Matrix<Velx,Sigmaxx>(shp, -timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Velx,Sigmaxy>(shp, -timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Velx,Sigmaxz>(shp, -timefacfac*normalvec_fluid(2), shp_tau);
            assembler.template Matrix<Vely,Sigmayx>(shp, -timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Vely,Sigmayy>(shp, -timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Vely,Sigmayz>(shp, -timefacfac*normalvec_fluid(2), shp_tau);
            assembler.template Matrix<Velz,Sigmazx>(shp, -timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Velz,Sigmazy>(shp, -timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Velz,Sigmazz>(shp, -timefacfac*normalvec_fluid(2), shp_tau);

            LINALG::Matrix<nsd,1> disctau_times_nf;
            disctau_times_nf.Multiply(tau,normalvec_fluid);
            //cout << "sigmaijnj : " << disctau_times_n << endl;
            assembler.template Vector<Velx>(shp, timefacfac*disctau_times_nf(0));
            assembler.template Vector<Vely>(shp, timefacfac*disctau_times_nf(1));
            assembler.template Vector<Velz>(shp, timefacfac*disctau_times_nf(2));


            if (monolithic_FSI)
            {
                  /*                 \
               - |  v^i , Dtau * n^f  |
                  \                 */
              patchassembler.template Matrix<Dispx,Sigmaxx>(Gds, shp_iface, -timefacfac*normalvec_solid(0), shp_tau);
              patchassembler.template Matrix<Dispx,Sigmaxy>(Gds, shp_iface, -timefacfac*normalvec_solid(1), shp_tau);
              patchassembler.template Matrix<Dispx,Sigmaxz>(Gds, shp_iface, -timefacfac*normalvec_solid(2), shp_tau);
              patchassembler.template Matrix<Dispy,Sigmayx>(Gds, shp_iface, -timefacfac*normalvec_solid(0), shp_tau);
              patchassembler.template Matrix<Dispy,Sigmayy>(Gds, shp_iface, -timefacfac*normalvec_solid(1), shp_tau);
              patchassembler.template Matrix<Dispy,Sigmayz>(Gds, shp_iface, -timefacfac*normalvec_solid(2), shp_tau);
              patchassembler.template Matrix<Dispz,Sigmazx>(Gds, shp_iface, -timefacfac*normalvec_solid(0), shp_tau);
              patchassembler.template Matrix<Dispz,Sigmazy>(Gds, shp_iface, -timefacfac*normalvec_solid(1), shp_tau);
              patchassembler.template Matrix<Dispz,Sigmazz>(Gds, shp_iface, -timefacfac*normalvec_solid(2), shp_tau);

              LINALG::Matrix<nsd,1> disctau_times_ns;
              disctau_times_ns.Multiply(tau,normalvec_solid);
              //cout << "sigmaijnj : " << disctau_times_n << endl;
              patchassembler.template Vector<Dispx>(rhsd, shp_iface, timefacfac*disctau_times_ns(0));
              patchassembler.template Vector<Dispy>(rhsd, shp_iface, timefacfac*disctau_times_ns(1));
              patchassembler.template Vector<Dispz>(rhsd, shp_iface, timefacfac*disctau_times_ns(2));
            }

            // here the interface force is integrated
            // this is done using test shape functions of the boundary mesh
            // hence, we can't use the local assembler here
            for (size_t inode = 0; inode < numnode_boundary; ++inode)
            {
              force_boundary(0,inode) += funct_boundary(inode) * (disctau_times_nf(0) * timefacfac);
              force_boundary(1,inode) += funct_boundary(inode) * (disctau_times_nf(1) * timefacfac);
              force_boundary(2,inode) += funct_boundary(inode) * (disctau_times_nf(2) * timefacfac);
            }

        } // end loop over gauss points

        // here we need to assemble into the global force vector of the boundary discretization
        // note that we assemble into a overlapping vector, hence we add only, if we are a xfem row element
        // this way, we can later add all contributions together when exporting to interface row elements
        if (ifaceForceContribution)
        {
          const Epetra_Map* dofcolmap = ih->cutterdis()->DofColMap();
          std::vector<int> gdofs(3);
          for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
          {
            ih->cutterdis()->Dof(nodes[inode],0,gdofs);
            (*iforcecol)[dofcolmap->LID(gdofs[0])] += force_boundary(0,inode);
            (*iforcecol)[dofcolmap->LID(gdofs[1])] += force_boundary(1,inode);
            (*iforcecol)[dofcolmap->LID(gdofs[2])] += force_boundary(2,inode);
          }
        }


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

    if (ih->ElementIntersected(ele->Id()))
    {
      SysmatBoundary4<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, etau, iforcecol, Gds, rhsd,
          timealgo, dt, theta, assembler, ifaceForceContribution, monolithic_FSI);
    }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUID::callSysmat4(
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
                Sysmat4<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex20:
                Sysmat4<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex27:
                Sysmat4<DRT::Element::hex27,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::tet4:
                Sysmat4<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::tet10:
                Sysmat4<DRT::Element::tet10,XFEM::standard_assembly>(
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
                Sysmat4<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex20:
                Sysmat4<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::hex27:
                Sysmat4<DRT::Element::hex27,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
                break;
            case DRT::Element::tet4:
                Sysmat4<DRT::Element::tet4,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, iforcecol, estif, eforce, Gds, rhsd,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution, monolithic_FSI, L2);
            case DRT::Element::tet10:
                Sysmat4<DRT::Element::tet10,XFEM::xfem_assembly>(
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
