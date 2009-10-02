/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat4.cpp

\brief element formulation for 3D combustion element

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
#include "combust3_utils.H"
#include "combust3_local_assembler.H"
#include "combust3_interpolation.H"
#include "../drt_f3/xfluid3_utils.H"
#include "../drt_f3/fluid3_stabilization.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_xfem/spacetime_boundary.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_mat/matlist.H"
#include "../drt_xfem/enrichment.H"

#include "../drt_lib/drt_element.H"


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
  template<> struct SizeFac<XFEM::xfem_assembly>     {static const std::size_t fac = 2;};


  //! fill a number of arrays with unknown values from the unknown vector given by the discretization
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            class M1, class V1, class M2, class V2>
  void fillElementUnknownsArrays(
          const XFEM::ElementDofManager& dofman,
          const DRT::ELEMENTS::Combust3::MyState& mystate,
          M1& evelnp,
          M1& eveln,
          M1& evelnm,
          M1& eaccn,
          V1& eprenp,
          V2& ephi,
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
      // put one here to create arrays of size 1, since they are not needed anyway
      // in the xfem assembly, the numparam is determined by the dofmanager
      const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

      const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
      const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
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
      // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
      // TODO: this is inefficient, but it is nice to have only fixed size matrices afterwards!
      for (size_t iparam=0; iparam<numnode; ++iparam)
          ephi(iparam) = mystate.phinp[iparam];
  }

/*!
  build domain integral entries for combustion problem
 */
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            size_t NUMDOF,
            size_t shpVecSize,
            size_t shpVecSizeStress>
  void BuildDomainIntegralsCombust(
      LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>&           assembler,
      const XFLUID::ApproxFunc<shpVecSize>&                     shp,
      const LINALG::Matrix<shpVecSizeStress,1>&  shp_tau,
      const double& fac,
      const double& timefac,
      const double& timefacfac,
      const double& visc,
      const LINALG::Matrix<3,1>& gpvelnp,
      const double&              pres,
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
      // inertia term (contribution to mass matrix)
      /*
                           /        \
                          |          |
                          |  v , Du  |
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

  // convection term, convective part
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
      // convection term, reactive part
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

  // viscous term
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

  // pressure term
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

  // solenoidality term - continuity equation
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
  /*
                  /    \
                 |      |
                 | v, f |             is this correct? henke 09/09
                 |      |
                  \    /
  */
  assembler.template Vector<Velx>(shp.d0, fac*rhsint(0));
  assembler.template Vector<Vely>(shp.d0, fac*rhsint(1));
  assembler.template Vector<Velz>(shp.d0, fac*rhsint(2));

  // Hellinger-Reissner terms
//  if (tauele_unknowns_present)
//  {
      /*
                       /                      \
                    - |  virt tau , eps(Dtau)  |
                       \                      /
      */
//
//      const double reciproke_viscfac = 1.0/(2.0*visc);
//      assembler.template Matrix<Sigmaxx,Sigmaxx>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
//      assembler.template Matrix<Sigmaxy,Sigmaxy>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
//      assembler.template Matrix<Sigmaxz,Sigmaxz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
//      assembler.template Matrix<Sigmayy,Sigmayy>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
//      assembler.template Matrix<Sigmayz,Sigmayz>(shp_tau, -reciproke_viscfac*timefacfac*2.0, shp_tau);
//      assembler.template Matrix<Sigmazz,Sigmazz>(shp_tau, -reciproke_viscfac*timefacfac, shp_tau);
//
//      assembler.template Vector<Sigmaxx>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,0));
//      assembler.template Vector<Sigmaxy>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,1)*2.0);
//      assembler.template Vector<Sigmaxz>(shp_tau,  reciproke_viscfac*timefacfac*tau(0,2)*2.0);
//      assembler.template Vector<Sigmayy>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,1));
//      assembler.template Vector<Sigmayz>(shp_tau,  reciproke_viscfac*timefacfac*tau(1,2)*2.0);
//      assembler.template Vector<Sigmazz>(shp_tau,  reciproke_viscfac*timefacfac*tau(2,2));
//
      /*             /                  \
                    | virt tau , eps(Du) |
                     \                  /
      */
//      assembler.template Matrix<Sigmaxx,Velx>(shp_tau,     timefacfac    , shp.dx);
//      assembler.template Matrix<Sigmaxy,Velx>(shp_tau,     timefacfac    , shp.dy);
//      assembler.template Matrix<Sigmaxy,Vely>(shp_tau,     timefacfac    , shp.dx);
//      assembler.template Matrix<Sigmaxz,Velx>(shp_tau,     timefacfac    , shp.dz);
//      assembler.template Matrix<Sigmaxz,Velz>(shp_tau,     timefacfac    , shp.dx);
//      assembler.template Matrix<Sigmayy,Vely>(shp_tau,     timefacfac    , shp.dy);
//      assembler.template Matrix<Sigmayz,Vely>(shp_tau,     timefacfac    , shp.dz);
//      assembler.template Matrix<Sigmayz,Velz>(shp_tau,     timefacfac    , shp.dy);
//      assembler.template Matrix<Sigmazz,Velz>(shp_tau,     timefacfac    , shp.dz);
//
//      assembler.template Vector<Sigmaxx>(shp_tau,    - timefacfac*vderxy(0, 0));
//      assembler.template Vector<Sigmaxy>(shp_tau,    - timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
//      assembler.template Vector<Sigmaxz>(shp_tau,    - timefacfac*(vderxy(0, 2) + vderxy(2, 0)));
//      assembler.template Vector<Sigmayy>(shp_tau,    - timefacfac*vderxy(1, 1));
//      assembler.template Vector<Sigmayz>(shp_tau,    - timefacfac*(vderxy(1, 2) + vderxy(1, 2)));
//      assembler.template Vector<Sigmazz>(shp_tau,    - timefacfac*vderxy(2, 2));
//
//
//      /* pressure-pressure coupling, rectangular part */
      /*
                     /                    \
                    |                      |
                  - | tr(virt tau^e) , p I |
                    |                      |
                     \                    /
      */
//      assembler.template Matrix<Sigmaxx,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
//      assembler.template Matrix<Sigmayy,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
//      assembler.template Matrix<Sigmazz,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
//
//      assembler.template Vector<Sigmaxx>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
//      assembler.template Vector<Sigmayy>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
//      assembler.template Vector<Sigmazz>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
//
//  }

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
  Calculate integrals in matrix and rhs for stationary two-phase flow problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class M2, class V2>
void SysmatTwoPhase(
    const DRT::ELEMENTS::Combust3*      ele,           ///< the element whose matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
    const M1&                           eveln,
    const M1&                           evelnm,
    const M1&                           eaccn,
    const V1&                           eprenp,
    const V2&                           ephi,
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
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler
)
{
    TEUCHOS_FUNC_TIME_MONITOR(" - evaluate - Sysmat - domain");

    // number of nodes for element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // space dimension for 3d fluid element
    const size_t nsd = 3;

    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // get node coordinates of the current element
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    // check if we really got a list of materials
    dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());

    // flag for higher order elements
    const bool higher_order_ele = XFLUID::secondDerivativesAvailable<DISTYPE>();

    //const int numparamvelx = getNumParam<ASSTYPE>(dofman, Velx, numnode);
    // different enrichments for pressure and velocity possible
    const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
//    std::cout << "element " << ele->Id() << std::endl;
//    std::cout << "velocity dofs " << numparamvelx << std::endl;
//    std::cout << "pressure dofs " << numparampres << std::endl;

    // stabilization parameter
    const double hk = XFLUID::HK<DISTYPE>(evelnp,xyze);
    const double mk = XFLUID::MK<DISTYPE>();

    // information about domain integration cells
    const GEO::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele));
//    cout << "element "<< ele->Id() << ": ";

    //std::cout << "number of IntCells " <<  domainIntCells.size() << std::endl;
    // loop over integration cells
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {
//      std::cout << "jetzt kommt die Schleife ueber alle Integrationszellen" << std::endl;
;
    	bool influidplus = cell->getDomainPlus();
//    	  if (influidplus)
//    	   std::cout << "domain + sys " << std::endl;
//    	  else
//    	   std::cout << "domain - sys " << std::endl;
    	int matid = 3;
//    XFEM::Enrichment::ApproachFrom   approachdirection = XFEM::Enrichment::approachFromPlus;
    	if (not influidplus)
    	{
    		matid = 4;
//    		std::cout << "domain - sys " << std::endl;
    		//approachdirection = XFEM::Enrichment::approachFromMinus;
    	}

        Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
        dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "Material law is not of type m_fluid.");
        const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
        // here: DYNAMIC VISCOSITY (mu)
        const double visc = mat->Viscosity();
// TODO: kann das so bleiben, das man direkt im Dat-File die dyn Viskosität vorgibt
        // ggf const double visc = kinvisc * dens
        const double dens = mat->Density();
        // Hence, the matrices have to be multiplied by the density.
        //std::cout << "Viscosity" << visc << std::endl;


//        const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape());
// special getXFEMGaussruleKinkEnr for kink enrichment is called as parabolic shape functions are obtained
// after multipying N and Psi
        const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussruleKinkEnr<DISTYPE>(ele, xyze, ele->Intersected(),cell->Shape());

        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop over Gauss points
        for(int iquad=0; iquad<intpoints.nquad; ++iquad)
        {

//        	std::cout << "jetzt kommt Gausspunktschleife" << std::endl;
//            // coordinates of the current integration point in cell coordinates \eta^domain
//            const LINALG::Matrix<nsd,1> pos_eta_domain(intpoints.qxg[iquad]);
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
//        	   std::cout << derxy(0,no) << "    " << derxy(1,no) << "    " << derxy(2,no) << std::endl;
////        	   std::cout << derxy(1,no) << std::endl;
////        	   std::cout << derxy(2,no) << std::endl;
//            }
            // compute second global derivative
            static LINALG::Matrix<6,numnode> derxy2;
            if (higher_order_ele)
            {
                static LINALG::Matrix<6,numnode> deriv2;
                DRT::UTILS::shape_function_3D_deriv2(deriv2,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
                DRT::UTILS::gder2<DISTYPE>(xjm, derxy, deriv2, xyze, derxy2);
                //TEST
//                for (std::size_t no=0; no<numnode; no++)
//                {
//            	   std::cout << derxy2(0,no) << "   " << derxy2(1,no) << "   " << derxy2(2,no) << "   " << derxy2(3,no) << "   " << derxy2(4,no) << "   " << derxy2(5,no) << std::endl;
////            	   std::cout << derxy2(1,no) << std::endl;
////            	   std::cout << derxy2(2,no) << std::endl;
////            	   std::cout << derxy2(3,no) << std::endl;
////            	   std::cout << derxy2(5,no) << std::endl;
////            	   std::cout << derxy2(5,no) << std::endl;
//                }
            }
            else
            {
                derxy2.Clear();
            }

            // the enrichment functions depend on the gauss point
            // therefore the comutation of the enrichment functions is called here
            // the gauss point is contained in funct!
            const XFEM::ElementEnrichmentValues enrvals(
                  *ele,
                  dofman,
                  ephi,
                  funct,
                  derxy,
                  derxy2);

            const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;


            //static Shp<shpVecSize> shp;

            static COMBUST::ApproxFunc<shpVecSize> shpvel;
            static COMBUST::ApproxFunc<shpVecSize> shppres;

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
                enrvals.ComputeKinkEnrichedNodalShapefunction(
                        Velx,
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
                enrvals.ComputeKinkEnrichedNodalShapefunction(
                        Pres,
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
            const LINALG::Matrix<nsd,1> gpvelnp = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shpvel.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpveln  = COMBUST::interpolateVectorFieldToIntPoint(eveln , shpvel.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpvelnm = COMBUST::interpolateVectorFieldToIntPoint(evelnm, shpvel.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpaccn  = COMBUST::interpolateVectorFieldToIntPoint(eaccn , shpvel.d0, numparamvelx);

            // commenting this section out leads to problems with instationary calculations   henke 01/09
//            if (ASSTYPE == XFEM::xfem_assembly and timealgo != timeint_stationary)
//            {
//              const bool valid_spacetime_cell_found = COMBUST::modifyOldTimeStepsValues<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, ivelcoln, ivelcolnm, iacccoln, gpveln, gpvelnm, gpaccn);
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
            const LINALG::Matrix<nsd,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
                gpveln, gpvelnm, gpaccn, timealgo, dt, theta);

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
//            cout << bodyforce << endl;
            ///////////////LINALG::SerialDenseVector bodyforce_(enr_edeadng_(i,j)*enr_funct_(j));

            // compute stabilization parameters (3 taus)
            double tau_stab_M  = 0.0;
            double tau_stab_Mp = 0.0;
            double tau_stab_C  = 0.0;
// TODO: was wird hier berechnet ist kin visc = visc/dens ok?
            if (ASSTYPE == XFEM::xfem_assembly)
            {
                XFLUID::computeStabilization(shpvel.dx, shpvel.dy, shpvel.dz, gpvelnp, numparamvelx, instationary, visc/dens, hk, mk, timefac, tau_stab_M, tau_stab_Mp, tau_stab_C);
//                std::cout << "xfem_assembly " << std::endl;
            }
            else
            	XFLUID::computeStabilization(shpvel.dx, shpvel.dy, shpvel.dz, gpvelnp, numparamvelx, instationary, visc/dens, hk, mk, timefac, tau_stab_M, tau_stab_Mp, tau_stab_C);

            //Modificate stabilization
            tau_stab_M /= dens;
            tau_stab_Mp /= dens;


            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            LINALG::Matrix<nsd,1> rhsint;
            LINALG::Matrix<nsd,1> bodyforce;
            bodyforce.Clear();
//            std::cout << "BodyForce" << std::endl;
            bodyforce(1) = -9.8;
// --------------- DAS GEHT AUCH, WENN MAN EINE VOLUMENLAST IM DAT-FILE VORGIBT!!!!!!!!!!!!!! -------
//            LINALG::SerialDenseMatrix edeadng = XFLUID::BodyForceTwoPhaseFlow<DISTYPE>(ele, 0.0);
//            for (std::size_t isd = 0; isd < nsd; isd++)
//            {
//            	for (std::size_t inode = 0; inode < numnode; inode++)
//            		bodyforce(isd) += edeadng(isd,inode) * funct(inode);
//            }
//            std::cout << "BodyForce " << std::endl;
//            std::cout << bodyforce(0) << std::endl;
//            std::cout << bodyforce(1) << std::endl;
//            std::cout << bodyforce(2) << std::endl;
            for (size_t isd = 0; isd < nsd; ++isd)
                rhsint(isd) = dens*histvec(isd) + bodyforce(isd)*timefac*dens;

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
            for (size_t isd = 0; isd < nsd; ++isd)
                res_old(isd) = -rhsint(isd)+timefac*(dens*conv_old(isd)+gradp(isd)-2.0*visc*visc_old(isd));

            if (instationary)
            {
            	for(size_t isd=0; isd<nsd; ++isd)
            		res_old(isd) += dens * gpvelnp(isd);
            }
            //res_old += gpvelnp;


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


//            std::cout << "----Berechnung der Matrixteile----" << std::endl;
            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

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
                assembler.template Matrix<Velx,Velx>(shpvel.d0, fac*dens, shpvel.d0);
                assembler.template Matrix<Vely,Vely>(shpvel.d0, fac*dens, shpvel.d0);
                assembler.template Matrix<Velz,Velz>(shpvel.d0, fac*dens, shpvel.d0);

                assembler.template Vector<Velx>(shpvel.d0, -fac*dens*gpvelnp(0));
                assembler.template Vector<Vely>(shpvel.d0, -fac*dens*gpvelnp(1));
                assembler.template Vector<Velz>(shpvel.d0, -fac*dens*gpvelnp(2));
            }

            /* convection, convective part */
            /*
                         /                          \
                        |         / n+1       \      |
                        | v , roh| u   o nabla | Du  |
                        |         \ (i)       /      |
                         \                          /
            */
            assembler.template Matrix<Velx,Velx>(shpvel.d0, timefacfac*dens, enr_conv_c_);
            assembler.template Matrix<Vely,Vely>(shpvel.d0, timefacfac*dens, enr_conv_c_);
            assembler.template Matrix<Velz,Velz>(shpvel.d0, timefacfac*dens, enr_conv_c_);

            assembler.template Vector<Velx>(shpvel.d0, -timefacfac*dens*(gpvelnp(0)*vderxy(0,0) // check order
                                                             +gpvelnp(1)*vderxy(0,1)
                                                             +gpvelnp(2)*vderxy(0,2)));
            assembler.template Vector<Vely>(shpvel.d0, -timefacfac*dens*(gpvelnp(0)*vderxy(1,0)
                                                             +gpvelnp(1)*vderxy(1,1)
                                                             +gpvelnp(2)*vderxy(1,2)));
            assembler.template Vector<Velz>(shpvel.d0, -timefacfac*dens*(gpvelnp(0)*vderxy(2,0)
                                                             +gpvelnp(1)*vderxy(2,1)
                                                             +gpvelnp(2)*vderxy(2,2)));

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
                assembler.template Matrix<Velx,Velx>(shpvel.d0, timefacfac*dens*vderxy(0,0), shpvel.d0);
                assembler.template Matrix<Velx,Vely>(shpvel.d0, timefacfac*dens*vderxy(0,1), shpvel.d0);
                assembler.template Matrix<Velx,Velz>(shpvel.d0, timefacfac*dens*vderxy(0,2), shpvel.d0);
                assembler.template Matrix<Vely,Velx>(shpvel.d0, timefacfac*dens*vderxy(1,0), shpvel.d0);
                assembler.template Matrix<Vely,Vely>(shpvel.d0, timefacfac*dens*vderxy(1,1), shpvel.d0);
                assembler.template Matrix<Vely,Velz>(shpvel.d0, timefacfac*dens*vderxy(1,2), shpvel.d0);
                assembler.template Matrix<Velz,Velx>(shpvel.d0, timefacfac*dens*vderxy(2,0), shpvel.d0);
                assembler.template Matrix<Velz,Vely>(shpvel.d0, timefacfac*dens*vderxy(2,1), shpvel.d0);
                assembler.template Matrix<Velz,Velz>(shpvel.d0, timefacfac*dens*vderxy(2,2), shpvel.d0);
            }

            /* Viskositaetsterm */
            /*
                          /                        \
                         |       / \         /  \   |
                         |  eps | v | , tau | Du |  |
                         |       \ /         \  /   |
                          \                        /
            */
            // visc is dynamic viscosity
            assembler.template Matrix<Velx,Velx>(shpvel.dx,   2.0*visc*timefacfac, shpvel.dx);
            assembler.template Matrix<Velx,Velx>(shpvel.dy,     visc*timefacfac, shpvel.dy);
            assembler.template Matrix<Velx,Vely>(shpvel.dy,     visc*timefacfac, shpvel.dx);
            assembler.template Matrix<Velx,Velx>(shpvel.dz,     visc*timefacfac, shpvel.dz);
            assembler.template Matrix<Velx,Velz>(shpvel.dz,     visc*timefacfac, shpvel.dx);

            assembler.template Matrix<Vely,Vely>(shpvel.dx,     visc*timefacfac, shpvel.dx);
            assembler.template Matrix<Vely,Velx>(shpvel.dx,     visc*timefacfac, shpvel.dy);
            assembler.template Matrix<Vely,Vely>(shpvel.dy,   2.0*visc*timefacfac, shpvel.dy);
            assembler.template Matrix<Vely,Vely>(shpvel.dz,     visc*timefacfac, shpvel.dz);
            assembler.template Matrix<Vely,Velz>(shpvel.dz,     visc*timefacfac, shpvel.dy);

            assembler.template Matrix<Velz,Velz>(shpvel.dx,     visc*timefacfac, shpvel.dx);
            assembler.template Matrix<Velz,Velx>(shpvel.dx,     visc*timefacfac, shpvel.dz);
            assembler.template Matrix<Velz,Velz>(shpvel.dy,     visc*timefacfac, shpvel.dy);
            assembler.template Matrix<Velz,Vely>(shpvel.dy,     visc*timefacfac, shpvel.dz);
            assembler.template Matrix<Velz,Velz>(shpvel.dz,   2.0*visc*timefacfac, shpvel.dz);

            assembler.template Vector<Velx>(shpvel.dx,     -visc*timefacfac*(vderxy(0, 0) + vderxy(0, 0)));
            assembler.template Vector<Velx>(shpvel.dy,     -visc*timefacfac*(vderxy(0, 1) + vderxy(1, 0)));
            assembler.template Vector<Velx>(shpvel.dz,     -visc*timefacfac*(vderxy(0, 2) + vderxy(2, 0)));

            assembler.template Vector<Vely>(shpvel.dx,     -visc*timefacfac*(vderxy(1, 0) + vderxy(0, 1)));
            assembler.template Vector<Vely>(shpvel.dy,     -visc*timefacfac*(vderxy(1, 1) + vderxy(1, 1)));
            assembler.template Vector<Vely>(shpvel.dz,     -visc*timefacfac*(vderxy(1, 2) + vderxy(2, 1)));

            assembler.template Vector<Velz>(shpvel.dx,     -visc*timefacfac*(vderxy(2, 0) + vderxy(0, 2)));
            assembler.template Vector<Velz>(shpvel.dy,     -visc*timefacfac*(vderxy(2, 1) + vderxy(1, 2)));
            assembler.template Vector<Velz>(shpvel.dz,     -visc*timefacfac*(vderxy(2, 2) + vderxy(2, 2)));

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
                	//fac*dens
                    assembler.template Matrix<Pres,Velx>(shppres.dx, timetauMp*dens, shpvel.d0);
                    assembler.template Matrix<Pres,Vely>(shppres.dy, timetauMp*dens, shpvel.d0);
                    assembler.template Matrix<Pres,Velz>(shppres.dz, timetauMp*dens, shpvel.d0);
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
                assembler.template Matrix<Pres,Velx>(shppres.dx, ttimetauMp*dens, enr_conv_c_);
                assembler.template Matrix<Pres,Vely>(shppres.dy, ttimetauMp*dens, enr_conv_c_);
                assembler.template Matrix<Pres,Velz>(shppres.dz, ttimetauMp*dens, enr_conv_c_);

                if (newton)
                {
                    /*  pressure stabilisation: convection, reactive part
                          /                                \
                         |              /          \   n+1  |
                         | grad q , roh| Du o nabla | u     |
                         |              \          /   (i)  |
                          \                                /
                    */
                    assembler.template Matrix<Pres,Velx>(shppres.dx, ttimetauMp*dens*vderxy(0,0), shpvel.d0);
                    assembler.template Matrix<Pres,Velx>(shppres.dy, ttimetauMp*dens*vderxy(1,0), shpvel.d0);
                    assembler.template Matrix<Pres,Velx>(shppres.dz, ttimetauMp*dens*vderxy(2,0), shpvel.d0);

                    assembler.template Matrix<Pres,Vely>(shppres.dx, ttimetauMp*dens*vderxy(0,1), shpvel.d0);
                    assembler.template Matrix<Pres,Vely>(shppres.dy, ttimetauMp*dens*vderxy(1,1), shpvel.d0);
                    assembler.template Matrix<Pres,Vely>(shppres.dz, ttimetauMp*dens*vderxy(2,1), shpvel.d0);

                    assembler.template Matrix<Pres,Velz>(shppres.dx, ttimetauMp*dens*vderxy(0,2), shpvel.d0);
                    assembler.template Matrix<Pres,Velz>(shppres.dy, ttimetauMp*dens*vderxy(1,2), shpvel.d0);
                    assembler.template Matrix<Pres,Velz>(shppres.dz, ttimetauMp*dens*vderxy(2,2), shpvel.d0);
                }

                /* pressure stabilisation: viscosity (-L_visc_u) */
                /*
                           /                             \
                          |                         /  \  |
                        - |  nabla q , nabla o tau | Du | |
                          |                         \  /  |
                           \                             /
                */
                assembler.template Matrix<Pres,Velx>(shppres.dx, -2.0*visc*ttimetauMp, enr_viscs2.xx);
                assembler.template Matrix<Pres,Vely>(shppres.dx, -2.0*visc*ttimetauMp, enr_viscs2.xy);
                assembler.template Matrix<Pres,Velz>(shppres.dx, -2.0*visc*ttimetauMp, enr_viscs2.xz);

                assembler.template Matrix<Pres,Velx>(shppres.dy, -2.0*visc*ttimetauMp, enr_viscs2.xy);
                assembler.template Matrix<Pres,Vely>(shppres.dy, -2.0*visc*ttimetauMp, enr_viscs2.yy);
                assembler.template Matrix<Pres,Velz>(shppres.dy, -2.0*visc*ttimetauMp, enr_viscs2.yz);

                assembler.template Matrix<Pres,Velx>(shppres.dz, -2.0*visc*ttimetauMp, enr_viscs2.xz);
                assembler.template Matrix<Pres,Vely>(shppres.dz, -2.0*visc*ttimetauMp, enr_viscs2.yz);
                assembler.template Matrix<Pres,Velz>(shppres.dz, -2.0*visc*ttimetauMp, enr_viscs2.zz);

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
                              /                               \
                             |              / n+1       \     |
                             | roh Du , roh| u   o nabla | v  |
                             |              \ (i)       /     |
                              \                               /
                    */
                	// timetauM * dens * dens, da einmal Dichte von Wichtungsfunktion und einmal von Testfunktion
                	// vgl. Paper bzw FEFluid-Skript bzw Blatt
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, timetauM*dens*dens, shpvel.d0);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, timetauM*dens*dens, shpvel.d0);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, timetauM*dens*dens, shpvel.d0);

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
                        assembler.template Matrix<Velx,Velx>(shpvel.dx, timetauM*dens*dens*gpvelnp(0), shpvel.d0);
                        assembler.template Matrix<Velx,Vely>(shpvel.dy, timetauM*dens*dens*gpvelnp(0), shpvel.d0);
                        assembler.template Matrix<Velx,Velz>(shpvel.dz, timetauM*dens*dens*gpvelnp(0), shpvel.d0);
                        assembler.template Matrix<Vely,Velx>(shpvel.dx, timetauM*dens*dens*gpvelnp(1), shpvel.d0);
                        assembler.template Matrix<Vely,Vely>(shpvel.dy, timetauM*dens*dens*gpvelnp(1), shpvel.d0);
                        assembler.template Matrix<Vely,Velz>(shpvel.dz, timetauM*dens*dens*gpvelnp(1), shpvel.d0);
                        assembler.template Matrix<Velz,Velx>(shpvel.dx, timetauM*dens*dens*gpvelnp(2), shpvel.d0);
                        assembler.template Matrix<Velz,Vely>(shpvel.dy, timetauM*dens*dens*gpvelnp(2), shpvel.d0);
                        assembler.template Matrix<Velz,Velz>(shpvel.dz, timetauM*dens*dens*gpvelnp(2), shpvel.d0);
                    }
                }

                const double ttimetauM  = timefac * timefac * tau_stab_M * fac;
                /* supg stabilisation: convective part ( L_conv_u) */
                /*
                     /                                                \
                    |     / n+1        \           / n+1        \      |
                    | roh| u    o nabla | v , roh | u    o nabla | Du  |
                    |     \ (i)        /           \ (i)        /      |
                     \                                                /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM*dens*dens, enr_conv_c_);
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM*dens*dens, enr_conv_c_);
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM*dens*dens, enr_conv_c_);
                /* supg stabilisation: pressure part  ( L_pres_p) */
                /*
                          /                                \
                         |      / n+1       \               |
                         | roh | u   o nabla | v , nabla Dp |
                         |      \ (i)       /               |
                          \                                /
                */
                assembler.template Matrix<Velx,Pres>(enr_conv_c_, ttimetauM*dens, shppres.dx);
                assembler.template Matrix<Vely,Pres>(enr_conv_c_, ttimetauM*dens, shppres.dy);
                assembler.template Matrix<Velz,Pres>(enr_conv_c_, ttimetauM*dens, shppres.dz);

                /* supg stabilisation: viscous part  (-L_visc_u) */
                /*
                      /                                           \
                     |               /  \       / n+1        \     |
                   - |  nabla o eps | Du |, roh| u    o nabla | v  |
                     |               \  /       \ (i)        /     |
                      \                                           /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.xx);
                assembler.template Matrix<Velx,Vely>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.xy);
                assembler.template Matrix<Velx,Velz>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.xz);

                assembler.template Matrix<Vely,Velx>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.yx);
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.yy);
                assembler.template Matrix<Vely,Velz>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.yz);

                assembler.template Matrix<Velz,Velx>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.zx);
                assembler.template Matrix<Velz,Vely>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.zy);
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, -2.0*visc*ttimetauM*dens, enr_viscs2.zz);

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
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, ttimetauM*dens*dens*vderxy(0,0), shpvel.d0);
                    assembler.template Matrix<Velx,Vely>(enr_conv_c_, ttimetauM*dens*dens*vderxy(0,1), shpvel.d0);
                    assembler.template Matrix<Velx,Velz>(enr_conv_c_, ttimetauM*dens*dens*vderxy(0,2), shpvel.d0);

                    assembler.template Matrix<Vely,Velx>(enr_conv_c_, ttimetauM*dens*dens*vderxy(1,0), shpvel.d0);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, ttimetauM*dens*dens*vderxy(1,1), shpvel.d0);
                    assembler.template Matrix<Vely,Velz>(enr_conv_c_, ttimetauM*dens*dens*vderxy(1,2), shpvel.d0);

                    assembler.template Matrix<Velz,Velx>(enr_conv_c_, ttimetauM*dens*dens*vderxy(2,0), shpvel.d0);
                    assembler.template Matrix<Velz,Vely>(enr_conv_c_, ttimetauM*dens*dens*vderxy(2,1), shpvel.d0);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, ttimetauM*dens*dens*vderxy(2,2), shpvel.d0);

                    /*
                             /                                                  |
                            |       / n+1         |   n+1      /           |    |
                            |  roh | u    o nabla | u    ,roh | Du o nabla | v  |
                            |      | (i)          |   (i)     |            |    |
                             \                                                  |
                    */
                    const double con0 = ttimetauM*(gpvelnp(0)*vderxy(0,0) + gpvelnp(1)*vderxy(0,1) + gpvelnp(2)*vderxy(0,2));
                    assembler.template Matrix<Velx,Velx>(shpvel.dx, con0*dens*dens, shpvel.d0);
                    assembler.template Matrix<Velx,Vely>(shpvel.dy, con0*dens*dens, shpvel.d0);
                    assembler.template Matrix<Velx,Velz>(shpvel.dz, con0*dens*dens, shpvel.d0);

                    const double con1 = ttimetauM*(gpvelnp(0)*vderxy(1,0) + gpvelnp(1)*vderxy(1,1) + gpvelnp(2)*vderxy(1,2));
                    assembler.template Matrix<Vely,Velx>(shpvel.dx, con1*dens*dens, shpvel.d0);
                    assembler.template Matrix<Vely,Vely>(shpvel.dy, con1*dens*dens, shpvel.d0);
                    assembler.template Matrix<Vely,Velz>(shpvel.dz, con1*dens*dens, shpvel.d0);

                    const double con2 = ttimetauM*(gpvelnp(0)*vderxy(2,0) + gpvelnp(1)*vderxy(2,1) + gpvelnp(2)*vderxy(2,2));
                    assembler.template Matrix<Velz,Velx>(shpvel.dx, con2*dens*dens, shpvel.d0);
                    assembler.template Matrix<Velz,Vely>(shpvel.dy, con2*dens*dens, shpvel.d0);
                    assembler.template Matrix<Velz,Velz>(shpvel.dz, con2*dens*dens, shpvel.d0);

                    /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
                    /*
                                    /                                  \
                                   |         n+1       /          \     |
                                   |  nabla p    ,roh | Du o nabla | v  |
                                   |         (i)       \          /     |
                                    \                                  /
                    */
                    assembler.template Matrix<Velx,Velx>(shpvel.dx, ttimetauM*dens*gradp(0), shpvel.d0);
                    assembler.template Matrix<Velx,Vely>(shpvel.dy, ttimetauM*dens*gradp(0), shpvel.d0);
                    assembler.template Matrix<Velx,Velz>(shpvel.dz, ttimetauM*dens*gradp(0), shpvel.d0);

                    assembler.template Matrix<Vely,Velx>(shpvel.dx, ttimetauM*dens*gradp(1), shpvel.d0);
                    assembler.template Matrix<Vely,Vely>(shpvel.dy, ttimetauM*dens*gradp(1), shpvel.d0);
                    assembler.template Matrix<Vely,Velz>(shpvel.dz, ttimetauM*dens*gradp(1), shpvel.d0);

                    assembler.template Matrix<Velz,Velx>(shpvel.dx, ttimetauM*dens*gradp(2), shpvel.d0);
                    assembler.template Matrix<Velz,Vely>(shpvel.dy, ttimetauM*dens*gradp(2), shpvel.d0);
                    assembler.template Matrix<Velz,Velz>(shpvel.dz, ttimetauM*dens*gradp(2), shpvel.d0);

                      /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                      /*
                              /                                            \
                             |               / n+1 \       /          \     |
                           - |  nabla o eps | u     |,roh | Du o nabla | v  |
                             |               \ (i) /       \          /     |
                              \                                            /
                      */
                    assembler.template Matrix<Velx,Velx>(shpvel.dx, -2.0*visc*ttimetauM*dens*visc_old(0), shpvel.d0);
                    assembler.template Matrix<Velx,Vely>(shpvel.dy, -2.0*visc*ttimetauM*dens*visc_old(0), shpvel.d0);
                    assembler.template Matrix<Velx,Velz>(shpvel.dz, -2.0*visc*ttimetauM*dens*visc_old(0), shpvel.d0);

                    assembler.template Matrix<Vely,Velx>(shpvel.dx, -2.0*visc*ttimetauM*dens*visc_old(1), shpvel.d0);
                    assembler.template Matrix<Vely,Vely>(shpvel.dy, -2.0*visc*ttimetauM*dens*visc_old(1), shpvel.d0);
                    assembler.template Matrix<Vely,Velz>(shpvel.dz, -2.0*visc*ttimetauM*dens*visc_old(1), shpvel.d0);

                    assembler.template Matrix<Velz,Velx>(shpvel.dx, -2.0*visc*ttimetauM*dens*visc_old(2), shpvel.d0);
                    assembler.template Matrix<Velz,Vely>(shpvel.dy, -2.0*visc*ttimetauM*dens*visc_old(2), shpvel.d0);
                    assembler.template Matrix<Velz,Velz>(shpvel.dz, -2.0*visc*ttimetauM*dens*visc_old(2), shpvel.d0);

                    /* supg stabilisation: bodyforce part, linearisation of test function */

                    /*
                                  /                                \
                                 |                 /          \     |
                               - |  rhsint   , roh| Du o nabla | v  |
                                 |                 \          /     |
                                  \                                /

                    */

                    assembler.template Matrix<Velx,Velx>(shpvel.dx, -timetauM*dens*rhsint(0), shpvel.d0);
                    assembler.template Matrix<Velx,Vely>(shpvel.dy, -timetauM*dens*rhsint(0), shpvel.d0);
                    assembler.template Matrix<Velx,Velz>(shpvel.dz, -timetauM*dens*rhsint(0), shpvel.d0);

                    assembler.template Matrix<Vely,Velx>(shpvel.dx, -timetauM*dens*rhsint(1), shpvel.d0);
                    assembler.template Matrix<Vely,Vely>(shpvel.dy, -timetauM*dens*rhsint(1), shpvel.d0);
                    assembler.template Matrix<Vely,Velz>(shpvel.dz, -timetauM*dens*rhsint(1), shpvel.d0);

                    assembler.template Matrix<Velz,Velx>(shpvel.dx, -timetauM*dens*rhsint(2), shpvel.d0);
                    assembler.template Matrix<Velz,Vely>(shpvel.dy, -timetauM*dens*rhsint(2), shpvel.d0);
                    assembler.template Matrix<Velz,Velz>(shpvel.dz, -timetauM*dens*rhsint(2), shpvel.d0);
                } // if newton

                // supg stabilisation
                assembler.template Vector<Velx>(enr_conv_c_, -timetauM*dens*res_old(0));
                assembler.template Vector<Vely>(enr_conv_c_, -timetauM*dens*res_old(1));
                assembler.template Vector<Velz>(enr_conv_c_, -timetauM*dens*res_old(2));
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
            } // endif cstab

        } // end loop over gauss points
    } // end loop over integration cells

    return;
}


/*!
  Calculate domain integrals in matrix and rhs for stationary combustion problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class M2, class V2>
void SysmatCombustDomain(
    const DRT::ELEMENTS::Combust3*      ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
    const M1&                           eveln,
    const M1&                           evelnm,
    const M1&                           eaccn,
    const V1&                           eprenp,
    const V2&                           ephi,
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
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler
)
{
    TEUCHOS_FUNC_TIME_MONITOR(" - evaluate - Sysmat - domain");

    // number space dimensions for 3d combustion element
    const size_t nsd = 3;

    // number of nodes of this element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // get node coordinates of the current element
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

#ifdef DEBUG
    // check if we really got a list of materials
    dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif

    // get material list for this element
    const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());

    // flag for higher order elements
    const bool higher_order_ele = XFLUID::secondDerivativesAvailable<DISTYPE>();

    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
#if 0
    // figure out whether we have stress unknowns at all
    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Sigmaxx, 0) > 0);
//    const bool velocity_unknowns_present = (getNumParam<ASSTYPE>(dofman, Velx, 1) > 0);
//    const bool pressure_unknowns_present = (getNumParam<ASSTYPE>(dofman, Pres, 1) > 0);
#endif
    // TODO remove as soon as precompiler flag above is turned off
    const bool tauele_unknowns_present = false;

    // get number of parameters (dofs) for each field
    // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
    const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
#if 0
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, numparam is determined by the dofmanager
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
#endif

    // stabilization parameter
    const double hk = XFLUID::HK<DISTYPE>(evelnp,xyze);
    const double mk = XFLUID::MK<DISTYPE>();

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
        // TODO: hack to compile - only one material is used! henke 06/09
        int matid = 3;
        // here I need a check on which side of the interface the cell is
        /*
        bool influidplus = cell->getDomainPlus();
        if(cell->indomainplus_) // cell is burnt
        {
          matid = 3; // plus == burnt material? check this!
        }
        else
        {
          matid = 4; // minus == unburnt material? check this!
        }*/
        // get material from list of materials
        Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
        // check if we really have a fluid material
        dsassert(matptr->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
        const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
        // get the kinematic viscosity \nu
//        const double kinvisc = mat->Viscosity();
        // get the density \rho
        const double dens = mat->Density();
        // get dynamic viscosity \mu
//        const double visc = kinvisc * dens;
        // TODO preliminary
        const double visc = mat->Viscosity();

        //------------------------------------------------------------------------------------------
        // evaluate the enrichment function for this integration cell
        //------------------------------------------------------------------------------------------
        const XFEM::ElementEnrichmentValues enrvals(*ele,dofman,*cell,ephi);

        //------------------------------------------------------------------------------------------
        // get Gaussian points for this integration cell
        //------------------------------------------------------------------------------------------
        const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Intersected(), cell->Shape(),false);
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

#if 0
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
#endif

            //--------------------------------------------------------------------------------------
            // procedures involving Jacobian matrix
            //--------------------------------------------------------------------------------------
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
            const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
            const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

            // (enriched) shape functions = approximation functions (P = N * \Psi)
            static XFLUID::ApproxFunc<shpVecSize> shp;
            // comment missing
            static LINALG::Matrix<shpVecSizeStress,1> shp_tau;
            if (ASSTYPE == XFEM::xfem_assembly)
            {
                // temporary arrays holding enriched shape functions (N * \Psi)
                static LINALG::Matrix<shpVecSize,1> enr_funct;
                static LINALG::Matrix<3,shpVecSize> enr_derxy;
//                static LINALG::Matrix<6,shpVecSize> enr_derxy2;

                // shape functions and derivatives for nodal parameters (dofs)
                enrvals.ComputeEnrichedNodalShapefunction(
                        Velx,
                        funct,
                        derxy,
                        enr_funct,
                        enr_derxy);

                // fill approximation functions for XFEM
                for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
                {
                  shp.d0(iparam) = enr_funct(iparam);
                  shp.dx(iparam) = enr_derxy(0,iparam);
                  shp.dy(iparam) = enr_derxy(1,iparam);
                  shp.dz(iparam) = enr_derxy(2,iparam);
//                  shp.dxdx(iparam) = enr_derxy2(0,iparam);
//                  shp.dxdy(iparam) = enr_derxy2(3,iparam);
//                  shp.dxdz(iparam) = enr_derxy2(4,iparam);
//                  shp.dydx(iparam) = shp.dxdy(iparam);
//                  shp.dydy(iparam) = enr_derxy2(1,iparam);
//                  shp.dydz(iparam) = enr_derxy2(5,iparam);
//                  shp.dzdx(iparam) = shp.dxdz(iparam);
//                  shp.dzdy(iparam) = shp.dydz(iparam);
//                  shp.dzdz(iparam) = enr_derxy2(2,iparam);
                }

#if 0
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
#endif
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
//                shp.dxdx(iparam) = derxy2(0,iparam);
//                shp.dxdy(iparam) = derxy2(3,iparam);
//                shp.dxdz(iparam) = derxy2(4,iparam);
//                shp.dydx(iparam) = shp.dxdy(iparam);
//                shp.dydy(iparam) = derxy2(1,iparam);
//                shp.dydz(iparam) = derxy2(5,iparam);
//                shp.dzdx(iparam) = shp.dxdz(iparam);
//                shp.dzdy(iparam) = shp.dydz(iparam);
//                shp.dzdz(iparam) = derxy2(2,iparam);
              }
#if 0
              if (tauele_unknowns_present)
              {
                dserror("no stress enrichments without xfem assembly");
              }
#endif
            }

            // get velocities and accelerations at integration point
            const LINALG::Matrix<nsd,1> gpvelnp = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpveln  = COMBUST::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpvelnm = COMBUST::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
            LINALG::Matrix<nsd,1> gpaccn  = COMBUST::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);

            // commenting this section out leads to problems with instationary calculations   henke 01/09
//            if (ASSTYPE == XFEM::xfem_assembly and timealgo != timeint_stationary)
//            {
//              const bool valid_spacetime_cell_found = COMBUST::modifyOldTimeStepsValues<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, ivelcoln, ivelcolnm, iacccoln, gpveln, gpvelnm, gpaccn);
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
            const LINALG::Matrix<nsd,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
                gpveln, gpvelnm, gpaccn, timealgo, dt, theta);

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
//              for (size_t iparam = 0; iparam < numparamvelx; ++iparam)
//              {
//                for (size_t isd = 0; isd < nsd; ++isd)
//                {
//                 vderxy2(isd,0) += evelnp(isd,iparam)*shp.dxdx(iparam);
//                  vderxy2(isd,1) += evelnp(isd,iparam)*shp.dydy(iparam);
//                  vderxy2(isd,2) += evelnp(isd,iparam)*shp.dzdz(iparam);
//                  vderxy2(isd,3) += evelnp(isd,iparam)*shp.dxdy(iparam);
//                  vderxy2(isd,4) += evelnp(isd,iparam)*shp.dxdz(iparam);
//                  vderxy2(isd,5) += evelnp(isd,iparam)*shp.dydz(iparam);
//               }
//              }
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
#if 0
            if (tauele_unknowns_present)
            {
              XFLUID::fill_tau(numparamtauxx, shp_tau, etau, tau);
            }
            else
#endif
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
            XFLUID::computeStabilization(shp.dx, shp.dy, shp.dz, gpvelnp, numparamvelx, instationary, visc/dens, hk, mk, timefac, tau_stab_M, tau_stab_Mp, tau_stab_C);

            //Modificate stabilization
            // TODO: is this correct?
            tau_stab_M /= dens;
            tau_stab_Mp /= dens;

            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            LINALG::Matrix<nsd,1> rhsint;
            LINALG::Matrix<nsd,1> bodyforce;
            bodyforce.Clear();
//            std::cout << "BodyForce" << std::endl;
            //bodyforce(0) = 1.0;
            //bodyforce(1) = -9.8;

// --------------- DAS GEHT AUCH, WENN MAN EINE VOLUMENLAST IM DAT-FILE VORGIBT!!!!!!!!!!!!!! -------
//            LINALG::SerialDenseMatrix edeadng = XFLUID::BodyForceTwoPhaseFlow<DISTYPE>(ele, 0.0);
//            for (std::size_t isd = 0; isd < nsd; isd++)
//            {
//            	for (std::size_t inode = 0; inode < numnode; inode++)
//            		bodyforce(isd) += edeadng(isd,inode) * funct(inode);
//            }
//            std::cout << "BodyForce " << std::endl;
//            std::cout << bodyforce(0) << std::endl;
//            std::cout << bodyforce(1) << std::endl;
//            std::cout << bodyforce(2) << std::endl;

            // TODO: modified by density
            for (size_t isd = 0; isd < nsd; ++isd)
                rhsint(isd) = histvec(isd)*dens + bodyforce(isd)*timefac*dens;

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
            // TODO: modified by density
            for (size_t isd = 0; isd < nsd; ++isd)
                res_old(isd) = -rhsint(isd)+timefac*(dens*conv_old(isd)+gradp(isd)-2.0*visc*visc_old(isd));

            // TODO: modified by density
            if (instationary)
            {
              for(size_t isd=0; isd<nsd; ++isd)
                res_old(isd) += dens * gpvelnp(isd);
//              res_old += gpvelnp;
	    }
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



            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

            BuildDomainIntegralsCombust<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
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
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2>
void SysmatCombustBoundary(
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
  TEUCHOS_FUNC_TIME_MONITOR(" - evaluate - Sysmat - boundary");

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
  // compute dynamic viscosity \mu
  const double visc_minus = kinvisc_minus * dens_minus;

  //------------------------------------------------------------------------------------------------
  // compute jump values
  //------------------------------------------------------------------------------------------------
  // velocity jump value
  const double ju = (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);
  // pressure jump value
  const double jp = -flamespeed*dens_minus*(1/dens_plus - 1/dens_minus);

  // Nitsche parameter velocity
  const double alphau = nitschevel;
  // Nitsche parameter pressure
  const double alphap = nitschepres;

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
    const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);

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
      static LINALG::Matrix<nsd,numnode> deriv;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      //--------------------------------------------------------------------------------------------
      // procedures involving Jacobian matrix for domain mapping
      //--------------------------------------------------------------------------------------------
      // get transposed of the jacobian matrix d x / d \xi
      // xjm(i,j) = deriv(i,k)*xyze(j,k)
      static LINALG::Matrix<nsd,nsd> xjm;
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
      xji.Invert(xjm);

      //--------------------------------------------------------------------------------------------
      // compute global derivates of shape functions at this Gaussian point
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<3,numnode> derxy;
      // derxy(i,j) = xji(i,k) * deriv(k,j)
      derxy.Multiply(xji,deriv);

      //--------------------------------------------------------------------------------------------
      // compute Jacobian matrix for mapping to boundary
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,numvertices> deriv_boundary;
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
      // compute total (time and spatial) integration factor (and coefficients of single terms)?
      const double timefacfac = timefac * fac;

      //--------------------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------------------
      const std::size_t shpVecSize = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
      static LINALG::Matrix<shpVecSize,1>       enrfunct_plus;
      static LINALG::Matrix<shpVecSize,1>       enrfunct_minus;

      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeEnrichedNodalShapefunction(Velx, funct, enrfunct_plus);
      enrvals_minus.ComputeEnrichedNodalShapefunction(Velx, funct, enrfunct_minus);

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFLUID::ApproxFunc<shpVecSize> shp_jump;
      static XFLUID::ApproxFunc<shpVecSize> shp_mean;
      static XFLUID::ApproxFunc<shpVecSize> shp_mean_visc;

      // fill approximation functions
      for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        shp_jump.d0(iparam) = enrfunct_plus(iparam) - enrfunct_minus(iparam);
        shp_mean.d0(iparam) = 0.5*(enrfunct_plus(iparam) + enrfunct_minus(iparam));
        shp_mean_visc.d0(iparam) = 0.5*(visc_plus*enrfunct_plus(iparam) + visc_minus*enrfunct_minus(iparam));
      }
      //--------------------------------------------------------------------------------------------
      // compute data at Gaussian point for rhs
      //--------------------------------------------------------------------------------------------

      // velocity jump
      static LINALG::Matrix<nsd,1> vjump(true);
      vjump.Clear();
      vjump = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);

      // mean velocity
      static LINALG::Matrix<nsd,1> vmean(true);
      vmean.Clear();
      vmean = COMBUST::interpolateVectorFieldToIntPoint(evelnp, shp_mean.d0, numparamvelx);

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

      // get pressure jump
      static double pjump;
      pjump = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pjump += shp_jump.d0(iparam)*eprenp(iparam);

      static double pmean;
      pmean = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pmean += shp_mean.d0(iparam)*eprenp(iparam);

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
//            if (timealgo != timeint_stationary)
//                for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
//                    for (std::size_t isd = 0; isd < nsd; ++isd)
//                        interface_gpvelnp(isd) += vel_boundary(isd,inode)*funct_boundary(inode);
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

      //---------------------------------    |                                     |
      // viscous adjoint consistency term  + |  < 2\mu epsilon( v ) > n, || Du ||  |
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
      // - |  < 2\mu epsilon( v ) > n, || u_i ||  |
      //   |                                      |

      assembler.template Vector<Velx>(shp_mean_visc.dx, -2.0* normal(0)           *timefacfac*vjump(0,0));
      assembler.template Vector<Velx>(shp_mean_visc.dy, -    (normal(0)+normal(1))*timefacfac*vjump(0,0));
      assembler.template Vector<Velx>(shp_mean_visc.dz, -    (normal(0)+normal(2))*timefacfac*vjump(0,0));

      assembler.template Vector<Vely>(shp_mean_visc.dx, -    (normal(1)+normal(0))*timefacfac*vjump(1,0));
      assembler.template Vector<Vely>(shp_mean_visc.dy, -2.0* normal(1)           *timefacfac*vjump(1,0));
      assembler.template Vector<Vely>(shp_mean_visc.dz, -    (normal(1)+normal(2))*timefacfac*vjump(1,0));

      assembler.template Vector<Velz>(shp_mean_visc.dx, -    (normal(2)+normal(0))*timefacfac*vjump(2,0));
      assembler.template Vector<Velz>(shp_mean_visc.dy, -    (normal(2)+normal(1))*timefacfac*vjump(2,0));
      assembler.template Vector<Velz>(shp_mean_visc.dz, -2.0* normal(2)           *timefacfac*vjump(2,0));

      //-------------------------------------    |                                    |
      // viscous adjoint consistency term RHS  + |  < 2\mu\rho epsilon( v ) > n, J_u  |
      //-------------------------------------    |                                    |

      assembler.template Vector<Velx>(shp_mean_visc.dx, -2.0* normal(0)           *timefacfac*normal(0)*ju);
      assembler.template Vector<Velx>(shp_mean_visc.dy, -    (normal(0)+normal(1))*timefacfac*normal(0)*ju);
      assembler.template Vector<Velx>(shp_mean_visc.dz, -    (normal(0)+normal(2))*timefacfac*normal(0)*ju);

      assembler.template Vector<Vely>(shp_mean_visc.dx, -    (normal(1)+normal(0))*timefacfac*normal(1)*ju);
      assembler.template Vector<Vely>(shp_mean_visc.dy, -2.0* normal(1)           *timefacfac*normal(1)*ju);
      assembler.template Vector<Vely>(shp_mean_visc.dz, -    (normal(1)+normal(2))*timefacfac*normal(1)*ju);

      assembler.template Vector<Velz>(shp_mean_visc.dx, -    (normal(2)+normal(0))*timefacfac*normal(2)*ju);
      assembler.template Vector<Velz>(shp_mean_visc.dy, -    (normal(2)+normal(1))*timefacfac*normal(2)*ju);
      assembler.template Vector<Velz>(shp_mean_visc.dz, -2.0* normal(2)           *timefacfac*normal(2)*ju);

      //--------------------------    |                     |
      // pressure consistency term  - |  || v ||, < Dp > n  |
      //--------------------------    |                     |

      assembler.template Matrix<Velx,Pres>(shp_jump.d0, -timefacfac*normal(0), shp_mean.d0);
      assembler.template Matrix<Vely,Pres>(shp_jump.d0, -timefacfac*normal(1), shp_mean.d0);
      assembler.template Matrix<Velz,Pres>(shp_jump.d0, -timefacfac*normal(2), shp_mean.d0);

      //   |                      |
      // + |  || v ||, < p_i > n  |
      //   |                      |

      assembler.template Vector<Velx>(shp_jump.d0, timefacfac*normal(0)*pmean);
      assembler.template Vector<Vely>(shp_jump.d0, timefacfac*normal(1)*pmean);
      assembler.template Vector<Velz>(shp_jump.d0, timefacfac*normal(2)*pmean);

      //----------------------------------    |                     |
      // pressure adjoint consistency term  + |  < q > n, || Du ||  |
      //----------------------------------    |                     |

      assembler.template Matrix<Pres,Velx>(shp_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<Pres,Vely>(shp_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<Pres,Velz>(shp_mean.d0, timefacfac*normal(2), shp_jump.d0);

      //   |                      |
      // - |  < q > n, || u_i ||  |
      //   |                      |

      assembler.template Vector<Pres>(shp_mean.d0, -timefacfac*(vmean(0,0)*normal(0)
                                                               +vmean(1,0)*normal(1)
                                                               +vmean(2,0)*normal(2)));

      //--------------------------------------    |                |
      // pressure adjoint consistency term RHS  + |  < q > n, J_u  |
      //--------------------------------------    |                |

      assembler.template Vector<Pres>(shp_mean.d0, timefacfac*ju);

      //------------------------------------------    |                |
      // pressure/flux jump (consistency) term RHS  + |  < v >, J_p*n  |
      //------------------------------------------    |                |

      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(0)*jp);
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(1)*jp);
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(2)*jp);

      //------------------------    |                              |
      // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
      //------------------------    |                              |

      assembler.template Matrix<Velx,Velx>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);

      //    |                               |
      //  - |  \alpha_u || v ||, || u_i ||  |
      //    |                               |

      assembler.template Vector<Velx>(shp_jump.d0, -alphau*timefacfac*vjump(0,0));
      assembler.template Vector<Vely>(shp_jump.d0, -alphau*timefacfac*vjump(1,0));
      assembler.template Vector<Velz>(shp_jump.d0, -alphau*timefacfac*vjump(2,0));

      //----------------------------    |                         |
      // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, J_u  |
      //----------------------------    |                         |

      assembler.template Vector<Velx>(shp_jump.d0, alphau*timefacfac*normal(0)*ju);
      assembler.template Vector<Vely>(shp_jump.d0, alphau*timefacfac*normal(1)*ju);
      assembler.template Vector<Velz>(shp_jump.d0, alphau*timefacfac*normal(2)*ju);

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

    } // loop Gaussian points
  } // loop boundary integration cells

  return;
}

/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void Sysmat(
        const DRT::ELEMENTS::Combust3*    ele,            ///< the element those matrix is calculated
        const Teuchos::RCP<COMBUST::InterfaceHandleCombust>  ih,   ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,         ///< dofmanager of the current element
        const DRT::ELEMENTS::Combust3::MyState&  mystate, ///< element state variables
        Epetra_SerialDenseMatrix&         estif,          ///< element matrix to calculate
        Epetra_SerialDenseVector&         eforce,         ///< element rhs to calculate
        Teuchos::RCP<const MAT::Material> material,       ///< fluid material
        const FLUID_TIMEINTTYPE           timealgo,       ///< time discretization type
        const double                      dt,             ///< delta t (time step size)
        const double                      theta,          ///< factor for one step theta scheme
        const bool                        newton,         ///< full Newton or fixed-point-like
        const bool                        pstab,          ///< flag for stabilisation
        const bool                        supg,           ///< flag for stabilisation
        const bool                        cstab,          ///< flag for stabilisation
        const bool                        instationary,   ///< switch between stationary and instationary formulation
        const INPAR::COMBUST::CombustionType combusttype, ///< switch for type of combusiton problem
        const double                      flamespeed,
        const double                      nitschevel,
        const double                      nitschepres
        )
{
    // initialize element stiffness matrix and force vector
    estif.Scale(0.0);
    eforce.Scale(0.0);

    const int NUMDOF = 4;

    // dead load in element nodes
    //////////////////////////////////////////////////// , LINALG::SerialDenseMatrix edeadng_(BodyForce(ele->Nodes(),time));

    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

    // split velocity and pressure (and stress)
    const int shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
    const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    LINALG::Matrix<shpVecSize,1> eprenp;
    LINALG::Matrix<3,shpVecSize> evelnp;
    LINALG::Matrix<3,shpVecSize> eveln;
    LINALG::Matrix<3,shpVecSize> evelnm;
    LINALG::Matrix<3,shpVecSize> eaccn;
    LINALG::Matrix<numnode,1> ephi;
    LINALG::Matrix<6,shpVecSizeStress> etau;

    fillElementUnknownsArrays<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp,
                                               ephi, etau);

    switch(combusttype)
    {
      case INPAR::COMBUST::combusttype_premixedcombustion:
      {
        SysmatCombustDomain<DISTYPE,ASSTYPE,NUMDOF>(
            ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etau,
            material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler);
        // boundary integrals are only added for intersected elements (fully enriched elements)
        if (ele->Intersected() == true)
        {
          SysmatCombustBoundary<DISTYPE,ASSTYPE,NUMDOF>(
            ele, ih, dofman, evelnp, eprenp, ephi, etau, material, timealgo, dt, theta, assembler,
            flamespeed,nitschevel,nitschepres);
        }
      }
      break;
      case INPAR::COMBUST::combusttype_twophaseflow:
      {
        SysmatTwoPhase<DISTYPE,ASSTYPE,NUMDOF>(
            ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etau,
            material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler);
      }
      break;
      default:
        dserror("unknown type of combustion problem");
    }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::callSysmat(
        const XFEM::AssemblyType             assembly_type,
        const DRT::ELEMENTS::Combust3*       ele,
        const Teuchos::RCP<COMBUST::InterfaceHandleCombust>&  ih,
        const XFEM::ElementDofManager&       eleDofManager,
        const DRT::ELEMENTS::Combust3::MyState&  mystate,   ///< element state variables
        Epetra_SerialDenseMatrix&            estif,
        Epetra_SerialDenseVector&            eforce,
        Teuchos::RCP<const MAT::Material>    material,
        const FLUID_TIMEINTTYPE              timealgo,      ///< time discretization type
        const double                         dt,            ///< delta t (time step size)
        const double                         theta,         ///< factor for one step theta scheme
        const bool                           newton ,
        const bool                           pstab  ,
        const bool                           supg   ,
        const bool                           cstab  ,
        const bool                           instationary,
        const INPAR::COMBUST::CombustionType combusttype,
        const double                         flamespeed,
        const double                         nitschevel,
        const double                         nitschepres
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::hex20:
                Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::hex27:
                Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::tet4:
                Sysmat<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::tet10:
                Sysmat<DRT::Element::tet10,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
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
                Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::hex20:
                Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::hex27:
                Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::tet4:
                Sysmat<DRT::Element::tet4,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            case DRT::Element::tet10:
                Sysmat<DRT::Element::tet10,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary,
                        combusttype,flamespeed,nitschevel,nitschepres);
                break;
            default:
                dserror("xfem_assembly Sysmat not templated yet");
        };
    }
}

#endif
#endif
