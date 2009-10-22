/*----------------------------------------------------------------------*/
/*!
\file xfluid3_sysmat_tp1.cpp

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
#include "xfluid3_interpolation.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/xfem_element_utils.H"
#include "../drt_fluid/time_integration_element.H"
#include "../drt_xfem/spacetime_boundary.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"


  using namespace XFEM::PHYSICS;

  //! size factor to allow fixed size arrays
  ///
  /// to allow fixed size arrays for a unknown number of unknowns, we make them bigger than necessary
  /// this factor is multiplied times numnode(distype) to get the size of many arrays
  template<XFEM::AssemblyType ASSTYPE>
  struct SizeFac {};
  /// specialization of SizeFac for XFEM::standard_assembly
  template<> struct SizeFac<XFEM::standard_assembly> {static const int fac = 1;};
  /// specialization of SizeFac for XFEM::xfem_assembly
  template<> struct SizeFac<XFEM::xfem_assembly>     {static const int fac = 3;};


  //! fill a number of arrays with unknown values from the unknown vector given by the discretization
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            class M1, class V1, class M2, class V2>
  void fillElementUnknownsArraysTP1(
          const XFEM::ElementDofManager& dofman,
          const DRT::ELEMENTS::XFluid3::MyState& mystate,
          M1& evelnp,
          M1& eveln,
          M1& evelnm,
          M1& eaccn,
          V1& eprenp,
          M2& etau,
          V2& ediscprenp
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
      const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Tauxx);
      const size_t numparamdiscpres = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::DiscPres);

      const size_t shpVecSize       = SizeFac<ASSTYPE>::fac*numnode;
      const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
      const size_t shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
      const DRT::Element::DiscretizationType discpresdistype = XFLUID::DiscPressureInterpolation3D<DISTYPE>::distype;
      const size_t shpVecSizeDiscPres = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

      if (numparamvelx > shpVecSize)
      {
        cout << "increase SizeFac for nodal unknowns" << endl;
      }
      if (numparamtauxx > shpVecSizeStress)
      {
        cout << "increase SizeFac for stress unknowns" << endl;
      }
      if (numparamdiscpres > shpVecSizeDiscPres)
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
      const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 0) > 0);
      if (tauele_unknowns_present)
      {
          const std::size_t numparamtauyy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyy, 1);
          const std::size_t numparamtauzz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauzz, 1);
          const std::size_t numparamtauxy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxy, 1);
          const std::size_t numparamtauxz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxz, 1);
          const std::size_t numparamtauyz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyz, 1);
          const std::vector<int>& tauxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauxx>());
          const std::vector<int>& tauyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauyy>());
          const std::vector<int>& tauzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauzz>());
          const std::vector<int>& tauxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauxy>());
          const std::vector<int>& tauxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauxz>());
          const std::vector<int>& tauyzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauyz>());
          for (std::size_t iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = mystate.velnp[tauxxdof[iparam]];
          for (std::size_t iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = mystate.velnp[tauyydof[iparam]];
          for (std::size_t iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = mystate.velnp[tauzzdof[iparam]];
          for (std::size_t iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = mystate.velnp[tauxydof[iparam]];
          for (std::size_t iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = mystate.velnp[tauxzdof[iparam]];
          for (std::size_t iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = mystate.velnp[tauyzdof[iparam]];
      }
      const bool discpres_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 0) > 0);
      if (discpres_unknowns_present)
      {
          const vector<int>& discpresdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::DiscPres>());
          for (std::size_t iparam=0; iparam<numparamdiscpres; ++iparam)   ediscprenp(iparam) = mystate.velnp[discpresdof[iparam]];
      }
  }

  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE,
            int NUMDOF,
            int shpVecSize,
            int shpVecSizeStress>
  void BuildStiffnessMatrixEntries(
      LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>&    assembler,
      const XFEM::ApproxFunc<shpVecSize>&                     shp,
      const LINALG::Matrix<shpVecSizeStress,1>&  shp_tau,
      const LINALG::Matrix<shpVecSizeStress,1>&  shp_discpres,
      const double fac,
      const double timefac,
      const double timefacfac,
      const double  visc,
      const LINALG::Matrix<3,1>& gpvelnp,
      const double             pres,
      const LINALG::Matrix<3,1>& gradp,
      const LINALG::Matrix<3,3>& vderxy,
      const LINALG::Matrix<3,1>& rhsint,
      const LINALG::Matrix<3,1>& res_old,
      const LINALG::Matrix<3,1>& visc_old,
      const LINALG::Matrix<3,3>& tau,
      const double&              discpres,
      const LINALG::Matrix<shpVecSize,1>& enr_conv_c_,
      const XFLUID::EnrViscs2<shpVecSize>& enr_viscs2,
      const bool tauele_unknowns_present,
      const bool instationary,
      const bool newton,
      const bool pstab,
      const bool supg,
      const bool cstab,
      const double tau_stab_Mp,
      const double tau_stab_M,
      const double tau_stab_C
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

  // Hellinger-Reissner terms
  if (tauele_unknowns_present)
  {
      dsassert(tauele_unknowns_present, "there should be disc pressure");
      /* Divergenzfreiheit*/
      /*
                     /               \
                    |                 |
                    | q^e , tr(eps^e) |
                    |                 |
                     \               /
      */
      const double presfac = -1.0;
      assembler.template Matrix<DiscPres,Tauxx>(shp_discpres, presfac*1.0/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<DiscPres,Tauyy>(shp_discpres, presfac*1.0/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<DiscPres,Tauzz>(shp_discpres, presfac*1.0/(2.0*visc)*timefacfac, shp_tau);

      const double trace_tau = (tau(0, 0) + tau(1, 1) + tau(2, 2));
      assembler.template Vector<DiscPres>(shp_discpres, -presfac*1.0/(2.0*visc)*timefacfac*trace_tau);
  }


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

      const double factor = 1.0;
      assembler.template Matrix<Tauxx,Tauxx>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauxy,Tauxy>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauxz,Tauxz>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauyx,Tauyx>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauyy,Tauyy>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauyz,Tauyz>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauzx,Tauzx>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauzy,Tauzy>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);
      assembler.template Matrix<Tauzz,Tauzz>(shp_tau, -    factor/(2.0*visc)*timefacfac, shp_tau);

      assembler.template Vector<Tauxx>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(0,0));
      assembler.template Vector<Tauxy>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(0,1));
      assembler.template Vector<Tauxz>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(0,2));
      assembler.template Vector<Tauyx>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(1,0));
      assembler.template Vector<Tauyy>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(1,1));
      assembler.template Vector<Tauyz>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(1,2));
      assembler.template Vector<Tauzx>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(2,0));
      assembler.template Vector<Tauzy>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(2,1));
      assembler.template Vector<Tauzz>(shp_tau,     factor/(2.0*visc)*timefacfac*tau(2,2));

                   /*                 \
                  | virt tau , eps(Du) |
                   \                 */

      const double rect_factor = 1.0;
      assembler.template Matrix<Tauxx,Velx>(shp_tau,     timefacfac*0.5*rect_factor, shp.dx);
      assembler.template Matrix<Tauxx,Velx>(shp_tau,     timefacfac*0.5*rect_factor, shp.dx);
      assembler.template Matrix<Tauxy,Velx>(shp_tau,     timefacfac*0.5*rect_factor, shp.dy);
      assembler.template Matrix<Tauxy,Vely>(shp_tau,     timefacfac*0.5*rect_factor, shp.dx);
      assembler.template Matrix<Tauxz,Velx>(shp_tau,     timefacfac*0.5*rect_factor, shp.dz);
      assembler.template Matrix<Tauxz,Velz>(shp_tau,     timefacfac*0.5*rect_factor, shp.dx);

      assembler.template Matrix<Tauyx,Vely>(shp_tau,     timefacfac*0.5*rect_factor, shp.dx);
      assembler.template Matrix<Tauyx,Velx>(shp_tau,     timefacfac*0.5*rect_factor, shp.dy);
      assembler.template Matrix<Tauyy,Vely>(shp_tau,     timefacfac*0.5*rect_factor, shp.dy);
      assembler.template Matrix<Tauyy,Vely>(shp_tau,     timefacfac*0.5*rect_factor, shp.dy);
      assembler.template Matrix<Tauyz,Vely>(shp_tau,     timefacfac*0.5*rect_factor, shp.dz);
      assembler.template Matrix<Tauyz,Velz>(shp_tau,     timefacfac*0.5*rect_factor, shp.dy);

      assembler.template Matrix<Tauzx,Velz>(shp_tau,     timefacfac*0.5*rect_factor, shp.dx);
      assembler.template Matrix<Tauzx,Velx>(shp_tau,     timefacfac*0.5*rect_factor, shp.dz);
      assembler.template Matrix<Tauzy,Velz>(shp_tau,     timefacfac*0.5*rect_factor, shp.dy);
      assembler.template Matrix<Tauzy,Vely>(shp_tau,     timefacfac*0.5*rect_factor, shp.dz);
      assembler.template Matrix<Tauzz,Velz>(shp_tau,     timefacfac*0.5*rect_factor, shp.dz);
      assembler.template Matrix<Tauzz,Velz>(shp_tau,     timefacfac*0.5*rect_factor, shp.dz);

      assembler.template Vector<Tauxx>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(0, 0));
      assembler.template Vector<Tauxx>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(0, 0));
      assembler.template Vector<Tauxy>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(0, 1));
      assembler.template Vector<Tauxy>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(1, 0));
      assembler.template Vector<Tauxz>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(0, 2));
      assembler.template Vector<Tauxz>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(2, 0));

      assembler.template Vector<Tauyx>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(1, 0));
      assembler.template Vector<Tauyx>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(0, 1));
      assembler.template Vector<Tauyy>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(1, 1));
      assembler.template Vector<Tauyy>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(1, 1));
      assembler.template Vector<Tauyz>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(1, 2));
      assembler.template Vector<Tauyz>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(2, 1));

      assembler.template Vector<Tauzx>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(2, 0));
      assembler.template Vector<Tauzx>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(0, 2));
      assembler.template Vector<Tauzy>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(2, 1));
      assembler.template Vector<Tauzy>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(1, 2));
      assembler.template Vector<Tauzz>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(2, 2));
      assembler.template Vector<Tauzz>(shp_tau, -    timefacfac*0.5*rect_factor*vderxy(2, 2));

      /* pressure-pressure coupling, quadratic part */
      /*
                     /                      \
                    |                        |
                    | tr(virt tau^e) , p^e I |
                    |                        |
                     \                      /
      */
      assembler.template Matrix<Tauxx,DiscPres>(shp_tau, 1.0/(2.0*visc)*timefacfac, shp_discpres);
      assembler.template Matrix<Tauyy,DiscPres>(shp_tau, 1.0/(2.0*visc)*timefacfac, shp_discpres);
      assembler.template Matrix<Tauzz,DiscPres>(shp_tau, 1.0/(2.0*visc)*timefacfac, shp_discpres);

      assembler.template Vector<Tauxx>(shp_tau, -1.0/(2.0*visc)*timefacfac*discpres);
      assembler.template Vector<Tauyy>(shp_tau, -1.0/(2.0*visc)*timefacfac*discpres);
      assembler.template Vector<Tauzz>(shp_tau, -1.0/(2.0*visc)*timefacfac*discpres);

      /* pressure-pressure coupling, rectangular part */
      /*
                     /                    \
                    |                      |
                  - | tr(virt tau^e) , p I |
                    |                      |
                     \                    /
      */
      assembler.template Matrix<Tauxx,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
      assembler.template Matrix<Tauyy,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);
      assembler.template Matrix<Tauzz,Pres>(shp_tau, -1.0/(2.0*visc)*timefacfac, shp.d0);

      assembler.template Vector<Tauxx>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
      assembler.template Vector<Tauyy>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);
      assembler.template Vector<Tauzz>(shp_tau, 1.0/(2.0*visc)*timefacfac*pres);

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
          class M1, class V1, class M2, class V2>
void SysmatDomainTP1(
    const DRT::Element*                 ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&      dofman,        ///< dofmanager of the current element
    const M1&                           evelnp,
    const M1&                           eveln,
    const M1&                           evelnm,
    const M1&                           eaccn,
    const V1&                           eprenp,
    const M2&                           etau,
    const V2&                           ediscpres,
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
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // space dimension for 3d fluid element
    const int nsd = 3;

    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // get node coordinates of the current element
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    // get older interface velocities and accelerations
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
    // kinematic viscosity \nu
    const double visc = actmat->Viscosity();
    const double dens = actmat->Density();
    // dynamic viscosity \mu
    const double dynvisc = visc * dens;

    // flag for higher order elements
    const bool higher_order_ele = DRT::UTILS::secondDerivativesZero<DISTYPE>();

    const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
    const DRT::Element::DiscretizationType discpresdistype = XFLUID::DiscPressureInterpolation3D<DISTYPE>::distype;

    // figure out whether we have stress unknowns at all
    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Tauxx, 0) > 0);
    const bool discpres_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, DiscPres, 0) > 0);
    if (tauele_unknowns_present != discpres_unknowns_present) dserror("you need disc. pressure unknowns, if you use stress unknowns");
    //cout << "discpres_unknowns_present: " << discpres_unknowns_present << endl;
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
    const int numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const int numparamtauxx = XFEM::getNumParam<ASSTYPE>(dofman, Tauxx, 1);
    const int numparamdiscpres = XFEM::getNumParam<ASSTYPE>(dofman, DiscPres, 1);

    // stabilization parameter
    const double hk = FLD::UTILS::HK<DISTYPE>(evelnp,xyze);
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
          else if (xlabelset.size() > 1) // multiple interface labels
          {
            compute = true;
          }
          else if (xlabelset.find(labelnp) == xlabelset.end()) // ???
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

        const DRT::UTILS::GaussRule3D gaussrule = XFLUID::getXFEMGaussrule<DISTYPE>(ele, xyze, ih->ElementIntersected(ele->Id()),cell->Shape(),false);

        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            static LINALG::Matrix<3,1> pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static LINALG::Matrix<3,1> posXiDomain;
            GEO::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = GEO::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);

            // shape functions and their first derivatives
            static LINALG::Matrix<numnode,1> funct;
            static LINALG::Matrix<3,numnode> deriv;
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
            // discontinouos pressure shape functions
            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement,1> funct_discpres;
            static LINALG::Matrix<3,DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement> deriv_discpres;
            if (ASSTYPE == XFEM::xfem_assembly)
            {
              if (discpres_unknowns_present)
              {
                DRT::UTILS::shape_function_3D(funct_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);
                DRT::UTILS::shape_function_3D_deriv1(deriv_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);
              }
              else
              {
                funct_discpres.Clear();
                deriv_discpres.Clear();
              }
            }

            // get transposed of the jacobian matrix d x / d \xi
            // xjm(i,j) = deriv(i,k)*xyze(j,k)
            static LINALG::Matrix<3,3> xjm;
            xjm.MultiplyNT(deriv,xyze);

            const double det = xjm.Determinant();
            const double fac = intpoints.qwgt[iquad]*det*detcell;

            if (det < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
            }

            // inverse of jacobian
            static LINALG::Matrix<3,3> xji;
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

            const int shpVecSize       = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
            const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
            const int shpVecSizeDiscPres = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

            static XFEM::ApproxFunc<shpVecSize> shp;

            static LINALG::Matrix<shpVecSizeStress,1>   shp_tau;
            static LINALG::Matrix<shpVecSizeDiscPres,1> shp_discpres;

            LINALG::SerialDenseVector enr_funct_discpres(numparamdiscpres);

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

                for (int iparam = 0; iparam != numparamvelx; ++iparam)
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
                    static LINALG::Matrix<shpVecSizeStress,1> enr_funct_stress;

                    // shape functions for element dofs
                    enrvals.ComputeEnrichedElementShapefunction(
                            Tauxx,
                            funct_stress,
                            enr_funct_stress);

                    for (int iparam = 0; iparam < numparamtauxx; ++iparam)
                    {
                      shp_tau(iparam) = enr_funct_stress(iparam);
                    }
                }
                else
                {
                  shp_tau.Clear();
                }

                if (discpres_unknowns_present)
                {
                    static LINALG::SerialDenseVector enr_funct_discpres(3*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement);

                    // shape functions for element dofs
                    enrvals.ComputeEnrichedElementShapefunction(
                            DiscPres,
                            funct_discpres,
                            enr_funct_discpres);

                    for (int iparam = 0; iparam < numparamdiscpres; ++iparam)
                    {
                      shp_discpres(iparam) = enr_funct_discpres(iparam);
                    }
                }
                else
                {
                  shp_discpres = 0.0;
                }
            }
            else
            {
              for (int iparam = 0; iparam != numnode; ++iparam)
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
            const LINALG::Matrix<3,1> gpvelnp = XFLUID::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
            LINALG::Matrix<3,1> gpveln  = XFLUID::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
            LINALG::Matrix<3,1> gpvelnm = XFLUID::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
            LINALG::Matrix<3,1> gpaccn  = XFLUID::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);

            double dtstar = -10000.0;
            if (timealgo != timeint_stationary)
            {
              const bool valid_spacetime_cell_found = XFLUID::modifyOldTimeStepsValues<DISTYPE>(ele, ih, xyze, posXiDomain, labelnp, dt, ivelcolnp, ivelcoln, ivelcolnm, iacccoln, gpveln, gpvelnm, gpaccn, dtstar);
              if (not valid_spacetime_cell_found)
                continue;
            }
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
            const LINALG::Matrix<3,1> histvec = FLD::TIMEINT_THETA_BDF2::GetOldPartOfRighthandside(
                gpveln, gpvelnm, gpaccn, timealgo, dt, theta);

            // get velocity (np,i) derivatives at integration point
            // vderxy = enr_derxy(j,k)*evelnp(i,k);
            static LINALG::Matrix<3,nsd> vderxy;
            vderxy.Clear();
            for (int iparam = 0; iparam < numparamvelx; ++iparam)
            {
              for (int isd = 0; isd < nsd; ++isd)
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
              for (int iparam = 0; iparam < numparamvelx; ++iparam)
              {
                for (int isd = 0; isd != nsd; ++isd)
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
            LINALG::Matrix<nsd,1> gradp;
            gradp.Clear();
            for (int iparam = 0; iparam != numparampres; ++iparam)
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
            for (int iparam = 0; iparam != numparampres; ++iparam)
              pres += shp.d0(iparam)*eprenp(iparam);

            // get discontinous pressure
            double discpres = 0.0;
            for (int iparam = 0; iparam < numparamdiscpres; ++iparam)
              discpres += shp_discpres(iparam)*ediscpres(iparam);

            // get viscous stress unknowns
            static LINALG::Matrix<3,3> tau;
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
            const double vel_norm = gpvelnp.Norm2();
            const double strle = FLD::UTILS::Streamlength(shp.dx, shp.dy, shp.dz, gpvelnp, vel_norm, numparamvelx);
            double tau_stab_M  = 0.0;
            double tau_stab_Mp = 0.0;
            double tau_stab_C  = 0.0;
            FLD::UTILS::computeStabilizationParams(gpvelnp, xji,
                instationary, dynvisc, dens, vel_norm, strle, hk, mk, timefac, dt, INPAR::FLUID::tautype_franca_barrenechea_valentin_wall,
                tau_stab_M, tau_stab_Mp, tau_stab_C);
            // correction due to definition of stabilization parameter in computeStabilizationParams()
            tau_stab_M  /= timefac;
            tau_stab_Mp /= timefac;
            tau_stab_C  /= timefac;



            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            static LINALG::Matrix<3,1> rhsint;
            static LINALG::Matrix<3,1> bodyforce;
            bodyforce.Clear();
            //bodyforce(0) = 1.0;
            for (int isd = 0; isd < nsd; ++isd)
                rhsint(isd) = histvec(isd) + bodyforce(isd)*timefac;

            /*----------------- get numerical representation of single operators ---*/

            /* Convective term  u_old * grad u_old: */
            static LINALG::Matrix<3,1> conv_old;
            //conv_old = vderxy(i, j)*gpvelnp(j);
            conv_old.Multiply(vderxy,gpvelnp);

            /* Viscous term  div epsilon(u_old) */
            static LINALG::Matrix<3,1> visc_old;
            visc_old(0) = vderxy2(0,0) + 0.5 * (vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
            visc_old(1) = vderxy2(1,1) + 0.5 * (vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
            visc_old(2) = vderxy2(2,2) + 0.5 * (vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));

            // evaluate residual once for all stabilisation right hand sides
            static LINALG::Matrix<3,1> res_old;
            for (int isd = 0; isd < nsd; ++isd)
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
            for (int iparam = 0; iparam != numparamvelx; ++iparam)
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

            for (int iparam = 0; iparam != numparamvelx; ++iparam)
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

            BuildStiffnessMatrixEntries<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress>(
                assembler, shp, shp_tau, shp_discpres, fac, timefac, timefacfac, visc,
                gpvelnp, pres, gradp, vderxy, rhsint, res_old, visc_old, tau, discpres,
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
          class M1, class M2, class V2>
void SysmatBoundaryTP1(
    const DRT::Element*               ele,           ///< the element those matrix is calculated
    const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const M2&                         etau,
    const V2&                         ediscpres,
    const Teuchos::RCP<const Epetra_Vector>& ivelcol,       ///< velocity for interface nodes
    const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
    const FLUID_TIMEINTTYPE           timealgo,      ///< time discretization type
    const double                      dt,            ///< delta t (time step size)
    const double                      theta,         ///< factor for one step theta scheme
    LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const bool                        ifaceForceContribution
)
{
    if (ASSTYPE != XFEM::xfem_assembly) dserror("works only with xfem assembly");

    const Epetra_BLAS blas;

    const int nsd = 3;

    // time integration constant
    const double timefac = FLD::TIMEINT_THETA_BDF2::ComputeTimeFac(timealgo, dt, theta);

    // number of nodes for element
    const int numnode_xele = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = XFEM::NumParam<numnode_xele,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
    //const int numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    //const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Tauxx, 1);
    const int numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Tauxx);
    const int numparamdiscpres = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::DiscPres);


    const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, Tauxx, 0) > 0);
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
        const int numnode_boundary = boundaryele->NumNode();

        // get current node coordinates
        LINALG::SerialDenseMatrix xyze_boundary(3,numnode_boundary);
        ih->fillBoundaryNodalPositionsNP(boundaryele, xyze_boundary);

        // get interface velocities at the boundary element nodes
        LINALG::SerialDenseMatrix vel_boundary(3,numnode_boundary);
        const DRT::Node*const* nodes = boundaryele->Nodes();
        for (int inode = 0; inode < numnode_boundary; ++inode)
        {
          const DRT::Node* node = nodes[inode];
          std::vector<int> lm = ih->cutterdis()->Dof(node);
          static std::vector<double> myvel(3);
          DRT::UTILS::ExtractMyValues(*ivelcol,myvel,lm);
          vel_boundary(0,inode) = myvel[0];
          vel_boundary(1,inode) = myvel[1];
          vel_boundary(2,inode) = myvel[2];
        }

        LINALG::SerialDenseMatrix force_boundary(3,numnode_boundary,true);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta^\boundary
            static LINALG::Matrix<2,1> pos_eta_boundary;
            pos_eta_boundary(0) = intpoints.qxg[iquad][0];
            pos_eta_boundary(1) = intpoints.qxg[iquad][1];

            // coordinates of the current integration point in element coordinates \xi^\boundary
            static LINALG::Matrix<2,1> posXiBoundary;
            mapEtaBToXiB(*cell, pos_eta_boundary, posXiBoundary);

            // coordinates of the current integration point in element coordinates \xi^\domain
            static LINALG::Matrix<3,1> posXiDomain;
            mapEtaBToXiD(*cell, pos_eta_boundary, posXiDomain);

            const double detcell = fabs(detEtaBToXiB(*cell, pos_eta_boundary)); //TODO: check normals
            if (detcell < 0.0)
            {
              cout << "detcel :  " << detcell << endl;
              dserror("negative detcell! should be a bug!");
            }

            // shape functions and their first derivatives
            LINALG::SerialDenseVector funct_boundary(DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()));
            DRT::UTILS::shape_function_2D(funct_boundary, posXiBoundary(0),posXiBoundary(1),boundaryele->Shape());
            LINALG::SerialDenseMatrix deriv_boundary(3, DRT::UTILS::getNumberOfElementNodes(boundaryele->Shape()));
            DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posXiBoundary(0),posXiBoundary(1),boundaryele->Shape());

            // shape functions and their first derivatives
            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1> funct;
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

            // stress shape function
            const DRT::Element::DiscretizationType stressdistype = XFLUID::StressInterpolation3D<DISTYPE>::distype;
            static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
            DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);

            // discontinouos pressure shape functions
            const DRT::Element::DiscretizationType discpresdistype = XFLUID::DiscPressureInterpolation3D<DISTYPE>::distype;
            static LINALG::SerialDenseVector funct_discpres(DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement);
            DRT::UTILS::shape_function_3D(funct_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);

            // position of the gausspoint in physical coordinates
            // gauss_pos_xyz = funct_boundary(j)*xyze_boundary(i,j);
            const LINALG::Matrix<3,1> gauss_pos_xyz = XFLUID::interpolateVectorFieldToIntPoint(xyze_boundary,funct_boundary,numnode_boundary);

            // get jacobian matrix d x / d \xi  (3x2)
            // dxyzdrs(i,j) = xyze_boundary(i,k)*deriv_boundary(j,k);
            static LINALG::Matrix<3,2> dxyzdrs;
            blas.GEMM('N','T',3,2,numnode_boundary,1.0,xyze_boundary.A(),xyze_boundary.LDA(),deriv_boundary.A(),deriv_boundary.LDA(),0.0,dxyzdrs.A(),dxyzdrs.M());

            // compute covariant metric tensor G for surface element (2x2)
            // metric = dxyzdrs(k,i)*dxyzdrs(k,j);
            static LINALG::Matrix<2,2> metric;
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

            const int shpVecSize         = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
            const int shpVecSizeStress   = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
            const int shpVecSizeDiscPres = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

            // temporary arrays
            static LINALG::Matrix<shpVecSize,1>          enr_funct;
            static LINALG::Matrix<shpVecSizeStress,1>    enr_funct_stress;
            static LINALG::Matrix<shpVecSizeDiscPres,1>  enr_funct_discpres;

            if (dofman.getUniqueEnrichments().size() > 1)
              dserror("for an intersected element, we assume only 1 enrichment for now!");
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
                    Tauxx,
                    funct_stress,
                    enr_funct_stress);
            // shape functions for element pressure dofs
            enrvals.ComputeEnrichedElementShapefunction(
                    DiscPres,
                    funct_discpres,
                    enr_funct_discpres);

            // perform integration for entire matrix and rhs
            static LINALG::Matrix<shpVecSize,1> shp;
            for (int iparam = 0; iparam < numparamvelx; ++iparam)
            {
              shp(iparam) = enr_funct(iparam);
            }
            static LINALG::Matrix<shpVecSizeStress,1> shp_tau;
            for (int iparam = 0; iparam < numparamtauxx; ++iparam)
            {
              shp_tau(iparam) = enr_funct_stress(iparam);
            }
            static LINALG::Matrix<shpVecSizeDiscPres,1> shp_discpres;
            for (int iparam = 0; iparam < numparamdiscpres; ++iparam)
            {
              shp_discpres(iparam) = enr_funct_discpres(iparam);
            }

            // get normal vector (in physical coordinates) to surface element at integration point
            static LINALG::Matrix<3,1> normalvec_solid;
            GEO::computeNormalToSurfaceElement(boundaryele->Shape(), xyze_boundary, posXiBoundary, normalvec_solid);
            static LINALG::Matrix<3,1> normalvec_fluid(true);
            normalvec_fluid.Update(-1.0,normalvec_solid,0.0);

            // get velocities (n+g,i) at integration point
            // gpvelnp = evelnp(i,j)*shp(j);
            static LINALG::Matrix<3,1> gpvelnp;
            gpvelnp.Clear();
            for (int iparam = 0; iparam < numparamvelx; ++iparam)
                for (int isd = 0; isd < nsd; ++isd)
                    gpvelnp(isd) += evelnp(isd,iparam)*shp(iparam);

            // get interface velocity
            static LINALG::Matrix<3,1> interface_gpvelnp;
            interface_gpvelnp.Clear();
            if (timealgo != timeint_stationary)
                for (int inode = 0; inode < numnode_boundary; ++inode)
                    for (int isd = 0; isd < 3; ++isd)
                        interface_gpvelnp(isd) += vel_boundary(isd,inode)*funct_boundary(inode);

            // get discontinous pressure
            //const double discpres(shp_discpres*ediscprenp);
            double discpres = 0.0;
            for (int iparam = 0; iparam < numparamdiscpres; ++iparam)
              discpres += shp_discpres(iparam)*ediscpres(iparam);

            // get viscous stress unknowns
            static LINALG::Matrix<3,3> tau;
            XFLUID::fill_tau(numparamtauxx, shp_tau, etau, tau);

            // integration factors and coefficients of single terms
            const double timefacfac = timefac * fac;

            //////////////////////////////////////
            // now build single stiffness terms //
            //////////////////////////////////////

               /*                      \
            - |  (virt tau) * n^f , Du  |
               \                      */

            assembler.template Matrix<Tauxx,Velx>(shp_tau, -timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Tauxy,Velx>(shp_tau, -timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Tauxz,Velx>(shp_tau, -timefacfac*normalvec_fluid(2), shp);
            assembler.template Matrix<Tauyx,Vely>(shp_tau, -timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Tauyy,Vely>(shp_tau, -timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Tauyz,Vely>(shp_tau, -timefacfac*normalvec_fluid(2), shp);
            assembler.template Matrix<Tauzx,Velz>(shp_tau, -timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<Tauzy,Velz>(shp_tau, -timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<Tauzz,Velz>(shp_tau, -timefacfac*normalvec_fluid(2), shp);

            assembler.template Vector<Tauxx>(shp_tau, timefacfac*normalvec_fluid(0)*gpvelnp(0));
            assembler.template Vector<Tauxy>(shp_tau, timefacfac*normalvec_fluid(1)*gpvelnp(0));
            assembler.template Vector<Tauxz>(shp_tau, timefacfac*normalvec_fluid(2)*gpvelnp(0));
            assembler.template Vector<Tauyx>(shp_tau, timefacfac*normalvec_fluid(0)*gpvelnp(1));
            assembler.template Vector<Tauyy>(shp_tau, timefacfac*normalvec_fluid(1)*gpvelnp(1));
            assembler.template Vector<Tauyz>(shp_tau, timefacfac*normalvec_fluid(2)*gpvelnp(1));
            assembler.template Vector<Tauzx>(shp_tau, timefacfac*normalvec_fluid(0)*gpvelnp(2));
            assembler.template Vector<Tauzy>(shp_tau, timefacfac*normalvec_fluid(1)*gpvelnp(2));
            assembler.template Vector<Tauzz>(shp_tau, timefacfac*normalvec_fluid(2)*gpvelnp(2));


               /*                            \
              |  (virt tau) * n^f , u^\iface  |
               \                            */

            assembler.template Vector<Tauxx>(shp_tau, -timefacfac*normalvec_fluid(0)*interface_gpvelnp(0));
            assembler.template Vector<Tauxy>(shp_tau, -timefacfac*normalvec_fluid(1)*interface_gpvelnp(0));
            assembler.template Vector<Tauxz>(shp_tau, -timefacfac*normalvec_fluid(2)*interface_gpvelnp(0));
            assembler.template Vector<Tauyx>(shp_tau, -timefacfac*normalvec_fluid(0)*interface_gpvelnp(1));
            assembler.template Vector<Tauyy>(shp_tau, -timefacfac*normalvec_fluid(1)*interface_gpvelnp(1));
            assembler.template Vector<Tauyz>(shp_tau, -timefacfac*normalvec_fluid(2)*interface_gpvelnp(1));
            assembler.template Vector<Tauzx>(shp_tau, -timefacfac*normalvec_fluid(0)*interface_gpvelnp(2));
            assembler.template Vector<Tauzy>(shp_tau, -timefacfac*normalvec_fluid(1)*interface_gpvelnp(2));
            assembler.template Vector<Tauzz>(shp_tau, -timefacfac*normalvec_fluid(2)*interface_gpvelnp(2));


               /*                      \
              |  (virt p^e) * n^f , Du  |
               \                      */

            assembler.template Matrix<DiscPres,Velx>(shp_discpres, timefacfac*normalvec_fluid(0), shp);
            assembler.template Matrix<DiscPres,Vely>(shp_discpres, timefacfac*normalvec_fluid(1), shp);
            assembler.template Matrix<DiscPres,Velz>(shp_discpres, timefacfac*normalvec_fluid(2), shp);

            assembler.template Vector<DiscPres>(shp_discpres, -timefacfac*normalvec_fluid(0)*gpvelnp(0));
            assembler.template Vector<DiscPres>(shp_discpres, -timefacfac*normalvec_fluid(1)*gpvelnp(1));
            assembler.template Vector<DiscPres>(shp_discpres, -timefacfac*normalvec_fluid(2)*gpvelnp(2));

               /*                            \
            - |  (virt p^e) * n^f , u^\iface  |
               \                            */

            assembler.template Vector<DiscPres>(shp_discpres, timefacfac*normalvec_fluid(0)*interface_gpvelnp(0));
            assembler.template Vector<DiscPres>(shp_discpres, timefacfac*normalvec_fluid(1)*interface_gpvelnp(1));
            assembler.template Vector<DiscPres>(shp_discpres, timefacfac*normalvec_fluid(2)*interface_gpvelnp(2));

               /*               \
            - |  v , Dtau * n^f  |
               \               */

            assembler.template Matrix<Velx,Tauxx>(shp, -timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Velx,Tauxy>(shp, -timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Velx,Tauxz>(shp, -timefacfac*normalvec_fluid(2), shp_tau);
            assembler.template Matrix<Vely,Tauyx>(shp, -timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Vely,Tauyy>(shp, -timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Vely,Tauyz>(shp, -timefacfac*normalvec_fluid(2), shp_tau);
            assembler.template Matrix<Velz,Tauzx>(shp, -timefacfac*normalvec_fluid(0), shp_tau);
            assembler.template Matrix<Velz,Tauzy>(shp, -timefacfac*normalvec_fluid(1), shp_tau);
            assembler.template Matrix<Velz,Tauzz>(shp, -timefacfac*normalvec_fluid(2), shp_tau);

            static LINALG::Matrix<3,1> disctau_times_n;
            disctau_times_n.Multiply(tau,normalvec_fluid);
            //cout << "sigmaijnj : " << disctau_times_n << endl;
            assembler.template Vector<Velx>(shp, timefacfac*disctau_times_n(0));
            assembler.template Vector<Vely>(shp, timefacfac*disctau_times_n(1));
            assembler.template Vector<Velz>(shp, timefacfac*disctau_times_n(2));

                         /*         \
                        |  v , Dp n  |
                         \         */

            assembler.template Matrix<Velx,DiscPres>(shp, timefacfac*normalvec_fluid(0), shp_discpres);
            assembler.template Matrix<Vely,DiscPres>(shp, timefacfac*normalvec_fluid(1), shp_discpres);
            assembler.template Matrix<Velz,DiscPres>(shp, timefacfac*normalvec_fluid(2), shp_discpres);

            assembler.template Vector<Velx>(shp, -timefacfac*normalvec_fluid(0)*discpres);
            assembler.template Vector<Vely>(shp, -timefacfac*normalvec_fluid(1)*discpres);
            assembler.template Vector<Velz>(shp, -timefacfac*normalvec_fluid(2)*discpres);


            // here the interface force is integrated
            // this is done using test shape functions of the boundary mesh
            // hence, we can't use the local assembler here
            for (int inode = 0; inode < numnode_boundary; ++inode)
            {
              force_boundary(0,inode) += funct_boundary(inode) * ((- discpres*normalvec_fluid(0) + disctau_times_n(0)) * timefacfac);
              force_boundary(1,inode) += funct_boundary(inode) * ((- discpres*normalvec_fluid(1) + disctau_times_n(1)) * timefacfac);
              force_boundary(2,inode) += funct_boundary(inode) * ((- discpres*normalvec_fluid(2) + disctau_times_n(2)) * timefacfac);
            }

        } // end loop over gauss points

        // here we need to assemble into the global force vector of the boundary discretization
        // note that we assemble into a overlapping vector, hence we add only, if we are a xfem row element
        // this way, we can later add all contributions together when exporting to interface row elements
        if (ifaceForceContribution)
        {
          const Epetra_Map* dofcolmap = ih->cutterdis()->DofColMap();
          std::vector<int> gdofs(3);
          for (int inode = 0; inode < numnode_boundary; ++inode)
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
void SysmatTP1(
        const DRT::Element*               ele,           ///< the element those matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const DRT::ELEMENTS::XFluid3::MyState&  mystate,  ///< element state variables
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
        const bool                        ifaceForceContribution
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
    const DRT::Element::DiscretizationType discpresdistype = XFLUID::DiscPressureInterpolation3D<DISTYPE>::distype;
    const int shpVecSizeStress = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    const int shpVecSizeDiscPres = SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;
    static LINALG::Matrix<shpVecSize,1> eprenp;
    static LINALG::Matrix<3,shpVecSize> evelnp;
    static LINALG::Matrix<3,shpVecSize> eveln;
    static LINALG::Matrix<3,shpVecSize> evelnm;
    static LINALG::Matrix<3,shpVecSize> eaccn;
    static LINALG::Matrix<6,shpVecSizeStress> etau;
    static LINALG::Matrix<shpVecSizeDiscPres,1> ediscpres;

    fillElementUnknownsArraysTP1<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp, etau, ediscpres);

    SysmatDomainTP1<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, etau, ediscpres,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, assembler);

    if (ASSTYPE == XFEM::xfem_assembly)
    {
      SysmatBoundaryTP1<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, etau, ediscpres, ivelcol, iforcecol,
          timealgo, dt, theta, assembler, ifaceForceContribution);
    }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUID::callSysmatTP1(
        const XFEM::AssemblyType          assembly_type,
        const DRT::ELEMENTS::XFluid3*     ele,
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,
        const XFEM::ElementDofManager&    eleDofManager,
        const DRT::ELEMENTS::XFluid3::MyState&  mystate,   ///< element state variables
        const Teuchos::RCP<const Epetra_Vector>& ivelcol,
        const Teuchos::RCP<Epetra_Vector>& iforcecol,     ///< reaction force due to given interface velocity
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
        const bool                        ifaceForceContribution
        )
{
    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                SysmatTP1<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
//            case DRT::Element::hex20:
//                SysmatTP1<DRT::Element::hex20,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
//                break;
//            case DRT::Element::hex27:
//                SysmatTP1<DRT::Element::hex27,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
//                break;
            case DRT::Element::tet4:
                SysmatTP1<DRT::Element::tet4,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
//            case DRT::Element::tet10:
//                SysmatTP1<DRT::Element::tet4,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
//                break;
            default:
                dserror("standard_assembly Sysmat not templated yet");
        };
    }
    else
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                SysmatTP1<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
                break;
//            case DRT::Element::hex20:
//                SysmatTP1<DRT::Element::hex20,XFEM::xfem_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
//                break;
//            case DRT::Element::hex27:
//                SysmatTP1<DRT::Element::hex27,XFEM::xfem_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
//                break;
//            case DRT::Element::tet10:
//                SysmatTP1<DRT::Element::tet4,XFEM::standard_assembly>(
//                        ele, ih, eleDofManager, mystate, ivelcol, iforcecol, estif, eforce,
//                        material, timealgo, dt, theta, newton, pstab, supg, cstab, instationary, ifaceForceContribution);
//                break;
            default:
                dserror("xfem_assembly Sysmat not templated yet");
        };
    }
}

#endif
#endif
