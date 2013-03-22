/*----------------------------------------------------------------------*/
/*!
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

#ifndef COMBUST3_SYSMAT_PREMIXED_STRESS_H
#define COMBUST3_SYSMAT_PREMIXED_STRESS_H


#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "combust3_stabilization.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_normals.H"


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
          size_t shpVecSizeStress,
          size_t shpVecSizeDiscPres>
void BuildDomainIntegralsStress(
    COMBUST::LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFunc<2,shpVecSize>& shp,                  // shape functions and derivatives
    const XFEM::ApproxFunc<0,shpVecSizeStress>& shp_epsilon,        // shape functions for stresses
    const XFEM::ApproxFunc<0,shpVecSizeDiscPres>& shp_discpres, // shape functions for disc. pressure
    const double& fac,                                          // factor for integration in space
    const double& timefac,                                      // factor for integration in time
    const double& timefacfac,                                   // overall integration factor
    const double& dynvisc,                                      // dynamic viscosity at n+1
    const double& pres,                                         // pressure at GP at n+1
    const LINALG::Matrix<3,3>& vderxy,                          // velocity derivative at GP at n+1
    const LINALG::Matrix<3,3>& epsilon,                         // element stress unknowns
    const double& discpres                                      // discontinuous pressure at GP at n+1
)
{
#ifdef COMBUST_EPSPRES_BASED_STRAIN_COUPLING
  // epsilon-epsilon coupling, discontinuous epsilon term
  /*
              /                              \
           - |  virt epsilon , 2\mu*Depsilon  |
              \                              /
   */
  const double twomu = 2.0*dynvisc;
  // remark: the symmetry properties are used here (XFEM::PHYSICS::Epsilonxy = XFEM::PHYSICS::Epsilonyx, XFEM::PHYSICS::Epsilonxz = XFEM::PHYSICS::Epsilonzx, XFEM::PHYSICS::Epsilonyz = XFEM::PHYSICS::Epsilonzy)
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac,     shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*2.0, shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*2.0, shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac,     shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*2.0, shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac,     shp_epsilon.d0);

  assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,0)    );
  assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,1)*2.0);
  assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,2)*2.0);
  assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,1)    );
  assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,2)*2.0);
  assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,2)    );

  // formulation without using symmetry of stress tensor
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonyx,XFEM::PHYSICS::Epsilonyx>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonzx,XFEM::PHYSICS::Epsilonzx>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonzy,XFEM::PHYSICS::Epsilonzy>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);

  //  assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,0));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,1));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,2));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,0));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,1));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,2));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,0));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,1));
  //  assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,2));

  // epsilon-epsilon coupling, continuous epsilon term
  /*               /                      \
                  | virt epsilon , 2\mu*Du |
                   \                      /
   */
  // remark: the symmetry properties are used here (XFEM::PHYSICS::Epsilonxy = XFEM::PHYSICS::Epsilonyx, XFEM::PHYSICS::Epsilonxz = XFEM::PHYSICS::Epsilonzx, XFEM::PHYSICS::Epsilonyz = XFEM::PHYSICS::Epsilonzy)
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac , shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac , shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac , shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac , shp.dz);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac , shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac , shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac , shp.dz);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac , shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac , shp.dz);

  assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac*vderxy(0,0));
  assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*(vderxy(0,1) + vderxy(1,0)));
  assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*(vderxy(0,2) + vderxy(2,0)));
  assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac*vderxy(1,1));
  assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*(vderxy(1,2) + vderxy(2,1)));
  assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac*vderxy(2,2));

  // formulation without using symmetry of stress tensor
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  //
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonyx,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonyx,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  //
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonzx,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonzx,XFEM::PHYSICS::Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonzy,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonzy,XFEM::PHYSICS::Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  //
  // assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,0));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,0));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,1));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,0));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,2));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,0));
  //
  // assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,0));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,1));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,1));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,1));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,2));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,1));
  //
  // assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,0));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,2));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,1));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,2));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,2));
  // assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,2));
#endif

#ifdef COMBUST_EPSPRES_BASED_PRESSURE_COUPLING
  // pressure-pressure coupling, discontinuous pressure term
  /*                 /                      \
                    | tr(virt epsilon) , p^e |
                     \                      /
   */
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::DiscPres>(shp_epsilon.d0, timefacfac, shp_discpres.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::DiscPres>(shp_epsilon.d0, timefacfac, shp_discpres.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::DiscPres>(shp_epsilon.d0, timefacfac, shp_discpres.d0);

  assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, -timefacfac*discpres);
  assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, -timefacfac*discpres);
  assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, -timefacfac*discpres);

  // pressure-pressure coupling, continuous pressure term
  /*                 /                    \
                  - | tr(virt epsilon) , p |
                     \                    /
   */
  assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Pres>(shp_epsilon.d0, -timefacfac, shp.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Pres>(shp_epsilon.d0, -timefacfac, shp.d0);
  assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Pres>(shp_epsilon.d0, -timefacfac, shp.d0);

  assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon.d0, timefacfac*pres);
  assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon.d0, timefacfac*pres);
  assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon.d0, timefacfac*pres);
#endif

#ifdef COMBUST_EPSPRES_BASED_ZERO_TRACE
  // enforce zero trace of discontinuous epsilon tensor (solenoidality of auxiliar field)
  /*                 /               \
                  - | q^e , tr(eps^e) |
                     \               /
   */
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Epsilonxx>(shp_discpres.d0, -timefacfac, shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Epsilonyy>(shp_discpres.d0, -timefacfac, shp_epsilon.d0);
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Epsilonzz>(shp_discpres.d0, -timefacfac, shp_epsilon.d0);

  // compute trace of discontinuous epsilon tensor
  const double trace_epsilon = (epsilon(0,0) + epsilon(1,1) + epsilon(2,2));
  assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres.d0, timefacfac*trace_epsilon);

  // enforce zero trace of continuous epsilon tensor (solenoidality of standard field)
  /*                 /                \
                    | q^e , tr(eps(u)) |
                     \                /
  */
  // new term; seems to have effect on third diagonal matrix entry
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velx>(shp_discpres.d0, timefacfac, shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Vely>(shp_discpres.d0, timefacfac, shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velz>(shp_discpres.d0, timefacfac, shp.dz);

  const double trace = vderxy(0,0) + vderxy(1,1) + vderxy(2,2);
  assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres.d0, -timefacfac*trace);
#endif

#ifdef  COMBUST_EPSPRES_BASED_FORTHTERM
  const double viscfac = 3.0/(2.0*dynvisc);
  // TODO the function of this forth term is not clear yet
  /*                 /                  \
                  - | q^e , 3/(2\mu) p^e |
                     \                  /
   */
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::DiscPres>(shp_discpres.d0, viscfac*timefacfac, shp_discpres.d0);
  assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres.d0, -viscfac*timefacfac*discpres);
  // new term
  /*                 /                \
                    | q^e , 3/(2\mu) p |
                     \                /
   */
  assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Pres>(shp_discpres.d0, -viscfac*timefacfac, shp.d0);
  assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres.d0, viscfac*timefacfac*pres);
#endif
}

template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          size_t shpVecSize,
          size_t shpVecSizeStress,
          size_t shpVecSizeDiscPres>
void BuildDomainIntegralsSigma(
    COMBUST::LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFunc<2,shpVecSize>& shp,                  // shape functions and derivatives
    const XFEM::ApproxFunc<0,shpVecSizeStress>& shp_sigma,      // shape functions for stresses
    const double& fac,                                          // factor for integration in space
    const double& timefac,                                      // factor for integration in time
    const double& timefacfac,                                   // overall integration factor
    const double& dynvisc,                                      // dynamic viscosity at n+1
    const double& pres,                                         // pressure at GP at n+1
    const LINALG::Matrix<3,3>& vderxy,                          // velocity derivative at GP at n+1
    const LINALG::Matrix<3,3>& sigma                            // element stress unknowns
)
{
#ifdef COMBUST_SIGMA_BASED_FIELD_COUPLING
  //
  /*         /                     \
          - |  virt tau , eps(Dtau) |
             \                     /
   */
  const double reciprocal_viscfac = 1.0/(2.0*dynvisc);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxx,XFEM::PHYSICS::Sigmaxx>(shp_sigma.d0, -reciprocal_viscfac*timefacfac,     shp_sigma.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxy,XFEM::PHYSICS::Sigmaxy>(shp_sigma.d0, -reciprocal_viscfac*timefacfac*2.0, shp_sigma.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxz,XFEM::PHYSICS::Sigmaxz>(shp_sigma.d0, -reciprocal_viscfac*timefacfac*2.0, shp_sigma.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmayy,XFEM::PHYSICS::Sigmayy>(shp_sigma.d0, -reciprocal_viscfac*timefacfac,     shp_sigma.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmayz,XFEM::PHYSICS::Sigmayz>(shp_sigma.d0, -reciprocal_viscfac*timefacfac*2.0, shp_sigma.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmazz,XFEM::PHYSICS::Sigmazz>(shp_sigma.d0, -reciprocal_viscfac*timefacfac,     shp_sigma.d0);

  assembler.template Vector<XFEM::PHYSICS::Sigmaxx>(shp_sigma.d0, reciprocal_viscfac*timefacfac*sigma(0,0)    );
  assembler.template Vector<XFEM::PHYSICS::Sigmaxy>(shp_sigma.d0, reciprocal_viscfac*timefacfac*sigma(0,1)*2.0);
  assembler.template Vector<XFEM::PHYSICS::Sigmaxz>(shp_sigma.d0, reciprocal_viscfac*timefacfac*sigma(0,2)*2.0);
  assembler.template Vector<XFEM::PHYSICS::Sigmayy>(shp_sigma.d0, reciprocal_viscfac*timefacfac*sigma(1,1)    );
  assembler.template Vector<XFEM::PHYSICS::Sigmayz>(shp_sigma.d0, reciprocal_viscfac*timefacfac*sigma(1,2)*2.0);
  assembler.template Vector<XFEM::PHYSICS::Sigmazz>(shp_sigma.d0, reciprocal_viscfac*timefacfac*sigma(2,2)    );

  // stress-sigma coupling, sigma term
  /*               /                    \
                  | virt sigma , eps(Du) |
                   \                    /
   */
  // remark: the symmetry properties are used here (XFEM::PHYSICS::Sigmaxy = XFEM::PHYSICS::Sigmayx, XFEM::PHYSICS::Sigmaxz = XFEM::PHYSICS::Sigmazx, XFEM::PHYSICS::Sigmayz = XFEM::PHYSICS::Sigmazy)
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxx,XFEM::PHYSICS::Velx>(shp_sigma.d0, timefacfac , shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxy,XFEM::PHYSICS::Velx>(shp_sigma.d0, timefacfac , shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxy,XFEM::PHYSICS::Vely>(shp_sigma.d0, timefacfac , shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxz,XFEM::PHYSICS::Velx>(shp_sigma.d0, timefacfac , shp.dz);
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxz,XFEM::PHYSICS::Velz>(shp_sigma.d0, timefacfac , shp.dx);
  assembler.template Matrix<XFEM::PHYSICS::Sigmayy,XFEM::PHYSICS::Vely>(shp_sigma.d0, timefacfac , shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::Sigmayz,XFEM::PHYSICS::Vely>(shp_sigma.d0, timefacfac , shp.dz);
  assembler.template Matrix<XFEM::PHYSICS::Sigmayz,XFEM::PHYSICS::Velz>(shp_sigma.d0, timefacfac , shp.dy);
  assembler.template Matrix<XFEM::PHYSICS::Sigmazz,XFEM::PHYSICS::Velz>(shp_sigma.d0, timefacfac , shp.dz);

  assembler.template Vector<XFEM::PHYSICS::Sigmaxx>(shp_sigma.d0, -timefacfac*vderxy(0,0));
  assembler.template Vector<XFEM::PHYSICS::Sigmaxy>(shp_sigma.d0, -timefacfac*(vderxy(0,1) + vderxy(1,0)));
  assembler.template Vector<XFEM::PHYSICS::Sigmaxz>(shp_sigma.d0, -timefacfac*(vderxy(0,2) + vderxy(2,0)));
  assembler.template Vector<XFEM::PHYSICS::Sigmayy>(shp_sigma.d0, -timefacfac*vderxy(1,1));
  assembler.template Vector<XFEM::PHYSICS::Sigmayz>(shp_sigma.d0, -timefacfac*(vderxy(1,2) + vderxy(2,1)));
  assembler.template Vector<XFEM::PHYSICS::Sigmazz>(shp_sigma.d0, -timefacfac*vderxy(2,2));
#endif

#ifdef COMBUST_SIGMA_BASED_TRACE
  // pressure-pressure coupling, continuous pressure term
  /*                 /                           \
                    | tr(virt sigma) , 1/(2\mu) p |
                     \                           /
   */
  assembler.template Matrix<XFEM::PHYSICS::Sigmaxx,XFEM::PHYSICS::Pres>(shp_sigma.d0, -reciprocal_viscfac*timefacfac, shp.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmayy,XFEM::PHYSICS::Pres>(shp_sigma.d0, -reciprocal_viscfac*timefacfac, shp.d0);
  assembler.template Matrix<XFEM::PHYSICS::Sigmazz,XFEM::PHYSICS::Pres>(shp_sigma.d0, -reciprocal_viscfac*timefacfac, shp.d0);

  assembler.template Vector<XFEM::PHYSICS::Sigmaxx>(shp_sigma.d0, reciprocal_viscfac*timefacfac*pres);
  assembler.template Vector<XFEM::PHYSICS::Sigmayy>(shp_sigma.d0, reciprocal_viscfac*timefacfac*pres);
  assembler.template Vector<XFEM::PHYSICS::Sigmazz>(shp_sigma.d0, reciprocal_viscfac*timefacfac*pres);
#endif
}

/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class V1, class M2, class V2, class V3>
void SysmatDomainStress(
    const DRT::ELEMENTS::Combust3*      ele,          ///< element whose matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih, ///< information about the interface
    const XFEM::ElementDofManager&      dofman,       ///< dofmanager of this element
    const M1&                           evelnp,       ///< nodal velocity n+1 values for this element
    const M1&                           eveln,        ///< nodal velocity n   values for this element
    const M1&                           evelnm,       ///< nodal velocity n-1 values for this element
    const M1&                           eaccn,        ///< nodal acceleration n values for this element
    const M1&                           eaccam,       ///< nodal acceleration n values for this element
    const V1&                           eprenp,       ///< nodal pressure n+1 values for this element
    const V2&                           ephi,         ///< nodal G-function values for this element
    const M2&                           eepsilon,     ///< element stresses
    const V3&                           ediscpres,    ///<
    Teuchos::RCP<const MAT::Material>   material,     ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,     ///< time discretization scheme
    const double                        dt,           ///< delta t (time step size)
    const double                        theta,        ///< parameter for one step theta scheme
    const double                        ga_alphaF,
    const double                        ga_alphaM,
    const double                        ga_gamma,
    const bool                          newton,       ///< full Newton or fixed-point-like
    const bool                          pstab,        ///< flag for stabilization
    const bool                          supg,         ///< flag for stabilization
    const bool                          graddiv,        ///< flag for stabilization
    const INPAR::FLUID::TauType         tautype,      ///< stabilization parameter definition
    const bool                          instationary, ///< switch between stationary and instationary formulation
    const bool                          genalpha,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>&   assembler
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

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;

  const bool velocity_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, 1) > 0);
  dsassert(velocity_unknowns_present == true, "there are no velocity unknowns!");
  const bool pressure_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Pres, 1) > 0);
  dsassert(pressure_unknowns_present == true, "there are no pressure unknowns!");

  // figure out whether we have stress unknowns at all
  const bool epsilonele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonxx, 0) > 0);
  const bool sigmaele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxx, 0) > 0);
  const bool discpres_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 0) > 0);
  if (ele->Bisected())
  {
#ifdef COMBUST_EPSPRES_BASED
    dsassert(epsilonele_unknowns_present == true, "there are no element stress unknowns!");
    dsassert(discpres_unknowns_present == true, "there are no discontinuous pressure unknowns!");
#endif
#ifdef COMBUST_SIGMA_BASED
    dsassert(sigmaele_unknowns_present == true, "there are no element stress unknowns!");
#endif
  }

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Epsilonxx);
// TODO @Florian clarify
//cout << "standard element: numparamepsilonxx " << numparamepsilonxx << endl;
//cout << "epsilonele_unknowns_present " << epsilonele_unknowns_present << endl;
  const size_t numparamdiscpres = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 1);
//cout << "standard element: numparamdiscpres " << numparamdiscpres << endl;
//cout << "discpres_unknowns_present " <<  discpres_unknowns_present << endl;
  const size_t numparamsigmaxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
//cout << "sigmaele_unknowns_present " << sigmaele_unknowns_present << endl;

  // stabilization parameter
  const double hk = COMBUST::UTILS::HK<DISTYPE>(evelnp,xyze);
  const double mk = COMBUST::UTILS::MK<DISTYPE>();

  // get domain integration cells for this element
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

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
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussrule<DISTYPE>(ele, xyze, ele->Bisected(), cell->Shape());
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

      // discontinuous stress shape functions
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
      if (ASSTYPE == XFEM::xfem_assembly)
      {
        if ((epsilonele_unknowns_present) or
            (sigmaele_unknowns_present))
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
      const size_t shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<2,shpVecSize> shp;
      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_epsilon;
      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_sigma;
      static XFEM::ApproxFunc<0,shpVecSizeDiscPres> shp_discpres;

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
            XFEM::PHYSICS::Velx,
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

        if (epsilonele_unknowns_present)
        {
          // shape functions for element dofs
          static LINALG::Matrix<shpVecSizeStress,1> enrfunct_stress;
          enrfunct_stress.Clear();
          enrvals.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Epsilonxx, funct_stress, enrfunct_stress);

          for (size_t iparam = 0; iparam < numparamepsilonxx; ++iparam)
          {
            shp_epsilon.d0(iparam) = enrfunct_stress(iparam);
          }
        }
        else if (sigmaele_unknowns_present)
        {
          // shape functions for element dofs
          static LINALG::Matrix<shpVecSizeStress,1> enrfunct_stress;
          enrfunct_stress.Clear();
          enrvals.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Sigmaxx, funct_stress, enrfunct_stress);

          for (size_t iparam = 0; iparam < numparamsigmaxx; ++iparam)
          {
            shp_sigma.d0(iparam) = enrfunct_stress(iparam);
          }
        }
        else
        {
          shp_epsilon.d0.Clear();
          shp_sigma.d0.Clear();
        }
        if (discpres_unknowns_present)
        {
          // shape functions for element dofs
          static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres;
          enrfunct_discpres.Clear();

          enrvals.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::DiscPres, funct_discpres, enrfunct_discpres);

          for (size_t iparam = 0; iparam < numparamdiscpres; ++iparam)
          {
            shp_discpres.d0(iparam) = enrfunct_discpres(iparam);
          }
        }
        else
        {
          shp_discpres.d0.Clear();
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
        if ((epsilonele_unknowns_present) or
            (sigmaele_unknowns_present))
        {
          dserror("no stress enrichments without xfem assembly");
        }
      }

      //------------------------------------------------------
      // get velocities and accelerations at integration point
      //------------------------------------------------------
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPoint(eveln , shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPoint(evelnm, shp.d0, numparamvelx);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPoint(eaccn , shp.d0, numparamvelx);
      // get acceleration at n+alpha_m for generalized alpha time integration scheme
      LINALG::Matrix<3,1> gpaccam(true);
      if (genalpha)
        gpaccam  = XFEM::interpolateVectorFieldToIntPoint(eaccam , shp.d0, numparamvelx);

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

      // get discontinous pressure
      double discpres = 0.0;
      if (discpres_unknowns_present)
      {
        for (size_t iparam = 0; iparam < numparamdiscpres; ++iparam)
          discpres += shp_discpres.d0(iparam)*ediscpres(iparam);
      }

      // get viscous stress unknowns
      static LINALG::Matrix<nsd,nsd> epsilon;
      if (epsilonele_unknowns_present)
      {
        XFEM::fill_tensor(numparamepsilonxx, shp_epsilon.d0, eepsilon, epsilon);
//cout << "epsilon " << epsilon << endl;
//cout << "discpres " << discpres << endl;
      }
      else
      {
        epsilon.Clear();
      }

      // get viscous stress unknowns
      static LINALG::Matrix<nsd,nsd> sigma;
      if (sigmaele_unknowns_present)
      {
        XFEM::fill_tensor(numparamsigmaxx, shp_sigma.d0, eepsilon, sigma);
//cout << "sigma " << sigma << endl;
      }
      else
      {
        sigma.Clear();
      }

      //------------------------------------------
      // compute stabilization parameters (3 taus)
      //------------------------------------------
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
          bodyforce(isd) += edead(isd,inode) * funct(inode);
      }
      // scale body force: body force * density^{n+1}
      bodyforce.Scale(densaf);

      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = bodyforce(isd);
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          rhsint(isd) = histvec(isd) + timefac*bodyforce(isd);
      }

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
      if (genalpha)
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = densaf*conv_old(isd) + gradp(isd) - 2.0*dynvisc*visc_old(isd) - bodyforce(isd);

        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpaccam(isd);
        }
      }
      else
      {
        for (size_t isd = 0; isd < nsd; ++isd)
          res_old(isd) = -rhsint(isd)+timefac*(densaf*conv_old(isd)+gradp(isd)-2.0*dynvisc*visc_old(isd));

        if (instationary)
        {
          for(size_t isd=0; isd < nsd; ++isd)
            res_old(isd) += densaf * gpvelnp(isd);
        }

      }

      // integration factors and coefficients of single terms
      const double timefacfac = timefac * fac;

      double rhsfac = fac;
      if (genalpha)
        rhsfac *= timefac / ga_alphaF;

      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      if (genalpha)
        dserror("SysmatDomainStress: Genalpha has not been tested for this, use at your own risk.");

      COMBUST::BuildStandardDomainIntegrals<DISTYPE,ASSTYPE,NUMDOF,shpVecSize>(
          assembler, shp, shp, fac, timefac, rhsfac, densaf, dynvisc,
          gpvelnp, pres, gpaccam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_conv_c, enr_viscs2,
          instationary, genalpha, newton, pstab, supg, graddiv,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);

      if (epsilonele_unknowns_present)
      {
      COMBUST::BuildDomainIntegralsStress<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress,shpVecSizeDiscPres>(
          assembler, shp, shp_epsilon, shp_discpres, fac, timefac, timefacfac,
          dynvisc, pres, vderxy, epsilon, discpres);
      }
      if (sigmaele_unknowns_present)
      {
        COMBUST::BuildDomainIntegralsSigma<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress,shpVecSizeDiscPres>(
            assembler, shp, shp_sigma, fac, timefac, timefacfac,
            dynvisc, pres, vderxy, epsilon);
      }

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
          class M1, class M2, class M3, class V1, class V2, class V3>
void SysmatBoundaryStress(
    const DRT::ELEMENTS::Combust3*    ele,                    ///< the element
    const COMBUST::InterfaceHandleCombust*&  ih, ///< interface information
    const XFEM::ElementDofManager&    dofman,                 ///< dofmanager of this element
    const M1&                         evelnp,                 ///< element velocity field
    const V1&                         eprenp,                 ///< element pressure field
    const V2&                         ephi,                   ///< element G-function field
    const M2&                         egradphi,               ///<
    const M3&                         eepsilon,               ///< discontinuous (deviatoric) stress
    const V3&                         ediscpres,              ///< discontinuous pressure
    Teuchos::RCP<const MAT::Material> material,               ///< fluid materials
    const INPAR::FLUID::TimeIntegrationScheme timealgo,       ///< time discretization type
    const double&                     dt,                     ///< delta t (time step size)
    const double&                     theta,                  ///< factor for one step theta scheme
    const double                      ga_alphaF,
    const double                      ga_alphaM,
    const double                      ga_gamma,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,      ///< local assembler
    const double                      flamespeed              ///< flame speed
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
  xyze.Clear();
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

  //-----------------------------
  // prepare interface conditions
  //-----------------------------
  // mass flux through flame front
  const double mflux = flamespeed*dens_minus;
#if 0
  // jump velocity vector ju_{n,t1,t2}
  LINALG::Matrix<nsd,1> ju_nt(true);
  ju_nt(0,0) = -mflux*(1.0/dens_minus - 1.0/dens_plus);
  ju_nt(1,0) = 0.0;
  ju_nt(2,0) = 0.0;
#else
  // velocity jump in normal direction
  const double ju_n = -mflux*(1.0/dens_minus - 1.0/dens_plus);
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTPRESJUMP
  // pressure jump value
  const double jp = 0.0; // (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);
#endif

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Epsilonxx);
  const size_t numparamdiscpres = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 1);

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
    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
    // TODO: finish -> template? Only tri3 boundary integration cells allowed?
    // choose an (arbitrary) number of Gaussian points for the output
    DRT::UTILS::GaussRule2D intrule2D = DRT::UTILS::intrule2D_undefined;
    if (cell->Shape() == DRT::Element::tri3)
    {
      // there should only be triangular boundary cells
      // TODO: are 6 Gaussian points enough?
      //const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
      intrule2D = DRT::UTILS::intrule_tri_37point;
    }
    else if (cell->Shape() == DRT::Element::quad4)
    {
      //dserror("triangular boundary integration cell expected");
      intrule2D = DRT::UTILS::intrule_quad_25point;
    }
    else
    {
      dserror("invalid type of boundary integration cell");
    }

    dserror("triangular boundary integration cell expected");
    // here, a triangular boundary integration cell is assumed (numvertices = 3)
    const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;

    const DRT::UTILS::IntegrationPoints2D intpoints(intrule2D);

    // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
    LINALG::SerialDenseMatrix cellXiDomaintmp = cell->CellNodalPosXiDomain();
    // get coordinates of vertices of boundary integration cell in physical coordinates \x^domain
    LINALG::SerialDenseMatrix cellXYZDomaintmp = cell->CellNodalPosXYZ();
    // transform to fixed size format
    const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);
    const LINALG::Matrix<nsd,numvertices> cellXYZDomain(cellXYZDomaintmp);

    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of this integration point in boundary cell coordinates \eta^boundary
      const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);

      // coordinates of this integration point in element coordinates \xi^domain
      LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

      //--------------------------------------------------------------------------------------------
      // compute normal vector (normed, in physical coordinates)
      // remark: for linear boundary integrastion cells this could be done before the loop over all
      //         Gaussian points
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> normal(true);
      normal.Clear();
      GEO::computeNormalToSurfaceElement(cell->Shape(), cellXYZDomain, posEtaBoundary, normal);

      // shape functions for fields with dofs at nodes
      static LINALG::Matrix<numnode,1> functele;
      functele.Clear();
      DRT::UTILS::shape_function_3D(functele,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      // TODO normal vector modified
#ifdef COLLAPSE_FLAME
      normal.Clear();
      for (unsigned i=0;i<numnode;i++)
      {
        normal(0) += functele(i)*xyze(0,i);
        normal(1) += functele(i)*xyze(1,i);
      }
      normal.Scale(1.0);
#endif
#ifdef FLAME_VORTEX
      normal.Clear();
      for (unsigned i=0;i<numnode;i++)
      {
        normal(0) += functele(i)*egradphi(0,i);
        normal(1) += functele(i)*egradphi(1,i);
      }
      normal.Scale(-1.0);
#endif
#ifdef COMBUST_2D
      normal(2) = 0.0;
#endif
      const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal.Scale(1.0/norm);

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

      // jump in normal stress
      const double jflux_n = mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);
      // jump of traction vector (stress tensor * normal)
      LINALG::Matrix<nsd,1> jflux(true);
      jflux(0) = jflux_n*normal(0);
      jflux(1) = jflux_n*normal(1);
      jflux(2) = jflux_n*normal(2);
      //cout << "jflux " << jflux << endl;
      // overwrite velocity jump condition (e.g. for specific test case)
      //jflux(0) = 0.0;
      //jflux(1) = 0.0;
      //jflux(2) = 0.0;

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
      // shape functions for element stress field
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
      funct_stress.Clear();
      DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
      // shape functions for discontinuous element pressure field
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement,1> funct_discpres;
      funct_discpres.Clear();
      DRT::UTILS::shape_function_3D(funct_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);

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

      //------------------------------------------------
      // compute Jacobian matrix for mapping to boundary
      //------------------------------------------------
      static LINALG::Matrix<nsd,numvertices> deriv_boundary;
      deriv_boundary.Clear();
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posEtaBoundary(0),posEtaBoundary(1),cell->Shape());

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

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const std::size_t shpVecSize         = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const std::size_t shpVecSizeStress   = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
      const std::size_t shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

      // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
      static LINALG::Matrix<shpVecSize,1>         enrfunct_plus;
      enrfunct_plus.Clear();
      static LINALG::Matrix<shpVecSize,1>         enrfunct_minus;
      enrfunct_minus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_plus;
      enrfunct_stress_plus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_minus;
      enrfunct_stress_minus.Clear();
      static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres_plus;
      enrfunct_discpres_plus.Clear();
      static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres_minus;
      enrfunct_discpres_minus.Clear();

      static LINALG::Matrix<3,shpVecSize> enrderxy_plus;
      enrderxy_plus.Clear();
      static LINALG::Matrix<3,shpVecSize> enrderxy_minus;
      enrderxy_minus.Clear();
      static LINALG::Matrix<6,shpVecSize> enrderxy2_plus;
      enrderxy2_plus.Clear();
      static LINALG::Matrix<6,shpVecSize> enrderxy2_minus;
      enrderxy2_minus.Clear();

      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_plus, enrderxy_plus, enrderxy2_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_minus, enrderxy_minus, enrderxy2_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Epsilonxx, funct_stress, enrfunct_stress_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Epsilonxx, funct_stress, enrfunct_stress_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::DiscPres, funct_discpres, enrfunct_discpres_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::DiscPres, funct_discpres, enrfunct_discpres_minus);

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<0,shpVecSize> shp_jump;
      static XFEM::ApproxFunc<0,shpVecSize> shp_mean;
      static XFEM::ApproxFunc<2,shpVecSize> shp_mean_visc;
      static XFEM::ApproxFunc<2,shpVecSize> shp_jump_visc;

      // fill approximation functions
      for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        shp_jump.d0(iparam) = enrfunct_minus(iparam) - enrfunct_plus(iparam);
        //cout << "shp_jump.d0: " << shp_jump.d0(iparam) << endl;
        shp_mean.d0(iparam) = 0.5*(enrfunct_plus(iparam) + enrfunct_minus(iparam));
        //cout << "shp_mean.d0: " << shp_mean.d0(iparam) << endl;
      }

      //static XFEM::ApproxFunc<0,shpVecSizeStress> shp_epsilon_jump;
      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_eps_visc_mean;
      for (std::size_t iparam = 0; iparam < numparamepsilonxx; ++iparam)
      {
        //shp_epsilon_jump.d0(iparam) = 2.0*dynvisc_minus*enrfunct_stress_minus(iparam) - 2.0*dynvisc_plus*enrfunct_stress_plus(iparam);
        // 1/2*(2*\mu^- * \eps^- + 2*\mu^+ * \eps^+)
        shp_eps_visc_mean.d0(iparam) = dynvisc_minus*enrfunct_stress_minus(iparam) + dynvisc_plus*enrfunct_stress_plus(iparam);
      }

      static XFEM::ApproxFunc<0,shpVecSizeDiscPres> shp_discpres_jump;
      static XFEM::ApproxFunc<0,shpVecSizeDiscPres> shp_discpres_mean;
      for (size_t iparam = 0; iparam < numparamdiscpres; ++iparam)
      {
        shp_discpres_jump.d0(iparam) = enrfunct_discpres_minus(iparam) - enrfunct_discpres_plus(iparam);
        shp_discpres_mean.d0(iparam) = 0.5*(enrfunct_discpres_minus(iparam) + enrfunct_discpres_plus(iparam));
      }

      //---------------------------------------
      // compute data at Gaussian point for rhs
      //---------------------------------------

      // velocity jump
      static LINALG::Matrix<nsd,1> veljump(true);
      veljump.Clear();
      veljump = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);

      // velocity minus
      static LINALG::Matrix<nsd,1> velminus(true);
      velminus.Clear();
      velminus = XFEM::interpolateVectorFieldToIntPoint(evelnp, enrfunct_minus, numparamvelx);

      // velocity plus
      static LINALG::Matrix<nsd,1> velplus(true);
      velplus.Clear();
      velplus = XFEM::interpolateVectorFieldToIntPoint(evelnp, enrfunct_plus, numparamvelx);

      // mean velocity
      static LINALG::Matrix<nsd,1> velmean(true);
      velmean.Clear();
      velmean = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_mean.d0, numparamvelx);

      // get discontinuous pressure jump
      static double pdiscjump;
      pdiscjump = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pdiscjump += shp_discpres_jump.d0(iparam)*ediscpres(iparam);

      // get discontinuous pressure average
      static double pdiscmean;
      pdiscmean = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pdiscmean += shp_discpres_mean.d0(iparam)*ediscpres(iparam);

//cout << "pdiscjump " << pdiscjump << endl;
//cout << "pdiscmean " << pdiscmean << endl;

      // get viscous stress unknowns
      //static LINALG::Matrix<nsd,nsd> epsilonjump;
      //epsilonjump.Clear();
      //XFEM::fill_tensor(numparamepsilonxx, shp_epsilon_jump.d0, eepsilon, epsilonjump);

      static LINALG::Matrix<nsd,nsd> epsilonmean;
      epsilonmean.Clear();
      XFEM::fill_tensor(numparamepsilonxx, shp_eps_visc_mean.d0, eepsilon, epsilonmean);

//cout << "epsilonjump " << epsilonjump << endl;
//cout << "epsilonmean " << epsilonmean << endl;
#ifdef COMBUST_EPSPRES_BASED_NEUMANN_INFLOW
//      static XFEM::ApproxFunc<0,shpVecSize> shp_inflow;
//      const double direction = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
//      if (direction < 0.0)
//      {
////        cout << "inflow right" << endl;
//        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
//        {
//          shp_inflow.d0(iparam) = enrfunct_minus(iparam);
//        }
//      }
//      else
//      {
////        cout << "inflow left" << endl;
//        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
//        {
//          shp_inflow.d0(iparam) = enrfunct_plus(iparam);
//        }
//      }
//      // inflow velocity
//      static LINALG::Matrix<nsd,1> velinflow(true);
//      velinflow.Clear();
//      velinflow = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_inflow.d0, numparamvelx);
//      //cout << velinflow << endl;
#endif

#ifdef COMBUST_EPSPRES_BASED_ADVECTIVEFLUX
      static XFEM::ApproxFunc<0,shpVecSize> shp_inflow;
      const double direction = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      if (direction < 0.0)
        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          shp_inflow.d0(iparam) = enrfunct_plus(iparam);
        }
      else
        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          shp_inflow.d0(iparam) = enrfunct_minus(iparam);
        }
      // inflow velocity
      static LINALG::Matrix<nsd,1> velinflow(true);
      velinflow.Clear();
      velinflow = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_inflow.d0, numparamvelx);
#endif
//--------------------------------------------------------------------------------------------------
      // build single boundary integral stiffness terms
      //
      // remarks: || x || stands for the jump operator ( || x || = x^- - x^+ )
      //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
      //             n    stands for the normal vector on the interface pointing from the positive (+)
      //                  domain into the negative (-) domain (n = n^+ = -n^-)
      //--------------------------------------------------------------------------------------------

#ifdef COMBUST_EPSPRES_BASED_CONSISTVISCJUMP
      //-------------------------    |                            |
      // viscous consistency term  + |  < v > , || Depsilon || n  |
      //-------------------------    |                            |

      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxx>(shp_mean.d0, timefacfac*normal(0), shp_epsilon_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxy>(shp_mean.d0, timefacfac*normal(1), shp_epsilon_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxz>(shp_mean.d0, timefacfac*normal(2), shp_epsilon_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyx>(shp_mean.d0, timefacfac*normal(0), shp_epsilon_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyy>(shp_mean.d0, timefacfac*normal(1), shp_epsilon_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyz>(shp_mean.d0, timefacfac*normal(2), shp_epsilon_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzx>(shp_mean.d0, timefacfac*normal(0), shp_epsilon_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzy>(shp_mean.d0, timefacfac*normal(1), shp_epsilon_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzz>(shp_mean.d0, timefacfac*normal(2), shp_epsilon_jump.d0);

      //   |                             |
      // - |  < v > , || epsilon_i || n  |
      //   |                             |

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*normal(0)*epsilonjump(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*normal(1)*epsilonjump(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*normal(2)*epsilonjump(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*normal(0)*epsilonjump(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*normal(1)*epsilonjump(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*normal(2)*epsilonjump(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*normal(0)*epsilonjump(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*normal(1)*epsilonjump(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*normal(2)*epsilonjump(2,2));
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTPRESJUMP
      //--------------------------    |                       |
      // pressure consistency term  - |  < v >, || Dp^e || n  |
      //--------------------------    |                       |

      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::DiscPres>(shp_mean.d0, -timefacfac*normal(0), shp_discpres_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::DiscPres>(shp_mean.d0, -timefacfac*normal(1), shp_discpres_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::DiscPres>(shp_mean.d0, -timefacfac*normal(2), shp_discpres_jump.d0);

      //   |                        |
      // + |  < v >, || p^e_i || n  |
      //   |                        |

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(0)*pdiscjump);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(1)*pdiscjump);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(2)*pdiscjump);

      //---------  ---------------------------------    |                |
      // pressure (flux) jump (consistency) term RHS  + |  < v >, j_p*n  |
      //-----------  -------------------------------    |                |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(0)*jp);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(1)*jp);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(2)*jp);
#endif

#ifdef COMBUST_EPSPRES_BASED_FLUXJUMP
      //---------------------------------------------------    |                 |
      // flux (traction/stress) jump (consistency) term RHS  - |  < v >, j_flux  |
      //---------------------------------------------------    |                 |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*jflux(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*jflux(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*jflux(2));

#if 0
      // j_sigma is a 3x3 Matrix -> j_sigma*normal
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(0)*jflux(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(1)*jflux(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(2)*jflux(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(0)*jflux(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(1)*jflux(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(2)*jflux(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(0)*jflux(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(1)*jflux(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(2)*jflux(2,2));
#endif
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTVISCMEAN
      //-------------------------    |                            |
      // viscous consistency term  + |  || v || , < Depsilon > n  |
      //-------------------------    |                            |
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxx>(shp_jump.d0, timefacfac*normal(0), shp_eps_visc_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxy>(shp_jump.d0, timefacfac*normal(1), shp_eps_visc_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxz>(shp_jump.d0, timefacfac*normal(2), shp_eps_visc_mean.d0);

      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyx>(shp_jump.d0, timefacfac*normal(0), shp_eps_visc_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyy>(shp_jump.d0, timefacfac*normal(1), shp_eps_visc_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyz>(shp_jump.d0, timefacfac*normal(2), shp_eps_visc_mean.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzx>(shp_jump.d0, timefacfac*normal(0), shp_eps_visc_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzy>(shp_jump.d0, timefacfac*normal(1), shp_eps_visc_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzz>(shp_jump.d0, timefacfac*normal(2), shp_eps_visc_mean.d0);

      //   |                             |
      // - |  || v || , < epsilon_i > n  |
      //   |                             |

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*normal(0)*epsilonmean(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*normal(1)*epsilonmean(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*normal(2)*epsilonmean(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*normal(0)*epsilonmean(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*normal(1)*epsilonmean(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*normal(2)*epsilonmean(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*normal(0)*epsilonmean(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*normal(1)*epsilonmean(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*normal(2)*epsilonmean(2,2));
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTPRESMEAN
      //--------------------------    |                       |
      // pressure consistency term  - |  || v ||, < Dp^e > n  |
      //--------------------------    |                       |

      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::DiscPres>(shp_jump.d0, -timefacfac*normal(0), shp_discpres_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::DiscPres>(shp_jump.d0, -timefacfac*normal(1), shp_discpres_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::DiscPres>(shp_jump.d0, -timefacfac*normal(2), shp_discpres_mean.d0);

      //   |                        |
      // + |  || v ||, < p^e_i > n  |
      //   |                        |

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, timefacfac*normal(0)*pdiscmean);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, timefacfac*normal(1)*pdiscmean);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, timefacfac*normal(2)*pdiscmean);
#endif

#ifdef COMBUST_EPSPRES_BASED_ADJCONSISTVISCJUMP
      //---------------------------------    |                             |
      // viscous adjoint consistency term  + |  < epsilon^e > n, || Du ||  |
      //---------------------------------    |                             |
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Velx>(shp_eps_visc_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Velx>(shp_eps_visc_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velx>(shp_eps_visc_mean.d0, timefacfac*normal(2), shp_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Epsilonyx,XFEM::PHYSICS::Vely>(shp_eps_visc_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Vely>(shp_eps_visc_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Vely>(shp_eps_visc_mean.d0, timefacfac*normal(2), shp_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Epsilonzx,XFEM::PHYSICS::Velz>(shp_eps_visc_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonzy,XFEM::PHYSICS::Velz>(shp_eps_visc_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Velz>(shp_eps_visc_mean.d0, timefacfac*normal(2), shp_jump.d0);

      //   |                              |
      // - |  < epsilon^e > n, || u_i ||  |
      //   |                              |

      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_eps_visc_mean.d0, -timefacfac*normal(0)*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_eps_visc_mean.d0, -timefacfac*normal(1)*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_eps_visc_mean.d0, -timefacfac*normal(2)*veljump(0));

      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_eps_visc_mean.d0, -timefacfac*normal(0)*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_eps_visc_mean.d0, -timefacfac*normal(1)*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_eps_visc_mean.d0, -timefacfac*normal(2)*veljump(1));

      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_eps_visc_mean.d0, -timefacfac*normal(0)*veljump(2));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_eps_visc_mean.d0, -timefacfac*normal(1)*veljump(2));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_eps_visc_mean.d0, -timefacfac*normal(2)*veljump(2));

      //-------------------------------------    |                        |
      // viscous adjoint consistency term RHS  + |  < epsilon^e > n, j_u  |
      //-------------------------------------    |                        |

#if 1
      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*ju(0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*ju(0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*ju(0));
#else
      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*normal(0)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*normal(0)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*normal(0)*ju_n);
#endif

#if 1
      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*ju(1));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*ju(1));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*ju(1));
#else
      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*normal(1)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*normal(1)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*normal(1)*ju_n);
#endif

#if 1
      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*ju(2));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*ju(2));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*ju(2));
#else
      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*normal(2)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*normal(2)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*normal(2)*ju_n);
#endif
#endif

#ifdef COMBUST_EPSPRES_BASED_ADJCONSISTPRESJUMP
      //----------------------------------      |                       |
      // pressure adjoint consistency term  +/- |  < q^e >, || Du || n  |
      //----------------------------------      |                       |
      // remark: the sign of this term is not clear, yet
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velx>(shp_discpres_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Vely>(shp_discpres_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velz>(shp_discpres_mean.d0, timefacfac*normal(2), shp_jump.d0);

      //     |                        |
      // -/+ |  < q^e >, || u_i || n  |
      //     |                        |
      // remark: the sign of this term is not clear, yet
      assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_mean.d0, -timefacfac*(veljump(0)*normal(0)
                                                                            +veljump(1)*normal(1)
                                                                            +veljump(2)*normal(2)));

      //--------------------------------------      |                  |
      // pressure adjoint consistency term RHS  +/- |  < q^e >, j_u n  |
      //--------------------------------------      |                  |
      // remark: the sign of this term is not clear, yet
//       assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_mean.d0, -timefacfac*ju);
       assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_mean.d0, timefacfac*(ju(0)*normal(0)
                                                                            +ju(1)*normal(1)
                                                                            +ju(2)*normal(2)));
#endif

#ifdef COMBUST_EPSPRES_BASED_NEUMANN_INFLOW
#if 0
      //-----------------    |                   |
      // inflow flux term  - |  v, u_n*|| Du ||  |
      //-----------------    |                   |
      //const double inflowvel = velinflow(0)*normal(0) + velinflow(1)*normal(1) + velinflow(2)*normal(2);
      const double inflowvel = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);

      //   |                    |
      // + |  v, u_n*|| u_i ||  |
      //   |                    |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_inflow.d0, timefacfac*inflowvel*velinflow(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_inflow.d0, timefacfac*inflowvel*velinflow(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_inflow.d0, timefacfac*inflowvel*velinflow(2));
#endif
      //-----------------    |                       |
      // inflow flux term  - |  < v >, u_n*|| Du ||  |
      //-----------------    |                       |
      const double inflowvel = -1.0*(velmean(0)*normal(0) + velmean(1)*normal(1) + velmean(2)*normal(2));
      //const double inflowvel = 1000*(veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2));
//      cout << "inflowvel" << inflowvel << endl;
//      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);

//      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);

      //   |                        |
      // + |  < v >, u_n*|| u_i ||  |
      //   |                        |
//      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel*veljump(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel*veljump(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel*veljump(2));

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*inflowvel*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*inflowvel*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*inflowvel*veljump(2));

//      assembler.template Vector<XFEM::PHYSICS::Velx>(enrfunct_minus, timefacfac*inflowvel*veljump(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(enrfunct_minus, timefacfac*inflowvel*veljump(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(enrfunct_minus, timefacfac*inflowvel*veljump(2));
      //   |                 |
      // + |  < v >, u_n*ju  |
      //   |                 |
//      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*inflowvel*ju(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*inflowvel*ju(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*inflowvel*ju(2));

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, timefacfac*inflowvel*ju(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, timefacfac*inflowvel*ju(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, timefacfac*inflowvel*ju(2));

//      assembler.template Vector<XFEM::PHYSICS::Velx>(enrfunct_minus, -timefacfac*inflowvel*ju(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(enrfunct_minus, -timefacfac*inflowvel*ju(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(enrfunct_minus, -timefacfac*inflowvel*ju(2));

#endif

#ifdef COMBUST_EPSPRES_BASED_ADVECTIVEFLUX
      //-----------------    |                   |
      // inflow flux term  + |  < v >, || Du ||  |
      //-----------------    |                   |
      const double inflowvel = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);

      //   |                    |
      // - |  < v >, || u_i ||  |
      //   |                    |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*inflowvel*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*inflowvel*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*inflowvel*veljump(2));

      //-----------------       |                |
      // inflow flux term RHS + |  < v >, j_u n  |
      //-----------------       |                |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel*ju*normal(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel*ju*normal(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel*ju*normal(2));

      //-----------------    |                     |
      // inflow flux term  - |  < q >, || Du || n  |
      //-----------------    |                     |
      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel*normal(2), shp_jump.d0);

      //   |                      |
      // + |  < q >, || u_i || n  |
      //   |                      |
      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, -timefacfac*inflowvel*(veljump(0)*normal(0)
                                                                        +veljump(1)*normal(1)
                                                                        +veljump(2)*normal(2)));

      //-----------------       |                |
      // inflow flux term RHS - |  < q >, j_u n  |
      //-----------------       |                |
      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, timefacfac*inflowvel*ju);
#endif
#ifdef COMBUST_STRESS_NITSCHEVEL
    const double alphau = 10.0;
    //------------------------    |                              |
    // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
    //------------------------    |                              |

    assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
    assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);

    //    |                               |
    //  - |  \alpha_u || v ||, || u_i ||  |
    //    |                               |

    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -alphau*timefacfac*veljump(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -alphau*timefacfac*veljump(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -alphau*timefacfac*veljump(2));

    //----------------------------    |                         |
    // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, j_u  |
    //----------------------------    |                         |

    //      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, alphau*timefacfac*ju(0));
    //      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, alphau*timefacfac*ju(1));
    //      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, alphau*timefacfac*ju(2));
    assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, alphau*timefacfac*ju(0));
    assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, alphau*timefacfac*ju(1));
    assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, alphau*timefacfac*ju(2));
#endif
    } // loop Gaussian points
  } // loop boundary integration cells

  return;
}

/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class M3, class V1, class V2, class V3>
void SysmatBoundarySigma(
    const DRT::ELEMENTS::Combust3*    ele,                    ///< the element
    const COMBUST::InterfaceHandleCombust*&  ih, ///< interface information
    const XFEM::ElementDofManager&    dofman,                 ///< dofmanager of this element
    const M1&                         evelnp,                 ///< element velocity field
    const V1&                         eprenp,                 ///< element pressure field
    const V2&                         ephi,                   ///< element G-function field
    const M2&                         egradphi,               ///<
    const M3&                         esigma,                 ///< discontinuous (deviatoric) stress
    const V3&                         ediscpres,              ///< discontinuous pressure
    Teuchos::RCP<const MAT::Material> material,               ///< fluid materials
    const INPAR::FLUID::TimeIntegrationScheme timealgo,       ///< time discretization type
    const double&                     dt,                     ///< delta t (time step size)
    const double&                     theta,                  ///< factor for one step theta scheme
    const double                      ga_alphaF,
    const double                      ga_alphaM,
    const double                      ga_gamma,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,      ///< local assembler
    const double                      flamespeed              ///< flame speed
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

  //-----------------------------
  // prepare interface conditions
  //-----------------------------
  // mass flux through flame front
  const double mflux = flamespeed*dens_minus;
#if 0
  // jump velocity vector ju_{n,t1,t2}
  LINALG::Matrix<nsd,1> ju_nt(true);
  ju_nt(0,0) = -mflux*(1.0/dens_minus - 1.0/dens_plus);
  ju_nt(1,0) = 0.0;
  ju_nt(2,0) = 0.0;
#else
  // velocity jump in normal direction
  const double ju_n = -mflux*(1.0/dens_minus - 1.0/dens_plus);
#endif

#ifdef COMBUST_SIGMA_BASED_CONSISTPRESJUMP
  // pressure jump value
  const double jp = 0.0; // (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);
#endif

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamsigmaxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);

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
    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
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

    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of this integration point in boundary cell coordinates \eta^boundary
      const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);

      // coordinates of this integration point in element coordinates \xi^domain
      LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

      //--------------------------------------------------------------------------------------------
      // compute normal vector (normed, in physical coordinates)
      // remark: for linear boundary integrastion cells this could be done before the loop over all
      //         Gaussian points
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> normal(true);
      normal.Clear();
      GEO::computeNormalToSurfaceElement(cell->Shape(), cellXYZDomain, posEtaBoundary, normal);

      // shape functions for fields with dofs at nodes
      static LINALG::Matrix<numnode,1> functele;
      functele.Clear();
      DRT::UTILS::shape_function_3D(functele,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

      // TODO normal vector modified
#ifdef COLLAPSE_FLAME
      normal.Clear();
      for (unsigned i=0;i<numnode;i++)
      {
        normal(0) += functele(i)*xyze(0,i);
        normal(1) += functele(i)*xyze(1,i);
      }
      normal.Scale(-1.0);
#endif
#ifdef FLAME_VORTEX
      normal.Clear();
      for (unsigned i=0;i<numnode;i++)
      {
        normal(0) += functele(i)*egradphi(0,i);
        normal(1) += functele(i)*egradphi(1,i);
      }
      normal.Scale(-1.0);
#endif
#ifdef COMBUST_2D
      normal(2) = 0.0;
#endif
      double norm = sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2));
      if (norm == 0.0) dserror("norm of normal vector is zero!");
      normal.Scale(1.0/norm);

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

      // jump in normal stress
      const double jflux_n = -mflux*mflux*(1.0/dens_minus - 1.0/dens_plus);
      // jump of traction vector (stress tensor * normal)
      LINALG::Matrix<nsd,1> jflux(true);
      jflux(0) = jflux_n*normal(0);
      jflux(1) = jflux_n*normal(1);
      jflux(2) = jflux_n*normal(2);
      //cout << "jflux " << jflux << endl;
      // overwrite velocity jump condition (e.g. for specific test case)
      //jflux(0) = 0.0;
      //jflux(1) = 0.0;
      //jflux(2) = 0.0;

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
      // shape functions for element stress field
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
      funct_stress.Clear();
      DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
      // shape functions for discontinuous element pressure field
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement,1> funct_discpres;
      funct_discpres.Clear();
      DRT::UTILS::shape_function_3D(funct_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);

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

      //------------------------------------------------
      // compute Jacobian matrix for mapping to boundary
      //------------------------------------------------
      static LINALG::Matrix<nsd,numvertices> deriv_boundary;
      deriv_boundary.Clear();
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, posEtaBoundary(0),posEtaBoundary(1),cell->Shape());

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

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const std::size_t shpVecSize         = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const std::size_t shpVecSizeStress   = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;

      // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
      static LINALG::Matrix<shpVecSize,1>         enrfunct_plus;
      enrfunct_plus.Clear();
      static LINALG::Matrix<shpVecSize,1>         enrfunct_minus;
      enrfunct_minus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_plus;
      enrfunct_stress_plus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_minus;
      enrfunct_stress_minus.Clear();

      static LINALG::Matrix<3,shpVecSize> enrderxy_plus;
      enrderxy_plus.Clear();
      static LINALG::Matrix<3,shpVecSize> enrderxy_minus;
      enrderxy_minus.Clear();
      static LINALG::Matrix<6,shpVecSize> enrderxy2_plus;
      enrderxy2_plus.Clear();
      static LINALG::Matrix<6,shpVecSize> enrderxy2_minus;
      enrderxy2_minus.Clear();

      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_plus, enrderxy_plus, enrderxy2_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_minus, enrderxy_minus, enrderxy2_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Sigmaxx, funct_stress, enrfunct_stress_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Sigmaxx, funct_stress, enrfunct_stress_minus);

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<0,shpVecSize> shp_jump;
      static XFEM::ApproxFunc<0,shpVecSize> shp_mean;

      // fill approximation functions
      for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        shp_mean.d0(iparam) = 0.5*(enrfunct_plus(iparam) + enrfunct_minus(iparam));
        //cout << "shp_mean.d0: " << shp_mean.d0(iparam) << endl;
      }

      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_sigma_mean;
      for (std::size_t iparam = 0; iparam < numparamsigmaxx; ++iparam)
      {
        //shp_sigma_jump.d0(iparam) = enrfunct_stress_minus(iparam) - enrfunct_stress_plus(iparam);
        shp_sigma_mean.d0(iparam) = 0.5*(enrfunct_stress_minus(iparam) + enrfunct_stress_plus(iparam));
      }

      //---------------------------------------
      // compute data at Gaussian point for rhs
      //---------------------------------------

      // velocity jump
      static LINALG::Matrix<nsd,1> veljump(true);
      veljump.Clear();
      veljump = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_jump.d0, numparamvelx);

      // velocity minus
      static LINALG::Matrix<nsd,1> velminus(true);
      velminus.Clear();
      velminus = XFEM::interpolateVectorFieldToIntPoint(evelnp, enrfunct_minus, numparamvelx);

      // velocity plus
      static LINALG::Matrix<nsd,1> velplus(true);
      velplus.Clear();
      velplus = XFEM::interpolateVectorFieldToIntPoint(evelnp, enrfunct_plus, numparamvelx);

      // mean velocity
      static LINALG::Matrix<nsd,1> velmean(true);
      velmean.Clear();
      velmean = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_mean.d0, numparamvelx);

      // get stress unknowns
      static LINALG::Matrix<nsd,nsd> sigmamean;
      sigmamean.Clear();
      XFEM::fill_tensor(numparamsigmaxx, shp_sigma_mean.d0, esigma, sigmamean);
//cout << "sigmamean " << sigmamean << endl;

#ifdef COMBUST_SIGMA_BASED_NEUMANN_INFLOW
//      static XFEM::ApproxFunc<0,shpVecSize> shp_inflow;
//      const double direction = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
//      if (direction < 0.0)
//      {
////        cout << "inflow right" << endl;
//        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
//        {
//          shp_inflow.d0(iparam) = enrfunct_minus(iparam);
//        }
//      }
//      else
//      {
////        cout << "inflow left" << endl;
//        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
//        {
//          shp_inflow.d0(iparam) = enrfunct_plus(iparam);
//        }
//      }
//      // inflow velocity
//      static LINALG::Matrix<nsd,1> velinflow(true);
//      velinflow.Clear();
//      velinflow = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_inflow.d0, numparamvelx);
//      //cout << velinflow << endl;
#endif

#ifdef COMBUST_SIGMA_BASED_ADVECTIVEFLUX
      static XFEM::ApproxFunc<0,shpVecSize> shp_inflow;
      const double direction = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      if (direction < 0.0)
        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          shp_inflow.d0(iparam) = enrfunct_plus(iparam);
        }
      else
        for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
        {
          shp_inflow.d0(iparam) = enrfunct_minus(iparam);
        }
      // inflow velocity
      static LINALG::Matrix<nsd,1> velinflow(true);
      velinflow.Clear();
      velinflow = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_inflow.d0, numparamvelx);
#endif
//--------------------------------------------------------------------------------------------------
      // build single boundary integral stiffness terms
      //
      // remarks: || x || stands for the jump operator ( || x || = x^- - x^+ )
      //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
      //             n    stands for the normal vector on the interface pointing from the positive (+)
      //                  domain into the negative (-) domain (n = n^+ = -n^-)
      //--------------------------------------------------------------------------------------------

#ifdef COMBUST_SIGMA_BASED_FLUXJUMP
      //---------------------------------------------------    |                 |
      // flux (traction/stress) jump (consistency) term RHS  - |  < v >, j_flux  |
      //---------------------------------------------------    |                 |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*jflux(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*jflux(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*jflux(2));

#if 0
      // j_sigma is a 3x3 Matrix -> j_sigma*normal
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(0)*jflux(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(1)*jflux(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*normal(2)*jflux(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(0)*jflux(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(1)*jflux(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*normal(2)*jflux(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(0)*jflux(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(1)*jflux(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*normal(2)*jflux(2,2));
#endif
#endif

#ifdef COMBUST_SIGMA_BASED_CONSISTVISCMEAN
      //-------------------------    |                          |
      // viscous consistency term  + |  || v || , < Dsigma > n  |
      //-------------------------    |                          |
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Sigmaxx>(shp_jump.d0, timefacfac*normal(0), shp_sigma_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Sigmaxy>(shp_jump.d0, timefacfac*normal(1), shp_sigma_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Sigmaxz>(shp_jump.d0, timefacfac*normal(2), shp_sigma_mean.d0);

      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Sigmayx>(shp_jump.d0, timefacfac*normal(0), shp_sigma_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Sigmayy>(shp_jump.d0, timefacfac*normal(1), shp_sigma_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Sigmayz>(shp_jump.d0, timefacfac*normal(2), shp_sigma_mean.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Sigmazx>(shp_jump.d0, timefacfac*normal(0), shp_sigma_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Sigmazy>(shp_jump.d0, timefacfac*normal(1), shp_sigma_mean.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Sigmazz>(shp_jump.d0, timefacfac*normal(2), shp_sigma_mean.d0);

      //   |                           |
      // - |  || v || , < sigma_i > n  |
      //   |                           |

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*normal(0)*sigmamean(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*normal(1)*sigmamean(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*normal(2)*sigmamean(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*normal(0)*sigmamean(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*normal(1)*sigmamean(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*normal(2)*sigmamean(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*normal(0)*sigmamean(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*normal(1)*sigmamean(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*normal(2)*sigmamean(2,2));
#endif

#ifdef COMBUST_SIGMA_BASED_ADJCONSISTVISCJUMP
      //---------------------------------    |                           |
      // viscous adjoint consistency term  + |  < sigma^e > n, || Du ||  |
      //---------------------------------    |                           |
      assembler.template Matrix<XFEM::PHYSICS::Sigmaxx,XFEM::PHYSICS::Velx>(shp_sigma_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Sigmaxy,XFEM::PHYSICS::Velx>(shp_sigma_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Sigmaxz,XFEM::PHYSICS::Velx>(shp_sigma_mean.d0, timefacfac*normal(2), shp_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Sigmayx,XFEM::PHYSICS::Vely>(shp_sigma_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Sigmayy,XFEM::PHYSICS::Vely>(shp_sigma_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Sigmayz,XFEM::PHYSICS::Vely>(shp_sigma_mean.d0, timefacfac*normal(2), shp_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Sigmazx,XFEM::PHYSICS::Velz>(shp_sigma_mean.d0, timefacfac*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Sigmazy,XFEM::PHYSICS::Velz>(shp_sigma_mean.d0, timefacfac*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Sigmazz,XFEM::PHYSICS::Velz>(shp_sigma_mean.d0, timefacfac*normal(2), shp_jump.d0);

      //   |                            |
      // - |  < sigma^e > n, || u_i ||  |
      //   |                            |

      assembler.template Vector<XFEM::PHYSICS::Sigmaxx>(shp_sigma_mean.d0, -timefacfac*normal(0)*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Sigmaxy>(shp_sigma_mean.d0, -timefacfac*normal(1)*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Sigmaxz>(shp_sigma_mean.d0, -timefacfac*normal(2)*veljump(0));

      assembler.template Vector<XFEM::PHYSICS::Sigmayx>(shp_sigma_mean.d0, -timefacfac*normal(0)*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Sigmayy>(shp_sigma_mean.d0, -timefacfac*normal(1)*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Sigmayz>(shp_sigma_mean.d0, -timefacfac*normal(2)*veljump(1));

      assembler.template Vector<XFEM::PHYSICS::Sigmazx>(shp_sigma_mean.d0, -timefacfac*normal(0)*veljump(2));
      assembler.template Vector<XFEM::PHYSICS::Sigmazy>(shp_sigma_mean.d0, -timefacfac*normal(1)*veljump(2));
      assembler.template Vector<XFEM::PHYSICS::Sigmazz>(shp_sigma_mean.d0, -timefacfac*normal(2)*veljump(2));

      //-------------------------------------    |                      |
      // viscous adjoint consistency term RHS  + |  < sigma^e > n, j_u  |
      //-------------------------------------    |                      |
#if 1
      assembler.template Vector<XFEM::PHYSICS::Sigmaxx>(shp_sigma_mean.d0, timefacfac*normal(0)*ju(0));
      assembler.template Vector<XFEM::PHYSICS::Sigmaxy>(shp_sigma_mean.d0, timefacfac*normal(1)*ju(0));
      assembler.template Vector<XFEM::PHYSICS::Sigmaxz>(shp_sigma_mean.d0, timefacfac*normal(2)*ju(0));
#else
      assembler.template Vector<XFEM::PHYSICS::Sigmaxx>(shp_sigma_mean.d0, timefacfac*normal(0)*normal(0)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Sigmaxy>(shp_sigma_mean.d0, timefacfac*normal(1)*normal(0)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Sigmaxz>(shp_sigma_mean.d0, timefacfac*normal(2)*normal(0)*ju_n);
#endif

#if 1
      assembler.template Vector<XFEM::PHYSICS::Sigmayx>(shp_sigma_mean.d0, timefacfac*normal(0)*ju(1));
      assembler.template Vector<XFEM::PHYSICS::Sigmayy>(shp_sigma_mean.d0, timefacfac*normal(1)*ju(1));
      assembler.template Vector<XFEM::PHYSICS::Sigmayz>(shp_sigma_mean.d0, timefacfac*normal(2)*ju(1));
#else
      assembler.template Vector<XFEM::PHYSICS::Sigmayx>(shp_sigma_mean.d0, timefacfac*normal(0)*normal(1)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Sigmayy>(shp_sigma_mean.d0, timefacfac*normal(1)*normal(1)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Sigmayz>(shp_sigma_mean.d0, timefacfac*normal(2)*normal(1)*ju_n);
#endif

#if 1
      assembler.template Vector<XFEM::PHYSICS::Sigmazx>(shp_sigma_mean.d0, timefacfac*normal(0)*ju(2));
      assembler.template Vector<XFEM::PHYSICS::Sigmazy>(shp_sigma_mean.d0, timefacfac*normal(1)*ju(2));
      assembler.template Vector<XFEM::PHYSICS::Sigmazz>(shp_sigma_mean.d0, timefacfac*normal(2)*ju(2));
#else
      assembler.template Vector<XFEM::PHYSICS::Sigmazx>(shp_sigma_mean.d0, timefacfac*normal(0)*normal(2)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Sigmazy>(shp_sigma_mean.d0, timefacfac*normal(1)*normal(2)*ju_n);
      assembler.template Vector<XFEM::PHYSICS::Sigmazz>(shp_sigma_mean.d0, timefacfac*normal(2)*normal(2)*ju_n);
#endif
#endif

#ifdef COMBUST_SIGMA_BASED_NEUMANN_INFLOW
#if 0
      //-----------------    |                   |
      // inflow flux term  - |  v, u_n*|| Du ||  |
      //-----------------    |                   |
      //const double inflowvel = velinflow(0)*normal(0) + velinflow(1)*normal(1) + velinflow(2)*normal(2);
      const double inflowvel = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);

      //   |                    |
      // + |  v, u_n*|| u_i ||  |
      //   |                    |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_inflow.d0, timefacfac*inflowvel*velinflow(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_inflow.d0, timefacfac*inflowvel*velinflow(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_inflow.d0, timefacfac*inflowvel*velinflow(2));
#endif
      //-----------------    |                       |
      // inflow flux term  - |  < v >, u_n*|| Du ||  |
      //-----------------    |                       |
      const double inflowvel = -1.0*(velmean(0)*normal(0) + velmean(1)*normal(1) + velmean(2)*normal(2));
      //const double inflowvel = 1000*(veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2));
//      cout << "inflowvel" << inflowvel << endl;
//      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);

//      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);

      //   |                        |
      // + |  < v >, u_n*|| u_i ||  |
      //   |                        |
//      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel*veljump(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel*veljump(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel*veljump(2));

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, -timefacfac*inflowvel*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, -timefacfac*inflowvel*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, -timefacfac*inflowvel*veljump(2));

//      assembler.template Vector<XFEM::PHYSICS::Velx>(enrfunct_minus, timefacfac*inflowvel*veljump(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(enrfunct_minus, timefacfac*inflowvel*veljump(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(enrfunct_minus, timefacfac*inflowvel*veljump(2));
      //   |                 |
      // + |  < v >, u_n*ju  |
      //   |                 |
//      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*inflowvel*ju(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*inflowvel*ju(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*inflowvel*ju(2));

      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_jump.d0, timefacfac*inflowvel*ju(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_jump.d0, timefacfac*inflowvel*ju(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_jump.d0, timefacfac*inflowvel*ju(2));

//      assembler.template Vector<XFEM::PHYSICS::Velx>(enrfunct_minus, -timefacfac*inflowvel*ju(0));
//      assembler.template Vector<XFEM::PHYSICS::Vely>(enrfunct_minus, -timefacfac*inflowvel*ju(1));
//      assembler.template Vector<XFEM::PHYSICS::Velz>(enrfunct_minus, -timefacfac*inflowvel*ju(2));

#endif

#ifdef COMBUST_SIGMA_BASED_ADVECTIVEFLUX
      //-----------------    |                   |
      // inflow flux term  + |  < v >, || Du ||  |
      //-----------------    |                   |
      const double inflowvel = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);

      //   |                    |
      // - |  < v >, || u_i ||  |
      //   |                    |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, -timefacfac*inflowvel*veljump(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, -timefacfac*inflowvel*veljump(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, -timefacfac*inflowvel*veljump(2));

      //-----------------       |                |
      // inflow flux term RHS + |  < v >, j_u n  |
      //-----------------       |                |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel*ju*normal(0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel*ju*normal(1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel*ju*normal(2));

      //-----------------    |                     |
      // inflow flux term  - |  < q >, || Du || n  |
      //-----------------    |                     |
      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velx>(shp_mean.d0, timefacfac*inflowvel*normal(0), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Vely>(shp_mean.d0, timefacfac*inflowvel*normal(1), shp_jump.d0);
      assembler.template Matrix<XFEM::PHYSICS::Pres,XFEM::PHYSICS::Velz>(shp_mean.d0, timefacfac*inflowvel*normal(2), shp_jump.d0);

      //   |                      |
      // + |  < q >, || u_i || n  |
      //   |                      |
      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, -timefacfac*inflowvel*(veljump(0)*normal(0)
                                                                        +veljump(1)*normal(1)
                                                                        +veljump(2)*normal(2)));

      //-----------------       |                |
      // inflow flux term RHS - |  < q >, j_u n  |
      //-----------------       |                |
      assembler.template Vector<XFEM::PHYSICS::Pres>(shp_mean.d0, timefacfac*inflowvel*ju);
#endif

    } // loop Gaussian points
  } // loop boundary integration cells

  return;
}

/*!
  Calculate boundary integrals in matrix and rhs for stationary combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class V1, class V2, class V3>
void SysmatBoundaryStressDoubleOneSided(
    const DRT::ELEMENTS::Combust3*    ele,           ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*&  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
    const M1&                         evelnp,
    const V1&                         eprenp,
    const V2&                         ephi,
    const M2&                         eepsilon,
    const V3&                         ediscpres,
    Teuchos::RCP<const MAT::Material> material,      ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,      ///< time discretization type
    const double&                     dt,            ///< delta t (time step size)
    const double&                     theta,         ///< factor for one step theta scheme
    const double                      ga_alphaF,
    const double                      ga_alphaM,
    const double                      ga_gamma,
    COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF>& assembler,
    const double                      flamespeed
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

  //--------------------------------------------------------------------------
  // get material parameters for all boundary integration cell of this element
  //--------------------------------------------------------------------------
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

  //--------------------
  // compute jump values
  //--------------------
  // velocity jump value
  const double juplus = 2.0; // -flamespeed*dens_minus*(1/dens_plus - 1/dens_minus);
  const double juminus = -4.0;
  //cout << "ju: " << ju << endl;
  // pressure jump value
  const double jpplus = 0.0; // (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);
  const double jpminus = 0.0;
  //cout << "jp: " << jp << endl;

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;

    // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert(numparamvelx == numparampres, "assumption violation");

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Epsilonxx);
  const size_t numparamdiscpres = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 1);

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
    //----------------------------------------------
    // get Gaussian points for this integration cell
    //----------------------------------------------
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

    //--------------------------------------
    // integration loop over Gaussian points
    //--------------------------------------
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      //---------------------------------------------
      // transform coordinates of this Gaussian point
      //---------------------------------------------
      // coordinates of this integration point in boundary cell coordinates \eta^boundary
      const LINALG::Matrix<2,1> posEtaBoundary(intpoints.qxg[iquad]);

      // coordinates of this integration point in element coordinates \xi^domain
      LINALG::Matrix<nsd,1> posXiDomain;
      GEO::mapEtaBToXiD(*cell, posEtaBoundary, posXiDomain);

      //--------------------------------------------------------------------------------------------
      // compute normal vector (in physical coordinates)
      // remark: for linear boundary integrastion cells this could be done before the loop over all
      //         Gaussian points
      //--------------------------------------------------------------------------------------------
      static LINALG::Matrix<nsd,1> normal(true);
      normal.Clear();
      GEO::computeNormalToSurfaceElement(cell->Shape(), cellXYZDomain, posEtaBoundary, normal);

      //------------------------------------------------
      // evaluate shape functions at this Gaussian point
      //------------------------------------------------
      // shape functions for fields with dofs at nodes
      static LINALG::Matrix<numnode,1> funct;
      funct.Clear();
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      // derivatives of shape functions for fields with dofs at nodes
      static LINALG::Matrix<nsd,numnode> deriv;
      deriv.Clear();
      DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      // shape functions for element stress field
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement,1> funct_stress;
      funct_stress.Clear();
      DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
      // shape functions for discontinuous element pressure field
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement,1> funct_discpres;
      funct_discpres.Clear();
      DRT::UTILS::shape_function_3D(funct_discpres,posXiDomain(0),posXiDomain(1),posXiDomain(2),discpresdistype);

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

      //------------------------------------------------
      // compute Jacobian matrix for mapping to boundary
      //------------------------------------------------
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
        dserror("negative Jacobian determinant detXitoEta: %f for element %d!", detXitoEta, ele->Id());
      }
//#endif

      //----------------------------
      // compute integration factors
      //----------------------------
      // compute spatial integration factor
      const double fac = intpoints.qwgt[iquad]*detXtoXi*detXitoEta;
      //cout << "boundary spatial factor " << fac << endl;
      // compute total (time and spatial) integration factor (and coefficients of single terms)?
      const double timefacfac = timefac * fac;
      //cout << "boundary time and spatial factor " << timefacfac << endl;

      //--------------------------------------------------------------------------------
      // rearrange (enriched) shape functions and derivatives as approximation functions
      //--------------------------------------------------------------------------------
      const std::size_t shpVecSize         = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const std::size_t shpVecSizeStress   = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
      const std::size_t shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

      // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
      static LINALG::Matrix<shpVecSize,1>         enrfunct_plus;
      enrfunct_plus.Clear();
      static LINALG::Matrix<shpVecSize,1>         enrfunct_minus;
      enrfunct_minus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_plus;
      enrfunct_stress_plus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_minus;
      enrfunct_stress_minus.Clear();
      static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres_plus;
      enrfunct_discpres_plus.Clear();
      static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres_minus;
      enrfunct_discpres_minus.Clear();

      static LINALG::Matrix<3,shpVecSize> enrderxy_plus;
      enrderxy_plus.Clear();
      static LINALG::Matrix<3,shpVecSize> enrderxy_minus;
      enrderxy_minus.Clear();
      static LINALG::Matrix<6,shpVecSize> enrderxy2_plus;
      enrderxy2_plus.Clear();
      static LINALG::Matrix<6,shpVecSize> enrderxy2_minus;
      enrderxy2_minus.Clear();

      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_plus, enrderxy_plus, enrderxy2_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(XFEM::PHYSICS::Velx, funct, derxy, derxy2,
          enrfunct_minus, enrderxy_minus, enrderxy2_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Epsilonxx, funct_stress, enrfunct_stress_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::Epsilonxx, funct_stress, enrfunct_stress_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::DiscPres, funct_discpres, enrfunct_discpres_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(XFEM::PHYSICS::DiscPres, funct_discpres, enrfunct_discpres_minus);

      // (enriched) shape functions = approximation functions (P = N * \Psi)
      static XFEM::ApproxFunc<0,shpVecSize> shp_plus;
      static XFEM::ApproxFunc<0,shpVecSize> shp_minus;

      // fill approximation functions
      for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        shp_plus.d0(iparam) = enrfunct_plus(iparam);
        shp_minus.d0(iparam) = enrfunct_minus(iparam);
      }

      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_epsilon_plus;
      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_epsilon_minus;
      for (std::size_t iparam = 0; iparam < numparamepsilonxx; ++iparam)
      {
        shp_epsilon_plus.d0(iparam) = enrfunct_stress_plus(iparam);
        shp_epsilon_minus.d0(iparam) = enrfunct_stress_minus(iparam);
      }

      static XFEM::ApproxFunc<0,shpVecSizeDiscPres> shp_discpres_plus;
      static XFEM::ApproxFunc<0,shpVecSizeDiscPres> shp_discpres_minus;
      for (size_t iparam = 0; iparam < numparamdiscpres; ++iparam)
      {
        shp_discpres_plus.d0(iparam) = enrfunct_discpres_plus(iparam);
        shp_discpres_minus.d0(iparam) = enrfunct_discpres_minus(iparam);
      }

      //---------------------------------------
      // compute data at Gaussian point for rhs
      //---------------------------------------

      // velocity jump
      static LINALG::Matrix<nsd,1> velplus(true);
      velplus.Clear();
      velplus = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_plus.d0, numparamvelx);

      // velocity jump
      static LINALG::Matrix<nsd,1> velminus(true);
      velminus.Clear();
      velminus = XFEM::interpolateVectorFieldToIntPoint(evelnp, shp_minus.d0, numparamvelx);

      // get discontinuous pressure jump
      static double pdiscplus;
      pdiscplus = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pdiscplus += shp_discpres_plus.d0(iparam)*ediscpres(iparam);

      // get discontinuous pressure average
      static double pdiscminus;
      pdiscminus = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pdiscminus += shp_discpres_minus.d0(iparam)*ediscpres(iparam);

      // get viscous stress unknowns
      static LINALG::Matrix<nsd,nsd> epsilonplus;
      epsilonplus.Clear();
      XFEM::fill_tensor(numparamepsilonxx, shp_epsilon_plus.d0, eepsilon, epsilonplus);
      static LINALG::Matrix<nsd,nsd> epsilonminus;
      epsilonminus.Clear();
      XFEM::fill_tensor(numparamepsilonxx, shp_epsilon_minus.d0, eepsilon, epsilonminus);

      //--------------------------------------------------------------------------------------------
      // build single boundary integral stiffness terms
      //
      // remarks: || x || stands for the jump operator ( || x || = x^+ - x^- )
      //           < x >  stands for the mean across the interface ( < x > = 1/2*(x^+ + x^-) )
      //             n    stands for the normal vector on the interface pointing into the burnt
      //                  domain (n = n^- = -n^+)
      //--------------------------------------------------------------------------------------------

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_TAU_PLUS
      //-------------------------    |                      |
      // viscous consistency term  + |  v^+ , Depsilon^+ n  |
      //-------------------------    |                      |
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxx>(shp_plus.d0, timefacfac*normal(0), shp_epsilon_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxy>(shp_plus.d0, timefacfac*normal(1), shp_epsilon_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxz>(shp_plus.d0, timefacfac*normal(2), shp_epsilon_plus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyx>(shp_plus.d0, timefacfac*normal(0), shp_epsilon_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyy>(shp_plus.d0, timefacfac*normal(1), shp_epsilon_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyz>(shp_plus.d0, timefacfac*normal(2), shp_epsilon_plus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzx>(shp_plus.d0, timefacfac*normal(0), shp_epsilon_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzy>(shp_plus.d0, timefacfac*normal(1), shp_epsilon_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzz>(shp_plus.d0, timefacfac*normal(2), shp_epsilon_plus.d0);

      //   |                       |
      // - |  v^+ , epsilon_i^+ n  |
      //   |                       |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_plus.d0, -timefacfac*normal(0)*epsilonplus(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_plus.d0, -timefacfac*normal(1)*epsilonplus(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_plus.d0, -timefacfac*normal(2)*epsilonplus(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_plus.d0, -timefacfac*normal(0)*epsilonplus(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_plus.d0, -timefacfac*normal(1)*epsilonplus(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_plus.d0, -timefacfac*normal(2)*epsilonplus(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_plus.d0, -timefacfac*normal(0)*epsilonplus(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_plus.d0, -timefacfac*normal(1)*epsilonplus(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_plus.d0, -timefacfac*normal(2)*epsilonplus(2,2));
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_TAU_MINUS
      //-------------------------    |                      |
      // viscous consistency term  - |  v^- , Depsilon^- n  |
      //-------------------------    |                      |
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxx>(shp_minus.d0, -timefacfac*normal(0), shp_epsilon_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxy>(shp_minus.d0, -timefacfac*normal(1), shp_epsilon_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::Epsilonxz>(shp_minus.d0, -timefacfac*normal(2), shp_epsilon_minus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyx>(shp_minus.d0, -timefacfac*normal(0), shp_epsilon_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyy>(shp_minus.d0, -timefacfac*normal(1), shp_epsilon_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::Epsilonyz>(shp_minus.d0, -timefacfac*normal(2), shp_epsilon_minus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzx>(shp_minus.d0, -timefacfac*normal(0), shp_epsilon_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzy>(shp_minus.d0, -timefacfac*normal(1), shp_epsilon_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::Epsilonzz>(shp_minus.d0, -timefacfac*normal(2), shp_epsilon_minus.d0);

      //   |                       |
      // + |  v^- , epsilon_i^- n  |
      //   |                       |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_minus.d0, timefacfac*normal(0)*epsilonminus(0,0));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_minus.d0, timefacfac*normal(1)*epsilonminus(0,1));
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_minus.d0, timefacfac*normal(2)*epsilonminus(0,2));

      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_minus.d0, timefacfac*normal(0)*epsilonminus(1,0));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_minus.d0, timefacfac*normal(1)*epsilonminus(1,1));
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_minus.d0, timefacfac*normal(2)*epsilonminus(1,2));

      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_minus.d0, timefacfac*normal(0)*epsilonminus(2,0));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_minus.d0, timefacfac*normal(1)*epsilonminus(2,1));
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_minus.d0, timefacfac*normal(2)*epsilonminus(2,2));
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_PRES_PLUS
#if 0
      //--------------------------    |                   |
      // pressure consistency term  - |  v^+ , Dp^{e+} n  |
      //--------------------------    |                   |
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::DiscPres>(shp_plus.d0, -timefacfac*normal(0), shp_discpres_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::DiscPres>(shp_plus.d0, -timefacfac*normal(1), shp_discpres_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::DiscPres>(shp_plus.d0, -timefacfac*normal(2), shp_discpres_plus.d0);

      //   |                    |
      // + |  v^+ , p^{e+}_i n  |
      //   |                    |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_plus.d0, timefacfac*normal(0)*pdiscplus);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_plus.d0, timefacfac*normal(1)*pdiscplus);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_plus.d0, timefacfac*normal(2)*pdiscplus);
#endif
#if 1
      //------------------------------------------    |                 |
      // pressure/flux jump (consistency) term RHS  + |  v^+ , j_p^+ n  |
      //------------------------------------------    |                 |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_plus.d0, timefacfac*normal(0)*jpplus);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_plus.d0, timefacfac*normal(1)*jpplus);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_plus.d0, timefacfac*normal(2)*jpplus);
#endif
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_CONSIST_PRES_MINUS
#if 0
      //--------------------------    |                   |
      // pressure consistency term  + |  v^- , Dp^{e-} n  |
      //--------------------------    |                   |
      assembler.template Matrix<XFEM::PHYSICS::Velx,XFEM::PHYSICS::DiscPres>(shp_minus.d0, timefacfac*normal(0), shp_discpres_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Vely,XFEM::PHYSICS::DiscPres>(shp_minus.d0, timefacfac*normal(1), shp_discpres_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Velz,XFEM::PHYSICS::DiscPres>(shp_minus.d0, timefacfac*normal(2), shp_discpres_minus.d0);

      //   |                    |
      // - |  v^- , p^{e-}_i n  |
      //   |                    |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_minus.d0, -timefacfac*normal(0)*pdiscminus);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_minus.d0, -timefacfac*normal(1)*pdiscminus);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_minus.d0, -timefacfac*normal(2)*pdiscminus);
#endif
#if 1
      //------------------------------------------    |                 |
      // pressure/flux jump (consistency) term RHS  + |  v^- , j_p^- n  |
      //------------------------------------------    |                 |
      assembler.template Vector<XFEM::PHYSICS::Velx>(shp_minus.d0, timefacfac*normal(0)*jpminus);
      assembler.template Vector<XFEM::PHYSICS::Vely>(shp_minus.d0, timefacfac*normal(1)*jpminus);
      assembler.template Vector<XFEM::PHYSICS::Velz>(shp_minus.d0, timefacfac*normal(2)*jpminus);
#endif
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_TAU_DIRICH_PLUS
      //---------------------------------    |                        |
      // viscous adjoint consistency term  - |  epsilon^{e+} n, Du^+  |
      //---------------------------------    |                        |
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Velx>(shp_epsilon_plus.d0, -timefacfac*normal(0), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Velx>(shp_epsilon_plus.d0, -timefacfac*normal(1), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velx>(shp_epsilon_plus.d0, -timefacfac*normal(2), shp_plus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Epsilonyx,XFEM::PHYSICS::Vely>(shp_epsilon_plus.d0, -timefacfac*normal(0), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Vely>(shp_epsilon_plus.d0, -timefacfac*normal(1), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Vely>(shp_epsilon_plus.d0, -timefacfac*normal(2), shp_plus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Epsilonzx,XFEM::PHYSICS::Velz>(shp_epsilon_plus.d0, -timefacfac*normal(0), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonzy,XFEM::PHYSICS::Velz>(shp_epsilon_plus.d0, -timefacfac*normal(1), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Velz>(shp_epsilon_plus.d0, -timefacfac*normal(2), shp_plus.d0);

      //   |                         |
      // + |  epsilon^{e+} n, u_i^+  |
      //   |                         |
      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon_plus.d0, timefacfac*normal(0)*velplus(0,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon_plus.d0, timefacfac*normal(1)*velplus(0,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon_plus.d0, timefacfac*normal(2)*velplus(0,0));

      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon_plus.d0, timefacfac*normal(0)*velplus(1,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon_plus.d0, timefacfac*normal(1)*velplus(1,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon_plus.d0, timefacfac*normal(2)*velplus(1,0));

      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon_plus.d0, timefacfac*normal(0)*velplus(2,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon_plus.d0, timefacfac*normal(1)*velplus(2,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon_plus.d0, timefacfac*normal(2)*velplus(2,0));

      //-------------------------------------    |                         |
      // viscous adjoint consistency term RHS  - |  epsilon^{e+} n, j_u^+  |
      //-------------------------------------    |                         |

      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon_plus.d0, -timefacfac*normal(0)*normal(0)*juplus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon_plus.d0, -timefacfac*normal(1)*normal(0)*juplus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon_plus.d0, -timefacfac*normal(2)*normal(0)*juplus);

      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon_plus.d0, -timefacfac*normal(0)*normal(1)*juplus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon_plus.d0, -timefacfac*normal(1)*normal(1)*juplus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon_plus.d0, -timefacfac*normal(2)*normal(1)*juplus);

      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon_plus.d0, -timefacfac*normal(0)*normal(2)*juplus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon_plus.d0, -timefacfac*normal(1)*normal(2)*juplus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon_plus.d0, -timefacfac*normal(2)*normal(2)*juplus);
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_TAU_DIRICH_MINUS
      //---------------------------------    |                        |
      // viscous adjoint consistency term  + |  epsilon^{e-} n, Du^-  |
      //---------------------------------    |                        |
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxx,XFEM::PHYSICS::Velx>(shp_epsilon_minus.d0, timefacfac*normal(0), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxy,XFEM::PHYSICS::Velx>(shp_epsilon_minus.d0, timefacfac*normal(1), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonxz,XFEM::PHYSICS::Velx>(shp_epsilon_minus.d0, timefacfac*normal(2), shp_minus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Epsilonyx,XFEM::PHYSICS::Vely>(shp_epsilon_minus.d0, timefacfac*normal(0), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonyy,XFEM::PHYSICS::Vely>(shp_epsilon_minus.d0, timefacfac*normal(1), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonyz,XFEM::PHYSICS::Vely>(shp_epsilon_minus.d0, timefacfac*normal(2), shp_minus.d0);

      assembler.template Matrix<XFEM::PHYSICS::Epsilonzx,XFEM::PHYSICS::Velz>(shp_epsilon_minus.d0, timefacfac*normal(0), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonzy,XFEM::PHYSICS::Velz>(shp_epsilon_minus.d0, timefacfac*normal(1), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::Epsilonzz,XFEM::PHYSICS::Velz>(shp_epsilon_minus.d0, timefacfac*normal(2), shp_minus.d0);

      //   |                         |
      // - |  epsilon^{e-} n, u_i^-  |
      //   |                         |
      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon_minus.d0, -timefacfac*normal(0)*velminus(0,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon_minus.d0, -timefacfac*normal(1)*velminus(0,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon_minus.d0, -timefacfac*normal(2)*velminus(0,0));

      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon_minus.d0, -timefacfac*normal(0)*velminus(1,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon_minus.d0, -timefacfac*normal(1)*velminus(1,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon_minus.d0, -timefacfac*normal(2)*velminus(1,0));

      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon_minus.d0, -timefacfac*normal(0)*velminus(2,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon_minus.d0, -timefacfac*normal(1)*velminus(2,0));
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon_minus.d0, -timefacfac*normal(2)*velminus(2,0));

      //-------------------------------------    |                         |
      // viscous adjoint consistency term RHS  + |  epsilon^{e-} n, j_u^-  |
      //-------------------------------------    |                         |

      assembler.template Vector<XFEM::PHYSICS::Epsilonxx>(shp_epsilon_minus.d0, timefacfac*normal(0)*normal(0)*juminus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonxy>(shp_epsilon_minus.d0, timefacfac*normal(1)*normal(0)*juminus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonxz>(shp_epsilon_minus.d0, timefacfac*normal(2)*normal(0)*juminus);

      assembler.template Vector<XFEM::PHYSICS::Epsilonyx>(shp_epsilon_minus.d0, timefacfac*normal(0)*normal(1)*juminus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonyy>(shp_epsilon_minus.d0, timefacfac*normal(1)*normal(1)*juminus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonyz>(shp_epsilon_minus.d0, timefacfac*normal(2)*normal(1)*juminus);

      assembler.template Vector<XFEM::PHYSICS::Epsilonzx>(shp_epsilon_minus.d0, timefacfac*normal(0)*normal(2)*juminus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonzy>(shp_epsilon_minus.d0, timefacfac*normal(1)*normal(2)*juminus);
      assembler.template Vector<XFEM::PHYSICS::Epsilonzz>(shp_epsilon_minus.d0, timefacfac*normal(2)*normal(2)*juminus);
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_PRES_DIRICH_PLUS
      //----------------------------------    |                  |
      // pressure adjoint consistency term  + |  q^{e+} n, Du^+  |
      //----------------------------------    |                  |
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velx>(shp_discpres_plus.d0, timefacfac*normal(0), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Vely>(shp_discpres_plus.d0, timefacfac*normal(1), shp_plus.d0);
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velz>(shp_discpres_plus.d0, timefacfac*normal(2), shp_plus.d0);
      //   |                   |
      // - |  q^{e+} n, u_i^+  |
      //   |                   |
      assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_plus.d0, -timefacfac*(velplus(0,0)*normal(0)
                                                                            +velplus(1,0)*normal(1)
                                                                            +velplus(2,0)*normal(2)));
      //--------------------------------------    |                   |
      // pressure adjoint consistency term RHS  + |  q^{e+} n, j_u^+  |
      //--------------------------------------    |                   |
      assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_plus.d0, timefacfac*juplus);
#endif

#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED_PRES_DIRICH_MINUS
      //----------------------------------    |                  |
      // pressure adjoint consistency term  - |  q^{e+} n, Du^+  |
      //----------------------------------    |                  |
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velx>(shp_discpres_minus.d0, -timefacfac*normal(0), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Vely>(shp_discpres_minus.d0, -timefacfac*normal(1), shp_minus.d0);
      assembler.template Matrix<XFEM::PHYSICS::DiscPres,XFEM::PHYSICS::Velz>(shp_discpres_minus.d0, -timefacfac*normal(2), shp_minus.d0);
      //   |                   |
      // + |  q^{e+} n, u_i^+  |
      //   |                   |
      assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_minus.d0, timefacfac*(velminus(0,0)*normal(0)
                                                                            +velminus(1,0)*normal(1)
                                                                            +velminus(2,0)*normal(2)));
      //--------------------------------------    |                   |
      // pressure adjoint consistency term RHS  - |  q^{e+} n, j_u^+  |
      //--------------------------------------    |                   |
      assembler.template Vector<XFEM::PHYSICS::DiscPres>(shp_discpres_minus.d0, -timefacfac*juminus);
#endif
    } // loop Gaussian points
  } // loop boundary integration cells

  return;
}

} // namespace COMBUST

#endif
