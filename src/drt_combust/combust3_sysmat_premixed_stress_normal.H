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

#ifndef COMBUST3_SYSMAT_PREMIXED_STRESS_NORMAL_H
#define COMBUST3_SYSMAT_PREMIXED_STRESS_NORMAL_H


#include "combust3_interpolation.H"
#include "combust_defines.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_normals.H"


#ifdef COMBUST_NORMAL_ENRICHMENT

namespace COMBUST
{

template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          size_t NUMDOF,
          size_t shpVecSizeVel,
          size_t shpVecSizePres,
          size_t shpVecSizeStress,
          size_t shpVecSizeDiscPres>
void BuildDomainIntegralsStressNormal(
    COMBUST::LocalAssembler<DISTYPE,ASSTYPE,NUMDOF>& assembler,
    const XFEM::ApproxFuncNormalVector<2,shpVecSizeVel>& shp,    // shape functions and derivatives
    const XFEM::ApproxFunc<1,shpVecSizePres>& shp_pres,         // shape functions and derivatives
    const XFEM::ApproxFunc<0,shpVecSizeStress>& shp_epsilon,    // shape functions for stresses
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
  // remark: the symmetry properties are used here (Epsilonxy = Epsilonyx, Epsilonxz = Epsilonzx, Epsilonyz = Epsilonzy)
  assembler.template Matrix<Epsilonxx,Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac,     shp_epsilon.d0);
  assembler.template Matrix<Epsilonxy,Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*2.0, shp_epsilon.d0);
  assembler.template Matrix<Epsilonxz,Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*2.0, shp_epsilon.d0);
  assembler.template Matrix<Epsilonyy,Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac,     shp_epsilon.d0);
  assembler.template Matrix<Epsilonyz,Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*2.0, shp_epsilon.d0);
  assembler.template Matrix<Epsilonzz,Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac,     shp_epsilon.d0);

  assembler.template Vector<Epsilonxx>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,0)    );
  assembler.template Vector<Epsilonxy>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,1)*2.0);
  assembler.template Vector<Epsilonxz>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,2)*2.0);
  assembler.template Vector<Epsilonyy>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,1)    );
  assembler.template Vector<Epsilonyz>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,2)*2.0);
  assembler.template Vector<Epsilonzz>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,2)    );

  // formulation without using symmetry of stress tensor
  //  assembler.template Matrix<Epsilonxx,Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonxy,Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonxz,Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonyx,Epsilonyx>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonyy,Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonyz,Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonzx,Epsilonzx>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonzy,Epsilonzy>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);
  //  assembler.template Matrix<Epsilonzz,Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac, shp_epsilon.d0);

  //  assembler.template Vector<Epsilonxx>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,0));
  //  assembler.template Vector<Epsilonxy>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,1));
  //  assembler.template Vector<Epsilonxz>(shp_epsilon.d0, twomu*timefacfac*epsilon(0,2));
  //  assembler.template Vector<Epsilonyx>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,0));
  //  assembler.template Vector<Epsilonyy>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,1));
  //  assembler.template Vector<Epsilonyz>(shp_epsilon.d0, twomu*timefacfac*epsilon(1,2));
  //  assembler.template Vector<Epsilonzx>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,0));
  //  assembler.template Vector<Epsilonzy>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,1));
  //  assembler.template Vector<Epsilonzz>(shp_epsilon.d0, twomu*timefacfac*epsilon(2,2));

  // epsilon-epsilon coupling, continuous epsilon term
  /*               /                      \
                  | virt epsilon , 2\mu*Du |
                   \                      /
   */
  // remark: the symmetry properties are used here (Epsilonxy = Epsilonyx, Epsilonxz = Epsilonzx, Epsilonyz = Epsilonzy)
  assembler.template MatrixSN<Epsilonxx,Velx,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.velx.dx);
  assembler.template MatrixSN<Epsilonxy,Velx,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.velx.dy);
  assembler.template MatrixSN<Epsilonxy,Vely,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.vely.dx);
  assembler.template MatrixSN<Epsilonxz,Velx,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.velx.dz);
  assembler.template MatrixSN<Epsilonxz,Velz,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.velz.dx);
  assembler.template MatrixSN<Epsilonyy,Vely,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.vely.dy);
  assembler.template MatrixSN<Epsilonyz,Vely,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.vely.dz);
  assembler.template MatrixSN<Epsilonyz,Velz,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.velz.dy);
  assembler.template MatrixSN<Epsilonzz,Velz,shpVecSizeStress>(shp_epsilon.d0, twomu*timefacfac , shp.velz.dz);

  assembler.template Vector<Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac*vderxy(0,0));
  assembler.template Vector<Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*(vderxy(0,1) + vderxy(1,0)));
  assembler.template Vector<Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*(vderxy(0,2) + vderxy(2,0)));
  assembler.template Vector<Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac*vderxy(1,1));
  assembler.template Vector<Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*(vderxy(1,2) + vderxy(2,1)));
  assembler.template Vector<Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac*vderxy(2,2));

  // formulation without using symmetry of stress tensor
  // assembler.template Matrix<Epsilonxx,Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<Epsilonxx,Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<Epsilonxy,Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<Epsilonxy,Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<Epsilonxz,Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<Epsilonxz,Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  //
  // assembler.template Matrix<Epsilonyx,Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<Epsilonyx,Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<Epsilonyy,Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<Epsilonyy,Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<Epsilonyz,Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<Epsilonyz,Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  //
  // assembler.template Matrix<Epsilonzx,Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dx);
  // assembler.template Matrix<Epsilonzx,Velx>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<Epsilonzy,Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dy);
  // assembler.template Matrix<Epsilonzy,Vely>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<Epsilonzz,Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  // assembler.template Matrix<Epsilonzz,Velz>(shp_epsilon.d0, twomu*timefacfac*0.5, shp.dz);
  //
  // assembler.template Vector<Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,0));
  // assembler.template Vector<Epsilonxx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,0));
  // assembler.template Vector<Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,1));
  // assembler.template Vector<Epsilonxy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,0));
  // assembler.template Vector<Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,2));
  // assembler.template Vector<Epsilonxz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,0));
  //
  // assembler.template Vector<Epsilonyx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,0));
  // assembler.template Vector<Epsilonyx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,1));
  // assembler.template Vector<Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,1));
  // assembler.template Vector<Epsilonyy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,1));
  // assembler.template Vector<Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,2));
  // assembler.template Vector<Epsilonyz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,1));
  //
  // assembler.template Vector<Epsilonzx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,0));
  // assembler.template Vector<Epsilonzx>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(0,2));
  // assembler.template Vector<Epsilonzy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,1));
  // assembler.template Vector<Epsilonzy>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(1,2));
  // assembler.template Vector<Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,2));
  // assembler.template Vector<Epsilonzz>(shp_epsilon.d0, -twomu*timefacfac*0.5*vderxy(2,2));
#endif

#ifdef COMBUST_EPSPRES_BASED_PRESSURE_COUPLING
  // pressure-pressure coupling, discontinuous pressure term
  /*                 /                      \
                    | tr(virt epsilon) , p^e |
                     \                      /
   */
  assembler.template Matrix<Epsilonxx,DiscPres>(shp_epsilon.d0, timefacfac, shp_discpres.d0);
  assembler.template Matrix<Epsilonyy,DiscPres>(shp_epsilon.d0, timefacfac, shp_discpres.d0);
  assembler.template Matrix<Epsilonzz,DiscPres>(shp_epsilon.d0, timefacfac, shp_discpres.d0);

  assembler.template Vector<Epsilonxx>(shp_epsilon.d0, -timefacfac*discpres);
  assembler.template Vector<Epsilonyy>(shp_epsilon.d0, -timefacfac*discpres);
  assembler.template Vector<Epsilonzz>(shp_epsilon.d0, -timefacfac*discpres);

  // pressure-pressure coupling, continuous pressure term
  /*                 /                    \
                  - | tr(virt epsilon) , p |
                     \                    /
   */
  assembler.template Matrix<Epsilonxx,Pres>(shp_epsilon.d0, -timefacfac, shp.d0);
  assembler.template Matrix<Epsilonyy,Pres>(shp_epsilon.d0, -timefacfac, shp.d0);
  assembler.template Matrix<Epsilonzz,Pres>(shp_epsilon.d0, -timefacfac, shp.d0);

  assembler.template Vector<Epsilonxx>(shp_epsilon.d0, timefacfac*pres);
  assembler.template Vector<Epsilonyy>(shp_epsilon.d0, timefacfac*pres);
  assembler.template Vector<Epsilonzz>(shp_epsilon.d0, timefacfac*pres);
#endif

#ifdef COMBUST_EPSPRES_BASED_ZERO_TRACE
  // enforce zero trace of discontinuous epsilon tensor (solenoidality of auxiliar field)
  /*                 /               \
                  - | q^e , tr(eps^e) |
                     \               /
   */
  assembler.template Matrix<DiscPres,Epsilonxx>(shp_discpres.d0, -timefacfac, shp_epsilon.d0);
  assembler.template Matrix<DiscPres,Epsilonyy>(shp_discpres.d0, -timefacfac, shp_epsilon.d0);
  assembler.template Matrix<DiscPres,Epsilonzz>(shp_discpres.d0, -timefacfac, shp_epsilon.d0);

  // compute trace of discontinuous epsilon tensor
  const double trace_epsilon = (epsilon(0,0) + epsilon(1,1) + epsilon(2,2));
  assembler.template Vector<DiscPres>(shp_discpres.d0, timefacfac*trace_epsilon);

  // enforce zero trace of continuous epsilon tensor (solenoidality of standard field)
  /*                 /                \
                    | q^e , tr(eps(u)) |
                     \                /
  */
  // new term; seems to have effect on third diagonal matrix entry
  assembler.template MatrixSN<DiscPres,Velx>(shp_discpres.d0, timefacfac, shp.velx.dx);
  assembler.template MatrixSN<DiscPres,Vely>(shp_discpres.d0, timefacfac, shp.vely.dy);
  assembler.template MatrixSN<DiscPres,Velz>(shp_discpres.d0, timefacfac, shp.velz.dz);

  const double trace = vderxy(0,0) + vderxy(1,1) + vderxy(2,2);
  assembler.template Vector<DiscPres>(shp_discpres.d0, -timefacfac*trace);
#endif

#ifdef  COMBUST_EPSPRES_BASED_FORTHTERM
  const double viscfac = 3.0/(2.0*dynvisc);
  // TODO the function of this forth term is not clear yet
  /*                 /                  \
                  - | q^e , 3/(2\mu) p^e |
                     \                  /
   */
  assembler.template Matrix<DiscPres,DiscPres>(shp_discpres.d0, viscfac*timefacfac, shp_discpres.d0);
  assembler.template Vector<DiscPres>(shp_discpres.d0, -viscfac*timefacfac*discpres);
  // new term
  /*                 /                \
                    | q^e , 3/(2\mu) p |
                     \                /
   */
  assembler.template Matrix<DiscPres,Pres>(shp_discpres.d0, -viscfac*timefacfac, shp_pres.d0);
  assembler.template Vector<DiscPres>(shp_discpres.d0, viscfac*timefacfac*pres);
#endif
}

/*!
  Calculate domain integrals in matrix and rhs for combustion problem formulation
 */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          int NUMDOF,
          class M1, class M2, class M3, class V1, class V2, class V3>
void SysmatDomainStressNormal(
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
    const M2&                           egradphi,     ///<
    const M3&                           eepsilon,     ///< element stresses
    const V3&                           ediscpres,    ///<
    Teuchos::RCP<const MAT::Material>   material,     ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,     ///< time discretization scheme
    const double                        dt,           ///< delta t (time step size)
    const double                        theta,        ///< parameter for one step theta scheme
    const bool                          newton,       ///< full Newton or fixed-point-like
    const bool                          pstab,        ///< flag for stabilization
    const bool                          supg,         ///< flag for stabilization
    const bool                          graddiv,        ///< flag for stabilization
    const INPAR::FLUID::TauType         tautype,      ///< stabilization parameter definition
    const bool                          instationary, ///< switch between stationary and instationary formulation
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
#ifdef COMBUST_NORMAL_ENRICHMENT
  const size_t numparamveln = XFEM::NumParam<0,ASSTYPE>::get(dofman, XFEM::PHYSICS::Veln);
#endif
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Epsilonxx);
//cout << "standard element: numparamepsilonxx " << numparamepsilonxx << endl;
  //cout << "epsilonele_unknowns_present " << epsilonele_unknowns_present << endl;
  const size_t numparamdiscpres = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 1);
//cout << "standard element: numparamdiscpres " << numparamdiscpres << endl;
//cout << "discpres_unknowns_present " <<  discpres_unknowns_present << endl;
  const size_t numparamsigmaxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Sigmaxx);
//cout << "sigmaele_unknowns_present " << sigmaele_unknowns_present << endl;

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
#ifndef COMBUST_NORMAL_ENRICHMENT
      const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
#else
      const size_t shpVecSizeVel    = COMBUST::SizeFacVel<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const size_t shpVecSizePres   = COMBUST::SizeFacPres<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
#endif
      const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
      const size_t shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

      // (enriched) shape functions = approximation functions (P = N * \Psi)
#ifdef COMBUST_NORMAL_ENRICHMENT
      // second derivatives needed for velocity
      static XFEM::ApproxFuncNormalVector<2,shpVecSizeVel> shp(false);
      shp.Clear();
      static XFEM::ApproxFunc<1,shpVecSizePres> shp_pres;
#else
      static XFEM::ApproxFunc<2,shpVecSize> shp;
#endif
      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_epsilon;
      static XFEM::ApproxFunc<0,shpVecSizeStress> shp_sigma;
      static XFEM::ApproxFunc<0,shpVecSizeDiscPres> shp_discpres;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // temporary arrays holding enriched shape functions (N * \Psi)
        static LINALG::Matrix<shpVecSizePres,1> enr_funct_pres;
        enr_funct_pres.Clear();
        static LINALG::Matrix<3,shpVecSizePres> enr_derxy_pres;
        enr_derxy_pres.Clear();
        static LINALG::Matrix<6,shpVecSizePres> enr_derxy2_pres;
        enr_derxy2_pres.Clear();

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            Pres,
            funct,
            derxy,
            derxy2,
            enr_funct_pres,
            enr_derxy_pres,
            enr_derxy2_pres);

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparampres; ++iparam)
        {
          shp_pres.d0(iparam) = enr_funct_pres(iparam);
          shp_pres.dx(iparam) = enr_derxy_pres(0,iparam);
          shp_pres.dy(iparam) = enr_derxy_pres(1,iparam);
          shp_pres.dz(iparam) = enr_derxy_pres(2,iparam);
        }

        //----------------------------------------
        // evaluate normal field at Gaussian point
        //----------------------------------------

        // fill approximation functions for XFEM
        for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
        {
          shp.velx.d0.s(iparam) = funct(iparam);
          shp.vely.d0.s(iparam) = funct(iparam);
          shp.velz.d0.s(iparam) = funct(iparam);

          shp.velx.dx.s(iparam) = derxy(0,iparam);
          shp.vely.dx.s(iparam) = derxy(0,iparam);
          shp.velz.dx.s(iparam) = derxy(0,iparam);

          shp.velx.dy.s(iparam) = derxy(1,iparam);
          shp.vely.dy.s(iparam) = derxy(1,iparam);
          shp.velz.dy.s(iparam) = derxy(1,iparam);

          shp.velx.dz.s(iparam) = derxy(2,iparam);
          shp.vely.dz.s(iparam) = derxy(2,iparam);
          shp.velz.dz.s(iparam) = derxy(2,iparam);

          shp.velx.dxdx.s(iparam) = derxy2(0,iparam);
          shp.vely.dxdx.s(iparam) = derxy2(0,iparam);
          shp.velz.dxdx.s(iparam) = derxy2(0,iparam);

          shp.velx.dxdy.s(iparam) = derxy2(3,iparam);
          shp.vely.dxdy.s(iparam) = derxy2(3,iparam);
          shp.velz.dxdy.s(iparam) = derxy2(3,iparam);

          shp.velx.dxdz.s(iparam) = derxy2(4,iparam);
          shp.vely.dxdz.s(iparam) = derxy2(4,iparam);
          shp.velz.dxdz.s(iparam) = derxy2(4,iparam);

          shp.velx.dydx.s(iparam) = shp.velx.dxdy.s(iparam);
          shp.vely.dydx.s(iparam) = shp.vely.dxdy.s(iparam);
          shp.velz.dydx.s(iparam) = shp.velz.dxdy.s(iparam);

          shp.velx.dydy.s(iparam) = derxy2(1,iparam);
          shp.vely.dydy.s(iparam) = derxy2(1,iparam);
          shp.velz.dydy.s(iparam) = derxy2(1,iparam);

          shp.velx.dydz.s(iparam) = derxy2(5,iparam);
          shp.vely.dydz.s(iparam) = derxy2(5,iparam);
          shp.velz.dydz.s(iparam) = derxy2(5,iparam);

          shp.velx.dzdx.s(iparam) = shp.velx.dxdz.s(iparam);
          shp.vely.dzdx.s(iparam) = shp.vely.dxdz.s(iparam);
          shp.velz.dzdx.s(iparam) = shp.velz.dxdz.s(iparam);

          shp.velx.dzdy.s(iparam) = shp.velx.dydz.s(iparam);
          shp.vely.dzdy.s(iparam) = shp.vely.dydz.s(iparam);
          shp.velz.dzdy.s(iparam) = shp.velz.dydz.s(iparam);

          shp.velx.dzdz.s(iparam) = derxy2(2,iparam);
          shp.vely.dzdz.s(iparam) = derxy2(2,iparam);
          shp.velz.dzdz.s(iparam) = derxy2(2,iparam);
        }

#ifdef COLLAPSE_FLAME
        // TODO normal vector modified
        LINALG::Matrix<3,1> normal(true);
        for (unsigned i=0;i<numnode;i++)
        {
          normal(0) += funct(i)*xyze(0,i);
          normal(1) += funct(i)*xyze(1,i);
        }
//        // get coordinates of cell vertices
//        const LINALG::SerialDenseMatrix& nodalPosXYZ(cell->CellNodalPosXYZ());
//        normal(0)= nodalPosXYZ(0,ivertex);
//        normal(1)= nodalPosXYZ(1,ivertex);
//        normal(2)=0.0;

        normal.Scale(-1.0);
        const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
        if (norm == 0.0) dserror("norm of normal vector is zero!");
        normal.Scale(1.0/norm);

        LINALG::Matrix<3,1> point(true);
        for (unsigned inode=0;inode<numnode;++inode)
        {
          point(0) += funct(inode)*xyze(0,inode);
          point(1) += funct(inode)*xyze(1,inode);
        }
#endif

        // shape functions and derivatives for nodal parameters (dofs)
        enrvals.ComputeNormalShapeFunction(
            funct,
            derxy,
            derxy2,
            egradphi,
#ifdef COLLAPSE_FLAME
            normal,
            point,
#endif
            shp);

//cout << std::setprecision(16) << std::scientific;
//cout << funct << endl;
//cout << derxy << endl;
//cout << derxy2 << endl;
//
        if (epsilonele_unknowns_present)
        {
          // shape functions for element dofs
          static LINALG::Matrix<shpVecSizeStress,1> enrfunct_stress;
          enrfunct_stress.Clear();
          enrvals.ComputeShiftedEnrichedElementShapefunction(Epsilonxx, funct_stress, enrfunct_stress);

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
          enrvals.ComputeShiftedEnrichedElementShapefunction(Sigmaxx, funct_stress, enrfunct_stress);

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

          enrvals.ComputeShiftedEnrichedElementShapefunction(DiscPres, funct_discpres, enrfunct_discpres);

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
          shp.velx.d0.s(iparam) = funct(iparam);
          shp.vely.d0.s(iparam) = funct(iparam);
          shp.velz.d0.s(iparam) = funct(iparam);

          shp.velx.dx.s(iparam) = derxy(0,iparam);
          shp.vely.dx.s(iparam) = derxy(0,iparam);
          shp.velz.dx.s(iparam) = derxy(0,iparam);

          shp.velx.dy.s(iparam) = derxy(1,iparam);
          shp.vely.dy.s(iparam) = derxy(1,iparam);
          shp.velz.dy.s(iparam) = derxy(1,iparam);

          shp.velx.dz.s(iparam) = derxy(2,iparam);
          shp.vely.dz.s(iparam) = derxy(2,iparam);
          shp.velz.dz.s(iparam) = derxy(2,iparam);

          shp.velx.dxdx.s(iparam) = derxy2(0,iparam);
          shp.vely.dxdx.s(iparam) = derxy2(0,iparam);
          shp.velz.dxdx.s(iparam) = derxy2(0,iparam);

          shp.velx.dxdy.s(iparam) = derxy2(3,iparam);
          shp.vely.dxdy.s(iparam) = derxy2(3,iparam);
          shp.velz.dxdy.s(iparam) = derxy2(3,iparam);

          shp.velx.dxdz.s(iparam) = derxy2(4,iparam);
          shp.vely.dxdz.s(iparam) = derxy2(4,iparam);
          shp.velz.dxdz.s(iparam) = derxy2(4,iparam);

          shp.velx.dydx.s(iparam) = shp.velx.dxdy.s(iparam);
          shp.vely.dydx.s(iparam) = shp.vely.dxdy.s(iparam);
          shp.velz.dydx.s(iparam) = shp.velz.dxdy.s(iparam);

          shp.velx.dydy.s(iparam) = derxy2(1,iparam);
          shp.vely.dydy.s(iparam) = derxy2(1,iparam);
          shp.velz.dydy.s(iparam) = derxy2(1,iparam);

          shp.velx.dydz.s(iparam) = derxy2(5,iparam);
          shp.vely.dydz.s(iparam) = derxy2(5,iparam);
          shp.velz.dydz.s(iparam) = derxy2(5,iparam);

          shp.velx.dzdx.s(iparam) = shp.velx.dxdz.s(iparam);
          shp.vely.dzdx.s(iparam) = shp.vely.dxdz.s(iparam);
          shp.velz.dzdx.s(iparam) = shp.velz.dxdz.s(iparam);

          shp.velx.dzdy.s(iparam) = shp.velx.dydz.s(iparam);
          shp.vely.dzdy.s(iparam) = shp.vely.dydz.s(iparam);
          shp.velz.dzdy.s(iparam) = shp.velz.dydz.s(iparam);

          shp.velx.dzdz.s(iparam) = derxy2(2,iparam);
          shp.vely.dzdz.s(iparam) = derxy2(2,iparam);
          shp.velz.dzdz.s(iparam) = derxy2(2,iparam);

          shp_pres.d0(iparam) = funct(iparam);
          shp_pres.dx(iparam) = derxy(0,iparam);
          shp_pres.dy(iparam) = derxy(1,iparam);
          shp_pres.dz(iparam) = derxy(2,iparam);
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
      // TODO ist das mit der Beschleunigung richtig so?
//cout << "evelnp "<< evelnp << endl;
//cout << "gpveln "<< gpveln << endl;
//cout << "gpvelnm "<< gpvelnm << endl;
//cout << "gpaccn "<< gpaccn << endl;
      const LINALG::Matrix<nsd,1> gpvelnp = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, shp, numparamvelx, ele, dofman);
      const LINALG::Matrix<nsd,1> gpveln  = XFEM::interpolateVectorFieldToIntPointNormal(eveln , shp, numparamvelx, ele, dofman);
      const LINALG::Matrix<nsd,1> gpvelnm = XFEM::interpolateVectorFieldToIntPointNormal(evelnm, shp, numparamvelx, ele, dofman);
      const LINALG::Matrix<nsd,1> gpaccn  = XFEM::interpolateVectorFieldToIntPointNormal(eaccn , shp, numparamvelx, ele, dofman);

//cout << "gpvelnp "<< gpvelnp << endl;
//cout << "gpveln "<< gpveln << endl;
//cout << "gpvelnm "<< gpvelnm << endl;
//cout << "gpaccn "<< gpaccn << endl;

      //----------------------------------------------
      // get velocity derivatives at integration point
      //----------------------------------------------
      static LINALG::Matrix<3,3> vderxy;
      vderxy.Clear();
      {
        const int* nodeids = ele->NodeIds();
        std::size_t velncounter = 0;

        // vderxy = enr_derxy(j,k)*evelnp(i,k);
        for (size_t inode = 0; inode < numparamvelx; ++inode)
        {
          // standard shape functions are identical for all vector components
          // e.g. shp.velx.dx.s == shp.vely.dx.s == shp.velz.dx.s
          vderxy(0,0) += evelnp(0,inode)*shp.velx.dx.s(inode);
          vderxy(0,1) += evelnp(0,inode)*shp.velx.dy.s(inode);
          vderxy(0,2) += evelnp(0,inode)*shp.velx.dz.s(inode);

          vderxy(1,0) += evelnp(1,inode)*shp.vely.dx.s(inode);
          vderxy(1,1) += evelnp(1,inode)*shp.vely.dy.s(inode);
          vderxy(1,2) += evelnp(1,inode)*shp.vely.dz.s(inode);

          vderxy(2,0) += evelnp(2,inode)*shp.velz.dx.s(inode);
          vderxy(2,1) += evelnp(2,inode)*shp.velz.dy.s(inode);
          vderxy(2,2) += evelnp(2,inode)*shp.velz.dz.s(inode);

          const int gid = nodeids[inode];
          const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);

          for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
              enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
          {
            if (enrfield->getField() == XFEM::PHYSICS::Veln)
            {
              vderxy(0,0) += evelnp(3,velncounter)*shp.velx.dx.n(velncounter);
              vderxy(0,1) += evelnp(3,velncounter)*shp.velx.dy.n(velncounter);
              vderxy(0,2) += evelnp(3,velncounter)*shp.velx.dz.n(velncounter);

              vderxy(1,0) += evelnp(3,velncounter)*shp.vely.dx.n(velncounter);
              vderxy(1,1) += evelnp(3,velncounter)*shp.vely.dy.n(velncounter);
              vderxy(1,2) += evelnp(3,velncounter)*shp.vely.dz.n(velncounter);

              vderxy(2,0) += evelnp(3,velncounter)*shp.velz.dx.n(velncounter);
              vderxy(2,1) += evelnp(3,velncounter)*shp.velz.dy.n(velncounter);
              vderxy(2,2) += evelnp(3,velncounter)*shp.velz.dz.n(velncounter);

              velncounter += 1;
            }
          }
        }
        if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
      }
//cout << "vderxy " << vderxy << endl;

      //-----------------------------------------------------
      // get second velocity derivatives at integration point
      //-----------------------------------------------------
      static LINALG::Matrix<3,6> vderxy2;
      vderxy2.Clear();
      //vderxy2 = evelnp(i,k)*enr_derxy2(j,k);
      if (higher_order_ele)
      {
        const int* nodeids = ele->NodeIds();
        std::size_t velncounter = 0;

        // vderxy = enr_derxy(j,k)*evelnp(i,k);
        for (size_t inode = 0; inode < numparamvelx; ++inode)
        {
          // standard shape functions are identical for all vector components
          // e.g. shp.velx.dxdx.s == shp.vely.dxdx.s == shp.velz.dxdx.s
          vderxy2(0,0) += evelnp(0,inode)*shp.velx.dxdx.s(inode);
          vderxy2(0,1) += evelnp(0,inode)*shp.velx.dydy.s(inode);
          vderxy2(0,2) += evelnp(0,inode)*shp.velx.dzdz.s(inode);
          vderxy2(0,3) += evelnp(0,inode)*shp.velx.dxdy.s(inode);
          vderxy2(0,4) += evelnp(0,inode)*shp.velx.dxdz.s(inode);
          vderxy2(0,5) += evelnp(0,inode)*shp.velx.dydz.s(inode);

          vderxy2(1,0) += evelnp(1,inode)*shp.vely.dxdx.s(inode);
          vderxy2(1,1) += evelnp(1,inode)*shp.vely.dydy.s(inode);
          vderxy2(1,2) += evelnp(1,inode)*shp.vely.dzdz.s(inode);
          vderxy2(1,3) += evelnp(1,inode)*shp.vely.dxdy.s(inode);
          vderxy2(1,4) += evelnp(1,inode)*shp.vely.dxdz.s(inode);
          vderxy2(1,5) += evelnp(1,inode)*shp.vely.dydz.s(inode);

          vderxy2(2,0) += evelnp(2,inode)*shp.velz.dxdx.s(inode);
          vderxy2(2,1) += evelnp(2,inode)*shp.velz.dydy.s(inode);
          vderxy2(2,2) += evelnp(2,inode)*shp.velz.dzdz.s(inode);
          vderxy2(2,3) += evelnp(2,inode)*shp.velz.dxdy.s(inode);
          vderxy2(2,4) += evelnp(2,inode)*shp.velz.dxdz.s(inode);
          vderxy2(2,5) += evelnp(2,inode)*shp.velz.dydz.s(inode);

          const int gid = nodeids[inode];
          const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);

          for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
              enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
          {
            if (enrfield->getField() == XFEM::PHYSICS::Veln)
            {
              vderxy2(0,0) += evelnp(3,velncounter)*shp.velx.dxdx.n(velncounter);
              vderxy2(0,1) += evelnp(3,velncounter)*shp.velx.dydy.n(velncounter);
              vderxy2(0,2) += evelnp(3,velncounter)*shp.velx.dzdz.n(velncounter);
              vderxy2(0,3) += evelnp(3,velncounter)*shp.velx.dxdy.n(velncounter);
              vderxy2(0,4) += evelnp(3,velncounter)*shp.velx.dxdz.n(velncounter);
              vderxy2(0,5) += evelnp(3,velncounter)*shp.velx.dydz.n(velncounter);

              vderxy2(1,0) += evelnp(3,velncounter)*shp.vely.dxdx.n(velncounter);
              vderxy2(1,1) += evelnp(3,velncounter)*shp.vely.dydy.n(velncounter);
              vderxy2(1,2) += evelnp(3,velncounter)*shp.vely.dzdz.n(velncounter);
              vderxy2(1,3) += evelnp(3,velncounter)*shp.vely.dxdy.n(velncounter);
              vderxy2(1,4) += evelnp(3,velncounter)*shp.vely.dxdz.n(velncounter);
              vderxy2(1,5) += evelnp(3,velncounter)*shp.vely.dydz.n(velncounter);

              vderxy2(2,0) += evelnp(3,velncounter)*shp.velz.dxdx.n(velncounter);
              vderxy2(2,1) += evelnp(3,velncounter)*shp.velz.dydy.n(velncounter);
              vderxy2(2,2) += evelnp(3,velncounter)*shp.velz.dzdz.n(velncounter);
              vderxy2(2,3) += evelnp(3,velncounter)*shp.velz.dxdy.n(velncounter);
              vderxy2(2,4) += evelnp(3,velncounter)*shp.velz.dxdz.n(velncounter);
              vderxy2(2,5) += evelnp(3,velncounter)*shp.velz.dydz.n(velncounter);

              velncounter += 1;
            }
          }
        }
        if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
      }

      //----------------------------------
      // get pressure at integration point
      //----------------------------------
      double pres = 0.0;
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
        pres += eprenp(iparam) * shp_pres.d0(iparam);

      //--------------------------------------------
      // get pressure gradients at integration point
      //--------------------------------------------
      static LINALG::Matrix<nsd,1> gradp;
      gradp.Clear();
      // gradp = enr_derxy(i,j)*eprenp(j);
      for (size_t iparam = 0; iparam != numparampres; ++iparam)
      {
        gradp(0) += eprenp(iparam) * shp_pres.dx(iparam);
        gradp(1) += eprenp(iparam) * shp_pres.dy(iparam);
        gradp(2) += eprenp(iparam) * shp_pres.dz(iparam);
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
      // compute streamlength based on element geometry (standard shape functions)
      //const double strle = FLD::UTILS::Streamlength(shp.velx.dx.s,shp.vely.dy.s,shp.velz.dz.s, gpvelnp, vel_norm, numparamvelx);
      // compute streamlength based on cell geometry (enriched shape functions)
      //const double strle = FLD::UTILS::Streamlength(shp, gpvelnp, vel_norm, numparamvelx);
      double tau_stab_Mu = 0.0;
      double tau_stab_Mp = 0.0;
      double tau_stab_C  = 0.0;
      // remark: stabilization parameters 'tau_stab_Mu' and 'tau_stab_Mp' are both computes based on
      //         'hk' instead of using the streamlength for 'tau_stab_Mu'
      FLD::UTILS::computeStabilizationParams(gpvelnp, xji,
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
      static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> enr_conv_c(false);
      enr_conv_c.Clear();

      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_conv_c.velx.d0.s(iparam) += shp.velx.dx.s(iparam)*gpvelnp(0);
        enr_conv_c.velx.d0.s(iparam) += shp.velx.dy.s(iparam)*gpvelnp(1);
        enr_conv_c.velx.d0.s(iparam) += shp.velx.dz.s(iparam)*gpvelnp(2);
      }
      enr_conv_c.vely.d0.s = enr_conv_c.velx.d0.s;
      enr_conv_c.velz.d0.s = enr_conv_c.velx.d0.s;
      for (size_t iparam = 0; iparam != numparamveln; ++iparam)
      {
        enr_conv_c.velx.d0.n(iparam) += shp.velx.dx.n(iparam)*gpvelnp(0);
        enr_conv_c.velx.d0.n(iparam) += shp.velx.dy.n(iparam)*gpvelnp(1);
        enr_conv_c.velx.d0.n(iparam) += shp.velx.dz.n(iparam)*gpvelnp(2);

        enr_conv_c.vely.d0.n(iparam) += shp.vely.dx.n(iparam)*gpvelnp(0);
        enr_conv_c.vely.d0.n(iparam) += shp.vely.dy.n(iparam)*gpvelnp(1);
        enr_conv_c.vely.d0.n(iparam) += shp.vely.dz.n(iparam)*gpvelnp(2);

        enr_conv_c.velz.d0.n(iparam) += shp.velz.dx.n(iparam)*gpvelnp(0);
        enr_conv_c.velz.d0.n(iparam) += shp.velz.dy.n(iparam)*gpvelnp(1);
        enr_conv_c.velz.d0.n(iparam) += shp.velz.dz.n(iparam)*gpvelnp(2);
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
      static XFEM::EnrViscs2N<shpVecSizeVel> enr_viscs2;
      //            enr_viscs2.Clear();
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enr_viscs2.xx.s(iparam) = 0.5 * (2.0 * shp.velx.dxdx.s(iparam) + shp.velx.dydy.s(iparam) + shp.velx.dzdz.s(iparam));
        enr_viscs2.xy.s(iparam) = 0.5 *  shp.vely.dxdy.s(iparam);
        enr_viscs2.xz.s(iparam) = 0.5 *  shp.velz.dxdz.s(iparam);
        enr_viscs2.yx.s(iparam) = 0.5 *  shp.velx.dydx.s(iparam);
        enr_viscs2.yy.s(iparam) = 0.5 * (shp.vely.dxdx.s(iparam) + 2.0 * shp.vely.dydy.s(iparam) + shp.vely.dzdz.s(iparam));
        enr_viscs2.yz.s(iparam) = 0.5 *  shp.vely.dydz.s(iparam);
        enr_viscs2.zx.s(iparam) = 0.5 *  shp.velx.dzdx.s(iparam);
        enr_viscs2.zy.s(iparam) = 0.5 *  shp.vely.dzdy.s(iparam);
        enr_viscs2.zz.s(iparam) = 0.5 * (shp.velz.dxdx.s(iparam) + shp.velz.dydy.s(iparam) + 2.0 * shp.velz.dzdz.s(iparam));
      }
      for (size_t iparam = 0; iparam != numparamveln; ++iparam)
      {
        enr_viscs2.xx.n(iparam) = 0.5 * (2.0 * shp.velx.dxdx.n(iparam) + shp.velx.dydy.n(iparam) + shp.velx.dzdz.n(iparam));
        enr_viscs2.xy.n(iparam) = 0.5 *  shp.vely.dxdy.n(iparam);
        enr_viscs2.xz.n(iparam) = 0.5 *  shp.velz.dxdz.n(iparam);
        enr_viscs2.yx.n(iparam) = 0.5 *  shp.velx.dydx.n(iparam);
        enr_viscs2.yy.n(iparam) = 0.5 * (shp.vely.dxdx.n(iparam) + 2.0 * shp.vely.dydy.n(iparam) + shp.vely.dzdz.n(iparam));
        enr_viscs2.yz.n(iparam) = 0.5 *  shp.vely.dydz.n(iparam);
        enr_viscs2.zx.n(iparam) = 0.5 *  shp.velx.dzdx.n(iparam);
        enr_viscs2.zy.n(iparam) = 0.5 *  shp.vely.dzdy.n(iparam);
        enr_viscs2.zz.n(iparam) = 0.5 * (shp.velz.dxdx.n(iparam) + shp.velz.dydy.n(iparam) + 2.0 * shp.velz.dzdz.n(iparam));
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

      //----------------------------------------
      // build stiffness matrix standard entries
      //----------------------------------------
      COMBUST::BuildStandardDomainIntegralsNormal<DISTYPE,ASSTYPE,NUMDOF,shpVecSizeVel,shpVecSizePres>(
          assembler, shp, shp_pres, fac, timefac, timefacfac, densaf, dynvisc,
          gpvelnp, pres, accam, gradp, vderxy, rhsint, res_old, visc_old, enr_conv_c, enr_viscs2,
          instationary, genalpha, newton, pstab, supg, graddiv,
          tau_stab_Mu, tau_stab_Mp, tau_stab_C);
#ifndef COMBUST_DECOUPLEDXFEM
      if (epsilonele_unknowns_present)
      {
      COMBUST::BuildDomainIntegralsStressNormal<DISTYPE,ASSTYPE,NUMDOF,shpVecSizeVel,shpVecSizePres,shpVecSizeStress,shpVecSizeDiscPres>(
          assembler, shp, shp_pres, shp_epsilon, shp_discpres, fac, timefac, timefacfac,
          dynvisc, pres, vderxy, epsilon, discpres);
      }
      if (sigmaele_unknowns_present)
      {
        dserror("uncomment this!");
//        COMBUST::BuildDomainIntegralsSigma<DISTYPE,ASSTYPE,NUMDOF,shpVecSize,shpVecSizeStress,shpVecSizeDiscPres>(
//            assembler, shp, shp_sigma, fac, timefac, timefacfac,
//            dynvisc, pres, vderxy, epsilon);
      }
#endif

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
void SysmatBoundaryStressNormal(
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

#ifdef COMBUST_EPSPRES_BASED_CONSISTPRESJUMP
  // pressure jump value
  const double jp = 0.0; // (flamespeed*flamespeed*dens_minus*dens_minus)*(1/dens_plus - 1/dens_minus);
#endif

  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;

  // get number of parameters (dofs) for each field
  // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparamveln = XFEM::NumParam<0,ASSTYPE>::get(dofman, XFEM::PHYSICS::Veln);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // put one here to create arrays of size 1, since they are not needed anyway
  // in the xfem assembly, numparam is determined by the dofmanager
  const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Epsilonxx);
  const size_t numparamdiscpres = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 1);

  // get domain integration cells for this element
  const GEO::BoundaryIntCells& boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

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
#ifndef COMBUST_NORMAL_ENRICHMENT
      const std::size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
#else
      const size_t shpVecSizeVel    = COMBUST::SizeFacVel<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
      const size_t shpVecSizePres   = COMBUST::SizeFacPres<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
#endif
      const std::size_t shpVecSizeStress   = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
      const std::size_t shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

      // temporary arrays holding enriched shape functions (N * \Psi) on either side of the interface
      static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> enrfunct_vel_plus(false);
      enrfunct_vel_plus.Clear();
      static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> enrfunct_vel_minus(false);
      enrfunct_vel_minus.Clear();

      static LINALG::Matrix<shpVecSizePres,1> enrfunct_pres_plus;
      enrfunct_pres_plus.Clear();
      static LINALG::Matrix<shpVecSizePres,1> enrfunct_pres_minus;
      enrfunct_pres_minus.Clear();

      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_plus;
      enrfunct_stress_plus.Clear();
      static LINALG::Matrix<shpVecSizeStress,1>   enrfunct_stress_minus;
      enrfunct_stress_minus.Clear();
      static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres_plus;
      enrfunct_discpres_plus.Clear();
      static LINALG::Matrix<shpVecSizeDiscPres,1> enrfunct_discpres_minus;
      enrfunct_discpres_minus.Clear();

      // shape functions for nodal parameters (dofs) on plus and minus side
      enrvals_plus.ComputeModifiedEnrichedNodalShapefunction(Pres, funct, enrfunct_pres_plus);
      enrvals_minus.ComputeModifiedEnrichedNodalShapefunction(Pres, funct, enrfunct_pres_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(Epsilonxx, funct_stress, enrfunct_stress_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(Epsilonxx, funct_stress, enrfunct_stress_minus);

      enrvals_plus.ComputeShiftedEnrichedElementShapefunction(DiscPres, funct_discpres, enrfunct_discpres_plus);
      enrvals_minus.ComputeShiftedEnrichedElementShapefunction(DiscPres, funct_discpres, enrfunct_discpres_minus);

      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enrfunct_vel_plus.velx.d0.s(iparam) = funct(iparam);
        enrfunct_vel_plus.vely.d0.s(iparam) = funct(iparam);
        enrfunct_vel_plus.velz.d0.s(iparam) = funct(iparam);
      }
      // fill approximation functions for XFEM
      for (size_t iparam = 0; iparam != numparamvelx; ++iparam)
      {
        enrfunct_vel_minus.velx.d0.s(iparam) = funct(iparam);
        enrfunct_vel_minus.vely.d0.s(iparam) = funct(iparam);
        enrfunct_vel_minus.velz.d0.s(iparam) = funct(iparam);
      }
      // shape functions and derivatives for nodal parameters (dofs)
      enrvals_minus.ComputeNormalShapeFunction(funct, egradphi,
#ifdef COLLAPSE_FLAME
          normal,
#endif
          enrfunct_vel_minus);
      enrvals_plus.ComputeNormalShapeFunction(funct, egradphi,
#ifdef COLLAPSE_FLAME
          normal,
#endif
          enrfunct_vel_plus);

      // (enriched) shape functions = approximation functions (P = N * \Psi)
#ifdef COMBUST_NORMAL_ENRICHMENT
      // second derivatives needed for velocity
      static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> shp_vel_jump(false);
      shp_vel_jump.Clear();
      static XFEM::ApproxFuncNormalVector<0,shpVecSizeVel> shp_vel_mean(false);
      shp_vel_mean.Clear();
      static XFEM::ApproxFunc<0,shpVecSizePres> shp_pres_jump;
      static XFEM::ApproxFunc<0,shpVecSizePres> shp_pres_mean;
#else
      static XFEM::ApproxFunc<0,shpVecSize> shp_jump;
      static XFEM::ApproxFunc<0,shpVecSize> shp_mean;
#endif

      for (std::size_t iparam = 0; iparam < numparamvelx; ++iparam)
      {
        // standard jump shape functions are zero by definition
        //shp_vel_jump.velx.d0.s(iparam) = 0.0;
        //shp_vel_jump.vely.d0.s(iparam) = 0.0;
        //shp_vel_jump.velz.d0.s(iparam) = 0.0;

        shp_vel_mean.velx.d0.s(iparam) = funct(iparam);
        shp_vel_mean.vely.d0.s(iparam) = funct(iparam);
        shp_vel_mean.velz.d0.s(iparam) = funct(iparam);
      }
      for (std::size_t iparam = 0; iparam < numparamveln; ++iparam)
      {
        shp_vel_jump.velx.d0.n(iparam) = enrfunct_vel_minus.velx.d0.n(iparam) - enrfunct_vel_plus.velx.d0.n(iparam);
        shp_vel_jump.vely.d0.n(iparam) = enrfunct_vel_minus.vely.d0.n(iparam) - enrfunct_vel_plus.vely.d0.n(iparam);
        shp_vel_jump.velz.d0.n(iparam) = enrfunct_vel_minus.velz.d0.n(iparam) - enrfunct_vel_plus.velz.d0.n(iparam);

        shp_vel_mean.velx.d0.n(iparam) = 0.5*(enrfunct_vel_minus.velx.d0.n(iparam) + enrfunct_vel_plus.velx.d0.n(iparam));
        shp_vel_mean.vely.d0.n(iparam) = 0.5*(enrfunct_vel_minus.vely.d0.n(iparam) + enrfunct_vel_plus.vely.d0.n(iparam));
        shp_vel_mean.velz.d0.n(iparam) = 0.5*(enrfunct_vel_minus.velz.d0.n(iparam) + enrfunct_vel_plus.velz.d0.n(iparam));
      }

      // fill approximation functions
      for (std::size_t iparam = 0; iparam < numparampres; ++iparam)
      {
        shp_pres_jump.d0(iparam) = enrfunct_pres_minus(iparam) - enrfunct_pres_plus(iparam);
        shp_pres_mean.d0(iparam) = 0.5*(enrfunct_pres_plus(iparam) + enrfunct_pres_minus(iparam));
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
      veljump = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, shp_vel_jump, numparamvelx, ele, dofman);

      // velocity minus
      static LINALG::Matrix<nsd,1> velminus(true);
      velminus.Clear();
      velminus = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, enrfunct_vel_minus, numparamvelx, ele, dofman);

      // velocity plus
      static LINALG::Matrix<nsd,1> velplus(true);
      velplus.Clear();
      velplus = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, enrfunct_vel_plus, numparamvelx, ele, dofman);

      // mean velocity
      static LINALG::Matrix<nsd,1> velmean(true);
      velmean.Clear();
      velmean = XFEM::interpolateVectorFieldToIntPointNormal(evelnp, shp_vel_mean, numparamvelx, ele, dofman);

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

      assembler.template Matrix<Velx,Epsilonxx>(shp_mean.d0, timefacfac*normal(0), shp_epsilon_jump.d0);
      assembler.template Matrix<Velx,Epsilonxy>(shp_mean.d0, timefacfac*normal(1), shp_epsilon_jump.d0);
      assembler.template Matrix<Velx,Epsilonxz>(shp_mean.d0, timefacfac*normal(2), shp_epsilon_jump.d0);

      assembler.template Matrix<Vely,Epsilonyx>(shp_mean.d0, timefacfac*normal(0), shp_epsilon_jump.d0);
      assembler.template Matrix<Vely,Epsilonyy>(shp_mean.d0, timefacfac*normal(1), shp_epsilon_jump.d0);
      assembler.template Matrix<Vely,Epsilonyz>(shp_mean.d0, timefacfac*normal(2), shp_epsilon_jump.d0);

      assembler.template Matrix<Velz,Epsilonzx>(shp_mean.d0, timefacfac*normal(0), shp_epsilon_jump.d0);
      assembler.template Matrix<Velz,Epsilonzy>(shp_mean.d0, timefacfac*normal(1), shp_epsilon_jump.d0);
      assembler.template Matrix<Velz,Epsilonzz>(shp_mean.d0, timefacfac*normal(2), shp_epsilon_jump.d0);

      //   |                             |
      // - |  < v > , || epsilon_i || n  |
      //   |                             |

      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac*normal(0)*epsilonjump(0,0));
      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac*normal(1)*epsilonjump(0,1));
      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac*normal(2)*epsilonjump(0,2));

      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac*normal(0)*epsilonjump(1,0));
      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac*normal(1)*epsilonjump(1,1));
      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac*normal(2)*epsilonjump(1,2));

      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac*normal(0)*epsilonjump(2,0));
      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac*normal(1)*epsilonjump(2,1));
      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac*normal(2)*epsilonjump(2,2));
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTPRESJUMP
      //--------------------------    |                       |
      // pressure consistency term  - |  < v >, || Dp^e || n  |
      //--------------------------    |                       |

      assembler.template Matrix<Velx,DiscPres>(shp_mean.d0, -timefacfac*normal(0), shp_discpres_jump.d0);
      assembler.template Matrix<Vely,DiscPres>(shp_mean.d0, -timefacfac*normal(1), shp_discpres_jump.d0);
      assembler.template Matrix<Velz,DiscPres>(shp_mean.d0, -timefacfac*normal(2), shp_discpres_jump.d0);

      //   |                        |
      // + |  < v >, || p^e_i || n  |
      //   |                        |

      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(0)*pdiscjump);
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(1)*pdiscjump);
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(2)*pdiscjump);

      //---------  ---------------------------------    |                |
      // pressure (flux) jump (consistency) term RHS  + |  < v >, j_p*n  |
      //-----------  -------------------------------    |                |
      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(0)*jp);
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(1)*jp);
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(2)*jp);
#endif

#ifdef COMBUST_EPSPRES_BASED_FLUXJUMP
      //---------------------------------------------------    |                 |
      // flux (traction/stress) jump (consistency) term RHS  - |  < v >, j_flux  |
      //---------------------------------------------------    |                 |
      assembler.template VectorN<Velx>(shp_vel_mean.velx.d0, -timefacfac*jflux(0));
      assembler.template VectorN<Vely>(shp_vel_mean.vely.d0, -timefacfac*jflux(1));
      assembler.template VectorN<Velz>(shp_vel_mean.velz.d0, -timefacfac*jflux(2));

#if 0
      // j_sigma is a 3x3 Matrix -> j_sigma*normal
      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(0)*jflux(0,0));
      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(1)*jflux(0,1));
      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*normal(2)*jflux(0,2));

      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(0)*jflux(1,0));
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(1)*jflux(1,1));
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*normal(2)*jflux(1,2));

      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(0)*jflux(2,0));
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(1)*jflux(2,1));
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*normal(2)*jflux(2,2));
#endif
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTVISCMEAN
      //-------------------------    |                            |
      // viscous consistency term  + |  || v || , < Depsilon > n  |
      //-------------------------    |                            |
      assembler.template MatrixNS<Velx,Epsilonxx>(shp_vel_jump.velx.d0, timefacfac*normal(0), shp_eps_visc_mean.d0);
      assembler.template MatrixNS<Velx,Epsilonxy>(shp_vel_jump.velx.d0, timefacfac*normal(1), shp_eps_visc_mean.d0);
      assembler.template MatrixNS<Velx,Epsilonxz>(shp_vel_jump.velx.d0, timefacfac*normal(2), shp_eps_visc_mean.d0);

      assembler.template MatrixNS<Vely,Epsilonyx>(shp_vel_jump.vely.d0, timefacfac*normal(0), shp_eps_visc_mean.d0);
      assembler.template MatrixNS<Vely,Epsilonyy>(shp_vel_jump.vely.d0, timefacfac*normal(1), shp_eps_visc_mean.d0);
      assembler.template MatrixNS<Vely,Epsilonyz>(shp_vel_jump.vely.d0, timefacfac*normal(2), shp_eps_visc_mean.d0);

      assembler.template MatrixNS<Velz,Epsilonzx>(shp_vel_jump.velz.d0, timefacfac*normal(0), shp_eps_visc_mean.d0);
      assembler.template MatrixNS<Velz,Epsilonzy>(shp_vel_jump.velz.d0, timefacfac*normal(1), shp_eps_visc_mean.d0);
      assembler.template MatrixNS<Velz,Epsilonzz>(shp_vel_jump.velz.d0, timefacfac*normal(2), shp_eps_visc_mean.d0);

      //   |                             |
      // - |  || v || , < epsilon_i > n  |
      //   |                             |
      assembler.template VectorN<Velx>(shp_vel_jump.velx.d0, -timefacfac*normal(0)*epsilonmean(0,0));
      assembler.template VectorN<Velx>(shp_vel_jump.velx.d0, -timefacfac*normal(1)*epsilonmean(0,1));
      assembler.template VectorN<Velx>(shp_vel_jump.velx.d0, -timefacfac*normal(2)*epsilonmean(0,2));

      assembler.template VectorN<Vely>(shp_vel_jump.vely.d0, -timefacfac*normal(0)*epsilonmean(1,0));
      assembler.template VectorN<Vely>(shp_vel_jump.vely.d0, -timefacfac*normal(1)*epsilonmean(1,1));
      assembler.template VectorN<Vely>(shp_vel_jump.vely.d0, -timefacfac*normal(2)*epsilonmean(1,2));

      assembler.template VectorN<Velz>(shp_vel_jump.velz.d0, -timefacfac*normal(0)*epsilonmean(2,0));
      assembler.template VectorN<Velz>(shp_vel_jump.velz.d0, -timefacfac*normal(1)*epsilonmean(2,1));
      assembler.template VectorN<Velz>(shp_vel_jump.velz.d0, -timefacfac*normal(2)*epsilonmean(2,2));
#endif

#ifdef COMBUST_EPSPRES_BASED_CONSISTPRESMEAN
      //--------------------------    |                       |
      // pressure consistency term  - |  || v ||, < Dp^e > n  |
      //--------------------------    |                       |

      assembler.template MatrixNS<Velx,DiscPres>(shp_vel_jump.velx.d0, -timefacfac*normal(0), shp_discpres_mean.d0);
      assembler.template MatrixNS<Vely,DiscPres>(shp_vel_jump.vely.d0, -timefacfac*normal(1), shp_discpres_mean.d0);
      assembler.template MatrixNS<Velz,DiscPres>(shp_vel_jump.velz.d0, -timefacfac*normal(2), shp_discpres_mean.d0);

      //   |                        |
      // + |  || v ||, < p^e_i > n  |
      //   |                        |

      assembler.template VectorN<Velx>(shp_vel_jump.velx.d0, timefacfac*normal(0)*pdiscmean);
      assembler.template VectorN<Vely>(shp_vel_jump.vely.d0, timefacfac*normal(1)*pdiscmean);
      assembler.template VectorN<Velz>(shp_vel_jump.velz.d0, timefacfac*normal(2)*pdiscmean);
#endif

#ifdef COMBUST_EPSPRES_BASED_ADJCONSISTVISCJUMP
      //---------------------------------    |                             |
      // viscous adjoint consistency term  + |  < epsilon^e > n, || Du ||  |
      //---------------------------------    |                             |
      assembler.template MatrixSN<Epsilonxx,Velx>(shp_eps_visc_mean.d0, timefacfac*normal(0), shp_vel_jump.velx.d0);
      assembler.template MatrixSN<Epsilonxy,Velx>(shp_eps_visc_mean.d0, timefacfac*normal(1), shp_vel_jump.velx.d0);
      assembler.template MatrixSN<Epsilonxz,Velx>(shp_eps_visc_mean.d0, timefacfac*normal(2), shp_vel_jump.velx.d0);

      assembler.template MatrixSN<Epsilonyx,Vely>(shp_eps_visc_mean.d0, timefacfac*normal(0), shp_vel_jump.vely.d0);
      assembler.template MatrixSN<Epsilonyy,Vely>(shp_eps_visc_mean.d0, timefacfac*normal(1), shp_vel_jump.vely.d0);
      assembler.template MatrixSN<Epsilonyz,Vely>(shp_eps_visc_mean.d0, timefacfac*normal(2), shp_vel_jump.vely.d0);

      assembler.template MatrixSN<Epsilonzx,Velz>(shp_eps_visc_mean.d0, timefacfac*normal(0), shp_vel_jump.velz.d0);
      assembler.template MatrixSN<Epsilonzy,Velz>(shp_eps_visc_mean.d0, timefacfac*normal(1), shp_vel_jump.velz.d0);
      assembler.template MatrixSN<Epsilonzz,Velz>(shp_eps_visc_mean.d0, timefacfac*normal(2), shp_vel_jump.velz.d0);

      //   |                              |
      // - |  < epsilon^e > n, || u_i ||  |
      //   |                              |

      assembler.template Vector<Epsilonxx>(shp_eps_visc_mean.d0, -timefacfac*normal(0)*veljump(0));
      assembler.template Vector<Epsilonxy>(shp_eps_visc_mean.d0, -timefacfac*normal(1)*veljump(0));
      assembler.template Vector<Epsilonxz>(shp_eps_visc_mean.d0, -timefacfac*normal(2)*veljump(0));

      assembler.template Vector<Epsilonyx>(shp_eps_visc_mean.d0, -timefacfac*normal(0)*veljump(1));
      assembler.template Vector<Epsilonyy>(shp_eps_visc_mean.d0, -timefacfac*normal(1)*veljump(1));
      assembler.template Vector<Epsilonyz>(shp_eps_visc_mean.d0, -timefacfac*normal(2)*veljump(1));

      assembler.template Vector<Epsilonzx>(shp_eps_visc_mean.d0, -timefacfac*normal(0)*veljump(2));
      assembler.template Vector<Epsilonzy>(shp_eps_visc_mean.d0, -timefacfac*normal(1)*veljump(2));
      assembler.template Vector<Epsilonzz>(shp_eps_visc_mean.d0, -timefacfac*normal(2)*veljump(2));

      //-------------------------------------    |                        |
      // viscous adjoint consistency term RHS  + |  < epsilon^e > n, j_u  |
      //-------------------------------------    |                        |

#if 1
      assembler.template Vector<Epsilonxx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*ju(0));
      assembler.template Vector<Epsilonxy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*ju(0));
      assembler.template Vector<Epsilonxz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*ju(0));
#else
      assembler.template Vector<Epsilonxx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*normal(0)*ju_n);
      assembler.template Vector<Epsilonxy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*normal(0)*ju_n);
      assembler.template Vector<Epsilonxz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*normal(0)*ju_n);
#endif

#if 1
      assembler.template Vector<Epsilonyx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*ju(1));
      assembler.template Vector<Epsilonyy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*ju(1));
      assembler.template Vector<Epsilonyz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*ju(1));
#else
      assembler.template Vector<Epsilonyx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*normal(1)*ju_n);
      assembler.template Vector<Epsilonyy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*normal(1)*ju_n);
      assembler.template Vector<Epsilonyz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*normal(1)*ju_n);
#endif

#if 1
      assembler.template Vector<Epsilonzx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*ju(2));
      assembler.template Vector<Epsilonzy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*ju(2));
      assembler.template Vector<Epsilonzz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*ju(2));
#else
      assembler.template Vector<Epsilonzx>(shp_eps_visc_mean.d0, timefacfac*normal(0)*normal(2)*ju_n);
      assembler.template Vector<Epsilonzy>(shp_eps_visc_mean.d0, timefacfac*normal(1)*normal(2)*ju_n);
      assembler.template Vector<Epsilonzz>(shp_eps_visc_mean.d0, timefacfac*normal(2)*normal(2)*ju_n);
#endif
#endif

#ifdef COMBUST_EPSPRES_BASED_ADJCONSISTPRESJUMP
      //----------------------------------      |                       |
      // pressure adjoint consistency term  +/- |  < q^e >, || Du || n  |
      //----------------------------------      |                       |
      // remark: the sign of this term is not clear, yet
      assembler.template MatrixSN<DiscPres,Velx>(shp_discpres_mean.d0, timefacfac*normal(0), shp_vel_jump.velx.d0);
      assembler.template MatrixSN<DiscPres,Vely>(shp_discpres_mean.d0, timefacfac*normal(1), shp_vel_jump.vely.d0);
      assembler.template MatrixSN<DiscPres,Velz>(shp_discpres_mean.d0, timefacfac*normal(2), shp_vel_jump.velz.d0);

      //     |                        |
      // -/+ |  < q^e >, || u_i || n  |
      //     |                        |
      // remark: the sign of this term is not clear, yet
      assembler.template Vector<DiscPres>(shp_discpres_mean.d0, -timefacfac*(veljump(0)*normal(0)
                                                                            +veljump(1)*normal(1)
                                                                            +veljump(2)*normal(2)));

      //--------------------------------------      |                  |
      // pressure adjoint consistency term RHS  +/- |  < q^e >, j_u n  |
      //--------------------------------------      |                  |
      // remark: the sign of this term is not clear, yet
//       assembler.template Vector<DiscPres>(shp_discpres_mean.d0, -timefacfac*ju);
       assembler.template Vector<DiscPres>(shp_discpres_mean.d0, timefacfac*(ju(0)*normal(0)
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
      assembler.template Matrix<Velx,Velx>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);
      assembler.template Matrix<Vely,Vely>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);
      assembler.template Matrix<Velz,Velz>(shp_inflow.d0, -timefacfac*inflowvel, shp_inflow.d0);

      //   |                    |
      // + |  v, u_n*|| u_i ||  |
      //   |                    |
      assembler.template Vector<Velx>(shp_inflow.d0, timefacfac*inflowvel*velinflow(0));
      assembler.template Vector<Vely>(shp_inflow.d0, timefacfac*inflowvel*velinflow(1));
      assembler.template Vector<Velz>(shp_inflow.d0, timefacfac*inflowvel*velinflow(2));
#endif
      //-----------------    |                       |
      // inflow flux term  - |  < v >, u_n*|| Du ||  |
      //-----------------    |                       |
      const double inflowvel = -1.0*(velmean(0)*normal(0) + velmean(1)*normal(1) + velmean(2)*normal(2));
      //const double inflowvel = 1000*(veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2));
//      cout << "inflowvel" << inflowvel << endl;
//      assembler.template Matrix<Velx,Velx>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<Vely,Vely>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<Velz,Velz>(shp_mean.d0, -timefacfac*inflowvel, shp_jump.d0);

      assembler.template Matrix<Velx,Velx>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_jump.d0, timefacfac*inflowvel, shp_jump.d0);

//      assembler.template Matrix<Velx,Velx>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<Vely,Vely>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);
//      assembler.template Matrix<Velz,Velz>(enrfunct_minus, -timefacfac*inflowvel, shp_jump.d0);

      //   |                        |
      // + |  < v >, u_n*|| u_i ||  |
      //   |                        |
//      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*inflowvel*veljump(0));
//      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*inflowvel*veljump(1));
//      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*inflowvel*veljump(2));

      assembler.template Vector<Velx>(shp_jump.d0, -timefacfac*inflowvel*veljump(0));
      assembler.template Vector<Vely>(shp_jump.d0, -timefacfac*inflowvel*veljump(1));
      assembler.template Vector<Velz>(shp_jump.d0, -timefacfac*inflowvel*veljump(2));

//      assembler.template Vector<Velx>(enrfunct_minus, timefacfac*inflowvel*veljump(0));
//      assembler.template Vector<Vely>(enrfunct_minus, timefacfac*inflowvel*veljump(1));
//      assembler.template Vector<Velz>(enrfunct_minus, timefacfac*inflowvel*veljump(2));
      //   |                 |
      // + |  < v >, u_n*ju  |
      //   |                 |
//      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac*inflowvel*ju(0));
//      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac*inflowvel*ju(1));
//      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac*inflowvel*ju(2));

      assembler.template Vector<Velx>(shp_jump.d0, timefacfac*inflowvel*ju(0));
      assembler.template Vector<Vely>(shp_jump.d0, timefacfac*inflowvel*ju(1));
      assembler.template Vector<Velz>(shp_jump.d0, timefacfac*inflowvel*ju(2));

//      assembler.template Vector<Velx>(enrfunct_minus, -timefacfac*inflowvel*ju(0));
//      assembler.template Vector<Vely>(enrfunct_minus, -timefacfac*inflowvel*ju(1));
//      assembler.template Vector<Velz>(enrfunct_minus, -timefacfac*inflowvel*ju(2));

#endif

#ifdef COMBUST_EPSPRES_BASED_ADVECTIVEFLUX
      //-----------------    |                   |
      // inflow flux term  + |  < v >, || Du ||  |
      //-----------------    |                   |
      const double inflowvel = veljump(0)*normal(0) + veljump(1)*normal(1) + veljump(2)*normal(2);
      assembler.template Matrix<Velx,Velx>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<Vely,Vely>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);
      assembler.template Matrix<Velz,Velz>(shp_mean.d0, timefacfac*inflowvel, shp_jump.d0);

      //   |                    |
      // - |  < v >, || u_i ||  |
      //   |                    |
      assembler.template Vector<Velx>(shp_mean.d0, -timefacfac*inflowvel*veljump(0));
      assembler.template Vector<Vely>(shp_mean.d0, -timefacfac*inflowvel*veljump(1));
      assembler.template Vector<Velz>(shp_mean.d0, -timefacfac*inflowvel*veljump(2));

      //-----------------       |                |
      // inflow flux term RHS + |  < v >, j_u n  |
      //-----------------       |                |
      assembler.template Vector<Velx>(shp_mean.d0, timefacfac*inflowvel*ju*normal(0));
      assembler.template Vector<Vely>(shp_mean.d0, timefacfac*inflowvel*ju*normal(1));
      assembler.template Vector<Velz>(shp_mean.d0, timefacfac*inflowvel*ju*normal(2));

      //-----------------    |                     |
      // inflow flux term  - |  < q >, || Du || n  |
      //-----------------    |                     |
      assembler.template Matrix<Pres,Velx>(shp_mean.d0, timefacfac*inflowvel*normal(0), shp_jump.d0);
      assembler.template Matrix<Pres,Vely>(shp_mean.d0, timefacfac*inflowvel*normal(1), shp_jump.d0);
      assembler.template Matrix<Pres,Velz>(shp_mean.d0, timefacfac*inflowvel*normal(2), shp_jump.d0);

      //   |                      |
      // + |  < q >, || u_i || n  |
      //   |                      |
      assembler.template Vector<Pres>(shp_mean.d0, -timefacfac*inflowvel*(veljump(0)*normal(0)
                                                                        +veljump(1)*normal(1)
                                                                        +veljump(2)*normal(2)));

      //-----------------       |                |
      // inflow flux term RHS - |  < q >, j_u n  |
      //-----------------       |                |
      assembler.template Vector<Pres>(shp_mean.d0, timefacfac*inflowvel*ju);
#endif
#ifdef COMBUST_STRESS_NITSCHEVEL
    const double alphau = 10.0;
    //------------------------    |                              |
    // Nitsche term (velocity)  + |  \alpha_u || v ||, || Du ||  |
    //------------------------    |                              |

    assembler.template Matrix<Velx,Velx>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
    assembler.template Matrix<Vely,Vely>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);
    assembler.template Matrix<Velz,Velz>(shp_jump.d0, alphau*timefacfac, shp_jump.d0);

    //    |                               |
    //  - |  \alpha_u || v ||, || u_i ||  |
    //    |                               |

    assembler.template Vector<Velx>(shp_jump.d0, -alphau*timefacfac*veljump(0));
    assembler.template Vector<Vely>(shp_jump.d0, -alphau*timefacfac*veljump(1));
    assembler.template Vector<Velz>(shp_jump.d0, -alphau*timefacfac*veljump(2));

    //----------------------------    |                         |
    // Nitsche term (velocity) RHS  + |  \alpha_u || v ||, j_u  |
    //----------------------------    |                         |

    //      assembler.template Vector<Velx>(shp_jump.d0, alphau*timefacfac*ju(0));
    //      assembler.template Vector<Vely>(shp_jump.d0, alphau*timefacfac*ju(1));
    //      assembler.template Vector<Velz>(shp_jump.d0, alphau*timefacfac*ju(2));
    assembler.template Vector<Velx>(shp_jump.d0, alphau*timefacfac*ju(0));
    assembler.template Vector<Vely>(shp_jump.d0, alphau*timefacfac*ju(1));
    assembler.template Vector<Velz>(shp_jump.d0, alphau*timefacfac*ju(2));
#endif
    } // loop Gaussian points
  } // loop boundary integration cells

  return;
}

} // namespace COMBUST

#endif // define COMBUST_NORMAL_ENRICHMENT

#endif
