/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat.cpp

\brief call system matrix formulation
       premixed combustion problem / two-phase flow problems

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
#include "combust3_sysmat_nitsche.H"
#include "combust3_sysmat_stress.H"
#include "combust3_sysmat_twophaseflow.H"
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


//! fill a number of (local) element arrays with unknown values from the (global) unknown vector given by the discretization
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M1, class V1, class M2, class V2, class V3>
void fillElementUnknownsArrays(
    const XFEM::ElementDofManager& dofman,
    const DRT::ELEMENTS::Combust3::MyState& mystate,
    M1& evelnp,
    M1& eveln,
    M1& evelnm,
    M1& eaccn,
    V1& eprenp,
    V2& ephi,
    M2& etau,
    V3& ediscpres
)
{
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Vely);
  const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velz);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz) and (numparamvelx == numparampres), "assumption violation");
  const size_t shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*numnode;
  if (numparamvelx > shpVecSize)
  {
    dserror("increase SizeFac for nodal unknowns");
  }

  const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
  const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
  const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
  const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());

  for (size_t iparam=0; iparam<numparamvelx; ++iparam)
  {
    evelnp(0,iparam) = mystate.velnp_[velxdof[iparam]];
    if (mystate.instationary_)
    {
      eveln( 0,iparam) = mystate.veln_[ velxdof[iparam]];
      evelnm(0,iparam) = mystate.velnm_[velxdof[iparam]];
      eaccn( 0,iparam) = mystate.accn_[ velxdof[iparam]];
    }
  }
  for (size_t iparam=0; iparam<numparamvely; ++iparam)
  {
    evelnp(1,iparam) = mystate.velnp_[velydof[iparam]];
    if (mystate.instationary_)
    {
      eveln( 1,iparam) = mystate.veln_[ velydof[iparam]];
      evelnm(1,iparam) = mystate.velnm_[velydof[iparam]];
      eaccn( 1,iparam) = mystate.accn_[ velydof[iparam]];
    }
  }
  for (size_t iparam=0; iparam<numparamvelz; ++iparam)
  {
    evelnp(2,iparam) = mystate.velnp_[velzdof[iparam]];
    if (mystate.instationary_)
    {
      eveln( 2,iparam) = mystate.veln_[ velzdof[iparam]];
      evelnm(2,iparam) = mystate.velnm_[velzdof[iparam]];
      eaccn( 2,iparam) = mystate.accn_[ velzdof[iparam]];
    }
  }
  for (size_t iparam=0; iparam<numparampres; ++iparam)
    eprenp(iparam) = mystate.velnp_[presdof[iparam]];

  const bool tauele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 0) > 0);
  if (tauele_unknowns_present)
  {
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const size_t numparamtauxx = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::Tauxx);
    const size_t numparamtauyy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyy, 1);
    const size_t numparamtauzz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauzz, 1);
    const size_t numparamtauxy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxy, 1);
    const size_t numparamtauxz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxz, 1);
    const size_t numparamtauyz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyz, 1);
    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
    const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    if (numparamtauxx > shpVecSizeStress)
    {
      dserror("increase SizeFac for stress unknowns");
    }
    const std::vector<int>& tauxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauxx>());
    const std::vector<int>& tauyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauyy>());
    const std::vector<int>& tauzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauzz>());
    const std::vector<int>& tauxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauxy>());
    const std::vector<int>& tauxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauxz>());
    const std::vector<int>& tauyzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Tauyz>());
    for (size_t iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = mystate.velnp_[tauxxdof[iparam]];
    for (size_t iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = mystate.velnp_[tauyydof[iparam]];
    for (size_t iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = mystate.velnp_[tauzzdof[iparam]];
    for (size_t iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = mystate.velnp_[tauxydof[iparam]];
    for (size_t iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = mystate.velnp_[tauxzdof[iparam]];
    for (size_t iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = mystate.velnp_[tauyzdof[iparam]];
  }
  const bool discpres_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 0) > 0);
  if (discpres_unknowns_present)
  {
    const size_t numparamdiscpres = XFEM::NumParam<1,ASSTYPE>::get(dofman, XFEM::PHYSICS::DiscPres);
    const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;
    const size_t shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;
    if (numparamdiscpres > shpVecSizeDiscPres)
    {
      dserror("increase SizeFac for stress unknowns");
    }
    const vector<int>& discpresdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::DiscPres>());
    for (std::size_t iparam=0; iparam<numparamdiscpres; ++iparam)   ediscpres(iparam) = mystate.velnp_[discpresdof[iparam]];
  }

  // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
  // TODO: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];
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
    const INPAR::FLUID::TauType       tautype,        ///< stabilization parameter definition
    const bool                        instationary,   ///< switch between stationary and instationary formulation
    const INPAR::COMBUST::CombustionType combusttype, ///< switch for type of combusiton problem
    const double                      flamespeed,     ///<
    const double                      nitschevel,     ///<
    const double                      nitschepres     ///<
)
{
  // initialize element stiffness matrix and force vector
  estif.Scale(0.0);
  eforce.Scale(0.0);

  const int NUMDOF = 4;

  LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

  // split velocity and pressure (and stress)
  const int shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;
  const int shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
  const int shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;
  LINALG::Matrix<shpVecSize,1> eprenp;
  LINALG::Matrix<3,shpVecSize> evelnp;
  LINALG::Matrix<3,shpVecSize> eveln;
  LINALG::Matrix<3,shpVecSize> evelnm;
  LINALG::Matrix<3,shpVecSize> eaccn;
  LINALG::Matrix<numnode,1> ephi;
  LINALG::Matrix<6,shpVecSizeStress> etau;
  LINALG::Matrix<shpVecSizeDiscPres,1> ediscpres;

  fillElementUnknownsArrays<DISTYPE,ASSTYPE>(dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp,
      ephi, etau, ediscpres);

  switch(combusttype)
  {
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
#ifdef COMBUST_NITSCHE
    COMBUST::SysmatDomainNitsche<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etau,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);
#endif
#ifdef COMBUST_STRESS_BASED
    COMBUST::SysmatDomainStress<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etau, ediscpres,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);
#endif

#ifndef COMBUST_DECOUPLEDXFEM
#ifdef COMBUST_NITSCHE
    // boundary integrals are only added for intersected elements (fully enriched elements)
    if (ele->Intersected() == true)
    {
      COMBUST::SysmatBoundaryNitsche<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, etau, material, timealgo, dt, theta, assembler,
          flamespeed,nitschevel,nitschepres);
    }
#endif
#ifdef COMBUST_STRESS_BASED
    // boundary integrals are only added for intersected elements (fully enriched elements)
    if (ele->Intersected() == true)
    {
#ifdef COMBUST_STRESS_BASED_DOUBLE_ONESIDED
      COMBUST::SysmatBoundaryStressDoubleOneSided<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, etau, ediscpres, material, timealgo, dt, theta, assembler,
          flamespeed);
#else
      COMBUST::SysmatBoundaryStress<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, etau, ediscpres, material, timealgo, dt, theta, assembler,
          flamespeed);
#endif
    }
#endif
#endif
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflow:
  {
    COMBUST::SysmatTwoPhaseFlow<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etau,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);
  }
  break;
  default:
    dserror("unknown type of combustion problem");
  }

  //cout << endl << "stiffness matrix of element: " << ele->Id() << " columns " << estif.N() << " rows " << estif.M() << endl << endl;
  int counter = 0;
  for (int row=0; row<estif.M(); ++row)
  {
    for (int col=0; col<estif.N(); ++col)
    {
      //    cout << estif(row,col);
      double diff = estif(row,col)-estif(col,row);
      if (!((diff>-1.0E-9) and (diff<+1.0E-9)))
      {
        //cout << counter << " difference of entry " << estif(row,col) << " is not 0.0, but " << diff << endl;
        //cout << "stiffness matrix entry " << estif(row,col) << " transpose " << estif(col,row) << endl;
      }
      counter++;
    }
  }
  //cout << "counter " << counter << endl;
  //dserror("STOP after first element matrix");
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
    const bool                           newton,
    const bool                           pstab,
    const bool                           supg,
    const bool                           cstab,
    const INPAR::FLUID::TauType          tautype,       ///< stabilization parameter definition
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
          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
          combusttype,flamespeed,nitschevel,nitschepres);
    break;
//    case DRT::Element::hex20:
//      Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::hex27:
//      Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet4:
//      Sysmat<DRT::Element::tet4,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet10:
//      Sysmat<DRT::Element::tet10,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
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
          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
          combusttype, flamespeed, nitschevel, nitschepres);
    break;
//    case DRT::Element::hex20:
//      Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::hex27:
//      Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet4:
//      Sysmat<DRT::Element::tet4,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet10:
//      Sysmat<DRT::Element::tet10,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
    default:
      dserror("xfem_assembly Sysmat not templated yet");
    };
  }
}


#endif
#endif
