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
#include "combust3_sysmat_premixed_nitsche.H"
#include "combust3_sysmat_premixed_nitsche_normal.H"
#include "combust3_sysmat_premixed_stress.H"
#include "combust3_sysmat_premixed_stress_normal.H"
#include "combust3_sysmat_twophaseflow.H"
#include "combust3_error_analysis.H"
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
    M2& etensor,
    V3& ediscpres
)
{
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

#ifndef COMBUST_NORMAL_ENRICHMENT
  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Vely);
  const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velz);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz) and (numparamvelx == numparampres), "assumption violation");
#else
  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Velx);
  const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Vely);
  const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Velz);
  const size_t numparamveln = XFEM::NumParam<0,ASSTYPE>::get(dofman,XFEM::PHYSICS::Veln);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Pres);
  dsassert((numparamvelx == 8) and (numparamvely == 8) and (numparamvelz == 8), "assumption violation");
#endif
#ifndef COMBUST_NORMAL_ENRICHMENT
  const size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*numnode;
  if (numparamvelx > shpVecSize)
  {
    dserror("increase SizeFac for nodal unknowns");
  }
#else
  const size_t shpVecSizeVel = COMBUST::SizeFacVel<ASSTYPE>::fac*numnode;
  if (numparamvelx > shpVecSizeVel)
  {
    dserror("increase SizeFac for nodal unknowns for velocity");
  }
  const size_t shpVecSizePres = COMBUST::SizeFacPres<ASSTYPE>::fac*numnode;
  if (numparampres > shpVecSizePres)
  {
    dserror("increase SizeFac for nodal unknowns for pressure");
  }
#endif

#ifndef COMBUST_NORMAL_ENRICHMENT
  const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
  const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
  const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
  const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());
#else
  const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
  const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
  const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
  const std::vector<int>& velndof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Veln>());
  const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());
#endif

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
#ifdef COMBUST_NORMAL_ENRICHMENT
  for (size_t iparam=0; iparam<numparamveln; ++iparam)
  {
    evelnp(3,iparam) = mystate.velnp_[velndof[iparam]];
    if (mystate.instationary_)
    {
      eveln( 3,iparam) = mystate.veln_[ velndof[iparam]];
      evelnm(3,iparam) = mystate.velnm_[velndof[iparam]];
      // TODO@Florian ist das richtig, mit der Beschleunigung?
      eaccn( 3,iparam) = mystate.accn_[ velndof[iparam]];
    }
  }
#endif
  for (size_t iparam=0; iparam<numparampres; ++iparam)
    eprenp(iparam) = mystate.velnp_[presdof[iparam]];

  const bool epsilonele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonxx, 0) > 0);
  if (epsilonele_unknowns_present)
  {
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const size_t numparamepsilonxx = XFEM::NumParam<1,ASSTYPE>::get(dofman,XFEM::PHYSICS::Epsilonxx);
    const size_t numparamepsilonyy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonyy, 1);
    const size_t numparamepsilonzz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonzz, 1);
    const size_t numparamepsilonxy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonxy, 1);
    const size_t numparamepsilonxz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonxz, 1);
    const size_t numparamepsilonyz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Epsilonyz, 1);
    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
    const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    if (numparamepsilonxx > shpVecSizeStress)
    {
      dserror("increase SizeFac for stress unknowns");
    }
    const std::vector<int>& epsilonxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Epsilonxx>());
    const std::vector<int>& epsilonyydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Epsilonyy>());
    const std::vector<int>& epsilonzzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Epsilonzz>());
    const std::vector<int>& epsilonxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Epsilonxy>());
    const std::vector<int>& epsilonxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Epsilonxz>());
    const std::vector<int>& epsilonyzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Epsilonyz>());
    for (size_t iparam=0; iparam<numparamepsilonxx; ++iparam)   etensor(0,iparam) = mystate.velnp_[epsilonxxdof[iparam]];
    for (size_t iparam=0; iparam<numparamepsilonyy; ++iparam)   etensor(1,iparam) = mystate.velnp_[epsilonyydof[iparam]];
    for (size_t iparam=0; iparam<numparamepsilonzz; ++iparam)   etensor(2,iparam) = mystate.velnp_[epsilonzzdof[iparam]];
    for (size_t iparam=0; iparam<numparamepsilonxy; ++iparam)   etensor(3,iparam) = mystate.velnp_[epsilonxydof[iparam]];
    for (size_t iparam=0; iparam<numparamepsilonxz; ++iparam)   etensor(4,iparam) = mystate.velnp_[epsilonxzdof[iparam]];
    for (size_t iparam=0; iparam<numparamepsilonyz; ++iparam)   etensor(5,iparam) = mystate.velnp_[epsilonyzdof[iparam]];
  }
  const bool sigmaele_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxx, 0) > 0);
  if (sigmaele_unknowns_present)
  {
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const size_t numparamsigmaxx = XFEM::NumParam<1,ASSTYPE>::get(dofman,XFEM::PHYSICS::Sigmaxx);
    const size_t numparamsigmayy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayy, 1);
    const size_t numparamsigmazz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmazz, 1);
    const size_t numparamsigmaxy = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxy, 1);
    const size_t numparamsigmaxz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmaxz, 1);
    const size_t numparamsigmayz = XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Sigmayz, 1);
    const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
    const size_t shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
    if (numparamsigmaxx > shpVecSizeStress)
    {
      dserror("increase SizeFac for stress unknowns");
    }
    const std::vector<int>& sigmaxxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxx>());
    const std::vector<int>& sigmayydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayy>());
    const std::vector<int>& sigmazzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmazz>());
    const std::vector<int>& sigmaxydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxy>());
    const std::vector<int>& sigmaxzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmaxz>());
    const std::vector<int>& sigmayzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Sigmayz>());
    for (size_t iparam=0; iparam<numparamsigmaxx; ++iparam)   etensor(0,iparam) = mystate.velnp_[sigmaxxdof[iparam]];
    for (size_t iparam=0; iparam<numparamsigmayy; ++iparam)   etensor(1,iparam) = mystate.velnp_[sigmayydof[iparam]];
    for (size_t iparam=0; iparam<numparamsigmazz; ++iparam)   etensor(2,iparam) = mystate.velnp_[sigmazzdof[iparam]];
    for (size_t iparam=0; iparam<numparamsigmaxy; ++iparam)   etensor(3,iparam) = mystate.velnp_[sigmaxydof[iparam]];
    for (size_t iparam=0; iparam<numparamsigmaxz; ++iparam)   etensor(4,iparam) = mystate.velnp_[sigmaxzdof[iparam]];
    for (size_t iparam=0; iparam<numparamsigmayz; ++iparam)   etensor(5,iparam) = mystate.velnp_[sigmayzdof[iparam]];
  }
  const bool discpres_unknowns_present = (XFEM::getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 0) > 0);
  if (discpres_unknowns_present)
  {
    const size_t numparamdiscpres = XFEM::NumParam<1,ASSTYPE>::get(dofman,XFEM::PHYSICS::DiscPres);
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
  // remark: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];
}
}

namespace COMBUST
{
//! fill a number of (local) element arrays
template <DRT::Element::DiscretizationType DISTYPE,
          class M>
void fillElementGradPhi(
    const DRT::ELEMENTS::Combust3::MyState& mystate,
    M& egradphi)
{
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  unsigned ipos;
  for (size_t iparam=0; iparam<numnode; ++iparam)
  {
    ipos = iparam*3;
    egradphi(0, iparam) = mystate.gradphinp_[ipos  ];
    egradphi(1, iparam) = mystate.gradphinp_[ipos+1];
    egradphi(2, iparam) = mystate.gradphinp_[ipos+2];
  }
}
}


/*------------------------------------------------------------------------------------------------*
 | get material parameters (constant within the domain integration cell)              henke 06/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::GetMaterialParams(
    Teuchos::RCP<const MAT::Material> material, // pointer to material (list)
    const bool indomplus, // boolean indicating side of the interface
    double&    dens,      // density
    double&    dynvisc    // dynamic viscosity
)
{
  //----------------------
  // get the material type
  //----------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // set default id in list of materials
  int matid = -1;
  // check on which side of the interface the cell is located
  if(indomplus) // cell belongs to burnt domain
  {
    matid = matlist->MatID(0); // burnt material (first material in material list)
  }
  else // cell belongs to unburnt domain
  {
    matid = matlist->MatID(1); // unburnt material (second material in material list)
  }
  // get material from list of materials
  Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
  INPAR::MAT::MaterialType mattype = matptr->MaterialType();

  // choose from different materials
  switch(mattype)
  {
  //--------------------------------------------------------
  // Newtonian fluid for incompressible flow (standard case)
  //--------------------------------------------------------
  case INPAR::MAT::m_fluid:
  {
    const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
    // get the dynamic viscosity \nu
    dynvisc = mat->Viscosity();
    // get the density \rho^{n+1}
    dens = mat->Density();
    break;
  }
  //------------------------------------------------
  // different types of materials (to be added here)
  //------------------------------------------------
  default:
    dserror("material type not supported");
  }

  // security check
  if (dens < 0 or dynvisc < 0)
    dserror("material parameters could not be determined");
  return;
}


/*------------------------------------------------------------------------------------------------*
 | get material parameters for both domains                                           henke 08/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::GetMaterialParams(
    Teuchos::RCP<const MAT::Material> material, // pointer to material (list)
    double&    dens_plus,    // density in "plus domain"
    double&    dynvisc_plus, // dynamic viscosity in "plus domain"
    double&    dens_minus,   // density in "minus domain"
    double&    dynvisc_minus // dynamic viscosity in "minus domain"
)
{
  //----------------------
  // get the material type
  //----------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // set default id in list of materials
  int matid = -1;

  // get material for both sides of the interface ("plus" and "minus" domain)
  for (int matcount=0;matcount<2;matcount++)
  {
    // get ID of material
    // matcount==0: material in burnt domain   (first material in material list)
    // matcount==1: material in unburnt domain (second material in material list)
    matid = matlist->MatID(matcount);

    // get material from list of materials
    Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
    INPAR::MAT::MaterialType mattype = matptr->MaterialType();

    // choose from different materials
    switch(mattype)
    {
    //--------------------------------------------------------
    // Newtonian fluid for incompressible flow (standard case)
    //--------------------------------------------------------
    case INPAR::MAT::m_fluid:
    {
      const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
      //--------------
      // plus material
      //--------------
      if (matcount==0)
      {
        // get the dynamic viscosity \nu
        dynvisc_plus = mat->Viscosity();
        // get the density \rho^{n+1}
        dens_plus = mat->Density();
      }
      //--------------
      // minus material
      //--------------
      if (matcount==1)
      {
        // get the dynamic viscosity \nu
        dynvisc_minus = mat->Viscosity();
        // get the density \rho^{n+1}
        dens_minus = mat->Density();
      }
      break;
    }
    //------------------------------------------------
    // different types of materials (to be added here)
    //------------------------------------------------
    default:
      dserror("material type not supported");
    }
  }

  // security check
  if ((dens_plus < 0  or dynvisc_plus < 0) and
      (dens_minus < 0 or dynvisc_minus < 0))
    dserror("material parameters could not be determined");

  return;
}

/*------------------------------------------------------------------------------------------------*
 | get surface tension coefficient                                               rasthofer  12/10 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::GetMaterialParams(
    Teuchos::RCP<const MAT::Material> material, // pointer to material (list)
    double&    surfacetensioncoeff       // surface tension coefficient
)
{
  //----------------------
  // get the material type
  //----------------------
#ifdef DEBUG
  // check if we really got a list of materials
  dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif
  // get material list for this element
  const MAT::MatList* matlist = static_cast<const MAT::MatList*>(material.get());
  // set default id in list of materials
  int matid = -1;
  matid = matlist->MatID(0);

//  // check on which side of the interface the cell is located
//  if(indomplus) // cell belongs to burnt domain
//  {
//    matid = matlist->MatID(0); // burnt material (first material in material list)
//  }
//  else // cell belongs to unburnt domain
//  {
//    matid = matlist->MatID(1); // unburnt material (second material in material list)
//  }
  // get material from list of materials
  Teuchos::RCP<const MAT::Material> matptr = matlist->MaterialById(matid);
  INPAR::MAT::MaterialType mattype = matptr->MaterialType();

  // choose from different materials
  switch(mattype)
  {
  //--------------------------------------------------------
  // Newtonian fluid for incompressible flow (standard case)
  //--------------------------------------------------------
  case INPAR::MAT::m_fluid:
  {
    const MAT::NewtonianFluid* mat = static_cast<const MAT::NewtonianFluid*>(matptr.get());
    // get the dynamic viscosity \nu
    surfacetensioncoeff = mat->Gamma();
//    // get the density \rho^{n+1}
//    dens = mat->Density();
    break;
  }
  //------------------------------------------------
  // different types of materials (to be added here)
  //------------------------------------------------
  default:
    dserror("material type not supported");
  }

  // security check
  if (surfacetensioncoeff < 0)
    dserror("material parameters could not be determined");
  return;
}

namespace COMBUST
{
/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void Sysmat(
    const DRT::ELEMENTS::Combust3*          ele,             ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih, ///< connection to the interface handler
    const XFEM::ElementDofManager&          dofman,          ///< dofmanager of the current element
    const DRT::ELEMENTS::Combust3::MyState& mystate,         ///< element state variables
    Epetra_SerialDenseMatrix&               estif,           ///< element matrix to calculate
    Epetra_SerialDenseVector&               eforce,          ///< element rhs to calculate
    Teuchos::RCP<const MAT::Material>       material,        ///< fluid material
    const INPAR::FLUID::TimeIntegrationScheme timealgo,      ///< time discretization type
    const double                            dt,              ///< delta t (time step size)
    const double                            theta,           ///< factor for one step theta scheme
    const bool                              newton,          ///< full Newton or fixed-point-like
    const bool                              pstab,           ///< flag for stabilisation
    const bool                              supg,            ///< flag for stabilisation
    const bool                              cstab,           ///< flag for stabilisation
    const INPAR::FLUID::TauType             tautype,         ///< stabilization parameter definition
    const bool                              instationary,    ///< switch between stationary and instationary formulation
    const INPAR::COMBUST::CombustionType    combusttype,     ///< switch for type of combusiton problem
    const double                            flamespeed,      ///<
    const double                            nitschevel,      ///<
    const double                            nitschepres,     ///<
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///<
    const bool                              connected_interface,
    const INPAR::COMBUST::VelocityJumpType  veljumptype,
    const INPAR::COMBUST::FluxJumpType      fluxjumptype,
    const bool                              smoothed_boundary_integration
)
{
  // initialize element stiffness matrix and force vector
  estif.Scale(0.0);
  eforce.Scale(0.0);

  const int NUMDOF = 4;

  COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

  // split velocity and pressure (and stress)
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const DRT::Element::DiscretizationType stressdistype = COMBUST::StressInterpolation3D<DISTYPE>::distype;
  const DRT::Element::DiscretizationType discpresdistype = COMBUST::DiscPressureInterpolation3D<DISTYPE>::distype;
  const int shpVecSizeStress = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<stressdistype>::numNodePerElement;
  const int shpVecSizeDiscPres = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<discpresdistype>::numNodePerElement;

#ifdef COMBUST_NORMAL_ENRICHMENT
  const size_t shpVecSizeVel = COMBUST::SizeFacVel<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const size_t shpVecSizePres = COMBUST::SizeFacPres<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  LINALG::Matrix<4,shpVecSizeVel> evelnp(true);
  LINALG::Matrix<4,shpVecSizeVel> eveln(true);
  LINALG::Matrix<4,shpVecSizeVel> evelnm(true);
  LINALG::Matrix<4,shpVecSizeVel> eaccn(true);
  LINALG::Matrix<shpVecSizePres,1> eprenp(true);
#else
  const int shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  LINALG::Matrix<3,shpVecSize> evelnp(true);
  LINALG::Matrix<3,shpVecSize> eveln(true);
  LINALG::Matrix<3,shpVecSize> evelnm(true);
  LINALG::Matrix<3,shpVecSize> eaccn(true);
  LINALG::Matrix<shpVecSize,1> eprenp(true);
#endif
  LINALG::Matrix<numnode,1> ephi(true);
  LINALG::Matrix<6,shpVecSizeStress> etensor(true);
  LINALG::Matrix<shpVecSizeDiscPres,1> ediscpres(true);

  COMBUST::fillElementUnknownsArrays<DISTYPE,ASSTYPE>(
      dofman, mystate, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etensor, ediscpres);

  switch(combusttype)
  {
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
#ifdef COMBUST_NITSCHE

    double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
    double ele_meas_minus = 0.0; // for different averages <> and {}
#ifndef COMBUST_NORMAL_ENRICHMENT
    COMBUST::SysmatDomainNitsche<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler,
        ele_meas_plus, ele_meas_minus);
#else
        LINALG::Matrix<3,numnode> egradphi(true);
     // boundary integrals are only added for intersected elements (fully enriched elements)
//    if (ele->Intersected() == true)
//    {
      COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi);
//    }
    COMBUST::SysmatDomainNitscheNormal<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, egradphi,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler,
        ele_meas_plus, ele_meas_minus);
#endif
#endif
#ifdef COMBUST_EPSPRES_BASED
#ifdef COMBUST_NORMAL_ENRICHMENT
    // get smoothed gradient of phi for surface tension applications
    LINALG::Matrix<3,numnode> egradphi(true);
     // boundary integrals are only added for intersected elements (fully enriched elements)
//    if (ele->Intersected() == true)
//    {
      COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi);
//    }

//cout << "phi gradient danach" << ele->Id() << " " << egradphi << endl;

    COMBUST::SysmatDomainStressNormal<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, egradphi, etensor, ediscpres,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);

#else
    COMBUST::SysmatDomainStress<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etensor, ediscpres,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);
#ifdef COMBUST_SIGMA_BASED
    // TODO: der aufruf ist doch der gleiche wie der darueber? -> COMBUST_SIGMA_BASED raus?
    COMBUST::SysmatDomainStress<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etensor, ediscpres,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler);
#endif
#endif
#endif

#ifndef COMBUST_DECOUPLEDXFEM
#ifdef COMBUST_NITSCHE
    // boundary integrals are added for intersected and touched elements (fully or partially enriched elements)
    if (ele->Intersected() == true || ele->Touched_Plus() == true )
    {
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      fillElementGradPhi<DISTYPE>(mystate, egradphi);
#ifndef COMBUST_NORMAL_ENRICHMENT
      COMBUST::SysmatBoundaryNitsche<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi, material, timealgo, dt, theta, assembler,
          flamespeed, nitschevel, nitschepres, ele_meas_plus, ele_meas_minus,
          surftensapprox, connected_interface, veljumptype,
          fluxjumptype, smoothed_boundary_integration);
#else
      COMBUST::SysmatBoundaryNitscheNormal<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi, material, timealgo, dt, theta, assembler,
          flamespeed, nitschevel, nitschepres, ele_meas_plus, ele_meas_minus,
          surftensapprox, connected_interface, veljumptype,
          fluxjumptype, smoothed_boundary_integration);
#endif
    }
#endif
    // boundary integrals are only added for intersected elements (fully enriched elements)
    if (ele->Intersected() == true)
    {
#ifdef COMBUST_STRESS_BASED
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi);

#ifdef COMBUST_EPSPRES_BASED
#ifdef COMBUST_NORMAL_ENRICHMENT
      COMBUST::SysmatBoundaryStressNormal<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi, etensor, ediscpres, material, timealgo, dt,
          theta, assembler, flamespeed);
#else
      COMBUST::SysmatBoundaryStress<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi, etensor, ediscpres, material, timealgo, dt,
          theta, assembler, flamespeed);
#endif
#ifdef COMBUST_SIGMA_BASED
      COMBUST::SysmatBoundarySigma<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi, etensor, ediscpres, material, timealgo, dt,
          theta, assembler, flamespeed);
#endif
#endif
#endif // COMBUST_STRESS_BASED
    }
#endif //COMBUST_DECOUPLEDXFEM
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflow:
  {
    double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
    double ele_meas_minus = 0.0; // for different averages <> and {}

    COMBUST::SysmatTwoPhaseFlow<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etensor,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler,
        ele_meas_plus, ele_meas_minus);
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  {
    double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
    double ele_meas_minus = 0.0; // for different averages <> and {}

    COMBUST::SysmatTwoPhaseFlow<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi, etensor,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler,
        ele_meas_plus, ele_meas_minus);

    // boundary integrals are added for intersected and touched elements (fully or partially enriched elements)
    if (ele->Intersected() == true || ele->Touched_Plus() == true )
    {
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi);

      COMBUST::SysmatBoundarySurfaceTension<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi, etensor,
          material, timealgo, dt, theta, assembler,
          flamespeed, nitschevel, nitschepres, ele_meas_plus, ele_meas_minus,
          surftensapprox, connected_interface, veljumptype,
          fluxjumptype, smoothed_boundary_integration);
    }
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
    double ele_meas_minus = 0.0; // for different averages <> and {}

    COMBUST::SysmatDomainNitsche<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelnp, eveln, evelnm, eaccn, eprenp, ephi,
        material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary, assembler,
        ele_meas_plus, ele_meas_minus);

    // boundary integrals are added for intersected and touched elements (fully or partially enriched elements)
    if (ele->Intersected() == true || ele->Touched_Plus() == true )
    {
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi);

      COMBUST::SysmatBoundaryNitsche<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelnp, eprenp, ephi, egradphi,
          material, timealgo, dt, theta, assembler,
          flamespeed, nitschevel, nitschepres, ele_meas_plus, ele_meas_minus,
          surftensapprox, connected_interface, veljumptype,
          fluxjumptype, smoothed_boundary_integration);
    }
  }
  break;
  default:
    dserror("unknown type of combustion problem");
  }

  //----------------------------------
  // symmetry check for element matrix
  // TODO: remove symmetry check
  //----------------------------------
//if(ele->Id()==2)
//{
//  //cout << endl << "stiffness matrix of element: " << ele->Id() << " columns " << estif.N() << " rows " << estif.M() << endl << endl;
//  bool sym = true;
//  int counter = 0;
//  for (int row=0; row<estif.M(); ++row)
//  {
//    for (int col=0; col<estif.N(); ++col)
//    {
//      //    cout << estif(row,col);
//      cout << " " << setw(4)<< std::setprecision(1) << estif(row,col);
//      double diff = estif(row,col)-estif(col,row);
//      if (!((diff>-1.0E-9) and (diff<+1.0E-9)))
//      {
//        cout << endl << counter << " difference of entry " << "Zeile "<< row << "Spalte " << col << " is not 0.0, but " << diff << endl;
//        sym = false;
//      }//std::setw(18) <<  << std::scientific
//      counter++;
//    }
//    cout << endl;
//  }
//  cout << "counter " << counter << endl;
//  if (sym==false) cout << "nicht symmetrisch"<<endl;
//  if (sym==true) cout << "symmetrisch"<<endl;
//  //dserror("STOP after middle element matrix");
//}



}
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::callSysmat(
    const XFEM::AssemblyType             assembly_type,
    const DRT::ELEMENTS::Combust3*       ele,
    const COMBUST::InterfaceHandleCombust* ih,
    const XFEM::ElementDofManager&       eleDofManager,
    const DRT::ELEMENTS::Combust3::MyState&  mystate,   ///< element state variables
    Epetra_SerialDenseMatrix&            estif,
    Epetra_SerialDenseVector&            eforce,
    Teuchos::RCP<const MAT::Material>    material,
    const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
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
    const double                         nitschepres,
    const INPAR::COMBUST::SurfaceTensionApprox  surftensapprox,
    const bool                             connected_interface,
    const INPAR::COMBUST::VelocityJumpType veljumptype,
    const INPAR::COMBUST::FluxJumpType     fluxjumptype,
    const bool                             smoothed_boundary_integration)
{
  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
          combusttype, flamespeed, nitschevel, nitschepres, surftensapprox,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration);
    break;
    case DRT::Element::hex20:
      COMBUST::Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration);
    break;
//    case DRT::Element::hex27:
//      COMBUST::Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox);
//    break;
//    case DRT::Element::tet4:
//      COMBUST::Sysmat<DRT::Element::tet4,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox);
//    break;
//    case DRT::Element::tet10:
//      COMBUST::Sysmat<DRT::Element::tet10,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox);
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
      COMBUST::Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
          combusttype, flamespeed, nitschevel, nitschepres, surftensapprox,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration);
    break;
    case DRT::Element::hex20:
      COMBUST::Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox,
            connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration);
    break;
//    case DRT::Element::hex27:
//      COMBUST::Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox);
//    break;
//    case DRT::Element::tet4:
//      COMBUST::Sysmat<DRT::Element::tet4,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox);
//    break;
//    case DRT::Element::tet10:
//      COMBUST::Sysmat<DRT::Element::tet10,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres,surftensapprox);
//    break;
    default:
      dserror("xfem_assembly Sysmat not templated yet");
    };
  }
}


namespace COMBUST
{
/*!
  Calculate Nitsche errors for Nitsche problem formulation
                                                                   schott Jun 15, 2010
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void NitscheErrors(
    ParameterList&                         eleparams,
    const INPAR::COMBUST::NitscheError&  NitscheErrorType,
    const DRT::ELEMENTS::Combust3*         ele,            ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih,   ///< connection to the interface handler
    const XFEM::ElementDofManager&                       dofman,         ///< dofmanager of the current element
    const DRT::ELEMENTS::Combust3::MyState&              mystate, ///< element state variables
    Teuchos::RCP<const MAT::Material>                    material,       ///< fluid material
    const bool                                           smoothed_boundary_integration
)
{
  const int NUMDOF = 4;

  // split velocity and pressure (and stress)
  const int shpVecSize       = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  LINALG::Matrix<shpVecSize,1> eprenp;
  LINALG::Matrix<3,shpVecSize> evelnp;
  LINALG::Matrix<numnode,1> ephi;


  //==============================================================================================================
  // fill velocity and pressure Arrays

  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Velx);
  const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Vely);
  const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Velz);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman,XFEM::PHYSICS::Pres);
  dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz) and (numparamvelx == numparampres), "assumption violation");

  if ((int)numparamvelx > shpVecSize)
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
  }
  for (size_t iparam=0; iparam<numparamvely; ++iparam)
  {
    evelnp(1,iparam) = mystate.velnp_[velydof[iparam]];
  }
  for (size_t iparam=0; iparam<numparamvelz; ++iparam)
  {
    evelnp(2,iparam) = mystate.velnp_[velzdof[iparam]];
  }
  for (size_t iparam=0; iparam<numparampres; ++iparam)
    eprenp(iparam) = mystate.velnp_[presdof[iparam]];

  // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
  // TODO: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];

  //=================================================================================================================
  double ele_meas_plus = 0.0;	// we need measure of element in plus domain and minus domain
  double ele_meas_minus = 0.0;	// for different averages <> and {}

  COMBUST::Nitsche_BuildDomainIntegratedErrors<DISTYPE,ASSTYPE,NUMDOF>(
      eleparams, NitscheErrorType, ele, ih, dofman, evelnp, eprenp, ephi, material, ele_meas_plus, ele_meas_minus);

  if (ele->Intersected() == true || ele->Touched_Plus() == true)
  {
    LINALG::Matrix<3,numnode> egradphi;
    egradphi.Clear();
    COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi);

    COMBUST::Nitsche_BuildBoundaryIntegratedErrors<DISTYPE,ASSTYPE,NUMDOF>(
        eleparams, NitscheErrorType, ele, ih, dofman, evelnp, eprenp, ephi, egradphi, material, ele_meas_plus, ele_meas_minus, smoothed_boundary_integration);
  }

  return;
}
}


/*----------------------------------------------------------------------*
 *----------------------------------------- schott Jun 15, 2010---------*/
void COMBUST::callNitscheErrors(
    ParameterList&                              eleparams,        ///< list of parameters
    const INPAR::COMBUST::NitscheError& NitscheErrorType, ///<
    const XFEM::AssemblyType                    assembly_type,    ///<
    const DRT::ELEMENTS::Combust3*              ele,              ///<
    const COMBUST::InterfaceHandleCombust*      ih,     ///<
    const XFEM::ElementDofManager&              eleDofManager,    ///<
    const DRT::ELEMENTS::Combust3::MyState&     mystate,          ///< element state variables
    Teuchos::RCP<const MAT::Material>           material,          ///<
    const bool                                  smoothed_boundary_integration
)
{
  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::NitscheErrors<DRT::Element::hex8,XFEM::standard_assembly>(
          eleparams, NitscheErrorType, ele, ih, eleDofManager, mystate, material,smoothed_boundary_integration);
    break;
//    case DRT::Element::hex20:
//      COMBUST::Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::hex27:
//      COMBUST::Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet4:
//      COMBUST::Sysmat<DRT::Element::tet4,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet10:
//      COMBUST::Sysmat<DRT::Element::tet10,XFEM::standard_assembly>(
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
      COMBUST::NitscheErrors<DRT::Element::hex8,XFEM::xfem_assembly>(
          eleparams, NitscheErrorType, ele, ih, eleDofManager, mystate, material,smoothed_boundary_integration);
    break;
//    case DRT::Element::hex20:
//      COMBUST::Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::hex27:
//      COMBUST::Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet4:
//      COMBUST::Sysmat<DRT::Element::tet4,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, cstab, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet10:
//      COMBUST::Sysmat<DRT::Element::tet10,XFEM::xfem_assembly>(
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
