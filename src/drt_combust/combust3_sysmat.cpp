/*----------------------------------------------------------------------*/
/*!
\file combust3_sysmat.cpp

\brief call system matrix formulation
       premixed combustion problem / two-phase flow problems

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/


#include "combust_interface.H"
#include "combust3_sysmat_premixed_nitsche.H"
#include "combust3_sysmat_twophaseflow.H"
#include "combust3_error_analysis.H"
#include "combust3_facemat_premixed_nitsche.H"
#include "combust3_facemat_twophaseflow.H"
#include "combust_defines.H"
#include "../drt_geometry/position_array.H"


namespace COMBUST
{
//! fill a number of (local) element arrays with unknown values from the (global) unknown vector given by the discretization
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M1, class V1, class M2, class V2, class V3>
void fillElementUnknownsArrays(
    const XFEM::ElementDofManager& dofman,
    const DRT::ELEMENTS::Combust3::MyState& mystate,
    M1& evelaf,
    M1& eveln,
    M1& evelnm,
    M1& eaccn,
    M1& eaccam,
    M1& fsevelaf,
    V1& epreaf,
    V2& ephi,
    M2& etensor,
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
  // since partial enrichments of single fields are allowed now,
  // this assumption is not valid any more
  //dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz) and (numparamvelx == numparampres), "assumption violation");

  const size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*numnode;
  if (numparamvelx > shpVecSize)
  {
    dserror("increase SizeFac for nodal unknowns");
  }

  const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
  const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
  const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
  const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());

  if (mystate.instationary_)
  {
    for (size_t iparam=0; iparam<numparamvelx; ++iparam)
    {
      eveln( 0,iparam) = mystate.veln_[ velxdof[iparam]];
      evelnm(0,iparam) = mystate.velnm_[velxdof[iparam]];
      eaccn( 0,iparam) = mystate.accn_[ velxdof[iparam]];
    }
    for (size_t iparam=0; iparam<numparamvely; ++iparam)
    {
      eveln( 1,iparam) = mystate.veln_[ velydof[iparam]];
      evelnm(1,iparam) = mystate.velnm_[velydof[iparam]];
      eaccn( 1,iparam) = mystate.accn_[ velydof[iparam]];
    }
    for (size_t iparam=0; iparam<numparamvelz; ++iparam)
    {
      eveln( 2,iparam) = mystate.veln_[ velzdof[iparam]];
      evelnm(2,iparam) = mystate.velnm_[velzdof[iparam]];
      eaccn( 2,iparam) = mystate.accn_[ velzdof[iparam]];
    }
    if(mystate.genalpha_)
    {
      for (size_t iparam=0; iparam<numparamvelx; ++iparam)
        evelaf(0,iparam) = mystate.velaf_[velxdof[iparam]];
      for (size_t iparam=0; iparam<numparamvely; ++iparam)
        evelaf(1,iparam) = mystate.velaf_[velydof[iparam]];
      for (size_t iparam=0; iparam<numparamvelz; ++iparam)
        evelaf(2,iparam) = mystate.velaf_[velzdof[iparam]];

      for (size_t iparam=0; iparam<numparamvelx; ++iparam)
        eaccam(0,iparam) = mystate.accam_[velxdof[iparam]];
      for (size_t iparam=0; iparam<numparamvely; ++iparam)
        eaccam(1,iparam) = mystate.accam_[velydof[iparam]];
      for (size_t iparam=0; iparam<numparamvelz; ++iparam)
        eaccam(2,iparam) = mystate.accam_[velzdof[iparam]];
    }
    else
    {
      for (size_t iparam=0; iparam<numparamvelx; ++iparam)
        evelaf(0,iparam) = mystate.velnp_[velxdof[iparam]];
      for (size_t iparam=0; iparam<numparamvely; ++iparam)
        evelaf(1,iparam) = mystate.velnp_[velydof[iparam]];
      for (size_t iparam=0; iparam<numparamvelz; ++iparam)
        evelaf(2,iparam) = mystate.velnp_[velzdof[iparam]];
    }
  }
  else
  {
    for (size_t iparam=0; iparam<numparamvelx; ++iparam)
      evelaf(0,iparam) = mystate.velnp_[velxdof[iparam]];
    for (size_t iparam=0; iparam<numparamvely; ++iparam)
      evelaf(1,iparam) = mystate.velnp_[velydof[iparam]];
    for (size_t iparam=0; iparam<numparamvelz; ++iparam)
      evelaf(2,iparam) = mystate.velnp_[velzdof[iparam]];
  }

  if (mystate.genalpha_)
  {
    for (size_t iparam=0; iparam<numparampres; ++iparam)
      epreaf(iparam) = mystate.velaf_[presdof[iparam]];
  }
  else
  {
    for (size_t iparam=0; iparam<numparampres; ++iparam)
      epreaf(iparam) = mystate.velnp_[presdof[iparam]];
  }

  if(mystate.avm3_)
  {
    for (size_t iparam=0; iparam<numparamvelx; ++iparam)
      fsevelaf( 0,iparam) = mystate.fsvelaf_[ velxdof[iparam]];
    for (size_t iparam=0; iparam<numparamvely; ++iparam)
      fsevelaf( 1,iparam) = mystate.fsvelaf_[ velydof[iparam]];
    for (size_t iparam=0; iparam<numparamvelz; ++iparam)
      fsevelaf( 2,iparam) = mystate.fsvelaf_[ velzdof[iparam]];
  }

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
    if (mystate.genalpha_)
    {
      for (size_t iparam=0; iparam<numparamepsilonxx; ++iparam)   etensor(0,iparam) = mystate.velaf_[epsilonxxdof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonyy; ++iparam)   etensor(1,iparam) = mystate.velaf_[epsilonyydof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonzz; ++iparam)   etensor(2,iparam) = mystate.velaf_[epsilonzzdof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonxy; ++iparam)   etensor(3,iparam) = mystate.velaf_[epsilonxydof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonxz; ++iparam)   etensor(4,iparam) = mystate.velaf_[epsilonxzdof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonyz; ++iparam)   etensor(5,iparam) = mystate.velaf_[epsilonyzdof[iparam]];
    }
    else
    {
      for (size_t iparam=0; iparam<numparamepsilonxx; ++iparam)   etensor(0,iparam) = mystate.velnp_[epsilonxxdof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonyy; ++iparam)   etensor(1,iparam) = mystate.velnp_[epsilonyydof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonzz; ++iparam)   etensor(2,iparam) = mystate.velnp_[epsilonzzdof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonxy; ++iparam)   etensor(3,iparam) = mystate.velnp_[epsilonxydof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonxz; ++iparam)   etensor(4,iparam) = mystate.velnp_[epsilonxzdof[iparam]];
      for (size_t iparam=0; iparam<numparamepsilonyz; ++iparam)   etensor(5,iparam) = mystate.velnp_[epsilonyzdof[iparam]];
    }
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
    if (mystate.genalpha_)
    {
      for (size_t iparam=0; iparam<numparamsigmaxx; ++iparam)   etensor(0,iparam) = mystate.velaf_[sigmaxxdof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmayy; ++iparam)   etensor(1,iparam) = mystate.velaf_[sigmayydof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmazz; ++iparam)   etensor(2,iparam) = mystate.velaf_[sigmazzdof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmaxy; ++iparam)   etensor(3,iparam) = mystate.velaf_[sigmaxydof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmaxz; ++iparam)   etensor(4,iparam) = mystate.velaf_[sigmaxzdof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmayz; ++iparam)   etensor(5,iparam) = mystate.velaf_[sigmayzdof[iparam]];
    }
    else
    {
      for (size_t iparam=0; iparam<numparamsigmaxx; ++iparam)   etensor(0,iparam) = mystate.velnp_[sigmaxxdof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmayy; ++iparam)   etensor(1,iparam) = mystate.velnp_[sigmayydof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmazz; ++iparam)   etensor(2,iparam) = mystate.velnp_[sigmazzdof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmaxy; ++iparam)   etensor(3,iparam) = mystate.velnp_[sigmaxydof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmaxz; ++iparam)   etensor(4,iparam) = mystate.velnp_[sigmaxzdof[iparam]];
      for (size_t iparam=0; iparam<numparamsigmayz; ++iparam)   etensor(5,iparam) = mystate.velnp_[sigmayzdof[iparam]];
    }
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
    const std::vector<int>& discpresdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::DiscPres>());
    for (std::size_t iparam=0; iparam<numparamdiscpres; ++iparam)   ediscpres(iparam) = mystate.velnp_[discpresdof[iparam]];
  }

  // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
  // remark: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];
}

//! fill a number of (local) element arrays with unknown values from the (global) unknown vector given by the discretization
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M, class V>
void fillElementUnknownsArrays(
    const XFEM::ElementDofManager& dofman,
    const DRT::ELEMENTS::Combust3::MyStateSurface& mystate,
    M& evelaf,
    V& ephi
)
{
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Vely);
  const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velz);
  dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz), "assumption violation");

  const size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*numnode;
  if (numparamvelx > shpVecSize)
    dserror("increase SizeFac for nodal unknowns");

  const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
  const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
  const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());

  if (mystate.instationary_)
  {
    if(mystate.genalpha_)
    {
      for (size_t iparam=0; iparam<numparamvelx; ++iparam)
        evelaf(0,iparam) = mystate.velaf_[velxdof[iparam]];
      for (size_t iparam=0; iparam<numparamvely; ++iparam)
        evelaf(1,iparam) = mystate.velaf_[velydof[iparam]];
      for (size_t iparam=0; iparam<numparamvelz; ++iparam)
        evelaf(2,iparam) = mystate.velaf_[velzdof[iparam]];
    }
    else
    {
      for (size_t iparam=0; iparam<numparamvelx; ++iparam)
        evelaf(0,iparam) = mystate.velnp_[velxdof[iparam]];
      for (size_t iparam=0; iparam<numparamvely; ++iparam)
        evelaf(1,iparam) = mystate.velnp_[velydof[iparam]];
      for (size_t iparam=0; iparam<numparamvelz; ++iparam)
        evelaf(2,iparam) = mystate.velnp_[velzdof[iparam]];
    }
  }
  else
  {
    for (size_t iparam=0; iparam<numparamvelx; ++iparam)
      evelaf(0,iparam) = mystate.velnp_[velxdof[iparam]];
    for (size_t iparam=0; iparam<numparamvely; ++iparam)
      evelaf(1,iparam) = mystate.velnp_[velydof[iparam]];
    for (size_t iparam=0; iparam<numparamvelz; ++iparam)
      evelaf(2,iparam) = mystate.velnp_[velzdof[iparam]];
  }

  // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
  // remark: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];
}


//! fill a number of (local) element arrays with unknown values from the (global) unknown vector given by the discretization
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE,
          class M, class V1, class V2>
void fillElementUnknownsArrays(
    const XFEM::ElementDofManager& dofman,
    const std::vector<double>& velnp,
    const std::vector<double>& velaf,
    const std::vector<double>& phinp,
    const bool instationary,
    const bool gen_alpha,
    M& evelaf,
    V1& epreaf,
    V2& ephi
)
{
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
  //const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
  const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);
  const size_t numparamvely = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Vely);
  const size_t numparamvelz = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velz);
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);
  // since partial enrichments of single fields are allowed now,
  // this assumption is not valid any more
  //dsassert((numparamvelx == numparamvely) and (numparamvelx == numparamvelz) and (numparamvelx == numparampres), "assumption violation");

  const size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*numnode;
  if (numparamvelx > shpVecSize)
  {
    dserror("increase SizeFac for nodal unknowns");
  }

  const std::vector<int>& velxdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velx>());
  const std::vector<int>& velydof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Vely>());
  const std::vector<int>& velzdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Velz>());
  const std::vector<int>& presdof(dofman.LocalDofPosPerField<XFEM::PHYSICS::Pres>());

  if(gen_alpha)
  {
    for (size_t iparam=0; iparam<numparamvelx; ++iparam)
      evelaf(0,iparam) = velaf[velxdof[iparam]];
    for (size_t iparam=0; iparam<numparamvely; ++iparam)
      evelaf(1,iparam) = velaf[velydof[iparam]];
    for (size_t iparam=0; iparam<numparamvelz; ++iparam)
      evelaf(2,iparam) = velaf[velzdof[iparam]];
  }
  else
  {
    for (size_t iparam=0; iparam<numparamvelx; ++iparam)
      evelaf(0,iparam) = velnp[velxdof[iparam]];
    for (size_t iparam=0; iparam<numparamvely; ++iparam)
      evelaf(1,iparam) = velnp[velydof[iparam]];
    for (size_t iparam=0; iparam<numparamvelz; ++iparam)
      evelaf(2,iparam) = velnp[velzdof[iparam]];
  }

  if (gen_alpha)
  {
    for (size_t iparam=0; iparam<numparampres; ++iparam)
      epreaf(iparam) = velaf[presdof[iparam]];
  }
  else
  {
    for (size_t iparam=0; iparam<numparampres; ++iparam)
      epreaf(iparam) = velnp[presdof[iparam]];
  }

  // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
  // remark: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam,0) = phinp[iparam];

  return;
}
}

namespace COMBUST
{
//! fill a number of (local) element arrays
template <DRT::Element::DiscretizationType DISTYPE,
          class M1,class M2, class V>
void fillElementGradPhi(
    const DRT::ELEMENTS::Combust3::MyState& mystate,
    M1& egradphi,
    M2& egradphi2,
    V& ecurv)
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

  unsigned ipos2;
  for (size_t iparam=0; iparam<numnode; ++iparam)
  {
    ipos2 = iparam*9;
    egradphi2(0, iparam) = mystate.gradphi2np_[ipos2  ];
    egradphi2(1, iparam) = mystate.gradphi2np_[ipos2+1];
    egradphi2(2, iparam) = mystate.gradphi2np_[ipos2+2];
    egradphi2(3, iparam) = mystate.gradphi2np_[ipos2+3];
    egradphi2(4, iparam) = mystate.gradphi2np_[ipos2+4];
    egradphi2(5, iparam) = mystate.gradphi2np_[ipos2+5];
    egradphi2(6, iparam) = mystate.gradphi2np_[ipos2+6];
    egradphi2(7, iparam) = mystate.gradphi2np_[ipos2+7];
    egradphi2(8, iparam) = mystate.gradphi2np_[ipos2+8];
  }

  for (size_t iparam=0; iparam<numnode; ++iparam)
  {
    ecurv(iparam) = mystate.curv_[iparam];
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


/*------------------------------------------------------------------------------------------------*
 | blend material parameters smoothly from original to target parameters              henke 10/11 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::BlendMaterial(
    const DRT::Element* ele,
    const bool          indomainplus,
    const double        time,
    double&             denstarget,    // target density
    double&             dynvisctarget, // target dynamic viscosity
    const double        densplus,      // plus density
    const double        dynviscplus,   // minus viscosity
    const double        densminus,     // plus density
    const double        dynviscminus   // minus viscosity
)
{
  std::vector<DRT::Condition*> cond;
  DRT::UTILS::FindElementConditions(ele, "BlendMaterial", cond);

  if (cond.size()>=1)
  {
    // get the domain to blend
    bool blenddomain = true;
    const std::string* sblenddomain = cond[0]->Get<std::string>("domain");
    if (sblenddomain->compare("plus") == 0)
      blenddomain = true;
    else if (sblenddomain->compare("minus") == 0)
      blenddomain = false;
    else
      dserror("Invalid domain in blending material");

    if (indomainplus == blenddomain)
    {
      // find out whether there is a time curve
      const std::vector<int>* curve = cond[0]->Get<std::vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];

      // initialisation
      double curvefac = 0.0;
      if (curvenum >= 0) // yes, we have a timecurve
      {
        // time factor for the intermediate step
        if(time >= 0.0)
          curvefac = DRT::Problem::Instance()->Funct(curvenum).EvaluateTime(time);
        else
          dserror("Negative time value in blending material: time = %f",time);
      }
      else // we do not have a timecurve --- timefactors are constant equal 1
        dserror("Blending of materials without time curve does not make sense. Provide a time curve (from 1 at t=0 to 0 at t=tend) in your input file!");

      double fac = -1.0;
      // get switch from the condition
      const std::vector<int>* onoff = cond[0]->Get<std::vector<int> > ("onoff");
      if (onoff)
        fac = (*onoff)[0]*curvefac;
      else
        dserror("could not read switch in blending material");

      if (blenddomain)
      {
        denstarget = densminus*fac + densplus*(1.0-fac);
        dynvisctarget = dynviscminus*fac + dynviscplus*(1.0-fac);
      }
      else
      {
        denstarget = densplus*fac + densminus*(1.0-fac);
        dynvisctarget = dynviscplus*fac + dynviscminus*(1.0-fac);
      }
    }
    else
    {
      if (indomainplus)
      {
        denstarget = densplus;
        dynvisctarget = dynviscplus;
      }
      else
      {
        denstarget = densminus;
        dynvisctarget = dynviscminus;
      }
    }
  }
  else
  {
    if (indomainplus)
    {
      denstarget = densplus;
      dynvisctarget = dynviscplus;
    }
    else
    {
      denstarget = densminus;
      dynvisctarget = dynviscminus;
    }
  }
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
    const double                            time,            ///< current time step
    const double                            dt,              ///< delta t (time step size)
    const double                            theta,           ///< factor for one step theta scheme
    const double                            ga_alphaF,
    const double                            ga_alphaM,
    const double                            ga_gamma,
    const bool                              newton,          ///< full Newton or fixed-point-like
    const bool                              pstab,           ///< flag for stabilisation
    const bool                              supg,            ///< flag for stabilisation
    const bool                              graddiv,           ///< flag for stabilisation
    const INPAR::FLUID::TauType             tautype,         ///< stabilization parameter definition
    const bool                              instationary,    ///< switch between stationary and instationary formulation
    const bool                              genalpha,
    const INPAR::COMBUST::CombustionType    combusttype,     ///< switch for type of combusiton problem
    const double                            flamespeed,      ///<
    const double                            marksteinlength, ///<
    const double                            nitschevel,      ///<
    const double                            nitschepres,     ///<
    const INPAR::COMBUST::SurfaceTensionApprox surftensapprox, ///< type of surface tension approximation
    const double                            variablesurftens,  ///< variable surface tesnion approximation
    const bool                              second_deriv, ///< for curvature calculation
    const bool                              connected_interface,
    const INPAR::COMBUST::VelocityJumpType  veljumptype,
    const INPAR::COMBUST::FluxJumpType      fluxjumptype,
    const bool                              smoothed_boundary_integration,
    const INPAR::COMBUST::WeightType        weighttype,
    const bool                              nitsche_convflux,
    const bool                              nitsche_convstab,
    const bool                              nitsche_convpenalty,
    const bool                              nitsche_mass,
    const INPAR::FLUID::TurbModelAction    turbmodel,
    const INPAR::FLUID::FineSubgridVisc    fssgv,
    const double                           Cs,
    const double                           Csgs,
    const double                           alpha,
    const bool                             CalcN,
    const double                           N,
    const INPAR::FLUID::RefVelocity        refvel,
    const INPAR::FLUID::RefLength          reflength,
    const double                           c_nu,
    const bool                             near_wall_limit,
    const bool                             B_gp,
    const bool                             mfs_is_conservative
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

  const int shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  LINALG::Matrix<3,shpVecSize> evelaf(true);
  LINALG::Matrix<3,shpVecSize> eveln(true);
  LINALG::Matrix<3,shpVecSize> evelnm(true);
  LINALG::Matrix<3,shpVecSize> eaccn(true);
  LINALG::Matrix<3,shpVecSize> eaccam(true);
  LINALG::Matrix<shpVecSize,1> epreaf(true);
  LINALG::Matrix<3,shpVecSize> fsevelaf(true);
  LINALG::Matrix<numnode,1> ephi(true);
  LINALG::Matrix<6,shpVecSizeStress> etensor(true);
  LINALG::Matrix<shpVecSizeDiscPres,1> ediscpres(true);

  COMBUST::fillElementUnknownsArrays<DISTYPE,ASSTYPE>(
      dofman, mystate, evelaf, eveln, evelnm, eaccn, eaccam, fsevelaf, epreaf, ephi, etensor, ediscpres);

  switch(combusttype)
  {
  case INPAR::COMBUST::combusttype_premixedcombustion:
  {
    double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
    double ele_meas_minus = 0.0; // for different averages <> and {}
      COMBUST::SysmatDomainNitsche<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelaf, eveln, evelnm, eaccn, eaccam, fsevelaf, epreaf, ephi,
          material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha, assembler,
          ele_meas_plus, ele_meas_minus,
          turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);

#ifndef COMBUST_DECOUPLEDXFEM
    // boundary integrals are added for intersected and touched elements (fully or partially enriched elements)
    if (ele->Intersected() )
    {
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      LINALG::Matrix<9,numnode> egradphi2;
      egradphi2.Clear();
      LINALG::Matrix<numnode,1> ecurv;
      ecurv.Clear();
      if (smoothed_boundary_integration)
        fillElementGradPhi<DISTYPE>(mystate, egradphi,egradphi2, ecurv);
      COMBUST::SysmatBoundaryNitsche<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelaf, epreaf, ephi, egradphi,egradphi2, ecurv, material, time,
          timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, assembler,
          flamespeed, marksteinlength, nitschevel, nitschepres, ele_meas_plus, ele_meas_minus,
          surftensapprox, variablesurftens, second_deriv, connected_interface, veljumptype,
          fluxjumptype, smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass);
    }
#endif //COMBUST_DECOUPLEDXFEM
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflow:
  {
    if (fssgv!=INPAR::FLUID::no_fssgv
        or turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
      dserror("Turbulence models not yet supported for kink enrichments");

    COMBUST::SysmatTwoPhaseFlow<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelaf, eveln, evelnm, eaccn, eaccam, epreaf, ephi, etensor,
        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary,
        genalpha, assembler);
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  {
    if (fssgv!=INPAR::FLUID::no_fssgv
        or turbmodel == INPAR::FLUID::multifractal_subgrid_scales)
      dserror("Turbulence models not yet supported for kink enrichments");

    COMBUST::SysmatTwoPhaseFlow<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelaf, eveln, evelnm, eaccn, eaccam, epreaf, ephi, etensor,
        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary,
        genalpha, assembler);

    // boundary integrals are added for intersected and touched elements (fully or partially enriched elements)
    if (ele->Intersected())
    {
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      LINALG::Matrix<9,numnode> egradphi2;
      egradphi2.Clear();
      LINALG::Matrix<numnode,1> ecurv;
      ecurv.Clear();
      if (smoothed_boundary_integration)
        COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi,egradphi2, ecurv);

      COMBUST::SysmatBoundarySurfaceTension<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelaf, epreaf, ephi, egradphi,egradphi2, ecurv, etensor,
          material, timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, assembler,
          flamespeed, nitschevel, nitschepres,
          surftensapprox, variablesurftens, second_deriv, connected_interface, smoothed_boundary_integration);
    }
  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
    double ele_meas_minus = 0.0; // for different averages <> and {}

    COMBUST::SysmatDomainNitsche<DISTYPE,ASSTYPE,NUMDOF>(
        ele, ih, dofman, evelaf, eveln, evelnm, eaccn, eaccam, fsevelaf, epreaf, ephi,
        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha, assembler,
        ele_meas_plus, ele_meas_minus,
        turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);

    // boundary integrals are added for intersected and touched elements (fully or partially enriched elements)
    if (ele->Intersected())
    {
      // get smoothed gradient of phi for surface tension applications
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      LINALG::Matrix<9,numnode> egradphi2;
      egradphi2.Clear();
      LINALG::Matrix<numnode,1> ecurv;
      ecurv.Clear();
      if (smoothed_boundary_integration)
        COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi, egradphi2, ecurv);

      COMBUST::SysmatBoundaryNitsche<DISTYPE,ASSTYPE,NUMDOF>(
          ele, ih, dofman, evelaf, epreaf, ephi, egradphi,egradphi2, ecurv,
          material, time, timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, assembler,
          flamespeed, marksteinlength, nitschevel, nitschepres, ele_meas_plus, ele_meas_minus,
          surftensapprox, variablesurftens, second_deriv, connected_interface, INPAR::COMBUST::vel_jump_none,
          INPAR::COMBUST::flux_jump_surface_tension, smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass);
    }
  }
  break;
  default:
    dserror("unknown type of combustion problem");
  }
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
    const double                         time,          ///< current time step
    const double                         dt,            ///< delta t (time step size)
    const double                         theta,         ///< factor for one step theta scheme
    const double                         ga_alphaF,
    const double                         ga_alphaM,
    const double                         ga_gamma,
    const bool                           newton,
    const bool                           pstab,
    const bool                           supg,
    const bool                           graddiv,
    const INPAR::FLUID::TauType          tautype,       ///< stabilization parameter definition
    const bool                           instationary,
    const bool                           genalpha,
    const INPAR::COMBUST::CombustionType combusttype,
    const double                         flamespeed,
    const double                         marksteinlength,
    const double                         nitschevel,
    const double                         nitschepres,
    const INPAR::COMBUST::SurfaceTensionApprox  surftensapprox,
    const double                           variablesurftens,
    const bool                             second_deriv,
    const bool                             connected_interface,
    const INPAR::COMBUST::VelocityJumpType veljumptype,
    const INPAR::COMBUST::FluxJumpType     fluxjumptype,
    const bool                             smoothed_boundary_integration,
    const INPAR::COMBUST::WeightType       weighttype,
    const bool                             nitsche_convflux,
    const bool                             nitsche_convstab,
    const bool                             nitsche_convpenalty,
    const bool                             nitsche_mass,
    const INPAR::FLUID::TurbModelAction    turbmodel,
    const INPAR::FLUID::FineSubgridVisc    fssgv,
    const double                           Cs,
    const double                           Csgs,
    const double                           alpha,
    const bool                             CalcN,
    const double                           N,
    const INPAR::FLUID::RefVelocity        refvel,
    const INPAR::FLUID::RefLength          reflength,
    const double                           c_nu,
    const bool                             near_wall_limit,
    const bool                             B_gp,
    const bool                             mfs_is_conservative
)
{
  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha,
          combusttype, flamespeed, marksteinlength, nitschevel, nitschepres, surftensapprox, variablesurftens, second_deriv,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass,
          turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);
    break;
    case DRT::Element::hex20:
      COMBUST::Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha,
          combusttype, flamespeed, marksteinlength, nitschevel, nitschepres,surftensapprox, variablesurftens, second_deriv,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass,
          turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);
    break;
    default:
      dserror("standard_assembly Sysmat not templated yet");
      break;
    };
  }
  else
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha,
          combusttype, flamespeed, marksteinlength, nitschevel, nitschepres, surftensapprox, variablesurftens, second_deriv,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass,
          turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);
    break;
    case DRT::Element::hex20:
      COMBUST::Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce,
          material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, newton, pstab, supg, graddiv, tautype, instationary, genalpha,
          combusttype, flamespeed, marksteinlength, nitschevel, nitschepres,surftensapprox, variablesurftens, second_deriv,
          connected_interface, veljumptype, fluxjumptype, smoothed_boundary_integration,
          weighttype,nitsche_convflux,nitsche_convstab,nitsche_convpenalty,nitsche_mass,
          turbmodel, fssgv, Cs, Csgs, alpha, CalcN, N, refvel, reflength, c_nu, near_wall_limit, B_gp, mfs_is_conservative);
    break;
    default:
      dserror("xfem_assembly Sysmat not templated yet");
      break;
    };
  }
}


namespace COMBUST
{
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void SysmatNeumannInflow(
    const DRT::ELEMENTS::Combust3*                  ele,
    const DRT::ELEMENTS::Combust3Surface*           elesurf,
    const XFEM::ElementDofManager&                  dofman,
    const DRT::ELEMENTS::Combust3::MyStateSurface&  mystate,       ///< element state variables
    Epetra_SerialDenseMatrix&                       estif,
    Epetra_SerialDenseVector&                       eforce,
    const GEO::BoundaryIntCells&                    surfintcells,
    Teuchos::RCP<const MAT::Material>               material,
    const INPAR::FLUID::TimeIntegrationScheme       timealgo,      ///< time discretization type
    const double                                    time,          ///< current time step
    const double                                    dt,            ///< delta t (time step size)
    const double                                    theta,         ///< factor for one step theta scheme
    const double                                    ga_gamma,
    const double                                    ga_alphaF,
    const double                                    ga_alphaM,
    const bool                                      newton,
    const bool                                      instationary,
    const bool                                      genalpha,
    const INPAR::COMBUST::CombustionType            combusttype
)
{
  //TEUCHOS_FUNC_TIME_MONITOR(" - evaluating - combustion sysmat - boundary");
  // initialize element stiffness matrix and force vector
  estif.Scale(0.0);
  eforce.Scale(0.0);

  const int NUMDOF = 4;
  COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

  // split velocity and pressure (and stress)
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  const int shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  LINALG::Matrix<3,shpVecSize> evelaf(true);

  LINALG::Matrix<numnode,1> ephi(true);

  COMBUST::fillElementUnknownsArrays<DISTYPE,ASSTYPE>(dofman, mystate, evelaf, ephi);

  switch(combusttype)
  {
  case INPAR::COMBUST::combusttype_premixedcombustion:
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    // time integration constant
    const double timefac = COMBUST::TIMEINT::ComputeTimeFac(timealgo, dt, theta, ga_alphaF, ga_alphaM, ga_gamma);

    // get node coordinates of the current element
    static LINALG::Matrix<3,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

    //---------------------------------------------------------------------------
    // get material parameters for all boundary integration cells of this element
    //---------------------------------------------------------------------------
#ifdef DEBUG
    // check if we really got a list of materials
    dsassert(material->MaterialType() == INPAR::MAT::m_matlist, "Material law is not of type m_matlist");
#endif

    // get number of parameters (dofs) for each field
    // remark: it is assumed that all fields are enriched -> equal for all velocity components and pressure
    const size_t numparamvelx = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Velx);

    //-------------------------------------
    // loop over boundary integration cells
    //-------------------------------------
    for (GEO::BoundaryIntCells::const_iterator cell = surfintcells.begin(); cell != surfintcells.end(); ++cell)
    {
      switch (cell->Shape())
      {
      case DRT::Element::tri3:
        COMBUST::Nitsche_SysmatNeumannInflow<DISTYPE,DRT::Element::tri3,ASSTYPE,NUMDOF>(
            assembler, ele, elesurf, dofman, *cell,
            DRT::UTILS::intrule_tri_37point, xyze,
            evelaf, ephi,
            numparamvelx, newton,
            material,
            time, timefac
        );
        break;
      case DRT::Element::quad4:
        COMBUST::Nitsche_SysmatNeumannInflow<DISTYPE,DRT::Element::quad4,ASSTYPE,NUMDOF>(
            assembler, ele, elesurf, dofman, *cell,
            DRT::UTILS::intrule_quad_25point, xyze,
            evelaf, ephi,
            numparamvelx, newton,
            material,
            time, timefac
        );
        break;
      default:
        dserror("invalid type of boundary integration cell");
      }
    } // loop boundary integration cells

  }
  break;
  case INPAR::COMBUST::combusttype_twophaseflow:
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  {
    dserror("Neumann inflow term only implemented for premixed combustion problems (jump enrichments)");
  }
  break;
  }

}
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::callSysmatNeumannInflow(
    const XFEM::AssemblyType                 assembly_type,
    const DRT::ELEMENTS::Combust3*           ele,
    const DRT::ELEMENTS::Combust3Surface*    elesurf,
    const XFEM::ElementDofManager&           dofman,
    const DRT::ELEMENTS::Combust3::MyStateSurface&  mystate,       ///< element state variables
    Epetra_SerialDenseMatrix&                estif,
    Epetra_SerialDenseVector&                eforce,
    const GEO::BoundaryIntCells&             surfintcells,
    Teuchos::RCP<const MAT::Material>    material,
    const INPAR::FLUID::TimeIntegrationScheme timealgo, ///< time discretization type
    const double                         time,          ///< current time step
    const double                         dt,            ///< delta t (time step size)
    const double                         theta,         ///< factor for one step theta scheme
    const double                         ga_gamma,
    const double                         ga_alphaF,
    const double                         ga_alphaM,
    const bool                           newton,
    const bool                           instationary,
    const bool                           genalpha,
    const INPAR::COMBUST::CombustionType combusttype
    )
{
  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::SysmatNeumannInflow<DRT::Element::hex8,XFEM::standard_assembly>(
          ele,
          elesurf,
          dofman,
          mystate,       ///< element state variables
          estif,
          eforce,
          surfintcells,
          material,
          timealgo, ///< time discretization type
          time,          ///< current time step
          dt,            ///< delta t (time step size)
          theta,         ///< factor for one step theta scheme
          ga_gamma,
          ga_alphaF,
          ga_alphaM,
          newton,
          instationary,
          genalpha,
          combusttype);
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
      COMBUST::SysmatNeumannInflow<DRT::Element::hex8,XFEM::xfem_assembly>(
          ele,
          elesurf,
          dofman,
          mystate,       ///< element state variables
          estif,
          eforce,
          surfintcells,
          material,
          timealgo, ///< time discretization type
          time,          ///< current time step
          dt,            ///< delta t (time step size)
          theta,         ///< factor for one step theta scheme
          ga_gamma,
          ga_alphaF,
          ga_alphaM,
          newton,
          instationary,
          genalpha,
          combusttype);
    break;
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
    Teuchos::ParameterList&                 eleparams,
    const INPAR::COMBUST::NitscheError&     NitscheErrorType,
    const DRT::ELEMENTS::Combust3*          ele,                          ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih,                           ///< connection to the interface handler
    const XFEM::ElementDofManager&          dofman,                       ///< dofmanager of the current element
    const DRT::ELEMENTS::Combust3::MyState& mystate,                      ///< element state variables
    Teuchos::RCP<const MAT::Material>       material,                     ///< fluid material
    const double                            time,                         ///< current time
    const bool                              smoothed_boundary_integration,
    const INPAR::COMBUST::CombustionType    combusttype
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
    evelnp(0,iparam) = mystate.velnp_[velxdof[iparam]];
  for (size_t iparam=0; iparam<numparamvely; ++iparam)
    evelnp(1,iparam) = mystate.velnp_[velydof[iparam]];
  for (size_t iparam=0; iparam<numparamvelz; ++iparam)
    evelnp(2,iparam) = mystate.velnp_[velzdof[iparam]];
  for (size_t iparam=0; iparam<numparampres; ++iparam)
    eprenp(iparam) = mystate.velnp_[presdof[iparam]];

  // copy element phi vector from std::vector (mystate) to LINALG::Matrix (ephi)
  // TODO: this is inefficient, but it is nice to have only fixed size matrices afterwards!
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];

  //=================================================================================================================
  double ele_meas_plus = 0.0;  // we need measure of element in plus domain and minus domain
  double ele_meas_minus = 0.0;  // for different averages <> and {}

  switch(combusttype)
  {
  case INPAR::COMBUST::combusttype_premixedcombustion:
  case INPAR::COMBUST::combusttype_twophaseflowjump:
  {
    COMBUST::Nitsche_BuildDomainIntegratedErrors<DISTYPE,ASSTYPE,NUMDOF>(
      eleparams, NitscheErrorType, ele, ih, dofman, evelnp, eprenp, ephi, material, time, ele_meas_plus, ele_meas_minus, false);

    if (ele->Intersected())
    {
      LINALG::Matrix<3,numnode> egradphi;
      egradphi.Clear();
      LINALG::Matrix<9,numnode> egradphi2;
      egradphi2.Clear();
      LINALG::Matrix<numnode,1> ecurv;
      ecurv.Clear();
      if (smoothed_boundary_integration)
        COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi, egradphi2, ecurv);

      COMBUST::Nitsche_BuildBoundaryIntegratedErrors<DISTYPE,ASSTYPE,NUMDOF>(
          eleparams, NitscheErrorType, ele, ih, dofman, evelnp, eprenp, ephi, egradphi, material, time, ele_meas_plus, ele_meas_minus, smoothed_boundary_integration, false);
    }

    break;
  }
  case INPAR::COMBUST::combusttype_twophaseflow:
  case INPAR::COMBUST::combusttype_twophaseflow_surf:
  {
    COMBUST::Nitsche_BuildDomainIntegratedErrors<DISTYPE,ASSTYPE,NUMDOF>(
      eleparams, NitscheErrorType, ele, ih, dofman, evelnp, eprenp, ephi, material, time, ele_meas_plus, ele_meas_minus, true);

    if (combusttype == INPAR::COMBUST::combusttype_twophaseflow_surf)
    {
      if (ele->Intersected())
      {
        LINALG::Matrix<3,numnode> egradphi;
        egradphi.Clear();
        LINALG::Matrix<9,numnode> egradphi2;
        egradphi2.Clear();
        LINALG::Matrix<numnode,1> ecurv;
        ecurv.Clear();
        if (smoothed_boundary_integration)
          COMBUST::fillElementGradPhi<DISTYPE>(mystate, egradphi,egradphi2, ecurv);

        COMBUST::Nitsche_BuildBoundaryIntegratedErrors<DISTYPE,ASSTYPE,NUMDOF>(
           eleparams, NitscheErrorType, ele, ih, dofman, evelnp, eprenp, ephi, egradphi, material, time, ele_meas_plus, ele_meas_minus, smoothed_boundary_integration, true);
      }
    }

    break;
  }
  default:
    dserror("unknown type of combustion problem");
  }

  return;
}
}


/*----------------------------------------------------------------------*
 *----------------------------------------- ----------------------------*/
void COMBUST::callNitscheErrors(
    Teuchos::ParameterList&                     eleparams,        ///< list of parameters
    const INPAR::COMBUST::NitscheError&         NitscheErrorType, ///<
    const XFEM::AssemblyType                    assembly_type,    ///<
    const DRT::ELEMENTS::Combust3*              ele,              ///<
    const COMBUST::InterfaceHandleCombust*      ih,               ///<
    const XFEM::ElementDofManager&              eleDofManager,    ///<
    const DRT::ELEMENTS::Combust3::MyState&     mystate,          ///< element state variables
    Teuchos::RCP<const MAT::Material>           material,         ///<
    const double                                time,             ///< current time
    const bool                                  smoothed_boundary_integration,
    const INPAR::COMBUST::CombustionType        combusttype
)
{
  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::NitscheErrors<DRT::Element::hex8,XFEM::standard_assembly>(
          eleparams, NitscheErrorType, ele, ih, eleDofManager, mystate, material, time, smoothed_boundary_integration, combusttype);
    break;
//    case DRT::Element::hex20:
//      COMBUST::Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::hex27:
//      COMBUST::Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet4:
//      COMBUST::Sysmat<DRT::Element::tet4,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet10:
//      COMBUST::Sysmat<DRT::Element::tet10,XFEM::standard_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
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
          eleparams, NitscheErrorType, ele, ih, eleDofManager, mystate, material,time, smoothed_boundary_integration, combusttype);
    break;
//    case DRT::Element::hex20:
//      COMBUST::Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::hex27:
//      COMBUST::Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet4:
//      COMBUST::Sysmat<DRT::Element::tet4,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
//    case DRT::Element::tet10:
//      COMBUST::Sysmat<DRT::Element::tet10,XFEM::xfem_assembly>(
//          ele, ih, eleDofManager, mystate, estif, eforce,
//          material, timealgo, dt, theta, newton, pstab, supg, graddiv, tautype, instationary,
//          combusttype, flamespeed, nitschevel, nitschepres);
//    break;
    default:
      dserror("xfem_assembly Sysmat not templated yet");
    };
  }
}



namespace COMBUST
{
/*!
  Integrate the pressure shapefunctions
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void IntegrateShape(
    const DRT::ELEMENTS::Combust3*          ele,             ///< the element those matrix is calculated
    const COMBUST::InterfaceHandleCombust*  ih,              ///< connection to the interface handler
    const XFEM::ElementDofManager&          dofman,          ///< dofmanager of the current element
    const DRT::ELEMENTS::Combust3::MyState& mystate,         ///< element state variables
    Epetra_SerialDenseMatrix&               estif,           ///< element matrix to calculate
    Epetra_SerialDenseVector&               eforce           ///< element rhs to calculate
)
{
  // Initialize the element force vector. The matrix is not used.
  eforce.Scale(0.0);

  // space dimension for 3d fluid element
  const size_t nsd = 3;

  const int NUMDOF = 4;
  COMBUST::LocalAssembler<DISTYPE, ASSTYPE, NUMDOF> assembler(dofman, estif, eforce);

  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // extract the element phi vector from mystate
  LINALG::Matrix<numnode,1> ephi(true);
  for (size_t iparam=0; iparam<numnode; ++iparam)
    ephi(iparam) = mystate.phinp_[iparam];

  // get node coordinates of the current element
  static LINALG::Matrix<nsd,numnode> xyze;
  GEO::fillInitialPositionArray<DISTYPE>(ele, xyze);

  // get the enrichments for pressure
  const size_t numparampres = XFEM::NumParam<numnode,ASSTYPE>::get(dofman, XFEM::PHYSICS::Pres);

  // information about domain integration cells
  const GEO::DomainIntCells&  domainIntCells(ih->ElementDomainIntCells(ele->Id()));

  // loop over integration cells
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {
    // special getXFEMGaussruleKinkEnr for kink enrichment is called as parabolic shape functions are obtained
    // after multipying N and Psi
    const DRT::UTILS::GaussRule3D gaussrule = XFEM::getXFEMGaussruleKinkEnr<DISTYPE>(ele, xyze, ele->Splited(),cell->Shape());

    // gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

    // integration loop over Gauss points
    for(int iquad=0; iquad<intpoints.nquad; ++iquad)
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

      // create dummy derivates since ElementEnrichmentValues needs them
      static LINALG::Matrix<3,numnode> dummyderxy;
      dummyderxy.Clear();

      static LINALG::Matrix<6,numnode> dummyderxy2;
      dummyderxy2.Clear();

      // the enrichment functions depend on the gauss point
      // therefore the comutation of the enrichment functions is called here
      // the gauss point is contained in funct!

      // kink enrichments are called with level-set values
      // jump enrichments are called with domain cells!
      const XFEM::ElementEnrichmentValues enrvals(*ele, dofman, ephi, *cell,
                                                  funct, dummyderxy, dummyderxy2);

      const size_t shpVecSize = COMBUST::SizeFac<ASSTYPE>::fac*DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      static XFEM::ApproxFunc<2,shpVecSize> shppres;

      static LINALG::Matrix<shpVecSize,1> enr_funct_pres;

      if (ASSTYPE == XFEM::xfem_assembly)
      {
        // shape function for nodal dofs pressure
        enrvals.ComputeModifiedEnrichedNodalShapefunction(
            XFEM::PHYSICS::Pres,
            funct,
            enr_funct_pres);

        for (size_t iparam = 0; iparam < numparampres; ++iparam)
          shppres.d0(iparam) = enr_funct_pres(iparam);
      }
      else // not xfem_assembly
      {
        for (size_t iparam = 0; iparam < numnode; ++iparam)
          shppres.d0(iparam) = funct(iparam);
      }

      // do the assembling
      assembler.template Vector<XFEM::PHYSICS::Pres>(shppres.d0, fac);

    } // end loop over gauss points
  } // end loop over integration cells

}
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void COMBUST::callIntegrateShape(
    const XFEM::AssemblyType             assembly_type,
    const DRT::ELEMENTS::Combust3*       ele,
    const COMBUST::InterfaceHandleCombust* ih,
    const XFEM::ElementDofManager&       eleDofManager,
    const DRT::ELEMENTS::Combust3::MyState& mystate,      ///< element state variables
    Epetra_SerialDenseMatrix&            estif,
    Epetra_SerialDenseVector&            eforce)
{
  if (assembly_type == XFEM::standard_assembly)
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::IntegrateShape<DRT::Element::hex8,XFEM::standard_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce);
    break;
    case DRT::Element::hex20:
      COMBUST::IntegrateShape<DRT::Element::hex20,XFEM::standard_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce);
    break;
    default:
      dserror("standard_assembly IntegrateShape not templated yet");
    };
  }
  else
  {
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
      COMBUST::IntegrateShape<DRT::Element::hex8,XFEM::xfem_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce);
    break;
    case DRT::Element::hex20:
      COMBUST::IntegrateShape<DRT::Element::hex20,XFEM::xfem_assembly>(
          ele, ih, eleDofManager, mystate, estif, eforce);
    break;
    default:
      dserror("xfem_assembly IntegrateShape not templated yet");
    };
  }
}


namespace COMBUST
{
template <DRT::Element::DiscretizationType DISTYPE,XFEM::AssemblyType M_ASSTYPE,XFEM::AssemblyType S_ASSTYPE>
void Facemat(
        const DRT::ELEMENTS::Combust3IntFace*      ele,             ///< internal face element
        const COMBUST::InterfaceHandleCombust*     ih,              ///< interface handler
        std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>&  estif_blocks,    ///< element matrix blocks
        std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>&  eforce_blocks,   ///< element vector blocks
        const INPAR::COMBUST::CombustionType       combusttype,     ///< problem type (combustion (jumps) or two-phase (kinks)
        const std::vector<double>&                 m_velnp,         ///< velocity/pressure unknowns at n+1 (master)
        const std::vector<double>&                 m_velaf,         ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (master)
        const std::vector<double>&                 m_phinp,         ///< phi at n+1 (master)
        const std::vector<double>&                 s_velnp,         ///< velocity/pressure unknowns at n+1 (slave)
        const std::vector<double>&                 s_velaf,         ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (slave)
        const std::vector<double>&                 s_phinp,         ///< phi at n+1 (slave)
        Teuchos::RCP<const MAT::Material>          material,        ///< material list
        const INPAR::FLUID::TimeIntegrationScheme  timealgo,        ///< time discretization type
        const double                               time,            ///< current time step
        const double                               dt,              ///< delta t (time step size)
        const double                               theta,           ///< factor for one step theta scheme
        const double                               ga_alphaF,       ///< factor for gen-alpha
        const double                               ga_alphaM,       ///< factor for gen-alpha
        const double                               ga_gamma,        ///< factor for gen-alpha
        const INPAR::FLUID::EOS_Pres               pres_stab,       ///< face term to evaluate: pressure stab
        const INPAR::FLUID::EOS_Conv_Stream        conv_stream_stab,///< face term to evaluate: conv stab stream
        const INPAR::FLUID::EOS_Conv_Cross         conv_cross_stab, ///< face term to evaluate: conv stab cross
        const INPAR::FLUID::EOS_Div                conti_stab,      ///< face term to evaluate: conti stab
        const INPAR::FLUID::EOS_ElementLength      hk_def,          ///< definition of characteristic element length
        const INPAR::FLUID::EOS_TauType            tau_def,         ///< definition of stab parameter
        const INPAR::FLUID::EOS_GP_Pattern         pattern,         ///< pattern for edge-based stabilization
        const bool                                 xfem_stab,       ///< flag to indicate wheter ghost penalties should be used
        const INPAR::XFEM::FaceType                facetype,        ///< flag to define faces which have to be integrated twice
        const std::map<XFEM::PHYSICS::Field,std::vector<int> >&                    lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
        const std::map<XFEM::PHYSICS::Field,std::vector<int> >&                    lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
        const XFEM::ElementDofManager&             m_eleDofManager, ///< element dof manager of master element
        const XFEM::ElementDofManager&             s_eleDofManager  ///< element dof manager of slave element
        )
{
  // initialize element stiffness matrix and force vector
    std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>::iterator it_stif;
    for (it_stif=estif_blocks.begin(); it_stif!=estif_blocks.end(); it_stif++)
      it_stif->second.Scale(0.0);
    std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>::iterator it_force;
    for (it_force=eforce_blocks.begin(); it_force!=eforce_blocks.end(); it_force++)
      it_force->second.Scale(0.0);

  const int NUMDOF = 4;

  // instationary formulation
  bool instationary = true;
  if (timealgo == INPAR::FLUID::timeint_stationary) instationary = false;
  // generalized alpha time integration scheme
  bool genalpha = false;
  if (timealgo == INPAR::FLUID::timeint_afgenalpha) genalpha = true;

  const unsigned int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  const unsigned int m_fac = COMBUST::SizeFac<M_ASSTYPE>::fac;
  const unsigned int s_fac = COMBUST::SizeFac<S_ASSTYPE>::fac;

  const unsigned int m_shpVecSize = m_fac * numnode;
  const unsigned int s_shpVecSize = s_fac * numnode;

  // here, we have to distinguish kink and jump enrichments
  switch(combusttype)
  {
    // jump cases
    case INPAR::COMBUST::combusttype_premixedcombustion:
    case INPAR::COMBUST::combusttype_twophaseflowjump:
    {
      switch (ele->ParentMasterElement()->Shape())
      {
        case DRT::Element::hex8:
        {
          if (ele->ParentSlaveElement()->Shape() == DRT::Element::hex8)
          {

            LINALG::Matrix<3,m_shpVecSize> evelaf_m(true);
            LINALG::Matrix<m_shpVecSize,1> epreaf_m(true);
            LINALG::Matrix<numnode,1> ephi_m(true);
            LINALG::Matrix<3,s_shpVecSize> evelaf_s(true);
            LINALG::Matrix<s_shpVecSize,1> epreaf_s(true);
            LINALG::Matrix<numnode,1> ephi_s(true);
            // master element
            COMBUST::fillElementUnknownsArrays<DRT::Element::hex8,M_ASSTYPE>(
                    m_eleDofManager,
                    m_velnp, m_velaf, m_phinp,
                    instationary, genalpha,
                    evelaf_m, epreaf_m, ephi_m);
            // slave element
            COMBUST::fillElementUnknownsArrays<DRT::Element::hex8,S_ASSTYPE>(
                    s_eleDofManager,
                    s_velnp, s_velaf, s_phinp,
                    instationary, genalpha,
                    evelaf_s, epreaf_s, ephi_s);
            COMBUST::LocalPatchAssembler<DRT::Element::hex8,DRT::Element::hex8,M_ASSTYPE,S_ASSTYPE,NUMDOF> assembler(
                    m_eleDofManager, s_eleDofManager,
                    estif_blocks, eforce_blocks, pattern,
                    lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch
                    );

            // evaluate face-based terms
            COMBUST::FaceMatNitsche<DISTYPE,DRT::Element::hex8,DRT::Element::hex8,M_ASSTYPE,S_ASSTYPE,NUMDOF>(
              ele, ih, assembler,
              evelaf_m, epreaf_m, ephi_m,
              evelaf_s, epreaf_s, ephi_s,
              material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, xfem_stab, facetype,
              pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def
              );
          }
          else
            dserror("Unknown element type");
          break;
        }
        default:
        {
          dserror("Unknown element type");
          break;
        }
      }
      break;
    }
    // kink cases
    case INPAR::COMBUST::combusttype_twophaseflow:
    case INPAR::COMBUST::combusttype_twophaseflow_surf:
    {
      switch (ele->ParentMasterElement()->Shape())
      {
        case DRT::Element::hex8:
        {
          if (ele->ParentSlaveElement()->Shape() == DRT::Element::hex8)
          {
            LINALG::Matrix<3,m_shpVecSize> evelaf_m(true);
            LINALG::Matrix<m_shpVecSize,1> epreaf_m(true);
            LINALG::Matrix<numnode,1> ephi_m(true);
            LINALG::Matrix<3,s_shpVecSize> evelaf_s(true);
            LINALG::Matrix<s_shpVecSize,1> epreaf_s(true);
            LINALG::Matrix<numnode,1> ephi_s(true);
            // master element
            COMBUST::fillElementUnknownsArrays<DRT::Element::hex8,M_ASSTYPE>(
                    m_eleDofManager,
                    m_velnp, m_velaf, m_phinp,
                    instationary, genalpha,
                    evelaf_m, epreaf_m, ephi_m);
            // slave element
            COMBUST::fillElementUnknownsArrays<DRT::Element::hex8,S_ASSTYPE>(
                    s_eleDofManager,
                    s_velnp, s_velaf, s_phinp,
                    instationary, genalpha,
                    evelaf_s, epreaf_s, ephi_s);
            COMBUST::LocalPatchAssembler<DRT::Element::hex8,DRT::Element::hex8,M_ASSTYPE,S_ASSTYPE,NUMDOF> assembler(
                    m_eleDofManager, s_eleDofManager,
                    estif_blocks, eforce_blocks, pattern,
                    lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch
                    );

            // evaluate face-based terms
            COMBUST::FaceMatTwoPhase<DISTYPE,DRT::Element::hex8,DRT::Element::hex8,M_ASSTYPE,S_ASSTYPE,NUMDOF>(
              ele, ih, assembler,
              evelaf_m, epreaf_m, ephi_m,
              evelaf_s, epreaf_s, ephi_s,
              material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma, xfem_stab, facetype,
              pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def
              );
          }
          else
            dserror("Unknown element type");
          break;
        }
        default:
        {
          dserror("Unknown element type");
          break;
        }
      }
      break;
    }
    default:
    {
      dserror("Unknown combust type!");
      break;
    }
  }

  return;
}
}

void COMBUST::callFacemat(
     const DRT::ELEMENTS::Combust3IntFace*      ele,             ///< internal face element
     const COMBUST::InterfaceHandleCombust*     ih,              ///< interface handler
     std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>&  estif_blocks,    ///< element matrix blocks
     std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>&  eforce_blocks,   ///< element vector blocks
     const INPAR::COMBUST::CombustionType       combusttype,     ///< problem type (combustion (jumps) or two-phase (kinks)
     const std::vector<double>&                 m_velnp,         ///< velocity/pressure unknowns at n+1 (master)
     const std::vector<double>&                 m_velaf,         ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (master)
     const std::vector<double>&                 m_phinp,         ///< phi at n+1 (master)
     const std::vector<double>&                 s_velnp,         ///< velocity/pressure unknowns at n+1 (slave)
     const std::vector<double>&                 s_velaf,         ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (slave)
     const std::vector<double>&                 s_phinp,         ///< phi at n+1 (slave)
     Teuchos::RCP<const MAT::Material>          material,        ///< material list
     const INPAR::FLUID::TimeIntegrationScheme  timealgo,        ///< time discretization type
     const double                               time,            ///< current time step
     const double                               dt,              ///< delta t (time step size)
     const double                               theta,           ///< factor for one step theta scheme
     const double                               ga_alphaF,       ///< factor for gen-alpha
     const double                               ga_alphaM,       ///< factor for gen-alpha
     const double                               ga_gamma,        ///< factor for gen-alpha
     const INPAR::FLUID::EOS_Pres               pres_stab,       ///< face term to evaluate: pressure stab
     const INPAR::FLUID::EOS_Conv_Stream        conv_stream_stab,///< face term to evaluate: conv stab stream
     const INPAR::FLUID::EOS_Conv_Cross         conv_cross_stab, ///< face term to evaluate: conv stab cross
     const INPAR::FLUID::EOS_Div                conti_stab,      ///< face term to evaluate: conti stab
     const INPAR::FLUID::EOS_ElementLength      hk_def,          ///< definition of characteristic element length
     const INPAR::FLUID::EOS_TauType            tau_def,         ///< definition of stab parameter
     const INPAR::FLUID::EOS_GP_Pattern         pattern,         ///< pattern for edge-based stabilization
     const bool                                 xfem_stab,       ///< flag to indicate wheter ghost penalties should be used
     const INPAR::XFEM::FaceType                facetype,        ///< flag to define faces which have to be integrated twice
     const XFEM::AssemblyType                   m_assembly_type, ///< assembly type of master element
     const XFEM::AssemblyType                   s_assembly_type, ///< assembly type of slave element
     const std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
     const std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
     const XFEM::ElementDofManager&             m_eleDofManager, ///< element dof manager of master element
     const XFEM::ElementDofManager&             s_eleDofManager  ///< element dof manager of slave element
     )
{
  switch (ele->Shape())
  {
    case DRT::Element::quad4:
    {
      switch (m_assembly_type)
      {
        case XFEM::standard_assembly:
        {
          if (s_assembly_type == XFEM::standard_assembly)
            COMBUST::Facemat<DRT::Element::quad4,XFEM::standard_assembly,XFEM::standard_assembly>(
                        ele, ih, estif_blocks, eforce_blocks, combusttype,
                        m_velnp, m_velaf, m_phinp,
                        s_velnp, s_velaf, s_phinp,
                        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma,
                        pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def, pattern, xfem_stab, facetype,
                        lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
                        m_eleDofManager, s_eleDofManager);
          else
            COMBUST::Facemat<DRT::Element::quad4,XFEM::standard_assembly,XFEM::xfem_assembly>(
                        ele, ih, estif_blocks, eforce_blocks, combusttype,
                        m_velnp, m_velaf, m_phinp,
                        s_velnp, s_velaf, s_phinp,
                        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma,
                        pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def, pattern, xfem_stab, facetype,
                        lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
                        m_eleDofManager, s_eleDofManager);
          break;
        }
        case XFEM::xfem_assembly:
        {
          if (s_assembly_type == XFEM::standard_assembly)
            COMBUST::Facemat<DRT::Element::quad4,XFEM::xfem_assembly,XFEM::standard_assembly>(
                        ele, ih, estif_blocks, eforce_blocks, combusttype,
                        m_velnp, m_velaf, m_phinp,
                        s_velnp, s_velaf, s_phinp,
                        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma,
                        pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def, pattern, xfem_stab, facetype,
                        lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
                        m_eleDofManager, s_eleDofManager);
          else
            COMBUST::Facemat<DRT::Element::quad4,XFEM::xfem_assembly,XFEM::xfem_assembly>(
                        ele, ih, estif_blocks, eforce_blocks, combusttype,
                        m_velnp, m_velaf, m_phinp,
                        s_velnp, s_velaf, s_phinp,
                        material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma,
                        pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def, pattern, xfem_stab, facetype,
                        lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
                        m_eleDofManager, s_eleDofManager);
          break;
        }
        default:
        {
          dserror("Unknown assembly type");
          break;
        }
      }
    }
    break;
    case DRT::Element::quad8:
    {
      dserror("Hex20 not yet considered for face integrals!");
    }
    break;
    default:
    {
      dserror("Unknown element type for face integrals!");
      break;
    }
  };

  return;
}

