/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_multiporo_reac.cpp

 \brief evaluation class containing routines for calculation of scalar transport
        within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_multiporo_reac.H"

#include "scatra_ele_parameter_timint.H"

#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/scatra_mat_multiporo.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/matlist_reactions.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::ScaTraEleCalcMultiPoroReac(const int numdofpernode,const int numscal,const std::string& disname)
: DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(numdofpernode,numscal,disname),
  efluxnp_(0),
  epresnp_(0),
  esatnp_(0),
  esolidpresnp_(true)
{
  // replace internal variable manager by internal variable manager for muliporo
  my::scatravarmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerMultiPoro<my::nsd_, my::nen_>(my::numscal_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype> * DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcMultiPoroReac* delete_me )
{
  static std::map<std::string,ScaTraEleCalcMultiPoroReac<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcMultiPoroReac<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcMultiPoroReac<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}

/*----------------------------------------------------------------------*
 | setup element evaluation                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::SetupCalc(
    DRT::Element*               ele,
    DRT::Discretization&        discretization
    )
{
  pororeac::SetupCalc(ele,discretization);

  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  if (material->MaterialType() == INPAR::MAT::m_matlist or
      material->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < my::numdofpernode_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numdofpernode_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      switch(singlemat->MaterialType())
      {
        case INPAR::MAT::m_scatra_multiporo:
        {
          const Teuchos::RCP<const MAT::ScatraMatMultiPoro>& poromat
            = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoro>(singlemat);

          VarManager()->SetPhaseID(k,poromat->PhaseID());
          break;
        }

        default:
        {
          dserror("Material type %i is not supported for multiphase flow through porous media!",singlemat->MaterialType());
          break;
        }
      }
    }
  }
  else
  {
    switch(material->MaterialType())
    {
      case INPAR::MAT::m_scatra_multiporo:
      {
        const Teuchos::RCP<const MAT::ScatraMatMultiPoro>& poromat
          = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoro>(material);

        VarManager()->SetPhaseID(0,poromat->PhaseID());
        break;
      }

      default:
      {
        dserror("Material type %i is not supported for multiphase flow through porous media!",material->MaterialType());
        break;
      }
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
)
{
  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

  // get number of dofset associated with pressure related dofs
  const int ndspres = params.get<int>("ndspres");

  // determine number of velocity related dofs per node (= number of phases)
  const int numphases = la[ndspres].lm_.size()/my::nen_;

  //resize state vectors based on number of phases
  efluxnp_.resize(numphases);
  epresnp_.resize(numphases);
  esatnp_.resize(numphases);

  std::string stateprefix = "flux";
  for(int curphase=0;curphase<numphases;curphase++)
  {
    std::stringstream statename;
    statename << stateprefix << curphase;

    // get convective (velocity - mesh displacement) velocity at nodes
    Teuchos::RCP<const Epetra_Vector> convel = discretization.GetState(ndsvel, statename.str());
    if(convel == Teuchos::null)
      dserror("Cannot get state vector %s",statename.str().c_str());

    // extract local values of convective velocity field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nsd_,my::nen_> >(*convel,efluxnp_[curphase],la[ndsvel].lm_);
  }

  // get additional state vector for ALE case: grid displacement
  if (my::scatrapara_->IsAle())
  {
    // get number of dofset associated with displacement related dofs
    const int ndsdisp = params.get<int>("ndsdisp");

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp==Teuchos::null)
      dserror("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp].lm_.size()/my::nen_;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(my::nsd_*my::nen_,-1);
    for (int inode=0; inode<my::nen_; ++inode)
      for (int idim=0; idim<my::nsd_; ++idim)
        lmdisp[inode*my::nsd_+idim] = la[ndsdisp].lm_[inode*numdispdofpernode+idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nsd_,my::nen_> >(*dispnp,my::edispnp_,lmdisp);

    // add nodal displacements to point coordinates
    my::UpdateNodeCoordinates();
  }
  else
  {
    my::edispnp_.Clear();
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist==Teuchos::null || phinp==Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'");

  //values of scatra field are always in first dofset
  const std::vector<int>&    lm = la[0].lm_;
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*hist,my::ehist_,lm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phinp,my::ephinp_,lm);

  if (my::scatraparatimint_->IsGenAlpha() and not my::scatraparatimint_->IsIncremental())
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    if (phin==Teuchos::null) dserror("Cannot get state vector 'phin'");
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phin,my::ephin_,lm);
  }

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> presnp = discretization.GetState(ndspres,"pressure");
  if (presnp==Teuchos::null)
    dserror("Cannot get state vector 'pressure'");

  // extract local values of pressure field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*presnp,epresnp_,la[ndspres].lm_);

  // get number of dofset associated with saturation related dofs
  const int ndssat = params.get<int>("ndssat");
  if((int)la[ndssat].lm_.size()!=numphases*my::nen_)
  {
    dserror("Number of DOFs of saturation vector unequal to number of phases given by the pressure vector!");
  }
  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> satnp = discretization.GetState(ndssat,"saturation");
  if (satnp==Teuchos::null)
    dserror("Cannot get state vector 'saturation'");
  // extract local values of saturation field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*satnp,esatnp_,la[ndssat].lm_);


  // get number of dofset associated with solid pressure related dofs
  const int nds_solid_pressure = params.get<int>("ndssolidpressure");
  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> solidprenp = discretization.GetState(nds_solid_pressure,"solid_pressure");
  if (solidprenp==Teuchos::null)
    dserror("Cannot get state vector 'solid_pressure'");
  // extract local values of solid pressure field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*solidprenp,esolidpresnp_,la[nds_solid_pressure].lm_);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::BodyForce(ele);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  my::OtherNodeBasedSourceTerms(lm,discretization,params);

  return;
}

/*----------------------------------------------------------------------*
 |  compute the solid pressure at gauss point  (protected)    vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::ComputePorePressure(
  )
{
  return VarManager()->SolidPressure();
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point

  )
{
  switch(material->MaterialType())
  {
    case INPAR::MAT::m_scatra_multiporo:
    {
      MatMultiPoro(material,k,densn,densnp,densam,visc,iquad);
      break;
    }

    default:
    {
      dserror("Material type %i is not supported for multiphase flow through porous media!",material->MaterialType());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::MatMultiPoro(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  if(iquad==-1)
    dserror("no gauss point given for evaluation of MatMultiPoro material. Check your input file.");

  const Teuchos::RCP<const MAT::ScatraMatMultiPoro>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoro>(material);

  //read the porosity from the diffusion manager and scale it with the saturation and the density
  const double porosity = poro::DiffManager()->GetPorosity(k)*VarManager()->Saturation(k)*actmat->Density();

  {
    // set diffusivity (scaled with porosity)
    poro::SetDiffusivity(actmat,k,porosity);

    // set densities (scaled with porosity)
    poro::SetDensities(porosity,densn,densnp,densam);
  }

  return;
} // ScaTraEleCalcMultiPoroReac<distype>::MatMultiPoro

/*------------------------------------------------------------------------------*
 | set internal variables                                           vuong 08/16 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::SetInternalVariablesForMatAndRHS()
{
  VarManager()->SetInternalVariablesMultiPoro(
      my::funct_,
      my::derxy_,
      my::ephinp_,
      my::ephin_,
      efluxnp_,
      epresnp_,
      esatnp_,
      esolidpresnp_,
      my::ehist_);

  return;
}

/*-------------------------------------------------------------------------------*
 |  Set advanced reaction terms and derivatives                      vuong 08/16 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::SetAdvancedReactionTerms(
    const int                                 k,           //!< index of current scalar
    const Teuchos::RCP<MAT::MatListReactions> matreaclist, //!< index of current scalar
    const double* gpcoord                                  //!< current Gauss-point coordinates
    )
{
  FillCouplingVector();

  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();

  remanager->AddToReaBodyForce(
      matreaclist->CalcReaBodyForceTerm(k,my::scatravarmanager_->Phinp(),couplingvalues_,gpcoord),
      k);

  matreaclist->CalcReaBodyForceDerivMatrix(
      k,
      remanager->GetReaBodyForceDerivVector(k),
      my::scatravarmanager_->Phinp(),
      couplingvalues_,
      gpcoord);

}

/*-------------------------------------------------------------------------------*
 |  fill the coupling vector                                         vuong 08/16 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::FillCouplingVector()
{
  // if it is empty rebuilt it
  if(couplingvalues_.empty())
  {
    //pressures
    const std::vector<double>& pressures = VarManager()->Pressure();
    const int numphases = pressures.size();
    for(int i =0;i<numphases;i++)
    {
      std::ostringstream temp;
      temp << i+1;
      couplingvalues_.push_back(std::pair<std::string,double>("p"+temp.str(),pressures[i]));
    }
    //saturation
    const std::vector<double>& saturations = VarManager()->Saturation();
    for(int i =0;i<numphases;i++)
    {
      std::ostringstream temp;
      temp << i+1;

      couplingvalues_.push_back(std::pair<std::string,double>("S"+temp.str(),saturations[i]));
    }
    //porosity
    couplingvalues_.push_back(std::pair<std::string,double>("porosity",poro::DiffManager()->GetPorosity(0)));
  }
  // directly copy values (rely on order for performance reasons)
  else
  {
    //pressures
    const std::vector<double>& pressures = VarManager()->Pressure();
    const int numphases = pressures.size();
    for(int i =0;i<numphases;i++)
    {
     // std::cout<<"pressure "<<i<<": "<<pressures[i]<<std::endl;
      couplingvalues_[i].second=pressures[i];
    }
    //saturation
    const std::vector<double>& saturations = VarManager()->Saturation();
    for(int i =0;i<numphases;i++)
    {
   //   std::cout<<"saturation "<<i<<": "<<saturations[i]<<std::endl;
      couplingvalues_[numphases+i].second=saturations[i];
    }
    //porosity
    couplingvalues_[2*numphases].second=poro::DiffManager()->GetPorosity(0);
    //std::cout<<"porosity: "<<poro::DiffManager()->GetPorosity(0)<<std::endl;
  }
}

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix in convective form    vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcMatConv(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  const double                  densnp,
  const LINALG::Matrix<my::nen_,1>& sgconv
  )
{
  //the only difference to the base class version is, that there is no scaling with the density
  pororeac::CalcMatConv(emat,k,timefacfac,1.0,sgconv);

  return;
} // ScaTraEleCalc<distype>::CalcMatConv


/*------------------------------------------------------------------------------------------*
 |  calculation of convective element matrix: add conservative contributions   vuong 08/16 |
 *------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcMatConvAddCons(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  const double                  vdiv,
  const double                  densnp
  )
{
  //the only difference to the base class version is, that there is no scaling with the density
  pororeac::CalcMatConvAddCons(emat,k,timefacfac,vdiv,1.0);

  return;
}

/*------------------------------------------------------------------- *
 | adaption of convective term for rhs                     vuong 08/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::RecomputeConvPhiForRhs(
  const int                     k,
  const LINALG::Matrix<my::nsd_,1>& sgvelint,
  const double                  densnp,
  const double                  densn,
  const double                  vdiv
  )
{
  //the only difference to the base class version is, that there is no scaling with the density
  pororeac::RecomputeConvPhiForRhs(k,sgvelint,1.0,1.0,vdiv);
  return;
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<DRT::Element::nurbs27>;
