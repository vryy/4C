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
          // set delta in the variablemanager
          VarManager()->SetDelta(poromat->Delta(), k);
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
        // set delta in the variablemanager
        VarManager()->SetDelta(poromat->Delta(), 0);
        break;
      }

      default:
      {
        dserror("Material type %i is not supported for multiphase flow through porous media!",material->MaterialType());
        break;
      }
    }
  }

  // set the fluid material in the element
  VarManager()->SetFluidPoromultiphaseMaterial(ele);

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

  // extract action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  VarManager()->SetAction(action);

  //---------------------------------------------------------------------------------------------
  //                                 STRUCTURE
  //---------------------------------------------------------------------------------------------

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
    for (unsigned inode=0; inode<my::nen_; ++inode)
      for (unsigned idim=0; idim<my::nsd_; ++idim)
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

  //---------------------------------------------------------------------------------------------
  //                                 SCATRA
  //---------------------------------------------------------------------------------------------

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

  //---------------------------------------------------------------------------------------------
  //                                 FLUID
  //---------------------------------------------------------------------------------------------

  // get number of dofset associated with pressure/fluid related dofs
  const int ndspres = params.get<int>("ndspres");

  // determine number of velocity related dofs per node (= number of phases)
  const int numphases = la[ndspres].lm_.size()/my::nen_;

  // this is a check if we have L2-based projection or evaluation at GP of fluid-quantities
  // if the fluid primary variable is present in the scatra discretization, we evaluate at gp
  // TODO: this is not very nice, is there a better way?
  //---------------------------------------------------------------------------------------------
  //                   CASE 1: no L2-projection --> fluid is handled by its own managers
  if(discretization.HasState(ndspres,"phinp_fluid"))
  {
    L2_projection_ = false;
    VarManager()->SetupPoroFluidManagers(ele,params,discretization,la,numphases);
    VarManager()->ExtractElementAndNodeValuesOfPoroFluid(ele,discretization,la,my::xyze_);
  }
  //---------------------------------------------------------------------------------------------
  //                   CASE 2: L2-projection --> fluid is handled by class variables of scatra
  else if(discretization.HasState(ndspres,"pressure"))
  {
    L2_projection_ = true;
    ExtractElementAndNodeValuesWithL2(ele, params, discretization, la, ndspres, numphases);
  }
  else
    dserror("Something went wrong here, scatra-dis has neither pressure nor fluid primary variables");

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
 | extract element based or nodal values (L2-projection)    vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::ExtractElementAndNodeValuesWithL2(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la,
    const int                     ndspres,
    const int                     numphases
    )
{
  //resize state vectors based on number of phases
  efluxnp_.resize(numphases);
  epresnp_.resize(numphases);
  esatnp_.resize(numphases);

  // get number of dofset associated with velocity related dofs
  const int ndsvel = params.get<int>("ndsvel");

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

  Teuchos::RCP<const Epetra_Vector> presnp = discretization.GetState(ndspres,"pressure");
  if (presnp==Teuchos::null)
    dserror("Cannot get state vector 'pressure'");

  // extract local values of pressure field from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*presnp,epresnp_,la[ndspres].lm_);

  // get number of dofset associated with saturation related dofs
  const int ndssat = params.get<int>("ndssat");
  if(la[ndssat].lm_.size()!= static_cast<unsigned>( numphases*my::nen_ ))
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
  double porosity = 0.0;
  // d_eff = d_0 * (porosity * saturation(k))^delta
  double d_eff = 0.0;
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    porosity = poro::DiffManager()->GetPorosity(k)*VarManager()->Saturation(k)*actmat->Density();
    d_eff = std::pow(poro::DiffManager()->GetPorosity(k)*VarManager()->Saturation(k),actmat->Delta());
  }

  {
    // set diffusivity (scaled with porosity)
    poro::SetDiffusivity(actmat,k,porosity*d_eff);

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

  if(!L2_projection_)
  {
  VarManager()->SetInternalVariablesMultiPoroWithoutL2(
      my::funct_,
      my::derxy_,
      my::deriv_,
      my::xjm_,
      pororeac::xyze0_,
      my::ephinp_,
      my::ephin_,
      my::ehist_);
  }
  else
  {
    VarManager()->SetInternalVariablesMultiPoroWithL2(
        my::funct_,
        my::derxy_,
        my::ephinp_,
        my::ephin_,
        efluxnp_,
        epresnp_,
        esatnp_,
        esolidpresnp_,
        my::ehist_);
  }

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
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();

  FillCouplingVectorAndAddVariables(k, matreaclist, remanager);

  const SCATRA::Action act = VarManager()->GetAction();

  remanager->AddToReaBodyForce(
      matreaclist->CalcReaBodyForceTerm(k,my::scatravarmanager_->Phinp(),couplingvalues_,gpcoord),
      k);

  std::vector<std::pair<std::string,double> > emptyconstants;

  if(act == SCATRA::Action::calc_mat_and_rhs ||
     act == SCATRA::Action::calc_initial_time_deriv)
  {
    matreaclist->CalcReaBodyForceDerivMatrix(
        k,
        remanager->GetReaBodyForceDerivVector(k),
        my::scatravarmanager_->Phinp(),
        couplingvalues_,
        gpcoord);
  }
  else if(act == SCATRA::Action::calc_scatra_mono_odblock_fluid)
  {
    matreaclist->CalcReaBodyForceDerivMatrixAddVariables(
        k,
        remanager->GetReaBodyForceDerivVectorAddVariables(k),
        my::scatravarmanager_->Phinp(),
        couplingvalues_,
        emptyconstants,
        gpcoord);
  }
  else if(act == SCATRA::Action::calc_scatra_mono_odblock_mesh)
  {
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnStruct())
    {
      matreaclist->CalcReaBodyForceDerivMatrixAddVariables(
          k,
          remanager->GetReaBodyForceDerivVectorAddVariables(k),
          my::scatravarmanager_->Phinp(),
          couplingvalues_,
          emptyconstants,
          gpcoord);
    }
  }
  else
    dserror("Wrong action type in VarManager(), action type is %d", act);

}

/*-------------------------------------------------------------------------------*
 |  fill the coupling vector and add variables to reactions          vuong 08/16 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::FillCouplingVectorAndAddVariables(
    const int k,
    const Teuchos::RCP<MAT::MatListReactions> matreaclist,
    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager)
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

    // initialize and add the variables to the reaction manager --> has to be done only once
    remanager->InitializeReaBodyForceDerivVectorAddVariables(my::numdofpernode_, couplingvalues_.size());
    for (int j = 0; j < my::numdofpernode_; j++)
      matreaclist->AddAdditionalVariables(j,couplingvalues_);
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
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    //the only difference to the base class version is, that there is no scaling with the density
    pororeac::CalcMatConv(emat,k,timefacfac,1.0,sgconv);
  }

  return;
} // ScaTraEleCalc<distype>::CalcMatConv

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix in convective form    vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcMatMass(
    Epetra_SerialDenseMatrix &     emat,
    const int &                    k,
    const double &                 fac,
    const double &                 densam
  )
{
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    //the only difference to the base class version is, that there is no scaling with the density
    pororeac::CalcMatMass(emat,k,fac,densam);
  }
  else
  {
    // If we have zero "densities" (porosity*saturation(k)), which mostly happens for tumor
    // cells, the whole equation will be equal to zero since it is scaled with the density
    // In that case also the mass fraction of the species (necrotic tumor cells) has to be zero
    // --> here we explicitly force it to be zero through a "Dirichlet" boundary condition
    for(unsigned vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+k;
      emat(fvi,fvi) += penalty_;
    }
  }
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

/*-------------------------------------------------------------------- *
 |  standard Galerkin convective term (OD mesh)       kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcConvODMesh(
  Epetra_SerialDenseMatrix&       emat,
  const int                       k,
  const int                       ndofpernodemesh,
  const double                    fac,
  const double                    rhsfac,
  const double                    densnp,
  const double                    J,
  const LINALG::Matrix<my::nsd_,1>&   gradphi,
  const LINALG::Matrix<my::nsd_,1>&   convelint
  )
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();

    const std::vector<LINALG::Matrix<my::nsd_,1> >& fluidgradphi = *(VarManager()->FluidVarManager()->GradPhinp());

    // current pressure gradient
    const LINALG::Matrix<my::nsd_,1> gradpres = VarManager()->PressureGradient(curphase);
    const double abspressgrad = VarManager()->AbsPressureGradient(curphase);

    // diffusion tensor
    LINALG::Matrix<my::nsd_,my::nsd_> difftensor(true);
    VarManager()->FluidPhaseManager()->PermeabilityTensor(curphase,difftensor);
    difftensor.Scale(VarManager()->FluidPhaseManager()->RelPermeability(curphase)/
        VarManager()->FluidPhaseManager()->DynViscosity(curphase, abspressgrad,2));

    // linearization of mesh motion
    //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
    // J denotes the determinant of the Jacobian of the mapping between current and parameter space, i.e. det(dx/ds)
    // in our case: rhsfac = J * dt * theta --> d(rhsfac)/dd = rhsfac * N_x
    for (unsigned vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+k;
      const double v =  rhsfac*my::funct_(vi)*(-1.0)*VarManager()->ConvPhi(k);

      for (unsigned ui=0; ui<my::nen_; ++ui)
      {
        for (unsigned idim=0; idim<my::nsd_; ++idim)
        {
          const int fui = ui*my::nsd_+idim;
          emat(fvi,fui) += v*my::derxy_(idim,ui);
        }
      }
    }

    //----------------------------------------------------------------
    // standard Galerkin terms  -- "shapederivatives" pressure gradient
    //----------------------------------------------------------------

    //gradient of pressure w.r.t. reference coordinates
    LINALG::Matrix<my::nsd_,1> refgradpres(true);
    refgradpres.Clear();

      //gradient of phi w.r.t. reference coordinates
    std::vector<LINALG::Matrix<my::nsd_,1> > reffluidgradphi(numphases,LINALG::Matrix<my::nsd_,1>(true));
    for (int idof=0; idof<numphases; ++idof)
      reffluidgradphi[idof].Multiply(my::xjm_,fluidgradphi[idof]);

      // compute the pressure gradient from the phi gradients
      for (int idof=0; idof<numphases; ++idof)
        refgradpres.Update(VarManager()->FluidPhaseManager()->PressureDeriv(curphase,idof),reffluidgradphi[idof],1.0);

      if(my::nsd_==3)
      {
        const double xjm_0_0   = my::xjm_(0, 0);
        const double xjm_0_1   = my::xjm_(0, 1);
        const double xjm_0_2   = my::xjm_(0, 2);
        const double xjm_1_0   = my::xjm_(1, 0);
        const double xjm_1_1   = my::xjm_(1, 1);
        const double xjm_1_2   = my::xjm_(1, 2);
        const double xjm_2_0   = my::xjm_(2, 0);
        const double xjm_2_1   = my::xjm_(2, 1);
        const double xjm_2_2   = my::xjm_(2, 2);

        {
          const double refgradpres_0   = refgradpres(0);
          const double refgradpres_1   = refgradpres(1);
          const double refgradpres_2   = refgradpres(2);

          const double gradphi_0   = gradphi(0);
          const double gradphi_1   = gradphi(1);
          const double gradphi_2   = gradphi(2);

          // TODO: anisotropic difftensor and
          //       non-constant viscosity (because of pressure gradient, probably not really necessary)
          const double vrhs = rhsfac*1.0/J*difftensor(0,0)*(-1.0);

          for (unsigned ui = 0; ui < my::nen_; ++ui)
          {
            const double v00 = + gradphi_1 * (
                                                refgradpres_0 * (my::deriv_(2, ui)*xjm_1_2 - my::deriv_(1, ui)*xjm_2_2)
                                              + refgradpres_1 * (my::deriv_(0, ui)*xjm_2_2 - my::deriv_(2, ui)*xjm_0_2)
                                              + refgradpres_2 * (my::deriv_(1, ui)*xjm_0_2 - my::deriv_(0, ui)*xjm_1_2)
                                                )
                               + gradphi_2 * (
                                                refgradpres_0 * (my::deriv_(1, ui)*xjm_2_1 - my::deriv_(2, ui)*xjm_1_1)
                                              + refgradpres_1 * (my::deriv_(2, ui)*xjm_0_1 - my::deriv_(0, ui)*xjm_2_1)
                                              + refgradpres_2 * (my::deriv_(0, ui)*xjm_1_1 - my::deriv_(1, ui)*xjm_0_1)
                                                );
            const double v01 = + gradphi_0 * (
                                                refgradpres_0 * (my::deriv_(1, ui)*xjm_2_2 - my::deriv_(2, ui)*xjm_1_2)
                                              + refgradpres_1 * (my::deriv_(2, ui)*xjm_0_2 - my::deriv_(0, ui)*xjm_2_2)
                                              + refgradpres_2 * (my::deriv_(0, ui)*xjm_1_2 - my::deriv_(1, ui)*xjm_0_2))
                               + gradphi_2 * (  refgradpres_0 * (my::deriv_(2, ui)*xjm_1_0 - my::deriv_(1, ui)*xjm_2_0)
                                              + refgradpres_1 * (my::deriv_(0, ui)*xjm_2_0 - my::deriv_(2, ui)*xjm_0_0)
                                              + refgradpres_2 * (my::deriv_(1, ui)*xjm_0_0 - my::deriv_(0, ui)*xjm_1_0)
                                                );
            const double v02 = + gradphi_0 * (
                                                refgradpres_0 * (my::deriv_(2, ui)*xjm_1_1 - my::deriv_(1, ui)*xjm_2_1)
                                              + refgradpres_1 * (my::deriv_(0, ui)*xjm_2_1 - my::deriv_(2, ui)*xjm_0_1)
                                              + refgradpres_2 * (my::deriv_(1, ui)*xjm_0_1 - my::deriv_(0, ui)*xjm_1_1)
                                                )
                               + gradphi_1 * (
                                                refgradpres_0 * (my::deriv_(1, ui)*xjm_2_0 - my::deriv_(2, ui)*xjm_1_0)
                                              + refgradpres_1 * (my::deriv_(2, ui)*xjm_0_0 - my::deriv_(0, ui)*xjm_2_0)
                                              + refgradpres_2 * (my::deriv_(0, ui)*xjm_1_0 - my::deriv_(1, ui)*xjm_0_0)
                                                );

            for (unsigned vi = 0; vi < my::nen_; ++vi)
            {
              const int fvi = vi*my::numdofpernode_+k;
              const double v = vrhs * my::funct_(vi);

              emat(fvi, ui * 3 + 0) += v * v00;
              emat(fvi, ui * 3 + 1) += v * v01;
              emat(fvi, ui * 3 + 2) += v * v02;
            }
          }
        }
      }
      else if(my::nsd_==2)
      {
        {
          const double refgradpres_0   = refgradpres(0);
          const double refgradpres_1   = refgradpres(1);

          const double gradphi_0   = gradphi(0);
          const double gradphi_1   = gradphi(1);

          // TODO: anisotropic difftensor and
          //       non-constant viscosity (because of pressure gradient, probably not really necessary)
          const double vrhs = rhsfac*1.0/J*difftensor(0,0)*(-1.0);

          for (unsigned ui = 0; ui < my::nen_; ++ui)
          {
            const double v00 = + gradphi_1 * (
                                                - refgradpres_0 * my::deriv_(1, ui)
                                                + refgradpres_1 * my::deriv_(0, ui)
                                                );
            const double v01 = + gradphi_0 * (
                                                  refgradpres_0 * my::deriv_(1, ui)
                                                - refgradpres_1 * my::deriv_(0, ui)
                                               )
                                                ;

            for (unsigned vi = 0; vi < my::nen_; ++vi)
            {
              const int fvi = vi*my::numdofpernode_+k;
              const double v = vrhs * my::funct_(vi);

              emat(fvi, ui * 2 + 0) += v * v00;
              emat(fvi, ui * 2 + 1) += v * v01;
            }
          }
        }
      }
      else
        dserror("shapederivatives not implemented for 1D!");

    //----------------------------------------------------------------
    // standard Galerkin terms  -- "shapederivatives" gradphi
    //----------------------------------------------------------------
    pororeac::ApplyShapeDerivsConv(emat,k,rhsfac,1.0,J,gradphi,convelint);
  }

  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin temporal term (OD mesh)         kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcLinMassODMesh(
  Epetra_SerialDenseMatrix&                  emat,
  const int                                  k,
  const int                                  ndofpernodemesh,
  const double                               rhsfac,
  const double                               fac,
  const double                               densam,
  const double                               densnp,
  const double                               phinp,
  const double                               hist,
  const double                               J,
  const LINALG::Matrix<1,my::nsd_*my::nen_>& dJ_dmesh
)
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    // call base class
    my::CalcLinMassODMesh(emat,
        k,
        ndofpernodemesh,
        rhsfac,
        fac,
        densam,
        densnp,
        phinp,
        hist,
        J,
        dJ_dmesh);

    double vtrans = 0.0;

    if (my::scatraparatimint_->IsGenAlpha())
      vtrans = rhsfac*densam*hist;
    else
    {
      // compute scalar at integration point
      vtrans = fac*densnp*phinp;
    }

    // linearization of porosity only if porosity is 'structure'-dependent
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnStruct())
      PorosityLinearizationStruct(emat, k, vtrans);
  }

  return;
}

/*-------------------------------------------------------------------- *
 | hist and source term (OD mesh)                     kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcHistAndSourceODMesh(
  Epetra_SerialDenseMatrix&                  emat,
  const int                                  k,
  const int                                  ndofpernodemesh,
  const double                               fac,
  const double                               rhsint,
  const double                               J,
  const LINALG::Matrix<1,my::nsd_*my::nen_>& dJ_dmesh,
  const double                               densnp
)
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    //const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();

    // call base class
    my::CalcHistAndSourceODMesh(emat,k,ndofpernodemesh,fac,rhsint,J,dJ_dmesh,densnp);

    double vrhs = -1.0*fac*rhsint;

    // linearization of porosity only if porosity is 'structure'-dependent
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnStruct())
      PorosityLinearizationStruct(emat, k, vrhs);

    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();
    if(remanager->Active() && VarManager()->FluidPhaseManager()->PorosityDependsOnStruct())
    {
      const std::vector<double> myderivs = remanager->GetReaBodyForceDerivVectorAddVariables(k);

      // porosity deriv at [2* numphases]: d reac / d d = d reac / d poro * d poro / d d
      // with
      // dporo/dd = dporo/dJ * dJ/dd = dporosity/dJ * J * N_x
      // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1

      const double poroderiv = myderivs[2*numphases]
                    *VarManager()->FluidPhaseManager()->JacobianDefGrad()
                    *VarManager()->FluidPhaseManager()->PorosityDerivWrtJacobianDefGrad();

      for (unsigned vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = vi*my::numdofpernode_+k;
        //TODO: gen-alpha
        const double v = my::funct_(vi)*poroderiv*my::scatraparatimint_->TimeFac()*fac*(-1.0)*densnp;

        for (unsigned ui=0; ui<my::nen_; ++ui)
        {
          for (unsigned idim=0; idim<my::nsd_; ++idim)
          {
            const int fui = ui*my::nsd_+idim;

            emat(fvi,fui) += v*my::derxy_(idim,ui);
          }
        }
      }
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 | diffusive term (OD mesh)                           kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcDiffODMesh(
      Epetra_SerialDenseMatrix&                   emat,
      const int                                   k,
      const int                                   ndofpernodemesh,
      const double                                fac,
      const double                                rhsfac,
      const double                                J,
      const LINALG::Matrix<my::nsd_,1>&           gradphi,
      const LINALG::Matrix<my::nsd_,1>&           convelint,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&  dJ_dmesh
  )
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    // call base class
    my::CalcDiffODMesh(emat,k,ndofpernodemesh,fac,rhsfac,J,gradphi,convelint,dJ_dmesh);

    // linearization of porosity only if it depends on deformation
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnStruct())
    {
      const double delta = VarManager()->GetDelta(k);

      // linearization of porosity
      //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
      // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
      // in our case: diffusivity is scaled with porosity^(delta + 1) --> scale it with 1.0/porosity^(delta + 1) here
      //              and build derivative d diff/d porosity = (delta + 1) * porosity^delta
      const double vrhs = rhsfac*my::diffmanager_->GetIsotropicDiff(k)
                         /std::pow(VarManager()->FluidPhaseManager()->Porosity(), delta + 1.0)
                         * (delta + 1.0) * std::pow(VarManager()->FluidPhaseManager()->Porosity(), delta )
                         *VarManager()->FluidPhaseManager()->JacobianDefGrad()
                         *VarManager()->FluidPhaseManager()->PorosityDerivWrtJacobianDefGrad();

      for (unsigned vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = vi*my::numdofpernode_+k;

        double laplawf(0.0);
        my::GetLaplacianWeakFormRHS(laplawf,gradphi,vi);
        const double v = vrhs * laplawf;
        for (unsigned ui=0; ui<my::nen_; ++ui)
        {
          for (unsigned idim=0; idim<my::nsd_; ++idim)
          {
            const int fui = ui*my::nsd_+idim;

            emat(fvi,fui) += v*my::derxy_(idim,ui);
          }
        }
      }
    }
  }
  return;
}

/*-------------------------------------------------------------------- *
 | reactive term (OD mesh)                            kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcReactODMesh(
      Epetra_SerialDenseMatrix&                   emat,
      const int                                   k,
      const int                                   ndofpernodemesh,
      const double                                rhsfac,
      const double                                rea_phi,
      const double                                J,
      const LINALG::Matrix<1,my::nsd_*my::nen_>&  dJ_dmesh
  )
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    // call base class
    my::CalcReactODMesh(emat,k,ndofpernodemesh,rhsfac,rea_phi,J,dJ_dmesh);

    if (my::reamanager_->Active())
    {
      // standard Galerkin term
      double vrhs = rhsfac*rea_phi;

      // linearization of porosity only if porosity is 'structure'-dependent
      if(VarManager()->FluidPhaseManager()->PorosityDependsOnStruct())
        PorosityLinearizationStruct(emat, k, vrhs);
    }
  }

  return;
}

/*--------------------------------------------------------------------------------- *
 | generic linearization of term scaled with porosity (OD struct)  kremheller 07/17 |
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::PorosityLinearizationStruct(
    Epetra_SerialDenseMatrix&                   emat,
    const int                                   k,
    double                                      prefac
)
{
  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  // in our case: prefac is scaled with density, i.e. porosity --> scale it with 1.0/porosity here

  prefac *= 1.0/(VarManager()->FluidPhaseManager()->Porosity())
              *VarManager()->FluidPhaseManager()->JacobianDefGrad()
              *VarManager()->FluidPhaseManager()->PorosityDerivWrtJacobianDefGrad();

  for (unsigned vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;
    const double v = my::funct_(vi)*prefac;

    for (unsigned ui=0; ui<my::nen_; ++ui)
    {
      for (unsigned idim=0; idim<my::nsd_; ++idim)
      {
        const int fui = ui*my::nsd_+idim;

        emat(fvi,fui) += v*my::derxy_(idim,ui);
      }
    }
  }

  return;
}

/*--------------------------------------------------------------------------------- *
 | generic linearization of term scaled with porosity (OD fluid)   kremheller 07/17 |
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::PorosityLinearizationFluid(
    Epetra_SerialDenseMatrix&                   emat,
    const int                                   k,
    const int                                   curphase,
    const int                                   numphases,
    double                                      prefac
)
{

  // linearization of porosity wrt fluid primary variables
  // in our case: prefac is scaled with density, i.e. porosity --> scale it with 1.0/porosity here

  prefac /= (VarManager()->FluidPhaseManager()->Porosity());

  for (unsigned vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;
    const double v = my::funct_(vi)*prefac;

    for (unsigned ui=0; ui<my::nen_; ++ui)
    {
      const double vfunct = v * my::funct_(ui);
      for (int idof=0; idof<numphases; ++idof)
      {
        const int fui = ui*numphases+idof;

        emat(fvi,fui) += vfunct*VarManager()->FluidPhaseManager()->PorosityDeriv(idof);
      }
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 | convective term (OD fluid)                         kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcMatConvODFluid(
  Epetra_SerialDenseMatrix&     emat,             //!< element matrix to be filled
  const int                     k,                //!< index of current scalar
  const int                     ndofpernodefluid, //!< number of dofs per node of fluid element
  const double                  rhsfac,           //!< domain-integration factor times time-integration factor
  const double                  densnp,           //!< density at time_(n+1)
  const LINALG::Matrix<my::nsd_,1>& gradphi           //!< scalar gradient
)
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();

    // current pressure gradient
    const LINALG::Matrix<my::nsd_,1> gradpres = VarManager()->PressureGradient(curphase);
    const double abspressgrad = VarManager()->AbsPressureGradient(curphase);

    // diffusion tensor
    LINALG::Matrix<my::nsd_,my::nsd_> difftensor(true);
    VarManager()->FluidPhaseManager()->PermeabilityTensor(curphase,difftensor);
    difftensor.Scale(VarManager()->FluidPhaseManager()->RelPermeability(curphase)/
        VarManager()->FluidPhaseManager()->DynViscosity(curphase, abspressgrad,2));

    // gradphi^T * difftensor
    LINALG::Matrix<1,my::nsd_> gradphiTdifftensor(true);
    gradphiTdifftensor.MultiplyTN(gradphi, difftensor);

    for (unsigned vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+k;
      const double v =  rhsfac*my::funct_(vi)*(-1.0);

      for (unsigned ui=0; ui<my::nen_; ++ui)
      {
        double laplawf(0.0);
        for (unsigned j = 0; j<my::nsd_; j++)
          laplawf += v*my::derxy_(j, ui)*gradphiTdifftensor(0, j);
        for (int idof=0; idof<numphases; ++idof)
        {
          const int fui = ui*numphases+idof;
          emat(fvi,fui) += laplawf*VarManager()->FluidPhaseManager()->PressureDeriv(curphase, idof);
        }
      }
    }

    //----------------------------------------------------------------
    // linearization of relative permeability w.r.t. dof
    //----------------------------------------------------------------
    if(not VarManager()->FluidPhaseManager()->HasConstantRelPermeability(curphase))
    {
      // diffusion tensor
      LINALG::Matrix<my::nsd_,my::nsd_> difftensor(true);
      VarManager()->FluidPhaseManager()->PermeabilityTensor(curphase,difftensor);
      difftensor.Scale(VarManager()->FluidPhaseManager()->RelPermeabilityDeriv(curphase)/
          VarManager()->FluidPhaseManager()->DynViscosity(curphase, abspressgrad,2));

      LINALG::Matrix<1,my::nsd_> gradphiTdifftensor(true);
      gradphiTdifftensor.MultiplyTN(gradphi, difftensor);

      for (unsigned vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = vi*my::numdofpernode_+k;
        const double v =  rhsfac*my::funct_(vi)*(-1.0);


        for (unsigned ui=0; ui<my::nen_; ++ui)
        {
          double laplawf(0.0);
          for (unsigned j = 0; j<my::nsd_; j++)
            laplawf += v*gradpres(j)*gradphiTdifftensor(0, j);
          for (int idof=0; idof<numphases; ++idof)
          {
            const int fui = ui*numphases+idof;
            emat(fvi,fui) += laplawf*my::funct_(ui)*VarManager()->FluidPhaseManager()->SaturationDeriv(curphase, idof);
          }
        }

      }
    }
    //----------------------------------------------------------------
    // linearization of dynamic viscosity w.r.t. dof --> TODO
    // however, I believe this is not necessary: FD-check does not fail
    //----------------------------------------------------------------
  }

  return;
}

/*-------------------------------------------------------------------------- *
 | convective term -- conservative contributions (OD fluid) kremheller 07/17 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcMatConvAddConsODFluid(
    Epetra_SerialDenseMatrix&     emat,             //!< element matrix to be filled
    const int                     k,                //!< index of current scalar
    const int                     ndofpernodefluid, //!< number of dofs per node of fluid element
    const double                  timefacfac,       //!< domain-integration factor times time-integration factor
    const double                  densnp,           //!< density at time_(n+1)
    const double                  phinp             //!< scalar at time_(n+1)
)
{

  dserror("CalcMatConvAddConsODFluid not yet available for scatre ele calc multiporo");
}

/*-------------------------------------------------------------------- *
 | temporal term (OD fluid)                           kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcLinMassODFluid(
    Epetra_SerialDenseMatrix&          emat,              //!< element matrix to be filled
    const int                          k,                 //!< index of current scalar
    const int                          ndofpernodemesh,   //!< number of dofs per node of fluid element // only a dummy variable
    const double                       rhsfac,            //!< time-integration factor for rhs times domain-integration factor
    const double                       fac,               //!< domain-integration factor
    const double                       densam,            //!< density at time_(n+am)
    const double                       densnp,            //!< density at time_(n+1)
    const double                       phinp,             //!< scalar at time_(n+1)
    const double                       hist               //!< history of time integartion
  )
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();

    double vtrans = 0.0;

    if (my::scatraparatimint_->IsGenAlpha())
      dserror("not implemented");
    else
    {
      // compute scalar at integration point
      vtrans = fac*densnp*phinp;
    }

    // linearization of saturation
    SaturationLinearization(emat, k, curphase, numphases, vtrans);

    // linearization of porosity only if porosity is pressure-dependent
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnFluid())
      PorosityLinearizationFluid(emat, k, curphase, numphases, vtrans);
  }

  return;
}

/*--------------------------------------------------------------------------------- *
 | generic linearization of term scaled with saturation (OD fluid) kremheller 07/17 |
 *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::SaturationLinearization(
    Epetra_SerialDenseMatrix&                   emat,
    const int                                   k,
    const int                                   curphase,
    const int                                   numphases,
    double                                      prefac
)
{
  // prefac is scaled with density, i.e. saturation, scale it with 1.0/saturation here
  prefac /= VarManager()->FluidPhaseManager()->Saturation(curphase);

  // dprefac / dphi_j = dprefac / dS_curphase * dS_curphase * d_phi_j = 1.0 * dS_curphase * d_phi_j
  for (unsigned vi=0; vi<my::nen_; ++vi)
  {
    const double v = prefac * my::funct_(vi);
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui=0; ui<my::nen_; ++ui)
    {
      const double vfunct = v*my::funct_(ui);
      for (int idof = 0; idof < numphases; ++idof)
      {
        const int fui = ui * numphases + idof;

        emat(fvi,fui) += vfunct*VarManager()->FluidPhaseManager()->SaturationDeriv(curphase, idof);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------- *
 |  standard Galerkin transient, old part of rhs and source term (OD fluid)     |
 |  + advanced reaction terms if reactionmanager is active     kremheller 07/17 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcHistAndSourceODFluid(
    Epetra_SerialDenseMatrix&          emat,
    const int                          k,
    const int                          ndofpernodemesh,
    const double                       fac,
    const double                       rhsint,
    const double                       densnp
  )
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();

    double vrhs = - fac*rhsint;

    // linearization of saturation
    SaturationLinearization(emat, k, curphase, numphases, vrhs);

    // linearization of porosity only if porosity is pressure-dependent
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnFluid())
      PorosityLinearizationFluid(emat, k, curphase, numphases, vrhs);

    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();
    if(remanager->Active() )
    {
      const std::vector<double> myderivs = remanager->GetReaBodyForceDerivVectorAddVariables(k);

      // derivatives after primary variables of fluid
      std::vector<double> phiderivs(numphases, 0.0);

      for(int i = 0; i < numphases; i++)
      {
        // porosity deriv at 2*numphases: d reac / d phi_i = d reac / d poro * d poro / d phi_i
        phiderivs[i] += myderivs[2*numphases]*VarManager()->FluidPhaseManager()->PorosityDeriv(i);
        for (int j = 0; j < numphases; j++)
        {
          // pressure derivs at       [0..numphases]: d reac / d phi_i = d reac / d pres_j * d pres_j / d phi_i
          // saturation derivs at [numph..2*numph-1]: d reac / d phi_i = d reac /  d sat_j *  d sat_j / d phi_i
          phiderivs[i] += myderivs[j+numphases]*VarManager()->FluidPhaseManager()->SaturationDeriv(j, i)
                        + myderivs[j          ]*VarManager()->FluidPhaseManager()->PressureDeriv(j, i);
        }
      }

      // fill matrix
      for (unsigned vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = vi*my::numdofpernode_+k;
        //TODO: gen-alpha?
        const double v = my::funct_(vi)*my::scatraparatimint_->TimeFac()*fac*(-1.0)*densnp;

        for (unsigned ui=0; ui<my::nen_; ++ui)
        {
          const double vfunct = v * my::funct_(ui);

          for (int idof=0; idof<numphases; ++idof)
          {
            const int fui = ui*numphases+idof;

            emat(fvi,fui) += vfunct * phiderivs[idof];
          }
        }
      }
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 | reactive term (OD fluid)                           kremheller 07/17 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcReactODFluid(
    Epetra_SerialDenseMatrix&          emat,
    const int                          k,
    const int                          ndofpernodemesh,
    const double                       rhsfac,
    const double                       rea_phi
  )
{
  if (my::reamanager_->Active() && fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();

    double vrhs = rhsfac*rea_phi;

    // linearization of saturation
    SaturationLinearization(emat, k, curphase, numphases, vrhs);

    // linearization of porosity only if porosity is pressure-dependent
    if(VarManager()->FluidPhaseManager()->PorosityDependsOnFluid())
      PorosityLinearizationFluid(emat, k, curphase, numphases, vrhs);
  }

  return;
}

/*------------------------------------------------------------------ *
 |  standard Galerkin diffusive term (OD fluid)     kremheller 07/17 |
 *----------------------------------------------------------------   */
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcMultiPoroReac<distype>::CalcDiffODFluid(
    Epetra_SerialDenseMatrix&           emat,             //!< element current to be filled
    const int                           k,                //!< index of current scalar
    const int                           ndofpernodemesh,  //!< number of dofs per node of ale element
    const double                        rhsfac,           //!< time-integration factor for rhs times domain-integration factor
    const LINALG::Matrix<my::nsd_,1>&   gradphi           //!< scalar gradient at Gauss point
  )
{
  // case of zero saturation
  if(fabs(VarManager()->Saturation(k)) > minimum_saturation_)
  {
    const int curphase = VarManager()->GetPhaseID(k);
    const int numphases = VarManager()->FluidPhaseManager()->NumPhases();
    const double delta = VarManager()->GetDelta(k);

    // linearization of saturation*porosity*d_eff = (saturation * porosity)^(delta+1) w.r.t saturation
    //-------------------------------------------------------------------------------------------------------------
    // in our case: diffusivity is scaled with saturation^(delta + 1) --> scale it with 1.0/saturation^(delta + 1) here
    //              and build derivative d diff/d saturation = (delta + 1) * saturation^delta
    const double vrhs_sat = rhsfac*my::diffmanager_->GetIsotropicDiff(k)
        /std::pow(VarManager()->FluidPhaseManager()->Saturation(curphase), delta + 1.0)
    * (delta + 1.0) * std::pow(VarManager()->FluidPhaseManager()->Saturation(curphase), delta );

    // linearization of saturation*porosity*d_eff = (saturation * porosity)^(delta+1) w.r.t porosity
    //-------------------------------------------------------------------------------------------------------------
    // in our case: diffusivity is scaled with porosity^(delta + 1) --> scale it with 1.0/porosity^(delta + 1) here
    //              and build derivative d diff/d porosity = (delta + 1) * porosity^delta
    const double vrhs_poro = rhsfac*my::diffmanager_->GetIsotropicDiff(k)
        /std::pow(VarManager()->FluidPhaseManager()->Porosity(), delta + 1.0)
    * (delta + 1.0) * std::pow(VarManager()->FluidPhaseManager()->Porosity(), delta );

    for (unsigned vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+k;

      double laplawf(0.0);
      my::GetLaplacianWeakFormRHS(laplawf,gradphi,vi);
      const double vsat = vrhs_sat * laplawf;
      const double vporo = vrhs_poro * laplawf;

      for (unsigned ui=0; ui<my::nen_; ++ui)
      {
        const double vfunctsat = vsat * my::funct_(ui);
        const double vfunctporo = vporo * my::funct_(ui);

        for (int idof=0; idof<numphases; ++idof)
        {
          const int fui = ui*numphases+idof;

          emat(fvi,fui) += vfunctsat*VarManager()->FluidPhaseManager()->SaturationDeriv(curphase,idof)
                         + vfunctporo*VarManager()->FluidPhaseManager()->PorosityDeriv(idof);
        }
      }
    }
  }
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
