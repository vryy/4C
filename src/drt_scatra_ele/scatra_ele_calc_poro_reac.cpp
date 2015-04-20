/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro_reac.cpp

 \brief

 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/structporo_reaction_ecm.H"
#include "../drt_mat/scatra_mat.H"

#include "scatra_ele_parameter.H"
#include "scatra_ele_calc_poro_reac.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode,numscal),
    DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal)
{
  // safety check
  if(not my::scatrapara_->TauGP())
    dserror("For poro reactions, tau needs to be evaluated by integration-point evaluations!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype> * DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcPoroReac<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcPoroReac<distype>(numdofpernode,numscal);
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
  //call poro base class to compute porosity
  poro::ComputePorosity(ele);

  //call advreac base class to evaluate porosity
  advreac::GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  switch(material->MaterialType())
  {
  case INPAR::MAT::m_scatra:
    MatScaTra(material,k,densn,densnp,densam,visc,iquad);
    break;
  case INPAR::MAT::m_scatra_poroECM:
    MatPoroECM(material,k,densn,densnp,densam,visc,iquad);
    break;
  default:
    dserror("Material type %i is not supported",material->MaterialType());
   break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::MatScaTra(
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
    dserror("no gauss point given for evaluation of scatra material. Check your input file.");

  const double porosity = poro::DiffManager()->GetPorosity();

  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  if(actmat == Teuchos::null)
    dserror("cast to ScatraMat failed");

  // set diffusivity (scaled with porosity)
  poro::SetDiffusivity(actmat,k,porosity);

  // set/calculate reaction coefficient
  // do not (!) scale with porosity, as reactive term is scaled with density (=porosity) before assembly
  //Todo
//  if (not advreac::iscoupled_)
//    poro::SetReaCoefficient(actmat,k,1.0);
//  else
    SetReactionTermsMatScatra(k,1.0);

  // set densities (scaled with porosity)
  poro::SetDensities(porosity,densn,densnp,densam);

  //NOTE: The reaction stuff happens in function SetAdvancedReactionTerms(...)
  return;
} // ScaTraEleCalcPoroReac<distype>::MatScaTra

/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::MatPoroECM(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  MatScaTra(material,k,densn,densnp,densam,visc,iquad);

  //Todo: clean up
  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding fluid element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(my::eid_);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", my::eid_);

  const Teuchos::RCP<const MAT::StructPoroReactionECM>& structmat
            = Teuchos::rcp_dynamic_cast<const MAT::StructPoroReactionECM>(structele->Material());
  if(structmat == Teuchos::null)
    dserror("invalid structure material for reactive poroelasticity model for ECM");

  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = advreac::ReaManager();

  const double porosity = poro::DiffManager()->GetPorosity();
  const double bodyforce = structmat->BodyForceTerm(porosity);
  const double bodyforce_old = reamanageradvreac->GetReaBodyForce(k);
  reamanageradvreac->SetReaBodyForce(bodyforce_old+bodyforce,k);

  return;
} // ScaTraEleCalcPoroReac<distype>::MatScaTra

/*-------------------------------------------------------------------------------*
 |  set body force, reaction coefficient and derivatives          vuong 06/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::SetReactionTermsMatScatra(
    const int                               k,
    const double                            porosity
                                    )
{
  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = advreac::ReaManager();
  if(reamanageradvreac==Teuchos::null) dserror("cast to ScaTraEleReaManagerAdvReac failed");

  reamanageradvreac->SetReaBodyForce( CalcReaBodyForceTerm(k,porosity) ,k);
  reamanageradvreac->SetReaCoeff( CalcReaCoeff(k,porosity) ,k);
  for (int j=0; j<my::numscal_ ;j++)
  {
    reamanageradvreac->SetReaBodyForceDerivMatrix( CalcReaBodyForceDerivMatrix(k,j,porosity) ,k,j );
    reamanageradvreac->SetReaCoeffDerivMatrix( CalcReaCoeffDerivMatrix(k,j,porosity) ,k,j );
  }
}

/*----------------------------------------------------------------------*
 |  Calculate K(c)                                           thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaCoeff(
    const int                        k,                       //!< id of current scalar
    const double                     porosity                 //!< current porosity
)
{
  double reactermK=0;

  for (int condnum = 0; condnum < advreac::numcond_ ; condnum++)
  {
    const std::vector<int>& stoich = advreac::stoich_[condnum]; //get stoichometrie
    const MAT::PAR::reaction_coupling couplingtype = advreac::couplingtype_[condnum]; //get coupling type
    const double& reaccoeff = advreac::reaccoeff_[condnum]; //get reaction coefficient
    const double& reacstart = advreac::reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] < 0)
    {
      double rcfac= CalcReaCoeffFac(stoich,couplingtype,k,porosity);

      if (advreac::reacstart_[condnum]>0 and rcfac>0)
        advreac::ReacStartForReaCoeff(k,condnum,reacstart,rcfac);

      reactermK += -reaccoeff*stoich[k]*rcfac; // scalar at integration point np
    }
  }
  return reactermK;
}

/*----------------------------------------------------------------------*
 |  helper for calculating K(c)                              thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaCoeffFac(
          const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
          const MAT::PAR::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
          const int                                 k,                       //!< id of current scalar
          const double                              porosity                //!< current porosity
)
{
  double rcfac=1.0;
  bool allpositive = true;

  for (int ii=0; ii < my::numscal_; ii++)
  {
    if (stoich[ii]<0)
    {
      allpositive = false;
      if (ii!=k)
      {
        switch (couplingtype)
        {
        case MAT::PAR::reac_coup_simple_multiplicative:
          rcfac *=my::funct_.Dot(my::ephinp_[ii])*porosity;
          break;
        //case ... :  //insert new Couplings here
        default:
          dserror("invalid reaction_coupling type");
          break;
        }
      }
//      //TODO: HACK!!!!!
//      else if(k==0)
//        rcfac *=10.1*(1-porosity)*porosity*porosity;
    }
  }
  if (allpositive)
    dserror("there must be at least one negative entry in each stoich list");

  return rcfac;
}

/*----------------------------------------------------------------------*
 |  calculate \frac{partial}{\partial c} K(c)                thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaCoeffDerivMatrix(
    const int                 k,                  //!< id of current scalar
    const int                 j,                   //!< concentration to be derived to
    const double              porosity                //!< current porosity
)
{
  double reacoeffderivmatrixKJ=0;

  for (int condnum = 0; condnum < advreac::numcond_; condnum++)
  {
    const std::vector<int>& stoich = advreac::stoich_[condnum]; //get stoichometrie
    const MAT::PAR::reaction_coupling& couplingtype = advreac::couplingtype_[condnum]; //get coupling type
    const double& reaccoeff = advreac::reaccoeff_[condnum]; //get reactioncoefficient
    const double& reacstart = advreac::reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] < 0)
    {
      double rcdmfac = CalcReaCoeffDerivFac(stoich,couplingtype,j,k,porosity);

      if (advreac::reacstart_[condnum]>0)
        advreac::ReacStartForReaCoeffDeriv(k,j,condnum,reacstart,rcdmfac,stoich,couplingtype);

      reacoeffderivmatrixKJ += -reaccoeff*stoich[k]*rcdmfac;
    } //end if(stoich[k] != 0)
  }
  return reacoeffderivmatrixKJ;
}

/*----------------------------------------------------------------------*
 |  helper for calculating \frac{partial}{\partial c} K(c)   thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaCoeffDerivFac(
          const std::vector<int>                  stoich,                  //!<stoichometrie of current condition
          const MAT::PAR::reaction_coupling  couplingtype,            //!<type of coupling the stoichometry coefficients
          const int                               toderive,                //!<concentration to be derived to
          const int                               k,                       //!< id of current scalar
          const double                            porosity                //!< current porosity
)
{
  double rcdmfac=1;

  if (stoich[toderive]<0 and toderive!=k)
  {
    for (int ii=0; ii < my::numscal_; ii++)
    {
      if (stoich[ii]<0)
      {
        switch (couplingtype)
        {
        case MAT::PAR::reac_coup_simple_multiplicative:
          if (ii!=k and ii!= toderive)
            rcdmfac *= my::funct_.Dot(my::ephinp_[ii])*porosity;
          break;
        //case ... :  //insert new Couplings here
        default:
          dserror("invalid reaction_coupling type");
          break;
        }
      }
    }
  }
  else
    rcdmfac = 0;

  return rcdmfac;
}

/*----------------------------------------------------------------------*
 |  calculate f(c)                                           thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaBodyForceTerm(
    const int                                 k,                       //!< id of current scalar
    const double                              porosity                //!< current porosity
)
{
  double bodyforcetermK=0.0;

  for (int condnum = 0; condnum < advreac::numcond_; condnum++)
  {
    const std::vector<int>& stoich = advreac::stoich_[condnum]; //get stoichometrie
    const MAT::PAR::reaction_coupling couplingtype = advreac::couplingtype_[condnum]; //get coupling type
    const double reaccoeff = advreac::reaccoeff_[condnum]; //get reactioncoefficient
    const double& reacstart = advreac::reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] > 0)
    {
      double bftfac = CalcReaBodyForceTermFac(stoich,couplingtype,porosity);// scalar at integration point np

      if (advreac::reacstart_[condnum]>0 and bftfac>0)
        advreac::ReacStartForReaBF(k,condnum,reacstart,bftfac);

      bodyforcetermK += reaccoeff*stoich[k]*bftfac;
    }
  }
  return bodyforcetermK;
}

/*----------------------------------------------------------------------*
 |  helper for calculating                                   thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaBodyForceTermFac(
          const std::vector<int>                      stoich,                 //!<stoichometrie of current condition
          const MAT::PAR::reaction_coupling      couplingtype,            //!<type of coupling the stoichometry coefficients
          const double                                porosity                //!< current porosity
)
{
  double bftfac=1.0;
  bool allpositive = true;

  for (int ii=0; ii < my::numscal_; ii++)
  {
    if (stoich[ii]<0)
    {
      allpositive = false;
      switch (couplingtype)
      {
      case MAT::PAR::reac_coup_simple_multiplicative:
        bftfac *=my::funct_.Dot(my::ephinp_[ii])*porosity;
        break;
      //case ... :  //insert new Couplings here
      default:
        dserror("invalid reaction_coupling type");
        break;
      }
    }
  }
  if (allpositive)
    dserror("there must be at least one negative entry in each stoich list");

  return bftfac;
}

/*----------------------------------------------------------------------*
 |  calculate \frac{partial}{\partial c} f(c)                thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaBodyForceDerivMatrix(
    const int                 k,                  //!< id of current scalar
    const int                 j,                   //!< concentration to be derived to
    const double              porosity                //!< current porosity
)
{
  double reabodyforcederivmatrixKJ=0;
  for (int condnum = 0; condnum < advreac::numcond_; condnum++)
  {
    //reading of conditions here, because of future implemantion of nonhomogeneous couplings
    const std::vector<int>& stoich = advreac::stoich_[condnum]; //get stoichometrie
    const MAT::PAR::reaction_coupling couplingtype = advreac::couplingtype_[condnum]; //get coupling type
    const double& reaccoeff = advreac::reaccoeff_[condnum]; //get reaction coefficient
    const double& reacstart = advreac::reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] > 0)
    {
      double bfdmfac = CalcReaBodyForceDerivFac(stoich,couplingtype,j,porosity);

      if (advreac::reacstart_[condnum]>0 and bfdmfac>0)
        advreac::ReacStartForReaBFDeriv(k,j,condnum,reacstart,bfdmfac,stoich,couplingtype);

      reabodyforcederivmatrixKJ += reaccoeff*stoich[k]*bfdmfac;
    }
  }
  return reabodyforcederivmatrixKJ;
}

/*-------------------------------------------------------------------------------*
 |  helper for calculating calculate \frac{partial}{\partial c} f(c)  thon 02/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::CalcReaBodyForceDerivFac(
        const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
        const MAT::PAR::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
        const int                                 toderive,                 //!<concentration to be derived to
        const double                              porosity                //!< current porosity
)
{
  double bfdmfac=1.0;
  if (stoich[toderive]<0)
  {
    for (int ii=0; ii < my::numscal_; ii++)
    {
      if (stoich[ii]<0)
      {
        switch (couplingtype)
        {
        case MAT::PAR::reac_coup_simple_multiplicative:
          if (ii!=toderive)
            bfdmfac *= my::funct_.Dot(my::ephinp_[ii])*porosity;
          break;
        //case ... :  //insert new Couplings here
        default:
          dserror("invalid reaction_coupling type");
          break;
        }
      }
    }
  }
  else
    bfdmfac = 0.0;

  return bfdmfac;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values                     vuong 04/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const std::vector<double> DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ExtractElementAndNodeValues(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm
)
{
  const Teuchos::RCP<Epetra_MultiVector> pre = params.get< Teuchos::RCP<Epetra_MultiVector> >("pressure field");
  LINALG::Matrix<1,my::nen_> eprenp;
  DRT::UTILS::ExtractMyNodeBasedValues(ele,eprenp,pre,1);

  //pressure values
  for (int i=0;i<my::nen_;++i)
  {
    my::eprenp_(i) = eprenp(0,i);
  }

  Teuchos::RCP<const Epetra_Vector> disp= discretization.GetState(1,"displacement");

  return poro::ExtractElementAndNodeValues(ele,params,discretization,lm);
}

// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs27>;

