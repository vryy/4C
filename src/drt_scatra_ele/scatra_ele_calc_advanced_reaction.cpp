/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.cpp

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_advanced_reaction.H"

#include "scatra_ele_parameter.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/biofilm.H"
#include "../drt_mat/scatra_growth_scd.H"
#include "../drt_mat/growth_scd.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype> * DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcAdvReac<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcAdvReac<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *  constructor---------------------------                   thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
  iscoupled_(false),
  isinit_(false),
  numcond_(0),
  stoich_(0,std::vector<int>(my::numscal_,0.0)),
  reaconst_(0,0.0),
  couplingtype_(0,DRT::ELEMENTS::simple_multiplicative),
  reacstart_(0,0.0)
{
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerAdvReac(my::numscal_));
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Materials(
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
    MatScaTra(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
    break;
  case INPAR::MAT::m_biofilm:
    MatBioFilm(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
    break;
  case INPAR::MAT::m_scatra_growth_scd:
    MatGrowthScd(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
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
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  my::MatScaTra(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);

  if (iscoupled_)
  {
    SetAdvancedReactionTerms(reamanager,k);
  }
}

/*----------------------------------------------------------------------*
 |  Material BioFilm                                         thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::MatBioFilm(
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
  const Teuchos::RCP<const MAT::Biofilm>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Biofilm>(material);

  // get constant diffusivity
  diffmanager->SetIsotropicDiff(actmat->Diffusivity(),k);

  // get substrate concentration at n+1 or n+alpha_F at integration point
  const double csnp = my::funct_.Dot(my::ephinp_[k]);

  // set reaction coefficient
  reamanager->SetReaCoeff(actmat->ComputeReactionCoeff(csnp),k);
  // set derivative of reaction coefficient
  reamanager->SetReaCoeffDerivMatrix(actmat->ComputeReactionCoeffDeriv(csnp),k,k);

  // set density at various time steps and density gradient factor to 1.0/0.0
  densn      = 1.0;
  densnp     = 1.0;
  densam     = 1.0;

  return;
}

/*----------------------------------------------------------------------*
 |  Material GrowthScd                                       vuong 01/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::MatGrowthScd(
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
  dsassert(my::numdofpernode_==1,"more than 1 dof per node for ScatraGrowthScd material");

  if(iquad < 0) dserror("ScatraGrowthScd material has to be evaluated at gauss point!");

  const Teuchos::RCP<const MAT::ScatraGrowthScd>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraGrowthScd>(material);

  //strategy to obtain theta from the structure at equivalent gauss-point
  //access structure discretization
   Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
   structdis = DRT::Problem::Instance()->GetDis("structure");
   //get corresponding structure element (it has the same global ID as the scatra element)
   DRT::Element* structele = structdis->gElement(my::eid_);
   if (structele == NULL)
     dserror("Structure element %i not on local processor", my::eid_);

   const Teuchos::RCP<const MAT::GrowthScd>& structmat
             = Teuchos::rcp_dynamic_cast<const MAT::GrowthScd>(structele->Material());
   if(structmat->MaterialType() != INPAR::MAT::m_growthscd)
     dserror("invalid structure material for scalar dependent growth");

   const double theta    = structmat->Gettheta_atgp(iquad);
   const double dtheta   = structmat->Getdtheta_atgp(iquad);
   const double thetaold = structmat->Getthetaold_atgp(iquad);
   const double detFe    = structmat->GetdetFe_atgp(iquad);

  // get constant diffusivity
  diffmanager->SetIsotropicDiff(actmat->Diffusivity(),k);

  // get substrate concentration at n+1 or n+alpha_F at integration point
  const double csnp = my::funct_.Dot(my::ephinp_[k]);

  // set reaction coefficient
  reamanager->SetReaCoeff(actmat->ComputeReactionCoeff(csnp,theta,dtheta,detFe),k);
  // set derivative of reaction coefficient
  reamanager->SetReaCoeffDerivMatrix(actmat->ComputeReactionCoeffDeriv(csnp,theta,thetaold,1.0),k,k);

  // set density at various time steps and density gradient factor to 1.0/0.0
  densn      = 1.0;
  densnp     = 1.0;
  densam     = 1.0;

  return;
}

/*----------------------------------------------------------------------*
 |  Get right hand side including reaction bodyforce term    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::GetRhsInt(
  double&      rhsint,   //!< rhs containing bodyforce at Gauss point
  const double densnp,  //!< density at t_(n+1)
  const int    k        //!< index of current scalar
  )
{
  // compute rhs containing bodyforce (divided by specific heat capacity) and,
  // for temperature equation, the time derivative of thermodynamic pressure,
  // if not constant, and for temperature equation of a reactive
  // equation system, the reaction-rate term


  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(my::reamanager_);

                                       //... + reaction terms not depending on phi(k) -> source term
  rhsint = my::bodyforce_[k].Dot(my::funct_) + densnp*reamanageradvreac->GetReaBodyForce(k);

  return;
}
/*----------------------------------------------------------------------*
 |  check if is coupled and if, read inputs--------------    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::IsCoupledAndRead(
  const DRT::Discretization&            ScaTraDiscretization //discretisation of the ScaTra field
  )
{
  //! note for the implementation of the homogenous scatra coupling:
  //! assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0
  //!
  //! if we assume the reaction is depending on the product of all
  //! reactants (this corresponds to couplingtype "simple_multiplicative"),
  //! the corresponding equations are: \partial_t A = -(4*1*B)*A  (negative since reactant)
  //!                              \partial_t B = -(4*2*A)*B  (negative since reactant)
  //!                              \partial_t C = + 4*3*A*B   (positive since product)
  //!
  //! this equation is in BACI achieved by the boundary condition:
  //! ------------------------DESIGN HOMOGENEOUS SCATRA COUPLING VOLUME CONDITIONS
  //! DVOL                           XYZ
  //! E XYZ - numscal 3 stoich -1 -2 3 reaccoeff 4.0 coupling simple_multiplicative
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to depend on c_i
  //! hence we have to calculate and set K(c)=(4*B;8*A;0) and f(c)=(0;0;12*A*B) and corresponding derivatives.



  //! note for the implementation of the homogenous scatra coupling reacstart feature:
  //! Assume concentration A is reproducing with reaction coefficient 1.0 and if the concentration
  //! exeeds some threshold 2.0 if starts to react A->3*B with reacion coefficient 4.0.
  //!
  //! the corresponding equations are:
  //!            \partial_t A = -(-1.0)*A - 4.0*(A - 2.0)_{+} (first termn postive, since equivalent as reactant with negative reaction coefficient)
  //!            \partial_t B = 3.0*4.0 (A - 2.0)_{+}   (positive since product)
  //!
  //! this equation is in BACI achieved by the boundary condition:
  //! ------------------------DESIGN HOMOGENEOUS SCATRA COUPLING VOLUME CONDITIONS
  //! DVOL                           XYZ
  //! E XYZ - numscal 2 stoich -1 0 reaccoeff -1.0 coupling simple_multiplicative
  //! E XYZ - numscal 2 stoich -1 3 reaccoeff 4.0 coupling simple_multiplicative reacstart 2.0
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to depend on c_i
  //! hence we have to calculate and set K(c)=(-A + 4*(A-2)_{+};0) and f(c)=(0;12*(A-2)_{+}) and corresponding derivatives.


  if (isinit_ == true)
    dserror("this function should be called just once...");

  std::vector<DRT::Condition*> HSTCConds;
  ScaTraDiscretization.GetCondition("HomoScaTraCoupling",HSTCConds);

  numcond_ = HSTCConds.size();

  if (numcond_ > 0)  //if there exists condition "DESIGN HOMOGENEOUS SCATRA COUPLING VOLUME CONDITIONS"
  {
    stoich_.resize(numcond_);
    couplingtype_.resize(numcond_);
    reaconst_.resize(numcond_);
    reacstart_.resize(numcond_);

    for (int i=0;i<numcond_;i++)
    {
    stoich_[i]=*HSTCConds[i]->GetMutable<std::vector<int> >("stoich");

    reaconst_[i]=HSTCConds[i]->GetDouble("reaccoeff");
//    to allow self-growing concentrations =>negativ reaccoeff needed in current implementation
//    if (reaconst_[i]<0)
//        dserror("reaccoeff of Condition %d is negativ",i);

    if ( *HSTCConds[i]->Get<std::string>("coupling") == "simple_multiplicative")
      couplingtype_[i] = DRT::ELEMENTS::simple_multiplicative;
    //else if  insert new couplings here
    else
      dserror("couplingtype OTHER is just a dummy and hence not implemented");

    reacstart_[i] =HSTCConds[i]->GetDouble("reacstart");

    }//end for
    return true; //1;
  } // end if
  else
  {
    return false;//0;
  }
}

/*----------------------------------------------------------------------*
 |  Calculate K(c)                                           thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaCoeff(
    const int                        k                       //!< id of current scalar
)
{
  double reactermK=0;
  //double reactermK=0;
  //std::cout <<"---------------IsCoupled-Loop!\n";
  for (int condnum = 0; condnum < numcond_ /*HSTCConds_.size()*/; condnum++)
  {
    const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
    const DRT::ELEMENTS::reaction_coupling couplingtype = couplingtype_[condnum]; //get coupling type
    const double reaccoeff = reaconst_[condnum]; //get reactioncoefficient    reaccoeff = HSTCConds_[condnum-1]->GetDouble("reaccoeff");

    if (stoich[k] < 0)
    {
      double rcfac= CalcReaCoeffFac(stoich,couplingtype,k);

      if (reacstart_[condnum]>0 and rcfac>0)
        ReacStartForReaCoeff(k,condnum,rcfac);

      reactermK += -reaccoeff*stoich[k]*rcfac; // scalar at integration point np
    }
  }
  return reactermK;
}

/*----------------------------------------------------------------------*
 |  helper for calculating K(c)                              thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaCoeffFac(
          const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
          const DRT::ELEMENTS::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
          const int                                 k                       //!< id of current scalar
)
{
  double rcfac=1;
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
        case DRT::ELEMENTS::simple_multiplicative:
          rcfac *=my::funct_.Dot(my::ephinp_[ii]);
          break;
        //case ... :  //insert new Couplings here
        }
      }
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
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaCoeffDerivMatrix(
    const int                 k,                  //!< id of current scalar
    const int                 j                   //!< concentration to be derived to
)
{
  double reacoeffderivmatrixKJ=0;

  for (int condnum = 0; condnum < numcond_; condnum++)
  {
    const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
    const DRT::ELEMENTS::reaction_coupling& couplingtype = couplingtype_[condnum]; //get coupling type
    const double& reaccoeff = reaconst_[condnum]; //get reactioncoefficient

    if (stoich[k] < 0)
    {
      double rcdmfac = CalcReaCoeffDerivFac(stoich,couplingtype,j,k);

      if (reacstart_[condnum]>0)
        ReacStartForReaCoeffDeriv(k,j,condnum,rcdmfac,stoich,couplingtype);

      reacoeffderivmatrixKJ += -reaccoeff*stoich[k]*rcdmfac;
    } //end if(stoich[k] != 0)
  }
  return reacoeffderivmatrixKJ;
}

/*----------------------------------------------------------------------*
 |  helper for calculating \frac{partial}{\partial c} K(c)   thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaCoeffDerivFac(
          const std::vector<int>                  stoich,                  //!<stoichometrie of current condition
          const DRT::ELEMENTS::reaction_coupling  couplingtype,            //!<type of coupling the stoichometry coefficients
          const int                               toderive,                //!<concentration to be derived to
          const int                               k                       //!< id of current scalar
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
        case DRT::ELEMENTS::simple_multiplicative:
          if (ii!=k and ii!= toderive)
            rcdmfac *= my::funct_.Dot(my::ephinp_[ii]);
          break;
        //case ... :  //insert new Couplings here
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
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaBodyForceTerm(
    const int                                 k                       //!< id of current scalar
)
{
  double bodyforcetermK=0;

  for (int condnum = 0; condnum < numcond_; condnum++)
  {
    const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
    const DRT::ELEMENTS::reaction_coupling couplingtype = couplingtype_[condnum]; //get coupling type
    const double reaccoeff = reaconst_[condnum]; //get reactioncoefficient    reaccoeff = HSTCConds_[condnum-1]->GetDouble("reaccoeff");

    if (stoich[k] > 0)
    {
      double bftfac = CalcReaBodyForceTermFac(stoich,couplingtype);// scalar at integration point np

      if (reacstart_[condnum]>0 and bftfac>0)
        ReacStartForReaBF(k,condnum,bftfac);

      bodyforcetermK += reaccoeff*stoich[k]*bftfac;
    }
  }
  return bodyforcetermK;
}

/*----------------------------------------------------------------------*
 |  helper for calculating                                   thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaBodyForceTermFac(
          const std::vector<int>                      stoich,                 //!<stoichometrie of current condition
          const DRT::ELEMENTS::reaction_coupling      couplingtype            //!<type of coupling the stoichometry coefficients
)
{
  double bftfac=1;
  bool allpositive = true;
  for (int ii=0; ii < my::numscal_; ii++)
  {
    if (stoich[ii]<0)
    {
      allpositive = false;
      switch (couplingtype)
      {
      case DRT::ELEMENTS::simple_multiplicative:
        bftfac *=my::funct_.Dot(my::ephinp_[ii]);
        break;
      //case ... :  //insert new Couplings here
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
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaBodyForceDerivMatrix(
    const int                 k,                  //!< id of current scalar
    const int                 j                   //!< concentration to be derived to
)
{
  double reabodyforcederivmatrixKJ=0;
  for (int condnum = 0; condnum < numcond_; condnum++)
  {
    //reading of conditions here, because of future implemantion of nonhomogeneous couplings
    const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
    const DRT::ELEMENTS::reaction_coupling couplingtype = couplingtype_[condnum]; //get coupling type
    const double reaccoeff = reaconst_[condnum]; //get reactioncoefficient

    if (stoich[k] > 0)
    {
      double bfdmfac = CalcReaBodyForceDerivFac(stoich,couplingtype,j);

      if (reacstart_[condnum]>0 and bfdmfac>0)
        ReacStartForReaBFDeriv(k,j,condnum,bfdmfac,stoich,couplingtype);

      reabodyforcederivmatrixKJ += reaccoeff*stoich[k]*bfdmfac;
    }
  }
  return reabodyforcederivmatrixKJ;
}

/*-------------------------------------------------------------------------------*
 |  helper for calculating calculate \frac{partial}{\partial c} f(c)  thon 02/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaBodyForceDerivFac(
        const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
        const DRT::ELEMENTS::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
        const int                                 toderive                 //!<concentration to be derived to
)
{
  double bfdmfac=1;
  if (stoich[toderive]<0)
  {
    for (int ii=0; ii < my::numscal_; ii++)
    {
      if (stoich[ii]<0)
      {
        switch (couplingtype)
        {
        case DRT::ELEMENTS::simple_multiplicative:
          if (ii!=toderive)
            bfdmfac *= my::funct_.Dot(my::ephinp_[ii]);
          break;
        //case ... :  //insert new Couplings here
        }
      }
    }
  }
  else
    bfdmfac = 0;

  return bfdmfac;
}

/*-------------------------------------------------------------------------------*
 |  calculate reaction coefficient for "delayed" reactions            thon 03/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ReacStartForReaCoeff(
  const int                               k,            //!< id of current scalar
  const int                               condnum,      //!< id of current condition
  double&                                 value         //!< current reaction value
  )
{
  double prod = value * my::funct_.Dot(my::ephinp_[k]); //for simple multiplikative only!

    if (prod > reacstart_[condnum] )
      value = value - reacstart_[condnum]/my::funct_.Dot(my::ephinp_[k]);
    else
      value = 0;

    return;
}

///*----------------------------------------------------------------------------------*
// |  calculate reaction coefficient derivative for "delayed" reactions    thon 03/14 |
// *----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ReacStartForReaCoeffDeriv(
        const int                               k,              //!< id of current scalar
        const int                               toderive,       //!<concentration to be derived to
        const int                               condnum,        //!< id of current condition
        double&                                 value,          //!< current reaction value
        const std::vector<int>                  stoich,         //!<stoichometrie of current condition
        const DRT::ELEMENTS::reaction_coupling  couplingtype    //!<type of coupling the stoichometry coefficients
  )
{
  double prod = CalcReaBodyForceTermFac(stoich,couplingtype);

  if ( prod > reacstart_[condnum])
  {
    if (k==toderive)
      value= value - (-reacstart_[condnum] / pow(my::funct_.Dot(my::ephinp_[k]),2) );
  }
  else
    value = 0;
  //std::cout<<__LINE__<<__FILE__<<"\t"<<value<<endl;
  return;
}

/*-------------------------------------------------------------------------------*
 |  calculate reactions body force term for "delayed" reactions       thon 03/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ReacStartForReaBF(
  const int                               k,            //!< id of current scalar
  const int                               condnum,      //!< id of current condition
  double&                                 value         //!< current reaction value
  )
{
  if (value > reacstart_[condnum] )
    value = value - reacstart_[condnum];
  else
    value = 0;

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  calculate reactions body force term derivative for "delayed" reactions      thon 03/14 |
 *-----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ReacStartForReaBFDeriv(
        const int                               k,              //!< id of current scalar
        const int                               toderive,       //!<concentration to be derived to
        const int                               condnum,        //!< id of current condition
        double&                                 value,          //!< current reaction value
        const std::vector<int>                  stoich,         //!<stoichometrie of current condition
        const DRT::ELEMENTS::reaction_coupling  couplingtype    //!<type of coupling the stoichometry coefficients
  )
{
  double prod = CalcReaBodyForceTermFac(stoich,couplingtype);

  if ( prod > reacstart_[condnum]) { }
  else
    value = 0;

  return;
}


/*-------------------------------------------------------------------------------*
 |  Evaluate including check for homogeneous ScaTra coupling          thon 02/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Evaluate(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra
  )
{

  if (not isinit_)
  {
    iscoupled_=IsCoupledAndRead(discretization);
    isinit_=true;
  }
  return my::Evaluate(
      ele,
      params,
      discretization,
      lm,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
}

/*--------------------------------------------------------------------------- *
 |  calculation of reactive element matrix for coupled reactions  thon 02/14  |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcMatReact(
  Epetra_SerialDenseMatrix&          emat,
  const int                          k,
  const double                       timefacfac,
  const double                       timetaufac,
  const double                       taufac,
  const double                       densnp,
  const double                       phinp,
  Teuchos::RCP<ScaTraEleReaManager>  reamanager,
  const LINALG::Matrix<my::nen_,1>&      conv,
  const LINALG::Matrix<my::nen_,1>&      sgconv,
  const LINALG::Matrix<my::nen_,1>&      diff
  )
{
  // -----------------first care for Term K(c)*(\partial_c c)=Id*K(c)--------------------------------------

  my::CalcMatReact(emat,k,timefacfac,timetaufac,taufac,densnp,phinp,reamanager,conv,sgconv,diff);


  // -----------------second care for Term (\partial_c K(c)) .* c - (\partial_c f_{reabody}(c))------------

  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(reamanager);

  for (int j=0; j<my::numscal_ ;j++)
  {
    const double fac_reac        = timefacfac*densnp*( reamanager->GetReaCoeffDerivMatrix(k,j)*phinp - reamanageradvreac->GetReaBodyForceDerivMatrix(k,j) );
    const double timetaufac_reac = timetaufac*densnp*reamanager->GetReaCoeffDerivMatrix(k,j)*phinp;

    //----------------------------------------------------------------
    // standard Galerkin reactive term
    //----------------------------------------------------------------
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac_reac*my::funct_(vi);
      const int fvi = vi*my::numdofpernode_+k;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+j;

        emat(fvi,fui) += v*my::funct_(ui);
      }
    }

    //----------------------------------------------------------------
    // stabilization of reactive term
    //----------------------------------------------------------------
    if(my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
    {
      double densreataufac = timetaufac_reac*densnp;
      // convective stabilization of reactive term (in convective form)
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double v = densreataufac*(conv(vi)+sgconv(vi)+my::scatrapara_->USFEMGLSFac()*1.0/my::scatraparatimint_->TimeFac()*my::funct_(vi));
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+j;

          emat(fvi,fui) += v*my::funct_(ui);
        }
      }

      if (my::use2ndderiv_)
      {
        // diffusive stabilization of reactive term
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = my::scatrapara_->USFEMGLSFac()*timetaufac_reac*diff(vi);
          const int fvi = vi*my::numdofpernode_+k;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+j;

            emat(fvi,fui) -= v*my::funct_(ui);
          }
        }
      }

      //----------------------------------------------------------------
      // reactive stabilization
      //----------------------------------------------------------------
      densreataufac = my::scatrapara_->USFEMGLSFac()*timetaufac_reac*densnp;

      // reactive stabilization of convective (in convective form) and reactive term
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double v = densreataufac*my::funct_(vi);
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+j;

          emat(fvi,fui) += v*(conv(ui)+reamanager->GetReaCoeff(k)*my::funct_(ui));
        }
      }

      if (my::use2ndderiv_)
      {
        // reactive stabilization of diffusive term
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = my::scatrapara_->USFEMGLSFac()*timetaufac_reac*my::funct_(vi);
          const int fvi = vi*my::numdofpernode_+k;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+j;

            emat(fvi,fui) -= v*diff(ui);
          }
        }
      }


      if (not my::scatraparatimint_->IsStationary())
      {
        // reactive stabilization of transient term
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = my::scatrapara_->USFEMGLSFac()*taufac*densnp*reamanager->GetReaCoeff(k)*densnp*my::funct_(vi);
          const int fvi = vi*my::numdofpernode_+k;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+j;

            emat(fvi,fui) += v*my::funct_(ui);
          }
        }

        if (my::use2ndderiv_ and reamanager->GetReaCoeff(k)!=0.0)
          dserror("Second order reactive stabilization is not fully implemented!! ");
      }
    }
  } //end for
  return;
}

/*-------------------------------------------------------------------------------*
 |  set body force, reaction coefficient and derivatives          vuong 06/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::SetAdvancedReactionTerms(
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,
    const int                               k,
    const double                            scale
                                    )
{
  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(reamanager);
  if(reamanageradvreac==Teuchos::null)
    dserror("cast to ScaTraEleReaManagerAdvReac failed");

  reamanageradvreac->SetReaBodyForce( CalcReaBodyForceTerm(k)*scale ,k);
  reamanageradvreac->SetReaCoeff( CalcReaCoeff(k)*scale ,k);
  for (int j=0; j<my::numscal_ ;j++)
  {
    reamanageradvreac->SetReaBodyForceDerivMatrix( CalcReaBodyForceDerivMatrix(k,j)*scale ,k,j );
    reamanager->SetReaCoeffDerivMatrix( CalcReaCoeffDerivMatrix(k,j)*scale ,k,j );
  }
}

// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs27>;
