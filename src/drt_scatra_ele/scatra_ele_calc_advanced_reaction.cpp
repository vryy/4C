/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.cpp

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_advanced_reaction.H"

#include "scatra_ele_parameter.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

//MATERIALS
#include "../drt_mat/biofilm.H"
#include "../drt_mat/scatra_growth_scd.H"
#include "../drt_mat/growth_scd.H"
#include "../drt_mat/growth_law.H"
#include "../drt_mat/matlist_reactions.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matlist.H"

  //! note for the implementation of the homogenous scatra coupling:
  //! assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0
  //!
  //! if we assume the reaction is depending on the product of all
  //! reactants (this corresponds to couplingtype "simple_multiplicative"),
  //! the corresponding equations are: \partial_t A = -(4*1*B)*A  (negative since reactant)
  //!                              \partial_t B = -(4*2*A)*B  (negative since reactant)
  //!                              \partial_t C = + 4*3*A*B   (positive since product)
  //!
  //! this equation is in BACI achieved by the MAT_scatra_reaction material:
  //! ----------------------------------------------------------MATERIALS
  //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 1 MATIDS 2 NUMREAC 1 REACIDS 3 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -2 3 REACOEFF 4.0 COUPLING simple_multiplicative
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to depend on c_i
  //! hence we have to calculate and set K(c)=(4*B;8*A;0) and f(c)=(0;0;12*A*B) and corresponding derivatives.



  //! note for the implementation of the homogenous scatra coupling reacstart feature:
  //! Assume concentration A is reproducing with reaction coefficient 1.0 and if the concentration
  //! exceeds some threshold 2.0 if starts to react A->3*B with reacion coefficient 4.0.
  //!
  //! the corresponding equations are:
  //!            \partial_t A = -(-1.0)*A - 4.0*(A - 2.0)_{+} (first termn postive, since equivalent as reactant with negative reaction coefficient)
  //!            \partial_t B = 3.0*4.0 (A - 2.0)_{+}   (positive since product)
  //!
  //! this equation is in BACI achieved by the boundary condition:
  //!   //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 1 MATIDS 2 NUMREAC 2 REACIDS 3 4 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 2 STOICH -1 0 REACOEFF -1.0 COUPLING simple_multiplicative
  //! MAT 4 MAT_scatra_reaction NUMSCAL 2 STOICH -1 3 REACOEFF 4.0 COUPLING simple_multiplicative REACSTART 2.0
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to depend on c_i
  //! hence we have to calculate and set K(c)=(-A + 4*(A-2)_{+};0) and f(c)=(0;12*(A-2)_{+}) and corresponding derivatives.

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
    numcond_(-1),
    stoich_(0),
    reaccoeff_(0),
    couplingtype_(MAT::PAR::reac_coup_none),
    reacstart_(0)
{
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerAdvReac(my::numscal_));
}

/*----------------------------------------------------------------------*
 |  get the material constants  (private)                      thon 09/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
//// get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();


  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {

    const Teuchos::RCP<const MAT::MatList>& actmat = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat,k,densn,densnp,densam,visc,iquad);
    }
  }
  else if (material->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {

    const Teuchos::RCP<const MAT::MatListReactions>& actmat = Teuchos::rcp_dynamic_cast<const MAT::MatListReactions>(material);
    if (actmat->NumMat() < my::numscal_) dserror("Not enough materials in MatList.");

    GetAdvancedReactionCoefficients(actmat); // read all reaction input from material and copy it into local variables

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat,k,densn,densnp,densam,visc,iquad);

      SetAdvancedReactionTerms(k,1.0); //every reaction calculation stuff happens in here!!
    }

  }
  else
  {
    Materials(material,0,densn,densnp,densam,visc,iquad);
  }
  return;
} //ScaTraEleCalc::GetMaterialParams

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
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  switch(material->MaterialType())
  {
  case INPAR::MAT::m_scatra:
    my::MatScaTra(material,k,densn,densnp,densam,visc,iquad);
    break;
  case INPAR::MAT::m_biofilm:
    MatBioFilm(material,k,densn,densnp,densam,visc,iquad);
    break;
  case INPAR::MAT::m_scatra_growth_scd:
    MatGrowthScd(material,k,densn,densnp,densam,visc,iquad);
    break;
  default:
    dserror("Material type %i is not supported",material->MaterialType());
    break;
  }
  return;
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
    double&                                 visc,         //!< fluid viscosity
    const int                               iquad         //!< id of current gauss point
  )
{
  const Teuchos::RCP<const MAT::Biofilm>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Biofilm>(material);

  // get constant diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->Diffusivity(),k);

  // get substrate concentration at n+1 or n+alpha_F at integration point
  const double csnp = my::funct_.Dot(my::ephinp_[k]);

  // set reaction coefficient
  my::reamanager_->SetReaCoeff(actmat->ComputeReactionCoeff(csnp),k);
  // set derivative of reaction coefficient
  my::reamanager_->SetReaCoeffDerivMatrix(actmat->ComputeReactionCoeffDeriv(csnp),k,k);

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
    double&                                 visc,         //!< fluid viscosity
    const int                               iquad         //!< id of current gauss point
  )
{
  dsassert(my::numdofpernode_==1,"more than 1 dof per node for ScatraGrowthScd material");

  if(iquad < 0) dserror("ScatraGrowthScd material has to be evaluated at gauss point!");

  const Teuchos::RCP<const MAT::ScatraGrowthScd>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraGrowthScd>(material);

  // get and save constant diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->Diffusivity(),k);

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
  if (structmat == Teuchos::null)
    dserror("dynamic cast of structure material GrowthScd failed.");
  if(structmat->MaterialType() != INPAR::MAT::m_growth_volumetric_scd)
    dserror("invalid structure material for scalar dependent growth");

  if (structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_linear or
      structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_exponential)
  {
    const double theta    = structmat->Gettheta_atgp(iquad);
    const double dtheta   = structmat->Getdtheta_atgp(iquad);
    const double thetaold = structmat->Getthetaold_atgp(iquad);
    const double detFe    = structmat->GetdetFe_atgp(iquad);

    // get substrate concentration at n+1 or n+alpha_F at integration point
    const double csnp = my::funct_.Dot(my::ephinp_[k]);

    // set reaction coefficient
    my::reamanager_->SetReaCoeff(actmat->ComputeReactionCoeff(csnp,theta,dtheta,detFe),k);
    // set derivative of reaction coefficient
    my::reamanager_->SetReaCoeffDerivMatrix(actmat->ComputeReactionCoeffDeriv(csnp,theta,thetaold,1.0),k,k);

    // set density at various time steps and density gradient factor to 1.0/0.0
    densn      = 1.0;
    densnp     = 1.0;
    densam     = 1.0;
  }
  else if (structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_ac or
           structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_ac_radial )
    {
    dserror("In the case of MAT_GrowthAC or MAT_GrowthACNormal one should not end up in here, "
        "since the growth does only change the scalars field size/volume. And this is already"
        " cared due to the conservative formulation you hopefully use!");
    }
  else
  {
    dserror("Your growth law is not a valid one!");
  }

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
                                       //... + reaction terms not depending on phi(k) -> source term
  rhsint = my::bodyforce_[k].Dot(my::funct_) + densnp*ReaManager()->GetReaBodyForce(k);

  return;
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

    for (int condnum = 0;condnum < numcond_;condnum++)
    {

      const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
      const double& reaccoeff = reaccoeff_[condnum]; //get reaction coefficient
      const MAT::PAR::reaction_coupling& couplingtype = couplingtype_[condnum]; //get coupling type
      const double& reacstart = reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] < 0)
    {
      double rcfac= CalcReaCoeffFac(stoich,couplingtype,k);

      if (reacstart>0 and rcfac>0)
        ReacStartForReaCoeff(k,condnum,reacstart,rcfac);

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
          const MAT::PAR::reaction_coupling         couplingtype,            //!<type of coupling the stoichiometry coefficients
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
        case MAT::PAR::reac_coup_simple_multiplicative:
          rcfac *=my::funct_.Dot(my::ephinp_[ii]);
          break;
        //case ... :  //insert new Couplings here
        case MAT::PAR::reac_coup_none:
          dserror("reac_coup_none is not a valid coupling");
          break;
        default:
          dserror("The couplingtype %i is not a valid coupling type.", couplingtype);
          break;
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

    for (int condnum = 0;condnum < numcond_;condnum++)
    {
      const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
      const double& reaccoeff = reaccoeff_[condnum]; //get reaction coefficient
      const MAT::PAR::reaction_coupling& couplingtype = couplingtype_[condnum]; //get coupling type
      const double& reacstart = reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] < 0)
    {
      double rcdmfac = CalcReaCoeffDerivFac(stoich,couplingtype,j,k);

      if (reacstart>0)
        ReacStartForReaCoeffDeriv(k,j,condnum,reacstart,rcdmfac,stoich,couplingtype);

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
          const MAT::PAR::reaction_coupling       couplingtype,            //!<type of coupling the stoichiometry coefficients
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
        case MAT::PAR::reac_coup_simple_multiplicative:
          if (ii!=k and ii!= toderive)
            rcdmfac *= my::funct_.Dot(my::ephinp_[ii]);
          break;
        //case ... :  //insert new Couplings here
        case MAT::PAR::reac_coup_none:
          dserror("reac_coup_none is not a valid coupling");
          break;
        default:
          dserror("The couplingtype %i is not a valid coupling type.", couplingtype);
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
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaBodyForceTerm(
    const int                                 k                       //!< id of current scalar
)
{
  double bodyforcetermK=0;

  for (int condnum = 0;condnum < numcond_;condnum++)
  {

    const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
    const double& reaccoeff = reaccoeff_[condnum]; //get reaction coefficient
    const MAT::PAR::reaction_coupling& couplingtype = couplingtype_[condnum]; //get coupling type
    const double& reacstart = reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] > 0)
    {
      double bftfac = CalcReaBodyForceTermFac(stoich,couplingtype);// scalar at integration point np

      if (reacstart>0 and bftfac>0)
        ReacStartForReaBF(k,condnum,reacstart,bftfac);

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
          const MAT::PAR::reaction_coupling           couplingtype            //!<type of coupling the stoichiometry coefficients
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
      case MAT::PAR::reac_coup_simple_multiplicative:
        bftfac *=my::funct_.Dot(my::ephinp_[ii]);
        break;
      //case ... :  //insert new Couplings here
      case MAT::PAR::reac_coup_none:
        dserror("reac_coup_none is not a valid coupling");
        break;
      default:
        dserror("The couplingtype %i is not a valid coupling type.", couplingtype);
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
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcReaBodyForceDerivMatrix(

    const int                 k,                  //!< id of current scalar
    const int                 j                   //!< concentration to be derived to
)
{
  double reabodyforcederivmatrixKJ=0;

    for (int condnum = 0;condnum < numcond_;condnum++)
    {
      const std::vector<int>& stoich = stoich_[condnum]; //get stoichometrie
      const double& reaccoeff = reaccoeff_[condnum]; //get reaction coefficient
      const MAT::PAR::reaction_coupling& couplingtype = couplingtype_[condnum]; //get coupling type
      const double& reacstart = reacstart_[condnum]; //get reactionstart coefficient

    if (stoich[k] > 0)
    {
      double bfdmfac = CalcReaBodyForceDerivFac(stoich,couplingtype,j);

      if (reacstart>0 and bfdmfac>0)
        ReacStartForReaBFDeriv(k,j,condnum,reacstart,bfdmfac,stoich,couplingtype);

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
        const MAT::PAR::reaction_coupling         couplingtype,            //!<type of coupling the stoichiometry coefficients
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
        case MAT::PAR::reac_coup_simple_multiplicative:
          if (ii!=toderive)
            bfdmfac *= my::funct_.Dot(my::ephinp_[ii]);
          break;
        //case ... :  //insert new Couplings here
        case MAT::PAR::reac_coup_none:
          dserror("reac_coup_none is not a valid coupling");
          break;
        default:
          dserror("The couplingtype %i is not a valid coupling type.", couplingtype);
          break;
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
  const double                            reacstart,      //!< value for reaction starting
  double&                                 value         //!< current reaction value
  )
{
  double prod = value * my::funct_.Dot(my::ephinp_[k]); //for simple multiplikative only!

    if (prod > reacstart )
      value = value - reacstart/my::funct_.Dot(my::ephinp_[k]);
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
        const double                            reacstart,      //!< value for reaction starting
        double&                                 value,          //!< current reaction value
        const std::vector<int>                  stoich,         //!<stoichometrie of current condition
        const MAT::PAR::reaction_coupling       couplingtype    //!<type of coupling the stoichiometry coefficients
  )
{
  double prod = CalcReaBodyForceTermFac(stoich,couplingtype);

  if ( prod > reacstart)
  {
    if (k==toderive)
      value= value - (-reacstart / pow(my::funct_.Dot(my::ephinp_[k]),2) );
  }
  else
    value = 0;
  return;
}

/*-------------------------------------------------------------------------------*
 |  calculate reactions body force term for "delayed" reactions       thon 03/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ReacStartForReaBF(
  const int                               k,            //!< id of current scalar
  const int                               condnum,      //!< id of current condition
  const double                            reacstart,      //!< value for reaction starting
  double&                                 value         //!< current reaction value
  )
{
  if (value > reacstart )
    value = value - reacstart;
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
        const double                            reacstart,      //!< value for reaction starting
        double&                                 value,          //!< current reaction value
        const std::vector<int>                  stoich,         //!<stoichometrie of current condition
        const MAT::PAR::reaction_coupling       couplingtype    //!<type of coupling the stoichiometry coefficients
  )
{
  double prod = CalcReaBodyForceTermFac(stoich,couplingtype);

  if ( prod > reacstart) { }
  else
    value = 0;

  return;
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
  const LINALG::Matrix<my::nen_,1>&      sgconv,
  const LINALG::Matrix<my::nen_,1>&      diff
  )
{
  // -----------------first care for Term K(c)*(\partial_c c)=Id*K(c)--------------------------------------

  my::CalcMatReact(emat,k,timefacfac,timetaufac,taufac,densnp,sgconv,diff);

  const double&                       phinp = my::scatravarmanager_->Phinp(k);
  const LINALG::Matrix<my::nen_,1>&   conv  = my::scatravarmanager_->Conv();

  // -----------------second care for Term (\partial_c K(c)) .* c - (\partial_c f_{reabody}(c))------------

  LINALG::Matrix<my::nen_,1> functint = my::funct_;
  if (not my::scatrapara_->MatGP())
    functint = funct_elementcenter_;

  for (int j=0; j<my::numscal_ ;j++)
  {
    const double fac_reac        = timefacfac*densnp*( my::reamanager_->GetReaCoeffDerivMatrix(k,j)*phinp - ReaManager()->GetReaBodyForceDerivMatrix(k,j) );
    const double timetaufac_reac = timetaufac*densnp*( my::reamanager_->GetReaCoeffDerivMatrix(k,j)*phinp - ReaManager()->GetReaBodyForceDerivMatrix(k,j) );

    //----------------------------------------------------------------
    // standard Galerkin reactive term
    //----------------------------------------------------------------
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac_reac*functint(vi);
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
        const double v = densreataufac*(conv(vi)+sgconv(vi)+my::scatrapara_->USFEMGLSFac()*1.0/my::scatraparatimint_->TimeFac()*functint(vi));
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
        const double v = densreataufac*functint(vi);
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+j;

          emat(fvi,fui) += v*(conv(ui)+my::reamanager_->GetReaCoeff(k)*my::funct_(ui));
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
          const double v = my::scatrapara_->USFEMGLSFac()*taufac*densnp*my::reamanager_->GetReaCoeff(k)*densnp*functint(vi);
          const int fvi = vi*my::numdofpernode_+k;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+j;

            emat(fvi,fui) += v*my::funct_(ui);
          }
        }

        if (my::use2ndderiv_ and my::reamanager_->GetReaCoeff(k)!=0.0)
          dserror("Second order reactive stabilization is not fully implemented!! ");
      }
    }
  } //end for
  return;
}

/*-----------------------------------------------------------------------------------------*
 |  get numcond, stoich list, reaction coefficient, couplingtpye from material  thon 09/14 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::GetAdvancedReactionCoefficients(
    const Teuchos::RCP<const MAT::Material> material //!< pointer to current material
  )
{
  const Teuchos::RCP<const MAT::MatListReactions>& actmat = Teuchos::rcp_dynamic_cast<const MAT::MatListReactions>(material);

  if (numcond_ == -1) //if not yet initialisied
  {
    if(actmat==Teuchos::null)
      dserror("cast to MatListReactions failed");

    numcond_= actmat->NumReac();

    stoich_.resize(numcond_);
    couplingtype_.resize(numcond_);
    reaccoeff_.resize(numcond_);
    reacstart_.resize(numcond_);
  }

  for (int i=0;i<numcond_;i++)
  {
    const int reacid = actmat->ReacID(i);
    const Teuchos::RCP<const MAT::ScatraReactionMat>& reacmat = Teuchos::rcp_dynamic_cast<const MAT::ScatraReactionMat>(actmat->MaterialById(reacid));

    stoich_[i] = *(reacmat->Stoich()); //get stoichometrie
    couplingtype_[i] = reacmat->Coupling(); //get coupling type
    reaccoeff_[i] = reacmat->ReacCoeff(); //get reaction coefficient
    reacstart_[i] = reacmat->ReacStart(); //get reaction start coefficient
  }
}

/*-------------------------------------------------------------------------------*
 |  set body force, reaction coefficient and derivatives              thon 09/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::SetAdvancedReactionTerms(
    const int                               k,
    const double                            scale
                                    )
{
  ReaManager()->SetReaBodyForce( CalcReaBodyForceTerm(k)*scale ,k);

  ReaManager()->SetReaCoeff( CalcReaCoeff(k)*scale ,k);

  for (int j=0; j<my::numscal_ ;j++)
  {
    ReaManager()->SetReaBodyForceDerivMatrix( CalcReaBodyForceDerivMatrix(k,j)*scale ,k,j );

    my::reamanager_->SetReaCoeffDerivMatrix( CalcReaCoeffDerivMatrix(k,j)*scale ,k,j );
  }

}

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at ele. center   jhoer 11/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::EvalShapeFuncAndDerivsAtEleCenter()
{
  const double vol = my::EvalShapeFuncAndDerivsAtEleCenter();

  //shape function at element center
  funct_elementcenter_ = my::funct_;

  return vol;

} //ScaTraImpl::EvalShapeFuncAndDerivsAtEleCenter


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
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs27>;
