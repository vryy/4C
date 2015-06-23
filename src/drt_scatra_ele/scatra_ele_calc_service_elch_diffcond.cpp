/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_diffcond.cpp

\brief evaluation of scatra elements for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"

#include "scatra_ele_calc_elch_diffcond.H"


/*-----------------------------------------------------------------------*
  |  Set scatra element parameter                             ehrl 01/14 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CheckElchElementParameter(
    DRT::Element*              ele
  )
{
  // get the material
  Teuchos::RCP<const MAT::Material> material = ele->Material();

  // 1) Check material specific options
  // 2) Check if numdofpernode, numscal is set correctly
  if (material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat>& actmat
          = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

    int numphase = actmat->NumPhase();

    // access mat_elchmat: container material for porous structures in elch
    if (numphase != 1) dserror("In the moment a single phase is only allowed.");

    // 1) loop over single phases
    for (int iphase=0; iphase < actmat->NumPhase();++iphase)
    {
      // access phase material
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlephase = actmat->PhaseById(phaseid);

      // dynmic cast: get access to mat_phase
      const Teuchos::RCP<const MAT::ElchPhase>& actphase
                = Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(singlephase);

      // Check if numdofpernode, numscal is set correctly
      int nummat = actphase->NumMat();
      // enough materials defined
      if (nummat != my::numscal_)
        dserror("The number of scalars defined in the material ElchMat does not correspond with "
                "the number of materials defined in the material MatPhase.");

      int numdofpernode = 0;
      if (ElchPara()->CurSolVar())
        numdofpernode = nummat+DRT::Problem::Instance()->NDim()+numphase;
      else
        numdofpernode = nummat+numphase;

      if(numdofpernode != my::numdofpernode_)
        dserror("The chosen element formulation (e.g. current as solution variable) "
                "does not correspond with the number of dof's defined in your material");

      // 2) loop over materials of the single phase
      for (int imat=0; imat < actphase->NumMat();++imat)
      {
        const int matid = actphase->MatID(imat);
        Teuchos::RCP<const MAT::Material> singlemat = actphase->MatById(matid);

        if(singlemat->MaterialType() == INPAR::MAT::m_newman)
        {
          // Material Newman is derived for a binary electrolyte utilizing the ENC to condense the non-reacting species
          if(my::numscal_>1)
            dserror("Material Newman is only valid for one scalar (binary electrolyte utilizing the ENC)");
        }
      }
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------------*
 | calculate element mass matrix and element residual for initial time derivative   fang 03/15 |
 *---------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcInitialTimeDerivative(
    DRT::Element*                 ele,              //!< current element
    Epetra_SerialDenseMatrix&     emat,             //!< element matrix
    Epetra_SerialDenseVector&     erhs,             //!< element residual
    Teuchos::ParameterList&       params,           //!< parameter list
    DRT::Discretization&          discretization,   //!< discretization
    DRT::Element::LocationArray&  la                //!< location array
    )
{
  // call base class routine
  myelch::CalcInitialTimeDerivative(
      ele,
      emat,
      erhs,
      params,
      discretization,
      la
      );

  // dummy mass matrix for the electric current dofs
  if(ElchPara()->CurSolVar())
  {
    // integration points and weights
    const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      for(int idim=0;idim<my::nsd_;++idim)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = fac*my::funct_(vi); // no density required here
          const int fvi = vi*my::numdofpernode_+my::numscal_+1+idim;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+my::numscal_+1+idim;

            emat(fvi,fui) += v*my::funct_(ui);
          }
        }
      }
    }
  }

  // In the moment the diffusion manager contains the porosity at the last Gauss point (previous call my::CalcInitialTimeDerivative())
  // Since the whole approach is valid only for constant porosities, we do not fill the diffusion manager again at the element center
  // The solution variable is the initial time derivative. Therefore, we have to correct emat by the initial porosity
  // Attention: this procedure is only valid for a constant porosity in the beginning
  emat.Scale(DiffManager()->GetPhasePoro(0));

  return;
}


/*----------------------------------------------------------------------*
 |  CorrectRHSFromCalcRHSLinMass                             ehrl 06/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CorrectRHSFromCalcRHSLinMass(
    Epetra_SerialDenseVector&     erhs,
    const int                     k,
    const double                  fac,
    const double                  densnp,
    const double                  phinp
  )
{
  // fac->-fac to change sign of rhs
  if (my::scatraparatimint_->IsIncremental())
     my::CalcRHSLinMass(erhs,k,0.0,-fac,0.0,DiffManager()->GetPhasePoro(0));
  else
    dserror("Must be incremental!");

  return;
}


/*-----------------------------------------------------------------------------------------*
 | extract element based or nodal values and return extracted values of phinp   fang 01/15 |
 *-----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const std::vector<double> DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
    )
{
  // call base class routine
  const std::vector<double> myphinp = myelch::ExtractElementAndNodeValues(ele,params,discretization,la);

  // get current density at element nodes if required
  if(ElchPara()->CurSolVar())
  {
    // get values for current at element nodes
    for (int ien=0;ien<my::nen_;++ien)
    {
      for(int idim=0; idim<my::nsd_; ++idim)
      {
        //current is stored after potential
        ecurnp_(idim,ien) = myphinp[ien*my::numdofpernode_+(my::numscal_+1)+idim];
      }
    }
  }
  else
    ecurnp_.Clear();

  return myphinp;
}


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     ae 05/15|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalculateFlux(
    LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac         //!< integration factor
)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
  case INPAR::SCATRA::flux_total_domain:
    // call base class routine
    myelectrode::CalculateFlux(q,fluxtype,k,fac);

    // no break statement here!
  case INPAR::SCATRA::flux_diffusive_domain:
    // diffusive flux contribution
//<<<<<<< .mine
//    q.Update(-dmedc->GetPhasePoroTort(0)*dmedc->GetIsotropicDiff(k),vmdc->GradPhi(k),1.0);
//=======
    q.Update(-DiffManager()->GetIsotropicDiff(k)*DiffManager()->GetPhasePoroTort(0),VarManager()->GradPhi(k),1.0);
    // flux due to ohmic overpotential
//<<<<<<< .mine
//    q.Update(-dmedc->GetTransNum(k)*vmdc->InvFVal(k)*dmedc->GetPhasePoroTort(0)*dmedc->GetCond(),vmdc->GradPot(),1.0);
//=======
    q.Update(-DiffManager()->GetTransNum(k)*VarManager()->InvFVal(k)*DiffManager()->GetCond()*DiffManager()->GetPhasePoroTort(0),VarManager()->GradPot(),1.0);
    // flux due to concentration overpotential
//<<<<<<< .mine
//    q.Update(-dmedc->GetTransNum(k)*vmdc->RTFFCVal(k)*dmedc->GetPhasePoroTort(0)*dmedc->GetCond()*dmedc->GetThermFac()*(myelch::elchpara_->NewmanConstA()+(myelch::elchpara_->NewmanConstB()*dmedc->GetTransNum(k)))*vmdc->ConIntInv(k),vmdc->GradPhi(k),1.0);
//=======
    q.Update(-DiffManager()->GetTransNum(k)*VarManager()->RTFFCVal(k)*DiffManager()->GetCond()*DiffManager()->GetPhasePoroTort(0)*DiffManager()->GetThermFac()*(ElchPara()->NewmanConstA()+(ElchPara()->NewmanConstB()*DiffManager()->GetTransNum(k)))*VarManager()->ConIntInv(k),VarManager()->GradPhi(k),1.0);
    break;
  default:
    dserror("received illegal flag inside flux evaluation for whole domain"); break;
  };

  return;
} // ScaTraCalc::CalculateFlux


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     ae 05/15|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalculateCurrent(
    LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const double                    fac         //!< integration factor
)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
  case INPAR::SCATRA::flux_total_domain:
  case INPAR::SCATRA::flux_diffusive_domain:
    // ohmic flux contribution
    q.Update(-DiffManager()->GetCond(),VarManager()->GradPot(),1.0);
    // diffusion overpotential flux contribution
    for (int k = 0; k<my::numscal_; ++k)
      q.Update(-VarManager()->RTF()/ElchPara()->NewmanConstC()*DiffManager()->GetCond()*DiffManager()->GetThermFac()*(ElchPara()->NewmanConstA()+(ElchPara()->NewmanConstB()*DiffManager()->GetTransNum(k)))*VarManager()->ConIntInv(k),VarManager()->GradPhi(k),1.0);

    break;
  default:
    dserror("received illegal flag inside flux evaluation for whole domain"); break;
  };


  return;
} // ScaTraCalc::CalculateCurrent

/*----------------------------------------------------------------------*
 | get conductivity                                          fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::GetConductivity(
    const enum INPAR::ELCH::EquPot   equpot,      //!< type of closing equation for electric potential
    double&                          sigma_all,   //!< conductivity of electrolyte solution
    std::vector<double>&             sigma,        //!< conductivity or a single ion + overall electrolyte solution
    bool                             effCond
    )
{
  // use precomputed conductivity
  sigma_all = DiffManager()->GetCond();

  if(effCond == true)
  {
    sigma_all = sigma_all * DiffManager()->GetPhasePoroTort(0);
  }

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::GetConductivity

/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution           gjb 10/08|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalErrorComparedToAnalytSolution(
  const DRT::Element*                   ele,
  Teuchos::ParameterList&               params,
  Epetra_SerialDenseVector&             errors
  )
{
  //at the moment, there is only one analytical test problem available!
  if (DRT::INPUT::get<SCATRA::Action>(params,"action") != SCATRA::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // in the ALE case add nodal displacements
  if (my::scatrapara_->IsAle()) dserror("No ALE for Kwok & Wu error calculation allowed.");

  // set constants for analytical solution
  const double t = my::scatraparatimint_->Time() + (1- my::scatraparatimint_->AlphaF())* my::scatraparatimint_->Dt(); //-(1-alphaF_)*dta_
  const double frt = ElchPara()->FRT();

  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // get material parameter (constants values)
  this->GetMaterialParams(ele,densn,densnp,densam,visc);

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcError errortype = DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag");
  switch(errortype)
  {
  case INPAR::SCATRA::calcerror_Kwok_Wu:
  {
    //   References:
    //   Kwok, Yue-Kuen and Wu, Charles C. K.
    //   "Fractional step algorithm for solving a multi-dimensional diffusion-migration equation"
    //   Numerical Methods for Partial Differential Equations
    //   1995, Vol 11, 389-397

    //   G. Bauer, V. Gravemeier, W.A. Wall,
    //   A 3D finite element approach for the coupled numerical simulation of
    //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339–1359.

    if (my::numscal_ != 1)
      dserror("Numscal_ != 1 for desired error calculation.");

    // working arrays
    double                      potint(0.0);
    LINALG::Matrix<1,1>         conint(true);
    LINALG::Matrix<my::nsd_,1>  xint(true);
    LINALG::Matrix<1,1>         c(true);
    double                      deltapot(0.0);
    LINALG::Matrix<1,1>         deltacon(true);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // get values of all transported scalars at integration point
      conint(0) = my::funct_.Dot(my::ephinp_[0]);

      // get el. potential solution at integration point
      potint = my::funct_.Dot(myelch::epotnp_);

      // get global coordinate of integration point
      xint.Multiply(my::xyze_,my::funct_);

      // compute various constants
//      const double d = frt*((DiffManager()->GetIsotropicDiff(0)*DiffManager()->GetValence(0)) - (DiffManager()->GetIsotropicDiff(1)*DiffManager()->GetValence(1)));
//      if (abs(d) == 0.0) dserror("division by zero");
      const double D = DiffManager()->GetIsotropicDiff(0);

      // compute analytical solution for cation and anion concentrations
      const double A0 = 2.0;
      const double m = 1.0;
      const double n = 2.0;
      const double k = 3.0;
      const double A_mnk = 1.0;
      double expterm;
      double c_0_0_0_t;

      if (my::nsd_==3)
      {
        expterm = exp((-D)*(m*m + n*n + k*k)*t*PI*PI);
        c(0) = A0 + (A_mnk*(cos(m*PI*xint(0))*cos(n*PI*xint(1))*cos(k*PI*xint(2)))*expterm);
        c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m + n*n + k*k)*t*PI*PI));
      }
      else if (my::nsd_==2)
      {
        expterm = exp((-D)*(m*m + n*n)*t*PI*PI);
        c(0) = A0 + (A_mnk*(cos(m*PI*xint(0))*cos(n*PI*xint(1)))*expterm);
        c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m + n*n)*t*PI*PI));
      }
      else if (my::nsd_==1)
      {
        expterm = exp((-D)*(m*m)*t*PI*PI);
        c(0) = A0 + (A_mnk*(cos(m*PI*xint(0)))*expterm);
        c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m)*t*PI*PI));
      }
      else
        dserror("Illegal number of space dimensions for analyt. solution: %d",my::nsd_);

      // compute analytical solution for el. potential
      //const double pot = ((DiffManager()->GetIsotropicDiff(1)-DiffManager()->GetIsotropicDiff(0))/d) * log(c(0)/c_0_0_0_t);
      const double pot = -1/frt * (ElchPara()->NewmanConstA() + (ElchPara()->NewmanConstB()*DiffManager()->GetTransNum(0))) / ElchPara()->NewmanConstC() * log(c(0)/c_0_0_0_t);

      // compute differences between analytical solution and numerical solution
      deltapot = potint - pot;
      deltacon.Update(1.0,conint,-1.0,c);

      // add square to L2 error
      errors[0] += deltacon(0)*deltacon(0)*fac; // cation concentration
      errors[1] += 0;                          // anion concentration
      errors[2] += deltapot*deltapot*fac; // electric potential in electrolyte solution

    } // end of loop over integration points
  } // Kwok and Wu
  break;
  default: dserror("Unknown analytical solution!"); break;
  } //switch(errortype)

  return;
} // CalErrorComparedToAnalytSolution


/*------------------------------------------------------------------------------*
 | set internal variables for diffusion-conduction formulation       fang 02/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchDiffCond(my::funct_,my::derxy_,my::ephinp_,my::ephin_,myelch::epotnp_,my::econvelnp_,my::ehist_,DiffManager(),ecurnp_);

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs27>;
