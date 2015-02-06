/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch.cpp

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
#include "scatra_ele_action.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"  // for time curve in body force
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_globalproblem.H"  // consistency check of formulation and material

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"

#include "scatra_ele_calc_elch.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateAction(
    DRT::ELEMENTS::Transport*   ele,
    Teuchos::ParameterList&     params,
    DRT::Discretization&        discretization,
    const SCATRA::Action&       action,
    const std::vector<int> &    lm,
    Epetra_SerialDenseMatrix&   elemat1_epetra,
    Epetra_SerialDenseMatrix&   elemat2_epetra,
    Epetra_SerialDenseVector&   elevec1_epetra,
    Epetra_SerialDenseVector&   elevec2_epetra,
    Epetra_SerialDenseVector&   elevec3_epetra
    )
{
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::check_scatra_element_parameter:
  {
    CheckElchElementParameter(ele);
    break;
  }
  case SCATRA::calc_initial_time_deriv:
  {
    // calculate matrix and rhs
    my::CalcInitialTimeDerivative(
      ele,
      elemat1_epetra,
      elevec1_epetra,
      params,
      discretization,
      lm
      );

    // do another loop over the integration points for the potential:
    // At this point time is not so critical since the CalcInitialTimeDerivative()
    // is called once in the beginning of the simulation!

    // we put a dummy mass matrix here in order to have a regular
    // matrix in the lower right block of the whole system-matrix
    // A identity matrix would cause problems with ML solver in the SIMPLE
    // schemes since ML needs to have off-diagonal entries for the aggregation!

    PrepMatAndRhsInitialTimeDerivative(elemat1_epetra,elevec1_epetra);

    break;
  }
  case SCATRA::integrate_shape_functions:
  {
    // calculate integral of shape functions
    const Epetra_IntSerialDenseVector dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
    my::IntegrateShapeFunctions(ele,elevec1_epetra,dofids);

    break;
  }
  case SCATRA::calc_flux_domain:
  {
    //--------------------------------------------------------------------------------
    // extract element based or nodal values
    //--------------------------------------------------------------------------------

    // get velocity values at the nodes
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,my::evelnp_,velocity,my::nsd_);
    const Teuchos::RCP<Epetra_MultiVector> convelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,my::econvelnp_,convelocity,my::nsd_);

    // need current values of transported scalar
    // -> extract local values from global vectors
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // fill all element arrays
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0;k<my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        my::ephinp_[k](i,0) = myphinp[k+(i*my::numdofpernode_)];
      }
      epotnp_(i) = myphinp[i*my::numdofpernode_+my::numscal_];
    } // for i

    //----------------------------------------------------------------------
    // calculation of element volume both for tau at ele. cent. and int. pt.
    //----------------------------------------------------------------------
    my::EvalShapeFuncAndDerivsAtEleCenter();

    //----------------------------------------------------------------------
    // get material and stabilization parameters (evaluation at element center)
    //----------------------------------------------------------------------
    // density at t_(n)
    double densn(1.0);
    // density at t_(n+1) or t_(n+alpha_F)
    double densnp(1.0);
    // density at t_(n+alpha_M)
    double densam(1.0);

    // fluid viscosity
    double visc(0.0);

    // material parameter at the element center
    if (not my::scatrapara_->MatGP())
      this->GetMaterialParams(ele,densn,densnp,densam,visc);

    //----------------------------------------------------------------------
    // integration loop for one element
    //----------------------------------------------------------------------
    // integration points and weights
    const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      //----------------------------------------------------------------------
      // get material parameters (evaluation at integration point)
      //----------------------------------------------------------------------
      if (my::scatrapara_->MatGP())
        this->GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

      // set internal variables
      SetInternalVariablesForMatAndRHS();

      // access control parameter for flux calculation
      INPAR::SCATRA::FluxType fluxtype = ElchPara()->WriteFlux();
      Teuchos::RCP<std::vector<int> > writefluxids = ElchPara()->WriteFluxIds();

      // do a loop for systems of transported scalars
      for (std::vector<int>::iterator it = writefluxids->begin(); it!=writefluxids->end(); ++it)
      {
        int k=0;
        // Actually, we compute here a weighted (and integrated) form of the fluxes!
        // On time integration level, these contributions are then used to calculate
        // an L2-projected representation of fluxes.
        // Thus, this method here DOES NOT YET provide flux values that are ready to use!!

        // allocate and initialize!
        LINALG::Matrix<my::nsd_,1> q(true);

        if((*it) != my::numdofpernode_)
        {
          k=(*it)-1;
          CalculateFlux(q,fluxtype,k,fac);
        }
        else if ((*it) == my::numdofpernode_)
        {
          k=my::numdofpernode_-1;
          CalculateCurrent(q,fluxtype,fac);
        }
        else
          dserror("Flux id, which should be calculated, does not exit in the dof set.");

        // integrate and assemble everything into the "flux" vector
        for (int vi=0; vi < my::nen_; vi++)
        {
          const int fvi = vi*my::numdofpernode_+k;
          elevec1_epetra[fvi] += fac*my::funct_(vi)*q(0);
          elevec2_epetra[fvi] += fac*my::funct_(vi)*q(1);
          if(my::nsd_<3)
            elevec3_epetra[fvi] = 0.0;
          else
            elevec3_epetra[fvi] += fac*my::funct_(vi)*q(2);
        } // vi
      }
    }

    break;
  }
  case SCATRA::calc_error:
  {
    // check if length suffices
    if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

    // need current solution
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // fill element arrays
    for (int i=0;i<my::nen_;++i)
    {
      // split for each transported scalar, insert into element arrays
      for (int k = 0; k< my::numscal_; ++k)
      {
        my::ephinp_[k](i) = myphinp[k+(i*my::numdofpernode_)];
      }
      // get values for el. potential at element nodes
      epotnp_(i) = myphinp[i*my::numdofpernode_+my::numscal_];
    } // for i

    CalErrorComparedToAnalytSolution(
      ele,
      params,
      elevec1_epetra);

    break;
  }
  case SCATRA::calc_elch_conductivity:
  {
    // extract local values from the global vector
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // fill element arrays
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0; k< my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        my::ephinp_[k](i,0) = myphinp[k+(i*my::numdofpernode_)];
      }
    } // for i

    // global element of processor 0 is printed to the screen
    if(discretization.Comm().MyPID()==0)
      std::cout << "Electrolyte conductivity evaluated at global element " << ele->Id() << ":" << std::endl;

    CalculateConductivity(ele,ElchPara()->EquPot(),elevec1_epetra);
    break;
  }
  case SCATRA::calc_elch_initial_potential:
  {
    // need initial field -> extract local values from the global vector
    Teuchos::RCP<const Epetra_Vector> phi0 = discretization.GetState("phi0");
    if (phi0==Teuchos::null) dserror("Cannot get state vector 'phi0'");
    std::vector<double> myphi0(lm.size());
    DRT::UTILS::ExtractMyValues(*phi0,myphi0,lm);

    // fill element arrays
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0; k< my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        my::ephinp_[k](i,0) = myphi0[k+(i*my::numdofpernode_)];
      }
      // get values for el. potential at element nodes
      epotnp_(i) = myphi0[i*my::numdofpernode_+my::numscal_];
    } // for i

    CalculateElectricPotentialField(ele,ElchPara()->EquPot(),elemat1_epetra,elevec1_epetra);

    break;
  }
  default:
  {
    my::EvaluateAction(ele,
                       params,
                       discretization,
                       action,
                       lm,
                       elemat1_epetra,
                       elemat2_epetra,
                       elevec1_epetra,
                       elevec2_epetra,
                       elevec3_epetra);
    break;
  }
  } // switch(action)

  return 0;
}


/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution           gjb 10/08|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalErrorComparedToAnalytSolution(
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
  const double t = my:: scatraparatimint_->Time();
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

//  if(diffcond_==true)
//  {
//    dserror("Analytical solution for Kwok and Wu is only valid for dilute electrolyte solutions!!\n"
//            "Compute corresponding transport properties on your on and activate it here");
//
//    diffus_[0] = 2.0e-3;
//    diffus_[1] = 4.0e-3;
//    valence_[0] = 1.0;
//    valence_[1] = -2.0;
//  }

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

    //if (numscal_ != 2)
    //  dserror("Numscal_ != 2 for desired error calculation.");

    // working arrays
    double                      potint(0.0);
    LINALG::Matrix<2,1>         conint(true);
    LINALG::Matrix<my::nsd_,1>  xint(true);
    LINALG::Matrix<2,1>         c(true);
    double                      deltapot(0.0);
    LINALG::Matrix<2,1>         deltacon(true);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // get values of all transported scalars at integration point
      for (int k=0; k<my::numscal_; ++k)
      {
        conint(k) = my::funct_.Dot(my::ephinp_[k]);
      }

      // get el. potential solution at integration point
      potint = my::funct_.Dot(epotnp_);

      // get global coordinate of integration point
      xint.Multiply(my::xyze_,my::funct_);

      // compute various constants
      const double d = frt*((DiffManager()->GetIsotropicDiff(0)*DiffManager()->GetValence(0)) - (DiffManager()->GetIsotropicDiff(1)*DiffManager()->GetValence(1)));
      if (abs(d) == 0.0) dserror("division by zero");
      const double D = frt*((DiffManager()->GetValence(0)*DiffManager()->GetIsotropicDiff(0)*DiffManager()->GetIsotropicDiff(1)) - (DiffManager()->GetValence(1)*DiffManager()->GetIsotropicDiff(1)*DiffManager()->GetIsotropicDiff(0)))/d;

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

      // compute analytical solution for anion concentration
      c(1) = (-DiffManager()->GetValence(0)/DiffManager()->GetValence(1))* c(0);
      // compute analytical solution for el. potential
      const double pot = ((DiffManager()->GetIsotropicDiff(1)-DiffManager()->GetIsotropicDiff(0))/d) * log(c(0)/c_0_0_0_t);

      // compute differences between analytical solution and numerical solution
      deltapot = potint - pot;
      deltacon.Update(1.0,conint,-1.0,c);

      // add square to L2 error
      errors[0] += deltacon(0)*deltacon(0)*fac; // cation concentration
      errors[1] += deltacon(1)*deltacon(1)*fac; // anion concentration
      errors[2] += deltapot*deltapot*fac; // electric potential in electrolyte solution

    } // end of loop over integration points
  } // Kwok and Wu
  break;
  case INPAR::SCATRA::calcerror_cylinder:
  {
    // two-ion system with Butler-Volmer kinetics between two concentric cylinders
    //   G. Bauer, V. Gravemeier, W.A. Wall,
    //   A 3D finite element approach for the coupled numerical simulation of
    //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339–1359.

    if (my::numscal_ != 2)
      dserror("Numscal_ != 2 for desired error calculation.");

    // working arrays
    LINALG::Matrix<2,1>     conint(true);
    LINALG::Matrix<my::nsd_,1>  xint(true);
    LINALG::Matrix<2,1>     c(true);
    LINALG::Matrix<2,1>     deltacon(true);

    // some constants that are needed
    const double c0_inner = 0.6147737641011396;
    const double c0_outer = 1.244249192148809;
    const double r_inner = 1.0;
    const double r_outer = 2.0;
    const double pot_inner = 2.758240847314454;
    const double b = log(r_outer/r_inner);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // get values of all transported scalars at integration point
      for (int k=0; k<my::numscal_; ++k)
      {
        conint(k) = my::funct_.Dot(my::ephinp_[k]);
      }

      // get el. potential solution at integration point
      const double potint = my::funct_.Dot(epotnp_);

      // get global coordinate of integration point
      xint.Multiply(my::xyze_,my::funct_);

      // evaluate analytical solution for cation concentration at radial position r
      if (my::nsd_==3)
      {
        const double r = sqrt(xint(0)*xint(0) + xint(1)*xint(1));
        c(0) = c0_inner + ((c0_outer- c0_inner)*(log(r) - log(r_inner))/b);
      }
      else
        dserror("Illegal number of space dimensions for analyt. solution: %d",my::nsd_);

      // compute analytical solution for anion concentration
      c(1) = (-DiffManager()->GetValence(0)/DiffManager()->GetValence(1))* c(0);
      // compute analytical solution for el. potential
      const double d = frt*((DiffManager()->GetIsotropicDiff(0)*DiffManager()->GetValence(0)) - (DiffManager()->GetIsotropicDiff(1)*DiffManager()->GetValence(1)));
      if (abs(d) == 0.0) dserror("division by zero");
      // reference value + ohmic resistance + concentration potential
      const double pot = pot_inner + log(c(0)/c0_inner); // + (((diffus_[1]-diffus_[0])/d) * log(c(0)/c0_inner));

      // compute differences between analytical solution and numerical solution
      double deltapot = potint - pot;
      deltacon.Update(1.0,conint,-1.0,c);

      // add square to L2 error
      errors[0] += deltacon(0)*deltacon(0)*fac; // cation concentration
      errors[1] += deltacon(1)*deltacon(1)*fac; // anion concentration
      errors[2] += deltapot*deltapot*fac; // electric potential in electrolyte solution

    } // end of loop over integration points
  } // concentric cylinders
  break;
  case INPAR::SCATRA::calcerror_electroneutrality:
  {
    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // get values of transported scalars at integration point
      // and compute electroneutrality
      double deviation(0.0);
      for (int k=0; k<my::numscal_; ++k)
      {
        const double conint_k = my::funct_.Dot(my::ephinp_[k]);
        deviation += DiffManager()->GetValence(k)*conint_k;
      }

    // add square to L2 error
    errors[0] += deviation*deviation*fac;
    } // loop over integration points
  }
  break;
  default: dserror("Unknown analytical solution!"); break;
  } //switch(errortype)

  return;
} // ScaTraImpl::CalErrorComparedToAnalytSolution


/*----------------------------------------------------------------------*
  |  Calculate conductivity (ELCH) (private)                   gjb 07/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalculateConductivity(
  const DRT::Element*               ele,
  const enum INPAR::ELCH::EquPot    equpot,
  Epetra_SerialDenseVector&         sigma
  )
{
  my::EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // material parameter at the element center
  this->GetMaterialParams(ele,densn,densnp,densam,visc);

  // compute the conductivity (1/(\Omega m) = 1 Siemens / m)
  double sigma_all(0.0);

  GetConductivity(equpot,sigma_all, sigma);

  // conductivity based on ALL ionic species (even eliminated ones!)
  sigma[my::numscal_] += sigma_all;

  return;

} //ScaTraEleCalcElch()


/*----------------------------------------------------------------------*
 | CalculateElectricPotentialField (protected)                gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalculateElectricPotentialField(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::EquPot    equpot,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs
  )
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //----------------------------------------------------------------------
    // get material and stabilization parameters (evaluation at gauss point)
    //----------------------------------------------------------------------
    // density at t_(n)
    double densn(1.0);
    // density at t_(n+1) or t_(n+alpha_F)
    double densnp(1.0);
    // density at t_(n+alpha_M)
    double densam(1.0);

    // fluid viscosity
    double visc(0.0);

    // material parameter at the element center
    this->GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

    // set internal variables
    SetInternalVariablesForMatAndRHS();

    CalcMatAndRhsElectricPotentialField(equpot,emat,erhs,fac,1.);
  } // integration loop

  return;

} //ScaTraImpl<distype>::CalculateElectricPotentialField


/*--------------------------------------------------------------------------------------------*
 | compute dummy element matrix and residual entries for electric potential dofs   fang 02/15 |
 *--------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::PrepMatAndRhsInitialTimeDerivative(
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra
)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // element integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac*my::funct_(vi); // no density required here
      const int fvi = vi*my::numdofpernode_+my::numscal_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+my::numscal_;

        elemat1_epetra(fvi,fui) += v*my::funct_(ui);
      }
    }
  }

  // set zero for the rhs entries associated with the electric potential
  for (int vi=0; vi<my::nen_; ++vi)
    elevec1_epetra[vi*my::numdofpernode_+my::numscal_] = 0.;

  return;
}


/*----------------------------------------------------------------------------------------*
 | finite difference check on element level (for debugging only) (protected)   fang 10/14 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::FDCheck(
  DRT::Element*                              ele,
  Epetra_SerialDenseMatrix&                  emat,
  Epetra_SerialDenseVector&                  erhs,
  Epetra_SerialDenseVector&                  subgrdiff
  )
{
  // screen output
  std::cout << "FINITE DIFFERENCE CHECK FOR ELEMENT " << ele->Id();

  // make a copy of state variables to undo perturbations later
  std::vector<LINALG::Matrix<my::nen_,1> > ephinp_original(my::numscal_);
  for(int k=0; k<my::numscal_; ++k)
    for(int i=0; i<my::nen_; ++i)
      ephinp_original[k](i,0) = my::ephinp_[k](i,0);
  LINALG::Matrix<my::nen_,1> epotnp_original(true);
  for(int i=0; i<my::nen_; ++i)
    epotnp_original(i) = epotnp_(i);

  // generalized-alpha time integration requires a copy of history variables as well
  std::vector<LINALG::Matrix<my::nen_,1> > ehist_original(my::numscal_);
  if(my::scatraparatimint_->IsGenAlpha())
  {
    for(int k=0; k<my::numscal_; ++k)
      for(int i=0; i<my::nen_; ++i)
        ehist_original[k](i,0)  = my::ehist_[k](i,0);
  }

  // initialize element matrix and vectors for perturbed state
  Epetra_SerialDenseMatrix emat_dummy(emat);
  Epetra_SerialDenseVector erhs_perturbed(erhs);
  Epetra_SerialDenseVector subgrdiff_dummy(subgrdiff);

  // initialize counter for failed finite difference checks
  unsigned counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  // loop over columns of element matrix by first looping over nodes and then over dofs at each node
  for(int inode=0; inode<my::nen_; ++inode)
  {
    for(int idof=0; idof<my::numdofpernode_; ++idof)
    {
      // number of current column of element matrix
      unsigned col = inode*my::numdofpernode_+idof;

      // clear element matrix and vectors for perturbed state
      emat_dummy.Scale(0.0);
      erhs_perturbed.Scale(0.0);
      subgrdiff_dummy.Scale(0.0);

      // fill state vectors with original state variables
      for(int k=0; k<my::numscal_; ++k)
        for(int i=0; i<my::nen_; ++i)
          my::ephinp_[k](i,0) = ephinp_original[k](i,0);
      for(int i=0; i<my::nen_; ++i)
        epotnp_(i) = epotnp_original(i);
      if(my::scatraparatimint_->IsGenAlpha())
        for(int k=0; k<my::numscal_; ++k)
          for(int i=0; i<my::nen_; ++i)
            my::ehist_[k](i,0)  = ehist_original[k](i,0);

      // impose perturbation
      if(my::scatraparatimint_->IsGenAlpha())
      {
        // perturbation on concentration
        if(idof < my::numscal_)
        {
          // perturbation of phi(n+alphaF), not of phi(n+1) => scale epsilon by factor alphaF
          my::ephinp_[idof](inode,0) += my::scatraparatimint_->AlphaF() * my::scatrapara_->FDCheckEps();

          // perturbation of phi(n+alphaF) by alphaF*epsilon corresponds to perturbation of phidtam (stored in ehist_)
          // by alphaM*epsilon/(gamma*dt); note: alphaF/timefac = alphaM/(gamma*dt)
          my::ehist_[idof](inode,0) += my::scatraparatimint_->AlphaF() / my::scatraparatimint_->TimeFac() * my::scatrapara_->FDCheckEps();
        }

        // perturbation on electric potential
        else
          epotnp_(inode) += my::scatraparatimint_->AlphaF() * my::scatrapara_->FDCheckEps();
      }
      else
      {
        // perturbation on concentration
        if(idof < my::numscal_)
          my::ephinp_[idof](inode,0) += my::scatrapara_->FDCheckEps();
        else
          epotnp_(inode) += my::scatrapara_->FDCheckEps();
      }

      // calculate element right-hand side vector for perturbed state
      Sysmat(ele,emat_dummy,erhs_perturbed,subgrdiff_dummy);

      // Now we compare the difference between the current entries in the element matrix
      // and their finite difference approximations according to
      // entries ?= (-erhs_perturbed + erhs_original) / epsilon

      // Note that the element right-hand side equals the negative element residual.
      // To account for errors due to numerical cancellation, we additionally consider
      // entries - erhs_original / epsilon ?= -erhs_perturbed / epsilon

      // Note that we still need to evaluate the first comparison as well. For small entries in the element
      // matrix, the second comparison might yield good agreement in spite of the entries being wrong!
      for(int row=0; row<my::numdofpernode_*my::nen_; ++row)
      {
        // get current entry in original element matrix
        const double entry = emat(row,col);

        // finite difference suggestion (first divide by epsilon and then subtract for better conditioning)
        const double fdval = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps() + erhs(row) / my::scatrapara_->FDCheckEps();

        // confirm accuracy of first comparison
        if(abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in first comparison
        const double abserr1 = entry - fdval;
        if(abs(abserr1) > abs(maxabserr))
          maxabserr = abserr1;
        double relerr1(0.);
        if(abs(entry) > 1.e-17)
          relerr1 = abserr1 / abs(entry);
        else if(abs(fdval) > 1.e-17)
          relerr1 = abserr1 / abs(fdval);
        if(abs(relerr1) > abs(maxrelerr))
          maxrelerr = relerr1;

        // evaluate first comparison
        if(abs(relerr1) > my::scatrapara_->FDCheckTol())
        {
          if(!counter)
            std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
          std::cout << "emat[" << row << "," << col << "]:  " << entry << "   ";
          std::cout << "finite difference suggestion:  " << fdval << "   ";
          std::cout << "absolute error:  " << abserr1 << "   ";
          std::cout << "relative error:  " << relerr1 << std::endl;

          counter++;
        }

        // first comparison OK
        else
        {
          // left-hand side in second comparison
          const double left  = entry - erhs(row) / my::scatrapara_->FDCheckEps();

          // right-hand side in second comparison
          const double right = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps();

          // confirm accuracy of second comparison
          if(abs(right) > 1.e-17 and abs(right) < 1.e-15)
            dserror("Finite difference check involves values too close to numerical zero!");

          // absolute and relative errors in second comparison
          const double abserr2 = left - right;
          if(abs(abserr2) > abs(maxabserr))
            maxabserr = abserr2;
          double relerr2(0.);
          if(abs(left) > 1.e-17)
            relerr2 = abserr2 / abs(left);
          else if(abs(right) > 1.e-17)
            relerr2 = abserr2 / abs(right);
          if(abs(relerr2) > abs(maxrelerr))
            maxrelerr = relerr2;

          // evaluate second comparison
          if(abs(relerr2) > my::scatrapara_->FDCheckTol())
          {
            if(!counter)
              std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
            std::cout << "emat[" << row << "," << col << "]-erhs[" << row << "]/eps:  " << left << "   ";
            std::cout << "-erhs_perturbed[" << row << "]/eps:  " << right << "   ";
            std::cout << "absolute error:  " << abserr2 << "   ";
            std::cout << "relative error:  " << relerr2 << std::endl;

            counter++;
          }
        }
      }
    }
  }

  // screen output in case finite difference check is passed
  if(!counter)
    std::cout << " --> PASSED WITH MAXIMUM ABSOLUTE ERROR " << maxabserr << " AND MAXIMUM RELATIVE ERROR " << maxrelerr << std::endl;

  // undo perturbations of state variables
  for(int k=0; k<my::numscal_; ++k)
    for(int i=0; i<my::nen_; ++i)
      my::ephinp_[k](i,0) = ephinp_original[k](i,0);
  for(int i=0; i<my::nen_; ++i)
    epotnp_(i) = epotnp_original(i);
  if(my::scatraparatimint_->IsGenAlpha())
    for(int k=0; k<my::numscal_; ++k)
      for(int i=0; i<my::nen_; ++i)
        my::ehist_[k](i,0)  = ehist_original[k](i,0);

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;
