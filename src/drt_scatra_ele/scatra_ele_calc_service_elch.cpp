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

#include "scatra_ele_calc_elch.H"

#include "scatra_ele.H"
#include "scatra_ele_action.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"  // for time curve in body force
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_globalproblem.H"  // consistency check of formulation and material

#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"


/*----------------------------------------------------------------------*
 * Action type: EvaluateService
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateService(
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
  // get element coordinates
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  // set element id
  my::eid_ = ele->Id();

  INPAR::ELCH::ElchType elchtype = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->ElchType();
  INPAR::ELCH::EquPot equpot = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->EquPot();

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");

  switch (action)
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
    // At this point time is not so critical since the CalcIntialTimeDerivative()
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
      this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

    // dynamic cast to elch-specific diffusion manager
    Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

    //----------------------------------------------------------------------
    // integration loop for one element
    //----------------------------------------------------------------------
    // integrations points and weights
    DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      //----------------------------------------------------------------------
      // get material parameters (evaluation at integration point)
      //----------------------------------------------------------------------
      if (my::scatrapara_->MatGP())
        this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc,iquad);

      // set internal variables
      varmanager_->SetInternalVariablesElch(my::funct_,my::derxy_,my::ephinp_,epotnp_,my::econvelnp_);

      SetFormulationSpecificInternalVariables(dme,varmanager_);

      // access control parameter for flux calculation
      INPAR::SCATRA::FluxType fluxtype = DRT::INPUT::get<INPAR::SCATRA::FluxType>(params, "fluxtype");

      // do a loop for systems of transported scalars
      for (int k = 0; k<my::numscal_; ++k)
        CalculateFlux(elevec1_epetra,elevec2_epetra,elevec3_epetra,fluxtype,k,fac,varmanager_,dme);

      if(elchpara_->EquPot()==INPAR::ELCH::equpot_divi)
        CalculateCurrent(elevec1_epetra,elevec2_epetra,elevec3_epetra,fluxtype,fac,varmanager_,dme);
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

    CalculateConductivity(ele,equpot,elevec1_epetra);
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

    CalculateElectricPotentialField(ele,equpot,elemat1_epetra,elevec1_epetra);

    break;
  }
  case SCATRA::calc_elch_electrode_kinetics:
  {
    dserror(" ");
    break;
  }
  default:
  {
    my::EvaluateService(ele,
                        params,
                        discretization,
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


/*-----------------------------------------------------------------------*
  |  Set scatra element parameter                             ehrl 01/14 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CheckElchElementParameter(
  DRT::ELEMENTS::Transport*  ele
  )
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

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
      if (elchpara_->CurSolVar()==true)
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
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

  //at the moment, there is only one analytical test problem available!
  if (DRT::INPUT::get<SCATRA::Action>(params,"action") != SCATRA::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // in the ALE case add nodal displacements
  if (my::scatrapara_->IsAle()) dserror("No ALE for Kwok & Wu error calculation allowed.");

  // set constants for analytical solution
  const double t = my:: scatraparatimint_->Time();
  const double frt = elchpara_->FRT();

  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // get material parameter (constants values)
  this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

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

  // integrations points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

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
      const double d = frt*((dme->GetIsotropicDiff(0)*dme->GetValence(0)) - (dme->GetIsotropicDiff(1)*dme->GetValence(1)));
      if (abs(d) == 0.0) dserror("division by zero");
      const double D = frt*((dme->GetValence(0)*dme->GetIsotropicDiff(0)*dme->GetIsotropicDiff(1)) - (dme->GetValence(1)*dme->GetIsotropicDiff(1)*dme->GetIsotropicDiff(0)))/d;

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
      c(1) = (-dme->GetValence(0)/dme->GetValence(1))* c(0);
      // compute analytical solution for el. potential
      const double pot = ((dme->GetIsotropicDiff(1)-dme->GetIsotropicDiff(0))/d) * log(c(0)/c_0_0_0_t);

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
      c(1) = (-dme->GetValence(0)/dme->GetValence(1))* c(0);
      // compute analytical solution for el. potential
      const double d = frt*((dme->GetIsotropicDiff(0)*dme->GetValence(0)) - (dme->GetIsotropicDiff(1)*dme->GetValence(1)));
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
        deviation += dme->GetValence(k)*conint_k;
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
  this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

  // compute the conductivity (1/(\Omega m) = 1 Siemens / m)
  double sigma_all(0.0);

  GetConductivity(equpot,sigma_all, sigma);

  // conductivity based on ALL ionic species (even eliminated ones!)
  sigma[my::numscal_] += sigma_all;

  return;

} //ScaTraEleCalcElch()


/*----------------------------------------------------------------------*
  |  CalculateElectricPotentialField (ELCH) (private)          gjb 04/10 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalculateElectricPotentialField(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::EquPot    equpot,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs
  )
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

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
    this->GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc,iquad);

    // set internal variables
    varmanager_->SetInternalVariablesElch(my::funct_,my::derxy_,my::ephinp_,epotnp_,my::econvelnp_);

    SetFormulationSpecificInternalVariables(dme,varmanager_);

    CalMatAndRhsElectricPotentialField(varmanager_,equpot,emat,erhs,fac,dme);
  } // integration loop

  return;

} //ScaTraImpl<distype>::CalculateElectricPotentialField


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
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
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;
