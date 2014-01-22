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
//TODO: SCATRA_ELE_CLEANING: Wie bekommen wir das sonst?
#include "../drt_lib/drt_discret.H"  // for time curve in body force
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on

#include "../drt_inpar/inpar_elch.H"
#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"

//#include "scatra_ele_parameter_timint.H"
//
//#include "../drt_lib/drt_utils.H"

//#include "../drt_lib/drt_discret.H"


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

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");

  switch (action)
  {
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

    // integrations points and weights
    DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    /*----------------------------------------------------------------------*/
    // element integration loop
    /*----------------------------------------------------------------------*/
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // loop starts at k=numscal_ !!
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

      // current as a solution variable
      if(cursolvar_)
      {
        for(int idim=0;idim<my::nsd_;++idim)
        {
          // loop starts at k=numscal_ !!
          for (int vi=0; vi<my::nen_; ++vi)
          {
            const double v = fac*my::funct_(vi); // no density required here
            const int fvi = vi*my::numdofpernode_+my::numscal_+1+idim;

            for (int ui=0; ui<my::nen_; ++ui)
            {
              const int fui = ui*my::numdofpernode_+my::numscal_+1+idim;

              elemat1_epetra(fvi,fui) += v*my::funct_(ui);
            }
          }
        }
      }
    }

    // set zero for the rhs of the potential
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+my::numscal_;

      elevec1_epetra[fvi] = 0.0; // zero out!

      //TODO: SCATRA_ELE_CLEANING: zero for current

    }

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

      // access control parameter for flux calculation
      INPAR::SCATRA::FluxType fluxtype = DRT::INPUT::get<INPAR::SCATRA::FluxType>(params, "fluxtype");

      // we always get an 3D flux vector for each node
      LINALG::Matrix<3,my::nen_> eflux(true);

      // do a loop for systems of transported scalars
      for (int k = 0; k<my::numscal_; ++k)
      {
        // calculate flux vectors for actual scalar
        eflux.Clear();
        CalculateFlux(eflux,ele,elchtype,fluxtype,k);
        // assembly
        for (int inode=0;inode<my::nen_;inode++)
        {
          const int fvi = inode*my::numdofpernode_+k;
          elevec1_epetra[fvi]+=eflux(0,inode);
          elevec2_epetra[fvi]+=eflux(1,inode);
          elevec3_epetra[fvi]+=eflux(2,inode);
        }
      } // loop over numscal

      break;
    }
  case SCATRA::calc_error:
  {
    // check if length suffices
    if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

    //TODO: SCATRA_ELE_CLEANING: DiffCondParameter zu Parameterliste hinzufügen
    // set specific parameter used in diffusion conduction formulation
    // this method need to be located inside ELCH
    //DiffCondParams(ele, params);

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
      elchtype,
      params,
      elevec1_epetra);

    break;
  }
  case SCATRA::calc_elch_conductivity:
  {
    //TODO: SCATRA_ELE_CLEANING: DiffCondParameter zu Parameterliste hinzufügen
    // set specific parameter used in diffusion conduction formulation
    // this method need to be located inside ELCH
    //DiffCondParams(ele, params);

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

    CalculateConductivity(ele,elchtype,elevec1_epetra);
    break;
  }
  case SCATRA::calc_elch_initial_potential:
  {
    //TODO: SCATRA_ELE_CLEANING: DiffCondParameter zu Parameterliste hinzufügen
    // set specific parameter used in diffusion conduction formulation
    // this method need to be located inside ELCH
    //DiffCondParams(ele, params);

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

    CalculateElectricPotentialField(ele,elchtype,elemat1_epetra,elevec1_epetra);

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


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     gjb 06/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalculateFlux(
LINALG::Matrix<3,my::nen_>&         flux,
const DRT::Element*             ele,
INPAR::ELCH::ElchType           elchtype,
const INPAR::SCATRA::FluxType   fluxtype,
const int                       k
)
{
  // set constants
  const double frt = elchpara_->FRT();

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

  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // get material parameters (evaluation at element center)
  if (not my::scatrapara_->MatGP()) my::GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

  // integration rule
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad< intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // get material parameters (evaluation at integration point)
    if (my::scatrapara_->MatGP()) my::GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

    // get velocity at integration point
    LINALG::Matrix<my::nsd_,1> velint(true);
    LINALG::Matrix<my::nsd_,1> convelint(true);
    velint.Multiply(my::evelnp_,my::funct_);
    convelint.Multiply(my::econvelnp_,my::funct_);

    // get scalar at integration point
    const double phi = my::funct_.Dot(my::ephinp_[k]);

    // get gradient of scalar at integration point
    LINALG::Matrix<my::nsd_,1> gradphi(true);
    gradphi.Multiply(my::derxy_,my::ephinp_[k]);

    // get gradient of electric potential at integration point
    LINALG::Matrix<my::nsd_,1> gradpot(true);
    gradpot.Multiply(my::derxy_,epotnp_);

    // allocate and initialize!
    LINALG::Matrix<my::nsd_,1> q(true);

    if(elchtype == INPAR::ELCH::elchtype_diffcond)
    {
      // Be careful: - Evaluation of phi only for actual scalar
      //             - loop over numscal instead of numdof
      dserror("domain fluy not implemented for diffusion-conduction formulation");
    }
    else
    {
      // add different flux contributions as specified by user input
      switch (fluxtype)
      {
      case INPAR::SCATRA::flux_total_domain:
        // convective flux contribution
        q.Update(densnp*phi,convelint);

        // no break statement here!
      case INPAR::SCATRA::flux_diffusive_domain:
        // diffusive flux contribution
        q.Update(-(my::diffmanager_->GetIsotropicDiff(k)),gradphi,1.0);

        q.Update(-my::diffmanager_->GetIsotropicDiff(k)*valence_[k]*frt*phi,gradpot,1.0);

        break;
      default:
        dserror("received illegal flag inside flux evaluation for whole domain"); break;
      };
      // q at integration point

      // integrate and assemble everything into the "flux" vector
      for (int vi=0; vi < my::nen_; vi++)
      {
        for (int idim=0; idim<my::nsd_ ;idim++)
        {
          flux(idim,vi) += fac*my::funct_(vi)*q(idim);
        } // idim
      } // vi
    }
  } // integration loop

  //set zeros for unused space dimensions
  for (int idim=my::nsd_; idim<3; idim++)
  {
    for (int vi=0;vi <my::nen_;vi++)
    {
      flux(idim,vi) = 0.0;
    }
  }

  return;
} // ScaTraCalc::CalculateFlux


/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution           gjb 10/08|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalErrorComparedToAnalytSolution(
  const DRT::Element*                   ele,
  const enum INPAR::ELCH::ElchType      elchtype,
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
  my::GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

  //TODO: SCATRA_ELE_CLEANING: DiffCondParameter zu Parameterliste hinzufügen
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
      const double d = frt*((my::diffmanager_->GetIsotropicDiff(0)*valence_[0]) - (my::diffmanager_->GetIsotropicDiff(1)*valence_[1]));
      if (abs(d) == 0.0) dserror("division by zero");
      const double D = frt*((valence_[0]*my::diffmanager_->GetIsotropicDiff(0)*my::diffmanager_->GetIsotropicDiff(1)) - (valence_[1]*my::diffmanager_->GetIsotropicDiff(1)*my::diffmanager_->GetIsotropicDiff(0)))/d;

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
      c(1) = (-valence_[0]/valence_[1])* c(0);
      // compute analytical solution for el. potential
      const double pot = ((my::diffmanager_->GetIsotropicDiff(1)-my::diffmanager_->GetIsotropicDiff(0))/d) * log(c(0)/c_0_0_0_t);

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
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

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
      c(1) = (-valence_[0]/valence_[1])* c(0);
      // compute analytical solution for el. potential
      const double d = frt*((my::diffmanager_->GetIsotropicDiff(0)*valence_[0]) - (my::diffmanager_->GetIsotropicDiff(1)*valence_[1]));
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
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // get values of transported scalars at integration point
      // and compute electroneutrality
      double deviation(0.0);
      for (int k=0; k<my::numscal_; ++k)
      {
        const double conint_k = my::funct_.Dot(my::ephinp_[k]);
        deviation += valence_[k]*conint_k;
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
  const enum INPAR::ELCH::ElchType  elchtype,
  Epetra_SerialDenseVector&         sigma
  )
{
  // calculate conductivity of electrolyte solution
  const double frt = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->FRT();

  my::EvalShapeFuncAndDerivsAtEleCenter();

  // get concentration of transported scalar k at integration point
  std::vector<double> conint(my::numscal_);
  for (int k = 0;k<my::numscal_;++k)
    conint[k] = my::funct_.Dot(my::ephinp_[k]);

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
  GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

  // compute the conductivity (1/(\Omega m) = 1 Siemens / m)
  double sigma_all(0.0);
  const double factor = frt*INPAR::ELCH::faraday_const; // = F^2/RT

  // Concentrated solution theory:
  // Conductivity given by a function is evaluated at bulk concentration
  if(elchtype == INPAR::ELCH::elchtype_diffcond)
  {
    // TODO: SCATRA_ELE_CLEANING: Gibt es eine besser Alternaive
    Teuchos::RCP<MAT::Material> material = ele->Material();
    const Teuchos::RCP<const MAT::ElchMat>& actmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);
    // loop over single phases
    for (int iphase=0; iphase < actmat->NumPhase();++iphase)
    {
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->PhaseById(phaseid);

      if(singlemat->MaterialType() == INPAR::MAT::m_elchphase)
      {
        const MAT::ElchPhase* actsinglemat = static_cast<const MAT::ElchPhase*>(singlemat.get());
        sigma_all = actsinglemat->ComputeConductivity(conint[0]);
      }
      else
        dserror("Conductivity has to be defined in m_elchphase! There is no material m_elchphase");
    }
  }
  // Dilute solution theory:
  // Conductivity is computed by
  // sigma = F^2/RT*Sum(z_k^2 D_k c_k)
  else
  {
    for(int k=0; k < my::numscal_; k++)
    {
      double sigma_k = factor*valence_[k]*my::diffmanager_->GetIsotropicDiff(k)*valence_[k]*conint[k];
      sigma[k] += sigma_k; // insert value for this ionic species
      sigma_all += sigma_k;

      // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
      if(elchtype==INPAR::ELCH::elchtype_enc_pde_elim)
      {
        sigma_all += factor*my::diffmanager_->GetIsotropicDiff(k)*valence_[k]*valence_[k]*(-conint[k]);
      }
    }
  }

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
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs
  )
{
  // calculate conductivity of electrolyte solution
  const double frt = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->FRT();

  const double faraday = INPAR::ELCH::faraday_const;

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
    GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc,iquad);

    // get concentration of transported scalar k at integration point
    std::vector<double> conint(true);
    for (int k = 0;k<my::numscal_;++k)
      conint[k] = my::funct_.Dot(my::ephinp_[k]);

    // loop over all transported scalars
    std::vector<LINALG::Matrix<my::nsd_,1> > gradphi(true);
    for (int k = 0; k < my::numscal_;++k)
      gradphi[k].Multiply(my::derxy_,my::ephinp_[k]);

    // get gradient of electric potential at integration point
    LINALG::Matrix<my::nsd_,1> gradpot(true);
    gradpot.Multiply(my::derxy_,epotnp_);

    //TODO: SCATRA_ELE_CLEANING: Grundsätzliches Konzept
    if(elchtype != INPAR::ELCH::elchtype_diffcond)
    {
      double sigmaint(0.0);
      for (int k=0; k<my::numscal_; ++k)
      {
        double sigma_k = frt*valence_[k]*my::diffmanager_->GetIsotropicDiff(k)*valence_[k]*conint[k];
        sigmaint += sigma_k;

        // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
        if(elchtype==INPAR::ELCH::elchtype_enc_pde_elim)
          sigmaint += frt*valence_[k]*my::diffmanager_->GetIsotropicDiff(k)*valence_[my::numscal_]*(-conint[k]);

        // diffusive terms on rhs
        const double vrhs = fac*my::diffmanager_->GetIsotropicDiff(k)*valence_[k];
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = vi*my::numdofpernode_+my::numscal_;
          double laplawf(0.0);
          my::GetLaplacianWeakFormRHS(laplawf,gradphi[k],vi);
          erhs[fvi] -= vrhs*laplawf;
          // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
          if(elchtype==INPAR::ELCH::elchtype_enc_pde_elim)
            erhs[fvi] -= -fac*valence_[k]*my::diffmanager_->GetIsotropicDiff(my::numscal_)*laplawf;
        }

        // provide something for conc. dofs: a standard mass matrix
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int    fvi = vi*my::numdofpernode_+k;
          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+k;
            emat(fvi,fui) += fac*my::funct_(vi)*my::funct_(ui);
          }
        }
      } // for k

      // ----------------------------------------matrix entries
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int    fvi = vi*my::numdofpernode_+my::numscal_;
        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+my::numscal_;
          double laplawf(0.0);
          my::GetLaplacianWeakForm(laplawf,ui,vi);
          emat(fvi,fui) += fac*sigmaint*laplawf;
        }

        double laplawf(0.0);
        my::GetLaplacianWeakFormRHS(laplawf,gradpot,vi);
        erhs[fvi] -= fac*sigmaint*laplawf;
      }
    }
    else
    {
      //TODO: SCATRA_ELE_CLEANING
      if(diffcondmat_==INPAR::ELCH::diffcondmat_diffcond)
        dserror("The function CalcInitialPotential is only implemented for Newman materials");

      for (int k=0; k<my::numscal_; ++k)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = vi*my::numdofpernode_+my::numscal_;
          double laplawf(0.0);
          my::GetLaplacianWeakFormRHS(laplawf,gradphi[k],vi);

          for (int iscal=0; iscal < my::numscal_; ++iscal)
          {
            erhs[fvi] -= fac*epstort_[0]/faraday/frt*cond_[0]*(therm_[0])*((a_+(b_*trans_[iscal]))/c_)/conint[iscal]*laplawf;
          }
        }

        // provide something for conc. dofs: a standard mass matrix
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int    fvi = vi*my::numdofpernode_+k;
          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+k;
            emat(fvi,fui) += fac*my::funct_(vi)*my::funct_(ui);
          }
        }
      } // for k

      // ----------------------------------------matrix entries
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int    fvi = vi*my::numdofpernode_+my::numscal_;
        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+my::numscal_;
          double laplawf(0.0);
          my::GetLaplacianWeakForm(laplawf,ui,vi);
          emat(fvi,fui) += fac*epstort_[0]/faraday*cond_[0]*laplawf;
        }

        double laplawf(0.0);
        my::GetLaplacianWeakFormRHS(laplawf,gradpot,vi);
        erhs[fvi] -= fac*epstort_[0]/faraday*cond_[0]*laplawf;
      }
    }
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
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;
