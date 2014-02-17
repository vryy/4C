/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_general.cpp

\brief Internal implementation of ScaTra element

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc.H"

#include "scatra_ele.H"
#include "scatra_ele_action.H"

#include "scatra_ele_parameter.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"

#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on


/*----------------------------------------------------------------------*
 * Action type: EvaluateService
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalc<distype>::EvaluateService(
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
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  if (scatrapara_->IsAle())
  {
    const Teuchos::RCP<Epetra_MultiVector> dispnp = params.get< Teuchos::RCP<Epetra_MultiVector> >("dispnp");
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_);
    // add nodal displacements to point coordinates
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,myknots_,weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  // set element id
  eid_ = ele->Id();

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch (action)
  {
  // calculate time derivative for time value t_0
  case SCATRA::calc_initial_time_deriv:
  {
    // calculate matrix and rhs
    CalcInitialTimeDerivative(
      ele,
      elemat1_epetra,
      elevec1_epetra,
      params,
      discretization,
      lm
      );
    break;
  }
  case SCATRA::integrate_shape_functions:
  {
    // calculate integral of shape functions
    const Epetra_IntSerialDenseVector dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
    IntegrateShapeFunctions(ele,elevec1_epetra,dofids);

    break;
  }
  case SCATRA::calc_flux_domain:
  {
    // get velocity values at the nodes
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,velocity,nsd_);
    const Teuchos::RCP<Epetra_MultiVector> convelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,convelocity,nsd_);

    // need current values of transported scalar
    // -> extract local values from global vectors
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // fill all element arrays
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
      }
    } // for i

    // access control parameter for flux calculation
    INPAR::SCATRA::FluxType fluxtype = DRT::INPUT::get<INPAR::SCATRA::FluxType>(params, "fluxtype");

    // we always get an 3D flux vector for each node
    LINALG::Matrix<3,nen_> eflux(true);

    // do a loop for systems of transported scalars
    for (int k = 0; k<numscal_; ++k)
    {
      // calculate flux vectors for actual scalar
      eflux.Clear();
      CalculateFlux(eflux,ele,fluxtype,k);
      // assembly
      for (int inode=0;inode<nen_;inode++)
      {
        const int fvi = inode*numdofpernode_+k;
        elevec1_epetra[fvi]+=eflux(0,inode);
        elevec2_epetra[fvi]+=eflux(1,inode);
        elevec3_epetra[fvi]+=eflux(2,inode);
      }
    } // loop over numscal

    break;
  }
  case SCATRA::calc_mean_scalars:
  {
    // NOTE: add integral values only for elements which are NOT ghosted!
    if (ele->Owner() == discretization.Comm().MyPID())
    {
      // get flag for inverting
      bool inverting = params.get<bool>("inverting");

      // need current scalar vector
      // -> extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

      // fill all element arrays
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
        }
      } // for i

      // calculate scalars and domain integral
      CalculateScalars(ele,elevec1_epetra,inverting);
    }
    break;
  }
  // calculated filtered fields for calculation of turbulent Prandtl number
  // required for dynamic Smagorinsky model in scatra
  case SCATRA::calc_scatra_box_filter:
  {
    if (nsd_ == 3)
    {
      // extract scalar and velocity values from global vectors

      // scalar field
      Teuchos::RCP<const Epetra_Vector> scalar = discretization.GetState("scalar");
      if (scalar == Teuchos::null)
        dserror("Cannot get scalar!");
      std::vector<double> myscalar(lm.size());
      DRT::UTILS::ExtractMyValues(*scalar,myscalar,lm);
      for (int i=0;i<nen_;++i)
          ephinp_[0](i,0) = myscalar[i];
      // velocity field
      const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity");
      DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,velocity,nsd_);

      // get thermodynamic pressure
      double thermpress = 0.0;
      if (scatrapara_->ScaTraType() == INPAR::SCATRA::scatratype_loma)
        thermpress = params.get<double>("thermpress");

      // initialize the contribution of this element to the patch volume to zero
      double volume_contribution = 0.0;
      // initialize the contributions of this element to the filtered scalar quantities
      double dens_hat = 0.0;
      double temp_hat = 0.0;
      double dens_temp_hat = 0.0;
      double phi2_hat=0.0;
      double phiexpression_hat=0.0;
      // get pointers for vector quantities
      Teuchos::RCP<std::vector<double> > vel_hat = params.get<Teuchos::RCP<std::vector<double> > >("vel_hat");
      Teuchos::RCP<std::vector<double> > densvel_hat = params.get<Teuchos::RCP<std::vector<double> > >("densvel_hat");
      Teuchos::RCP<std::vector<double> > densveltemp_hat = params.get<Teuchos::RCP<std::vector<double> > >("densveltemp_hat");
      Teuchos::RCP<std::vector<double> > densstraintemp_hat = params.get<Teuchos::RCP<std::vector<double> > >("densstraintemp_hat");
      Teuchos::RCP<std::vector<double> > phi_hat = params.get<Teuchos::RCP<std::vector<double> > >("phi_hat");
      Teuchos::RCP<std::vector<std::vector<double> > > alphaijsc_hat = params.get<Teuchos::RCP<std::vector<std::vector<double> > > >("alphaijsc_hat");
      // integrate the convolution with the box filter function for this element
      // the results are assembled onto the *_hat arrays
      switch (distype)
      {
      case DRT::Element::hex8:
      {
        scatra_apply_box_filter(
            thermpress,
            dens_hat,
            temp_hat,
            dens_temp_hat,
            phi2_hat,
            phiexpression_hat,
            vel_hat,
            densvel_hat,
            densveltemp_hat,
            densstraintemp_hat,
            phi_hat,
            alphaijsc_hat,
            volume_contribution,
            ele);
        break;
      }
      default:
      {
        dserror("Unknown element type for box filter application\n");
      }
      }

      // hand down the volume contribution to the time integration algorithm
      params.set<double>("volume_contribution",volume_contribution);
      // as well as the filtered scalar quantities
      params.set<double>("dens_hat",dens_hat);
      params.set<double>("temp_hat",temp_hat);
      params.set<double>("dens_temp_hat",dens_temp_hat);
      params.set<double>("phi2_hat",phi2_hat);
      params.set<double>("phiexpression_hat",phiexpression_hat);

    }// end if (nsd == 3)
    else dserror("action 'calc_scatra_box_filter' is 3D specific action");

    break;
  }
  // calculate turbulent prandtl number of dynamic Smagorinsky model
  case SCATRA::calc_turbulent_prandtl_number:
  {
    if (nsd_ == 3)
    {
      // get required quantities, set in dynamic Smagorinsky class
      Teuchos::RCP<Epetra_MultiVector> col_filtered_vel =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("col_filtered_vel");
      Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("col_filtered_dens_vel");
      Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_vel_temp =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("col_filtered_dens_vel_temp");
      Teuchos::RCP<Epetra_MultiVector> col_filtered_dens_rateofstrain_temp =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("col_filtered_dens_rateofstrain_temp");
      Teuchos::RCP<Epetra_Vector> col_filtered_temp =
        params.get<Teuchos::RCP<Epetra_Vector> >("col_filtered_temp");
      Teuchos::RCP<Epetra_Vector> col_filtered_dens =
        params.get<Teuchos::RCP<Epetra_Vector> >("col_filtered_dens");
      Teuchos::RCP<Epetra_Vector> col_filtered_dens_temp =
        params.get<Teuchos::RCP<Epetra_Vector> >("col_filtered_dens_temp");

      // initialize variables to calculate
      double LkMk   = 0.0;
      double MkMk   = 0.0;
      double xcenter  = 0.0;
      double ycenter  = 0.0;
      double zcenter  = 0.0;

      // calculate LkMk and MkMk
      switch (distype)
      {
      case DRT::Element::hex8:
      {
        scatra_calc_smag_const_LkMk_and_MkMk(
            col_filtered_vel                   ,
            col_filtered_dens_vel              ,
            col_filtered_dens_vel_temp         ,
            col_filtered_dens_rateofstrain_temp,
            col_filtered_temp                  ,
            col_filtered_dens                  ,
            col_filtered_dens_temp             ,
            LkMk                               ,
            MkMk                               ,
            xcenter                            ,
            ycenter                            ,
            zcenter                            ,
            ele                                );
        break;
      }
      default:
      {
        dserror("Unknown element type for box filter application\n");
      }
      }

      // set Prt without averaging (only clipping)
      // calculate inverse of turbulent Prandtl number times (C_s*Delta)^2
      double inv_Prt = 0.0;
      if (abs(MkMk) < 1E-16)
      {
        //std::cout << "warning: abs(MkMk) < 1E-16 -> set inverse of turbulent Prandtl number to zero!"  << std::endl;
        inv_Prt = 0.0;
      }
      else  inv_Prt = LkMk / MkMk;
      if (inv_Prt<0.0)
        inv_Prt = 0.0;

      // set all values in parameter list
      params.set<double>("LkMk",LkMk);
      params.set<double>("MkMk",MkMk);
      params.set<double>("xcenter",xcenter);
      params.set<double>("ycenter",ycenter);
      params.set<double>("zcenter",zcenter);
      params.set<double>("ele_Prt",inv_Prt);
    }
    else dserror("action 'calc_turbulent_prandtl_number' is a 3D specific action");

    break;
  }
  case SCATRA::calc_vreman_scatra:
  {
    if (nsd_ == 3)
    {

      Teuchos::RCP<Epetra_MultiVector> col_filtered_phi =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("col_filtered_phi");
      Teuchos::RCP<Epetra_Vector> col_filtered_phi2 =
        params.get<Teuchos::RCP<Epetra_Vector> >("col_filtered_phi2");
      Teuchos::RCP<Epetra_Vector> col_filtered_phiexpression =
        params.get<Teuchos::RCP<Epetra_Vector> >("col_filtered_phiexpression");
      Teuchos::RCP<Epetra_MultiVector> col_filtered_alphaijsc =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("col_filtered_alphaijsc");

      // initialize variables to calculate
      double dt_numerator=0.0;
      double dt_denominator=0.0;

      // calculate LkMk and MkMk
      switch (distype)
      {
      case DRT::Element::hex8:
      {
        scatra_calc_vreman_dt(
            col_filtered_phi                   ,
            col_filtered_phi2              ,
            col_filtered_phiexpression         ,
            col_filtered_alphaijsc,
            dt_numerator,
            dt_denominator,
            ele                                );
        break;
      }
      default:
      {
        dserror("Unknown element type for vreman scatra application\n");
      }
      break;
      }

      elevec1_epetra(0)=dt_numerator;
      elevec1_epetra(1)=dt_denominator;
    }
    else dserror("action 'calc_vreman_scatra' is a 3D specific action");


    break;
  }
  // calculate normalized subgrid-diffusivity matrix
  case SCATRA::calc_subgrid_diffusivity_matrix:
  {
    // calculate mass matrix and rhs
    CalcSubgrDiffMatrix(
      ele,
      elemat1_epetra);

    break;
  }
  // calculate mean Cai of multifractal subgrid-scale modeling approach
  case SCATRA::calc_mean_Cai:
  {
    // get nodel velocites
    const Teuchos::RCP<Epetra_MultiVector> convelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,convelocity,nsd_);
    // get phi for material parameters
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null)
      dserror("Cannot get state vector 'phinp'");
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
    // fill all element arrays
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
      }
    }

    if (scatrapara_->TurbModel() != INPAR::FLUID::multifractal_subgrid_scales)
      dserror("Multifractal_Subgrid_Scales expected");

    double Cai = 0.0;
    double vol = 0.0;
    // calculate Cai and volume, do not include elements of potential inflow section
    if (scatrapara_->AdaptCsgsPhi() and scatrapara_->Nwl() and (not SCATRA::InflowElement(ele)))
    {
      // use one-point Gauss rule to do calculations at the element center
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);
      vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints,0);

      // adopt integrations points and weights for gauss point evaluation of B
      if (scatrapara_->BD_Gp())
      {
        DRT::UTILS::IntPointsAndWeights<nsd_> gauss_intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
        intpoints = gauss_intpoints;
      }

      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

        // density at t_(n)
        double densn(1.0);
        // density at t_(n+1) or t_(n+alpha_F)
        double densnp(1.0);
        // density at t_(n+alpha_M)
        double densam(1.0);

        // diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma

        diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManager(numscal_));

        // fluid viscosity
        double visc(0.0);

        // get material
        GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);

        // get velocity at integration point
        LINALG::Matrix<nsd_,1> convelint(true);
        convelint.Multiply(econvelnp_,funct_);

        // calculate characteristic element length
        double hk = CalcRefLength(vol,convelint);

        // estimate norm of strain rate
        double strainnorm = GetStrainRate(econvelnp_);
        strainnorm /= sqrt(2.0);

        // get Re from strain rate
        double Re_ele_str = strainnorm * hk * hk * densnp / visc;
        if (Re_ele_str < 0.0) dserror("Something went wrong!");
        // ensure positive values
        if (Re_ele_str < 1.0)
          Re_ele_str = 1.0;

        // calculate corrected Cai
        //           -3/16
        //  =(1 - (Re)   )
        //
        Cai += (1.0-pow(Re_ele_str,-3.0/16.0)) * fac;
      }
    }

    // hand down the Cai and volume contribution to the time integration algorithm
    params.set<double>("Cai_int",Cai);
    params.set<double>("ele_vol",vol);

    break;
  }
  // calculate dissipation introduced by stabilization and turbulence models
  case SCATRA::calc_dissipation:
  {
    CalcDissipation(params,
                    ele,
                    discretization,
                    lm);
    break;
  }
  case SCATRA::calc_integr_objf:
  {
    // don't want ghosted elements
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      double integralvalue = params.get<double>("integr_objf");
      const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
      double fac = 0.0;
      // integration loop
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        fac += EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);
      }
      double aux(1.0);
      GetMaterialParams(ele,aux,aux,aux,diffmanager_,reamanager_,aux);
      double reacoeff = reamanager_->GetReaCoeff(0);

      // integralvalue = 0.5 * ( diffus_[0]*diffus_[0] + reacoeff_[0]*reacoeff_[0] ) * fac;
      integralvalue = 0.5 * ( reacoeff * reacoeff ) * fac;
      params.set<double>("integr_objf",integralvalue);
    }
    break;
  }
  case SCATRA::calc_integr_pat_rhsvec:
  {
    // extract local values from the global vectors w and phi
    Teuchos::RCP<const Epetra_Vector> rhsvec = discretization.GetState("rhsnodebasedvals");
    if (rhsvec==Teuchos::null)
      dserror("Cannot get state vector 'rhsnodebasedvals' ");
    std::vector<double> myrhsvec(lm.size());
    DRT::UTILS::ExtractMyValues(*rhsvec,myrhsvec,lm);

    Epetra_SerialDenseVector rhs(nen_);
    for(int i=0; i<nen_; ++i)
      rhs(i) = myrhsvec[i];

    Epetra_SerialDenseMatrix mass(nen_,nen_);
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = fac*funct_(vi); // no density required here

        for (int ui=0; ui<nen_; ++ui)
        {
          mass(vi,ui) += v*funct_(ui);
        }
      }
    }

    elevec1_epetra.Multiply('N','N',1.0,mass,rhs,0.0);

    break;
  }
  case SCATRA::calc_integr_grad_reac:
  {
    // don't want ghosted elements
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      double integr_grad_reac = params.get<double>("integr_grad_reac");
      double integralvalue = 0.0;
      bool sign = params.get<bool>("signum_mu");
      double sign_fac = 1.0;
      if(sign) sign_fac = -1.0;

      // extract local values from the global vectors w and phi
      Teuchos::RCP<const Epetra_Vector> w = discretization.GetState("w");
      Teuchos::RCP<const Epetra_Vector> psi = discretization.GetState("psi");
      Teuchos::RCP<const Epetra_Vector> phi = discretization.GetState("phi");
      if (w==Teuchos::null || phi==Teuchos::null || psi==Teuchos::null)
        dserror("Cannot get state vector 'w' and/or 'phi' and/or 'psi'");
      std::vector<double> myw(lm.size());
      std::vector<double> myphi(lm.size());
      std::vector<double> mypsi(lm.size());
      DRT::UTILS::ExtractMyValues(*w,myw,lm);
      DRT::UTILS::ExtractMyValues(*phi,myphi,lm);
      DRT::UTILS::ExtractMyValues(*psi,mypsi,lm);

      // calculate integral
      // integrations points and weights
      const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
      double fac = 0.0;
      // integration loop
      double temp = 0.0;
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);
        temp += fac;
        double sumphi = 0.0;
        double sumw   = 0.0;
        for (int i=0; i<nen_; ++i)
        {
          sumphi += funct_(i) * myphi[i];
          sumw   += funct_(i) * (sign_fac * myw[i] + mypsi[i]);
        }
        integralvalue += fac * sumphi * sumw;
      }
      double aux(1.0);
      GetMaterialParams(ele,aux,aux,aux,diffmanager_,reamanager_,aux);
      double reacoeff = reamanager_->GetReaCoeff(0);
      integralvalue *= reacoeff;
      integralvalue += temp * params.get<double>("regular_mu") * reacoeff;
      integr_grad_reac += integralvalue;
      params.set<double>("integr_grad_reac",integr_grad_reac);
    }
    break;
  }
  default:
  {
    dserror("Not acting on this action. Forgot implementation?");
    break;
  }
  } // switch(action)

  return 0;
}


/*----------------------------------------------------------------------------*
 | calculate mass matrix + rhs for initial time derivative calc.     gjb 03/12|
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcInitialTimeDerivative(
  DRT::ELEMENTS::Transport*             ele,
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  Teuchos::ParameterList&               params,
  DRT::Discretization&                  discretization,
  const std::vector<int>&               lm
  )
{
  // dummy matrix + vectors required for Evaluate() call (zero size)
  Epetra_SerialDenseMatrix  elemat2_epetra=Teuchos::null;
  Epetra_SerialDenseVector  elevec2_epetra=Teuchos::null;
  Epetra_SerialDenseVector  elevec3_epetra=Teuchos::null;

  Evaluate(
      ele,
      params,
      discretization,
      lm,
      emat,
      elemat2_epetra,
      erhs,
      elevec2_epetra,
      elevec3_epetra
  );

  // undo the matrix from the standard call, only a mass matrix is needed here created below
  emat.Scale(0.0);

  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0);

  //------------------------------------------------------------------------------------
  // get material parameters and stabilization parameters (evaluation at element center)
  //------------------------------------------------------------------------------------
    // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // the stabilisation parameters (one per transported scalar)
  std::vector<double> tau(numscal_,0.0);

  if (not scatrapara_->MatGP() or not scatrapara_->TauGP())
  {

    GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);

    if (not scatrapara_->TauGP())
    {
      // get velocity at element center
      LINALG::Matrix<nsd_,1> velint(true);
      LINALG::Matrix<nsd_,1> convelint(true);
      velint.Multiply(evelnp_,funct_);
      convelint.Multiply(econvelnp_,funct_);

      for (int k = 0;k<numscal_;++k) // loop of each transported scalar
      {
        // calculation of stabilization parameter at element center
        CalcTau(tau[k],diffmanager_->GetIsotropicDiff(k),reamanager_->GetReaCoeff(k),densnp,convelint,vol,k);
      }
    }
  }

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  /*----------------------------------------------------------------------*/
  // element integration loop
  /*----------------------------------------------------------------------*/
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (scatrapara_->MatGP()) GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);

    //------------ get values of variables at integration point
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // get phi at integration point for all scalars
      double phiint = funct_.Dot(ephinp_[k]);

      // get velocity at integration point
      LINALG::Matrix<nsd_,1> velint(true);
      LINALG::Matrix<nsd_,1> convelint(true);
      velint.Multiply(evelnp_,funct_);
      convelint.Multiply(econvelnp_,funct_);

      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      LINALG::Matrix<nen_,1> conv(true);
      conv.MultiplyTN(derxy_,convelint);

      // scalar at integration point at time step n+1
      const double phinp = funct_.Dot(ephinp_[k]);

      // velocity divergence required for conservative form
      double vdiv(0.0);
      if (scatrapara_->IsConservative()) GetDivergence(vdiv,evelnp_);

      // diffusive part used in stabilization terms
      LINALG::Matrix<nen_,1> diff(true);
      // diffusive term using current scalar value for higher-order elements
      if (use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        GetLaplacianStrongForm(diff);
        diff.Scale(diffmanager_->GetIsotropicDiff(k));
      }

      // calculation of stabilization parameter at integration point
      if (scatrapara_->TauGP()) CalcTau(tau[k],diffmanager_->GetIsotropicDiff(k),reamanager_->GetReaCoeff(k),densnp,convelint,vol,k);

      const double fac_tau = fac*tau[k];

      //----------------------------------------------------------------
      // element matrix: transient term
      //----------------------------------------------------------------
      // transient term
      CalcMatMass(emat,k,fac,densam);

      //----------------------------------------------------------------
      // element matrix: stabilization of transient term
      //----------------------------------------------------------------
      if(scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      {
        // subgrid-scale velocity (dummy)
        LINALG::Matrix<nen_,1> sgconv(true);
        CalcMatMassStab(emat,k,fac_tau,densam,densnp,conv,sgconv,diff);

        // remove convective stabilization of inertia term
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;
          erhs(fvi)+=fac_tau*densnp*conv(vi)*densnp*phiint;
        }
      }

      if (scatraparatimint_->IsIncremental())
      {
         // zero for quantities only required for gen-alpha which is not the case here
         // fac->-fac to change sign of rhs
         CalcRHSLinMass(erhs,k,0.0,-fac,0.0,densnp,phinp,0.0);
      }
      else
        dserror("Must be incremental!");

    } // loop over each scalar k
  } // integration loop

  // correct scaling of rhs (after subtraction!!!!)
  erhs.Scale(1.0/scatraparatimint_->TimeFac());

  return;
} // ScaTraEleCalc::CalcInitialTimeDerivative()


/*----------------------------------------------------------------------*
 |  Integrate shape functions over domain (private)           gjb 07/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::IntegrateShapeFunctions(
  const DRT::Element*             ele,
  Epetra_SerialDenseVector&       elevec1,
  const Epetra_IntSerialDenseVector& dofids
  )
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // safety check
  if (dofids.M() < numdofpernode_)
    dserror("Dofids vector is too short. Received not enough flags");

  // TODO (ehrl): order is not effective since integrating shape functions result always in the same result
  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,gpid);

    // compute integral of shape functions (only for dofid)
    for (int k=0;k<numdofpernode_;k++)
    {
      if (dofids[k] >= 0)
      {
        for (int node=0;node<nen_;node++)
        {
          elevec1[node*numdofpernode_+k] += funct_(node) * fac;
        }
      }
    }

  } //loop over integration points

  return;

} //ScaTraEleCalc::IntegrateShapeFunction


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     gjb 06/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalculateFlux(
LINALG::Matrix<3,nen_>&         flux,
const DRT::Element*             ele,
const INPAR::SCATRA::FluxType   fluxtype,
const int                       k
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

  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  // get material parameters (evaluation at element center)
  if (not scatrapara_->MatGP()) GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);

  // integration rule
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad< intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // get material parameters (evaluation at integration point)
    if (scatrapara_->MatGP()) GetMaterialParams(ele,densn,densnp,densam,diffmanager_,reamanager_,visc);

    // get velocity at integration point
    LINALG::Matrix<nsd_,1> velint(true);
    LINALG::Matrix<nsd_,1> convelint(true);
    velint.Multiply(evelnp_,funct_);
    convelint.Multiply(econvelnp_,funct_);

    // get scalar at integration point
    const double phi = funct_.Dot(ephinp_[k]);

    // get gradient of scalar at integration point
    LINALG::Matrix<nsd_,1> gradphi(true);
    gradphi.Multiply(derxy_,ephinp_[k]);

    // allocate and initialize!
    LINALG::Matrix<nsd_,1> q(true);

    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
    case INPAR::SCATRA::flux_total_domain:
      // convective flux contribution
      q.Update(densnp*phi,convelint);

      // no break statement here!
    case INPAR::SCATRA::flux_diffusive_domain:
      // diffusive flux contribution
      q.Update(-(diffmanager_->GetIsotropicDiff(k)),gradphi,1.0);

      break;
    default:
      dserror("received illegal flag inside flux evaluation for whole domain"); break;
    };
    // q at integration point

    // integrate and assemble everything into the "flux" vector
    for (int vi=0; vi < nen_; vi++)
    {
      for (int idim=0; idim<nsd_ ;idim++)
      {
        flux(idim,vi) += fac*funct_(vi)*q(idim);
      } // idim
    } // vi

  } // integration loop

  //set zeros for unused space dimensions
  for (int idim=nsd_; idim<3; idim++)
  {
    for (int vi=0; vi < nen_; vi++)
    {
      flux(idim,vi) = 0.0;
    }
  }

  return;
} // ScaTraCalc::CalculateFlux


/*----------------------------------------------------------------------*
|  calculate scalar(s) and domain integral                     vg 09/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalculateScalars(
const DRT::Element*             ele,
Epetra_SerialDenseVector&       scalars,
const bool                      inverting
  )
{
  // integrations points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // calculate integrals of (inverted) scalar(s) and domain
    if (inverting)
    {
      for (int i=0; i<nen_; i++)
      {
        const double fac_funct_i = fac*funct_(i);
        for (int k = 0; k < numscal_; k++)
        {
          if (std::abs(ephinp_[k](i,0))> EPS14)
            scalars[k] += fac_funct_i/ephinp_[k](i,0);
          else
            dserror("Division by zero");
        }
        // for domain volume
        scalars[numscal_] += fac_funct_i;
      }
    }
    else
    {
      for (int i=0; i<nen_; i++)
      {
        const double fac_funct_i = fac*funct_(i);
        for (int k = 0; k < numscal_; k++)
        {
          scalars[k] += fac_funct_i*ephinp_[k](i,0);
        }
        // for domain volume
        scalars[numscal_] += fac_funct_i;
      }
    }
  } // loop over integration points

return;
} // ScaTraEleCalc::CalculateScalars


/*----------------------------------------------------------------------*
 | calculate normalized subgrid-diffusivity matrix              vg 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalc<distype>::CalcSubgrDiffMatrix(
  const DRT::Element*           ele,
  Epetra_SerialDenseMatrix&     emat
  )
{
  /*----------------------------------------------------------------------*/
  // integration loop for one element
  /*----------------------------------------------------------------------*/
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    for (int k=0;k<numscal_;++k)
    {
      // set diffusion coeff to 1.0
      diffmanager_->SetIsotropicDiff(1.0,k);

      // calculation of diffusive element matrix
      double timefacfac = scatraparatimint_->TimeFac() * fac;
      CalcMatDiff(emat,k,timefacfac,diffmanager_);

      /*subtract SUPG term */
      //emat(fvi,fui) -= taufac*conv(vi)*conv(ui);
    }
  } // integration loop

  return;
} // ScaTraImpl::CalcSubgrDiffMatrix


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalc<DRT::Element::nurbs27>;
