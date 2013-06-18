/*----------------------------------------------------------------------*/
/*!
  \file scatra_ele_impl.cpp

  \brief Internal implementation of scalar transport elements

  <pre>
  Maintainer: Andreas Ehrl
  ehrl@lnm.mw.tum.de
  http://www.lnm.mw.tum.de
  089 - 289-15252
  </pre>
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#include "scatra_ele_action.H"
#include "scatra_ele_boundary_impl.H"

#include "scatra_ele_impl.H"
#include "scatra_ele_impl_reinit.H"

#include "../drt_lib/drt_globalproblem.H"  // for time curve in body force
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"
//#include <Teuchos_StandardParameterEntryValidators.hpp>  // included by inpar files
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"

// print error file for function EvaluateErrorComparedToAnalyticalSol()
#include "../drt_io/io_control.H"

#include "../drt_mat/diffcond.H"
#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"
#include "../drt_mat/newman.H"

//TODO(ehrl): Dirichlet for eliminated species

//#define DEBUG_BATTERY


/*----------------------------------------------------------------------*
  | calculate matrix and rhs for electrochemistry problem      gjb 10/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::DiffCondParams(
    DRT::Element*    ele,
    Teuchos::ParameterList&   params)
{
  // electrochemical formulation is evaluated based on the material
  Teuchos::RCP<MAT::Material> material = ele->Material();
  if (material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

    if(actmat->SpecById(actmat->SpecID(0))->MaterialType()==INPAR::MAT::m_diffcond or
       actmat->SpecById(actmat->SpecID(0))->MaterialType()==INPAR::MAT::m_newman)
    {
      diffcond_=true;

      Teuchos::ParameterList& diffcondparams = params.sublist("DIFFCOND");
      chemdiffcoupltransp_ = DRT::INPUT::IntegralValue<int>(diffcondparams,"CHEMICALDIFF_COUPLING_TRANSP");
      chemdiffcouplcurr_ = DRT::INPUT::IntegralValue<int>(diffcondparams,"CHEMICALDIFF_COUPLING_CURR");
      diffbased_ = DRT::INPUT::IntegralValue<int>(diffcondparams,"DIFFBASED");

      // current as solution variable: definition in parameter list and material have to be consistent
      cursolvar_ = actmat->Current();
      if(cursolvar_ != DRT::INPUT::IntegralValue<int>(diffcondparams,"CURRENT_SOLUTION_VAR"))
        dserror("current does not fit");

      equpot_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams,"EQUPOT");
      constparams_ = DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(diffcondparams,"CONSTPARAMS");

      frt_ = params.get<double>("frt");

      if(not DRT::INPUT::IntegralValue<int>(diffcondparams,"DIFFCOND_FORMULATION"))
        dserror("flag DIFFCOND_FORMULATION in the sublist DIFFCOND "
            "need to be activated in order to use the formulation");
    }
    else
      dserror("In the moment material elchmat supports only matrials diffcond and newman!!");

    if(actmat->SpecById(actmat->SpecID(0))->MaterialType()==INPAR::MAT::m_newman)
    {
      newman_=true;

      if(constparams_!=false)
        dserror("Newman only for variable parameter");

      if(actmat->NumSpec()!=1)
        dserror("Newman material is only valid for binary electrolytes!! \n"
            "Reduce the number of materials for ionic species to one.");

      if(chemdiffcoupltransp_==true)
        dserror("The coupling term is included in the concentration dependent diffusion coefficient!! \n"
            "Therefore the coupling term need to be turned off.");

      // this term is included in Newman formulation (experimentally it can be switched off)
      if(chemdiffcouplcurr_==false)
        dserror("The coupling term in current equation need to be activated for Newmann material");

      if(diffbased_==true)
        dserror("The coupling term in current equation need to be based on the transference number "
            "and the conductivity for Newmann material");

      if(equpot_==INPAR::ELCH::equpot_enc)
        dserror("ENC cannot be used for Newman materials");
    }
  }
  else
    diffcond_=false;


  return;
}

/*----------------------------------------------------------------------*
  |  get the material constants  (private)                      gjb 10/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetMaterialParamsDiffCond(
  Teuchos::RCP<MAT::Material> material
  )
{
  const Teuchos::RCP<const MAT::ElchMat>& actmat
        = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

  bool materialNewman = false;

  // access ionic species
  if (actmat->NumSpec() < numscal_) dserror("Not enough materials in MatList.");
  // acces phases
  if (actmat->NumPhase() != 1) dserror("In the moment a single phase is only allowed.");

  for (int k = 0;k<numscal_;++k)
  {
    const int specid = actmat->SpecID(k);
    Teuchos::RCP<const MAT::Material> singlemat = actmat->SpecById(specid);

    // constant physical properties: valence, diffusion coefficient and transference number
    // TODO(ehrl): Is it possible to combine both materials??
    if(singlemat->MaterialType() == INPAR::MAT::m_diffcond)
    {
      const MAT::DiffCond* actsinglemat = static_cast<const MAT::DiffCond*>(singlemat.get());
      diffus_[k] = actsinglemat->Diffusivity();
      valence_[k] = actsinglemat->Valence();
      trans_[k] = actsinglemat->Transference();

      diffusvalence_[k] = valence_[k]*diffus_[k];

      if(k==numscal_-1 and constparams_==false)
      {
        double sum=0.0;
        for(int ispec=0; ispec<numscal_;++ispec)
        {
          sum += valence_[ispec]*valence_[ispec]*diffus_[ispec]*conint_[ispec];
        }
        double denomin = sum*sum;

        for(int kk = 0;kk<numscal_;++kk)
        {
          trans_[kk] = (valence_[kk]*valence_[kk]*diffus_[kk]*conint_[kk])/sum;
          //test += trans_[kk];
          for(int ispec=0; ispec<numscal_;++ispec)
          {
            if(kk==ispec)
              ((transderiv_[kk])[ispec])=(valence_[kk]*valence_[kk]*diffus_[kk]*(sum-valence_[kk]*valence_[kk]*diffus_[kk]*conint_[kk]))/denomin;
            else
              ((transderiv_[kk])[ispec])=(-valence_[kk]*valence_[kk]*diffus_[kk]*conint_[kk]*valence_[ispec]*valence_[ispec]*diffus_[ispec])/denomin;
          }
        }
      }
    }
    // concentration depending properties
    else if (singlemat->MaterialType() == INPAR::MAT::m_newman)
    {
      materialNewman = true;
      const MAT::Newman* actsinglemat = static_cast<const MAT::Newman*>(singlemat.get());
      valence_[k] = actsinglemat->Valence();

      // diffusion coefficient
      diffus_[k] = actsinglemat->ComputeDiffusionCoefficient(conint_[k]);
      diffusderiv_[k] = actsinglemat->ComputeFirstDerivDiffCoeff(conint_[k]);

      // transference number
      trans_[k] = actsinglemat->ComputeTransferenceNumber(conint_[k]);
      (transderiv_[k])[k] = actsinglemat->ComputeFirstDerivTrans(conint_[k]);

      /// dilute solution theory (diffusion potential in current equation):
      ///    a          b
      ///   |--|  |----------|
      ///   z_1 + (z_2 - z_1) t_1
      /// ------------------------ (RT/F kappa 1/c_k grad c_k)
      ///      z_1 z_2
      ///     |________|
      ///         c
      a_ = actsinglemat->A();
      b_ = actsinglemat->B();
      c_ = actsinglemat->C();

      // diffusvalence
      diffusvalence_[k] = valence_[k]*diffus_[k];
    }

    // check whether there is negative (physical) diffusivity
    if (diffus_[k] < -EPS15) dserror("negative (physical) diffusivity of species %i",k);
  }  // end for-loop over species k

  // loop over single phases
  for (int iphase=0; iphase < actmat->NumPhase();++iphase)
  {
    const int phaseid = actmat->PhaseID(iphase);
    Teuchos::RCP<const MAT::Material> singlemat = actmat->PhaseById(phaseid);

    if(singlemat->MaterialType() == INPAR::MAT::m_elchphase)
    {
      const MAT::ElchPhase* actsinglemat = static_cast<const MAT::ElchPhase*>(singlemat.get());

      //TODO(ehrl): conductivity and first derivative can maximally depend on one concentration
      // since time curve is used as input routine
      cond_[iphase] = actsinglemat->ComputeConductivity(conint_[0]);
      condderiv_[iphase] = actsinglemat->ComputeFirstDerivCond(conint_[0]);
      eps_[iphase] = actsinglemat->Epsilon();
      tort_[iphase] = actsinglemat->Tortuosity();
      epstort_[iphase]=eps_[iphase]*tort_[iphase];

      if(constparams_==false and materialNewman == false)
      {
        cond_[0] = 0.0;
        const double faraday = INPAR::SCATRA::faraday_const;
        for(int ispec = 0;ispec<numscal_;++ispec)
        {
         cond_[0] += frt_*faraday*valence_[ispec]*valence_[ispec]*diffus_[ispec]*conint_[ispec];
         condderiv_[ispec]= frt_*faraday*valence_[ispec]*valence_[ispec]*diffus_[ispec];
        }
      }
    }
  }

#ifdef DEBUG_BATTERY
  for (int k = 0;k<numscal_;++k)
  {
    std::cout << "concentration " << k << ":   " << conint_[k] << std::endl;
    std::cout << "diffusion coefficient " << k << ":   " << diffus_[k] << std::endl;
    std::cout << "derivation diffusion coefficient " << k << ":   " << diffusderiv_[k] << std::endl;
    std::cout << "conductivity " << k << ":   " << cond_[k] << std::endl;
    std::cout << "derivation conductivity " << k << ":   " << condderiv_[k] << std::endl;
    std::cout << "transference number " << k << ":   " << trans_[k] << std::endl;
    std::cout << "derivation transference number " << k << ":   " << (transderiv_[k])[k] << std::endl;
  }
#endif

}


/*----------------------------------------------------------------------*
  | calculate matrix and rhs for electrochemistry problem      gjb 10/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatElchBat(
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const double                          frt,
  const double                          timefac,
  const double                          alphaF,
  const double                          fac,
  const enum INPAR::SCATRA::ScaTraType  scatratype
  )
{
  const double faraday = INPAR::SCATRA::faraday_const;

  // get gradient of electric potential at integration point
  gradpot_.Multiply(derxy_,epotnp_);

  // current density at Gauss point
  curint_.Multiply(ecurnp_,funct_);

#if 0
  // derivation of current density at Gauss point
  LINALG::Matrix<nsd_,nsd_> iderxy(true);
  iderxy.MultiplyNT(ecurnp_,derxy_);
  double divi = 0.0;
  for (int idim=0;idim<nsd_;++idim)
    divi += iderxy(idim,idim);
#endif

  // migration term (convective part without z_k D_k): -F/RT\grad{\Phi}\grad
  migconv_.MultiplyTN(1.0,derxy_,gradpot_);

  for (int idim = 0; idim <nsd_; ++idim)
  {
    //get vdiv at time n+1 for np_genalpha,
    LINALG::Matrix<nsd_,nsd_> curderxy(true);
    curderxy.MultiplyNT(ecurnp_,derxy_);
    curdiv_ += curderxy(idim, idim);
  }

  // loop over all transported scalars
  for (int k = 0; k < numscal_;++k)
    gradphicoupling_[k].Multiply(derxy_,ephinp_[k]);


  // time parameter

  double timefacfac   = 0.0;
  //double timetaufac   = 0.0;
  double rhsfac       = 0.0;
  //double rhstaufac    = 0.0;


  // perform time-integration specific actions
  if (is_stationary_)
  {
    // do not include any timefac for stationary calculations!
    timefacfac  = fac;
    rhsfac      = fac;
  }
  else
  {
    timefacfac  = timefac * fac;

    if (is_genalpha_)
    {
      rhsfac    = timefacfac/alphaF;
      //rhstaufac = timetaufac/alphaF;
    }
    else
    {
      rhsfac    = timefacfac;
    }
  } // if(is_genalpha_)

#if 0
  // DEBUG output
  std::cout<<std::endl<<"values at GP:"<<std::endl;
  std::cout<<"factor F/RT = "<<frt<<std::endl;
  for (int k=0;k<numscal_;++k)
  {std::cout<<"conint_["<<k<<"] = "<<conint_[k]<<std::endl;}
  for (int k=0;k<nsd_;++k)
  {std::cout<<"gradpot_["<<k<<"] = "<<gradpot_(k)<<std::endl;}
#endif

  // loop over all transported scalars
  for (int k = 0; k < numscal_;++k)
  {
    // compute gradient of scalar k at integration point
    gradphi_.Multiply(derxy_,ephinp_[k]);
    //std::cout << "grad_phi_k:  " <<gradphi_(k) << std::endl;

    // factor D_k * z_k
    const double diffus_valence_k = diffusvalence_[k];

    //double diff_ephinp_k(0.0);
    //double migrea_k(0.0);
    if (use2ndderiv_) // only necessary for higher order elements
    {
      diff_.Clear();
      migrea_.Clear();

      // diffusive part:  diffus_k * ( N,xx  +  N,yy +  N,zz )
      diff_.Update(diffus_[k],laplace_);

      // get Laplacian of electric potential at integration point
      double lappot = laplace_.Dot(epotnp_);
      // reactive part of migration term
      migrea_.Update(-frt*diffus_valence_k*lappot,funct_);

      //diff_ephinp_k = diff_.Dot(ephinp_[k]);   // diffusion
      //migrea_k      = migrea_.Dot(ephinp_[k]); // reactive part of migration term
    }
    else
    {
      diff_.Clear();
      migrea_.Clear();
    }

    // further short cuts and definitions
    const double conv_ephinp_k = conv_.Dot(ephinp_[k]);

    double rhsint       = rhs_[k]; // source/sink terms at int. point
    //double residual     = 0.0;

    // perform time-integration specific actions
    if (is_stationary_)
    {

    }
    else
    {
      if (is_genalpha_)
      {
        // note: in hist_ we receive the time derivative phidtam at time t_{n+alpha_M} !!
        //residual  = hist_[k] + conv_ephinp_k - diff_ephinp_k - rhsint;

        rhsint   *= (timefac/alphaF);  // not nice, but necessary !

        // rhs contribution due to incremental formulation (phidtam)
        // Standard Galerkin term
        const double vtrans = rhsfac*eps_[0]*hist_[k];
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vtrans*funct_(vi);
        }
      }
      else
      {
        rhsint = eps_[0]*hist_[k] + (rhs_[k]*timefac); // contributions from t_n and \theta*dt*bodyforce(t_{n+1})
        //residual  = conint_[k] + timefac*(conv_ephinp_k - diff_ephinp_k) - rhsint;

        // rhs contribution due to incremental formulation (phinp)
        // Standard Galerkin term
        const double vtrans = fac*conint_[k];
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vtrans*funct_(vi);
        }
      } // if(is_genalpha_)

    //----------------------------------------------------------------
    // 2) element matrix
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // 2a) element matrix: governing equation for concentration
    //----------------------------------------------------------------

  //------------------------------------------------------------
  // k ionic species
  //------------------------------------------------------------

      // instationary terms
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;
        const double fac_funct_vi = fac*funct_(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          /* Standard Galerkin term: */
          // (w, c)
          emat(fvi, fui) += fac_funct_vi*eps_[0]*funct_(ui) ;

        } // for ui
      } // for vi
    } // if (is_stationary_)
    // end time-integration specific actions

#ifdef PRINT_ELCH_DEBUG
    std::cout<<"tau["<<k<<"]    = "<<tau_[k]<<std::endl;
    std::cout<<"taufac["<<k<<"] = "<<taufac<<std::endl;
    if (tau_[k] != 0.0)
      std::cout<<"residual["<<k<<"] = "<< residual<<std::endl;
    std::cout<<"conv_eff_k    = "<<conv_eff_k<<std::endl;
    std::cout<<"conv_ephinp_k  = "<<conv_ephinp_k<<std::endl;
    std::cout<<"Dkzk_mig_ephinp_k = "<<Dkzk_mig_ephinp_k<<std::endl;
    std::cout<<"diff_ephinp_k = "<<diff_ephinp_k<<std::endl;
    std::cout<<"migrea_k      = "<<migrea_k <<std::endl;
    std::cout<<std::endl;
#endif

    //----------------------------------------------------------------------------------------------
    // stationary terms (transport equation)
    //----------------------------------------------------------------------------------------------

    // standard: convection and self-diffusion term
    for (int vi=0; vi<nen_; ++vi)
    {
      const int    fvi = vi*numdofpernode_+k;
      const double timefacfac_funct_vi = timefacfac*funct_(vi);

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        // matrix entries
        double matvalconc = 0.0;

        // convective term (transport equation)
        //
        // (w, u grad Dc)
        //
        matvalconc += timefacfac_funct_vi*eps_[0]*conv_(ui) ;

        // ionic diffusive term (transport equation)
        //
        // (grad w, D grad Dc)
        //
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi); // compute once, reuse below!
        matvalconc += timefacfac*epstort_[0]*diffus_[k]*laplawf;

        // try to access the element matrix not too often. Can be costly
        emat(fvi,fui)                        += matvalconc;
      } // for ui
    } // for vi

    //additional term for Newman material: concentration depending diffusion coefficient
    if(newman_==true)
    {
      //linearization of diffusion coefficient in the ionic diffusion term (transport equation)
      //
      // (grad w, D(D(c)) grad c)
      //
      if(constparams_==false)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            double laplawf=0.0;
            GetLaplacianWeakFormRHS(laplawf,derxy_,gradphicoupling_[k],vi);

            emat(vi*numdofpernode_+k,ui*numdofpernode_+k)
              += timefacfac*epstort_[0]*diffusderiv_[k]*laplawf*funct_(ui);
          }
        }
      }
    }

    //----------------------------------------------------------------------------------------------
    // Coupling term: diffusion potential (transport equation)
    //----------------------------------------------------------------------------------------------
    if(chemdiffcoupltransp_==true)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          //the coupling term cancel out for a 2 equation system
          //(current not a solution variable)
          if(cursolvar_==true)
          {
            // compute once, reuse below!
            double laplawf(0.0);
            GetLaplacianWeakForm(laplawf, derxy_,ui,vi);

            for (int iscal=0; iscal < numscal_; ++iscal)
            {
              // formulation a): plain ionic diffusion coefficients without using ENC
              //
              // (grad w, sum(D_i grad Dc))
              //
              emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                += -timefacfac*epstort_[0]*trans_[k]/valence_[k]*valence_[iscal]*diffus_[iscal]*laplawf;

              // formulation b):  plain ionic diffusion coefficients simplified by ENC

              if(constparams_==false)
              {
                //linearization of transference number in the coupling term (transport equation)
                //
                // (grad w, Dt_k(c)/z_k (Sum_i z_i D_i grad c_i))
                //
                for(int iscal2=0; iscal2<numscal_; ++iscal2)
                {
                  double laplawf2=0.0;
                  GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal2],vi);

                  emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                    += -timefacfac*epstort_[0]*((transderiv_[k])[iscal])*funct_(ui)*valence_[iscal2]/valence_[k]*diffus_[iscal2]*laplawf2;
                } // for(iscal2)
              } // if(constparams_)
            } // for(iscal)
          } // if(cursolvar_)
        } // for ui
      } // for vi
    } // end i(chemdiffcoupltransp)


    //----------------------------------------------------------------------------------------------
    // electrical conduction or migration term (transport equation)
    //----------------------------------------------------------------------------------------------

    if(cursolvar_ == true)
    {
      // current term (with current as a solution variable)
      for (int vi=0; vi<nen_; ++vi)
      {
        const int    fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          for (int idim=0; idim<nsd_; ++idim)
          {
            const int fui = ui*numdofpernode_+(numscal_+1)+idim;

            //linearization of conduction term (transport equation)
            //
            // (grad w, t_k/(F z_k) Di)
            //
            emat(fvi,fui) += -timefacfac*derxy_(idim,vi)*trans_[k]/(faraday*valence_[k])*funct_(ui);

            if(constparams_==false)
            {
              //linearization of transference number in conduction term (transport equation)
              //
              // (grad w, Dt_k(c)/(F z_k) i)
              //
              for (int iscal=0; iscal<numscal_; ++iscal)
                emat(fvi,ui*numdofpernode_+iscal)
                  += -timefacfac*derxy_(idim,vi)*((transderiv_[k])[iscal])*funct_(ui)/(faraday*valence_[k])*curint_(idim);
            }
          }
        }
      }
    }
    else
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          double laplawf(0.0);
          GetLaplacianWeakForm(laplawf, derxy_,ui,vi); // compute once, reuse below!

          //linearization of conduction term (transport equation)
          //
          // (grad w, t_k kappa/(F z_k) D(grad phi))
          //
          emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) +=
            timefacfac*epstort_[0]*trans_[k]*cond_[0]/faraday/valence_[k]*laplawf;

          if(constparams_==false)
          {
            double laplawf2(0.0);
            GetLaplacianWeakFormRHS(laplawf2,derxy_,gradpot_,vi);

            for(int iscal=0; iscal<numscal_;++iscal)
            {
              //Todo (ehrl): Convergence in the first two steps is quite bad (no current, ENC)

              //linearization of the conductivity in the conduction term (transport equation)
              //
              // (grad w, t_k D(kappa(c))/(F z_k) grad phi)
              //
              emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                 += timefacfac*epstort_[0]*trans_[k]/faraday/valence_[k]*condderiv_[iscal]*funct_(ui)*laplawf2;

              //linearization of the transference number in the conduction term (transport equation)
              //
              // (grad w, D(t_k(c)) kappa/(F z_k) grad phi)
              //
              emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                 += timefacfac*epstort_[0]*((transderiv_[k])[iscal])*funct_(ui)/faraday/valence_[k]*cond_[0]*laplawf2;
            }
          }

          // inconsitstence for ..
          if(newman_==true)
          {
            double laplawf(0.0);
            GetLaplacianWeakForm(laplawf, derxy_,ui,vi);

            for(int iscal=0;iscal<numscal_;++iscal)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                += timefacfac*epstort_[0]/frt_/faraday/valence_[k]*trans_[k]*cond_[0]*((a_+(b_*trans_[iscal]))/c_)/conint_[iscal]*laplawf;

              //TODO (ehrl): Linearization only for one species (otherwise you need two for-dim loops)
              {
                double laplawf2(0.0);
                GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal],vi);

                // Linearization wrt ln c
                emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                  += -timefacfac*epstort_[0]/frt_/faraday/valence_[k]*trans_[k]*cond_[0]*(a_+(b_*trans_[iscal]))/c_/conint_[iscal]/conint_[iscal]*laplawf2*funct_(ui);

                // Linearization wrt kappa
                emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                  += timefacfac*epstort_[0]/frt_/faraday/valence_[k]*trans_[k]*(a_+(b_*trans_[iscal]))/c_/conint_[iscal]*laplawf2*condderiv_[iscal]*funct_(ui);

                // Linearization wrt transference number
                emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                  += timefacfac*epstort_[0]/frt_/faraday/valence_[k]*cond_[0]/c_/conint_[iscal]*laplawf2*(a_+b_*trans_[iscal])*(transderiv_[iscal])[iscal]*funct_(ui);

                // Linearization wrt transference number
                emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                  += timefacfac*epstort_[0]/frt_/faraday/valence_[k]*cond_[0]/c_/conint_[iscal]*laplawf2*trans_[iscal]*b_*(transderiv_[iscal])[iscal]*funct_(ui);
              }

            }
          }

          if(diffbased_==false and equpot_==INPAR::ELCH::equpot_divi and constparams_==true)
          {
            double laplawf(0.0);
            GetLaplacianWeakForm(laplawf, derxy_,ui,vi);

            for(int iscal=0;iscal<numscal_;++iscal)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+iscal)
                += timefacfac*epstort_[0]/frt_/faraday/valence_[k]*trans_[k]*cond_[0]*trans_[iscal]/valence_[iscal]/conint_[iscal]*laplawf;
            }
          }
        } // for ui
      } // for vi
    } // end if(cursolvar_)

    //----------------------------------------------------------------------------------------------
    // equation for potential (electroneutrality)
    //----------------------------------------------------------------------------------------------

    if(equpot_==INPAR::ELCH::equpot_enc)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;
        const double alphaF_valence_k_fac_funct_vi = alphaF*valence_[k]*fac*funct_(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          //linearization of the transference number in the conduction term (transport equation)
          //
          // (w, sum(z_k c_k))
          //
          emat(pvi, fui) += alphaF_valence_k_fac_funct_vi*funct_(ui);
        } // for ui
      }
    }


    //-----------------------------------------------------------------------
    // 3) element right hand side vector (neg. residual of nonlinear problem)
    //-----------------------------------------------------------------------

    //----------------------------------------------------------------
    // 3a) rhs: governing equation for concentration
    //----------------------------------------------------------------
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      // RHS source term (contains old part of rhs for OST / BDF2)
      erhs[fvi] += fac*eps_[0]*funct_(vi)*rhsint ;

      // convective term
      erhs[fvi] -= rhsfac*eps_[0]*funct_(vi)*conv_ephinp_k;

      {
        // diffusive term
        double laplawf(0.0);
        GetLaplacianWeakFormRHS(laplawf,derxy_,gradphicoupling_[k],vi);
        erhs[fvi] -= rhsfac*epstort_[0]*diffus_[k]*laplawf;
      }

      if(chemdiffcoupltransp_==true)
      {
        for (int iscal=0; iscal < numscal_; ++iscal)
        {
          // diffusive term second
          double laplawf(0.0);
          GetLaplacianWeakFormRHS(laplawf,derxy_,gradphicoupling_[iscal],vi); // compute once, reuse below!

          //this term cancel out if current is not a solution variable
          if(cursolvar_==true)
          {
            erhs[fvi] -= - rhsfac*epstort_[0]*trans_[k]*valence_[iscal]/valence_[k]*diffus_[iscal]*laplawf;
          }
        }
      }

      if(cursolvar_ == true)
      // current density term
      {
        double laplawf=0.0;
        GetLaplacianWeakFormRHS(laplawf,derxy_,curint_,vi);
        erhs[fvi]-= -rhsfac*trans_[k]/faraday/valence_[k]*laplawf;

        //TODO(ehrl): concentration dependence of transference number
      }
      else
      {
        // diffusive term
        double laplawf=0.0;
        GetLaplacianWeakFormRHS(laplawf,derxy_,gradpot_,vi);
        erhs[fvi]-= rhsfac*epstort_[0]*trans_[k]*cond_[0]/faraday/valence_[k]*laplawf;

        if(newman_==true)
        {
          for(int iscal=0; iscal<numscal_;++iscal)
          {
            // diffusive term second
            double laplawf2(0.0);
            GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal],vi); // compute once, reuse below!

            erhs[fvi]
              -= rhsfac*epstort_[0]/frt_/faraday/valence_[k]*trans_[k]*cond_[0]*((a_+(b_*trans_[iscal]))/c_)/conint_[iscal]*laplawf2;
          }
        }

        if(diffbased_==false and equpot_==INPAR::ELCH::equpot_divi and constparams_==true)
        {
          for(int iscal=0; iscal<numscal_;++iscal)
          {
            // diffusive term second
            double laplawf2(0.0);
            GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal],vi); // compute once, reuse below!

            erhs[fvi]
              -= rhsfac*epstort_[0]/frt_/faraday/valence_[k]*trans_[k]*cond_[0]*trans_[iscal]/valence_[iscal]/conint_[iscal]*laplawf2;
           }
         }
      }
    } // for vi

    //----------------------------------------------------------------
    // 3c) rhs: governing equation for potential (electroneutrality)
    //----------------------------------------------------------------

    if(equpot_==INPAR::ELCH::equpot_enc)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        // electroneutrality condition
        // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
        erhs[vi*numdofpernode_+numscal_] -= valence_[k]*fac*funct_(vi)*conint_[k];
      } // for vi
    // RHS vector finished
    }
  } // loop over scalars


//-------------------------------------------------------------------------
// governing equation for current
//-------------------------------------------------------------------------
  if(cursolvar_==true)
  {
    // (v, i)
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int ui=0; ui<nen_; ++ui)
      {
        for (int idim=0; idim<nsd_; ++idim)
        {
          const int fvi = vi*numdofpernode_+(numscal_+1)+idim;
          const int fui = ui*numdofpernode_+(numscal_+1)+idim;

          emat(fvi,fui) += timefacfac/faraday*funct_(vi)*funct_(ui);
        }
      }
    }

    // (v, kappa grad phi)
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int ui=0; ui<nen_; ++ui)
      {
        for (int idim=0; idim<nsd_; ++idim)
        {
          const int fvi = vi*numdofpernode_+(numscal_+1)+idim;
          const int fui = ui*numdofpernode_+numscal_;

          emat(fvi,fui) += timefacfac/faraday*epstort_[0]*funct_(vi)*cond_[0]*derxy_(idim,ui);

          //linearization of conductivity in the ohmic resistance term (current equation)
          //
          // (w, D(kappa(c)) grad phi)
          //
          if(constparams_==false)
          {
            for(int k=0;k<numscal_;++k)
              emat(fvi,ui*numdofpernode_+k) += timefacfac/faraday*epstort_[0]*funct_(vi)*condderiv_[k]*funct_(ui)*gradpot_(idim);
          }
        }
      }
    }

    // Coupling term
    if(chemdiffcouplcurr_==true)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          // diffusive term
          // (grad w, D grad c)
          for (int idim = 0; idim < nsd_; ++idim)
          {
            for (int iscal=0; iscal < numscal_; ++iscal)
            {
              if(newman_==false)
              {
                if(diffbased_==true)
                  emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                    += timefacfac/faraday*epstort_[0]*funct_(vi)*faraday*valence_[iscal]*diffus_[iscal]*derxy_(idim,ui);
                else
                  emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                    += timefacfac/faraday*epstort_[0]*funct_(vi)/frt_*cond_[0]*trans_[iscal]/valence_[iscal]/conint_[iscal]*derxy_(idim,ui);
              }
              else
              {
                emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                  += timefacfac/faraday*epstort_[0]*funct_(vi)/frt_*cond_[0]*((a_+(b_*trans_[iscal]))/c_)/conint_[iscal]*derxy_(idim,ui);

// linearization of coupling term in the current equation??
#if 0
                emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                  += timefacfac*funct_(vi)/frt_*cond_[0]*trans_[iscal]/valence_[iscal]/((-1)*conint_[iscal]*conint_[iscal])*funct_(ui)*gradphicoupling_[iscal](idim);

                if(constparams_==false)
                {
                  for(int iscal2=0; iscal2<numscal_;++iscal2)
                  {
                    //Check if necessary??
                    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                      += timefacfac*funct_(vi)/frt_*condderiv_[iscal]*funct_(ui)*trans_[iscal2]*funct_(ui)/valence_[iscal2]/conint_[iscal2]*gradphicoupling_[iscal2](idim);


                    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                      += timefacfac*funct_(vi)/frt_*cond_[0]*(transderiv_[iscal2])[iscal]*funct_(ui)/valence_[iscal2]/conint_[iscal2]*gradphicoupling_[iscal2](idim);
                  }
                }
#endif
              }
            }
          }
        } // for ui
      } // for vi
    }

    //-------------------------------------------------------------------------
    // 2c) element matrix: governing equation for potential (div i)
    //-------------------------------------------------------------------------
    if(equpot_==INPAR::ELCH::equpot_divi)
    {
      for (int vi=0; vi<nen_; ++vi)
       {
         for (int ui=0; ui<nen_; ++ui)
         {
           for (int idim = 0; idim <nsd_; ++idim)
           {
             const int fvi = numdofpernode_*vi+numscal_;
             const int fui = numdofpernode_*ui+(numscal_+1)+idim;
             /* current continuity term */
             /*
                  /               \
                 |                 |
                 | w, nabla o Di   |
                 |                 |
                  \               /
             */
             //emat(fvi,fui) += timefacfac*funct_(vi);*derxy_(idim,ui);

             /* current continuity term */
             /*
                  /               \
                 |                 |
                 | grad phi,  Di   |
                 |                 |
                  \               /
             */
             // version a: (grad phi,  Di)
             emat(fvi,fui) -= timefacfac/faraday*derxy_(idim,vi)*funct_(ui);
             // version b: (phi, div Di) -> not partially integrated
             //emat(fvi,fui) += timefacfac*funct_(vi)*derxy_(idim,ui);
           } // end for(idim)
         } // end for(ui)
       } // end for(vi)
    }
  }
  // if(cursolvar_==true)
  else
  {
    if(equpot_==INPAR::ELCH::equpot_divi)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          // diffusive term
          // (grad w, D grad c)
          double laplawf(0.0);
          GetLaplacianWeakForm(laplawf, derxy_,ui,vi); // compute once, reuse below!

          emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_) += timefacfac*epstort_[0]/faraday*cond_[0]*laplawf;

          if(constparams_==false)
          {
            double laplawf2(0.0);
            GetLaplacianWeakFormRHS(laplawf2,derxy_,gradpot_,vi);

            for(int iscal=0;iscal<numscal_;++iscal)
            {
              emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                += timefacfac*epstort_[0]/faraday*condderiv_[iscal]*funct_(ui)*laplawf2;
            }
          }

          for (int iscal=0; iscal < numscal_; ++iscal)
          {
            if(newman_==false)
            {
              if(diffbased_==true)
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                  += timefacfac/faraday*epstort_[0]*faraday*valence_[iscal]*diffus_[iscal]*laplawf;
              else
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                  += timefacfac/faraday*epstort_[0]/frt_*cond_[0]*trans_[iscal]/valence_[iscal]/conint_[iscal]*laplawf;
            }
            else
            {
              emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                += timefacfac*epstort_[0]/faraday/frt_*cond_[0]*(a_+(b_*trans_[iscal]))/c_/conint_[iscal]*laplawf;

              //TODO (ehrl): Linearization only for one species (otherwise you need two for-dim loops)
              {
                double laplawf2(0.0);
                GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal],vi);

                // Linearization wrt ln c
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                  += -timefacfac*epstort_[0]/faraday/frt_*cond_[0]*(a_+(b_*trans_[iscal]))/c_/conint_[iscal]/conint_[iscal]*laplawf2*funct_(ui);

                // Linearization wrt kappa
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                  += timefacfac*epstort_[0]/faraday/frt_*(a_+(b_*trans_[iscal]))/c_/conint_[iscal]*laplawf2*condderiv_[iscal]*funct_(ui);

                // Linearization wrt transference number
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+iscal)
                  += timefacfac*epstort_[0]/faraday/frt_*cond_[0]/c_/conint_[iscal]*laplawf2*b_*(transderiv_[iscal])[iscal]*funct_(ui);
              }
            }
          }
        } // for ui
      } // for vi

      // rhs of
      for (int vi=0; vi<nen_; ++vi)
      {
        // diffusive term
        // (grad w, D grad c)
        {
          double laplawf(0.0);
          GetLaplacianWeakFormRHS(laplawf, derxy_,gradpot_,vi); // compute once, reuse below!

          erhs[vi*numdofpernode_+numscal_] -= rhsfac*epstort_[0]/faraday*cond_[0]*laplawf;
        }

        for (int iscal=0; iscal < numscal_; ++iscal)
        {
          if(newman_== false)
          {
            // diffusive term second
            double laplawf2(0.0);
            GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal],vi); // compute once, reuse below!

            if(diffbased_==true)
              erhs[vi*numdofpernode_+numscal_] -= rhsfac/faraday*epstort_[0]*faraday*valence_[iscal]*diffus_[iscal]*laplawf2;
            else
              erhs[vi*numdofpernode_+numscal_] -= rhsfac/faraday*epstort_[0]/frt_*cond_[0]*trans_[iscal]/valence_[iscal]/conint_[iscal]*laplawf2;
          }
          else
          {
            // diffusive term second
            double laplawf2(0.0);
            GetLaplacianWeakFormRHS(laplawf2,derxy_,gradphicoupling_[iscal],vi); // compute once, reuse below!

            erhs[vi*numdofpernode_+numscal_]
              -= rhsfac*epstort_[0]/faraday/frt_*cond_[0]*((a_+(b_*trans_[iscal]))/c_)/conint_[iscal]*laplawf2;
          }
        }
      }
    }
  }


  //----------------------------------------------------------------
  // 2d) Stabilization terms
  //----------------------------------------------------------------

  // no stabilization terms

  //----------------------------------------------------------------
  // 3b) rhs: governing equation for current
  //----------------------------------------------------------------
  if(cursolvar_ == true)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        // (v, i)
        erhs[vi*numdofpernode_+(numscal_+1)+idim]-=rhsfac/faraday*funct_(vi)*curint_(idim);

        // (v, kappa grad phi)
        erhs[vi*numdofpernode_+(numscal_+1)+idim]-=rhsfac/faraday*epstort_[0]*funct_(vi)*cond_[0]*gradpot_(idim);
      }

      if(equpot_==INPAR::ELCH::equpot_divi)
      {
        // (v, i)
        //erhs[vi*numdofpernode_+(numdofpernode_-1)]-=rhsfac*funct_(vi)*curdiv_;
        {
          double laplawf=0.0;
          // version a: (grad phi,  Di)
          GetLaplacianWeakFormRHS(laplawf,derxy_,curint_,vi);
          erhs[vi*numdofpernode_+numscal_]-= -rhsfac/faraday*laplawf;
          // version b: (phi, div Di) -> not partially integrated
          //erhs[vi*numdofpernode_+numscal_] -= rhsfac*funct_(vi)*divi;
        }
      }

      // Coupling term
      if(chemdiffcouplcurr_==true)
      {
        // diffusive term
        // (grad w, D grad c)
        for (int idim = 0; idim < nsd_; ++idim)
        {
          for (int iscal=0; iscal < numscal_; ++iscal)
          {
            if(newman_== false)
            {
              if(diffbased_==true)
                erhs[vi*numdofpernode_+(numscal_+1)+idim] -= rhsfac/faraday*epstort_[0]*funct_(vi)*faraday*valence_[iscal]*diffus_[iscal]*gradphicoupling_[iscal](idim);
              else
              {
                erhs[vi*numdofpernode_+(numscal_+1)+idim]
                     -= rhsfac/faraday*epstort_[0]*funct_(vi)/frt_*cond_[0]*trans_[iscal]/valence_[iscal]/conint_[iscal]*gradphicoupling_[iscal](idim);
              }
            }
            else
            {
              erhs[vi*numdofpernode_+(numscal_+1)+idim]
                -= rhsfac/faraday*epstort_[0]*funct_(vi)/frt_*cond_[0]*((a_+(b_*trans_[iscal]))/c_)/conint_[iscal]*gradphicoupling_[iscal](idim);
              //erhs[vi*numdofpernode_+(numscal_+1)+idim]
              //  -= rhsfac*funct_(vi)*faraday*(2.0e-5-4.0e-5)*gradphicoupling_[iscal](idim);
            }
          }
        }
      }
    }
  }

  //----------------------------------------------------------------
  // 3c) rhs: governing equation for potential
  //----------------------------------------------------------------
  //if(eid_==1)
    //PrintEleMatToExcel(emat,erhs);
  //---------------------------------------------------------------
  // 3d) rhs: Stabilization terms
  //----------------------------------------------------------------

  return;
} // ScaTraImpl::CalMatElch

/*----------------------------------------------------------------------*
  |  Calculate conductivity (ELCH) (private)                   gjb 07/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::PrintEleMatToExcel(
    Epetra_SerialDenseMatrix&             emat,
    Epetra_SerialDenseVector&             erhs
  )
{
  std::ostringstream oss;
  oss << eid_;
  std::string value = oss.str();
  //ostringstream temp;

  {
    const std::string output = "Matrix"+value+".csv";

    std::ofstream f;
    f.open(output.c_str());

    int nodeh = 0;
    int nodev = 0;
    for (int ii=0; ii<nen_*numdofpernode_; ++ii)
    {
     for (int jj=0; jj<nen_*numdofpernode_; ++jj)
     {
         f << emat(ii,jj)<< ",";
         if(jj== numdofpernode_*nodeh+numdofpernode_-1)
         {
           f << ",";
           nodeh +=1;
         }
         if(jj== nen_*numdofpernode_-1)
           f << std::endl;
     } // end for(ui)
     nodeh =0;
     if(ii== numdofpernode_*nodev+numdofpernode_-1)
     {
       f << std::endl;
       nodev +=1;
     }
    } // end for(vi)

    f.flush();
    f.close();
  }

  {
    //ostringstream temp;
    const std::string output = "residual"+value+".csv";

    std::ofstream f;
    f.open(output.c_str());
    int nodev = 0;
    for (int ii=0; ii<nen_*numdofpernode_; ++ii)
    {
     f << erhs(ii)<< std::endl;

     if(ii== numdofpernode_*nodev+numdofpernode_-1)
     {
       f << std::endl;
       nodev +=1;
     }
    } // end for(vi)

    f.flush();
    f.close();
  }
  if(eid_==9)
    //dserror("element matrix was printed to a csv file");
  return;
}


template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tet4>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tet10>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::ScaTraImpl<DRT::Element::nurbs27>;

