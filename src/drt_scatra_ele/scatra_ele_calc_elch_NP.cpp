/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.cpp

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_elch_NP.H"
#include "scatra_ele.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/ion.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"
#include "../drt_inpar/inpar_elch.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchNP<distype> * DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcElchNP<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcElchNP<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::ScaTraEleCalcElchNP(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalcElch<distype>::ScaTraEleCalcElch(numdofpernode,numscal)
{
  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs                           ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalMatAndRhs(
  Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
  Epetra_SerialDenseMatrix&               emat,         //!< element matrix to calculate
  Epetra_SerialDenseVector&               erhs,         //!< element rhs to calculate+
  const int                               k,            //!< index of current scalar
  const double                            fac,          //!< domain-integration factor
  const double                            timefacfac,   //!< domain-integration factor times time-integration factor
  const double                            rhsfac,       //!< time-integration factor for rhs times domain-integration factor
  Teuchos::RCP<ScaTraEleDiffManagerElch>& dme,          //!< diffusion manager
  double&                                 rhsint,       //!< rhs at Gauss point
  const double                            hist          //!< history
  )
{
#if 0
  //----------------------------------------------------------------
  // 1) element matrix: instationary terms
  //----------------------------------------------------------------
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;
    const double fac_funct_vi = fac*my::funct_(vi);

//          // compute effective convective stabilization operator
//          double conv_eff_vi = conv_(vi);
//          if (migrationstab_)
//          {
//            conv_eff_vi += dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv_(vi);
//          }

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;

      /* Standard Galerkin term: */
      emat(fvi, fui) += fac_funct_vi*my::funct_(ui) ;

//            /* 1) convective stabilization of transient term*/
//            emat(fvi, fui) += taufac*conv_eff_vi*my::funct_(ui);

      /* 2) diffusive stabilization */
      // not implemented. Only stabilization of SUPG type

      /* 3) reactive stabilization (reactive part of migration term) */
      // not implemented. Only stabilization of SUPG type

    } // for ui
  } // for vi

  //----------------------------------------------------------------
  // 2) element matrix: stationary terms
  //----------------------------------------------------------------
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int    fvi = vi*my::numdofpernode_+k;

    // compute effective convective stabilization operator
    double conv_eff_vi = conv(vi);
    if (migrationstab_)
    {
      conv_eff_vi += dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi);
    }

    const double timefacfac_funct_vi = timefacfac*my::funct_(vi);
    const double timefacfac_diffus_valence_k_mig_vi = timefacfac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi);

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      // matrix entries
      double matvalconc = 0.0;
      double matvalpot = 0.0;

      // convective term
      matvalconc += timefacfac_funct_vi*conv(ui) ;

      // addition to convective term for conservative form
      if (my::scatrapara_->IsConservative())
      {
        // convective term using current scalar value
        matvalconc += timefacfac_funct_vi*vdiv*my::funct_(ui);
      }

      // diffusive term
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf,ui,vi); // compute once, reuse below!
      matvalconc += timefacfac*dme->GetIsotropicDiff(k)*laplawf;

      // migration term
      // a) derivative w.r.t. concentration c_k
      matvalconc -= timefacfac_diffus_valence_k_mig_vi*my::funct_(ui);
      // b) derivative w.r.t. electric potential
      matvalpot += frt*timefacfac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k]*laplawf;

      // TODO (ehrl)
      // Including stabilization yields in different results for the uncharged particle and
      // the binary electrolyte solution
      // -> Check calculation procedure of the method

      //----------------------------------------------------------------
      // Stabilization terms
      //----------------------------------------------------------------

      /* 0) transient stabilization */
      // not implemented. Only stabilization of SUPG type

      /* 1) convective stabilization */

      /* convective term */

      // I) linearization of residual part of stabilization term

      // effective convective stabilization of convective term
      // derivative of convective term in residual w.r.t. concentration c_k
      matvalconc += timetaufac*conv_eff_vi*conv(ui);

      // migration convective stabilization of convective term
      double val_ui;
      my::GetLaplacianWeakFormRHS(val_ui,gradphi[k],ui);
      if (migrationinresidual_)
      {
        // a) derivative w.r.t. concentration_k
        matvalconc += timetaufac*conv_eff_vi*dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(ui);

        // b) derivative w.r.t. electric potential
        matvalpot -= timetaufac*conv_eff_vi*dme->GetIsotropicDiff(k)*dme->GetValence(k)*frt*val_ui;

        // note: higher-order and instationary parts of residuum part are linearized elsewhere!
      }

      // II) linearization of convective stabilization operator part of stabilization term
      if (migrationstab_)
      {
        // a) derivative w.r.t. concentration_k
        //    not necessary -> zero

        // b) derivative w.r.t. electric potential
        double laplacewf(0.0);
        my::GetLaplacianWeakForm(laplacewf,ui,vi);
        matvalpot -= timetaufac*residual*dme->GetIsotropicDiff(k)*dme->GetValence(k)*frt*laplacewf;
      }

      // III) linearization of tau part of stabilization term
      if (migrationintau_)
      {
        // derivative of tau (only effective for Taylor_Hughes_Zarins) w.r.t. electric potential
        const double tauderiv_ui = ((tauderpot[k])(ui,0));
        matvalpot += timefacfac*tauderiv_ui*conv_eff_vi*residual;
      }

      // try to access the element matrix not too often. Can be costly
      emat(fvi,fui)                        += matvalconc;
      emat(fvi,ui*my::numdofpernode_+my::numscal_) += matvalpot;

    } // for ui
  } // for vi

  //-------------------------------------------------------------------------
  // 2b) element matrix: stationary terms (governing equation for potential)
  //-------------------------------------------------------------------------
  // what's the governing equation for the electric potential field?
  // we provide a lot of different options here:
  switch (elchpara_->ElchType())
  {
  case INPAR::ELCH::elchtype_enc:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;
      const double alphaF_valence_k_fac_funct_vi = my::scatraparatimint_->AlphaF()*dme->GetValence(k)*fac*my::funct_(vi);

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;

        // electroneutrality condition (only derivative w.r.t. concentration c_k)
        emat(pvi, fui) += alphaF_valence_k_fac_funct_vi*my::funct_(ui);
      } // for ui
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_enc_pde:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;
      const double timefacfac_diffus_valence_k_mig_vi = timefacfac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi);

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;

        double laplawf(0.0);
        my::GetLaplacianWeakForm(laplawf,ui,vi);

        // use 2nd order pde derived from electroneutrality condition (k=1,...,m)
        // a) derivative w.r.t. concentration c_k
        emat(pvi, fui) -= dme->GetValence(k)*(timefacfac_diffus_valence_k_mig_vi*my::funct_(ui));
        emat(pvi, fui) += dme->GetValence(k)*(timefacfac*dme->GetIsotropicDiff(k)*laplawf);
        // b) derivative w.r.t. electric potential
        emat(pvi, ui*my::numdofpernode_+my::numscal_) += dme->GetValence(k)*(frt*timefacfac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k]*laplawf);
      } // for ui
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_enc_pde_elim:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;
      const double timefacfac_diffus_valence_k_mig_vi = timefacfac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi);
      const double timefacfac_diffus_valence_m_mig_vi = timefacfac*dme->GetIsotropicDiff(my::numscal_)*dme->GetValence(my::numscal_)*migconv(vi);

      for (int ui=0; ui<my::nen_; ++ui)
      {
        // matrix entries
        double matvalconc = 0.0;
        double matvalpot = 0.0;

        double laplawf(0.0);
        my::GetLaplacianWeakForm(laplawf,ui,vi);

        // use 2nd order pde derived from electroneutrality condition (k=1,...,m-1)
        // a) derivative w.r.t. concentration c_k
        matvalconc -= (timefacfac_diffus_valence_k_mig_vi*my::funct_(ui));
        matvalconc += (timefacfac*dme->GetIsotropicDiff(k)*laplawf);
        // b) derivative w.r.t. electric potential
        matvalpot += (frt*timefacfac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k]*laplawf);

        // care for eliminated species with index m
        //(diffus_ and valence_ vector were extended in GetMaterialParams()!)
        // a) derivative w.r.t. concentration c_k
        matvalconc += (timefacfac_diffus_valence_m_mig_vi*my::funct_(ui));
        matvalconc -= (timefacfac*dme->GetIsotropicDiff(my::numscal_)*laplawf);
        // b) derivative w.r.t. electric potential
        matvalpot -= (frt*timefacfac*dme->GetIsotropicDiff(my::numscal_)*dme->GetValence(my::numscal_)*conint[k]*laplawf);

        // try to access the element matrix not too often. Can be costly
        const int fui = ui*my::numdofpernode_+k;
        emat(pvi,fui) += dme->GetValence(k)*matvalconc;
        const int pui = ui*my::numdofpernode_+my::numscal_;
        emat(pvi,pui) += dme->GetValence(k)*matvalpot;

      } // for ui
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_poisson:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;
      const double alphaF_valence_k_fac_funct_vi = my::scatraparatimint_->AlphaF()*dme->GetValence(k)*fac*my::funct_(vi);

      for (int ui=0; ui<my::nen_; ++ui)
      {
        // we have a loop over k around. So prevent that the potential
        // term is added more than once!!
        if (k==0)
        {
          const int pui = ui*my::numdofpernode_+my::numscal_;
          double laplawf(0.0);
          my::GetLaplacianWeakForm(laplawf,ui,vi);

          const double epsbyF = epsilon/faraday;
          emat(pvi,pui) += my::scatraparatimint_->AlphaF()*fac*epsbyF*laplawf;
        }
        const int fui = ui*my::numdofpernode_+k;
        // electroneutrality condition (only derivative w.r.t. concentration c_k)
        emat(pvi, fui) += alphaF_valence_k_fac_funct_vi*my::funct_(ui);
      } // for ui
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_laplace:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        // we have a loop over k around. So prevent that the potential
        // term is added more than once!!
        if (k==0)
        {
          const int pui = ui*my::numdofpernode_+my::numscal_;
          double laplawf(0.0);
          my::GetLaplacianWeakForm(laplawf,ui,vi);
          emat(pvi,pui) += my::scatraparatimint_->AlphaF()*fac*laplawf;
        }
      } // for ui
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_diffcond:
    break;
  default:
  {
    dserror ("How did you reach this point?");
    break;
  }
  } // end switch(elchparam_->ElchType())


  if (my::use2ndderiv_)
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+k;

      // compute effective convective stabilization operator
      double conv_eff_vi = conv(vi);
      if (migrationstab_)
      {
        conv_eff_vi += dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi);
      }

      const double timetaufac_conv_eff_vi = timetaufac*conv_eff_vi;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;

        // 1) convective stabilization

        // diffusive term
        // derivative w.r.t. concentration c_k
        emat(fvi,fui) -= timetaufac_conv_eff_vi*diff(ui) ;

      } // for ui

      // reactive part of migration term
      if (migrationinresidual_)
      {
        const double timetaufac_conv_eff_vi_conint_k_frt_valence_k =timetaufac_conv_eff_vi*conint[k]*frt*dme->GetValence(k);
        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+k;

          // a) derivative w.r.t. concentration_k
          emat(fvi,fui) += timetaufac_conv_eff_vi*migrea(ui) ;
          // note: migrea_ already contains frt*diffus_valence!!!

          // b) derivative w.r.t. electric potential
          emat(fvi, ui*my::numdofpernode_+my::numscal_) -= timetaufac_conv_eff_vi_conint_k_frt_valence_k*diff(ui);
          // note: diff_ already includes factor D_k

        } // for ui
      }

      // 2) diffusive stabilization
      // not implemented. Only stabilization of SUPG type

      // 3) reactive stabilization (reactive part of migration term)
      // not implemented. Only stabilization of SUPG type

    } // for vi
  } // use2ndderiv

  //-----------------------------------------------------------------------
  // 3) element right hand side vector (neg. residual of nonlinear problem)
  //-----------------------------------------------------------------------
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    //----------------------------------------------------------------
    // standard Galerkin terms (ion transport equations)
    //----------------------------------------------------------------

    // RHS source term (contains old part of rhs for OST / BDF2)
    erhs[fvi] += fac*my::funct_(vi)*rhsint ;

    // nonlinear migration term
    erhs[fvi] += rhsfac*conint[k]*dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi);

    // convective term
    erhs[fvi] -= rhsfac*my::funct_(vi)*conv_ephinp_k;

    // addition to convective term for conservative form
    // (not included in residual)
    if (my::scatrapara_->IsConservative())
    {
      // convective term in conservative form
      erhs[fvi] -= rhsfac*my::funct_(vi)*conint[k]*vdiv;
    }

    // diffusive term
    double laplawf(0.0);
    my::GetLaplacianWeakFormRHS(laplawf,gradphi[k],vi);
    erhs[fvi] -= rhsfac*dme->GetIsotropicDiff(k)*laplawf;


    //----------------------------------------------------------------
    // Stabilization terms
    //----------------------------------------------------------------

    // 0) transient stabilization
    //    not implemented. Only stabilization of SUPG type

    // 1) convective stabilization

    erhs[fvi] -= rhstaufac*conv(vi)*residual;
    if (migrationstab_)
    {
      erhs[fvi] -=  rhstaufac*dme->GetIsotropicDiff(k)*dme->GetValence(k)*migconv(vi)*residual;
    }

    // 2) diffusive stabilization
    //    not implemented. Only stabilization of SUPG type

    // 3) reactive stabilization (reactive part of migration term)
    //    not implemented. Only stabilization of SUPG type
  } // for vi

    //----------------------------------------------------------------
    // standard Galerkin terms (equation for electric potential)
    //----------------------------------------------------------------
    // what's the governing equation for the electric potential field ?
  switch (elchpara_->ElchType())
  {
  case INPAR::ELCH::elchtype_enc:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;

      // electroneutrality condition
      // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
      erhs[pvi] -= dme->GetValence(k)*fac*my::funct_(vi)*conint[k];
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_enc_pde:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;

      double laplawf(0.0);
      my::GetLaplacianWeakFormRHS(laplawf,gradphi[k],vi);

      // use 2nd order pde derived from electroneutrality condition (k=1,...,m)
      erhs[pvi] += rhsfac*dme->GetValence(k)*((dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k]*migconv(vi))-(dme->GetIsotropicDiff(k)*laplawf));
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_enc_pde_elim:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;

      double laplawf(0.0);
      my::GetLaplacianWeakFormRHS(laplawf,gradphi[k],vi);

      // use 2nd order pde derived from electroneutrality condition (k=0,...,m-1)
      erhs[pvi] += rhsfac*dme->GetValence(k)*((dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k]*migconv(vi))-(dme->GetIsotropicDiff(k)*laplawf));
      // care for eliminated species with index m
      //(diffus_ and valence_ vector were extended in GetMaterialParams()!)
      erhs[pvi] -= rhsfac*dme->GetValence(k)*((dme->GetIsotropicDiff(my::numscal_)*dme->GetValence(my::numscal_)*conint[k]*migconv(vi))-(dme->GetIsotropicDiff(my::numscal_)*laplawf));
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_poisson:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;

      // we have a loop over k around. So prevent that the potential
      // term is added more than once!!
      if (k==0)
      {
        double laplawf(0.0);
        my::GetLaplacianWeakFormRHS(laplawf,gradpot,vi);
        const double epsbyF = epsilon/faraday;
        erhs[pvi] -= fac*epsbyF*laplawf;
      }

      // electroneutrality condition
      // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
      erhs[pvi] -= dme->GetValence(k)*fac*my::funct_(vi)*conint[k];
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_laplace:
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int pvi = vi*my::numdofpernode_+my::numscal_;

      // we have a loop over k around. So prevent that the potential
      // term is added more than once!!
      if (k==0)
      {
        double laplawf(0.0);
        my::GetLaplacianWeakFormRHS(laplawf,gradpot,vi);
        erhs[pvi] -= fac*laplawf;
      }
    } // for vi
    break;
  }
  case INPAR::ELCH::elchtype_diffcond:
    break;
  default:
  {
    dserror ("How did you reach this point?");
    break;
  }
  } // end switch (elchparam_->ElchType())
#endif
  return;
}

/*----------------------------------------------------------------------*
|  Correct sysmat for fluxes accros DC                       ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CorrectionForFluxAccrosDC(
  DRT::Discretization&        discretization,
  const std::vector<int>&     lm,
  Epetra_SerialDenseMatrix&   emat,
  Epetra_SerialDenseVector&   erhs)
{
#if 0
  if((scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim) or
     (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde))
  {
    double val(0.0);
    const DRT::Node* const* nodes = ele->Nodes();
    std::string condname = "Dirichlet";

    for (int vi=0; vi<nen_; ++vi)
    {
      std::vector<DRT::Condition*> dirichcond0;
      nodes[vi]->GetCondition(condname,dirichcond0);

      // there is at least one Dirichlet condition on this node
      if (dirichcond0.size() > 0)
      {
        //std::cout<<"Ele Id = "<<ele->Id()<<"  Found one Dirichlet node for vi="<<vi<<std::endl;
        const std::vector<int>*    onoff = dirichcond0[0]->Get<std::vector<int> >   ("onoff");
        for (int k=0; k<numscal_; ++k)
        {
          if ((*onoff)[k])
          {
            //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
            //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<valence_[k]<<std::endl;
            const int fvi = vi*numdofpernode_+k;
            // We use the fact, that the rhs vector value for boundary nodes
            // is equivalent to the integrated negative normal flux
            // due to diffusion and migration
            val = erhs[fvi];
            erhs[vi*numdofpernode_+numscal_] += valence_[k]*(-val);
            // corresponding linearization
            for (int ui=0; ui<nen_; ++ui)
            {
              val = emat(vi*numdofpernode_+k,ui*numdofpernode_+k);
              emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k)+=valence_[k]*(-val);
              val = emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_);
              emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_)+=valence_[k]*(-val);
            }
          }
        } // for k
      } // if Dirichlet at node vi
    } // for vi
  }  // elim
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  get the material constants  (private)                     ehrl 01/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>  reamanager,   //!< reaction manager
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  if(material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
    }
  }
  else
    dserror("");

  return;
} //ScaTraEleCalc::GetMaterialParams


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  if(material->MaterialType() == INPAR::MAT::m_ion)
    myelch::MatIon(material,k,densn,densnp,densam,Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(diffmanager),reamanager,visc,iquad);
  else dserror("Material type is not supported");

  return;
}



// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs27>;



