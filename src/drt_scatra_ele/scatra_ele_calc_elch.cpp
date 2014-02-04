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

#include "scatra_ele_calc_elch.H"
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
DRT::ELEMENTS::ScaTraEleCalcElch<distype> * DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcElch<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcElch<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElch<distype>::ScaTraEleCalcElch(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    diffcondmat_(INPAR::ELCH::diffcondmat_undefined),
    epotnp_(my::numscal_),
    ecurnp_(my::numscal_),
    //TODO: ELCH: cleaning
    migrationstab_(false),
    migrationintau_(false),
    migrationinresidual_(true),
    transelim_(0.0),
    diffuselimderiv_(my::numscal_,0.0),
    diffuselim_(0.0),
    eps_(1,1.0),
    tort_(1,1.0),
    epstort_(1,1.0)
{
  // initialization of diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElch(my::numscal_));

  // set appropriate parameter list
  my::scatrapara_ = DRT::ELEMENTS::ScaTraEleParameterElch::Instance();
  elchpara_ = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_);

  // flag: current solution variable
  cursolvar_ = elchpara_->CurSolVar();

  // flag if elch system is closed by ENC or div i
  equpot_= elchpara_->EquPot();

  return;
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Evaluate(
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
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  //get element coordinates
  GEO::fillInitialPositionArray<distype,my::nsd_,LINALG::Matrix<my::nsd_,my::nen_> >(ele,my::xyze_);

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,my::myknots_,my::weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  // get myphi vector for potential and current
  const std::vector<double> myphinp = my::ExtractElementAndNodeValues(ele,params,discretization,lm);

  // get additional values for el. potential at element nodes
  for (int ien=0;ien<my::nen_;++ien)
    epotnp_(ien) = myphinp[ien*my::numdofpernode_+my::numscal_];

  //TODO:SCATRA_ELE_CLEANING: Generalized-alpha time integration scheme: Auswertung von strom richtig?
  //TODO: SCATRA_ELE_CLEANING: Wie soll man das steuern? material, Eingabeparameter, beides?
  if(cursolvar_)
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

  //--------------------------------------------------------------------------------
  // prepare turbulence models
  //--------------------------------------------------------------------------------

  int nlayer = 0;
  my::ExtractTurbulenceApproach(ele,params,discretization,lm,nlayer);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  Sysmat(
    ele,
    elemat1_epetra,
    elevec1_epetra,
    elevec2_epetra);

//  if(my::eid_==10)
//  {
//  std::cout << "matrix: " << std::endl;
//  std::cout <<  elemat1_epetra << std::endl;
//  std::cout << "vector: " << std::endl;
//  std::cout <<  elevec1_epetra << std::endl;
//  }

  // for certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries
  CorrectionForFluxAccrosDC(discretization,lm,elemat1_epetra,elevec1_epetra);

#if 0
  // for debugging of matrix entries
  if((ele->Id()==2) and (time < 1.3 and time > 1.1))
  {
    FDcheck(
      ele,
      elemat1_epetra,
      elevec1_epetra,
      elevec2_epetra,
      time,
      dt,
      timefac,
      alphaF,
      whichassgd,
      whichfssgd,
      assgd,
      fssgd,
      turbmodel_,
      Cs,
      tpn,
      frt,
      scatratype);
  }
#endif

  return 0;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 ehrl  08/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Sysmat(
  DRT::Element*                         ele, ///< the element those matrix is calculated
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs, ///< element rhs to calculate
  Epetra_SerialDenseVector&             subgrdiff ///< subgrid-diff.-scaling vector
  )
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

  // definition of constants
  const double faraday = INPAR::ELCH::faraday_const;
  const double epsilon = INPAR::ELCH::epsilon;
  const double frt=elchpara_->FRT();

  // pre-calculation of divisions which are regularly used in following
  // (a division is much more expensive than a multiplication)
  // constant parameter 1/F
  const double invf = 1/faraday;
  // constant parameter RT/F
  const double rtf = 1/frt;
  // constant parameter RT/F^2/constC
  // attentio: this is a newman specific parameter
  const double rtffc=rtf*invf/elchpara_->NewmanConstC();

  // flag if term for diffusion overpotential is based on the diffusion coefficents or the transference numbers
  bool diffbased = elchpara_->DiffusionCoeffBased();

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
    GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc);

  //! the stabilisation parameters (one per transported scalar)
  std::vector<double> tau(my::numscal_);
  std::vector<LINALG::Matrix<my::nen_,1> > tauderpot(my::numscal_);

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
      GetMaterialParams(ele,densn,densnp,densam,my::diffmanager_,my::reamanager_,visc,iquad);

    // pre-calculation of divisions which are regularly used in following
    // (a division is much more expensive than a multiplication)
    std::vector <double> invfval(my::numscal_);
    std::vector <double> rtffcval(my::numscal_);
    for(int k=0; k<my::numscal_; ++k)
    {
      invfval[k]=invf/dme->GetValence(k);
      rtffcval[k]=rtffc/dme->GetValence(k);
    }

    // get concentration of transported scalar k at integration point
    // evaluation of all concentrations is necessary at this point since
    // -> homogeneous reactions of scalar k may depend on all concentrations
    // -> concentration depending material parameters for the diffusion-confection formulation
    // -> avoiding of possible errors (concentration was always defined as a vector where only on
    //    entry was filled)
    std::vector<double> conint(my::numscal_);
    std::vector<double> conintinv(my::numscal_);
    for (int k = 0;k<my::numscal_;++k)
    {
      conint[k] = my::funct_.Dot(my::ephinp_[k]);
      // pre-calculation of divisions which are regularly used in following
      // (a division is much more expensive than a multiplication)
      conintinv[k] = 1/conint[k];
    }

    // loop over all transported scalars
    std::vector<LINALG::Matrix<my::nsd_,1> > gradphi(my::numscal_);
    for (int k = 0; k < my::numscal_;++k)
      gradphi[k].Multiply(my::derxy_,my::ephinp_[k]);


    // get gradient of electric potential at integration point
    LINALG::Matrix<my::nsd_,1> gradpot(true);
    gradpot.Multiply(my::derxy_,epotnp_);

    // current density at Gauss point
    LINALG::Matrix<my::nsd_,1> curint;
    curint.Multiply(ecurnp_,my::funct_);

    double curdiv(true);
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      //get vdiv at time n+1 for np_genalpha,
      LINALG::Matrix<my::nsd_,my::nsd_> curderxy(true);
      curderxy.MultiplyNT(ecurnp_,my::derxy_);
      curdiv += curderxy(idim, idim);
    }

    // get velocity at integration point
    LINALG::Matrix<my::nsd_,1> velint(true);
    LINALG::Matrix<my::nsd_,1> convelint(true);
    velint.Multiply(my::evelnp_,my::funct_);
    convelint.Multiply(my::econvelnp_,my::funct_);

    // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
    LINALG::Matrix<my::nen_,1> conv(true);
    conv.MultiplyTN(my::derxy_,convelint);

    // dummy variable
    LINALG::Matrix<my::nen_,1> sgconv(true);

    // velocity divergence required for conservative form
    double vdiv(0.0);
    if (my::scatrapara_->IsConservative())
      this->GetDivergence(vdiv,my::evelnp_);

    LINALG::Matrix<my::nen_,1> migconv;
    if(elchpara_->NernstPlanck())
    {
      // migration term (convective part without z_k D_k): -F/RT\grad{\Phi}\grad
      migconv.MultiplyTN(-frt,my::derxy_,gradpot);
    }

    // TODO: check ELCH: time integration
    // stabilization parameter and integration factors
    const double timefacfac = my::scatraparatimint_->TimeFac()*fac;
    //TODO: ELCH stabilization
    const double timetaufac = 0.0;
    double rhsfac    = my::scatraparatimint_->TimeFacRhs()*fac;

    // loop all scalars
    for (int k=0;k<my::numscal_;++k) // deal with a system of transported scalars
    {
      const double taufac = tau[k]*fac;  // corresponding stabilization parameter
      //double rhsint    = rhs;
      double rhstaufac = my::scatraparatimint_->TimeFacRhsTau() * taufac;
      double residual     = 0.0;

      // scalar at integration point at time step n
      //const double phin = my::funct_.Dot(my::ephin_[k]);

      // convective term using current scalar value
      //double conv_phi(0.0);
      //conv_phi = convelint.Dot(gradphi[k]);

      // further short cuts and definitions
      const double conv_ephinp_k = conv.Dot(my::ephinp_[k]);

      // diffusive part used in stabilization terms
      //double diff_phi(0.0);
      // reactive part of migration term used in stabilization terms
      LINALG::Matrix<my::nen_,1> migrea;
      // diffusive term
      LINALG::Matrix<my::nen_,1> diff(true);

      // diffusive term using current scalar value for higher-order elements
      if (my::use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        my::GetLaplacianStrongForm(diff);
        diff.Scale(dme->GetIsotropicDiff(k));
        //diff_phi = diff.Dot(my::ephinp_[k]);

        LINALG::Matrix<my::nen_,1> laplace;
        my::GetLaplacianStrongForm(laplace);

        // get Laplacian of electric potential at integration point
        double lappot = laplace.Dot(epotnp_);
        // reactive part of migration term
        migrea.Update(-frt*dme->GetIsotropicDiff(k)*dme->GetValence(k)*lappot,my::funct_);

        //diff_ephinp_k = diff_.Dot(ephinp_[k]);   // diffusion
        //migrea_k      = migrea_.Dot(ephinp_[k]); // reactive part of migration term
      }

      // reactive part of the form: (reaction coefficient)*phi
      //double rea_phi(0.0);
      //rea_phi = densnp*phinp*my::reamanager_->GetReaCoeff(k);

      // get history data (or acceleration)
      double hist(0.0);
      hist = my::funct_.Dot(my::ehist_[k]);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhs(0.0);
      my::GetRhs(rhs,densnp,k);

      double rhsint = rhs;
      // TODO: ELCH: Check implementation of time integration schemes carefully
      // perform time-integration specific actions
      if (not my::scatraparatimint_->IsStationary())
      {
        //TODO: ELCH: Zeitintegration genalpha Konzentrationsabhängige Parameter alter Zeitschritt
        if (my::scatraparatimint_->IsGenAlpha())
        {
          // note: in hist_ we receive the time derivative phidtam at time t_{n+alpha_M} !!
          //residual  = hist_[k] + conv_ephinp_k - diff_ephinp_k - rhsint;

          rhsint   *= my::scatraparatimint_->TimeFacRhs(); //(timefac/alphaF);  // not nice, but necessary !

          // rhs contribution due to incremental formulation (phidtam)
          // Standard Galerkin term
          const double vtrans = rhsfac*eps_[0]*hist;
          for (int vi=0; vi<my::nen_; ++vi)
          {
            const int fvi = vi*my::numdofpernode_+k;

            erhs[fvi] -= vtrans*my::funct_(vi);
          }
        }
        else
        {
          // TODO: do I need this term
          rhsint = eps_[0]*hist + (rhs*my::scatraparatimint_->TimeFac()); // contributions from t_n and \theta*dt*bodyforce(t_{n+1})
          //residual  = conint_[k] + timefac*(conv_ephinp_k - diff_ephinp_k) - rhsint;

          // rhs contribution due to incremental formulation (phinp)
          // Standard Galerkin term
          const double vtrans = fac*eps_[0]*conint[k];
          for (int vi=0; vi<my::nen_; ++vi)
          {
            const int fvi = vi*my::numdofpernode_+k;

            erhs[fvi] -= vtrans*my::funct_(vi);
          }
        } // if(is_genalpha_)
      } // if (is_stationary_)

      // if(Nernst-Planck formulation or diffusion-conduction formulation)
      if(elchpara_->NernstPlanck())
      {
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
      }
      // else(Nernst-Planck formulation or diffusion-conduction formulation)
      else
      {
        // specific constants for the Newman-material:
        // switch between a dilute solution theory like formulation and the classical concentrated solution theory
        double newman_const_a = elchpara_->NewmanConstA();
        double newman_const_b = elchpara_->NewmanConstB();
        //Newmans-specific costant newman_const_c is included in pre-calculation of divisions

        //TODO: SCATRA_ELE_CLEANING: generalized alpha scheme and history vector??

        //----------------------------------------------------------------
        // 1) element matrix: instationary terms
        //----------------------------------------------------------------

        if (not my::scatraparatimint_->IsStationary())
          my::CalcMatMass(emat,k,fac,densam*eps_[0],densnp);

        //----------------------------------------------------------------
        // 2) element matrix: stationary terms
        //----------------------------------------------------------------

        // 2a)  element matrix: convective term
        my::CalcMatConv(emat,k,timefacfac,densnp*eps_[0],conv,sgconv);

        // 2b)  element matrix: diffusion term
        //      - constant diffusion coefficient
        my::CalcMatDiff(emat,k,timefacfac*epstort_[0],dme);

        //      - concentation depending diffusion coefficient
        //        (additional term for Newman material)
        if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
        {
          //linearization of diffusion coefficient in the ionic diffusion term (transport equation)
          //
          // (grad w, D(D(c)) grad c)
          //
          for (int vi=0; vi<my::nen_; ++vi)
          {
            for (int ui=0; ui<my::nen_; ++ui)
            {
              double laplawfrhs_gradphi=0.0;
              my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[k],vi);

              emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k)
                += timefacfac*epstort_[0]*dme->GetDerivIsoDiffCoef(k,k)*laplawfrhs_gradphi*my::funct_(ui);
            }
          }
        }

        //----------------------------------------------------------------------------------------------
        // electrical conduction term (transport equation)
        //----------------------------------------------------------------------------------------------

        // equation for current is inserted in the mass transport equation
        //
        // mass transport equation:
        //            |     diffusion term      | |     conduction term    |
        //            |                         | |                        |
        //   dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
        //
        // equation for current:
        //   i = - kappa nabla phi + RT/F kappa (thermfactor) f(t_k) nabla ln c_k
        //
        if (not cursolvar_)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            for (int ui=0; ui<my::nen_; ++ui)
            {
              double laplawf(0.0);
              my::GetLaplacianWeakForm(laplawf,ui,vi); // compute once, reuse below!

              // linearization of conduction term depending on the potential
              //
              // (grad w, t_k kappa/(F z_k) D(grad phi))
              //
              emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=
                timefacfac*epstort_[0]*dme->GetTransNum(k)*dme->GetCond()*invfval[k]*laplawf;

              double laplawfrhs_gradpot(0.0);
              my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);

              for(int iscal=0; iscal<my::numscal_;++iscal)
              {
                //linearization of the conductivity in the conduction term depending on the potential
                //
                // (grad w, t_k D(kappa(c))/(F z_k) grad phi)
                //
                emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                   += timefacfac*epstort_[0]*dme->GetTransNum(k)*invfval[k]*dme->GetDerivCond(iscal)*my::funct_(ui)*laplawfrhs_gradpot;

                //linearization of the transference number in the conduction term depending on the potential
                //
                // (grad w, D(t_k(c)) kappa/(F z_k) grad phi)
                //
                emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                   += timefacfac*epstort_[0]*(dme->GetDerivTransNum(k,iscal))*my::funct_(ui)*invfval[k]*dme->GetCond()*laplawfrhs_gradpot;
              }

              // inconsitstence for ..
              if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
              {
                double laplawf(0.0);
                my::GetLaplacianWeakForm(laplawf,ui,vi);

                for(int iscal=0;iscal<my::numscal_;++iscal)
                {
                  // linearization of conduction term depending on the concentration
                  //
                  // (grad w, RT/(z_k F^2) kappa thermfac f(t_+) D(grad ln c_k))
                  //
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffcval[k]*dme->GetTransNum(k)*dme->GetCond()*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawf;

                  // linearization of conduction term depending on the concentration is implemented
                  // only for one species
                  // otherwise you would need a second loop over the all scalars
                  double laplawfrhs_gradc(0.0);
                  my::GetLaplacianWeakFormRHS(laplawfrhs_gradc,gradphi[iscal],vi);

                  // Linearization wrt ln c
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += -timefacfac*epstort_[0]*rtffcval[k]*dme->GetTransNum(k)*dme->GetCond()*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*conintinv[iscal]*laplawfrhs_gradc*my::funct_(ui);

                  // Linearization wrt kappa
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffcval[k]*dme->GetTransNum(k)*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradc*dme->GetDerivCond(iscal)*my::funct_(ui);

                  // Linearization wrt transference number 1
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffcval[k]*dme->GetCond()*conintinv[iscal]*laplawfrhs_gradc*dme->GetThermFac()*(newman_const_a+newman_const_b*dme->GetTransNum(iscal))*dme->GetDerivTransNum(iscal,iscal)*my::funct_(ui);

                  // Linearization wrt transference number 2
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffcval[k]*dme->GetCond()*dme->GetThermFac()*conintinv[iscal]*laplawfrhs_gradc*dme->GetTransNum(iscal)*newman_const_b*dme->GetDerivTransNum(iscal,iscal)*my::funct_(ui);

                  // Linearization wrt thermodynamic factor
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffcval[k]*dme->GetTransNum(k)*dme->GetCond()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradc*dme->GetDerivThermFac(iscal)*my::funct_(ui);
                }
              }
            } // for ui
          } // for vi
        } // end if(cursolvar_)
        // equation for current is solved independently
        //
        // mass transport equation:
        //            |     diffusion term      | |     conduction term    |
        //            |                         | |                        |
        //   dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
        //
        // equation for current:
        //   i = - kappa nabla phi + RT/F kappa (thermfactor) f(t_k) nabla ln c_k
        //
        else
        {
          // current term (with current as a solution variable)
          for (int vi=0; vi<my::nen_; ++vi)
          {
            const int    fvi = vi*my::numdofpernode_+k;

            for (int ui=0; ui<my::nen_; ++ui)
            {
              for (int idim=0; idim<my::nsd_; ++idim)
              {
                const int fui = ui*my::numdofpernode_+(my::numscal_+1)+idim;

                //linearization of conduction term depending on current flow
                //
                // (grad w, t_k/(F z_k) Di)
                //
                emat(fvi,fui) += -timefacfac*my::derxy_(idim,vi)*dme->GetTransNum(k)*invfval[k]*my::funct_(ui);

                //linearization of transference number in conduction term depending on current flow
                //
                // (grad w, Dt_k(c)/(F z_k) i)
                //
                for (int iscal=0; iscal<my::numscal_; ++iscal)
                  emat(fvi,ui*my::numdofpernode_+iscal)
                    += -timefacfac*my::derxy_(idim,vi)*(dme->GetDerivTransNum(k,iscal))*my::funct_(ui)*invfval[k]*curint(idim);
              }

              // this coupling term cancels out for a 2 equation system
              if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
              {
                // compute once, reuse below!
                double laplawf(0.0);
                my::GetLaplacianWeakForm(laplawf,ui,vi);

                for (int iscal=0; iscal < my::numscal_; ++iscal)
                {
                  // formulation a): plain ionic diffusion coefficients without using ENC
                  //
                  // (grad w, sum(D_i grad Dc))
                  //
                  emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                    += -timefacfac*epstort_[0]*dme->GetTransNum(k)/dme->GetValence(k)*dme->GetValence(iscal)*dme->GetIsotropicDiff(iscal)*laplawf;

                  //linearization of transference number in the coupling term (transport equation)
                  //
                  // (grad w, Dt_k(c)/z_k (Sum_i z_i D_i grad c_i))
                  //
                  for(int iscal2=0; iscal2<my::numscal_; ++iscal2)
                  {
                    double laplawfrhs_gradphi=0.0;
                    this->GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[iscal2],vi);

                    emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                      += -timefacfac*epstort_[0]*(dme->GetDerivTransNum(k,iscal))*my::funct_(ui)*dme->GetValence(iscal2)/dme->GetValence(k)*dme->GetIsotropicDiff(iscal2)*laplawfrhs_gradphi;
                  } // for(iscal2)

                  // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
                  //                  -> not implemented yet

                  // formulation c):  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i) replaced by sum (t_i/c_i nabla c_i)
                  //                  -> not activated
                  //                  -> linearization is missing
                  // emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
                  //    += timefacfac*epstort_[0]/frt/faraday/dme->GetValence(k)*dme->GetTransNum(k)*dme->GetCond()*dme->GetTransNum(iscal)/dme->GetValence(iscal]/conint[iscal]*laplawf;
                } // end for(iscal)
              } // end if(diffcondmat)
            } // end for ui
          } // end for vi
        } // end if(not cursolvar_)

        //-----------------------------------------------------------------------
        // 3) element right hand side vector (neg. residual of nonlinear problem)
        //-----------------------------------------------------------------------

        //----------------------------------------------------------------
        // 3a) rhs: governing equation for concentration
        //----------------------------------------------------------------
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const int fvi = vi*my::numdofpernode_+k;

          // RHS source term (contains old part of rhs for OST / BDF2)
          erhs[fvi] += fac*eps_[0]*my::funct_(vi)*rhsint;

          // convective term
          erhs[fvi] -= rhsfac*eps_[0]*my::funct_(vi)*conv_ephinp_k;

          // diffusive term
          double laplawfrhs_gradphi(0.0);
          this->GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[k],vi);
          erhs[fvi] -= rhsfac*epstort_[0]*dme->GetIsotropicDiff(k)*laplawfrhs_gradphi;

          if(not cursolvar_)
          {
            // diffusive term
            double laplawfrhs_gradpot=0.0;
            my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);
            erhs[fvi]-= rhsfac*epstort_[0]*dme->GetTransNum(k)*dme->GetCond()*invfval[k]*laplawfrhs_gradpot;

            if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
            {
              for(int iscal=0; iscal<my::numscal_;++iscal)
              {
                // diffusive term second
                double laplawfrhs_gradphi(0.0);
                this->GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[iscal],vi); // compute once, reuse below!

                // formulation a): plain ionic diffusion coefficients without using ENC
                //
                // (grad w, sum(D_i grad Dc))
                erhs[fvi]
                  -= rhsfac*epstort_[0]*rtffcval[k]*dme->GetTransNum(k)*dme->GetCond()*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradphi;

                // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
                //                  -> not implemented yet
              }
            }
            // if(diffcondmat_==INPAR::ELCH::diffcondmat_ion): all terms cancel out
          }
          else
          {
            double laplawfrhs_cur=0.0;
            my::GetLaplacianWeakFormRHS(laplawfrhs_cur,curint,vi);
            erhs[fvi]-= -rhsfac*dme->GetTransNum(k)*invfval[k]*laplawfrhs_cur;

            // this coupling term cancels out for a 2 equation system
            // (current not a solution variable)
            if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
            {
              for (int iscal=0; iscal <my::numscal_; ++iscal)
              {
                // diffusive term second
                double laplawfrhs_gradphi(0.0);
                this->GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[iscal],vi); // compute once, reuse below!

                // formulation a:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i)
                erhs[fvi] -= - rhsfac*epstort_[0]*dme->GetTransNum(k)*dme->GetValence(iscal)/dme->GetValence(k)*dme->GetIsotropicDiff(iscal)*laplawfrhs_gradphi;

                // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
                //                  -> not implemented yet
                // formulation c:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i) replaced by sum (t_i/c_i nabla c_i)
                //                  -> not activated
                // erhs[fvi]
                //   -= rhsfac*epstort_[0]/frt/faraday/dme->GetValence(k)*dme->GetTransNum(k)*dme->GetCond()*dme->GetTransNum(iscal)/dme->GetValence(iscal]/conint[iscal]*laplawf2;
              }
            }
          }
        } // for vi
      } // end(Nernst-Planck formulation or diffusion-conduction formulation)
    }  // end loop over scalar

    // if(Nernst-Planck formulation or diffusion-conduction formulation)
    if(not elchpara_->NernstPlanck())
    {
      // specific constants for the Newman-material:
      // switch between a dilute solution theory like formulation and the classical concentrated solution theory
      double newman_const_a = elchpara_->NewmanConstA();
      double newman_const_b = elchpara_->NewmanConstB();
      //Newmans-specific costant newman_const_c is included in pre-calculation of divisions

      //-------------------------------------------------------------------------
      // governing equation for current
      //-------------------------------------------------------------------------

      // equation for current is inserted in the mass transport equation
      //
      // mass transport equation:
      //            |     diffusion term      | |     conduction term    |
      //            |                         | |                        |
      //   dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k /z_k (i/F))
      //
      // equation for current:
      //   i/F = - 1/F kappa nabla phi + RT/F^2 kappa (thermfactor) f(t_k) nabla ln c_k
      //
      if(not cursolvar_)
      {
        if(equpot_==INPAR::ELCH::equpot_divi)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            for (int ui=0; ui<my::nen_; ++ui)
            {
              double laplawf(0.0);
              my::GetLaplacianWeakForm(laplawf,ui,vi); // compute once, reuse below!

              // linearization of the ohmic term
              //
              // (grad w, 1/F kappa D(grad pot))
              //
              emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += timefacfac*epstort_[0]*invf*dme->GetCond()*laplawf;

              for(int iscal=0;iscal<my::numscal_;++iscal)
              {
                double laplawfrhs_gradpot(0.0);
                my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);

                // linearization of the ohmic term wrt conductivity
                //
                // (grad w, 1/F kappa D(grad pot))
                //
                emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                  += timefacfac*epstort_[0]*invf*dme->GetDerivCond(iscal)*my::funct_(ui)*laplawfrhs_gradpot;

                // linearization of the diffusion overpotential term
                //
                // (grad w, RT/F^2 kappa (thermfactor) f(t_k) D nabla ln c_k)
                //
                // thermodynamic factor only implemented for Newman
                if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
                {
                  // linearization of the diffusion overpotential term
                   //
                   // (grad w, RT/F^2 kappa (thermfactor) f(t_k) 1/c_k D nabla c_k)
                   //
                  emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffc*dme->GetCond()*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawf;

                  // linearization of conduction term depending on the concentration is implemented
                  // only for one species
                  // otherwise you would need a second loop over the all scalars
                  double laplawfrhs_gradphi(0.0);
                  my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[iscal],vi);

                  // Linearization wrt ln c
                  emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                    += -timefacfac*epstort_[0]*rtffc*dme->GetCond()*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*conintinv[iscal]*laplawfrhs_gradphi*my::funct_(ui);

                  // Linearization wrt kappa
                  emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffc*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradphi*dme->GetDerivCond(iscal)*my::funct_(ui);

                  // Linearization wrt transference number
                  emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffc*dme->GetCond()*dme->GetThermFac()*conintinv[iscal]*laplawfrhs_gradphi*newman_const_b*dme->GetDerivTransNum(iscal,iscal)*my::funct_(ui);

                  // Linearization wrt thermodynamic factor
                  emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                    += timefacfac*epstort_[0]*rtffc*dme->GetCond()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradphi*dme->GetDerivThermFac(iscal)*my::funct_(ui);
                }
                else if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
                {
                  if(diffbased==true)
                    emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                      += timefacfac*epstort_[0]*dme->GetValence(iscal)*dme->GetIsotropicDiff(iscal)*laplawf;
                  else
                    emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
                      += timefacfac*epstort_[0]*rtf*invf/dme->GetValence(iscal)*dme->GetCond()*dme->GetTransNum(iscal)*conintinv[iscal]*laplawf;
                }
                else
                  dserror("Diffusion-Conduction material is not specified");
              }
            }// for ui
          }// for vi
        }
        else if(equpot_==INPAR::ELCH::equpot_enc)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            for (int ui=0; ui<my::nen_; ++ui)
            {
              for (int iscal=0; iscal < my::numscal_; ++iscal)
              {
                //linearization of the transference number in the conduction term (transport equation)
                //
                // (w, sum(z_k c_k))
                //
                emat(vi*my::numdofpernode_+my::numscal_, ui*my::numdofpernode_+iscal) += my::scatraparatimint_->AlphaF()*dme->GetValence(iscal)*fac*my::funct_(vi)*my::funct_(ui);
              }
            }
          }
        }
        else
          dserror("");
      } //end if(not cursolvar_)
      // equation for current is solved independently
      //
      // mass transport equation:
      //            |     diffusion term      | |     conduction term    |
      //            |                         | |                        |
      //   dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
      //
      // equation for current:
      //   i = - kappa nabla phi + RT/F kappa (thermfactor) f(t_k) nabla ln c_k
      //
      else
      {
        // (v, i)
        for (int vi=0; vi<my::nen_; ++vi)
        {
          for (int ui=0; ui<my::nen_; ++ui)
          {
            for (int idim=0; idim<my::nsd_; ++idim)
            {
              const int fvi = vi*my::numdofpernode_+(my::numscal_+1)+idim;
              const int fui = ui*my::numdofpernode_+(my::numscal_+1)+idim;

              emat(fvi,fui) += timefacfac*invf*my::funct_(vi)*my::funct_(ui);
            }
          }
        }

        // (v, kappa grad phi)
        for (int vi=0; vi<my::nen_; ++vi)
        {
          for (int ui=0; ui<my::nen_; ++ui)
          {
            for (int idim=0; idim<my::nsd_; ++idim)
            {
              const int fvi = vi*my::numdofpernode_+(my::numscal_+1)+idim;
              const int fui = ui*my::numdofpernode_+my::numscal_;

              emat(fvi,fui) += timefacfac*invf*epstort_[0]*my::funct_(vi)*dme->GetCond()*my::derxy_(idim,ui);

              //linearization of conductivity in the ohmic resistance term (current equation)
              //
              // (w, D(kappa(c)) grad phi)
              //
              for(int k=0;k<my::numscal_;++k)
                emat(fvi,ui*my::numdofpernode_+k) += timefacfac*invf*epstort_[0]*my::funct_(vi)*dme->GetDerivCond(k)*my::funct_(ui)*gradpot(idim);
            }
          }
        }

        for (int vi=0; vi<my::nen_; ++vi)
        {
          for (int ui=0; ui<my::nen_; ++ui)
          {
            // diffusive term
            // (grad w, D grad c)
            for (int idim = 0; idim < my::nsd_; ++idim)
            {
              for (int iscal=0; iscal < my::numscal_; ++iscal)
              {
                if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
                {
                  emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+iscal)
                    += timefacfac*rtffc*epstort_[0]*my::funct_(vi)*dme->GetCond()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*my::derxy_(idim,ui);

                  //  // linearization of coupling term in the current equation is still missing
                  //  emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                  //    += timefacfac*my::funct_(vi)/frt_*dme->GetCond()*dme->GetTransNum(iscal)/dme->GetValence(iscal]/((-1)*conint_[iscal]*conint_[iscal])*my::funct_(ui)*gradphicoupling_[iscal](idim);
                  //
                  //  for(int iscal2=0; iscal2<numscal_;++iscal2)
                  //  {
                  //    //Check if necessary??
                  //    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                  //       += timefacfac*my::funct_(vi)/frt_*dme->GetDerivCond(iscal)*my::funct_(ui)*trans_[iscal2]*my::funct_(ui)/dme->GetValence(iscal2]/conint_[iscal2]*gradphicoupling_[iscal2](idim);
                  //
                  //
                  //    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+iscal)
                  //      += timefacfac*my::funct_(vi)/frt_*dme->GetCond()*(transderiv_[iscal2])[iscal]*my::funct_(ui)/dme->GetValence(iscal2]/conint_[iscal2]*gradphicoupling_[iscal2](idim);
                  //  }
                }
                else if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
                {
                  if(diffbased==true)
                    emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+iscal)
                      += timefacfac*epstort_[0]*my::funct_(vi)*dme->GetValence(iscal)*dme->GetIsotropicDiff(iscal)*my::derxy_(idim,ui);
                  else
                    emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+iscal)
                      += timefacfac*epstort_[0]*invfval[iscal]*rtf*my::funct_(vi)*dme->GetCond()*dme->GetTransNum(iscal)*conintinv[iscal]*my::derxy_(idim,ui);
                }
                else
                  dserror("Diffusion-Conduction material is not specified");
              }
            }
          } // for ui
        } // for vi

        //-------------------------------------------------------------------------
        // 2c) element matrix: governing equation for potential (div i)
        //-------------------------------------------------------------------------
        if(equpot_==INPAR::ELCH::equpot_divi)
        {
          for (int vi=0; vi<my::nen_; ++vi)
           {
             for (int ui=0; ui<my::nen_; ++ui)
             {
               for (int idim = 0; idim <my::nsd_; ++idim)
               {
                 const int fvi = my::numdofpernode_*vi+my::numscal_;
                 const int fui = my::numdofpernode_*ui+(my::numscal_+1)+idim;
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
                 emat(fvi,fui) -= timefacfac*invf*my::derxy_(idim,vi)*my::funct_(ui);
                 // version b: (phi, div Di) -> not partially integrated
                 //emat(fvi,fui) += timefacfac*funct_(vi)*derxy_(idim,ui);
               } // end for(idim)
             } // end for(ui)
           } // end for(vi)
        }
        else if(equpot_==INPAR::ELCH::equpot_enc)
        {
          for (int vi=0; vi<my::nen_; ++vi)
           {
             for (int ui=0; ui<my::nen_; ++ui)
             {
               for (int iscal=0; iscal < my::numscal_; ++iscal)
               {
                 //linearization of the transference number in the conduction term (transport equation)
                 //
                 // (w, sum(z_k c_k))
                 //
                 emat(vi*my::numdofpernode_+my::numscal_, ui*my::numdofpernode_+iscal) += my::scatraparatimint_->AlphaF()*dme->GetValence(iscal)*fac*my::funct_(vi)*my::funct_(ui);
               }
             }
           }
         }
         else
           dserror("");
      }

      //----------------------------------------------------------------
      // 2d) Stabilization terms
      //----------------------------------------------------------------

      // no stabilization terms

      //----------------------------------------------------------------
      // 3b) rhs: governing equation for current
      //----------------------------------------------------------------
      if(not cursolvar_)
      {
        // add rhs
        if(equpot_==INPAR::ELCH::equpot_divi)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            double laplawfrhs_gradpot(0.0);
            this->GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi); // compute once, reuse below!

            erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*epstort_[0]*invf*dme->GetCond()*laplawfrhs_gradpot;

            for (int iscal=0; iscal < my::numscal_; ++iscal)
            {
              if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
              {
                // diffusive term second
                double laplawf2(0.0);
                this->GetLaplacianWeakFormRHS(laplawf2,gradphi[iscal],vi); // compute once, reuse below!

                if(diffbased==true)
                  erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*epstort_[0]*dme->GetValence(iscal)*dme->GetIsotropicDiff(iscal)*laplawf2;
                else
                  erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*epstort_[0]*rtf*invf/dme->GetValence(iscal)*dme->GetCond()*dme->GetTransNum(iscal)*conintinv[iscal]*laplawf2;
              }
              // thermodynamic factor only implemented for Newman
              else if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
              {
                // diffusive term second
                double laplawf2(0.0);
                this->GetLaplacianWeakFormRHS(laplawf2,gradphi[iscal],vi); // compute once, reuse below!

                erhs[vi*my::numdofpernode_+my::numscal_]
                  -= rhsfac*epstort_[0]*rtffc*dme->GetCond()*dme->GetThermFac()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*laplawf2;
              }
              else
                dserror("Diffusion-Conduction material is not specified");
            }
          }
        }
        else if(equpot_==INPAR::ELCH::equpot_enc)
        {
          for (int vi=0; vi<my::nen_; ++vi)
          {
            // electroneutrality condition
            // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
            for (int iscal=0; iscal < my::numscal_; ++iscal)
              erhs[vi*my::numdofpernode_+my::numscal_] -= dme->GetValence(iscal)*fac*my::funct_(vi)*conint[iscal];
          }
        }
        else
          dserror("");
      }
      else
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          for (int idim = 0; idim <my::nsd_; ++idim)
          {
            // (v, i)
            erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]-=rhsfac/faraday*my::funct_(vi)*curint(idim);

            // (v, kappa grad phi)
            erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]-=rhsfac/faraday*epstort_[0]*my::funct_(vi)*dme->GetCond()*gradpot(idim);
          }

          if(equpot_==INPAR::ELCH::equpot_divi)
          {
            // (v, i)
            //erhs[vi*my::numdofpernode_+(my::numdofpernode_-1)]-=rhsfac*my::funct_(vi)*curdiv_;
            {
              double laplawf=0.0;
              // version a: (grad phi,  Di)
              my::GetLaplacianWeakFormRHS(laplawf,curint,vi);
              erhs[vi*my::numdofpernode_+my::numscal_]-= -rhsfac/faraday*laplawf;
              // version b: (phi, div Di) -> not partially integrated
              //erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*my::funct_(vi)*divi;
            }
          }
          else if(equpot_==INPAR::ELCH::equpot_enc)
          {
            // electroneutrality condition
            // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
            for (int iscal=0; iscal < my::numscal_; ++iscal)
              erhs[vi*my::numdofpernode_+my::numscal_] -= dme->GetValence(iscal)*fac*my::funct_(vi)*conint[iscal];
          }
          else
            dserror("");

          // diffusive term
          // (grad w, D grad c)
          for (int idim = 0; idim < my::nsd_; ++idim)
          {
            for (int iscal=0; iscal < my::numscal_; ++iscal)
            {
              if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
              {
                erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]
                  -= rhsfac*epstort_[0]*my::funct_(vi)*rtffc*dme->GetCond()*(newman_const_a+(newman_const_b*dme->GetTransNum(iscal)))*conintinv[iscal]*gradphi[iscal](idim);
                //erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]
                //  -= rhsfac*my::funct_(vi)*faraday*(2.0e-5-4.0e-5)*gradphicoupling_[iscal](idim);
              }
              else if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
              {
                if(diffbased==true)
                  erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim] -= rhsfac/faraday*epstort_[0]*my::funct_(vi)*faraday*dme->GetValence(iscal)*dme->GetIsotropicDiff(iscal)*gradphi[iscal](idim);
                else
                {
                  erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]
                       -= rhsfac/faraday*epstort_[0]*my::funct_(vi)/frt*dme->GetCond()*dme->GetTransNum(iscal)/dme->GetValence(iscal)*conintinv[iscal]*gradphi[iscal](idim);
                }
              }
              else
                dserror("Diffusion-Conduction material is not specified");
            }
          }
        }
      }
    } // end(Nernst-Planck formulation or diffusion-conduction formulation)
  }

  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 ehrl  08/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CorrectionForFluxAccrosDC(
  DRT::Discretization&        discretization,
  const std::vector<int>&     lm,
  Epetra_SerialDenseMatrix&   emat,
  Epetra_SerialDenseVector&   erhs)
{
  // get dirichlet toggle from the discretization
  Teuchos::RCP<const Epetra_Vector> dctoggle = discretization.GetState("dctoggle");
    std::vector<double> mydctoggle(lm.size());
    DRT::UTILS::ExtractMyValues(*dctoggle,mydctoggle,lm);

  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

  double val=0.0;
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int k=0; k<my::numscal_; ++k)
    {
      //
      if (mydctoggle[vi*my::numdofpernode_+k] == 1)
      {
        const int fvi = vi*my::numdofpernode_+k;
        // We use the fact, that the rhs vector value for boundary nodes
        // is equivalent to the integrated negative normal flux
        // due to diffusion and migration

        // scaling of div i results in a matrix with better condition number
        val = erhs[fvi];
        erhs[vi*my::numdofpernode_+my::numscal_] += dme->GetValence(k)*(-val);
        // corresponding linearization
        for (int ui=0; ui<my::nen_; ++ui)
        {
          val = emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)+=dme->GetValence(k)*(-val);
          val = emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_)+=dme->GetValence(k)*(-val);
        }
      }

      // Dirichlet conditions on the potential are only allowed for the newman material
      // since additional information about the reacting species is required. This is fulfilled naturally for
      // the Newman material since only one species is allowed in this case.
      // Newman material models binary electrolytes where the second species is condensed via the ENC!
      if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
      {
        if (mydctoggle[vi*my::numdofpernode_+my::numscal_]==1)
        {
          // reacting species 0:
          // Newman material: reacting species is always the first species since there is only one species
          // other materials: one have to find a way to define the reacting species
          int k =0;

          const int fvi = vi*my::numdofpernode_+my::numscal_;
          // We use the fact, that the rhs vector value for boundary nodes
          // is equivalent to the integrated negative normal flux
          // due to diffusion and migration

          // scaling of div i results in a matrix with better condition number
          val = erhs[fvi];
          erhs[vi*my::numdofpernode_+k] += 1.0/dme->GetValence(k)*(-val);
          // corresponding linearization
          for (int ui=0; ui<my::nen_; ++ui)
          {
            val = emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k);
            emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += 1.0/dme->GetValence(k)*(-val);
            val = emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_);
            emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=1.0/dme->GetValence(k)*(-val);
          }
        }
      }
    }  // for k
  }

//   // for concentrated solution theory (using div i as closing term for the potential)
//   // additional flux terms / currents across Dirichlet boundaries
//   if(elchpara_->ElchType()==INPAR::ELCH::elchtype_diffcond and
//       diffcondmat_==INPAR::ELCH::diffcondmat_newman and
//       equpot_==INPAR::ELCH::equpot_divi)
//   {
//     //const double faraday = INPAR::SCATRA::faraday_const;
//     double val(0.0);
//     const DRT::Node* const* nodes = ele->Nodes();
//     std::string condname = "Dirichlet";
//
//     for (int vi=0; vi<my::nen_; ++vi)
//     {
//       std::vector<DRT::Condition*> dirichcond0;
//       nodes[vi]->GetCondition(condname,dirichcond0);
//
//       // there is at least one Dirichlet condition on this node
//       if (dirichcond0.size() > 0)
//       {
//         //std::cout<<"Ele Id = "<<ele->Id()<<"  Found one Dirichlet node for vi="<<vi<<std::endl;
//         const std::vector<int>*    onoff = dirichcond0[0]->Get<std::vector<int> >   ("onoff");
//         for (int k=0; k<my::numscal_; ++k)
//         {
//           if ((*onoff)[k])
//           {
//             //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
//             //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<valence_[k]<<std::endl;
//             const int fvi = vi*my::numdofpernode_+k;
//             // We use the fact, that the rhs vector value for boundary nodes
//             // is equivalent to the integrated negative normal flux
//             // due to diffusion and migration
//
//             // scaling of div i results in a matrix with better condition number
//             val = elevec1_epetra[fvi];
//             elevec1_epetra[vi*my::numdofpernode_+my::numscal_] += dme_->GetValence(k)*(-val);
//             // corresponding linearization
//             for (int ui=0; ui<my::nen_; ++ui)
//             {
//               val = elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
//               elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)+=dme_->GetValence(k)*(-val);
//               val = elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
//               elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_)+=dme_->GetValence(k)*(-val);
//             }
//           }
//         } // for k
//         // dirichlet condition for potential
//         if ((*onoff)[my::numscal_])
//         {
//           //reacting species 0
//           int k =0;
//
//           //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
//           //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<dme_->GetValence(k)<<std::endl;
//           const int fvi = vi*my::numdofpernode_+my::numscal_;
//           // We use the fact, that the rhs vector value for boundary nodes
//           // is equivalent to the integrated negative normal flux
//           // due to diffusion and migration
//
//           // scaling of div i results in a matrix with better condition number
//           val = elevec1_epetra[fvi];
//           elevec1_epetra[vi*my::numdofpernode_+k] += 1.0/dme_->GetValence(k)*(-val);
//           // corresponding linearization
//           for (int ui=0; ui<my::nen_; ++ui)
//           {
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += 1.0/dme_->GetValence(k)*(-val);
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=1.0/dme_->GetValence(k)*(-val);
//           }
//         }
//       } // if Dirichlet at node vi
//     } // for vi
//
//     // Nernst boundary conditions have to be handled like Dirichlet conditions!!!
//     std::string condname2 = "ElectrodeKinetics";
//     for (int vi=0; vi<my::nen_; ++vi)
//     {
//       std::vector<DRT::Condition*> elctrodeKinetics;
//       nodes[vi]->GetCondition(condname2,elctrodeKinetics);
//
//       // there is at least one Dirichlet condition on this node
//       if (elctrodeKinetics.size() == 1)
//       {
//         const int  kinetics = elctrodeKinetics[0]->GetInt("kinetic model");
//
//         if (kinetics==INPAR::SCATRA::nernst)
//         {
//           //reacting species 0
//           int k = 0;
//
//           //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
//           //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<dme_->GetValence(k)<<std::endl;
//           const int fvi = vi*my::numdofpernode_+my::numscal_;
//           // We use the fact, that the rhs vector value for boundary nodes
//           // is equivalent to the integrated negative normal flux
//           // due to diffusion and migration
//
//           // scaling of div i results in a matrix with better condition number
//           val = elevec1_epetra[fvi];
//           elevec1_epetra[vi*my::numdofpernode_+k] += 1.0/dme_->GetValence(k)*(-val);
//           // corresponding linearization
//           for (int ui=0; ui<my::nen_; ++ui)
//           {
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += 1.0/dme_->GetValence(k)*(-val);
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=1.0/dme_->GetValence(k)*(-val);
//           }
//         }
//       } // if Dirichlet at node vi
//     } // for vi
//   }

  return;
}


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                     ehrl 01/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::GetMaterialParams(
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

  if (material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat>& actmat
          = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

    // TODO: ELCH add to check validy
    // access mat_elchmat: container material for porous structures in elch
    //if (actmat->NumPhase() != 1) dserror("In the moment a single phase is only allowed.");

    // 1) loop over single phases
    for (int iphase=0; iphase < actmat->NumPhase();++iphase)
    {
      // access phase material
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlephase = actmat->PhaseById(phaseid);
      // get parameters (porosity and tortuosity) from mat_phase
      Materials(singlephase,iphase,densn,densnp,densam,diffmanager,reamanager,visc,iquad);

      // dynmic cast: get access to mat_phase
      const Teuchos::RCP<const MAT::ElchPhase>& actphase
                = Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(singlephase);

      // TODO: ELCH add to check validy
      // enough materials defined
      //if (actphase->NumMat() < my::numscal_) dserror("Not enough materials in MatList.");

      // 2) loop over materials of the single phase
      for (int imat=0; imat < actphase->NumMat();++imat)
      {
        const int matid = actphase->MatID(imat);
        Teuchos::RCP<const MAT::Material> singlemat = actphase->MatById(matid);

        Materials(singlemat,imat,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
      }
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
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Materials(
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
  if(material->MaterialType() == INPAR::MAT::m_elchphase)
    MatPhase(material,k,densn,densnp,densam,Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(diffmanager),reamanager,visc,iquad);
  else if(material->MaterialType() == INPAR::MAT::m_newman)
    MatNewman(material,k,densn,densnp,densam,Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(diffmanager),reamanager,visc,iquad);
  else if(material->MaterialType() == INPAR::MAT::m_ion)
    MatIon(material,k,densn,densnp,densam,Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(diffmanager),reamanager,visc,iquad);
  else dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material Phase                                           ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::MatPhase(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               iphase,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  const MAT::ElchPhase* actsinglemat = static_cast<const MAT::ElchPhase*>(material.get());

  // get porosity
  eps_[0] = actsinglemat->Epsilon();

  // get tortuosity
  tort_[0] = actsinglemat->Tortuosity();

  epstort_[0]=eps_[0]*tort_[0];

  return;
}


/*----------------------------------------------------------------------*
 |  Material Newman                                          ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::MatNewman(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  //materialNewman = true;
  const MAT::Newman* actmat = static_cast<const MAT::Newman*>(material.get());

  // TODO: ELCH add to check validy
  // Material Newman is derived for a binary electrolyte utilizing the ENC to condense the non-reacting species
  // -> k=0
  //if(my::numscal_>1)
  //  dserror("Material Newman is only valid for one scalar (binary electrolyte utilizing the ENC)");

  // set material
  diffcondmat_ = INPAR::ELCH::diffcondmat_newman;

  // concentration of species k at the integration point
  const double conint = my::funct_.Dot(my::ephinp_[k]);

  // valence of ionic species
  diffmanager->SetValence(actmat->Valence(),k);

  // concentration depending diffusion coefficient
  diffmanager->SetIsotropicDiff(actmat->ComputeDiffusionCoefficient(conint),k);
  // derivation of concentration depending diffusion coefficient wrt all ionic species
  diffmanager->SetDerivIsoDiffCoef(actmat->ComputeFirstDerivDiffCoeff(conint),k,k);

  // concentration depending transference number
  diffmanager->SetTransNum(actmat->ComputeTransferenceNumber(conint),k);
  // derivation of concentration depending transference number wrt all ionic species
  diffmanager->SetDerivTransNum(actmat->ComputeFirstDerivTrans(conint),k,k);

  // thermodynamic factor of electrolyte solution
  diffmanager->SetThermFac(actmat->ComputeThermFac(conint));
  // derivative of conductivity with respect to concentrations
  diffmanager->SetDerivThermFac(actmat->ComputeFirstDerivThermFac(conint),0);

  // conductivity and first derivative can maximally depend on one concentration
  // since time curve is used as input routine
  // conductivity of electrolyte solution
  diffmanager->SetCond(actmat->ComputeConductivity(conint));
  // derivative of conductivity with respect to concentrations
  diffmanager->SetDerivCond(actmat->ComputeFirstDerivCond(conint),0);

  return;
}


/*----------------------------------------------------------------------*
 |  Material ION                                             ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::MatIon(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  diffcondmat_ = INPAR::ELCH::diffcondmat_ion;

  const MAT::Ion* actmat = static_cast<const MAT::Ion*>(material.get());

  // valence of ionic species
  diffmanager->SetValence(actmat->Valence(),k);

  // concentration depending diffusion coefficient
  diffmanager->SetIsotropicDiff(actmat->Diffusivity(),k);

  // Loop over materials is finished - now all material parameter are set
  if(k==(my::numscal_-1))
  {
    // Material data of eliminated ion species is read from the LAST ion material
    // in the matlist!
    if(elchpara_->ElchType()==INPAR::ELCH::elchtype_enc_pde_elim)
    {
      diffmanager->IncreaseLengthVector(k, my::numscal_);

      // valence of ionic species
      diffmanager->SetValence(actmat->ElimValence(),my::numscal_);

      // concentration depending diffusion coefficient
      diffmanager->SetIsotropicDiff(actmat->ElimDiffusivity(),my::numscal_);
    }

    // calculation of conductivity and transference number based on diffusion coefficient and valence
    // for mat_ion
    if(elchpara_->ElchType()==INPAR::ELCH::elchtype_diffcond)
    {
      // dynamic cast to elch-specific diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

      std::vector<double> conint(my::numscal_);
      for (int k = 0;k<my::numscal_;++k)
        conint[k] = my::funct_.Dot(my::ephinp_[k]);

      dme->CalcConductivity(my::numscal_,INPAR::ELCH::faraday_const*elchpara_->FRT(),conint);
      dme->CalcTransNum(my::numscal_,conint);
    }
  }

  return;
}


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



