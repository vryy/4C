/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_diffcond.cpp

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_elch_diffcond.H"
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
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype> * DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcElchDiffCond<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcElchDiffCond<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::ScaTraEleCalcElchDiffCond(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalcElch<distype>::ScaTraEleCalcElch(numdofpernode,numscal),
    diffcondmat_(INPAR::ELCH::diffcondmat_undefined),
    cursolvar_(false),
    equpot_(INPAR::ELCH::equpot_undefined),
    eps_(1,1.0),
    tort_(1,1.0),
    epstort_(1,1.0)
{
  // initialization of diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCond(my::numscal_));

  // set appropriate parameter list
  myelch::elchpara_ = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_);

  // flag: current solution variable
  cursolvar_ = myelch::elchpara_->CurSolVar();

  // flag if elch system is closed by ENC or div i
  equpot_= myelch::elchpara_->EquPot();

  // initialize internal variable manager
  myelch::varmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerElchDiffCond<my::nsd_, my::nen_>(my::numscal_,my::nsd_,myelch::elchpara_));

  return;
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Evaluate(
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
    myelch::epotnp_(ien) = myphinp[ien*my::numdofpernode_+my::numscal_];

  //TODO:SCATRA_ELE_CLEANING: Generalized-alpha time integration scheme: Auswertung von strom richtig?
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

  myelch::Sysmat(
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
|  calculate system matrix and rhs                           ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalMatAndRhs(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&                 emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&                 erhs,         //!< element rhs to calculate+
    const int                                 k,            //!< index of current scalar
    const double                              fac,          //!< domain-integration factor
    const double                              timefacfac,   //!< domain-integration factor times time-integration factor
    const double                              rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,          //!< diffusion manager
    double&                                   rhsint,       //!< rhs at Gauss point
    const double                              hist          //!< history
  )
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc = Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCond>(dme);

  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCond <my::nsd_,my::nen_> > vmdc
    = Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchDiffCond <my::nsd_,my::nen_> >(vm);

  //----------------------------------------------------------------
  // 1) element matrix: instationary terms
  //----------------------------------------------------------------

  if (not my::scatraparatimint_->IsStationary())
    my::CalcMatMass(emat,k,fac,eps_[0]);

  //----------------------------------------------------------------
  // 2) element matrix: stationary terms of ion-transport equation
  //----------------------------------------------------------------

  // 2a)  element matrix: convective term
  my::CalcMatConv(emat,k,timefacfac,eps_[0],vmdc->Conv(),vmdc->SGConv());

  // 2b)  element matrix: diffusion term
  //      i)  constant diffusion coefficient
  my::CalcMatDiff(emat,k,timefacfac*epstort_[0],dmedc);

  //      ii) concentation depending diffusion coefficient
  //          (additional term for Newman material)
  if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
    CalcMatDiffCoeffLin(emat,k,timefacfac*epstort_[0],dmedc,vmdc->GradPhi(k));

  // 2c) electrical conduction term (transport equation)
  //
  //     mass transport equation:
  //
  //               |     diffusion term      | |     conduction term    |
  //               |                         | |                        |
  //      dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
  //
  //    equation for current (based on Newman):
  //
  //          | ohmic overpot.   |         concentration overpotential           |
  //          |                  |                                               |
  //      i = - kappa nabla phi  + RT/F kappa (thermfactor) f(t_k) nabla ln c_k

  // equation for current is inserted in the mass transport equation
  if (not cursolvar_)
  {
    //    i)  conduction term + ohmic overpotential
    //        (w_k, - t_k kappa nabla phi /(z_k F))
    CalcMatCondOhm(emat,k,timefacfac,vmdc->InvFVal(k),dmedc,vmdc->GradPot());

    //    ii) conduction term + concentration overpotential
    //        (w_k, - t_k RT/F kappa (thermfactor) f(t_k) nabla ln c_k /(z_k F))
    if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
      CalcMatCondConc(emat,k,timefacfac,vmdc->RTFFCVal(k),myelch::elchpara_->NewmanConstA(),myelch::elchpara_->NewmanConstB(),dmedc,vmdc->GradPhi(k),vmdc->ConIntInv());
  }
  // equation for current is solved independently
  else
  {
    // current term (with current as a solution variable)
    CalcMatCond(emat,k,timefacfac,vmdc->InvFVal(k),dmedc,vmdc->CurInt());

    // this coupling term cancels out for a 2 equation system
    if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
      CalcMatCondDiff(emat,k,timefacfac,vmdc->InvFVal(k),dmedc,vmdc->GradPhi());
  } // end if(not cursolvar_)

  //----------------------------------------------------------------
  // 3)   governing equation for the electric potential field and current
  //----------------------------------------------------------------
  // see function CalMatAndRhsOutsideScalarLoop()

  //-----------------------------------------------------------------------
  // 4) element right hand side vector (neg. residual of nonlinear problem)
  //-----------------------------------------------------------------------

  //adaption of rhs with respect to time integration
  //my::ComputeRhsInt(rhsint,eps_[0],eps_[0],hist);

  double rhs(0.0);
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
      const double vtrans = fac*eps_[0]*vmdc->ConInt(k);
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const int fvi = vi*my::numdofpernode_+k;

        erhs[fvi] -= vtrans*my::funct_(vi);
      }
    } // if(is_genalpha_)
  } // if (is_stationary_)

  // add RHS and history contribution
  my::CalcRHSHistAndSource(erhs,k,fac*eps_[0],rhsint);

  // 3a)  element rhsx: convective term
  my::CalcRHSConv(erhs,k,rhsfac*eps_[0],vmdc->ConvPhi(k));

  // 3b)  element rhs: diffusion term
  my::CalcRHSDiff(erhs,k,rhsfac*epstort_[0],dmedc,vmdc->GradPhi(k));

  // 3c) electrical conduction term (transport equation)
  //     equation for current is inserted in the mass transport equation
  //
  //     mass transport equation:
  //
  //               |     diffusion term      | |     conduction term    |
  //               |                         | |                        |
  //      dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
  //
  //    equation for current:
  //
  //          | ohmic overpot.   |         concentration overpotential           |
  //          |                  |                                               |
  //      i = - kappa nabla phi  + RT/F kappa (thermfactor) f(t_k) nabla ln c_k

  if(not cursolvar_)
  {
    CalcRhsCondOhm(erhs,k,rhsfac,vmdc->InvFVal(k),dmedc,vmdc->GradPot());

    // if(diffcondmat_==INPAR::ELCH::diffcondmat_ion): all terms cancel out
    if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
      CalcRhsCondConc(erhs,k,rhsfac,vmdc->RTFFCVal(k),myelch::elchpara_->NewmanConstA(),myelch::elchpara_->NewmanConstB(),dmedc,vmdc->GradPhi(k),vmdc->ConIntInv());
  }
  // equation for current is solved independently
  else
  {
    CalcRhsCond(erhs,k,rhsfac,vmdc->InvFVal(k),dmedc,vmdc->CurInt());

    if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
      CalcRhsCondDiff(erhs,k,rhsfac,dmedc,vmdc->GradPhi());
  }

  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs                           ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalMatAndRhsOutsideScalarLoop(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&                 emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&                 erhs,         //!< element rhs to calculate
    const double                             fac,          //!< domain-integration factor
    const double                             timefacfac,   //!< domain-integration factor times time-integration factor
    const double                             rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme           //!< diffusion manager
  )
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc = Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCond>(dme);

  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCond <my::nsd_,my::nen_> > vmdc
    = Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchDiffCond <my::nsd_,my::nen_> >(vm);

  //----------------------------------------------------------------
  // 3)   governing equation for the electric potential field
  //----------------------------------------------------------------
  //
  // mass transport equation:
  //            |     diffusion term      | |     conduction term    |
  //            |                         | |                        |
  //   dc_k/dt - nabla dot (D_k nabla c_k) + nabla dot (t_k i/(z_k F)
  //
  // equation for current:
  //   i = - kappa nabla phi + RT/F kappa (thermfactor) f(t_k) nabla ln c_k
  //
  // equation for potential:
  //
  // a) nabla cdot i = 0
  // b) ENC

  // // equation for current is NOT solved independently
  if (not cursolvar_)
  {
    // equation for current is inserted in the mass transport equation
    // 3a)  nabla dot i = 0
    if(equpot_==INPAR::ELCH::equpot_divi)
    {
      //  i)  ohmic overpotential (implemented after the scalar loop)
      //      (w_k, - kappa nabla phi)
      CalcMatPotEquDiviOhm(emat,timefacfac,vmdc->InvF(),dmedc,vmdc->GradPot());

      //
      CalcRhsPotEquDiviOhm(erhs,rhsfac,vmdc->InvF(),dmedc,vmdc->GradPot());

      //  i)  concentration  overpotential
      //      (w_k, RT/F kappa (thermfactor) f(t_k) nabla ln c_k)
      for (int k=0;k<my::numscal_;++k)
      {
        //
        CalcMatPotEquDiviConc(emat,k,timefacfac,vmdc->RTFFC(),vmdc->RTF(),vmdc->InvF(),myelch::elchpara_->NewmanConstA(),myelch::elchpara_->NewmanConstB(),dmedc,vmdc->GradPhi(k),vmdc->ConIntInv(k));

        //
        CalcRhsPotEquDiviConc(erhs,k,rhsfac,vmdc->RTF(),vmdc->InvFVal(),vmdc->RTFFC(),myelch::elchpara_->NewmanConstA(),myelch::elchpara_->NewmanConstB(),dmedc,vmdc->GradPhi(k),vmdc->ConIntInv(k));
      }
    }
    // 3b)  ENC
    else if(equpot_==INPAR::ELCH::equpot_enc)
    {
      for (int k=0;k<my::numscal_;++k)
      {
        //
        CalcMatPotEquENC(emat,k,fac,my::scatraparatimint_->AlphaF(),dmedc);

        //
        CalcRhsPotEquENC(erhs,k,fac,dmedc,vmdc->ConInt(k));
      }
    }
    else
      dserror("(div i, ENC) are the options available in the Diffusion-Conduction framework");
  }
  // equation for current is solved independently
  else
  {
    // (xsi_i,Di)
    CalcMatCurEquCur(emat, timefacfac,vmdc->InvF());

    // (nabla xsi, -D(kappa phi))
    CalcMatCurEquOhm(emat,timefacfac,vmdc->InvF(),dmedc,vmdc->GradPot());

    //
    CalcMatCurEquConc(emat,timefacfac,vmdc->RTF(),vmdc->RTFFC(),vmdc->InvFVal(),myelch::elchpara_->NewmanConstA(),myelch::elchpara_->NewmanConstB(),dmedc,vmdc->GradPhi(),vmdc->ConIntInv());

    //
    CalcRhsCurEquCur(erhs,rhsfac,vmdc->InvF(),dmedc,vmdc->CurInt());

    //
    CalcRhsCurEquOhm(erhs,rhsfac,vmdc->InvF(),dmedc,vmdc->GradPot());

    //
    CalcRhsCurEquConc(erhs,rhsfac,vmdc->RTF(),vmdc->InvFVal(),vmdc->RTFFC(),myelch::elchpara_->NewmanConstA(),myelch::elchpara_->NewmanConstB(),dmedc,vmdc->GradPhi(),vmdc->ConIntInv());

    if(equpot_==INPAR::ELCH::equpot_divi)
    {
      CalcMatPotEquDivi(emat,timefacfac,vmdc->InvF());

      CalcRhsPotEquDivi(erhs,rhsfac,vmdc->InvF(),vmdc->CurInt());
    }
    else if(equpot_==INPAR::ELCH::equpot_enc)
    {
      for (int k=0; k < my::numscal_; ++k)
      {
        //
        CalcMatPotEquENC(emat,k,fac,my::scatraparatimint_->AlphaF(),dmedc);

        //
        CalcRhsPotEquENC(erhs,k,fac,dmedc,vmdc->ConInt(k));
      }
    }
    else
      dserror("(div i, ENC) are the options available in the Diffusion-Conduction framework");
  }

  return;
}

/*----------------------------------------------------------------------------------*
|  CalcMat: Linearization of diffusion coefficient in diffusion term     ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatDiffCoeffLin(
  Epetra_SerialDenseMatrix&                       emat,         //!< element matrix to be filled
  const int                                       k,            //!< index of current scalar
  const double                                    timefacfac,   //!< domain-integration factor times time-integration factor
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,          //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&               gradphi       //!< gradient of concentration at GP
)
{
  //linearization of diffusion coefficient in the ionic diffusion term (transport equation)
  //
  // (nabla w, D(D(c)) nabla c)
  //
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double laplawfrhs_gradphi=0.0;
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi,vi);

      emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k)
        += timefacfac*epstort_[0]*dmedc->GetDerivIsoDiffCoef(k,k)*laplawfrhs_gradphi*my::funct_(ui);
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Conduction term with inserted current - ohmic overpotential  ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCondOhm(
  Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
  const int                               k,            //!< index of current scalar
  const double                            timefacfac,   //!< domain-integration factor times time-integration factor
  const double                            invfval,      //!< 1/(F z_k)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,          //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradpot       //!< gradient of potenial at GP
)
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
        timefacfac*epstort_[0]*dmedc->GetTransNum(k)*dmedc->GetCond()*invfval*laplawf;

      double laplawfrhs_gradpot(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);

      for(int iscal=0; iscal<my::numscal_;++iscal)
      {
        //linearization of the conductivity in the conduction term depending on the potential
        //
        // (grad w, t_k D(kappa(c))/(F z_k) grad phi)
        //
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
           += timefacfac*epstort_[0]*dmedc->GetTransNum(k)*invfval*dmedc->GetDerivCond(iscal)*my::funct_(ui)*laplawfrhs_gradpot;

        //linearization of the transference number in the conduction term depending on the potential
        //
        // (grad w, D(t_k(c)) kappa/(F z_k) grad phi)
        //
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
           += timefacfac*epstort_[0]*(dmedc->GetDerivTransNum(k,iscal))*my::funct_(ui)*invfval*dmedc->GetCond()*laplawfrhs_gradpot;
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Conduction term with inserted current - conc. overpotential  ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCondConc(
  Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
  const int                               k,              //!< index of current scalar
  const double                            timefacfac,     //!< domain-integration factor times time-integration factor
  const double                            rtffcval,       //!< RT/F^2/Newman_const_c/z_k
  const double                            newman_const_a, //!< Newman constant a
  const double                            newman_const_b, //!< Newman constant b
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,            //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradphi,        //!< gradient of concentration at GP
  const std::vector<double>&               conintinv       //!< inverted concentration at GP
)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf,ui,vi);

      for(int iscal=0;iscal<my::numscal_;++iscal)
      {
        // TODO: check linearization: iscal in machen Linearisierungen kommt mir komisch vor

        // linearization of conduction term depending on the concentration
        //
        // (grad w, RT/(z_k F^2) kappa thermfac f(t_+) D(grad ln c_k))
        //
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += timefacfac*epstort_[0]*rtffcval*dmedc->GetTransNum(k)*dmedc->GetCond()*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(iscal)))*conintinv[iscal]*laplawf;

        // linearization of conduction term depending on the concentration is implemented
        // only for one species
        // otherwise you would need a second loop over the all scalars
        double laplawfrhs_gradc(0.0);
        my::GetLaplacianWeakFormRHS(laplawfrhs_gradc,gradphi,vi);

        // Linearization wrt ln c
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += -timefacfac*epstort_[0]*rtffcval*dmedc->GetTransNum(k)*dmedc->GetCond()*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(iscal)))*conintinv[iscal]*conintinv[iscal]*laplawfrhs_gradc*my::funct_(ui);

        // Linearization wrt kappa
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += timefacfac*epstort_[0]*rtffcval*dmedc->GetTransNum(k)*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradc*dmedc->GetDerivCond(iscal)*my::funct_(ui);

        // Linearization wrt transference number 1
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += timefacfac*epstort_[0]*rtffcval*dmedc->GetCond()*conintinv[iscal]*laplawfrhs_gradc*dmedc->GetThermFac()*(newman_const_a+newman_const_b*dmedc->GetTransNum(iscal))*dmedc->GetDerivTransNum(iscal,iscal)*my::funct_(ui);

        // Linearization wrt transference number 2
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += timefacfac*epstort_[0]*rtffcval*dmedc->GetCond()*dmedc->GetThermFac()*conintinv[iscal]*laplawfrhs_gradc*dmedc->GetTransNum(iscal)*newman_const_b*dmedc->GetDerivTransNum(iscal,iscal)*my::funct_(ui);

        // Linearization wrt thermodynamic factor
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += timefacfac*epstort_[0]*rtffcval*dmedc->GetTransNum(k)*dmedc->GetCond()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradc*dmedc->GetDerivThermFac(iscal)*my::funct_(ui);
      }
    } // for ui
  } // for vi

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Conduction term without inserted current                     ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCond(
  Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
  const int                               k,            //!< index of current scalar
  const double                            timefacfac,   //!< domain-integration factor times time-integration factor
  const double                            invfval,      //!< 1/(F z_k)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,          //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        curint        //!< current at GP
)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      for (int idim=0; idim<my::nsd_; ++idim)
      {
        const int fui = ui*my::numdofpernode_+(my::numscal_+1)+idim;

        //linearization of conduction term depending on current flow
        //
        // (grad w, t_k/(F z_k) Di)
        //
        emat(fvi,fui) += -timefacfac*my::derxy_(idim,vi)*dmedc->GetTransNum(k)*invfval*my::funct_(ui);

        //linearization of transference number in conduction term depending on current flow
        //
        // (grad w, Dt_k(c)/(F z_k) i)
        //
        for (int iscal=0; iscal<my::numscal_; ++iscal)
          emat(fvi,ui*my::numdofpernode_+iscal)
            += -timefacfac*my::derxy_(idim,vi)*(dmedc->GetDerivTransNum(k,iscal))*my::funct_(ui)*invfval*curint(idim);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Additional diffusion term without inserted current           ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCondDiff(
  Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
  const int                               k,            //!< index of current scalar
  const double                            timefacfac,   //!< domain-integration factor times time-integration factor
  const double                            invfval,      //!< 1/(F z_k)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,          //!< diffusion manager
  const std::vector<LINALG::Matrix<my::nsd_,1> >&   gradphi       //!< gradient of concentration at GP
)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      // compute once, reuse below!
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf,ui,vi);

      for (int iscal=0; iscal < my::numscal_; ++iscal)
      {
        // formulation a): plain ionic diffusion coefficients without using ENC
        //
        // (grad w, t_k/z_k*sum_i(D_i grad Dc))
        //
        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += -timefacfac*epstort_[0]*dmedc->GetTransNum(k)/dmedc->GetValence(k)*dmedc->GetValence(iscal)*dmedc->GetIsotropicDiff(iscal)*laplawf;
      }

      //linearization of transference number in the coupling term (transport equation)
      //
      // (grad w, Dt_k(c)/z_k (Sum_i z_i D_i grad c_i))
      //
      for(int iscal=0; iscal<my::numscal_; ++iscal)
      {
        double term_vi = 0.0;
        for(int iscal2=0; iscal2<my::numscal_; ++iscal2)
        {
          double laplawfrhs_gradphi=0.0;
          my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[iscal2],vi);

          term_vi += dmedc->GetValence(iscal2)*dmedc->GetIsotropicDiff(iscal2)*laplawfrhs_gradphi;
        }

        emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
          += -timefacfac*epstort_[0]*(dmedc->GetDerivTransNum(k,iscal))*my::funct_(ui)/dmedc->GetValence(k)*term_vi;
      } // for(iscal)

      // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
      //                  -> not implemented yet

      // formulation c):  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i) replaced by sum (t_i/c_i nabla c_i)
      //                  -> not activated
      //                  -> linearization is missing
      // emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+iscal)
      //    += timefacfac*epstort_[0]/frt/faraday/dmedc->GetValence(k)*dmedc->GetTransNum(k)*dmedc->GetCond()*dmedc->GetTransNum(iscal)/dmedc->GetValence(iscal]/conint[iscal]*laplawf;
    } // end for ui
  } // end for vi

  return;
}

/*--------------------------------------------------------------------------------------*
|  CalcMat: Potential equation div i inserted current - ohmic overpotential  ehrl  02/14|
*---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatPotEquDiviOhm(
  Epetra_SerialDenseMatrix&               emat,         //!< element matrix to be filled
  const double                            timefacfac,   //!< domain-integration factor times time-integration factor
  const double                            invf,         //!< 1/F
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,          //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradpot       //!< gradient of potenial at GP
)
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
      emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += timefacfac*epstort_[0]*invf*dmedc->GetCond()*laplawf;

      for(int iscal=0;iscal<my::numscal_;++iscal)
      {
        double laplawfrhs_gradpot(0.0);
        my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);

        // linearization of the ohmic term wrt conductivity
        //
        // (grad w, 1/F kappa D(grad pot))
        //
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal)
          += timefacfac*epstort_[0]*invf*dmedc->GetDerivCond(iscal)*my::funct_(ui)*laplawfrhs_gradpot;
      }
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------*
|  CalcMat: Potential equation div i inserted current - conc. overpotential   ehrl  02/14|
*----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatPotEquDiviConc(
  Epetra_SerialDenseMatrix&                       emat,           //!< element matrix to be filled
  const int                                       k,              //!< index of current scalar
  const double                                    timefacfac,     //!< domain-integration factor times time-integration factor
  const double                                    rtffc,          //!< RT/(F^2 Newman_const_c)
  const double                                    rtf,            //!< RT/F
  const double                                    invf,           //!< 1/F
  const double                                    newman_const_a, //!< Newman constant a
  const double                                    newman_const_b, //!< Newman constant b
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>& dmedc,          //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&               gradphi,        //!< gradient of concentration at GP
  const double                                    conintinv       //!< inverted concentration at GP
)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf,ui,vi);

      if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
      {
        // linearization of the diffusion overpotential term
         //
         // (grad w, RT/F^2 kappa (thermfactor) f(t_k) 1/c_k D nabla c_k)
         //
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
          += timefacfac*epstort_[0]*rtffc*dmedc->GetCond()*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv*laplawf;

        // linearization of conduction term depending on the concentration is implemented
        // only for one species
        // otherwise you would need a second loop over the all scalars
        double laplawfrhs_gradphi(0.0);
        my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi,vi);

        // Linearization wrt ln c
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
          += -timefacfac*epstort_[0]*rtffc*dmedc->GetCond()*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv*conintinv*laplawfrhs_gradphi*my::funct_(ui);

        // Linearization wrt kappa
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
          += timefacfac*epstort_[0]*rtffc*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv*laplawfrhs_gradphi*dmedc->GetDerivCond(k)*my::funct_(ui);

        // Linearization wrt transference number
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
          += timefacfac*epstort_[0]*rtffc*dmedc->GetCond()*dmedc->GetThermFac()*conintinv*laplawfrhs_gradphi*newman_const_b*dmedc->GetDerivTransNum(k,k)*my::funct_(ui);

        // Linearization wrt thermodynamic factor
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
          += timefacfac*epstort_[0]*rtffc*dmedc->GetCond()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv*laplawfrhs_gradphi*dmedc->GetDerivThermFac(k)*my::funct_(ui);
      }
      else if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
      {
        if(myelch::elchpara_->DiffusionCoeffBased()==true)
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
            += timefacfac*epstort_[0]*dmedc->GetValence(k)*dmedc->GetIsotropicDiff(k)*laplawf;
        else
        {
          // TODO: Linearization of transference number, conductivity, ... is still missing
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)
            += timefacfac*epstort_[0]*rtf*invf/dmedc->GetValence(k)*dmedc->GetCond()*dmedc->GetTransNum(k)*conintinv*laplawf;
        }
      }
      else
        dserror("Diffusion-Conduction material is not specified");
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Potential equation ENC                                       ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatPotEquENC(
  Epetra_SerialDenseMatrix&               emat,   //!< element matrix to be filled
  const int                               k,      //!< index of current scalar
  const double                            fac,    //!< domain-integration factor
  const double                            alphaf, //!< time factor for ENC
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc     //!< diffusion manager

)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      //linearization of the transference number in the conduction term (transport equation)
      //
      // (w, sum(z_k c_k))
      //
      emat(vi*my::numdofpernode_+my::numscal_, ui*my::numdofpernode_+k) += alphaf*dmedc->GetValence(k)*fac*my::funct_(vi)*my::funct_(ui);
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Potential equation div i without inserted current            ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatPotEquDivi(
  Epetra_SerialDenseMatrix&   emat,       //!< element matrix to be filled
  const double                timefacfac, //!< domain-integration factor times time-integration factor
  const double               invf        //!< 1/F
)
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

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Current equation current                                     ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCurEquCur(
  Epetra_SerialDenseMatrix&   emat,       //!< element matrix to be filled
  const double                timefacfac, //!< domain-integration factor times time-integration factor
  const double                invf        //!< 1/F
)
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

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Current equation ohmic overpotential                         ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCurEquOhm(
  Epetra_SerialDenseMatrix&               emat,       //!< element matrix to be filled
  const double                            timefacfac, //!< domain-integration factor times time-integration factor
  const double                           invf,       //!< 1/F
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,        //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradpot     //!< gradient of potenial at GP
)
{
  // (v, kappa grad phi)
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      for (int idim=0; idim<my::nsd_; ++idim)
      {
        const int fvi = vi*my::numdofpernode_+(my::numscal_+1)+idim;
        const int fui = ui*my::numdofpernode_+my::numscal_;

        emat(fvi,fui) += timefacfac*invf*epstort_[0]*my::funct_(vi)*dmedc->GetCond()*my::derxy_(idim,ui);

        //linearization of conductivity in the ohmic resistance term (current equation)
        //
        // (w, D(kappa(c)) grad phi)
        //
        for(int k=0;k<my::numscal_;++k)
          emat(fvi,ui*my::numdofpernode_+k) += timefacfac*invf*epstort_[0]*my::funct_(vi)*dmedc->GetDerivCond(k)*my::funct_(ui)*gradpot(idim);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Current equation concentration overpotential                 ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcMatCurEquConc(
  Epetra_SerialDenseMatrix&               emat,           //!< element matrix to be filled
  const double                            timefacfac,     //!< domain-integration factor times time-integration factor
  const double                            rtf,            //!< RT/F
  const double                            rtffc,          //!< RT/(F^2 Newman_const_c)
  const std::vector <double>&              invfval,        //!< 1/(F z_k)
  const double                           newman_const_a, //!< Newman constant a
  const double                          newman_const_b, //!< Newman constant b
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,            //!< diffusion manager
  const std::vector< LINALG::Matrix<my::nsd_,1> >& gradphi, //!< gradient of concentration at GP
  const std::vector<double>&               conintinv       //!< inverted concentration at GP
)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      // diffusive term
      // (grad w, D grad c)
      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        for (int k=0; k < my::numscal_; ++k)
        {
          if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
          {
            emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+k)
              += timefacfac*rtffc*epstort_[0]*my::funct_(vi)*dmedc->GetCond()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv[k]*my::derxy_(idim,ui);

            // TODO: Linearization of coupling terms is still missing

            //  //linearization of coupling term in the current equation is still missing
            //  emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+k)
            //    += timefacfac*my::funct_(vi)/frt_*dmedc->GetCond()*dmedc->GetTransNum(k)/dmedc->GetValence(k]/((-1)*conint_[k]*conint_[k])*my::funct_(ui)*gradphicoupling_[k](idim);
            //
            //  for(int k2=0; k2<numscal_;++k2)
            //  {
            //    //Check if necessary??
            //    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+k)
            //       += timefacfac*my::funct_(vi)/frt_*dmedc->GetDerivCond(k)*my::funct_(ui)*trans_[k2]*my::funct_(ui)/dmedc->GetValence(k2]/conint_[k2]*gradphicoupling_[k2](idim);
            //
            //
            //    emat(vi*numdofpernode_+(numscal_+1)+idim,ui*numdofpernode_+k)
            //      += timefacfac*my::funct_(vi)/frt_*dmedc->GetCond()*(transderiv_[k2])[k]*my::funct_(ui)/dmedc->GetValence(k2]/conint_[k2]*gradphicoupling_[k2](idim);
            //  }
          }
          else if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
          {
            if(myelch::elchpara_->DiffusionCoeffBased()==true)
              emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+k)
                += timefacfac*epstort_[0]*my::funct_(vi)*dmedc->GetValence(k)*dmedc->GetIsotropicDiff(k)*my::derxy_(idim,ui);
            else
            {
              // linarization wrt nabla c_k
              emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+k)
                += timefacfac*epstort_[0]*invfval[k]*rtf*my::funct_(vi)*dmedc->GetCond()*dmedc->GetTransNum(k)*conintinv[k]*my::derxy_(idim,ui);

              // linearization wrt 1/c_k
              emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+k)
                 += -timefacfac*epstort_[0]*rtf*invfval[k]*dmedc->GetCond()*dmedc->GetTransNum(k)*conintinv[k]*conintinv[k]*my::funct_(vi)*(gradphi[k])(idim)*my::funct_(ui);

              // linearization wrt kappa
              double term_vi = 0.0;
              for (int iscal=0; iscal < my::numscal_; ++iscal)
              {
                term_vi += invfval[iscal]*dmedc->GetTransNum(iscal)*conintinv[iscal]*my::funct_(vi)*(gradphi[iscal])(idim)*my::funct_(ui);
              }
              emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+k)
                += timefacfac*epstort_[0]*dmedc->GetDerivCond(k)*rtf*term_vi;

              // linearization wrt transference number
              for (int iscal=0; iscal < my::numscal_; ++iscal)
              {
                emat(vi*my::numdofpernode_+(my::numscal_+1)+idim,ui*my::numdofpernode_+iscal)
                  += timefacfac*epstort_[0]*rtf*invfval[k]*dmedc->GetCond()*(dmedc->GetDerivTransNum(k,iscal))*conintinv[k]*my::funct_(vi)*(gradphi[k])(idim)*my::funct_(ui);
              }
            }
          }
          else
            dserror("Diffusion-Conduction material is not specified");
        }
      }
    } // for ui
  } // for vi

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Conduction term with inserted current - ohmic overpotential    ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCondOhm(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const int                               k,        //!< index of current scalar
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                            invfval,  //!< 1/(F z_k)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> & dmedc,      //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradpot   //!< gradient of potenial at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    // diffusive term
    double laplawfrhs_gradpot=0.0;
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);
    erhs[vi*my::numdofpernode_+k]-= rhsfac*epstort_[0]*dmedc->GetTransNum(k)*dmedc->GetCond()*invfval*laplawfrhs_gradpot;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Conduction term with inserted current - conc. overpotential    ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCondConc(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const int                               k,        //!< index of current scalar
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                            rtffcval, //!< RT/(F^2 Newman_const_c z_k)
  const double                            newman_const_a, //!< Newman constant a
  const double                            newman_const_b, //!< Newman constant b
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradphi,  //!< gradient of concentration at GP
  const std::vector<double>&               conintinv //!< inverted concentration at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for(int iscal=0; iscal<my::numscal_;++iscal)
    {
      // diffusive term second
      double laplawfrhs_gradphi(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi,vi); // compute once, reuse below!

      // formulation a): plain ionic diffusion coefficients without using ENC
      //
      // (grad w, sum(D_i grad Dc))
      erhs[vi*my::numdofpernode_+k]
        -= rhsfac*epstort_[0]*rtffcval*dmedc->GetTransNum(k)*dmedc->GetCond()*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(iscal)))*conintinv[iscal]*laplawfrhs_gradphi;

      // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
      //                  -> not implemented yet
    }
  }

  return;
}

/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Conduction term without inserted current                       ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCond(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const int                               k,        //!< index of current scalar
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                            invfval,  //!< 1/(F z_k)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        curint    //!< current at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    double laplawfrhs_cur=0.0;
    my::GetLaplacianWeakFormRHS(laplawfrhs_cur,curint,vi);

    erhs[vi*my::numdofpernode_+k]-= -rhsfac*dmedc->GetTransNum(k)*invfval*laplawfrhs_cur;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Additional diffusion term without inserted current             ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCondDiff(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const int                               k,        //!< index of current scalar
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
  const std::vector< LINALG::Matrix<my::nsd_,1> >&  gradphi   //!< gradient of concentration at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int iscal=0; iscal <my::numscal_; ++iscal)
    {
      // diffusive term second
      double laplawfrhs_gradphi(0.0);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi[iscal],vi); // compute once, reuse below!

      // formulation a:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i)
      erhs[vi*my::numdofpernode_+k] -= - rhsfac*epstort_[0]*dmedc->GetTransNum(k)*dmedc->GetValence(iscal)/dmedc->GetValence(k)*dmedc->GetIsotropicDiff(iscal)*laplawfrhs_gradphi;

      // formulation b):  plain ionic diffusion coefficients: one species eleminated by ENC
      //                  -> not implemented yet
      // formulation c:  plain ionic diffusion coefficients: sum (z_i D_i nabla c_i) replaced by sum (t_i/c_i nabla c_i)
      //                  -> not activated
      // erhs[fvi]
      //   -= rhsfac*epstort_[0]/frt/faraday/dmedc->GetValence(k)*dmedc->GetTransNum(k)*dmedc->GetCond()*dmedc->GetTransNum(iscal)/dmedc->GetValence(iscal]/conint[iscal]*laplawf2;
    }
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 | CalcRhs: Potential equation div i inserted current - ohmic overpotential  ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsPotEquDiviOhm(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                            invf,     //!< 1/F
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> & dmedc,      //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradpot   //!< gradient of potenial at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    double laplawfrhs_gradpot(0.0);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi); // compute once, reuse below!

    erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*epstort_[0]*invf*dmedc->GetCond()*laplawfrhs_gradpot;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |CalcRhs: Potential equation div i inserted current - conc. overpotential  ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsPotEquDiviConc(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const int                               k,        //!< index of current scalar
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                            rtf,      //!< RT/F
  const std::vector<double>&               invfval,  //!< 1/(F z_k)
  const double                            rtffc, //!< RT/(F^2 Newman_const_c)
  const double                            newman_const_a, //!< Newman constant a
  const double                            newman_const_b, //!< Newman constant b
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradphi,  //!< gradient of concentration at GP
  const double                            conintinv //!< inverted concentration at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
    {
      // diffusive term second
      double laplawf2(0.0);
      my::GetLaplacianWeakFormRHS(laplawf2,gradphi,vi); // compute once, reuse below!

      if(myelch::elchpara_->DiffusionCoeffBased()==true)
        erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*epstort_[0]*dmedc->GetValence(k)*dmedc->GetIsotropicDiff(k)*laplawf2;
      else
        erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*epstort_[0]*rtf*invfval[k]*dmedc->GetCond()*dmedc->GetTransNum(k)*conintinv*laplawf2;
    }
    // thermodynamic factor only implemented for Newman
    else if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
    {
      // diffusive term second
      double laplawf2(0.0);
      my::GetLaplacianWeakFormRHS(laplawf2,gradphi,vi); // compute once, reuse below!

      erhs[vi*my::numdofpernode_+my::numscal_]
        -= rhsfac*epstort_[0]*rtffc*dmedc->GetCond()*dmedc->GetThermFac()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv*laplawf2;
    }
    else
      dserror("Diffusion-Conduction material is not specified");
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Potential equation ENC                                         ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsPotEquENC(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const int                               k,        //!< index of current scalar
  const double                            fac,   //!< time-integration factor for rhs times domain-integration factor
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
  const double                            conint    //!< concentration at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    // electroneutrality condition
    // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
    erhs[vi*my::numdofpernode_+my::numscal_] -= dmedc->GetValence(k)*fac*my::funct_(vi)*conint;
  }

  return;
}

/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Potential equation divi without inserted current               ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsPotEquDivi(
  Epetra_SerialDenseVector&   erhs,     //!< element vector to be filled
  const double                rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                invf,     //!< 1/F
  const LINALG::Matrix<my::nsd_,1>&  curint    //!< current at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    double laplawf=0.0;
    // version a: (grad phi,  Di)
    my::GetLaplacianWeakFormRHS(laplawf,curint,vi);
    erhs[vi*my::numdofpernode_+my::numscal_]-= -rhsfac*invf*laplawf;
    // version b: (phi, div Di) -> not partially integrated
    //erhs[vi*my::numdofpernode_+my::numscal_] -= rhsfac*my::funct_(vi)*divi;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Current equation - current                                     ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCurEquCur(
    Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
    const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
    const double                            invf,     //!< 1/F
    Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
    const LINALG::Matrix<my::nsd_,1> &       curint    //!< current at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      // (v, i)
      erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]-=rhsfac*invf*my::funct_(vi)*curint(idim);
    }
  }

  return;
}

/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Current equation - ohmic overpotential                         ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCurEquOhm(
  Epetra_SerialDenseVector&               erhs,     //!< element vector to be filled
  const double                            rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                            invf,     //!< 1/F
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&  dmedc,      //!< diffusion manager
  const LINALG::Matrix<my::nsd_,1>&        gradpot   //!< gradient of potenial at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      // (v, kappa grad phi)
      erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]-=rhsfac*invf*epstort_[0]*my::funct_(vi)*dmedc->GetCond()*gradpot(idim);
    }
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Current equation - concentration overpotential                 ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcRhsCurEquConc(
  Epetra_SerialDenseVector&                 erhs,     //!< element vector to be filled
  const double                              rhsfac,   //!< time-integration factor for rhs times domain-integration factor
  const double                              rtf,      //!< RT/F
  const std::vector<double>&                 invfval,  //!< 1/(F z_k)
  const double                              rtffc,    //!< RT/(F^2 Newman_const_c)
  const double                              newman_const_a, //!< Newman constant a
  const double                              newman_const_b, //!< Newman constant b
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&    dmedc,      //!< diffusion manager
  const std::vector<LINALG::Matrix<my::nsd_,1> >&  gradphi,  //!< vector of gradient of concentration at GP
  const std::vector<double>&                       conintinv //!< inverted concentration at GP
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int idim = 0; idim <my::nsd_; ++idim)
    {
      for (int k=0; k < my::numscal_; ++k)
      {
        if(diffcondmat_==INPAR::ELCH::diffcondmat_newman)
        {
          erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]
            -= rhsfac*epstort_[0]*my::funct_(vi)*rtffc*dmedc->GetCond()*(newman_const_a+(newman_const_b*dmedc->GetTransNum(k)))*conintinv[k]*gradphi[k](idim);
        }
        else if(diffcondmat_==INPAR::ELCH::diffcondmat_ion)
        {
          if(myelch::elchpara_->DiffusionCoeffBased()==true)
            erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim] -= rhsfac*epstort_[0]*my::funct_(vi)*dmedc->GetValence(k)*dmedc->GetIsotropicDiff(k)*gradphi[k](idim);
          else
          {
            erhs[vi*my::numdofpernode_+(my::numscal_+1)+idim]
                 -= rhsfac*epstort_[0]*my::funct_(vi)*rtf*dmedc->GetCond()*invfval[k]*dmedc->GetTransNum(k)*conintinv[k]*gradphi[k](idim);
          }
        }
        else
          dserror("Diffusion-Conduction material is not specified");
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
|  Correct sysmat for fluxes accros DC                       ehrl  02/14|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CorrectionForFluxAccrosDC(
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
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc = Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCond>(my::diffmanager_);

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
        erhs[vi*my::numdofpernode_+my::numscal_] += dmedc->GetValence(k)*(-val);
        // corresponding linearization
        for (int ui=0; ui<my::nen_; ++ui)
        {
          val = emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)+=dmedc->GetValence(k)*(-val);
          val = emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_)+=dmedc->GetValence(k)*(-val);
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
          erhs[vi*my::numdofpernode_+k] += 1.0/dmedc->GetValence(k)*(-val);
          // corresponding linearization
          for (int ui=0; ui<my::nen_; ++ui)
          {
            val = emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k);
            emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += 1.0/dmedc->GetValence(k)*(-val);
            val = emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_);
            emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=1.0/dmedc->GetValence(k)*(-val);
          }
        }
      }
    }  // for k
  }

//   // for concentrated solution theory (using div i as closing term for the potential)
//   // additional flux terms / currents across Dirichlet boundaries
//   if(myelch::elchpara_->ElchType()==INPAR::ELCH::elchtype_diffcond and
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
//             elevec1_epetra[vi*my::numdofpernode_+my::numscal_] += dmedc_->GetValence(k)*(-val);
//             // corresponding linearization
//             for (int ui=0; ui<my::nen_; ++ui)
//             {
//               val = elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
//               elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k)+=dmedc_->GetValence(k)*(-val);
//               val = elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
//               elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_)+=dmedc_->GetValence(k)*(-val);
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
//           //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<dmedc_->GetValence(k)<<std::endl;
//           const int fvi = vi*my::numdofpernode_+my::numscal_;
//           // We use the fact, that the rhs vector value for boundary nodes
//           // is equivalent to the integrated negative normal flux
//           // due to diffusion and migration
//
//           // scaling of div i results in a matrix with better condition number
//           val = elevec1_epetra[fvi];
//           elevec1_epetra[vi*my::numdofpernode_+k] += 1.0/dmedc_->GetValence(k)*(-val);
//           // corresponding linearization
//           for (int ui=0; ui<my::nen_; ++ui)
//           {
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += 1.0/dmedc_->GetValence(k)*(-val);
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=1.0/dmedc_->GetValence(k)*(-val);
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
//           //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<dmedc_->GetValence(k)<<std::endl;
//           const int fvi = vi*my::numdofpernode_+my::numscal_;
//           // We use the fact, that the rhs vector value for boundary nodes
//           // is equivalent to the integrated negative normal flux
//           // due to diffusion and migration
//
//           // scaling of div i results in a matrix with better condition number
//           val = elevec1_epetra[fvi];
//           elevec1_epetra[vi*my::numdofpernode_+k] += 1.0/dmedc_->GetValence(k)*(-val);
//           // corresponding linearization
//           for (int ui=0; ui<my::nen_; ++ui)
//           {
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += 1.0/dmedc_->GetValence(k)*(-val);
//             val = elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_);
//             elemat1_epetra(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_) +=1.0/dmedc_->GetValence(k)*(-val);
//           }
//         }
//       } // if Dirichlet at node vi
//     } // for vi
//   }

  return;
}


/*----------------------------------------------------------------------*
 |  set formulation-specific variables                        ehrl 01/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::SetFormulationSpecificInternalVariables(
  Teuchos::RCP<ScaTraEleDiffManagerElch>&                                   dme,
  Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >&  vm
)
{
  // dynamic cast to elch diffusion-conduction-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc = Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCond>(dme);

  // dynamic cast to elch diffusion conduction-specific diffusion manager
  Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCond <my::nsd_,my::nen_> > vmdc =
      Teuchos::rcp_static_cast< ScaTraEleInternalVariableManagerElchDiffCond <my::nsd_,my::nen_> >(vm);

  vmdc->SetInternalVariablesElchDiffCond(dmedc,my::funct_,my::derxy_,ecurnp_);

  return;
}


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                     ehrl 01/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::GetMaterialParams(
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
    //TODO: ELCH: dynamic cast ist relative teuer
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
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::Materials(
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
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc = Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCond>(diffmanager);

  if(material->MaterialType() == INPAR::MAT::m_elchphase)
    MatPhase(material,k,densn,densnp,densam,dmedc,reamanager,visc,iquad);
  else if(material->MaterialType() == INPAR::MAT::m_newman)
  {
    MatNewman(material,k,densn,densnp,densam,dmedc,reamanager,visc,iquad);
    diffcondmat_ = INPAR::ELCH::diffcondmat_newman;
  }
  else if(material->MaterialType() == INPAR::MAT::m_ion)
  {
    myelch::MatIon(material,k,densn,densnp,densam,dmedc,reamanager,visc,iquad);

    // set flag for material type
    diffcondmat_ = INPAR::ELCH::diffcondmat_ion;

    // Loop over materials is finished - now all material parameter are set
    // calculation of conductivity and transference number based on diffusion coefficient and valence
    // for mat_ion
    if(k==(my::numscal_-1))
    {
      std::vector<double> conint(my::numscal_);
      for (int k = 0;k<my::numscal_;++k)
        conint[k] = my::funct_.Dot(my::ephinp_[k]);

      dmedc->CalcConductivity(my::numscal_,INPAR::ELCH::faraday_const*myelch::elchpara_->FRT(),conint);
      dmedc->CalcTransNum(my::numscal_,conint);
    }
  }
  else dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material Phase                                           ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::MatPhase(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               iphase,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>  dmedc,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
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
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::MatNewman(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>  dmedc,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
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

  // concentration of species k at the integration point
  const double conint = my::funct_.Dot(my::ephinp_[k]);

  // valence of ionic species
  dmedc->SetValence(actmat->Valence(),k);

  // concentration depending diffusion coefficient
  dmedc->SetIsotropicDiff(actmat->ComputeDiffusionCoefficient(conint),k);
  // derivation of concentration depending diffusion coefficient wrt all ionic species
  dmedc->SetDerivIsoDiffCoef(actmat->ComputeFirstDerivDiffCoeff(conint),k,k);

  // concentration depending transference number
  dmedc->SetTransNum(actmat->ComputeTransferenceNumber(conint),k);
  // derivation of concentration depending transference number wrt all ionic species
  dmedc->SetDerivTransNum(actmat->ComputeFirstDerivTrans(conint),k,k);

  // thermodynamic factor of electrolyte solution
  dmedc->SetThermFac(actmat->ComputeThermFac(conint));
  // derivative of conductivity with respect to concentrations
  dmedc->SetDerivThermFac(actmat->ComputeFirstDerivThermFac(conint),0);

  // conductivity and first derivative can maximally depend on one concentration
  // since time curve is used as input routine
  // conductivity of electrolyte solution
  dmedc->SetCond(actmat->ComputeConductivity(conint));
  // derivative of conductivity with respect to concentrations
  dmedc->SetDerivCond(actmat->ComputeFirstDerivCond(conint),0);

  return;
}


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs27>;



