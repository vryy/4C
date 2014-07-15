/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_lsreinit.cpp

\brief evaluation of scatra elements for reinitialization equation

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_lsreinit.H"

#include "scatra_ele.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_parameter_lsreinit.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on

#define USE_PHIN_FOR_VEL
//#define MODIFIED_EQ

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype> * DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcLsReinit<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcLsReinit<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::ScaTraEleCalcLsReinit(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    ephizero_(my::numscal_)  // size of vector
{
  // set appropriate parameter list
  my::scatrapara_ = DRT::ELEMENTS::ScaTraEleParameterLsReinit::Instance();
  // set appropriate diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerLsReinit<my::nsd_>(my::numscal_));
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Evaluate(
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

  // clear all unused variables
//  my::ephin_.clear();
  my::edispnp_.Clear();
  my::weights_.Clear();
  my::myknots_.clear();
  my::evelnp_.Clear();
  my::eaccnp_.Clear();
  my::eprenp_.Clear();
  my::bodyforce_.clear();

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  Teuchos::RCP<const Epetra_Vector> phizero = discretization.GetState("phizero");
  if (hist==Teuchos::null || phinp==Teuchos::null || phin==Teuchos::null || phizero==Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'/'phin' and/or 'phizero'");
  std::vector<double> myhist(lm.size());
  std::vector<double> myphinp(lm.size());
  std::vector<double> myphin(lm.size());
  std::vector<double> myphizero(lm.size());
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
  DRT::UTILS::ExtractMyValues(*phin,myphin,lm);
  DRT::UTILS::ExtractMyValues(*phizero,myphizero,lm);

  // fill all element arrays
  for (int i=0;i<my::nen_;++i)
  {
    for (int k = 0; k< my::numscal_; ++k)
    {
      // split for each transported scalar, insert into element arrays
      my::ephinp_[k](i,0) = myphinp[k+(i*my::numdofpernode_)];
      my::ephin_[k](i,0) = myphin[k+(i*my::numdofpernode_)];
      ephizero_[k](i,0) = myphizero[k+(i*my::numdofpernode_)];
      // the history vector contains information of time step t_n
      my::ehist_[k](i,0) = myhist[k+(i*my::numdofpernode_)];
    } // for k
  } // for i

  if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->UseProjectedVel())
  {
    // get velocity at nodes (pre-computed via L2 projection)
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("reinitialization velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,my::econvelnp_,velocity,my::nsd_);
  }

  // calculate element coefficient matrix and rhs
  Sysmat(
    elemat1_epetra,
    elevec1_epetra);

  return 0;
}


#if 0
/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)             rasthofer 12/13 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Sysmat(
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs ///< element rhs to calculate
  )
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0);

  //----------------------------------------------------------------------
  // calculation of characteristic element length
  //----------------------------------------------------------------------

  // get gradient of initial phi at element center
  LINALG::Matrix<my::nsd_,1> gradphizero(true);
  gradphizero.Multiply(my::derxy_,ephizero_[0]);

  // get characteristic element length
  const double charelelength = CalcCharEleLengthReinit(vol,gradphizero);

  //----------------------------------------------------------------------
  // calculation of stabilization parameter at element center
  //----------------------------------------------------------------------

  // the stabilization parameter
  double tau = 0.0;

  if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
  {
    if (not my::scatrapara_->TauGP())
    {
      // gradient of current scalar value at element center
      LINALG::Matrix<my::nsd_,1> gradphi(true);
      gradphi.Multiply(my::derxy_,my::ephinp_[0]);
      // get norm
      const double gradphi_norm = gradphi.Norm2();
      // get sign function
      double signphi = 0.0;
      // initial phi at element center
      double phizero = 0.0;
      phizero = my::funct_.Dot(ephizero_[0]);
      // current phi at element center
      double phi = 0.0;
      phi = my::funct_.Dot(my::ephinp_[0]);
      SignFunction(signphi,charelelength,phizero,gradphizero,phi,gradphi);

      // get velocity at element center
//      LINALG::Matrix<my::nsd_,1> convelint(true);
//      if (gradphi_norm>1e-8)
//        convelint.Update(signphi/gradphi_norm,gradphi);
      //TODO:
      // get velocity at integration point
      LINALG::Matrix<my::nsd_,1> convelint(true);
      convelint.Multiply(my::evelnp_,my::funct_);
      // otherwise gradphi is almost zero and we keep a zero velocity

      // calculation of stabilization parameter at element center
      my::CalcTau(tau,0.0,0.0,1.0,convelint,vol,0);
    }
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //--------------------------------------------------------------------
    // set all Gauss point quantities
    //--------------------------------------------------------------------

    // gradient of current scalar value Gauss point
    LINALG::Matrix<my::nsd_,1> gradphi(true);
    gradphi.Multiply(my::derxy_,my::ephinp_[0]);
    // get norm
    const double gradphi_norm = gradphi.Norm2();
    // get sign function
    double signphi = 0.0;
    // initial phi at Gauss point
    double phizero = 0.0;
    phizero = my::funct_.Dot(ephizero_[0]);
    gradphizero.Clear();
    gradphizero.Multiply(my::derxy_,ephizero_[0]);
//    std::cout << phizero << std::endl;
    // current phi at Gauss point
    double phi = 0.0;
    phi = my::funct_.Dot(my::ephinp_[0]);
    SignFunction(signphi,charelelength,phizero,gradphizero,phi,gradphi);

    // get velocity at element center
    LINALG::Matrix<my::nsd_,1> convelint(true);
    //if (gradphi_norm>1e-8)
      //convelint.Update(signphi/gradphi_norm,gradphi);
//    LINALG::Matrix<my::nsd_,1> mygradphin(true);
//    mygradphin.Multiply(my::derxy_,my::ephin_[0]);
//    double mynorm = mygradphin.Norm2();
//    if (mynorm>1e-8)
//    convelint.Update(signphi/mynorm,mygradphin);
    // otherwise gradphi is almost zero and we keep a zero velocity
    //TODO:
    // get velocity at integration point
//    LINALG::Matrix<my::nsd_,1> convelint(true);
    convelint.Multiply(my::evelnp_,my::funct_);

    // convective part in convective form: u_x*N,x+ u_y*N,y
    LINALG::Matrix<my::nen_,1> conv(true);
    conv.MultiplyTN(my::derxy_,convelint);

    // convective term using current scalar value
    double conv_phi(0.0);
    conv_phi = convelint.Dot(gradphi);

    // get history data (or acceleration)
    double hist(0.0);
    // TODO:
    // use history vector of global level
    //hist = my::funct_.Dot(my::ehist_[0]);
    // recompute history
    // as long as the correction is not applied as a corrector step both
    // ways are equivalent
    // if we use a correction step than we loose the link used in the hist calculation
#if 1
    LINALG::Matrix<my::nsd_,1> gradphin(true);
    gradphin.Multiply(my::derxy_,my::ephin_[0]);
    // get norm
    const double gradphin_norm = gradphin.Norm2();
    double phin = 0.0;
    phin = my::funct_.Dot(my::ephin_[0]);
    double oldsign = 0.0;
    SignFunction(oldsign,charelelength,phizero,gradphizero,phin,gradphin);
    LINALG::Matrix<my::nsd_,1> convelintold(true);
    if (gradphin_norm>1e-8)
         convelintold.Update(oldsign/gradphin_norm,gradphin);
    hist = phin - my::scatraparatimint_->Dt() * (1.0 - my::scatraparatimint_->TimeFac()/my::scatraparatimint_->Dt()) * (convelintold.Dot(gradphin)-oldsign);
#endif

    //--------------------------------------------------------------------
    // calculation of stabilization parameter at integration point
    //--------------------------------------------------------------------

    // subgrid-scale velocity vector in gausspoint
    LINALG::Matrix<my::nsd_,1> sgvelint(true);

    if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
    {
      if (my::scatrapara_->TauGP())
        // calculation of stabilization parameter at integration point
        my::CalcTau(tau,0.0,0.0,1.0,convelint,vol,0);
    }

    // residual of convection-diffusion-reaction eq
    double scatrares(0.0);
    // residual-based subgrid-scale scalar (just a dummy here)
    double sgphi(0.0);

    // compute residual of scalar transport equation and
    // subgrid-scale part of scalar
    my::CalcResidualAndSubgrScalar(0,scatrares,sgphi,1.0,1.0,hist,conv_phi,0.0,0.0,signphi,tau);

    //----------------------------------------------------------------
    // evaluation of matrix and rhs
    //----------------------------------------------------------------

    // stabilization parameter and integration factors
    const double taufac     = tau*fac;
    const double timefacfac = my::scatraparatimint_->TimeFac()*fac;
    const double timetaufac = my::scatraparatimint_->TimeFac()*taufac;

    //----------------------------------------------------------------
    // 1) element matrix: instationary terms
    //----------------------------------------------------------------

    my::CalcMatMass(emat,0,fac,1.0,1.0);

    // diffusive part used in stabilization terms (dummy here)
    LINALG::Matrix<my::nen_,1> diff(true);
    // subgrid-scale velocity (dummy)
    LINALG::Matrix<my::nen_,1> sgconv(true);
    if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      my::CalcMatMassStab(emat,0,taufac,1.0,1.0,conv,sgconv,diff);

    //----------------------------------------------------------------
    // 2) element matrix: convective term in convective form
    //----------------------------------------------------------------

    my::CalcMatConv(emat,0,timefacfac,1.0,conv,sgconv);

    // convective stabilization of convective term (in convective form)
    // transient stabilization of convective term (in convective form)
    if(my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      my::CalcMatTransConvDiffStab(emat,0,timetaufac,1.0,conv,sgconv,diff);

    //----------------------------------------------------------------
    // 3) element right hand side
    //----------------------------------------------------------------

    double rhsint    = signphi;
    double rhsfac    = my::scatraparatimint_->TimeFacRhs() * fac;
    double rhstaufac = my::scatraparatimint_->TimeFacRhsTau() * taufac;

    // linearization of transient term
    my::CalcRHSLinMass(erhs,0,rhsfac,fac,1.0,1.0,hist);

    // the order of the following three functions is important
    // and must not be changed
    my::ComputeRhsInt(rhsint,1.0,1.0,hist);
    double rea_phi(0.0); // dummy
    my::RecomputeScatraResForRhs(scatrares,0,convelint,gradphi,diff,1.0,1.0,conv_phi,rea_phi,rhsint);
    // note: the third function is not required here, since we neither have a subgrid velocity
    //       nor a conservative form

    // standard Galerkin transient, old part of rhs and bodyforce term
    my::CalcRHSHistAndSource(erhs,0,fac,rhsint);

    // linearization of convective term
    my::CalcRHSConv(erhs,0,rhsfac,conv_phi);

    // linearization of stabilization terms
    if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      my::CalcRHSTransConvDiffStab(erhs,0,rhstaufac,1.0,scatrares,conv,sgconv,diff);

  } // end: loop all Gauss points

  return;
}
#endif


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)             rasthofer 12/13 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Sysmat(
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs ///< element rhs to calculate
  )
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = my::EvalShapeFuncAndDerivsAtEleCenter();

  //----------------------------------------------------------------------
  // calculation of characteristic element length
  //----------------------------------------------------------------------

  // get gradient of initial phi at element center
  LINALG::Matrix<my::nsd_,1> gradphizero(true);
  gradphizero.Multiply(my::derxy_,ephizero_[0]);

  // get characteristic element length
  const double charelelength = CalcCharEleLengthReinit(vol,gradphizero);

  //----------------------------------------------------------------------
  // prepare diffusion manager
  //----------------------------------------------------------------------

  // set diffusion coefficient of scalar 0 to 0.0
  if (not my::scatrapara_->MatGP())
    my::diffmanager_->SetIsotropicDiff(0.0,0);

  //----------------------------------------------------------------------
  // calculation of stabilization parameter at element center
  //----------------------------------------------------------------------

  // the stabilization parameter
  double tau = 0.0;

  if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
  {
    if (not my::scatrapara_->TauGP())
    {
      // get velocity at element center
      LINALG::Matrix<my::nsd_,1> convelint(true);

      // switch type for velocity field
      if (not dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->UseProjectedVel())
      {
#ifndef USE_PHIN_FOR_VEL
        // gradient of current scalar value at element center
        LINALG::Matrix<my::nsd_,1> gradphinp(true);
        gradphinp.Multiply(my::derxy_,my::ephinp_[0]);
        // get norm
        const double gradphinp_norm = gradphinp.Norm2();
        // get sign function
        double signphi = 0.0;
        // initial phi at element center
        double phizero = 0.0;
        phizero = my::funct_.Dot(ephizero_[0]);
        // current phi at element center
        double phinp = 0.0;
        phinp = my::funct_.Dot(my::ephinp_[0]);
        SignFunction(signphi,charelelength,phizero,gradphizero,phinp,gradphinp);

        if (gradphinp_norm>1e-8)
          convelint.Update(signphi/gradphinp_norm,gradphinp);
        // otherwise gradphi is almost zero and we keep a zero velocity
#else
        // gradient of scalar value at t_n at element center
        LINALG::Matrix<my::nsd_,1> gradphin(true);
        gradphin.Multiply(my::derxy_,my::ephin_[0]);
        // get norm
        const double gradphin_norm = gradphin.Norm2();
        // get sign function
        double signphi = 0.0;
        // initial phi at element center
        double phizero = 0.0;
        phizero = my::funct_.Dot(ephizero_[0]);
        // phi at element center
        double phin = 0.0;
        phin = my::funct_.Dot(my::ephin_[0]);
        SignFunction(signphi,charelelength,phizero,gradphizero,phin,gradphin);

        if (gradphin_norm>1e-8)
          convelint.Update(signphi/gradphin_norm,gradphin);
        // otherwise gradphi is almost zero and we keep a zero velocity
#endif
      }
      else
      {
        convelint.Multiply(my::econvelnp_,my::funct_);
      }

      // calculation of stabilization parameter at element center
      // here, second argument is isoptropic diffusion, which is zero!
      my::CalcTau(tau,0.0,0.0,1.0,convelint,vol);
    }
  }

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

    // set diffusion coefficient of scalar 0 to 0.0
    if (my::scatrapara_->MatGP())
      my::diffmanager_->SetIsotropicDiff(0.0,0);

    //--------------------------------------------------------------------
    // set all Gauss point quantities
    //--------------------------------------------------------------------

    // gradient of current scalar value at integration point
    LINALG::Matrix<my::nsd_,1> gradphinp(true);
    gradphinp.Multiply(my::derxy_,my::ephinp_[0]);
    // scalar at integration point at time step n+1
    const double phinp = my::funct_.Dot(my::ephinp_[0]);

    // initial phi at integration point
    double phizero = 0.0;
    phizero = my::funct_.Dot(ephizero_[0]);
    // and corresponding gradient
    gradphizero.Clear();
    gradphizero.Multiply(my::derxy_,ephizero_[0]);

    // scalar at integration point at time step n
    const double phin = my::funct_.Dot(my::ephin_[0]);

#if 0 // OLD but running
    // gradient of current scalar value Gauss point
    LINALG::Matrix<my::nsd_,1> gradphi(true);
    gradphi.Multiply(my::derxy_,my::ephinp_[0]);
    // get norm
    const double gradphi_norm = gradphi.Norm2();
    // get sign function
    double signphi = 0.0;
    // initial phi at Gauss point
    double phizero = 0.0;
    phizero = my::funct_.Dot(ephizero_[0]);
    gradphizero.Clear();
    gradphizero.Multiply(my::derxy_,ephizero_[0]);
//    std::cout << phizero << std::endl;
    // current phi at Gauss point
    double phi = 0.0;
    phi = my::funct_.Dot(my::ephinp_[0]);
    SignFunction(signphi,charelelength,phizero,gradphizero,phi,gradphi);
    if (my::eid_==10) std::cout << signphi << std::endl;

    // get velocity at element center
    LINALG::Matrix<my::nsd_,1> convelint(true);
    //if (gradphi_norm>1e-8)
      //convelint.Update(signphi/gradphi_norm,gradphi);
    LINALG::Matrix<my::nsd_,1> mygradphin(true);
    mygradphin.Multiply(my::derxy_,my::ephin_[0]);
    double mynorm = mygradphin.Norm2();
    if (mynorm>1e-8)
    convelint.Update(signphi/mynorm,mygradphin);
    // otherwise gradphi is almost zero and we keep a zero velocity
    //TODO:
    // get velocity at integration point
//    LINALG::Matrix<my::nsd_,1> convelint(true);
//    convelint.Multiply(my::econvelnp_,my::funct_);
#endif

    // get velocity at element center
    LINALG::Matrix<my::nsd_,1> convelint(true);

    // get sign function
    double signphi = 0.0;
#ifndef USE_PHIN_FOR_VEL
    SignFunction(signphi,charelelength,phizero,gradphizero,phinp,gradphinp);
#else
    // gradient of scalar value at t_n at integration point
    LINALG::Matrix<my::nsd_,1> gradphin(true);
    gradphin.Multiply(my::derxy_,my::ephin_[0]);

    SignFunction(signphi,charelelength,phizero,gradphizero,phin,gradphin);
#endif

    // switch type for velocity field
    if (not dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->UseProjectedVel())
    {
#ifndef USE_PHIN_FOR_VEL
      // get norm
      const double gradphinp_norm = gradphinp.Norm2();

      if (gradphinp_norm>1e-8)
        convelint.Update(signphi/gradphinp_norm,gradphinp);
      // otherwise gradphi is almost zero and we keep a zero velocity
#else
      // get norm
      const double gradphin_norm = gradphin.Norm2();

      if (gradphin_norm>1e-8)
        convelint.Update(signphi/gradphin_norm,gradphin);
      // otherwise gradphi is almost zero and we keep a zero velocity
#endif
    }
    else
    {
      convelint.Multiply(my::econvelnp_,my::funct_);
    }

    // convective part in convective form: u_x*N,x+ u_y*N,y
    LINALG::Matrix<my::nen_,1> conv(true);
    conv.MultiplyTN(my::derxy_,convelint);

    // convective term using current scalar value
    double conv_phi(0.0);
    conv_phi = convelint.Dot(gradphinp);

    // diffusive part used in stabilization terms
    double diff_phi(0.0);
    LINALG::Matrix<my::nen_,1> diff(true);
    // diffusive term using current scalar value for higher-order elements
    if (my::use2ndderiv_)
    {
      // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
      my::GetLaplacianStrongForm(diff);
      diff.Scale(0.0);
      diff_phi = diff.Dot(my::ephinp_[0]);
    }

    // get history data (or acceleration)
    double hist(0.0);
    // use history vector of global level
    hist = my::funct_.Dot(my::ehist_[0]);
    // recompute history
    // as long as the correction is not applied as a corrector step both
    // ways are equivalent
    // if we use a correction step than we loose the link used in the hist calculation
#if 0
#ifndef USE_PHIN_FOR_VEL
    LINALG::Matrix<my::nsd_,1> gradphin(true);
    gradphin.Multiply(my::derxy_,my::ephin_[0]);
    // get norm
    const double gradphin_norm = gradphin.Norm2();
    double phin = 0.0;
    phin = my::funct_.Dot(my::ephin_[0]);
#endif
    double oldsign = 0.0;
    SignFunction(oldsign,charelelength,phizero,gradphizero,phin,gradphin);
    LINALG::Matrix<my::nsd_,1> convelintold(true);
    const double gradphin_norm = gradphin.Norm2();
    if (gradphin_norm>1e-8)
         convelintold.Update(oldsign/gradphin_norm,gradphin);
    hist = phin - my::scatraparatimint_->Dt() * (1.0 - my::scatraparatimint_->TimeFac()/my::scatraparatimint_->Dt()) * (convelintold.Dot(gradphin)-oldsign);
#endif

    //--------------------------------------------------------------------
    // calculation of stabilization parameter at integration point
    //--------------------------------------------------------------------

    // subgrid-scale velocity vector in gausspoint
    //LINALG::Matrix<my::nsd_,1> sgvelint(true);

    if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
    {
      if (my::scatrapara_->TauGP())
        // calculation of stabilization parameter at integration point
        // here, second argument is isoptropic diffusion, which is zero!
        my::CalcTau(tau,0.0,0.0,1.0,convelint,vol);
    }

    //--------------------------------------------------------------------
    // calculation of artificial diffusion
    //--------------------------------------------------------------------

    if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->ArtDiff()
        != INPAR::SCATRA::artdiff_none)
    {
      // residual of reinit eq
      double scatrares = 0.0;
      // subgrid scalar (dummy)
      double sgphi = 0.0;
      my::CalcResidualAndSubgrScalar(0,scatrares,sgphi,1.0,1.0,phinp,hist,conv_phi,diff_phi,0.0,signphi,tau);

      // compute artificial diffusion
      // diffusion coefficient has been explicitly set to zero
      // additionally stored in subgrid diffusion coefficient
      my::CalcArtificialDiff(vol,0,my::diffmanager_,1.0,convelint,gradphinp,conv_phi,scatrares,tau);

      if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->ArtDiff()
          == INPAR::SCATRA::artdiff_crosswind)
        Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerLsReinit<my::nsd_> >(my::diffmanager_)->SetVelocityForCrossWindDiff(convelint);

#ifdef MODIFIED_EQ
      // recompute tau to get adaption to artificial diffusion
      if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      {
        if (my::scatrapara_->TauGP())
          // calculation of stabilization parameter at integration point
          // here, second argument is isoptropic diffusion, which is zero!
          my::CalcTau(tau,my::diffmanager_->GetIsotropicDiff(0),0.0,1.0,convelint,vol,0);
      }

      // recompute diff operator
      if (my::use2ndderiv_)
      {
        // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
        my::GetLaplacianStrongForm(diff);
        diff.Scale(my::diffmanager_->GetIsotropicDiff(0));
        diff_phi = diff.Dot(my::ephinp_[0]);
      }
#endif
    }

    // residual of convection-diffusion-reaction eq
    double scatrares(0.0);
    // residual-based subgrid-scale scalar (just a dummy here)
    double sgphi(0.0);
    // compute residual of scalar transport equation and
    // subgrid-scale part of scalar
    my::CalcResidualAndSubgrScalar(0,scatrares,sgphi,1.0,1.0,phinp,hist,conv_phi,diff_phi,0.0,signphi,tau);

    //----------------------------------------------------------------
    // evaluation of matrix and rhs
    //----------------------------------------------------------------

    // stabilization parameter and integration factors
    const double taufac     = tau*fac;
    const double timefacfac = my::scatraparatimint_->TimeFac()*fac;
    const double dtfac      = my::scatraparatimint_->Dt()*fac;
    const double timetaufac = my::scatraparatimint_->TimeFac()*taufac;

    //----------------------------------------------------------------
    // 1) element matrix: instationary terms
    //----------------------------------------------------------------

    my::CalcMatMass(emat,0,fac,1.0);

    // subgrid-scale velocity (dummy)
    LINALG::Matrix<my::nen_,1> sgconv(true);
    if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->LinForm() == INPAR::SCATRA::newton)
    {
      if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
        my::CalcMatMassStab(emat,0,taufac,1.0,1.0,conv,sgconv,diff);
    }

    //----------------------------------------------------------------
    // 2) element matrix: convective term in convective form
    //----------------------------------------------------------------

    if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->LinForm() == INPAR::SCATRA::newton)
      my::CalcMatConv(emat,0,timefacfac,1.0,conv,sgconv);

    // convective stabilization of convective term (in convective form)
    // transient stabilization of convective term (in convective form)
    if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->LinForm() == INPAR::SCATRA::newton)
    {
      if(my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
        my::CalcMatTransConvDiffStab(emat,0,timetaufac,1.0,conv,sgconv,diff);
    }


    // calculation of diffusive element matrix
    if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->ArtDiff()
        != INPAR::SCATRA::artdiff_none)
#ifndef MODIFIED_EQ
      CalcMatDiff(emat,0,dtfac,my::diffmanager_); //implicit treatment
#else
      CalcMatDiff(emat,0,timefacfac,my::diffmanager_);
#endif

    //----------------------------------------------------------------
    // 3) element right hand side
    //----------------------------------------------------------------

    double rhsint    = signphi;
    double rhsfac    = my::scatraparatimint_->TimeFacRhs() * fac;
    double rhstaufac = my::scatraparatimint_->TimeFacRhsTau() * taufac;

    // linearization of transient term
    my::CalcRHSLinMass(erhs,0,rhsfac,fac,1.0,1.0,phinp,hist);

    // the order of the following three functions is important
    // and must not be changed
    my::ComputeRhsInt(rhsint,1.0,1.0,hist);
    double rea_phi(0.0); // dummy
    my::RecomputeScatraResForRhs(scatrares,0,convelint,gradphinp,diff,1.0,1.0,conv_phi,rea_phi,phin,my::reamanager_,rhsint);
    // note: the third function is not required here, since we neither have a subgrid velocity
    //       nor a conservative form

    // standard Galerkin transient, old part of rhs and bodyforce term
    my::CalcRHSHistAndSource(erhs,0,fac,rhsint);

    // linearization of convective term
    my::CalcRHSConv(erhs,0,rhsfac,conv_phi);

    // linearization of diffusive term
    if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->ArtDiff()
        != INPAR::SCATRA::artdiff_none)
#ifndef MODIFIED_EQ
      CalcRHSDiff(erhs,0,dtfac,my::diffmanager_,gradphinp);  //implicit treatment
#else
      CalcRHSDiff(erhs,0,rhsfac,my::diffmanager_,gradphinp);
#endif

    // linearization of stabilization terms
    if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      my::CalcRHSTransConvDiffStab(erhs,0,rhstaufac,1.0,scatrares,conv,sgconv,diff);

  } // end: loop all Gauss points

  return;
}


/*----------------------------------------------------------------------*
 |  sign function                                       rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::SignFunction(
  double&                           sign_phi,
  const double                      charelelength,
  const double                      phizero,
  const LINALG::Matrix<my::nsd_,1>& gradphizero,
  const double                      phi,
  const LINALG::Matrix<my::nsd_,1>& gradphi
  )
{
  // compute interface thickness
  const double epsilon = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->InterfaceThicknessFac() * charelelength;

  if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->SignType()
      == INPAR::SCATRA::signtype_nonsmoothed)
  {
    if      (phizero < 0.0) sign_phi = -1.0;
    else if (phizero > 0.0) sign_phi = +1.0;
    else                    sign_phi =  0.0;
  }
  else if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->SignType()
           == INPAR::SCATRA::signtype_SussmanSmerekaOsher1994)
  {
    sign_phi = phizero/sqrt(phizero*phizero + epsilon*epsilon);
  }
  else if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->SignType()
           == INPAR::SCATRA::signtype_PengEtAl1999)
  {
    const double grad_norm_phi = gradphi.Norm2();
    sign_phi = phi/sqrt(phi*phi + epsilon*epsilon * grad_norm_phi*grad_norm_phi);
  }
  else if (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->SignType()
           == INPAR::SCATRA::signtype_SussmanFatemi1999)
  {
    if(fabs(epsilon) < 1e-15) dserror("divide by zero in evaluate for smoothed sign function");

    if (phizero < -epsilon)       sign_phi = -1.0;
    else if (phizero > epsilon)   sign_phi = +1.0;
    else                          sign_phi = phizero/epsilon + 1.0/PI*sin(PI*phizero/epsilon);
  }
  else dserror("unknown type of sign function!");
  return;
}


/*----------------------------------------------------------------------*
 | derivative of sign function                          rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::DerivSignFunction(
  double&                           deriv_sign,
  const double                      charelelength,
  const double                      phizero
  )
{
  // compute interface thickness
  const double epsilon = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->InterfaceThicknessFac() * charelelength;

  if (phizero < -epsilon)       deriv_sign = 0.0;
  else if (phizero > epsilon)   deriv_sign = 0.0;
  else                          deriv_sign = 1.0/(2.0*epsilon)*(1.0 + cos(PI*phizero/epsilon));

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of characteristic element length        rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::CalcCharEleLengthReinit(
  const double                      vol,
  const LINALG::Matrix<my::nsd_,1>& gradphizero
  )
{
  // define and initialize length
  double h = 0.0;

  //---------------------------------------------------------------------
  // select from various definitions for characteristic element length
  //---------------------------------------------------------------------
  switch (dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterLsReinit*>(my::scatrapara_)->CharEleLengthReinit())
  {
    // a) streamlength due to Kees et al. (2011)
    case INPAR::SCATRA::streamlength_reinit:
    {
      // get norm of gradient of phi
      double gradphi_norm = gradphizero.Norm2();
      LINALG::Matrix<my::nsd_,1> gradphi_scaled(true);
      if (gradphi_norm>=1e-8) gradphi_scaled.Update(1.0/gradphi_norm,gradphizero);
      else
      {
        //TODO: clearify this
        dserror("gradphi_norm=0: cannot compute characteristic element length");
        gradphi_scaled.Update(1.0,gradphizero);
        //gradphi_scaled.Clear();
        //gradphi_scaled(0,0) = 1.0;
      }

      // computation of covariant metric tensor
      double G;
      double Gnormgradphi(0.0);
      for (int nn=0;nn<my::nsd_;++nn)
      {
        for (int rr=0;rr<my::nsd_;++rr)
        {
          G = my::xij_(nn,0)*my::xij_(rr,0);
          for(int tt=1;tt<my::nsd_;tt++)
          {
            G += my::xij_(nn,tt)*my::xij_(rr,tt);
          }
          Gnormgradphi+=gradphi_scaled(nn,0)*G*gradphi_scaled(rr,0);
        }
      }

      h = 1.0/std::sqrt(Gnormgradphi);

    }
    break;

    // b) cubic/square root of element volume/area or element length (3-/2-/1-D)
    case INPAR::SCATRA::root_of_volume_reinit:
    {
      // cast dimension to a double varibale -> pow()
      const double dim = double (my::nsd_);
      h = std::pow(vol,1.0/dim);
    }
    break;

    default: dserror("unknown characteristic element length\n");
    break;
  }

  return h;
}


/*------------------------------------------------------------------- *
 | calculation of diffusive element matrix            rasthofer 12/13 |
 | here we consider both isotropic and crosswind artificial diffusion |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::CalcMatDiff(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  Teuchos::RCP<ScaTraEleDiffManager>  diffmanager
  )
{
  // flag for anisotropic diffusion
  const bool crosswind = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerLsReinit<my::nsd_> >(diffmanager)->HaveCrossWindDiff();

  // set scalar diffusion: isotropic or factor for anisotropic tensor in case of
  // crosswind diffusion
  double diffus = diffmanager->GetIsotropicDiff(k);

  // diffusive factor
  const double fac_diffus = timefacfac * diffus;

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;
      double laplawf(0.0);

      // in case of anisotropic or crosswind diffusion, multiply 'derxy_' with diffusion tensor
      // inside 'GetLaplacianWeakForm'
      if (crosswind)
        my::GetLaplacianWeakForm(laplawf,
                                 Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerLsReinit<my::nsd_> >(diffmanager)->GetCrosswindTensor(),
                                 ui,vi);
      else
        my::GetLaplacianWeakForm(laplawf,ui,vi);

      emat(fvi,fui) += fac_diffus*laplawf;
    }
  }
  return;
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side rasthofer 12/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::CalcRHSDiff(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac,
  Teuchos::RCP<ScaTraEleDiffManager>  diffmanager,
  const LINALG::Matrix<my::nsd_,1>& gradphi
  )
{
  // flag for anisotropic diffusion
  const bool crosswind = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerLsReinit<my::nsd_> >(diffmanager)->HaveCrossWindDiff();

  // set scalar diffusion: isotropic or factor for anisotropic tensor in case of
  // crosswind diffusion
  double vrhs = rhsfac*diffmanager->GetIsotropicDiff(k);

  LINALG::Matrix<my::nsd_,1> gradphirhs(true);
  if (crosswind)
    // in case of anisotropic or crosswind diffusion, multiply 'gradphi' with diffusion tensor
    gradphirhs.Multiply(Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerLsReinit<my::nsd_> >(diffmanager)->GetCrosswindTensor(),gradphi);
  else
    gradphirhs.Update(1.0,gradphi,0.0);

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    double laplawf(0.0);

    this->GetLaplacianWeakFormRHS(laplawf,gradphirhs,vi);

    erhs[fvi] -= vrhs*laplawf;
  }

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::nurbs27>;



