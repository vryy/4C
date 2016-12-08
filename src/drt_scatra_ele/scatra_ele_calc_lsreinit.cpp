/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_lsreinit.cpp

\brief evaluation of scatra elements for reinitialization equation

\level 2

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_lsreinit.H"

#include "scatra_ele.H"
#include "scatra_ele_parameter_lsreinit.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/integrationcell_coordtrafo.H"


#define USE_PHIN_FOR_VEL
//#define MODIFIED_EQ

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype> * DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcLsReinit* delete_me )
{
  static std::map<std::string,ScaTraEleCalcLsReinit<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcLsReinit<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcLsReinit<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::ScaTraEleCalcLsReinit(
    const int numdofpernode,
    const int numscal,
    const std::string& disname
    ) :
    DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
    ephizero_(my::numscal_),  // size of vector
    lsreinitparams_(DRT::ELEMENTS::ScaTraEleParameterLsReinit::Instance(disname)) // parameter class
{
  // set appropriate diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerLsReinit<my::nsd_>(my::numscal_));
  // set appropriate internal variable manager
  my::scatravarmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerLsReinit<my::nsd_,my::nen_>(my::numscal_));

  // safety checks
  if(my::scatrapara_->RBSubGrVel())
    dserror("CalcSubgrVelocityLevelSet not available anymore");
  if(lsreinitparams_->ArtDiff() != INPAR::SCATRA::artdiff_none)
  {
    if(not my::scatrapara_->MatGP() or not my::scatrapara_->TauGP())
     dserror("Evaluation of material and stabilization parameters need to be done at the integration points for reinitialization due to artificial diff!");
  }

  return;
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::Evaluate(
  DRT::Element*                 ele,
  Teuchos::ParameterList&       params,
  DRT::Discretization&          discretization,
  DRT::Element::LocationArray&  la,
  Epetra_SerialDenseMatrix&     elemat1_epetra,
  Epetra_SerialDenseMatrix&     elemat2_epetra,
  Epetra_SerialDenseVector&     elevec1_epetra,
  Epetra_SerialDenseVector&     elevec2_epetra,
  Epetra_SerialDenseVector&     elevec3_epetra
  )
{
  // setup
  if(SetupCalc(ele,discretization) == -1)
    return 0;

  //(for now) only first dof set considered
  const std::vector<int> &    lm = la[0].lm_;

  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==Teuchos::null)
    dserror("Cannot get state vector 'phinp'");
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phinp,my::ephinp_,lm);

  // distinguish reinitalization
  switch(lsreinitparams_->ReinitType())
  {
    case INPAR::SCATRA::reinitaction_sussman:
    {
      // extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
      Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
      Teuchos::RCP<const Epetra_Vector> phizero = discretization.GetState("phizero");
      if (hist==Teuchos::null || phinp==Teuchos::null || phin==Teuchos::null || phizero==Teuchos::null)
        dserror("Cannot get state vector 'hist' and/or 'phin' and/or 'phizero'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phin,my::ephin_,lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phizero,ephizero_,lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*hist,my::ehist_,lm);

      if(lsreinitparams_->UseProjectedVel())
      {
        // get velocity at nodes (pre-computed via L2 projection)
        const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("reinitialization velocity field");
        DRT::UTILS::ExtractMyNodeBasedValues(ele,my::econvelnp_,velocity,my::nsd_);
      }

      // calculate element coefficient matrix and rhs
      SysmatHyperbolic(
        elemat1_epetra,
        elevec1_epetra);
      break;
    }
    case INPAR::SCATRA::reinitaction_ellipticeq:
    {
      // extract boundary integration cells for elliptic reinitialization
      // define empty list
      GEO::BoundaryIntCells boundaryIntCells = GEO::BoundaryIntCells();

      Teuchos::RCP<std::map<int,GEO::BoundaryIntCells > > allcells =
              params.get<Teuchos::RCP<std::map<int,GEO::BoundaryIntCells > > >("boundary cells", Teuchos::null);

      std::map<int,GEO::BoundaryIntCells>::const_iterator tmp = allcells->find(my::eid_);
      if (tmp != allcells->end())
        boundaryIntCells = tmp->second;

      if(lsreinitparams_->Project())
      {
        const Teuchos::RCP<Epetra_MultiVector> gradphi = params.get< Teuchos::RCP<Epetra_MultiVector> >("gradphi");
        DRT::UTILS::ExtractMyNodeBasedValues(ele,my::econvelnp_,gradphi,my::nsd_);
      }

      // calculate element coefficient matrix and rhs
      SysmatElliptic(
          elemat1_epetra,
          elevec1_epetra,
          boundaryIntCells);
      break;
    }
    default:
      dserror("Unknown reinitialization equation");
      break;
  }

  return 0;
}


#if 0
/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)             rasthofer 12/13 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::SysmatHyperbolic(
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
  // integration points and weights
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

    // compute residual of scalar transport equation and
    // subgrid-scale part of scalar
    my::CalcStrongResidual(0,scatrares,1.0,1.0,hist,conv_phi,0.0,signphi,tau);

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
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::SysmatHyperbolic(
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
      if (not lsreinitparams_->UseProjectedVel())
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
  // integration points and weights
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

    // also store values in variable manager
    VarManager()->SetPhinp(0,phinp);
    VarManager()->SetPhin(0,phin);
    VarManager()->SetGradPhi(0,gradphinp);

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
    if (not lsreinitparams_->UseProjectedVel())
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

    //set changed values in variable manager
    VarManager()->SetConv(conv);
    VarManager()->SetConVel(convelint);
    VarManager()->SetConvPhi(0,conv_phi);

    LINALG::Matrix<my::nen_,1> diff(true);
    // diffusive term using current scalar value for higher-order elements
    if (my::use2ndderiv_)
    {
      // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
      my::GetLaplacianStrongForm(diff);
      diff.Scale(0.0);
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
    //set changed values in variable manager
    VarManager()->SetHist(0,hist);

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

    if (lsreinitparams_->ArtDiff()
        != INPAR::SCATRA::artdiff_none)
    {
      // residual of reinit eq
      double scatrares = 0.0;

      my::CalcStrongResidual(0,scatrares,1.0,1.0,0.0,signphi,tau);

      // compute artificial diffusion
      // diffusion coefficient has been explicitly set to zero
      // additionally stored in subgrid diffusion coefficient
      my::CalcArtificialDiff(vol,0,1.0,convelint,gradphinp,conv_phi,scatrares,tau);

      if (lsreinitparams_->ArtDiff()
          == INPAR::SCATRA::artdiff_crosswind)
        DiffManager()->SetVelocityForCrossWindDiff(convelint);

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
    // compute residual of scalar transport equation and
    // subgrid-scale part of scalar
    my::CalcStrongResidual(0,scatrares,1.0,1.0,0.0,signphi,tau);

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
    if (lsreinitparams_->LinForm() == INPAR::SCATRA::newton)
    {
      if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
        my::CalcMatMassStab(emat,0,taufac,1.0,1.0,sgconv,diff);
    }

    //----------------------------------------------------------------
    // 2) element matrix: convective term in convective form
    //----------------------------------------------------------------

    if (lsreinitparams_->LinForm() == INPAR::SCATRA::newton)
      my::CalcMatConv(emat,0,timefacfac,1.0,sgconv);

    // convective stabilization of convective term (in convective form)
    // transient stabilization of convective term (in convective form)
    if (lsreinitparams_->LinForm() == INPAR::SCATRA::newton)
    {
      if(my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
        my::CalcMatTransConvDiffStab(emat,0,timetaufac,1.0,sgconv,diff);
    }


    // calculation of diffusive element matrix
    if (lsreinitparams_->ArtDiff()
        != INPAR::SCATRA::artdiff_none)
#ifndef MODIFIED_EQ
      CalcMatDiff(emat,0,dtfac); //implicit treatment
#else
      CalcMatDiff(emat,0,timefacfac);
#endif

    //----------------------------------------------------------------
    // 3) element right hand side
    //----------------------------------------------------------------

    double rhsint    = signphi;
    double rhsfac    = my::scatraparatimint_->TimeFacRhs() * fac;
    double rhstaufac = my::scatraparatimint_->TimeFacRhsTau() * taufac;

    // linearization of transient term
    my::CalcRHSLinMass(erhs,0,rhsfac,fac,1.0,1.0);

    // the order of the following three functions is important
    // and must not be changed
    my::ComputeRhsInt(rhsint,1.0,1.0,hist);
    double rea_phi(0.0); // dummy
    my::RecomputeScatraResForRhs(scatrares,0,diff,1.0,1.0,rea_phi,rhsint);
    // note: the third function is not required here, since we neither have a subgrid velocity
    //       nor a conservative form

    // standard Galerkin transient, old part of rhs and bodyforce term
    my::CalcRHSHistAndSource(erhs,0,fac,rhsint);

    // linearization of convective term
    my::CalcRHSConv(erhs,0,rhsfac);

    // linearization of diffusive term
    if (lsreinitparams_->ArtDiff()
        != INPAR::SCATRA::artdiff_none)
#ifndef MODIFIED_EQ
      CalcRHSDiff(erhs,0,dtfac,gradphinp);  //implicit treatment
#else
      CalcRHSDiff(erhs,0,rhsfac,gradphinp);
#endif

    // linearization of stabilization terms
    if (my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
      my::CalcRHSTransConvDiffStab(erhs,0,rhstaufac,1.0,scatrares,sgconv,diff);

  } // end: loop all Gauss points

  return;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)             rasthofer 09/14 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::SysmatElliptic(
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs, ///< element rhs to calculate
  const GEO::BoundaryIntCells&          boucell ///< interface for penalty term
  )
{
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

    // gradient of current scalar value at integration point
    LINALG::Matrix<my::nsd_,1> gradphinp(true);
    if (lsreinitparams_->Project())
      gradphinp.Multiply(my::econvelnp_,my::funct_);
    else
      gradphinp.Multiply(my::derxy_,my::ephinp_[0]);

    // TODO: remove
//    if (std::abs(gradphinp(2,0))>1.0e-8)
//    {
//        std::cout << gradphinp << std::setprecision(8) << my::ephinp_[0] << std::endl;
//        dserror("ENDE");
//    }

    double normgradphi = gradphinp.Norm2();

    //----------------------------------------------------------------------
    // prepare diffusion manager
    //----------------------------------------------------------------------

    // calculate nonlinear diffusivity
    double diff = 0.0;
    switch (lsreinitparams_->DiffFct())
    {
      case INPAR::SCATRA::hyperbolic:
      {
        // the basic form: goes to -infinity, if norm(grad)=1
        // yields directly desired signed-distance field
        if (normgradphi>1.0e-8)
          diff = 1.0 - (1.0/normgradphi);
        else
          diff = 1.0 - (1.0/1.0e-8);

        break;
      }
      case INPAR::SCATRA::hyperbolic_smoothed_positive:
      {
        // version as suggested by Basting and Kuzmin 2013
        // returns to positive values for norm(grad)<0.5
        if (normgradphi>1.0)
          diff = 1.0 - (1.0/normgradphi);
        else
          diff = 2.0 * normgradphi * normgradphi - 3.0 * normgradphi + 1.0;

        break;
      }
      case INPAR::SCATRA::hyperbolic_clipped_05:
      {
          if (normgradphi>(2.0/3.0))
            diff = 1.0 - (1.0/normgradphi);
          else
            diff = -0.5;

          break;
      }
      case INPAR::SCATRA::hyperbolic_clipped_1:
      {
        if (normgradphi>0.5)
          diff = 1.0 - (1.0/normgradphi);
        else
          diff = -1.0;

        break;
      }
      default:
      {
        dserror("Unknown diffusivity function!");
        break;
      }
    }

    // set diffusivity of scalar 0 to 1.0 for lhs term
    my::diffmanager_->SetIsotropicDiff(1.0,0);

    //----------------------------------------------------------------
    // evaluation of matrix and rhs
    //----------------------------------------------------------------

    //----------------------------------------------------------------
    // 1) element matrix: diffusion matrix
    //----------------------------------------------------------------

    my::CalcMatDiff(emat,0,fac);

    //----------------------------------------------------------------
    // 2) element right hand side
    //----------------------------------------------------------------

    // set diffusivity for rhs term
    my::diffmanager_->SetIsotropicDiff((-diff+1.0),0);
    // set gradphi for rhs term
    my::scatravarmanager_->SetGradPhi(0,gradphinp);

    my::CalcRHSDiff(erhs,0,-fac);

  } // end: loop all Gauss points

  //----------------------------------------------------------------
  // 3) evaluation of penalty term at initial interface
  //----------------------------------------------------------------

  EvaluateInterfaceTerm(emat,erhs,boucell);

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
  const double epsilon = lsreinitparams_->InterfaceThicknessFac() * charelelength;

  if (lsreinitparams_->SignType()
      == INPAR::SCATRA::signtype_nonsmoothed)
  {
    if      (phizero < 0.0) sign_phi = -1.0;
    else if (phizero > 0.0) sign_phi = +1.0;
    else                    sign_phi =  0.0;
  }
  else if (lsreinitparams_->SignType()
           == INPAR::SCATRA::signtype_SussmanSmerekaOsher1994)
  {
    sign_phi = phizero/sqrt(phizero*phizero + epsilon*epsilon);
  }
  else if (lsreinitparams_->SignType()
           == INPAR::SCATRA::signtype_PengEtAl1999)
  {
    const double grad_norm_phi = gradphi.Norm2();
    sign_phi = phi/sqrt(phi*phi + epsilon*epsilon * grad_norm_phi*grad_norm_phi);
  }
  else if (lsreinitparams_->SignType()
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
  const double epsilon = lsreinitparams_->InterfaceThicknessFac() * charelelength;

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
  switch (lsreinitparams_->CharEleLengthReinit())
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
  const double                  timefacfac
  )
{
  // flag for anisotropic diffusion
  const bool crosswind = DiffManager()->HaveCrossWindDiff();

  // set scalar diffusion: isotropic or factor for anisotropic tensor in case of
  // crosswind diffusion
  double diffus = DiffManager()->GetIsotropicDiff(k);

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
                                 DiffManager()->GetCrosswindTensor(),
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
  const LINALG::Matrix<my::nsd_,1>& gradphi
  )
{
  // flag for anisotropic diffusion
  const bool crosswind = DiffManager()->HaveCrossWindDiff();

  // set scalar diffusion: isotropic or factor for anisotropic tensor in case of
  // crosswind diffusion
  double vrhs = rhsfac*DiffManager()->GetIsotropicDiff(k);

  LINALG::Matrix<my::nsd_,1> gradphirhs(true);
  if (crosswind)
    // in case of anisotropic or crosswind diffusion, multiply 'gradphi' with diffusion tensor
    gradphirhs.Multiply(DiffManager()->GetCrosswindTensor(),gradphi);
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


/*-------------------------------------------------------------------- *
 | calculation of interface penalty term for elliptic reinitialization |
 |                                                     rasthofer 09/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::EvaluateInterfaceTerm(
  Epetra_SerialDenseMatrix&             emat,     //!< element matrix to calculate
  Epetra_SerialDenseVector&             erhs,     //!< element vector to calculate
  const GEO::BoundaryIntCells&          boucell   //!< interface for penalty term
)
{
  //------------------------------------------------------------------------------------------------
  // loop over boundary integration cells
  //------------------------------------------------------------------------------------------------
  for (GEO::BoundaryIntCells::const_iterator cell = boucell.begin(); cell != boucell.end(); ++cell)
  {
    // get shape of boundary cell
    DRT::Element::DiscretizationType celldistype = cell->Shape();

    GEO::BoundaryIntCell actcell = *cell;

    switch (celldistype)
    {
      case DRT::Element::tri3:
      {
        CalcPenaltyTerm<DRT::Element::tri3>(emat,erhs,actcell);
        break;
      }
      case DRT::Element::quad4:
      {
        CalcPenaltyTerm<DRT::Element::quad4>(emat,erhs,actcell);
        break;
      }
      default:
        dserror("cell distype not implemented yet");
        break;
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 | calculation of interface penalty term for elliptic reinitialization |
 | gauss loop                                          rasthofer 09/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType celldistype>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype>::CalcPenaltyTerm(
  Epetra_SerialDenseMatrix&             emat,     //!< element matrix to calculate
  Epetra_SerialDenseVector&             erhs,     //!< element vector to calculate
  const GEO::BoundaryIntCell            cell      //!< interface cell
)
{
  // get number of vertices of cell
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement;
  const size_t nsd = 3;
  if (my::nsd_!=3) dserror("Extend for other dimensions");
  const size_t nsd_cell = 2; //my::nsd_-1;
    // get coordinates of vertices of boundary integration cell in element coordinates \xi^domain
    LINALG::SerialDenseMatrix cellXiDomaintmp = cell.CellNodalPosXiDomain();
    //cellXiDomaintmp.Print(std::cout);
    // transform to fixed size format
    const LINALG::Matrix<nsd,numvertices> cellXiDomain(cellXiDomaintmp);

    //----------------------------------------------------------------------------------------------
    // integration loop over Gaussian points
    //----------------------------------------------------------------------------------------------
    // integrations points and weights
    DRT::UTILS::IntegrationPoints2D intpoints(SCATRA::CellTypeToOptGaussRule<celldistype>::rule);

    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // new transformation for boundary integrals
      // 1. define a coupled transformation x_3D(xi_3D(eta_2D)): transformation from 2D->3D
      // 2. compute the corresponding Jacobian J_eta2D->x_3D and
      // 3. the corresponding surface integral factor sqrt(det(J_eta2D->x_3D^T * J_eta2D->x_3D))
      // 4. approximate integral with Gauss rule in eta coordinate system
      // 5. evaluate the transformed integrand f(x(xi(eta)))

      const LINALG::Matrix<nsd_cell,1> gpinEta2D(intpoints.qxg[iquad]);

      // Jacobian for coupled transformation
      // get derivatives dxi_3D/deta_2D
      static LINALG::Matrix<nsd_cell,numvertices> deriv_eta2D;
      DRT::UTILS::shape_function_2D_deriv1(deriv_eta2D,gpinEta2D(0,0),gpinEta2D(1,0),celldistype);

      // calculate dxi3Ddeta2D
      static LINALG::Matrix<nsd,nsd_cell> dXi3Ddeta2D;
      dXi3Ddeta2D.Clear();

      for (int i = 0; i < (int)nsd; i++)   // dimensions
        for (int j = 0; j < (int)nsd_cell; j++) // derivatives
          for (int k = 0; k < (int)numvertices; k++)
          {
            dXi3Ddeta2D(i,j) += cellXiDomain(i,k)*deriv_eta2D(j,k);
          }

      // transform Gauss point to xi3D space (element parameter space)
      static LINALG::Matrix<nsd,1> gpinXi3D;
      gpinXi3D.Clear();

      // coordinates of this integration point in element coordinates \xi^domain
      GEO::mapEtaBToXiD(cell, gpinEta2D, gpinXi3D);

      static LINALG::Matrix<nsd,my::nen_> deriv_xi3D;
      DRT::UTILS::shape_function_3D_deriv1(deriv_xi3D,gpinXi3D(0,0), gpinXi3D(1,0), gpinXi3D(2,0), distype);

      // calculate dx3Ddxi3D
      static LINALG::Matrix<nsd,nsd> dX3DdXi3D;
      dX3DdXi3D.Clear();
      for (int i = 0; i < (int)nsd; i++)   // dimensions
        for (int j = 0; j < (int)nsd; j++) // derivatives
          for (int k = 0; k < my::nen_; k++)
            dX3DdXi3D(i,j) += my::xyze_(i,k)*deriv_xi3D(j,k);

      // get the coupled Jacobian dx3Ddeta2D
      static LINALG::Matrix<3,2> dx3Ddeta2D;
      dx3Ddeta2D.Clear();
      for (int i = 0; i < (int)nsd; i++)   // dimensions
        for (int j = 0; j < (int)nsd_cell; j++) // derivatives
          for (int k = 0; k < (int)nsd; k++)
            dx3Ddeta2D(i,j) += dX3DdXi3D(i,k) * dXi3Ddeta2D(k,j);

      // get deformation factor
      static LINALG::Matrix<nsd_cell,nsd_cell> Jac_tmp; // J^T*J
      Jac_tmp.Clear();
      Jac_tmp.MultiplyTN(dx3Ddeta2D,dx3Ddeta2D);

      if(Jac_tmp.Determinant() == 0.0) dserror("deformation factor for boundary integration is zero");
      const double deform_factor = sqrt(Jac_tmp.Determinant()); // sqrt(det(J^T*J))

      const double fac = intpoints.qwgt[iquad]*deform_factor;

      LINALG::Matrix<nsd_cell,1> posEtaBoundary;
      posEtaBoundary.Clear();
      for (int i= 0; i < (int)nsd_cell; i++)
        posEtaBoundary(i,0) = gpinEta2D(i,0);

      LINALG::Matrix<nsd,1> posXiDomain;
      posXiDomain.Clear();
      for (int i= 0; i< (int)nsd; i++)
        posXiDomain(i,0) = gpinXi3D(i,0);

      //--------------------------------------------------------------------------------------------
      // evaluate shape functions and their first derivatives at this Gaussian point
      //--------------------------------------------------------------------------------------------
      my::funct_.Clear();
      DRT::UTILS::shape_function_3D(my::funct_,posXiDomain(0),posXiDomain(1),posXiDomain(2),distype);

      //--------------------------------------------------------------------------------------------
      // evaluate mat and rhs
      //--------------------------------------------------------------------------------------------
      // caution density of original function is replaced by penalty parameter here
      my::CalcMatMass(emat,0,fac,lsreinitparams_->PenaltyPara());

      }// loop Gaussian points

  return;
}


// template classes

#include "scatra_ele_calc_fwd.hpp"
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri6>;
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
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::nurbs27>;
