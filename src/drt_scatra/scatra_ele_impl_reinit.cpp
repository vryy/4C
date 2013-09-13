/*----------------------------------------------------------------------*/
/*!
  \file scatra_ele_impl_reinit.cpp

  \brief Internal implementation of scalar transport elements

  <pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
  </pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_impl_reinit.H"
#include "scatra_ele_action.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/biofilm.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/yoghurt.H"
#include "../drt_mat/matlist.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ReInitImpl<distype> * DRT::ELEMENTS::ReInitImpl<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create
  )
{
  static ReInitImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ReInitImpl<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ReInitImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0,0,false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ReInitImpl<distype>::ReInitImpl(const int numdofpernode, const int numscal)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    is_ale_(false),           // bool set
    is_reactive_(false),      // bool set
    is_coupled_(false),       // bool set
    diffreastafac_(0.0),
    is_stationary_(false),    // bool set
    is_genalpha_(false),      // bool set
    is_incremental_(false),   // bool set
    is_conservative_(false),  // bool set
    betterconsistency_(false), // bool set
    // whichtau_ not initialized
    gradphi_(true),     // initialized to zero
    ephin_(numscal_),   // size of vector
    ephinp_(numscal_),  // size of vector
    ephiam_(numscal_),  // size of vector
    hist_(numscal_),    // size of vector
    ehist_(numscal_),   // size of vector
    ephi0_Reinit_Reference_(numscal_),// size of vector
    ephi0_penalty_(numscal_),         // size of vector
    fsphinp_(numscal_), // size of vector
    evelnp_(true),      // initialized to zero
    econvelnp_(true),   // initialized to zero
    efsvel_(true),      // initialized to zero
    eaccnp_(true),      // initialized to zero
    edispnp_(true),     // initialized to zero
    velint_(true),      // initialized to zero
    convelint_(true),   // initialized to zero
    sgvelint_(true),    // initialized to zero
    fsvelint_(true),    // initialized to zero
    mfsgvelint_(true),  // initialized to zero
    migvelint_(true),   // initialized to zero
    conv_(true),        // initialized to zero
    sgconv_(true),      // initialized to zero
    vdiv_(0.0),         // set double
    eprenp_(true),      // initialized to zero
    densn_(numscal_),        // size of vector
    densnp_(numscal_),       // size of vector
    densam_(numscal_),       // size of vector
    densgradfac_(numscal_),  // size of vector
    diffus_(numscal_),       // size of vector
    sgdiff_(numscal_),       // size of vector
    reacoeff_(numscal_),     // size of vector
    reacoeffderiv_(numscal_),// size of vector
    reacoeffderivmatrix_(numscal_,std::vector<double>(numscal_)),// size of matrix + initialized to zero
    shc_(0.0),      // set double
    visc_(0.0),     // set double
    diff_(true),    // initialized to zero
    migconv_(true), // initialized to zero
    migrea_(true),  // initialized to zero
    xsi_(true),     // initialized to zero
    xyze_(true),    // initialized to zero
    funct_(true),   // initialized to zero
    deriv_(true),   // initialized to zero
    deriv2_(true),  // initialized to zero
    derxy_(true),   // initialized to zero
    derxy2_(true),  // initialized to zero
    xjm_(true),     // initialized to zero
    xij_(true),     // initialized to zero
    xder2_(true),   // initialized to zero
    laplace_(true), // initialized to zero
    rhs_(numdofpernode_),       // size of vector
    reatemprhs_(numdofpernode_),// size of vector
    bodyforce_(numdofpernode_), // size of vector
    scatrares_(numscal_),  // size of vector
    conv_phi_(numscal_),   // size of vector
    diff_phi_(numscal_),   // size of vector
    rea_phi_(numscal_),    // size of vector
    tau_(numscal_),             // size of vector
    tauderpot_(numscal_),       // size of vector
    efluxreconstr_(numscal_)   // size of vector
{
  return;
}


//! compute largest element diameter for reinitialization pseudo time step size
template<DRT::Element::DiscretizationType distype, class M1>
double getEleDiameter(const M1& xyze)
{
  double elediam = 0.0;

  // number of nodes of this element
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // check all possible connections between nodes of an element
  for(size_t i_start=0; i_start< numnode-2; ++i_start)
  {
    for(size_t i_end= i_start+1; i_end < numnode-1; ++i_end)
    {
      LINALG::Matrix<3,1> direction;
      direction.Clear();
      direction(0) = xyze(0, i_start) - xyze(0, i_end);
      direction(1) = xyze(1, i_start) - xyze(1, i_end);
      direction(2) = xyze(2, i_start) - xyze(2, i_end);

      // update elediam
      if (direction.Norm2() > elediam) elediam=direction.Norm2();
    }
  }

  return elediam;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ReInitImpl<distype>::Evaluate(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra
  )
{
  // --------mandatory are performed here at first ------------
  // get node coordinates (we do this for all actions!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // get additional state vector for ALE case: grid displacement
  is_ale_ = params.get<bool>("isale",false);
  if (is_ale_)
  {
    const RCP<Epetra_MultiVector> dispnp = params.get< RCP<Epetra_MultiVector> >("dispnp",Teuchos::null);
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_);
    // add nodal displacements to point coordinates
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Set parameter SCATRATYPE in your input file!");

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch (action)
  {
  case SCATRA::reinitialize_levelset:
  {

    bool reinitswitch = params.get<bool>("reinitswitch",false);
    if(reinitswitch == false) dserror("action reinitialize_levelset should be called only with reinitswitch=true");

    const INPAR::SCATRA::PenaltyMethod reinit_penalty_method     = params.get<INPAR::SCATRA::PenaltyMethod>("reinit_penalty_method");
    const double reinit_epsilon_bandwidth                        = params.get<double>("reinit_epsilon_bandwidth",false);
    const double reinit_penalty_interface                        = params.get<double>("reinit_penalty_interface", false);
    const INPAR::SCATRA::SmoothedSignType smoothedSignType       = params.get<INPAR::SCATRA::SmoothedSignType>("reinit_smoothed_sign_type");
    const double reinit_pseudo_timestepsize_factor               = params.get<double>("reinit_pseudotimestepfactor",false);
    const INPAR::SCATRA::ReinitializationStrategy reinitstrategy = params.get<INPAR::SCATRA::ReinitializationStrategy>("reinit_strategy");


    // extract local values from the global vectors
    RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    RCP<const Epetra_Vector> phin  = discretization.GetState("phin");
    RCP<const Epetra_Vector> phi0_Reinit_Reference = discretization.GetState("phistart");

    if (phinp==Teuchos::null || phin==Teuchos::null || phi0_Reinit_Reference==Teuchos::null)
      dserror("Cannot get state vector 'phinp' or 'phi0_Reinit_Reference'");
    std::vector<double> myphinp(lm.size());
    std::vector<double> myphin(lm.size());
    std::vector<double> myphi0(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
    if (reinitswitch==true) DRT::UTILS::ExtractMyValues(*phin,myphin,lm);
    if (reinitswitch==true) DRT::UTILS::ExtractMyValues(*phi0_Reinit_Reference,myphi0,lm);


    // fill all element arrays
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephinp_[k](i,0)                 = myphinp[k+(i*numdofpernode_)];
        // the phi reference vector contains information of pseudo time step tau=0
        ephi0_Reinit_Reference_[k](i,0) = myphin[k+(i*numdofpernode_)]; //TODO: Ursula: hier steht phin, das ist aber glaube ich falsch
                                                                        // und muesste phin0 sein: das aendert aber die Ergebnisse aller Tests
                                                                        // ich schau mir das dann im Zusammenhang mit dem ls-algo genau an
        ephin_[k](i,0)                  = myphin[k+(i*numdofpernode_)];
        ephi0_penalty_[k](i,0)          = myphi0[k+(i*numdofpernode_)];
      }
    } // for i


    if (reinitstrategy == INPAR::SCATRA::reinitstrategy_pdebased_characteristic_galerkin)
    {
      // calculate element coefficient matrix and rhs
      Sysmat_Reinit_TG(
        ele,
        elemat1_epetra,
        elevec1_epetra,
        reinitswitch,
        reinit_pseudo_timestepsize_factor,
        smoothedSignType,
        reinit_epsilon_bandwidth,
        reinit_penalty_method,
        reinit_penalty_interface);
    }
    else if (reinitstrategy == INPAR::SCATRA::reinitstrategy_pdebased_linear_convection)
    {
      const double theta_reinit = params.get<double>("theta_reinit",false);
      const bool reinit_shock_capturing = params.get<int>("reinit_shock_capturing",false);
      const double reinit_shock_capturing_diffusivity = params.get<double>("reinit_shock_capturing_diffusivity",false);

      if(theta_reinit != 1.0) dserror(" correct implementation of hist_vector!!!");
      const double meshsize = getEleDiameter<distype>(xyze_);
      const double dt = reinit_pseudo_timestepsize_factor * meshsize;


      // set flag for including reactive terms to false initially
      // flag will be set to true below when reactive material is included
      is_reactive_ = false;

      // get control parameters
      is_stationary_  = params.get<bool>("using stationary formulation");
      is_genalpha_    = params.get<bool>("using generalized-alpha time integration");
      is_incremental_ = params.get<bool>("incremental solver");

      // get time factor and alpha_F if required
      // one-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
      double timefac = 1.0;
      timefac =  theta_reinit * dt;

      // set flag for conservative form
      const INPAR::SCATRA::ConvForm convform =
        DRT::INPUT::get<INPAR::SCATRA::ConvForm>(params, "form of convective term");
      is_conservative_ = false;
      if (convform ==INPAR::SCATRA::convform_conservative) is_conservative_ = true;



      // set parameters for stabilization
      Teuchos::ParameterList& stablist = params.sublist("STABILIZATION");

      // get definition for stabilization parameter tau
      whichtau_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(stablist,"DEFINITION_TAU");

      // set correct stationary definition for stabilization parameter automatically
      // and ensure that exact stabilization parameter is only used in stationary case
      if (is_stationary_)
      {
        if (whichtau_ == INPAR::SCATRA::tau_taylor_hughes_zarins)
          whichtau_ = INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt;
        else if (whichtau_ == INPAR::SCATRA::tau_franca_valentin)
          whichtau_ = INPAR::SCATRA::tau_franca_valentin_wo_dt;
        else if (whichtau_ == INPAR::SCATRA::tau_shakib_hughes_codina)
          whichtau_ = INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt;
        else if (whichtau_ == INPAR::SCATRA::tau_codina)
          whichtau_ = INPAR::SCATRA::tau_codina_wo_dt;
        else if (whichtau_ == INPAR::SCATRA::tau_franca_madureira_valentin)
          whichtau_ = INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt;
      }
      else
      {
        if (whichtau_ == INPAR::SCATRA::tau_exact_1d)
          dserror("exact stabilization parameter only available for stationary case");
      }

      // get characteristic element length for stabilization parameter definition
      charelelength_ = DRT::INPUT::IntegralValue<INPAR::SCATRA::CharEleLength>(stablist,"CHARELELENGTH");

      // set (sign) factor for diffusive and reactive stabilization terms
      // (factor is zero for SUPG) and overwrite tau definition when there
      // is no stabilization
      const INPAR::SCATRA::StabType stabinp = DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(stablist,"STABTYPE");
      switch(stabinp)
      {
      case INPAR::SCATRA::stabtype_no_stabilization:
        whichtau_ = INPAR::SCATRA::tau_zero;
        break;
      case INPAR::SCATRA::stabtype_SUPG:
        diffreastafac_ = 0.0;
        break;
      case INPAR::SCATRA::stabtype_GLS:
        diffreastafac_ = 1.0;
        break;
      case INPAR::SCATRA::stabtype_USFEM:
        diffreastafac_ = -1.0;
        break;
      default:
        dserror("unknown definition for stabilization parameter");
      }

      // set flags for potential evaluation of tau and material law at int. point
      const INPAR::SCATRA::EvalTau tauloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(stablist,"EVALUATION_TAU");
      tau_gp_ = (tauloc == INPAR::SCATRA::evaltau_integration_point); // set true/false
      const INPAR::SCATRA::EvalMat matloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(stablist,"EVALUATION_MAT");
      mat_gp_ = (matloc == INPAR::SCATRA::evalmat_integration_point); // set true/false

      // get velocity at nodes
      const RCP<Epetra_MultiVector> reinit_velocity = params.get< RCP<Epetra_MultiVector> >("reinit velocity field",Teuchos::null);
      DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,reinit_velocity,nsd_);
      const RCP<Epetra_MultiVector> reinit_convelocity = params.get< RCP<Epetra_MultiVector> >("reinit convective velocity field",Teuchos::null);
      DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,reinit_convelocity,nsd_);

      // calculate element coefficient matrix and rhs
      Sysmat_Reinit_OST(
        ele,
        elemat1_epetra,
        elevec1_epetra,
        dt,
        timefac,
        meshsize,
        reinitswitch,
        reinit_pseudo_timestepsize_factor,
        smoothedSignType,
        reinit_penalty_method,
        reinit_penalty_interface,
        reinit_epsilon_bandwidth,
        reinit_shock_capturing,
        reinit_shock_capturing_diffusivity,
        scatratype);
    }
    else dserror("reinitstrategy not a known type");

    break;
  }
  case SCATRA::calc_TG_mat_and_rhs:
  {
    // get timealgo
    const INPAR::SCATRA::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::SCATRA::TimeIntegrationScheme>(params, "timealgo");

    // get current time-step length
    const double dt   = params.get<double>("time-step length");

    // get velocity at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",Teuchos::null);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,velocity,nsd_);
    const RCP<Epetra_MultiVector> convelocity = params.get< RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,convelocity,nsd_);

    // extract local values from the global vectors
    RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    RCP<const Epetra_Vector> phin  = discretization.GetState("phin");

    if (phinp==Teuchos::null || phin==Teuchos::null)
      dserror("Cannot get state vector 'phinp' or 'phin_'");
    std::vector<double> myphinp(lm.size());
    std::vector<double> myphin(lm.size());

    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
    DRT::UTILS::ExtractMyValues(*phin,myphin,lm);


    // fill all element arrays
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
        ephin_[k](i,0) = myphin[k+(i*numdofpernode_)];
      }
    } // for i


    // calculate element coefficient matrix and rhs
    if(scatratype==INPAR::SCATRA::scatratype_levelset)
    {
      Sysmat_Transport_TG(
          ele,
          elemat1_epetra,
          elevec1_epetra,
          dt,
          timealgo);
    }
    else dserror("Sysmat_Taylor_Galerkin only for level set problems available");

    break;
  }
  case SCATRA::calc_error_reinit:
  {
    // add error only for elements which are not ghosted
    if(ele->Owner() == discretization.Comm().MyPID())
    {

      // need current solution
      RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

      // extract local values from the global vector
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);


      // fill element arrays
      for (int i=0;i<nen_;++i)
      {
        // split for each transported scalar, insert into element arrays
        for (int k = 0; k< numscal_; ++k)
        {
          ephinp_[k](i) = myphinp[k+(i*numdofpernode_)];
        }
      } // for i

      CalErrorsReinitialization(
        ele,
        params);
    }
    break;
  }
  default:
  {
    dserror("Not acting on this action. Forgot implementation?");
  }
  } // switch(action)

  // work is done
  return 0;
}



/*----------------------------------------------------------------------*
  |  get the material constants  (private)                      gjb 10/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::GetMaterialParams(
  const DRT::Element*  ele,
  const enum INPAR::SCATRA::ScaTraType  scatratype
  )
{
// get the material
  RCP<MAT::Material> material = ele->Material();
/*
// get diffusivity / diffusivities
  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());
    if (actmat->NumMat() < numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<numscal_;++k)
    {
      // set reaction coeff. and temperature rhs for reactive equation system to zero
      reacoeff_[k]   = 0.0;
      reatemprhs_[k] = 0.0;

      // set specific heat capacity at constant pressure to 1.0
      shc_ = 1.0;

      // set density at various time steps and density gradient factor to 1.0/0.0
      densn_[k]       = 1.0;
      densnp_[k]      = 1.0;
      densam_[k]      = 1.0;
      densgradfac_[k] = 0.0;

      const int matid = actmat->MatID(k);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == INPAR::MAT::m_ion)
      {
        const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(singlemat.get());
        valence_[k] = actsinglemat->Valence();
        diffus_[k] = actsinglemat->Diffusivity();
        diffusvalence_[k] = valence_[k]*diffus_[k];

        // Material data of eliminated ion species is read from the LAST ion material
        // in the matlist!
        if ((scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim) and (k==(numscal_-1)))
        {
          if (diffus_.size() == (unsigned) numscal_)
          {
            // For storing additional data, we increase the vector for
            // diffusivity and valences by one!
            std::cout<<"k = "<<k<<"   Did push back for diffus_ and valence_!"<<std::endl;
            diffus_.push_back(actsinglemat->ElimDiffusivity());
            valence_.push_back(actsinglemat->ElimValence());
            diffusvalence_.push_back(valence_[numscal_]*diffus_[numscal_]);
            // we also enlarge some other vectors by one
            tau_.push_back(0.0);
            LINALG::Matrix<nen_,1> mat(true);
            tauderpot_.push_back(mat);
          }
          else if (diffus_.size() == (unsigned) (numscal_+1))
          {
            diffus_[numscal_]  = actsinglemat->ElimDiffusivity();
            valence_[numscal_] = actsinglemat->ElimValence();
            diffusvalence_[numscal_] = valence_[numscal_]*diffus_[numscal_];
          }
          else
            dserror("Something is wrong with eliminated ion species data");
          //if (ele->Id()==0)
          //  std::cout<<"data: "<<diffus_[numscal_]<<"   "<<valence_[numscal_]<<std::endl;
          // data check:
          if (abs(diffus_[numscal_])< EPS13) dserror("No diffusivity for eliminated species read!");
          if (abs(valence_[numscal_])< EPS13) dserror("No valence for eliminated species read!");
        }
      }
      else if (singlemat->MaterialType() == INPAR::MAT::m_arrhenius_spec)
      {
        const MAT::ArrheniusSpec* actsinglemat = static_cast<const MAT::ArrheniusSpec*>(singlemat.get());

        // compute temperature
        const double tempnp = funct_.Dot(ephinp_[numscal_-1]);

        // compute diffusivity according to Sutherland law
        diffus_[k] = actsinglemat->ComputeDiffusivity(tempnp);

        // compute reaction coefficient for species equation
        reacoeff_[k] = actsinglemat->ComputeReactionCoeff(tempnp);
        reacoeffderiv_[k] = reacoeff_[k];
        // set reaction flag to true
        is_reactive_ = true;
      }
      else if (singlemat->MaterialType() == INPAR::MAT::m_arrhenius_temp)
      {
        if (k != numscal_-1) dserror("Temperature equation always needs to be the last variable for reactive equation system!");

        const MAT::ArrheniusTemp* actsinglemat = static_cast<const MAT::ArrheniusTemp*>(singlemat.get());

        // get specific heat capacity at constant pressure
        shc_ = actsinglemat->Shc();

        // compute species mass fraction and temperature
        const double spmf   = funct_.Dot(ephinp_[0]);
        const double tempnp = funct_.Dot(ephinp_[k]);

        // compute diffusivity according to Sutherland law
        diffus_[k] = actsinglemat->ComputeDiffusivity(tempnp);

        // compute density based on temperature and thermodynamic pressure
        densnp_[k] = actsinglemat->ComputeDensity(tempnp,thermpressnp_);

        if (is_genalpha_)
        {
          // compute density at n+alpha_M
          const double tempam = funct_.Dot(ephiam_[k]);
          densam_[k] = actsinglemat->ComputeDensity(tempam,thermpressam_);

          if (not is_incremental_)
          {
            // compute density at n (thermodynamic pressure approximated at n+alpha_M)
            const double tempn = funct_.Dot(ephin_[k]);
            densn_[k] = actsinglemat->ComputeDensity(tempn,thermpressam_);
          }
          else densn_[k] = 1.0;
        }
        else densam_[k] = densnp_[k];

        // factor for density gradient
        densgradfac_[k] = -densnp_[k]/tempnp;

        // compute sum of reaction rates for temperature equation divided by specific
        // heat capacity -> will be considered a right-hand side contribution
        reatemprhs_[k] = actsinglemat->ComputeReactionRHS(spmf,tempnp)/shc_;

        // set reaction flag to true
        is_reactive_ = true;
      }
      else if (singlemat->MaterialType() == INPAR::MAT::m_scatra)
      {
        const MAT::ScatraMat* actsinglemat = static_cast<const MAT::ScatraMat*>(singlemat.get());
        diffus_[k] = actsinglemat->Diffusivity();

        // in case of reaction with constant coefficient, read coefficient and
        // set reaction flag to true
        reacoeff_[k] = actsinglemat->ReaCoeff();
        if (reacoeff_[k] > EPS14) is_reactive_ = true;
        if (reacoeff_[k] < -EPS14)
          dserror("Reaction coefficient for species %d is not positive: %f",k, reacoeff_[k]);
        reacoeffderiv_[k] = reacoeff_[k];
      }
      else if (singlemat->MaterialType() == INPAR::MAT::m_biofilm)
      {
        const MAT::Biofilm* actsinglemat = static_cast<const MAT::Biofilm*>(singlemat.get());
        diffus_[k] = actsinglemat->Diffusivity();
        // double rearate_k = actsinglemat->ReaRate();
        // double satcoeff_k = actsinglemat->SatCoeff();

        // set reaction flag to true
        is_reactive_ = true;

        // get substrate concentration at n+1 or n+alpha_F at integration point
        const double csnp = funct_.Dot(ephinp_[k]);
        //const double conp = funct_.Dot(ephinp_[1]);

        // compute reaction coefficient for species equation
        reacoeff_[k] = actsinglemat->ComputeReactionCoeff(csnp);
        reacoeffderiv_[k] = actsinglemat->ComputeReactionCoeffDeriv(csnp);
      }
      else dserror("material type not allowed");

      // check whether there is negative (physical) diffusivity
      if (diffus_[k] < -EPS15) dserror("negative (physical) diffusivity");
    }
  }
  else*/ if (material->MaterialType() == INPAR::MAT::m_scatra)
  {
    const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

    // get constant diffusivity
    diffus_[0] = actmat->Diffusivity();

    // in case of reaction with (non-zero) constant coefficient:
    // read coefficient and set reaction flag to true
    reacoeff_[0] = actmat->ReaCoeff();
    if (reacoeff_[0] > EPS14) is_reactive_ = true;
    if (reacoeff_[0] < -EPS14)
      dserror("Reaction coefficient is not positive: %f",0, reacoeff_[0]);

    reacoeffderiv_[0] = reacoeff_[0];

    // set specific heat capacity at constant pressure to 1.0
    shc_ = 1.0;

    // set temperature rhs for reactive equation system to zero
    reatemprhs_[0] = 0.0;

    // set density at various time steps and density gradient factor to 1.0/0.0
    densn_[0]       = 1.0;
    densnp_[0]      = 1.0;
    densam_[0]      = 1.0;
    densgradfac_[0] = 0.0;
  }
/*  else if (material->MaterialType() == INPAR::MAT::m_ion)
  {
    const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for single ion material");

    // set reaction coeff. and temperature rhs for reactive equation system to zero
    reacoeff_[0]   = 0.0;
    reatemprhs_[0] = 0.0;
    // set specific heat capacity at constant pressure to 1.0
    shc_ = 1.0;
    // set density at various time steps and density gradient factor to 1.0/0.0
    densn_[0]       = 1.0;
    densnp_[0]      = 1.0;
    densam_[0]      = 1.0;
    densgradfac_[0] = 0.0;

    // get constant diffusivity
    diffus_[0] = actsinglemat->Diffusivity();
    valence_[0] = 0.0; // remains unused -> we only do convection-diffusion in this case!
    diffusvalence_[0] = 0.0; // remains unused
  }
  else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
  {
    const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for mixture-fraction material");

    // compute mixture fraction at n+1 or n+alpha_F
    const double mixfracnp = funct_.Dot(ephinp_[0]);

    // compute dynamic diffusivity at n+1 or n+alpha_F based on mixture fraction
    diffus_[0] = actmat->ComputeDiffusivity(mixfracnp);

    // compute density at n+1 or n+alpha_F based on mixture fraction
    densnp_[0] = actmat->ComputeDensity(mixfracnp);

    // set specific heat capacity at constant pressure to 1.0
    shc_ = 1.0;

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      const double mixfracam = funct_.Dot(ephiam_[0]);
      densam_[0] = actmat->ComputeDensity(mixfracam);

      if (not is_incremental_)
      {
        // compute density at n
        const double mixfracn = funct_.Dot(ephin_[0]);
        densn_[0] = actmat->ComputeDensity(mixfracn);
      }
      else densn_[0] = 1.0;
    }
    else densam_[0] = densnp_[0];

    // factor for density gradient
    densgradfac_[0] = -densnp_[0]*densnp_[0]*actmat->EosFacA();

    // set reaction coeff. and temperature rhs for reactive equation system to zero
    reacoeff_[0] = 0.0;
    reacoeffderiv_[0] = 0.0;
    reatemprhs_[0] = 0.0;

    // get also fluid viscosity if subgrid-scale velocity is to be included
    // or multifractal subgrid-scales are used
    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      visc_ = actmat->ComputeViscosity(mixfracnp);
  }
  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
  {
    const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for Sutherland material");

    // get specific heat capacity at constant pressure
    shc_ = actmat->Shc();

    // compute temperature at n+1 or n+alpha_F
    const double tempnp = funct_.Dot(ephinp_[0]);

    // compute diffusivity according to Sutherland law
    diffus_[0] = actmat->ComputeDiffusivity(tempnp);

    // compute density at n+1 or n+alpha_F based on temperature
    // and thermodynamic pressure
    densnp_[0] = actmat->ComputeDensity(tempnp,thermpressnp_);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      const double tempam = funct_.Dot(ephiam_[0]);
      densam_[0] = actmat->ComputeDensity(tempam,thermpressam_);

      if (not is_incremental_)
      {
        // compute density at n (thermodynamic pressure approximated at n+alpha_M)
        const double tempn = funct_.Dot(ephin_[0]);
        densn_[0] = actmat->ComputeDensity(tempn,thermpressam_);
      }
      else densn_[0] = 1.0;
    }
    else densam_[0] = densnp_[0];

    // factor for density gradient
    densgradfac_[0] = -densnp_[0]/tempnp;

    // set reaction coeff. and temperature rhs for reactive equation system to zero
    reacoeff_[0] = 0.0;
    reacoeffderiv_[0] = 0.0;
    reatemprhs_[0] = 0.0;

    // get also fluid viscosity if subgrid-scale velocity is to be included
    // or multifractal subgrid-scales are used
    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      visc_ = actmat->ComputeViscosity(tempnp);
  }
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
  {
    const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for progress-variable material");

    // get progress variable at n+1 or n+alpha_F
    const double provarnp = funct_.Dot(ephinp_[0]);

    // get specific heat capacity at constant pressure and
    // compute temperature based on progress variable
    shc_ = actmat->ComputeShc(provarnp);
    const double tempnp = actmat->ComputeTemperature(provarnp);

    // compute density at n+1 or n+alpha_F
    densnp_[0] = actmat->ComputeDensity(provarnp);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      const double provaram = funct_.Dot(ephiam_[0]);
      densam_[0] = actmat->ComputeDensity(provaram);

      if (not is_incremental_)
      {
        // compute density at n
        const double provarn = funct_.Dot(ephin_[0]);
        densn_[0] = actmat->ComputeDensity(provarn);
      }
      else densn_[0] = 1.0;
    }
    else densam_[0] = densnp_[0];

    // factor for density gradient
    densgradfac_[0] = -densnp_[0]*actmat->ComputeFactor(provarnp);

    // compute diffusivity according to Sutherland law
    diffus_[0] = actmat->ComputeDiffusivity(tempnp);

    // compute reaction coefficient for progress variable
    reacoeff_[0] = actmat->ComputeReactionCoeff(tempnp);
    reacoeffderiv_[0] = reacoeff_[0];
    // compute right-hand side contribution for progress variable
    // -> equal to reaction coefficient
    reatemprhs_[0] = reacoeff_[0];

    // set reaction flag to true
    is_reactive_ = true;

    // get also fluid viscosity if subgrid-scale velocity is to be included
    // or multifractal subgrid-scales are used
    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      visc_ = actmat->ComputeViscosity(tempnp);
  }
  else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
  {
    const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for progress-variable material");

    // get progress variable at n+1 or n+alpha_F
    const double provarnp = funct_.Dot(ephinp_[0]);

    // get specific heat capacity at constant pressure and
    // compute temperature based on progress variable
    shc_ = actmat->ComputeShc(provarnp);
    const double tempnp = actmat->ComputeTemperature(provarnp);

    // compute density at n+1 or n+alpha_F
    densnp_[0] = actmat->ComputeDensity(provarnp);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      const double provaram = funct_.Dot(ephiam_[0]);
      densam_[0] = actmat->ComputeDensity(provaram);

      if (not is_incremental_)
      {
        // compute density at n
        const double provarn = funct_.Dot(ephin_[0]);
        densn_[0] = actmat->ComputeDensity(provarn);
      }
      else densn_[0] = 1.0;
    }
    else densam_[0] = densnp_[0];

    // factor for density gradient
    densgradfac_[0] = -densnp_[0]*actmat->ComputeFactor(provarnp);

    // compute diffusivity according to Sutherland law
    diffus_[0] = actmat->ComputeDiffusivity(tempnp);

    // compute reaction coefficient for progress variable
    reacoeff_[0] = actmat->ComputeReactionCoeff(provarnp);
    reacoeffderiv_[0] = reacoeff_[0];
    // compute right-hand side contribution for progress variable
    // -> equal to reaction coefficient
    reatemprhs_[0] = reacoeff_[0];

    // set reaction flag to true
    is_reactive_ = true;

    // get also fluid viscosity if subgrid-scale velocity is to be included
    // or multifractal subgrid-scales are used
    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
      visc_ = actmat->ComputeViscosity(tempnp);
  }
  else if (material->MaterialType() == INPAR::MAT::m_biofilm)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for BIOFILM material");

    const MAT::Biofilm* actmat = static_cast<const MAT::Biofilm*>(material.get());
    diffus_[0] = actmat->Diffusivity();
    // double rearate_k = actmat->ReaRate();
    // double satcoeff_k = actmat->SatCoeff();

    // set reaction flag to true
    is_reactive_ = true;

    // get substrate concentration at n+1 or n+alpha_F at integration point
    const double csnp = funct_.Dot(ephinp_[0]);
    //const double conp = funct_.Dot(ephinp_[1]);

    // compute reaction coefficient for species equation
    reacoeff_[0] = actmat->ComputeReactionCoeff(csnp);
    reacoeffderiv_[0] = actmat->ComputeReactionCoeffDeriv(csnp);

    // set specific heat capacity at constant pressure to 1.0
    shc_ = 1.0;

    // set temperature rhs for reactive equation system to zero
    reatemprhs_[0] = 0.0;

    // set density at various time steps and density gradient factor to 1.0/0.0
    densn_[0]       = 1.0;
    densnp_[0]      = 1.0;
    densam_[0]      = 1.0;
    densgradfac_[0] = 0.0;
  }
  else if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for isotropic Fourier material");

    const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());

    // get constant diffusivity (conductivity divided by specific heat capacity)
    diffus_[0] = actmat->Conductivity()/actmat->Capacity();

    // set density at various time steps and density gradient factor to 1.0/0.0
    densn_[0]       = 1.0;
    densnp_[0]      = 1.0;
    densam_[0]      = 1.0;
    densgradfac_[0] = 0.0;

    // set specific heat capacity at constant volume
    // (value divided by density here for its intended use on right-hand side)
    shc_ = actmat->Capacity()/densnp_[0];

    // set reaction coeff. and temperature rhs for reactive equation system to zero
    reacoeff_[0]      = 0.0;
    reacoeffderiv_[0] = 0.0;
    reatemprhs_[0]    = 0.0;
  }
  else if (material->MaterialType() == INPAR::MAT::m_thermostvenant)
  {
    dsassert(numdofpernode_==1,"more than 1 dof per node for thermo St. Venant-Kirchhoff material");

    const MAT::ThermoStVenantKirchhoff* actmat = static_cast<const MAT::ThermoStVenantKirchhoff*>(material.get());

    // get constant diffusivity (conductivity divided by specific heat capacity)
    diffus_[0] = actmat->Conductivity()/actmat->Capacity();

    // set density at various time steps and set density gradient factor to 1.0/0.0
    densnp_[0]      = actmat->Density();
    densam_[0]      = densnp_[0];
    densn_[0]       = densnp_[0];
    densgradfac_[0] = 0.0;

    // set specific heat capacity at constant volume
    // (value divided by density here for its intended use on right-hand side)
    shc_ = actmat->Capacity()/densnp_[0];

    // compute reaction coefficient
    // (divided by density due to later multiplication by density in CalMatAndRHS)
    const double stmodulus = actmat->STModulus();
    reacoeff_[0] = -vdiv_*stmodulus/(actmat->Capacity()*densnp_[0]);

    // set reaction flag to true, check whether reaction coefficient is positive
    // and set derivative of reaction coefficient
    if (reacoeff_[0] > EPS14 or reacoeff_[0] < -EPS14) is_reactive_ = true;
    reacoeffderiv_[0] = reacoeff_[0];

    // set temperature rhs for reactive equation system to zero
    reatemprhs_[0] = 0.0;

    // set temporal derivative of thermodynamic pressure to zero for
    // the present structure-based scalar transport
    thermpressdt_ = 0.0;
  }
  else if (material->MaterialType() == INPAR::MAT::m_yoghurt)
  {
    const MAT::Yoghurt* actmat = static_cast<const MAT::Yoghurt*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for Yoghurt material");

    // get specific heat capacity at constant pressure
    shc_ = actmat->Shc();

    // compute diffusivity
    diffus_[0] = actmat->ComputeDiffusivity();

    // get constant density
    densnp_[0] = actmat->Density();
    densam_[0] = densnp_[0];
    densn_[0] = densnp_[0];

    // set reaction coeff. and temperature rhs for reactive equation system to zero
    reacoeff_[0] = 0.0;
    reacoeffderiv_[0] = 0.0;
    reatemprhs_[0] = 0.0;

    // get also fluid viscosity if subgrid-scale velocity is to be included
    // or multifractal subgrid-scales are used
    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      // compute temperature at n+1 or n+alpha_F
      const double tempnp = funct_.Dot(ephinp_[0]);

      // compute rate of strain
      double rateofstrain = -1.0e30;
      dserror("check this point.");
      //rateofstrain = GetStrainRate(evelnp_);

      // compute viscosity for Yoghurt-like flows according to Afonso et al. (2003)
      visc_ = actmat->ComputeViscosity(rateofstrain,tempnp);
    }
  }*/
  else dserror("Material type is not supported");

// check whether there is negative (physical) diffusivity
  if (diffus_[0] < -EPS15) dserror("negative (physical) diffusivity");

  return;
} //ReInitImpl::GetMaterialParams

/*----------------------------------------------------------------------*
  |  calculate stabilization parameter  (private)              gjb 06/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalTau(
  DRT::Element*                         ele,
  double                                diffus,
  const double                          dt,
  const double                          timefac,
  const double                          vol,
  const int                             k,
  const double                          frt,
  const bool                            migrationintau
  )
{
  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();
  // reset
  tauderpot_[k].Clear();

  //----------------------------------------------------------------------
  // computation of stabilization parameters depending on definition used
  //----------------------------------------------------------------------
  switch (whichtau_)
  {
  case INPAR::SCATRA::tau_taylor_hughes_zarins:
  case INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt:
  {
    /*

    literature:
    1) C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
    of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
    (1998) 155-196.
    2) V. Gravemeier, W.A. Wall, An algebraic variational multiscale-
    multigrid method for large-eddy simulation of turbulent variable-
    density flow at low Mach number, J. Comput. Phys. 229 (2010)
    6047-6070.
    -> version for variable-density scalar transport equation as
    implemented here, which corresponds to constant-density
    version as given in the previous publication when density
    is constant

    1
    +-                                               -+ - -
    |        2                                        |   2
    | c_1*rho                                  2      |
    tau = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
    |     2                                           |
    |   dt                                            |
    +-                                               -+

    with the constants and covariant metric tensor defined as follows:

    C   = 1.0 (not explicitly defined here),
    c_1 = 4.0,
    c_2 = 1.0 (not explicitly defined here),
    c_3 = 12.0/m_k (36.0 for linear and 144.0 for quadratic elements)

    +-           -+   +-           -+   +-           -+
    |             |   |             |   |             |
    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
    G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
    ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
    |    i     j  |   |    i     j  |   |    i     j  |
    +-           -+   +-           -+   +-           -+

    +----
    \
    G : G =   +   G   * G
    /     ij    ij
    +----
    i,j
    +----
    \
    rho*u*G*rho*u  =   +   rho*u * G  *rho*u
    /        i   ij      j
    +----
    i,j
    */
    // effective velocity at element center:
    // (weighted) convective velocity + individual migration velocity
    LINALG::Matrix<nsd_,1> veleff(convelint_,false);

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to be zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_[k];
    if (whichtau_ == INPAR::SCATRA::tau_taylor_hughes_zarins) sigma_tot += 1.0/dt;

    // computation of various values derived from covariant metric tensor
    double G;
    double normG(0.0);
    double Gnormu(0.0);
    const double dens_sqr = densnp_[k]*densnp_[k];
    for (int nn=0;nn<nsd_;++nn)
    {
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xij_(nn,0)*xij_(rr,0);
        for(int tt=1;tt<nsd_;tt++)
        {
          G += xij_(nn,tt)*xij_(rr,tt);
        }
        normG+=G*G;
        Gnormu+=dens_sqr*veleff(nn,0)*G*veleff(rr,0);
      }
    }

    // definition of constants as described above
    const double c1 = 4.0;
    const double c3 = 12.0/mk;

    // compute diffusive part
    const double Gdiff = c3*diffus*diffus*normG;

    // computation of stabilization parameter tau
    tau_[k] = 1.0/(sqrt(c1*dens_sqr*DSQR(sigma_tot) + Gnormu + Gdiff));

  }
  break;
  case INPAR::SCATRA::tau_franca_valentin:
  {
    /*

    literature:
    L.P. Franca, F. Valentin, On an improved unusual stabilized
    finite element method for the advective-reactive-diffusive
    equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.


    xi1,xi2 ^
    |      /
    |     /
    |    /
    1 +---+
    |
    |
    |
    +--------------> re1,re2
    1

    */
    // get Euclidean norm of (weighted) velocity at element center
    double vel_norm;
    // dserror("FrancaValentin with migrationintau not available at the moment");
    /*
    // get Euclidean norm of effective velocity at element center:
    // (weighted) convective velocity + individual migration velocity
    LINALG::Matrix<nsd_,1> veleff(velint_,false);

    veleff.Update(diffusvalence_[k],migvelint_,1.0);
    vel_norm = veleff.Norm2();

    #ifdef VISUALIZE_ELEMENT_DATA
    veleff.Update(diffusvalence_[k],migvelint_,0.0);
    double vel_norm_mig = veleff.Norm2();
    double migepe2 = mk * vel_norm_mig * h / diffus;

    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
    if (!actele) dserror("cast to Transport* failed");
    std::vector<double> v(1,migepe2);
    std::ostringstream temp;
    temp << k;
    std::string name = "Pe_mig_"+temp.str();
    actele->AddToData(name,v);
    name = "hk_"+temp.str();
    v[0] = h;
    actele->AddToData(name,v);
    #endif
    }
    else*/
    vel_norm = convelint_.Norm2();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to be zero in GetMaterialParams for non-reactive material)
    const double sigma_tot = 1.0/timefac + reacoeff_[k];

    // calculate characteristic element length
    const double h = CalcCharEleLength(vol,vel_norm);

    // various parameter computations:
    // relating convective to viscous part
    if (diffus < EPS14) dserror("Invalid diffusion coefficent");
    const double epe = mk * densnp_[k] * vel_norm * h / diffus;
    // relating viscous to reactive part
    const double epe1 = 2.0*diffus/(mk*densnp_[k]*sigma_tot*DSQR(h));

    // respective "switching" parameters
    const double xi  = std::max(epe,1.0);
    const double xi1 = std::max(epe1,1.0);

    tau_[k] = DSQR(h)/(DSQR(h)*densnp_[k]*sigma_tot*xi1 + 2.0*diffus*xi/mk);

#ifdef VISUALIZE_ELEMENT_DATA
    // visualize resultant Pe number
    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
    if (!actele) dserror("cast to Transport* failed");
    std::vector<double> v(1,epe);
    std::ostringstream temp;
    temp << k;
    std::string name = "Pe_"+temp.str();
    actele->AddToData(name,v);
#endif
  }
  break;
  case INPAR::SCATRA::tau_franca_valentin_wo_dt:
  {
    /*

    stabilization parameter as above without inclusion of dt-part

    */
    // get Euclidean norm of (weighted) velocity at element center
    double vel_norm;
    // dserror("FrancaValentin with migrationintau not available at the moment");
    /*
    // get Euclidean norm of effective velocity at element center:
    // (weighted) convective velocity + individual migration velocity
    LINALG::Matrix<nsd_,1> veleff(velint_,false);

    veleff.Update(diffusvalence_[k],migvelint_,1.0);
    vel_norm = veleff.Norm2();

    #ifdef VISUALIZE_ELEMENT_DATA
    veleff.Update(diffusvalence_[k],migvelint_,0.0);
    double vel_norm_mig = veleff.Norm2();
    double migepe2 = mk * vel_norm_mig * h / diffus;

    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
    if (!actele) dserror("cast to Transport* failed");
    std::vector<double> v(1,migepe2);
    std::ostringstream temp;
    temp << k;
    std::string name = "Pe_mig_"+temp.str();
    actele->AddToData(name,v);
    name = "hk_"+temp.str();
    v[0] = h;
    actele->AddToData(name,v);
    #endif
    }
    else*/
    vel_norm = convelint_.Norm2();

    // calculate characteristic element length
    const double h = CalcCharEleLength(vol,vel_norm);

    // various parameter computations for case without dt:
    // relating convective to viscous part
    if (diffus < EPS14) dserror("Invalid diffusion coefficent");
    const double epe = mk * densnp_[k] * vel_norm * h / diffus;
    // relating viscous to reactive part
    double epe1 = 0.0;
    if (is_reactive_) epe1 = 2.0*diffus/(mk*densnp_[k]*reacoeff_[k]*DSQR(h));

    // respective "switching" parameters
    const double xi  = std::max(epe,1.0);
    const double xi1 = std::max(epe1,1.0);

    tau_[k] = DSQR(h)/(DSQR(h)*densnp_[k]*reacoeff_[k]*xi1 + 2.0*diffus*xi/mk);

#ifdef VISUALIZE_ELEMENT_DATA
    // visualize resultant Pe number
    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
    if (!actele) dserror("cast to Transport* failed");
    std::vector<double> v(1,epe);
    std::ostringstream temp;
    temp << k;
    std::string name = "Pe_"+temp.str();
    actele->AddToData(name,v);
#endif
  }
  break;
  case INPAR::SCATRA::tau_shakib_hughes_codina:
  case INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt:
  {
    /*

    literature:
    1) F. Shakib, Finite element analysis of the compressible Euler and
    Navier-Stokes equations, PhD thesis, Division of Applied Mechanics,
    Stanford University, Stanford, CA, USA, 1989.
    2) F. Shakib, T.J.R. Hughes, A new finite element formulation for
    computational fluid dynamics: IX. Fourier analysis of space-time
    Galerkin/least-squares algorithms, Comput. Methods Appl. Mech.
    Engrg. 87 (1991) 35-58.
    3) R. Codina, Stabilized finite element approximation of transient
    incompressible flows using orthogonal subscales, Comput. Methods
    Appl. Mech. Engrg. 191 (2002) 4295-4321.

    All those proposed definitions were for non-reactive incompressible
    flow; they are adapted to potentially reactive scalar transport
    equations with potential density variations here.

    constants defined as in Shakib (1989) / Shakib and Hughes (1991),
    merely slightly different with respect to c_3:

    c_1 = 4.0,
    c_2 = 4.0,
    c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

    Codina (2002) proposed present version without dt and explicit
    definition of constants
    (condition for constants as defined here: c_2 <= sqrt(c_3)).

    */
    // get Euclidean norm of velocity
    const double vel_norm = convelint_.Norm2();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to be zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_[k];
    if (whichtau_ == INPAR::SCATRA::tau_shakib_hughes_codina) sigma_tot += 1.0/dt;

    // calculate characteristic element length
    const double h = CalcCharEleLength(vol,vel_norm);

    // definition of constants as described above
    const double c1 = 4.0;
    const double c2 = 4.0;
    const double c3 = 4.0/(mk*mk);
    // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

    tau_[k] = 1.0/(sqrt(c1*DSQR(densnp_[k])*DSQR(sigma_tot)
                        + c2*DSQR(densnp_[k])*DSQR(vel_norm)/DSQR(h)
                        + c3*DSQR(diffus)/(DSQR(h)*DSQR(h))));
  }
  break;
  case INPAR::SCATRA::tau_codina:
  case INPAR::SCATRA::tau_codina_wo_dt:
  {
    /*

    literature:
    R. Codina, Comparison of some finite element methods for solving
    the diffusion-convection-reaction equation, Comput. Methods
    Appl. Mech. Engrg. 156 (1998) 185-210.

    constants:
    c_1 = 1.0,
    c_2 = 2.0,
    c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)

    Codina (1998) proposed present version without dt.

    */
    // get Euclidean norm of velocity
    const double vel_norm = convelint_.Norm2();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to be zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_[k];
    if (whichtau_ == INPAR::SCATRA::tau_codina) sigma_tot += 1.0/dt;

    // calculate characteristic element length
    const double h = CalcCharEleLength(vol,vel_norm);

    // definition of constants as described above
    const double c1 = 1.0;
    const double c2 = 2.0;
    const double c3 = 4.0/mk;

    tau_[k] = 1.0/(c1*densnp_[k]*sigma_tot
                   + c2*densnp_[k]*vel_norm/h
                   + c3*diffus/(h*h));
  }
  break;
  case INPAR::SCATRA::tau_franca_madureira_valentin:
  case INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt:
  {
    /*

    This stabilization parameter is only intended to be used for
    reactive-diffusive problems such as structure-based scalar
    transport problems in case of potentially dominating reaction.

    literature:
    L.P. Franca, A.L. Madureira, F. Valentin, Towards multiscale
    functions: enriching finite element spaces with local but not
    bubble-like functions, Comput. Methods Appl. Mech. Engrg. 194
    (2005) 3006-3021.

    */
    // get Euclidean norm of velocity at element center
//    double vel_norm;
//    vel_norm = convelint_.Norm2();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to be zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_[k];
    if (whichtau_ == INPAR::SCATRA::tau_franca_madureira_valentin)
      sigma_tot += 1.0/timefac;

    // calculate characteristic element length
    // -> currently: cubic/square root of element volume/area or
    //    element length (3-/2-/1-D)
    // cast dimension to a double variable -> pow()
    const double dim = (double) nsd_;
    const double h = std::pow(vol,1/dim);


    // parameter relating reactive to diffusive part
    const double epe = 2.0*diffus/(mk*densnp_[k]*sigma_tot*DSQR(h));

    // respective "switching" parameter
    const double xi = std::max(epe,1.0);

    // constant c_u as suggested in Badia and Codina (2010), method A
    // is set to be 1.0 here as in Franca et al. (2005)
    // alternative: 4.0 as suggested in Badia and Codina (2010) for
    // Darcy flow
    const double c_u = 1.0;

    tau_[k] = DSQR(h)/(c_u*DSQR(h)*densnp_[k]*sigma_tot*xi + (2.0*diffus/mk));
  }
  break;
  case INPAR::SCATRA::tau_exact_1d:
  {
    // get number of dimensions (convert from int to double)
    const double dim = (double) nsd_;

    // get characteristic element length
    double h = std::pow(vol,(1.0/dim)); // equals streamlength in 1D

    // get Euclidean norm of (weighted) velocity at element center
    double vel_norm(0.0);
    vel_norm = convelint_.Norm2();

    if (diffus < EPS14) dserror("Invalid diffusion coefficent");
    double epe = 0.5 * densnp_[k] * vel_norm * h / diffus;

    const double pp = exp(epe);
    const double pm = exp(-epe);
    double xi = 0.0;
    if (epe >= 700.0)
      tau_[k] = 0.5*h/vel_norm;
    else if (epe < 700.0 and epe > EPS15)
    {
      xi = (((pp+pm)/(pp-pm))-(1.0/epe)); // xi = coth(epe) - 1/epe
      // compute optimal stabilization parameter
      tau_[k] = 0.5*h*xi/vel_norm;

#if 0
      std::cout<<"epe = "<<epe<<std::endl;
      std::cout<<"xi_opt  = "<<xi<<std::endl;
      std::cout<<"vel_norm  = "<<vel_norm<<std::endl;
      std::cout<<"tau_opt = "<<tau_[k]<<std::endl<<std::endl;
#endif
    }
    else tau_[k] = 0.0;
  }
  break;
  case INPAR::SCATRA::tau_zero:
  {
    // set tau's to zero (-> no stabilization effect)
    tau_[k] = 0.0;
  }
  break;
  default: dserror("unknown definition for stabilization parameter tau\n");
  } //switch (whichtau_)

#if 0
  std::cout<<"diffus  for k "<<k <<" is = "<<diffus<<std::endl;
#endif
#ifdef VISUALIZE_ELEMENT_DATA
  // visualize stabilization parameter
  DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
  if (!actele) dserror("cast to Transport* failed");
  std::vector<double> v(1,tau_[k]);
  std::ostringstream temp;
  temp << k;
  std::string name = "tau_"+ temp.str();
  actele->AddToData(name,v);
#endif

  return;
} //ReInitImpl::CalTau


/*----------------------------------------------------------------------*
  |  calculation of characteristic element length               vg 01/11 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ReInitImpl<distype>::CalcCharEleLength(
  const double  vol,
  const double  vel_norm
  )
{
  // define and initialize streamlength
  double h = 0.0;
  
  //---------------------------------------------------------------------
  // select from various definitions for characteristic element length
  //---------------------------------------------------------------------
  switch (charelelength_)
  {
    // a) streamlength due to Tezduyar et al. (1992) -> default
    // normed velocity vector
    case INPAR::SCATRA::streamlength:
    {
      LINALG::Matrix<nsd_,1> velino(true);
      if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,convelint_);
      else
      {
        velino.Clear();
        velino(0,0) = 1.0;
      }

      // get streamlength using the normed velocity at element centre
      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_,velino);
      const double val = tmp.Norm1();
      h = 2.0/val; // h=streamlength
    }
    break;

    // b) volume-equivalent diameter (warning: 3-D formula!)
    case INPAR::SCATRA::volume_equivalent_diameter:
    {
      h = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
    }
    break;

    // c) cubic/square root of element volume/area or element length (3-/2-/1-D)
    case INPAR::SCATRA::root_of_volume:
    {
      // cast dimension to a double varibale -> pow()
      const double dim = double (nsd_);
      h = std::pow(vol,1/dim);
    }
    break;

    default: dserror("unknown characteristic element length\n");
    break;
  } //switch (charelelength_)

  return h;
} // ReInitImpl<distype>::CalcCharEleLength()



/*----------------------------------------------------------------------*
  | evaluate shape functions and derivatives at int. point     gjb 08/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ReInitImpl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
  const int                                    iquad,      ///< id of current Gauss point
  const int                                    eleid       ///< the element id
  )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  if (use2ndderiv_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
  }

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */

  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad]*det;

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  // compute second global derivatives (if needed)
  if (use2ndderiv_)
  {
    // get global second derivatives
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else
    derxy2_.Clear();

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

} //ScaTraImpl::CalcSubgrVelocity


/*----------------------------------------------------------------------*
 |  calculate the Laplacian for all shape functions (strong form)       |
 |                                                  (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ReInitImpl<distype>::GetLaplacianStrongForm(
  LINALG::Matrix<nen_,1>& diff,
  const LINALG::Matrix<numderiv2_,nen_>& deriv2
  )
{
  diff.Clear();
  // compute N,xx  +  N,yy +  N,zz for each shape function at integration point
  for (int i=0; i<nen_; ++i)
  {
    for (int j = 0; j<nsd_; ++j)
    {
      diff(i) += deriv2(j,i);
    }
  }
  return;
}; // ReInitImpl<distype>::GetLaplacianStrongForm

/*----------------------------------------------------------------------*
 |  calculate the Laplacian (weak form)             (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ReInitImpl<distype>::GetLaplacianWeakForm(
  double& val,
  const LINALG::Matrix<nsd_,nen_>& derxy,
  const int vi,
  const int ui)
{
  val = 0.0;
  for (int j = 0; j<nsd_; j++)
  {
    val += derxy(j, vi)*derxy(j, ui);
  }
  return;
}; // ReInitImpl<distype>::GetLaplacianWeakForm

/*----------------------------------------------------------------------*
 |  calculate rhs of Laplacian (weak form)          (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ReInitImpl<distype>::GetLaplacianWeakFormRHS(
    double& val,
    const LINALG::Matrix<nsd_,nen_>& derxy,
    const LINALG::Matrix<nsd_,1>&   gradphi,
    const int vi)
{
  val = 0.0;
  for (int j = 0; j<nsd_; j++)
  {
    val += derxy(j,vi)*gradphi(j);
  }
  return;
}; // ReInitImpl<distype>::GetLaplacianWeakFormRHS


/*----------------------------------------------------------------------*
 |  calculate divergence of vector field (e.g., velocity)  (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ReInitImpl<distype>::GetDivergence(
  double&                          vdiv,
  const LINALG::Matrix<nsd_,nen_>& evel,
  const LINALG::Matrix<nsd_,nen_>& derxy)
{
  LINALG::Matrix<nsd_,nsd_> vderxy;
  vderxy.MultiplyNT(evel,derxy);

  vdiv = 0.0;
  // compute vel x,x  + vel y,y +  vel z,z at integration point
  for (int j = 0; j<nsd_; ++j)
  {
    vdiv += vderxy(j,j);
  }
  return;
};


/*----------------------------------------------------------------------*
  |  calculate system matrix and rhs for TaylorGalerkin scheme          |
  |  just for the pure transport equation                   schott 05/11|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::Sysmat_Transport_TG(
    DRT::Element*                         ele,               ///< the element those matrix is calculated
    Epetra_SerialDenseMatrix&             emat,              ///< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,              ///< element rhs to calculate
    const double                          dt,                ///< current time-step length
    const enum INPAR::SCATRA::TimeIntegrationScheme timealgo ///< time algorithm
)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,ele->Id());

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  // Assemble element rhs and vector for domain!!! integrals
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // get velocity at integration point
    velint_.Multiply(evelnp_,funct_);
    convelint_.Multiply(econvelnp_,funct_);

    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // REMARK:
      // the bodyforce vector is evaluated at each Gaussian point as a nonlinear function
      // a node-wise rhs_-vector for the reinitialization bodyforce (smoothed sign-function) leads not the desired results

      // compute matrix and rhs
      if(timealgo == INPAR::SCATRA::timeint_tg2)
      {
        // implicit characteristic galerkin scheme 2nd order
        CalMatAndRHS_Transport_ICG2(emat, erhs, fac, k, ele, dt );
      }
      else if(timealgo == INPAR::SCATRA::timeint_tg3)
      {
        // explicit Taylor-Galerkin scheme 3rd order
        CalMatAndRHS_Transport_TG3(emat, erhs, fac, k, ele, dt );
      }
      else dserror("no characteristic/Taylor Galerkin method chosen here");

    } // loop over each scalar
  } // integration loop

  return;
}



/*-------------------------------------------------------------------------------*
  |  evaluate element matrix and rhs for transport equation based                 |
  |  on an implicit characteristic Galerkin scheme 2nd order          schott 12/10|
  *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalMatAndRHS_Transport_ICG2(
    Epetra_SerialDenseMatrix&  emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector&  erhs,     //!< element rhs to calculate
    const double               fac,      //!< integration factor
    const int                  k,        //!< index of current scalar
    DRT::Element*              ele,      //!< the element we are dealing with
    const double               dt        //!< current time-step length
  )
{
  //==========================================================
  // evaluate element vectors and gradients
  //==========================================================
  // get element vectors
  // dist_npi   distance vector for reinitialization at current timestep np and old increment i
  // dist_n     distance vector for reinitialization at old timestep n
  // phi_0      reference phi vector for smoothed sign function -> needed for directed transport along characteristics
  const double dist_n   = funct_.Dot(ephin_[k]);   // d^n
  const double dist_npi = funct_.Dot(ephinp_[k]);  // d^(n+1)_i


  // get gradients and norms
  // dist_npi   distance vector for reinitialization at current timestep np and old increment i
  // dist_n     distance vector for reinitialization at old timestep n
  // phi_0      reference phi vector for smoothed sign function -> needed for directed transport along characteristics

  LINALG::Matrix<nsd_,1> grad_dist_n(true);
  grad_dist_n.Multiply(derxy_,ephin_[k]);

  LINALG::Matrix<nsd_,1> grad_dist_npi(true);
  grad_dist_npi.Multiply(derxy_,ephinp_[k]);

  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------

  //     |           |
  //     | w, D(psi) |
  //     |           |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      emat(fvi,fui) += funct_(vi)*fac*funct_(ui);
    }
  }


  //               |                           |
  //   1/4   dt^2  | u*grad(v), u*grad(D(psi)) |
  //               |                           |

  LINALG::Matrix<nen_,1> uGradDphi(true);
  uGradDphi.MultiplyTN(derxy_,convelint_);


  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      emat(fvi,fui) += uGradDphi(vi,0)*uGradDphi(ui,0)* (fac*dt*dt/4.0);
    }
  }


  //----------  --------------    |          n+1     n    |
  //  rhs                       - | w, u*(phi   - phi  )  |
  //--------------------------    |          i            |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= funct_(vi)*fac*(dist_npi-dist_n);
  }


  //----------  --------------                |              n   |
  //  rhs                                 -dt | w, u*grad(phi )  |
  //--------------------------                |                  |

  LINALG::Matrix<1,1> uGradPhi(true);
  uGradPhi.MultiplyTN(convelint_,grad_dist_n);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= funct_(vi)*dt*fac*uGradPhi(0,0);
  }

  //               |                       n+1     n  |
  //    -dt*dt/4.0 |  grad(w)*u, u*grad(psi   + psi ) |
  //               |                       i          |

  LINALG::Matrix<nsd_,1> sum_phi(true);
  sum_phi.Update(1.0,grad_dist_npi, 1.0, grad_dist_n);

  LINALG::Matrix<1,1> uGradSumPhi(true);
  uGradSumPhi.MultiplyTN(convelint_,sum_phi);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= uGradDphi(vi,0)*dt*dt/4.0*fac*uGradSumPhi(0,0);
  }

  return;
} //ReInitImpl::CalMatAndRHS_ICG2



/*-------------------------------------------------------------------------------*
  |  evaluate element matrix and rhs for transport equation with TG3  schott 12/10|
  *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalMatAndRHS_Transport_TG3(
    Epetra_SerialDenseMatrix&  emat,     //!< element matrix to calculate
    Epetra_SerialDenseVector&  erhs,     //!< element rhs to calculate
    const double               fac,      //!< integration factor
    const int                  k,        //!< index of current scalar
    DRT::Element*              ele,      //!< the element we are dealing with
    const double               dt        //!< current time-step length
  )
{

  //==========================================================
  // evaluate element vectors and gradients
  //==========================================================
  // get element vectors
  // dist_npi   distance vector for reinitialization at current timestep np and old increment i
  // dist_n     distance vector for reinitialization at old timestep n
  // phi_0      reference phi vector for smoothed sign function -> needed for directed transport along characteristics
  const double dist_n   = funct_.Dot(ephin_[k]);   // d^n
  const double dist_npi = funct_.Dot(ephinp_[k]);  // d^(n+1)_i


  // get gradients and norms
  // dist_npi   distance vector for reinitialization at current timestep np and old increment i
  // dist_n     distance vector for reinitialization at old timestep n
  // phi_0      reference phi vector for smoothed sign function -> needed for directed transport along characteristics

  LINALG::Matrix<nsd_,1> grad_dist_n(true);
  grad_dist_n.Multiply(derxy_,ephin_[k]);

  LINALG::Matrix<nsd_,1> grad_dist_npi(true);
  grad_dist_npi.Multiply(derxy_,ephinp_[k]);

  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------

  //     |           |
  //     | w, D(phi) |
  //     |           |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      emat(fvi,fui) += funct_(vi)*fac*funct_(ui);
    }
  }

  //               |                           |
  //    1/6 *dt^2  | a*grad(w), a*grad(D(phi)) |
  //               |                           |

  // a*grad(w) bzw. a*grad(D(phi))
  LINALG::Matrix<nen_,1> aGradD(true);
  aGradD.MultiplyTN(derxy_,convelint_);

  // a*grad(phi_n)
  double a_phi_n = convelint_.Dot(grad_dist_n);

  // a*grad(phi_npi)
  double a_phi_npi = convelint_.Dot(grad_dist_npi);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      emat(fvi,fui) += aGradD(vi)*aGradD(ui)* (fac*dt*dt/6.0);
    }
  }

  //----------  --------------    |                n    |
  //  rhs                     +dt | a*grad(w) , phi     |
  //--------------------------    |                     |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] += dt*fac*aGradD(vi)*dist_n;
  }

  //----------  --------------    |        n+1     n|
  //  rhs                       - | w , phi    -phi |
  //--------------------------    |        i        |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= fac*funct_(vi)*(dist_npi-dist_n);
  }


  //----------  --------------                |                      n   |
  //  rhs                         -1/2*dt*dt* | a*grad(w), a*grad(phi )  |
  //--------------------------                |                          |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= aGradD(vi)*dt*dt*0.5*fac*a_phi_n;
  }

  //----------  --------------                |                      n+1    n    |
  //  rhs                         -1/6*dt*dt* | a*grad(w), a*grad(phi   -phi  )  |
  //--------------------------                |                      i           |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= aGradD(vi)*dt*dt/6.0*fac*(a_phi_npi-a_phi_n);
  }

  return;
} //ReInitImpl::CalMatAndRHS_TG3


/*-------------------------------------------------------------------------------*
  |  evaluate a smoothed sign function for reinitialization          schott 12/10|
  *-------------------------------------------------------------------------------*/
double EvaluateSmoothedSign(
  double                             phi_0,
  double                             grad_norm_phi_0,
  const double                       epsilon_bandwidth,
  const double                       mesh_size,
  INPAR::SCATRA::SmoothedSignType    smoothedSignType
  )
{
  //==========================================================
  // evaluate smoothed sign function for reference phi-field
  //==========================================================

  double sign_phi_0 = 0.0;

  if(smoothedSignType == INPAR::SCATRA::signtype_nonsmoothed)
  {
    if      (phi_0 < 0.0) sign_phi_0 = -1.0;
    else if (phi_0 > 0.0) sign_phi_0 = +1.0;
    else                  sign_phi_0 =  0.0;
  }
  else if (smoothedSignType == INPAR::SCATRA::signtype_LinEtAl2005)
  {
    double alpha = epsilon_bandwidth * mesh_size;
    sign_phi_0 = phi_0/sqrt(phi_0*phi_0 + alpha*alpha * grad_norm_phi_0*grad_norm_phi_0);
  }
  else if (smoothedSignType == INPAR::SCATRA::signtype_LinEtAl_normalized)
  {
    double alpha = epsilon_bandwidth * mesh_size;
    sign_phi_0 = phi_0/sqrt(phi_0*phi_0 + alpha*alpha);
  }
  else if (smoothedSignType == INPAR::SCATRA::signtype_Nagrath2005)
  {
    double alpha = epsilon_bandwidth * mesh_size;
    if(fabs(alpha) < 1e-15) dserror("divide by zero in evaluate for smoothed sign function");

    if (phi_0 < -alpha)       sign_phi_0 = -1.0;
    else if (phi_0 > alpha)   sign_phi_0 = +1.0;
    else                      sign_phi_0 = (1.0 + phi_0/alpha + 1.0/PI*sin(PI*phi_0/alpha))-1.0;
  }
  else dserror("unknown type of smoothed sign function!");


  return sign_phi_0;
}

/*-------------------------------------------------------------------------------*
  |  evaluate the derivative of a smoothed sign function             schott 12/10|
  *-------------------------------------------------------------------------------*/
double EvaluateSmoothedSignDeriv(
  double                             phi_0,
  const double                       epsilon_bandwidth,
  const double                       mesh_size,
  INPAR::SCATRA::SmoothedSignType    smoothedSignType
  )
{
  //==========================================================
  // evaluate smoothed sign function for reference phi-field
  //==========================================================

  double deriv_smoothed_Heavyside = 0.0;

  if(smoothedSignType == INPAR::SCATRA::signtype_nonsmoothed)
  {
    deriv_smoothed_Heavyside = 0.0;
  }
  else if (smoothedSignType == INPAR::SCATRA::signtype_LinEtAl2005)
  {
    // double alpha = epsilon_bandwidth * mesh_size;
    // sign_phi_0 = phi_0/sqrt(phi_0*phi_0 + alpha*alpha * grad_norm_phi_0*grad_norm_phi_0);
    dserror("derivative of smoothed Heavyside function not implemented yet");
  }
  else if (smoothedSignType == INPAR::SCATRA::signtype_LinEtAl_normalized)
  {
    // double alpha = epsilon_bandwidth * mesh_size;
    // sign_phi_0 = phi_0/sqrt(phi_0*phi_0 + alpha*alpha);
    dserror("derivative of smoothed Heavyside function not implemented yet");
  }
  else if (smoothedSignType == INPAR::SCATRA::signtype_Nagrath2005)
  {
    double alpha = epsilon_bandwidth * mesh_size;
    if(fabs(alpha) < 1e-15) dserror("divide by zero in evaluate for smoothed sign function");

    if (phi_0 < -alpha)       deriv_smoothed_Heavyside = 0.0;
    else if (phi_0 > alpha)   deriv_smoothed_Heavyside = 0.0;
    else                      deriv_smoothed_Heavyside = 1.0/(2.0*alpha)*(1.0 + cos(PI*phi_0/alpha));
  }
  else dserror("unknown type of smoothed sign function!");


  return deriv_smoothed_Heavyside;
}

/*----------------------------------------------------------------------*
  | evaluate shape functions and derivatives at int. point  schott 01/11 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ReInitImpl<distype>::EvalShapeFuncAndDerivsAtIntPointREINIT(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
  const int                                    iquad,      ///< id of current Gauss point
  const int                                    eleid       ///< the element id
  )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  if (not DRT::NURBS::IsNurbs(distype))
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
    if (use2ndderivReinitialization_)
    {
      // get the second derivatives of standard element at current GP
      DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
    }
  }
  else // nurbs elements are always somewhat special...
  {
    dserror("what to do in case of nurbs?");
  }

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */

  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad]*det;

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  // compute second global derivatives (if needed)
  if (use2ndderivReinitialization_)
  {
    // get global second derivatives
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else
    derxy2_.Clear();

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

}

/*---------------------------------------------------------------------*
  |  calculate error for reinitialization                   schott 12/10|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalErrorsReinitialization(
    const DRT::Element*          ele,   //!< the element
    Teuchos::ParameterList&      params //!< parameter list
  )
{
  // evaluate the error only within the transition region of the smoothed heavyside function
#define L1_NORM_TRANSITION_REGION

  //==================================== REINITIALIZATION error calculation  ========================
  // gradient norm of phi || grad(phi) -1.0 ||_L1(Omega)                                 schott 12/10
  //=================================================================================================

  // get element params
  double eleL1gradienterr  = params.get<double>("L1 integrated gradient error"); //squared errors
  double elevolume         = params.get<double>("volume");                       // non-squared volume

  int k = 0; // we assume only one scalar

  // get Gaussian points for integrated L2-norm and volume calculation
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_reinit(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  // caculate element wise errors and volume
  for (int iquad=0; iquad<intpoints_reinit.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints_reinit,iquad,ele->Id());

#ifdef L1_NORM_TRANSITION_REGION
    const double mesh_size = getEleDiameter<distype>(xyze_);

    double phi_ref = funct_.Dot(ephinp_[k]);
    double deriv_smoothed_Heavyside = EvaluateSmoothedSignDeriv(phi_ref, 3,mesh_size, INPAR::SCATRA::signtype_Nagrath2005);
#endif

    LINALG::Matrix<nsd_,1> gradphi;
    gradphi.Clear();
    gradphi.Multiply(derxy_,ephinp_[k]);
    double gradphi_norm = gradphi.Norm2();
    double gradphi_norm_err = gradphi_norm - 1.0;

    // integrate L1-error ( || grad(phi)-1.0 || )_L1(Omega_ele)

#ifdef L1_NORM_TRANSITION_REGION
    eleL1gradienterr += fabs(gradphi_norm_err) * fac* deriv_smoothed_Heavyside;
    // integrate volume
    elevolume        += fac * deriv_smoothed_Heavyside;
#else
    eleL1gradienterr += fabs(gradphi_norm_err) * fac;

    // integrate volume
    elevolume        += fac;
#endif
  } // Gaussian points for correction factor

  // set new element params
  params.set<double>("L1 integrated gradient error", eleL1gradienterr);
  params.set<double>("volume", elevolume);


  return;
}



/*-------------------------------------------------------------------------------*
  |  evaluate element matrix and rhs for reinitialization equation based          |
  |  on an implicit characteristic Galerkin scheme 2nd order          schott 12/10|
  *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalMatAndRHS_REINIT_ICG2(
    Epetra_SerialDenseMatrix&             emat,                     //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                     //!< element rhs to calculate
    const double                          fac,                      //!< integration factor
    const int                             k,                        //!< index of current scalar
    DRT::Element*                         ele,                      //!< the element we are dealing with
    double                                pseudo_timestep_size,     //!< pseudo time step size
    double                                mesh_size,                //!< meshsize
    const INPAR::SCATRA::PenaltyMethod    penalty_method,           //!< penalty method to keep the interface the interface position
    const double                          penalty_interface_reinit, //!< penalty parameter to keep the interface position
    const double                          epsilon_bandwidth,        //!< epsilon bandwith for the smoothed sign function
    const INPAR::SCATRA::SmoothedSignType smoothedSignType          //!< type of smoothing the sign function
  )
{

//==========================================================
// evaluate element vectors and gradients
//==========================================================
// get element vectors
// dist_npi   distance vector for reinitialization at current timestep np and old increment i
// dist_n     distance vector for reinitialization at old timestep n
// phi_0      reference phi vector for smoothed sign function -> needed for directed transport along characteristics
  const double dist_n   = funct_.Dot(ephin_[k]);   // d^n
  const double dist_npi = funct_.Dot(ephinp_[k]);  // d^(n+1)_i
  const double phi_0    = funct_.Dot(ephi0_Reinit_Reference_[k]);

// get gradients and norms
// grad_dist_npi   distance vector for reinitialization at current timestep np and old increment i
// grad_dist_n     distance vector for reinitialization at old timestep n
// grad_phi_0      reference phi vector for smoothed sign function -> needed for directed transport along characteristics

  LINALG::Matrix<nsd_,1> grad_dist_n(true);
  grad_dist_n.Multiply(derxy_,ephin_[k]);

  LINALG::Matrix<nsd_,1> grad_dist_npi(true);
  grad_dist_npi.Multiply(derxy_,ephinp_[k]);

  LINALG::Matrix<nsd_,1> grad_phi_0(true);
  grad_phi_0.Multiply(derxy_,ephi0_Reinit_Reference_[k]);

  double grad_norm_dist_n   = grad_dist_n.Norm2();
  double grad_norm_phi_0    = grad_phi_0.Norm2();


// get 2nd order derivatives
  LINALG::Matrix<numderiv2_,1> second_dist_n(true);
  second_dist_n.Multiply(derxy2_,ephin_[k]);

  LINALG::Matrix<numderiv2_,1> second_dist_npi(true);
  second_dist_npi.Multiply(derxy2_,ephinp_[k]);


  double sign_phi_0 = EvaluateSmoothedSign(phi_0, grad_norm_phi_0, epsilon_bandwidth, mesh_size, smoothedSignType);


  if(penalty_method == INPAR::SCATRA::penalty_method_akkerman)
  {
    // get derivative of smoothed sign_function
    //evaluate phinp and phi_0
    double phinp = funct_.Dot(ephinp_[k]);
    double phi_ref = funct_.Dot(ephi0_Reinit_Reference_[k]);

    double deriv_smoothed_Heavyside = EvaluateSmoothedSignDeriv(phi_ref, epsilon_bandwidth,mesh_size, smoothedSignType);

    const double densfac_penalty = pseudo_timestep_size*fac*densnp_[k]*penalty_interface_reinit*deriv_smoothed_Heavyside;
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densfac_penalty*funct_(vi);
      const int fvi = vi*numdofpernode_+k;
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;
        emat(fvi,fui) += v*funct_(ui);
      }
    }

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;
      erhs[fvi] -= funct_(vi)*densfac_penalty*( phinp - phi_ref) ;
    }
  }


  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------

  //     |           |
  //     | w, D(psi) |
  //     |           |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += funct_(vi)*fac*funct_(ui);
    }
  }


  //               |                       |
  //   1/4 dtau^2  | grad(w), grad(D(psi)) |
  //               |                       |

  LINALG::Matrix<nen_,nen_>  derxyMultderxy(true);
  derxyMultderxy.MultiplyTN(derxy_,derxy_);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      emat(fvi,fui) += derxyMultderxy(vi,ui)* (fac*pseudo_timestep_size*pseudo_timestep_size/4.0);
    }
  }


  //----------  --------------    |                       m         |
  //  rhs                    dtau | w, so(1.0- || grad(psi ) ||)    |
  //--------------------------    |                                 |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] += funct_(vi)*pseudo_timestep_size*fac*sign_phi_0*(1.0-grad_norm_dist_n);
  }


  //----------  --------------                |                 m   |
  //  rhs                    - 1.0/2.0*dtau^2 | grad(w),grad(psi )  |
  //--------------------------                |                     |

  LINALG::Matrix<nen_,1> derxyMultGradn(true);
  derxyMultGradn.MultiplyTN(derxy_,grad_dist_n);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= derxyMultGradn(vi)*pseudo_timestep_size*pseudo_timestep_size*fac/2.0;
  }


  // Assemble rhs for linear part of weak formulation for nonlinear iteration


  // get difference between psi^m+1 - psi^m
  double Delta_Psi = dist_npi - dist_n;

  LINALG::Matrix<nsd_,1> Delta_Grad_psi(true);
  Delta_Grad_psi.Update(1.0,grad_dist_npi, -1.0, grad_dist_n);


  //     |         m+1     m  |
  //     | -w, (psi   - psi ) |
  //     |         i          |

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= funct_(vi)*fac*Delta_Psi;
  }


  //                    |                   m+1     m  |
  //    1/4*delta_tau^2 | -grad(w), grad(psi   - psi ) |
  //                    |                   i          |

  LINALG::Matrix<nen_,1> Grad_w_Grad_Dpsi(true);
  Grad_w_Grad_Dpsi.MultiplyTN(derxy_,Delta_Grad_psi);


  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= Grad_w_Grad_Dpsi(vi,0)*fac*pseudo_timestep_size*pseudo_timestep_size/4.0;
  }

  return;
} //ReInitImpl::CalMatAndRHS_REINIT_ICG2




/*---------------------------------------------------------------------------------------------*
  |  calculate element matrix and rhs vector for the intersection penalty method   schott 12/10|
  *--------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalMatAndRHS_REINIT_Penalty(
    Epetra_SerialDenseMatrix&             emat,                     //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                     //!< element rhs to calculate
    const int                             k,                        //!< index of current scalar
    DRT::Element*                         ele,                      //!< the element we are dealing with
    const double                          penalty_interface_reinit  //!< penalty parameter to keep the interface position
  )
{

  //=======================================================================================================
  // get intersection point for each edge of the element
  // loop over edges
  //{
  // for each edge: check the G-function values, if there is an intersection point
  // if there is an intersection point: get the right local coordinates (interpolation of the local coordinates of the element)
  // evaluate the shape functions at this point
  // multiply the shape function with the phi-values
  // assemble G^T*G*phi = 0 into the sysmat
  //}
  //=======================================================================================================

  const size_t numnode = ele->NumNode();

  std::vector<Teuchos::RCP<DRT::Element> > linesVec = ele->Lines();
  typedef std::vector<Teuchos::RCP<DRT::Element> >::iterator lines_iterator;

  for(lines_iterator line = linesVec.begin(); line != linesVec.end(); line++)
  {
    if(distype==DRT::Element::hex8)
    {
      RCP<DRT::Element> pt_to_line = *line;
      const DRT::Node* const* vecOfPtsToNode = pt_to_line->Nodes();

      const int numberOfNodesOfLine = pt_to_line->NumNode();

      if (numberOfNodesOfLine != 2) dserror("not exact 2 nodes on this line");

      // get phi_value of current node
      int nodeID_start = vecOfPtsToNode[0]->Id();
      int nodeID_end   = vecOfPtsToNode[1]->Id();

      const int* ptToNodeIds_adj = ele->NodeIds();

      int ID_param_space_start = -1;
      int ID_param_space_end   = -1;

      for (size_t inode=0; inode < numnode; inode++)
      {
        // get local number of node actnode in ele_adj
        if(nodeID_start == ptToNodeIds_adj[inode]) ID_param_space_start = inode;
        if(nodeID_end   == ptToNodeIds_adj[inode]) ID_param_space_end   = inode;
      }

      if(ID_param_space_start == -1) dserror("node of line not a node of element!?!?!?");
      if(ID_param_space_end   == -1) dserror("node of line not a node of element!?!?!?");


      // get node xi coordinates
      LINALG::Matrix<3,1> node_Xicoordinates_start;
      node_Xicoordinates_start.Clear();
      node_Xicoordinates_start = DRT::UTILS::getNodeCoordinates(ID_param_space_start, distype);

      LINALG::Matrix<3,1> node_Xicoordinates_end;
      node_Xicoordinates_end.Clear();
      node_Xicoordinates_end = DRT::UTILS::getNodeCoordinates(ID_param_space_end, distype);


      // get the intersection point (linear interpolation)
      double interp_alpha = 0.0;
      double phi_end   = ephi0_penalty_[k](ID_param_space_end);
      double phi_start = ephi0_penalty_[k](ID_param_space_start);

      // if an intersection point is given
      if (phi_end * phi_start <= 0.0)
      {

        double phi_diff = phi_start-phi_end;
        if(fabs(phi_diff)< 1e-12)
        {
          // maybe a complete edge is zero -> do nothing for this element
          std::cout << "!!! WARNING: one element edge is zero in element " << ele->Id() << "-> check this penalty case!!! (do nothing at the moment)" << std::endl;
          return;
        }

        interp_alpha = -phi_end/(phi_diff);

        LINALG::Matrix<3,1> intersection_local(true);
        intersection_local.Update(interp_alpha,node_Xicoordinates_start, (1.0-interp_alpha), node_Xicoordinates_end);


        // evaluate the shape functions at the intersection point
        LINALG::Matrix<nen_,1> funct_intersection(true);
        DRT::UTILS::shape_function_3D(funct_intersection,intersection_local(0),intersection_local(1),intersection_local(2),distype);

        // assemble shape functions in sysmat
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) += funct_intersection(ui)*funct_intersection(vi)*penalty_interface_reinit;
          }
        }

        // assemble shape functions in rhs
        LINALG::Matrix<1,1> distnpi_intersection(true);

        distnpi_intersection.MultiplyTN(funct_intersection,ephinp_[k]);

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= penalty_interface_reinit*funct_intersection(vi,0)*distnpi_intersection(0,0);
        }

      } // end if intersection
    } // end of hex8
    else dserror("penalty not implemented for this type of element");
  } // end loop over lines

  return;
}




/*-------------------------------------------------------------------------------*
  |  evaluate element matrix and rhs for reinitialization (private)   schott 12/10|
  *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::CalMatAndRHS_REINIT_OST(
    Epetra_SerialDenseMatrix&             emat,                        //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                        //!< element rhs to calculate
    const double                          fac,                         //!< integration factor
    const int                             k,                           //!< index of current scalar
    DRT::Element*                         ele,                         //!< the element we are dealing with
    double                                pseudo_timestep_size_factor, //!< pseudo time step factor
    double                                meshsize,                    //!< meshsize
    const INPAR::SCATRA::PenaltyMethod    penalty_method,              //!< penalty method to keep the interface the interface position
    const double                          penalty_interface_reinit,    //!< penalty parameter to keep the interface position
    const double                          epsilon_bandwidth,           //!< epsilon bandwith for the smoothed sign function
    const INPAR::SCATRA::SmoothedSignType smoothedSignType,            //!< type of smoothing the sign function
    const bool                            shock_capturing,             //!< shock capturing switched on/off
    const double                          shock_capturing_diffusivity, //!< shock capturing diffusivity
    double                                timefac                      //!< time discretization factor
  )
{


  //----------------------------------------------------------------
  // 1) element matrix: stationary terms
  //----------------------------------------------------------------
  // stabilization parameter and integration factors
  const double taufac     = tau_[k]*fac;
  const double timefacfac = timefac*fac;
  const double timetaufac = timefac*taufac;
  const double fac_diffus = timefacfac*diffus_[k];



  //----------------------------------------------------------------
  // standard Galerkin terms
  //----------------------------------------------------------------
  // convective term in convective form
  const double densfac = timefacfac*densnp_[k];
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = densfac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*(conv_(ui)+sgconv_(ui));
    }
  }

  // addition to convective term for conservative form
  if (is_conservative_)
  {
    dserror("no conservative form for reinitialization equation implemented");
  }

  // diffusive term
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
      emat(fvi,fui) += fac_diffus*laplawf;
    }
  }


  if(shock_capturing)
  {
    const double fac_shock_capt = timefacfac*shock_capturing_diffusivity;

    // diffusive shock capturing term
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
        emat(fvi,fui) += fac_shock_capt*laplawf;
      }
    }


    // gradient of current scalar value
    gradphi_.Multiply(derxy_,ephinp_[k]);
    // diffusive shock capturing term
    double vrhs_shock_capt = fac_shock_capt;
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      double laplawf(0.0);
      GetLaplacianWeakFormRHS(laplawf,derxy_,gradphi_,vi);
      erhs[fvi] -= vrhs_shock_capt*laplawf;
    }
  }

  if(penalty_method == INPAR::SCATRA::penalty_method_akkerman)
  {

    // get derivative of smoothed sign_function
    //evaluate phinp and phi_0
    double phinp = funct_.Dot(ephinp_[k]);
    double phi_ref = funct_.Dot(ephi0_Reinit_Reference_[k]);

    double deriv_smoothed_Heavyside = EvaluateSmoothedSignDeriv(phi_ref, epsilon_bandwidth,meshsize, smoothedSignType);

    const double densfac_penalty = timefacfac*densnp_[k]*penalty_interface_reinit*deriv_smoothed_Heavyside;
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densfac_penalty*funct_(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*funct_(ui);
      }
    }

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      erhs[fvi] -= funct_(vi)*densfac_penalty*( phinp - phi_ref) ;
    }
  }

  //----------------------------------------------------------------
  // convective stabilization term
  //----------------------------------------------------------------
  // convective stabilization of convective term (in convective form)
  const double dens2taufac = timetaufac*densnp_[k]*densnp_[k];
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = dens2taufac*(conv_(vi)+sgconv_(vi));
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;

      emat(fvi,fui) += v*conv_(ui);
    }
  }

  //----------------------------------------------------------------
  // stabilization terms for higher-order elements
  //----------------------------------------------------------------
  if (use2ndderiv_)
  {
    // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
    GetLaplacianStrongForm(diff_, derxy2_);
    diff_.Scale(diffus_[k]);

    const double denstaufac = timetaufac*densnp_[k];
    // convective stabilization of diffusive term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = denstaufac*(conv_(vi)+sgconv_(vi));
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) -= v*diff_(ui);
      }
    }

    const double densdifftaufac = diffreastafac_*denstaufac;
    // diffusive stabilization of convective term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densdifftaufac*diff_(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) -= v*conv_(ui);
      }
    }

    const double difftaufac = diffreastafac_*timetaufac;
    // diffusive stabilization of diffusive term
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = difftaufac*diff_(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*diff_(ui);
      }
    }
  }

  //----------------------------------------------------------------
  // 2) element matrix: instationary terms
  //----------------------------------------------------------------
  if (not is_stationary_)
  {
    const double densamfac = fac*densam_[k];
    //----------------------------------------------------------------
    // standard Galerkin transient term
    //----------------------------------------------------------------
    // transient term
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densamfac*funct_(vi);
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*funct_(ui);
      }
    }

    const double densamnptaufac = taufac*densam_[k]*densnp_[k];
    //----------------------------------------------------------------
    // stabilization of transient term
    //----------------------------------------------------------------
    // convective stabilization of transient term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densamnptaufac*(conv_(vi)+sgconv_(vi));
      const int fvi = vi*numdofpernode_+k;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        emat(fvi,fui) += v*funct_(ui);
      }
    }

    if (use2ndderiv_)
    {
      const double densamreataufac = diffreastafac_*taufac*densam_[k];
      // diffusive stabilization of transient term
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = densamreataufac*diff_(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) -= v*funct_(ui);
        }
      }
    }
  }


  //----------------------------------------------------------------
  // 4) element right hand side
  //----------------------------------------------------------------
  //----------------------------------------------------------------
  // computation of bodyforce (and potentially history) term,
  // residual, integration factors and standard Galerkin transient
  // term (if required) on right hand side depending on respective
  // (non-)incremental stationary or time-integration scheme
  //----------------------------------------------------------------
  double rhsint    = rhs_[k];
  double residual  = 0.0;
  double rhsfac    = 0.0;
  double rhstaufac = 0.0;
  double conv_phi  = 0.0;
  double diff_phi  = 0.0;
  double rea_phi   = 0.0;
  if (is_incremental_ and is_genalpha_)
  {
    dserror("generalized alpha implementation not yet available");
  }
  else if (not is_incremental_ and is_genalpha_)
  {
    dserror("generalized alpha implementation not yet available");
  }
  else if (is_incremental_ and not is_genalpha_)
  {
    // gradient of current scalar value
    gradphi_.Multiply(derxy_,ephinp_[k]);

    // convective term using current scalar value
    conv_phi = convelint_.Dot(gradphi_);

    // diffusive term using current scalar value for higher-order elements
    if (use2ndderiv_) diff_phi = diff_.Dot(ephinp_[k]);

    // reactive term using current scalar value
    if (is_reactive_)
    {
      dserror("no reactive terms in reinitializaton equation");
    }

    if (not is_stationary_)
    {
      // compute scalar at integration point
      double dens_phi = funct_.Dot(ephinp_[k]);

      rhsint  *= timefac;
      rhsint  += densnp_[k]*hist_[k];
      residual = densnp_[k]*dens_phi + timefac*(densnp_[k]*conv_phi - diff_phi + rea_phi) - rhsint;
      rhsfac   = timefacfac;

      const double vtrans = fac*densnp_[k]*dens_phi;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] -= vtrans*funct_(vi);
      }
    }
    else
    {
      residual = densnp_[k]*conv_phi - diff_phi + rea_phi - rhsint;
      rhsfac   = fac;
    }
    rhstaufac = taufac;

    // addition to convective term due to subgrid-scale velocity
    // (not included in residual)
    double sgconv_phi = sgvelint_.Dot(gradphi_);
    conv_phi += sgconv_phi;

    // addition to convective term for conservative form
    // (not included in residual)
    if (is_conservative_)
    {
      dserror("no conservative form implemented for reinitialization equation");
    }

    // multiply convective term by density
    conv_phi *= densnp_[k];
  }
  else
  {
    if (not is_stationary_)
    {
      rhsint *= timefac;
      rhsint += densnp_[k]*hist_[k];
    }
    residual  = -rhsint;
    rhstaufac = taufac;
  }

  //----------------------------------------------------------------
  // standard Galerkin bodyforce term
  //----------------------------------------------------------------
  double vrhs = fac*rhsint;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] += vrhs*funct_(vi);
  }

  //----------------------------------------------------------------
  // standard Galerkin terms on right hand side
  //----------------------------------------------------------------
  // convective term
  vrhs = rhsfac*conv_phi;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= vrhs*funct_(vi);
  }

  // diffusive term
  vrhs = rhsfac*diffus_[k];
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf,derxy_,gradphi_,vi);
    erhs[fvi] -= vrhs*laplawf;
  }


  //----------------------------------------------------------------
  // stabilization terms
  //----------------------------------------------------------------
  // convective rhs stabilization (in convective form)
  vrhs = rhstaufac*residual*densnp_[k];
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= vrhs*(conv_(vi)+sgconv_(vi));
  }

  // diffusive rhs stabilization
  if (use2ndderiv_)
  {
    vrhs = rhstaufac*residual;
    // diffusive stabilization of convective temporal rhs term (in convective form)
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;
      erhs[fvi] += diffreastafac_*vrhs*diff_(vi);
    }
  }

  return;
} //ReInitImpl::CalMatAndRHS_LinearAdvection_REINITIALIZATION




/*----------------------------------------------------------------------*
  | calculate system matrix and rhs for the first order                 |
  | reinitialization equation based on a                                 |
  | implicit characteristic Galerkin scheme                  schott 12/10|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::Sysmat_Reinit_TG(
    DRT::Element*                         ele,                               //!< the element those matrix is calculated
    Epetra_SerialDenseMatrix&             emat,                              //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,                              //!< element rhs to calculate
    const bool                            reinitswitch,                      //!< flag if scatra object is a reinitializer
    const double                          reinit_pseudo_timestepsize_factor, //!< pseudo timestep factor
    const INPAR::SCATRA::SmoothedSignType smoothedSignType,                  //!< time of smoothed sign function
    const double                          epsilon_bandwidth,                 //!< epsilon bandwith
    const INPAR::SCATRA::PenaltyMethod    penalty_method,                    //!< penalty method to keep the interface position
    const double                          penalty_interface_reinit           //!< penalty parameter to keep the interface position
  )
{

  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,ele->Id());


  const double mesh_size = getEleDiameter<distype>(xyze_);

  const double pseudo_timestep_size = mesh_size * reinit_pseudo_timestepsize_factor;

  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_reinit(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  //==================================== new implementation of REINITIALIZATION================================
  // reinitialization according to sussman 1994, nagrath 2005                                     schott 12/10
  //===========================================================================================================


  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------

  // Assemble element rhs and vector for domain!!! integrals
  for (int iquad=0; iquad<intpoints_reinit.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPointREINIT(intpoints_reinit,iquad,ele->Id());

    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      CalMatAndRHS_REINIT_ICG2(emat,
          erhs,
          fac,
          k,
          ele,
          pseudo_timestep_size,
          mesh_size,
          penalty_method,
          penalty_interface_reinit,
          epsilon_bandwidth,
          smoothedSignType);

    } // loop over each scalar
  } // integration loop


  if(penalty_method == INPAR::SCATRA::penalty_method_intersection_points)
  {
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      CalMatAndRHS_REINIT_Penalty(emat,
          erhs,
          k,
          ele,
          penalty_interface_reinit);
    }
  }

  return;
}



/*----------------------------------------------------------------------*
  |  calculate system matrix and rhs (public)                schott 05/11|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ReInitImpl<distype>::Sysmat_Reinit_OST(
    DRT::Element*                         ele,           //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,          //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,          //!< element rhs to calculate
    const double                          dt,            //!< current time-step length
    const double                          timefac,       //!< time discretization factor
    const double                          meshsize,      //!< meshsize
    const bool                            reinitswitch,  //!< flag if scatra object is a reinitialization object
    const double                          reinit_pseudo_timestepsize_factor, //!< pseudo time step factor
    const INPAR::SCATRA::SmoothedSignType smoothedSignType,                  //!< type of smoothed sign function
    const INPAR::SCATRA::PenaltyMethod    penalty_method,                    //!< penalty method to keep the interface position
    const double                          penalty_interface_reinit,          //!< penalty parameter to keep the interface position
    const double                          epsilon_bandwidth,                 //!< epsilon bandwith
    const bool                            shock_capturing,                   //!< shockcapturing operator turned on/off
    const double                          shock_capturing_diffusivity,       //!< shockcapturing diffusivity
    const enum INPAR::SCATRA::ScaTraType  scatratype                         //!< type of scalar transport problem
  )
{
#define REINIT_LINEAR_ADVECTION_PHIGRADIENT
//#define REINIT_LINEAR_ADVECTION_RECONSTRUCTED_NORMALS
//#define DONT_EVALUATE_SMALL_GRADIENTS

  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,ele->Id());

  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------
  if (not mat_gp_ or not tau_gp_) GetMaterialParams(ele, scatratype);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  //=========================================================================================
  const double phi_gradient_TOL = 1e-003;
  bool do_evaluate = true;

  // decide if element is evaluated or not!!!
  // get phi_gradient at midpoint
  LINALG::Matrix<nsd_,1> midpoint(true); // midpoint is (0.0, 0.0, 0.0)
  DRT::UTILS::shape_function_deriv1<distype>(midpoint,deriv_);

  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  grad_phi_0.Clear();
  grad_phi_0.Multiply(derxy_,ephi0_Reinit_Reference_[0]);
  if(numscal_>1) dserror("evaluate check implemented for one scalar only");

  if(fabs(grad_phi_0.Norm2()) < phi_gradient_TOL)
  {
    do_evaluate=false;
    std::cout << "only mass matrix assembled in element " << ele->Id() << "Too small gradients for reinitialization" << std::endl;
  }

  //=========================================================================================

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(ele, scatratype);

    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // get gradients and norms
      grad_phi_0.Clear();
      grad_phi_0.Multiply(derxy_,ephi0_Reinit_Reference_[k]);

      // evaluate phi and gradnormphi
      const double phi_0    = funct_.Dot(ephi0_Reinit_Reference_[k]);

#ifdef REINIT_LINEAR_ADVECTION_PHIGRADIENT
      // use orginial phi-gradients for computation of reinit velocity

      double grad_norm_phi_0    = grad_phi_0.Norm2();
      if(fabs(grad_norm_phi_0)>1e-12) convelint_.Update(1.0/grad_norm_phi_0, grad_phi_0, 0.0);
      else                            convelint_.Clear();
#endif
#ifdef REINIT_LINEAR_ADVECTION_RECONSTRUCTED_NORMALS

      // use reconstructed phi-gradients for computation of reinit velocity
      // get velocity at integration point
      velint_.Multiply(evelnp_,funct_);
      convelint_.Multiply(econvelnp_,funct_);
#endif

      double smoothedSign = EvaluateSmoothedSign(phi_0, grad_norm_phi_0, epsilon_bandwidth,meshsize, smoothedSignType);


      if(do_evaluate==false)
      {
#ifdef DONT_EVALUATE_SMALL_GRADIENTS
      // do not advect in elements with small gradients, but assemble the mass matrix
      smoothedSign = 0.0;
#endif
      }


      // evaluate signum function and scale the normalized direction stored in velint_;
      velint_.Scale(smoothedSign);
      convelint_.Scale(smoothedSign);

      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      conv_.MultiplyTN(derxy_,convelint_);

      // velocity divergence required for conservative form
      if (is_conservative_) GetDivergence(vdiv_,evelnp_,derxy_);

      // ensure that subgrid-scale velocity and subgrid-scale convective part
      // are zero if not computed below
      sgvelint_.Clear();
      sgconv_.Clear();

      //--------------------------------------------------------------------
      // calculation of (fine-scale) subgrid diffusivity, subgrid-scale
      // velocity and stabilization parameter(s) at integration point
      //--------------------------------------------------------------------
      if (true) // tau_gp_
      {
        // calculation of stabilization parameter at integration point
                CalTau(ele,diffus_[k],dt,timefac,vol,k,0.0,false);
      }

      hist_[k] = funct_.Dot(ephin_[k]);

      // set the rhs_ for reinitialization problems
      rhs_[k] = densnp_[k]*smoothedSign;


      // compute matrix and rhs
      CalMatAndRHS_REINIT_OST(emat,
                              erhs,
                              fac,
                              k,
                              ele,
                              reinit_pseudo_timestepsize_factor,
                              meshsize,
                              penalty_method,
                              penalty_interface_reinit,
                              epsilon_bandwidth,
                              smoothedSignType,
                              shock_capturing,
                              shock_capturing_diffusivity,
                              timefac);
    } // loop over each scalar
  } // integration loop

  for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
  {
    if(penalty_method == INPAR::SCATRA::penalty_method_intersection_points)
    {
      CalMatAndRHS_REINIT_Penalty(emat,
                                  erhs,
                                  k,
                                  ele,
                                  penalty_interface_reinit);
    }
  }

  return;
}



/// here we tell the compiler which instances have to be built !!!
template class DRT::ELEMENTS::ReInitImpl<DRT::Element::hex8>;


