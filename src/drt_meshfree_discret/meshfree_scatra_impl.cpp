/*----------------------------------------------------------------------*/
/*!
  \file meshfree_scatra_impl.cpp

  \brief Internal implementation of meshfree scalar transport cells

  <pre>
  Maintainer: Keijo Nissen
  nissen@lnm.mw.tum.de
  http://www.lnm.mw.tum.de
  089 - 289-15253
  </pre>
*/
/*----------------------------------------------------------------------*/

#include "meshfree_scatra_impl.H"           // class declarations
#include "drt_meshfree_discret.H"              // for cast to get knots
#include "drt_meshfree_cell.H"              // for cast to get knots
#include "drt_meshfree_cell_utils.H"        // to get Gauss points in real space
#include "../drt_scatra/scatra_ele_action.H"// for enum of scatra actions
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H" // basis function evaluation
#include "../drt_mat/scatra_mat.H"          // in GetMaterialParams(): type ScatraMat
#include "../drt_lib/drt_globalproblem.H"   // in BodyForce(): DRT::Problem::Instance()
#include "../drt_lib/drt_utils.H"           // in Evaluate(): ExtractMyValues()
#include "../drt_lib/drt_condition_utils.H" // in BodyForce(): FindElementConditions()

/*==========================================================================*
 * class MeshfreeScaTraImplInterface                                        *
 *==========================================================================*/

/*--------------------------------------------------------------------------*
 |  internal implementation for meshfree scatra cells    (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
DRT::ELEMENTS::MeshfreeScaTraImplInterface* DRT::ELEMENTS::MeshfreeScaTraImplInterface::Impl(
  const DRT::Element* ele,
  const enum INPAR::SCATRA::ScaTraType& scatratype)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));

  switch (ele->Shape()){
  case DRT::Element::hex8 : {return MeshfreeScaTraImpl<DRT::Element::hex8 >::Instance(numdofpernode);}
  case DRT::Element::tet4 : {return MeshfreeScaTraImpl<DRT::Element::tet4 >::Instance(numdofpernode);}
  case DRT::Element::quad4: {return MeshfreeScaTraImpl<DRT::Element::quad4>::Instance(numdofpernode);}
  case DRT::Element::tri3 : {return MeshfreeScaTraImpl<DRT::Element::tri3 >::Instance(numdofpernode);}
  case DRT::Element::line2: {return MeshfreeScaTraImpl<DRT::Element::line2>::Instance(numdofpernode);}
  default: dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(ele->Shape()).c_str());}
  return NULL;
}

/*==========================================================================*
 * class MeshfreeScaTraImpl                                                 *
 *==========================================================================*/


/*--------------------------------------------------------------------------*
 |  static members                                      (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const int DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::nsd_;
//template <DRT::Element::DiscretizationType distype>
//const int DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::nek_;

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::MeshfreeScaTraImpl(const int numdofpernode, const int numscal)
: numdofpernode_(numdofpernode),
  numscal_(numscal),
  nen_(),
  phi_(),
  hist_(numdofpernode_),
  gradphi_(nsd_),
  ephin_(numscal_),
  ephinp_(numscal_),
  ephiam_(numscal_),
  ehist_(numdofpernode_),
  velint_(nsd_),
  convelint_(nsd_),
  vdiv_(0.0),
  evelnp_(),
  econvelnp_(),
  eaccnp_(),
  eprenp_(),
  edispnp_(true),
  funct_(),
  deriv_(),
  bodyforce_(numdofpernode_), // size of vector
  rhs_(numdofpernode_),
  reatemprhs_(numdofpernode_),
  scatrares_(numscal_),
  conv_phi_(numscal_),
  diff_phi_(numscal_),
  rea_phi_(numscal_),
  laplace_(true),
  efluxreconstr_(numscal_),
  densn_(numscal_),
  densnp_(numscal_),
  densam_(numscal_),
  densgradfac_(numscal_),
  diffus_(numscal_),
  reacoeff_(numscal_),
  reacoeffderiv_(numscal_),
  shc_(0.0),
  visc_(),
  conv_(),
  thermpressnp_(0.0),
  thermpressam_(0.0),
  thermpressdt_(0.0),
  is_elch_(numdofpernode-numscal==1),
  is_reactive_(false) // flag
{
  return;
}


/*--------------------------------------------------------------------------*
 |  singleton access method                              (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeScaTraImpl<distype> * DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create
  )
{
  static MeshfreeScaTraImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new MeshfreeScaTraImpl<distype>(numdofpernode,numscal);
    }
  }
  else // delete
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*--------------------------------------------------------------------------*
 |  called upon singleton destruction                    (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0,0,false);
}

/*--------------------------------------------------------------------------*
 |  evaluate meshfree scatra cell                        (public) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::Evaluate(
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

  // get number of nodes
  nen_ = ele->NumNode();

  // get discretization
  discret_ = dynamic_cast<DRT::MESHFREE::MeshfreeDiscretization*>(&(discretization));

  // set size of all vectors of SerialDense element arrays
  funct_.LightSize(nen_);
  deriv_.LightShape(nsd_,nen_);
  evelnp_.LightShape(nsd_,nen_);
  econvelnp_.LightShape(nsd_,nen_);
  conv_.LightSize(nen_);
  for (int k=0;k<numscal_;++k){
    // set size of all vectors of SerialDenseVectors
    ephin_[k].LightSize(nen_); // without initialisation
    ephinp_[k].LightSize(nen_); // without initialisation
    ehist_[k].LightSize(nen_); // without initialisation
    bodyforce_[k].LightSize(nen_);
    // set size of all vectors of SerialDenseMatrices
  }
  bodyforce_[numdofpernode_-1].LightSize(nen_);


  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");

  // set parameters for stabilization
  Teuchos::ParameterList& stablist = params.sublist("STABILIZATION");

  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Set parameter SCATRATYPE in your input file!");

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch (action)
  {
  case SCATRA::calc_mat_and_rhs:
  {
    // get control parameters
    is_stationary_  = params.get<bool>("using stationary formulation");
    is_genalpha_    = params.get<bool>("using generalized-alpha time integration");
    is_incremental_ = params.get<bool>("incremental solver");

    // get current time and time-step length
    const double time = params.get<double>("total time");
    const double dt   = params.get<double>("time-step length");

    // get time factor and alpha_F if required
    // one-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    if (not is_stationary_)
    {
      timefac = params.get<double>("time factor");
      if (is_genalpha_)
      {
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
      }
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // set flag for conservative form
    is_conservative_ = false;
    const INPAR::SCATRA::ConvForm convform =
      DRT::INPUT::get<INPAR::SCATRA::ConvForm>(params, "form of convective term");
    if (convform ==INPAR::SCATRA::convform_conservative) is_conservative_ = true;

    // set flag for material law and locality at int. point ? why name stablist ?
    const INPAR::SCATRA::EvalMat matloc = DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(stablist,"EVALUATION_MAT");
    mat_gp_ = (matloc == INPAR::SCATRA::evalmat_integration_point); // set true/false

    // get velocity at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evelnp_,velocity,nsd_);
    const RCP<Epetra_MultiVector> convelocity = params.get< RCP<Epetra_MultiVector> >("convective velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,econvelnp_,convelocity,nsd_);

    // extract local values from the global vectors
    RCP<const Epetra_Vector> hist = discretization.GetState("hist");
    RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (hist==Teuchos::null || phinp==Teuchos::null)
      dserror("Cannot get state vector 'hist' and/or 'phinp'");
    std::vector<double> myhist(lm.size());
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // fill all element arrays
    for (int i=0;i<nen_;++i)
    {
      int temp = i*numdofpernode_;
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephinp_[k][i] = myphinp[k+temp];
        // the history vectors contains information of time step t_n
        ehist_[k][i] = myhist[k+temp];
      }
      if (is_elch_)
      {
        // the history vectors contains information of time step t_n
        ehist_[numdofpernode_][i] = myhist[numdofpernode_+temp];
      }
    } // for i

    if (is_genalpha_ and not is_incremental_)
    {
      // extract additional local values from global vector
      RCP<const Epetra_Vector> phin = discretization.GetState("phin");
      if (phin==Teuchos::null) dserror("Cannot get state vector 'phin'");
      std::vector<double> myphin(lm.size());
      DRT::UTILS::ExtractMyValues(*phin,myphin,lm);

      // fill element array
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephin_[k][i] = myphin[k+(i*numdofpernode_)];
        }
      } // for i
    }

    // calculate element coefficient matrix and rhs
    Sysmat(
      ele,
      elemat1_epetra,
      elevec1_epetra,
      time,
      dt,
      timefac,
      alphaF,
      scatratype);
#if 0
    // for debugging of matrix entries
    if(ele->Id()==2) // and (time < 3 or time > 99.0))
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
        scatratype);
    }
#endif
    break;
  }
  case SCATRA::calc_domain_and_bodyforce:
  {
    // NOTE: add integral values only for elements which are NOT ghosted!

    if (ele->Owner() == discretization.Comm().MyPID())
    {
      const double time = params.get<double>("total time");

      // calculate domain and bodyforce integral
      CalculateDomainAndBodyforce(elevec1_epetra,ele,time);
    }

    break;
  }
  case SCATRA::integrate_shape_functions:
  {
    // calculate integral of shape functions
    const Epetra_IntSerialDenseVector dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
    IntegrateShapeFunctions(ele,elevec1_epetra,dofids);

    break;
  }
  case SCATRA::get_material_parameters:
  {
    // get the material
    RCP<MAT::Material> material = ele->Material();
    params.set("thermodynamic pressure",0.0);

    break;
  }
  default:
  {
    dserror("Not acting on action No. %i. Forgot implementation?",action);
  }
  }
  // work is done
  return 0;
}

/*--------------------------------------------------------------------------*
 |  initiates calculation of system matrix and rhs      (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::Sysmat(
  DRT::Element*                         ele,
  Epetra_SerialDenseMatrix&             sys_mat,
  Epetra_SerialDenseVector&             residual,
  const double &                        time,
  const double &                        dt,
  const double &                        timefac,
  const double &                        alphaF,
  const enum INPAR::SCATRA::ScaTraType& scatratype
  )
{
  //----------------------------------------------------------------------
  // cast element pointer to cell pointer to get access to knot information
  //----------------------------------------------------------------------
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const *>(ele);
  if (cell==NULL)
    dserror("dynamic_cast of element to meshfree cell failed!");

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele,time);

  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------
  if (not mat_gp_) GetMaterialParams(cell,scatratype);

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);

  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the meshfree basis functions.");

    // DEBUG
//    std::cout << "funct_ = ";
//    for (int i=0; i<nen_; i++)
//      std::cout << funct_(i) << " ";
//    std::cout << std::endl;
//
//    std::cout << "deriv_ = ";
//    for (int i=0; i<nsd_; i++){
//      for (int j=0; j<nen_; j++)
//        std::cout << deriv_(i,j) << " ";
//      std::cout << std::endl;
//    }
    // END DEBUG

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(cell,scatratype);

    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // get velocity at integration point
      velint_.Multiply('N','N',1.0,evelnp_,funct_,0.0);
      convelint_.Multiply('N','N',1.0,econvelnp_,funct_,0.0);

      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      conv_.Multiply('T','N',1.0,deriv_,convelint_,0.0);

      // gradient of current scalar value
      gradphi_.Multiply('N','N',1.0,deriv_,ephinp_[k],0.0);

      // velocity divergence required for conservative form
      if (is_conservative_) GetDivergence(vdiv_,evelnp_,deriv_);

      // get history data (or acceleration)
      hist_[k] = funct_.Dot(ehist_[k]);

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term

      // TODO: this only works for non-negative basis functions
      //       create:         LINALG::SerialDenseVector::Sum()
      //       or even better: LINALG::SerialDenseVector::Sum(double)
      rhs_[k] = bodyforce_[k].Dot(funct_)/shc_; // normal bodyforce
      rhs_[k] += thermpressdt_/shc_;            // time derivative of thermodynamic pressure (loma)
      rhs_[k] += densnp_[k]*reatemprhs_[k];     // temperature reaction-rate term

      CalMatAndRHS(sys_mat,residual,fac,timefac,dt,alphaF,k);
    } // loop over each scalar
  }// integration loop

  return;


} //MeshfreeScaTraImpl::Sysmat


/*--------------------------------------------------------------------------*
 |  get the body force                                 (private)  nis Feb12 |
 |                                                                          |
 |  this funtions needs to set bodyforce_ every time even for static or     |
 |  even zero forces, since it is not an element, but an impl bodyforce     |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::BodyForce(
  DRT::Element const * const & ele,
  const double&                time
  )
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique Neumann condition
  switch(nsd_){
  case 3: DRT::UTILS::FindElementConditions(ele, "VolumeNeumann" , myneumcond); break;
  case 2: DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond); break;
  case 1: DRT::UTILS::FindElementConditions(ele, "LineNeumann"   , myneumcond); break;
  default: dserror("Illegal number of spatial dimensions: %d",nsd_);
  }

  if (myneumcond.size()>1)
    dserror("More than one Neumann condition on one node!");

  if (myneumcond.size()==1)
  {
    // check for potential time curve
    const std::vector<int>* curve  = myneumcond[0]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];

    // initialization of time-curve factor
    double curvefac(0.0);

    // compute potential time curve or set time-curve factor to one
    if (curvenum >= 0)
    {
      // time factor (negative time indicating error)
      if (time >= 0.0)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      else dserror("Negative time in bodyforce calculation: time = %f",time);
    }
    else curvefac = 1.0;

    // get values and switches from the condition
    const std::vector<int>*    onoff = myneumcond[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = myneumcond[0]->Get<std::vector<double> >("val"  );

    // set this condition to the bodyforce array
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      for (int jnode=0; jnode<nen_; jnode++)
      {
        (bodyforce_[idof])(jnode) = (*onoff)[idof]*(*val)[idof]*curvefac;
      }
    }
  }
  else
  {
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      // no bodyforce - set all entries to zero
      bodyforce_[idof].Zero();
    }
  }

  return;
} //MeshfreeScaTraImpl::BodyForce


/*--------------------------------------------------------------------------*
 |  get the material constants                         (private)  nis Feb12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::GetMaterialParams(
  const DRT::Element*  ele,
  const enum INPAR::SCATRA::ScaTraType& scatratype
  )
{
// get the material
  RCP<MAT::Material> material = ele->Material();

// get diffusivity / diffusivities
  if (material->MaterialType() == INPAR::MAT::m_scatra)
  {
    const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

    dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

    // get constant diffusivity
    diffus_[0] = actmat->Diffusivity();

    // in case of reaction with (non-zero) constant coefficient:
    // read coefficient and set reaction flag to true
    reacoeff_[0] = actmat->ReaCoeff();
    if (reacoeff_[0] >  1.e-14) is_reactive_ = true;
    if (reacoeff_[0] < -1.e-14)
      dserror("Reaction coefficient for species %d is not positive: %f",0, reacoeff_[0]);

    // set specific heat capacity at constant pressure to 1.0
    shc_ = 1.0;

    // set density at various time steps and density gradient factor to 1.0/0.0
    densn_[0]       = 1.0;
    densnp_[0]      = 1.0;
    densam_[0]      = 1.0;
    densgradfac_[0] = 0.0;

    reacoeffderiv_[0] = reacoeff_[0];
  }
  else dserror("Material type is not supported");

// check whether there is negative (physical) diffusivity
  if (diffus_[0] < -1.e-15) dserror("negative (physical) diffusivity");

  return;
} //MeshfreeScaTraImpl::GetMaterialParams

/*----------------------------------------------------------------------*
  |  evaluate element matrix and rhs (private)                   vg 02/09|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::CalMatAndRHS(
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const double&                         fac,
  const double&                         timefac,
  const double&                         dt,
  const double&                         alphaF,
  const int&                            dofindex
  )
{
  //----------------------------------------------------------------
  // 1) element matrix: stationary terms
  //----------------------------------------------------------------

  // integration factors
  const double timefacfac = timefac*fac;
  const double fac_diffus = timefacfac*diffus_[dofindex];

//  std::cout << "fac = " << fac << "; timefac = " << timefac << "; timefacfac = " << timefacfac << "; fac_diffus = " << fac_diffus << std::endl;

  //----------------------------------------------------------------
  // standard Galerkin terms
  //----------------------------------------------------------------

  // convective term in convective form
  const double densfac = timefacfac*densnp_[dofindex];
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = densfac*funct_(vi);
    const int fvi = vi*numdofpernode_+dofindex;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+dofindex;

      emat(fvi,fui) += v*conv_(ui);
    }
  }

  // addition to convective term for conservative form
  if (is_conservative_)
  {
    // gradient of current scalar value
    deriv_.Multiply(false,ephinp_[dofindex],gradphi_); // gradphi_ = deriv_*ephinp_[dofindex]

    // convective term using current scalar value
    const double cons_conv_phi = convelint_.Dot(gradphi_);

    const double consfac = timefacfac*(densnp_[dofindex]*vdiv_+densgradfac_[dofindex]*cons_conv_phi);
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = consfac*funct_(vi);
      const int fvi = vi*numdofpernode_+dofindex;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+dofindex;

        emat(fvi,fui) += v*funct_(ui);
      }
    }
  }

  // diffusive term
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+dofindex;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+dofindex;
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf, deriv_,ui,vi);
      emat(fvi,fui) += fac_diffus*laplawf;
    }
  }

  //----------------------------------------------------------------
  // 2) element matrix: instationary terms
  //----------------------------------------------------------------
  if (not is_stationary_)
  {
    const double densamfac = fac*densam_[dofindex];
    //----------------------------------------------------------------
    // standard Galerkin transient term
    //----------------------------------------------------------------
    // transient term
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = densamfac*funct_(vi);
      const int fvi = vi*numdofpernode_+dofindex;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+dofindex;

        emat(fvi,fui) += v*funct_(ui);
      }
    }

  }

  //----------------------------------------------------------------
  // 3) element matrix: reactive terms
  //----------------------------------------------------------------
  if (is_reactive_)
  {
    const double fac_reac        = timefacfac*densnp_[dofindex]*reacoeffderiv_[dofindex];
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = fac_reac*funct_(vi);
      const int fvi = vi*numdofpernode_+dofindex;

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+dofindex;

        emat(fvi,fui) += v*funct_(ui);
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
  double rhsint    = rhs_[dofindex];
  double rhsfac    = 0.0;

  if (is_incremental_ and is_genalpha_)
  {
    // gradient of current scalar value
    deriv_.Multiply(false,ephinp_[dofindex],gradphi_); // gradphi_=deriv_*ephinp_[dofindex]

    rhsfac    = timefacfac/alphaF;
    rhsint   *= (timefac/alphaF);

    const double vtrans = rhsfac*densam_[dofindex]*hist_[dofindex];
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+dofindex;

      erhs[fvi] -= vtrans*funct_(vi);
    }

    // addition to convective term for conservative form
    // (not included in residual)
    if (is_conservative_)
    {
      // scalar at integration point at time step n
      phi_ = funct_.Dot(ephinp_[dofindex]);

      // convective term in conservative form
      conv_phi_[dofindex] += phi_*(vdiv_+(densgradfac_[dofindex]/densnp_[dofindex])*conv_phi_[dofindex]);
    }

    // multiply convective term by density
    conv_phi_[dofindex] *= densnp_[dofindex];
  }
  else if (not is_incremental_ and is_genalpha_)
  {
   // gradient of scalar value at n
    deriv_.Multiply(false,ephin_[dofindex],gradphi_); // gradphi_=deriv_*ephin_[dofindex]

    // convective term using scalar value at n
    conv_phi_[dofindex] = convelint_.Dot(gradphi_);

    // reactive term using scalar value at n
    if (is_reactive_)
    {
      // scalar at integration point
      phi_ = funct_.Dot(ephin_[dofindex]);

      rea_phi_[dofindex] = densnp_[dofindex]*reacoeff_[dofindex]*phi_;
    }

    rhsint   += densam_[dofindex]*hist_[dofindex]*(alphaF/timefac);
    scatrares_[dofindex] = (1.0-alphaF) * (densn_[dofindex]*conv_phi_[dofindex]
                                           - diff_phi_[dofindex] + rea_phi_[dofindex]) - rhsint;
    rhsfac    = timefacfac*(1.0-alphaF)/alphaF;
    rhsint   *= (timefac/alphaF);

    // addition to convective term for conservative form
    // (not included in residual)
    if (is_conservative_)
    {
      // scalar at integration point at time step n
      phi_ = funct_.Dot(ephin_[dofindex]);

      // convective term in conservative form
      // caution: velocity divergence is for n+1 and not for n!
      // -> hopefully, this inconsistency is of small amount
      conv_phi_[dofindex] += phi_*(vdiv_+(densgradfac_[dofindex]/densn_[dofindex])*conv_phi_[dofindex]);
    }

    // multiply convective term by density
    conv_phi_[dofindex] *= densn_[dofindex];
  }
  else if (is_incremental_ and not is_genalpha_)
  {
    // gradient of current scalar value
    gradphi_.Multiply('N','N',1.0,deriv_,ephinp_[dofindex],0.0); // gradphi_=deriv_*ephinp_[dofindex]

    if (not is_stationary_)
    {
      scatrares_[dofindex] *= dt;
      rhsint               *= timefac;
      rhsint               += densnp_[dofindex]*hist_[dofindex];
      rhsfac                = timefacfac;

      // compute scalar at integration point
      phi_ = funct_.Dot(ephinp_[dofindex]);

      const double vtrans = fac*densnp_[dofindex]*phi_;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        erhs[fvi] -= vtrans*funct_(vi);
      }
    }
    else rhsfac   = fac;

    // addition to convective term for conservative form
    // (not included in residual)
    if (is_conservative_)
    {
      // scalar at integration point at time step n
      phi_ = funct_.Dot(ephinp_[dofindex]);

      // convective term in conservative form
      conv_phi_[dofindex] += phi_*(vdiv_+(densgradfac_[dofindex]/densnp_[dofindex])*conv_phi_[dofindex]);
    }

    // multiply convective term by density
    conv_phi_[dofindex] *= densnp_[dofindex];
  }
  else
  {
    if (not is_stationary_)
    {
      rhsint *= timefac;
      rhsint += densnp_[dofindex]*hist_[dofindex];
    }
    scatrares_[dofindex] = -rhsint;
  }

  //----------------------------------------------------------------
  // standard Galerkin bodyforce term
  //----------------------------------------------------------------
  double vrhs = fac*rhsint;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+dofindex;

    erhs[fvi] += vrhs*funct_(vi);
  }

  //----------------------------------------------------------------
  // standard Galerkin terms on right hand side
  //----------------------------------------------------------------

  // convective term
  vrhs = rhsfac*conv_phi_[dofindex];
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+dofindex;

    erhs[fvi] -= vrhs*funct_(vi);
  }

  // diffusive term
  vrhs = rhsfac*diffus_[dofindex];
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+dofindex;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf,deriv_,gradphi_,vi);
    erhs[fvi] -= vrhs*laplawf;
  }

  //----------------------------------------------------------------
  // reactive terms (standard Galerkin) on rhs
  //----------------------------------------------------------------

  // standard Galerkin term
  if (is_reactive_)
  {
    vrhs = rhsfac*rea_phi_[dofindex];
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+dofindex;

      erhs[fvi] -= vrhs*funct_(vi);
    }
  }

  return;
} //MeshfreeScaTraImpl::CalMatAndRHS

/*--------------------------------------------------------------------------*
 |  calc mass matrix + rhs for initial time derivative  (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::InitialTimeDerivative(
  DRT::Element*                         ele,
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const enum INPAR::SCATRA::ScaTraType& scatratype
  )
{
  // dead load in element nodes at initial point in time
  const double time = 0.0;

  //----------------------------------------------------------------------
  // cast element pointer to cell pointer to get access to knot information
  //----------------------------------------------------------------------
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const * >(ele);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele,time);

  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------
  if (not mat_gp_) GetMaterialParams(cell,scatratype);

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight (here called fac instead of cgw)

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=1; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(ele,scatratype);

    //------------ get values of variables at integration point
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // get bodyforce in gausspoint (divided by specific heat capacity)
      // (For temperature equation, time derivative of thermodynamic pressure
      //  is added, if not constant.)
      rhs_[k] = bodyforce_[k].Dot(funct_)/shc_; // normal bodyforce
      rhs_[k] += thermpressdt_/shc_;            // time derivative of thermodynamic pressure (loma)

      // get velocity at integration point
      evelnp_.Multiply(false,funct_,velint_);       // velint_ = evelnp_ * funct_
      econvelnp_.Multiply(false,funct_,convelint_); // convelint_ = econvelnp_ * funct_

      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      deriv_.Multiply(false,convelint_,conv_);      // conv_ = deriv_ * convelint_;

      // velocity divergence required for conservative form
      if (is_conservative_) GetDivergence(vdiv_,evelnp_,deriv_);

      // diffusive integration factor
      const double fac_diffus = fac*diffus_[k];

      // get value of current scalar
      phi_ = funct_.Dot(ephinp_[k]);

      // gradient of current scalar value
      deriv_.Multiply(false,ephinp_[k],gradphi_);    // gradphi_ = deriv_ * ephip_;

      // convective part in convective form times initial scalar field
      double conv_ephi0_k = conv_.Dot(ephinp_[k]);

      // addition to convective term for conservative form
      // -> spatial variation of density not yet accounted for
      if (is_conservative_)
        conv_ephi0_k += phi_*(vdiv_+(densgradfac_[k]/densnp_[k])*conv_ephi0_k);

      //----------------------------------------------------------------
      // element matrix: transient term
      //----------------------------------------------------------------
      // transient term
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = fac*funct_(vi)*densnp_[k];
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) += v*funct_(ui);
        }
      }

      //----------------------------------------------------------------
      // element right hand side: convective term in convective form
      //----------------------------------------------------------------
      double vrhs = fac*densnp_[k]*conv_ephi0_k;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] -= vrhs*funct_(vi);
      }

      //----------------------------------------------------------------
      // element right hand side: diffusive term
      //----------------------------------------------------------------
      vrhs = fac_diffus;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        double laplawf(0.0);
        GetLaplacianWeakFormRHS(laplawf,deriv_,gradphi_,vi);
        erhs[fvi] -= vrhs*laplawf;
      }

      //----------------------------------------------------------------
      // element right hand side: reactive term
      //----------------------------------------------------------------
      if (is_reactive_)
      {
        vrhs = fac*densnp_[k]*reacoeff_[k]*phi_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vrhs*funct_(vi);
        }
      }

      //----------------------------------------------------------------
      // element right hand side: bodyforce term
      //----------------------------------------------------------------
      vrhs = fac*rhs_[k];
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] += vrhs*funct_(vi);
      }
    } // loop over each scalar k

  } // integration loop

  return;
} // MeshfreeScaTraImpl::InitialTimeDerivative

/*--------------------------------------------------------------------------*
 |  calculate mass matrix + rhs for time deriv. reinit. (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::TimeDerivativeReinit(
  DRT::Element*                         ele,
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const double&                         dt,
  const double&                         timefac,
  const enum INPAR::SCATRA::ScaTraType& scatratype
  )
{
  //----------------------------------------------------------------------
  // cast element pointer to cell pointer to get access to knot information
  //----------------------------------------------------------------------
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const *>(ele);

  //------------------------------------------------------------------------------------
  // get material parameters and stabilization parameters (evaluation at element center)
  //------------------------------------------------------------------------------------
  if (not mat_gp_) GetMaterialParams(ele,scatratype);

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight (here not called cgw)

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(ele,scatratype);

    //------------ get values of variables at integration point
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      // get bodyforce in gausspoint (divided by specific heat capacity)
      // (For temperature equation, time derivative of thermodynamic pressure
      //  is added, if not constant.)
      rhs_[k] = bodyforce_[k].Dot(funct_) / shc_;

      // get velocity at integration point
      evelnp_.Multiply(false,funct_,velint_);       // velint_ = evelnp_ * funct_
      econvelnp_.Multiply(false,funct_,convelint_); // convelint_ = econvelnp_ * funct_

      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      deriv_.Multiply(false,convelint_,conv_);      // conv_ = deriv_ * convelint_;

      // velocity divergence required for conservative form
      if (is_conservative_) GetDivergence(vdiv_,evelnp_,deriv_);

      // diffusive integration factor
      const double fac_diffus = fac*diffus_[k];

      // get value of current scalar
      phi_ = funct_.Dot(ephinp_[k]);

      // gradient of current scalar value
      deriv_.Multiply(false,ephinp_[k],gradphi_);    // gradphi_ = deriv_ * ephinp_[dofindex]

      // convective part in convective form times initial scalar field
      double conv_ephi0_k = conv_.Dot(ephinp_[k]);

      // addition to convective term for conservative form
      // -> spatial variation of density not yet accounted for
      if (is_conservative_)
        conv_ephi0_k += phi_*(vdiv_+(densgradfac_[k]/densnp_[k])*conv_ephi0_k);

      //----------------------------------------------------------------
      // element matrix: transient term
      //----------------------------------------------------------------
      // transient term
      for (int vi=0; vi<nen_; ++vi)
      {
        const double v = fac*funct_(vi)*densnp_[k];
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) += v*funct_(ui);
        }
      }

      //----------------------------------------------------------------
      // element right hand side: convective term in convective form
      //----------------------------------------------------------------
      double vrhs = fac*densnp_[k]*conv_ephi0_k;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] -= vrhs*funct_(vi);
      }

      //----------------------------------------------------------------
      // element right hand side: diffusive term
      //----------------------------------------------------------------
      vrhs = fac_diffus;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        double laplawf(0.0);
        GetLaplacianWeakFormRHS(laplawf,deriv_,gradphi_,vi);
        erhs[fvi] -= vrhs*laplawf;
      }

      //----------------------------------------------------------------
      // element right hand side: reactive term
      //----------------------------------------------------------------
      if (is_reactive_)
      {
        vrhs = fac*densnp_[k]*reacoeff_[k]*phi_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vrhs*funct_(vi);
        }
      }

      //----------------------------------------------------------------
      // element right hand side: bodyforce term
      //----------------------------------------------------------------
      vrhs = fac*rhs_[k];
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] += vrhs*funct_(vi);
      }
    } // loop over each scalar k

  } // integration loop

  return;
} // MeshfreeScaTraImpl::TimeDerivativeReinit

/*--------------------------------------------------------------------------*
 |  calc weighted mass flux (no reactive flux so far)   (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::CalculateFlux(
  LINALG::SerialDenseMatrix &     flux,
  const DRT::Element*             ele,
  const INPAR::SCATRA::FluxType&  fluxtype,
  const int&                      dofindex,
  const enum INPAR::SCATRA::ScaTraType& scatratype
  )
{
  /*------------------------------------------------------------------------
  //
  // Actually, we compute here a weighted (and integrated) form of the fluxes!
  // On time integration level, these contributions are then used to calculate
  // an L2-projected representation of fluxes.
  // Thus, this method here DOES NOT YET provide flux values that are ready to
  // use!!
  //
  //   /                                \
  //   |                /   \           |
  //   | w, -D * nabla | phi | + u*phi  |
  //   |                \   /           |
  //   \                      [optional]/
  //
  //------------------------------------------------------------------------*/

  // cast element pointer to cell pointer to get access to knot information
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const * >(ele);

  // get material parameters (evaluation at element center)
  if (not mat_gp_) GetMaterialParams(ele,scatratype);

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

    // get material parameters (evaluation at integration point)
    if (mat_gp_) GetMaterialParams(ele,scatratype);

    // get velocity at integration point
    evelnp_.Multiply(false,funct_,velint_);       // velint_ = evelnp_ * funct_
    econvelnp_.Multiply(false,funct_,convelint_); // convelint_ = econvelnp_ * funct_

    // get scalar at integration point
    phi_ = funct_.Dot(ephinp_[dofindex]);

    // gradient of current scalar value
    deriv_.Multiply(false,ephinp_[dofindex],gradphi_);    // gradphi_ = deriv_ * ephinp_[dofindex];

    // allocate and initialize!
    LINALG::SerialDenseVector q(nsd_);

    // add different flux contributions as specified by user input
    switch (fluxtype)
    {
    case INPAR::SCATRA::flux_total_domain:

      // convective flux contribution
      q = convelint_;
      q.Scale(densnp_[dofindex]*phi_);

      // no break statement here!
    case INPAR::SCATRA::flux_diffusive_domain:
      // diffusive flux contribution
      q.Update(-diffus_[dofindex],gradphi_,1.0);

      break;
    default:
      dserror("received illegal flag inside flux evaluation for whole domain");
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
} // MeshfreeScaTraImpl::CalculateFlux

/*--------------------------------------------------------------------------*
 |  calculate scalar(s) and domain integral             (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::CalculateScalars(
  const DRT::Element*             ele,
  const std::vector<double>&       ephinp,
  Epetra_SerialDenseVector&       scalars,
  const bool                      inverting
  )
{
  //----------------------------------------------------------------------
  // cast element pointer to cell pointer to get access to knot information
  //----------------------------------------------------------------------
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const * >(ele);

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

    // calculate integrals of (inverted) scalar(s) and domain
    if (inverting)
    {
      for (int i=0; i<nen_; i++)
      {
        const double fac_funct_i = fac*funct_(i);
        for (int k = 0; k < numscal_; k++)
        {
          scalars[k] += fac_funct_i/ephinp[i*numdofpernode_+k];
        }
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
          scalars[k] += fac_funct_i*ephinp[i*numdofpernode_+k];
        }
        scalars[numscal_] += fac_funct_i;
      }
    }
  } // loop over integration points

  return;
} // ScaTraImpl::CalculateScalars


/*--------------------------------------------------------------------------*
 |  calculate domain integral                           (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::CalculateDomainAndBodyforce(
  Epetra_SerialDenseVector&  scalars,
  const DRT::Element*        ele,
  const double&              time
  )
{
  //----------------------------------------------------------------------
  // cast element pointer to cell pointer to get access to knot information
  //----------------------------------------------------------------------
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const * >(ele);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------

  BodyForce(ele,time);

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

    // get bodyforce in gausspoint
    rhs_[0] = bodyforce_[0].Dot(funct_);

    // calculate integrals of domain and bodyforce
    for (int i=0; i<nen_; i++)
    {
      scalars[0] += fac*funct_(i);
    }
    scalars[1] += fac*rhs_[0];

  } // loop over integration points

  return;
} // ScaTraImpl::CalculateDomain


/*--------------------------------------------------------------------------*
 |  Integrate shape functions over domain               (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::IntegrateShapeFunctions(
  DRT::Element const *               ele,
  Epetra_SerialDenseVector&          elevec1,
  const Epetra_IntSerialDenseVector& dofids
  )
{
  //----------------------------------------------------------------------
  // cast element pointer to cell pointer to get access to knot information
  //----------------------------------------------------------------------
  DRT::MESHFREE::Cell const * cell = dynamic_cast<DRT::MESHFREE::Cell const * >(ele);

  // safety check
  if (dofids.M() < numdofpernode_)
    dserror("Dofids vector is too short. Received not enough flags");

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  LINALG::SerialDenseMatrix gxyz; // read: Gauss xyz-coordinate
  LINALG::SerialDenseVector gw;   // read: Gauss xyz-coordinate
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(cell->Knots(), gxyz, gw);
  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  DRT::Node const * const * const nodes = cell->Nodes(); // node pointer
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate
  double fac;     // current Gauss weight

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz[iquad]; // read: current gauss xyz-coordinate
    fac = gw[iquad];     // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nodes[i]->X();
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int error = discret_->solutionfunct_->GetMeshfreeBasisFunction(funct_,deriv_,distng,nsd_);
    if (error) dserror("Something went wrong when calculating the max-ent basis functions.");

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

} //MeshfreecaTraImpl<distype>::IntegrateShapeFunction

/*--------------------------------------------------------------------------*
 |  Do a finite difference check for a given element id.                    |
 |  Meant for debugging only!                           (private) nis Mar12 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeScaTraImpl<distype>::FDcheck(
  DRT::Element*                         ele,
  Epetra_SerialDenseMatrix&             sys_mat,
  Epetra_SerialDenseVector&             residual,
  const double&                         time,
  const double&                         dt,
  const double&                         timefac,
  const double&                         alphaF,
  const enum INPAR::SCATRA::ScaTraType& scatratype
  )
{
  // magnitude of dof perturbation
  const double epsilon=1e-6; // 1.e-8 seems already too small!

  // make a copy of all input parameters potentially modified by Sysmat
  // call --- they are not intended to be modified

  // alloc the vectors that will store the original, non-perturbed values
  std::vector<LINALG::SerialDenseVector> origephinp(numscal_);
  std::vector<LINALG::SerialDenseVector> origehist(numscal_);

  // copy original concentrations and potentials to these storage arrays
  for (int i=0;i<nen_;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      origephinp[k](i) = ephinp_[k](i);
      origehist[k](i)  = ehist_[k](i);
    }
  } // for i

  // allocate arrays to compute element matrices and vectors at perturbed positions
  Epetra_SerialDenseMatrix  checkmat1(sys_mat);
  Epetra_SerialDenseVector  checkvec1(residual);

  // echo to screen
  printf("+-------------------------------------------+\n");
  printf("| FINITE DIFFERENCE CHECK FOR ELEMENT %5d |\n",ele->Id());
  printf("+-------------------------------------------+\n");
  printf("\n");

  // loop columns of matrix by looping nodes and then dof per nodes
  // loop nodes
  for(int nn=0;nn<nen_;++nn)
  {
    printf("-------------------------------------\n");
    printf("-------------------------------------\n");
    printf("NODE of element local id %d\n",nn);
    // loop dofs
    for(int rr=0;rr<numdofpernode_;++rr)
    {
      // number of the matrix column to check
      int dof=nn*(numdofpernode_)+rr;

      // clear element matrices and vectors to assemble
      checkmat1.Scale(0.0);
      checkvec1.Scale(0.0);

      // first put the non-perturbed values to the working arrays
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          ephinp_[k](i) = origephinp[k](i);
          ehist_[k](i)  = origehist[k](i);
        }
      } // for i

      printf("concentration dof %d (%d)\n",rr,nn);

      if (is_genalpha_)
      {
        // perturbation of phi(n+1) in phi(n+alphaF) => additional factor alphaF
        ephinp_[rr](nn)+=(alphaF*epsilon);

        // perturbation of solution variable phi(n+1) for gen.alpha
        // leads to perturbation of phidtam (stored in ehist_)
        // with epsilon*alphaM/(gamma*dt)
        const double factor = alphaF/timefac; // = alphaM/(gamma*dt)
        ehist_[rr](nn)+=(factor*epsilon);

      }
      else
      {
        ephinp_[rr](nn,0)+=epsilon;
      }

      // calculate the right hand side for the perturbed vector
      Sysmat(
        ele,
        checkmat1,
        checkvec1,
        time,
        dt,
        timefac,
        alphaF,
        scatratype);

      // compare the difference between linaer approximation and
      // (nonlinear) right hand side evaluation

      // note that it makes more sense to compare these quantities
      // than to compare the matrix entry to the difference of the
      // the right hand sides --- the latter causes numerical problems
      // do to deletion //gammi

      // however, matrix entries delivered from the element are compared
      // with the finite-difference suggestion, too. It works surprisingly well
      // for epsilon set to 1e-6 (all displayed digits nearly correct)
      // and allows a more obvious comparison!
      // when matrix entries are small, lin. and nonlin. approximation
      // look identical, although the matrix entry may be rubbish!
      // gjb

      for(int mm=0;mm<(numdofpernode_*nen_);++mm)
      {
        double val   =-residual(mm)/epsilon;
        double lin   =-residual(mm)/epsilon+sys_mat(mm,dof);
        double nonlin=-checkvec1(mm)/epsilon;

        double norm=abs(lin);
        if(norm<1e-12)
        {
          norm=1e-12;
          std::cout<<"warning norm of lin is set to 10e-12"<<std::endl;
        }

        // output to screen
        {
          printf("relerr  %+12.5e   ",(lin-nonlin)/norm);
          printf("abserr  %+12.5e   ",lin-nonlin);
          printf("orig. value  %+12.5e   ",val);
          printf("lin. approx. %+12.5e   ",lin);
          printf("nonlin. funct.  %+12.5e   ",nonlin);
          printf("matrix[%d,%d]  %+12.5e   ",mm,dof,sys_mat(mm,dof));
          // finite difference approximation (FIRST divide by epsilon and THEN subtract!)
          // ill-conditioned operation has to be done as late as possible!
          printf("FD suggestion  %+12.5e ",((residual(mm)/epsilon)-(checkvec1(mm)/epsilon)) );
          printf("\n");
        }
      }
    }
  } // loop nodes

  // undo changes in state variables
  for (int i=0;i<nen_;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      ephinp_[k](i) = origephinp[k](i);
      ehist_[k](i)  = origehist[k](i);
    }
  } // for i

  return;
}//MeshfreeScaTraImpl::FDcheck
