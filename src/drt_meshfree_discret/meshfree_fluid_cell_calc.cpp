/*!---------------------------------------------------------------------------

\file meshfree_fluid_cell_calc.cpp

\brief main file containing routines for calculation of meshfree fluid cell

\maintainer Keijo Nissen

\level 3

*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell_calc.H"             // class declarations

#include "meshfree_fluid_cell.H"                  //
#include "drt_meshfree_utils.H"                    //
#include "drt_meshfree_node.H"                    //
#include "drt_meshfree_cell.H"                    //
#include "drt_meshfree_cell_utils.H"              // to get Gauss points in real space
#include "drt_meshfree_discret.H"                 // for cast to get points
#include "../drt_fluid_ele/fluid_ele_parameter.H"
#include "../drt_mat/newtonianfluid.H"

#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::MeshfreeFluidCellCalc():
  eid_(-1.0),
  is_inflow_ele_(false),
  kxyz_(nsd_,nek_),
  gxyz_(nsd_,ngp_),
  gw_(ngp_),
  sfunct_(Teuchos::rcp(new LINALG::SerialDenseVector())),
  sderiv_(Teuchos::rcp(new LINALG::SerialDenseMatrix())),
  wfunct_(Teuchos::rcp(new LINALG::SerialDenseVector())),
  wderiv_(Teuchos::rcp(new LINALG::SerialDenseMatrix())),
  vderxy_(false),
  bodyforce_(false),
  histmom_(false),
  velint_(false),
  gridvelint_(false),
  convvelint_(false),
  accint_(true),        // need to be set true
  gradp_(false),
  conv_c_(false),
  vdiv_(0.0),
  rhsmom_(false),
  conv_old_(false),
  det_(0.0),
  fac_(0.0),
  visc_(0.0),
  reacoeff_(0.0),
  densaf_(1.0),         // initialized to 1.0 (filled in Fluid::GetMaterialParams)
  densam_(1.0),         // initialized to 1.0 (filled in Fluid::GetMaterialParams)
  densn_(1.0)           // initialized to 1.0 (filled in Fluid::GetMaterialParams)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
}

/*--------------------------------------------------------------------------*
 |  generic evaluation of the element                    (public) nis Jan13 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::Evaluate(
  DRT::ELEMENTS::MeshfreeFluid* cell,
  DRT::Discretization &         discretization,
  const std::vector<int> &      lm,
  Teuchos::ParameterList&       params,
  Teuchos::RCP<MAT::Material> & mat,
  Epetra_SerialDenseMatrix&     elemat1_epetra,
  Epetra_SerialDenseMatrix&     elemat2_epetra,
  Epetra_SerialDenseVector&     elevec1_epetra,
  Epetra_SerialDenseVector&     elevec2_epetra,
  Epetra_SerialDenseVector&     elevec3_epetra,
  bool                          offdiag)
{
  //----------------------------------------------------------------------
  // cast discretization pointer to derived meshfree type
  //----------------------------------------------------------------------
  discret_ = dynamic_cast<DRT::MESHFREE::MeshfreeDiscretization*>( &(discretization) );
  if (discret_==NULL)
    dserror("dynamic_cast of discretization to meshfree discretization failed!");

  // ---------------------------------------------------------------------
  // construct views on element matrices and vectors
  // ---------------------------------------------------------------------

  LINALG::SerialDenseMatrix elemat1(View,elemat1_epetra.A(),elemat1_epetra.LDA(),elemat1_epetra.M(),elemat1_epetra.N(),false);
  LINALG::SerialDenseMatrix elemat2; // not used and error when constructing view; probably not initialised
  LINALG::SerialDenseVector elevec1(View,elevec1_epetra.Values(),elevec1_epetra.Length());
  // elevec2 and elevec3 are currently not in use

  // ---------------------------------------------------------------------
  // set size of all vectors of SerialDense element arrays
  // ---------------------------------------------------------------------

  // get number of nodes
  nen_ = cell->NumNode();

  // resize matrices and vectors
  conv_c_.LightSize(nen_);

  // get global node coordinates
  nxyz_.LightShape(nsd_,nen_);
  double const * cnxyz;
  for (int j=0; j<nen_; j++){
    cnxyz =  cell->Nodes()[j]->X();
    for (int k=0; k<nsd_; k++){
      nxyz_(k,j) = cnxyz[k];
    }
  }

  // get global point coordinates
  double const * ckxyz;
  for (int j=0; j<nek_; j++){
    ckxyz =  cell->Points()[j]->X();
    for (int k=0; k<nsd_; k++){
      kxyz_(k,j) = ckxyz[k];
    }
  }

  // ---------------------------------------------------------------------
  // get all general state vectors: velocity/pressure and history
  // velocity/pressure are at time n+alpha_F/n+alpha_M for generalized-alpha
  // scheme and at time n+1/n for all other schemes acceleration values are at
  // time n+alpha_M for generalized-alpha scheme and at time n+1 for all other
  // schemes
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::SerialDenseMatrix evelaf(nsd_,nen_,true);// need to be set to true??
  LINALG::SerialDenseVector epreaf(nen_,true);// need to be set to true??
  ExtractValuesFromGlobalVector(discretization,lm, &evelaf, &epreaf,"velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::SerialDenseMatrix evelnp(nsd_,nen_,true);// need to be set to true??
  LINALG::SerialDenseVector eprenp(nen_,true);// need to be set to true??
  if (fldparatimint_->IsGenalphaNP())
    ExtractValuesFromGlobalVector(discretization,lm, &evelnp, &eprenp,"velnp");

  LINALG::SerialDenseMatrix emhist(nsd_,nen_,true);// need to be set to true??
  ExtractValuesFromGlobalVector(discretization,lm, &emhist, NULL,"hist");

  LINALG::SerialDenseMatrix eaccam(nsd_,nen_,true);// need to be set to true??
  ExtractValuesFromGlobalVector(discretization,lm,  &eaccam, NULL,"accam");

  if (!fldparatimint_->IsGenalpha())
    eaccam.Zero();

  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::SerialDenseMatrix edispnp(nsd_,nen_,true);// need to be set to true??
  LINALG::SerialDenseMatrix egridv(nsd_,nen_,true);// need to be set to true??

  if (cell->IsAle())
  {
    // stationary formulation does not support ALE formulation
    if (fldparatimint_->IsStationary())
      dserror("No ALE support within stationary fluid solver.");

    // ALE in Oseen and Stokes problems questionable
    switch (fldpara_->PhysicalType())
    {
    case INPAR::FLUID::oseen:
    case INPAR::FLUID::stokes:
    {
      dserror("ALE with Oseen or Stokes seems to be a tricky combination. Think deep before removing dserror!");
      break;
    }
    default:
    {
      ExtractValuesFromGlobalVector(discretization,lm, &edispnp, NULL,"dispnp");
      ExtractValuesFromGlobalVector(discretization,lm, &egridv, NULL,"gridv");
    }
    }
  }

  // identify elements of inflow section
  InflowElement(cell);

  // set element id
  eid_ = cell->Id();

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = Evaluate(
    cell,
    params,
    elemat1,
    elemat2,
    elevec1,
    evelaf,
    epreaf,
    eprenp,
    evelnp,
    emhist,
    eaccam,
    edispnp,
    egridv,
    mat,
    cell->IsAle(),
    offdiag);

  return result;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::Evaluate(
  const DRT::ELEMENTS::MeshfreeFluid* cell,
  const Teuchos::ParameterList &      params,
  LINALG::SerialDenseMatrix &         elemat1,
  LINALG::SerialDenseMatrix &         elemat2,
  LINALG::SerialDenseVector &         elevec1,
  const LINALG::SerialDenseMatrix &   evelaf,
  const LINALG::SerialDenseVector &   epreaf,
  const LINALG::SerialDenseVector &   eprenp,
  const LINALG::SerialDenseMatrix &   evelnp,
  const LINALG::SerialDenseMatrix &   emhist,
  const LINALG::SerialDenseMatrix &   eaccam,
  const LINALG::SerialDenseMatrix &   edispnp,
  const LINALG::SerialDenseMatrix &   egridv,
  Teuchos::RCP<const MAT::Material>   mat,
  const bool                          isale,
  const bool                          offdiag
  )
{
  if (offdiag)
    dserror("No-off-diagonal matrix evaluation in standard fluid implementation!!");

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(
    cell,
    evelaf,
    evelnp,
    epreaf,
    eprenp,
    eaccam,
    emhist,
    edispnp,
    egridv,
    elemat1,
    elemat2,  // -> emesh
    elevec1,
    mat,
    isale);

  return 0;
}

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::Sysmat(
  const DRT::ELEMENTS::MeshfreeFluid* cell,
  const LINALG::SerialDenseMatrix&    evelaf,
  const LINALG::SerialDenseMatrix&    evelnp,
  const LINALG::SerialDenseVector&    epreaf,
  const LINALG::SerialDenseVector&    eprenp,
  const LINALG::SerialDenseMatrix&    eaccam, // ??
  const LINALG::SerialDenseMatrix&    emhist,
  const LINALG::SerialDenseMatrix&    edispnp,
  const LINALG::SerialDenseMatrix&    egridv,
  LINALG::SerialDenseMatrix&          estif,
  LINALG::SerialDenseMatrix&          emesh,
  LINALG::SerialDenseVector&          eforce,
  Teuchos::RCP<const MAT::Material>   material,
  const bool                          isale
  )
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of global matrices
  LINALG::SerialDenseMatrix estif_u(nen_*nsd_,nen_*nsd_,true);
  LINALG::SerialDenseMatrix estif_p_v(nen_*nsd_,nen_,true);
  LINALG::SerialDenseMatrix estif_q_u(nen_, nen_*nsd_,true);
  LINALG::SerialDenseMatrix ppmat(nen_,nen_,true);

  // definition of global RHS force vectors
  LINALG::SerialDenseVector preforce(nen_,true);
  LINALG::SerialDenseMatrix velforce(nsd_,nen_,true);

  // definition of velocity-based momentum residual vectors
  LINALG::SerialDenseMatrix lin_resM_Du(nsd_*nsd_,nen_,false);
  LINALG::Matrix<nsd_,1>    resM_Du(false);

  // add displacement when fluid nodes and points move in the ALE case
  if (isale)
  {
    kxyz_ += edispnp;
    nxyz_ += edispnp;
  }

  // if polynomial pressure projection: reset variables
  if (fldpara_->PPP())
  {
    D_ = 0;
    Eu_.LightResize(nen_);
    Eu_.Zero();
    Fu_.LightReshape(nsd_,nen_);
    Fu_.Zero();
    Ev_.LightResize(nen_);
    Ev_.Zero();
    Fv_.LightReshape(nsd_,nen_);
    Fv_.Zero();
  }

  //------------------------------------------------------------------------
  // potential evaluation of material parameters at element center
  //------------------------------------------------------------------------
  if (not fldpara_->MatGp())
  {
    // get material parameters at element center
    GetMaterialParams(material,evelaf);
  }

  //----------------------------------------------------------------------
  // get integrations points and weights in xyz-system
  //----------------------------------------------------------------------
  int ngp = DRT::MESHFREE::CellGaussPointInterface::Impl(distype)->GetCellGaussPointsAtX(kxyz_, gxyz_, gw_);

  LINALG::SerialDenseMatrix distng(nsd_,nen_); // matrix for distance between node and Gauss point
  double const * cgxyz; // read: current Gauss xyz-coordinate
  double const * cnxyz; // read: current node xyz-coordinate

  //------------------------------------------------------------------------
  //  loop over integration points for current cell
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<ngp; ++iquad)
  {
    //----------------------------------------------------------------------
    // get basis function values at current Gauss point
    //----------------------------------------------------------------------

    // get xyz-coordinates and weight of current Gauss point
    cgxyz = gxyz_[iquad]; // read: current gauss xyz-coordinate
    fac_ = gw_[iquad];    // read: current gauss weight

    // coordinates of the current integration point
    for (int i=0; i<nen_; ++i){
      // get current node xyz-coordinate
      cnxyz = nxyz_[i];
      for (int j=0; j<nsd_; ++j){
        // get distance between
        distng(j,i) = cnxyz[j] - cgxyz[j];
      }
    }

    // calculate basis functions and derivatives via max-ent optimization
    int err = discret_->GetSolutionApprox()->GetMeshfreeBasisFunction(nsd_,Teuchos::rcpFromRef(distng),sfunct_,sderiv_);
    if (err>0)
    {
      std::cout << "When computing the solution basis functions at gauss point " << iquad << " of cell " << cell->Id() << ":" << std::endl;
      DRT::MESHFREE::OutputMeshfreeError(err);
    }

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including derivatives and grid velocity)
    //  2) pressure (including derivatives)
    //  3) body-force vector
    //  4) "history" vector for momentum equation
    //----------------------------------------------------------------------

    // construct SerialDense-view on velint_
    LINALG::SerialDenseVector velint_sdm(View, velint_.A(), nsd_);
    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_sdm.Multiply('N','N',1.0,evelaf,*sfunct_,0.0);

    // construct SerialDense-view on vderxy_
    LINALG::SerialDenseMatrix vderxy_sdm(View, vderxy_.A(), nsd_, nsd_, nsd_);
    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_sdm.Multiply('N','T',1.0,evelaf,*sderiv_,0.0);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    double press(true);
    if(fldparatimint_->IsGenalphaNP())
      press = sfunct_->Dot(eprenp);
    else
      press = sfunct_->Dot(epreaf);

    // construct SerialDense-view on gradp_
    LINALG::SerialDenseVector gradp_sdm(View, gradp_.A(), nsd_);
    // get pressure gradient at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    if(fldparatimint_->IsGenalphaNP())
      gradp_sdm.Multiply('N','N',1.0,*sderiv_,eprenp,0.0);
    else
      gradp_sdm.Multiply('N','N',1.0,*sderiv_,epreaf,0.0);

    // construct SerialDense-view on histmom_
    LINALG::SerialDenseVector histmom_sdm(View, histmom_.A(), nsd_);
    // get momentum history data at integration point
    // (only required for one-step-theta and BDF2 time-integration schemes)
    histmom_sdm.Multiply('N','N',1.0,emhist,*sfunct_,0.0);

    //----------------------------------------------------------------------
    // potential evaluation of material parameters
    //----------------------------------------------------------------------

    // get material parameters at integration point
    if (fldpara_->MatGp())
    {
      GetMaterialParams(material,evelaf);
    }

    //----------------------------------------------------------------------
    //  evaluation of various velocities at integration point
    //----------------------------------------------------------------------

    // get convective velocity at integration point
    switch (fldpara_->PhysicalType())
    {
    case INPAR::FLUID::incompressible:
    {
      convvelint_ = velint_;
      break;
    }
    case INPAR::FLUID::oseen:
    {
      const int advefuncno = fldpara_->OseenFieldFuncNo();
      const double time = fldparatimint_->Time();
      for(int idim=0;idim<nsd_;++idim)
        convvelint_(idim) = DRT::Problem::Instance()->Funct(advefuncno-1).Evaluate(idim,cgxyz,time);
      break;
    }
    case INPAR::FLUID::stokes:
    {
      convvelint_.Clear();
      break;
    }
    default:
      dserror("Physical type not implemented in meshfree fluid.");
    }

    // (ALE case handled implicitly here using the (potential
    //  mesh-movement-dependent) convective velocity, avoiding
    //  various ALE terms used to be calculated before)
    if (isale)
    {
      dserror("Meshfree fluid not tested for ALE, yet. Remove dserror at own risk.");
      // construct SerialDense-view on gridvelint_
      LINALG::SerialDenseVector gridvelint_sdm(View, gridvelint_.A(), nsd_);
      gridvelint_sdm.Multiply('N','N',1.0,egridv,*sfunct_,0.0);
      convvelint_.Update(-1.0,gridvelint_,1.0);
    }

    // compute convective term from previous iteration
    conv_old_.MultiplyNN(1.0,vderxy_,convvelint_);
    // construct SerialDense-view on convvelint_
    LINALG::SerialDenseVector convvelint_sdm(View, convvelint_.A(), nsd_);
    // compute convective term of convective operator
    conv_c_.Multiply('T','N',1.0,*sderiv_,convvelint_sdm,0.0);

    // compute divergence of velocity from previous iteration
    vdiv_ = 0.0;
    if (not fldparatimint_->IsGenalphaNP())
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        vdiv_ += vderxy_(idim, idim);
      }
    }
    else
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        //get vdiv at time n+1 for np_genalpha,
        LINALG::SerialDenseMatrix vderxy(nsd_,nsd_);
        vderxy.Multiply('N','T',1.0,evelnp,*sderiv_,0.0); // not vderxy_ because of 'evelnp'
        vdiv_ += vderxy(idim, idim);
      }
    }

    // calculate weighting basis functions and derivatives via max-ent optimization
    Teuchos::RCP<LINALG::SerialDenseVector> upsilon = Teuchos::rcp(new LINALG::SerialDenseVector(Copy,convvelint_.A(),nsd_));
    upsilon->Scale(-1.0/visc_);
    err = discret_->GetWeightingApprox()->GetMeshfreeBasisFunction(nsd_,Teuchos::rcpFromRef(distng),wfunct_,wderiv_,upsilon,sfunct_,sderiv_);
    if (err>0)
    {
      std::cout << "When computing the weighting basis functions at gauss point " << iquad << " of cell " << cell->Id() << ":" << std::endl;
      DRT::MESHFREE::OutputMeshfreeError(err);
    }

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac    = fldparatimint_->TimeFac()    * fac_; // ??
    const double timefacfacpre = fldparatimint_->TimeFacPre() * fac_; // ??
    const double rhsfac        = fldparatimint_->TimeFacRhs() * fac_; // ??

    //----------------------------------------------------------------------
    // compute residual of momentum equation
    // -> different for generalized-alpha and other time-integration schemes
    //----------------------------------------------------------------------
    if (fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
      dserror("Boussinesq approximation has not been implemented for meshfree fluid yet!");

    // compute bodyforce at current gauss point
    BodyForce(cell,cgxyz);
    rhsmom_.Update(densaf_,bodyforce_,0.0);

    if (fldparatimint_->IsGenalpha())
    {
      // construct SerialDense-view on accint_
      LINALG::SerialDenseVector accint_sdm(View, accint_.A(), nsd_);
      // get acceleration at time n+alpha_M at integration point
      accint_sdm.Multiply('N','N',1.0,eaccam,*sfunct_,0.0);
    }
    else
    {
      if (not fldparatimint_->IsStationary())
      {
        rhsmom_.Update((densn_/fldparatimint_->Dt()/fldparatimint_->Theta()),histmom_,1.0);
      }
    }

    // set velocity-based momentum residual vectors to zero
    lin_resM_Du.Zero();
    resM_Du.Clear();

    // compute first version of velocity-based momentum residual containing
    // inertia term, convection term (convective and reactive part),
    // reaction term and cross-stress term
    LinGalMomResU(lin_resM_Du,
                  timefacfac);

    //----------------------------------------------------------------------
    // computation of standard Galerkin and stabilization contributions to
    // element matrix and right-hand-side vector
    //----------------------------------------------------------------------
    // 1) standard Galerkin inertia, convection and reaction terms
    //    (convective and reactive part for convection term)
    //    as well as first part of cross-stress term on left-hand side
    InertiaConvectionReactionGalPart(estif_u,
                                     velforce,
                                     lin_resM_Du,
                                     resM_Du,
                                     rhsfac);

    // 2) standard Galerkin viscous term
    //    (including viscous stress computation,
    //     excluding viscous part for low-Mach-number flow)
    LINALG::Matrix<nsd_,nsd_> viscstress(true);
    ViscousGalPart(estif_u,
                   velforce,
                   viscstress,
                   timefacfac,
                   rhsfac);

    if (not fldpara_->PPP())
    {
      // 3) standard Galerkin pressure term
      PressureGalPart(estif_p_v,
                      velforce,
                      timefacfac,
                      timefacfacpre,
                      rhsfac,
                      press);

      // 4) standard Galerkin continuity term
      ContinuityGalPart(estif_q_u,
                        preforce,
                        timefacfac,
                        timefacfacpre,
                        rhsfac);
    }

    // 5) standard Galerkin bodyforce term on right-hand side
    BodyForceRhsTerm(velforce,
                     rhsfac);

    // 6) additional standard Galerkin terms due to conservative formulation
    if (fldpara_->IsConservative())
    {
      ConservativeFormulation(estif_u,
                              velforce,
                              timefacfac,
                              rhsfac);
    }

    // 7) inf-sup stabilization if desired
    //  Dohrmann-Bochev polynomial pressure projection term with constant
    //  basis a = 1. See Dohrmann and Bochev 2004, IJNMF for notation.
    if (fldpara_->PPP())
    {
      PressureProjection(ppmat);
    }

  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  // if polynomial pressure projection: finish ppmat
  if (fldpara_->PPP())
  {
    if (fldparatimint_->IsGenalphaNP())
      PressureProjectionFinalize(estif_p_v, estif_q_u, ppmat, velforce, preforce, evelnp, eprenp);
    else
      PressureProjectionFinalize(estif_p_v, estif_q_u, ppmat, velforce, preforce, evelaf, epreaf);
  }

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------

  // add pressure part to right-hand-side vector
  for (int vi=0; vi<nen_; ++vi)
  {
    eforce(numdofpernode_*vi+nsd_)+=preforce(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      eforce(numdofpernode_*vi+idim)+=velforce(idim,vi);
    }
  }

  //------------------------------------------------------------------------
  //  add contributions to element matrix
  //------------------------------------------------------------------------

  // add pressure-pressure part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int fuippp = numdofpernode_*ui+nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int numdof_vi_p_nsd = numdofpernode_*vi+nsd_;

      estif(numdof_vi_p_nsd,fuippp)+=ppmat(vi,ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int numdof_vi = numdofpernode_*vi;
      const int nsd_vi = nsd_*vi;

      for (int jdim=0; jdim < nsd_;++jdim)
      {
        const int numdof_ui_jdim = numdof_ui+jdim;
        const int nsd_ui_jdim = nsd_ui+jdim;

        for (int idim=0; idim <nsd_; ++idim)
        {
          estif(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui_nsd = numdofpernode_*ui + nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int nsd_vi = nsd_*vi;
      const int numdof_vi = numdofpernode_*vi;

      for (int idim=0; idim <nsd_; ++idim)
      {
        estif(numdof_vi+idim, numdof_ui_nsd) += estif_p_v(nsd_vi+idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
        estif(numdofpernode_*vi+nsd_, numdof_ui_jdim) += estif_q_u(vi, nsd_ui_jdim);
    }
  }

  return;
} // end DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::Sysmat


/*----------------------------------------------------------------------*
 |  compute body force at element nodes (private)              vg 10/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::BodyForce(
  const DRT::ELEMENTS::MeshfreeFluid *      cell,
  const double *                            cgxyz)
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique Neumann condition
  if (nsd_==3)
    DRT::UTILS::FindElementConditions(cell, "VolumeNeumann", myneumcond);
  else if (nsd_==2)
    DRT::UTILS::FindElementConditions(cell, "SurfaceNeumann", myneumcond);
  else
    dserror("Body force for 1D problem not yet implemented!");

  if (myneumcond.size()>1)
    dserror("More than one Neumann condition on one node!");

  if (myneumcond.size()==1)
  {
    // initialization of time-curve factor
    const double time = fldparatimint_->Time();

    // get values and switches from the condition
    const std::vector<int>*    onoff     = myneumcond[0]->Get<std::vector<int   > >("onoff");
    const std::vector<double>* val       = myneumcond[0]->Get<std::vector<double> >("val"  );
    const std::vector<int>*    functions = myneumcond[0]->Get<std::vector<int   > >("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the ebofoaf array
    for (int isd=0;isd<nsd_;isd++)
    {
      // get factor given by spatial function
      if (functions) functnum = (*functions)[isd];
      else functnum = -1;

      double num = (*onoff)[isd]*(*val)[isd];

      // evaluate function at the position of the current node
      if (functnum>0)
        functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,cgxyz,time);
      else
        functionfac = 1.0;

      // get usual body force
      bodyforce_(isd) = num*functionfac;
    } // isd
  } // if neumanncond
  else
    bodyforce_.Clear();
} // end of void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::BodyForce

/*----------------------------------------------------------------------*
 |  compute material parameters                   (protected) nis Jan13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::GetMaterialParams(
  Teuchos::RCP<const MAT::Material> material,
  const LINALG::SerialDenseMatrix&  evelaf
)
{
  // initially set density values and values with respect to continuity rhs
  densam_        = 1.0;
  densaf_        = 1.0;
  densn_         = 1.0;

  if (material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

    // get constant dynamic viscosity
    visc_ = actmat->Viscosity();

    densaf_ = actmat->Density();
    densam_ = densaf_;
    densn_  = densaf_;
  }
  else
    dserror("Material type is not supported");

  return;
} // MeshfreeFluidCellCalc::GetMaterialParams

/*----------------------------------------------------------------------------*
 | calculate prefactor of instationary/convective/reaction terms of momentum  |
 | equations later multiplied with 'Du'                             nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::LinGalMomResU(
  LINALG::SerialDenseMatrix & lin_resM_Du,
  const double                timefacfac)
{
  /*
    instationary  convection, reactive        convection, convective      reaction
       +-----+  +-----------------------+   +------------------------+   +-------+
       |     |  |                       |   |                        |   |       |

                 /       n+1       \         /                \  n+1
       rho*Du + |   rho*u   o nabla | Du  + |   rho*Du o nabla | u     +  sigma*Du
                 \      (i)        /         \                /   (i)

  */

  int idim_nsd_p_idim[nsd_];

  for (int idim = 0; idim <nsd_; ++idim)
  {
    idim_nsd_p_idim[idim]=idim*nsd_+idim;
  }

  // instationary part
  if (not fldparatimint_->IsStationary())
  {
    const double fac_densam=fac_*densam_;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double v=fac_densam*(*sfunct_)(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }
  }

  const double timefacfac_densaf=timefacfac*densaf_;

  // convection, reactive
  for (int ui=0; ui<nen_; ++ui)
  {
    const double v=timefacfac_densaf*conv_c_(ui);

    for (int idim = 0; idim <nsd_; ++idim)
    {
      lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
    }
  }

  // convection, convective (only for Newton)
  if(fldpara_->IsNewton() and (not (fldpara_->PhysicalType()==INPAR::FLUID::stokes)))
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const double temp=timefacfac_densaf*(*sfunct_)(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        const int idim_nsd=idim*nsd_;

        for(int jdim=0;jdim<nsd_;++jdim)
        {
          lin_resM_Du(idim_nsd+jdim,ui)+=temp*vderxy_(idim,jdim);
        }
      }
    }
  }


  // reaction part
  if (fldpara_->Reaction())
  {
    const double fac_reac=timefacfac*reacoeff_;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double v=fac_reac*(*sfunct_)(ui);

      for (int idim = 0; idim <nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim],ui)+=v;
      }
    }
  }

  return;
} // end of MeshfreeFluidCellCalc::LinGalMomResU

/*----------------------------------------------------------------------------*
 |  calculate instationary/convective/reaction terms of mom eqs     nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::InertiaConvectionReactionGalPart(
  LINALG::SerialDenseMatrix & estif_u,
  LINALG::SerialDenseMatrix & velforce,
  const LINALG::SerialDenseMatrix & lin_resM_Du,
  LINALG::Matrix<nsd_,1> &    resM_Du,
  const double                rhsfac)
{
  //------------------------------------------------------------------------
  // multiply instationary/convective/reaction-aggregate with test function
  //------------------------------------------------------------------------

  if (fldpara_->IsNewton())
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        const int idim_nsd=idim*nsd_;

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = nsd_*vi;

          const int fvi_p_idim = fvi+idim;

          for (int jdim= 0; jdim<nsd_;++jdim)
          {
            estif_u(fvi_p_idim,fui+jdim) += (*wfunct_)(vi)*lin_resM_Du(idim_nsd+jdim,ui);
          } // end for (jdim)
        } // end for (idim)
      } //vi
    } // ui
  }
  else
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi;

        for (int idim = 0; idim <nsd_; ++idim)
        {
          estif_u(fvi+idim,fui+idim) += (*wfunct_)(vi)*lin_resM_Du(idim*nsd_+idim,ui);
        } // end for (idim)
      } //vi
    } // ui
  }

  //------------------------------------------------------------------------
  // compute instationary/convective/reaction aggregate for rhs
  //------------------------------------------------------------------------

  // inertia terms of rhs
  if (not fldparatimint_->IsStationary())
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      if (fldparatimint_->IsGenalpha()) // ??
        resM_Du(idim)+=rhsfac*densam_*accint_(idim);
      else
        resM_Du(idim)+=fac_*densaf_*velint_(idim);
    }
  } // end if (not stationary)

  // convective terms of rhs
  for (int idim = 0; idim <nsd_; ++idim)
  {
    resM_Du(idim)+=rhsfac*densaf_*conv_old_(idim);
  } // end for(idim)

  // reactive terms of rhs
  if (fldpara_->Reaction())
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      resM_Du(idim) += rhsfac*reacoeff_*velint_(idim);
    }
  } // end if (reaction_)

  //------------------------------------------------------------------------
  // multiply aggregate for rhs with test function and add to velforce
  //------------------------------------------------------------------------

  for (int vi=0; vi<nen_; ++vi)
  {
    for(int idim = 0; idim <nsd_; ++idim)
    {
      velforce(idim,vi)-=resM_Du(idim)*(*wfunct_)(vi);
    }
  }

  return;
} // end of MeshfreeFluidCellCalc::InertiaConvectionReactionGalPart(...)


/*----------------------------------------------------------------------------*
 |  calculate viscous term of the momentum equations                nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::ViscousGalPart(
  LINALG::SerialDenseMatrix & estif_u,
  LINALG::SerialDenseMatrix & velforce,
  LINALG::Matrix<nsd_,nsd_> & viscstress,
  const double                timefacfac,
  const double                rhsfac)
{
  const double visc_timefacfac = visc_*timefacfac;

  /* viscosity term */
  /*
                   /                        \
                  |       /  \         / \   |
            2 mu  |  eps | Du | , eps | v |  |
                  |       \  /         \ /   |
                   \                        /
  */
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = nsd_*vi;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = nsd_*ui;
      const double * sderivs = (*sderiv_)[ui];

      for (int jdim=0; jdim<nsd_;++jdim)
      {
        const double temp=visc_timefacfac*(*wderiv_)(jdim,vi);

        for (int idim=0; idim <nsd_; ++idim)
        {
          const int fvi_p_idim = fvi+idim;

          estif_u(fvi_p_idim,fui+jdim) += temp*sderivs[idim];
          estif_u(fvi_p_idim,fui+idim) += temp*sderivs[jdim];
        } // end for (idim)
      } // ui
    } // end for (idim)
  } // vi


//  // solving weak form of: - \visc \laplace u + \grad p = f
//  for (int vi=0; vi<nen_; ++vi)
//  {
//    const int fvi = nsd_*vi;
//
//    for (int ui=0; ui<nen_; ++ui)
//    {
//      const int fui = nsd_*ui;
//
//      for (int jdim=0; jdim<nsd_;++jdim)
//      {
//
//        for (int idim=0; idim <nsd_; ++idim)
//        {
//
//          estif_u(fvi+jdim,fui+jdim) += visc_timefacfac*(*wderiv_)(idim,vi)*(*sderiv_)(idim, ui);
//
//        } // end for (idim)
//      } // ui
//    } // end for (idim)
//  } // vi

  for (int jdim=0; jdim<nsd_; ++jdim)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      viscstress(idim,jdim) = vderxy_(jdim,idim)+vderxy_(idim,jdim);
    }
  }
  viscstress.Scale(visc_*rhsfac);

  // computation of right-hand-side viscosity term
  for (int vi=0; vi<nen_; ++vi)
  {
    const double * wderivs = (*wderiv_)[vi];

    for (int idim=0; idim<nsd_; ++idim)
    {
      for (int jdim=0; jdim<nsd_; ++jdim)
      {
        /* viscosity term on right-hand side */
        velforce(idim,vi) -= viscstress(idim,jdim)*wderivs[jdim];
      }
    }
  }

  return;
} // end of MeshfreeFluidCellCalc::ViscousGalPart

/*----------------------------------------------------------------------------*
 |  calculate pressure term of the momentum equations               nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::PressureGalPart(
  LINALG::SerialDenseMatrix & estif_p_v,
  LINALG::SerialDenseMatrix & velforce,
  const double                timefacfac,
  const double                timefacfacpre,
  const double                rhsfac,
  const double                press)
{
  // pressure term on system matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const double v = -timefacfacpre*(*sfunct_)(ui);

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = nsd_*vi;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        /* pressure term */
        /*
           /                \
          |                  |
          | -Dp , nabla o v  |
          |                  |
           \                /
        */
        estif_p_v(fvi + idim,ui) += v*(*wderiv_)(idim, vi);
      }
    }
  }

  // pressure term on right-hand side
  velforce.Update(press*rhsfac,*wderiv_,1.0);

  return;
} // end of MeshfreeFluidCellCalc::PressureGalPart

/*----------------------------------------------------------------------------*
 |  calculate continuity equations                                  nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::ContinuityGalPart(
  LINALG::SerialDenseMatrix & estif_q_u,
  LINALG::SerialDenseVector & preforce,
  const double                timefacfac,
  const double                timefacfacpre,
  const double                rhsfac)
{
  // continuity term on system matrix
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = timefacfacpre*(*wfunct_)(vi);

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        /* continuity term */
        /*
             /                \
            |                  |
            | nabla o Du  , q  |
            |                  |
             \                /
        */
        estif_q_u(vi,fui+idim) += v*(*sderiv_)(idim,ui);
      }
    }
  }  // end for(idim)

  // continuity term on right-hand side
  preforce.Update(-rhsfac * vdiv_,*wfunct_,1.0);

  return;
} // end of MeshfreeFluidCellCalc::ContinuityGalPart

/*----------------------------------------------------------------------------*
 |  add integration point contribution to pressure matrices         nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::PressureProjection(
  LINALG::SerialDenseMatrix & ppmat
  )
{
  // mass matrix of pressure basis functions - used as temp for M
  ppmat.Multiply('N','T',fac_, *wfunct_, *sfunct_, 1.0);

  // "mass matrix" of projection modes - here element volume_equivalent_diameter_pc
  D_ += fac_;

  // prolongator(?) of projection with solution basis functions
  Eu_.Update(fac_, *sfunct_, 1.0);

  // integrated discrete velocity divergence with solution basis functions
  Fu_.Update(fac_, *sderiv_, 1.0);

  // prolongator(?) of projection with weighting basis functions
  Ev_.Update(fac_, *wfunct_, 1.0);

  // integrated discrete velocity divergence with weighting basis functions
  Fv_.Update(fac_, *wderiv_, 1.0);

}

/*----------------------------------------------------------------------------*
 |  finalize pressure projection and compute rhs-contribution       nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::PressureProjectionFinalize(
  LINALG::SerialDenseMatrix & estif_p_v,
  LINALG::SerialDenseMatrix & estif_q_u,
  LINALG::SerialDenseMatrix & ppmat,
  LINALG::SerialDenseMatrix & velforce,
  LINALG::SerialDenseVector & preforce,
  const LINALG::SerialDenseMatrix & evel,
  const LINALG::SerialDenseVector & epre
  )
{
  //------------------------------------------------------------------------
  // pressure projection matrices
  //------------------------------------------------------------------------

  // compute stabilisation matrix as difference of consistent and projection
  // pressure mass matrices
  ppmat.Multiply('N','T',-1.0/D_,Ev_,Eu_,1.0);
  ppmat.Scale(1.0/visc_);

  // compute pressure-Galerkin and continuity-Galerkin term from F_
  LINALG::SerialDenseVector Fu_vec(View,Fu_.A(),nsd_*nen_);
  LINALG::SerialDenseVector Fv_vec(View,Fv_.A(),nsd_*nen_);
  estif_q_u.Multiply('N','T', 1.0/D_,Ev_,Fu_vec,0.0);
  estif_p_v.Multiply('N','T',-1.0/D_,Fv_vec,Eu_,0.0);

  //------------------------------------------------------------------------
  // rhs-contributions of pressure projection matrices
  //------------------------------------------------------------------------

  // compute velocity rhs-contribution from pressure-Galerkin term
  LINALG::SerialDenseVector temp_v(nsd_*nen_,false);
  temp_v.Multiply('N','N',1.0,estif_p_v,epre,0.0);

  // add velocity rhs-contributions
  LINALG::SerialDenseMatrix temp_v_mat(View,temp_v.A(),nsd_,nsd_,nen_);
  velforce.Update(-fldparatimint_->TimeFacRhs(),temp_v_mat,1.0);

  // compute pressure rhs-contribution from pressure-projection term
  LINALG::SerialDenseVector temp_p(nen_,false);
  temp_p.Multiply('N','N',1.0,ppmat,epre,0.0);
  // compute pressure rhs-contribution from continuity-Galerkin term
  LINALG::SerialDenseVector evel_vec(View,evel.A(),nsd_*nen_);
  temp_p.Multiply('N','N',1.0,estif_q_u,evel_vec,1.0);

  // add pressure rhs-contributions
  preforce.Update(-fldparatimint_->TimeFacRhs(),temp_p,1.0);

  //------------------------------------------------------------------------
  // final scale of pressure projection matrices
  //------------------------------------------------------------------------

  // scale matrices with pressure time factor
  estif_p_v.Scale(fldparatimint_->TimeFacPre());
  estif_q_u.Scale(fldparatimint_->TimeFacPre());
  ppmat.Scale(fldparatimint_->TimeFacPre());
}


/*----------------------------------------------------------------------------*
 |  calculate                                                       nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::BodyForceRhsTerm(
  LINALG::SerialDenseMatrix & velforce,
  const double                rhsfac)
{
  for (int idim = 0; idim <nsd_; ++idim)
  {
    const double scaled_rhsmom=rhsfac*rhsmom_(idim);

    for (int vi=0; vi<nen_; ++vi)
    {
      velforce(idim,vi)+=scaled_rhsmom*(*wfunct_)(vi);
    }
  }  // end for(idim)

  return;
} // end of MeshfreeFluidCellCalc::BodyForceRhsTerm


/*----------------------------------------------------------------------------*
 |  ??                                                              nis Jan13 |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MeshfreeFluidCellCalc<distype>::ConservativeFormulation(
    LINALG::SerialDenseMatrix & estif_u,
    LINALG::SerialDenseMatrix & velforce,
    const double                timefacfac,
    const double                rhsfac)
{
  //----------------------------------------------------------------------
  // computation of additions to convection term (convective and
  // reactive part) for conservative form of convection term including
  // right-hand-side contribution
  //----------------------------------------------------------------------

  // compute prefactor
  double v = timefacfac*densaf_*vdiv_;

  for (int idim = 0; idim <nsd_; ++idim)
  {
    //----------------------------------------------------------------------
    // matrix additions
    //----------------------------------------------------------------------

    /* convection, convective part (conservative addition) */
    /*
      /                                                 \
      |      /              n+1    n+1           \      |
      |  Du | rho*nabla o u    +  u   *nabla rho | , v  |
      |      \             (i)     (i)          /       |
      \                                                 /
    */

    for (int ui=0; ui<nen_; ++ui)
    {
      const int    fui   = nsd_*ui + idim;
      const double v1 = v*(*sfunct_)(ui);

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi + idim;
        estif_u(fvi  , fui  ) += (*wfunct_)(vi)*v1;
      }
    }

    /*  convection, reactive part (conservative addition) */
    /*
      /                              \
      |  n+1  /               \      |
      | u    | rho*nabla o Du | , v  |
      |  (i)  \              /       |
      \                             /
    */

    if (fldpara_->IsNewton())
    {
      const double v_idim = timefacfac*densaf_*velint_(idim);

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi + idim;
        const double v1_idim = v_idim*(*wfunct_)(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = nsd_*ui;

          for(int jdim=0; jdim<nsd_;++jdim)
            estif_u(fvi,  fui+jdim  ) += v1_idim*(*sderiv_)(jdim, ui);
        }
      }
    }

    //----------------------------------------------------------------------
    // right hand side additions
    //----------------------------------------------------------------------

    /* convection (conservative addition) on right-hand side */
    double temp = -rhsfac*densaf_*velint_(idim)*vdiv_;

    for (int vi=0; vi<nen_; ++vi)
      velforce(idim, vi    ) += temp*(*wfunct_)(vi);

  }  // end for(idim)

  return;
}

// template classes
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::MeshfreeFluidCellCalc<DRT::Element::tri3>;
