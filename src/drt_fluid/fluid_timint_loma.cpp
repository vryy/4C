/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_loma.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/sutherland.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLoma::TimIntLoma(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      thermpressaf_(1.0),
      thermpressam_(1.0),
      thermpressdtaf_(0.0),
      thermpressdtam_(0.0)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::Init()
{
  // conservative formulation currently not supported in low-Mach-number case
  // when using generalized-alpha time-integration scheme
  if (convform_ == "conservative")
    dserror(
        "conservative formulation currently not supported for low-Mach-number flow within "
        "generalized-alpha time-integration scheme");

  // ---------------------------------------------------------------------
  // set density variable to 1.0 and get gas constant for low-Mach-number
  // flow and get constant density variable for incompressible flow
  // ---------------------------------------------------------------------

  // get gas constant
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_sutherland);
  if (id == -1)
    dserror("Could not find sutherland material");
  else
  {
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::Sutherland* actmat = static_cast<const MAT::PAR::Sutherland*>(mat);
    // we need the kinematic viscosity here
    gasconstant_ = actmat->gasconst_;
  }

  // potential check here -> currently not executed
  // if (gasconstant_ < EPS15) dserror("received zero or negative gas constant");

  // set some Loma-specific parameters
  SetElementCustomParameter();
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntLoma::~TimIntLoma() { return; }

/*----------------------------------------------------------------------*
 | set fields for scatra - fluid coupling, esp.                         |
 | set fields for low-Mach-number flow within iteration loop   vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
    Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
    Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
    Teuchos::RCP<DRT::Discretization> scatradis)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  //--------------------------------------------------------------------------
  // Filling the scaaf-vector and scaam-vector at time n+alpha_F/n+1 and
  // n+alpha_M/n, respectively, with scalar at pressure dofs
  // Additionally, filling the scaam-vector at time n+alpha_M/n with
  // velocity at time n at velocity dofs for OST/BDF2
  // Filling the accam-vector at time n+alpha_M/n+1, respectively, with
  // scalar time derivative values at pressure dofs
  //--------------------------------------------------------------------------
  // get velocity values at time n in scaam-vector as copy from veln-vector
  scaam_->Update(1.0, *veln_, 0.0);

  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0, lscatranode);
    const int globalscatradofid = scatradis->Dof(0, lscatranode, numscatradof - 1);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0) dserror("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(0, lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int numdof = discret_->NumDof(0, lnode);
    const int globaldofid = discret_->Dof(0, lnode, numdof - 1);
    const int localdofid = scaam_->Map().LID(globaldofid);
    if (localdofid < 0) dserror("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into scaaf_");

    value = (*scalaram)[localscatradofid];
    err = scaam_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into scaam_");

    if (scalardtam != Teuchos::null)
    {
      value = (*scalardtam)[localscatradofid];
    }
    else
    {
      value = 0.0;  // for safety reasons: set zeros in accam_
    }
    err = accam_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into accam_");

    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (fsscalaraf != Teuchos::null)
        value = (*fsscalaraf)[localscatradofid];
      else
        dserror("Expected fine-scale scalar!");

      err = fsscaaf_->ReplaceMyValue(localdofid, 0, value);
      if (err != 0) dserror("error while inserting value into fsscaaf_");
    }
  }

  //--------------------------------------------------------------------------
  // get thermodynamic pressure at n+alpha_F/n+1 and n+alpha_M/n and
  // time derivative of thermodyn. press. at n+alpha_F/n+1 and n+alpha_M/n+1
  //--------------------------------------------------------------------------
  thermpressaf_ = thermpressaf;
  thermpressam_ = thermpressam;
  thermpressdtaf_ = thermpressdtaf;
  thermpressdtam_ = thermpressdtam;

  return;
}  // TimIntLoma::SetLomaIterScalarFields


/*----------------------------------------------------------------------*
 | set scalar fields     vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
    const double thermpressnp, Teuchos::RCP<const Epetra_Vector> scatraresidual,
    Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar)
{
  FluidImplicitTimeInt::SetScalarFields(
      scalarnp, thermpressnp, scatraresidual, scatradis, whichscalar);
  //--------------------------------------------------------------------------
  // get thermodynamic pressure at n+1
  //--------------------------------------------------------------------------
  thermpressaf_ = thermpressnp;


  return;

}  // TimIntLoma::SetScalarFields

// -------------------------------------------------------------------
// set loma parameters                               rasthofer 03/2012
// -------------------------------------------------------------------
void FLD::TimIntLoma::SetElementCustomParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_loma_parameter);

  // set parameters to update material with subgrid-scale temperature
  // potential inclusion of additional subgrid-scale terms in continuity equation
  eleparams.sublist("LOMA") = params_->sublist("LOMA");
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") =
      params_->sublist("RESIDUAL-BASED STABILIZATION");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") =
      params_->sublist("MULTIFRACTAL SUBGRID SCALES");

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| print info about turbulence model (loma-specific)            bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::PrintTurbulenceModel()
{
  FluidImplicitTimeInt::PrintTurbulenceModel();

  if (physicaltype_ == INPAR::FLUID::loma and turbmodel_ == INPAR::FLUID::smagorinsky)
  {
    if (DRT::INPUT::IntegralValue<int>(params_->sublist("SUBGRID VISCOSITY"), "C_INCLUDE_CI"))
    {
      if (params_->sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") > 0.0)
      {
        std::cout << "with Yoshizawa constant Ci= ";
        std::cout << params_->sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") << "\n";
      }
      else
        dserror("Ci expected!");
    }
    else
      std::cout << "Yoshizawa constant Ci not included";

    std::cout << &std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in AssembleMatAndRHS                   bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1", thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n", thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1", thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1", thermpressdtam_);
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in ApplyNonlinearBoundaryConditions    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetCustomEleParamsApplyNonlinearBoundaryConditions(
    Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1", thermpressaf_);
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in LinearRelaxationSolve               bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetCustomEleParamsLinearRelaxationSolve(Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1", thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n", thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1", thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1", thermpressdtam_);
  return;
}

/*----------------------------------------------------------------------*
| call the statistics manager including thermpress parameters  bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::CallStatisticsManager()
{
  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
  // compute equation-of-state factor
  const double eosfac = thermpressaf_ / gasconstant_;
  statisticsmanager_->DoTimeSample(
      step_, eosfac, thermpressaf_, thermpressam_, thermpressdtaf_, thermpressdtam_);
  return;
}

/*----------------------------------------------------------------------*
 | prepare AVM3-based scale separation                         vg 10/08 |
 | overloaded in TimIntRedModels and TimIntLoma               bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::AVM3Preparation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // necessary here, because some application time integrations add something to the residual
  // before the Neumann loads are added
  residual_->PutScalar(0.0);

  eleparams.set("thermpress at n+alpha_F/n+1", thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n", thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1", thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1", thermpressdtam_);

  AVM3AssembleMatAndRHS(eleparams);

  // get scale-separation matrix
  AVM3GetScaleSeparationMatrix();

  // perform initial separation to initialize fsvelaf_
  // required for loma
  if (physicaltype_ == INPAR::FLUID::loma)
  {
    UpdateVelafGenAlpha();
    Sep_Multiply();
  }

  return;
}  // TimIntLoma::AVM3Preparation
