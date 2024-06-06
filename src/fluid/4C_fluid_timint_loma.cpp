/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_loma.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_turbulence_statistic_manager.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_sutherland.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLoma::TimIntLoma(const Teuchos::RCP<Discret::Discretization>& actdis,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
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
    FOUR_C_THROW(
        "conservative formulation currently not supported for low-Mach-number flow within "
        "generalized-alpha time-integration scheme");

  // ---------------------------------------------------------------------
  // set density variable to 1.0 and get gas constant for low-Mach-number
  // flow and get constant density variable for incompressible flow
  // ---------------------------------------------------------------------

  // get gas constant
  int id = Global::Problem::Instance()->Materials()->FirstIdByType(Core::Materials::m_sutherland);
  if (id == -1)
    FOUR_C_THROW("Could not find sutherland material");
  else
  {
    const Core::Mat::PAR::Parameter* mat =
        Global::Problem::Instance()->Materials()->ParameterById(id);
    const Mat::PAR::Sutherland* actmat = static_cast<const Mat::PAR::Sutherland*>(mat);
    // we need the kinematic viscosity here
    gasconstant_ = actmat->gasconst_;
  }

  // potential check here -> currently not executed
  // if (gasconstant_ < 1e-15) FOUR_C_THROW("received zero or negative gas constant");

  // set some Loma-specific parameters
  set_element_custom_parameter();
  return;
}



/*----------------------------------------------------------------------*
 | set fields for scatra - fluid coupling, esp.                         |
 | set fields for low-Mach-number flow within iteration loop   vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::set_loma_iter_scalar_fields(Teuchos::RCP<const Epetra_Vector> scalaraf,
    Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
    Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
    Teuchos::RCP<Discret::Discretization> scatradis)
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
    Core::Nodes::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0, lscatranode);
    const int globalscatradofid = scatradis->Dof(0, lscatranode, numscatradof - 1);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0) FOUR_C_THROW("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(0, lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int numdof = discret_->NumDof(0, lnode);
    const int globaldofid = discret_->Dof(0, lnode, numdof - 1);
    const int localdofid = scaam_->Map().LID(globaldofid);
    if (localdofid < 0) FOUR_C_THROW("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) FOUR_C_THROW("error while inserting value into scaaf_");

    value = (*scalaram)[localscatradofid];
    err = scaam_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) FOUR_C_THROW("error while inserting value into scaam_");

    if (scalardtam != Teuchos::null)
    {
      value = (*scalardtam)[localscatradofid];
    }
    else
    {
      value = 0.0;  // for safety reasons: set zeros in accam_
    }
    err = accam_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) FOUR_C_THROW("error while inserting value into accam_");

    if (turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales)
    {
      if (fsscalaraf != Teuchos::null)
        value = (*fsscalaraf)[localscatradofid];
      else
        FOUR_C_THROW("Expected fine-scale scalar!");

      err = fsscaaf_->ReplaceMyValue(localdofid, 0, value);
      if (err != 0) FOUR_C_THROW("error while inserting value into fsscaaf_");
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
}  // TimIntLoma::set_loma_iter_scalar_fields


/*----------------------------------------------------------------------*
 | set scalar fields     vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
    const double thermpressnp, Teuchos::RCP<const Epetra_Vector> scatraresidual,
    Teuchos::RCP<Discret::Discretization> scatradis, const int whichscalar)
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
void FLD::TimIntLoma::set_element_custom_parameter()
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
void FLD::TimIntLoma::print_turbulence_model()
{
  FluidImplicitTimeInt::print_turbulence_model();

  if (physicaltype_ == Inpar::FLUID::loma and turbmodel_ == Inpar::FLUID::smagorinsky)
  {
    if (Core::UTILS::IntegralValue<int>(params_->sublist("SUBGRID VISCOSITY"), "C_INCLUDE_CI"))
    {
      if (params_->sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") > 0.0)
      {
        std::cout << "with Yoshizawa constant Ci= ";
        std::cout << params_->sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA") << "\n";
      }
      else
        FOUR_C_THROW("Ci expected!");
    }
    else
      std::cout << "Yoshizawa constant Ci not included";

    std::cout << &std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in assemble_mat_and_rhs                   bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::set_custom_ele_params_assemble_mat_and_rhs(Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1", thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n", thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1", thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1", thermpressdtam_);
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in apply_nonlinear_boundary_conditions    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::set_custom_ele_params_apply_nonlinear_boundary_conditions(
    Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1", thermpressaf_);
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in linear_relaxation_solve               bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::set_custom_ele_params_linear_relaxation_solve(
    Teuchos::ParameterList& eleparams)
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
void FLD::TimIntLoma::call_statistics_manager()
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
void FLD::TimIntLoma::av_m3_preparation()
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

  av_m3_assemble_mat_and_rhs(eleparams);

  // get scale-separation matrix
  av_m3_get_scale_separation_matrix();

  // perform initial separation to initialize fsvelaf_
  // required for loma
  if (physicaltype_ == Inpar::FLUID::loma)
  {
    UpdateVelafGenAlpha();
    Sep_Multiply();
  }

  return;
}  // TimIntLoma::av_m3_preparation

FOUR_C_NAMESPACE_CLOSE
