/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_loma.cpp
\brief TimIntLoma

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_loma.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/turbulence_statistic_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/sutherland.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntLoma::TimIntLoma(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      thermpressaf_(1.0),
      thermpressam_(1.0),
      thermpressdtaf_(0.0),
      thermpressdtam_(0.0)
{

  // conservative formulation currently not supported in low-Mach-number case
  // when using generalized-alpha time-integration scheme
  if (convform_ == "conservative")
     dserror("conservative formulation currently not supported for low-Mach-number flow within generalized-alpha time-integration scheme");

  // ---------------------------------------------------------------------
  // set density variable to 1.0 and get gas constant for low-Mach-number
  // flow and get constant density variable for incompressible flow
  // ---------------------------------------------------------------------

  // get gas constant
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_sutherland);
  if (id==-1)
    dserror("Could not find sutherland material");
  else
  {
    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::Sutherland* actmat = static_cast<const MAT::PAR::Sutherland*>(mat);
    // we need the kinematic viscosity here
    gasconstant_ = actmat->gasconst_;
  }

  // potential check here -> currently not executed
  //if (gasconstant_ < EPS15) dserror("received zero or negative gas constant");

  //set some Loma-specific parameters
  SetElementCustomParameter();
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntLoma::~TimIntLoma()
{
  return;
}

/*----------------------------------------------------------------------*
 | set fields for low-Mach-number flow within iteration loop   vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetIterLomaFields(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<const Epetra_Vector> scalaram,
   Teuchos::RCP<const Epetra_Vector> scalardtam,
   Teuchos::RCP<const Epetra_Vector> fsscalaraf,
   const double             thermpressaf,
   const double             thermpressam,
   const double             thermpressdtaf,
   const double             thermpressdtam,
   Teuchos::RCP<DRT::Discretization> scatradis)
{
  FluidImplicitTimeInt::SetIterLomaFields(scalaraf,scalaram,scalardtam,fsscalaraf,
      thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,scatradis);

  //--------------------------------------------------------------------------
  // get thermodynamic pressure at n+alpha_F/n+1 and n+alpha_M/n and
  // time derivative of thermodyn. press. at n+alpha_F/n+1 and n+alpha_M/n+1
  //--------------------------------------------------------------------------
  thermpressaf_   = thermpressaf;
  thermpressam_   = thermpressam;
  thermpressdtaf_ = thermpressdtaf;
  thermpressdtam_ = thermpressdtam;

  return;

} // TimIntLoma::SetIterLomaFields

/*----------------------------------------------------------------------*
 | set fields for low-Mach-number flow at end of time step     vg 09/09 |
 *----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetTimeLomaFields(
   Teuchos::RCP<const Epetra_Vector> scalarnp,
   const double             thermpressnp,
   Teuchos::RCP<const Epetra_Vector> scatraresidual,
   Teuchos::RCP<DRT::Discretization> scatradis,
   const int                whichscalar)
{

  FluidImplicitTimeInt::SetTimeLomaFields(scalarnp,thermpressnp,scatraresidual,scatradis,whichscalar);
  //--------------------------------------------------------------------------
  // get thermodynamic pressure at n+1
  //--------------------------------------------------------------------------
  thermpressaf_ = thermpressnp;


  return;

} // TimIntLoma::SetTimeLomaFields

// -------------------------------------------------------------------
// set loma parameters                               rasthofer 03/2012
// -------------------------------------------------------------------
void FLD::TimIntLoma::SetElementCustomParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_loma_parameter);

  // set parameters to update material with subgrid-scale temperature
  // potential inclusion of additional subgrid-scale terms in continuity equation
  eleparams.sublist("LOMA") = params_->sublist("LOMA");
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
| print info about turbulence model (loma-specific)            bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::PrintTurbulenceModel()
{
  FluidImplicitTimeInt::PrintTurbulenceModel();

  if (physicaltype_==INPAR::FLUID::loma and turbmodel_ == INPAR::FLUID::smagorinsky)
  {
    if (DRT::INPUT::IntegralValue<int>(params_->sublist("SUBGRID VISCOSITY"),"C_INCLUDE_CI"))
    {
      if (params_->sublist("SUBGRID VISCOSITY").get<double>("C_YOSHIZAWA")>0.0)
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
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1",thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in ApplyNonlinearBoundaryConditions    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetCustomEleParamsApplyNonlinearBoundaryConditions(Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  return;
}

/*----------------------------------------------------------------------*
| set custom ele params in LinearRelaxationSolve               bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntLoma::SetCustomEleParamsLinearRelaxationSolve(Teuchos::ParameterList& eleparams)
{
  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1",thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);
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
  const double eosfac = thermpressaf_/gasconstant_;
  statisticsmanager_->DoTimeSample(step_,eosfac,
                                   thermpressaf_,thermpressam_,
                                   thermpressdtaf_,thermpressdtam_);
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

  //necessary here, because some application time integrations add something to the residual
  //before the Neumann loads are added
  residual_->PutScalar(0.0);

  eleparams.set("thermpress at n+alpha_F/n+1",thermpressaf_);
  eleparams.set("thermpress at n+alpha_M/n",thermpressam_);
  eleparams.set("thermpressderiv at n+alpha_F/n+1",thermpressdtaf_);
  eleparams.set("thermpressderiv at n+alpha_M/n+1",thermpressdtam_);

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
}// TimIntLoma::AVM3Preparation
