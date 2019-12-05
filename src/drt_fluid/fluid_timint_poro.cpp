/*-----------------------------------------------------------*/
/*! \file

\brief base class of time integration schemes for porous fluid

\maintainer Johannes Kremheller

\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_poro.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"

#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntPoro::TimIntPoro(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntPoro::Init()
{
  Teuchos::ParameterList* stabparams;
  stabparams = &(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if (stabparams->get<std::string>("STABTYPE") == "residual_based")
    if (stabparams->get<std::string>("TDS") == "time_dependent")
      dserror(
          "TDS is not implemented for Poro yet. An error will occur in "
          "FluidImplicitTimeInt::TimeUpdate().");

  if (not alefluid_) dserror("poro fluid has to be an ale fluid!");

  // set some poro-specific parameters
  SetElementCustomParameter();
  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntPoro::~TimIntPoro() { return; }

/*----------------------------------------------------------------------*
| set params in constructor                                    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetElementGeneralFluidParameter()
{
  // set some poro-specific parameters only in specific poro cases
  FluidImplicitTimeInt::SetElementGeneralFluidParameter();
  return;
}

/*----------------------------------------------------------------------*
| set params in constructor                                    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetElementTurbulenceParameters()
{
  // set some poro-specific parameters only in specific poro cases
  FluidImplicitTimeInt::SetElementTurbulenceParameters();
  return;
}

void FLD::TimIntPoro::AssembleMatAndRHS()
{
  FluidImplicitTimeInt::AssembleMatAndRHS();
  PoroIntUpdate();

  return;
}

void FLD::TimIntPoro::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_, step);
  reader.ReadVector(gridv_, "gridv");

  return;
}

// -------------------------------------------------------------------
// set poro parameters                               vuong  11/2012
// -------------------------------------------------------------------
void FLD::TimIntPoro::SetElementCustomParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_poro_parameter);

  // set general element parameters
  eleparams.set("form of convective term", convform_);
  eleparams.set<int>("Linearisation", newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // set poro specific element parameters
  eleparams.set<bool>("conti partial integration", params_->get<bool>("conti partial integration"));
  eleparams.set<int>("Transient Terms Poro Fluid", params_->get<int>("Transient Terms Poro Fluid"));

  // parameter for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") =
      params_->sublist("RESIDUAL-BASED STABILIZATION");
  eleparams.sublist("EDGE-BASED STABILIZATION") = params_->sublist("EDGE-BASED STABILIZATION");

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}

/*----------------------------------------------------------------------*
 |  set initial field for porosity                                      |
 *----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetInitialPorosityField(
    const INPAR::POROELAST::InitialField init, const int startfuncno)
{
  std::cout << "FLD::TimIntPoro::SetInitialPorosityField()" << std::endl;

  switch (init)
  {
    case INPAR::POROELAST::initfield_field_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(lnode);

        int numdofs = nodedofset.size();
        double initialval =
            DRT::Problem::Instance()->Funct(startfuncno - 1).Evaluate(0, lnode->X(), time_);

        // check whether there are invalid values of porosity
        if (initialval < EPS15) dserror("zero or negative initial porosity");
        if (initialval > 1.0) dserror("initial porosity greater than 1");
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          int err = initporosityfield_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) dserror("dof not on proc");
        }
      }

      break;
    }
    default:
      dserror("Unknown option for initial field: %d", init);
      break;
  }  // switch(init)

  return;
}  // TimIntPoro::SetInitialPorosityField

/*----------------------------------------------------------------------*
| add some functionality to UpdateIterIncrementally            bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::UpdateIterIncrementally(
    Teuchos::RCP<const Epetra_Vector> vel)  //!< input residual velocities

{
  FluidImplicitTimeInt::UpdateIterIncrementally(vel);
  // set the new solution we just got
  if (vel != Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(*(discret_->DofRowMap(0)), true);

    // only one step theta
    // new end-point accelerations
    aux->Update(1.0 / (theta_ * dta_), *velnp_, -1.0 / (theta_ * dta_), *(*veln_)(0), 0.0);
    aux->Update(-(1.0 - theta_) / theta_, *(*accn_)(0), 1.0);
    // put only to free/non-DBC DOFs
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(accnp_), aux);
    *accnp_ = *aux;
  }

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 | overloading function                                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntPoro::Output()
{
  FluidImplicitTimeInt::Output();
  // output of solution
  if (step_ % upres_ == 0)
  {
    Teuchos::RCP<Epetra_Vector> convel = Teuchos::rcp(new Epetra_Vector(*velnp_));
    convel->Update(-1.0, *gridv_, 1.0);
    output_->WriteVector("convel", convel);
    output_->WriteVector("gridv", gridv_);
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ > 0 && step_ % uprestart_ == 0)
  {
    output_->WriteVector("gridv", gridv_);
  }
  return;
}  // TimIntPoro::Output

/*----------------------------------------------------------------------*
| set params in constructor                                    bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  eleparams.set<int>("Physical Type", physicaltype_);

  // just for poroelasticity
  discret_->SetState("dispn", dispn_);
  discret_->SetState("accnp", accnp_);
  discret_->SetState("accn", accn_);
  discret_->SetState("gridvn", gridvn_);

  eleparams.set("total time", time_);
  eleparams.set("delta time", dta_);

  return;
}

/*----------------------------------------------------------------------*
| update sysmat after AssembleMatAndRHS                        bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::PoroIntUpdate()
{
  sysmat_->UnComplete();

  std::string condname = "PoroPartInt";
  std::vector<DRT::Condition*> poroPartInt;
  discret_->GetCondition(condname, poroPartInt);
  if (poroPartInt.size())
  {
    Teuchos::ParameterList eleparams;

    // set action for elements
    eleparams.set<int>("action", FLD::poro_boundary);
    eleparams.set("total time", time_);
    eleparams.set("delta time", dta_);
    eleparams.set<POROELAST::coupltype>("coupling", POROELAST::fluidfluid);
    eleparams.set<int>("Physical Type", physicaltype_);

    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
    discret_->SetState("velnp", velnp_);
    discret_->SetState("scaaf", scaaf_);
    discret_->EvaluateCondition(
        eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null, condname);
    discret_->ClearState();
  }

  condname = "PoroPresInt";
  std::vector<DRT::Condition*> poroPresInt;
  discret_->GetCondition(condname, poroPresInt);
  if (poroPresInt.size())
  {
    Teuchos::ParameterList eleparams;

    // set action for elements
    eleparams.set<int>("action", FLD::poro_prescoupl);
    eleparams.set<POROELAST::coupltype>("coupling", POROELAST::fluidfluid);
    eleparams.set<int>("Physical Type", physicaltype_);

    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
    discret_->SetState("gridv", gridv_);
    discret_->SetState("velnp", velnp_);
    discret_->EvaluateCondition(
        eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null, condname);
    discret_->ClearState();
  }
  sysmat_->Complete();

  return;
}

/*----------------------------------------------------------------------*
| Calculate acceleration for poro                              bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntPoro::TimIntCalculateAcceleration()
{
  // for poro problems, there is a time derivative of the porosity/pressure
  // in the continuity equation. Therefore, we potentially need time
  // derivatives of the pressure and thus do not split the state vectors
  CalculateAcceleration(velnp_, veln_, velnm_, accn_, accnp_);

  return;
}
