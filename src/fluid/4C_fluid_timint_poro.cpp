/*-----------------------------------------------------------*/
/*! \file

\brief base class of time integration schemes for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint_poro.hpp"

#include "4C_fluid_ele_action.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

FLD::TimIntPoro::TimIntPoro(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid)
{
}

void FLD::TimIntPoro::Init()
{
  Teuchos::ParameterList* stabparams;
  stabparams = &(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if (stabparams->get<std::string>("STABTYPE") == "residual_based")
  {
    if (stabparams->get<std::string>("TDS") == "time_dependent")
    {
      FOUR_C_THROW(
          "TDS is not implemented for Poro yet. An error will occur in "
          "FluidImplicitTimeInt::TimeUpdate().");
    }
  }

  if (not alefluid_) FOUR_C_THROW("poro fluid has to be an ale fluid!");

  // set some poro-specific parameters
  SetElementCustomParameter();
}

void FLD::TimIntPoro::SetElementGeneralFluidParameter()
{
  // set some poro-specific parameters only in specific poro cases
  FluidImplicitTimeInt::SetElementGeneralFluidParameter();
}

void FLD::TimIntPoro::SetElementTurbulenceParameters()
{
  // set some poro-specific parameters only in specific poro cases
  FluidImplicitTimeInt::SetElementTurbulenceParameters();
}

void FLD::TimIntPoro::AssembleMatAndRHS()
{
  FluidImplicitTimeInt::AssembleMatAndRHS();
  PoroIntUpdate();
}

void FLD::TimIntPoro::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_, GLOBAL::Problem::Instance()->InputControlFile(), step);
  reader.ReadVector(gridv_, "gridv");
}

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
  eleparams.set<bool>("convective term", params_->get<bool>("convective term"));

  // parameter for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") =
      params_->sublist("RESIDUAL-BASED STABILIZATION");
  eleparams.sublist("EDGE-BASED STABILIZATION") = params_->sublist("EDGE-BASED STABILIZATION");

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

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
        double initialval = GLOBAL::Problem::Instance()
                                ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                                .Evaluate(lnode->X().data(), time_, 0);

        // check whether there are invalid values of porosity
        if (initialval < 1e-15) FOUR_C_THROW("zero or negative initial porosity");
        if (initialval > 1.0) FOUR_C_THROW("initial porosity greater than 1");
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          int err = init_porosity_field_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }

      break;
    }
    default:
      FOUR_C_THROW("Unknown option for initial field: %d", init);
      break;
  }
}

void FLD::TimIntPoro::UpdateIterIncrementally(
    Teuchos::RCP<const Epetra_Vector> vel)  //!< input residual velocities

{
  FluidImplicitTimeInt::UpdateIterIncrementally(vel);
  // set the new solution we just got
  if (vel != Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.
    Teuchos::RCP<Epetra_Vector> aux = CORE::LINALG::CreateVector(*(discret_->DofRowMap(0)), true);

    // only one step theta
    // new end-point accelerations
    aux->Update(1.0 / (theta_ * dta_), *velnp_, -1.0 / (theta_ * dta_), *(*veln_)(0), 0.0);
    aux->Update(-(1.0 - theta_) / theta_, *(*accn_)(0), 1.0);
    // put only to free/non-DBC DOFs
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(accnp_), aux);
    *accnp_ = *aux;
  }
}

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
}

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
}

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
    eleparams.set<POROELAST::Coupltype>("coupling", POROELAST::fluidfluid);
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
    eleparams.set<POROELAST::Coupltype>("coupling", POROELAST::fluidfluid);
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
}

void FLD::TimIntPoro::TimIntCalculateAcceleration()
{
  // for poro problems, there is a time derivative of the porosity/pressure
  // in the continuity equation. Therefore, we potentially need time
  // derivatives of the pressure and thus do not split the state vectors
  CalculateAcceleration(velnp_, veln_, velnm_, accn_, accnp_);
}

FOUR_C_NAMESPACE_CLOSE
