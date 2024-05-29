/*----------------------------------------------------------------------*/
/*! \file
 \brief time integration schemes for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_scatra_timint_poromulti.hpp"

#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_poromultiphase_scatra_utils.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                             vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMulti::ScaTraTimIntPoroMulti(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output), L2_projection_(false)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::Init() { return; }

/*----------------------------------------------------------------------*
 | set solution fields on given dof sets                    vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::set_l2_flux_of_multi_fluid(
    Teuchos::RCP<const Epetra_MultiVector> multiflux)
{
  // set L2-projection to true
  L2_projection_ = true;

  // safety check
  if (NdsVel() >= discret_->NumDofSets()) FOUR_C_THROW("Too few dofsets on scatra discretization!");

  if (multiflux->NumVectors() % nsd_ != 0)
    FOUR_C_THROW("Unexpected length of flux vector: %i", multiflux->NumVectors());

  const int totalnumdof = multiflux->NumVectors() / nsd_;

  std::string stateprefix = "flux";

  for (int curphase = 0; curphase < totalnumdof; ++curphase)
  {
    // initialize velocity vectors
    Teuchos::RCP<Epetra_Vector> phaseflux =
        CORE::LINALG::CreateVector(*discret_->dof_row_map(NdsVel()), true);

    std::stringstream statename;
    statename << stateprefix << curphase;

    // loop all nodes on the processor
    for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
    {
      // get the processor local node
      CORE::Nodes::Node* lnode = discret_->lRowNode(lnodeid);

      // get dofs associated with current node
      std::vector<int> nodedofs = discret_->Dof(NdsVel(), lnode);

      if ((int)nodedofs.size() != nsd_)
        FOUR_C_THROW(
            "Expected number of DOFs to be equal to the number of space dimensions for flux "
            "state!");

      for (int index = 0; index < nsd_; ++index)
      {
        // get global and local dof IDs
        const int gid = nodedofs[index];
        const int lid = phaseflux->Map().LID(gid);
        if (lid < 0) FOUR_C_THROW("Local ID not found in map for given global ID!");

        const double value = (*(*multiflux)(curphase * nsd_ + index))[lnodeid];

        int err = phaseflux->ReplaceMyValue(lid, 0, value);
        if (err != 0) FOUR_C_THROW("error while inserting a value into convel");
      }
    }

    // provide scatra discretization with convective velocity
    discret_->set_state(NdsVel(), statename.str(), phaseflux);
  }
}  // ScaTraTimIntImpl::SetSolutionFields

/*----------------------------------------------------------------------*
 | set solution fields on given dof sets              kremheller  07/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::set_solution_field_of_multi_fluid(
    Teuchos::RCP<const Epetra_Vector> phinp_fluid, Teuchos::RCP<const Epetra_Vector> phin_fluid)
{
  if (NdsPressure() >= discret_->NumDofSets())
    FOUR_C_THROW("Too few dofsets on scatra discretization!");

  // provide scatra discretization with fluid primary variable field
  discret_->set_state(NdsPressure(), "phinp_fluid", phinp_fluid);
  discret_->set_state(NdsPressure(), "phin_fluid", phin_fluid);
}

/*----------------------------------------------------------------------*
 | add parameters depending on the problem                  vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::add_problem_specific_parameters_and_vectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  // provide pressure field
  params.set<bool>("L2-projection", L2_projection_);
}

/*----------------------------------------------------------------------*
 |  write current state to BINIO                           vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::output_state()
{
  // solution
  output_->WriteVector("phinp", phinp_);

  // displacement field
  if (isale_)
  {
    const auto dispnp = discret_->GetState(NdsDisp(), "dispnp");
    if (dispnp == Teuchos::null)
      FOUR_C_THROW("Cannot extract displacement field from discretization");

    const auto dispnp_multi = convert_dof_vector_to_componentwise_node_vector(dispnp, NdsDisp());
    output_->WriteVector("dispnp", dispnp_multi, IO::nodevector);
  }

  if (has_external_force_)
  {
    const int nds_vel = NdsVel();

    const auto external_force = discret_->GetState(nds_vel, "external_force");
    const auto output_external_force =
        convert_dof_vector_to_componentwise_node_vector(external_force, NdsVel());
    output_->WriteVector("external_force", output_external_force, IO::nodevector);

    const auto mobility = discret_->GetState(nds_vel, "intrinsic_mobility");
    const auto output_intrinsic_mobility =
        convert_dof_vector_to_componentwise_node_vector(mobility, NdsVel());
    output_->WriteVector("intrinsic_mobility", output_intrinsic_mobility, IO::nodevector);
  }
}  // ScaTraTimIntImpl::output_state

/*----------------------------------------------------------------------*
 | problem specific output                             kremheller 10/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::output_problem_specific()
{
  // oxygen partial pressure (if desired)
  output_oxygen_partial_pressure();

  return;
}

/*----------------------------------------------------------------------*
 | output of oxygen partial pressure                   kremheller 10/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::output_oxygen_partial_pressure()
{
  // extract conditions for oxygen partial pressure
  std::vector<CORE::Conditions::Condition*> conditions;
  discret_->GetCondition("PoroMultiphaseScatraOxyPartPressCalcCond", conditions);

  // perform all following operations only if there is at least one condition for oxygen partial
  // pressure
  if (conditions.size() > 0)
  {
    const Teuchos::RCP<Epetra_Vector> oxypartpress =
        Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(), true));

    // this condition is supposed to be for output of oxygen partial pressure over whole domain
    // it does not make sense to have more than one condition
    if (conditions.size() != 1)
      FOUR_C_THROW(
          "Should have only one PoroMultiphaseScatraOxyPartPressCalcCond per discretization");

    // extract nodal cloud from condition
    const std::vector<int>* nodegids = conditions[0]->GetNodes();

    // output
    double Pb = 0.0;

    // read input from condition
    const auto oxyscalar = conditions[0]->parameters().Get<int>("SCALARID") - 1;
    const auto CaO2_max = conditions[0]->parameters().Get<double>("CaO2_max");
    const auto Pb50 = conditions[0]->parameters().Get<double>("Pb50");
    const auto n = conditions[0]->parameters().Get<double>("n");
    const auto alpha_eff = conditions[0]->parameters().Get<double>("alpha_bl_eff");
    const auto rho_oxy = conditions[0]->parameters().Get<double>("rho_oxy");
    const auto rho_bl = conditions[0]->parameters().Get<double>("rho_bl");

    // loop over all nodes
    for (unsigned inode = 0; inode < nodegids->size(); ++inode)
    {
      // extract global ID of current node
      const int nodegid((*nodegids)[inode]);
      // process only nodes stored by current processor
      if (discret_->HaveGlobalNode(nodegid))
      {
        // extract current node
        const CORE::Nodes::Node* const node = discret_->gNode(nodegid);

        // process only nodes owned by current processor
        if (node->Owner() == discret_->Comm().MyPID())
        {
          // get dof
          int myoxydof = discret_->Dof(0, node, oxyscalar);
          const int lidoxydof = discret_->dof_row_map()->LID(myoxydof);
          if (lidoxydof < 0) FOUR_C_THROW("Couldn't extract local ID of oxygen dof!");
          // compute CaO2
          const double CaO2 = (*phinp_)[lidoxydof] * rho_bl / rho_oxy;
          // compute Pb
          POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
              Pb, CaO2, CaO2_max, Pb50, n, alpha_eff);
          // replace value
          oxypartpress->ReplaceGlobalValue(node->Id(), 0, Pb);
        }
      }
    }
    output_->WriteVector("oxypartpress", oxypartpress, IO::nodevector);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiOST::ScaTraTimIntPoroMultiOST(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}



/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiOST::Update()
{
  TimIntOneStepTheta::Update();
  ScaTraTimIntPoroMulti::Update();

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiBDF2::ScaTraTimIntPoroMultiBDF2(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntBDF2(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiBDF2::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntBDF2::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}



/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiBDF2::Update()
{
  TimIntBDF2::Update();
  ScaTraTimIntPoroMulti::Update();

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiGenAlpha::ScaTraTimIntPoroMultiGenAlpha(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}



/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiGenAlpha::Update()
{
  TimIntGenAlpha::Update();
  ScaTraTimIntPoroMulti::Update();

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiStationary::ScaTraTimIntPoroMultiStationary(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntStationary(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiStationary::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}



/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                         vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiStationary::Update()
{
  TimIntStationary::Update();
  ScaTraTimIntPoroMulti::Update();

  return;
}

FOUR_C_NAMESPACE_CLOSE
