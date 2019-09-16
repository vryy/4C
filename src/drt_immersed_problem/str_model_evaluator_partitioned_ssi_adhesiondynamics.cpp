/*----------------------------------------------------------------------*/
/*! \file

\brief Model evaluator for partitioned ssi adhesion dynamics

\level 3

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/


#include "str_model_evaluator_partitioned_ssi_adhesiondynamics.H"

#include "../drt_structure_new/str_timint_implicit.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "Epetra_Comm.h"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::PartitionedSSIAdhesionDynamics()
    : PartitionedSSI(Teuchos::null),
      adhesion_dynamics_(false),
      penalty_(0.0),
      num_bound_species_(-1),
      adhesion_force_ptr_(Teuchos::null),
      stiff_adhesion_force_ptr_(Teuchos::null),
      disnp_ptr_(Teuchos::null),
      fixpoint_coord_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::Setup()
{
  CheckInit();

  // call setup of base class
  STR::MODELEVALUATOR::PartitionedSSI::Setup();

  //-------- check if adhesion dynamics is turned on and setup pointers if yes
  // get global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // get parameter list for cell migration
  const Teuchos::ParameterList& cellmigrationparams = problem->CellMigrationParams();

  if (cellmigrationparams.get<std::string>("ADHESION_DYNAMICS") == "yes" and
      cellmigrationparams.sublist("ADHESION MODULE").get<std::string>("COUPMETHOD") == "Penalty")
    adhesion_dynamics_ = true;

  if (adhesion_dynamics_)
  {
    // setup the pointers for displacement and stiffness
    disnp_ptr_ = GState().GetMutableDisNp();
    stiff_adhesion_force_ptr_ =
        Teuchos::rcp(new LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));

    adhesion_force_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(), true));

    // set penalty parameter
    penalty_ = cellmigrationparams.sublist("ADHESION MODULE").get<double>("PENALTY");

    // set numdof of bound integrin
    num_bound_species_ =
        problem->CellMigrationParams().sublist("ADHESION MODULE").get<int>("NUM_BOUNDSPECIES");
  }

  // set flag
  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  if (adhesion_dynamics_)
  {
    // update the structural displacement vector
    disnp_ptr_ = GState().GetDisNp();

    // Zero out force and stiffness contributions
    adhesion_force_ptr_->PutScalar(0.0);
    stiff_adhesion_force_ptr_->Zero();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::GetCurrentSolutionPtr() const
{
  CheckInit();
  return GState().GetDisNp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::GetLastTimeStepSolutionPtr() const
{
  CheckInit();
  return GState().GetDisN();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::EvaluateForce()
{
  CheckInitSetup();

  // --- evaluate adhesion contributions -------------------------------
  if (adhesion_dynamics_)
  {
    if (fixpoint_coord_ptr_ == Teuchos::null)
    {
      // do nothing because in the predictor adhesion dynamic specific
      // vectors are not yet initialized.
    }
    else
    {
      Discret().ClearState();
      Discret().SetState(0, "adh_fixpoint_coords", fixpoint_coord_ptr_);

      // set states
      Discret().SetState(0, "displacement", disnp_ptr_);

      // add action for bond traction evaluation
      Teuchos::ParameterList params;
      params.set<std::string>("action", "calc_cell_adhesion_forces");

      // add penalty parameter
      params.set<double>("penalty_parameter", penalty_);

      // add number of different bound species
      params.set<int>("NUM_BOUNDSPECIES", num_bound_species_);

      // add type of evaluation
      params.set<std::string>("eval_type", "EvaluateForce");

      // evaluate least squares traction
      Discret().EvaluateCondition(params, stiff_adhesion_force_ptr_, Teuchos::null,
          adhesion_force_ptr_, Teuchos::null, Teuchos::null, "CellFocalAdhesion", -1);
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::EvaluateStiff()
{
  CheckInitSetup();

  // --- evaluate adhesion contributions -----------------------------
  if (adhesion_dynamics_)
  {
    // set states
    Discret().SetState(0, "displacement", disnp_ptr_);

    // add action for bond traction evaluation
    Teuchos::ParameterList params;
    params.set<std::string>("action", "calc_cell_adhesion_forces");

    // add penalty parameter
    Discret().SetState(0, "adh_fixpoint_coords", fixpoint_coord_ptr_);

    // add penalty parameter
    params.set<double>("penalty_parameter", penalty_);

    // add number of different bound species
    params.set<int>("NUM_BOUNDSPECIES", num_bound_species_);

    // add type of evaluation
    params.set<std::string>("eval_type", "EvaluateStiff");

    // evaluate least squares traction
    Discret().EvaluateCondition(params, stiff_adhesion_force_ptr_, Teuchos::null,
        adhesion_force_ptr_, Teuchos::null, Teuchos::null, "CellFocalAdhesion", -1);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::EvaluateForceStiff()
{
  CheckInitSetup();

  // --- evaluate adhesion contributions -----------------------------
  if (adhesion_dynamics_)
  {
    // set states
    Discret().SetState(0, "displacement", disnp_ptr_);

    // add action for bond traction evaluation
    Teuchos::ParameterList params;
    params.set<std::string>("action", "calc_cell_adhesion_forces");

    // add penalty parameter
    Discret().SetState(0, "adh_fixpoint_coords", fixpoint_coord_ptr_);

    // add penalty parameter
    params.set<double>("penalty_parameter", penalty_);

    // add number of different bound species
    params.set<int>("NUM_BOUNDSPECIES", num_bound_species_);

    // add type of evaluation
    params.set<std::string>("eval_type", "EvaluateForceStiff");

    // evaluate least squares traction
    Discret().EvaluateCondition(params, stiff_adhesion_force_ptr_, Teuchos::null,
        adhesion_force_ptr_, Teuchos::null, Teuchos::null, "CellFocalAdhesion", -1);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  if (adhesion_dynamics_) LINALG::AssembleMyVector(1.0, f, timefac_np, *adhesion_force_ptr_);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  if (adhesion_dynamics_)
  {
    stiff_adhesion_force_ptr_->Complete();

    Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr = GState().ExtractDisplBlock(jac);
    jac_dd_ptr->Add(*stiff_adhesion_force_ptr_, false, timefac_np, 1.0);
    // no need to keep it
    stiff_adhesion_force_ptr_->Zero();
  }

  return true;
}
