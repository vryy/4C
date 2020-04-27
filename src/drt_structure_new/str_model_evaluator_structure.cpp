/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all structure terms


\level 3

*/
/*-----------------------------------------------------------*/

#include "str_model_evaluator_structure.H"
#include "str_model_evaluator_data.H"
#include "str_timint_implicit.H"
#include "str_predict_generic.H"
#include "str_utils.H"
#include "str_integrator.H"
#include "str_dbc.H"
#include "str_timint_basedataio_runtime_vtk_output.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_discret.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/discretization_runtime_vtu_writer.H"

#include "str_discretization_runtime_vtu_output_params.H"
#include "../drt_beam3/beam_discretization_runtime_vtu_writer.H"
#include "../drt_beam3/beam_discretization_runtime_vtu_output_params.H"

#include "../drt_lib/prestress_service.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Structure::Structure()
    : dt_ele_ptr_(NULL),
      masslin_type_(INPAR::STR::ml_none),
      stiff_ptr_(NULL),
      stiff_ptc_ptr_(Teuchos::null),
      dis_incr_ptr_(Teuchos::null),
      vtu_writer_ptr_(Teuchos::null),
      beam_vtu_writer_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // get the global state content
  {
    // structural element evaluation time
    dt_ele_ptr_ = &(GState().GetMutableElementEvaluationTime());
  }

  // displ-displ block
  stiff_ptr_ = dynamic_cast<LINALG::SparseMatrix*>(GState().CreateStructuralStiffnessMatrixBlock());

  // modified stiffness pointer for storing element based scaling operator (PTC)
  stiff_ptc_ptr_ =
      Teuchos::rcp(new LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));

  if (stiff_ptr_ == NULL) dserror("Dynamic cast to LINALG::SparseMatrix failed!");

  // get the structural dynamic content
  {
    // setup important evaluation booleans
    masslin_type_ = TimInt().GetDataSDyn().GetMassLinType();
  }
  // setup new variables
  {
    dis_incr_ptr_ = Teuchos::rcp(new Epetra_Vector(DisNp().Map(), true));
  }



  if (GInOutput().GetRuntimeVtkOutputParams() != Teuchos::null)
  {
    if (GInOutput().GetRuntimeVtkOutputParams()->OutputStructure()) InitOutputRuntimeVtkStructure();

    if (GInOutput().GetRuntimeVtkOutputParams()->OutputBeams()) InitOutputRuntimeVtkBeams();
  }

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  /* --- reset external forces
   * Please note, that PutScalar is safer (but maybe slower) than
   * Scale(0.0), because of possible NaN and inf values! */
  FextNp().PutScalar(0.0);

  /* --- reset internal forces
   * Please note, that PutScalar is safer (but maybe slower) than
   * Scale(0.0), because of possible NaN and inf values! */
  FintNp().PutScalar(0.0);

  // reset stiffness matrix
  Stiff().Zero();

  // reset modified stiffness matrix
  StiffPTC().Zero();


  // set evaluation time back to zero
  *dt_ele_ptr_ = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::EvaluateForce()
{
  CheckInitSetup();
  bool ok = true;
  // ---------------------------------------
  // (1) EXTERNAL FORCES
  // ---------------------------------------
  ok = ApplyForceExternal();

  // ---------------------------------------
  // (2) INTERNAL FORCES
  // ---------------------------------------
  // ordinary internal force
  ok = (ok ? ApplyForceInternal() : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::EvaluateStiff()
{
  CheckInitSetup();
  bool ok = true;

  /* We use the same routines as for the ApplyForceStiff case, but we
   * do not update the global force vector, which is used for the
   * solution process in the NOX library.
   * This is meaningful, since the computational overhead, which is
   * generated by evaluating the right hand side is negligible */
  // *********** time measurement ***********
  double dtcpu = GState().GetTimer()->WallTime();
  // *********** time measurement ***********
  // ---------------------------------------------------------------------
  // (1) EXTRERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  ok = ApplyForceStiffExternal();

  // ---------------------------------------------------------------------
  // (2) INTERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  // ordinary internal force
  ok = (ok ? ApplyForceStiffInternal() : false);

  // *********** time measurement ***********
  *dt_ele_ptr_ += GState().GetTimer()->WallTime() - dtcpu;
  // *********** time measurement ***********

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::EvaluateForceStiff()
{
  CheckInitSetup();
  bool ok = true;

  // *********** time measurement ***********
  double dtcpu = GState().GetTimer()->WallTime();
  // *********** time measurement ***********
  // ---------------------------------------------------------------------
  // (1) EXTRERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  ok = ApplyForceStiffExternal();

  // ---------------------------------------------------------------------
  // (2) INTERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  // ordinary internal force
  ok = (ok ? ApplyForceStiffInternal() : false);

  // *********** time measurement ***********
  *dt_ele_ptr_ += GState().GetTimer()->WallTime() - dtcpu;
  // *********** time measurement ***********

  // that's it
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::AssembleForce(Epetra_Vector& f, const double& timefac_np) const
{
  LINALG::AssembleMyVector(1.0, f, -timefac_np, FextNp());
  LINALG::AssembleMyVector(1.0, f, timefac_np, FintNp());

  // add the scaled force contributions of the old time step
  // structural dofs of the right-hand-side vector at t_{n+timefac_n} (read-only)
  Teuchos::RCP<const Epetra_Vector> fstructold_ptr = GState().GetFstructureOld();
  LINALG::AssembleMyVector(1.0, f, 1.0, *fstructold_ptr);

  // add the visco and mass contributions
  Int().AddViscoMassContributions(f);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  int err = Stiff().Scale(timefac_np);
  GState().AssignModelBlock(jac, Stiff(), Type(), DRT::UTILS::block_displ_displ);

  // add the visco and mass contributions
  Int().AddViscoMassContributions(jac);

  return (err == 0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::InitializeInertiaAndDamping()
{
  CheckInitSetup();

  // currently a fixed number of matrix and vector pointers are supported
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  // create vector with zero entries
  Teuchos::RCP<const Epetra_Vector> zeros = Int().GetDbc().GetZerosPtr();

  // set vector values needed by elements
  // --> initially zero !!!
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", zeros);
  Discret().SetState(0, "displacement", zeros);

  // set action type and evaluation matrix and vector pointers
  StaticContributions(&eval_mat[0], &eval_vec[0]);
  MaterialDampingContributions(&eval_mat[0]);
  InertialContributions(&eval_mat[0], &eval_vec[0]);

  // evaluate
  EvaluateInternal(&eval_mat[0], &eval_vec[0]);

  // complete stiffness and mass matrix
  FillComplete();

  // assemble the rayleigh damping matrix
  RayleighDampingMatrix();

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceInternal()
{
  CheckInitSetup();

  // currently a fixed number of matrix and vector pointers are supported
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Discret().SetState(0, "velocity", GState().GetVelNp());

  // set action type and evaluation matrix and vector pointers
  StaticContributions(&eval_vec[0]);
  MaterialDampingContributions(&eval_mat[0]);
  InertialContributions(&eval_vec[0]);

  // evaluate ...
  EvaluateInternal(&eval_mat[0], &eval_vec[0]);

  // evaluate inertia and visco forces
  InertialAndViscousForces();

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceExternal()
{
  CheckInitSetup();

  // Set to default value, because it is unnecessary for the
  // EvaluateNeumann routine.
  EvalData().SetActionType(DRT::ELEMENTS::none);
  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisN());
  if (EvalData().GetDampingType() == INPAR::STR::damp_material)
    Discret().SetState(0, "velocity", GState().GetVelN());
  Discret().SetState(0, "displacement new", GState().GetDisNp());
  EvaluateNeumann(GState().GetMutableFextNp(), Teuchos::null);

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffExternal()
{
  CheckInitSetup();

  if (PreApplyForceStiffExternal(FextNp(), *stiff_ptr_)) return true;

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisN());

  if (EvalData().GetDampingType() == INPAR::STR::damp_material)
    Discret().SetState(0, "velocity", GState().GetVelN());

  // get load vector
  if (!TimInt().GetDataSDyn().GetLoadLin())
    EvaluateNeumann(GState().GetMutableFextNp(), Teuchos::null);
  else
  {
    Discret().SetState(0, "displacement new", GState().GetDisNp());
    /* Add the linearization of the external force to the stiffness
     * matrix. */
    EvaluateNeumann(GState().GetMutableFextNp(), Teuchos::rcpFromRef(*stiff_ptr_));
  }

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::PreApplyForceStiffExternal(
    Epetra_Vector& fextnp, LINALG::SparseMatrix& stiff) const
{
  CheckInitSetup();

  const STR::TIMINT::Implicit* impl_ptr = dynamic_cast<const STR::TIMINT::Implicit*>(&TimInt());
  if (impl_ptr) return impl_ptr->Predictor().PreApplyForceExternal(fextnp);

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffInternal()
{
  CheckInitSetup();
  // currently a fixed number of matrix and vector pointers are supported
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Discret().SetState(0, "velocity", GState().GetVelNp());

  // set action types and evaluate matrices/vectors
  StaticContributions(&eval_mat[0], &eval_vec[0]);
  MaterialDampingContributions(&eval_mat[0]);
  if (masslin_type_ != INPAR::STR::ml_none) InertialContributions(&eval_mat[0], &eval_vec[0]);

  // evaluate
  EvaluateInternal(&eval_mat[0], &eval_vec[0]);

  // complete stiffness and mass matrix
  FillComplete();

  // evaluate inertial and viscous forces
  InertialAndViscousForces();

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::StaticContributions(
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
  // set default matrix
  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptr_);
  // set default force vector
  eval_vec[0] = GState().GetMutableFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::StaticContributions(Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_internalforce);
  // set default force vector
  eval_vec[0] = GState().GetMutableFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::MaterialDampingContributions(
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat)
{
  if (EvalData().GetDampingType() != INPAR::STR::damp_material) return;

  // action for elements
  // (reset the action type to be independent of the calling order)
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
  // set the discretization state
  Discret().SetState(0, "velocity", GState().GetVelNp());
  // reset damping matrix
  Damp().Zero();
  // add the stiffness matrix as well (also for the ApplyForce case!)
  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptr_);
  // set damping matrix
  eval_mat[1] = GState().GetMutableDampMatrix();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InertialContributions(
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  CheckInitSetup();

  if (TimInt().GetDataSDynPtr()->NeglectInertia()) return;

  // overwrite element action
  if (TimInt().GetDataSDyn().IsMassLumping())
    EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstifflmass);
  else
    EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiffmass);

  // set the discretization state
  Discret().SetState(0, "velocity", GState().GetVelNp());
  Discret().SetState(0, "acceleration", GState().GetAccNp());
  // reset the mass matrix
  Mass().Zero();
  // set mass matrix
  eval_mat[1] = GState().GetMutableMassMatrix();
  // set inertial vector if necessary
  eval_vec[1] = GetInertialForce();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InertialContributions(Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  CheckInitSetup();

  if (masslin_type_ == INPAR::STR::ml_none or TimInt().GetDataSDynPtr()->NeglectInertia()) return;

  // overwrite element action
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_internalinertiaforce);
  // set the discretization state
  Discret().SetState(0, "velocity", GState().GetVelNp());
  Discret().SetState(0, "acceleration", GState().GetAccNp());

  // set inertial vector if necessary
  eval_vec[1] = GetInertialForce();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InertialAndViscousForces()
{
  CheckInitSetup();

  if (masslin_type_ == INPAR::STR::ml_none and !TimInt().GetDataSDynPtr()->NeglectInertia())
  {
    // calculate the inertial force at t_{n+1}
    Mass().Multiply(false, *GState().GetAccNp(), FinertialNp());
  }

  // calculate the viscous/damping force at t_{n+1}
  if (EvalData().GetDampingType() != INPAR::STR::damp_none)
  {
    if (not Damp().Filled()) Damp().Complete();
    Damp().Multiply(false, *GState().GetVelNp(), FviscoNp());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::FillComplete()
{
  if (not stiff_ptr_->Filled()) stiff_ptr_->Complete();

  if (not Mass().Filled()) Mass().Complete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RayleighDampingMatrix()
{
  if (EvalData().GetDampingType() != INPAR::STR::damp_rayleigh) return;

  const double& dampk = TimInt().GetDataSDyn().GetDampingStiffnessFactor();
  const double& dampm = TimInt().GetDataSDyn().GetDampingMassFactor();

  // damping matrix with initial stiffness
  Damp().Add(Stiff(), false, dampk, 0.0);
  Damp().Add(Mass(), false, dampm, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::MODELEVALUATOR::Structure::GetInertialForce()
{
  switch (masslin_type_)
  {
    case INPAR::STR::ml_rotations:
    case INPAR::STR::ml_standard:
    {
      FinertialNp().PutScalar(0.0);
      // set inertial force
      return GState().GetMutableFinertialNp();
      break;
    }
    case INPAR::STR::ml_none:
      // do nothing
      break;
    default:
      dserror("Unknown mass linearization type!");
      exit(EXIT_FAILURE);
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InitOutputRuntimeVtkStructure()
{
  CheckInit();

  vtu_writer_ptr_ = Teuchos::rcp(new DiscretizationRuntimeVtuWriter());

  // get the parameter container object
  const STR::TIMINT::ParamsRuntimeVtkOutput& vtu_output_params =
      *GInOutput().GetRuntimeVtkOutputParams();

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  if (GInOutput().GetRuntimeVtkOutputParams()->OutputEveryIteration())
    num_timesteps_in_simulation_upper_bound *= 1000;

  // initialize the writer object
  vtu_writer_ptr_->Initialize(
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(
          const_cast<STR::MODELEVALUATOR::Structure*>(this)->DiscretPtr(), true),
      num_timesteps_in_simulation_upper_bound, GState().GetTimeN(),
      vtu_output_params.WriteBinaryOutput());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteTimeStepOutputRuntimeVtkStructure() const
{
  CheckInitSetup();

  // export displacement state to column format
  const DRT::Discretization& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  LINALG::Export(*GState().GetDisN(), *disn_col);

  if (not GInOutput().GetRuntimeVtkOutputParams()->OutputEveryIteration())
    WriteOutputRuntimeVtkStructure(disn_col, GState().GetStepN(), GState().GetTimeN());
  else
    WriteOutputRuntimeVtkStructure(disn_col, 10000 * GState().GetStepN(), GState().GetTimeN());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteIterationOutputRuntimeVtkStructure() const
{
  CheckInitSetup();

  // export displacement state to column format
  const DRT::Discretization& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disnp_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  LINALG::Export(*GState().GetDisNp(), *disnp_col);

  const int augmented_timestep_number_incl_iteration_count =
      10000 * GState().GetStepN() + 1 * EvalData().GetNlnIter();

  const double augmented_time_incl_iteration_count =
      GState().GetTimeN() + 1e-8 * EvalData().GetNlnIter();

  WriteOutputRuntimeVtkStructure(disnp_col, augmented_timestep_number_incl_iteration_count,
      augmented_time_incl_iteration_count);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteOutputRuntimeVtkStructure(
    const Teuchos::RCP<Epetra_Vector>& displacement_state_vector, int timestep_number,
    double time) const
{
  CheckInitSetup();

  // get the parameter container object
  const DRT::ELEMENTS::StructureRuntimeVtuOutputParams& structure_vtu_output_params =
      *GInOutput().GetRuntimeVtkOutputParams()->GetStructureParams();

  // reset time and time step of the writer object
  vtu_writer_ptr_->ResetTimeAndTimeStep(time, timestep_number);

  // append all desired output data to the writer object's storage

  // append displacement if desired
  if (structure_vtu_output_params.OutputDisplacementState())
    vtu_writer_ptr_->AppendDofBasedResultDataVector(
        displacement_state_vector, 3, 0, "displacement");

  // append element owner if desired
  if (structure_vtu_output_params.OutputElementOwner())
    vtu_writer_ptr_->AppendElementOwner("element_owner");

  // append element GIDs if desired
  if (structure_vtu_output_params.OutputElementGID())
    vtu_writer_ptr_->AppendElementGID("element_gid");

  // finalize everything and write all required files to filesystem
  vtu_writer_ptr_->WriteFiles();
  vtu_writer_ptr_->WriteCollectionFileOfAllWrittenFiles();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InitOutputRuntimeVtkBeams()
{
  beam_vtu_writer_ptr_ = Teuchos::rcp(new BeamDiscretizationRuntimeVtuWriter());

  // get the parameter container object
  const DRT::ELEMENTS::BeamRuntimeVtuOutputParams& beam_vtu_output_params =
      *GInOutput().GetRuntimeVtkOutputParams()->GetBeamParams();

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  if (GInOutput().GetRuntimeVtkOutputParams()->OutputEveryIteration())
    num_timesteps_in_simulation_upper_bound *= 1000;

  // export displacement state to column format
  const DRT::Discretization& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  LINALG::Export(*GState().GetDisN(), *disn_col);

  // get bounding box object only if periodic boundaries are active
  Teuchos::RCP<GEO::MESHFREE::BoundingBox> bounding_box_ptr =
      TimInt().GetDataSDynPtr()->GetPeriodicBoundingBox();

  // initialize the writer object with current displacement state
  beam_vtu_writer_ptr_->Initialize(
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(
          const_cast<STR::MODELEVALUATOR::Structure*>(this)->DiscretPtr(), true),
      disn_col, beam_vtu_output_params.UseAbsolutePositions(), bounding_box_ptr,
      num_timesteps_in_simulation_upper_bound, GState().GetTimeN(),
      GInOutput().GetRuntimeVtkOutputParams()->WriteBinaryOutput());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteTimeStepOutputRuntimeVtkBeams() const
{
  CheckInitSetup();

  // export displacement state to column format
  const DRT::Discretization& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  LINALG::Export(*GState().GetDisN(), *disn_col);

  if (not GInOutput().GetRuntimeVtkOutputParams()->OutputEveryIteration())
    WriteOutputRuntimeVtkBeams(disn_col, GState().GetStepN(), GState().GetTimeN());
  else
    WriteOutputRuntimeVtkBeams(disn_col, 10000 * GState().GetStepN(), GState().GetTimeN());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteIterationOutputRuntimeVtkBeams() const
{
  CheckInitSetup();

  // export displacement state to column format
  const DRT::Discretization& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disnp_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  LINALG::Export(*GState().GetDisNp(), *disnp_col);

  const int augmented_timestep_number_incl_iteration_count =
      10000 * GState().GetStepN() + 1 * EvalData().GetNlnIter();

  const double augmented_time_incl_iteration_count =
      GState().GetTimeN() + 1e-8 * EvalData().GetNlnIter();

  WriteOutputRuntimeVtkBeams(disnp_col, augmented_timestep_number_incl_iteration_count,
      augmented_time_incl_iteration_count);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteOutputRuntimeVtkBeams(
    const Teuchos::RCP<Epetra_Vector>& displacement_state_vector, int timestep_number,
    double time) const
{
  CheckInitSetup();

  // get the parameter container object
  const DRT::ELEMENTS::BeamRuntimeVtuOutputParams& beam_vtu_output_params =
      *GInOutput().GetRuntimeVtkOutputParams()->GetBeamParams();

  // set geometry
  beam_vtu_writer_ptr_->SetGeometryFromBeamDiscretization(displacement_state_vector);

  // reset time and time step of the writer object
  beam_vtu_writer_ptr_->ResetTimeAndTimeStep(time, timestep_number);

  // append all desired output data to the writer object's storage
  beam_vtu_writer_ptr_->AppendElementOwningProcessor();

  // append beam radius
  beam_vtu_writer_ptr_->AppendElementCircularCrossSectionRadius();

  // append displacement if desired
  if (beam_vtu_output_params.IsWriteInternalEnergyElement())
    beam_vtu_writer_ptr_->AppendElementInternalEnergy();

  // append displacement if desired
  if (beam_vtu_output_params.IsWriteKineticEnergyElement())
    beam_vtu_writer_ptr_->AppendElementKineticEnergy();

  // append displacement if desired
  if (beam_vtu_output_params.OutputDisplacementState())
    beam_vtu_writer_ptr_->AppendDisplacementField(displacement_state_vector);

  // append triads if desired
  if (beam_vtu_output_params.IsWriteTriadVisualizationPoints())
    beam_vtu_writer_ptr_->AppendTriadField(displacement_state_vector);

  // append material cross-section strain resultants if desired
  if (beam_vtu_output_params.IsWriteMaterialStrainsGaussPoints())
    beam_vtu_writer_ptr_->AppendGaussPointMaterialCrossSectionStrainResultants();

  // append material cross-section stress resultants if desired
  if (beam_vtu_output_params.IsWriteMaterialStressesGaussPoints())
    beam_vtu_writer_ptr_->AppendGaussPointMaterialCrossSectionStressResultants();

  // append filament id and type if desired
  if (beam_vtu_output_params.IsWriteElementFilamentCondition())
    beam_vtu_writer_ptr_->AppendElementFilamentIdAndType();

  // append filament id and type if desired
  if (beam_vtu_output_params.IsWriteOrientationParamter())
    beam_vtu_writer_ptr_->AppendElementOrientationParamater(displacement_state_vector);

  // export displacement state to column format
  if (beam_vtu_output_params.IsWriteRVECrosssectionForces())
    beam_vtu_writer_ptr_->AppendRVECrosssectionForces(displacement_state_vector);

  // export beam element IDs
  if (beam_vtu_output_params.IsWriteElementGID()) beam_vtu_writer_ptr_->AppendElementGID();

  // finalize everything and write all required VTU files to filesystem
  beam_vtu_writer_ptr_->WriteFiles();

  // write collection file
  beam_vtu_writer_ptr_->WriteCollectionFileOfAllWrittenFiles();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternal(
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  PreEvaluateInternal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());

  EvaluateInternal(p, eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternal(Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  if (p.numParams() > 1)
    dserror(
        "Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  if (not p.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>(
          "interface"))
    dserror("The given parameter has the wrong type!");

  // FixMe as soon as possible: write data to the parameter list.
  // this is about to go, once the old time integration is deleted
  ParamsInterface2ParameterList(EvalDataPtr(), p);

  Discret().Evaluate(p, eval_mat[0], eval_mat[1], eval_vec[0], eval_vec[1], eval_vec[2]);
  Discret().ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternalSpecifiedElements(Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec,
    const Epetra_Map* ele_map_to_be_evaluated)
{
  if (p.numParams() > 1)
    dserror(
        "Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  if (not p.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>(
          "interface"))
    dserror("The given parameter has the wrong type!");

  // write data to the parameter list.
  // this is about to go, once the old time integration is deleted
  ParamsInterface2ParameterList(EvalDataPtr(), p);

  DRT::UTILS::Evaluate(*(Teuchos::rcp_dynamic_cast<DRT::Discretization>(DiscretPtr(), true)), p,
      *eval_mat, *eval_vec, ele_map_to_be_evaluated);

  Discret().ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateNeumann(const Teuchos::RCP<Epetra_Vector>& eval_vec,
    const Teuchos::RCP<LINALG::SparseOperator>& eval_mat)
{
  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());
  EvaluateNeumann(p, eval_vec, eval_mat);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateNeumann(Teuchos::ParameterList& p,
    const Teuchos::RCP<Epetra_Vector>& eval_vec,
    const Teuchos::RCP<LINALG::SparseOperator>& eval_mat)
{
  if (p.numParams() > 1)
    dserror(
        "Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  if (not p.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>(
          "interface"))
    dserror("The given parameter has the wrong type!");
  Discret().EvaluateNeumann(p, eval_vec, eval_mat);
  Discret().ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // write forces
  iowriter.WriteVector("fstructure_old", GState().GetFstructureOld());

  if (forced_writerestart) return;

  iowriter.WriteVector("displacement", GState().GetDisN());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  // read structural force vector
  ioreader.ReadVector(GState().GetMutableFstructureOld(), "fstructure_old");
  // read displacement field
  Teuchos::RCP<Epetra_Vector>& disnp = GState().GetMutableDisNp();
  ioreader.ReadVector(disnp, "displacement");
  GState().GetMutableMultiDis()->UpdateSteps(*disnp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Predict(const INPAR::STR::PredEnum& pred_type)
{
  // set the element action
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_predict);
  EvalData().SetPredictorType(pred_type);

  // set the matrix and vector pointers to Teuchos::null
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RunPostComputeX(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  CheckInitSetup();
  Reset(xnew);
  /* set the class internal displacement increment vector. Check if it is
   * meaningful/necessary in some cases, like incremental strains etc. */
  dis_incr_ptr_ = GState().ExtractDisplEntries(dir);
  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);
  Discret().SetState(0, "displacement", GState().GetMutableDisNp());
  // set the element action
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_recover);
  // set the matrix and vector pointers to Teuchos::null
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RunPostIterate(const NOX::Solver::Generic& solver)
{
  CheckInitSetup();

  if (vtu_writer_ptr_ != Teuchos::null and
      GInOutput().GetRuntimeVtkOutputParams()->OutputEveryIteration())
    WriteIterationOutputRuntimeVtkStructure();

  // write special output for beams if desired
  if (beam_vtu_writer_ptr_ != Teuchos::null and
      GInOutput().GetRuntimeVtkOutputParams()->OutputEveryIteration())
    WriteIterationOutputRuntimeVtkBeams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepState(const double& timefac_n)
{
  CheckInitSetup();
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  GState().GetMutableMultiDis()->UpdateSteps(DisNp());

  // new velocities at t_{n+1} -> t_{n}
  //    V_{n} := V_{n+1}
  GState().GetMutableMultiVel()->UpdateSteps(*GState().GetVelNp());

  // new at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  GState().GetMutableMultiAcc()->UpdateSteps(*GState().GetAccNp());

  // store the old external force
  GState().GetMutableFextN()->Scale(1.0, FextNp());

  // new at t_{n+1} -> t_{n+timefac_n}
  //    F^{struct}_{n+timefac_n} := timefac_n * F^{struct}_{n+1}
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetMutableFstructureOld();
  fstructold_ptr->Update(timefac_n, FintNp(), 1.0);
  fstructold_ptr->Update(-timefac_n, FextNp(), 1.0);

  // set the displacement increment back to zero
  dis_incr_ptr_->Scale(0.0);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateJacobianContributionsFromElementLevelForPTC()
{
  CheckInitSetup();
  // currently a fixed number of matrix and vector pointers are supported
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_addjacPTC);

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisNp());

  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptc_ptr_);

  // evaluate
  EvaluateInternal(&eval_mat[0], &eval_vec[0]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::AssembleJacobianContributionsFromElementLevelForPTC(
    Teuchos::RCP<LINALG::SparseMatrix>& modjac, const double& timefac_n)
{
  GState().AssignModelBlock(*modjac, StiffPTC(), Type(), DRT::UTILS::block_displ_displ);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepElement()
{
  CheckInitSetup();
  // other parameters that might be needed by the elements
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);

  const INPAR::STR::PreStress prestress_type = TimInt().GetDataSDyn().GetPreStressType();
  const double prestress_time = TimInt().GetDataSDyn().GetPreStressTime();
  bool isDuringPrestressing =
      ::UTILS::PRESTRESS::IsActive(GState().GetTimeN(), prestress_type, prestress_time);

  if (isDuringPrestressing)
  {
    switch (prestress_type)
    {
      case INPAR::STR::PreStress::mulf:
        if (Discret().Comm().MyPID() == 0) IO::cout << "====== Entering MULF update" << IO::endl;

        // Choose special update action for elements in case of MULF
        EvalData().SetActionType(DRT::ELEMENTS::struct_update_prestress);
        break;
      case INPAR::STR::PreStress::id:
        dserror(
            "Inverse design prestressing is only implemented in the old time integration "
            "framework.");
        break;
      default:
        dserror(
            "The type of prestressing algorithm is unknown in the new time "
            "integration framework!");
        break;
    }
  }
  else
  {
    // Call the normal element  update routine
    EvalData().SetActionType(DRT::ELEMENTS::struct_calc_update_istep);
  }


  // go to elements
  Discret().ClearState();
  Discret().SetState("displacement", GState().GetDisN());

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};
  EvaluateInternal(eval_mat, eval_vec);

  // Check for prestressing
  if (isDuringPrestressing)
  {
    switch (prestress_type)
    {
      case INPAR::STR::PreStress::mulf:
        // This is a MULF step, hence we do not update the displacements at the end of the timestep.
        // This is achieved by resetting the displacements, velocities and accelerations.
        GState().GetMutableDisN()->PutScalar(0.0);
        GState().GetMutableVelN()->PutScalar(0.0);
        GState().GetMutableAccN()->PutScalar(0.0);
        break;
      default:
        break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineStressStrain()
{
  CheckInitSetup();

  if (GInOutput().GetStressOutputType() == INPAR::STR::stress_none and
      GInOutput().GetCouplingStressOutputType() == INPAR::STR::stress_none and
      GInOutput().GetStrainOutputType() == INPAR::STR::strain_none and
      GInOutput().GetPlasticStrainOutputType() == INPAR::STR::strain_none)
    return;

  // set all parameters in the evaluation data container
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_stress);
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);
  EvalData().SetStressData(Teuchos::rcp(new std::vector<char>()));
  EvalData().SetCouplingStressData(Teuchos::rcp(new std::vector<char>()));
  EvalData().SetStrainData(Teuchos::rcp(new std::vector<char>()));
  EvalData().SetPlasticStrainData(Teuchos::rcp(new std::vector<char>()));

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineStrainEnergy(
    const Epetra_Vector& disnp, const bool global)
{
  CheckInitSetup();

  // set required parameters in the evaluation data container
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_energy);
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);

  // set state vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", Teuchos::rcpFromRef(disnp));
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  PreEvaluateInternal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());

  // evaluate energy contributions on element level (row elements only)
  EvaluateInternalSpecifiedElements(p, eval_mat, eval_vec, Discret().ElementRowMap());

  if (global)
  {
    double my_int_energy = EvalData().GetEnergyData(STR::internal_energy);
    double gsum = 0.0;
    Discret().Comm().SumAll(&my_int_energy, &gsum, 1);

    EvalData().SetValueForEnergyType(gsum, STR::internal_energy);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineEnergy()
{
  DetermineEnergy(*GState().GetDisNp(), GState().GetVelNp().get(), false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineEnergy(
    const Epetra_Vector& disnp, const Epetra_Vector* velnp, const bool global)
{
  DetermineStrainEnergy(disnp, global);

  // global calculation of kinetic energy
  if (masslin_type_ == INPAR::STR::ml_none and velnp != NULL)
  {
    double kinetic_energy_times2 = 0.0;

    Teuchos::RCP<Epetra_Vector> linear_momentum =
        LINALG::CreateVector(*GState().DofRowMapView(), true);

    Mass().Multiply(false, *velnp, *linear_momentum);

    linear_momentum->Dot(*velnp, &kinetic_energy_times2);

    // only add the result on one processor because we sum over all procs later
    if (global or GState().GetMyRank() == 0)
    {
      EvalData().AddContributionToEnergyType(0.5 * kinetic_energy_times2, STR::kinetic_energy);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineOptionalQuantity()
{
  CheckInitSetup();

  switch (GInOutput().GetOptQuantityOutputType())
  {
    case INPAR::STR::optquantity_none:
    {
      // do nothing and return
      return;
    }
    case INPAR::STR::optquantity_membranethickness:
    {
      // evaluate thickness of membrane finite elements
      EvalData().SetActionType(DRT::ELEMENTS::struct_calc_thickness);
      break;
    }
    default:
      dserror("Type of optional quantity not implemented yet!");
  }

  // set all parameters in the evaluation data container
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);
  EvalData().SetOptQuantityData(Teuchos::rcp(new std::vector<char>()));

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::DetermineElementVolumes(
    const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols)
{
  // set action in params-interface
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_mass_volume);

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());

  // set vector values needed by elements
  Discret().ClearState();
  Teuchos::RCP<const Epetra_Vector> disnp = GState().ExtractDisplEntries(x);
  Discret().SetState(0, "displacement", disnp);

  DRT::Discretization& discret = dynamic_cast<DRT::Discretization&>(Discret());

  // start evaluation
  const Epetra_Map* relemap = Discret().ElementRowMap();
  ele_vols = Teuchos::rcp(new Epetra_Vector(*relemap, true));
  const unsigned my_num_reles = relemap->NumMyElements();

  DRT::Element::LocationArray la(Discret().NumDofSets());
  LINALG::SerialDenseVector ele_vol(6, true);

  LINALG::SerialDenseMatrix empty_dummy_mat;
  LINALG::SerialDenseVector empty_dummy_vec;

  STR::ELEMENTS::EvalErrorFlag ele_eval_error = STR::ELEMENTS::ele_error_none;
  for (unsigned elid = 0; elid < my_num_reles; ++elid)
  {
    DRT::Element* rele = Discret().lRowElement(elid);
    rele->LocationVector(discret, la, false);

    EvalData().SetActionType(DRT::ELEMENTS::analyse_jacobian_determinant);
    rele->Evaluate(p, discret, la, empty_dummy_mat, empty_dummy_mat, ele_vol, empty_dummy_vec,
        empty_dummy_vec);

    if (not EvalData().IsEleEvalError())
    {
      EvalData().SetActionType(DRT::ELEMENTS::struct_calc_mass_volume);
      rele->Evaluate(p, discret, la, empty_dummy_mat, empty_dummy_mat, ele_vol, empty_dummy_vec,
          empty_dummy_vec);
    }

    // set a negative value, if the evaluation failed
    if (EvalData().IsEleEvalError())
    {
      ele_eval_error = EvalData().GetEleEvalErrorFlag();
      // reset for the next element
      EvalData().SetEleEvalErrorFlag(STR::ELEMENTS::ele_error_none);
      ele_vol(2) = -1.0;
    }

    const int rele_lid = relemap->LID(rele->Id());
    (*ele_vols)[rele_lid] = ele_vol(2);
    ele_vol.Zero();
  }

  Discret().ClearState();
  EvalData().SetEleEvalErrorFlag(ele_eval_error);

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

  // write output every iteration for debug purposes
  if (GInOutput().IsOutputEveryIter())
  {
    iowriter.WriteVector("displacement", GState().GetDisNp());
    /* for visualization of vel and acc do not forget to comment in
     * corresponding lines in StructureEnsightWriter */
    if (GInOutput().IsWriteVelAcc())
    {
      iowriter.WriteVector("velocity", GState().GetVelNp());
      iowriter.WriteVector("acceleration", GState().GetAccNp());
    }
    return;
  }
  else
  {
    // write default output...
    iowriter.WriteVector("displacement", GState().GetDisN());

    /* for visualization of vel and acc do not forget to comment in
     * corresponding lines in StructureEnsightWriter */
    if (GInOutput().IsWriteVelAcc())
    {
      iowriter.WriteVector("velocity", GState().GetVelN());
      iowriter.WriteVector("acceleration", GState().GetAccN());
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RuntimePreOutputStepState() { CheckInitSetup(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RuntimeOutputStepState() const
{
  CheckInitSetup();

  if (vtu_writer_ptr_ != Teuchos::null) WriteTimeStepOutputRuntimeVtkStructure();

  // write special output for beams if desired
  if (beam_vtu_writer_ptr_ != Teuchos::null) WriteTimeStepOutputRuntimeVtkBeams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ResetStepState()
{
  CheckInitSetup();

  // reset disp, vel, acc state vector
  GStatePtr()->GetMutableDisNp()->Update(1.0, (*GStatePtr()->GetDisN()), 0.0);
  GStatePtr()->GetMutableVelNp()->Update(1.0, (*GStatePtr()->GetVelN()), 0.0);
  GStatePtr()->GetMutableAccNp()->Update(1.0, (*GStatePtr()->GetAccN()), 0.0);

  // other parameters that might be needed by the elements
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_reset_istep);

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};
  EvaluateInternal(eval_mat, eval_vec);

  DiscretPtr()->ClearState();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Structure::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::GetCurrentSolutionPtr() const
{
  CheckInit();
  return GState().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::GetLastTimeStepSolutionPtr() const
{
  CheckInit();
  return GState().GetDisN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::PostOutput()
{
  CheckInitSetup();
  // empty

  return;
}  // PostOutput()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FintNp()
{
  CheckInit();
  if (GState().GetMutableFintNp().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FintNp() const
{
  CheckInit();
  if (GState().GetFintNp().is_null()) dserror("NULL pointer!");

  return *GState().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FextNp()
{
  CheckInit();
  if (GState().GetMutableFextNp().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableFextNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FextNp() const
{
  CheckInit();
  if (GState().GetFextNp().is_null()) dserror("NULL pointer!");

  return *GState().GetFextNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FinertialNp()
{
  CheckInit();
  if (GState().GetMutableFinertialNp().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableFinertialNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FinertialNp() const
{
  CheckInit();
  if (GState().GetFinertialNp().is_null()) dserror("NULL pointer!");

  return *GState().GetFinertialNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FviscoNp()
{
  CheckInit();
  if (GState().GetMutableFviscoNp().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableFviscoNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FviscoNp() const
{
  CheckInit();
  if (GState().GetFviscoNp().is_null()) dserror("NULL pointer!");

  return *GState().GetFviscoNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::DisNp()
{
  CheckInit();
  if (GState().GetMutableDisNp().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::DisNp() const
{
  CheckInit();
  if (GState().GetDisNp().is_null()) dserror("NULL pointer!");

  return *GState().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::SparseMatrix& STR::MODELEVALUATOR::Structure::Stiff() const
{
  CheckInit();
  if (not stiff_ptr_) dserror("NULL pointer!");

  return *stiff_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::SparseMatrix& STR::MODELEVALUATOR::Structure::StiffPTC() const
{
  CheckInit();
  if (stiff_ptc_ptr_ == Teuchos::null) dserror("NULL pointer!");

  return *stiff_ptc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Mass()
{
  CheckInit();
  if (GState().GetMutableMassMatrix().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableMassMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Mass() const
{
  CheckInit();
  if (GState().GetMassMatrix().is_null()) dserror("NULL pointer!");

  return *GState().GetMassMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Damp()
{
  CheckInit();
  if (GState().GetMutableDampMatrix().is_null()) dserror("NULL pointer!");

  return *GState().GetMutableDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Damp() const
{
  CheckInit();
  if (GState().GetDampMatrix().is_null()) dserror("NULL pointer!");

  return *GState().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ParamsInterface2ParameterList(
    Teuchos::RCP<STR::MODELEVALUATOR::Data> interface_ptr, Teuchos::ParameterList& params)
{
  if (interface_ptr == Teuchos::null) dserror("ParamsInterface pointer not set");

  params.set<double>("delta time", interface_ptr->GetDeltaTime());
  params.set<double>("total time", interface_ptr->GetTotalTime());
  params.set<double>("timintfac_dis", interface_ptr->GetTimIntFactorDisp());
  params.set<double>("timintfac_vel", interface_ptr->GetTimIntFactorVel());

  DRT::ELEMENTS::ActionType act = interface_ptr->GetActionType();
  std::string action;
  switch (act)
  {
    case DRT::ELEMENTS::struct_calc_linstiff:
      action = "calc_struct_linstiff";
      break;
    case DRT::ELEMENTS::struct_calc_nlnstiff:
      action = "calc_struct_nlnstiff";
      break;
    case DRT::ELEMENTS::struct_calc_internalforce:
      action = "calc_struct_internalforce";
      break;
    case DRT::ELEMENTS::struct_calc_linstiffmass:
      action = "calc_struct_linstiffmass";
      break;
    case DRT::ELEMENTS::struct_calc_nlnstiffmass:
      action = "calc_struct_nlnstiffmass";
      break;
    case DRT::ELEMENTS::struct_calc_nlnstifflmass:
      action = "calc_struct_nlnstifflmass";
      break;
    case DRT::ELEMENTS::struct_calc_nlnstiff_gemm:
      action = "calc_struct_nlnstiff_gemm";
      break;
    case DRT::ELEMENTS::struct_calc_stress:
      action = "calc_struct_stress";
      break;
    case DRT::ELEMENTS::struct_calc_thickness:
      action = "calc_struct_thickness";
      break;
    case DRT::ELEMENTS::struct_calc_eleload:
      action = "calc_struct_eleload";
      break;
    case DRT::ELEMENTS::struct_calc_fsiload:
      action = "calc_struct_fsiload";
      break;
    case DRT::ELEMENTS::struct_calc_update_istep:
      action = "calc_struct_update_istep";
      break;
    case DRT::ELEMENTS::struct_calc_reset_istep:
      action = "calc_struct_reset_istep";
      break;
    case DRT::ELEMENTS::struct_calc_store_istep:
      action = "calc_struct_store_istep";
      break;
    case DRT::ELEMENTS::struct_calc_recover_istep:
      action = "calc_struct_recover_istep";
      break;
    case DRT::ELEMENTS::struct_calc_reset_all:
      action = "calc_struct_reset_all";
      break;
    case DRT::ELEMENTS::struct_calc_energy:
      action = "calc_struct_energy";
      break;
    case DRT::ELEMENTS::struct_calc_errornorms:
      action = "calc_struct_errornorms";
      break;
    case DRT::ELEMENTS::multi_init_eas:
      action = "multi_eas_init";
      break;
    case DRT::ELEMENTS::multi_set_eas:
      action = "multi_eas_set";
      break;
    case DRT::ELEMENTS::multi_readrestart:
      action = "multi_readrestart";
      break;
    case DRT::ELEMENTS::multi_calc_dens:
      action = "multi_calc_dens";
      break;
    case DRT::ELEMENTS::struct_postprocess_stress:
      action = "postprocess_stress";
      break;
    case DRT::ELEMENTS::struct_postprocess_thickness:
      action = "postprocess_thickness";
      break;
    case DRT::ELEMENTS::struct_update_prestress:
      action = "calc_struct_prestress_update";
      break;
    case DRT::ELEMENTS::inversedesign_update:
      action = "calc_struct_inversedesign_update";
      break;
    case DRT::ELEMENTS::inversedesign_switch:
      action = "calc_struct_inversedesign_switch";
      break;
    case DRT::ELEMENTS::struct_calc_global_gpstresses_map:
      action = "calc_global_gpstresses_map";
      break;
    case DRT::ELEMENTS::struct_interpolate_velocity_to_point:
      action = "interpolate_velocity_to_given_point";
      break;
    case DRT::ELEMENTS::struct_calc_mass_volume:
      action = "calc_struct_mass_volume";
      break;
    case DRT::ELEMENTS::struct_calc_recover:
      action = "calc_struct_recover";
      break;
    case DRT::ELEMENTS::struct_calc_predict:
      action = "calc_struct_predict";
      break;
    default:
      action = "unknown";
      break;
  }
  params.set<std::string>("action", action);

  params.set<Teuchos::RCP<std::vector<char>>>("stress", interface_ptr->MutableStressDataPtr());
  params.set<Teuchos::RCP<std::vector<char>>>("strain", interface_ptr->MutableStrainDataPtr());
  params.set<Teuchos::RCP<std::vector<char>>>(
      "plstrain", interface_ptr->MutablePlasticStrainDataPtr());
  params.set<Teuchos::RCP<std::vector<char>>>(
      "optquantity", interface_ptr->MutableOptQuantityDataPtr());
  params.set<int>("iostress", (int)interface_ptr->GetStressOutputType());
  params.set<int>("iostrain", (int)interface_ptr->GetStrainOutputType());
  params.set<int>("ioplstrain", (int)interface_ptr->GetPlasticStrainOutputType());
  params.set<int>("iooptquantity", (int)interface_ptr->GetOptQuantityOutputType());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::CreateBackupState(const Epetra_Vector& dir)
{
  CheckInitSetup();

  // set all parameters in the evaluation data container
  EvalData().SetActionType(DRT::ELEMENTS::struct_create_backup);

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Teuchos::RCP<const Epetra_Vector> dir_displ = GState().ExtractDisplEntries(dir);
  Discret().SetState(0, "residual displacement", dir_displ);

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RecoverFromBackupState()
{
  CheckInitSetup();

  // set all parameters in the evaluation data container
  EvalData().SetActionType(DRT::ELEMENTS::struct_recover_from_backup);

  // set vector values needed by elements
  Discret().ClearState();

  // set dummy evaluation vectors and matrices
  Teuchos::RCP<Epetra_Vector> eval_vec[3] = {Teuchos::null, Teuchos::null, Teuchos::null};
  Teuchos::RCP<LINALG::SparseOperator> eval_mat[2] = {Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat, eval_vec);
}
