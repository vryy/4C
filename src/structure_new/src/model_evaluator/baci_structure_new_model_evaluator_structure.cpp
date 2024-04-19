/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all structure terms


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_model_evaluator_structure.hpp"

#include "baci_beam3_base.hpp"
#include "baci_beam3_discretization_runtime_output_params.hpp"
#include "baci_beam3_discretization_runtime_vtu_writer.hpp"
#include "baci_discretization_fem_general_utils_gauss_point_postprocess.hpp"
#include "baci_global_data.hpp"
#include "baci_io.hpp"
#include "baci_io_discretization_visualization_writer_mesh.hpp"
#include "baci_io_pstream.hpp"
#include "baci_io_visualization_parameters.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_utils_discret.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_linalg_sparseoperator.hpp"
#include "baci_linalg_utils_sparse_algebra_assemble.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_structure_new_dbc.hpp"
#include "baci_structure_new_discretization_runtime_output_params.hpp"
#include "baci_structure_new_gauss_point_data_output_manager.hpp"
#include "baci_structure_new_integrator.hpp"
#include "baci_structure_new_model_evaluator_data.hpp"
#include "baci_structure_new_predict_generic.hpp"
#include "baci_structure_new_timint_basedataio_runtime_vtk_output.hpp"
#include "baci_structure_new_timint_implicit.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Structure::Structure()
    : dt_ele_ptr_(nullptr),
      masslin_type_(INPAR::STR::ml_none),
      stiff_ptr_(nullptr),
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
  FOUR_C_ASSERT(IsInit(), "Init() has not been called, yet!");

  // get the global state content
  {
    // structural element evaluation time
    dt_ele_ptr_ = &(GState().GetElementEvaluationTime());
  }

  // displ-displ block
  stiff_ptr_ =
      dynamic_cast<CORE::LINALG::SparseMatrix*>(GState().CreateStructuralStiffnessMatrixBlock());

  // modified stiffness pointer for storing element based scaling operator (PTC)
  stiff_ptc_ptr_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));

  FOUR_C_ASSERT(stiff_ptr_ != nullptr, "Dynamic cast to CORE::LINALG::SparseMatrix failed!");

  // get the structural dynamic content
  {
    // setup important evaluation booleans
    masslin_type_ = TimInt().GetDataSDyn().GetMassLinType();
  }
  // setup new variables
  {
    dis_incr_ptr_ = Teuchos::rcp(new Epetra_Vector(DisNp().Map(), true));
  }

  // setup output writers
  {
    if (GInOutput().GetRuntimeOutputParams() != Teuchos::null)
    {
      visualization_params_ = IO::VisualizationParametersFactory(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
          *GLOBAL::Problem::Instance()->OutputControlFile(), GState().GetTimeN());

      // We only want to create the respective writers if they are actually needed. Therefore, we
      // get the global number ob beam and non-beam elements here. Based on that number we know
      // which output writers need to be initialized.
      const auto discretization =
          Teuchos::rcp_dynamic_cast<const DRT::Discretization>(DiscretPtr(), true);
      int number_my_solid_elements = std::count_if(discretization->MyRowElementRange().begin(),
          discretization->MyRowElementRange().end(),
          [](const auto* row_element)
          { return dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(row_element) == nullptr; });
      int number_my_beam_elements = discretization->NumMyRowElements() - number_my_solid_elements;
      int number_global_solid_elements = 0;
      int number_global_beam_elements = 0;
      discretization->Comm().MaxAll(&number_my_solid_elements, &number_global_solid_elements, 1);
      discretization->Comm().MaxAll(&number_my_beam_elements, &number_global_beam_elements, 1);

      if (GInOutput().GetRuntimeOutputParams()->OutputStructure() &&
          number_global_solid_elements > 0)
        InitOutputRuntimeStructure();

      if (GInOutput().GetRuntimeOutputParams()->OutputBeams() && number_global_beam_elements > 0)
        InitOutputRuntimeBeams();
    }
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
   * putScalar(0.0), because of possible NaN and inf values! */
  FextNp().PutScalar(0.0);

  /* --- reset internal forces
   * Please note, that PutScalar is safer (but maybe slower) than
   * putScalar(0.0), because of possible NaN and inf values! */
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
  double dtcpu = GState().GetTimer()->wallTime();
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
  *dt_ele_ptr_ += GState().GetTimer()->wallTime() - dtcpu;
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
  double dtcpu = GState().GetTimer()->wallTime();
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
  *dt_ele_ptr_ += GState().GetTimer()->wallTime() - dtcpu;
  // *********** time measurement ***********

  // that's it
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::AssembleForce(Epetra_Vector& f, const double& timefac_np) const
{
  CORE::LINALG::AssembleMyVector(1.0, f, -timefac_np, FextNp());
  CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, FintNp());

  // add the scaled force contributions of the old time step
  // structural dofs of the right-hand-side vector at t_{n+timefac_n} (read-only)
  Teuchos::RCP<const Epetra_Vector> fstructold_ptr = GState().GetFstructureOld();
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fstructold_ptr);

  // add the visco and mass contributions
  Int().AddViscoMassContributions(f);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::AssembleJacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  int err = Stiff().Scale(timefac_np);
  GState().AssignModelBlock(jac, Stiff(), Type(), STR::MatBlockType::displ_displ);

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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  // create vector with zero entries
  Teuchos::RCP<const Epetra_Vector> zeros = Int().GetDbc().GetZerosPtr();

  // set vector values needed by elements
  // --> initially zero !!!
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", zeros);
  Discret().SetState(0, "displacement", zeros);

  // set action type and evaluation matrix and vector pointers
  StaticContributions(eval_mat.data(), eval_vec.data());
  MaterialDampingContributions(eval_mat.data());
  InertialContributions(eval_mat.data(), eval_vec.data());

  // evaluate
  EvaluateInternal(eval_mat.data(), eval_vec.data());

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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Discret().SetState(0, "velocity", GState().GetVelNp());

  // set action type and evaluation matrix and vector pointers
  StaticContributions(eval_vec.data());
  MaterialDampingContributions(eval_mat.data());
  InertialContributions(eval_vec.data());

  // evaluate ...
  EvaluateInternal(eval_mat.data(), eval_vec.data());

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
  EvaluateNeumann(GState().GetFextNp(), Teuchos::null);

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
    EvaluateNeumann(GState().GetFextNp(), Teuchos::null);
  else
  {
    Discret().SetState(0, "displacement new", GState().GetDisNp());
    /* Add the linearization of the external force to the stiffness
     * matrix. */
    EvaluateNeumann(GState().GetFextNp(), Teuchos::rcpFromRef(*stiff_ptr_));
  }

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::PreApplyForceStiffExternal(
    Epetra_Vector& fextnp, CORE::LINALG::SparseMatrix& stiff) const
{
  CheckInitSetup();

  const auto* impl_ptr = dynamic_cast<const STR::TIMINT::Implicit*>(&TimInt());
  if (impl_ptr) return impl_ptr->Predictor().PreApplyForceExternal(fextnp);

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::ApplyForceStiffInternal()
{
  CheckInitSetup();
  // currently a fixed number of matrix and vector pointers are supported
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);
  Discret().SetState(0, "displacement", GState().GetDisNp());
  Discret().SetState(0, "velocity", GState().GetVelNp());

  // set action types and evaluate matrices/vectors
  StaticContributions(eval_mat.data(), eval_vec.data());
  MaterialDampingContributions(eval_mat.data());
  if (masslin_type_ != INPAR::STR::ml_none) InertialContributions(eval_mat.data(), eval_vec.data());

  // evaluate
  EvaluateInternal(eval_mat.data(), eval_vec.data());

  // complete stiffness and mass matrix
  FillComplete();

  // evaluate inertial and viscous forces
  InertialAndViscousForces();

  return EvalErrorCheck();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::StaticContributions(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
  // set default matrix
  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptr_);
  // set default force vector
  eval_vec[0] = GState().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::StaticContributions(Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_internalforce);
  // set default force vector
  eval_vec[0] = GState().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::MaterialDampingContributions(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat)
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
  eval_mat[1] = GState().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InertialContributions(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
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
  eval_mat[1] = GState().GetMassMatrix();
  // set inertial vector if necessary
  eval_vec[1] = GetInertialForce();
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
      return GState().GetFinertialNp();
      break;
    }
    case INPAR::STR::ml_none:
      // do nothing
      break;
    default:
      FOUR_C_THROW("Unknown mass linearization type!");
      exit(EXIT_FAILURE);
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InitOutputRuntimeStructure()
{
  CheckInit();
  const auto discretization = Teuchos::rcp_dynamic_cast<const DRT::Discretization>(
      const_cast<STR::MODELEVALUATOR::Structure*>(this)->DiscretPtr(), true);
  vtu_writer_ptr_ = Teuchos::rcp(
      new IO::DiscretizationVisualizationWriterMesh(discretization, visualization_params_));

  if (GInOutput().GetRuntimeOutputParams()->GetStructureParams()->GaussPointDataOutput() !=
      INPAR::STR::GaussPointDataOutputType::none)
  {
    InitOutputRuntimeStructureGaussPointData();
  }
}

void STR::MODELEVALUATOR::Structure::InitOutputRuntimeStructureGaussPointData()
{
  // Set all parameters in the evaluation data container.
  EvalData().SetActionType(DRT::ELEMENTS::struct_init_gauss_point_data_output);
  EvalData().SetGaussPointDataOutputManagerPtr(Teuchos::rcp(new GaussPointDataOutputManager(
      GInOutput().GetRuntimeOutputParams()->GetStructureParams()->GaussPointDataOutput())));
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);

  // Set vector values needed by elements.
  Discret().ClearState();

  // Set dummy evaluation vectors and matrices.
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat.data(), eval_vec.data());

  EvalData().GetGaussPointDataOutputManagerPtr()->DistributeQuantities(Discret().Comm());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteTimeStepOutputRuntimeStructure() const
{
  CheckInitSetup();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisN(), *disn_col);
  Teuchos::RCP<Epetra_Vector> veln_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetVelN(), *veln_col);

  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, GState().GetTimeN(), GState().GetStepN());
  WriteOutputRuntimeStructure(disn_col, veln_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteIterationOutputRuntimeStructure() const
{
  CheckInitSetup();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disnp_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisNp(), *disnp_col);
  Teuchos::RCP<Epetra_Vector> velnp_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetVelNp(), *velnp_col);

  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, GState().GetTimeN(), GState().GetStepN(), EvalData().GetNlnIter());
  WriteOutputRuntimeStructure(disnp_col, velnp_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteOutputRuntimeStructure(
    const Teuchos::RCP<Epetra_Vector>& displacement_state_vector,
    const Teuchos::RCP<Epetra_Vector>& velocity_state_vector, int timestep_number,
    double time) const
{
  CheckInitSetup();

  // get the parameter container object
  const DRT::ELEMENTS::StructureRuntimeOutputParams& structure_output_params =
      *GInOutput().GetRuntimeOutputParams()->GetStructureParams();

  // reset time and time step of the writer object
  vtu_writer_ptr_->Reset();

  // append all desired output data to the writer object's storage

  // append displacement if desired
  if (structure_output_params.OutputDisplacementState())
    vtu_writer_ptr_->AppendDofBasedResultDataVector(
        displacement_state_vector, 3, 0, "displacement");

  // append velocity if desired
  if (structure_output_params.OutputVelocityState())
    vtu_writer_ptr_->AppendDofBasedResultDataVector(velocity_state_vector, 3, 0, "velocity");

  // append element owner if desired
  if (structure_output_params.OutputElementOwner())
    vtu_writer_ptr_->AppendElementOwner("element_owner");

  // append element GIDs if desired
  if (structure_output_params.OutputElementGID()) vtu_writer_ptr_->AppendElementGID("element_gid");

  // append element ghosting information if desired
  if (structure_output_params.OutputElementGhosting())
    vtu_writer_ptr_->AppendElementGhostingInformation();

  // append node GIDs if desired
  if (structure_output_params.OutputNodeGID()) vtu_writer_ptr_->AppendNodeGID("node_gid");

  // append stress if desired
  if (structure_output_params.OutputStressStrain() and
      GInOutput().GetStressOutputType() != INPAR::STR::stress_none)
  {
    std::string name_nodal = "";
    std::string name_element = "";

    if (GInOutput().GetStressOutputType() == INPAR::STR::stress_2pk)
    {
      name_nodal = "nodal_2PK_stresses_xyz";
      name_element = "element_2PK_stresses_xyz";
    }
    else if (GInOutput().GetStressOutputType() == INPAR::STR::stress_cauchy)
    {
      name_nodal = "nodal_cauchy_stresses_xyz";
      name_element = "element_cauchy_stresses_xyz";
    }

    // Write nodal stress data.
    vtu_writer_ptr_->AppendNodeBasedResultDataVector(
        EvalData().GetStressDataNodePostprocessed(), 6, name_nodal);

    // Write element stress data.
    vtu_writer_ptr_->AppendElementBasedResultDataVector(
        EvalData().GetStressDataElementPostprocessed(), 6, name_element);
  }

  // append strain if desired.
  if (structure_output_params.OutputStressStrain() and
      GInOutput().GetStrainOutputType() != INPAR::STR::strain_none)
  {
    std::string name_nodal = "";
    std::string name_element = "";

    if (GInOutput().GetStrainOutputType() == INPAR::STR::strain_gl)
    {
      name_nodal = "nodal_GL_strains_xyz";
      name_element = "element_GL_strains_xyz";
    }
    else if (GInOutput().GetStrainOutputType() == INPAR::STR::strain_ea)
    {
      name_nodal = "nodal_EA_strains_xyz";
      name_element = "element_EA_strains_xyz";
    }
    else if (GInOutput().GetStrainOutputType() == INPAR::STR::strain_log)
    {
      name_nodal = "nodal_LOG_strains_xyz";
      name_element = "element_LOG_strains_xyz";
    }

    // Write nodal strain data.
    vtu_writer_ptr_->AppendNodeBasedResultDataVector(
        EvalData().GetStrainDataNodePostprocessed(), 6, name_nodal);

    // Write element strain data.
    vtu_writer_ptr_->AppendElementBasedResultDataVector(
        EvalData().GetStrainDataElementPostprocessed(), 6, name_element);
  }

  // Add gauss point data if desired
  if (structure_output_params.GaussPointDataOutput() != INPAR::STR::GaussPointDataOutputType::none)
  {
    const GaussPointDataOutputManager& elementDataManager =
        *EvalData().GetGaussPointDataOutputManagerPtr();
    for (const auto& nameAndSize : elementDataManager.GetQuantities())
    {
      const std::string& name = nameAndSize.first;
      const int size = nameAndSize.second;

      switch (elementDataManager.GetOutputType())
      {
        case INPAR::STR::GaussPointDataOutputType::element_center:
        {
          Teuchos::RCP<Epetra_MultiVector> data =
              elementDataManager.GetElementCenterData().at(name);
          vtu_writer_ptr_->AppendElementBasedResultDataVector(data, size, name);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::gauss_points:
        {
          const std::vector<Teuchos::RCP<Epetra_MultiVector>>& data_list =
              elementDataManager.GetGaussPointData().at(name);
          for (std::size_t gp = 0; gp < data_list.size(); ++gp)
          {
            const std::string name_with_gp = name + "_gp_" + std::to_string(gp);
            vtu_writer_ptr_->AppendElementBasedResultDataVector(data_list[gp], size, name_with_gp);
          }
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::nodes:
        {
          Teuchos::RCP<Epetra_MultiVector> data = elementDataManager.GetNodalData().at(name);
          vtu_writer_ptr_->AppendNodeBasedResultDataVector(data, size, name);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::none:
          FOUR_C_THROW("Gauss point data output type is none");
        default:
          FOUR_C_THROW("Gauss point data output type is not implemented yet");
      }
    }
  }

  // finalize everything and write all required files to filesystem
  vtu_writer_ptr_->WriteToDisk(time, timestep_number);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::OutputRuntimeStructurePostprocessStressStrain()
{
  CheckInitSetup();

  if (not(GInOutput().GetStressOutputType() == INPAR::STR::stress_none and
          GInOutput().GetStrainOutputType() == INPAR::STR::strain_none))
  {
    // Set all parameters in the evaluation data container.
    EvalData().SetActionType(DRT::ELEMENTS::struct_calc_stress);
    EvalData().SetTotalTime(GState().GetTimeNp());
    EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);
    EvalData().SetStressData(Teuchos::rcp(new std::vector<char>()));
    EvalData().SetCouplingStressData(Teuchos::rcp(new std::vector<char>()));
    EvalData().SetStrainData(Teuchos::rcp(new std::vector<char>()));
    EvalData().SetPlasticStrainData(Teuchos::rcp(new std::vector<char>()));

    // Set vector values needed by elements.
    Discret().ClearState();
    Discret().SetState(0, "displacement", GState().GetDisNp());
    Discret().SetState(0, "residual displacement", dis_incr_ptr_);

    // GState().GetDisNp()->Print(std::cout);

    // Set dummy evaluation vectors and matrices.
    std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
        Teuchos::null, Teuchos::null, Teuchos::null};
    std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
        Teuchos::null, Teuchos::null};

    EvaluateInternalSpecifiedElements(eval_mat.data(), eval_vec.data(), Discret().ElementRowMap());

    auto DoPostprocessingOnElement = [](const DRT::Element& ele)
    {
      // If it is not a beam element, we post-process it.
      return dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(&ele) == nullptr;
    };

    auto EvaluateGaussPointData = [&](const std::vector<char>& raw_data)
    {
      // Get the values at the Gauss-points.
      std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> mapdata{};
      std::vector<char>::size_type position = 0;
      for (int i = 0; i < DiscretPtr()->ElementRowMap()->NumMyElements(); ++i)
      {
        if (DoPostprocessingOnElement(*Discret().lRowElement(i)))
        {
          Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstress =
              Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
          CORE::COMM::ParObject::ExtractfromPack(position, raw_data, *gpstress);
          mapdata[DiscretPtr()->ElementRowMap()->GID(i)] = gpstress;
        }
      }
      return mapdata;
    };

    auto PostprocessGaussPointDataToNodes =
        [&](const std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>& map_data,
            Epetra_MultiVector& assembled_data)
    {
      DiscretPtr()->Evaluate(
          [&](DRT::Element& ele)
          {
            if (DoPostprocessingOnElement(ele))
              CORE::FE::ExtrapolateGaussPointQuantityToNodes(
                  ele, *map_data.at(ele.Id()), Discret(), assembled_data);
          });
    };

    auto PostprocessGaussPointDataToElementCenter =
        [&](const std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>& map_data,
            Epetra_MultiVector& assembled_data)
    {
      DiscretPtr()->Evaluate(
          [&](DRT::Element& ele)
          {
            if (DoPostprocessingOnElement(ele))
              CORE::FE::EvaluateGaussPointQuantityAtElementCenter(
                  ele, *map_data.at(ele.Id()), assembled_data);
          });
    };

    const auto* discret = dynamic_cast<const DRT::Discretization*>(&Discret());

    // Postprocess the result vectors.
    if (GInOutput().GetStressOutputType() != INPAR::STR::stress_none)
    {
      std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> gp_stress_data =
          EvaluateGaussPointData(*EvalData().GetStressData());

      CORE::COMM::Exporter ex(
          *(Discret().ElementRowMap()), *(discret->ElementColMap()), Discret().Comm());
      ex.Export(gp_stress_data);

      EvalData().GetStressDataNodePostprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret->NodeColMap(), 6, true));
      EvalData().GetStressDataElementPostprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret->ElementRowMap(), 6, true));


      Epetra_MultiVector row_nodal_data(*discret->NodeRowMap(), 6, true);
      PostprocessGaussPointDataToNodes(gp_stress_data, row_nodal_data);
      CORE::LINALG::Export(row_nodal_data, *EvalData().GetStressDataNodePostprocessed());

      PostprocessGaussPointDataToElementCenter(
          gp_stress_data, *EvalData().GetStressDataElementPostprocessed());
    }
    if (GInOutput().GetStrainOutputType() != INPAR::STR::strain_none)
    {
      std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> gp_strain_data =
          EvaluateGaussPointData(*EvalData().GetStrainData());

      CORE::COMM::Exporter ex(
          *(Discret().ElementRowMap()), *(discret->ElementColMap()), Discret().Comm());
      ex.Export(gp_strain_data);

      EvalData().GetStrainDataNodePostprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret->NodeColMap(), 6, true));
      EvalData().GetStrainDataElementPostprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret->ElementRowMap(), 6, true));

      Epetra_MultiVector row_nodal_data(*discret->NodeRowMap(), 6, true);
      PostprocessGaussPointDataToNodes(gp_strain_data, row_nodal_data);
      CORE::LINALG::Export(row_nodal_data, *EvalData().GetStrainDataNodePostprocessed());

      PostprocessGaussPointDataToElementCenter(
          gp_strain_data, *EvalData().GetStrainDataElementPostprocessed());
    }
  }
}

void STR::MODELEVALUATOR::Structure::OutputRuntimeStructureGaussPointData()
{
  const DRT::ELEMENTS::StructureRuntimeOutputParams& structure_output_params =
      *GInOutput().GetRuntimeOutputParams()->GetStructureParams();
  if (structure_output_params.GaussPointDataOutput() != INPAR::STR::GaussPointDataOutputType::none)
  {
    CheckInitSetup();

    EvalData().SetActionType(DRT::ELEMENTS::struct_gauss_point_data_output);
    EvalData().SetTotalTime(GState().GetTimeNp());
    EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);

    const auto* discret = dynamic_cast<const DRT::Discretization*>(&Discret());
    EvalData().GaussPointDataOutputManagerPtr()->PrepareData(
        *discret->NodeColMap(), *discret->ElementRowMap());

    Discret().ClearState();
    Discret().SetState(0, "displacement", GState().GetDisNp());
    Discret().SetState(0, "residual displacement", dis_incr_ptr_);

    std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
        Teuchos::null, Teuchos::null, Teuchos::null};
    std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
        Teuchos::null, Teuchos::null};

    EvaluateInternal(eval_mat.data(), eval_vec.data());

    EvalData().GaussPointDataOutputManagerPtr()->PostEvaluate();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::InitOutputRuntimeBeams()
{
  beam_vtu_writer_ptr_ = Teuchos::rcp(
      new BeamDiscretizationRuntimeOutputWriter(visualization_params_, DisNp().Comm()));

  // get the parameter container object
  const DRT::ELEMENTS::BeamRuntimeOutputParams& beam_output_params =
      *GInOutput().GetRuntimeOutputParams()->GetBeamParams();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisN(), *disn_col);

  // get bounding box object only if periodic boundaries are active
  Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> bounding_box_ptr =
      TimInt().GetDataSDynPtr()->GetPeriodicBoundingBox();

  // initialize the writer object with current displacement state
  beam_vtu_writer_ptr_->Initialize(const_cast<STR::MODELEVALUATOR::Structure*>(this)->DiscretPtr(),
      beam_output_params.UseAbsolutePositions(),
      beam_output_params.GetNumberVisualizationSubsegments(), bounding_box_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteTimeStepOutputRuntimeBeams() const
{
  CheckInitSetup();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisN(), *disn_col);

  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, GState().GetTimeN(), GState().GetStepN());
  WriteOutputRuntimeBeams(disn_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteIterationOutputRuntimeBeams() const
{
  CheckInitSetup();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disnp_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisNp(), *disnp_col);

  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, GState().GetTimeN(), GState().GetStepN(), EvalData().GetNlnIter());
  WriteOutputRuntimeBeams(disnp_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::WriteOutputRuntimeBeams(
    const Teuchos::RCP<Epetra_Vector>& displacement_state_vector, int timestep_number,
    double time) const
{
  CheckInitSetup();

  // get the parameter container object
  const DRT::ELEMENTS::BeamRuntimeOutputParams& beam_output_params =
      *GInOutput().GetRuntimeOutputParams()->GetBeamParams();

  // set geometry
  beam_vtu_writer_ptr_->SetGeometryFromBeamDiscretization(displacement_state_vector);

  // append all desired output data to the writer object's storage
  beam_vtu_writer_ptr_->AppendElementOwningProcessor();

  // append beam radius
  beam_vtu_writer_ptr_->AppendElementCircularCrossSectionRadius();

  // append displacement if desired
  if (beam_output_params.IsWriteInternalEnergyElement())
    beam_vtu_writer_ptr_->AppendElementInternalEnergy();

  // append displacement if desired
  if (beam_output_params.IsWriteKineticEnergyElement())
    beam_vtu_writer_ptr_->AppendElementKineticEnergy();

  // append displacement if desired
  if (beam_output_params.OutputDisplacementState())
    beam_vtu_writer_ptr_->AppendDisplacementField(displacement_state_vector);

  // append triads if desired
  if (beam_output_params.IsWriteTriadVisualizationPoints())
    beam_vtu_writer_ptr_->AppendTriadField(displacement_state_vector);

  // append material cross-section strain resultants if desired
  if (beam_output_params.IsWriteMaterialStrainsGaussPoints())
    beam_vtu_writer_ptr_->AppendGaussPointMaterialCrossSectionStrainResultants();

  // append material cross-section strain resultants if desired
  if (beam_output_params.IsWriteMaterialStrainsContinuous())
    beam_vtu_writer_ptr_->AppendGaussPointMaterialCrossSectionStrainResultantsContinuous();

  // append material cross-section stress resultants if desired
  if (beam_output_params.IsWriteMaterialStressesGaussPoints())
    beam_vtu_writer_ptr_->AppendGaussPointMaterialCrossSectionStressResultants();

  // append material cross-section stress resultants if desired
  if (beam_output_params.IsWriteMaterialStressContinuous())
    beam_vtu_writer_ptr_->AppendGaussPointMaterialCrossSectionStressResultantsContinuous();

  // append filament id and type if desired
  if (beam_output_params.IsWriteElementFilamentCondition())
    beam_vtu_writer_ptr_->AppendElementFilamentIdAndType();

  // append filament id and type if desired
  if (beam_output_params.IsWriteOrientationParamter())
    beam_vtu_writer_ptr_->AppendElementOrientationParamater(displacement_state_vector);

  // append reference length if desired.
  if (beam_output_params.IsWriteRefLength()) beam_vtu_writer_ptr_->AppendRefLength();

  // export displacement state to column format
  if (beam_output_params.IsWriteRVECrosssectionForces())
    beam_vtu_writer_ptr_->AppendRVECrosssectionForces(displacement_state_vector);

  // export beam element IDs
  if (beam_output_params.IsWriteElementGID()) beam_vtu_writer_ptr_->AppendElementGID();

  // Ghosting information
  if (beam_output_params.IsWriteElementGhosting())
    beam_vtu_writer_ptr_->AppendElementGhostingInformation();

  // finalize everything and write all required VTU files to filesystem
  beam_vtu_writer_ptr_->WriteToDisk(time, timestep_number);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternal(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  PreEvaluateInternal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());

  EvaluateInternal(p, eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternal(Teuchos::ParameterList& p,
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  if (p.numParams() > 1)
  {
    FOUR_C_THROW(
        "Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  }
  if (not p.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>(
          "interface"))
    FOUR_C_THROW("The given parameter has the wrong type!");

  // FixMe as soon as possible: write data to the parameter list.
  // this is about to go, once the old time integration is deleted
  ParamsInterface2ParameterList(EvalDataPtr(), p);

  Discret().Evaluate(p, eval_mat[0], eval_mat[1], eval_vec[0], eval_vec[1], eval_vec[2]);
  Discret().ClearState();
}

void STR::MODELEVALUATOR::Structure::EvaluateInternalSpecifiedElements(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec,
    const Epetra_Map* ele_map_to_be_evaluated)
{
  PreEvaluateInternal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());

  EvaluateInternalSpecifiedElements(p, eval_mat, eval_vec, ele_map_to_be_evaluated);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateInternalSpecifiedElements(Teuchos::ParameterList& p,
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec,
    const Epetra_Map* ele_map_to_be_evaluated)
{
  if (p.numParams() > 1)
  {
    FOUR_C_THROW(
        "Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  }
  if (not p.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>(
          "interface"))
    FOUR_C_THROW("The given parameter has the wrong type!");

  // write data to the parameter list.
  // this is about to go, once the old time integration is deleted
  ParamsInterface2ParameterList(EvalDataPtr(), p);

  DRT::UTILS::Evaluate(*DiscretPtr(), p, *eval_mat, *eval_vec, ele_map_to_be_evaluated);

  Discret().ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateNeumann(const Teuchos::RCP<Epetra_Vector>& eval_vec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat)
{
  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());
  EvaluateNeumann(p, eval_vec, eval_mat);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateNeumann(Teuchos::ParameterList& p,
    const Teuchos::RCP<Epetra_Vector>& eval_vec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat)
{
  if (p.numParams() > 1)
  {
    FOUR_C_THROW(
        "Please use the STR::ELEMENTS::Interface and its derived "
        "classes to set and get parameters.");
  }
  if (not p.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>(
          "interface"))
    FOUR_C_THROW("The given parameter has the wrong type!");
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
  iowriter.WriteVector("fint", GState().GetFintN());

  if (forced_writerestart) return;

  iowriter.WriteVector("displacement", GState().GetDisN());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ReadRestart(IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();
  // read structural force vector
  ioreader.ReadVector(GState().GetFstructureOld(), "fstructure_old");
  ioreader.ReadVector(GState().GetFintN(), "fint");
  // read displacement field
  Teuchos::RCP<Epetra_Vector>& disnp = GState().GetDisNp();
  ioreader.ReadVector(disnp, "displacement");
  GState().GetMultiDis()->UpdateSteps(*disnp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Predict(const INPAR::STR::PredEnum& pred_type)
{
  // set the element action
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_predict);
  EvalData().SetPredictorType(pred_type);

  // set the matrix and vector pointers to Teuchos::null
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RunPreComputeX(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  CheckInitSetup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RunRecover()
{
  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "residual displacement", dis_incr_ptr_);
  Discret().SetState(0, "displacement", GState().GetDisNp());
  // set the element action
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_recover);
  // set the matrix and vector pointers to Teuchos::null
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat.data(), eval_vec.data());
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
  RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RunPostIterate(const ::NOX::Solver::Generic& solver)
{
  CheckInitSetup();

  if (vtu_writer_ptr_ != Teuchos::null and
      GInOutput().GetRuntimeOutputParams()->OutputEveryIteration())
  {
    OutputRuntimeStructurePostprocessStressStrain();
    OutputRuntimeStructureGaussPointData();
    WriteIterationOutputRuntimeStructure();
  }

  // write special output for beams if desired
  if (beam_vtu_writer_ptr_ != Teuchos::null and
      GInOutput().GetRuntimeOutputParams()->OutputEveryIteration())
    WriteIterationOutputRuntimeBeams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepState(const double& timefac_n)
{
  CheckInitSetup();
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  GState().GetMultiDis()->UpdateSteps(DisNp());

  // new velocities at t_{n+1} -> t_{n}
  //    V_{n} := V_{n+1}
  GState().GetMultiVel()->UpdateSteps(*GState().GetVelNp());

  // new at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  GState().GetMultiAcc()->UpdateSteps(*GState().GetAccNp());

  // store the old external force
  GState().GetFextN()->Scale(1.0, FextNp());

  // store the old reaction force
  GState().GetFreactN()->Scale(1.0, *GState().GetFreactNp());

  // store the old internal force
  GState().GetFintN()->Scale(1.0, FintNp());

  // new at t_{n+1} -> t_{n+timefac_n}
  //    F^{struct}_{n+timefac_n} := timefac_n * F^{struct}_{n+1}
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetFstructureOld();
  fstructold_ptr->Update(timefac_n, FintNp(), 1.0);
  fstructold_ptr->Update(-timefac_n, FextNp(), 1.0);

  // set the displacement increment back to zero
  dis_incr_ptr_->PutScalar(0.0);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::EvaluateJacobianContributionsFromElementLevelForPTC()
{
  CheckInitSetup();
  // currently a fixed number of matrix and vector pointers are supported
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_addjacPTC);

  // set vector values needed by elements
  Discret().ClearState();
  Discret().SetState(0, "displacement", GState().GetDisNp());

  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptc_ptr_);

  // evaluate
  EvaluateInternal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::AssembleJacobianContributionsFromElementLevelForPTC(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n)
{
  GState().AssignModelBlock(*modjac, StiffPTC(), Type(), STR::MatBlockType::displ_displ);
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
  bool isDuringPrestressing = prestress_type != INPAR::STR::PreStress::none &&
                              GState().GetTimeN() <= prestress_time + 1.0e-15;

  if (isDuringPrestressing && prestress_type == INPAR::STR::PreStress::mulf)
  {
    if (Discret().Comm().MyPID() == 0)
      IO::cout << "====== Entering PRESTRESSING update" << IO::endl;

    // Choose special update action for elements in case of MULF
    EvalData().SetActionType(DRT::ELEMENTS::struct_update_prestress);
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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};
  EvaluateInternal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateResidual()
{
  CheckInitSetup();
  dis_incr_ptr_->Update(-1.0, *GState().GetDisN(), 1.0, *GState().GetDisNp(), 0.0);
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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternalSpecifiedElements(eval_mat.data(), eval_vec.data(), Discret().ElementRowMap());
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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  PreEvaluateInternal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalDataPtr());

  // evaluate energy contributions on element level (row elements only)
  EvaluateInternalSpecifiedElements(p, eval_mat.data(), eval_vec.data(), Discret().ElementRowMap());

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
  if (masslin_type_ == INPAR::STR::ml_none and velnp != nullptr)
  {
    double kinetic_energy_times2 = 0.0;

    Teuchos::RCP<Epetra_Vector> linear_momentum =
        CORE::LINALG::CreateVector(*GState().DofRowMapView(), true);

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
      FOUR_C_THROW("Type of optional quantity not implemented yet!");
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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat.data(), eval_vec.data());
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

  auto& discret = dynamic_cast<DRT::Discretization&>(Discret());

  // start evaluation
  const Epetra_Map* relemap = Discret().ElementRowMap();
  ele_vols = Teuchos::rcp(new Epetra_Vector(*relemap, true));
  const unsigned my_num_reles = relemap->NumMyElements();

  DRT::Element::LocationArray la(Discret().NumDofSets());
  CORE::LINALG::SerialDenseVector ele_vol(6, true);

  CORE::LINALG::SerialDenseMatrix empty_dummy_mat;
  CORE::LINALG::SerialDenseVector empty_dummy_vec;

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
    ele_vol.putScalar(0.0);
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
void STR::MODELEVALUATOR::Structure::RuntimePreOutputStepState()
{
  CheckInitSetup();

  if (vtu_writer_ptr_ != Teuchos::null)
  {
    OutputRuntimeStructurePostprocessStressStrain();
    OutputRuntimeStructureGaussPointData();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RuntimeOutputStepState() const
{
  CheckInitSetup();

  if (vtu_writer_ptr_ != Teuchos::null) WriteTimeStepOutputRuntimeStructure();

  // write special output for beams if desired
  if (beam_vtu_writer_ptr_ != Teuchos::null) WriteTimeStepOutputRuntimeBeams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ResetStepState()
{
  CheckInitSetup();

  // reset disp, vel, acc state vector
  GStatePtr()->GetDisNp()->Update(1.0, (*GStatePtr()->GetDisN()), 0.0);
  GStatePtr()->GetVelNp()->Update(1.0, (*GStatePtr()->GetVelN()), 0.0);
  GStatePtr()->GetAccNp()->Update(1.0, (*GStatePtr()->GetAccN()), 0.0);

  // other parameters that might be needed by the elements
  EvalData().SetTotalTime(GState().GetTimeNp());
  EvalData().SetDeltaTime((*GState().GetDeltaTime())[0]);
  // action for elements
  EvalData().SetActionType(DRT::ELEMENTS::struct_calc_reset_istep);

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};
  EvaluateInternal(eval_mat.data(), eval_vec.data());

  DiscretPtr()->ClearState();
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
}  // PostOutput()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FintNp()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFintNp().is_null(), "nullptr!");

  return *GState().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FintNp() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFintNp().is_null(), "nullptr!");

  return *GState().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FintN() const
{
  CheckInit();
  if (GState().GetFintN().is_null()) FOUR_C_THROW("NULL pointer!");

  return *GState().GetFintN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FextNp()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFextNp().is_null(), "nullptr!");

  return *GState().GetFextNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FextNp() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFextNp().is_null(), "nullptr!");

  return *GState().GetFextNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FextN() const
{
  CheckInit();
  if (GState().GetFextN().is_null()) FOUR_C_THROW("NULL pointer!");

  return *GState().GetFextN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FinertialNp()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFinertialNp().is_null(), "nullptr!");

  return *GState().GetFinertialNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FinertialNp() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFinertialNp().is_null(), "nullptr!");

  return *GState().GetFinertialNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::FviscoNp()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFviscoNp().is_null(), "nullptr!");

  return *GState().GetFviscoNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::FviscoNp() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetFviscoNp().is_null(), "nullptr!");

  return *GState().GetFviscoNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::DisNp()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetDisNp().is_null(), "nullptr!");

  return *GState().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::DisNp() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetDisNp().is_null(), "nullptr!");

  return *GState().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix& STR::MODELEVALUATOR::Structure::Stiff() const
{
  CheckInit();
  FOUR_C_ASSERT(stiff_ptr_, "nullptr!");

  return *stiff_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix& STR::MODELEVALUATOR::Structure::StiffPTC() const
{
  CheckInit();
  FOUR_C_ASSERT(stiff_ptc_ptr_ != Teuchos::null, "nullptr!");

  return *stiff_ptc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Mass()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetMassMatrix().is_null(), "nullptr!");

  return *GState().GetMassMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Mass() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetMassMatrix().is_null(), "nullptr!");

  return *GState().GetMassMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Damp()
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetDampMatrix().is_null(), "nullptr!");

  return *GState().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::Damp() const
{
  CheckInit();
  FOUR_C_ASSERT(!GState().GetDampMatrix().is_null(), "nullptr!");

  return *GState().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ParamsInterface2ParameterList(
    Teuchos::RCP<STR::MODELEVALUATOR::Data> interface_ptr, Teuchos::ParameterList& params)
{
  FOUR_C_ASSERT(interface_ptr != Teuchos::null, "ParamsInterface pointer not set");

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
    case DRT::ELEMENTS::struct_calc_energy:
      action = "calc_struct_energy";
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
    case DRT::ELEMENTS::struct_postprocess_thickness:
      action = "postprocess_thickness";
      break;
    case DRT::ELEMENTS::struct_update_prestress:
      action = "calc_struct_prestress_update";
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
    case DRT::ELEMENTS::struct_init_gauss_point_data_output:
      action = "struct_init_gauss_point_data_output";
      break;
    case DRT::ELEMENTS::struct_gauss_point_data_output:
      action = "struct_gauss_point_data_output";
      break;
    case DRT::ELEMENTS::struct_poro_calc_fluidcoupling:
      action = "struct_poro_calc_fluidcoupling";
      break;
    case DRT::ELEMENTS::struct_poro_calc_scatracoupling:
      action = "struct_poro_calc_scatracoupling";
      break;
    default:
      action = "unknown";
      break;
  }
  params.set<std::string>("action", action);

  params.set<Teuchos::RCP<std::vector<char>>>("stress", interface_ptr->StressDataPtr());
  params.set<Teuchos::RCP<std::vector<char>>>("strain", interface_ptr->StrainDataPtr());
  params.set<Teuchos::RCP<std::vector<char>>>("plstrain", interface_ptr->PlasticStrainDataPtr());
  params.set<Teuchos::RCP<std::vector<char>>>("optquantity", interface_ptr->OptQuantityDataPtr());
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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat.data(), eval_vec.data());
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
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  EvaluateInternal(eval_mat.data(), eval_vec.data());
}

FOUR_C_NAMESPACE_CLOSE
