/*-----------------------------------------------------------*/
/*! \file

\brief Evaluation and assembly of all structure terms


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_model_evaluator_structure.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beam3_discretization_runtime_output_params.hpp"
#include "4C_beam3_discretization_runtime_vtu_writer.hpp"
#include "4C_discretization_fem_general_utils_gauss_point_postprocess.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_discretization_visualization_writer_mesh.hpp"
#include "4C_io_pstream.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_utils_discret.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_discretization_runtime_output_params.hpp"
#include "4C_structure_new_integrator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_predict_generic.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"
#include "4C_structure_new_timint_implicit.hpp"
#include "4C_utils_exceptions.hpp"

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
  FOUR_C_ASSERT(is_init(), "Init() has not been called, yet!");

  // get the global state content
  {
    // structural element evaluation time
    dt_ele_ptr_ = &(g_state().get_element_evaluation_time());
  }

  // displ-displ block
  stiff_ptr_ = dynamic_cast<CORE::LINALG::SparseMatrix*>(
      g_state().create_structural_stiffness_matrix_block());

  // modified stiffness pointer for storing element based scaling operator (PTC)
  stiff_ptc_ptr_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*g_state().DofRowMapView(), 81, true, true));

  FOUR_C_ASSERT(stiff_ptr_ != nullptr, "Dynamic cast to CORE::LINALG::SparseMatrix failed!");

  // get the structural dynamic content
  {
    // setup important evaluation booleans
    masslin_type_ = TimInt().GetDataSDyn().GetMassLinType();
  }
  // setup new variables
  {
    dis_incr_ptr_ = Teuchos::rcp(new Epetra_Vector(dis_np().Map(), true));
  }

  // setup output writers
  {
    if (g_in_output().get_runtime_output_params() != Teuchos::null)
    {
      visualization_params_ = IO::VisualizationParametersFactory(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
          *GLOBAL::Problem::Instance()->OutputControlFile(), g_state().GetTimeN());

      // We only want to create the respective writers if they are actually needed. Therefore, we
      // get the global number ob beam and non-beam elements here. Based on that number we know
      // which output writers need to be initialized.
      const auto discretization =
          Teuchos::rcp_dynamic_cast<const DRT::Discretization>(discret_ptr(), true);
      int number_my_solid_elements = std::count_if(discretization->MyRowElementRange().begin(),
          discretization->MyRowElementRange().end(),
          [](const auto* row_element)
          { return dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(row_element) == nullptr; });
      int number_my_beam_elements = discretization->NumMyRowElements() - number_my_solid_elements;
      int number_global_solid_elements = 0;
      int number_global_beam_elements = 0;
      discretization->Comm().MaxAll(&number_my_solid_elements, &number_global_solid_elements, 1);
      discretization->Comm().MaxAll(&number_my_beam_elements, &number_global_beam_elements, 1);

      if (g_in_output().get_runtime_output_params()->OutputStructure() &&
          number_global_solid_elements > 0)
        init_output_runtime_structure();

      if (g_in_output().get_runtime_output_params()->OutputBeams() &&
          number_global_beam_elements > 0)
        init_output_runtime_beams();
    }
  }

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Reset(const Epetra_Vector& x)
{
  check_init_setup();

  /* --- reset external forces
   * Please note, that PutScalar is safer (but maybe slower) than
   * putScalar(0.0), because of possible NaN and inf values! */
  fext_np().PutScalar(0.0);

  /* --- reset internal forces
   * Please note, that PutScalar is safer (but maybe slower) than
   * putScalar(0.0), because of possible NaN and inf values! */
  fint_np().PutScalar(0.0);

  // reset stiffness matrix
  stiff().Zero();

  // reset modified stiffness matrix
  stiff_ptc().Zero();


  // set evaluation time back to zero
  *dt_ele_ptr_ = 0.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::evaluate_force()
{
  check_init_setup();
  bool ok = true;
  // ---------------------------------------
  // (1) EXTERNAL FORCES
  // ---------------------------------------
  ok = apply_force_external();

  // ---------------------------------------
  // (2) INTERNAL FORCES
  // ---------------------------------------
  // ordinary internal force
  ok = (ok ? apply_force_internal() : false);

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::evaluate_stiff()
{
  check_init_setup();
  bool ok = true;

  /* We use the same routines as for the ApplyForceStiff case, but we
   * do not update the global force vector, which is used for the
   * solution process in the NOX library.
   * This is meaningful, since the computational overhead, which is
   * generated by evaluating the right hand side is negligible */
  // *********** time measurement ***********
  double dtcpu = g_state().GetTimer()->wallTime();
  // *********** time measurement ***********
  // ---------------------------------------------------------------------
  // (1) EXTRERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  ok = apply_force_stiff_external();

  // ---------------------------------------------------------------------
  // (2) INTERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  // ordinary internal force
  ok = (ok ? apply_force_stiff_internal() : false);

  // *********** time measurement ***********
  *dt_ele_ptr_ += g_state().GetTimer()->wallTime() - dtcpu;
  // *********** time measurement ***********

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::evaluate_force_stiff()
{
  check_init_setup();
  bool ok = true;

  // *********** time measurement ***********
  double dtcpu = g_state().GetTimer()->wallTime();
  // *********** time measurement ***********
  // ---------------------------------------------------------------------
  // (1) EXTRERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  ok = apply_force_stiff_external();

  // ---------------------------------------------------------------------
  // (2) INTERNAL FORCES and STIFFNESS ENTRIES
  // ---------------------------------------------------------------------
  // ordinary internal force
  ok = (ok ? apply_force_stiff_internal() : false);

  // *********** time measurement ***********
  *dt_ele_ptr_ += g_state().GetTimer()->wallTime() - dtcpu;
  // *********** time measurement ***********

  // that's it
  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  CORE::LINALG::AssembleMyVector(1.0, f, -timefac_np, fext_np());
  CORE::LINALG::AssembleMyVector(1.0, f, timefac_np, fint_np());

  // add the scaled force contributions of the old time step
  // structural dofs of the right-hand-side vector at t_{n+timefac_n} (read-only)
  Teuchos::RCP<const Epetra_Vector> fstructold_ptr = GState().GetFstructureOld();
  CORE::LINALG::AssembleMyVector(1.0, f, 1.0, *fstructold_ptr);

  // add the visco and mass contributions
  Int().add_visco_mass_contributions(f);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::assemble_jacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  int err = stiff().Scale(timefac_np);
  GState().AssignModelBlock(jac, stiff(), Type(), STR::MatBlockType::displ_displ);

  // add the visco and mass contributions
  Int().add_visco_mass_contributions(jac);

  return (err == 0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::initialize_inertia_and_damping()
{
  check_init_setup();

  // currently a fixed number of matrix and vector pointers are supported
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  // create vector with zero entries
  Teuchos::RCP<const Epetra_Vector> zeros = Int().GetDbc().GetZerosPtr();

  // set vector values needed by elements
  // --> initially zero !!!
  discret().ClearState();
  discret().set_state(0, "residual displacement", zeros);
  discret().set_state(0, "displacement", zeros);

  // set action type and evaluation matrix and vector pointers
  static_contributions(eval_mat.data(), eval_vec.data());
  material_damping_contributions(eval_mat.data());
  inertial_contributions(eval_mat.data(), eval_vec.data());

  // evaluate
  evaluate_internal(eval_mat.data(), eval_vec.data());

  // complete stiffness and mass matrix
  fill_complete();

  // assemble the rayleigh damping matrix
  rayleigh_damping_matrix();

  return eval_error_check();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::apply_force_internal()
{
  check_init_setup();

  // currently a fixed number of matrix and vector pointers are supported
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "residual displacement", dis_incr_ptr_);
  discret().set_state(0, "displacement", g_state().GetDisNp());
  discret().set_state(0, "velocity", g_state().GetVelNp());

  // set action type and evaluation matrix and vector pointers
  static_contributions(eval_vec.data());
  material_damping_contributions(eval_mat.data());
  inertial_contributions(eval_vec.data());

  // evaluate ...
  evaluate_internal(eval_mat.data(), eval_vec.data());

  // evaluate inertia and visco forces
  inertial_and_viscous_forces();

  return eval_error_check();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::apply_force_external()
{
  check_init_setup();

  // Set to default value, because it is unnecessary for the
  // evaluate_neumann routine.
  eval_data().SetActionType(DRT::ELEMENTS::none);
  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", g_state().GetDisN());
  if (eval_data().GetDampingType() == INPAR::STR::damp_material)
    discret().set_state(0, "velocity", g_state().GetVelN());
  discret().set_state(0, "displacement new", g_state().GetDisNp());
  evaluate_neumann(g_state().GetFextNp(), Teuchos::null);

  return eval_error_check();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::apply_force_stiff_external()
{
  check_init_setup();

  if (pre_apply_force_stiff_external(fext_np(), *stiff_ptr_)) return true;

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", g_state().GetDisN());

  if (eval_data().GetDampingType() == INPAR::STR::damp_material)
    discret().set_state(0, "velocity", g_state().GetVelN());

  // get load vector
  if (!TimInt().GetDataSDyn().GetLoadLin())
    evaluate_neumann(g_state().GetFextNp(), Teuchos::null);
  else
  {
    discret().set_state(0, "displacement new", g_state().GetDisNp());
    /* Add the linearization of the external force to the stiffness
     * matrix. */
    evaluate_neumann(g_state().GetFextNp(), Teuchos::rcpFromRef(*stiff_ptr_));
  }

  return eval_error_check();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::pre_apply_force_stiff_external(
    Epetra_Vector& fextnp, CORE::LINALG::SparseMatrix& stiff) const
{
  check_init_setup();

  const auto* impl_ptr = dynamic_cast<const STR::TIMINT::Implicit*>(&TimInt());
  if (impl_ptr) return impl_ptr->Predictor().pre_apply_force_external(fextnp);

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::apply_force_stiff_internal()
{
  check_init_setup();
  // currently a fixed number of matrix and vector pointers are supported
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "residual displacement", dis_incr_ptr_);
  discret().set_state(0, "displacement", g_state().GetDisNp());
  discret().set_state(0, "velocity", g_state().GetVelNp());

  // set action types and evaluate matrices/vectors
  static_contributions(eval_mat.data(), eval_vec.data());
  material_damping_contributions(eval_mat.data());
  if (masslin_type_ != INPAR::STR::ml_none)
    inertial_contributions(eval_mat.data(), eval_vec.data());

  // evaluate
  evaluate_internal(eval_mat.data(), eval_vec.data());

  // complete stiffness and mass matrix
  fill_complete();

  // evaluate inertial and viscous forces
  inertial_and_viscous_forces();

  return eval_error_check();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::static_contributions(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  // action for elements
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
  // set default matrix
  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptr_);
  // set default force vector
  eval_vec[0] = g_state().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::static_contributions(Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  // action for elements
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_internalforce);
  // set default force vector
  eval_vec[0] = g_state().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::material_damping_contributions(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat)
{
  if (eval_data().GetDampingType() != INPAR::STR::damp_material) return;

  // action for elements
  // (reset the action type to be independent of the calling order)
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiff);
  // set the discretization state
  discret().set_state(0, "velocity", g_state().GetVelNp());
  // reset damping matrix
  damp().Zero();
  // add the stiffness matrix as well (also for the ApplyForce case!)
  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptr_);
  // set damping matrix
  eval_mat[1] = g_state().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::inertial_contributions(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  check_init_setup();

  if (TimInt().GetDataSDynPtr()->NeglectInertia()) return;

  // overwrite element action
  if (TimInt().GetDataSDyn().IsMassLumping())
    eval_data().SetActionType(DRT::ELEMENTS::struct_calc_nlnstifflmass);
  else
    eval_data().SetActionType(DRT::ELEMENTS::struct_calc_nlnstiffmass);

  // set the discretization state
  discret().set_state(0, "velocity", g_state().GetVelNp());
  discret().set_state(0, "acceleration", g_state().GetAccNp());
  // reset the mass matrix
  mass().Zero();
  // set mass matrix
  eval_mat[1] = g_state().GetMassMatrix();
  // set inertial vector if necessary
  eval_vec[1] = get_inertial_force();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::inertial_contributions(Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  check_init_setup();

  if (masslin_type_ == INPAR::STR::ml_none or TimInt().GetDataSDynPtr()->NeglectInertia()) return;

  // overwrite element action
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_internalinertiaforce);
  // set the discretization state
  discret().set_state(0, "velocity", g_state().GetVelNp());
  discret().set_state(0, "acceleration", g_state().GetAccNp());

  // set inertial vector if necessary
  eval_vec[1] = get_inertial_force();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::inertial_and_viscous_forces()
{
  check_init_setup();

  if (masslin_type_ == INPAR::STR::ml_none and !TimInt().GetDataSDynPtr()->NeglectInertia())
  {
    // calculate the inertial force at t_{n+1}
    mass().Multiply(false, *g_state().GetAccNp(), finertial_np());
  }

  // calculate the viscous/damping force at t_{n+1}
  if (eval_data().GetDampingType() != INPAR::STR::damp_none)
  {
    if (not damp().Filled()) damp().Complete();
    damp().Multiply(false, *g_state().GetVelNp(), fvisco_np());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::fill_complete()
{
  if (not stiff_ptr_->Filled()) stiff_ptr_->Complete();

  if (not mass().Filled()) mass().Complete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::rayleigh_damping_matrix()
{
  if (eval_data().GetDampingType() != INPAR::STR::damp_rayleigh) return;

  const double& dampk = TimInt().GetDataSDyn().get_damping_stiffness_factor();
  const double& dampm = TimInt().GetDataSDyn().get_damping_mass_factor();

  // damping matrix with initial stiffness
  damp().Add(stiff(), false, dampk, 0.0);
  damp().Add(mass(), false, dampm, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::MODELEVALUATOR::Structure::get_inertial_force()
{
  switch (masslin_type_)
  {
    case INPAR::STR::ml_rotations:
    case INPAR::STR::ml_standard:
    {
      finertial_np().PutScalar(0.0);
      // set inertial force
      return g_state().GetFinertialNp();
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
void STR::MODELEVALUATOR::Structure::init_output_runtime_structure()
{
  check_init();
  const auto discretization = Teuchos::rcp_dynamic_cast<const DRT::Discretization>(
      const_cast<STR::MODELEVALUATOR::Structure*>(this)->discret_ptr(), true);
  vtu_writer_ptr_ = Teuchos::rcp(
      new IO::DiscretizationVisualizationWriterMesh(discretization, visualization_params_,
          [](const DRT::Element* element)
          {
            // Skip beam elements which live in the same discretization but use a different output
            // mechanism
            return !dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element);
          }));

  if (g_in_output().get_runtime_output_params()->GetStructureParams()->gauss_point_data_output() !=
      INPAR::STR::GaussPointDataOutputType::none)
  {
    init_output_runtime_structure_gauss_point_data();
  }
}

void STR::MODELEVALUATOR::Structure::init_output_runtime_structure_gauss_point_data()
{
  // Set all parameters in the evaluation data container.
  eval_data().SetActionType(DRT::ELEMENTS::struct_init_gauss_point_data_output);
  eval_data().set_gauss_point_data_output_manager_ptr(Teuchos::rcp(new GaussPointDataOutputManager(
      g_in_output().get_runtime_output_params()->GetStructureParams()->gauss_point_data_output())));
  eval_data().SetTotalTime(g_state().GetTimeNp());
  eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);

  // Set vector values needed by elements.
  discret().ClearState();

  // Set dummy evaluation vectors and matrices.
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal(eval_mat.data(), eval_vec.data());

  eval_data().get_gauss_point_data_output_manager_ptr()->distribute_quantities(discret().Comm());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_time_step_output_runtime_structure() const
{
  check_init_setup();

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
  write_output_runtime_structure(disn_col, veln_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_iteration_output_runtime_structure() const
{
  check_init_setup();

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
  write_output_runtime_structure(disnp_col, velnp_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_output_runtime_structure(
    const Teuchos::RCP<Epetra_Vector>& displacement_state_vector,
    const Teuchos::RCP<Epetra_Vector>& velocity_state_vector, int timestep_number,
    double time) const
{
  check_init_setup();

  // get the parameter container object
  const DRT::ELEMENTS::StructureRuntimeOutputParams& structure_output_params =
      *GInOutput().get_runtime_output_params()->GetStructureParams();

  // reset time and time step of the writer object
  vtu_writer_ptr_->Reset();

  // append all desired output data to the writer object's storage

  // append displacement if desired
  if (structure_output_params.output_displacement_state())
    vtu_writer_ptr_->append_dof_based_result_data_vector(
        displacement_state_vector, 3, 0, "displacement");

  // append velocity if desired
  if (structure_output_params.OutputVelocityState())
    vtu_writer_ptr_->append_dof_based_result_data_vector(velocity_state_vector, 3, 0, "velocity");

  // append element owner if desired
  if (structure_output_params.OutputElementOwner())
    vtu_writer_ptr_->AppendElementOwner("element_owner");

  // append element GIDs if desired
  if (structure_output_params.OutputElementGID()) vtu_writer_ptr_->AppendElementGID("element_gid");

  // append element ghosting information if desired
  if (structure_output_params.output_element_ghosting())
    vtu_writer_ptr_->append_element_ghosting_information();

  // append node GIDs if desired
  if (structure_output_params.OutputNodeGID()) vtu_writer_ptr_->AppendNodeGID("node_gid");

  // append stress if desired
  if (structure_output_params.output_stress_strain() and
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
    vtu_writer_ptr_->append_node_based_result_data_vector(
        EvalData().get_stress_data_node_postprocessed(), 6, name_nodal);

    // Write element stress data.
    vtu_writer_ptr_->append_element_based_result_data_vector(
        EvalData().get_stress_data_element_postprocessed(), 6, name_element);
  }

  // append strain if desired.
  if (structure_output_params.output_stress_strain() and
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
    vtu_writer_ptr_->append_node_based_result_data_vector(
        EvalData().get_strain_data_node_postprocessed(), 6, name_nodal);

    // Write element strain data.
    vtu_writer_ptr_->append_element_based_result_data_vector(
        EvalData().get_strain_data_element_postprocessed(), 6, name_element);
  }

  // Add gauss point data if desired
  if (structure_output_params.gauss_point_data_output() !=
      INPAR::STR::GaussPointDataOutputType::none)
  {
    const GaussPointDataOutputManager& elementDataManager =
        *EvalData().get_gauss_point_data_output_manager_ptr();
    for (const auto& nameAndSize : elementDataManager.GetQuantities())
    {
      const std::string& name = nameAndSize.first;
      const int size = nameAndSize.second;

      switch (elementDataManager.GetOutputType())
      {
        case INPAR::STR::GaussPointDataOutputType::element_center:
        {
          Teuchos::RCP<Epetra_MultiVector> data =
              elementDataManager.get_element_center_data().at(name);
          vtu_writer_ptr_->append_element_based_result_data_vector(data, size, name);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::gauss_points:
        {
          const std::vector<Teuchos::RCP<Epetra_MultiVector>>& data_list =
              elementDataManager.GetGaussPointData().at(name);
          for (std::size_t gp = 0; gp < data_list.size(); ++gp)
          {
            const std::string name_with_gp = name + "_gp_" + std::to_string(gp);
            vtu_writer_ptr_->append_element_based_result_data_vector(
                data_list[gp], size, name_with_gp);
          }
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::nodes:
        {
          Teuchos::RCP<Epetra_MultiVector> data = elementDataManager.GetNodalData().at(name);
          vtu_writer_ptr_->append_node_based_result_data_vector(data, size, name);
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
void STR::MODELEVALUATOR::Structure::output_runtime_structure_postprocess_stress_strain()
{
  check_init_setup();

  if (not(g_in_output().GetStressOutputType() == INPAR::STR::stress_none and
          g_in_output().GetStrainOutputType() == INPAR::STR::strain_none))
  {
    // Set all parameters in the evaluation data container.
    eval_data().SetActionType(DRT::ELEMENTS::struct_calc_stress);
    eval_data().SetTotalTime(g_state().GetTimeNp());
    eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);
    eval_data().SetStressData(Teuchos::rcp(new std::vector<char>()));
    eval_data().set_coupling_stress_data(Teuchos::rcp(new std::vector<char>()));
    eval_data().SetStrainData(Teuchos::rcp(new std::vector<char>()));
    eval_data().set_plastic_strain_data(Teuchos::rcp(new std::vector<char>()));

    // Set vector values needed by elements.
    discret().ClearState();
    discret().set_state(0, "displacement", g_state().GetDisNp());
    discret().set_state(0, "residual displacement", dis_incr_ptr_);

    // GState().GetDisNp()->Print(std::cout);

    // Set dummy evaluation vectors and matrices.
    std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
        Teuchos::null, Teuchos::null, Teuchos::null};
    std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
        Teuchos::null, Teuchos::null};

    evaluate_internal_specified_elements(
        eval_mat.data(), eval_vec.data(), discret().ElementRowMap());

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
      for (int i = 0; i < discret_ptr()->ElementRowMap()->NumMyElements(); ++i)
      {
        if (DoPostprocessingOnElement(*discret().lRowElement(i)))
        {
          Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstress =
              Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
          CORE::COMM::ParObject::ExtractfromPack(position, raw_data, *gpstress);
          mapdata[discret_ptr()->ElementRowMap()->GID(i)] = gpstress;
        }
      }
      return mapdata;
    };

    auto PostprocessGaussPointDataToNodes =
        [&](const std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>& map_data,
            Epetra_MultiVector& assembled_data)
    {
      discret_ptr()->Evaluate(
          [&](DRT::Element& ele)
          {
            if (DoPostprocessingOnElement(ele))
              CORE::FE::ExtrapolateGaussPointQuantityToNodes(
                  ele, *map_data.at(ele.Id()), discret(), assembled_data);
          });
    };

    auto PostprocessGaussPointDataToElementCenter =
        [&](const std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>& map_data,
            Epetra_MultiVector& assembled_data)
    {
      discret_ptr()->Evaluate(
          [&](DRT::Element& ele)
          {
            if (DoPostprocessingOnElement(ele))
              CORE::FE::EvaluateGaussPointQuantityAtElementCenter(
                  ele, *map_data.at(ele.Id()), assembled_data);
          });
    };

    // Postprocess the result vectors.
    if (g_in_output().GetStressOutputType() != INPAR::STR::stress_none)
    {
      std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> gp_stress_data =
          EvaluateGaussPointData(*eval_data().GetStressData());

      CORE::COMM::Exporter ex(
          *(discret().ElementRowMap()), *(discret().ElementColMap()), discret().Comm());
      ex.Export(gp_stress_data);

      eval_data().get_stress_data_node_postprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret().NodeColMap(), 6, true));
      eval_data().get_stress_data_element_postprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret().ElementRowMap(), 6, true));


      Epetra_MultiVector row_nodal_data(*discret().NodeRowMap(), 6, true);
      PostprocessGaussPointDataToNodes(gp_stress_data, row_nodal_data);
      CORE::LINALG::Export(row_nodal_data, *eval_data().get_stress_data_node_postprocessed());

      PostprocessGaussPointDataToElementCenter(
          gp_stress_data, *eval_data().get_stress_data_element_postprocessed());
    }
    if (g_in_output().GetStrainOutputType() != INPAR::STR::strain_none)
    {
      std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>> gp_strain_data =
          EvaluateGaussPointData(*eval_data().GetStrainData());

      CORE::COMM::Exporter ex(
          *(discret().ElementRowMap()), *(discret().ElementColMap()), discret().Comm());
      ex.Export(gp_strain_data);

      eval_data().get_strain_data_node_postprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret().NodeColMap(), 6, true));
      eval_data().get_strain_data_element_postprocessed() =
          Teuchos::rcp(new Epetra_MultiVector(*discret().ElementRowMap(), 6, true));

      Epetra_MultiVector row_nodal_data(*discret().NodeRowMap(), 6, true);
      PostprocessGaussPointDataToNodes(gp_strain_data, row_nodal_data);
      CORE::LINALG::Export(row_nodal_data, *eval_data().get_strain_data_node_postprocessed());

      PostprocessGaussPointDataToElementCenter(
          gp_strain_data, *eval_data().get_strain_data_element_postprocessed());
    }
  }
}

void STR::MODELEVALUATOR::Structure::output_runtime_structure_gauss_point_data()
{
  const DRT::ELEMENTS::StructureRuntimeOutputParams& structure_output_params =
      *g_in_output().get_runtime_output_params()->GetStructureParams();
  if (structure_output_params.gauss_point_data_output() !=
      INPAR::STR::GaussPointDataOutputType::none)
  {
    check_init_setup();

    eval_data().SetActionType(DRT::ELEMENTS::struct_gauss_point_data_output);
    eval_data().SetTotalTime(g_state().GetTimeNp());
    eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);

    eval_data().gauss_point_data_output_manager_ptr()->PrepareData(
        *discret().NodeColMap(), *discret().ElementRowMap());

    discret().ClearState();
    discret().set_state(0, "displacement", g_state().GetDisNp());
    discret().set_state(0, "residual displacement", dis_incr_ptr_);

    std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
        Teuchos::null, Teuchos::null, Teuchos::null};
    std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
        Teuchos::null, Teuchos::null};

    evaluate_internal(eval_mat.data(), eval_vec.data());

    eval_data().gauss_point_data_output_manager_ptr()->post_evaluate();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::init_output_runtime_beams()
{
  beam_vtu_writer_ptr_ = Teuchos::rcp(
      new BeamDiscretizationRuntimeOutputWriter(visualization_params_, dis_np().Comm()));

  // get the parameter container object
  const DRT::ELEMENTS::BeamRuntimeOutputParams& beam_output_params =
      *g_in_output().get_runtime_output_params()->GetBeamParams();

  // export displacement state to column format
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret().DofColMap(), true));
  CORE::LINALG::Export(*g_state().GetDisN(), *disn_col);

  // get bounding box object only if periodic boundaries are active
  Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> bounding_box_ptr =
      TimInt().GetDataSDynPtr()->get_periodic_bounding_box();

  // initialize the writer object with current displacement state
  beam_vtu_writer_ptr_->Initialize(const_cast<STR::MODELEVALUATOR::Structure*>(this)->discret_ptr(),
      beam_output_params.use_absolute_positions(),
      beam_output_params.get_number_visualization_subsegments(), bounding_box_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_time_step_output_runtime_beams() const
{
  check_init_setup();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disn_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisN(), *disn_col);

  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, GState().GetTimeN(), GState().GetStepN());
  write_output_runtime_beams(disn_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_iteration_output_runtime_beams() const
{
  check_init_setup();

  // export displacement state to column format
  const auto& discret = dynamic_cast<const DRT::Discretization&>(Discret());
  Teuchos::RCP<Epetra_Vector> disnp_col =
      Teuchos::rcp(new Epetra_Vector(*discret.DofColMap(), true));
  CORE::LINALG::Export(*GState().GetDisNp(), *disnp_col);

  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, GState().GetTimeN(), GState().GetStepN(), EvalData().GetNlnIter());
  write_output_runtime_beams(disnp_col, output_step, output_time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_output_runtime_beams(
    const Teuchos::RCP<Epetra_Vector>& displacement_state_vector, int timestep_number,
    double time) const
{
  check_init_setup();

  // get the parameter container object
  const DRT::ELEMENTS::BeamRuntimeOutputParams& beam_output_params =
      *GInOutput().get_runtime_output_params()->GetBeamParams();

  // set geometry
  beam_vtu_writer_ptr_->set_geometry_from_beam_discretization(displacement_state_vector);

  // append all desired output data to the writer object's storage
  beam_vtu_writer_ptr_->append_element_owning_processor();

  // append beam radius
  beam_vtu_writer_ptr_->append_element_circular_cross_section_radius();

  // append displacement if desired
  if (beam_output_params.is_write_internal_energy_element())
    beam_vtu_writer_ptr_->append_element_internal_energy();

  // append displacement if desired
  if (beam_output_params.is_write_kinetic_energy_element())
    beam_vtu_writer_ptr_->append_element_kinetic_energy();

  // append displacement if desired
  if (beam_output_params.output_displacement_state())
    beam_vtu_writer_ptr_->append_displacement_field(displacement_state_vector);

  // append triads if desired
  if (beam_output_params.is_write_triad_visualization_points())
    beam_vtu_writer_ptr_->AppendTriadField(displacement_state_vector);

  // append material cross-section strain resultants if desired
  if (beam_output_params.is_write_material_strains_gauss_points())
    beam_vtu_writer_ptr_->append_gauss_point_material_cross_section_strain_resultants();

  // append material cross-section strain resultants if desired
  if (beam_output_params.is_write_material_strains_continuous())
    beam_vtu_writer_ptr_->append_gauss_point_material_cross_section_strain_resultants_continuous();

  // append material cross-section stress resultants if desired
  if (beam_output_params.is_write_material_stresses_gauss_points())
    beam_vtu_writer_ptr_->append_gauss_point_material_cross_section_stress_resultants();

  // append material cross-section stress resultants if desired
  if (beam_output_params.is_write_material_stress_continuous())
    beam_vtu_writer_ptr_->append_gauss_point_material_cross_section_stress_resultants_continuous();

  // append filament id and type if desired
  if (beam_output_params.is_write_element_filament_condition())
    beam_vtu_writer_ptr_->append_element_filament_id_and_type();

  // append filament id and type if desired
  if (beam_output_params.is_write_orientation_paramter())
    beam_vtu_writer_ptr_->append_element_orientation_paramater(displacement_state_vector);

  // append reference length if desired.
  if (beam_output_params.IsWriteRefLength()) beam_vtu_writer_ptr_->AppendRefLength();

  // export displacement state to column format
  if (beam_output_params.is_write_rve_crosssection_forces())
    beam_vtu_writer_ptr_->append_rve_crosssection_forces(displacement_state_vector);

  // export beam element IDs
  if (beam_output_params.IsWriteElementGID()) beam_vtu_writer_ptr_->AppendElementGID();

  // Ghosting information
  if (beam_output_params.is_write_element_ghosting())
    beam_vtu_writer_ptr_->append_element_ghosting_information();

  // finalize everything and write all required VTU files to filesystem
  beam_vtu_writer_ptr_->WriteToDisk(time, timestep_number);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::evaluate_internal(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec)
{
  pre_evaluate_internal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", eval_data_ptr());

  evaluate_internal(p, eval_mat, eval_vec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::evaluate_internal(Teuchos::ParameterList& p,
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
  params_interface2_parameter_list(eval_data_ptr(), p);

  discret().Evaluate(p, eval_mat[0], eval_mat[1], eval_vec[0], eval_vec[1], eval_vec[2]);
  discret().ClearState();
}

void STR::MODELEVALUATOR::Structure::evaluate_internal_specified_elements(
    Teuchos::RCP<CORE::LINALG::SparseOperator>* eval_mat, Teuchos::RCP<Epetra_Vector>* eval_vec,
    const Epetra_Map* ele_map_to_be_evaluated)
{
  pre_evaluate_internal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", eval_data_ptr());

  evaluate_internal_specified_elements(p, eval_mat, eval_vec, ele_map_to_be_evaluated);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::evaluate_internal_specified_elements(Teuchos::ParameterList& p,
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
  params_interface2_parameter_list(eval_data_ptr(), p);

  DRT::UTILS::Evaluate(*discret_ptr(), p, *eval_mat, *eval_vec, ele_map_to_be_evaluated);

  discret().ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::evaluate_neumann(const Teuchos::RCP<Epetra_Vector>& eval_vec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& eval_mat)
{
  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", eval_data_ptr());
  evaluate_neumann(p, eval_vec, eval_mat);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::evaluate_neumann(Teuchos::ParameterList& p,
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
  discret().evaluate_neumann(p, eval_vec, eval_mat);
  discret().ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::write_restart(
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
void STR::MODELEVALUATOR::Structure::read_restart(IO::DiscretizationReader& ioreader)
{
  check_init_setup();
  // read structural force vector
  ioreader.ReadVector(g_state().GetFstructureOld(), "fstructure_old");
  ioreader.ReadVector(g_state().GetFintN(), "fint");
  // read displacement field
  Teuchos::RCP<Epetra_Vector>& disnp = g_state().GetDisNp();
  ioreader.ReadVector(disnp, "displacement");
  g_state().GetMultiDis()->UpdateSteps(*disnp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::Predict(const INPAR::STR::PredEnum& pred_type)
{
  // set the element action
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_predict);
  eval_data().SetPredictorType(pred_type);

  // set the matrix and vector pointers to Teuchos::null
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::run_pre_compute_x(
    const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
{
  check_init_setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::RunRecover()
{
  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "residual displacement", dis_incr_ptr_);
  discret().set_state(0, "displacement", g_state().GetDisNp());
  // set the element action
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_recover);
  // set the matrix and vector pointers to Teuchos::null
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::run_post_compute_x(
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  check_init_setup();
  Reset(xnew);
  /* set the class internal displacement increment vector. Check if it is
   * meaningful/necessary in some cases, like incremental strains etc. */
  dis_incr_ptr_ = g_state().ExtractDisplEntries(dir);
  RunRecover();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::run_post_iterate(const ::NOX::Solver::Generic& solver)
{
  check_init_setup();

  if (vtu_writer_ptr_ != Teuchos::null and
      g_in_output().get_runtime_output_params()->output_every_iteration())
  {
    output_runtime_structure_postprocess_stress_strain();
    output_runtime_structure_gauss_point_data();
    write_iteration_output_runtime_structure();
  }

  // write special output for beams if desired
  if (beam_vtu_writer_ptr_ != Teuchos::null and
      g_in_output().get_runtime_output_params()->output_every_iteration())
    write_iteration_output_runtime_beams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepState(const double& timefac_n)
{
  check_init_setup();
  // update state
  // new displacements at t_{n+1} -> t_n
  //    D_{n} := D_{n+1}
  g_state().GetMultiDis()->UpdateSteps(dis_np());

  // new velocities at t_{n+1} -> t_{n}
  //    V_{n} := V_{n+1}
  g_state().GetMultiVel()->UpdateSteps(*g_state().GetVelNp());

  // new at t_{n+1} -> t_n
  //    A_{n} := A_{n+1}
  g_state().GetMultiAcc()->UpdateSteps(*g_state().GetAccNp());

  // store the old external force
  g_state().GetFextN()->Scale(1.0, fext_np());

  // store the old reaction force
  g_state().GetFreactN()->Scale(1.0, *g_state().GetFreactNp());

  // store the old internal force
  g_state().GetFintN()->Scale(1.0, fint_np());

  // new at t_{n+1} -> t_{n+timefac_n}
  //    F^{struct}_{n+timefac_n} := timefac_n * F^{struct}_{n+1}
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr = g_state().GetFstructureOld();
  fstructold_ptr->Update(timefac_n, fint_np(), 1.0);
  fstructold_ptr->Update(-timefac_n, fext_np(), 1.0);

  // set the displacement increment back to zero
  dis_incr_ptr_->PutScalar(0.0);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::evaluate_jacobian_contributions_from_element_level_for_ptc()
{
  check_init_setup();
  // currently a fixed number of matrix and vector pointers are supported
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_addjacPTC);

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", g_state().GetDisNp());

  eval_mat[0] = Teuchos::rcpFromRef(*stiff_ptc_ptr_);

  // evaluate
  evaluate_internal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::assemble_jacobian_contributions_from_element_level_for_ptc(
    Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac, const double& timefac_n)
{
  g_state().AssignModelBlock(*modjac, stiff_ptc(), Type(), STR::MatBlockType::displ_displ);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::UpdateStepElement()
{
  check_init_setup();
  // other parameters that might be needed by the elements
  eval_data().SetTotalTime(g_state().GetTimeNp());
  eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);

  const INPAR::STR::PreStress prestress_type = TimInt().GetDataSDyn().GetPreStressType();
  const double prestress_time = TimInt().GetDataSDyn().GetPreStressTime();
  bool isDuringPrestressing = prestress_type != INPAR::STR::PreStress::none &&
                              g_state().GetTimeN() <= prestress_time + 1.0e-15;

  if (isDuringPrestressing && prestress_type == INPAR::STR::PreStress::mulf)
  {
    if (discret().Comm().MyPID() == 0)
      IO::cout << "====== Entering PRESTRESSING update" << IO::endl;

    // Choose special update action for elements in case of MULF
    eval_data().SetActionType(DRT::ELEMENTS::struct_update_prestress);
  }
  else
  {
    // Call the normal element  update routine
    eval_data().SetActionType(DRT::ELEMENTS::struct_calc_update_istep);
  }


  // go to elements
  discret().ClearState();
  discret().set_state("displacement", g_state().GetDisN());

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};
  evaluate_internal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::update_residual()
{
  check_init_setup();
  dis_incr_ptr_->Update(-1.0, *g_state().GetDisN(), 1.0, *g_state().GetDisNp(), 0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::determine_stress_strain()
{
  check_init_setup();

  if (g_in_output().GetStressOutputType() == INPAR::STR::stress_none and
      g_in_output().get_coupling_stress_output_type() == INPAR::STR::stress_none and
      g_in_output().GetStrainOutputType() == INPAR::STR::strain_none and
      g_in_output().get_plastic_strain_output_type() == INPAR::STR::strain_none)
    return;

  // set all parameters in the evaluation data container
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_stress);
  eval_data().SetTotalTime(g_state().GetTimeNp());
  eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);
  eval_data().SetStressData(Teuchos::rcp(new std::vector<char>()));
  eval_data().set_coupling_stress_data(Teuchos::rcp(new std::vector<char>()));
  eval_data().SetStrainData(Teuchos::rcp(new std::vector<char>()));
  eval_data().set_plastic_strain_data(Teuchos::rcp(new std::vector<char>()));

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", g_state().GetDisNp());
  discret().set_state(0, "residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal_specified_elements(eval_mat.data(), eval_vec.data(), discret().ElementRowMap());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::determine_strain_energy(
    const Epetra_Vector& disnp, const bool global)
{
  check_init_setup();

  // set required parameters in the evaluation data container
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_energy);
  eval_data().SetTotalTime(g_state().GetTimeNp());
  eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);

  // set state vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", Teuchos::rcpFromRef(disnp));
  discret().set_state(0, "residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  pre_evaluate_internal();

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", eval_data_ptr());

  // evaluate energy contributions on element level (row elements only)
  evaluate_internal_specified_elements(
      p, eval_mat.data(), eval_vec.data(), discret().ElementRowMap());

  if (global)
  {
    double my_int_energy = eval_data().GetEnergyData(STR::internal_energy);
    double gsum = 0.0;
    discret().Comm().SumAll(&my_int_energy, &gsum, 1);

    eval_data().set_value_for_energy_type(gsum, STR::internal_energy);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineEnergy()
{
  DetermineEnergy(*g_state().GetDisNp(), g_state().GetVelNp().get(), false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::DetermineEnergy(
    const Epetra_Vector& disnp, const Epetra_Vector* velnp, const bool global)
{
  determine_strain_energy(disnp, global);

  // global calculation of kinetic energy
  if (masslin_type_ == INPAR::STR::ml_none and velnp != nullptr)
  {
    double kinetic_energy_times2 = 0.0;

    Teuchos::RCP<Epetra_Vector> linear_momentum =
        CORE::LINALG::CreateVector(*g_state().DofRowMapView(), true);

    mass().Multiply(false, *velnp, *linear_momentum);

    linear_momentum->Dot(*velnp, &kinetic_energy_times2);

    // only add the result on one processor because we sum over all procs later
    if (global or g_state().GetMyRank() == 0)
    {
      eval_data().add_contribution_to_energy_type(0.5 * kinetic_energy_times2, STR::kinetic_energy);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::determine_optional_quantity()
{
  check_init_setup();

  switch (g_in_output().get_opt_quantity_output_type())
  {
    case INPAR::STR::optquantity_none:
    {
      // do nothing and return
      return;
    }
    case INPAR::STR::optquantity_membranethickness:
    {
      // evaluate thickness of membrane finite elements
      eval_data().SetActionType(DRT::ELEMENTS::struct_calc_thickness);
      break;
    }
    default:
      FOUR_C_THROW("Type of optional quantity not implemented yet!");
  }

  // set all parameters in the evaluation data container
  eval_data().SetTotalTime(g_state().GetTimeNp());
  eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);
  eval_data().SetOptQuantityData(Teuchos::rcp(new std::vector<char>()));

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", g_state().GetDisNp());
  discret().set_state(0, "residual displacement", dis_incr_ptr_);

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Structure::determine_element_volumes(
    const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols)
{
  // set action in params-interface
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_mass_volume);

  Teuchos::ParameterList p;
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", eval_data_ptr());

  // set vector values needed by elements
  discret().ClearState();
  Teuchos::RCP<const Epetra_Vector> disnp = g_state().ExtractDisplEntries(x);
  discret().set_state(0, "displacement", disnp);

  // start evaluation
  const Epetra_Map* relemap = discret().ElementRowMap();
  ele_vols = Teuchos::rcp(new Epetra_Vector(*relemap, true));
  const unsigned my_num_reles = relemap->NumMyElements();

  DRT::Element::LocationArray la(discret().NumDofSets());
  CORE::LINALG::SerialDenseVector ele_vol(6, true);

  CORE::LINALG::SerialDenseMatrix empty_dummy_mat;
  CORE::LINALG::SerialDenseVector empty_dummy_vec;

  STR::ELEMENTS::EvalErrorFlag ele_eval_error = STR::ELEMENTS::ele_error_none;
  for (unsigned elid = 0; elid < my_num_reles; ++elid)
  {
    DRT::Element* rele = discret().lRowElement(elid);
    rele->LocationVector(discret(), la, false);

    eval_data().SetActionType(DRT::ELEMENTS::analyse_jacobian_determinant);
    rele->Evaluate(p, discret(), la, empty_dummy_mat, empty_dummy_mat, ele_vol, empty_dummy_vec,
        empty_dummy_vec);

    if (not eval_data().IsEleEvalError())
    {
      eval_data().SetActionType(DRT::ELEMENTS::struct_calc_mass_volume);
      rele->Evaluate(p, discret(), la, empty_dummy_mat, empty_dummy_mat, ele_vol, empty_dummy_vec,
          empty_dummy_vec);
    }

    // set a negative value, if the evaluation failed
    if (eval_data().IsEleEvalError())
    {
      ele_eval_error = eval_data().GetEleEvalErrorFlag();
      // reset for the next element
      eval_data().SetEleEvalErrorFlag(STR::ELEMENTS::ele_error_none);
      ele_vol(2) = -1.0;
    }

    const int rele_lid = relemap->LID(rele->Id());
    (*ele_vols)[rele_lid] = ele_vol(2);
    ele_vol.putScalar(0.0);
  }

  discret().ClearState();
  eval_data().SetEleEvalErrorFlag(ele_eval_error);

  return eval_error_check();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();

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
void STR::MODELEVALUATOR::Structure::runtime_pre_output_step_state()
{
  check_init_setup();

  if (vtu_writer_ptr_ != Teuchos::null)
  {
    output_runtime_structure_postprocess_stress_strain();
    output_runtime_structure_gauss_point_data();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::runtime_output_step_state() const
{
  check_init_setup();

  if (vtu_writer_ptr_ != Teuchos::null) write_time_step_output_runtime_structure();

  // write special output for beams if desired
  if (beam_vtu_writer_ptr_ != Teuchos::null) write_time_step_output_runtime_beams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::ResetStepState()
{
  check_init_setup();

  // reset disp, vel, acc state vector
  g_state_ptr()->GetDisNp()->Update(1.0, (*g_state_ptr()->GetDisN()), 0.0);
  g_state_ptr()->GetVelNp()->Update(1.0, (*g_state_ptr()->GetVelN()), 0.0);
  g_state_ptr()->GetAccNp()->Update(1.0, (*g_state_ptr()->GetAccN()), 0.0);

  // other parameters that might be needed by the elements
  eval_data().SetTotalTime(g_state().GetTimeNp());
  eval_data().SetDeltaTime((*g_state().GetDeltaTime())[0]);
  // action for elements
  eval_data().SetActionType(DRT::ELEMENTS::struct_calc_reset_istep);

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};
  evaluate_internal(eval_mat.data(), eval_vec.data());

  discret_ptr()->ClearState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Structure::get_block_dof_row_map_ptr() const
{
  check_init_setup();
  return GState().dof_row_map();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::get_current_solution_ptr() const
{
  check_init();
  return GState().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Structure::get_last_time_step_solution_ptr()
    const
{
  check_init();
  return GState().GetDisN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::PostOutput()
{
  check_init_setup();
  // empty
}  // PostOutput()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::fint_np()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetFintNp().is_null(), "nullptr!");

  return *g_state().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::fint_np() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetFintNp().is_null(), "nullptr!");

  return *GState().GetFintNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::fint_n() const
{
  check_init();
  if (GState().GetFintN().is_null()) FOUR_C_THROW("NULL pointer!");

  return *GState().GetFintN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::fext_np()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetFextNp().is_null(), "nullptr!");

  return *g_state().GetFextNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::fext_np() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetFextNp().is_null(), "nullptr!");

  return *GState().GetFextNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::fext_n() const
{
  check_init();
  if (GState().GetFextN().is_null()) FOUR_C_THROW("NULL pointer!");

  return *GState().GetFextN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::finertial_np()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetFinertialNp().is_null(), "nullptr!");

  return *g_state().GetFinertialNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::finertial_np() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetFinertialNp().is_null(), "nullptr!");

  return *GState().GetFinertialNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::fvisco_np()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetFviscoNp().is_null(), "nullptr!");

  return *g_state().GetFviscoNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::fvisco_np() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetFviscoNp().is_null(), "nullptr!");

  return *GState().GetFviscoNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& STR::MODELEVALUATOR::Structure::dis_np()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetDisNp().is_null(), "nullptr!");

  return *g_state().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& STR::MODELEVALUATOR::Structure::dis_np() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetDisNp().is_null(), "nullptr!");

  return *GState().GetDisNp();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix& STR::MODELEVALUATOR::Structure::stiff() const
{
  check_init();
  FOUR_C_ASSERT(stiff_ptr_, "nullptr!");

  return *stiff_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix& STR::MODELEVALUATOR::Structure::stiff_ptc() const
{
  check_init();
  FOUR_C_ASSERT(stiff_ptc_ptr_ != Teuchos::null, "nullptr!");

  return *stiff_ptc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::mass()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetMassMatrix().is_null(), "nullptr!");

  return *g_state().GetMassMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::mass() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetMassMatrix().is_null(), "nullptr!");

  return *GState().GetMassMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::damp()
{
  check_init();
  FOUR_C_ASSERT(!g_state().GetDampMatrix().is_null(), "nullptr!");

  return *g_state().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CORE::LINALG::SparseOperator& STR::MODELEVALUATOR::Structure::damp() const
{
  check_init();
  FOUR_C_ASSERT(!GState().GetDampMatrix().is_null(), "nullptr!");

  return *GState().GetDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::params_interface2_parameter_list(
    Teuchos::RCP<STR::MODELEVALUATOR::Data> interface_ptr, Teuchos::ParameterList& params)
{
  FOUR_C_ASSERT(interface_ptr != Teuchos::null, "params_interface pointer not set");

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
  params.set<Teuchos::RCP<std::vector<char>>>("plstrain", interface_ptr->plastic_strain_data_ptr());
  params.set<Teuchos::RCP<std::vector<char>>>("optquantity", interface_ptr->OptQuantityDataPtr());
  params.set<int>("iostress", (int)interface_ptr->GetStressOutputType());
  params.set<int>("iostrain", (int)interface_ptr->GetStrainOutputType());
  params.set<int>("ioplstrain", (int)interface_ptr->get_plastic_strain_output_type());
  params.set<int>("iooptquantity", (int)interface_ptr->get_opt_quantity_output_type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::CreateBackupState(const Epetra_Vector& dir)
{
  check_init_setup();

  // set all parameters in the evaluation data container
  eval_data().SetActionType(DRT::ELEMENTS::struct_create_backup);

  // set vector values needed by elements
  discret().ClearState();
  discret().set_state(0, "displacement", g_state().GetDisNp());
  Teuchos::RCP<const Epetra_Vector> dir_displ = g_state().ExtractDisplEntries(dir);
  discret().set_state(0, "residual displacement", dir_displ);

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal(eval_mat.data(), eval_vec.data());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Structure::recover_from_backup_state()
{
  check_init_setup();

  // set all parameters in the evaluation data container
  eval_data().SetActionType(DRT::ELEMENTS::struct_recover_from_backup);

  // set vector values needed by elements
  discret().ClearState();

  // set dummy evaluation vectors and matrices
  std::array<Teuchos::RCP<Epetra_Vector>, 3> eval_vec = {
      Teuchos::null, Teuchos::null, Teuchos::null};
  std::array<Teuchos::RCP<CORE::LINALG::SparseOperator>, 2> eval_mat = {
      Teuchos::null, Teuchos::null};

  evaluate_internal(eval_mat.data(), eval_vec.data());
}

FOUR_C_NAMESPACE_CLOSE
