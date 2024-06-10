/*--------------------------------------------------------------------------*/
/*! \file

\brief base class for all elastohydrodynamic lubrication (lubrication structure interaction)
algorithms

\level 3

*/
/*--------------------------------------------------------------------------*/

#include "4C_ehl_base.hpp"

#include "4C_adapter_coupling_ehl_mortar.hpp"
#include "4C_adapter_lubrication.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_node.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_matchingoctree.hpp"
#include "4C_ehl_partitioned.hpp"
#include "4C_ehl_utils.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_lubrication_timint_implicit.hpp"
#include "4C_mat_lubrication_mat.hpp"

#include <Epetra_MultiVector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                     (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Base::Base(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string lubrication_disname)
    : AlgorithmBase(comm, globaltimeparams),
      structure_(Teuchos::null),
      lubrication_(Teuchos::null),
      fieldcoupling_(Core::UTILS::IntegralValue<Inpar::EHL::FieldCoupling>(
          Global::Problem::Instance()->elasto_hydro_dynamic_params(), "FIELDCOUPLING")),
      dry_contact_(Core::UTILS::IntegralValue<bool>(
          Global::Problem::Instance()->elasto_hydro_dynamic_params(), "DRY_CONTACT_MODEL"))
{
  Global::Problem* problem = Global::Problem::Instance();

  // get the solver number used for Lubrication solver
  const int linsolvernumber = lubricationparams.get<int>("LINEAR_SOLVER");

  // 2.- Setup discretizations and coupling.
  setup_discretizations(comm, struct_disname, lubrication_disname);

  setup_field_coupling(struct_disname, lubrication_disname);

  // 3.- Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<Core::FE::Discretization> structdis =
      Global::Problem::Instance()->GetDis(struct_disname);

  // set moving grid
  bool isale = true;

  // determine which time params to use to build the single fields
  // in case of time stepping time params have to be read from single field sections
  // in case of equal timestep size for all fields the time params are controlled solely
  // by the problem section (e.g. ehl or cell dynamic)
  const Teuchos::ParameterList* structtimeparams = &globaltimeparams;
  const Teuchos::ParameterList* lubricationtimeparams = &globaltimeparams;
  if (Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->elasto_hydro_dynamic_params(), "DIFFTIMESTEPSIZE"))
  {
    structtimeparams = &structparams;
    lubricationtimeparams = &lubricationparams;
  }

  Teuchos::RCP<Adapter::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new Adapter::StructureBaseAlgorithm(
          *structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<Adapter::Structure>(structure->structure_field());
  structure_->Setup();
  lubrication_ = Teuchos::rcp(new Adapter::LubricationBaseAlgorithm());
  lubrication_->Setup(*lubricationtimeparams, lubricationparams,
      problem->SolverParams(linsolvernumber), lubrication_disname, isale);
  mortaradapter_->store_dirichlet_status(structure_field()->GetDBCMapExtractor());

  // Structure displacement at the lubricated interface
  Teuchos::RCP<Epetra_Vector> disp = Core::LinAlg::CreateVector(*(structdis->dof_row_map()), true);

  mortaradapter_->Integrate(disp, Dt());
  // the film thickness initialization for very first time step
  heightold_ = mortaradapter_->Nodal_Gap();
}

/*----------------------------------------------------------------------*
 | read restart information for given time step   (public) wirtz 12/15  |
 *----------------------------------------------------------------------*/
void EHL::Base::read_restart(int restart)
{
  if (restart)
  {
    lubrication_->LubricationField()->read_restart(restart);
    structure_->read_restart(restart);
    SetTimeStep(structure_->TimeOld(), restart);

    mortaradapter_->Interface()->set_state(Mortar::state_old_displacement, *structure_->Dispn());
    mortaradapter_->Interface()->set_state(Mortar::state_new_displacement, *structure_->Dispn());
    mortaradapter_->Interface()->evaluate_nodal_normals();
    mortaradapter_->Interface()->export_nodal_normals();
    mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::n_old);
    mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::dm);
    mortaradapter_->Integrate(structure_->Dispnp(), Dt());
    heightold_ = mortaradapter_->Nodal_Gap();

    Core::IO::DiscretizationReader reader(lubrication_->LubricationField()->discretization(),
        Global::Problem::Instance()->InputControlFile(), restart);
    mortaradapter_->read_restart(reader);
  }
}

/*----------------------------------------------------------------------*
 | calculate velocities by a FD approximation               wirtz 12/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> EHL::Base::calc_velocity(Teuchos::RCP<const Epetra_Vector> dispnp)
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector(*(structure_->Dispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / Dt(), *dispnp, -1. / Dt());

  return vel;
}  // calc_velocity()

/*----------------------------------------------------------------------*
 | test results (if necessary)                     (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::TestResults(const Epetra_Comm& comm)
{
  Global::Problem* problem = Global::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(lubrication_->create_lubrication_field_test());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::setup_discretizations(const Epetra_Comm& comm, const std::string struct_disname,
    const std::string lubrication_disname)
{
  // Scheme   : the structure discretization is received from the input. Then, an ale-lubrication
  // disc. is cloned.

  Global::Problem* problem = Global::Problem::Instance();

  // 1.-Initialization.
  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<Core::FE::Discretization> lubricationdis = problem->GetDis(lubrication_disname);
  if (!structdis->Filled()) structdis->fill_complete();
  if (!lubricationdis->Filled()) lubricationdis->fill_complete();

  // first call fill_complete for single discretizations.
  // This way the physical dofs are numbered successively
  structdis->fill_complete();
  lubricationdis->fill_complete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_lubrication = lubricationdis->NumDof(0, lubricationdis->lRowNode(0));
  const int ndofperelement_lubrication = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;

  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux_lubrication =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
          ndofpernode_lubrication, ndofperelement_lubrication, 0, true));
  if (structdis->AddDofSet(dofsetaux_lubrication) != 1)
    FOUR_C_THROW("unexpected dof sets in structure field");

  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux_struct =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
          ndofpernode_struct, ndofperelement_struct, 0, true));
  if (lubricationdis->AddDofSet(dofsetaux_struct) != 1)
    FOUR_C_THROW("unexpected dof sets in lubrication field");

  // call assign_degrees_of_freedom also for auxiliary dofsets
  // note: the order of fill_complete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. lubrication dofs
  // 3. structure auxiliary dofs
  // 4. lubrication auxiliary dofs
  structdis->fill_complete(true, false, false);
  lubricationdis->fill_complete(true, false, false);
}

/*----------------------------------------------------------------------*
 | set structure solution on lubrication field              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::set_struct_solution(Teuchos::RCP<const Epetra_Vector> disp)
{
  //---------------------------------------------------------
  // 1. Update the Mortar Coupling
  //---------------------------------------------------------

  //  //Extract the structure displacement at the lubricated interface
  //  Teuchos::RCP<Epetra_Vector> idisp = Core::LinAlg::CreateVector(*(mergedrowmapextr_->Map(0)),
  //  true);//Structure displacement at the lubricated interface
  //  mergedrowmapextr_->ExtractVector(disp,0,idisp);
  // Reevalute the mortar martices D and M
  mortaradapter_->Integrate(disp, Dt());

  // Displace the mesh of the lubrication field in accordance with the slave-side interface
  set_mesh_disp(disp);

  // Calculate the average tangential fractions of the structure velocities at the interface and
  // provide them to the lubrication field
  set_average_velocity_field();

  // Calculate the relative tangential fractions of the structure velocities at the interface and
  // provide them to the lubrication field
  set_relative_velocity_field();

  // Provide the gap at the interface
  set_height_field();

  // provide the heightdot (time derivative of the gap)
  set_height_dot();

  // Create DBC map for unprojectable nodes
  setup_unprojectable_dbc();

  return;
}

/*----------------------------------------------------------------------*
 | calc tractions, resulting from fluid (pressure and viscous) seitz 01/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> EHL::Base::EvaluateFluidForce(
    Teuchos::RCP<const Epetra_Vector> pressure)
{
  // safety: unprojectable nodes to zero pressure
  if (inf_gap_toggle_lub_ != Teuchos::null)
    for (int i = 0; i < lubrication_->LubricationField()->Prenp()->Map().NumMyElements(); ++i)
    {
      if (abs(inf_gap_toggle_lub_->operator[](inf_gap_toggle_lub_->Map().LID(
                  lubrication_->LubricationField()->Prenp()->Map().GID(i))) -
              1) < 1.e-2)
        lubrication_->LubricationField()->Prenp()->operator[](i) = 0.;
    }

  // Forces on the interfaces due to the fluid traction
  Teuchos::RCP<Epetra_Vector> slaveiforce =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMortarMatrixD()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> masteriforce =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMortarMatrixM()->DomainMap()));

  stritraction_D_ = Teuchos::rcp(new Epetra_Vector(*ada_strDisp_to_lubDisp_->MasterDofMap()));
  stritraction_M_ = Teuchos::rcp(new Epetra_Vector(*ada_strDisp_to_lubDisp_->MasterDofMap()));

  // add pressure force
  add_pressure_force(slaveiforce, masteriforce);
  // add poiseuille flow force
  add_poiseuille_force(slaveiforce, masteriforce);
  // add couette flow force
  add_couette_force(slaveiforce, masteriforce);

  // External force vector (global)
  Teuchos::RCP<Epetra_Vector> strforce =
      Teuchos::rcp(new Epetra_Vector(*(structure_->dof_row_map())));

  // Insert both interface forces into the global force vector
  slaverowmapextr_->InsertVector(slaveiforce, 0, strforce);
  masterrowmapextr_->InsertVector(masteriforce, 0, strforce);

  return strforce;
}

/*----------------------------------------------------------------------*
 | set tractions, resulting from lubrication pressure       seitz 01/18 |
 *----------------------------------------------------------------------*/
void EHL::Base::set_lubrication_solution(Teuchos::RCP<const Epetra_Vector> pressure)
{
  // Provide the structure field with the force vector
  // Note that the mid-point values (gen-alpha) of the interface forces are evaluated in
  // STR::TimIntGenAlpha::evaluate_force_residual()
  structure_->SetForceInterface(EvaluateFluidForce(pressure));
}

void EHL::Base::add_pressure_force(
    Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce)
{
  Teuchos::RCP<Epetra_Vector> stritraction;

  Teuchos::RCP<Epetra_Vector> p_full =
      Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->dof_row_map(1)));
  if (lubrimaptransform_->Apply(*lubrication_->LubricationField()->Prenp(), *p_full))
    FOUR_C_THROW("apply failed");
  Teuchos::RCP<Epetra_Vector> p_exp =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  p_exp = ada_strDisp_to_lubDisp_->SlaveToMaster(p_full);
  stritraction = Teuchos::rcp(new Epetra_Vector(*mortaradapter_->Normals()));
  stritraction->Multiply(-1., *mortaradapter_->Normals(), *p_exp, 0.);

  // Get the Mortar D and M Matrix
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortard = mortaradapter_->GetMortarMatrixD();
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortarm = mortaradapter_->GetMortarMatrixM();

  // f_slave = D^T*t
  int err = mortard->Multiply(true, *stritraction, *slaveiforce);
  if (err != 0) FOUR_C_THROW("error while calculating slave side interface force");
  if (stritraction_D_->Update(1., *stritraction, 1.)) FOUR_C_THROW("Update failed");

  // f_master = -M^T*t
  err = mortarm->Multiply(true, *stritraction, *masteriforce);
  if (err != 0) FOUR_C_THROW("error while calculating master side interface force");
  masteriforce->Scale(-1.0);
  if (stritraction_M_->Update(-1., *stritraction, 1.)) FOUR_C_THROW("update failed");
}

void EHL::Base::add_poiseuille_force(
    Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce)
{
  // poiseuille flow forces
  Teuchos::RCP<Epetra_Vector> p_int =
      ada_strDisp_to_lubPres_->SlaveToMaster(lubrication_->LubricationField()->Prenp());
  Teuchos::RCP<Epetra_Vector> p_int_full =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  Core::LinAlg::Export(*p_int, *p_int_full);

  Teuchos::RCP<Epetra_Vector> nodal_gap =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *nodal_gap))
    FOUR_C_THROW("multiply failed");

  Core::LinAlg::SparseMatrix m(*mortaradapter_->SurfGradMatrix());

  m.LeftScale(*nodal_gap);
  m.Scale(-.5);

  Teuchos::RCP<Epetra_Vector> poiseuille_force =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  m.Apply(*p_int_full, *poiseuille_force);

  Teuchos::RCP<Epetra_Vector> slave_psl =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMortarMatrixD()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> master_psl =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMortarMatrixM()->DomainMap()));

  // f_slave = D^T*t
  if (mortaradapter_->GetMortarMatrixD()->Multiply(true, *poiseuille_force, *slave_psl))
    FOUR_C_THROW("Multiply failed");
  if (stritraction_D_->Update(1., *poiseuille_force, 1.)) FOUR_C_THROW("Update failed");

  // f_master = +M^T*t // attention: no minus sign here: poiseuille points in same direction on
  // slave and master side
  if (mortaradapter_->GetMortarMatrixM()->Multiply(true, *poiseuille_force, *master_psl))
    FOUR_C_THROW("Multiply failed");
  if (stritraction_M_->Update(1., *poiseuille_force, 1.)) FOUR_C_THROW("update failed");

  // add the contribution
  if (slaveiforce->Update(1., *slave_psl, 1.)) FOUR_C_THROW("Update failed");
  if (masteriforce->Update(1., *master_psl, 1.)) FOUR_C_THROW("Update failed");
}


void EHL::Base::add_couette_force(
    Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce)
{
  const int ndim = Global::Problem::Instance()->NDim();
  const Teuchos::RCP<const Epetra_Vector> relVel = mortaradapter_->RelTangVel();
  Teuchos::RCP<Epetra_Vector> height =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *height))
    FOUR_C_THROW("multiply failed");
  Teuchos::RCP<Epetra_Vector> h_inv =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (h_inv->Reciprocal(*height)) FOUR_C_THROW("Reciprocal failed");
  Teuchos::RCP<Epetra_Vector> hinv_relV =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  hinv_relV->Multiply(1., *h_inv, *relVel, 0.);

  Core::FE::Discretization& lub_dis = *lubrication_->LubricationField()->discretization();
  Teuchos::RCP<Epetra_Vector> visc_vec =
      Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->dof_row_map(1)));
  for (int i = 0; i < lub_dis.NodeRowMap()->NumMyElements(); ++i)
  {
    Core::Nodes::Node* lnode = lub_dis.lRowNode(i);
    if (!lnode) FOUR_C_THROW("node not found");
    const double p = lubrication_->LubricationField()->Prenp()->operator[](
        lubrication_->LubricationField()->Prenp()->Map().LID(lub_dis.Dof(0, lnode, 0)));

    Teuchos::RCP<Core::Mat::Material> mat = lnode->Elements()[0]->Material(0);
    if (mat.is_null()) FOUR_C_THROW("null pointer");
    Teuchos::RCP<Mat::LubricationMat> lmat =
        Teuchos::rcp_dynamic_cast<Mat::LubricationMat>(mat, true);
    const double visc = lmat->ComputeViscosity(p);

    for (int d = 0; d < ndim; ++d) visc_vec->ReplaceGlobalValue(lub_dis.Dof(1, lnode, d), 0, visc);
  }
  Teuchos::RCP<Epetra_Vector> visc_vec_str = ada_strDisp_to_lubDisp_->SlaveToMaster(visc_vec);
  Teuchos::RCP<Epetra_Vector> couette_force =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  couette_force->Multiply(-1., *visc_vec_str, *hinv_relV, 0.);

  Teuchos::RCP<Epetra_Vector> slave_cou =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMortarMatrixD()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> master_cou =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMortarMatrixM()->DomainMap()));
  // f_slave = D^T*t
  if (mortaradapter_->GetMortarMatrixD()->Multiply(true, *couette_force, *slave_cou))
    FOUR_C_THROW("Multiply failed");
  if (stritraction_D_->Update(1., *couette_force, 1.)) FOUR_C_THROW("Update failed");

  // f_master = -M^T*t
  if (mortaradapter_->GetMortarMatrixM()->Multiply(true, *couette_force, *master_cou))
    FOUR_C_THROW("Multiply failed");
  if (stritraction_M_->Update(-1., *couette_force, 1.)) FOUR_C_THROW("update failed");

  // add the contribution
  if (slaveiforce->Update(1., *slave_cou, 1.)) FOUR_C_THROW("Update failed");
  if (masteriforce->Update(-1., *master_cou, 1.)) FOUR_C_THROW("Update failed");
}

/*----------------------------------------------------------------------*
 | set structure velocity fields on lubrication field       seitz 12/17 |
 *----------------------------------------------------------------------*/
void EHL::Base::set_average_velocity_field()
{
  Teuchos::RCP<Epetra_Vector> avVelLub =
      ada_strDisp_to_lubDisp_->MasterToSlave(mortaradapter_->AvTangVel());
  lubrication_->LubricationField()->set_average_velocity_field(1, avVelLub);
}

/*----------------------------------------------------------------------*
 | set structure relative velocity fields on lub. field     faraji 02/19 |
 *----------------------------------------------------------------------*/
void EHL::Base::set_relative_velocity_field()
{
  Teuchos::RCP<Epetra_Vector> relVelLub =
      ada_strDisp_to_lubDisp_->MasterToSlave(mortaradapter_->RelTangVel());
  lubrication_->LubricationField()->set_relative_velocity_field(1, relVelLub);
}

/*----------------------------------------------------------------------*
 | set film height on lubrication field                      wirtz 01/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::set_height_field()
{
  //  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortardinv = mortaradapter_->GetDinvMatrix();
  Teuchos::RCP<Epetra_Vector> discretegap =
      Core::LinAlg::CreateVector(*(slaverowmapextr_->Map(0)), true);

  // get the weighted gap and store it in slave dof map (for each node, the scalar value is stored
  // in the 0th dof)
  int err = slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *discretegap);
  if (err != 0) FOUR_C_THROW("error while transforming map of weighted gap");

  // store discrete gap in lubrication disp dof map (its the film height)
  Teuchos::RCP<Epetra_Vector> height = ada_strDisp_to_lubDisp_->MasterToSlave(discretegap);

  // provide film height to lubrication discretization
  lubrication_->LubricationField()->set_height_field(1, height);
}

/*----------------------------------------------------------------------*
 | set time derivative of film height on lubrication field   Faraji 03/18|
 *----------------------------------------------------------------------*/
void EHL::Base::set_height_dot()
{
  Teuchos::RCP<Epetra_Vector> heightdot =
      Teuchos::rcp(new Epetra_Vector(*(mortaradapter_->Nodal_Gap())));
  Teuchos::RCP<const Epetra_Vector> heightnp = mortaradapter_->Nodal_Gap();

  heightdot->Update(-1.0 / Dt(), *heightold_, 1.0 / Dt());

  Teuchos::RCP<Epetra_Vector> discretegap =
      Core::LinAlg::CreateVector(*(slaverowmapextr_->Map(0)), true);
  // get the weighted heightdot and store it in slave dof map (for each node, the scalar value is
  // stored in the 0th dof)
  int err = slavemaptransform_->Multiply(false, *heightdot, *discretegap);
  if (err != 0) FOUR_C_THROW("error while transforming map of weighted gap");
  // store discrete heightDot in lubrication disp dof map (its the film height time derivative)
  Teuchos::RCP<Epetra_Vector> heightdotSet = ada_strDisp_to_lubDisp_->MasterToSlave(discretegap);

  // provide film height time derivative to lubrication discretization
  lubrication_->LubricationField()->SetHeightDotField(1, heightdotSet);
}

/*----------------------------------------------------------------------*
 | set structure mesh displacement on lubrication field     wirtz 03/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::set_mesh_disp(Teuchos::RCP<const Epetra_Vector> disp)
{
  // Extract the structure displacement at the slave-side interface
  Teuchos::RCP<Epetra_Vector> slaveidisp = Core::LinAlg::CreateVector(
      *(slaverowmapextr_->Map(0)), true);  // Structure displacement at the lubricated interface
  slaverowmapextr_->ExtractVector(disp, 0, slaveidisp);

  // Transfer the displacement vector onto the lubrication field
  Teuchos::RCP<Epetra_Vector> lubridisp = ada_strDisp_to_lubDisp_->MasterToSlave(slaveidisp);

  // Provide the lubrication discretization with the displacement
  lubrication_->LubricationField()->ApplyMeshMovement(lubridisp, 1);
}


/*----------------------------------------------------------------------*
 | Create DBC toggle for unprojectable nodes                seitz 01/18 |
 *----------------------------------------------------------------------*/
void EHL::Base::setup_unprojectable_dbc()
{
  if (not Core::UTILS::IntegralValue<int>(
          ((Global::Problem::Instance()->elasto_hydro_dynamic_params())), "UNPROJ_ZERO_DBC"))
    return;

  Teuchos::RCP<Epetra_FEVector> inf_gap_toggle =
      Teuchos::rcp(new Epetra_FEVector(*mortaradapter_->SlaveDofMap(), true));
  for (int i = 0; i < mortaradapter_->Interface()->SlaveRowNodes()->NumMyElements(); ++i)
  {
    Core::Nodes::Node* node = mortaradapter_->Interface()->Discret().gNode(
        mortaradapter_->Interface()->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("gnode returned nullptr");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("dynamic cast failed");
    if (cnode->Data().Getg() > 1.e11)
    {
      for (int e = 0; e < cnode->NumElement(); ++e)
      {
        Core::Elements::Element* ele = cnode->Elements()[e];
        for (int nn = 0; nn < ele->num_node(); ++nn)
        {
          CONTACT::Node* cnn = dynamic_cast<CONTACT::Node*>(ele->Nodes()[nn]);
          if (!cnn) FOUR_C_THROW("cast failed");
          for (int j = 0; j < 3; ++j)
          {
            const int row = cnn->Dofs()[j];
            const double one = 1.;
            inf_gap_toggle->SumIntoGlobalValues(1, &row, &one, 0);
          }
        }
      }
    }
  }
  if (inf_gap_toggle->GlobalAssemble(Epetra_Max, false) != 0)
    FOUR_C_THROW("global_assemble failed");
  for (int i = 0; i < inf_gap_toggle->Map().NumMyElements(); ++i)
    if (inf_gap_toggle->operator()(0)->operator[](i) > 0.5)
      inf_gap_toggle->operator()(0)->operator[](i) = 1.;

  Teuchos::RCP<Epetra_Vector> exp =
      Teuchos::rcp(new Epetra_Vector(*ada_strDisp_to_lubPres_->MasterDofMap()));
  Core::LinAlg::Export(*inf_gap_toggle, *exp);
  inf_gap_toggle_lub_ = ada_strDisp_to_lubPres_->MasterToSlave(exp);

  static Teuchos::RCP<Epetra_Vector> old_toggle = Teuchos::null;
  if (old_toggle != Teuchos::null)
  {
    for (int i = 0; i < inf_gap_toggle_lub_->Map().NumMyElements(); ++i)
      if (abs(inf_gap_toggle_lub_->operator[](i) - old_toggle->operator[](i)) > 1.e-12)
      {
        if (!Comm().MyPID())
          std::cout << "dbc of unprojectable nodes changed boundary condition" << std::endl;
        break;
      }
  }
  else
  {
    double d = 0.;
    inf_gap_toggle_lub_->MaxValue(&d);

    if (!Comm().MyPID())
      std::cout << "dbc of unprojectable nodes changed boundary condition" << std::endl;
  }
  old_toggle = Teuchos::rcp(new Epetra_Vector(*inf_gap_toggle_lub_));

  lubrication_->LubricationField()->InfGapToggle() =
      Teuchos::rcp(new Epetra_Vector(*inf_gap_toggle_lub_));
}

/*----------------------------------------------------------------------*
 | setup adapters for EHL on boundary                       wirtz 01/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::setup_field_coupling(
    const std::string struct_disname, const std::string lubrication_disname)
{
  Global::Problem* problem = Global::Problem::Instance();
  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<Core::FE::Discretization> lubricationdis = problem->GetDis(lubrication_disname);

  if (structdis.is_null()) FOUR_C_THROW("structure dis does not exist");
  if (lubricationdis.is_null()) FOUR_C_THROW("lubrication dis does not exist");

  const int ndim = Global::Problem::Instance()->NDim();

  //------------------------------------------------------------------
  // 1. Mortar coupling: Slave-side structure <-> Master-side structure
  //------------------------------------------------------------------

  // A mortar coupling adapter, using the "EHL Coupling Condition" is set up. The Coupling is
  // between the slave- and the master-side interface of the structure. Dofs, which, need to be
  // transfered from the master-side to the lubrication field, need to be mortar-projected to the
  // slave-side interface and then transfered by a matching-node coupling,  and vice versa. The
  // matching node coupling is defined below.

  std::vector<int> coupleddof(ndim, 1);
  mortaradapter_ = Teuchos::rcp(new Adapter::CouplingEhlMortar(Global::Problem::Instance()->NDim(),
      Global::Problem::Instance()->mortar_coupling_params(),
      Global::Problem::Instance()->contact_dynamic_params(),
      Global::Problem::Instance()->spatial_approximation_type()));
  mortaradapter_->Setup(structdis, structdis, coupleddof, "EHLCoupling");

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(
          mortaradapter_->Interface()->interface_params(), "STRATEGY") !=
      Inpar::CONTACT::solution_ehl)
    FOUR_C_THROW("you need to set ---CONTACT DYNAMIC: STRATEGY   Ehl");

  Teuchos::RCP<Epetra_Vector> idisp = Core::LinAlg::CreateVector(
      *(structdis->dof_row_map()), true);  // Structure displacement at the lubricated interface
  mortaradapter_->Interface()->Initialize();
  mortaradapter_->Interface()->set_state(Mortar::state_old_displacement, *idisp);
  mortaradapter_->Interface()->set_state(Mortar::state_new_displacement, *idisp);
  mortaradapter_->Interface()->evaluate_nodal_normals();
  mortaradapter_->Interface()->export_nodal_normals();
  mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::n_old);
  mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::dm);
  mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::activeold);
  mortaradapter_->Integrate(idisp, Dt());

  // Maps of the interface dofs
  Teuchos::RCP<Epetra_Map> masterdofrowmap = mortaradapter_->Interface()->MasterRowDofs();
  Teuchos::RCP<Epetra_Map> slavedofrowmap = mortaradapter_->Interface()->SlaveRowDofs();
  Teuchos::RCP<Epetra_Map> mergeddofrowmap =
      Core::LinAlg::MergeMap(masterdofrowmap, slavedofrowmap, false);

  // Map extractors with the structure dofs as full maps and local interface maps
  slaverowmapextr_ = Teuchos::rcp(
      new Core::LinAlg::MapExtractor(*(structdis->dof_row_map()), slavedofrowmap, false));
  masterrowmapextr_ = Teuchos::rcp(
      new Core::LinAlg::MapExtractor(*(structdis->dof_row_map()), masterdofrowmap, false));
  mergedrowmapextr_ = Teuchos::rcp(
      new Core::LinAlg::MapExtractor(*(structdis->dof_row_map()), mergeddofrowmap, false));


  //----------------------------------------------------------
  // 2. build coupling adapters
  //----------------------------------------------------------
  Teuchos::RCP<Epetra_Map> strucnodes = mortaradapter_->Interface()->SlaveRowNodes();
  const Epetra_Map* lubrinodes = lubricationdis->NodeRowMap();
  ada_strDisp_to_lubDisp_ = Teuchos::rcp(new Core::Adapter::Coupling);
  ada_strDisp_to_lubDisp_->setup_coupling(
      *structdis, *lubricationdis, *strucnodes, *lubrinodes, ndim, true, 1.e-8, 0, 1);

  ada_lubPres_to_lubDisp_ = Teuchos::rcp(new Core::Adapter::Coupling);
  ada_lubPres_to_lubDisp_->setup_coupling(*lubricationdis, *lubricationdis,
      *lubricationdis->NodeRowMap(), *lubricationdis->NodeRowMap(), 1, true, 1.e-8, 0, 1);

  ada_strDisp_to_lubPres_ = Teuchos::rcp(new Core::Adapter::Coupling);
  ada_strDisp_to_lubPres_->setup_coupling(mortaradapter_->Interface()->Discret(), *lubricationdis,
      *mortaradapter_->Interface()->SlaveRowNodes(), *lubricationdis->NodeRowMap(), 1, true, 1.e-3);

  // Setup of transformation matrix: slave node map <-> slave disp dof map
  slavemaptransform_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*slavedofrowmap, 81, false, false));
  for (int i = 0; i < mortaradapter_->Interface()->SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = mortaradapter_->Interface()->SlaveRowNodes()->GID(i);
    Core::Nodes::Node* node = structdis->gNode(gid);
    std::vector<int> dofs = structdis->Dof(0, node);
    // slavemaptransform_->Assemble(1.0,dofs[0],gid);
    for (unsigned int idim = 0; idim < dofs.size(); idim++)
    {
      int row = dofs[idim];
      slavemaptransform_->Assemble(1.0, row, gid);
    }
  }
  slavemaptransform_->Complete(*(mortaradapter_->Interface()->SlaveRowNodes()), *slavedofrowmap);

  // Setup of transformation matrix: lubrication pre dof map <-> lubrication disp dof map
  lubrimaptransform_ = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(*(lubricationdis->dof_row_map(1)), 81, false, false));
  for (int inode = 0; inode < lubricationdis->NumMyRowNodes(); ++inode)
  {
    Core::Nodes::Node* node = lubricationdis->lRowNode(inode);
    std::vector<int> nodepredof = lubricationdis->Dof(0, node);
    std::vector<int> nodedispdofs = lubricationdis->Dof(1, node);
    for (unsigned int idim = 0; idim < nodedispdofs.size(); idim++)
      lubrimaptransform_->Assemble(1.0, nodedispdofs[idim], nodepredof[0]);
  }
  lubrimaptransform_->Complete(
      *(lubricationdis->dof_row_map(0)), *(lubricationdis->dof_row_map(1)));
}


/*----------------------------------------------------------------------*
 | update (protected)                                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Base::update()
{
  heightold_ = mortaradapter_->Nodal_Gap();
  mortaradapter_->Interface()->set_state(
      Mortar::state_old_displacement, *structure_field()->Dispnp());
  mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::n_old);
  mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::dm);
  mortaradapter_->Interface()->StoreToOld(Mortar::StrategyBase::activeold);
  structure_field()->Update();
  lubrication_->LubricationField()->Update();

  return;
}

/*----------------------------------------------------------------------*
 | output (protected)                                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Base::output(bool forced_writerestart)
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.

  //===========================
  // output for structurefield:
  //===========================
  //  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());
  structure_field()->Output(forced_writerestart);

  // Additional output on structure field
  structure_field()->disc_writer()->WriteVector("fluid_force",
      EvaluateFluidForce(lubrication_->LubricationField()->Prenp()), Core::IO::dofvector);

  if (dry_contact_)
  {
    Teuchos::RCP<Epetra_Vector> active_toggle, slip_toggle;
    mortaradapter_->create_active_slip_toggle(&active_toggle, &slip_toggle);
    for (int i = 0; i < active_toggle->Map().NumMyElements(); ++i)
      slip_toggle->operator[](i) += active_toggle->operator[](i);
    Teuchos::RCP<Epetra_Vector> active =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->NodeRowMap()));
    Teuchos::RCP<Epetra_Vector> slip =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->NodeRowMap()));
    Core::LinAlg::Export(*active_toggle, *active);
    Core::LinAlg::Export(*slip_toggle, *slip);
    structure_field()->disc_writer()->WriteVector("active", active, Core::IO::dofvector);
    structure_field()->disc_writer()->WriteVector("slip", slip, Core::IO::dofvector);
  }
  if (dry_contact_)
  {
    Teuchos::RCP<Epetra_Vector> n, t;
    mortaradapter_->CreateForceVec(n, t);
    Teuchos::RCP<Epetra_Vector> ne =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->dof_row_map()));
    Teuchos::RCP<Epetra_Vector> te =
        Teuchos::rcp(new Epetra_Vector(*structure_field()->discretization()->dof_row_map()));
    Core::LinAlg::Export(*n, *ne);
    Core::LinAlg::Export(*t, *te);
    structure_field()->disc_writer()->WriteVector("normal_contact", ne, Core::IO::dofvector);
    structure_field()->disc_writer()->WriteVector("tangential_contact", te, Core::IO::dofvector);
  }

  //=============================
  // output for lubricationfield:
  //=============================
  set_mesh_disp(structure_field()->Dispnp());
  lubrication_->LubricationField()->Output(forced_writerestart);

  // ============================
  // output for mortar interface
  // ============================
  mortaradapter_->write_restart(*lubrication_->LubricationField()->DiscWriter());

  // Addtitional output on the lubrication field
  {
    Teuchos::RCP<Epetra_Vector> discretegap =
        Core::LinAlg::CreateVector(*(slaverowmapextr_->Map(0)), true);

    // get the weighted gap and store it in slave dof map (for each node, the scalar value is stored
    // in the 0th dof)
    int err = slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *discretegap);
    if (err != 0) FOUR_C_THROW("error while transforming map of weighted gap");

    // store discrete gap in lubrication disp dof map (its the film height)
    Teuchos::RCP<Epetra_Vector> height = ada_strDisp_to_lubDisp_->MasterToSlave(discretegap);

    Teuchos::RCP<Epetra_Vector> height_ex =
        Teuchos::rcp(new Epetra_Vector(*ada_lubPres_to_lubDisp_->SlaveDofMap()));
    Core::LinAlg::Export(*height, *height_ex);
    Teuchos::RCP<Epetra_Vector> h1 = ada_lubPres_to_lubDisp_->SlaveToMaster(height_ex);
    lubrication_->LubricationField()->DiscWriter()->WriteVector("height", h1, Core::IO::dofvector);

    if (inf_gap_toggle_lub_ != Teuchos::null)
      lubrication_->LubricationField()->DiscWriter()->WriteVector(
          "no_gap_DBC", inf_gap_toggle_lub_, Core::IO::dofvector);

    // output for viscosity

    const int ndim = Global::Problem::Instance()->NDim();
    Teuchos::RCP<Epetra_Vector> visc_vec =
        Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->dof_row_map(1)));
    for (int i = 0;
         i < lubrication_->LubricationField()->discretization()->NodeRowMap()->NumMyElements(); ++i)
    {
      Core::Nodes::Node* lnode = lubrication_->LubricationField()->discretization()->lRowNode(i);
      if (!lnode) FOUR_C_THROW("node not found");
      const double p = lubrication_->LubricationField()->Prenp()->operator[](
          lubrication_->LubricationField()->Prenp()->Map().LID(
              lubrication_->LubricationField()->discretization()->Dof(0, lnode, 0)));
      Teuchos::RCP<Core::Mat::Material> mat = lnode->Elements()[0]->Material(0);
      if (mat.is_null()) FOUR_C_THROW("null pointer");
      Teuchos::RCP<Mat::LubricationMat> lmat =
          Teuchos::rcp_dynamic_cast<Mat::LubricationMat>(mat, true);
      const double visc = lmat->ComputeViscosity(p);

      for (int d = 0; d < ndim; ++d)
        visc_vec->ReplaceGlobalValue(
            lubrication_->LubricationField()->discretization()->Dof(1, lnode, d), 0, visc);
    }

    Teuchos::RCP<Epetra_Vector> visc_vec_ex =
        Teuchos::rcp(new Epetra_Vector(*ada_lubPres_to_lubDisp_->SlaveDofMap()));

    Core::LinAlg::Export(*visc_vec, *visc_vec_ex);

    Teuchos::RCP<Epetra_Vector> v1 = ada_lubPres_to_lubDisp_->SlaveToMaster(visc_vec_ex);

    lubrication_->LubricationField()->DiscWriter()->WriteVector(
        "viscosity", v1, Core::IO::dofvector);
  }

  // reset states
  structure_field()->discretization()->ClearState(true);
  lubrication_->LubricationField()->discretization()->ClearState(true);
}  // Output()

FOUR_C_NAMESPACE_CLOSE
