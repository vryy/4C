/*--------------------------------------------------------------------------*/
/*!
\file ehl_base.cpp

\brief base class for all elastohydrodynamic lubrication (lubrication structure interaction)
algorithms

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/

#include "ehl_base.H"

#include "ehl_partitioned.H"
#include "ehl_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_lubrication.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_ehl_mortar.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lubrication/lubrication_timint_implicit.H"
#include "../drt_mat/lubrication_mat.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_matchingoctree.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"

#include "Epetra_MultiVector.h"

/*----------------------------------------------------------------------*
 | constructor                                     (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
EHL::Base::Base(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& lubricationparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string lubrication_disname)
    : AlgorithmBase(comm, globaltimeparams),
      structure_(Teuchos::null),
      lubrication_(Teuchos::null),
      fieldcoupling_(DRT::INPUT::IntegralValue<INPAR::EHL::FieldCoupling>(
          DRT::Problem::Instance()->ElastoHydroDynamicParams(), "FIELDCOUPLING")),
      dry_contact_(DRT::INPUT::IntegralValue<bool>(
          DRT::Problem::Instance()->ElastoHydroDynamicParams(), "DRY_CONTACT_MODEL"))
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // get the solver number used for Lubrication solver
  const int linsolvernumber = lubricationparams.get<int>("LINEAR_SOLVER");

  // 2.- Setup discretizations and coupling.
  SetupDiscretizations(comm, struct_disname, lubrication_disname);

  SetupFieldCoupling(struct_disname, lubrication_disname);

  // 3.- Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis(struct_disname);

  // set moving grid
  bool isale = true;

  // determine which time params to use to build the single fields
  // in case of time stepping time params have to be read from single field sections
  // in case of equal timestep size for all fields the time params are controlled solely
  // by the problem section (e.g. ehl or cell dynamic)
  const Teuchos::ParameterList* structtimeparams = &globaltimeparams;
  const Teuchos::ParameterList* lubricationtimeparams = &globaltimeparams;
  if (DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->ElastoHydroDynamicParams(), "DIFFTIMESTEPSIZE"))
  {
    structtimeparams = &structparams;
    lubricationtimeparams = &lubricationparams;
  }

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
          *structtimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::Structure>(structure->StructureField());
  structure_->Setup();
  lubrication_ = Teuchos::rcp(new ADAPTER::LubricationBaseAlgorithm());
  lubrication_->Setup(*lubricationtimeparams, lubricationparams,
      problem->SolverParams(linsolvernumber), lubrication_disname, isale);
  mortaradapter_->StoreDirichletStatus(StructureField()->GetDBCMapExtractor());

  // Structure displacement at the lubricated interface
  Teuchos::RCP<Epetra_Vector> disp = LINALG::CreateVector(*(structdis->DofRowMap()), true);

  mortaradapter_->Integrate(disp, Dt());
  // the film thickness initialization for very first time step
  heightold_ = mortaradapter_->Nodal_Gap();
}

/*----------------------------------------------------------------------*
 | read restart information for given time step   (public) wirtz 12/15  |
 *----------------------------------------------------------------------*/
void EHL::Base::ReadRestart(int restart)
{
  if (restart)
  {
    lubrication_->LubricationField()->ReadRestart(restart);
    structure_->ReadRestart(restart);
    SetTimeStep(structure_->TimeOld(), restart);

    mortaradapter_->Interface()->SetState(MORTAR::state_old_displacement, *structure_->Dispn());
    mortaradapter_->Interface()->SetState(MORTAR::state_new_displacement, *structure_->Dispn());
    mortaradapter_->Interface()->EvaluateNodalNormals();
    mortaradapter_->Interface()->ExportNodalNormals();
    mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::n_old);
    mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::dm);
    mortaradapter_->Integrate(structure_->Dispnp(), Dt());
    heightold_ = mortaradapter_->Nodal_Gap();

    IO::DiscretizationReader reader(lubrication_->LubricationField()->Discretization(), restart);
    mortaradapter_->ReadRestart(reader);
  }

  return;
}

/*----------------------------------------------------------------------*
 | calculate velocities by a FD approximation               wirtz 12/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> EHL::Base::CalcVelocity(Teuchos::RCP<const Epetra_Vector> dispnp)
{
  Teuchos::RCP<Epetra_Vector> vel = Teuchos::null;
  // copy D_n onto V_n+1
  vel = Teuchos::rcp(new Epetra_Vector(*(structure_->Dispn())));
  // calculate velocity with timestep Dt()
  //  V_n+1^k = (D_n+1^k - D_n) / Dt
  vel->Update(1. / Dt(), *dispnp, -1. / Dt());

  return vel;
}  // CalcVelocity()

/*----------------------------------------------------------------------*
 | read restart information for given time        (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::ReadRestartfromTime(double restarttime)
{
  if (restarttime > 0.0)
  {
    int restartstructure = EHL::Utils::CheckTimeStepping(structure_->Dt(), restarttime);
    int restartlubrication =
        EHL::Utils::CheckTimeStepping(lubrication_->LubricationField()->Dt(), restarttime);

    lubrication_->LubricationField()->ReadRestart(restartlubrication);
    structure_->ReadRestart(restartstructure);
    SetTimeStep(structure_->TimeOld(), restartstructure);
  }

  return;
}

/*----------------------------------------------------------------------*
 | test results (if necessary)                     (public) wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(lubrication_->CreateLubricationFieldTest());
  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*
 | setup discretizations and dofsets                        wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname,
    const std::string lubrication_disname)
{
  // Scheme   : the structure discretization is received from the input. Then, an ale-lubrication
  // disc. is cloned.

  DRT::Problem* problem = DRT::Problem::Instance();

  // 1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> lubricationdis = problem->GetDis(lubrication_disname);
  if (!structdis->Filled()) structdis->FillComplete();
  if (!lubricationdis->Filled()) lubricationdis->FillComplete();

  // first call FillComplete for single discretizations.
  // This way the physical dofs are numbered successively
  structdis->FillComplete();
  lubricationdis->FillComplete();

  // build auxiliary dofsets, i.e. pseudo dofs on each discretization
  const int ndofpernode_lubrication = lubricationdis->NumDof(0, lubricationdis->lRowNode(0));
  const int ndofperelement_lubrication = 0;
  const int ndofpernode_struct = structdis->NumDof(0, structdis->lRowNode(0));
  const int ndofperelement_struct = 0;

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux_lubrication =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(
          ndofpernode_lubrication, ndofperelement_lubrication, 0, true));
  if (structdis->AddDofSet(dofsetaux_lubrication) != 1)
    dserror("unexpected dof sets in structure field");

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux_struct = Teuchos::rcp(
      new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct, ndofperelement_struct, 0, true));
  if (lubricationdis->AddDofSet(dofsetaux_struct) != 1)
    dserror("unexpected dof sets in lubrication field");

  // call AssignDegreesOfFreedom also for auxiliary dofsets
  // note: the order of FillComplete() calls determines the gid numbering!
  // 1. structure dofs
  // 2. lubrication dofs
  // 3. structure auxiliary dofs
  // 4. lubrication auxiliary dofs
  structdis->FillComplete(true, false, false);
  lubricationdis->FillComplete(true, false, false);
}

/*----------------------------------------------------------------------*
 | set structure solution on lubrication field              wirtz 12/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetStructSolution(Teuchos::RCP<const Epetra_Vector> disp)
{
  //---------------------------------------------------------
  // 1. Update the Mortar Coupling
  //---------------------------------------------------------

  //  //Extract the structure displacement at the lubricated interface
  //  Teuchos::RCP<Epetra_Vector> idisp = LINALG::CreateVector(*(mergedrowmapextr_->Map(0)),
  //  true);//Structure displacement at the lubricated interface
  //  mergedrowmapextr_->ExtractVector(disp,0,idisp);
  // Reevalute the mortar martices D and M
  mortaradapter_->Integrate(disp, Dt());

  // Displace the mesh of the lubrication field in accordance with the slave-side interface
  SetMeshDisp(disp);

  // Calculate the average tangential fractions of the structure velocities at the interface and
  // provide them to the lubrication field
  SetAverageVelocityField();

  // Calculate the relative tangential fractions of the structure velocities at the interface and
  // provide them to the lubrication field
  SetRelativeVelocityField();

  // Provide the gap at the interface
  SetHeightField();

  // provide the heightdot (time derivative of the gap)
  SetHeightDot();

  // Create DBC map for unprojectable nodes
  SetupUnprojectableDBC();

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
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetDMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> masteriforce =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMMatrix()->DomainMap()));

  stritraction_D_ = Teuchos::rcp(new Epetra_Vector(*ada_strDisp_to_lubDisp_->MasterDofMap()));
  stritraction_M_ = Teuchos::rcp(new Epetra_Vector(*ada_strDisp_to_lubDisp_->MasterDofMap()));

  // add pressure force
  AddPressureForce(slaveiforce, masteriforce);
  // add poiseuille flow force
  AddPoiseuilleForce(slaveiforce, masteriforce);
  // add couette flow force
  AddCouetteForce(slaveiforce, masteriforce);

  // External force vector (global)
  Teuchos::RCP<Epetra_Vector> strforce =
      Teuchos::rcp(new Epetra_Vector(*(structure_->DofRowMap())));

  // Insert both interface forces into the global force vector
  slaverowmapextr_->InsertVector(slaveiforce, 0, strforce);
  masterrowmapextr_->InsertVector(masteriforce, 0, strforce);

  return strforce;
}

/*----------------------------------------------------------------------*
 | set tractions, resulting from lubrication pressure       seitz 01/18 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetLubricationSolution(Teuchos::RCP<const Epetra_Vector> pressure)
{
  // Provide the structure field with the force vector
  // Note that the mid-point values (gen-alpha) of the interface forces are evaluated in
  // STR::TimIntGenAlpha::EvaluateForceResidual()
  structure_->SetForceInterface(EvaluateFluidForce(pressure));
}

void EHL::Base::AddPressureForce(
    Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce)
{
  Teuchos::RCP<Epetra_Vector> stritraction;

  Teuchos::RCP<Epetra_Vector> p_full =
      Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->DofRowMap(1)));
  if (lubrimaptransform_->Apply(*lubrication_->LubricationField()->Prenp(), *p_full))
    dserror("apply failed");
  Teuchos::RCP<Epetra_Vector> p_exp =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  p_exp = ada_strDisp_to_lubDisp_->SlaveToMaster(p_full);
  stritraction = Teuchos::rcp(new Epetra_Vector(*mortaradapter_->Normals()));
  stritraction->Multiply(-1., *mortaradapter_->Normals(), *p_exp, 0.);

  // Get the Mortar D and M Matrix
  const Teuchos::RCP<LINALG::SparseMatrix> mortard = mortaradapter_->GetDMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix> mortarm = mortaradapter_->GetMMatrix();

  // f_slave = D^T*t
  int err = mortard->Multiply(true, *stritraction, *slaveiforce);
  if (err != 0) dserror("error while calculating slave side interface force");
  if (stritraction_D_->Update(1., *stritraction, 1.)) dserror("Update failed");

  // f_master = -M^T*t
  err = mortarm->Multiply(true, *stritraction, *masteriforce);
  if (err != 0) dserror("error while calculating master side interface force");
  masteriforce->Scale(-1.0);
  if (stritraction_M_->Update(-1., *stritraction, 1.)) dserror("update failed");
}

void EHL::Base::AddPoiseuilleForce(
    Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce)
{
  // poiseuille flow forces
  Teuchos::RCP<Epetra_Vector> p_int =
      ada_strDisp_to_lubPres_->SlaveToMaster(lubrication_->LubricationField()->Prenp());
  Teuchos::RCP<Epetra_Vector> p_int_full =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  LINALG::Export(*p_int, *p_int_full);

  Teuchos::RCP<Epetra_Vector> nodal_gap =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *nodal_gap))
    dserror("multiply failed");

  LINALG::SparseMatrix m(*mortaradapter_->SurfGradMatrix());

  m.LeftScale(*nodal_gap);
  m.Scale(-.5);

  Teuchos::RCP<Epetra_Vector> poiseuille_force =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  m.Apply(*p_int_full, *poiseuille_force);

  Teuchos::RCP<Epetra_Vector> slave_psl =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetDMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> master_psl =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMMatrix()->DomainMap()));

  // f_slave = D^T*t
  if (mortaradapter_->GetDMatrix()->Multiply(true, *poiseuille_force, *slave_psl))
    dserror("Multiply failed");
  if (stritraction_D_->Update(1., *poiseuille_force, 1.)) dserror("Update failed");

  // f_master = +M^T*t // attention: no minus sign here: poiseuille points in same direction on
  // slave and master side
  if (mortaradapter_->GetMMatrix()->Multiply(true, *poiseuille_force, *master_psl))
    dserror("Multiply failed");
  if (stritraction_M_->Update(1., *poiseuille_force, 1.)) dserror("update failed");

  // add the contribution
  if (slaveiforce->Update(1., *slave_psl, 1.)) dserror("Update failed");
  if (masteriforce->Update(1., *master_psl, 1.)) dserror("Update failed");
}


void EHL::Base::AddCouetteForce(
    Teuchos::RCP<Epetra_Vector> slaveiforce, Teuchos::RCP<Epetra_Vector> masteriforce)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  const Teuchos::RCP<const Epetra_Vector> relVel = mortaradapter_->RelTangVel();
  Teuchos::RCP<Epetra_Vector> height =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *height))
    dserror("multiply failed");
  Teuchos::RCP<Epetra_Vector> h_inv =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  if (h_inv->Reciprocal(*height)) dserror("Reciprocal failed");
  Teuchos::RCP<Epetra_Vector> hinv_relV =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  hinv_relV->Multiply(1., *h_inv, *relVel, 0.);

  DRT::Discretization& lub_dis = *lubrication_->LubricationField()->Discretization();
  Teuchos::RCP<Epetra_Vector> visc_vec =
      Teuchos::rcp(new Epetra_Vector(*lubrication_->LubricationField()->DofRowMap(1)));
  for (int i = 0; i < lub_dis.NodeRowMap()->NumMyElements(); ++i)
  {
    DRT::Node* lnode = lub_dis.lRowNode(i);
    if (!lnode) dserror("node not found");
    const double p = lubrication_->LubricationField()->Prenp()->operator[](
        lubrication_->LubricationField()->Prenp()->Map().LID(lub_dis.Dof(0, lnode, 0)));

    Teuchos::RCP<MAT::Material> mat = lnode->Elements()[0]->Material(0);
    if (mat.is_null()) dserror("null pointer");
    Teuchos::RCP<MAT::LubricationMat> lmat =
        Teuchos::rcp_dynamic_cast<MAT::LubricationMat>(mat, true);
    const double visc = lmat->ComputeViscosity(p);

    for (int d = 0; d < ndim; ++d) visc_vec->ReplaceGlobalValue(lub_dis.Dof(1, lnode, d), 0, visc);
  }
  Teuchos::RCP<Epetra_Vector> visc_vec_str = ada_strDisp_to_lubDisp_->SlaveToMaster(visc_vec);
  Teuchos::RCP<Epetra_Vector> couette_force =
      Teuchos::rcp(new Epetra_Vector(*mortaradapter_->SlaveDofMap()));
  couette_force->Multiply(-1., *visc_vec_str, *hinv_relV, 0.);

  Teuchos::RCP<Epetra_Vector> slave_cou =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetDMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> master_cou =
      Teuchos::rcp(new Epetra_Vector(mortaradapter_->GetMMatrix()->DomainMap()));
  // f_slave = D^T*t
  if (mortaradapter_->GetDMatrix()->Multiply(true, *couette_force, *slave_cou))
    dserror("Multiply failed");
  if (stritraction_D_->Update(1., *couette_force, 1.)) dserror("Update failed");

  // f_master = -M^T*t
  if (mortaradapter_->GetMMatrix()->Multiply(true, *couette_force, *master_cou))
    dserror("Multiply failed");
  if (stritraction_M_->Update(-1., *couette_force, 1.)) dserror("update failed");

  // add the contribution
  if (slaveiforce->Update(1., *slave_cou, 1.)) dserror("Update failed");
  if (masteriforce->Update(-1., *master_cou, 1.)) dserror("Update failed");
}

/*----------------------------------------------------------------------*
 | set structure velocity fields on lubrication field       seitz 12/17 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetAverageVelocityField()
{
  Teuchos::RCP<Epetra_Vector> avVelLub =
      ada_strDisp_to_lubDisp_->MasterToSlave(mortaradapter_->AvTangVel());
  lubrication_->LubricationField()->SetAverageVelocityField(1, avVelLub);
}

/*----------------------------------------------------------------------*
 | set structure relative velocity fields on lub. field     faraji 02/19 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetRelativeVelocityField()
{
  Teuchos::RCP<Epetra_Vector> relVelLub =
      ada_strDisp_to_lubDisp_->MasterToSlave(mortaradapter_->RelTangVel());
  lubrication_->LubricationField()->SetRelativeVelocityField(1, relVelLub);
}

/*----------------------------------------------------------------------*
 | set film height on lubrication field                      wirtz 01/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetHeightField()
{
  //  const Teuchos::RCP<LINALG::SparseMatrix> mortardinv = mortaradapter_->GetDinvMatrix();
  Teuchos::RCP<Epetra_Vector> discretegap = LINALG::CreateVector(*(slaverowmapextr_->Map(0)), true);

  // get the weighted gap and store it in slave dof map (for each node, the scalar value is stored
  // in the 0th dof)
  int err = slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *discretegap);
  if (err != 0) dserror("error while transforming map of weighted gap");

  // store discrete gap in lubrication disp dof map (its the film height)
  Teuchos::RCP<Epetra_Vector> height = ada_strDisp_to_lubDisp_->MasterToSlave(discretegap);

  // provide film height to lubrication discretization
  lubrication_->LubricationField()->SetHeightField(1, height);
}

/*----------------------------------------------------------------------*
 | set time derivative of film height on lubrication field   Faraji 03/18|
 *----------------------------------------------------------------------*/
void EHL::Base::SetHeightDot()
{
  Teuchos::RCP<Epetra_Vector> heightdot =
      Teuchos::rcp(new Epetra_Vector(*(mortaradapter_->Nodal_Gap())));
  Teuchos::RCP<const Epetra_Vector> heightnp = mortaradapter_->Nodal_Gap();

  heightdot->Update(-1.0 / Dt(), *heightold_, 1.0 / Dt());

  Teuchos::RCP<Epetra_Vector> discretegap = LINALG::CreateVector(*(slaverowmapextr_->Map(0)), true);
  // get the weighted heightdot and store it in slave dof map (for each node, the scalar value is
  // stored in the 0th dof)
  int err = slavemaptransform_->Multiply(false, *heightdot, *discretegap);
  if (err != 0) dserror("error while transforming map of weighted gap");
  // store discrete heightDot in lubrication disp dof map (its the film height time derivative)
  Teuchos::RCP<Epetra_Vector> heightdotSet = ada_strDisp_to_lubDisp_->MasterToSlave(discretegap);

  // provide film height time derivative to lubrication discretization
  lubrication_->LubricationField()->SetHeightDotField(1, heightdotSet);
}

/*----------------------------------------------------------------------*
 | set structure mesh displacement on lubrication field     wirtz 03/15 |
 *----------------------------------------------------------------------*/
void EHL::Base::SetMeshDisp(Teuchos::RCP<const Epetra_Vector> disp)
{
  // Extract the structure displacement at the slave-side interface
  Teuchos::RCP<Epetra_Vector> slaveidisp = LINALG::CreateVector(
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
void EHL::Base::SetupUnprojectableDBC()
{
  if (not DRT::INPUT::IntegralValue<int>(
          ((DRT::Problem::Instance()->ElastoHydroDynamicParams())), "UNPROJ_ZERO_DBC"))
    return;

  Teuchos::RCP<Epetra_FEVector> inf_gap_toggle =
      Teuchos::rcp(new Epetra_FEVector(*mortaradapter_->SlaveDofMap(), true));
  for (int i = 0; i < mortaradapter_->Interface()->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = mortaradapter_->Interface()->Discret().gNode(
        mortaradapter_->Interface()->SlaveRowNodes()->GID(i));
    if (!node) dserror("gnode returned NULL");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("dynamic cast failed");
    if (cnode->CoData().Getg() > 1.e11)
    {
      for (int e = 0; e < cnode->NumElement(); ++e)
      {
        DRT::Element* ele = cnode->Elements()[e];
        for (int nn = 0; nn < ele->NumNode(); ++nn)
        {
          CONTACT::CoNode* cnn = dynamic_cast<CONTACT::CoNode*>(ele->Nodes()[nn]);
          if (!cnn) dserror("cast failed");
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
  if (inf_gap_toggle->GlobalAssemble(Epetra_Max, false) != 0) dserror("global_assemble failed");
  for (int i = 0; i < inf_gap_toggle->Map().NumMyElements(); ++i)
    if (inf_gap_toggle->operator()(0)->operator[](i) > 0.5)
      inf_gap_toggle->operator()(0)->operator[](i) = 1.;

  Teuchos::RCP<Epetra_Vector> exp =
      Teuchos::rcp(new Epetra_Vector(*ada_strDisp_to_lubPres_->MasterDofMap()));
  LINALG::Export(*inf_gap_toggle, *exp);
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
void EHL::Base::SetupFieldCoupling(
    const std::string struct_disname, const std::string lubrication_disname)
{
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);
  Teuchos::RCP<DRT::Discretization> lubricationdis = problem->GetDis(lubrication_disname);

  if (structdis.is_null()) dserror("structure dis does not exist");
  if (lubricationdis.is_null()) dserror("lubrication dis does not exist");

  const int ndim = DRT::Problem::Instance()->NDim();

  //------------------------------------------------------------------
  // 1. Mortar coupling: Slave-side structure <-> Master-side structure
  //------------------------------------------------------------------

  // A mortar coupling adapter, using the "EHL Coupling Condition" is set up. The Coupling is
  // between the slave- and the master-side interface of the structure. Dofs, which, need to be
  // transfered from the master-side to the lubrication field, need to be mortar-projected to the
  // slave-side interface and then transfered by a matching-node coupling,  and vice versa. The
  // matching node coupling is defined below.

  std::vector<int> coupleddof(ndim, 1);
  mortaradapter_ = Teuchos::rcp(new ADAPTER::CouplingEhlMortar);
  mortaradapter_->Setup(structdis, structdis, coupleddof, "EHLCoupling");

  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(
          mortaradapter_->Interface()->IParams(), "STRATEGY") != INPAR::CONTACT::solution_ehl)
    dserror("you need to set ---CONTACT DYNAMIC: STRATEGY   Ehl");

  Teuchos::RCP<Epetra_Vector> idisp = LINALG::CreateVector(
      *(structdis->DofRowMap()), true);  // Structure displacement at the lubricated interface
  mortaradapter_->Interface()->Initialize();
  mortaradapter_->Interface()->SetState(MORTAR::state_old_displacement, *idisp);
  mortaradapter_->Interface()->SetState(MORTAR::state_new_displacement, *idisp);
  mortaradapter_->Interface()->EvaluateNodalNormals();
  mortaradapter_->Interface()->ExportNodalNormals();
  mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::n_old);
  mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::dm);
  mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::activeold);
  mortaradapter_->Integrate(idisp, Dt());

  // Maps of the interface dofs
  Teuchos::RCP<Epetra_Map> masterdofrowmap = mortaradapter_->Interface()->MasterRowDofs();
  Teuchos::RCP<Epetra_Map> slavedofrowmap = mortaradapter_->Interface()->SlaveRowDofs();
  Teuchos::RCP<Epetra_Map> mergeddofrowmap =
      LINALG::MergeMap(masterdofrowmap, slavedofrowmap, false);

  // Map extractors with the structure dofs as full maps and local interface maps
  slaverowmapextr_ =
      Teuchos::rcp(new LINALG::MapExtractor(*(structdis->DofRowMap()), slavedofrowmap, false));
  masterrowmapextr_ =
      Teuchos::rcp(new LINALG::MapExtractor(*(structdis->DofRowMap()), masterdofrowmap, false));
  mergedrowmapextr_ =
      Teuchos::rcp(new LINALG::MapExtractor(*(structdis->DofRowMap()), mergeddofrowmap, false));


  //----------------------------------------------------------
  // 2. build coupling adapters
  //----------------------------------------------------------
  Teuchos::RCP<Epetra_Map> strucnodes = mortaradapter_->Interface()->SlaveRowNodes();
  const Epetra_Map* lubrinodes = lubricationdis->NodeRowMap();
  ada_strDisp_to_lubDisp_ = Teuchos::rcp(new ADAPTER::Coupling);
  ada_strDisp_to_lubDisp_->SetupCoupling(
      *structdis, *lubricationdis, *strucnodes, *lubrinodes, ndim, true, 1.e-8, 0, 1);

  ada_lubPres_to_lubDisp_ = Teuchos::rcp(new ADAPTER::Coupling);
  ada_lubPres_to_lubDisp_->SetupCoupling(*lubricationdis, *lubricationdis,
      *lubricationdis->NodeRowMap(), *lubricationdis->NodeRowMap(), 1, true, 1.e-8, 0, 1);

  ada_strDisp_to_lubPres_ = Teuchos::rcp(new ADAPTER::Coupling);
  ada_strDisp_to_lubPres_->SetupCoupling(mortaradapter_->Interface()->Discret(), *lubricationdis,
      *mortaradapter_->Interface()->SlaveRowNodes(), *lubricationdis->NodeRowMap(), 1, true, 1.e-3);

  // Setup of transformation matrix: slave node map <-> slave disp dof map
  slavemaptransform_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap, 81, false, false));
  for (int i = 0; i < mortaradapter_->Interface()->SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = mortaradapter_->Interface()->SlaveRowNodes()->GID(i);
    DRT::Node* node = structdis->gNode(gid);
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
  lubrimaptransform_ =
      Teuchos::rcp(new LINALG::SparseMatrix(*(lubricationdis->DofRowMap(1)), 81, false, false));
  for (int inode = 0; inode < lubricationdis->NumMyRowNodes(); ++inode)
  {
    DRT::Node* node = lubricationdis->lRowNode(inode);
    std::vector<int> nodepredof = lubricationdis->Dof(0, node);
    std::vector<int> nodedispdofs = lubricationdis->Dof(1, node);
    for (unsigned int idim = 0; idim < nodedispdofs.size(); idim++)
      lubrimaptransform_->Assemble(1.0, nodedispdofs[idim], nodepredof[0]);
  }
  lubrimaptransform_->Complete(*(lubricationdis->DofRowMap(0)), *(lubricationdis->DofRowMap(1)));
}


/*----------------------------------------------------------------------*
 | update (protected)                                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Base::Update()
{
  heightold_ = mortaradapter_->Nodal_Gap();
  mortaradapter_->Interface()->SetState(
      MORTAR::state_old_displacement, *StructureField()->Dispnp());
  mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::n_old);
  mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::dm);
  mortaradapter_->Interface()->StoreToOld(MORTAR::StrategyBase::activeold);
  StructureField()->Update();
  lubrication_->LubricationField()->Update();

  return;
}

/*----------------------------------------------------------------------*
 | output (protected)                                       wirtz 01/16 |
 *----------------------------------------------------------------------*/
void EHL::Base::Output(bool forced_writerestart)
{
  // Note: The order in the output is important here!

  // In here control file entries are written. And these entries define the
  // order in which the filters handle the Discretizations, which in turn
  // defines the dof number ordering of the Discretizations.

  //===========================
  // output for structurefield:
  //===========================
  //  ApplyLubricationCouplingState(lubrication_->LubricationField()->Prenp());
  StructureField()->Output(forced_writerestart);

  // Additional output on structure field
  StructureField()->DiscWriter()->WriteVector(
      "fluid_force", EvaluateFluidForce(lubrication_->LubricationField()->Prenp()), IO::dofvector);

  if (dry_contact_)
  {
    Teuchos::RCP<Epetra_Vector> active_toggle, slip_toggle;
    mortaradapter_->CreateActiveSlipToggle(&active_toggle, &slip_toggle);
    for (int i = 0; i < active_toggle->Map().NumMyElements(); ++i)
      slip_toggle->operator[](i) += active_toggle->operator[](i);
    Teuchos::RCP<Epetra_Vector> active =
        Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->NodeRowMap()));
    Teuchos::RCP<Epetra_Vector> slip =
        Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->NodeRowMap()));
    LINALG::Export(*active_toggle, *active);
    LINALG::Export(*slip_toggle, *slip);
    StructureField()->DiscWriter()->WriteVector("active", active, IO::dofvector);
    StructureField()->DiscWriter()->WriteVector("slip", slip, IO::dofvector);
  }
  if (dry_contact_)
  {
    Teuchos::RCP<Epetra_Vector> n, t;
    mortaradapter_->CreateForceVec(n, t);
    Teuchos::RCP<Epetra_Vector> ne =
        Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap()));
    Teuchos::RCP<Epetra_Vector> te =
        Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap()));
    LINALG::Export(*n, *ne);
    LINALG::Export(*t, *te);
    StructureField()->DiscWriter()->WriteVector("normal_contact", ne, IO::dofvector);
    StructureField()->DiscWriter()->WriteVector("tangential_contact", te, IO::dofvector);
  }

  //=============================
  // output for lubricationfield:
  //=============================
  SetMeshDisp(StructureField()->Dispnp());
  lubrication_->LubricationField()->Output(forced_writerestart);

  // ============================
  // output for mortar interface
  // ============================
  mortaradapter_->WriteRestart(*lubrication_->LubricationField()->DiscWriter());

  // Addtitional output on the lubrication field
  {
    Teuchos::RCP<Epetra_Vector> discretegap =
        LINALG::CreateVector(*(slaverowmapextr_->Map(0)), true);

    // get the weighted gap and store it in slave dof map (for each node, the scalar value is stored
    // in the 0th dof)
    int err = slavemaptransform_->Multiply(false, *mortaradapter_->Nodal_Gap(), *discretegap);
    if (err != 0) dserror("error while transforming map of weighted gap");

    // store discrete gap in lubrication disp dof map (its the film height)
    Teuchos::RCP<Epetra_Vector> height = ada_strDisp_to_lubDisp_->MasterToSlave(discretegap);

    Teuchos::RCP<Epetra_Vector> height_ex =
        Teuchos::rcp(new Epetra_Vector(*ada_lubPres_to_lubDisp_->SlaveDofMap()));
    LINALG::Export(*height, *height_ex);
    Teuchos::RCP<Epetra_Vector> h1 = ada_lubPres_to_lubDisp_->SlaveToMaster(height_ex);
    lubrication_->LubricationField()->DiscWriter()->WriteVector("height", h1, IO::dofvector);

    if (inf_gap_toggle_lub_ != Teuchos::null)
      lubrication_->LubricationField()->DiscWriter()->WriteVector(
          "no_gap_DBC", inf_gap_toggle_lub_, IO::dofvector);
  }

  // reset states
  StructureField()->Discretization()->ClearState(true);
  lubrication_->LubricationField()->Discretization()->ClearState(true);
}  // Output()
