/*----------------------------------------------------------------------*/
/*! \file

 \brief  Basis of all porous media algorithms

 \level 2


 *-----------------------------------------------------------------------*/

#include "poroelast_base.H"

// needed for PrintNewton
#include <cstddef>

#include "poroelast_defines.H"
#include "poroelast_utils.H"

#include "adapter_coupling.H"
#include "adapter_coupling_volmortar.H"
#include "adapter_fld_base_algorithm.H"
#include "adapter_fld_poro.H"
#include "adapter_str_fpsiwrapper.H"

// new structural time integration
#include "adapter_str_structure_new.H"
#include "adapter_str_factory.H"

// contact
#include "contact_poro_lagrange_strategy.H"
#include "contact_meshtying_contact_bridge.H"

#include "io_control.H"

#include "lib_assemblestrategy.H"
#include "lib_condition_utils.H"
#include "lib_dofset_gidbased_wrapper.H"
#include "lib_globalproblem.H"

#include "linalg_utils_sparse_algebra_assemble.H"
#include "linalg_utils_sparse_algebra_create.H"
#include "solver_linalg_solver.H"

// contact
#include "mortar_manager_base.H"

#include "structure_aux.H"


POROELAST::PoroBase::PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : AlgorithmBase(comm, timeparams),
      is_part_of_multifield_problem_(false),
      matchinggrid_(DRT::INPUT::IntegralValue<bool>(
          DRT::Problem::Instance()->PoroelastDynamicParams(), "MATCHINGGRID")),
      oldstructimint_(DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
                          DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") ==
                      INPAR::STR::int_old)
{
  if (DRT::Problem::Instance()->GetProblemType() != ProblemType::poroelast)
    is_part_of_multifield_problem_ = true;

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  if (!matchinggrid_)
  {
    Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("porofluid");
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_ = Teuchos::rcp(new ADAPTER::MortarVolCoupl());

    // build material strategy
    Teuchos::RCP<UTILS::PoroMaterialStrategy> materialstrategy =
        Teuchos::rcp(new UTILS::PoroMaterialStrategy());

    // setup projection matrices
    volcoupl_->Init(structdis, fluiddis, nullptr, nullptr, nullptr, nullptr, materialstrategy);
    volcoupl_->Redistribute();
    volcoupl_->Setup();
  }

  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // create the structural time integrator (Init() called inside)
  // clean up as soon as old time integration is unused!
  if (oldstructimint_)
  {
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
            timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
    structure_ =
        Teuchos::rcp_dynamic_cast<ADAPTER::FPSIStructureWrapper>(structure->StructureField());
    structure_->Setup();
  }
  else
  {
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr =
        ADAPTER::STR::BuildStructureAlgorithm(sdyn);
    adapterbase_ptr->Init(timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
    adapterbase_ptr->Setup();
    structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::FPSIStructureWrapper>(
        adapterbase_ptr->StructureField());
  }

  if (structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FPSIStructureWrapper failed");

  // ask base algorithm for the fluid time integrator
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fluiddynparams = problem->FluidDynamicParams();
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams, fluiddynparams, "porofluid", true));
  fluid_ = Teuchos::rcp_dynamic_cast<ADAPTER::FluidPoro>(fluid->FluidField());

  if (fluid_ == Teuchos::null)
    dserror("cast from ADAPTER::FluidBaseAlgorithm to ADAPTER::FluidPoro failed");

  // as this is a two way coupled problem, every discretization needs to know the other one.
  // For this we use DofSetProxies and coupling objects which are setup here
  SetupCoupling();

  if (submeshes_) ReplaceDofSets();

  // extractor for constraints on structure phase
  //
  // when using constraints applied via Lagrange-Multipliers there is a
  // difference between StructureField()->DofRowMap() and StructureField()->DofRowMap(0).
  // StructureField()->DofRowMap(0) returns the DofRowMap
  // known to the discretization (without lagrange multipliers)
  // while StructureField()->DofRowMap() returns the DofRowMap known to
  // the constraint manager (with lagrange multipliers)
  cond_splitter_ = Teuchos::rcp(
      new LINALG::MapExtractor(*StructureField()->DofRowMap(), StructureField()->DofRowMap(0)));

  // look for special poro conditions and set flags
  CheckForPoroConditions();

  // do some checks
  {
    // access the problem-specific parameter lists
    const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();

    std::vector<DRT::Condition*> porocoupl;
    FluidField()->Discretization()->GetCondition("PoroCoupling", porocoupl);
    if (porocoupl.size() == 0)
      dserror("no Poro Coupling Condition defined for porous media problem. Fix your input file!");

    // check time integration algo -> currently only one-step-theta scheme supported
    auto structtimealgo = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP");
    auto fluidtimealgo =
        DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR");

    if (not((structtimealgo == INPAR::STR::dyna_onesteptheta and
                fluidtimealgo == INPAR::FLUID::timeint_one_step_theta) or
            (structtimealgo == INPAR::STR::dyna_statics and
                fluidtimealgo == INPAR::FLUID::timeint_stationary) or
            (structtimealgo == INPAR::STR::dyna_genalpha and
                (fluidtimealgo == INPAR::FLUID::timeint_afgenalpha or
                    fluidtimealgo == INPAR::FLUID::timeint_npgenalpha))))
    {
      dserror(
          "porous media problem is limited in functionality (only one-step-theta scheme, "
          "stationary and (af)genalpha case possible)");
    }

    if (fluidtimealgo == INPAR::FLUID::timeint_npgenalpha)
    {
      dserror(
          "npgenalpha time integration for porous fluid is possibly not valid. Either check the "
          "theory or use afgenalpha instead!");
    }

    if (structtimealgo == INPAR::STR::dyna_onesteptheta and
        fluidtimealgo == INPAR::FLUID::timeint_one_step_theta)
    {
      double theta_struct = sdyn.sublist("ONESTEPTHETA").get<double>("THETA");
      double theta_fluid = fdyn.get<double>("THETA");

      if (theta_struct != theta_fluid)
      {
        dserror(
            "porous media problem is limited in functionality. Only one-step-theta scheme with "
            "equal theta for both fields possible. Fix your input file.");
      }
    }

    std::string damping = sdyn.get<std::string>("DAMPING");
    if (damping != "Material")
    {
      dserror(
          "Material damping has to be used for porous media! Set DAMPING to 'Material' in the "
          "STRUCTURAL DYNAMIC section.");
    }

    // access the problem-specific parameter lists
    const Teuchos::ParameterList& pedyn = DRT::Problem::Instance()->PoroelastDynamicParams();
    auto physicaltype =
        DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(pedyn, "PHYSICAL_TYPE");
    if (porosity_dof_ and physicaltype != INPAR::FLUID::poro_p1)
    {
      dserror(
          "Poro P1 elements need a special fluid. Set 'PHYSICAL_TYPE' to 'Poro_P1' in the FLUID "
          "DYNAMIC section!");
    }

    auto transientfluid =
        DRT::INPUT::IntegralValue<INPAR::POROELAST::TransientEquationsOfPoroFluid>(
            pedyn, "TRANSIENT_TERMS");

    if (fluidtimealgo == INPAR::FLUID::timeint_stationary)
    {
      if (transientfluid != INPAR::POROELAST::transient_none)
      {
        dserror(
            "Invalid option for stationary fluid! Set 'TRANSIENT_TERMS' in section POROELASTICITY "
            "DYNAMIC to 'none'!");
      }
    }
    else
    {
      if (transientfluid == INPAR::POROELAST::transient_none)
      {
        dserror(
            "Invalid option for stationary fluid! Set 'TRANSIENT_TERMS' in section POROELASTICITY "
            "DYNAMIC to valid parameter!");
      }
    }

    if (transientfluid == INPAR::POROELAST::transient_momentum_only)
    {
      dserror(
          "Option 'momentum' for parameter 'TRANSIENT_TERMS' in section POROELASTICITY DYNAMIC is "
          "not working properly! There is probably a bug in the linearization ....");
    }
  }
}

void POROELAST::PoroBase::ReadRestart(const int step)
{
  if (step)
  {
    if (not oldstructimint_) structure_->Setup();

    // apply current velocity and pressures to structure
    SetFluidSolution();
    // apply current structural displacements to fluid
    SetStructSolution();

    FluidField()->ReadRestart(step);
    StructureField()->ReadRestart(step);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (submeshes_) ReplaceDofSets();

    // apply current velocity and pressures to structure
    SetFluidSolution();
    // apply current structural displacements to fluid
    SetStructSolution();

    // second ReadRestart needed due to the coupling variables
    FluidField()->ReadRestart(step);
    StructureField()->ReadRestart(step);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (submeshes_) ReplaceDofSets();

    // set the current time in the algorithm (taken from fluid field)
    SetTimeStep(FluidField()->Time(), step);

    // Material pointers to other field were deleted during ReadRestart().
    // They need to be reset.
    if (matchinggrid_)
    {
      POROELAST::UTILS::SetMaterialPointersMatchingGrid(
          StructureField()->Discretization(), FluidField()->Discretization());
    }
    else
    {
      // build material strategy
      Teuchos::RCP<UTILS::PoroMaterialStrategy> materialstrategy =
          Teuchos::rcp(new UTILS::PoroMaterialStrategy());

      volcoupl_->AssignMaterials(
          StructureField()->Discretization(), FluidField()->Discretization(), materialstrategy);
    }
  }
}

void POROELAST::PoroBase::PrepareTimeStep()
{
  // counter and print header
  IncrementTimeAndStep();
  if (!is_part_of_multifield_problem_) PrintHeader();

  // set fluid velocities and pressures onto the structure
  SetFluidSolution();

  // call the predictor
  StructureField()->PrepareTimeStep();

  // set structure displacements onto the fluid
  SetStructSolution();

  // call the predictor
  FluidField()->PrepareTimeStep();
}

void POROELAST::PoroBase::Update()
{
  StructureField()->Update();
  FluidField()->Update();
  // clean up as soon as old time integration is unused!
  if (oldstructimint_)
  {
    if (StructureField()->MeshtyingContactBridge() != Teuchos::null)
    {
      if (StructureField()->MeshtyingContactBridge()->HaveContact() && !nit_contact_)
      {
        (static_cast<CONTACT::PoroLagrangeStrategy&>(
             StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy()))
            .UpdatePoroContact();
      }
    }
  }
}

void POROELAST::PoroBase::PrepareOutput(bool force_prepare_timestep)
{
  StructureField()->PrepareOutput(force_prepare_timestep);
}

void POROELAST::PoroBase::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(StructureField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
}

Teuchos::RCP<Epetra_Vector> POROELAST::PoroBase::StructureToFluidField(
    Teuchos::RCP<const Epetra_Vector> iv)
{
  if (matchinggrid_)
  {
    if (submeshes_)
      return coupling_fluid_structure_->MasterToSlave(psi_extractor_->ExtractCondVector(iv));
    else
      return coupling_fluid_structure_->MasterToSlave(iv);
  }
  else
  {
    Teuchos::RCP<const Epetra_Vector> mv = volcoupl_->ApplyVectorMapping21(iv);

    Teuchos::RCP<Epetra_Vector> sv =
        LINALG::CreateVector(*(FluidField()->VelPresSplitter()->OtherMap()));

    std::copy(mv->Values(),
        mv->Values() + (static_cast<ptrdiff_t>(mv->MyLength() * mv->NumVectors())), sv->Values());
    return sv;
  }
}

void POROELAST::PoroBase::SetStructSolution()
{
  Teuchos::RCP<const Epetra_Vector> dispnp;
  // apply current displacements and velocities to the fluid field
  if (StructureField()->HaveConstraint())
  {
    // displacement vector without lagrange-multipliers
    dispnp = cond_splitter_->ExtractCondVector(StructureField()->Dispnp());
  }
  else
    dispnp = StructureField()->Dispnp();

  Teuchos::RCP<const Epetra_Vector> velnp = StructureField()->Velnp();

  // transfer the current structure displacement to the fluid field
  Teuchos::RCP<Epetra_Vector> structdisp = StructureToFluidField(dispnp);
  FluidField()->ApplyMeshDisplacement(structdisp);

  // transfer the current structure velocity to the fluid field
  Teuchos::RCP<Epetra_Vector> structvel = StructureToFluidField(velnp);
  FluidField()->ApplyMeshVelocity(structvel);
}

void POROELAST::PoroBase::SetFluidSolution()
{
  if (matchinggrid_)
  {
    StructureField()->Discretization()->SetState(1, "fluidvel", FluidField()->Velnp());
  }
  else
  {
    StructureField()->Discretization()->SetState(
        1, "fluidvel", volcoupl_->ApplyVectorMapping12(FluidField()->Velnp()));
  }
}

void POROELAST::PoroBase::TimeLoop()
{
  while (NotFinished())
  {
    // solve one time step
    DoTimeStep();
  }
}

void POROELAST::PoroBase::Output(bool forced_writerestart)
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField()->StatisticsAndOutput();
  StructureField()->Output(forced_writerestart);
}

void POROELAST::PoroBase::SetupCoupling()
{
  // get discretizations
  Teuchos::RCP<DRT::Discretization> structdis = StructureField()->Discretization();
  Teuchos::RCP<DRT::Discretization> fluiddis = FluidField()->Discretization();

  // if one discretization is a subset of the other, they will differ in node number (and element
  // number) we assume matching grids for the overlapping part here
  const Epetra_Map* structnoderowmap = structdis->NodeRowMap();
  const Epetra_Map* fluidnoderowmap = fluiddis->NodeRowMap();

  const int numglobalstructnodes = structnoderowmap->NumGlobalElements();
  const int numglobalfluidnodes = fluidnoderowmap->NumGlobalElements();

  if (matchinggrid_)
  {
    // check for submeshes
    submeshes_ = (numglobalstructnodes != numglobalfluidnodes);
  }
  else
    submeshes_ = false;

  const int ndim = DRT::Problem::Instance()->NDim();
  const int numglobalstructdofs = structdis->DofRowMap()->NumGlobalElements();
  if (numglobalstructdofs == numglobalstructnodes * ndim)
    porosity_dof_ = false;
  else
  {
    porosity_dof_ = true;
    porosity_splitter_ = POROELAST::UTILS::BuildPoroSplitter(StructureField()->Discretization());
  }

  coupling_fluid_structure_ = Teuchos::rcp(new ADAPTER::Coupling());
  int ndof = ndim;

  // if the porosity is a primary variable, we get one more dof
  if (porosity_dof_) ndof++;

  if (matchinggrid_)
  {
    if (submeshes_)
    {
      // for submeshes we only couple a part of the structure disc. with the fluid disc.
      // we use the fact, that we have matching grids and matching gids
      // The node matching search tree is used to find matching structure and fluid nodes.
      // Note, that the structure discretization must be the bigger one (because it is the
      // masterdis).
      coupling_fluid_structure_->SetupCoupling(
          *structdis, *fluiddis, *fluidnoderowmap, *fluidnoderowmap, ndof, false);
    }
    else
    {
      // matching grid case: we rely on that the cloning strategy build the fluid node map with
      // equal node gids as the structure and also identical parallel distribution. Hence, we do not
      // use the node search tree here and use the same fluid node map also as permuted map.
      coupling_fluid_structure_->SetupCoupling(
          *structdis, *fluiddis, *structnoderowmap, *fluidnoderowmap, *fluidnoderowmap, ndof);
    }

    FluidField()->SetMeshMap(coupling_fluid_structure_->SlaveDofMap());

    if (submeshes_)
      psi_extractor_ = Teuchos::rcp(new LINALG::MapExtractor(
          *StructureField()->DofRowMap(), coupling_fluid_structure_->MasterDofMap()));
  }
  else
  {
    FluidField()->SetMeshMap(FluidField()->VelPresSplitter()->OtherMap());
  }
}

void POROELAST::PoroBase::ReplaceDofSets()
{
  // the problem is two way coupled, thus each discretization must know the other discretization

  // get discretizations
  Teuchos::RCP<DRT::Discretization> structdis = StructureField()->Discretization();
  Teuchos::RCP<DRT::Discretization> fluiddis = FluidField()->Discretization();

  /* When coupling porous media with a pure structure we will have two discretizations
   * of different size. In this case we need a special proxy, which can handle submeshes.
   */
  if (submeshes_)
  {
    Teuchos::RCP<DRT::DofSetGIDBasedWrapper> structsubdofset =
        Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(structdis, structdis->GetDofSetProxy()));
    Teuchos::RCP<DRT::DofSetGIDBasedWrapper> fluidsubdofset =
        Teuchos::rcp(new DRT::DofSetGIDBasedWrapper(fluiddis, fluiddis->GetDofSetProxy()));

    fluiddis->ReplaceDofSet(1, structsubdofset);
    structdis->ReplaceDofSet(1, fluidsubdofset);
  }
  else
  {
    // build a proxy of the structure discretization for the fluid field
    Teuchos::RCP<DRT::DofSetInterface> structdofsetproxy = structdis->GetDofSetProxy();
    // build a proxy of the fluid discretization for the structure field
    Teuchos::RCP<DRT::DofSetInterface> fluiddofsetproxy = fluiddis->GetDofSetProxy();

    fluiddis->ReplaceDofSet(1, structdofsetproxy);
    structdis->ReplaceDofSet(1, fluiddofsetproxy);
  }

  fluiddis->FillComplete(true, true, true);
  structdis->FillComplete(true, true, true);
}

void POROELAST::PoroBase::CheckForPoroConditions()
{
  std::vector<DRT::Condition*> nopencond;
  FluidField()->Discretization()->GetCondition("NoPenetration", nopencond);
  nopen_handle_ = Teuchos::rcp(new POROELAST::NoPenetrationConditionHandle(nopencond));

  part_int_cond_ = false;
  std::vector<DRT::Condition*> poroPartInt;
  FluidField()->Discretization()->GetCondition("PoroPartInt", poroPartInt);
  if (poroPartInt.size()) part_int_cond_ = true;

  pres_int_cond_ = false;
  std::vector<DRT::Condition*> poroPresInt;
  FluidField()->Discretization()->GetCondition("PoroPresInt", poroPresInt);
  if (poroPresInt.size()) pres_int_cond_ = true;
}

void POROELAST::NoPenetrationConditionHandle::BuidNoPenetrationMap(
    const Epetra_Comm& comm, Teuchos::RCP<const Epetra_Map> dofRowMap)
{
  std::vector<int> condIDs;
  std::set<int>::iterator it;
  for (it = cond_ids_->begin(); it != cond_ids_->end(); it++)
  {
    condIDs.push_back(*it);
  }
  Teuchos::RCP<Epetra_Map> nopendofmap =
      Teuchos::rcp(new Epetra_Map(-1, int(condIDs.size()), condIDs.data(), 0, comm));

  nopenetration_ = Teuchos::rcp(new LINALG::MapExtractor(*dofRowMap, nopendofmap));
}

void POROELAST::NoPenetrationConditionHandle::ApplyCondRHS(
    Teuchos::RCP<Epetra_Vector> iterinc, Teuchos::RCP<Epetra_Vector> rhs)
{
  if (has_cond_)
  {
    const Teuchos::RCP<const Epetra_Map>& nopenetrationmap = nopenetration_->Map(1);
    LINALG::ApplyDirichlettoSystem(iterinc, rhs, cond_rhs_, *nopenetrationmap);
  }
}

void POROELAST::NoPenetrationConditionHandle::Clear(POROELAST::coupltype coupltype)
{
  if (has_cond_)
  {
    cond_rhs_->PutScalar(0.0);
    cond_ids_->clear();
    switch (coupltype)
    {
      case POROELAST::fluidfluid:
        fluid_fluid_constraint_matrix_->Zero();
        cond_dofs_->PutScalar(0.0);
        break;
      case POROELAST::fluidstructure:
        fluid_structure_constraint_matrix_->Zero();
        structure_vel_constraint_matrix_->Zero();
        break;
      default:
        cond_dofs_->PutScalar(0.0);
        fluid_fluid_constraint_matrix_->Zero();
        fluid_structure_constraint_matrix_->Zero();
        structure_vel_constraint_matrix_->Zero();
        break;
    }
  }
}

void POROELAST::NoPenetrationConditionHandle::Setup(
    Teuchos::RCP<const Epetra_Map> dofRowMap, const Epetra_Map* dofRowMapFluid)
{
  if (has_cond_)
  {
    cond_rhs_ = Teuchos::rcp(new Epetra_Vector(*dofRowMap, true));

    cond_dofs_ = Teuchos::rcp(new Epetra_Vector(*dofRowMapFluid, true));

    fluid_fluid_constraint_matrix_ =
        Teuchos::rcp(new LINALG::SparseMatrix(*dofRowMapFluid, 81, true, true));

    fluid_structure_constraint_matrix_ =
        Teuchos::rcp(new LINALG::SparseMatrix(*dofRowMapFluid, 81, true, true));

    structure_vel_constraint_matrix_ =
        Teuchos::rcp(new LINALG::SparseMatrix(*dofRowMapFluid, 81, true, true));
  }
}

Teuchos::RCP<LINALG::SparseMatrix> POROELAST::NoPenetrationConditionHandle::ConstraintMatrix(
    POROELAST::coupltype coupltype)
{
  if (has_cond_)
  {
    if (coupltype == POROELAST::fluidfluid)
      return fluid_fluid_constraint_matrix_;
    else if (coupltype == POROELAST::fluidstructure)
      return fluid_structure_constraint_matrix_;
  }
  return Teuchos::null;
}

Teuchos::RCP<LINALG::SparseMatrix>
POROELAST::NoPenetrationConditionHandle::StructVelConstraintMatrix(POROELAST::coupltype coupltype)
{
  if (has_cond_)
  {
    if (coupltype == POROELAST::fluidfluid)
      return Teuchos::null;
    else if (coupltype == POROELAST::fluidstructure)
      return structure_vel_constraint_matrix_;
  }
  return Teuchos::null;
}
