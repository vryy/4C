/*----------------------------------------------------------------------*/
/*! \file

\brief Algorithm for the calculation of biofilm growth.
       It consists of:
       - an inner timeloop (resolving fsi and scatra (in both fluid and structure)
       at fluid-dynamic time-scale
       - an outer timeloop (resolving only the biofilm growth)
       at biological time-scale

\level 3


 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_fs3i_biofilm_fsi.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fem_geometry_update_reference_config.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fs3i_biofilm_fsi_utils.hpp"
#include "4C_fsi_monolithicfluidsplit.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

// #define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::BiofilmFSI::BiofilmFSI(const Epetra_Comm& comm) : PartFS3I1Wc(comm), comm_(comm)
{
  // has to stay empty
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::Init()
{
  if (comm_.MyPID() == 0)
    std::cout << "\n WARNING ! The implementation of BiofilmFSI is not well tested,\n"
                 " buggy, and introduction of just Init(...) and Setup() in commit\n"
                 " to revision 22366 led to differing results slightly above the\n"
                 " convergence tolerance. Rework on this problem type is necessary!\n\n"
              << std::endl;

  // call Init() in base class
  FS3I::PartFS3I1Wc::Init();

  //---------------------------------------------------------------------
  // set up struct ale
  //---------------------------------------------------------------------

  // this algorithm needs an ale discretization also for the structure in order to be able to handle
  // the growth
  Global::Problem* problem = Global::Problem::Instance();
  problem->GetDis("structale")->fill_complete();

  // create struct ale elements if not yet existing
  Teuchos::RCP<Discret::Discretization> structaledis = problem->GetDis("structale");
  Teuchos::RCP<Discret::Discretization> structdis = problem->GetDis("structure");

  // time measurement
  Teuchos::Time time("biofilm_fsi_Init", true);

  if (structaledis->NumGlobalNodes() == 0)
  {
    Teuchos::RCP<Core::FE::DiscretizationCreator<ALE::UTILS::AleCloneStrategy>> alecreator =
        Teuchos::rcp(new Core::FE::DiscretizationCreator<ALE::UTILS::AleCloneStrategy>());
    alecreator->create_matching_discretization(structdis, structaledis, 11);
    structaledis->fill_complete();
  }
  if (comm_.MyPID() == 0)
  {
    std::cout << "Created discretization " << (structaledis->Name())
              << " as a clone of discretization " << (structdis->Name()) << " in...."
              << time.totalElapsedTime(true) << " secs\n\n";
  }

  // ask base algorithm for the ale time integrator
  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  Teuchos::RCP<Adapter::AleBaseAlgorithm> ale =
      Teuchos::rcp(new Adapter::AleBaseAlgorithm(fsidyn, structaledis));
  ale_ = Teuchos::rcp_dynamic_cast<Adapter::AleFsiWrapper>(ale->ale_field());
  if (ale_ == Teuchos::null)
    FOUR_C_THROW("cast from Adapter::Ale to Adapter::AleFsiWrapper failed");


  //---------------------------------------------------------------------
  // getting and initializing problem-specific parameters
  //---------------------------------------------------------------------

  const Teuchos::ParameterList& biofilmcontrol =
      Global::Problem::Instance()->biofilm_control_params();

  // make sure that initial time derivative of concentration is not calculated
  // automatically (i.e. field-wise)
  const Teuchos::ParameterList& scatradyn =
      Global::Problem::Instance()->scalar_transport_dynamic_params();
  if (Core::UTILS::IntegralValue<int>(scatradyn, "SKIPINITDER") == false)
    FOUR_C_THROW(
        "Initial time derivative of phi must not be calculated automatically -> set SKIPINITDER to "
        "false");

  // fsi parameters
  dt_fsi_ = fsidyn.get<double>("TIMESTEP");
  nstep_fsi_ = fsidyn.get<int>("NUMSTEP");
  maxtime_fsi_ = fsidyn.get<double>("MAXTIME");
  step_fsi_ = 0;
  time_fsi_ = 0.;

  // growth parameters
  dt_bio_ = biofilmcontrol.get<double>("BIOTIMESTEP");
  nstep_bio_ = biofilmcontrol.get<int>("BIONUMSTEP");
  fluxcoef_ = biofilmcontrol.get<double>("FLUXCOEF");
  normforceposcoef_ = biofilmcontrol.get<double>("NORMFORCEPOSCOEF");
  normforcenegcoef_ = biofilmcontrol.get<double>("NORMFORCENEGCOEF");
  tangoneforcecoef_ = biofilmcontrol.get<double>("TANGONEFORCECOEF");
  tangtwoforcecoef_ = biofilmcontrol.get<double>("TANGTWOFORCECOEF");
  step_bio_ = 0;
  time_bio_ = 0.;

  // total time
  time_ = 0.;

  // safety checks
  if (volume_fieldcouplings_[0] == Inpar::FS3I::coupling_nonmatch or
      volume_fieldcouplings_[1] == Inpar::FS3I::coupling_nonmatch)
    FOUR_C_THROW("Mortar volume coupling is yet not implemented for biofilm-fs3i.");
  if (!problem->GetDis("scatra1")->GetCondition("ScaTraFluxCalc") or
      !problem->GetDis("scatra2")->GetCondition("ScaTraFluxCalc"))
    FOUR_C_THROW(
        "Fluid-scatra and solid-scatra discretizations must have boundary conditions for flux "
        "calculation at FSI interface!");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::Setup()
{
  // call Setup() in base class
  FS3I::PartFS3I1Wc::Setup();

  Teuchos::RCP<Discret::Discretization> structaledis =
      Global::Problem::Instance()->GetDis("structale");

  // create fluid-ALE Dirichlet Map Extractor for FSI step
  ale_->SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_std);

  // create fluid-ALE Dirichlet Map Extractor for growth step
  ale_->SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_biofilm, ale_->Interface());

  // create fluid-ALE Dirichlet Map Extractor for growth step
  fsi_->ale_field()->SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_std, Teuchos::null);

  // create fluid-ALE Dirichlet Map Extractor for FSI step
  fsi_->ale_field()->SetupDBCMapEx(
      ALE::UTILS::MapExtractor::dbc_set_biofilm, fsi_->ale_field()->Interface());

  //---------------------------------------------------------------------
  // set up couplings
  //---------------------------------------------------------------------

  const std::string condname = "FSICoupling";
  const int ndim = Global::Problem::Instance()->NDim();

  // set up ale-fluid couplings
  icoupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoupfa_->setup_condition_coupling(*(fsi_->fluid_field()->discretization()),
      (fsi_->fluid_field()->Interface()->FSICondMap()), *(fsi_->ale_field()->discretization()),
      (fsi_->ale_field()->Interface()->FSICondMap()), condname, ndim);
  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fsi_->fluid_field()->discretization()->NodeRowMap();
  const Epetra_Map* fluidalenodemap = fsi_->ale_field()->discretization()->NodeRowMap();
  coupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupfa_->setup_coupling(*(fsi_->fluid_field()->discretization()),
      *(fsi_->ale_field()->discretization()), *fluidnodemap, *fluidalenodemap, ndim);

  // set up structale-structure couplings
  icoupsa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoupsa_->setup_condition_coupling(*(fsi_->structure_field()->discretization()),
      fsi_->structure_field()->Interface()->FSICondMap(), *structaledis,
      ale_->Interface()->FSICondMap(), condname, ndim);
  // the structure-structale coupling always matches
  const Epetra_Map* structurenodemap = fsi_->structure_field()->discretization()->NodeRowMap();
  const Epetra_Map* structalenodemap = structaledis->NodeRowMap();
  coupsa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupsa_->setup_coupling(*(fsi_->structure_field()->discretization()), *structaledis,
      *structurenodemap, *structalenodemap, ndim);

  /// do we need this? What's for???
  fsi_->fluid_field()->SetMeshMap(coupfa_->MasterDofMap());

  idispn_ = fsi_->fluid_field()->extract_interface_veln();
  idispnp_ = fsi_->fluid_field()->extract_interface_veln();
  iveln_ = fsi_->fluid_field()->extract_interface_veln();

  struidispn_ = fsi_->structure_field()->extract_interface_dispn();
  struidispnp_ = fsi_->structure_field()->extract_interface_dispn();
  struiveln_ = fsi_->structure_field()->extract_interface_dispn();

  struct_growth_disp_ = AleToStructField(ale_->WriteAccessDispnp());
  fluid_growth_disp_ = ale_to_fluid_field(fsi_->ale_field()->WriteAccessDispnp());
  scatra_struct_growth_disp_ = Teuchos::rcp(new Epetra_MultiVector(
      *(scatravec_[1]->ScaTraField()->discretization())->NodeRowMap(), 3, true));
  scatra_fluid_growth_disp_ = Teuchos::rcp(new Epetra_MultiVector(
      *(scatravec_[0]->ScaTraField()->discretization())->NodeRowMap(), 3, true));

  idispn_->PutScalar(0.0);
  idispnp_->PutScalar(0.0);
  iveln_->PutScalar(0.0);

  struidispn_->PutScalar(0.0);
  struidispnp_->PutScalar(0.0);
  struiveln_->PutScalar(0.0);

  struct_growth_disp_->PutScalar(0.0);
  fluid_growth_disp_->PutScalar(0.0);
  scatra_struct_growth_disp_->PutScalar(0.0);
  scatra_fluid_growth_disp_->PutScalar(0.0);

  norminflux_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  normtraction_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  tangtractionone_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  tangtractiontwo_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::Timeloop()
{
  check_is_init();
  check_is_setup();


  const Teuchos::ParameterList& biofilmcontrol =
      Global::Problem::Instance()->biofilm_control_params();
  const int biofilmgrowth = Core::UTILS::IntegralValue<int>(biofilmcontrol, "BIOFILMGROWTH");
  const int outputgmsh_ = Core::UTILS::IntegralValue<int>(biofilmcontrol, "OUTPUT_GMSH");

  std::cout << std::endl << "--------------SIMULATION PARAMETERS-----------------" << std::endl;
  std::cout << "FSI TIMESTEP = " << dt_fsi_ << "; FSI NUMSTEP = " << nstep_fsi_ << std::endl;
  std::cout << "BIO TIMESTEP = " << dt_bio_ << "; BIO NUMSTEP = " << nstep_bio_ << std::endl;
  std::cout << "FLUXCOEF = " << fluxcoef_ << ";" << std::endl;
  std::cout << "NORMFORCEPOSCOEF = " << normforceposcoef_
            << "; NORMFORCENEGCOEF = " << normforcenegcoef_ << std::endl;
  std::cout << "TANGONEFORCECOEF = " << tangoneforcecoef_
            << "; TANGTWOFORCECOEF = " << tangtwoforcecoef_ << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;

  if (biofilmgrowth)
  {
    // outer loop for biofilm growth
    while (step_bio_ <= nstep_bio_)
    {
      // update step and time
      step_bio_++;
      time_bio_ += dt_bio_;
      time_ = time_bio_ + time_fsi_;

      if (step_bio_ == 1 && outputgmsh_)
      {
        StructGmshOutput();
        FluidGmshOutput();
      }

      // inner loop for fsi and scatra
      InnerTimeloop();

      // gmsh output only if requested
      if (outputgmsh_)
      {
        StructGmshOutput();
        FluidGmshOutput();
      }

      if (Comm().MyPID() == 0)
      {
        std::cout << "\n***********************\n     GROWTH STEP \n***********************\n";
        printf(" growth step = %3d   \n", step_bio_);
        printf(" Total time = %3f   \n", time_);
      }

      // compute interface displacement and velocity
      compute_interface_vectors(idispnp_, iveln_, struidispnp_, struiveln_);

      // do all the settings and solve the fluid on a deforming mesh
      FluidAleSolve();

      // do all the settings and solve the structure on a deforming mesh
      StructAleSolve();

      fsi_->output();
      ScatraOutput();

      // reset step and state vectors
      fsi_->structure_field()->Reset();
      // fluid reset can be bypassed, in this way the next step starts from a solution closer to the
      // final one fsi_->fluid_field()->Reset(false, false, step_bio);
      fsi_->ale_field()->Reset();

      fsi_->ale_field()->create_system_matrix(fsi_->ale_field()->Interface());
    }
  }

  if (!biofilmgrowth)
  {
    InnerTimeloop();

    // gmsh output only if requested
    if (outputgmsh_)
    {
      StructGmshOutput();
      FluidGmshOutput();
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::InnerTimeloop()
{
  // initialize time and step each time we enter the innerloop
  double t = 0.;
  step_fsi_ = 0;
  // initialize fluxes and tractions each time we enter the innerloop
  norminflux_->PutScalar(0.0);
  normtraction_->PutScalar(0.0);
  tangtractionone_->PutScalar(0.0);
  tangtractiontwo_->PutScalar(0.0);

  // output of initial state
  //  ScatraOutput();

  fsi_->PrepareTimeloop();

  // Calculation of growth can be based both on values averaged during the inner timeloop
  // (in this case for the time being it takes in account also the initial transient state!),
  // or only on the last values coming from the fsi-scatra simulation
  const Teuchos::ParameterList& biofilmcontrol =
      Global::Problem::Instance()->biofilm_control_params();
  const int avgrowth = Core::UTILS::IntegralValue<int>(biofilmcontrol, "AVGROWTH");
  // in case of averaged values we need temporary variables
  Teuchos::RCP<Epetra_Vector> normtempinflux_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  Teuchos::RCP<Epetra_Vector> normtemptraction_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  Teuchos::RCP<Epetra_Vector> tangtemptractionone_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  Teuchos::RCP<Epetra_Vector> tangtemptractiontwo_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->structure_field()->discretization()->NodeRowMap())));
  normtempinflux_->PutScalar(0.0);
  normtemptraction_->PutScalar(0.0);
  tangtemptractionone_->PutScalar(0.0);
  tangtemptractiontwo_->PutScalar(0.0);

  while (step_fsi_ < nstep_fsi_ and t + 1e-10 * dt_fsi_ < maxtime_fsi_)
  {
    step_fsi_++;
    t += dt_fsi_;

    fsi_->prepare_time_step();
    fsi_->TimeStep(fsi_);

    constexpr bool force_prepare = false;
    fsi_->prepare_output(force_prepare);
    fsi_->update();

    SetFSISolution();

    if (Comm().MyPID() == 0)
    {
      std::cout << "\n***********************\n GAS TRANSPORT SOLVER \n***********************\n";
    }

    // first scatra field is associated with fluid, second scatra field is
    // associated with structure

    bool stopnonliniter = false;
    int itnum = 0;

    prepare_time_step();

    while (stopnonliniter == false)
    {
      scatra_evaluate_solve_iter_update();
      itnum++;
      if (scatra_convergence_check(itnum)) break;
    }

    // calculation of the flux at the interface based on normal influx values before time shift of
    // results is performed in Update
    Teuchos::RCP<Epetra_MultiVector> strufluxn =
        scatravec_[1]->ScaTraField()->CalcFluxAtBoundary(false);

    UpdateScatraFields();

    // this is necessary because we want to write all the steps except the last one
    // the last one will be written only after the calculation of the growth
    // in this way also the displacement due to growth is written
    if (step_fsi_ < nstep_fsi_ and t + 1e-10 * dt_fsi_ < maxtime_fsi_)
    {
      fsi_->output();
      ScatraOutput();
    }

    // access structure discretization
    Teuchos::RCP<Discret::Discretization> strudis = fsi_->structure_field()->discretization();

    // recovery of forces at the interface nodes based on lagrange multipliers values
    // lambda_ is defined only at the interface, while lambdafull on the entire fluid/structure
    // field.
    Teuchos::RCP<Epetra_Vector> lambda_;
    Teuchos::RCP<Epetra_Vector> lambdafull;

    // at the purpose to compute lambdafull, it is necessary to know which coupling algorithm is
    // used however the imposition of a Dirichlet condition on the interface produce wrong lambda_
    // when structuresplit is used
    const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
    int coupling = Teuchos::getIntegralValue<int>(fsidyn, "COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit)
    {
      Teuchos::RCP<Epetra_Vector> lambdafluid = fsi_->GetLambda();
      lambda_ = fsi_->fluid_to_struct(lambdafluid);
    }
    else if (coupling == fsi_iter_monolithicstructuresplit)
    {
      lambda_ = fsi_->GetLambda();
    }

    lambdafull = fsi_->structure_field()->Interface()->InsertFSICondVector(lambda_);

    // calculate interface normals in deformed configuration
    Teuchos::RCP<Epetra_Vector> nodalnormals =
        Teuchos::rcp(new Epetra_Vector(*(strudis->dof_row_map())));

    Teuchos::ParameterList eleparams;
    eleparams.set("action", "calc_cur_nodal_normals");
    strudis->ClearState();
    strudis->set_state("displacement", fsi_->structure_field()->Dispnp());
    strudis->evaluate_condition(eleparams, Teuchos::null, Teuchos::null, nodalnormals,
        Teuchos::null, Teuchos::null, "FSICoupling");
    strudis->ClearState();

    const Epetra_Map* dofrowmap = strudis->dof_row_map();
    const Epetra_Map* noderowmap = strudis->NodeRowMap();
    Teuchos::RCP<Epetra_MultiVector> lambdanode =
        Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));

    // lagrange multipliers defined on a nodemap are necessary
    for (int lnodeid = 0; lnodeid < strudis->NumMyRowNodes(); lnodeid++)
    {
      // get the processor local node
      Core::Nodes::Node* lnode = strudis->lRowNode(lnodeid);
      // get the dofs of the node
      //      std::vector<int> dofs= strudis->Dof(lnode);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = strudis->Dof(0, lnode);

      for (unsigned index = 0; index < nodedofset.size(); ++index)
      {
        // get global id of the dof
        int gid = nodedofset[index];
        // get local id of the dof
        int lid = dofrowmap->LID(gid);

        //        int lnodeid=lnode->LID();

        double lambdai = (*lambdafull)[lid];
        int lnodeid = noderowmap->LID(lnode->Id());
        ((*lambdanode)(index))->ReplaceMyValues(1, &lambdai, &lnodeid);
      }
    }
    // loop over all local interface nodes of structure discretization
    Teuchos::RCP<Epetra_Map> condnodemap =
        Core::Conditions::ConditionNodeRowMap(*strudis, "FSICoupling");
    for (int nodei = 0; nodei < condnodemap->NumMyElements(); nodei++)
    {
      // Here we rely on the fact that the structure scatra discretization is a clone of the
      // structure mesh

      // get the processor's local node with the same lnodeid
      int gnodeid = condnodemap->GID(nodei);
      Core::Nodes::Node* strulnode = strudis->gNode(gnodeid);
      // get the degrees of freedom associated with this node
      std::vector<int> strunodedofs = strudis->Dof(0, strulnode);
      // determine number of space dimensions
      const int numdim = ((int)strunodedofs.size());

      std::vector<int> doflids(numdim);
      double temp = 0.;
      std::vector<double> unitnormal(3);
      for (int i = 0; i < numdim; ++i)
      {
        doflids[i] = strudis->dof_row_map()->LID(strunodedofs[i]);
        unitnormal[i] = (*nodalnormals)[doflids[i]];
        temp += unitnormal[i] * unitnormal[i];
      }
      double unitnormalabsval = sqrt(temp);
      int lnodeid = strudis->NodeRowMap()->LID(gnodeid);

      // compute average unit nodal normal
      std::vector<double> Values(numdim);
      for (int j = 0; j < numdim; ++j)
      {
        unitnormal[j] /= unitnormalabsval;
      }

      // compute tangents
      std::vector<double> unittangentone(3);
      std::vector<double> unittangenttwo(3);

      // take care of special case
      double TOL = 1e-11;
      if (abs(unitnormal[0]) < TOL && abs(unitnormal[1]) < TOL)
      {
        unittangentone[0] = 1.0;
        unittangentone[1] = 0.0;
        unittangentone[2] = 0.0;

        unittangenttwo[0] = 0.0;
        unittangenttwo[1] = 1.0;
        unittangenttwo[2] = 0.0;
      }
      else
      {
        // first unit tangent
        unittangentone[0] = -unitnormal[1];
        unittangentone[1] = unitnormal[0];
        unittangentone[2] = 0.0;
        temp = 0.;
        for (int i = 0; i < numdim; ++i)
        {
          temp += unittangentone[i] * unittangentone[i];
        }
        double unittangentoneabsval = sqrt(temp);
        for (int j = 0; j < numdim; ++j)
        {
          unittangentone[j] /= unittangentoneabsval;
        }

        // second unit tangent
        unittangenttwo[0] = -unitnormal[0] * unitnormal[2];
        unittangenttwo[1] = -unitnormal[1] * unitnormal[2];
        unittangenttwo[2] = unitnormal[0] * unitnormal[0] + unitnormal[1] * unitnormal[1];
        temp = 0.;
        for (int i = 0; i < numdim; ++i)
        {
          temp += unittangenttwo[i] * unittangenttwo[i];
        }
        double unittangenttwoabsval = sqrt(temp);
        for (int j = 0; j < numdim; ++j)
        {
          unittangenttwo[j] /= unittangenttwoabsval;
        }
      }

      double tempflux = 0.0;
      double tempnormtrac = 0.0;
      double temptangtracone = 0.0;
      double temptangtractwo = 0.0;
      for (int index = 0; index < numdim; ++index)
      {
        double fluxcomp = (*((*strufluxn)(index)))[lnodeid];
        tempflux += fluxcomp * unitnormal[index];
        // for the calculation of the growth and erosion both the tangential and the normal
        // components of the forces acting on the interface are important.
        // Since probably they will have a different effect on the biofilm growth,
        // they are calculated separately and different coefficients can be used.
        double traccomp = (*((*lambdanode)(index)))[lnodeid];
        tempnormtrac += traccomp * unitnormal[index];
        temptangtracone += traccomp * unittangentone[index];
        temptangtractwo += traccomp * unittangenttwo[index];
      }

      if (avgrowth)
      {
        (*((*normtempinflux_)(0)))[lnodeid] += tempflux;
        (*((*normtemptraction_)(0)))[lnodeid] += abs(tempnormtrac);
        (*((*tangtemptractionone_)(0)))[lnodeid] += abs(temptangtracone);
        (*((*tangtemptractiontwo_)(0)))[lnodeid] += abs(temptangtractwo);
      }
      else
      {
        (*((*norminflux_)(0)))[lnodeid] = tempflux;
        (*((*normtraction_)(0)))[lnodeid] = abs(tempnormtrac);
        (*((*tangtractionone_)(0)))[lnodeid] = abs(temptangtracone);
        (*((*tangtractiontwo_)(0)))[lnodeid] = abs(temptangtractwo);
      }
    }
  }

  // here is the averaging of variables needed for biofilm growth, in case the average way was
  // chosen
  if (avgrowth)
  {
    Teuchos::RCP<Discret::Discretization> strudis = fsi_->structure_field()->discretization();

    // loop over all local interface nodes of structure discretization
    Teuchos::RCP<Epetra_Map> condnodemap =
        Core::Conditions::ConditionNodeRowMap(*strudis, "FSICoupling");
    for (int i = 0; i < condnodemap->NumMyElements(); i++)
    {
      // get the processor's local node with the same lnodeid
      int gnodeid = condnodemap->GID(i);
      int lnodeid = strudis->NodeRowMap()->LID(gnodeid);

      (*((*norminflux_)(0)))[lnodeid] = (*((*normtempinflux_)(0)))[lnodeid] / step_fsi_;
      (*((*normtraction_)(0)))[lnodeid] = (*((*normtemptraction_)(0)))[lnodeid] / step_fsi_;
      (*((*tangtractionone_)(0)))[lnodeid] = (*((*tangtemptractionone_)(0)))[lnodeid] / step_fsi_;
      (*((*tangtractiontwo_)(0)))[lnodeid] = (*((*tangtemptractiontwo_)(0)))[lnodeid] / step_fsi_;
    }
  }

  time_fsi_ += t;
}

/*----------------------------------------------------------------------*
 | write FSI solutions into scatra discretisation             Thon 11/14|
 *----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::SetFSISolution()
{
  set_mesh_disp();
  set_velocity_fields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::compute_interface_vectors(Teuchos::RCP<Epetra_Vector> idispnp,
    Teuchos::RCP<Epetra_Vector> iveln, Teuchos::RCP<Epetra_Vector> struidispnp,
    Teuchos::RCP<Epetra_Vector> struiveln)
{
  // initialize structure interface displacement at time t^{n+1}
  // shouldn't that be zeroed?
  struidispnp->PutScalar(0.0);

  // select biofilm growth boundaries
  std::string biogrcondname = "BioGrCoupling";

  // set action for elements: compute normal vectors at nodes (for reference configuration)
  Teuchos::RCP<Discret::Discretization> strudis = fsi_->structure_field()->discretization();
  Teuchos::RCP<Epetra_Vector> nodalnormals =
      Teuchos::rcp(new Epetra_Vector(*(strudis->dof_row_map())));
  Teuchos::ParameterList eleparams;
  eleparams.set("action", "calc_ref_nodal_normals");
  strudis->evaluate_condition(eleparams, Teuchos::null, Teuchos::null, nodalnormals, Teuchos::null,
      Teuchos::null, biogrcondname);

  // select row map with nodes from condition
  Teuchos::RCP<Epetra_Map> condnodemap =
      Core::Conditions::ConditionNodeRowMap(*strudis, biogrcondname);

  // loop all conditioned nodes
  for (int i = 0; i < condnodemap->NumMyElements(); ++i)
  {
    int nodegid = condnodemap->GID(i);
    if (strudis->HaveGlobalNode(nodegid) == false) FOUR_C_THROW("node not found on this proc");
    Core::Nodes::Node* actnode = strudis->gNode(nodegid);
    std::vector<int> globaldofs = strudis->Dof(0, actnode);
    const int numdim = (int)(globaldofs.size());

    // extract averaged nodal normal and compute its absolute value
    std::vector<double> unitnormal(numdim);
    double temp = 0.;
    for (int j = 0; j < numdim; ++j)
    {
      unitnormal[j] = (*nodalnormals)[strudis->dof_row_map()->LID(globaldofs[j])];
      temp += unitnormal[j] * unitnormal[j];
    }
    double unitnormalabsval = sqrt(temp);
    int lnodeid = strudis->NodeRowMap()->LID(nodegid);
    double influx = (*norminflux_)[lnodeid];
    double normforces = (*normtraction_)[lnodeid];
    double tangoneforce = (*tangtractionone_)[lnodeid];
    double tangtwoforce = (*tangtractiontwo_)[lnodeid];

    // compute average unit nodal normal and "interface velocity"
    std::vector<double> Values(numdim, 0);

    for (int j = 0; j < numdim; ++j)
    {
      unitnormal[j] /= unitnormalabsval;

      // Traction and compression probably have different effect on biofilm growth -->
      // different coefficients can be used
      double normforcecoef_;
      if (normforces > 0)
        normforcecoef_ = normforceposcoef_;
      else
        normforcecoef_ = normforcenegcoef_;

      // explanation of signs present in the following phenomenological laws
      // influx<0     --> growth  -->  - fluxcoef_
      // normforces>0 --> erosion -->  - normforcecoef_
      // tangforces>0 --> erosion -->  - tangforcecoef_
      // for pseudo-3D problems the second tangent should not be taken in account
      Values[j] = -fluxcoef_ * influx * unitnormal[j] -
                  normforcecoef_ * normforces * unitnormal[j] -
                  tangoneforcecoef_ * tangoneforce * unitnormal[j] -
                  tangtwoforcecoef_ * tangtwoforce * unitnormal[j];
    }

    int error = struiveln_->ReplaceGlobalValues(numdim, Values.data(), globaldofs.data());
    if (error > 0) FOUR_C_THROW("Could not insert values into vector struiveln_: error %d", error);
  }

  struidispnp->Update(dt_bio_, *struiveln_, 0.0);

  Teuchos::RCP<Epetra_Vector> fluididisp = fsi_->struct_to_fluid(struidispnp);
  idispnp->Update(1.0, *fluididisp, 0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::FluidAleSolve()
{
  Teuchos::RCP<Discret::Discretization> fluidaledis =
      fsi_->ale_field()->write_access_discretization();

  // if we have values at the fluid interface we need to apply them
  if (idispnp_ != Teuchos::null)
  {
    fsi_->ale_field()->apply_interface_displacements(FluidToAle(idispnp_));
  }

  fsi_->ale_field()->create_system_matrix(Teuchos::null);
  fsi_->ale_field()->Evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_biofilm);
  int error = fsi_->ale_field()->Solve();
  if (error == 1) FOUR_C_THROW("Could not solve fluid ALE in biofilm FS3I!");
  fsi_->ale_field()->UpdateIter();

  // change nodes reference position of the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp =
      ale_to_fluid_field(fsi_->ale_field()->WriteAccessDispnp());
  Teuchos::RCP<Discret::Discretization> fluiddis = fsi_->fluid_field()->discretization();
  Core::Geo::update_reference_config_with_disp(fluiddis, fluiddisp);



  // change nodes reference position also for the fluid ale field
  Teuchos::RCP<Epetra_Vector> fluidaledisp = fsi_->ale_field()->WriteAccessDispnp();
  Core::Geo::update_reference_config_with_disp(fluidaledis, fluidaledisp);

  // change nodes reference position also for scatra fluid field
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[0];
  Teuchos::RCP<Discret::Discretization> scatradis = scatra->ScaTraField()->discretization();
  FS3I::BioFilm::UTILS::ScatraChangeConfig(scatradis, fluiddis, fluiddisp);

  // set the total displacement due to growth for output reasons
  // fluid
  fluid_growth_disp_->Update(1.0, *fluiddisp, 1.0);
  fsi_->fluid_field()->SetFldGrDisp(fluid_growth_disp_);
  // fluid scatra
  VecToScatravec(scatradis, fluid_growth_disp_, scatra_fluid_growth_disp_);
  scatra->ScaTraField()->SetScFldGrDisp(scatra_fluid_growth_disp_);

  // computation of fluid solution
  // fluid_->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::StructAleSolve()
{
  Teuchos::RCP<Discret::Discretization> structaledis = ale_->write_access_discretization();

  // if we have values at the structure interface we need to apply them
  if (struidispnp_ != Teuchos::null)
  {
    ale_->apply_interface_displacements(StructToAle(struidispnp_));
  }

  ale_->create_system_matrix(Teuchos::null);
  ale_->Evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_biofilm);
  int error = ale_->Solve();
  if (error == 1) FOUR_C_THROW("Could not solve fluid ALE in biofilm FS3I!");
  ale_->UpdateIter();

  // change nodes reference position of the structure field
  Teuchos::RCP<Epetra_Vector> structdisp = AleToStructField(ale_->WriteAccessDispnp());
  Teuchos::RCP<Discret::Discretization> structdis = fsi_->structure_field()->discretization();
  Core::Geo::update_reference_config_with_disp(structdis, structdisp);
  structdis->fill_complete(false, true, true);

  // change nodes reference position also for the struct ale field
  Core::Geo::update_reference_config_with_disp(structaledis, ale_->WriteAccessDispnp());

  // change nodes reference position also for scatra structure field
  Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> struscatra = scatravec_[1];
  Teuchos::RCP<Discret::Discretization> struscatradis = struscatra->ScaTraField()->discretization();
  FS3I::BioFilm::UTILS::ScatraChangeConfig(struscatradis, structdis, structdisp);

  // set the total displacement due to growth for output reasons
  // structure
  struct_growth_disp_->Update(1.0, *structdisp, 1.0);
  fsi_->structure_field()->SetStrGrDisp(struct_growth_disp_);
  // structure scatra
  VecToScatravec(struscatradis, struct_growth_disp_, scatra_struct_growth_disp_);
  struscatra->ScaTraField()->SetScStrGrDisp(scatra_struct_growth_disp_);

  // computation of structure solution
  // structure_->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::ale_to_fluid_field(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::AleToStructField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::AleToStructField(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupsa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::StructToAle(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupsa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::VecToScatravec(Teuchos::RCP<Discret::Discretization> scatradis,
    Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<Epetra_MultiVector> scatravec)
{
  // define error variable
  int err(0);

  // loop over all local nodes of scatra discretization
  for (int lnodeid = 0; lnodeid < scatradis->NumMyRowNodes(); lnodeid++)
  {
    // determine number of space dimensions
    const int numdim = Global::Problem::Instance()->NDim();

    for (int index = 0; index < numdim; ++index)
    {
      double vecval = (*vec)[index + numdim * lnodeid];

      // insert value into node-based vector
      err = scatravec->ReplaceMyValue(lnodeid, index, vecval);

      if (err != 0) FOUR_C_THROW("Error while inserting value into vector scatravec!");
    }

    // for 1- and 2-D problems: set all unused vector components to zero
    for (int index = numdim; index < 3; ++index)
    {
      err = scatravec->ReplaceMyValue(lnodeid, index, 0.0);
      if (err != 0) FOUR_C_THROW("Error while inserting value into vector scatravec!");
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::StructGmshOutput()
{
  const Teuchos::RCP<Discret::Discretization> structdis = fsi_->structure_field()->discretization();
  const Teuchos::RCP<Discret::Discretization> structaledis = ale_->write_access_discretization();
  Teuchos::RCP<Discret::Discretization> struscatradis =
      scatravec_[1]->ScaTraField()->discretization();

  const std::string filename = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles("struct",
      structdis->Writer()->Output()->FileName(), step_bio_, 701, false, structdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  Teuchos::RCP<const Epetra_Vector> structdisp = fsi_->structure_field()->Dispn();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "struct displacement \" {" << std::endl;
    // draw vector field 'struct displacement' for every element
    Core::IO::Gmsh::VectorFieldDofBasedToGmsh(structdis, structdisp, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<const Epetra_Vector> structaledisp = ale_->Dispnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "struct ale displacement \" {" << std::endl;
    // draw vector field 'struct ale displacement' for every element
    Core::IO::Gmsh::VectorFieldDofBasedToGmsh(structaledis, structaledisp, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<const Epetra_Vector> structphi = scatravec_[1]->ScaTraField()->Phinp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "struct phi \" {" << std::endl;
    // draw vector field 'struct phi' for every element
    Core::IO::Gmsh::ScalarFieldToGmsh(struscatradis, structphi, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::FluidGmshOutput()
{
  const Teuchos::RCP<Discret::Discretization> fluiddis = fsi_->fluid_field()->discretization();
  const Teuchos::RCP<Discret::Discretization> fluidaledis =
      fsi_->ale_field()->write_access_discretization();
  Teuchos::RCP<Discret::Discretization> fluidscatradis =
      scatravec_[0]->ScaTraField()->discretization();

  const std::string filenamefluid = Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles("fluid",
      fluiddis->Writer()->Output()->FileName(), step_bio_, 701, false, fluiddis->Comm().MyPID());
  std::ofstream gmshfilecontent(filenamefluid.c_str());

  Teuchos::RCP<const Epetra_Vector> fluidvel = fsi_->fluid_field()->Velnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "fluid velocity \" {" << std::endl;
    // draw vector field 'fluid velocity' for every element
    Core::IO::Gmsh::VectorFieldDofBasedToGmsh(fluiddis, fluidvel, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<Epetra_Vector> fluidaledisp = fsi_->ale_field()->WriteAccessDispnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "fluid ale displacement \" {" << std::endl;
    // draw vector field 'fluid ale displacement' for every element
    Core::IO::Gmsh::VectorFieldDofBasedToGmsh(fluidaledis, fluidaledisp, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<Epetra_Vector> fluidphi = scatravec_[0]->ScaTraField()->Phinp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "fluid phi \" {" << std::endl;
    // draw vector field 'fluid phi' for every element
    Core::IO::Gmsh::ScalarFieldToGmsh(fluidscatradis, fluidphi, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();

  return;
}

FOUR_C_NAMESPACE_CLOSE
