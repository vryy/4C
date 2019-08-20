/*----------------------------------------------------------------------*/
/*! \file

\brief Algorithm for the calculation of biofilm growth.
       It consists of:
       - an inner timeloop (resolving fsi and scatra (in both fluid and structure)
       at fluid-dynamic time-scale
       - an outer timeloop (resolving only the biofilm growth)
       at biological time-scale

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15249

 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "biofilm_fsi.H"
#include "biofilm_fsi_utils.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_utils_materials.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_scatra/scatra_timint_implicit.H"

//#define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::BiofilmFSI::BiofilmFSI(const Epetra_Comm& comm) : PartFS3I_1WC(comm), comm_(comm)
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
  FS3I::PartFS3I_1WC::Init();

  //---------------------------------------------------------------------
  // set up struct ale
  //---------------------------------------------------------------------

  // this algorithm needs an ale discretization also for the structure in order to be able to handle
  // the growth
  DRT::Problem* problem = DRT::Problem::Instance();
  problem->GetDis("structale")->FillComplete();

  // create struct ale elements if not yet existing
  Teuchos::RCP<DRT::Discretization> structaledis = problem->GetDis("structale");
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");

  // access the communicator for time measurement
  Epetra_Time time(comm_);

  if (structaledis->NumGlobalNodes() == 0)
  {
    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<ALE::UTILS::AleCloneStrategy>> alecreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<ALE::UTILS::AleCloneStrategy>());
    alecreator->CreateMatchingDiscretization(structdis, structaledis, 11);
    structaledis->FillComplete();
  }
  if (comm_.MyPID() == 0)
  {
    std::cout << "Created discretization " << (structaledis->Name())
              << " as a clone of discretization " << (structdis->Name()) << " in...."
              << time.ElapsedTime() << " secs\n\n";
  }

  // ask base algorithm for the ale time integrator
  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale =
      Teuchos::rcp(new ADAPTER::AleBaseAlgorithm(fsidyn, structaledis));
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFsiWrapper>(ale->AleField());
  if (ale_ == Teuchos::null) dserror("cast from ADAPTER::Ale to ADAPTER::AleFsiWrapper failed");


  //---------------------------------------------------------------------
  // getting and initializing problem-specific parameters
  //---------------------------------------------------------------------

  const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();

  // make sure that initial time derivative of concentration is not calculated
  // automatically (i.e. field-wise)
  const Teuchos::ParameterList& scatradyn =
      DRT::Problem::Instance()->ScalarTransportDynamicParams();
  if (DRT::INPUT::IntegralValue<int>(scatradyn, "SKIPINITDER") == false)
    dserror(
        "Initial time derivative of phi must not be calculated automatically -> set SKIPINITDER to "
        "false");

  // fsi parameters
  dt_fsi = fsidyn.get<double>("TIMESTEP");
  nstep_fsi = fsidyn.get<int>("NUMSTEP");
  maxtime_fsi = fsidyn.get<double>("MAXTIME");
  step_fsi = 0;
  time_fsi = 0.;

  // growth parameters
  dt_bio = biofilmcontrol.get<double>("BIOTIMESTEP");
  nstep_bio = biofilmcontrol.get<int>("BIONUMSTEP");
  fluxcoef_ = biofilmcontrol.get<double>("FLUXCOEF");
  normforceposcoef_ = biofilmcontrol.get<double>("NORMFORCEPOSCOEF");
  normforcenegcoef_ = biofilmcontrol.get<double>("NORMFORCENEGCOEF");
  tangoneforcecoef_ = biofilmcontrol.get<double>("TANGONEFORCECOEF");
  tangtwoforcecoef_ = biofilmcontrol.get<double>("TANGTWOFORCECOEF");
  step_bio = 0;
  time_bio = 0.;

  // total time
  time_ = 0.;

  // safety checks
  if (volume_fieldcouplings_[0] == INPAR::FS3I::coupling_nonmatch or
      volume_fieldcouplings_[1] == INPAR::FS3I::coupling_nonmatch)
    dserror("Mortar volume coupling is yet not implemented for biofilm-fs3i.");
  if (!problem->GetDis("scatra1")->GetCondition("ScaTraFluxCalc") or
      !problem->GetDis("scatra2")->GetCondition("ScaTraFluxCalc"))
    dserror(
        "Fluid-scatra and solid-scatra discretizations must have boundary conditions for flux "
        "calculation at FSI interface!");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::Setup()
{
  // call Setup() in base class
  FS3I::PartFS3I_1WC::Setup();

  Teuchos::RCP<DRT::Discretization> structaledis = DRT::Problem::Instance()->GetDis("structale");

  // create fluid-ALE Dirichlet Map Extractor for FSI step
  ale_->SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_std);

  // create fluid-ALE Dirichlet Map Extractor for growth step
  ale_->SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_biofilm, ale_->Interface());

  // create fluid-ALE Dirichlet Map Extractor for growth step
  fsi_->AleField()->SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_std, Teuchos::null);

  // create fluid-ALE Dirichlet Map Extractor for FSI step
  fsi_->AleField()->SetupDBCMapEx(
      ALE::UTILS::MapExtractor::dbc_set_biofilm, fsi_->AleField()->Interface());

  //---------------------------------------------------------------------
  // set up couplings
  //---------------------------------------------------------------------

  const std::string condname = "FSICoupling";
  const int ndim = DRT::Problem::Instance()->NDim();

  // set up ale-fluid couplings
  icoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupfa_->SetupConditionCoupling(*(fsi_->FluidField()->Discretization()),
      (fsi_->FluidField()->Interface()->FSICondMap()), *(fsi_->AleField()->Discretization()),
      (fsi_->AleField()->Interface()->FSICondMap()), condname, ndim);
  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fsi_->FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* fluidalenodemap = fsi_->AleField()->Discretization()->NodeRowMap();
  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_->SetupCoupling(*(fsi_->FluidField()->Discretization()),
      *(fsi_->AleField()->Discretization()), *fluidnodemap, *fluidalenodemap, ndim);

  // set up structale-structure couplings
  icoupsa_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupsa_->SetupConditionCoupling(*(fsi_->StructureField()->Discretization()),
      fsi_->StructureField()->Interface()->FSICondMap(), *structaledis,
      ale_->Interface()->FSICondMap(), condname, ndim);
  // the structure-structale coupling always matches
  const Epetra_Map* structurenodemap = fsi_->StructureField()->Discretization()->NodeRowMap();
  const Epetra_Map* structalenodemap = structaledis->NodeRowMap();
  coupsa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsa_->SetupCoupling(*(fsi_->StructureField()->Discretization()), *structaledis,
      *structurenodemap, *structalenodemap, ndim);

  /// do we need this? What's for???
  fsi_->FluidField()->SetMeshMap(coupfa_->MasterDofMap());

  idispn_ = fsi_->FluidField()->ExtractInterfaceVeln();
  idispnp_ = fsi_->FluidField()->ExtractInterfaceVeln();
  iveln_ = fsi_->FluidField()->ExtractInterfaceVeln();

  struidispn_ = fsi_->StructureField()->ExtractInterfaceDispn();
  struidispnp_ = fsi_->StructureField()->ExtractInterfaceDispn();
  struiveln_ = fsi_->StructureField()->ExtractInterfaceDispn();

  struct_growth_disp = AleToStructField(ale_->WriteAccessDispnp());
  fluid_growth_disp = AleToFluidField(fsi_->AleField()->WriteAccessDispnp());
  scatra_struct_growth_disp = Teuchos::rcp(new Epetra_MultiVector(
      *(scatravec_[1]->ScaTraField()->Discretization())->NodeRowMap(), 3, true));
  scatra_fluid_growth_disp = Teuchos::rcp(new Epetra_MultiVector(
      *(scatravec_[0]->ScaTraField()->Discretization())->NodeRowMap(), 3, true));

  idispn_->PutScalar(0.0);
  idispnp_->PutScalar(0.0);
  iveln_->PutScalar(0.0);

  struidispn_->PutScalar(0.0);
  struidispnp_->PutScalar(0.0);
  struiveln_->PutScalar(0.0);

  struct_growth_disp->PutScalar(0.0);
  fluid_growth_disp->PutScalar(0.0);
  scatra_struct_growth_disp->PutScalar(0.0);
  scatra_fluid_growth_disp->PutScalar(0.0);

  norminflux_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  normtraction_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  tangtractionone_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  tangtractiontwo_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::Timeloop()
{
  CheckIsInit();
  CheckIsSetup();


  const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();
  const int biofilmgrowth = DRT::INPUT::IntegralValue<int>(biofilmcontrol, "BIOFILMGROWTH");
  const int outputgmsh_ = DRT::INPUT::IntegralValue<int>(biofilmcontrol, "OUTPUT_GMSH");

  std::cout << std::endl << "--------------SIMULATION PARAMETERS-----------------" << std::endl;
  std::cout << "FSI TIMESTEP = " << dt_fsi << "; FSI NUMSTEP = " << nstep_fsi << std::endl;
  std::cout << "BIO TIMESTEP = " << dt_bio << "; BIO NUMSTEP = " << nstep_bio << std::endl;
  std::cout << "FLUXCOEF = " << fluxcoef_ << ";" << std::endl;
  std::cout << "NORMFORCEPOSCOEF = " << normforceposcoef_
            << "; NORMFORCENEGCOEF = " << normforcenegcoef_ << std::endl;
  std::cout << "TANGONEFORCECOEF = " << tangoneforcecoef_
            << "; TANGTWOFORCECOEF = " << tangtwoforcecoef_ << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;

  if (biofilmgrowth)
  {
    // outer loop for biofilm growth
    while (step_bio <= nstep_bio)
    {
      // update step and time
      step_bio++;
      time_bio += dt_bio;
      time_ = time_bio + time_fsi;

      if (step_bio == 1 && outputgmsh_)
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
        printf(" growth step = %3d   \n", step_bio);
        printf(" Total time = %3f   \n", time_);
      }

      // compute interface displacement and velocity
      ComputeInterfaceVectors(idispnp_, iveln_, struidispnp_, struiveln_);

      // do all the settings and solve the fluid on a deforming mesh
      FluidAleSolve();

      // do all the settings and solve the structure on a deforming mesh
      StructAleSolve();

      fsi_->Output();
      ScatraOutput();

      // reset step and state vectors
      fsi_->StructureField()->Reset();
      // fluid reset can be bypassed, in this way the next step starts from a solution closer to the
      // final one fsi_->FluidField()->Reset(false, false, step_bio);
      fsi_->AleField()->Reset();

      fsi_->AleField()->CreateSystemMatrix(fsi_->AleField()->Interface());
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
  step_fsi = 0;
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
  const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();
  const int avgrowth = DRT::INPUT::IntegralValue<int>(biofilmcontrol, "AVGROWTH");
  // in case of averaged values we need temporary variables
  Teuchos::RCP<Epetra_Vector> normtempinflux_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  Teuchos::RCP<Epetra_Vector> normtemptraction_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  Teuchos::RCP<Epetra_Vector> tangtemptractionone_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  Teuchos::RCP<Epetra_Vector> tangtemptractiontwo_ =
      Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
  normtempinflux_->PutScalar(0.0);
  normtemptraction_->PutScalar(0.0);
  tangtemptractionone_->PutScalar(0.0);
  tangtemptractiontwo_->PutScalar(0.0);

  while (step_fsi < nstep_fsi and t + 1e-10 * dt_fsi < maxtime_fsi)
  {
    step_fsi++;
    t += dt_fsi;

    fsi_->PrepareTimeStep();
    fsi_->TimeStep(fsi_);
    fsi_->PrepareOutput();
    fsi_->Update();

    SetFSISolution();

    if (Comm().MyPID() == 0)
    {
      std::cout << "\n***********************\n GAS TRANSPORT SOLVER \n***********************\n";
    }

    // first scatra field is associated with fluid, second scatra field is
    // associated with structure

    bool stopnonliniter = false;
    int itnum = 0;

    PrepareTimeStep();

    while (stopnonliniter == false)
    {
      ScatraEvaluateSolveIterUpdate();
      itnum++;
      if (ScatraConvergenceCheck(itnum)) break;
    }

    // calculation of the flux at the interface based on normal influx values before time shift of
    // results is performed in Update
    Teuchos::RCP<Epetra_MultiVector> strufluxn =
        scatravec_[1]->ScaTraField()->CalcFluxAtBoundary(false);

    UpdateScatraFields();

    // this is necessary because we want to write all the steps except the last one
    // the last one will be written only after the calculation of the growth
    // in this way also the displacement due to growth is written
    if (step_fsi < nstep_fsi and t + 1e-10 * dt_fsi < maxtime_fsi)
    {
      fsi_->Output();
      ScatraOutput();
    }

    // access structure discretization
    Teuchos::RCP<DRT::Discretization> strudis = fsi_->StructureField()->Discretization();

    // recovery of forces at the interface nodes based on lagrange multipliers values
    // lambda_ is defined only at the interface, while lambdafull on the entire fluid/structure
    // field.
    Teuchos::RCP<Epetra_Vector> lambda_;
    Teuchos::RCP<Epetra_Vector> lambdafull;

    // at the purpose to compute lambdafull, it is necessary to know which coupling algorithm is
    // used however the imposition of a Dirichlet condition on the interface produce wrong lambda_
    // when structuresplit is used
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = Teuchos::getIntegralValue<int>(fsidyn, "COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit)
    {
      Teuchos::RCP<Epetra_Vector> lambdafluid = fsi_->GetLambda();
      lambda_ = fsi_->FluidToStruct(lambdafluid);
    }
    else if (coupling == fsi_iter_monolithicstructuresplit)
    {
      lambda_ = fsi_->GetLambda();
    }

    lambdafull = fsi_->StructureField()->Interface()->InsertFSICondVector(lambda_);

    // calculate interface normals in deformed configuration
    Teuchos::RCP<Epetra_Vector> nodalnormals =
        Teuchos::rcp(new Epetra_Vector(*(strudis->DofRowMap())));

    Teuchos::ParameterList eleparams;
    eleparams.set("action", "calc_cur_nodal_normals");
    strudis->ClearState();
    strudis->SetState("displacement", fsi_->StructureField()->Dispnp());
    strudis->EvaluateCondition(eleparams, Teuchos::null, Teuchos::null, nodalnormals, Teuchos::null,
        Teuchos::null, "FSICoupling");
    strudis->ClearState();

    const Epetra_Map* dofrowmap = strudis->DofRowMap();
    const Epetra_Map* noderowmap = strudis->NodeRowMap();
    Teuchos::RCP<Epetra_MultiVector> lambdanode =
        Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));

    // lagrange multipliers defined on a nodemap are necessary
    for (int lnodeid = 0; lnodeid < strudis->NumMyRowNodes(); lnodeid++)
    {
      // get the processor local node
      DRT::Node* lnode = strudis->lRowNode(lnodeid);
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
    Teuchos::RCP<Epetra_Map> condnodemap = DRT::UTILS::ConditionNodeRowMap(*strudis, "FSICoupling");
    for (int nodei = 0; nodei < condnodemap->NumMyElements(); nodei++)
    {
      // Here we rely on the fact that the structure scatra discretization is a clone of the
      // structure mesh

      // get the processor's local node with the same lnodeid
      int gnodeid = condnodemap->GID(nodei);
      DRT::Node* strulnode = strudis->gNode(gnodeid);
      // get the degrees of freedom associated with this node
      std::vector<int> strunodedofs = strudis->Dof(0, strulnode);
      // determine number of space dimensions
      const int numdim = ((int)strunodedofs.size());

      std::vector<int> doflids(numdim);
      double temp = 0.;
      std::vector<double> unitnormal(3);
      for (int i = 0; i < numdim; ++i)
      {
        doflids[i] = strudis->DofRowMap()->LID(strunodedofs[i]);
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
    Teuchos::RCP<DRT::Discretization> strudis = fsi_->StructureField()->Discretization();

    // loop over all local interface nodes of structure discretization
    Teuchos::RCP<Epetra_Map> condnodemap = DRT::UTILS::ConditionNodeRowMap(*strudis, "FSICoupling");
    for (int i = 0; i < condnodemap->NumMyElements(); i++)
    {
      // get the processor's local node with the same lnodeid
      int gnodeid = condnodemap->GID(i);
      int lnodeid = strudis->NodeRowMap()->LID(gnodeid);

      (*((*norminflux_)(0)))[lnodeid] = (*((*normtempinflux_)(0)))[lnodeid] / step_fsi;
      (*((*normtraction_)(0)))[lnodeid] = (*((*normtemptraction_)(0)))[lnodeid] / step_fsi;
      (*((*tangtractionone_)(0)))[lnodeid] = (*((*tangtemptractionone_)(0)))[lnodeid] / step_fsi;
      (*((*tangtractiontwo_)(0)))[lnodeid] = (*((*tangtemptractiontwo_)(0)))[lnodeid] / step_fsi;
    }
  }

  time_fsi += t;
}

/*----------------------------------------------------------------------*
 | write FSI solutions into scatra discretisation             Thon 11/14|
 *----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::SetFSISolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::ComputeInterfaceVectors(Teuchos::RCP<Epetra_Vector> idispnp,
    Teuchos::RCP<Epetra_Vector> iveln, Teuchos::RCP<Epetra_Vector> struidispnp,
    Teuchos::RCP<Epetra_Vector> struiveln)
{
  // initialize structure interface displacement at time t^{n+1}
  // shouldn't that be zeroed?
  struidispnp->PutScalar(0.0);

  // select biofilm growth boundaries
  std::string biogrcondname = "BioGrCoupling";

  // set action for elements: compute normal vectors at nodes (for reference configuration)
  Teuchos::RCP<DRT::Discretization> strudis = fsi_->StructureField()->Discretization();
  Teuchos::RCP<Epetra_Vector> nodalnormals =
      Teuchos::rcp(new Epetra_Vector(*(strudis->DofRowMap())));
  Teuchos::ParameterList eleparams;
  eleparams.set("action", "calc_ref_nodal_normals");
  strudis->EvaluateCondition(eleparams, Teuchos::null, Teuchos::null, nodalnormals, Teuchos::null,
      Teuchos::null, biogrcondname);

  // select row map with nodes from condition
  Teuchos::RCP<Epetra_Map> condnodemap = DRT::UTILS::ConditionNodeRowMap(*strudis, biogrcondname);

  // loop all conditioned nodes
  for (int i = 0; i < condnodemap->NumMyElements(); ++i)
  {
    int nodegid = condnodemap->GID(i);
    if (strudis->HaveGlobalNode(nodegid) == false) dserror("node not found on this proc");
    DRT::Node* actnode = strudis->gNode(nodegid);
    std::vector<int> globaldofs = strudis->Dof(0, actnode);
    const int numdim = (int)(globaldofs.size());

    // extract averaged nodal normal and compute its absolute value
    std::vector<double> unitnormal(numdim);
    double temp = 0.;
    for (int j = 0; j < numdim; ++j)
    {
      unitnormal[j] = (*nodalnormals)[strudis->DofRowMap()->LID(globaldofs[j])];
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

    int error = struiveln_->ReplaceGlobalValues(numdim, &Values[0], &globaldofs[0]);
    if (error > 0) dserror("Could not insert values into vector struiveln_: error %d", error);
  }

  struidispnp->Update(dt_bio, *struiveln_, 0.0);

  Teuchos::RCP<Epetra_Vector> fluididisp = fsi_->StructToFluid(struidispnp);
  idispnp->Update(1.0, *fluididisp, 0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::FluidAleSolve()
{
  Teuchos::RCP<DRT::Discretization> fluidaledis = fsi_->AleField()->WriteAccessDiscretization();

  // if we have values at the fluid interface we need to apply them
  if (idispnp_ != Teuchos::null)
  {
    fsi_->AleField()->ApplyInterfaceDisplacements(FluidToAle(idispnp_));
  }

  fsi_->AleField()->CreateSystemMatrix(Teuchos::null);
  fsi_->AleField()->Evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_biofilm);
  int error = fsi_->AleField()->Solve();
  if (error == 1) dserror("Could not solve fluid ALE in biofilm FS3I!");
  fsi_->AleField()->UpdateIter();

  // change nodes reference position of the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(fsi_->AleField()->WriteAccessDispnp());
  Teuchos::RCP<DRT::Discretization> fluiddis = fsi_->FluidField()->Discretization();
  DRT::UTILS::UpdateMaterialConfigWithDispVector(fluiddis, fluiddisp);



  // change nodes reference position also for the fluid ale field
  Teuchos::RCP<Epetra_Vector> fluidaledisp = fsi_->AleField()->WriteAccessDispnp();
  DRT::UTILS::UpdateMaterialConfigWithDispVector(fluidaledis, fluidaledisp);

  // change nodes reference position also for scatra fluid field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[0];
  Teuchos::RCP<DRT::Discretization> scatradis = scatra->ScaTraField()->Discretization();
  FS3I::BIOFILM::UTILS::ScatraChangeConfig(scatradis, fluiddis, fluiddisp);

  // set the total displacement due to growth for output reasons
  // fluid
  fluid_growth_disp->Update(1.0, *fluiddisp, 1.0);
  fsi_->FluidField()->SetFldGrDisp(fluid_growth_disp);
  // fluid scatra
  VecToScatravec(scatradis, fluid_growth_disp, scatra_fluid_growth_disp);
  scatra->ScaTraField()->SetScFldGrDisp(scatra_fluid_growth_disp);

  // computation of fluid solution
  // fluid_->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::StructAleSolve()
{
  Teuchos::RCP<DRT::Discretization> structaledis = ale_->WriteAccessDiscretization();

  // if we have values at the structure interface we need to apply them
  if (struidispnp_ != Teuchos::null)
  {
    ale_->ApplyInterfaceDisplacements(StructToAle(struidispnp_));
  }

  ale_->CreateSystemMatrix(Teuchos::null);
  ale_->Evaluate(Teuchos::null, ALE::UTILS::MapExtractor::dbc_set_biofilm);
  int error = ale_->Solve();
  if (error == 1) dserror("Could not solve fluid ALE in biofilm FS3I!");
  ale_->UpdateIter();

  // change nodes reference position of the structure field
  Teuchos::RCP<Epetra_Vector> structdisp = AleToStructField(ale_->WriteAccessDispnp());
  Teuchos::RCP<DRT::Discretization> structdis = fsi_->StructureField()->Discretization();
  DRT::UTILS::UpdateMaterialConfigWithDispVector(structdis, structdisp);
  structdis->FillComplete(false, true, true);

  // change nodes reference position also for the struct ale field
  DRT::UTILS::UpdateMaterialConfigWithDispVector(structaledis, ale_->WriteAccessDispnp());

  // change nodes reference position also for scatra structure field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> struscatra = scatravec_[1];
  Teuchos::RCP<DRT::Discretization> struscatradis = struscatra->ScaTraField()->Discretization();
  FS3I::BIOFILM::UTILS::ScatraChangeConfig(struscatradis, structdis, structdisp);

  // set the total displacement due to growth for output reasons
  // structure
  struct_growth_disp->Update(1.0, *structdisp, 1.0);
  fsi_->StructureField()->SetStrGrDisp(struct_growth_disp);
  // structure scatra
  VecToScatravec(struscatradis, struct_growth_disp, scatra_struct_growth_disp);
  struscatra->ScaTraField()->SetScStrGrDisp(scatra_struct_growth_disp);

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
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
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
void FS3I::BiofilmFSI::VecToScatravec(Teuchos::RCP<DRT::Discretization> scatradis,
    Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<Epetra_MultiVector> scatravec)
{
  // define error variable
  int err(0);

  // loop over all local nodes of scatra discretization
  for (int lnodeid = 0; lnodeid < scatradis->NumMyRowNodes(); lnodeid++)
  {
    // determine number of space dimensions
    const int numdim = DRT::Problem::Instance()->NDim();

    for (int index = 0; index < numdim; ++index)
    {
      double vecval = (*vec)[index + numdim * lnodeid];

      // insert value into node-based vector
      err = scatravec->ReplaceMyValue(lnodeid, index, vecval);

      if (err != 0) dserror("Error while inserting value into vector scatravec!");
    }

    // for 1- and 2-D problems: set all unused vector components to zero
    for (int index = numdim; index < 3; ++index)
    {
      err = scatravec->ReplaceMyValue(lnodeid, index, 0.0);
      if (err != 0) dserror("Error while inserting value into vector scatravec!");
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::StructGmshOutput()
{
  const Teuchos::RCP<DRT::Discretization> structdis = fsi_->StructureField()->Discretization();
  const Teuchos::RCP<DRT::Discretization> structaledis = ale_->WriteAccessDiscretization();
  Teuchos::RCP<DRT::Discretization> struscatradis = scatravec_[1]->ScaTraField()->Discretization();

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(
      "struct", step_bio, 701, false, structdis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  Teuchos::RCP<const Epetra_Vector> structdisp = fsi_->StructureField()->Dispn();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "struct displacement \" {" << std::endl;
    // draw vector field 'struct displacement' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(structdis, structdisp, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<const Epetra_Vector> structaledisp = ale_->Dispnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "struct ale displacement \" {" << std::endl;
    // draw vector field 'struct ale displacement' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(structaledis, structaledisp, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<const Epetra_Vector> structphi = scatravec_[1]->ScaTraField()->Phinp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "struct phi \" {" << std::endl;
    // draw vector field 'struct phi' for every element
    IO::GMSH::ScalarFieldToGmsh(struscatradis, structphi, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::FluidGmshOutput()
{
  const Teuchos::RCP<DRT::Discretization> fluiddis = fsi_->FluidField()->Discretization();
  const Teuchos::RCP<DRT::Discretization> fluidaledis =
      fsi_->AleField()->WriteAccessDiscretization();
  Teuchos::RCP<DRT::Discretization> fluidscatradis = scatravec_[0]->ScaTraField()->Discretization();

  const std::string filenamefluid = IO::GMSH::GetNewFileNameAndDeleteOldFiles(
      "fluid", step_bio, 701, false, fluiddis->Comm().MyPID());
  std::ofstream gmshfilecontent(filenamefluid.c_str());

  Teuchos::RCP<const Epetra_Vector> fluidvel = fsi_->FluidField()->Velnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "fluid velocity \" {" << std::endl;
    // draw vector field 'fluid velocity' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(fluiddis, fluidvel, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<Epetra_Vector> fluidaledisp = fsi_->AleField()->WriteAccessDispnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "fluid ale displacement \" {" << std::endl;
    // draw vector field 'fluid ale displacement' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(fluidaledis, fluidaledisp, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  Teuchos::RCP<Epetra_Vector> fluidphi = scatravec_[0]->ScaTraField()->Phinp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "fluid phi \" {" << std::endl;
    // draw vector field 'fluid phi' for every element
    IO::GMSH::ScalarFieldToGmsh(fluidscatradis, fluidphi, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();

  return;
}
