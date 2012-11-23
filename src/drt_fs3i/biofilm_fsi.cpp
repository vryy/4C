#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fsi/fsi_dyn.H"

#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "../drt_fsi/fs_monolithic.H"

#include "../drt_fsi/fsi_nox_aitken.H"
#include "../drt_fsi/fsi_nox_group.H"
#include "../drt_fsi/fsi_nox_newton.H"
#include "../drt_fsi/fsi_statustest.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_scatra/scatra_utils.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_ale/ale_utils_clonestrategy.H"

#include "../drt_structure/stru_aux.H"

#include <mpi.h>
#include <Epetra_SerialComm.h>

#include "biofilm_fsi.H"
#include "../drt_adapter/ad_str_bio.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_lib/drt_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_adapter/ad_fld_wrapper.H"

//#define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::BiofilmFSI::BiofilmFSI(const Epetra_Comm& comm)
  :PartFS3I_1WC(comm),
   comm_(comm)
{

  DRT::Problem* problem = DRT::Problem::Instance();
  problem->GetDis("structale")->FillComplete();

  //---------------------------------------------------------------------
  // create ale elements if not yet existing
  //---------------------------------------------------------------------

  RCP<DRT::Discretization> structaledis = problem->GetDis("structale");
  if (structaledis->NumGlobalNodes()==0)
  {
    RCP<DRT::Discretization> structdis = problem->GetDis("structure");
    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<ALE::UTILS::AleCloneStrategy> > alecreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<ALE::UTILS::AleCloneStrategy>() );

    alecreator->CreateMatchingDiscretization(structdis,structaledis,11);
  }

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(fsidyn,1));

  ale_ = ale->AleFieldrcp();

  //---------------------------------------------------------------------
  // set up couplings
  //---------------------------------------------------------------------

  const string condname = "FSICoupling";

  // set up ale-fluid couplings
	const int  ndim = DRT::Problem::Instance()->NDim();
  icoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupfa_->SetupConditionCoupling(*(fsi_->FluidField().Discretization()),
                                   (fsi_->FluidField().Interface()->FSICondMap()),
                                   *(fsi_->AleField().Discretization()),
                                   (fsi_->AleField().Interface()->FSICondMap()),
                                   condname,
                                   ndim);

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fsi_->FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* fluidalenodemap   = fsi_->AleField().Discretization()->NodeRowMap();

  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_->SetupCoupling(*(fsi_->FluidField().Discretization()),
                         *(fsi_->AleField().Discretization()),
                         *fluidnodemap,
                         *fluidalenodemap,
                         ndim);


  // set up ale-structure couplings
  icoupsa_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupsa_->SetupConditionCoupling(*(fsi_->StructureField()->Discretization()),
                                   fsi_->StructureField()->Interface()->FSICondMap(),
                                   *structaledis,
                                   ale_->Interface()->FSICondMap(),
                                   condname,
                                   ndim);

  // the structure-ale coupling always matches
  const Epetra_Map* structurenodemap = fsi_->StructureField()->Discretization()->NodeRowMap();
  const Epetra_Map* structalenodemap   = structaledis->NodeRowMap();

  coupsa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsa_->SetupCoupling(*(fsi_->StructureField()->Discretization()),
                         *structaledis,
                         *structurenodemap,
                         *structalenodemap,
                         ndim);

  /// do we need this. What's for???
  fsi_->FluidField().SetMeshMap(coupfa_->MasterDofMap());

  //---------------------------------------------------------------------
  // getting and initializing problem-specific parameters
  //---------------------------------------------------------------------

  const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();

  // make sure that initial time derivative of concentration is not calculated
  // automatically (i.e. field-wise)
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  if (DRT::INPUT::IntegralValue<int>(scatradyn,"SKIPINITDER")==false)
    dserror("Initial time derivative of phi must not be calculated automatically -> set SKIPINITDER to false");

  //fsi parameters
  dt_fsi = fsidyn.get<double>("TIMESTEP");
  nstep_fsi = fsidyn.get<int>("NUMSTEP");
  maxtime_fsi = fsidyn.get<double>("MAXTIME");
  step_fsi = 0;
  time_fsi = 0.;

  //surface growth parameters
  dt_bio= biofilmcontrol.get<double>("BIOTIMESTEP");
  nstep_bio= biofilmcontrol.get<int>("BIONUMSTEP");
  grownvolume_ = biofilmcontrol.get<double>("GROWNVOLUME");
  step_bio=0;
  time_bio = 0.;

  //total time
  time_ = 0.;

  idispn_= fsi_->FluidField().ExtractInterfaceVeln();
  idispnp_= fsi_->FluidField().ExtractInterfaceVeln();
  iveln_= fsi_->FluidField().ExtractInterfaceVeln();

  struidispn_= fsi_->StructureField()->ExtractInterfaceDispn();
  struidispnp_= fsi_->StructureField()->ExtractInterfaceDispn();
  struiveln_= fsi_->StructureField()->ExtractInterfaceDispn();

  struct_growth_disp= AleToStructField(ale_->ExtractDispnp());
  fluid_growth_disp= AleToFluidField(fsi_->AleField().ExtractDispnp());

  idispn_->PutScalar(0.0);
  idispnp_->PutScalar(0.0);
  iveln_->PutScalar(0.0);

  struidispn_->PutScalar(0.0);
  struidispnp_->PutScalar(0.0);
  struiveln_->PutScalar(0.0);

  struct_growth_disp->PutScalar(0.0);
  fluid_growth_disp->PutScalar(0.0);

  norminflux_ = Teuchos::rcp(new Epetra_Vector(*(fsi_->StructureField()->Discretization()->NodeRowMap())));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::Timeloop()
{

#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();
  const int surfgrowth = DRT::INPUT::IntegralValue<int>(biofilmcontrol,"SURFACEGROWTH");

  if (surfgrowth)
  {
    //outer loop for surface growth
    while (step_bio < nstep_bio)
    {
      //fsi_->SetupSystem();

      // inner loop for fsi and scatra
      InnerTimeloop();

      if (Comm().MyPID()==0)
      {
        cout<<"\n***********************\n     GROWTH STEP \n***********************\n";
        printf(" surface growth step = %3d   \n",step_bio);
        printf(" Total time = %3f   \n",time_);
      }

      // compute interface displacement and velocity
      ComputeInterfaceVectors(idispnp_,iveln_,struidispnp_,struiveln_);

      if (idispnp_!=Teuchos::null)
	    {
	      // if we have values at the interface we need to apply them
	      fsi_->AleField().ApplyInterfaceDisplacements(FluidToAle(idispnp_));
	    }

      // do all the settings and solve the fluid on a deforming mesh
      FluidAleSolve(idispnp_);

      if (struidispnp_!=Teuchos::null)
      {
        // if we have values at the interface we need to apply them
        ale_->ApplyInterfaceDisplacements(StructToAle(struidispnp_));
      }

      // do all the settings and solve the structure on a deforming mesh
      StructAleSolve(struidispnp_);

      step_bio++;
      time_bio+=dt_bio;
      time_ = time_bio + time_fsi;

      fsi_->StructureField()->Reset();
      fsi_->FluidField().Reset(true, true, step_bio);
      fsi_->AleField().Reset();

      fsi_->AleField().BuildSystemMatrix(false);

      StructGmshOutput();
      FluidGmshOutput();
      //fsi_->StructureField()->DiscWriter()->WriteMesh(step_bio, time_bio);
      //fsi_->FluidField().DiscWriter()->WriteMesh(step_bio, time_bio);

      GrowthOutput();
    }
  }

  if (!surfgrowth) InnerTimeloop();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::InnerTimeloop()
{
//  // output of initial state
//  ScatraOutput();

  fsi_->PrepareTimeloop();

  double t=0.;
  step_fsi=0;

  vector<std::string> condnames(1);
  condnames[0] = "FSICoupling";
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[0];
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> struscatra = scatravec_[1];

  flux = scatra->ScaTraField().CalcFluxAtBoundary(condnames,false);
  struflux = struscatra->ScaTraField().CalcFluxAtBoundary(condnames,false);

  flux->PutScalar(0.0);
  struflux->PutScalar(0.0);

  norminflux_->PutScalar(0.0);

  while (step_fsi < nstep_fsi and t+1e-10*dt_fsi < maxtime_fsi)
  {
    step_fsi++;
    t+=dt_fsi;

    DoFSIStep();
    SetFSISolution();
    DoScatraStep();

    // calculation of the flux at the interface based on normal influx values

    Teuchos::RCP<Epetra_MultiVector> strufluxn = struscatra->ScaTraField().CalcFluxAtBoundary(condnames,false);
    Teuchos::RCP<DRT::Discretization> strudis = fsi_->StructureField()->Discretization();

    // calculate interface normals in deformed configuration
    Teuchos::RCP<Epetra_Vector> nodalnormals = Teuchos::rcp(new Epetra_Vector(*(strudis->DofRowMap())));
    std::string condname = "FSICoupling";
    Teuchos::ParameterList eleparams;
    eleparams.set("action","calc_cur_nodal_normals");
    strudis->ClearState();
    strudis->SetState("displacement",fsi_->StructureField()->Dispnp());
    strudis->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,nodalnormals,Teuchos::null,Teuchos::null,condname);
    strudis->ClearState();

    // loop over all local interface nodes of structure discretization
    Teuchos::RCP<Epetra_Map> condnodemap = DRT::UTILS::ConditionNodeRowMap(*strudis, condname);

    for (int i=0; i < condnodemap->NumMyElements(); i++)
    {
      // Here we rely on the fact that the scatra discretization
      // is a clone of the fluid mesh. => a scatra node has the same
      // local (and global) ID as its corresponding fluid node!

      // get the processor's local node with the same lnodeid
      int gnodeid = condnodemap->GID(i);
      DRT::Node* strulnode = strudis->gNode(gnodeid);
      // get the degrees of freedom associated with this node
      vector<int> strunodedofs = strudis->Dof(strulnode);

      const int numdim = ((int) strunodedofs.size());
      // number of dof per node in ScaTra

      std::vector<int> doflids(numdim);
      double temp = 0.;
      std::vector<double> unitnormal(3);
      for (int i=0; i<numdim; ++i)
      {
        doflids[i] = strudis->DofRowMap()->LID(strunodedofs[i]);
        unitnormal[i] = (*nodalnormals)[doflids[i]];
        temp += unitnormal[i]*unitnormal[i];

      }
      double absval = sqrt(temp);
      int lnodeid = strudis->NodeRowMap()->LID(gnodeid);

      // compute average unit nodal normal
      std::vector<double> Values(numdim);
      for(int j=0; j<numdim; ++j)
      {
        unitnormal[j] /= absval;
      }
      double tempflux = 0.0;
      for(int index=0;index<numdim;++index)
      {
        double fluxcomp = (*((*strufluxn)(index)))[lnodeid];
        tempflux += fluxcomp*unitnormal[index];

      }
      (*((*norminflux_)(0)))[lnodeid] = tempflux;
    }

  }

  time_fsi+=t;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::ComputeInterfaceVectors(RCP<Epetra_Vector> idispnp,
                                               RCP<Epetra_Vector> iveln,
                                               RCP<Epetra_Vector> struidispnp,
                                               RCP<Epetra_Vector> struiveln)
{
	Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[0];
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> struscatra = scatravec_[1];

	// access discretizations
	RCP<DRT::Discretization> fluiddis = fsi_->FluidField().Discretization();
	RCP<DRT::Discretization> scatradis = scatra->ScaTraField().Discretization();
  RCP<DRT::Discretization> strudis = fsi_->StructureField()->Discretization();
  RCP<DRT::Discretization> struscatradis = struscatra->ScaTraField().Discretization();

  // it would be better to introduce a special condition "growth - surface" to include nodes
  // that are not part of the fsi interface

  Teuchos::RCP<Epetra_Vector> nodalnormals = Teuchos::rcp(new Epetra_Vector(*(strudis->DofRowMap())));

  std::string condname = "FSICoupling";

  Teuchos::ParameterList eleparams;

  // set action for elements
  eleparams.set("action","calc_ref_nodal_normals");
  strudis->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,nodalnormals,Teuchos::null,Teuchos::null,condname);

  Teuchos::RCP<Epetra_Map> condnodemap = DRT::UTILS::ConditionNodeRowMap(*strudis, condname);

  // shouldn't that be zeroed?
  struidispnp->PutScalar(0.0);

  // loop all conditioned nodes
  for (int i=0; i<condnodemap->NumMyElements(); ++i)
  {
    int nodegid = condnodemap->GID(i);
    if (strudis->HaveGlobalNode(nodegid)==false) dserror("node not found on this proc");
    DRT::Node* actnode = strudis->gNode(nodegid);
    std::vector<int> globaldofs = strudis->Dof(actnode);
    const int numdim = (int)(globaldofs.size());

    // extract averaged nodal normal and compute its absolute value
    std::vector<double> unitnormal(numdim);
    double temp = 0.;
    for (int j=0; j<numdim; ++j)
    {
      unitnormal[j] = (*nodalnormals)[strudis->DofRowMap()->LID(globaldofs[j])];
      temp += unitnormal[j]*unitnormal[j];

    }
    double absval = sqrt(temp);
    int lnodeid = strudis->NodeRowMap()->LID(nodegid);
    double influx = (*norminflux_)[lnodeid];

    // compute average unit nodal normal and "interface velocity"
    std::vector<double> Values(numdim);
    for(int j=0; j<numdim; ++j)
    {

    //introduced a tolerance otherwise NAN will be produced in case of zero absval
      double TOL = 1e-6;
      if (absval > TOL)
      {
        unitnormal[j] /= absval;
        Values[j] = - grownvolume_ * influx * unitnormal[j];
      }
    }

    int error = struiveln_->ReplaceGlobalValues(numdim,&Values[0],&globaldofs[0]);
    if (error > 0) dserror("Could not insert values into vector struiveln_: error %d",error);
  }

  struidispnp->Update(dt_bio,*struiveln_,0.0);

  Teuchos::RCP<Epetra_Vector> fluididisp = fsi_->StructToFluid(struidispnp);
  idispnp->Update(1.0, *fluididisp, 0.0);

	return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::FluidAleSolve(
    Teuchos::RCP<Epetra_Vector> idisp)
{
  fsi_->AleField().SetupDBCMapEx(1);
  fsi_->AleField().BuildSystemMatrix();
  fsi_->AleField().Solve();

  //change nodes reference position
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(fsi_->AleField().ExtractDispnp());
  RCP<DRT::Discretization> fluiddis = fsi_->FluidField().Discretization();

  fluid_growth_disp->Update(1,*fluiddisp,1.0);

  ChangeConfig(fluiddis, fluiddisp);

  //change nodes reference position also for ale field
  Teuchos::RCP<Epetra_Vector> fluidaledisp = fsi_->AleField().ExtractDispnp();
  RCP<DRT::Discretization> fluidaledis = fsi_->AleField().Discretization();
  ChangeConfig(fluidaledis, fluidaledisp);

  //change nodes reference position also for scatra fluid field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[0];
  RCP<DRT::Discretization> scatradis = scatra->ScaTraField().Discretization();
  ScatraChangeConfig(scatradis, fluiddis, fluiddisp);

  fsi_->AleField().SetupDBCMapEx(0);

  // computation of structure solution
  //structure_->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::StructAleSolve(
    Teuchos::RCP<Epetra_Vector> idisp)
{
  ale_->SetupDBCMapEx(1);
  ale_->BuildSystemMatrix();
  ale_->Solve();

  //change nodes reference position
  Teuchos::RCP<Epetra_Vector> structdisp = AleToStructField(ale_->ExtractDispnp());
  RCP<DRT::Discretization> structdis = fsi_->StructureField()->Discretization();

  struct_growth_disp->Update(1,*structdisp,1.0);

  ChangeConfig(structdis, structdisp);

  structdis->FillComplete(false, true, true);

  //change nodes reference position also for ale field
  RCP<DRT::Discretization> structaledis = ale_->Discretization();
  ChangeConfig(structaledis, ale_->ExtractDispnp());

  //change nodes reference position also for scatra structure field
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> struscatra = scatravec_[1];
  RCP<DRT::Discretization> struscatradis = struscatra->ScaTraField().Discretization();
  ScatraChangeConfig(struscatradis, structdis, structdisp);

  ale_->SetupDBCMapEx(0);

  // computation of structure solution
  //structure_->Solve();

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
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::AleToStructField(Teuchos::RCP<const Epetra_Vector> iv) const
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
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmFSI::StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupsa_->MasterToSlave(iv);
}


// is it needed?
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//void FS3I::BiofilmFSI::UpdateAndOutput()
//{
//  for (unsigned i=0; i<scatravec_.size(); ++i)
//  {
//    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
//    scatra->ScaTraField().Update();
//
//    // perform time shift of interface displacement
//    idispn_->Update(1.0, *idispnp_ , 0.0);
//    struidispn_->Update(1.0, *struidispnp_ , 0.0);
//
//    // in order to do not have artificial velocity
//    idispnp_=idispn_;
//    struidispnp_=struidispn_;
//
//    scatra->ScaTraField().Output();
//  }
//}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::ChangeConfig(RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Vector> disp)
{

  const int numnode = (dis->NodeRowMap())->NumMyElements();
  const Epetra_Vector& gvector =*disp;

  // loop over all nodes
  for (int i = 0; i < numnode; ++i)
  {
    // get current node
    int gid = (dis->NodeRowMap())->GID(i);

    DRT::Node* mynode = dis->gNode(gid);

    vector<int> globaldofs = dis->Dof(mynode);
    std::vector<double> nvector(globaldofs.size());

    // extract local values from the global vector
    //const size_t ldim = globaldofs.size();
    //nvector.resize(ldim);

    for (size_t i=0; i<3; ++i)
    {
      const int lid = gvector.Map().LID(globaldofs[i]);

      if (lid<0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",gvector.Comm().MyPID(),globaldofs[i]);
      nvector[i] += gvector[lid];

    }

    mynode->ChangePos(nvector);

  }

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::ScatraChangeConfig(RCP<DRT::Discretization> scatradis, RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Vector> disp)
{
  const int numnode = (scatradis->NodeRowMap())->NumMyElements();
  const Epetra_Vector& gvector =*disp;

  // loop over all nodes
  for (int index = 0; index < numnode; ++index)
  {
    // get current node
    int gid = (scatradis->NodeRowMap())->GID(index);
    DRT::Node* mynode = scatradis->gNode(gid);

    // get local fluid/structure node with the same local node id
    DRT::Node* lnode = dis->lRowNode(index);

    // get degrees of freedom associated with this fluid/structure node
    // (first dofset always considered, allowing for using multiple
    //  dof sets, e.g., for structure-based scalar transport)
    vector<int> nodedofs = dis->Dof(0,lnode);

    std::vector<double> nvector(nodedofs.size());


    // extract local values from the global vector
//    const size_t ldim = globaldofs.size();
    //nvector.resize(ldim);

    for (size_t i=0; i<3; ++i)
    {

      const int lid = gvector.Map().LID(nodedofs[i]);

      if (lid<0)
      dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",gvector.Comm().MyPID(),nodedofs[i]);
      nvector[i] += gvector[lid];

    }

    mynode->ChangePos(nvector);

  }

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::StructGmshOutput()
{
  // write the structure mesh
  const Teuchos::RCP<DRT::Discretization> structdis = fsi_->StructureField()->Discretization();
  const Teuchos::RCP<DRT::Discretization> structaledis = ale_->Discretization();
  RCP<DRT::Discretization> struscatradis = scatravec_[1]->ScaTraField().Discretization();

  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("struct", step_bio, 701, true, structdis->Comm().MyPID());

  std::ofstream gmshfilecontent(filename.c_str());

  Teuchos::RCP<Epetra_Vector> structdisp = fsi_->StructureField()->ExtractDispn();

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "struct displacement \" {" << endl;
    // draw vector field 'struct displacement' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(structdis,structdisp,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  Teuchos::RCP<Epetra_Vector> strucaletdisp = ale_->ExtractDispnp();

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "struct ale displacement \" {" << endl;
    // draw vector field 'struct ale displacement' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(structaledis,strucaletdisp,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "struct phi \" {" << endl;
    // draw vector field 'struct phi' for every element
    IO::GMSH::ScalarFieldToGmsh(struscatradis,scatravec_[1]->ScaTraField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }
  gmshfilecontent.close();
  std::cout << " done" << endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::FluidGmshOutput()
{

  //write the fluid mesh
  const Teuchos::RCP<DRT::Discretization> fluiddis = fsi_->FluidField().Discretization();
  const Teuchos::RCP<DRT::Discretization> fluidaledis = fsi_->AleField().Discretization();
  RCP<DRT::Discretization> struscatradis = scatravec_[0]->ScaTraField().Discretization();

  const std::string filenamefluid = IO::GMSH::GetNewFileNameAndDeleteOldFiles("fluid", step_bio, 701, true, fluiddis->Comm().MyPID());

  std::ofstream gmshfilecontent(filenamefluid.c_str());

  Teuchos::RCP<const Epetra_Vector> fluidvel = fsi_->FluidField().Velnp();
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "fluid velocity \" {" << endl;
    // draw vector field 'fluid velocity' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(fluiddis,fluidvel,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  Teuchos::RCP<Epetra_Vector> fluidaledisp = fsi_->AleField().ExtractDispnp();

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "fluid ale displacement \" {" << endl;
    // draw vector field 'fluid ale displacement' for every element
    IO::GMSH::VectorFieldDofBasedToGmsh(fluidaledis,fluidaledisp,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "fluid phi \" {" << endl;
    // draw vector field 'fluid phi' for every element
    IO::GMSH::ScalarFieldToGmsh(struscatradis,scatravec_[0]->ScaTraField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();
  std::cout << " done" << endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmFSI::GrowthOutput()
{

//  RCP<DRT::Discretization> actdis;
//
//  actdis = fsi_->FluidField().Discretization();

  fsi_->FluidField().DiscWriter()->WriteVector("growth_displ", fluid_growth_disp);

//  actdis = fsi_->StructureField()->Discretization();

  fsi_->StructureField()->DiscWriter()->WriteVector("growth_displ", struct_growth_disp);

  //introduce output also for scatra
}
