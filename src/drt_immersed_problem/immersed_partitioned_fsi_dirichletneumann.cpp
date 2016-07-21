/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann.cpp

\brief partitioned immersed fsi algorithm for neumann-(dirichlet-neumann) like coupling

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_dirichletneumann.H"
#include "immersed_partitioned_fsi.H"
#include "fsi_partitioned_immersed.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_fld_fluid_immersed.H"
#include "../drt_adapter/ad_fld_fluid_ale_immersed.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_inpar/inpar_immersed.H"

#include "../linalg/linalg_utils.H"

// search tree related
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

// time monitoring
#include <Teuchos_TimeMonitor.hpp>

// relaxation
#include "../drt_fsi/fsi_nox_aitken.H"


IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ImmersedPartitionedFSIDirichletNeumann(const Epetra_Comm& comm)
  : ImmersedBase(),
    FSI::PartitionedImmersed(comm)
{
  // important variables for parallel simulations
  myrank_  = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  immersedstructure_=Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(StructureField());

  if(globalproblem_->FSIDynamicParams().get<std::string>("COUPALGO") == "iter_stagg_fixed_rel_param")
  {
    coupalgo_= INPAR::IMMERSED::fixed;
    if(myrank_==0)
      std::cout<<"\n Using FIXED relaxation parameter. "<<std::endl;
  }
  else if(globalproblem_->FSIDynamicParams().get<std::string>("COUPALGO") == "iter_stagg_AITKEN_rel_param")
  {
    coupalgo_ = INPAR::IMMERSED::aitken;
    if(myrank_==0)
      std::cout<<"\n Using AITKEN relaxation parameter. "<<std::endl;
  }
  else
    dserror("Unknown definition of COUPALGO in FSI DYNAMIC section for Immersed FSI.");


  // get integration rule for fluid elements cut by structural boundary
  int num_gp_fluid_bound = globalproblem_->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  if(num_gp_fluid_bound == 8)
    degree_gp_fluid_bound_ = 3;
  else if (num_gp_fluid_bound == 64)
    degree_gp_fluid_bound_ = 7;
  else if (num_gp_fluid_bound == 125)
    degree_gp_fluid_bound_ = 9;
  else if (num_gp_fluid_bound == 343)
    degree_gp_fluid_bound_ = 13;
  else if (num_gp_fluid_bound == 729)
    degree_gp_fluid_bound_ = 17;
  else if (num_gp_fluid_bound == 1000)
    degree_gp_fluid_bound_ = 19;
  else
    dserror("Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, 729 and 1000).");

  // get coupling variable
  displacementcoupling_ = globalproblem_->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE") == "Displacement";
  if(displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned FSI scheme :  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned FSI scheme :  Force "<<std::endl;

  // get pointer to discretizations
  fluiddis_  = globalproblem_->GetDis("fluid");
  structdis_ = globalproblem_->GetDis("structure");

  // decide whether multiple structural bodies or not
  std::vector<DRT::Condition*> conditions;
  structdis_->GetCondition("ImmersedSearchbox",conditions);
  if((int)conditions.size()>0)
  {
    if(myrank_==0)
      std::cout<<" MULTIBODY SIMULATION   Number of bodies: "<<(int)conditions.size()<<std::endl;
    multibodysimulation_ = true;
  }
  else
    multibodysimulation_ = false;

  // vector of fluid stresses interpolated to structural bdry int points and integrated over structural surface
  struct_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(immersedstructure_->DofRowMap()),true));

  // vector with fluid velocities interpolated from structure
  fluid_artificial_velocity_ = Teuchos::rcp(new Epetra_Vector(*(MBFluidField()->FluidField()->DofRowMap()),true));
  //fluid_artificial_velocity_old_ = Teuchos::rcp(new Epetra_Vector(*(MBFluidField()->FluidField()->DofRowMap()),true));


  // build 3D search tree for fluid domain
  fluid_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // find positions of the background fluid discretization
  for (int lid = 0; lid < fluiddis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = fluiddis_->lColNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    currpositions_fluid_[node->Id()] = currpos;
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*fluiddis_,currpositions_fluid_);
  fluid_SearchTree_->initializeTree(rootBox,*fluiddis_,GEO::TreeType(GEO::OCTTREE));

  if(myrank_==0)
    std::cout<<"\n Build Fluid SearchTree ... "<<std::endl;


  // construct 3D search tree for structural domain
  // initialized in SetupStructuralDiscretization()
  structure_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  correct_boundary_velocities_= (DRT::INPUT::IntegralValue<int>(globalproblem_->ImmersedMethodParams(), "CORRECT_BOUNDARY_VELOCITIES"));

  output_evry_nlniter_= (DRT::INPUT::IntegralValue<int>(globalproblem_->ImmersedMethodParams(), "OUTPUT_EVRY_NLNITER"));

  comm.Barrier();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  // reinitialize the transfer vectors
  //fluid_artificial_velocity_old_->Update(1.0,*fluid_artificial_velocity_,0.0);
  fluid_artificial_velocity_->PutScalar(0.0);
  struct_bdry_traction_->PutScalar(0.0);

  // reset element and node information about immersed method
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::reset_immersed_ele);
  params.set<int>("intpoints_fluid_bound",degree_gp_fluid_bound_);
  EvaluateSubsetElements(params,
                         fluiddis_,
                         curr_subset_of_fluiddis_,
                         (int)FLD::reset_immersed_ele);

  if (displacementcoupling_)
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL FluidOp
    ////////////////////
    SetStatesFluidOP();
    PrepareFluidOp();
    FluidOp(fluid_artificial_velocity_, fillFlag);

    ////////////////////
    // CALL StructOp
    ////////////////////
    Teuchos::RCP<Epetra_Vector> idispnp =
    StructOp(immersedstructure_->Interface()->ExtractIMMERSEDCondVector(struct_bdry_traction_), fillFlag);

  int err = CalcResidual(F,idispnp,idispn);

  if(err != 0)
    dserror("Vector update of FSI-residual returned err=%d",err);

  }
  else if(!displacementcoupling_) // FORCE COUPLING
  {
    // get the initial guess
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL StructOp
    ////////////////////
    StructOp(iforcen, fillFlag);

    ////////////////////
    // CALL FluidOp
    ////////////////////
    SetStatesFluidOP();
    PrepareFluidOp();
    FluidOp(fluid_artificial_velocity_, fillFlag);

    int err = CalcResidual(F,struct_bdry_traction_,iforcen);

    if(err != 0)
      dserror("Vector update of FSI-residual returned err=%d",err);

  }

  // write output after every solve of fluid and structure
  // current limitations:
  // max 100 partitioned iterations and max 100 timesteps in total
  if(output_evry_nlniter_)
  {
    int iter = ((FSI::Partitioned::IterationCounter())[0]);
    Teuchos::rcp_dynamic_cast<ADAPTER::FluidAleImmersed>(MBFluidField())->Output((Step()*100)+(iter-1),Time()-Dt()*((100-iter)/100.0));
    StructureField()->PrepareOutput();
    Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(StructureField())->Output(false,(Step()*100)+(iter-1),Time()-Dt()*((100-iter)/100.0));
  }

  // perform n steps max; then set converged
  bool nlnsolver_continue = globalproblem_->ImmersedMethodParams().get<std::string>("DIVERCONT") == "continue";
  int  itemax = globalproblem_->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<int>("ITEMAX");
  if((FSI::Partitioned::IterationCounter())[0] == itemax and nlnsolver_continue)
  {
    if(myrank_==0)
      std::cout<<"\n  Continue with next time step after ITEMAX = "<<(FSI::Partitioned::IterationCounter())[0]<<" iterations. \n"<<std::endl;

    // !!! EXPERIMENTAL !!!
    // set F to zero to tell NOX that this timestep is converged
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp( new Epetra_Vector(F.Map(),true));
    F.Update(1.0,*zeros,0.0);
    // !!! EXPERIMENTAL !!!

    // clear states after time step was set converged
    immersedstructure_->Discretization()->ClearState();
    MBFluidField()->Discretization()->ClearState();
  }

  if(globalproblem_->ImmersedMethodParams().get<std::string>("TIMESTATS")=="everyiter")
  {
    Teuchos::TimeMonitor::summarize();
    Teuchos::TimeMonitor::zeroOutTimers();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FluidOp(Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity,
                               const FillType fillFlag)
{
  // print
  FSI::Partitioned::FluidOp(fluid_artificial_velocity,fillFlag);
  Teuchos::ParameterList params;
  params.set<int>("intpoints_fluid_bound",degree_gp_fluid_bound_);

  if (fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    // normal fluid solve

    const int itemax = MBFluidField()->Itemax();

    // calc the fluid velocity from the structural displacements
    DRT::AssembleStrategy fluid_vol_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        fluid_artificial_velocity,  // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );

    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate Dirichlet Values from immersed elements which overlap the "<<MBFluidField()->Discretization()->Name()<<" nodes ..."<<std::endl;
    }

      EvaluateImmersed(params,
                                        MBFluidField()->Discretization(),
                                        &fluid_vol_strategy,
                                        &curr_subset_of_fluiddis_,
                                        structure_SearchTree_,
                                        &currpositions_struct_,
                                        (int)FLD::interpolate_velocity_to_given_point_immersed,
                                        false);

    BuildImmersedDirichMap(MBFluidField()->Discretization(), dbcmap_immersed_, MBFluidField()->FluidField()->GetDBCMapExtractor()->CondMap());
    AddDirichCond();

    // apply immersed dirichlets
    DoImmersedDirichletCond(MBFluidField()->FluidField()->WriteAccessVelnp(),fluid_artificial_velocity, dbcmap_immersed_);
    double normofvelocities;
    MBFluidField()->FluidField()->ExtractVelocityPart(fluid_artificial_velocity)->Norm2(&normofvelocities);


    if(myrank_ == 0)
    {
      std::cout<<"###   Norm of Dirichlet Values:   "<<std::setprecision(7)<<normofvelocities<<std::endl;
      std::cout<<"###   Norm of transferred divergence of structural velocity:   "<<std::setprecision(10)<<"not used"<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }

    // solve fluid
    SolveFluid();

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration)
    RemoveDirichCond();

    //****************************************************************************
    // Correct velocity in nodes of fluid elements cut by the structural surface
    // (fluid is solved second time with different dirichlet values)
    //****************************************************************************
    if (correct_boundary_velocities_)
    {

      SetStatesVelocityCorrection();

      if(myrank_ == 0)
      {
        std::cout<<"\nCorrection step "<<std::endl;
        std::cout<<"################################################################################################"<<std::endl;
        std::cout<<"###   Correct Velocity in fluid boundary elements "<<std::endl;
      }

      // calculate new dirichlet velocities for fluid elements cut by structure
      EvaluateImmersed(params,
                                        MBFluidField()->Discretization(),
                                        &fluid_vol_strategy,
                                        &curr_subset_of_fluiddis_,
                                        structure_SearchTree_,
                                        &currpositions_struct_,
                                        (int)FLD::correct_immersed_fluid_bound_vel,
                                        true);

      // Build new dirich map
      BuildImmersedDirichMap(MBFluidField()->Discretization(), dbcmap_immersed_, MBFluidField()->FluidField()->GetDBCMapExtractor()->CondMap());
      AddDirichCond();

      // apply new dirichlets after velocity correction
      DoImmersedDirichletCond(MBFluidField()->FluidField()->WriteAccessVelnp(),fluid_artificial_velocity_, dbcmap_immersed_);
      double normofnewvelocities;
      MBFluidField()->FluidField()->ExtractVelocityPart(fluid_artificial_velocity)->Norm2(&normofnewvelocities);

      if(myrank_ == 0)
      {
        std::cout<<"###   Norm of new Dirichlet Values:   "<<std::setprecision(7)<<normofnewvelocities<<std::endl;
        std::cout<<"################################################################################################"<<std::endl;
      }

      // solve fluid again with new dirichlet values
      SolveFluid();

      // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration step)
      RemoveDirichCond();
    } // correct_boundary_velocities_ finished

    // set max number of Newton iterations
    MBFluidField()->SetItemax(itemax);

    // calculate new fluid tractions interpolated to structural surface
    CalcFluidTractionsOnStructure();

  }

  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::StructOp(Teuchos::RCP<Epetra_Vector> struct_bdry_traction,
                                                           const FillType fillFlag)
{
  FSI::Partitioned::StructOp(struct_bdry_traction,fillFlag);

  if(fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
    return Teuchos::null;
  }
  else
  {
    // prescribe neumann values at structural boundary dofs
    ApplyInterfaceForces(struct_bdry_traction);

    // solve
    SolveStruct();

    return ExtractInterfaceDispnp();
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::InitialGuess()
{

  if(displacementcoupling_)
    return immersedstructure_->PredictImmersedInterfaceDispnp();
  else
  {
    return immersedstructure_->Interface()->ExtractIMMERSEDCondVector(struct_bdry_traction_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& dirichmap,const Teuchos::RCP<const Epetra_Map>& dirichmap_original)
{
  const Epetra_Map* elecolmap = dis->ElementColMap();
  std::vector<int> mydirichdofs(0);

  for(int i=0; i<elecolmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(dis->gElement(elecolmap->GID(i)));
    if(immersedele->HasProjectedDirichlet())
    {
      DRT::Node** nodes = immersedele->Nodes();
      for (int inode=0; inode<(immersedele->NumNode()); inode++)
      {
        if(static_cast<IMMERSED::ImmersedNode* >(nodes[inode])->IsMatched() and nodes[inode]->Owner()==myrank_)
        {
          std::vector<int> dofs = dis->Dof(nodes[inode]);

          for (int dim=0;dim<3;++dim)
          {
            if(dirichmap_original->LID(dofs[dim]) == -1) // if not already in original dirich map
              mydirichdofs.push_back(dofs[dim]);
          }

          // include also pressure dof if node does not belong to a boundary background element
          //if((nodes[inode]->IsBoundaryImmersed())==0)
          //  mydirichdofs.push_back(dofs[3]);
        }
      }
    }
  }

  int nummydirichvals = mydirichdofs.size();
  dirichmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap)
{
  int mynumvals = dbcmap->NumMyElements();
  double* myvals = dirichvals->Values();

  for(int i=0;i<mynumvals;++i)
  {
    int gid = dbcmap->GID(i);

#ifdef DEBUG
    int err = -2;
    int lid = dirichvals->Map().LID(gid);
    err = statevector -> ReplaceGlobalValue(gid,0,myvals[lid]);
    if(err==-1)
      dserror("VectorIndex >= NumVectors()");
    else if (err==1)
        dserror("GlobalRow not associated with calling processor");
    else if (err != -1 and err != 1 and err != 0)
      dserror("Trouble using ReplaceGlobalValue on fluid state vector. ErrorCode = %d",err);
#else
    int lid = dirichvals->Map().LID(gid);
    statevector -> ReplaceGlobalValue(gid,0,myvals[lid]);
#endif

  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetupStructuralDiscretization()
{
  // ghost structure on each proc (for search algorithm)
  if(numproc_ > 1)
  {
    // fill complete inside
    CreateGhosting(structdis_);
  }
  else
  {
    // fill complete to incorporate changes due to ghosting and build geometries
    structdis_->FillComplete();
  }

  // find positions of the immersed structural discretization
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_struct;
  for (int lid = 0; lid < structdis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = structdis_->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    my_currpositions_struct[node->Id()] = currpos;
  }
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[numproc_];
  for(int i=0;i<numproc_;i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_struct,currpositions_struct_,numproc_,&procs[0],Comm());

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox2 = GEO::getXAABBofDis(*structdis_,currpositions_struct_);
  structure_SearchTree_->initializeTree(rootBox2,*structdis_,GEO::TreeType(GEO::OCTTREE));

  if(myrank_==0)
    std::cout<<"\n Build Structure SearchTree ... "<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetStatesFluidOP()
{
  // for FluidOP
  structdis_->SetState(0,"displacement",immersedstructure_->Dispnp());
  structdis_->SetState(0,"velocity",immersedstructure_->Velnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetStatesVelocityCorrection()
{
  fluiddis_->SetState(0,"velnp",MBFluidField()->FluidField()->Velnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetStatesStructOP()
{
  fluiddis_->SetState(0,"velnp",Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->FluidField()->Velnp());
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SolveFluid()
{
  MBFluidField()->NonlinearSolve(Teuchos::null,Teuchos::null);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SolveStruct()
{
  immersedstructure_->Solve();

return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::PrepareFluidOp()
{
  TEUCHOS_FUNC_TIME_MONITOR("IMMERSED::PrepareFluidOp()");

  double structsearchradiusfac = globalproblem_->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

// DEBUG option:
// mark elements in which structural boundary is immersed
#ifdef DEBUG
//  Teuchos::ParameterList params;
//  params.set<std::string>("action","mark_immersed_elements");
//
//  DRT::AssembleStrategy dummy_strategy(
//      0,              // struct dofset for row
//      0,              // struct dofset for column
//      Teuchos::null,  // matrix 1
//      Teuchos::null,  //
//      Teuchos::null,  // vector 1
//      Teuchos::null,  //
//      Teuchos::null   //
//  );
//
//  EvaluateInterpolationCondition( immersedstructure_->Discretization(), params, dummy_strategy,"IMMERSEDCoupling", -1 );
# endif

  //
  // determine subset of fluid discretization which is potentially underlying the structural discretization
  //
  // get state
  Teuchos::RCP<const Epetra_Vector> displacements = immersedstructure_->Dispnp();
  //std::cout<<*displacements<<std::endl;
  // find current positions for immersed structural discretization
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_struct;
  for (int lid = 0; lid < structdis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = structdis_->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;
    std::vector<int> dofstoextract(3);
    std::vector<double> mydisp(3);

    // get the current displacement
    structdis_->Dof(node,0,dofstoextract);
    DRT::UTILS::ExtractMyValues(*displacements,mydisp,dofstoextract);

    currpos(0) = node->X()[0]+mydisp.at(0);
    currpos(1) = node->X()[1]+mydisp.at(1);
    currpos(2) = node->X()[2]+mydisp.at(2);

    my_currpositions_struct[node->Id()] = currpos;
  }
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[numproc_];
  for(int i=0;i<numproc_;i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_struct,currpositions_struct_,numproc_,&procs[0],Comm());

//  DEBUG output
//  std::cout<<"Proc "<<myrank_<<": my_curr_pos.size()="<<my_currpositions_struct.size()<<std::endl;
//  std::cout<<"Proc "<<myrank_<<":    curr_pos.size()="<<currpositions_struct_.size()<<std::endl;

  if (multibodysimulation_ == false)
  {
    // get bounding box of current configuration of structural dis
    const LINALG::Matrix<3,2> structBox = GEO::getXAABBofDis(*structdis_,currpositions_struct_);
    double max_radius = sqrt(pow(structBox(0,0)-structBox(0,1),2)+pow(structBox(1,0)-structBox(1,1),2)+pow(structBox(2,0)-structBox(2,1),2));
    // search for background elements within a certain radius around the center of the immersed bounding box
    LINALG::Matrix<3,1> boundingboxcenter;
    boundingboxcenter(0) = structBox(0,0)+(structBox(0,1)-structBox(0,0))*0.5;
    boundingboxcenter(1) = structBox(1,0)+(structBox(1,1)-structBox(1,0))*0.5;
    boundingboxcenter(2) = structBox(2,0)+(structBox(2,1)-structBox(2,0))*0.5;

# ifdef DEBUG
    std::cout<<"Bounding Box of Structure: "<<structBox<<" on PROC "<<myrank_<<std::endl;
    std::cout<<"Bounding Box Center of Structure: "<<boundingboxcenter<<" on PROC "<<myrank_<<std::endl;
    std::cout<<"Search Radius Around Center: "<<structsearchradiusfac*max_radius<<" on PROC "<<myrank_<<std::endl;
    std::cout<<"Length of Dispnp()="<<displacements->MyLength()<<" on PROC "<<myrank_<<std::endl;
    std::cout<<"Size of currpositions_struct_="<<currpositions_struct_.size()<<" on PROC "<<myrank_<<std::endl;
    std::cout<<"My DofRowMap Size="<<structdis_->DofRowMap()->NumMyElements()<<"  My DofColMap Size="<<structdis_->DofColMap()->NumMyElements()<<" on PROC "<<myrank_<<std::endl;
    std::cout<<"Dis Structure NumColEles: "<<structdis_->NumMyColElements()<<" on PROC "<<myrank_<<std::endl;
# endif

    UpdateCurrentPositionsFluidNodes();
    SearchPotentiallyCoveredBackgrdElements(&curr_subset_of_fluiddis_,fluid_SearchTree_,*fluiddis_,currpositions_fluid_,boundingboxcenter,structsearchradiusfac*max_radius,0);

    if(curr_subset_of_fluiddis_.empty() == false)
      std::cout<<"\nPrepareFluidOp returns "<<curr_subset_of_fluiddis_.begin()->second.size()<<" background elements on Proc "<<myrank_<<std::endl;
  }
  else
  {
    // get searchbox conditions on bodies
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    std::vector<DRT::Condition*> conditions;
    structdis_->GetCondition("ImmersedSearchbox",conditions);

    // build element list
    std::map<int, std::set<int> > elementList;
    std::set<int> settoinsert;
    for(int i=0;i<(int)conditions.size();++i)
    {
      for (curr=conditions[i]->Geometry().begin(); curr!=conditions[i]->Geometry().end(); ++curr)
      {
        settoinsert.insert(curr->second->Id());
      }
      elementList.insert(std::pair<int, std::set<int> >(i,settoinsert));
      settoinsert.clear();
    }

    // get bounding boxes of the bodies
    std::vector<LINALG::Matrix<3, 2> > structboxes = GEO::computeXAABBForLabeledStructures(*structdis_,currpositions_struct_,elementList);

    double max_radius;

    // search for background elements within a certain radius around the center of the immersed bounding box
    for(int i=0;i<(int)structboxes.size();++i)
    {
      max_radius=sqrt(pow(structboxes[i](0,0)-structboxes[i](0,1),2)+pow(structboxes[i](1,0)-structboxes[i](1,1),2)+pow(structboxes[i](2,0)-structboxes[i](2,1),2));

      LINALG::Matrix<3,1> boundingboxcenter;
      boundingboxcenter(0) = structboxes[i](0,0)+( structboxes[i](0,1)- structboxes[i](0,0))*0.5;
      boundingboxcenter(1) = structboxes[i](1,0)+( structboxes[i](1,1)- structboxes[i](1,0))*0.5;
      boundingboxcenter(2) = structboxes[i](2,0)+( structboxes[i](2,1)- structboxes[i](2,0))*0.5;

      std::map<int, std::set<int> > tempmap = fluid_SearchTree_->searchElementsInRadius(*fluiddis_,currpositions_fluid_,boundingboxcenter,0.5*max_radius,0);
      curr_subset_of_fluiddis_.insert(std::pair<int, std::set<int> >(i,(tempmap.begin()->second)));
    }

    for(int i=0;i<(int)structboxes.size();++i)
      std::cout<<"\nPrepareFluidOp returns "<<curr_subset_of_fluiddis_.at(i).size()<<" background elements for body "<<i<<std::endl;
  }

//  // DEBUG option:
//  // output of sets:
//
//  std::cout<<"\n Structural bounding box diagonal = "<<max_radius<<std::endl;
//  std::cout<<"\n bounding box center coordinate = "<<boundingboxcenter<<std::endl;
//  if(curr_subset_of_fluiddis_.empty() == false)
//  {
//    int counter = 0;
//    for(std::map<int, std::set<int> >::const_iterator closele = curr_subset_of_fluiddis_.begin(); closele != curr_subset_of_fluiddis_.end(); closele++)
//    {
//      for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
//      {
//        counter ++;
//        std::cout<<"PROC"<<myrank_<<" key "<<closele->first<<" eleIter "<<counter<<" : "<<*eleIter<<std::endl;
//      }
//    }
//  }
//  else
//    std::cout<<"curr_subset_of_fluiddis_ empty on PROC "<<myrank_<<std::endl;
//
//  //___________________________
//
//  if(currpositions_struct_.empty() == false)
//  {
//    int counter = 0;
//    for(std::map<int,LINALG::Matrix<3,1> >::const_iterator closele = currpositions_struct_.begin(); closele != currpositions_struct_.end(); closele++)
//    {
//      counter ++;
//      if(myrank_==1)
//      std::cout<<"PROC"<<myrank_<<" key: "<<closele->first<<"  position: "<<counter<<closele->second(0)<<" "<<closele->second(1)<<" "<<closele->second(2)<<std::endl;
//    }
//  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ExtractInterfaceDispnp()
{
  return immersedstructure_->ExtractImmersedInterfaceDispnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> full_traction_vec)
{
  immersedstructure_->ApplyImmersedInterfaceForces(Teuchos::null,full_traction_vec);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::AddDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->AddDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::RemoveDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->RemoveDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CalcResidual(Epetra_Vector &F, const Teuchos::RCP<Epetra_Vector> newstate , const Teuchos::RCP<Epetra_Vector> oldstate)
{
  int err;
  if(!displacementcoupling_)
     err = F.Update(1.0,*(immersedstructure_->Interface()->ExtractIMMERSEDCondVector(newstate)), -1.0, *oldstate, 0.0);
  else
     err = F.Update(1.0,*newstate, -1.0, *oldstate, 0.0);

  return err;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CalcFluidTractionsOnStructure()
{
  // calc new tractions on immersed boundary
  // after leaving the FluidOp the vector struct_bdry_traction_
  // contains the current tractions on the immersed boundary
  Teuchos::ParameterList params;
  params.set<std::string>("action","calc_fluid_traction");
  params.set<std::string>("backgrddisname","fluid");
  params.set<std::string>("immerseddisname","structure");

  SetStatesStructOP();

  DRT::AssembleStrategy struct_bdry_strategy(
      0,              // struct dofset for row
      0,              // struct dofset for column
      Teuchos::null,  // matrix 1
      Teuchos::null,  //
      struct_bdry_traction_ ,  // vector 1
      Teuchos::null,  //
      Teuchos::null   //
  );
  if(myrank_ == 0)
  {
    std::cout<<"################################################################################################"<<std::endl;
    std::cout<<"###   Interpolate fluid stresses to structural surface and calculate tractions                  "<<std::endl;

  }
  EvaluateInterpolationCondition( immersedstructure_->Discretization(), params, struct_bdry_strategy, "IMMERSEDCoupling", -1 );

  double normorstructbdrytraction;
  struct_bdry_traction_->Norm2(&normorstructbdrytraction);
  if(myrank_ == 0)
  {
    std::cout<<"###   Norm of Boundary Traction:   "<<normorstructbdrytraction<<std::endl;
    std::cout<<"################################################################################################"<<std::endl;
  }

}
