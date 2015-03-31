/*!----------------------------------------------------------------------
\file immersed_partitioned_cellmigration.cpp

\brief partitioned immersed cell migration algorithm

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_cellmigration.H"

#include "../linalg/linalg_utils.H"

#include "../drt_poroelast/poroelast_monolithic.H"

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H" //< type of poro structure

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_inpar/inpar_immersed.H"

#include "../drt_fsi/fsi_nox_aitken.H"

#include "../drt_structure/stru_aux.H"

#include <Teuchos_TimeMonitor.hpp>


IMMERSED::ImmersedPartitionedCellMigration::ImmersedPartitionedCellMigration(const Epetra_Comm& comm)
  : ImmersedPartitioned(comm)
{
  // important variables for parallel simulations
  myrank_  = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // get pointer to discretizations
  backgroundfluiddis_     = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_            = globalproblem_->GetDis("cell");

  // decide whether multiple structural bodies or not
  std::vector<DRT::Condition*> conditions;
  immerseddis_->GetCondition("ImmersedSearchbox",conditions);
  if((int)conditions.size()>0)
  {
    if(myrank_==0)
      std::cout<<" MULTI CELL MIGRATION SIMULATION   Number of cells: "<<(int)conditions.size()<<std::endl;
    multicellmigration_ = true;
  }
  else
    multicellmigration_ = false;

  // setup the relaxation parameters
  SetupRelaxation();

  // build field structure
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> cellstructure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(globalproblem_->CellMigrationParams(), const_cast<Teuchos::ParameterList&>(globalproblem_->StructuralDynamicParams()), immerseddis_));
    cellstructure_ = Teuchos::rcp_dynamic_cast< ::ADAPTER::FSIStructureWrapper>(cellstructure->StructureField());
    std::cout<<"Created Field Cell Sturucture ..."<<std::endl;

  // create instance of poroelast subproblem
  poroelast_subproblem_ = Teuchos::rcp(new POROELAST::Monolithic(comm, globalproblem_->PoroelastDynamicParams()));
  // setup of poro monolithic algorithm
  SetupPoroAlgorithm();
  if(myrank_==0)
    std::cout<<"Created Field Poroelast ..."<<std::endl;

  // vector of fluid stresses interpolated to cell bdry int points and integrated over cell surface
  cell_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // vector with fluid velocities interpolated from structure
  porofluid_artificial_velocity_ = Teuchos::rcp(new Epetra_Vector(*(poroelast_subproblem_->FluidField()->DofRowMap()),true));

  // construct 3D search tree for fluid domain
  // initialized in SetupSBackgroundDiscretization()
  fluid_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // construct 3D search tree for cell domain
  // initialized in SetupImmersedDiscretization()
  cell_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // get coupling variable
  displacementcoupling_ = globalproblem_->ImmersedMethodParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE") == "Displacement";
  if(displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Force "<<std::endl;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  // reinitialize the transfer vectors
  porofluid_artificial_velocity_->PutScalar(0.0);
  cell_bdry_traction_->PutScalar(0.0);

  // reset element and node information about immersed method
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::reset_immersed_ele);
  params.set<int>("Physical Type", poroelast_subproblem_->FluidField()->PhysicalType());
  backgroundfluiddis_->Evaluate(params);


  if (displacementcoupling_) // DISPLACEMENT COUPLING
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp();
    BackgroundOp(porofluid_artificial_velocity_, fillFlag);



    // reset element and node information about immersed method
    params.set<int>("action",FLD::reset_immersed_ele);
    params.set<int>("reset_isimmersed",0);
    params.set<int>("reset_hasprojecteddirichlet",0);
    poroelast_subproblem_->FluidField()->Discretization()->Evaluate(params);

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    Teuchos::RCP<Epetra_Vector> idispnp =
    ImmersedOp(cell_bdry_traction_, fillFlag);

    int err = F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
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
    ImmersedOp(cell_bdry_traction_, fillFlag);

    ////////////////////
    // CALL FluidOp
    ////////////////////
    PrepareBackgroundOp(); // determine elements cut by the boundary
    BackgroundOp(porofluid_artificial_velocity_, fillFlag);

    // reset element and node information about immersed method
    params.set<int>("action",FLD::reset_immersed_ele);
    params.set<int>("reset_isimmersed",0);
    params.set<int>("reset_hasprojecteddirichlet",0);
    backgroundfluiddis_->Evaluate(params);

    ///////////////////////////////////////////////////
    // CALL StructOp again to recalculate tractions (only tractions due to StructOp( , User))
    //////////////////////////////////////////////////
    const Teuchos::RCP<Epetra_Vector> iforcenp = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
    ImmersedOp(iforcenp, User);

    int err = F.Update(1.0, *(cellstructure_->Interface()->ExtractFSICondVector(iforcenp)), -1.0, *iforcen, 0.0);

    if(err != 0)
      dserror("Vector update of FSI-residual returned err=%d",err);

  }

  // perform n steps max; then set converged
  bool nlnsolver_continue = globalproblem_->ImmersedMethodParams().get<std::string>("DIVERCONT") == "continue";
  int  itemax = globalproblem_->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<int>("ITEMAX");
  if((IterationCounter())[0] == itemax and nlnsolver_continue)
  {
    if(Comm().MyPID()==0)
      std::cout<<"\n  Continue with next time step after ITEMAX = "<<(IterationCounter())[0]<<" iterations. \n"<<std::endl;

    // !!! EXPERIMENTAL !!!
    // set F to zero to tell NOX that this timestep is converged
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp( new Epetra_Vector(F.Map(),true));
    F.Update(1.0,*zeros,0.0);
    // !!! EXPERIMENTAL !!!

    immerseddis_->ClearState();
    backgroundfluiddis_->ClearState();
  }

  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
                                                              const FillType fillFlag)
{
  // print
  IMMERSED::ImmersedPartitioned::BackgroundOp(backgrd_dirichlet_values,fillFlag);

  if (fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    // calc the fluid velocity from the cell displacements
    DRT::AssembleStrategy fluid_vol_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        backgrd_dirichlet_values ,  // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );


    // relaxation
    if(displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::fixed)
    {
      std::cout<<"Apply relaxation parameter omega = "<<velrelax_<<" on full artificial velocity field."<<std::endl;
      backgrd_dirichlet_values->Scale(velrelax_);
    }
    else if(displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::aitken)
    {
      double omega_aitken = Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(aitkenfactory_)->GetAitken()->getOmega();
      std::cout<<"Apply AITKEN relaxation parameter omega = "<<omega_aitken<<" on full artificial velocity field."<<std::endl;

      backgrd_dirichlet_values->Scale(omega_aitken);
    }
    else if(displacementcoupling_ and relaxvelglobally_ == false)
      std::cout<<"Selective velocity relaxation not implemented yet. Parameter VEL_RELAX = "<<velrelax_<<" has no effect."<<std::endl;


    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate Dirichlet Values from immersed elements which overlap the "<<backgroundfluiddis_->Name()<<" nodes ..."<<std::endl;
    }

    if(curr_subset_of_backgrounddis_.empty()==false)
      EvaluateWithInternalCommunication(backgroundfluiddis_,&fluid_vol_strategy, curr_subset_of_backgrounddis_, cell_SearchTree_, currpositions_cell_);

    // dofsetnum of backgroundfluid is 0 for backgroundfluiddis in poroelast_subproblem_
    BuildImmersedDirichMap(backgroundfluiddis_, dbcmap_immersed_, poroelast_subproblem_->FluidField()->GetDBCMapExtractor()->CondMap(),0);
    poroelast_subproblem_->FluidField()->AddDirichCond(dbcmap_immersed_);

    // apply immersed dirichlets to porofluid field
    DoImmersedDirichletCond(poroelast_subproblem_->FluidField()->WriteAccessVelnp(),backgrd_dirichlet_values, dbcmap_immersed_);
    // rebuild the combined dbcmap
    poroelast_subproblem_->BuildCombinedDBCMap();

    double normofvelocities;
    poroelast_subproblem_->FluidField()->ExtractVelocityPart(backgrd_dirichlet_values)->Norm2(&normofvelocities);
//
//    // apply pressure neumann based on divergence of structural field above fluid nodes
//    Teuchos::RCP<const Epetra_Vector> fluid_artificial_veldiv = MBFluidField()->FluidField()->ExtractPressurePart(fluid_artificial_velocity);
//    double normofdivergence;
//    fluid_artificial_veldiv->Norm2(&normofdivergence);

    if(myrank_ == 0)
    {
      std::cout<<"###   Norm of Dirichlet values:   "<<std::setprecision(7)<<normofvelocities<<std::endl;
      std::cout<<"###   Norm of transferred divergence of cell velocity:   "<<std::setprecision(10)<<"not used"<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }


    // solve poro
    //MBFluidField()->FluidField()->AddContributionToNeumannLoads(fluid_artificial_veldiv);
    poroelast_subproblem_->Solve();

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration)
    poroelast_subproblem_->FluidField()->RemoveDirichCond(dbcmap_immersed_);
    // rebuild the combined dbcmap
    poroelast_subproblem_->BuildCombinedDBCMap();

  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
                                                       const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::ImmersedOp(bdry_traction,fillFlag);

  if(fillFlag==User)
  {
    Teuchos::ParameterList params;
    params.set<std::string>("action","calc_fluid_traction");
    params.set<std::string>("backgrddisname","porofluid");
    params.set<std::string>("immerseddisname","cell");

    backgroundfluiddis_->SetState(0,"velnp",poroelast_subproblem_->FluidField()->Velnp());
    backgroundfluiddis_->SetState(0,"veln", poroelast_subproblem_->FluidField()->Veln());

    DRT::AssembleStrategy struct_bdry_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        bdry_traction , // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );
    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Calculate Boundary Tractions on Structure for Initial Guess or Convergence Check ...      "<<std::endl;

    }
    EvaluateInterpolationCondition( immerseddis_, params, struct_bdry_strategy,"FSICoupling", -1 );

    double normorstructbdrytraction;
    bdry_traction->Norm2(&normorstructbdrytraction);
    if(myrank_ == 0)
    {
      std::cout<<"###   Norm of Boundary Traction:   "<<normorstructbdrytraction<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }

    return Teuchos::null;
  }
  else
  {
    Teuchos::ParameterList params;
    params.set<std::string>("action","calc_fluid_traction");
    params.set<std::string>("backgrddisname","porofluid");
    params.set<std::string>("immerseddisname","cell");

    backgroundfluiddis_->SetState(0,"velnp",poroelast_subproblem_->FluidField()->Velnp());
    backgroundfluiddis_->SetState(0,"veln", poroelast_subproblem_->FluidField()->Veln());

    DRT::AssembleStrategy struct_bdry_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        bdry_traction , // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );
    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate fluid stresses to structural surface and calculate tractions                  "<<std::endl;

    }
    EvaluateInterpolationCondition( immerseddis_, params, struct_bdry_strategy, "FSICoupling", -1 );

    double normorstructbdrytraction;
    bdry_traction->Norm2(&normorstructbdrytraction);
    if(myrank_ == 0)
    {
      std::cout<<"###   Norm of Boundary Traction:   "<<normorstructbdrytraction<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }

    // relaxation
    if(!displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::fixed)
    {
      std::cout<<"Apply relaxation parameter omega = "<<velrelax_<<" on whole structural force interface."<<std::endl;
      bdry_traction->Scale(velrelax_);
    }
    else if(!displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::aitken)
    {
      double omega_aitken = Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(aitkenfactory_)->GetAitken()->getOmega();
      std::cout<<"Apply AITKEN relaxation parameter omega = "<<omega_aitken<<" on whole structural force interface."<<std::endl;

      bdry_traction->Scale(omega_aitken);
    }
    else if(!displacementcoupling_ and relaxvelglobally_ == false)
      std::cout<<"Selective velocity relaxation not implemented yet. Parameter VEL_RELAX = "<<velrelax_<<" has no effect."<<std::endl;

    // prescribe neumann values at structural boundary dofs
    cellstructure_->ApplyInterfaceForces(cellstructure_->Interface()->ExtractFSICondVector(*bdry_traction));

    // solve
    cellstructure_->Solve();

    return cellstructure_->ExtractInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if(myrank_==0)
    std::cout<<"Cell Predictor: "<<std::endl;
  cellstructure_->PrepareTimeStep();
  if(myrank_==0)
    std::cout<<"Poro Predictor: "<<std::endl;
  poroelast_subproblem_->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::InitialGuess()
{
  if(displacementcoupling_)
    return cellstructure_->PredictInterfaceDispnp();
  else
  {
    immerseddis_->SetState(0,"displacement",cellstructure_->Dispnp());
    cell_bdry_traction_old_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
    ImmersedOp(cell_bdry_traction_old_, User);

    return cellstructure_->Interface()->ExtractFSICondVector(cell_bdry_traction_old_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
                                                                        Teuchos::RCP<Epetra_Map>& dirichmap,
                                                                        const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
                                                                        int dofsetnum)
{
  const Epetra_Map* elerowmap = dis->ElementRowMap();
  std::vector<int> mydirichdofs(0);

  for(int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(dis->gElement(elerowmap->GID(i)));
    if(immersedele->HasProjectedDirichlet())
    {
      DRT::Node** nodes = immersedele->Nodes();
      for (int inode=0; inode<(immersedele->NumNode()); inode++)
      {
        if(static_cast<IMMERSED::ImmersedNode* >(nodes[inode])->IsMatched())
        {
          std::vector<int> dofs = dis->Dof(dofsetnum,nodes[inode]);

          for (int dim=0;dim<3;++dim)
          {
            if(dirichmap_original->LID(dofs[dim]) == -1) // if not already in original dirich map
              mydirichdofs.push_back(dofs[dim]);
          }

          // include also pressure dof if node does not belong to a boundary background element
          //if((nodes[inode]->IsImmersedBoundary())==0)
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
void IMMERSED::ImmersedPartitionedCellMigration::DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap)
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
void IMMERSED::ImmersedPartitionedCellMigration::SetupImmersedDiscretization()
{
  // ghost structure on each proc (for search algorithm)
  if(numproc_ > 1)
  {
    // fill complete inside
    CreateGhosting(immerseddis_);
  }
  else
  {
    // fill complete to incorporate changes due to ghosting and build geometries
    immerseddis_->FillComplete();
  }

  // find positions of the immersed discretization on proc
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_cell;
  for (int lid = 0; lid < immerseddis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = immerseddis_->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    my_currpositions_cell[node->Id()] = currpos;
  }
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[numproc_];
  for(int i=0;i<numproc_;i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_cell,currpositions_cell_,numproc_,&procs[0],Comm());

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox2 = GEO::getXAABBofDis(*immerseddis_,currpositions_cell_);
  cell_SearchTree_->initializeTree(rootBox2,*immerseddis_,GEO::TreeType(GEO::OCTTREE));

  if(myrank_==0)
    std::cout<<"\n Build Cell SearchTree ... "<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::SetupBackgroundDiscretization()
{
  // find positions of the background fluid discretization
  for (int lid = 0; lid < backgroundfluiddis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = backgroundfluiddis_->lColNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    currpositions_ECM_[node->Id()] = currpos;
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*backgroundfluiddis_,currpositions_ECM_);
  fluid_SearchTree_->initializeTree(rootBox,*backgroundfluiddis_,GEO::TreeType(GEO::OCTTREE));

  if(myrank_==0)
    std::cout<<"\n Build Fluid SearchTree ... "<<std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareImmersedOp()
{
  // not used yet
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareBackgroundOp()
{
  TEUCHOS_FUNC_TIME_MONITOR("IMMERSED::PrepareBackgroundOp()");

  immerseddis_->SetState(0,"displacement",cellstructure_->Dispnp());
  immerseddis_->SetState(0,"displacement_old",cellstructure_->Dispn());
  immerseddis_->SetState(0,"velocity_old",cellstructure_->Veln());
  immerseddis_->SetState(0,"velocity",cellstructure_->Velnp());
  backgroundfluiddis_ ->SetState(0,"veln",poroelast_subproblem_->FluidField()->Veln());

  double structsearchradiusfac = DRT::Problem::Instance()->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

  //
  // determine subset of fluid discretization which is potentially underlying the immersed discretization
  //
  // get state
  Teuchos::RCP<const Epetra_Vector> celldisplacements = cellstructure_->Dispnp();

  // find current positions for immersed structural discretization
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_cell;
  for (int lid = 0; lid < immerseddis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = immerseddis_->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;
    std::vector<int> dofstoextract(3);
    std::vector<double> mydisp(3);

    // get the current displacement
    immerseddis_->Dof(node,0,dofstoextract);
    DRT::UTILS::ExtractMyValues(*celldisplacements,mydisp,dofstoextract);

    currpos(0) = node->X()[0]+mydisp.at(0);
    currpos(1) = node->X()[1]+mydisp.at(1);
    currpos(2) = node->X()[2]+mydisp.at(2);

    my_currpositions_cell[node->Id()] = currpos;
  }
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[Comm().NumProc()];
  for(int i=0;i<Comm().NumProc();i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_cell,currpositions_cell_,Comm().NumProc(),&procs[0],Comm());

//  DEBUG output
//  std::cout<<"Proc "<<Comm().MyPID()<<": my_curr_pos.size()="<<my_currpositions_cell.size()<<std::endl;
//  std::cout<<"Proc "<<Comm().MyPID()<<":    curr_pos.size()="<<currpositions_struct_.size()<<std::endl;

  if (multicellmigration_ == false)
  {
    // get bounding box of current configuration of structural dis
    const LINALG::Matrix<3,2> structBox = GEO::getXAABBofDis(*immerseddis_,currpositions_cell_);
    double max_radius = sqrt(pow(structBox(0,0)-structBox(0,1),2)+pow(structBox(1,0)-structBox(1,1),2)+pow(structBox(2,0)-structBox(2,1),2));
    // search for background elements within a certain radius around the center of the immersed bounding box
    LINALG::Matrix<3,1> boundingboxcenter;
    boundingboxcenter(0) = structBox(0,0)+(structBox(0,1)-structBox(0,0))*0.5;
    boundingboxcenter(1) = structBox(1,0)+(structBox(1,1)-structBox(1,0))*0.5;
    boundingboxcenter(2) = structBox(2,0)+(structBox(2,1)-structBox(2,0))*0.5;

# ifdef DEBUG
    std::cout<<"Bounding Box of Structure: "<<structBox<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Bounding Box Center of Structure: "<<boundingboxcenter<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Search Radius Around Center: "<<structsearchradiusfac*max_radius<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Length of Dispnp()="<<celldisplacements->MyLength()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Size of currpositions_struct_="<<currpositions_cell_.size()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"My DofRowMap Size="<<immerseddis_->DofRowMap()->NumMyElements()<<"  My DofColMap Size="<<immerseddis_->DofColMap()->NumMyElements()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Dis Structure NumColEles: "<<immerseddis_->NumMyColElements()<<" on PROC "<<Comm().MyPID()<<std::endl;
# endif

    curr_subset_of_backgrounddis_ = fluid_SearchTree_->searchElementsInRadius(*backgroundfluiddis_,currpositions_ECM_,boundingboxcenter,structsearchradiusfac*max_radius,0);

    if(curr_subset_of_backgrounddis_.empty() == false)
      std::cout<<"\nPrepareBackgroundOp returns "<<curr_subset_of_backgrounddis_.begin()->second.size()<<" background elements on Proc "<<Comm().MyPID()<<std::endl;
  }
  else
  {
    dserror("multi cell migration not implemented yet ... ");
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::SetupPoroAlgorithm()
{
  poroelast_subproblem_->SetupSystem();
  poroelast_subproblem_->SetupSolver();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareOutput()
{
  poroelast_subproblem_->PrepareOutput();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::Output()
{
  cellstructure_->Output();
  poroelast_subproblem_->Output();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::Update()
{
  cellstructure_->Update();
  poroelast_subproblem_->Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::SetupRelaxation()
{
  // initialize some relaxation related member variables
    relaxforceglobally_ = globalproblem_->ImmersedMethodParams().get<std::string>("APPLY_FORCE_RELAX") == "globally";
    relaxvelglobally_   = globalproblem_->ImmersedMethodParams().get<std::string>("APPLY_VEL_RELAX")   == "globally";
    forcerelax_         = globalproblem_->ImmersedMethodParams().get<double>("FORCE_RELAX");
    velrelax_           = globalproblem_->ImmersedMethodParams().get<double>("VEL_RELAX");

    if(globalproblem_->CellMigrationParams().get<std::string>("COUPALGO") == "iter_stagg_fixed_rel_param")
    {
      coupalgo_= INPAR::IMMERSED::fixed;
      if(myrank_==0)
        std::cout<<"\n"
                 " Relax Force globally = "<<relaxforceglobally_<<" with relaxation parameter = "<<forcerelax_<<"\n"
                 " Relax Vel   globally = "<<relaxvelglobally_  <<" with relaxation parameter = "<<velrelax_<<"\n"<<std::endl;
    }
    else if(globalproblem_->CellMigrationParams().get<std::string>("COUPALGO") == "iter_stagg_AITKEN_rel_param")
    {
      coupalgo_ = INPAR::IMMERSED::aitken;
      if(myrank_==0)
        std::cout<<"\n Using AITKEN relaxation parameter. "<<std::endl;
    }
    else
      dserror("Unknown definition of COUPALGO in FSI DYNAMIC section for Immersed FSI.");

  return;
}


