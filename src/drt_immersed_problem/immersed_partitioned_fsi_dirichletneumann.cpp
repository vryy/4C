/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann.cpp

\brief partitioned immersed fsi algorithm for neumann-neumann like coupling (volume force coupling)

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_dirichletneumann.H"
#include "immersed_partitioned_fsi.H"
#include "fsi_partitioned_immersed.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_assemblestrategy.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_immersed.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_so3/so_surface.H"
#include "../drt_fluid/fluidimplicitintegration.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_immersed.H"

#include "../linalg/linalg_utils.H"

// search tree related
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

// time monitoring
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

// relaxation
#include "../drt_fsi/fsi_nox_aitken.H"


IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ImmersedPartitionedFSIDirichletNeumann(const Epetra_Comm& comm)
  : ImmersedBase(),
    FSI::PartitionedImmersed(comm)
{
  int myrank = comm.MyPID();

  double fsirelax = DRT::Problem::Instance()->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<double>("RELAX");
  if (fsirelax != 1.0 and myrank==0)
  {
    std::cout<<"!!!! WARNING !!! \n"
               "Relaxation parameter set in FSI DYNAMIC/PARTITIONED SOLVER section is not equal 1.0.\n"
               "This is not possible in Immersed FSI. Adapt the relaxation parameters in IMMERSED METHOD section, instead.\n"
               "NOX relaxes only the interface increment which is not used in Immersed FSI, since the whole artificial fluid \n"
               "volume gets dirichlet values."<<std::endl;
    dserror("Invalid parameter");
  }

  // initialize some relaxation related member variables
  relaxforceglobally_ = DRT::Problem::Instance()->ImmersedMethodParams().get<std::string>("APPLY_FORCE_RELAX") == "globally";
  relaxvelglobally_   = DRT::Problem::Instance()->ImmersedMethodParams().get<std::string>("APPLY_VEL_RELAX")   == "globally";
  forcerelax_ = DRT::Problem::Instance()->ImmersedMethodParams().get<double>("FORCE_RELAX");
  velrelax_ = DRT::Problem::Instance()->ImmersedMethodParams().get<double>("VEL_RELAX");

  if(DRT::Problem::Instance()->FSIDynamicParams().get<std::string>("COUPALGO") == "iter_stagg_fixed_rel_param")
  {
    coupalgo_= INPAR::IMMERSED::fixed;
    if(myrank==0)
      std::cout<<"\n"
               " Relax Force globally = "<<relaxforceglobally_<<" with relaxation parameter = "<<forcerelax_<<"\n"
               " Relax Vel   globally = "<<relaxvelglobally_  <<" with relaxation parameter = "<<velrelax_<<"\n"<<std::endl;
  }
  else if(DRT::Problem::Instance()->FSIDynamicParams().get<std::string>("COUPALGO") == "iter_stagg_AITKEN_rel_param")
  {
    coupalgo_ = INPAR::IMMERSED::aitken;
    if(myrank==0)
      std::cout<<"\n Using AITKEN relaxation parameter. "<<std::endl;
  }
  else
    dserror("Unknown definition of COUPALGO in FSI DYNAMIC section for Immersed FSI.");


  // get coupling variable
  displacementcoupling_ = DRT::Problem::Instance()->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE") == "Displacement";
  if(displacementcoupling_ and myrank==0)
    std::cout<<" Coupling variable for partitioned FSI scheme :  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank==0)
    std::cout<<" Coupling variable for partitioned FSI scheme :  Force "<<std::endl;

  // get pointer to discretizations
  fluiddis_  = DRT::Problem::Instance()->GetDis("fluid");
  structdis_ = DRT::Problem::Instance()->GetDis("structure");

  // decide whether multiple structural bodies or not
  std::vector<DRT::Condition*> conditions;
  structdis_->GetCondition("ImmersedSearchbox",conditions);
  if((int)conditions.size()>0)
  {
    if(myrank==0)
      std::cout<<" MULTIBODY SIMULATION   Number of bodies: "<<(int)conditions.size()<<std::endl;
    multibodysimulation_ = true;
  }
  else
    multibodysimulation_ = false;

  // vector of fluid stresses interpolated to structural bdry int points and integrated over structural surface
  struct_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(StructureField()->DofRowMap()),true));

  // vector with fluid velocities interpolated from structure
  fluid_artificial_velocity_ = Teuchos::rcp(new Epetra_Vector(*(MBFluidField()->FluidField()->DofRowMap()),true));


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

  if(myrank==0)
    std::cout<<"\n Build Fluid SearchTree ... "<<std::endl;


  // construct 3D search tree for structural domain
  // initialized in SetupStructuralDiscretization()
  structure_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));


  comm.Barrier();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  // reinitialize the transfer vectors
  fluid_artificial_velocity_->PutScalar(0.0);
  struct_bdry_traction_->PutScalar(0.0);

  // reset element and node information about immersed method
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::reset_immersed_ele);
  MBFluidField()->Discretization()->Evaluate(params);


  if (displacementcoupling_)
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL FluidOp
    ////////////////////
    PrepareFluidOp();
    FluidOp(fluid_artificial_velocity_, fillFlag);

    // reset element and node information about immersed method
    params.set<int>("action",FLD::reset_immersed_ele);
    params.set<int>("reset_isimmersed",0);
    params.set<int>("reset_hasprojecteddirichlet",0);
    MBFluidField()->Discretization()->Evaluate(params);

    ////////////////////
    // CALL StructOp
    ////////////////////
    Teuchos::RCP<Epetra_Vector> idispnp =
    StructOp(struct_bdry_traction_, fillFlag);

// DEBUG option:
// Allows to write a certain nonlinear iteration step
// no matter if it is really converged or not.
//#ifdef DEBUG
//    //////////////////////////////
//    // DEBUG (dummy zero residual)
//    /////////////////////////////
//    int steptowrite = 2;
//
//    if ((FSI::Partitioned::IterationCounter())[0] == steptowrite)
//    {
//      Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(F.Map(),true));
//      F.Update(1.0,*zeros,0.0);
//    }
//    else
//    {
//      Teuchos::RCP<Epetra_Vector> ones = Teuchos::rcp(new Epetra_Vector(F.Map(),true));
//      ones->PutScalar(1.0);
//      F.Update(1.0,*ones,0.0);
//    }
//#else
  int err = F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);

  if(err != 0)
    dserror("Vector update of FSI-residual returned err=%d",err);
//
//#endif

//    double norm;
//    F.Norm2(&norm);
//    std::cout<<"\n PARTITIONED FSI RESIDUAL (L2-Norm) = "<<norm<<std::endl;

  }
  else if(!displacementcoupling_) // FORCE COUPLING
  {
    // get the initial guess
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));


    ////////////////////
    // CALL StructOp
    ////////////////////
    StructOp(struct_bdry_traction_, fillFlag);

    ////////////////////
    // CALL FluidOp
    ////////////////////
    PrepareFluidOp(); // determine elements cut by the boundary
    FluidOp(fluid_artificial_velocity_, fillFlag);

    // reset element and node information about immersed method
    params.set<int>("action",FLD::reset_immersed_ele);
    params.set<int>("reset_isimmersed",0);
    params.set<int>("reset_hasprojecteddirichlet",0);
    MBFluidField()->Discretization()->Evaluate(params);

    ///////////////////////////////////////////////////
    // CALL StructOp again to recalculate tractions (only tractions due to StructOp( , User))
    //////////////////////////////////////////////////
    const Teuchos::RCP<Epetra_Vector> iforcenp = Teuchos::rcp(new Epetra_Vector(*(StructureField()->DofRowMap()),true));
    StructOp(iforcenp, User);

    int err = F.Update(1.0, *(StructureField()->Interface()->ExtractFSICondVector(iforcenp)), -1.0, *iforcen, 0.0);

    if(err != 0)
      dserror("Vector update of FSI-residual returned err=%d",err);

  }


  // perform n steps max; then set converged
  bool nlnsolver_continue = DRT::Problem::Instance()->ImmersedMethodParams().get<std::string>("DIVERCONT") == "continue";
  int  itemax = DRT::Problem::Instance()->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<int>("ITEMAX");
  if((FSI::Partitioned::IterationCounter())[0] == itemax and nlnsolver_continue)
  {
    if(Comm().MyPID()==0)
      std::cout<<"\n  Continue with next time step after "<<(FSI::Partitioned::IterationCounter())[0]<<" iterations. \n"<<std::endl;

    // !!! EXPERIMENTAL !!!
    // set F to zero to tell NOX that this timestep is converged
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp( new Epetra_Vector(F.Map(),true));
    F.Update(1.0,*zeros,0.0);
    // !!! EXPERIMENTAL !!!

    StructureField()->Discretization()->ClearState();
    MBFluidField()->Discretization()->ClearState();
  }

  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FluidOp(Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity,
                               const FillType fillFlag)
{
  // print
  FSI::Partitioned::FluidOp(fluid_artificial_velocity,fillFlag);

  if (fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    // normal fluid solve

    // A rather simple hack. We need something better!
    const int itemax = MBFluidField()->Itemax();
    if (fillFlag==MF_Res and mfresitemax_ > 0)
      MBFluidField()->SetItemax(mfresitemax_ + 1);

    // calc the fluid velocity from the structural displacements
    DRT::AssembleStrategy fluid_vol_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        fluid_artificial_velocity ,  // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );


    // relaxation
    if(displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::fixed)
    {
      std::cout<<"Apply relaxation parameter omega = "<<velrelax_<<" on full artificial velocity field."<<std::endl;
      fluid_artificial_velocity->Scale(velrelax_);
    }
    else if(displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::aitken)
    {
      double omega_aitken = Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(aitkenfactory_)->GetAitken()->getOmega();
      std::cout<<"Apply AITKEN relaxation parameter omega = "<<omega_aitken<<" on full artificial velocity field."<<std::endl;

      fluid_artificial_velocity->Scale(omega_aitken);
    }
    else if(displacementcoupling_ and relaxvelglobally_ == false)
      std::cout<<"Selective velocity relaxation not implemented yet. Parameter VEL_RELAX = "<<velrelax_<<" has no effect."<<std::endl;


    if(StructureField()->Discretization()->Comm().MyPID() == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate Dirichlet Values from immersed elements which overlap the "<<MBFluidField()->Discretization()->Name()<<" nodes ..."<<std::endl;
    }

    if(curr_subset_of_fluiddis_.empty()==false)
      EvaluateWithInternalCommunication(MBFluidField()->Discretization(),&fluid_vol_strategy, curr_subset_of_fluiddis_, structure_SearchTree_, currpositions_struct_);

    BuildImmersedDirichMap(MBFluidField()->Discretization(), dbcmap_immersed_, MBFluidField()->FluidField()->GetDBCMapExtractor()->CondMap());
    Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->AddDirichCond(dbcmap_immersed_);

    // apply immersed dirichlets
    DoImmersedDirichletCond(MBFluidField()->FluidField()->WriteAccessVelnp(),fluid_artificial_velocity_, dbcmap_immersed_);
    double normofvelocities;
    MBFluidField()->FluidField()->ExtractVelocityPart(fluid_artificial_velocity)->Norm2(&normofvelocities);

    // apply pressure neumann based on divergence of structural field above fluid nodes
    Teuchos::RCP<const Epetra_Vector> fluid_artificial_veldiv = MBFluidField()->FluidField()->ExtractPressurePart(fluid_artificial_velocity);
    double normofdivergence;
    fluid_artificial_veldiv->Norm2(&normofdivergence);

    if(StructureField()->Discretization()->Comm().MyPID() == 0)
    {
      std::cout<<"###   Norm of Dirichlet Values:   "<<std::setprecision(7)<<normofvelocities<<std::endl;
      std::cout<<"###   Norm of transferred divergence of structural velocity:   "<<std::setprecision(10)<<normofdivergence<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }


    // solve fluid
    MBFluidField()->FluidField()->AddContributionToNeumannLoads(fluid_artificial_veldiv);
    MBFluidField()->NonlinearSolve(Teuchos::null,Teuchos::null);

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration)
    Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->RemoveDirichCond(dbcmap_immersed_);

    MBFluidField()->SetItemax(itemax);
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
    Teuchos::ParameterList params;
    params.set<std::string>("action","calc_fluid_traction");

    MBFluidField()->Discretization()->SetState(0,"velnp",Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->FluidField()->Velnp());
    MBFluidField()->Discretization()->SetState(0,"veln",Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->FluidField()->Veln());

    DRT::AssembleStrategy struct_bdry_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        struct_bdry_traction ,  // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );
    if(StructureField()->Discretization()->Comm().MyPID() == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Calculate Boundary Tractions on Structure for Initial Guess or Convergence Check ...      "<<std::endl;

    }
    EvaluateInterpolationCondition( StructureField()->Discretization(), params, struct_bdry_strategy,"FSICoupling", -1 );

    double normorstructbdrytraction;
    struct_bdry_traction->Norm2(&normorstructbdrytraction);
    if(StructureField()->Discretization()->Comm().MyPID() == 0)
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

    MBFluidField()->Discretization()->SetState(0,"velnp",Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->FluidField()->Velnp());
    MBFluidField()->Discretization()->SetState(0,"veln",Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->FluidField()->Veln());

    DRT::AssembleStrategy struct_bdry_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        struct_bdry_traction ,  // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );
    if(StructureField()->Discretization()->Comm().MyPID() == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate fluid stresses to structural surface and calculate tractions                  "<<std::endl;

    }
    EvaluateInterpolationCondition( StructureField()->Discretization(), params, struct_bdry_strategy, "FSICoupling", -1 );

    double normorstructbdrytraction;
    struct_bdry_traction->Norm2(&normorstructbdrytraction);
    if(StructureField()->Discretization()->Comm().MyPID() == 0)
    {
      std::cout<<"###   Norm of Boundary Traction:   "<<normorstructbdrytraction<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }

    // relaxation
    if(!displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::fixed)
    {
      std::cout<<"Apply relaxation parameter omega = "<<velrelax_<<" on whole structural force interface."<<std::endl;
      struct_bdry_traction->Scale(velrelax_);
    }
    else if(!displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::aitken)
    {
      double omega_aitken = Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(aitkenfactory_)->GetAitken()->getOmega();
      std::cout<<"Apply AITKEN relaxation parameter omega = "<<omega_aitken<<" on whole structural force interface."<<std::endl;

      struct_bdry_traction->Scale(omega_aitken);
    }
    else if(!displacementcoupling_ and relaxvelglobally_ == false)
      std::cout<<"Selective velocity relaxation not implemented yet. Parameter VEL_RELAX = "<<velrelax_<<" has no effect."<<std::endl;

    // prescribe neumann values at structural boundary dofs
    StructureField()->ApplyInterfaceForces(StructureField()->Interface()->ExtractFSICondVector(*struct_bdry_traction));

    // solve
    StructureField()->Solve();

    return StructureField()->ExtractInterfaceDispnp();
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::InitialGuess()
{

  if(displacementcoupling_)
    return StructureField()->PredictInterfaceDispnp();
  else
  {
    StructureField()->Discretization()->SetState(0,"displacement",StructureField()->Dispnp());
    struct_bdry_traction_old_ = Teuchos::rcp(new Epetra_Vector(*(StructureField()->DofRowMap()),true));
    StructOp(struct_bdry_traction_old_, User);

    return StructureField()->Interface()->ExtractFSICondVector(struct_bdry_traction_old_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& dirichmap,const Teuchos::RCP<const Epetra_Map>& dirichmap_original)
{
  const Epetra_Map* elerowmap = dis->ElementRowMap();
  std::vector<int> mydirichdofs(0);

  for(int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    DRT::ELEMENTS::FluidImmersed* immersedele = static_cast<DRT::ELEMENTS::FluidImmersed*>(dis->gElement(elerowmap->GID(i)));
    if(immersedele->HasProjectedDirichlet())
    {
      DRT::Node** nodes = immersedele->Nodes();
      for (int inode=0; inode<(immersedele->NumNode()); inode++)
      {
        if(static_cast<IMMERSED::ImmersedNode* >(nodes[inode])->IsMatched())
        {
          std::vector<int> dofs = dis->Dof(nodes[inode]);

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
  int numproc = Comm().NumProc();

  // ghost structure on each proc (for search algorithm)
  if(numproc > 1)
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
  int procs[numproc];
  for(int i=0;i<numproc;i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_struct,currpositions_struct_,numproc,&procs[0],Comm());

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox2 = GEO::getXAABBofDis(*structdis_,currpositions_struct_);
  structure_SearchTree_->initializeTree(rootBox2,*structdis_,GEO::TreeType(GEO::OCTTREE));

  if(Comm().MyPID()==0)
    std::cout<<"\n Build Structure SearchTree ... "<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::PrepareFluidOp()
{
  TEUCHOS_FUNC_TIME_MONITOR("IMMERSED::PrepareFluidOp()");

  structdis_->SetState(0,"displacement",StructureField()->Dispnp());
  structdis_->SetState(0,"displacement_old",StructureField()->Dispn());
  structdis_->SetState(0,"velocity_old",StructureField()->Veln());
  structdis_->SetState(0,"velocity",StructureField()->Velnp());
  fluiddis_ ->SetState(0,"veln",Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed >(MBFluidField())->FluidField()->Veln());

  double structsearchradiusfac = DRT::Problem::Instance()->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

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
//  EvaluateInterpolationCondition( StructureField()->Discretization(), params, dummy_strategy,"FSICoupling", -1 );
# endif

  //
  // determine subset of fluid discretization which is potentially underlying the structural discretization
  //
  // get state
  Teuchos::RCP<const Epetra_Vector> displacements = StructureField()->Dispnp();
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
  int procs[Comm().NumProc()];
  for(int i=0;i<Comm().NumProc();i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_struct,currpositions_struct_,Comm().NumProc(),&procs[0],Comm());

//  DEBUG output
//  std::cout<<"Proc "<<Comm().MyPID()<<": my_curr_pos.size()="<<my_currpositions_struct.size()<<std::endl;
//  std::cout<<"Proc "<<Comm().MyPID()<<":    curr_pos.size()="<<currpositions_struct_.size()<<std::endl;

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
    std::cout<<"Bounding Box of Structure: "<<structBox<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Bounding Box Center of Structure: "<<boundingboxcenter<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Search Radius Around Center: "<<structsearchradiusfac*max_radius<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Length of Dispnp()="<<displacements->MyLength()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Size of currpositions_struct_="<<currpositions_struct_.size()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"My DofRowMap Size="<<structdis_->DofRowMap()->NumMyElements()<<"  My DofColMap Size="<<structdis_->DofColMap()->NumMyElements()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Dis Structure NumColEles: "<<structdis_->NumMyColElements()<<" on PROC "<<Comm().MyPID()<<std::endl;
# endif

    curr_subset_of_fluiddis_ = fluid_SearchTree_->searchElementsInRadius(*fluiddis_,currpositions_fluid_,boundingboxcenter,structsearchradiusfac*max_radius,0);

    if(curr_subset_of_fluiddis_.empty() == false)
      std::cout<<"\nPrepareFluidOp returns "<<curr_subset_of_fluiddis_.begin()->second.size()<<" background elements on Proc "<<Comm().MyPID()<<std::endl;
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
//        std::cout<<"PROC"<<Comm().MyPID()<<" key "<<closele->first<<" eleIter "<<counter<<" : "<<*eleIter<<std::endl;
//      }
//    }
//  }
//  else
//    std::cout<<"curr_subset_of_fluiddis_ empty on PROC "<<Comm().MyPID()<<std::endl;
//
//  //___________________________
//
//  if(currpositions_struct_.empty() == false)
//  {
//    int counter = 0;
//    for(std::map<int,LINALG::Matrix<3,1> >::const_iterator closele = currpositions_struct_.begin(); closele != currpositions_struct_.end(); closele++)
//    {
//      counter ++;
//      if(Comm().MyPID()==1)
//      std::cout<<"PROC"<<Comm().MyPID()<<" key: "<<closele->first<<"  position: "<<counter<<closele->second(0)<<" "<<closele->second(1)<<" "<<closele->second(2)<<std::endl;
//    }
//  }

  return;
}

