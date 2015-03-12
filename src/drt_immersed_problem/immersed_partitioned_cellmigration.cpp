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
#include "../drt_adapter/ad_str_fsiwrapper.H"


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

  // build field structure
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> cellstructure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(globalproblem_->CellMigrationParams(), const_cast<Teuchos::ParameterList&>(globalproblem_->StructuralDynamicParams()), immerseddis_));
    cellstructure_ = Teuchos::rcp_dynamic_cast< ::ADAPTER::FSIStructureWrapper>(cellstructure->StructureField());
    std::cout<<"Created Field Cell Sturucture ..."<<std::endl;

  // create instance of poroelast subproblem
  poroelast_subproblem_ = Teuchos::rcp(new POROELAST::Monolithic(comm, globalproblem_->PoroelastDynamicParams()));
  std::cout<<"Created Field Poroelast ..."<<std::endl;

  // construct 3D search tree for fluid domain
  // initialized in SetupSBackgroundDiscretization()
  fluid_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // construct 3D search tree for cell domain
  // initialized in SetupImmersedDiscretization()
  cell_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // get coupling variable
  displacementcoupling_ = DRT::Problem::Instance()->ImmersedMethodParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE") == "Displacement";
  if(displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Force "<<std::endl;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::CouplingOp()
{
  std::cout<<"Reached CouplingOP() ... "<<std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::BackgroundOp()
{
  std::cout<<"Reached BackgroundOp() ... "<<std::endl;
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::ImmersedOp()
{
  std::cout<<"Reached ImmersedOp() ... "<<std::endl;
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::InitialGuess()
{
  std::cout<<"Reached InitialGuess() ... "<<std::endl;
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::BuildImmersedDirichMap()
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::DoImmersedDirichletCond()
{
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

  // find positions of the immersed structural discretization
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_struct;
  for (int lid = 0; lid < immerseddis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = immerseddis_->lRowNode(lid);
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
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_struct,currpositions_cell_,numproc_,&procs[0],Comm());

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
  std::cout<<"Reached PrepareImmersedOp() ... "<<std::endl;
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareOutput()
{
  std::cout<<"Reached PrepareOutput() ... "<<std::endl;
  return;
}
