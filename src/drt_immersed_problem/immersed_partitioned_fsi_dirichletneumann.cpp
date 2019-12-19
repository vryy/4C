/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed fsi algorithm for neumann-(dirichlet-neumann) like coupling

\level 2

\maintainer Jonas Eichinger

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

#include "../linalg/linalg_utils_densematrix_communication.H"

// search tree related
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"

// time monitoring
#include <Teuchos_TimeMonitor.hpp>


IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ImmersedPartitionedFSIDirichletNeumann(
    const Epetra_Comm& comm)
    : ImmersedBase(),
      FSI::PartitionedImmersed(comm),
      struct_bdry_traction_(Teuchos::null),
      fluid_artificial_velocity_(Teuchos::null),
      dbcmap_immersed_(Teuchos::null),
      fluid_SearchTree_(Teuchos::null),
      structure_SearchTree_(Teuchos::null),
      myrank_(comm.MyPID()),
      numproc_(comm.NumProc()),
      globalproblem_(NULL),
      displacementcoupling_(false),
      multibodysimulation_(false),
      output_evry_nlniter_(false),
      is_relaxation_(false),
      isALE_(false),
      correct_boundary_velocities_(0),
      degree_gp_fluid_bound_(0),
      artificial_velocity_isvalid_(false),
      boundary_traction_isvalid_(false),
      immersed_info_isvalid_(false),
      fluiddis_(Teuchos::null),
      structdis_(Teuchos::null),
      immersedstructure_(Teuchos::null)

{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedFSIDirichletNeumann::Init(const Teuchos::ParameterList& params)
{
  // reset the setup flag
  SetIsSetup(false);

  // do all init stuff here

  // set isinit_ flag true
  SetIsInit(true);

  return 0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ReadRestart(int step)
{
  FSI::Partitioned::ReadRestart(step);

  SetupStructuralDiscretization();

  if (not displacementcoupling_)
  {
    CalcArtificialVelocity();
    CalcFluidTractionsOnStructure();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // call setup of base class
  FSI::PartitionedImmersed::Setup();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // get pointer to discretizations
  fluiddis_ = globalproblem_->GetDis("fluid");
  structdis_ = globalproblem_->GetDis("structure");

  // cast to specialized adapter
  immersedstructure_ =
      Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(StructureField());

  // get coupling variable
  displacementcoupling_ = globalproblem_->FSIDynamicParams()
                              .sublist("PARTITIONED SOLVER")
                              .get<std::string>("COUPVARIABLE") == "Displacement";
  if (displacementcoupling_ and myrank_ == 0)
    std::cout << "\n Coupling variable for partitioned FSI scheme :  Displacements " << std::endl;
  else if (!displacementcoupling_ and myrank_ == 0)
    std::cout << "\n Coupling variable for partitioned FSI scheme :  Force " << std::endl;

  // set switch for interface velocity correction
  correct_boundary_velocities_ = (DRT::INPUT::IntegralValue<int>(
      globalproblem_->ImmersedMethodParams(), "CORRECT_BOUNDARY_VELOCITIES"));

  // set switch for output in every nln. iteration (for debugging)
  output_evry_nlniter_ = (DRT::INPUT::IntegralValue<int>(
      globalproblem_->ImmersedMethodParams(), "OUTPUT_EVRY_NLNITER"));

  // print acceleration method
  if (globalproblem_->FSIDynamicParams().get<std::string>("COUPALGO") ==
      "iter_stagg_fixed_rel_param")
  {
    is_relaxation_ = false;
    if (myrank_ == 0) std::cout << "\n Using FIXED relaxation parameter. " << std::endl;
  }
  else if (globalproblem_->FSIDynamicParams().get<std::string>("COUPALGO") ==
           "iter_stagg_AITKEN_rel_param")
  {
    is_relaxation_ = true;
    if (myrank_ == 0) std::cout << "\n Using AITKEN relaxation parameter. " << std::endl;
  }
  else
    dserror("Unknown definition of COUPALGO in FSI DYNAMIC section for Immersed FSI.");

  // check for unfeasible combination
  if (correct_boundary_velocities_ and displacementcoupling_ and is_relaxation_)
    dserror(
        "Interface velocity correction is not possible with displacement coupled Immersed FSI in "
        "combination with relaxation.");

  // get integration rule for fluid elements cut by structural boundary
  int num_gp_fluid_bound = globalproblem_->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  if (num_gp_fluid_bound == 8)
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
    dserror(
        "Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, "
        "729, and 1000).");

  // Decide whether multiple structural bodies or not.
  // Bodies need to be labeled with "ImmersedSearchbox" condition.
  std::vector<DRT::Condition*> conditions;
  structdis_->GetCondition("ImmersedSearchbox", conditions);
  if ((int)conditions.size() > 0)
  {
    if (myrank_ == 0)
      std::cout << " MULTIBODY SIMULATION   Number of bodies: " << (int)conditions.size()
                << std::endl;
    multibodysimulation_ = true;
  }
  else
    multibodysimulation_ = false;

  // vector of fluid stresses interpolated to structural bdry. int. points and integrated over
  // structural surface
  struct_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(immersedstructure_->DofRowMap()), true));

  // vector with fluid velocities interpolated from structure
  fluid_artificial_velocity_ =
      Teuchos::rcp(new Epetra_Vector(*(MBFluidField()->FluidField()->DofRowMap()), true));

  // build 3D search tree for fluid domain
  fluid_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // find positions of the background fluid discretization
  for (int lid = 0; lid < fluiddis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = fluiddis_->lColNode(lid);
    LINALG::Matrix<3, 1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    currpositions_fluid_[node->Id()] = currpos;
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3, 2> rootBox = GEO::getXAABBofPositions(currpositions_fluid_);
  fluid_SearchTree_->initializeTree(rootBox, *fluiddis_, GEO::TreeType(GEO::OCTTREE));

  if (myrank_ == 0) std::cout << "\n Build Fluid SearchTree ... " << std::endl;

  // construct 3D search tree for structural domain
  // initialized in SetupStructuralDiscretization()
  structure_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // Validation flag for velocity in artificial domain. After each structure solve the velocity
  // becomes invalid and needs to be projected again.
  artificial_velocity_isvalid_ = false;

  // Validation flag for bdry. traction on structure. After each fluid solve the traction becomes
  // invalid and needs to be integrated again.
  boundary_traction_isvalid_ = false;

  // Validation flag for immersed information. After each structure solve, we have to assess
  // again, which fluid elements are covered by the structure and which fluid eles are cut by
  // the interface.
  // todo: NOTE: There is little inconsistency right now in this method.
  //             Fluid elements are labeled isboundarimmersed_ if:
  //             1) at least one but not all nodes are covered, or
  //             2) if a gp of the structural surface lies in a fluid element.
  //
  //             The first criterion is checked in CalcArtificialVelocity(). This is
  //             done at the same time as struct vel. projection is performed.
  //             The second criterion is checked in CalcFluidTractionsOnStructure().
  //             This is done at the same time as the bdry. traction is integrated.
  //
  //             Since, the fluid field decides which gps of the fluid are compressible,
  //             based on IsImmersed() and IsBoundaryImmersed() information, it might
  //             happen, that after performing CalcArtificialVelocity(), the nodal criterion
  //             is updated correctly, but since the structure has moved, the gp criterion
  //             is invalid. It is only updated after the next struct solve. So the fluid
  //             might be solved with incorrect compressible gps.
  //
  //             To fix this, one would have to split the information update and the projections.
  //
  //             However:
  //             1) This would be more expensive.
  //             2) I assume the error we make by the procedure explained above is very small, since
  //                the deformations between iteration steps should be small. Especially, when we
  //                are close to convergence, the difference in structure deformation should be so
  //                small, that no new gps should have changed their states inbetween the last two
  //                structure solves.
  immersed_info_isvalid_ = false;

  // wait for all processors to arrive here
  Comm().Barrier();

  // set flag issetup true
  SetIsSetup(true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FSIOp(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  CheckIsInit();
  CheckIsSetup();

  if (displacementcoupling_)
  {
    // get the current artificial velocity state
    const Teuchos::RCP<Epetra_Vector> artificial_velocity_n = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL FluidOp
    ////////////////////
    FluidOp(artificial_velocity_n, fillFlag);  //!< solve the fluid

    ////////////////////
    // CALL StructOp
    ////////////////////
    CalcFluidTractionsOnStructure();  //!< calculate new fluid tractions interpolated to structural
                                      //!< surface
    StructOp(immersedstructure_->Interface()->ExtractIMMERSEDCondVector(struct_bdry_traction_),
        fillFlag);               //!< solve the structure
    ResetImmersedInformation();  //!< structure moved; immersed info are invalid -> reset
    const Teuchos::RCP<Epetra_Vector> artificial_velocity_np =
        CalcArtificialVelocity();  //!< calc new projected velocities and update immersed
                                   //!< information

    int err = CalcResidual(F, artificial_velocity_np, artificial_velocity_n);
    if (err != 0) dserror("Vector update of FSI-residual returned err=%d", err);
  }
  else if (!displacementcoupling_)  // FORCE COUPLING
  {
    // get the current interface force state
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL StructOp
    ////////////////////
    StructOp(iforcen, fillFlag);  //!< solve the structure
    ResetImmersedInformation();   //!< structure moved; immersed info are invalid -> reset

    ////////////////////
    // CALL FluidOp
    ////////////////////
    CalcArtificialVelocity();  //!< calc the new velocity in the artificial fluid domain, immersed
                               //!< info are set inside
    FluidOp(fluid_artificial_velocity_, fillFlag);  //!< solve the fluid
    CalcFluidTractionsOnStructure();  //!< calculate new fluid tractions integrated over structural
                                      //!< surface

    int err = CalcResidual(F, struct_bdry_traction_, iforcen);
    if (err != 0) dserror("Vector update of FSI-residual returned err=%d", err);
  }

  // write output after every solve of fluid and structure
  // current limitations:
  // max 100 partitioned iterations and max 100 timesteps in total
  if (output_evry_nlniter_)
  {
    int iter = ((FSI::Partitioned::IterationCounter())[0]);
    Teuchos::rcp_dynamic_cast<ADAPTER::FluidAleImmersed>(MBFluidField())
        ->Output((Step() * 100) + (iter - 1), Time() - Dt() * ((100 - iter) / 100.0));
    StructureField()->PrepareOutput();
    Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(StructureField())
        ->Output(false, (Step() * 100) + (iter - 1), Time() - Dt() * ((100 - iter) / 100.0));
  }

  // perform n steps max; then set converged
  bool nlnsolver_continue =
      globalproblem_->ImmersedMethodParams().get<std::string>("DIVERCONT") == "continue";
  int itemax = globalproblem_->FSIDynamicParams().sublist("PARTITIONED SOLVER").get<int>("ITEMAX");
  if ((FSI::Partitioned::IterationCounter())[0] == itemax and nlnsolver_continue)
  {
    if (myrank_ == 0)
      std::cout << "\n  Continue with next time step after ITEMAX = "
                << (FSI::Partitioned::IterationCounter())[0] << " iterations. \n"
                << std::endl;

    // !!! EXPERIMENTAL !!!
    // set F to zero to tell NOX that this timestep is converged
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(F.Map(), true));
    F.Update(1.0, *zeros, 0.0);
    // !!! EXPERIMENTAL !!!

    // clear states after time step was set converged
    immersedstructure_->Discretization()->ClearState();
    MBFluidField()->Discretization()->ClearState();
  }

  if (globalproblem_->ImmersedMethodParams().get<std::string>("TIMESTATS") == "everyiter")
  {
    Teuchos::TimeMonitor::summarize();
    Teuchos::TimeMonitor::zeroOutTimers();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumann::FluidOp(
    Teuchos::RCP<Epetra_Vector> fluid_artificial_velocity, const FillType fillFlag)
{
  // print
  FSI::Partitioned::FluidOp(fluid_artificial_velocity, fillFlag);

  if (fillFlag == User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    // get maximum number of Newton iterations
    const int itemax = MBFluidField()->Itemax();

    // apply the given artificial velocity to the fluid field
    ApplyImmersedDirichlet(fluid_artificial_velocity);

    // solve fluid
    SolveFluid();

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration)
    RemoveDirichCond();

    // correct the quality of the interface solution
    CorrectInterfaceVelocity();

    // set max number of Newton iterations
    MBFluidField()->SetItemax(itemax);

    // we just invalidated the boundary tractions
    boundary_traction_isvalid_ = false;

  }  // fillflag is not User

  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumann::StructOp(
    Teuchos::RCP<Epetra_Vector> struct_bdry_traction, const FillType fillFlag)
{
  FSI::Partitioned::StructOp(struct_bdry_traction, fillFlag);

  if (fillFlag == User)
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

    // structure moved; we just invalidated the artificial velocity
    artificial_velocity_isvalid_ = false;

    // we also invalidated the immersed info
    immersed_info_isvalid_ = false;


    return ExtractInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFSIDirichletNeumann::InitialGuess()
{
  if (myrank_ == 0) std::cout << "\n Do Initial Guess." << std::endl;

  if (displacementcoupling_)
    return CalcArtificialVelocity();
  else
    return immersedstructure_->Interface()->ExtractIMMERSEDCondVector(struct_bdry_traction_);

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::BuildImmersedDirichMap(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& dirichmap,
    const Teuchos::RCP<const Epetra_Map>& dirichmap_original)
{
  const Epetra_Map* elecolmap = dis->ElementColMap();
  std::vector<int> mydirichdofs(0);

  for (int i = 0; i < elecolmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    DRT::ELEMENTS::FluidImmersedBase* immersedele =
        dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(dis->gElement(elecolmap->GID(i)));
    if (immersedele->HasProjectedDirichlet())
    {
      DRT::Node** nodes = immersedele->Nodes();
      for (int inode = 0; inode < (immersedele->NumNode()); inode++)
      {
        if (static_cast<IMMERSED::ImmersedNode*>(nodes[inode])->IsMatched() and
            nodes[inode]->Owner() == myrank_)
        {
          std::vector<int> dofs = dis->Dof(nodes[inode]);

          for (int dim = 0; dim < 3; ++dim)
          {
            // if not already in original dirich map
            if (dirichmap_original->LID(dofs[dim]) == -1) mydirichdofs.push_back(dofs[dim]);
          }
          // include also pressure dof if node does not belong to a boundary background element
          // if((nodes[inode]->IsBoundaryImmersed())==0)
          //  mydirichdofs.push_back(dofs[3]);
        }
      }
    }
  }

  int nummydirichvals = mydirichdofs.size();
  dirichmap = Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichdofs[0]), 0, dis->Comm()));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::DoImmersedDirichletCond(
    Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals,
    Teuchos::RCP<Epetra_Map> dbcmap)
{
  int mynumvals = dbcmap->NumMyElements();
  double* myvals = dirichvals->Values();

  for (int i = 0; i < mynumvals; ++i)
  {
    int gid = dbcmap->GID(i);

#ifdef DEBUG
    int err = -2;
    int lid = dirichvals->Map().LID(gid);
    err = statevector->ReplaceGlobalValue(gid, 0, myvals[lid]);
    if (err == -1)
      dserror("VectorIndex >= NumVectors()");
    else if (err == 1)
      dserror("GlobalRow not associated with calling processor");
    else if (err != -1 and err != 1 and err != 0)
      dserror("Trouble using ReplaceGlobalValue on fluid state vector. ErrorCode = %d", err);
#else
    int lid = dirichvals->Map().LID(gid);
    statevector->ReplaceGlobalValue(gid, 0, myvals[lid]);
#endif
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetupStructuralDiscretization()
{
  // find positions of the immersed structural discretization
  std::map<int, LINALG::Matrix<3, 1>> my_currpositions_struct;
  for (int lid = 0; lid < structdis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = structdis_->lRowNode(lid);
    LINALG::Matrix<3, 1> currpos;

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
  for (int i = 0; i < numproc_; i++) procs[i] = i;
  LINALG::Gather<int, LINALG::Matrix<3, 1>>(
      my_currpositions_struct, currpositions_struct_, numproc_, &procs[0], Comm());

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3, 2> rootBox2 = GEO::getXAABBofDis(*structdis_, currpositions_struct_);
  structure_SearchTree_->initializeTree(rootBox2, *structdis_, GEO::TreeType(GEO::OCTTREE));

  if (myrank_ == 0) std::cout << "\n Build Structure SearchTree ... " << std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetStatesFluidOP()
{
  // for FluidOP
  structdis_->SetState(0, "displacement", immersedstructure_->Dispnp());
  structdis_->SetState(0, "velocity", immersedstructure_->Velnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetStatesVelocityCorrection()
{
  structdis_->SetState(0, "displacement", immersedstructure_->Dispnp());
  structdis_->SetState(0, "velocity", immersedstructure_->Velnp());
  fluiddis_->SetState(0, "velnp", MBFluidField()->FluidField()->Velnp());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SetStatesStructOP()
{
  structdis_->SetState(0, "displacement", immersedstructure_->Dispnp());
  fluiddis_->SetState(0, "velnp",
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed>(MBFluidField())->FluidField()->Velnp());
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::SolveFluid()
{
  MBFluidField()->NonlinearSolve(Teuchos::null, Teuchos::null);

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

  // search radius factor around center of structure bounding box (fac*diagonal of bounding box)
  const double structsearchradiusfac =
      globalproblem_->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

  // determine subset of fluid discretization which is potentially underlying the structural
  // discretization
  //
  // get state
  Teuchos::RCP<const Epetra_Vector> displacements = immersedstructure_->Dispnp();

  // find current positions for immersed structural discretization
  std::map<int, LINALG::Matrix<3, 1>> my_currpositions_struct;
  for (int lid = 0; lid < structdis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = structdis_->lRowNode(lid);
    LINALG::Matrix<3, 1> currpos;
    std::vector<int> dofstoextract(3);
    std::vector<double> mydisp(3);

    // get the current displacement
    structdis_->Dof(node, 0, dofstoextract);
    DRT::UTILS::ExtractMyValues(*displacements, mydisp, dofstoextract);

    currpos(0) = node->X()[0] + mydisp.at(0);
    currpos(1) = node->X()[1] + mydisp.at(1);
    currpos(2) = node->X()[2] + mydisp.at(2);

    my_currpositions_struct[node->Id()] = currpos;
  }

  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[numproc_];
  for (int i = 0; i < numproc_; i++) procs[i] = i;
  LINALG::Gather<int, LINALG::Matrix<3, 1>>(
      my_currpositions_struct, currpositions_struct_, numproc_, &procs[0], Comm());

  // take special care in case of multibody simulations
  if (multibodysimulation_ == false)
  {
    // get bounding box of current configuration of structural dis
    const LINALG::Matrix<3, 2> structBox = GEO::getXAABBofDis(*structdis_, currpositions_struct_);
    double max_radius =
        sqrt(pow(structBox(0, 0) - structBox(0, 1), 2) + pow(structBox(1, 0) - structBox(1, 1), 2) +
             pow(structBox(2, 0) - structBox(2, 1), 2));
    // search for background elements within a certain radius around the center of the immersed
    // bounding box
    LINALG::Matrix<3, 1> boundingboxcenter;
    boundingboxcenter(0) = structBox(0, 0) + (structBox(0, 1) - structBox(0, 0)) * 0.5;
    boundingboxcenter(1) = structBox(1, 0) + (structBox(1, 1) - structBox(1, 0)) * 0.5;
    boundingboxcenter(2) = structBox(2, 0) + (structBox(2, 1) - structBox(2, 0)) * 0.5;

#ifdef DEBUG
    std::cout << "Bounding Box of Structure: " << structBox << " on PROC " << myrank_ << std::endl;
    std::cout << "Bounding Box Center of Structure: " << boundingboxcenter << " on PROC " << myrank_
              << std::endl;
    std::cout << "Search Radius Around Center: " << structsearchradiusfac * max_radius
              << " on PROC " << myrank_ << std::endl;
    std::cout << "Length of Dispnp()=" << displacements->MyLength() << " on PROC " << myrank_
              << std::endl;
    std::cout << "Size of currpositions_struct_=" << currpositions_struct_.size() << " on PROC "
              << myrank_ << std::endl;
    std::cout << "My DofRowMap Size=" << structdis_->DofRowMap()->NumMyElements()
              << "  My DofColMap Size=" << structdis_->DofColMap()->NumMyElements() << " on PROC "
              << myrank_ << std::endl;
    std::cout << "Dis Structure NumColEles: " << structdis_->NumMyColElements() << " on PROC "
              << myrank_ << std::endl;
#endif

    if (isALE_) UpdateCurrentPositionsFluidNodes();

    SearchPotentiallyCoveredBackgrdElements(&curr_subset_of_fluiddis_, fluid_SearchTree_,
        *fluiddis_, currpositions_fluid_, boundingboxcenter, structsearchradiusfac * max_radius, 0);

    if (curr_subset_of_fluiddis_.empty() == false)
      std::cout << "\nPrepareFluidOp returns " << curr_subset_of_fluiddis_.begin()->second.size()
                << " background elements on Proc " << myrank_ << std::endl;
  }
  else
  {
    // get searchbox conditions on bodies
    std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
    std::vector<DRT::Condition*> conditions;
    structdis_->GetCondition("ImmersedSearchbox", conditions);

    // build element list
    std::map<int, std::set<int>> elementList;
    std::set<int> settoinsert;
    for (int i = 0; i < (int)conditions.size(); ++i)
    {
      for (curr = conditions[i]->Geometry().begin(); curr != conditions[i]->Geometry().end();
           ++curr)
      {
        settoinsert.insert(curr->second->Id());
      }
      elementList.insert(std::pair<int, std::set<int>>(i, settoinsert));
      settoinsert.clear();
    }

    // get bounding boxes of the bodies
    std::vector<LINALG::Matrix<3, 2>> structboxes =
        GEO::computeXAABBForLabeledStructures(*structdis_, currpositions_struct_, elementList);

    double max_radius;

    // search for background elements within a certain radius around the center of the immersed
    // bounding box
    for (int i = 0; i < (int)structboxes.size(); ++i)
    {
      max_radius = sqrt(pow(structboxes[i](0, 0) - structboxes[i](0, 1), 2) +
                        pow(structboxes[i](1, 0) - structboxes[i](1, 1), 2) +
                        pow(structboxes[i](2, 0) - structboxes[i](2, 1), 2));

      LINALG::Matrix<3, 1> boundingboxcenter;
      boundingboxcenter(0) =
          structboxes[i](0, 0) + (structboxes[i](0, 1) - structboxes[i](0, 0)) * 0.5;
      boundingboxcenter(1) =
          structboxes[i](1, 0) + (structboxes[i](1, 1) - structboxes[i](1, 0)) * 0.5;
      boundingboxcenter(2) =
          structboxes[i](2, 0) + (structboxes[i](2, 1) - structboxes[i](2, 0)) * 0.5;

      std::map<int, std::set<int>> tempmap = fluid_SearchTree_->searchElementsInRadius(
          *fluiddis_, currpositions_fluid_, boundingboxcenter, 0.5 * max_radius, 0);
      curr_subset_of_fluiddis_.insert(std::pair<int, std::set<int>>(i, (tempmap.begin()->second)));
    }

    for (int i = 0; i < (int)structboxes.size(); ++i)
      std::cout << "\nPrepareFluidOp returns " << curr_subset_of_fluiddis_.at(i).size()
                << " background elements for body " << i << std::endl;
  }  //

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ExtractInterfaceDispnp()
{
  return immersedstructure_->ExtractImmersedInterfaceDispnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ApplyInterfaceForces(
    Teuchos::RCP<Epetra_Vector> full_traction_vec)
{
  double normorstructbdrytraction;
  full_traction_vec->Norm2(&normorstructbdrytraction);
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Norm of Boundary Traction:   " << normorstructbdrytraction << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }
  immersedstructure_->ApplyImmersedInterfaceForces(Teuchos::null, full_traction_vec);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::AddDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed>(MBFluidField())
      ->AddDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::RemoveDirichCond()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed>(MBFluidField())
      ->RemoveDirichCond(dbcmap_immersed_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CalcResidual(Epetra_Vector& F,
    const Teuchos::RCP<Epetra_Vector> newstate, const Teuchos::RCP<Epetra_Vector> oldstate)
{
  int err = -1234;

  if (!displacementcoupling_)
    err = F.Update(1.0, *(immersedstructure_->Interface()->ExtractIMMERSEDCondVector(newstate)),
        -1.0, *oldstate, 0.0);
  else
    err = F.Update(1.0, *newstate, -1.0, *oldstate, 0.0);

  return err;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CalcFluidTractionsOnStructure()
{
  // sanity check
  if (boundary_traction_isvalid_)
    dserror(
        "Boundary traction from fluid onto immersed structure is still valid!\n"
        "If you really need to calc them anew, invalidate flag boundary_traction_isvalid_ at the "
        "proper position.");

  // reinitialize the transfer vector
  struct_bdry_traction_->Scale(0.0);

  // declare and fill parameter list
  Teuchos::ParameterList params;
  params.set<std::string>("action", "calc_fluid_traction");
  params.set<std::string>("backgrddisname", "fluid");
  params.set<std::string>("immerseddisname", "structure");

  // set the states needed for evaluation
  SetStatesStructOP();

  DRT::AssembleStrategy struct_bdry_strategy(0,  // struct dofset for row
      0,                                         // struct dofset for column
      Teuchos::null,                             // matrix 1
      Teuchos::null,                             //
      struct_bdry_traction_,                     // vector 1
      Teuchos::null,                             //
      Teuchos::null                              //
  );
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Interpolate fluid stresses to structural surface and calculate traction "
                 "...              "
              << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }
  EvaluateInterpolationCondition(
      immersedstructure_->Discretization(), params, struct_bdry_strategy, "IMMERSEDCoupling", -1);

  // we just validate the boundary tractions
  boundary_traction_isvalid_ = true;

  // we just validated the immersed info again.
  // technically this is not entirely true.
  // see remark in constructor of this class.
  // here we additionally validated the IsBoundaryImmersed
  // information based on the struct. bdry. int. points.
  immersed_info_isvalid_ = true;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CalcArtificialVelocity()
{
  if (not artificial_velocity_isvalid_)
  {
    // reinitialize the transfer vector
    fluid_artificial_velocity_->Scale(0.0);

    // declare parameter list
    Teuchos::ParameterList params;

    // provide number of integration points in fluid elements cut by boundary
    params.set<int>("intpoints_fluid_bound", degree_gp_fluid_bound_);
    // provide name of immersed discretization
    params.set<std::string>("immerseddisname", "structure");

    // set the states needed for evaluation
    SetStatesFluidOP();
    // update search trees, etc. ...
    PrepareFluidOp();

    // calc the fluid velocity from the structural displacements
    DRT::AssembleStrategy fluid_vol_strategy(0,  // struct dofset for row
        0,                                       // struct dofset for column
        Teuchos::null,                           // matrix 1
        Teuchos::null,                           //
        fluid_artificial_velocity_,              // vector 1
        Teuchos::null,                           //
        Teuchos::null                            //
    );

    if (myrank_ == 0)
    {
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Interpolate Dirichlet Values from immersed elements which overlap the "
                << MBFluidField()->Discretization()->Name() << " nodes ..." << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
    }

    EvaluateImmersed(params, MBFluidField()->Discretization(), &fluid_vol_strategy,
        &curr_subset_of_fluiddis_, structure_SearchTree_, &currpositions_struct_,
        FLD::interpolate_velocity_to_given_point_immersed, false);

    // we just validated the artificial velocity
    artificial_velocity_isvalid_ = true;

    // we just validated the immersed info again.
    // technically this is not entirely true.
    // see remark in constructor of this class.
    immersed_info_isvalid_ = true;
  }

  return fluid_artificial_velocity_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ApplyImmersedDirichlet(
    Teuchos::RCP<Epetra_Vector> artificial_velocity)
{
  BuildImmersedDirichMap(MBFluidField()->Discretization(), dbcmap_immersed_,
      MBFluidField()->FluidField()->GetDBCMapExtractor()->CondMap());
  AddDirichCond();

  // apply immersed dirichlets
  DoImmersedDirichletCond(
      MBFluidField()->FluidField()->WriteAccessVelnp(), artificial_velocity, dbcmap_immersed_);
  double normofvelocities = -1234.0;
  MBFluidField()->FluidField()->ExtractVelocityPart(artificial_velocity)->Norm2(&normofvelocities);

  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Norm of Dirichlet Values:   " << std::setprecision(7) << normofvelocities
              << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::CorrectInterfaceVelocity()
{
  //***********************************************************************************
  // Correct velocity at nodes of fluid elements being cut by the structural surface
  // (fluid is solved a second time with different Dirichlet values)
  //***********************************************************************************
  if (correct_boundary_velocities_)
  {
    // declare parameter list
    Teuchos::ParameterList params;

    // provide number of integration points in fluid elements cut by boundary
    params.set<int>("intpoints_fluid_bound", degree_gp_fluid_bound_);

    // calc the fluid velocity from the structural displacements
    DRT::AssembleStrategy fluid_vol_strategy(0,  // struct dofset for row
        0,                                       // struct dofset for column
        Teuchos::null,                           // matrix 1
        Teuchos::null,                           //
        fluid_artificial_velocity_,              // vector 1
        Teuchos::null,                           //
        Teuchos::null                            //
    );

    SetStatesVelocityCorrection();

    if (myrank_ == 0)
    {
      std::cout << "\nCorrection step " << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Correct Velocity in fluid boundary elements " << std::endl;
    }

    // calculate new dirichlet velocities for fluid elements cut by structure
    EvaluateImmersed(params, MBFluidField()->Discretization(), &fluid_vol_strategy,
        &curr_subset_of_fluiddis_, structure_SearchTree_, &currpositions_struct_,
        (int)FLD::correct_immersed_fluid_bound_vel, true);

    // Build new dirich map
    BuildImmersedDirichMap(MBFluidField()->Discretization(), dbcmap_immersed_,
        MBFluidField()->FluidField()->GetDBCMapExtractor()->CondMap());
    AddDirichCond();

    // apply new dirichlets after velocity correction
    DoImmersedDirichletCond(MBFluidField()->FluidField()->WriteAccessVelnp(),
        fluid_artificial_velocity_, dbcmap_immersed_);
    double normofnewvelocities;
    MBFluidField()
        ->FluidField()
        ->ExtractVelocityPart(fluid_artificial_velocity_)
        ->Norm2(&normofnewvelocities);

    if (myrank_ == 0)
    {
      std::cout << "###   Norm of new Dirichlet Values:   " << std::setprecision(7)
                << normofnewvelocities << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
    }

    // solve fluid again with new dirichlet values
    SolveFluid();

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration step)
    RemoveDirichCond();

  }  // correct_boundary_velocities_ finished

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumann::ResetImmersedInformation()
{
  if (immersed_info_isvalid_)
    dserror(
        "Immersed information are valid! Reconsider your call to ResetImmersedInformation().\n"
        "Did you forget to invalidate the flag immersed_info_isvalid_?");

  if (myrank_ == 0) std::cout << "\nReset Immersed Information ...\n" << std::endl;

  // reset element and node information about immersed method
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::reset_immersed_ele);
  params.set<int>("intpoints_fluid_bound", degree_gp_fluid_bound_);
  EvaluateSubsetElements(params, fluiddis_, curr_subset_of_fluiddis_, (int)FLD::reset_immersed_ele);

  return;
}
