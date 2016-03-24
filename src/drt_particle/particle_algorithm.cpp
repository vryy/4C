/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.cpp

\brief Algorithm to control particle simulations

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_adapter/ad_str_structure.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_inpar/inpar_meshfree.H"
#include "../drt_inpar/inpar_cavitation.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_gmsh.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
PARTICLE::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : AlgorithmBase(comm,params),
  BinningStrategy(comm),
  particles_(Teuchos::null),
  bincolmap_(Teuchos::null),
  structure_(Teuchos::null),
  particlewalldis_(Teuchos::null),
  moving_walls_(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"MOVING_WALLS"))
{
  const Teuchos::ParameterList& meshfreeparams = DRT::Problem::Instance()->MeshfreeParams();
  // safety check
  INPAR::MESHFREE::meshfreetype meshfreetype = DRT::INPUT::IntegralValue<INPAR::MESHFREE::meshfreetype>(meshfreeparams,"TYPE");
  if (meshfreetype!=INPAR::MESHFREE::particle)
    dserror("MESHFREE -> TYPE must be Particle in input file.");

  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  gravity_acc_.PutScalar(0.0);
  // get acceleration vector due to gravity for particles
  std::istringstream accstream(Teuchos::getNumericStringParameter(particleparams,"GRAVITY_ACCELERATION"));
  for(int dim=0; dim<3; dim++)
  {
    double value = 0.0;
    if(accstream >> value)
      gravity_acc_(dim) = value;
  }

  if(particle_dim_ == INPAR::PARTICLE::particle_2Dz)
  {
    gravity_acc_(2) = 0.0;
    if(myrank_ == 0)
      IO::cout << "gravity in z-direction ignored as this is a pseudo-2D problem" << IO::endl;
  }

  // initial setup of particle discretization
  particledis_ = DRT::Problem::Instance()->GetDis("particle");
  // new dofs are numbered from zero, minnodgid is ignored and it does not register in static_dofsets_
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset = Teuchos::rcp(new DRT::IndependentDofSet(true));
  particledis_->ReplaceDofSet(independentdofset);

  if(sparse_binning_ && DRT::Problem::Instance()->ProblemType() == prb_particle)
    dserror("sparse bin scheme cannot be used for pure particle problems");

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the particle algorithm                      ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // transfer particles into their correct bins
    if(Step()%100 == 0 and Comm().NumProc() != 1)
      DynamicLoadBalancing();

    // counter and print header
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // calculate stresses, strains, energies
    PrepareOutput();

    // transfer particles into their correct bins
    TransferParticles(true);

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Init(bool restarted)
{
  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);

  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));

  Teuchos::RCP<Epetra_Map> binrowmap;
  if(not restarted)
  {
    CreateBins(particledis_);
    binrowmap = DistributeBinsToProcs();
  }
  else
  {
    binrowmap = Teuchos::rcp(new Epetra_Map(*particledis_->ElementRowMap()));
  }

  // setup pbcs after bins have been created
  BuildParticlePeriodicBC();

  if(binrowmap->NumGlobalElements() > particlerowmap->NumGlobalElements() / 4.0)
    IO::cout << "\n\n\n WARNING: Reduction of number of bins recommended!! Increase cutoff radius. \n\n\n" << IO::endl;

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBinsRoundRobin(homelessparticles);

  // ghost bins and particles according to the bins --> final FillComplete() call included
  SetupGhosting(binrowmap);

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // add fully redundant discret for particle walls with identical dofs to full structural discret

    // get input lists
    const Teuchos::ParameterList& particledyn = DRT::Problem::Instance()->ParticleParams();
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    // access the structural discretization
    Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(particledyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
    structure_ = structure->StructureField();

    SetupParticleWalls(structdis);

    // assign wall elements to bins initially once for fixed walls (additionally rebuild pointers after ghosting)
    if(!moving_walls_)
      AssignWallElesToBins();

    // create time integrator based on structural time integration
    Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm(particledyn, particledis_));
    particles_ = particles->ParticleField();

    // set particle algorithm into time integration
    particles_->SetParticleAlgorithm(Teuchos::rcp(this,false));
    particles_->Init();

    // determine consistent initial acceleration for the particles
    CalculateAndApplyForcesToParticles();
    particles_->DetermineMassDampConsistAccel();
  }
  else
  {
    // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
    bool rebuildwallpointer = true;
    if(moving_walls_)
      rebuildwallpointer = false;
    BuildElementToBinPointers(rebuildwallpointer);
  }

  // some output
  if (myrank_ == 0)
    IO::cout << "after ghosting of particles" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  // apply dirichlet boundary conditions
  particles_->PrepareTimeStep();

  if(structure_ != Teuchos::null)
    structure_->PrepareTimeStep();

  // do rough safety check if bin size is appropriate
  BinSizeSafetyCheck(particles_->Dt());

  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Integrate()
{
  CalculateAndApplyForcesToParticles();

  Teuchos::RCP<Epetra_Vector> walldisn = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> walldisnp = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> wallvelnp = Teuchos::null;

  // solve for structural (wall) problem
  if(moving_walls_)
  {
    structure_->Solve();

    // extract displacement and velocity from full structural field to obtain wall states
    walldisn = wallextractor_->ExtractCondVector(structure_->Dispn());
    walldisnp = wallextractor_->ExtractCondVector(structure_->Dispnp());
    wallvelnp = wallextractor_->ExtractCondVector(structure_->Velnp());
  }
  else
  {
    walldisn = LINALG::CreateVector(*particlewalldis_->DofRowMap(), true);
    walldisnp = LINALG::CreateVector(*particlewalldis_->DofRowMap(), true);
    wallvelnp = LINALG::CreateVector(*particlewalldis_->DofRowMap(), true);
  }

  particlewalldis_->SetState("walldisn", walldisn);
  particlewalldis_->SetState("walldisnp", walldisnp);
  particlewalldis_->SetState("wallvelnp", wallvelnp);

  // assign wall elements dynamically to bins
  if(moving_walls_)
    AssignWallElesToBins();

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::Integrate");
  particles_->IntegrateStep();

  return;
}


/*----------------------------------------------------------------------*
 | calculate forces on particle and apply it               ghamm 02/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::CalculateAndApplyForcesToParticles()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::CalculateAndApplyForcesToParticles");

  // vector to be filled with forces
  Teuchos::RCP<Epetra_Vector> particleforces = LINALG::CreateVector(*particledis_->DofRowMap(),true);

  // mass of particles
  Teuchos::RCP<const Epetra_Vector> mass_p = particles_->Mass();

  // all row particles are evaluated
  const int numrownodes = particledis_->NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    /*------------------------------------------------------------------*/
    //// gravity forces = mass_p * g
    for(int dim=0; dim<3; ++dim)
    {
      (*particleforces)[i*3+dim] = (*mass_p)[i] * gravity_acc_(dim);
    }
    /*------------------------------------------------------------------*/
  }

  // apply forces to particles
  particles_->SetForceInterface(particleforces);

  return;
}


/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Update()
{
  if(structure_ != Teuchos::null)
    structure_->Update();

  // write state vectors from n+1 to n
  particles_->Update();

  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ReadRestart(int restart)
{
  // 1st) loop over bins and remove initial particle info
  const int numrowbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numrowbin; ++ibin)
  {
    DRT::Element* actele = particledis_->lColElement(ibin);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele)->DeleteNodes();
  }

  // 2nd) initial particles need to be removed from particledis_
  particledis_->DeleteNodes();

  // read in particles for restart
  {
    IO::DiscretizationReader reader(particledis_, restart);
    reader.ReadNodesOnly(restart);
  }

  // Init() is needed to obtain connectivity -> includes FillComplete())
  Init(true);

  // now, correct map layouts are available and states can be read
  particles_->ReadRestart(restart);
  SetTimeStep(particles_->TimeOld(),restart);

  // read restart for walls
  if(structure_ != Teuchos::null)
    structure_->ReadRestart(restart);

  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::Algorithm::DistributeBinsToProcs()
{
  // initial dummy distribution
  const int numbin = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];
  Teuchos::RCP<Epetra_Map> roweles = Teuchos::rcp(new Epetra_Map(numbin,0,Comm()));

  const int maxband = 26;
  Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*roweles,maxband,false));

  // fill all local entries into the graph
  {
    for (int lid=0; lid<roweles->NumMyElements(); ++lid)
    {
      const int binId = roweles->GID(lid);

      std::vector<int> neighbors;
      GetBinConnectivity(binId,neighbors);

      int err = graph->InsertGlobalIndices(binId,(int)neighbors.size(),&neighbors[0]);
      if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,binId);
    }
  }

  // finish graph
  graph->FillComplete();
  graph->OptimizeStorage();

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "PARTITION");

  Epetra_CrsGraph *balanced_graph = NULL;
  try {
    balanced_graph =
      Isorropia::Epetra::createBalancedCopy(*graph, paramlist);

  }
  catch(std::exception& exc) {
    std::cout << "Isorropia::createBalancedCopy threw "
         << "exception '" << exc.what() << "' on proc "
         << myrank_ << std::endl;
    dserror("Error within Isorropia (graph balancing)");
  }

  // obtain the row map
    Teuchos::RCP<Epetra_CrsGraph> rcp_balanced_graph = Teuchos::rcp(balanced_graph);
  rcp_balanced_graph->FillComplete();
  rcp_balanced_graph->OptimizeStorage();
  roweles = Teuchos::rcp(new Epetra_Map(-1,
      rcp_balanced_graph->RowMap().NumMyElements(),
      rcp_balanced_graph->RowMap().MyGlobalElements(),0,Comm()));

  // fill bins into discret
  for(int i=0; i<roweles->NumMyElements(); i++)
  {
    const int gid = roweles->GID(i);
    Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
    particledis_->AddElement(bin);
  }

  // return binrowmap
  return roweles;
}


/*----------------------------------------------------------------------*
| dynamic load balancing for bin distribution               ghamm 08/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::DynamicLoadBalancing()
{
  const Epetra_Map* oldrowmap = particledis_->ElementRowMap();

  Teuchos::RCP<const Epetra_CrsGraph> constgraph = CreateGraph();

  // Now we're going to create a Epetra_Vector with vertex weights to
  // be used in the repartitioning operation.
  Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*oldrowmap, false);
  // weights must be at least one for zoltan
  double* vals = vweights->Values();
  for(int i=0; i<oldrowmap->NumMyElements(); ++i)
  {
    const int numnode = particledis_->lRowElement(i)->NumNode();
    vals[i] = 1.0 + numnode*3 + numnode*numnode;
  }

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
   Teuchos::rcp(new Isorropia::Epetra::CostDescriber);

  costs->setVertexWeights(vweights);

  Teuchos::ParameterList paramlist;
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "REPARTITION");
  // Now create the partitioner object

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
   Teuchos::rcp(new Isorropia::Epetra::Partitioner(constgraph, costs, paramlist));

  // Extract repartitioned map
  Teuchos::RCP<Epetra_Map> newelerowmap = partitioner->createNewMap();

  //--------------------------------------------------------------------
  // rebuild of the system with the new map

  // export elements to new layout
  particledis_->ExportRowElements(*newelerowmap);

  // export row nodes to new layout
  {
    // create a set of row particle IDs for each proc
    std::set<int> particles;
    for (int lid=0; lid<newelerowmap->NumMyElements(); ++lid)
    {
      DRT::Element* bin = particledis_->gElement(newelerowmap->GID(lid));
      const int* particleids = bin->NodeIds();
      for(int iparticle=0; iparticle<bin->NumNode(); ++iparticle)
        particles.insert(particleids[iparticle]);
    }

    // copy particlegids to a vector and create particlerowmap
    std::vector<int> rowparticles(particles.begin(),particles.end());
    Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(-1,(int)rowparticles.size(),&rowparticles[0],0,Comm()));

    // place all nodes on the correct processor
    particledis_->ExportRowNodes(*particlerowmap);
  }

  // ghost bins and particles according to the bins --> final FillComplete() call included
  SetupGhosting(newelerowmap);

  BuildElementToBinPointers(true);

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  return;
}


/*----------------------------------------------------------------------*
| dynamic load balancing for bin distribution               ghamm 08/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_CrsGraph> PARTICLE::Algorithm::CreateGraph()
{
  const Epetra_Map* oldrowmap = particledis_->ElementRowMap();

  const int maxband = 26;
  Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*oldrowmap,maxband,false));

  // fill all local entries into the graph
  {
    for (int lid=0; lid<oldrowmap->NumMyElements(); ++lid)
    {
      const int binId = oldrowmap->GID(lid);

      std::vector<int> neighbors;
      GetBinConnectivity(binId,neighbors);

      int err = graph->InsertGlobalIndices(binId,(int)neighbors.size(),&neighbors[0]);
      if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,binId);
    }
  }

  return graph;
}


/*----------------------------------------------------------------------*
 | rough safety check for proper bin size                  ghamm 12/15  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BinSizeSafetyCheck(const double dt)
{
  // rough safety check whether bin size is large enough for proper contact detection
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  INPAR::PARTICLE::ContactStrategy contact_strategy =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");
  if(contact_strategy != INPAR::PARTICLE::None)
  {
    double extrema[2] = {0.0, 0.0};
    particles_->Veln()->MinValue(&extrema[0]);
    particles_->Veln()->MaxValue(&extrema[1]);
    const double maxvel = std::max(-extrema[0], extrema[1]);
    double maxrad = 0.0;
    particles_->Radius()->MaxValue(&maxrad);
    if(maxrad + maxvel*dt > 0.5*cutoff_radius_)
    {
      int lid = -1;
      // find entry with max velocity
      for(int i=0; i<particles_->Veln()->MyLength(); ++i)
      {
        if((*particles_->Veln())[i] < extrema[0]+1.0e-12 || (*particles_->Veln())[i] > extrema[1]-1.0e-12)
        {
          lid = i;
          break;
        }
      }

#ifdef DEBUG
      particles_->Veln()->Print(std::cout);
#endif
      if(lid != -1)
      {
        const int gid = particles_->Veln()->Map().GID(lid);
        dserror("Particle %i (gid) travels more than one bin per time step (%f > %f). Increase bin size or reduce step size."
            "Max radius is: %f", (int)gid/3, 2.0*(maxrad + maxvel*dt), cutoff_radius_, maxrad);
      }
    }
  }
}


/*----------------------------------------------------------------------*
| fill particles into their correct bin on according proc   ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::FillParticlesIntoBinsRoundRobin(std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>& homelessparticles)
{
  //--------------------------------------------------------------------
  // -> 2) round robin loop

  const int numproc = particledis_->Comm().NumProc();
  const int myrank = particledis_->Comm().MyPID();       // me
  const int torank = (myrank + 1) % numproc;             // to
  const int fromrank = (myrank + numproc - 1) % numproc; // from

  DRT::Exporter exporter(particledis_->Comm());

  for (int irobin = 0; irobin < numproc; ++irobin)
  {
    std::vector<char> sdata;
    std::vector<char> rdata;

    // ---- pack data for sending -----
    {
      DRT::PackBuffer data;
      for (std::set<Teuchos::RCP<DRT::Node> >::const_iterator currparticle=homelessparticles.begin(); currparticle != homelessparticles.end(); ++currparticle)
      {
//        cout << " Id:" << (*currparticle)->Id() << " was packed on proc: " << myrank_ << endl;
        (*currparticle)->Pack(data);
      }
      data.StartPacking();
      for (std::set<Teuchos::RCP<DRT::Node> >::const_iterator currparticle=homelessparticles.begin(); currparticle != homelessparticles.end(); ++currparticle)
      {
        (*currparticle)->Pack(data);
        particledis_->DeleteNode((*currparticle)->Id());
      }
      std::swap(sdata, data());
    }


    // ---- send ----
    MPI_Request request;
    exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);


    // ---- receive ----
    int length = rdata.size();
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234 or from != fromrank)
      dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);


    // ---- unpack ----
    {
      // Put received nodes either into discretization or into list of homeless particles
      homelessparticles.clear();
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::vector<char> data;
        DRT::ParObject::ExtractfromPack(index,rdata,data);
        // this Teuchos::rcp holds the memory of the node
        Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data),true);
        Teuchos::RCP<DRT::Node> node = Teuchos::rcp_dynamic_cast<DRT::Node>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        // process received particle
        const double* currpos = node->X();
        PlaceNodeCorrectly(node, currpos, homelessparticles);
      }
    }


    // wait for all communication to finish
    exporter.Wait(request);
    particledis_->Comm().Barrier(); // I feel better this way ;-)
  } // end for irobin

  if(homelessparticles.size())
  {
    std::cout << " There are " << homelessparticles.size() << " particles which have left the computational domain on rank " << myrank << std::endl;
    // erase everything that is left
    homelessparticles.clear();
  }

  return;
}


/*----------------------------------------------------------------------*
| fill particles into their correct bin on according proc   ghamm 03/16 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::FillParticlesIntoBinsRemoteIdList(std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>& homelessparticles)
{
  const int numproc = particledis_->Comm().NumProc();
  if(numproc == 1)
  {
    if(homelessparticles.size())
    {
      std::cout << " There are " << homelessparticles.size() << " particles which have left the computational domain on rank " << myrank_ << std::endl;
      std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>::const_iterator hlp;
      for(hlp = homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
      {
        particledis_->DeleteNode((*hlp)->Id());
      }
      homelessparticles.clear();
    }

    return;
  }

  // parallel case
  // ---- find new host procs for particles -----
  const int fullsize = (int)homelessparticles.size();
  std::vector<int> targetbinIdlist;
  targetbinIdlist.reserve(fullsize);
  std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>::const_iterator hlp;
  for(hlp=homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
  {
    const int binId = ConvertPosToGid((*hlp)->X());
    targetbinIdlist.push_back(binId);
  }

  // get proc which will be the future host of homeless particles
  std::vector<int> pidlist(fullsize);
  {
    // only unique id lists are accepted in RemoteIDList
    // 1) make gid list unique
    std::set<int> unique_targetbinIdlist(targetbinIdlist.begin(), targetbinIdlist.end());
    std::vector<int> uniquevec_targetbinIdlist(unique_targetbinIdlist.begin(), unique_targetbinIdlist.end());
    const int uniquesize = (int)unique_targetbinIdlist.size();

    // 2) communication
    std::vector<int> unique_pidlist(uniquesize);
    int err = particledis_->ElementRowMap()->RemoteIDList(uniquesize,uniquevec_targetbinIdlist.data(),unique_pidlist.data(),NULL);
    if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d",err);

    // 3) build full pid list via lookup table
    std::map<int,int> lookuptable;
    for(int s=0; s<uniquesize; ++s)
      lookuptable.insert(lookuptable.end(), std::pair<int,int>(uniquevec_targetbinIdlist[s], unique_pidlist[s]));
    for(int s=0; s<fullsize; ++s)
      pidlist[s] = lookuptable[targetbinIdlist[s]];
  }

  // ---- pack data for sending -----
  std::map<int, std::vector<char> > sdata;
  std::vector<int> targetprocs(numproc,0);
  int counter=0;
  int iter=0;
  for(hlp = homelessparticles.begin(); hlp != homelessparticles.end(); ++hlp)
  {
    Teuchos::RCP<DRT::Node> iterhomelessparticle= *hlp;

    // ---- pack data for sending -----
    const int targetproc = pidlist[iter];
    if(targetproc != -1)
    {
      DRT::PackBuffer data;
      iterhomelessparticle->Pack(data);
      data.StartPacking();
      iterhomelessparticle->Pack(data);
      particledis_->DeleteNode(iterhomelessparticle->Id());
      sdata[targetproc].insert(sdata[targetproc].end(),data().begin(),data().end());
      targetprocs[targetproc] = 1;
    }
    else
    {
      ++counter;
      particledis_->DeleteNode(iterhomelessparticle->Id());
    }
    ++iter;
  }
  if(counter)
    std::cout << " There are " << counter << " particles which have left the computational domain on rank " << myrank_ << std::endl;
  homelessparticles.clear();

  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc,0);
  particledis_->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);


  // ---- send ----
  DRT::Exporter exporter(particledis_->Comm());
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for(std::map<int, std::vector<char> >::const_iterator p=sdata.begin(); p!=sdata.end(); ++p)
  {
    exporter.ISend(myrank_, p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // ---- receive ----
  for(int rec=0; rec<summedtargets[myrank_]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // ---- unpack ----
    {
      // Put received nodes into discretization
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::vector<char> data;
        DRT::ParObject::ExtractfromPack(index,rdata,data);
        // this Teuchos::rcp holds the memory of the node
        Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data),true);
        Teuchos::RCP<DRT::Node> node = Teuchos::rcp_dynamic_cast<DRT::Node>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        // process received particle
        const double* currpos = node->X();
        PlaceNodeCorrectly(node, currpos, homelessparticles);
        if(homelessparticles.size())
          dserror("particle (id: %i) was sent to proc %i but corresponding bin (gid: %i) is missing", node->Id(), myrank_, ConvertPosToGid(currpos));
      }
    }
  }

  // wait for all communications to finish
  {
    for (int i=0; i<length; ++i)
      exporter.Wait(request[i]);
  }

  particledis_->Comm().Barrier(); // I feel better this way ;-)

  return;
}


/*----------------------------------------------------------------------*
| node is placed into the correct row bin                   ghamm 09/12 |
 *----------------------------------------------------------------------*/
bool PARTICLE::Algorithm::PlaceNodeCorrectly
(Teuchos::RCP<DRT::Node> node,
  const double* currpos,
  std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>& homelessparticles
  )
{
//  std::cout << "on proc: " << myrank_ << " node with ID: " << node->Id() << " and owner: " << node->Owner() << " arrived in PlaceNodeCorrectly" << std::endl;
  const int binId = ConvertPosToGid(currpos);

  // check whether the current node belongs into a bin on this proc
  const bool found = particledis_->HaveGlobalElement(binId);

  // either fill particle into correct bin on this proc or mark it as homeless
  if(found == true)
  {
    DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>( particledis_->gElement(binId) );
#ifdef DEBUG
    if(currbin == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    // check whether it is a row bin
    if(currbin->Owner() == myrank_) // row bin
    {
//      std::cout << "on proc: " << myrank_ << " for node " << node->Id() << " a row bin was found" << std::endl;
      // node already exists (either row or ghost)
      if( particledis_->HaveGlobalNode(node->Id()) == true)
      {
        DRT::Node* existingnode = particledis_->gNode(node->Id());
        // existing node is a row node, this means that node is equal existingnode
        if(existingnode->Owner() == myrank_)
        {
//          std::cout << "on proc: " << myrank_ << " existingnode row node " << existingnode->Id() << " (ID from outside node: " << node->Id() << ") is added to element: " << currbin->Id() << std::endl;

          // assign node to the correct bin
          currbin->AddNode(existingnode);
        }
        else // ghost node becomes row node and node from outside is trashed
        {
//          std::cout << "on proc: " << myrank_ << " existingnode ghost node " << existingnode->Id() << " (ID from outside node: " << node->Id() << ") is added to element: " << currbin->Id() << " after setting ownership" << std::endl;

          // change owner of the node to this proc
          existingnode->SetOwner(myrank_);

          // received node is no longer needed, X() of former ghost node has to be updated for output reasons
          {
            std::vector<double> update(3,0.0);
            const double* refposparticle = existingnode->X();
            for(int dim=0; dim<3; dim++)
              update[dim] = currpos[dim] - refposparticle[dim];
            // change X() of existing node
            existingnode->ChangePos(update);
          }

          // assign node to the correct bin
          currbin->AddNode(existingnode);
        }
      }
      else // fill newly received node into discretization
      {
        // change owner of the node to this proc and add it to the discretization
        node->SetOwner(myrank_);
        particledis_->AddNode(node);
//        std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to the discretization and assigned to element: " << currbin->Id() << std::endl;
        // assign node to the correct bin
        currbin->AddNode(node.get());
      }

      return true;
    }
    else // ghost bin
    {
      homelessparticles.insert(node);
//      std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to homeless because of becoming a future ghost node" << std::endl;
      return false;
    }
  }
  else // bin not found on this proc
  {
    homelessparticles.insert(node);
//    std::cout << "on proc: " << myrank_ << " node " << node->Id() << " is added to homeless because bin is not on this proc " << std::endl;
    return false;
  }

}


/*----------------------------------------------------------------------*
| setup ghosting of bins and particles                      ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap)
{
  // 1st step: ghosting of bins
  {
    // gather bins of rowmap and all its neighbors (row + ghost)
    std::set<int> bins;
    for (int lid=0;lid<binrowmap->NumMyElements();++lid)
    {
      const int gid = binrowmap->GID(lid);
      int ijk[3] = {-1,-1,-1};
      ConvertGidToijk(gid, ijk);

      // get all neighboring cells, including the element itself: one layer ghosting
      for(int i=-1;i<2;i++)
      {
        for(int j=-1;j<2;j++)
        {
          for(int k=-1;k<2;k++)
          {
            const int ijk_neighbor[3] = {ijk[0]+i, ijk[1]+j, ijk[2]+k};

            const int neighborgid = ConvertijkToGid(&ijk_neighbor[0]);
            if(neighborgid != -1)
            {
              bins.insert(neighborgid);
            }
          } // end for int k
        } // end for int j
      } // end for int i
    } // end for lid

    // remove non-existing ghost bins from original bin set
    if(sparse_binning_)
    {
      // create copy of column bins
      std::set<int> ghostbins(bins);
      // find ghost bins and check for existence
      for (int lid=0;lid<binrowmap->NumMyElements();++lid)
      {
        const int gid = binrowmap->GID(lid);
        std::set<int>::iterator iter = ghostbins.find(gid);
        if(iter != ghostbins.end())
          ghostbins.erase(iter);
      }
      // only ghost bins remain
      std::vector<int> ghostbins_vec(ghostbins.begin(),ghostbins.end());
      const int size = (int)ghostbins.size();
      std::vector<int> pidlist(size);
      const int err = binrowmap->RemoteIDList(size,ghostbins_vec.data(),pidlist.data(),NULL);
      if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d",err);

      for(int i=0; i<size; ++i)
      {
        if(pidlist[i] == -1)
        {
          std::set<int>::iterator iter = bins.find(ghostbins_vec[i]);
          if(iter == bins.end())
            dserror("bin id is missing in bin set");
          // erase non-existing id
          bins.erase(iter);
        }
      }
    }

    // copy bingids to a vector and create bincolmap
    std::vector<int> bincolmap(bins.begin(),bins.end());
    bincolmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)bincolmap.size(),&bincolmap[0],0,Comm()));

    if(bincolmap_->NumGlobalElements() == 1 && bincolmap_->Comm().NumProc() > 1)
      dserror("one bin cannot be run in parallel -> reduce CUTOFF_RADIUS");

    // create ghosting for bins (each knowing its particle ids)
    particledis_->ExtendedGhosting(*bincolmap_,true,false,true,false);
  }


#ifdef DEBUG
    // check whether each proc has only particles that are within bins on this proc
    for(int k=0; k<particledis_->NumMyColElements(); k++)
    {
      int binid = particledis_->lColElement(k)->Id();
      DRT::Node** particles = particledis_->lColElement(k)->Nodes();

      for(int iparticle=0; iparticle<particledis_->lColElement(k)->NumNode(); iparticle++)
      {
        int ijk[3] = {-1,-1,-1};
        for(int dim=0; dim<3; ++dim)
        {
          ijk[dim] = (int)((particles[iparticle]->X()[dim]-XAABB_(dim,0)) / bin_size_[dim]);
        }

        int gidofbin = ConvertijkToGid(&ijk[0]);
        if(gidofbin != binid)
          dserror("after ghosting: particle which should be in bin no. %i is in %i",gidofbin,binid);
      }
    }
#endif


  return;
}


/*----------------------------------------------------------------------*
 | particles are checked and transferred if necessary       ghamm 10/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::TransferParticles(const bool updatestates, const bool ghosting)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles");

  // leave here in case nothing to do
  if(particles_->Radius()->GlobalLength() == 0)
    return;

  // set of homeless particles
  std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less> homelessparticles;

  // current positions of particles
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();

  // apply periodic boundary conditions for particles
  if(havepbc_)
  {
    // offset delta for pbc direction
    // delta equal zero otherwise
    std::vector<double> delta(3);
    for(int dim=0; dim<3; ++dim)
      delta[dim] = pbcbounds_[dim][1] - pbcbounds_[dim][0];

    for(int i=0; i<disnp->MyLength(); i++)
    {
      const int dim = i%3;
      if((*disnp)[i] < pbcbounds_[dim][0])
      {
        (*disnp)[i] += delta[dim];
         continue;
      }
      if((*disnp)[i] > pbcbounds_[dim][1])
        (*disnp)[i] -= delta[dim];
    }
  }

  std::set<int> examinedbins;
  // check in each bin whether particles have moved out
  // first run over particles and then process whole bin in which particle is located
  // until all particles have been checked
  const int numrownodes = particledis_->NodeRowMap()->NumMyElements();
  for(int i=0; i<numrownodes; ++i)
  {
    DRT::Node *currparticle = particledis_->lRowNode(i);

#ifdef DEBUG
    if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
#endif

    DRT::Element** currele = currparticle->Elements();
    DRT::Element* currbin = currele[0];
    // as checked above, there is only one element in currele array
    const int binId = currbin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.count(binId) == 1 )
    {
      continue;
    }
    // else: bin is examined for the first time --> new entry in examinedbins
    else
    {
      examinedbins.insert(binId);
    }

    // now all particles in this bin are processed
#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::Node** particles = currbin->Nodes();
    std::vector<int> tobemoved(0);
    for(int iparticle=0; iparticle<currbin->NumNode(); iparticle++)
    {
      DRT::Node* currnode = particles[iparticle];
      // get the first gid of a node and convert it into a LID
      const int gid = particledis_->Dof(currnode, 0);
      const int lid = disnp->Map().LID(gid);

      double currpos[3];
      for(int dim=0; dim<3; ++dim)
        currpos[dim] = (*disnp)[lid+dim];

      // update reference configuration of particle for correct output and correct placement via MPI
      {
        std::vector<double> update(3,0.0);
        const double* refposparticle = currnode->X();
        for(int dim=0; dim<3; ++dim)
        {
          update[dim] = currpos[dim] - refposparticle[dim];
        }
        // change X() of current particle
        currnode->ChangePos(update);
//        std::cout << "particle (Id: " << currnode->Id() << " ) position is updated to" << currnode->X()[0] << "  "<< currnode->X()[1] << "  "<< currnode->X()[2] << "  " << std::endl;
      }

      const int gidofbin = ConvertPosToGid(currpos);
      if(gidofbin != binId) // particle has left current bin
      {
        // gather all node Ids that will be removed and remove them afterwards
        // (looping over nodes and deleting at the same time is detrimental)
        tobemoved.push_back(currnode->Id());
        // find new bin for particle
        /*bool placed = */PlaceNodeCorrectly(Teuchos::rcp(currnode,false), currpos, homelessparticles);
      }

    } // end for inode

    // finally remove nodes from their old bin
    for(size_t iter=0; iter<tobemoved.size(); iter++)
    {
      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currbin)->DeleteNode(tobemoved[iter]);
    }

  } // end for ibin

#ifdef DEBUG
  if(homelessparticles.size())
    std::cout << "There are " << homelessparticles.size() << " homeless particles on proc" << myrank_ << std::endl;
#endif

  // homeless particles are sent to their new processors where they are inserted into their correct bin
  FillParticlesIntoBinsRemoteIdList(homelessparticles);

  // check whether all procs have a filled particledis_,
  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
  particledis_->CheckFilledGlobally();

  // new ghosting if necessary
  if (ghosting)
    particledis_->ExtendedGhosting(*bincolmap_,true,false,true,false);
  else
    particledis_->FillComplete(true, false, true);

  // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
  bool rebuildwallpointer = true;
  if(moving_walls_)
    rebuildwallpointer = false;
  BuildElementToBinPointers(rebuildwallpointer);

  // update state vectors in time integrator to the new layout
  if(updatestates)
  {
    particles_->UpdateStatesAfterParticleTransfer();
    UpdateStates();
  }

  return;
}


/*----------------------------------------------------------------------*
| particle walls are added from the structural discret      ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupParticleWalls(Teuchos::RCP<DRT::Discretization> basediscret)
{
  //--------------------------------------------------------------------
  // 1st step: build fully redundant discretization with wall elements
  //--------------------------------------------------------------------

  // declare struct objects in wall condition
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > > structgelements; // col map of structure elements
  std::map<int, std::map<int, DRT::Node*> > dummy2;  // dummy map
  std::map<int, std::map<int, DRT::Node*> > structgnodes; // col map of structure nodes

  //initialize struct objects in wall condition
  DRT::UTILS::FindConditionObjects(*basediscret, dummy2, structgnodes, structgelements,"ParticleWall");

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit;

  // initialize new particle wall discretizations
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(basediscret->Comm().Clone());
  const std::string discret_name = "particlewalls";
  Teuchos::RCP<DRT::Discretization> particlewalldis = Teuchos::rcp(new DRT::Discretization(discret_name,com));

  // number of dofs is important for transparent dof set
  // only zeros are applied to the wall displacements when fluid domain is basediscret
  // -> number of dofs is irrelevant when reading data for wall discret in this case
  // future implementation using ALE needs to be handled like a structure
  std::stringstream elename;
  if(structure_ != Teuchos::null)
    elename << "BELE3_" << 3;
  else
    elename << "BELE3_" << 4;

  std::vector<int> nodeids;
  std::vector<int> eleids;
  // loop over all particle wall nodes and elements and fill new discretization
  for(std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit=structgelements.begin(); meit!=structgelements.end(); ++meit)
  {
    // care about particle wall nodes (only row nodes in case of moving walls, otherwise everything)
    std::map<int, DRT::Node*> wallgnodes = structgnodes[meit->first];
    for (std::map<int, DRT::Node* >::iterator nit=wallgnodes.begin(); nit != wallgnodes.end(); ++nit)
    {
      DRT::Node* currnode = (*nit).second;
      if ( currnode->Owner() == myrank_ || sparse_binning_ )
      {
        nodeids.push_back(currnode->Id());
        particlewalldis->AddNode(Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
      }
    }

    // care about particle wall eles (only row elements in case of moving walls, otherwise everything)
    std::map<int, Teuchos::RCP<DRT::Element> > structelementsinterf = structgelements[meit->first];
    for (std::map<int, Teuchos::RCP<DRT::Element> >::iterator eit=structelementsinterf.begin(); eit != structelementsinterf.end(); ++eit)
    {
      Teuchos::RCP<DRT::Element> currele = eit->second;
      if ( currele->Owner() == myrank_ || sparse_binning_ )
      {
        eleids.push_back(currele->Id() );
        // structural surface elements cannot be distributed --> Bele3 element is used
        Teuchos::RCP<DRT::Element> wallele = DRT::UTILS::Factory(elename.str(),"Polynomial", currele->Id(), currele->Owner());
        wallele->SetNodeIds(currele->NumNode(), currele->NodeIds());
        particlewalldis->AddElement( wallele );
      }
    }
  }

  // extended ghosting of wall elements for sparse bin scheme
  if(sparse_binning_)
  {
    // fill complete in order to obtain element col map
    particlewalldis->FillComplete(false, false, false);

    std::map<int, std::set<int> > rowelesinbin;
    DistributeElesToBins(particlewalldis, rowelesinbin);

    // get extended column map for wall elements
    std::map<int, std::set<int> > dummy;
     Teuchos::RCP<Epetra_Map> wallelecolmap =
         ExtendGhosting(particlewalldis->ElementColMap(), rowelesinbin, dummy, bincolmap_);

     // extend ghosting (add nodes/elements) according to the new column layout
     particlewalldis->ExtendedGhosting(*wallelecolmap,false, false, false, true);
  }
  else  // ... or otherwise do fully redundant storage of wall elements in case of moving boundaries
  {
    // row node map of walls
    Teuchos::RCP<Epetra_Map> wallnoderowmap = Teuchos::rcp(new Epetra_Map(-1,nodeids.size(),&nodeids[0],0,particlewalldis->Comm()));

    // fully overlapping node map
    Teuchos::RCP<Epetra_Map> wallrednodecolmap = LINALG::AllreduceEMap(*wallnoderowmap);

    // row ele map of walls
    Teuchos::RCP<Epetra_Map> wallelerowmap = Teuchos::rcp(new Epetra_Map(-1,eleids.size(),&eleids[0],0,particlewalldis->Comm()));
    // fully overlapping ele map
    Teuchos::RCP<Epetra_Map> wallredelecolmap = LINALG::AllreduceEMap(*wallelerowmap);

    // do the fully overlapping ghosting of the wall elements to have everything redundant
    particlewalldis->ExportColumnNodes(*wallrednodecolmap);
    particlewalldis->ExportColumnElements(*wallredelecolmap);
  }

  // find out if we are in parallel; needed for TransparentDofSet
  bool parallel = (particlewalldis->Comm().NumProc() == 1) ? false : true;

  // dofs of the original discretization are used to set same dofs for the new particle wall discretization
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentDofSet(basediscret,parallel));
  particlewalldis->ReplaceDofSet(newdofset);
  newdofset=Teuchos::null;

  // final fill complete to reorganize everything in the discretization
  particlewalldis->FillComplete(true, false, false);
  particlewalldis_ = particlewalldis;

  // some output to screen and initialization of binary output
  if(myrank_ == 0)
    std::cout << "after adding particle walls" << std::endl;
  DRT::UTILS::PrintParallelDistribution(*particlewalldis_);
  if(moving_walls_)
    wallextractor_ = Teuchos::rcp(new LINALG::MapExtractor(*(structure_->Discretization()->DofRowMap()),Teuchos::rcp(particlewalldis_->DofRowMap(), false)));

  return;
}


/*----------------------------------------------------------------------*
| particle walls are added from the structural discret      ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::AssignWallElesToBins()
{
  // loop over all bins and remove assigned wall elements
  const int numcolbins = particledis_->ElementColMap()->NumMyElements();
  for(int binlid=0; binlid<numcolbins; ++binlid)
  {
    DRT::Element *currentbin = particledis_->lColElement(binlid);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currentbin)->RemoveAssociatedWallEles();
  }

  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  if(moving_walls_)
  {
    Teuchos::RCP<const Epetra_Vector> walldisn = particlewalldis_->GetState("walldisn");
    for (int lid=0; lid<particlewalldis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = particlewalldis_->lColNode(lid);
      std::vector<int> lm_node;
      lm_node.reserve(3);
      particlewalldis_->Dof(node,lm_node);

      // nodal displacements
      static LINALG::Matrix<3,1> node_disn;
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*walldisn, node_disn, lm_node);

      LINALG::Matrix<3,1> currpos;
      const double* X = node->X();
      for(int dim=0; dim<3; ++dim)
        currpos(dim) = X[dim] + node_disn(dim);
      currentpositions[node->Id()] = currpos;
    }
  }
  else
  {
    for(int lid=0; lid<particlewalldis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = particlewalldis_->lColNode(lid);

      LINALG::Matrix<3,1> currpos;
      const double* X = node->X();
      for(int dim=0; dim<3; ++dim)
        currpos(dim) = X[dim];
      currentpositions[node->Id()] = currpos;
    }
  }

  double bincircumcircle = 0.0;
  for(int dim=0; dim<3; ++dim)
  {
    bincircumcircle += std::pow(bin_size_[dim]/2.0,2.0);
  }
  bincircumcircle = sqrt(bincircumcircle);

  // minimal bin size
  double min_bin_size = bin_size_[0];
  for(int dim=1; dim<3; ++dim)
    min_bin_size = std::min(min_bin_size, bin_size_[dim]);

  // find bins for all wall elements
  const int numcolwalleles = particlewalldis_->NumMyColElements();
  for(int lid=0; lid<numcolwalleles; ++lid)
  {
    DRT::Element* wallele = particlewalldis_->lColElement(lid);
    const int *nodeids = wallele->NodeIds();
    const int numnode = wallele->NumNode();
    // variable to store bin ids in which this wall element is located
    std::set<int> binIds;

    // do a positive search and get all bins enclosed in the bounding box of each wall element
    {
      // initialize ijk_range with ijk of first node of wall element
      int ijk[3];
      ConvertPosToijk(currentpositions[nodeids[0]], ijk);

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; ++j)
      {
        int ijk[3];
        ConvertPosToijk(currentpositions[nodeids[j]], ijk);

        for(int dim=0; dim<3; ++dim)
        {
          if(ijk[dim]<ijk_range[dim*2])
            ijk_range[dim*2] = ijk[dim];
          if(ijk[dim]>ijk_range[dim*2+1])
            ijk_range[dim*2+1] = ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range and fill them into binIds
      GidsInijkRange(&ijk_range[0], binIds, true);
    }

    // if no bins on this proc were found, next wall element can be processed
    if(binIds.empty())
      continue;

    // do a negative search and remove bins that are too far away from the wall element
    {
      std::set<int> binfaraway;
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
      {
        const LINALG::Matrix<3,1> bincentroid = GetBinCentroid(*biniter);

        // search for the closest object, more exactly it's coordinates
        LINALG::Matrix<3,1> minDistCoords;
        GEO::nearest3DObjectOnElement(wallele, currentpositions, bincentroid, minDistCoords);

        LINALG::Matrix<3,1> distance;
        distance.Update(1.0, bincentroid, -1.0, minDistCoords);
        double dist = distance.Norm2();

        // if distance is larger than radius of circumcircle of bin --> too far away
        if(dist > bincircumcircle)
        {
          binfaraway.insert(*biniter);
        }
        // if distance is smaller than half the minimum bin size --> very close
        else if(dist <= min_bin_size*0.5)
        {
          continue;
        }
        // if distance is between half the minimum bin size and radius of the circumcircle --> further checks
        else
        {
          std::vector<LINALG::Matrix<3,1> > bincorners;
          GetBinCorners(*biniter, bincorners);

          // in case wall element is axis aligned, it might not be detected as inside because projection points
          // are located on the edges of the bin --> Remedy: bin centroid is tested as well
          bincorners.push_back(bincentroid);

          bool projpointinsidebin = false;
          // all corners of the close bin are projected onto the wall element: if at least one projection
          // point is inside the bin, it won't be removed from the list
          for(size_t corner=0; corner<bincorners.size(); ++corner)
          {
            // search for the closest object, more exactly it's coordinates
            LINALG::Matrix<3,1> minDistCoords;
            GEO::nearest3DObjectOnElement(wallele, currentpositions, bincorners[corner], minDistCoords);

            const int gid = ConvertPosToGid(minDistCoords);
            if(gid == *biniter)
            {
              projpointinsidebin = true;
              break;
            }
          }
          if(projpointinsidebin == false)
            binfaraway.insert(*biniter);
        }
      }
      for(std::set<int>::const_iterator biniter=binfaraway.begin(); biniter!=binfaraway.end(); ++biniter)
        binIds.erase(*biniter);
    }

    // assign wall element to remaining bins
    {
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->gElement(*biniter))->AddAssociatedWallEle(wallele->Id(), wallele);
    }

  } // end lid

  return;
}


/*----------------------------------------------------------------------*
| build connectivity from particle wall elements to bins    ghamm 04/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BuildElementToBinPointers(bool wallpointer)
{
  if(wallpointer == true)
  {
    // loop over column bins and fill wall elements
    const int numcolbin = particledis_->NumMyColElements();
    for (int ibin=0; ibin<numcolbin; ++ibin)
    {
      DRT::Element* actele = particledis_->lColElement(ibin);
      DRT::MESHFREE::MeshfreeMultiBin* actbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
      const int numwallele = actbin->NumAssociatedWallEle();
      const int* walleleids = actbin->AssociatedWallEleIds();
      std::vector<DRT::Element*> wallelements(numwallele);
      for(int iwall=0; iwall<numwallele; ++iwall)
      {
        const int wallid = walleleids[iwall];
        wallelements[iwall] = particlewalldis_->gElement(wallid);
      }
      actbin->BuildWallElePointers(&wallelements[0]);
    }
  }
  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors based on an        ghamm 11/12 |
| underlying discretization                                             |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::Algorithm::DistributeBinsToProcsBasedOnUnderlyingDiscret(
  Teuchos::RCP<DRT::Discretization> underlyingdis,
  std::map<int, std::set<int> >& rowelesinbin)
{

  //--------------------------------------------------------------------
  // 1st step: exploiting bounding box idea for scatra elements and bins
  //--------------------------------------------------------------------

  DistributeElesToBins(underlyingdis, rowelesinbin);

  //--------------------------------------------------------------------
  // 2nd step: decide which proc will be owner of each bin
  //--------------------------------------------------------------------

  std::vector<int> rowbins;
  {
    // NOTE: This part of the setup can be the bottleneck because vectors of all bins
    // are needed on each proc (memory issue!!); std::map could perhaps help when gathering
    // num fluid nodes in each bin, then block wise communication after copying data to vector

    const int numbins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];
    std::vector<int> mynumeles_per_bin(numbins,0);

    std::map<int, std::set<int> >::const_iterator iter;
    for(iter=rowelesinbin.begin(); iter!=rowelesinbin.end(); ++ iter)
    {
      mynumeles_per_bin[iter->first] = iter->second.size();
    }

    // find maximum number of eles in each bin over all procs (init with -1)
    std::vector<int> maxnumeles_per_bin(numbins,-1);
    underlyingdis->Comm().MaxAll(&mynumeles_per_bin[0], &maxnumeles_per_bin[0], numbins);

    // it is possible that several procs have the same number of eles in a bin
    // only proc which has maximum number of eles in a bin writes its rank
    std::vector<int> myrank_per_bin(numbins,-1);
    for(int i=0; i<numbins; ++i)
    {
      if(mynumeles_per_bin[i] == maxnumeles_per_bin[i])
        myrank_per_bin[i] = myrank_;
    }

    mynumeles_per_bin.clear();
    if(sparse_binning_ == false)
      maxnumeles_per_bin.clear();

    // find maximum myrank for each bin over all procs (init with -1)
    std::vector<int> maxmyrank_per_bin(numbins,-1);
    underlyingdis->Comm().MaxAll(&myrank_per_bin[0], &maxmyrank_per_bin[0], numbins);

    // distribute bins to proc with highest rank
    for(int gid=0; gid<numbins; ++gid)
    {
      if(myrank_ == maxmyrank_per_bin[gid])
      {
        // do not add empty bins in case of sparse bin distribution
        if(sparse_binning_ && maxnumeles_per_bin[gid] == 0)
          continue;

        Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
        particledis_->AddElement(bin);
        rowbins.push_back(gid);
      }
    }

    maxnumeles_per_bin.clear();
    myrank_per_bin.clear();
    maxmyrank_per_bin.clear();
  }

  // return binrowmap (without having called FillComplete on particledis_ so far)
  return Teuchos::rcp(new Epetra_Map(-1,(int)rowbins.size(),&rowbins[0],0,Comm()));
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(particles_->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
  return;
}


/*----------------------------------------------------------------------*
 | calculate stresses, strains, energies                   ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::PrepareOutput()
{
  particles_->PrepareOutput();
  if(structure_ != Teuchos::null)
    structure_->PrepareOutput();

  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Output(bool forced_writerestart /*= false*/)
{
  // INFO regarding output: Bins are not written to file because they cannot
  // be post-processed anyway (no nodes and connectivity available)
  particles_->OutputStep(forced_writerestart);
  if(structure_ != Teuchos::null)
  {
    structure_->Output();
    // add missing restart information if necessary
    if(forced_writerestart)
      structure_->Output(forced_writerestart);
  }

//  const std::string filename = IO::GMSH::GetFileName("particle_data", Step(), true, Comm().MyPID());
//  std::ofstream gmshfilecontent(filename.c_str());
//
//  // velocity
//  {
//    gmshfilecontent << "View \" " << "velocity" << " \" {\n";
//    LINALG::Matrix<3,1> vectorvalue(true);
//
//    for(int n=0; n<particledis_->NumMyRowNodes(); n++)
//    {
//      DRT::Node* actnode = particledis_->lRowNode(n);
//      // get the first gid of a node and convert it into a LID
//      int gid = particledis_->Dof(actnode, 0);
//      int lid = particles_->Dispnp()->Map().LID(gid);
//      Teuchos::RCP<const Epetra_Vector> disnp = particles_->Dispnp();
//      Teuchos::RCP<const Epetra_Vector> velnp = particles_->Velnp();
//      LINALG::Matrix<3,1> posXYZDomain(true);
//      for(int dim=0; dim < 3; dim++)
//      {
//        posXYZDomain(dim) = (*disnp)[lid+dim];
//        vectorvalue(dim) = (*velnp)[lid+dim];
//      }
//
//      // write data to Gmsh file
//      IO::GMSH::VectorToStream(posXYZDomain, vectorvalue, gmshfilecontent);
//    }
//
//    gmshfilecontent << "};\n";
//  }
//
//  // density
//  {
//    gmshfilecontent << "View \" " << "density" << " \" {\n";
//
//    for(int n=0; n<particledis_->NumMyRowNodes(); n++)
//    {
//      DRT::Node* actnode = particledis_->lRowNode(n);
//      // get the first gid of a node and convert it into a LID
//      int gid = particledis_->Dof(actnode, 0);
//      int lid = particles_->Dispnp()->Map().LID(gid);
//      Teuchos::RCP<const Epetra_Vector> disnp = particles_->Dispnp();
//      LINALG::Matrix<3,1> posXYZDomain(true);
//      for (int dim=0; dim<3; dim++)
//      {
//        posXYZDomain(dim) = (*disnp)[lid+dim];
//      }
//
//      const double density = particles_->ParticleDensity();
//
//      // write data to Gmsh file
//      IO::GMSH::ScalarToStream(posXYZDomain, density, gmshfilecontent);
//    }
//
//    gmshfilecontent << "};\n";
//  }
//
//  // radius
//  {
//    gmshfilecontent << "View \" " << "radius" << " \" {\n";
//
//    for(int n=0; n<particledis_->NumMyRowNodes(); n++)
//    {
//      DRT::Node* actnode = particledis_->lRowNode(n);
//      // get the first gid of a node and convert it into a LID
//      int gid = particledis_->Dof(actnode, 0);
//      int lid = particles_->Dispnp()->Map().LID(gid);
//      Teuchos::RCP<const Epetra_Vector> disnp = particles_->Dispnp();
//      LINALG::Matrix<3,1> posXYZDomain(true);
//      for (int dim=0; dim<3; dim++)
//      {
//        posXYZDomain(dim) = (*disnp)[lid+dim];
//      }
//
//      double radius = (*particles_->Radius())[n];
//
//      // write data to Gmsh file
//      IO::GMSH::ScalarToStream(posXYZDomain, radius, gmshfilecontent);
//    }
//
//    gmshfilecontent << "};\n";
//  }
//
//  gmshfilecontent.close();

  return;
}
