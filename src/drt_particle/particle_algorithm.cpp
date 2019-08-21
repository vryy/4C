/*----------------------------------------------------------------------*/
/*! \file

\brief Algorithm to control particle simulations

\level 2

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"
#include "particle_timint_strategy.H"
#include "particle_node.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include "../drt_binstrategy/drt_meshfree_multibin.H"

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
PARTICLE::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : AlgorithmBase(comm, params),
      ParticleHandler(comm),
      particles_(Teuchos::null),
      dofimporter_(Teuchos::null),
      nodeimporter_(Teuchos::null),
      writeresultsevery_(0),
      particlewalldis_(Teuchos::null),
      particlewallelecolmap_standardghosting_(Teuchos::null),
      moving_walls_((bool)DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->ParticleParamsOld(), "MOVING_WALLS")),
      rep_strategy_(DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::RepartitionStrategy>(
          DRT::Problem::Instance()->ParticleParamsOld(), "REPARTITIONSTRATEGY")),
      dis_at_last_redistr_(Teuchos::null),
      particleInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::ParticleInteractions>(
          DRT::Problem::Instance()->ParticleParamsOld(), "PARTICLE_INTERACTION")),
      extendedGhosting_(DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::ExtendedGhosting>(
          DRT::Problem::Instance()->ParticleParamsOld(), "EXTENDED_GHOSTING")),
      particleMat_(0),
      bin_wallcontent_(BINSTRATEGY::UTILS::BELE3)
{
  // get particle params list
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParamsOld();

  if (moving_walls_) dserror("moving walls do not supported in old particle framework");

  gravity_acc_.PutScalar(0.0);
  // get acceleration vector due to gravity for particles
  std::istringstream accstream(
      Teuchos::getNumericStringParameter(particleparams, "GRAVITY_ACCELERATION"));
  for (int dim = 0; dim < 3; dim++)
  {
    double value = 0.0;
    if (accstream >> value) gravity_acc_(dim) = value;
  }

  if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dz)
  {
    gravity_acc_(2) = 0.0;
    if (MyRank() == 0)
      IO::cout << "gravity in z-direction ignored as this is a pseudo-2D problem" << IO::endl;
  }

  // initial setup of particle discretization
  BinStrategy()->BinDiscret() = DRT::Problem::Instance()->GetDis("particle");
  // new dofs are numbered from zero, minnodgid is ignored and it does not register in
  // static_dofsets_
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset =
      Teuchos::rcp(new DRT::IndependentDofSet(true));
  BinStrategy()->BinDiscret()->ReplaceDofSet(independentdofset);

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
    // redistribute load in parallel
    DynamicLoadBalancing();

    // counter and print header
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

  return;
}

/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupSystem() { return; }

/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Init(bool restarted)
{
  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  BinStrategy()->BinDiscret()->FillComplete(false, false, false);

  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap =
      Teuchos::rcp(new Epetra_Map(*BinStrategy()->BinDiscret()->NodeRowMap()));

  Teuchos::RCP<Epetra_Map> binrowmap;
  if (not restarted)
  {
    BinStrategy()->CreateBinsBasedOnCutoffAndXAABB(BinStrategy()->BinDiscret());
    // setup pbcs after bins have been created
    BinStrategy()->BuildPeriodicBC();
    binrowmap = DistributeBinsToProcs();
  }
  else
  {
    // setup pbcs after bins have been created
    BinStrategy()->BuildPeriodicBC();
    binrowmap = Teuchos::rcp(new Epetra_Map(*BinStrategy()->BinDiscret()->ElementRowMap()));
  }

  if (binrowmap->NumGlobalElements() > particlerowmap->NumGlobalElements() / 4.0 && MyRank() == 0)
    IO::cout << "\n\n\n CAREFUL: Reduction of number of bins recommended! Performance might be "
                "deteriorated. Increase cutoff radius. \n\n\n"
             << IO::endl;

  //--------------------------------------------------------------------
  // -> 1) create a list of homeless particles that are not in a bin on this proc
  std::list<Teuchos::RCP<DRT::Node>> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = BinStrategy()->BinDiscret()->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node, false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBinsRoundRobin(homelessparticles);

  // ghost bins and particles according to the bins --> final FillComplete() call included
  SetupGhosting(binrowmap);

  // determine boundary bins (physical boundary as well as boundary to other procs)
  BinStrategy()->DetermineBoundaryRowBins();

  // determine one layer ghosting around boundary bins determined in previous step
  BinStrategy()->DetermineBoundaryColBinsIds();

  // setup importer for dof and node based vectors
  if (particleInteractionType_ != INPAR::PARTICLEOLD::None) SetupImporter();

  // the following has only to be done once --> skip in case of restart
  if (not restarted)
  {
    // set up the links to the materials for easy access
    // make sure that a particle material is defined in the dat-file
    InitMaterials();

    // get input parameters for particles
    const Teuchos::ParameterList& particledyn = DRT::Problem::Instance()->ParticleParamsOld();

    // access structure and build particle walls
    AccessStructure();

    // create time integrator based on structural time integration
    Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm(particledyn, BinStrategy()->BinDiscret()));
    particles_ = particles->ParticleField();

    writeresultsevery_ = particledyn.get<int>("RESULTSEVRY");

    // set particle algorithm into time integration
    particles_->SetParticleAlgorithm(Teuchos::rcp(this, false));

    particles_->Init();

    // in case random noise is added to the particle position, particle transfer is necessary
    double amplitude = particledyn.get<double>("RANDOM_AMPLITUDE");
    if (amplitude)
    {
      SetParticleNodePos();
      TransferParticles(true, true);
    }

    // determine consistent initial acceleration for the particles
    particles_->UpdateExtActions();

    particles_->DetermineMassDampConsistAccel();
  }
  else
  {
    // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
    BuildElementToBinPointers(not moving_walls_);
  }

  // store current displacement state as displacements of last redistribution
  if (rep_strategy_ == INPAR::PARTICLEOLD::repstr_adaptive)
    dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*particles_->Dispn()));

  // some output
  if (MyRank() == 0) IO::cout << "after ghosting of particles" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*BinStrategy()->BinDiscret());

  return;
}

/*----------------------------------------------------------------------*
 | build connectivity from particle wall elements to bins   ghamm 04/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BuildElementToBinPointers(bool wallpointer)
{
  if (wallpointer == true)
  {
    // loop over column bins and fill wall elements
    const int numcolbin = BinStrategy()->BinDiscret()->NumMyColElements();
    for (int ibin = 0; ibin < numcolbin; ++ibin)
    {
      DRT::Element* actele = BinStrategy()->BinDiscret()->lColElement(ibin);
      DRT::MESHFREE::MeshfreeMultiBin* actbin =
          dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
      const int numwallele = actbin->NumAssociatedEle(bin_wallcontent_);
      const int* walleleids = actbin->AssociatedEleIds(bin_wallcontent_);
      std::vector<DRT::Element*> wallelements(numwallele);
      for (int iwall = 0; iwall < numwallele; ++iwall)
      {
        const int wallid = walleleids[iwall];
        wallelements[iwall] = particlewalldis_->gElement(wallid);
      }
      actbin->BuildElePointers(bin_wallcontent_, &wallelements[0]);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | assign wall elements and gids to bins                   sfuchs 08/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::AssignWallElesAndGidsToBins()
{
  // loop over wall elements
  for (std::map<int, std::set<int>>::iterator eleiter = relwallgidtobinids_.begin();
       eleiter != relwallgidtobinids_.end(); ++eleiter)
  {
    // global id of current wall element
    int wallgid = eleiter->first;

    // bins related to current wall element
    std::set<int> binIds = eleiter->second;

    // assign wall element to bins
    for (std::set<int>::const_iterator biniter = binIds.begin(); biniter != binIds.end(); ++biniter)
    {
      if (BinStrategy()->BinDiscret()->HaveGlobalElement(*biniter))
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(
            BinStrategy()->BinDiscret()->gElement(*biniter))
            ->AddAssociatedEle(bin_wallcontent_, wallgid, particlewalldis_->gElement(wallgid));
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | set up pointers to material bundles                     catta 09/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::InitMaterials()
{
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);

  if (id < 0)
    id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids);
  else if (DRT::Problem::Instance()->Materials()->FirstIdByType(
               INPAR::MAT::m_particlemat_ellipsoids) >= 0)
    dserror("Cannot have materials for spherical and ellipsoidal particles at the same time!");

  // check
  if (id == -1) dserror("Could not find particle material or material type");

  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  particleMat_.push_back(static_cast<const MAT::PAR::ParticleMat* const>(mat));

  return;
}

/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::PrepareTimeStep(bool print_header)
{
  IncrementTimeAndStep();
  if (print_header) PrintHeader();

  // apply dirichlet boundary conditions
  particles_->PrepareTimeStep();

  return;
}

/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Integrate()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::Integrate");

  particles_->UpdateExtActions();

  // set states of wall discretization
  SetUpWallDiscret();

  // particle time integration
  particles_->IntegrateStep();

  // check bin size for current time step
  BinSizeSafetyCheck();

  return;
}

/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Update()
{
  // write state vectors from n+1 to n
  particles_->Update();

  return;
}

/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ReadRestart(int restart)
{
  RemoveAllParticles();

  // read in particles for restart
  {
    IO::DiscretizationReader reader(BinStrategy()->BinDiscret(), restart);
    reader.ReadNodesOnly(restart);
  }

  // Init() is needed to obtain connectivity -> includes FillComplete())
  Init(true);

  // now, correct map layouts are available and states can be read
  particles_->ReadRestart(restart);
  SetTimeStep(particles_->TimeOld(), restart);

  return;
}

/*----------------------------------------------------------------------*
| dynamic load balancing for bin distribution               ghamm 08/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::DynamicLoadBalancing()
{
  // decide whether it is time for dynamic load balancing
  if (Step() % 100 != 0 or Comm().NumProc() == 1) return;

  // transfer particles into their correct bins
  if (rep_strategy_ == INPAR::PARTICLEOLD::repstr_adaptive) TransferParticles(true);

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::DynamicLoadBalancing()");

  const Epetra_Map* oldrowmap = BinStrategy()->BinDiscret()->ElementRowMap();

  Teuchos::RCP<const Epetra_CrsGraph> constgraph = CreateGraph();

  // Now we're going to create a Epetra_Vector with vertex weights to
  // be used in the repartitioning operation.
  Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*oldrowmap, false);
  // weights must be at least one for zoltan
  double* vals = vweights->Values();
  for (int i = 0; i < oldrowmap->NumMyElements(); ++i)
  {
    const int numnode = BinStrategy()->BinDiscret()->lRowElement(i)->NumNode();
    vals[i] = 1.0 + numnode * 3 + numnode * numnode;
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
  BinStrategy()->BinDiscret()->ExportRowElements(*newelerowmap);

  // export row nodes to new layout
  {
    // create a set of row particle IDs for each proc
    std::set<int> particles;
    for (int lid = 0; lid < newelerowmap->NumMyElements(); ++lid)
    {
      DRT::Element* bin = BinStrategy()->BinDiscret()->gElement(newelerowmap->GID(lid));
      const int* particleids = bin->NodeIds();
      for (int iparticle = 0; iparticle < bin->NumNode(); ++iparticle)
        particles.insert(particleids[iparticle]);
    }

    // copy particle gids to a vector and create particle row map
    std::vector<int> rowparticles(particles.begin(), particles.end());
    Teuchos::RCP<Epetra_Map> particlerowmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)rowparticles.size(), &rowparticles[0], 0, Comm()));

    // place all nodes on the correct processor
    BinStrategy()->BinDiscret()->ExportRowNodes(*particlerowmap);
  }

  // ghost bins and particles according to the bins --> final FillComplete() call included
  SetupGhosting(newelerowmap);

  // determine boundary bins (physical boundary as well as boundary to other procs)
  BinStrategy()->DetermineBoundaryRowBins();

  // determine one layer ghosting around boundary bins determined in previous step
  BinStrategy()->DetermineBoundaryColBinsIds();

  // update walls and connectivity
  if (particlewalldis_ != Teuchos::null)
  {
    if (not moving_walls_)
      BinStrategy()->ExtendEleGhosting(particlewalldis_, particlewallelecolmap_standardghosting_,
          BinColMap(), true, false, false);

    BuildElementToBinPointers(true);
  }

  // setup importer for dof and node based vectors
  if (particleInteractionType_ != INPAR::PARTICLEOLD::None) SetupImporter();

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  if (rep_strategy_ == INPAR::PARTICLEOLD::repstr_adaptive)
  {
    // update vector according to the new distribution of particles
    Teuchos::RCP<Epetra_Vector> temp = dis_at_last_redistr_;
    dis_at_last_redistr_ = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*temp, *dis_at_last_redistr_);
  }

  DRT::UTILS::PrintParallelDistribution(*BinStrategy()->BinDiscret());

  return;
}

/*----------------------------------------------------------------------*
 | safety check for proper bin size                        sfuchs 02/18 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BinSizeSafetyCheck()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::BinSizeSafetyCheck");

  if (particleInteractionType_ == INPAR::PARTICLEOLD::None) return;

  // note: checking two important criteria concerning bin size to maintain proper particle
  // interaction 1) the particle interaction distance may not be larger than one bin size 2) the
  // maximum particle movement in one time step may not be larger than one bin size

  // particle interaction distance
  double interaction_distance = ParticleInteractionDistance();

  if (interaction_distance > BinStrategy()->CutoffRadius())
    dserror("The particle interaction distance is larger than one bin size (%f > %f)!",
        interaction_distance, BinStrategy()->CutoffRadius());

  // displacement increment in this time step
  Teuchos::RCP<Epetra_Vector> dis_incr =
      LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);

  LINALG::Matrix<3, 1> lastpos;
  LINALG::Matrix<3, 1> currpos;

  for (int i = 0; i < BinStrategy()->BinDiscret()->NumMyRowNodes(); ++i)
  {
    lastpos.Clear();
    currpos.Clear();

    // get a pointer at i-th row node
    DRT::Node* node = BinStrategy()->BinDiscret()->lRowNode(i);

    // get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = BinStrategy()->BinDiscret()->Dof(node);

    for (int dim = 0; dim < 3; ++dim)
    {
      int doflid = (*particles_->Dispn()).Map().LID(dofnode[dim]);
      lastpos(dim) = (*particles_->Dispn())[doflid];
      currpos(dim) = (*particles_->Dispnp())[doflid];

      // consider periodic boundary conditions in current spatial direction
      if (BinStrategy()->HavePBC(dim))
      {
        // periodic length in current spatial direction
        double pbc_length = BinStrategy()->PBCDelta(dim);

        // shift position back if necessary
        if (lastpos(dim) - currpos(dim) > 0.5 * pbc_length)
          currpos(dim) += pbc_length;
        else if (lastpos(dim) - currpos(dim) < -0.5 * pbc_length)
          currpos(dim) -= pbc_length;
      }

      (*dis_incr)[doflid] = currpos(dim) - lastpos(dim);
    }
  }

  // get maximum displacement increment over all procs
  double extrema[2] = {0.0, 0.0};
  dis_incr->MinValue(&extrema[0]);
  dis_incr->MaxValue(&extrema[1]);
  const double gmaxdisincr = std::max(-extrema[0], extrema[1]);

#ifdef DEBUG
  if (MyRank() == 0)
    std::cout << "maximum particle movement in this time step " << gmaxdisincr << std::endl;
#endif

  if (gmaxdisincr > BinStrategy()->CutoffRadius())
  {
    // find entry with maximum particle movement
    int doflid = -1;
    for (int i = 0; i < dis_incr->MyLength(); ++i)
    {
      if ((*dis_incr)[i] < extrema[0] + 1.0e-12 || (*dis_incr)[i] > extrema[1] - 1.0e-12)
      {
        doflid = i;
        break;
      }
    }

    const int gid = dis_incr->Map().GID(doflid);

    dserror("Particle %i (gid) travels more than one bin per time step (%f > %f)!", (int)gid / 3,
        gmaxdisincr, BinStrategy()->CutoffRadius());
  }

  return;
}

/*----------------------------------------------------------------------*
 | update connectivity                                     sfuchs 08/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::UpdateConnectivity()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::UpdateConnectivity");

  // relate wall gids to bins ids
  if (moving_walls_) RelateWallGidsToBinIds();

  // set current position to particle nodes
  SetParticleNodePos();

  // transfer particles into their correct bins adaptively
  if (rep_strategy_ == INPAR::PARTICLEOLD::repstr_adaptive)
  {
    // check if repartitioning is necessary
    if (CheckAdaptiveRepartition())
    {
      // transfer particles into their correct bins
      TransferParticles(true);

      // setup importer for dof and node based vectors
      if (particleInteractionType_ != INPAR::PARTICLEOLD::None) SetupImporter();

      // update displacement state of particles after redistribution
      dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*particles_->Dispn()));
    }
  }
  // transfer particles into their correct bins every time step
  else if (rep_strategy_ == INPAR::PARTICLEOLD::repstr_everydt)
  {
    // transfer particles into their correct bins
    TransferParticles(true);

    // setup importer for dof and node based vectors
    if (particleInteractionType_ != INPAR::PARTICLEOLD::None) SetupImporter();
  }
  // default
  else
    dserror("Unknown particle repartitioning strategy!");

  // assign wall elements and gids to bins (after repartitioning of binning discretization)
  if (moving_walls_) AssignWallElesAndGidsToBins();

  return;
}

/*----------------------------------------------------------------------*
 | set position of particle nodes                          sfuchs 02/18 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetParticleNodePos()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::SetParticleNodePos");

  // Get current displacements
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();
  if (disnp == Teuchos::null) dserror("Attention: disnp is not set!");

  std::vector<int> examinedbins(BinStrategy()->BinDiscret()->NumMyRowElements(), 0);
  // first run over particles and then process whole bin in which particle is located
  // until all particles have been checked
  for (int i = 0; i < BinStrategy()->BinDiscret()->NumMyRowNodes(); ++i)
  {
    DRT::Node* currparticle = BinStrategy()->BinDiscret()->lRowNode(i);

#ifdef DEBUG
    if (currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
#endif

    DRT::MESHFREE::MeshfreeMultiBin* currbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currparticle->Elements()[0]);
    // as checked above, there is only one element in currele array
    const int binId = currbin->Id();
    const int rlid = BinStrategy()->BinDiscret()->ElementRowMap()->LID(binId);

    // if a bin has already been examined --> continue with next particle
    if (examinedbins[rlid]) continue;
    // else: bin is examined for the first time --> new entry in examinedbins
    else
      examinedbins[rlid] = 1;

    DRT::Node** particles = currbin->Nodes();
    for (int iparticle = 0; iparticle < currbin->NumNode(); ++iparticle)
    {
      // get current node
      DRT::Node* currnode = particles[iparticle];

      // get the first gid of a node and convert it into a lid
      const int gid = BinStrategy()->BinDiscret()->Dof(currnode, 0);
      const int lid = disnp->Map().LID(gid);
      if (lid < 0) dserror("displacement for node %d not stored on this proc: %d", gid, MyRank());

      // get current position
      LINALG::Matrix<3, 1> pos_mat(true);
      for (int dim = 0; dim < 3; ++dim) pos_mat(dim) = (*disnp)[lid + dim];

      // shift in case of periodic boundary conditions
      BinStrategy()->PeriodicBoundaryShift3D(pos_mat);

      for (int dim = 0; dim < 3; ++dim) (*disnp)[lid + dim] = pos_mat(dim);

      // transform to array
      std::vector<double> pos(3, 0.0);
      for (int dim = 0; dim < 3; ++dim) pos[dim] = pos_mat(dim);

      // change X() of current particle
      currnode->SetPos(pos);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | check movement of particles since last redistribution   sfuchs 02/18 |
 *----------------------------------------------------------------------*/
bool PARTICLE::Algorithm::CheckAdaptiveRepartition()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::CheckAdaptiveRepartition");

#ifdef DEBUG
  // safety check
  if (not dis_at_last_redistr_->Map().SameAs(*BinStrategy()->BinDiscret()->DofRowMap()))
    dserror(
        "Current particle dof row map and map of disp vector after last redistribution are not the "
        "same!");
#endif

  // displacement increment from particle position at last redistribution to current position
  Teuchos::RCP<Epetra_Vector> dis_incr =
      LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);

  LINALG::Matrix<3, 1> lastpos;
  LINALG::Matrix<3, 1> currpos;

  for (int i = 0; i < BinStrategy()->BinDiscret()->NumMyRowNodes(); ++i)
  {
    lastpos.Clear();
    currpos.Clear();

    // get a pointer at i-th row node
    DRT::Node* node = BinStrategy()->BinDiscret()->lRowNode(i);

    // get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = BinStrategy()->BinDiscret()->Dof(node);

    for (int dim = 0; dim < 3; ++dim)
    {
      int doflid = dis_at_last_redistr_->Map().LID(dofnode[dim]);
      lastpos(dim) = (*dis_at_last_redistr_)[doflid];
      currpos(dim) = (*particles_->Dispn())[doflid];

      // consider periodic boundary conditions in current spatial direction
      if (BinStrategy()->HavePBC(dim))
      {
        // periodic length in current spatial direction
        double pbc_length = BinStrategy()->PBCDelta(dim);

        // shift position back if necessary
        if (lastpos(dim) - currpos(dim) > 0.5 * pbc_length)
          currpos(dim) += pbc_length;
        else if (lastpos(dim) - currpos(dim) < -0.5 * pbc_length)
          currpos(dim) -= pbc_length;
      }

      (*dis_incr)[doflid] = currpos(dim) - lastpos(dim);
    }
  }

  // get maximum displacement increment since last redistribution over all procs
  double extrema[2] = {0.0, 0.0};
  dis_incr->MinValue(&extrema[0]);
  dis_incr->MaxValue(&extrema[1]);
  const double gmaxdisincr = std::max(-extrema[0], extrema[1]);

#ifdef DEBUG
  if (MyRank() == 0)
    std::cout << "maximum particle movement since last redistribution " << gmaxdisincr << std::endl;
#endif

  // particle interaction distance
  double interaction_distance = ParticleInteractionDistance();

  // repartition of particles necessary
  // note: it is assumed that (in a worst case scenario) two particles approach each other with
  // maximum displacement
  bool repartition = ((interaction_distance + 2.0 * gmaxdisincr) > (BinStrategy()->CutoffRadius()));

  return repartition;
}

/*----------------------------------------------------------------------*
 | particles are checked and transferred if necessary       ghamm 10/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::list<int>> PARTICLE::Algorithm::TransferParticles(
    const bool updatestates, const bool ghosting)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles");

  // leave here in case nothing to do
  if (particles_->Radiusn()->GlobalLength() == 0) return Teuchos::rcp(new std::list<int>(0));

  // transfer particles to new bins
  Teuchos::RCP<std::list<int>> deletedparticles =
      PARTICLE::ParticleHandler::TransferParticles(ghosting);

  // check whether all procs have a filled bindis_,
  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
  BinStrategy()->BinDiscret()->CheckFilledGlobally();

  // new ghosting if necessary
  if (ghosting)
    BinStrategy()->BinDiscret()->ExtendedGhosting(*ExtendedBinColMap(), true, false, true, false);
  else
    BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
  BuildElementToBinPointers(not moving_walls_);

  // update state vectors in time integrator to the new layout
  if (updatestates)
  {
    TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles::UpdateStates");
    particles_->UpdateStatesAfterParticleTransfer();
    UpdateStates();
  }

  return deletedparticles;
}

/*----------------------------------------------------------------------*
 | setup importer for dof and node based vectors           sfuchs 02/18 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupImporter()
{
  // setup importer for dof based vectors
  dofimporter_ = Teuchos::rcp(new Epetra_Import(
      *BinStrategy()->BinDiscret()->DofColMap(), *BinStrategy()->BinDiscret()->DofRowMap()));

  // setup importer for node based vectors
  nodeimporter_ = Teuchos::rcp(new Epetra_Import(
      *BinStrategy()->BinDiscret()->NodeColMap(), *BinStrategy()->BinDiscret()->NodeRowMap()));

  return;
}

/*----------------------------------------------------------------------*
 | extended bin column map                                 sfuchs 06/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::Algorithm::ExtendedBinColMap()
{
  std::set<int> colbins;

  // insert standard one layer ghosting
  for (int lid = 0; lid < BinColMap()->NumMyElements(); ++lid)
    colbins.insert(BinColMap()->GID(lid));

  // insert an additional (second) ghost layer
  if (extendedGhosting_ == INPAR::PARTICLEOLD::AddLayerGhosting) AddLayerGhosting(colbins);

  // insert bins in proximity of ghosted wall elements
  if (extendedGhosting_ == INPAR::PARTICLEOLD::WallElementGhosting) WallElementGhosting(colbins);

  std::vector<int> colbinsvec(colbins.begin(), colbins.end());

  return Teuchos::rcp(
      new Epetra_Map(-1, (int)colbinsvec.size(), &colbinsvec[0], 0, BinColMap()->Comm()));
}

/*----------------------------------------------------------------------*
 | ghosting of an additional second bin layer              sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::AddLayerGhosting(std::set<int>& colbins)
{
  // get boundary column bins
  std::set<int> bdrycolbins = BinStrategy()->BoundaryColBinsIds();

  std::vector<int> binvec(27);
  std::set<int>::const_iterator biniter;
  for (biniter = bdrycolbins.begin(); biniter != bdrycolbins.end(); ++biniter)
  {
    binvec.clear();
    BinStrategy()->GetNeighborAndOwnBinIds(*biniter, binvec);
    colbins.insert(binvec.begin(), binvec.end());
  }

  return;
}

/*----------------------------------------------------------------------*
 | ghosting of bins in proximity of wall elements          sfuchs 08/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::WallElementGhosting(std::set<int>& colbins)
{
  // insert bins that are touched by a relevant wall element
  std::map<int, std::set<int>>::const_iterator iter;
  for (iter = relwallgidtobinids_.begin(); iter != relwallgidtobinids_.end(); ++iter)
  {
    std::set<int>::const_iterator biniter;
    std::vector<int> binvec(27);
    for (biniter = iter->second.begin(); biniter != iter->second.end(); ++biniter)
    {
      binvec.clear();
      BinStrategy()->GetNeighborAndOwnBinIds(*biniter, binvec);
      colbins.insert(binvec.begin(), binvec.end());
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | set wall states                                         sfuchs 02/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetWallStates(
    Teuchos::RCP<const Epetra_Vector> walldispn,   ///< wall displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> walldispnp,  ///< wall displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> wallvelnp    ///< wall velocities at \f$t_{n+1}\f$
)
{
  walldispn_ = walldispn;
  walldispnp_ = walldispnp;
  wallvelnp_ = wallvelnp;

  return;
}

/*------------------------------------------------------------------------*
 | set up wall discretizations                               sfuchs 04/17 |
 *------------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetUpWallDiscret()
{
  if (particlewalldis_ != Teuchos::null)
  {
    // supply states for discretization
    particlewalldis_->SetState("walldisn", walldispn_);
    particlewalldis_->SetState("walldisnp", walldispnp_);
    particlewalldis_->SetState("wallvelnp", wallvelnp_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | access structure and setup particle wall                 ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::AccessStructure()
{
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  // initialize structure if necessary
  int numstructnode = structdis->NumMyColElements();
  int eleexist = 0;
  structdis->Comm().MaxAll(&numstructnode, &eleexist, 1);

  if (eleexist)
  {
    // fill discretization
    if (not structdis->Filled()) structdis->FillComplete();

    // add fully redundant discretization for particle walls with identical dofs to full structural
    // discret
    SetupParticleWalls(structdis, "BELE3_3");

    // assign wall elements to bins initially once for fixed walls (additionally rebuild pointers
    // after ghosting)
    if (not moving_walls_)
    {
      RelateWallGidsToBinIds();
      AssignWallElesAndGidsToBins();
    }
  }

  // safety check
  else if (moving_walls_)
    dserror("Moving walls indicated in input file despite empty structure discretization!");

  return;
}

/*----------------------------------------------------------------------*
 | particle walls are added from the structural discret     ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupParticleWalls(
    Teuchos::RCP<DRT::Discretization> basediscret, const std::string elename)
{
  // elename = "BELE3_3" or "BELE3_4"
  // number of dofs is important for transparent dof set
  // only zeros are applied to the wall displacements when fluid domain is basediscret
  // -> number of dofs is irrelevant when reading data for wall discret in this case
  // future implementation using ALE needs to be handled like a structure

  //--------------------------------------------------------------------
  // 1st step: build fully redundant discretization with wall elements
  //--------------------------------------------------------------------

  // declare struct objects in wall condition
  std::map<int, Teuchos::RCP<DRT::Element>> structgelements;  // col map of structure elements
  std::map<int, DRT::Node*> dummy2;                           // dummy map
  std::map<int, DRT::Node*> structgnodes;                     // col map of structure nodes

  // initialize struct objects in wall condition
  DRT::UTILS::FindConditionObjects(
      *basediscret, dummy2, structgnodes, structgelements, "ParticleWall");

  // initialize new particle wall discretizations
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(basediscret->Comm().Clone());
  const std::string discret_name = "particlewalls";
  Teuchos::RCP<DRT::Discretization> particlewalldis =
      Teuchos::rcp(new DRT::Discretization(discret_name, com));

  std::vector<int> nodeids;
  std::vector<int> eleids;

  {
    // care about particle wall nodes:
    // fill everything in case of static walls and only row nodes in case of moving walls
    for (std::map<int, DRT::Node*>::iterator nit = structgnodes.begin(); nit != structgnodes.end();
         ++nit)
    {
      DRT::Node* currnode = (*nit).second;
      if (not moving_walls_ || ((currnode->Owner() == MyRank()) && moving_walls_))
      {
        nodeids.push_back(currnode->Id());
        particlewalldis->AddNode(
            Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
      }
    }

    // care about particle wall eles:
    // fill everything in case of static walls and only row elements in case of moving walls
    for (std::map<int, Teuchos::RCP<DRT::Element>>::iterator eit = structgelements.begin();
         eit != structgelements.end(); ++eit)
    {
      Teuchos::RCP<DRT::Element> currele = eit->second;
      if (not moving_walls_ || ((currele->Owner() == MyRank()) && moving_walls_))
      {
        eleids.push_back(currele->Id());
        // structural surface elements cannot be distributed --> Bele3 element is used
        Teuchos::RCP<DRT::Element> wallele =
            DRT::UTILS::Factory(elename, "Polynomial", currele->Id(), currele->Owner());
        wallele->SetNodeIds(currele->NumNode(), currele->NodeIds());
        particlewalldis->AddElement(wallele);
      }
    }
  }

  // extended ghosting of wall elements for static walls
  if (not moving_walls_)
  {
    // fill complete in order to obtain element col map
    particlewalldis->FillComplete(false, false, false);

    particlewallelecolmap_standardghosting_ =
        Teuchos::rcp(new Epetra_Map(*particlewalldis->ElementColMap()));

    // extend ghosting to the bin col map
    BinStrategy()->ExtendEleGhosting(
        particlewalldis, particlewallelecolmap_standardghosting_, BinColMap(), false, false, false);
  }
  else  // ... or otherwise do fully redundant storage of wall elements in case of moving boundaries
  {
    // row node map of walls
    Teuchos::RCP<Epetra_Map> wallnoderowmap =
        Teuchos::rcp(new Epetra_Map(-1, nodeids.size(), &nodeids[0], 0, particlewalldis->Comm()));

    // fully overlapping node map
    Teuchos::RCP<Epetra_Map> wallrednodecolmap = LINALG::AllreduceEMap(*wallnoderowmap);

    // row ele map of walls
    Teuchos::RCP<Epetra_Map> wallelerowmap =
        Teuchos::rcp(new Epetra_Map(-1, eleids.size(), &eleids[0], 0, particlewalldis->Comm()));
    // fully overlapping ele map
    Teuchos::RCP<Epetra_Map> wallredelecolmap = LINALG::AllreduceEMap(*wallelerowmap);

    // do the fully overlapping ghosting of the wall elements to have everything redundant
    particlewalldis->ExportColumnNodes(*wallrednodecolmap);
    particlewalldis->ExportColumnElements(*wallredelecolmap);
  }

  // find out if we are in parallel; needed for TransparentDofSet
  bool parallel = (particlewalldis->Comm().NumProc() == 1) ? false : true;

  // dofs of the original discretization are used to set same dofs for the new particle wall
  // discretization
  Teuchos::RCP<DRT::DofSet> newdofset =
      Teuchos::rcp(new DRT::TransparentDofSet(basediscret, parallel));
  particlewalldis->ReplaceDofSet(newdofset);
  newdofset = Teuchos::null;

  // final fill complete to reorganize everything in the discretization
  particlewalldis->FillComplete(true, false, false);
  particlewalldis_ = particlewalldis;

  // some output to screen and initialization of binary output
  if (MyRank() == 0) std::cout << "after adding particle walls" << std::endl;
  DRT::UTILS::PrintParallelDistribution(*particlewalldis_);

  if (moving_walls_)
    wallextractor_ = Teuchos::rcp(new LINALG::MapExtractor(
        *(basediscret->DofRowMap()), Teuchos::rcp(particlewalldis_->DofRowMap(), false)));

  // write initial wall discret for visualization
  particlewalldis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(particlewalldis_)));
  Teuchos::RCP<IO::DiscretizationWriter> output = particlewalldis_->Writer();
  output->WriteMesh(Step(), Time());
  output->NewStep(Step(), Time());
  output->WriteVector("displacement", LINALG::CreateVector(*particlewalldis_->DofRowMap(), true));

  // initialize wall states
  walldispn_ = LINALG::CreateVector(*particlewalldis_->DofRowMap(), true);
  walldispnp_ = LINALG::CreateVector(*particlewalldis_->DofRowMap(), true);
  wallvelnp_ = LINALG::CreateVector(*particlewalldis_->DofRowMap(), true);

  return;
}

/*----------------------------------------------------------------------*
 | relate wall gids to bin ids                              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::RelateWallGidsToBinIds()
{
  // clear map relating wall gids to bin ids
  relwallgidtobinids_.clear();

  // remove assigned wall elements
  BinStrategy()->RemoveSpecificElesFromBins(bin_wallcontent_);

  std::map<int, LINALG::Matrix<3, 1>> currentpositions;
  if (moving_walls_)
  {
    Teuchos::RCP<const Epetra_Vector> walldisn = particlewalldis_->GetState("walldisn");
    for (int lid = 0; lid < particlewalldis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = particlewalldis_->lColNode(lid);
      std::vector<int> lm_node;
      lm_node.reserve(3);
      particlewalldis_->Dof(node, lm_node);

      // nodal displacements
      static LINALG::Matrix<3, 1> node_disn;
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3, 1>>(*walldisn, node_disn, lm_node);

      LINALG::Matrix<3, 1> currpos;
      const double* X = node->X();
      for (int dim = 0; dim < 3; ++dim) currpos(dim) = X[dim] + node_disn(dim);
      currentpositions[node->Id()] = currpos;
    }
  }
  else
  {
    for (int lid = 0; lid < particlewalldis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = particlewalldis_->lColNode(lid);

      LINALG::Matrix<3, 1> currpos;
      const double* X = node->X();
      for (int dim = 0; dim < 3; ++dim) currpos(dim) = X[dim];
      currentpositions[node->Id()] = currpos;
    }
  }

  // circumcircle of bin
  double bincircumcircle = 0.0;
  for (int dim = 0; dim < 3; ++dim)
  {
    bincircumcircle += std::pow(BinStrategy()->BinSize()[dim] / 2.0, 2.0);
  }
  bincircumcircle = sqrt(bincircumcircle);

  // minimal bin size
  double min_bin_size = BinStrategy()->BinSize()[0];
  for (int dim = 1; dim < 3; ++dim)
    min_bin_size = std::min(min_bin_size, BinStrategy()->BinSize()[dim]);

  // find bins for all wall elements
  const int numcolwalleles = particlewalldis_->NumMyColElements();
  for (int lid = 0; lid < numcolwalleles; ++lid)
  {
    DRT::Element* wallele = particlewalldis_->lColElement(lid);
    const int* nodeids = wallele->NodeIds();
    const int numnode = wallele->NumNode();

    // variable to store bin ids in which this wall element is located
    std::set<int> binIds;

    // do a positive search and get all bins enclosed in the bounding box of each wall element
    {
      // initialize ijk_range with ijk of first node of wall element
      int ijk[3];
      BinStrategy()->ConvertPosToijk(currentpositions[nodeids[0]], ijk);

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j = 1; j < numnode; ++j)
      {
        int ijk[3];
        BinStrategy()->ConvertPosToijk(currentpositions[nodeids[j]], ijk);

        for (int dim = 0; dim < 3; ++dim)
        {
          if (ijk[dim] < ijk_range[dim * 2]) ijk_range[dim * 2] = ijk[dim];
          if (ijk[dim] > ijk_range[dim * 2 + 1]) ijk_range[dim * 2 + 1] = ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range and fill them into binIds
      BinStrategy()->GidsInijkRange(&ijk_range[0], binIds, false);
    }  // do a positive search

    // if no bins on this proc were found, next wall element can be processed
    if (binIds.empty()) continue;

    // do a negative search and remove bins that are too far away from the wall element
    {
      std::set<int> binfaraway;
      for (std::set<int>::const_iterator biniter = binIds.begin(); biniter != binIds.end();
           ++biniter)
      {
        const LINALG::Matrix<3, 1> bincentroid = BinStrategy()->GetBinCentroid(*biniter);

        // search for the closest object, more exactly it's coordinates
        LINALG::Matrix<3, 1> minDistCoords;
        GEO::nearest3DObjectOnElement(wallele, currentpositions, bincentroid, minDistCoords);

        LINALG::Matrix<3, 1> distance;
        distance.Update(1.0, bincentroid, -1.0, minDistCoords);
        double dist = distance.Norm2();

        // if distance is larger than radius of circumcircle of bin --> too far away
        if (dist > bincircumcircle)
        {
          binfaraway.insert(*biniter);
        }
        // if distance is smaller than half the minimum bin size --> very close
        else if (dist <= min_bin_size * 0.5)
        {
          continue;
        }
        // if distance is between half the minimum bin size and radius of the circumcircle -->
        // further checks
        else
        {
          std::vector<LINALG::Matrix<3, 1>> bincorners;
          BinStrategy()->GetBinCorners(*biniter, bincorners);

          // in case wall element is axis aligned, it might not be detected as inside because
          // projection points are located on the edges of the bin --> Remedy: bin centroid is
          // tested as well
          bincorners.push_back(bincentroid);

          bool projpointinsidebin = false;
          // all corners of the close bin are projected onto the wall element: if at least one
          // projection point is inside the bin, it won't be removed from the list
          for (size_t corner = 0; corner < bincorners.size(); ++corner)
          {
            // search for the closest object, more exactly it's coordinates
            LINALG::Matrix<3, 1> minDistCoords;
            GEO::nearest3DObjectOnElement(
                wallele, currentpositions, bincorners[corner], minDistCoords);

            const int gid = BinStrategy()->ConvertPosToGid(minDistCoords);
            if (gid == *biniter)
            {
              projpointinsidebin = true;
              break;
            }
          }
          if (projpointinsidebin == false) binfaraway.insert(*biniter);
        }
      }

      // erase bins that are far away from a wall element
      for (std::set<int>::const_iterator biniter = binfaraway.begin(); biniter != binfaraway.end();
           ++biniter)
        binIds.erase(*biniter);
    }  // do a negative search

    // if none of found bins is in BinColMap, next wall element can be processed
    {
      bool noBinInBinColMap = true;
      for (std::set<int>::const_iterator biniter = binIds.begin(); biniter != binIds.end();
           ++biniter)
        if (BinColMap()->LID(*biniter) >= 0)
        {
          noBinInBinColMap = false;
          continue;
        }

      if (noBinInBinColMap) continue;
    }

    // insert found bins into map relating wall gids to bin ids
    relwallgidtobinids_[wallele->Id()].insert(binIds.begin(), binIds.end());

  }  // end lid

  return;
}

/*----------------------------------------------------------------------*
 | single fields are tested                                 ghamm 09/12 |
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

  return;
}

/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 10/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Output(bool forced_writerestart /*= false*/)
{
  // INFO regarding output: Bins are not written to file because they cannot
  // be post-processed anyway (no nodes and connectivity available)
  particles_->OutputStep(forced_writerestart);

  if (writeresultsevery_ and (Step() % writeresultsevery_ == 0))
  {
    // visualize bins according to specification in input file
    BinStrategy()->WriteBinOutput(Step(), Time());

    // write moving walls to file
    if (moving_walls_)
    {
      particlewalldis_->Writer()->NewStep(Step(), Time());
      particlewalldis_->Writer()->WriteVector("displacement", walldispnp_);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | get neighbouring particles and walls                    ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetNeighbouringItems(DRT::Node* particle,
    std::list<DRT::Node*>& neighboursLinf_p,
    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w) const
{
  if (particle->NumElement() != 1) dserror("More than one element for this particle");

  DRT::Element** CurrentBin = particle->Elements();

  GetNeighbouringItems(CurrentBin[0]->Id(), neighboursLinf_p, &neighboursLinf_w);

  return;
}

/*----------------------------------------------------------------------*
 | get neighbouring particles and walls (bin version)      katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetNeighbouringItems(const int binId,
    std::list<DRT::Node*>& neighboursLinf_p,
    boost::unordered_map<int, DRT::Element*>* neighboursLinf_w) const
{
  std::vector<int> binIds;
  binIds.reserve(27);

  BinStrategy().GetNeighborAndOwnBinIds(binId, binIds);

  GetBinContent(neighboursLinf_p, neighboursLinf_w, binIds);

  return;
}

/*----------------------------------------------------------------------*
 | get particles and wall elements in given bins           ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetBinContent(std::list<DRT::Node*>& bin_p,
    boost::unordered_map<int, DRT::Element*>* bin_w, std::vector<int>& binIds) const
{
  // loop over all bins
  for (std::vector<int>::const_iterator bin = binIds.begin(); bin != binIds.end(); ++bin)
  {
    // extract bins from discretization after checking on existence
    const int lid = BinStrategy().BinDiscret()->ElementColMap()->LID(*bin);
    if (lid < 0) continue;

#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(
        BinStrategy().BinDiscret()->lColElement(lid));
    if (test == NULL)
      dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* neighboringbin =
        static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy().BinDiscret()->lColElement(lid));

    // gather particles
    DRT::Node** nodes = neighboringbin->Nodes();
    bin_p.insert(bin_p.end(), nodes, nodes + neighboringbin->NumNode());

    // gather wall elements
    if (bin_w != NULL)
    {
      DRT::Element** walleles = neighboringbin->AssociatedEles(bin_wallcontent_);
      const int numwalls = neighboringbin->NumAssociatedEle(bin_wallcontent_);
      const int* ids = neighboringbin->AssociatedEleIds(bin_wallcontent_);
      for (int iwall = 0; iwall < numwalls; ++iwall)
      {
        if (ids[iwall] != walleles[iwall]->Id())
        {
          std::cout << "Comm().MyPID(): " << Comm().MyPID() << std::endl;
          std::cout << "ids[iwall]: " << ids[iwall] << std::endl;
          std::cout << "walleles[iwall]->Id(): " << walleles[iwall]->Id() << std::endl << std::endl;
          dserror(
              "Ids of wall elements are different from Ids stored in multibin! Accidential access "
              "of random storage?");
        }

        (*bin_w)[walleles[iwall]->Id()] = walleles[iwall];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | determine particle interaction distance                 sfuchs 02/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::Algorithm::ParticleInteractionDistance()
{
  // determine maximum particle interaction distance for different particle interaction types
  // TODO: differentiate the case of particle contact in combination with adhesive forces (active
  // for distances > 2*radius)

  double maxrad = 0.0;
  particles_->Radiusn()->MaxValue(&maxrad);

  return 2.0 * maxrad;
}

/*------------------------------------------------------------------------*
 | return gravity acceleration                               meier 05/17  |
 *------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> PARTICLE::Algorithm::GetGravityAcc(const double time) { return gravity_acc_; }
