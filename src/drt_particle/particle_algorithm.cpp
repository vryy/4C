/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.cpp

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
#include "particle_utils.H"
#include "particle_node.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/extparticle_mat.H"
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
#include "particle_sph_rendering.H"

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
PARTICLE::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : AlgorithmBase(comm,params), ParticleHandler(comm),
  particles_(Teuchos::null),
  writeresultsevery_(0),
  particlewalldis_(Teuchos::null),
  particlewallelecolmap_standardghosting_(Teuchos::null),
  moving_walls_((bool)DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"MOVING_WALLS")),
  transfer_every_(DRT::Problem::Instance()->ParticleParams().get<int>("TRANSFER_EVERY")),
  particleInteractionType_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleInteractions>(DRT::Problem::Instance()->ParticleParams(),"PARTICLE_INTERACTION")),
  extendedGhosting_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ExtendedGhosting>(DRT::Problem::Instance()->ParticleParams(),"EXTENDED_GHOSTING")),
  particleMat_(NULL),
  particleMat2_(NULL),
  extParticleMat_(NULL),
  extParticleMat2_(NULL),
  bin_wallcontent_(BINSTRATEGY::UTILS::BELE3),
  rendering_(Teuchos::null)
{

  //Check number of space dimensions chosen for SPH weight functions
  if(particleInteractionType_==INPAR::PARTICLE::SPH)
  {
    if(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")!=INPAR::PARTICLE::particle_3D)
      dserror("The general Particle SPH Interactions framework (binning strategy etc.) does so far only cover 3D problems (DIMENSION  3D).\n"
              "However, if you want to treat quasi-2D or -1D problems, set the input parameter WEIGHT_FUNCTION_DIM to WF_2D or WF_1D, respectively.\n"
              "In Particle SPH Interactions, the definition of the weight functions have to be adapted according to the number of space dimensions considered!");

    if(MyRank() == 0)
    {
      switch(DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunctionDim>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION_DIM"))
      {
        case INPAR::PARTICLE::WF_3D :
          IO::cout << "Welcome to Particle SPH Interactions in 3D!" << IO::endl;
        break;
        case INPAR::PARTICLE::WF_2D :
          IO::cout << "Welcome to Particle SPH Interactions in 2D!" << IO::endl;
        break;
        case INPAR::PARTICLE::WF_1D :
          IO::cout << "Welcome to Particle SPH Interactions in 1D!" << IO::endl;
        break;
        default :  //nothing to do here
        break;
      }
    }

    INPAR::PARTICLE::DynamicType timinttype = DRT::INPUT::IntegralValue<INPAR::PARTICLE::DynamicType>(DRT::Problem::Instance()->ParticleParams(),"DYNAMICTYP");
    if(timinttype!=INPAR::PARTICLE::dyna_kickdrift and DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"TRANSPORT_VELOCITY")==true)
      dserror("Modified particle convection velocities based on a TRANSPORT_VELOCITY field only possible for KickDrift time integration scheme!");
  }
  else
  {
    if(DRT::Problem::Instance()->ParticleParams().get<double>("GRAVITY_RAMP_TIME")>=0.0)
      dserror("Ramp time for smooth increase of gravity force only possible for SPH applications so far!");
  }

  if (moving_walls_ == true and DRT::Problem::Instance()->ProblemType() != prb_pasi)
    dserror("MOVING_WALLS flag is activated!\n"
        "Set parameter PROBLEMTYP to 'Particle_Structure_Interaction' in ---PROBLEM TYP section\n"
        "and set the particle structure interaction parameters in ---PASI DYNAMIC section\n");

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

  if(BinStrategy()->ParticleDim() == INPAR::PARTICLE::particle_2Dz)
  {
    gravity_acc_(2) = 0.0;
    if(MyRank() == 0)
      IO::cout << "gravity in z-direction ignored as this is a pseudo-2D problem" << IO::endl;
  }

  // initial setup of particle discretization
  BinStrategy()->BinDiscret() = DRT::Problem::Instance()->GetDis("particle");
  // new dofs are numbered from zero, minnodgid is ignored and it does not register in static_dofsets_
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset = Teuchos::rcp(new DRT::IndependentDofSet(true));
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
  BinStrategy()->BinDiscret()->FillComplete(false,false,false);

  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*BinStrategy()->BinDiscret()->NodeRowMap()));

  Teuchos::RCP<Epetra_Map> binrowmap;
  if(not restarted)
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

  if(binrowmap->NumGlobalElements() > particlerowmap->NumGlobalElements() / 4.0 && MyRank() == 0)
    IO::cout << "\n\n\n CAREFUL: Reduction of number of bins recommended! Performance might be deteriorated. Increase cutoff radius. \n\n\n" << IO::endl;

  //--------------------------------------------------------------------
  // -> 1) create a list of homeless particles that are not in a bin on this proc
  std::list<Teuchos::RCP<DRT::Node> > homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = BinStrategy()->BinDiscret()->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBinsRoundRobin(homelessparticles);

  // ghost bins and particles according to the bins --> final FillComplete() call included
  SetupGhosting(binrowmap);

  // determine boundary bins (physical boundary as well as boundary to other procs)
  BinStrategy()->DetermineBoundaryRowBins();

  // determine one layer ghosting around boundary bins determined in previous step
  BinStrategy()->DetermineBoundaryColBinsIds();

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // set up the links to the materials for easy access
    // make sure that a particle material is defined in the dat-file
    InitMaterials();

    // get input parameters for particles
    const Teuchos::ParameterList& particledyn = DRT::Problem::Instance()->ParticleParams();

    // access structure and build particle walls
    AccessStructure();

    // create time integrator based on structural time integration
    Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm(particledyn, BinStrategy()->BinDiscret()));
    particles_ = particles->ParticleField();

    writeresultsevery_ = particledyn.get<int>("RESULTSEVRY");

    // set particle algorithm into time integration
    particles_->SetParticleAlgorithm(Teuchos::rcp(this,false));

    particles_->Init();

    // in case random noise is added to the particle position, particle transfer is necessary
    double amplitude = DRT::Problem::Instance()->ParticleParams().get<double>("RANDOM_AMPLITUDE");
    if(amplitude)
      TransferParticles(true, true);

    // determine consistent initial acceleration for the particles
    particles_->UpdateExtActions();

    particles_->DetermineMassDampConsistAccel();
  }
  else
  {
    // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
    BuildElementToBinPointers(not moving_walls_);
  }

  // some output
  if (MyRank() == 0)
    IO::cout << "after ghosting of particles" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*BinStrategy()->BinDiscret());

  return;
}

/*----------------------------------------------------------------------*
 | build connectivity from particle wall elements to bins   ghamm 04/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BuildElementToBinPointers(bool wallpointer)
{
  if(wallpointer == true)
  {
    // loop over column bins and fill wall elements
    const int numcolbin = BinStrategy()->BinDiscret()->NumMyColElements();
    for (int ibin=0; ibin<numcolbin; ++ibin)
    {
      DRT::Element* actele = BinStrategy()->BinDiscret()->lColElement(ibin);
      DRT::MESHFREE::MeshfreeMultiBin* actbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
      const int numwallele = actbin->NumAssociatedEle(bin_wallcontent_);
      const int* walleleids = actbin->AssociatedEleIds(bin_wallcontent_);
      std::vector<DRT::Element*> wallelements(numwallele);
      for(int iwall=0; iwall<numwallele; ++iwall)
      {
        const int wallid = walleleids[iwall];
        wallelements[iwall] = particlewalldis_->gElement(wallid);
      }
      actbin->BuildElePointers(bin_wallcontent_,&wallelements[0]);
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
  for ( std::map<int, std::set<int> >::iterator eleiter = relwallgidtobinids_.begin(); eleiter!=relwallgidtobinids_.end(); ++eleiter )
  {
    // global id of current wall element
    int wallgid = eleiter->first;

    // bins related to current wall element
    std::set<int> binIds = eleiter->second;

    // assign wall element to bins
    for ( std::set<int>::const_iterator biniter = binIds.begin(); biniter != binIds.end(); ++biniter )
    {
      if ( BinStrategy()->BinDiscret()->HaveGlobalElement(*biniter) )
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy()->BinDiscret()->gElement(*biniter))->AddAssociatedEle(bin_wallcontent_, wallgid, particlewalldis_->gElement(wallgid));
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | set up pointers to material bundles                     catta 09/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::InitMaterials()
{
  int id = -1;
  switch (particleInteractionType_)
  {
  case INPAR::PARTICLE::SPH :
  {
    int testid = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_extparticlemat);
    if(testid!=1)
      dserror("In SPH simulations, the first material ID has always to be 1!");

    id = 1;
    break;
  }
  default :
  {
    id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    if(id < 0)
      id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids);
    else if(DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat_ellipsoids) >= 0)
      dserror("Cannot have materials for spherical and ellipsoidal particles at the same time!");
    break;
  }
  }
  // check
  if (id==-1)
    dserror("Could not find particle material or material type");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  particleMat_ = static_cast<const MAT::PAR::ParticleMat*>(mat);
  if (particleInteractionType_ == INPAR::PARTICLE::SPH)
    extParticleMat_ = static_cast<const MAT::PAR::ExtParticleMat*>(mat);

  if(particleInteractionType_==INPAR::PARTICLE::SPH)
  {
    const INPAR::PARTICLE::FreeSurfaceType freeSurfaceType=DRT::INPUT::IntegralValue<INPAR::PARTICLE::FreeSurfaceType>(DRT::Problem::Instance()->ParticleParams(),"FREE_SURFACE_TYPE");
    if(freeSurfaceType==INPAR::PARTICLE::TwoPhase)
    {
      int id2 = 2;
      const MAT::PAR::Parameter* mat2 = DRT::Problem::Instance()->Materials()->ParameterById(id2);
      particleMat2_ = static_cast<const MAT::PAR::ParticleMat*>(mat2);
      extParticleMat2_ = static_cast<const MAT::PAR::ExtParticleMat*>(mat2);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::PrepareTimeStep(bool print_header)
{
  IncrementTimeAndStep();
  if (print_header)
    PrintHeader();

  // apply dirichlet boundary conditions
  particles_->PrepareTimeStep();

  // do rough safety check if bin size is appropriate
  BinSizeSafetyCheck(particles_->Dt());

  return;
}

/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Integrate()
{
  particles_->UpdateExtActions();

  SetUpWallDiscret();

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::Integrate");

  particles_->IntegrateStep();

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
  SetTimeStep(particles_->TimeOld(),restart);

  return;
}

/*----------------------------------------------------------------------*
| dynamic load balancing for bin distribution               ghamm 08/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::DynamicLoadBalancing()
{
  // decide whether it is time for dynamic load balancing
  if(Step()%100 != 0 or Comm().NumProc() == 1)
    return;

  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::DynamicLoadBalancing()");

  const Epetra_Map* oldrowmap = BinStrategy()->BinDiscret()->ElementRowMap();

  Teuchos::RCP<const Epetra_CrsGraph> constgraph = CreateGraph();

  // Now we're going to create a Epetra_Vector with vertex weights to
  // be used in the repartitioning operation.
  Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*oldrowmap, false);
  // weights must be at least one for zoltan
  double* vals = vweights->Values();
  for(int i=0; i<oldrowmap->NumMyElements(); ++i)
  {
    const int numnode = BinStrategy()->BinDiscret()->lRowElement(i)->NumNode();
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
  BinStrategy()->BinDiscret()->ExportRowElements(*newelerowmap);

  // export row nodes to new layout
  {
    // create a set of row particle IDs for each proc
    std::set<int> particles;
    for (int lid=0; lid<newelerowmap->NumMyElements(); ++lid)
    {
      DRT::Element* bin = BinStrategy()->BinDiscret()->gElement(newelerowmap->GID(lid));
      const int* particleids = bin->NodeIds();
      for(int iparticle=0; iparticle<bin->NumNode(); ++iparticle)
        particles.insert(particleids[iparticle]);
    }

    // copy particlegids to a vector and create particlerowmap
    std::vector<int> rowparticles(particles.begin(),particles.end());
    Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(-1,(int)rowparticles.size(),&rowparticles[0],0,Comm()));

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
  if(particlewalldis_ != Teuchos::null)
  {
    if( not moving_walls_)
      BinStrategy()->ExtendEleGhosting(particlewalldis_, particlewallelecolmap_standardghosting_, BinColMap(), true, false, false);

    BuildElementToBinPointers(true);
  }

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  DRT::UTILS::PrintParallelDistribution(*BinStrategy()->BinDiscret());

  return;
}

/*----------------------------------------------------------------------*
 | rough safety check for proper bin size                  ghamm 12/15  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BinSizeSafetyCheck(const double dt)
{
  // rough safety check whether bin size is large enough for proper contact detection
  if(particleInteractionType_ != INPAR::PARTICLE::None)
  {
    double extrema[2] = {0.0, 0.0};
    particles_->Veln()->MinValue(&extrema[0]);
    particles_->Veln()->MaxValue(&extrema[1]);
    const double maxvel = std::max(-extrema[0], extrema[1]);
    double maxrad = 0.0;

    switch (particleInteractionType_)
    {
    case INPAR::PARTICLE::SPH :
    {
      particles_->Radiusnp()->MaxValue(&maxrad);
      break;
    }
    default :
    {
      particles_->Radiusn()->MaxValue(&maxrad);
      break;
    }
    }

    //Determine effective interaction distance (maximal particle distance at which interaction forces are active) for different particle applications
    //TODO: Differentiate the case of particle contact in combination with adhesive forces (active for distances > 2*radius)
    double half_interaction_distance=0.0;
    if(particleInteractionType_==INPAR::PARTICLE::SPH)
      half_interaction_distance=maxrad/2.0;
    else
      half_interaction_distance=maxrad;

    //Here we changed the factor transfer_every_ applied in the previous revision (r23144) to (transfer_every_-1), since
    //the call of the method UpdateConnectivity() has been shifted to the location directly between displacement update and evaluation of
    //contact/interaction forces. Thus, in case UpdateConnectivity() is called in every time step (transfer_every_=1), the bin edge length
    //has at least to be as large as the interaction distance!
    if(half_interaction_distance + (transfer_every_-1)*maxvel*dt > 0.5*BinStrategy()->CutoffRadius())
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
            "Max radius is: %f", (int)gid/3, 2.0*(half_interaction_distance + maxvel*dt), BinStrategy()->CutoffRadius(), half_interaction_distance);
        std::cout << "Particle %i (gid) travels more than one bin per time step!!!!!!!" << std::endl;
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | particles are checked and transferred if necessary       ghamm 10/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<std::list<int> > PARTICLE::Algorithm::TransferParticles(const bool updatestates,
                                            const bool ghosting)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles");

  // leave here in case nothing to do
  if(particles_->Radiusn()->GlobalLength() == 0)
    return Teuchos::rcp(new std::list<int>(0));

  // Get current displacements
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();

  // transfer particles to new bins
  Teuchos::RCP<std::list<int> > deletedparticles =
      PARTICLE::ParticleHandler::TransferParticles( disnp, ghosting);

  // check whether all procs have a filled bindis_,
  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
  BinStrategy()->BinDiscret()->CheckFilledGlobally();

  // new ghosting if necessary
  if (ghosting)
    BinStrategy()->BinDiscret()->ExtendedGhosting(*ExtendedBinColMap(),true,false,true,false);
  else
    BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
  BuildElementToBinPointers(not moving_walls_);

  // update state vectors in time integrator to the new layout
  if(updatestates)
  {
    TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles::UpdateStates");
    particles_->UpdateStatesAfterParticleTransfer();
    UpdateStates();
  }

  return deletedparticles;
}

/*----------------------------------------------------------------------*
 | extended bin column map                                 sfuchs 06/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::Algorithm::ExtendedBinColMap()
{
  std::set<int> colbins;

  // insert standard one layer ghosting
  for (int lid=0; lid<BinColMap()->NumMyElements(); ++lid)
    colbins.insert(BinColMap()->GID(lid));

  // insert an additional (second) ghost layer
  if ( extendedGhosting_ == INPAR::PARTICLE::AddLayerGhosting )
    AddLayerGhosting(colbins);

  // insert bins in the proximity of ghosted boundary particles
  if ( extendedGhosting_ == INPAR::PARTICLE::BdryParticleGhosting )
    BdryParticleGhosting(colbins);

  // insert bins in proximity of ghosted wall elements
  if ( extendedGhosting_ == INPAR::PARTICLE::WallElementGhosting )
    WallElementGhosting(colbins);

  std::vector<int> colbinsvec(colbins.begin(),colbins.end());

  return Teuchos::rcp(new Epetra_Map(-1,(int)colbinsvec.size(),&colbinsvec[0],0,BinColMap()->Comm()));
}

/*----------------------------------------------------------------------*
 | ghosting of an additional second bin layer              sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::AddLayerGhosting(
    std::set<int>& colbins
    )
{
  // get boundary column bins
  std::set<int> bdrycolbins = BinStrategy()->BoundaryColBinsIds();

  std::vector<int> binvec(27);
  std::set< int >::const_iterator biniter;
  for ( biniter = bdrycolbins.begin(); biniter != bdrycolbins.end() ; ++biniter )
  {
    binvec.clear();
    BinStrategy()->GetNeighborAndOwnBinIds( *biniter, binvec );
    colbins.insert( binvec.begin(), binvec.end() );
  }

  return;
}

/*----------------------------------------------------------------------*
 | ghosting of bins in proximity of boundary particles     sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::BdryParticleGhosting(
    std::set<int>& colbins
    )
{
  // map relating extended column bins to receiving processor ids
  std::map< int, std::set<int> > towhomwhat;

  // set to store receiving processors ids
  std::set<int> receivingprocids;

  // get boundary row bins
  std::list<DRT::Element*> const bdryrowbins = BinStrategy()->BoundaryRowBins();
  if ( bdryrowbins.size() == 0 )
    dserror("Boundary row bins unknown, call function DetermineBoundaryRowBins() first!");

  // loop over boundary row bins
  std::list<DRT::Element*>::const_iterator iter;
  for ( iter = bdryrowbins.begin(); iter != bdryrowbins.end(); ++iter )
  {
    // current boundary row bin
    DRT::Element* bdryrowbin = *iter;

    bool containsbdryparticle = false;

    // loop over all particles in current boundary row bin
    DRT::Node** bdryrowbinparticles = bdryrowbin->Nodes();
    for ( int i = 0; i < bdryrowbin->NumNode(); ++i )
    {
      // get current particle
      PARTICLE::ParticleNode* particleNode_i = dynamic_cast<PARTICLE::ParticleNode*>(bdryrowbinparticles[i]);
      if ( particleNode_i == NULL )
        dserror("Dynamic cast to ParticleNode failed");

      // current particle is a boundary particle
      if ( particleNode_i->Is_bdry_particle() )
      {
        containsbdryparticle = true;
        break;
      }
    }

    // current boundary row bin contains no boundary particles
    if ( not containsbdryparticle )
      continue;

    // clear set containing receiving processors ids
    receivingprocids.clear();

    // get neighboring bins of current boundary row bin
    std::vector<int> neighborbins;
    neighborbins.reserve(26);
    BinStrategy()->GetNeighborBinIds( bdryrowbin->Id(), neighborbins );

    // loop over neighboring bins
    std::vector<int>::const_iterator biniter;
    for ( biniter = neighborbins.begin(); biniter != neighborbins.end(); ++biniter )
    {
      // get current bin and owner
      DRT::Element* neighboringbin = BinStrategy()->BinDiscret()->gElement(*biniter);
      int neighboringbinowner = neighboringbin->Owner();

      // neighboring bin on same processor or does not contain any particles
      if ( ( neighboringbinowner == MyRank() ) or ( neighboringbin->NumNode() == 0 ) )
        continue;

      // store owner of receiving neighboring bin
      receivingprocids.insert(neighboringbinowner);
    }

    // insert current neighboring bins in map
    std::set<int>::iterator iter;
    for ( iter = receivingprocids.begin(); iter != receivingprocids.end(); ++iter )
      towhomwhat[*iter].insert( neighborbins.begin(), neighborbins.end() );
  }

  // ---- send ---- ( we do not need to pack anything)
  DRT::Exporter exporter(BinStrategy()->BinDiscret()->Comm());
  const int length = towhomwhat.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  const int numproc = BinStrategy()->BinDiscret()->Comm().NumProc();
  std::vector<int> targetprocs( numproc, 0 );

  std::map< int, std::set<int> >::const_iterator p;
  for ( p = towhomwhat.begin(); p != towhomwhat.end(); ++p )
  {
    targetprocs[p->first] = 1;
    std::vector<int> bins( (p->second).begin(), (p->second).end() );
    exporter.ISend( MyRank(), p->first, &((bins)[0]), (int)(bins).size(), 1234, request[tag] );
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // ---- prepare receiving procs -----
  std::vector<int> summedtargets( numproc, 0) ;
  BinStrategy()->BinDiscret()->Comm().SumAll( targetprocs.data(), summedtargets.data(), numproc );

  // ---- receive ----- (we do not need to unpack anything)
  for( int rec = 0; rec < summedtargets[MyRank()]; ++rec )
  {
    std::vector<int> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny( from, tag ,rdata, length );
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", MyRank(), from);

    // insert received bins
    colbins.insert( rdata.begin(), rdata.end() );
  }

  // wait for all communication to finish
  for ( int i = 0; i < length; ++i )
    exporter.Wait( request[i] );

  // should be no time operation (if we have done everything correctly)
  BinStrategy()->BinDiscret()->Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | ghosting of bins in proximity of wall elements          sfuchs 08/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::WallElementGhosting(
    std::set<int>& colbins
    )
{
  // insert bins that are touched by a relevant wall element
  std::map<int, std::set<int> >::const_iterator iter;
  for( iter = relwallgidtobinids_.begin(); iter != relwallgidtobinids_.end(); ++iter )
  {
    std::set<int>::const_iterator biniter;
    std::vector<int> binvec(27);
    for( biniter = iter->second.begin(); biniter != iter->second.end(); ++biniter )
    {
      binvec.clear();
      BinStrategy()->GetNeighborAndOwnBinIds( *biniter, binvec );
      colbins.insert( binvec.begin(), binvec.end() );
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
  if(particlewalldis_ != Teuchos::null)
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

  if(eleexist)
  {
    // fill discretization
    if (not structdis->Filled())
      structdis->FillComplete();

    // add fully redundant discretization for particle walls with identical dofs to full structural discret
    SetupParticleWalls(structdis,"BELE3_3");

    // assign wall elements to bins initially once for fixed walls (additionally rebuild pointers after ghosting)
    if ( not moving_walls_ )
    {
      RelateWallGidsToBinIds();
      AssignWallElesAndGidsToBins();
    }
  }

  // safety check
  else if(moving_walls_)
    dserror("Moving walls indicated in input file despite empty structure discretization!");

  return;
}

/*----------------------------------------------------------------------*
 | particle walls are added from the structural discret     ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupParticleWalls(Teuchos::RCP<DRT::Discretization> basediscret,const std::string elename)
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
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > > structgelements; // col map of structure elements
  std::map<int, std::map<int, DRT::Node*> > dummy2;  // dummy map
  std::map<int, std::map<int, DRT::Node*> > structgnodes; // col map of structure nodes

  //initialize struct objects in wall condition
  DRT::UTILS::FindConditionObjects(*basediscret, dummy2, structgnodes, structgelements,"ParticleWall");

  // initialize new particle wall discretizations
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(basediscret->Comm().Clone());
  const std::string discret_name = "particlewalls";
  Teuchos::RCP<DRT::Discretization> particlewalldis = Teuchos::rcp(new DRT::Discretization(discret_name,com));

  std::vector<int> nodeids;
  std::vector<int> eleids;
  // loop over all particle wall nodes and elements and fill new discretization
  for(std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit=structgelements.begin(); meit!=structgelements.end(); ++meit)
  {
    // care about particle wall nodes:
    // fill everything in case of static walls and only row nodes in case of moving walls
    std::map<int, DRT::Node*> wallgnodes = structgnodes[meit->first];
    for (std::map<int, DRT::Node* >::iterator nit=wallgnodes.begin(); nit != wallgnodes.end(); ++nit)
    {
      DRT::Node* currnode = (*nit).second;
      if ( not moving_walls_ || ((currnode->Owner() == MyRank()) && moving_walls_) )
      {
        nodeids.push_back(currnode->Id());
        particlewalldis->AddNode(Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
      }
    }

    // care about particle wall eles:
    // fill everything in case of static walls and only row elements in case of moving walls
    std::map<int, Teuchos::RCP<DRT::Element> > structelementsinterf = structgelements[meit->first];
    for (std::map<int, Teuchos::RCP<DRT::Element> >::iterator eit=structelementsinterf.begin(); eit != structelementsinterf.end(); ++eit)
    {
      Teuchos::RCP<DRT::Element> currele = eit->second;
      if ( not moving_walls_ || ((currele->Owner() == MyRank()) && moving_walls_) )
      {
        eleids.push_back(currele->Id() );
        // structural surface elements cannot be distributed --> Bele3 element is used
        Teuchos::RCP<DRT::Element> wallele = DRT::UTILS::Factory(elename,"Polynomial", currele->Id(), currele->Owner());
        wallele->SetNodeIds(currele->NumNode(), currele->NodeIds());
        particlewalldis->AddElement( wallele );
      }
    }
  }

  // extended ghosting of wall elements for static walls
  if(not moving_walls_)
  {
    // fill complete in order to obtain element col map
    particlewalldis->FillComplete(false, false, false);

    particlewallelecolmap_standardghosting_ = Teuchos::rcp(new Epetra_Map(*particlewalldis->ElementColMap()));

    // extend ghosting to the bin col map
    BinStrategy()->ExtendEleGhosting(particlewalldis, particlewallelecolmap_standardghosting_, BinColMap(), false, false, false);
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
  if(MyRank() == 0)
    std::cout << "after adding particle walls" << std::endl;
  DRT::UTILS::PrintParallelDistribution(*particlewalldis_);

  if(moving_walls_)
    wallextractor_ = Teuchos::rcp(new LINALG::MapExtractor(*(basediscret->DofRowMap()),Teuchos::rcp(particlewalldis_->DofRowMap(), false)));

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

  // circumcircle of bin
  double bincircumcircle = 0.0;
  for(int dim=0; dim<3; ++dim)
  {
    bincircumcircle += std::pow(BinStrategy()->BinSize()[dim]/2.0,2.0);
  }
  bincircumcircle = sqrt(bincircumcircle);

  // minimal bin size
  double min_bin_size = BinStrategy()->BinSize()[0];
  for(int dim=1; dim<3; ++dim)
    min_bin_size = std::min(min_bin_size, BinStrategy()->BinSize()[dim]);

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
      BinStrategy()->ConvertPosToijk(currentpositions[nodeids[0]], ijk);

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; ++j)
      {
        int ijk[3];
        BinStrategy()->ConvertPosToijk(currentpositions[nodeids[j]], ijk);

        for(int dim=0; dim<3; ++dim)
        {
          if(ijk[dim]<ijk_range[dim*2])
            ijk_range[dim*2] = ijk[dim];
          if(ijk[dim]>ijk_range[dim*2+1])
            ijk_range[dim*2+1] = ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range and fill them into binIds
      BinStrategy()->GidsInijkRange(&ijk_range[0], binIds, false);
    } // do a positive search

    // if no bins on this proc were found, next wall element can be processed
    if ( binIds.empty() )
      continue;

    // do a negative search and remove bins that are too far away from the wall element
    {
      std::set<int> binfaraway;
      for ( std::set<int>::const_iterator biniter = binIds.begin(); biniter != binIds.end(); ++biniter )
      {
        const LINALG::Matrix<3,1> bincentroid = BinStrategy()->GetBinCentroid(*biniter);

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
          BinStrategy()->GetBinCorners(*biniter, bincorners);

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

            const int gid = BinStrategy()->ConvertPosToGid(minDistCoords);
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

      // erase bins that are far away from a wall element
      for ( std::set<int>::const_iterator biniter = binfaraway.begin(); biniter != binfaraway.end(); ++biniter )
        binIds.erase(*biniter);
    } // do a negative search

    // if none of found bins is in BinColMap, next wall element can be processed
    {
      bool noBinInBinColMap = true;
      for ( std::set<int>::const_iterator biniter = binIds.begin(); biniter != binIds.end(); ++biniter )
        if ( BinColMap()->LID(*biniter) >= 0 )
        {
          noBinInBinColMap = false;
          continue;
        }

      if ( noBinInBinColMap )
        continue;
    }

    // insert found bins into map relating wall gids to bin ids
    relwallgidtobinids_[wallele->Id()].insert(binIds.begin(), binIds.end() );

  } // end lid

  return;
}

/*----------------------------------------------------------------------*
 | single fields are tested                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(particles_->CreateFieldTest());

  if (GetRendering() != Teuchos::null)
    DRT::Problem::Instance()->AddFieldTest(GetRendering()->CreateFieldTest());

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

  if (writeresultsevery_ and (Step()%writeresultsevery_ == 0))
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
 | update connectivity                                    sfuchs 08/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::UpdateConnectivity()
{
  // relate wall gids to bins ids
  if ( moving_walls_ )
    RelateWallGidsToBinIds();

  // transfer particles into their correct bins
  if ( Step()%transfer_every_ == 0 )
    TransferParticles(true);

  // assign wall elements and gids to bins (after ghosting)
  if ( moving_walls_ )
    AssignWallElesAndGidsToBins();

  return;
}

/*----------------------------------------------------------------------*
 | get neighbouring particles and walls                    ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetNeighbouringItems(
    DRT::Node* particle,
    std::list<DRT::Node*>& neighboursLinf_p,
    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w) const
{
  if (particle->NumElement() != 1)
    dserror("More than one element for this particle");

  DRT::Element** CurrentBin = particle->Elements();

  GetNeighbouringItems(CurrentBin[0]->Id(),neighboursLinf_p,&neighboursLinf_w);

  return;
}

/*----------------------------------------------------------------------*
 | get neighbouring particles and walls (bin version)      katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetNeighbouringItems(
    const int binId,
    std::list<DRT::Node*>& neighboursLinf_p,
    boost::unordered_map<int, DRT::Element*>* neighboursLinf_w) const
{
  std::vector<int> binIds;
  binIds.reserve(27);

  BinStrategy().GetNeighborAndOwnBinIds(binId,binIds);

  GetBinContent(neighboursLinf_p, neighboursLinf_w, binIds);

  return;
}

/*----------------------------------------------------------------------*
 | get particles and wall elements in given bins           ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetBinContent(
    std::list<DRT::Node*>& bin_p,
    boost::unordered_map<int, DRT::Element*>* bin_w,
    std::vector<int> &binIds) const
{
  // loop over all bins
  for(std::vector<int>::const_iterator bin=binIds.begin(); bin!=binIds.end(); ++bin)
  {
    // extract bins from discretization after checking on existence
    const int lid = BinStrategy().BinDiscret()->ElementColMap()->LID(*bin);
    if(lid<0)
      continue;

#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy().BinDiscret()->lColElement(lid));
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* neighboringbin =
        static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy().BinDiscret()->lColElement(lid));

    // gather particles
    DRT::Node** nodes = neighboringbin->Nodes();
    bin_p.insert(bin_p.end(), nodes, nodes+neighboringbin->NumNode());

    // gather wall elements
    if (bin_w != NULL)
    {
      DRT::Element** walleles = neighboringbin->AssociatedEles(bin_wallcontent_);
      const int numwalls = neighboringbin->NumAssociatedEle(bin_wallcontent_);
      const int* ids=neighboringbin->AssociatedEleIds(bin_wallcontent_);
      for(int iwall=0;iwall<numwalls; ++iwall)
      {
        if(ids[iwall]!=walleles[iwall]->Id())
        {
          std::cout << "Comm().MyPID(): " << Comm().MyPID() << std::endl;
          std::cout << "ids[iwall]: " << ids[iwall] << std::endl;
          std::cout << "walleles[iwall]->Id(): " << walleles[iwall]->Id() << std::endl << std::endl;
          dserror("Ids of wall elements are different from Ids stored in multibin! Accidential access of random storage?");
        }

        (*bin_w)[walleles[iwall]->Id()] = walleles[iwall];
      }
    }
  }
  return;
}

/*------------------------------------------------------------------------*
 | return gravity acceleration                               meier 05/17  |
 *------------------------------------------------------------------------*/
LINALG::Matrix<3,1> PARTICLE::Algorithm::GetGravityAcc(const double time)
{
  double fac=1.0;
  double gravity_ramp_time=DRT::Problem::Instance()->ParticleParams().get<double>("GRAVITY_RAMP_TIME");
  if(gravity_ramp_time>0 and time >=0)
  {
    if(time<gravity_ramp_time)
       fac=0.5*(1-cos(time*PI/gravity_ramp_time));
  }

  LINALG::Matrix<3,1> scaled_gravity_acc(gravity_acc_);
  scaled_gravity_acc.Scale(fac);

  return scaled_gravity_acc;
}
