/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.cpp

\brief Algorithm to control particle simulations

\level 2

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"
#include "particle_utils.H"
#include "particleMeshFree_rendering.H"
#include "particleMeshFree_interaction.H"
#include "../drt_adapter/adapter_particle.H"
#include "binning_strategy.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_inpar/inpar_meshfree.H"
#include "../drt_inpar/inpar_particle.H"

#include "../drt_mat/particle_mat.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_mat/matpar_bundle.H"

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

#include "particle_heatSource.H"
#include "particle_node.H"

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
  particleMat_(NULL),
  extParticleMat_(NULL),
  bin_wallcontent_(BINSTRATEGY::UTILS::BELE3),
  rendering_(Teuchos::null)
{
  const Teuchos::ParameterList& meshfreeparams = DRT::Problem::Instance()->MeshfreeParams();
  // safety check
  INPAR::MESHFREE::meshfreetype meshfreetype = DRT::INPUT::IntegralValue<INPAR::MESHFREE::meshfreetype>(meshfreeparams,"TYPE");
  if (meshfreetype!=INPAR::MESHFREE::particle)
    dserror("MESHFREE -> TYPE must be Particle in input file.");

  if (moving_walls_ == true and DRT::Problem::Instance()->ProblemType() != prb_pasi)
    dserror("MOVING_WALLS flag is activated!\n"
        "Set parameter PROBLEMTYP to 'Particle_Structure_Interaction' in ---PROBLEM TYP section.\n"
        "Set parameter COUPALGO to 'partitioned_onewaycoup' in ---PASI DYNAMIC section\n"
        "and set the particle structure interaction time parameters.");

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

    // adaptions for Normal_DEM_thermo
    NormDemThermoAdapt();

    // calculate stresses, strains, energies
    PrepareOutput();

    // transfer particles and heat sources into their correct bins
    UpdateConnectivity();

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
  BinStrategy()->BinDiscret()->FillComplete(false,false,false);

  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*BinStrategy()->BinDiscret()->NodeRowMap()));

  Teuchos::RCP<Epetra_Map> binrowmap;
  if(not restarted)
  {
    BinStrategy()->CreateBins(BinStrategy()->BinDiscret());
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

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // set up the links to the materials for easy access
    // make sure that a particle material is defined in the dat-file
    InitMaterials();

    // get input parameters for particles
    const Teuchos::ParameterList& particledyn = DRT::Problem::Instance()->ParticleParams();

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

    // set up Heat Sources in a map
    SetUpHeatSources();
    UpdateHeatSourcesConnectivity(false);

    // access structure and build particle walls
    AccessStructure();
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
  if (MyRank() == 0)
    IO::cout << "after ghosting of particles" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*BinStrategy()->BinDiscret());

  // update connectivity
  //UpdateHeatSourcesConnectivity(true);

  if (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->ParticleParams(),"RENDERING"))
  {
    rendering_ = Teuchos::rcp(new PARTICLE::Rendering(Teuchos::rcp(this,false)));
  }
}

/*----------------------------------------------------------------------*
| build connectivity from particle wall elements to bins    ghamm 04/13 |
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
 | set up pointers to material bundles                     catta 09/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::InitMaterials()
{
  int id = -1;
  switch (particleInteractionType_)
  {
  case INPAR::PARTICLE::MeshFree :
  case INPAR::PARTICLE::Normal_DEM_thermo :
  {
    id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_extparticlemat);
    break;
  }
  default :
  {
    id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    break;
  }
  }
  // check
  if (id==-1)
    dserror("Could not find particle material or material type");
  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  particleMat_ = static_cast<const MAT::PAR::ParticleMat*>(mat);
  if (particleInteractionType_ == INPAR::PARTICLE::MeshFree || particleInteractionType_ == INPAR::PARTICLE::Normal_DEM_thermo)
    extParticleMat_ = static_cast<const MAT::PAR::ExtParticleMat*>(mat);
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

  UpdateHeatSourcesConnectivity(true);

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

  // update walls and connectivity
  if(particlewalldis_ != Teuchos::null)
  {
    if( not moving_walls_)
      BinStrategy()->ExtendGhosting(particlewalldis_, particlewallelecolmap_standardghosting_, BinColMap(), true, false, false, false);

    BuildElementToBinPointers(true);
  }

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  // restart heat source map
  UpdateHeatSourcesConnectivity(true);

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
    case INPAR::PARTICLE::MeshFree :
    case INPAR::PARTICLE::Normal_DEM_thermo :
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

    if(maxrad + transfer_every_*maxvel*dt > 0.5*BinStrategy()->CutoffRadius())
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
            "Max radius is: %f", (int)gid/3, 2.0*(maxrad + maxvel*dt), BinStrategy()->CutoffRadius(), maxrad);
      }
    }
  }
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
    BinStrategy()->BinDiscret()->ExtendedGhosting(*BinColMap(),true,false,true,false);
  else
    BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
  bool rebuildwallpointer = true;
  if(moving_walls_)
    rebuildwallpointer = false;
  BuildElementToBinPointers(rebuildwallpointer);

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
 | set structural states                                   sfuchs 02/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetStructStates(
    Teuchos::RCP<const Epetra_Vector> structdispn,   ///< structural displacements at t_n
    Teuchos::RCP<const Epetra_Vector> structdispnp,  ///< structural displacements at t_n+1
    Teuchos::RCP<const Epetra_Vector> structvelnp    ///< structural velocities at t_n
    )
{
  structdispn_ = structdispn;
  structdispnp_ = structdispnp;
  structvelnp_ = structvelnp;

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
    if(!moving_walls_)
      AssignWallElesToBins();
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

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit;

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
    BinStrategy()->ExtendGhosting(particlewalldis, particlewallelecolmap_standardghosting_, BinColMap(), false, false, false, false);
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

  return;
}


/*----------------------------------------------------------------------*
| particle walls are added from the structural discret      ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::AssignWallElesToBins()
{
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
      BinStrategy()->GidsInijkRange(&ijk_range[0], binIds, true);
    }

    // if no bins on this proc were found, next wall element can be processed
    if(binIds.empty())
      continue;

    // do a negative search and remove bins that are too far away from the wall element
    {
      std::set<int> binfaraway;
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
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
      for(std::set<int>::const_iterator biniter=binfaraway.begin(); biniter!=binfaraway.end(); ++biniter)
        binIds.erase(*biniter);
    }

    // assign wall element to remaining bins
    {
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy()->BinDiscret()->gElement(*biniter))->AddAssociatedEle(bin_wallcontent_,wallele->Id(), wallele);
    }

  } // end lid

  return;
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

  return;
}

/*----------------------------------------------------------------------*
 | adaptions for Normal_DEM_thermo                        sfuchs 02/17  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::NormDemThermoAdapt()
{
  if (particleInteractionType_ == INPAR::PARTICLE::Normal_DEM_thermo)
  {
    // compute thermodynamic expansion
    ThermalExpansion();

    // particle dismembering
    ParticleDismemberer();
  }

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

  if(moving_walls_ and writeresultsevery_ and (Step()%writeresultsevery_ == 0))
  {
    Teuchos::RCP<Epetra_Vector> walldisnp = wallextractor_->ExtractCondVector(structdispnp_);
    particlewalldis_->Writer()->NewStep(Step(), Time());
    particlewalldis_->Writer()->WriteVector("displacement", walldisnp);
  }

//  const std::string filename = IO::GMSH::GetFileName("particle_data", Step(), true, Comm().MyPID());
//  std::ofstream gmshfilecontent(filename.c_str());
//
//  // velocity
//  {
//    gmshfilecontent << "View \" " << "velocity" << " \" {\n";
//    LINALG::Matrix<3,1> vectorvalue(true);
//
//    for(int n=0; n<BinStrategy()->BinDiscret()->NumMyRowNodes(); n++)
//    {
//      DRT::Node* actnode = BinStrategy()->BinDiscret()->lRowNode(n);
//      // get the first gid of a node and convert it into a LID
//      int gid = BinStrategy()->BinDiscret()->Dof(actnode, 0);
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
//    for(int n=0; n<BinStrategy()->BinDiscret()->NumMyRowNodes(); n++)
//    {
//      DRT::Node* actnode = BinStrategy()->BinDiscret()->lRowNode(n);
//      // get the first gid of a node and convert it into a LID
//      int gid = BinStrategy()->BinDiscret()->Dof(actnode, 0);
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
//    for(int n=0; n<BinStrategy()->BinDiscret()->NumMyRowNodes(); n++)
//    {
//      DRT::Node* actnode = BinStrategy()->BinDiscret()->lRowNode(n);
//      // get the first gid of a node and convert it into a LID
//      int gid = BinStrategy()->BinDiscret()->Dof(actnode, 0);
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

/*----------------------------------------------------------------------*
 | set up heat sources                                      catta 06/16 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetUpHeatSources()
{
  // extract heat source conditions
  std::vector<DRT::Condition*> conds;
  BinStrategy()->BinDiscret()->GetCondition("ParticleHeatSource", conds);

  for (size_t iHS=0; iHS<conds.size(); iHS++)
  {
    // extract condition
    std::vector<double> HSZone_minVer = *conds[iHS]->Get<std::vector<double> >("vertex0");
    std::vector<double> HSZone_maxVer = *conds[iHS]->Get<std::vector<double> >("vertex1");
    const double HSQDot = conds[iHS]->GetDouble("HSQDot");
    const double HSTstart = conds[iHS]->GetDouble("HSTstart");
    const double HSTend = conds[iHS]->GetDouble("HSTend");

    // vertex sort
    double buffer;
    for (int idim=0; idim<3; idim++)
    {
      if (HSZone_maxVer[idim]<HSZone_minVer[idim])
        {
          buffer = HSZone_maxVer[idim];
          HSZone_maxVer[idim] = HSZone_minVer[idim];
          HSZone_minVer[idim] = buffer;
        }
    }

    // create the heat source class
    heatSources_.push_back(Teuchos::rcp(new HeatSource(false,
                                                      iHS,
                                                      HSZone_minVer,
                                                      HSZone_maxVer,
                                                      HSQDot,
                                                      HSTstart,
                                                      HSTend)));
  }
}


/*----------------------------------------------------------------------*
 | update of the map bins->heat Sources                     catta 06/16 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::UpdateHeatSourcesConnectivity(bool trg_forceRestart)
{
  // clear the map in case of force restart
  if (trg_forceRestart)
    bins2heatSources_.clear();

  // heat source activation
  for (std::list<Teuchos::RCP<HeatSource> >::const_iterator iHS = heatSources_.begin(); iHS != heatSources_.end(); ++iHS)
  {
    // force restart heatSources status
    if (trg_forceRestart)
      (*iHS)->active_ = false;

    // map assignment
    if ((*iHS)->Tstart_<=Time() && (*iHS)->Tend_>=Time() && (*iHS)->active_ == false)
    {
      // find the bins
      int minVerZone_ijk[3];
      int maxVerZone_ijk[3];
      BinStrategy()->ConvertPosToijk(&((*iHS)->minVerZone_[0]), minVerZone_ijk);
      BinStrategy()->ConvertPosToijk(&((*iHS)->maxVerZone_[0]), maxVerZone_ijk);
      const int ijk_range[] = {
          minVerZone_ijk[0],maxVerZone_ijk[0],
          minVerZone_ijk[1],maxVerZone_ijk[1],
          minVerZone_ijk[2],maxVerZone_ijk[2]};
      std::set<int>  binIds;
      BinStrategy()->GidsInijkRange(&ijk_range[0], binIds, false);

      if (binIds.empty()) dserror("Weird! Heat Source %i found but could not be assigned to bins. Is it outside of bins?",(*iHS)->id_);

      // create/update the map
      for (std::set<int>::const_iterator iBin = binIds.begin(); iBin != binIds.end(); ++iBin)
          if(BinStrategy()->BinDiscret()->ElementRowMap()->LID(*iBin) >= 0)
            bins2heatSources_[*iBin].push_back((*iHS));

      (*iHS)->active_ = true;
    }
  }

  // heat source deactivation
  for (std::list<Teuchos::RCP<HeatSource> >::const_iterator iHS = heatSources_.begin(); iHS != heatSources_.end(); ++iHS)
  {
    if (((*iHS)->Tend_<Time() || (*iHS)->Tstart_>Time()) && (*iHS)->active_ == true)
    {
      // remove elements from the map
      for (std::map<int,std::list<Teuchos::RCP<HeatSource> > >::iterator iBin = bins2heatSources_.begin(); iBin != bins2heatSources_.end(); ++iBin)
      {
        for (std::list<Teuchos::RCP<HeatSource> >::iterator iHSb=iBin->second.begin(); iHSb != iBin->second.end(); ++iHSb)
        {
          if (*iHS == *iHSb)
          {
            iHSb = iBin->second.erase(iHSb);
            break;
          }
        }
      }
      (*iHS)->active_ = false;
    }
  }
}


/*----------------------------------------------------------------------*
 | update connectivity                                     catta 06/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::UpdateConnectivity()
{
  if(Step()%transfer_every_ == 0)
  {
    // transfer particles into their correct bins
    TransferParticles(true);
  }
  // update heat sources
  UpdateHeatSourcesConnectivity(false);
}


/*----------------------------------------------------------------------*
 | get neighbouring particles and walls                    ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetNeighbouringItems(
    DRT::Node* particle,
    std::list<DRT::Node*>& neighboursLinf_p,
    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs)
{
  if (particle->NumElement() != 1)
    dserror("More than one element for this particle");

  DRT::Element** CurrentBin = particle->Elements();

  GetNeighbouringItems(CurrentBin[0]->Id(),neighboursLinf_p,neighboursLinf_w, neighboursLinf_hs);
}


/*----------------------------------------------------------------------*
 | get neighbouring particles and walls (bin version)      katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetNeighbouringItems(
    const int binId,
    std::list<DRT::Node*>& neighboursLinf_p,
    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs)
{
  std::vector<int> binIds;
  binIds.reserve(27);

  BinStrategy()->GetNeighborAndOwnBinIds(binId,binIds);

  GetBinContent(neighboursLinf_p, neighboursLinf_w, neighboursLinf_hs, binIds);
}


/*----------------------------------------------------------------------*
 | get particles and wall elements in given bins           ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetBinContent(
    std::list<DRT::Node*>& bin_p,
    boost::unordered_map<int, DRT::Element*>& bin_w,
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > bin_hs,
  std::vector<int> &binIds)
{
  // loop over all bins
  for(std::vector<int>::const_iterator bin=binIds.begin(); bin!=binIds.end(); ++bin)
  {
    // extract bins from discretization after checking on existence
    const int lid = BinStrategy()->BinDiscret()->ElementColMap()->LID(*bin);
    if(lid<0)
      continue;

#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy()->BinDiscret()->lColElement(lid));
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* neighboringbin =
        static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy()->BinDiscret()->lColElement(lid));

    // gather particles
    DRT::Node** nodes = neighboringbin->Nodes();
    bin_p.insert(bin_p.end(), nodes, nodes+neighboringbin->NumNode());

    // gather wall elements
    DRT::Element** walleles = neighboringbin->AssociatedEles(bin_wallcontent_);
    const int numwalls = neighboringbin->NumAssociatedEle(bin_wallcontent_);
    for(int iwall=0;iwall<numwalls; ++iwall)
    {
      bin_w[walleles[iwall]->Id()] = walleles[iwall];
    }
    // gather heat sources
    // it is a set so that there are no repetitions of the heat source
    if (bin_hs != Teuchos::null)
    {
      std::list<Teuchos::RCP<HeatSource> >::const_iterator iHS;
      for (iHS = bins2heatSources_[*bin].begin(); iHS != bins2heatSources_[*bin].end(); ++iHS)
        (*bin_hs)[(*iHS)->Id()] = (*iHS);
    }
  }
}

/*----------------------------------------------------------------------*
 | particleDismemberer                                     catta 07/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ParticleDismemberer()
{
  // extract the material parameters
  const double dismemberRadius = extParticleMat_->dismemberRadius_;
  const double specEnthalpyTL = extParticleMat_->SpecEnthalpyTL();

  // extract the state vectors
  Teuchos::RCP<Epetra_Vector> dispn = particles_->WriteAccessDispnp();
  Teuchos::RCP<Epetra_Vector> mass = particles_->WriteAccessMass();
  Teuchos::RCP<Epetra_Vector> radiusn = particles_->WriteAccessRadiusnp();
  Teuchos::RCP<Epetra_Vector> specEnthalpyn = particles_->WriteAccessSpecEnthalpynp();
  Teuchos::RCP<const Epetra_Vector> specEnthalpy = particles_->SpecEnthalpyn();

  // with this snapshotting the addition of nodes does not affect the main loop
  const int maxLidNode_old = BinStrategy()->BinDiscret()->NodeRowMap()->NumMyElements();

  std::vector<int> listOrganizer(maxLidNode_old);
  std::list<homelessParticleTemp > newParticleList;

  for (int lidNode_old = 0; lidNode_old < maxLidNode_old; ++lidNode_old)
  {
    DRT::Node* currParticle_old = BinStrategy()->BinDiscret()->lRowNode(lidNode_old);
    const int lidDof_old = BinStrategy()->BinDiscret()->DofRowMap()->LID(BinStrategy()->BinDiscret()->Dof(currParticle_old, 0));

    // checks
    if (lidNode_old == -1)
      dserror("Invalid lidNode\n");
    if (lidDof_old == -1)
      dserror("Invalid lidDof\n");
    if (dismemberRadius > ((*radiusn)[lidNode_old]/3))
      continue;  // we do not compress particles without creating new particles

    // check if in this time step we completed the transition. This check is based on the temperature history
    if (Step() > 0 && (*specEnthalpy)[lidNode_old] <= specEnthalpyTL && (*specEnthalpyn)[lidNode_old]>specEnthalpyTL)
    {
      // --------------------------------------------------------------------
      // position of the new particles temporarily stocked to compute the gid
      // --------------------------------------------------------------------
      const double x_step = 2 * dismemberRadius ;
      const int semiLengthInParticlesx = ComputeSemiLengthInParticlesForParticleDismemberer((*radiusn)[lidNode_old], x_step/2.0);
      const double y_step = sqrt(3) * dismemberRadius;
      const int semiLengthInParticlesy = ComputeSemiLengthInParticlesForParticleDismemberer((*radiusn)[lidNode_old], y_step/2.0);
      const double z_step = 2 * sqrt(2) * (1/(sqrt(3))) * dismemberRadius;
      const int semiLengthInParticlesz = ComputeSemiLengthInParticlesForParticleDismemberer((*radiusn)[lidNode_old], z_step/2.0);
      LINALG::Matrix<3,1> newRelativeParticlePosition(true);
      for (int ix=-semiLengthInParticlesx; ix<=semiLengthInParticlesx; ++ix)
      {
        for (int iy=-semiLengthInParticlesy; iy<=semiLengthInParticlesy; ++iy)
        {
          for (int iz=-semiLengthInParticlesz; iz<=semiLengthInParticlesz; ++iz)
          {
            // x position
            newRelativeParticlePosition(0) = ix * x_step;
            // line offset if it is odd in the y direction
            if (std::abs(iy)%2 == 1)
            {
              newRelativeParticlePosition(0) += x_step/2.0;
            }

            // y position
            newRelativeParticlePosition(1) = iy * y_step;
            // plane offset if it is odd in the z direction
            if (std::abs(iz)%2 == 1)
            {
              newRelativeParticlePosition(0) -= x_step/2.0;
              newRelativeParticlePosition(1) += y_step/3.0;
            }
            // z position
            newRelativeParticlePosition(2) = iz * z_step;
            // is it inside the old radius?
            if ((newRelativeParticlePosition.Norm2()+dismemberRadius <= (*radiusn)[lidNode_old]) && !(ix == 0 && iy == 0 && iz == 0))
            {
              homelessParticleTemp hpt;
              std::vector<double> newParticlePosition(3);
              for (int ii=0; ii<3; ++ii)
              {
                newParticlePosition.at(ii) = newRelativeParticlePosition(ii) + (*dispn)[lidDof_old+ ii];
              }
              hpt.pos = newParticlePosition;
              hpt.lidNode_old = lidNode_old;
              hpt.lidDof_old = lidDof_old;
              newParticleList.push_back(hpt);
              ++listOrganizer[lidNode_old];
            }
          }
        }
      }
    }
  }

  // -----------------------------------------------
  // communication: determination of the correct gid
  // -----------------------------------------------
  const int nextMaxNode_gid = BinStrategy()->BinDiscret()->NodeRowMap()->MaxAllGID();
  const int numproc = BinStrategy()->BinDiscret()->Comm().NumProc();
  std::vector<int> myentries(numproc,0);
  std::vector<int> globentries(numproc,0);
  const int MyNewParticleListSize0 = newParticleList.size();
  myentries[BinStrategy()->BinDiscret()->Comm().MyPID()] = MyNewParticleListSize0;
  BinStrategy()->BinDiscret()->Comm().SumAll(&myentries[0], &globentries[0], numproc);  // parallel communication
  int MyOffset = 0;
  for(int ii = 0; ii < BinStrategy()->BinDiscret()->Comm().MyPID(); ++ii)
  {
   MyOffset += globentries[ii];
  }
  // finally the offset for the gid
  const int finalOffset = nextMaxNode_gid + MyOffset;

  // --------------------------------------------------
  // place the node in the bins and update connectivity
  // --------------------------------------------------
  int newParticleIDcounter = finalOffset;
  for (std::list<homelessParticleTemp >::const_iterator iNodeList = newParticleList.begin(); iNodeList != newParticleList.end(); ++iNodeList)
  {
    ++newParticleIDcounter;

    std::list<Teuchos::RCP<DRT::Node> > homelessparticles;
    Teuchos::RCP<DRT::Node> newParticle = Teuchos::rcp(new PARTICLE::ParticleNode(
        newParticleIDcounter, &((*iNodeList).pos[0]), MyRank()));

    PlaceNodeCorrectly(newParticle, newParticle->X(), homelessparticles);
    // rare case when after dismembering some particles fall into another bin
    if(homelessparticles.size() != 0)
    {
      BinStrategy()->BinDiscret()->AddNode(newParticle);
      // assign node to an arbitrary row bin -> correct placement will follow in the timeloop in TransferParticles
      DRT::MESHFREE::MeshfreeMultiBin* firstbinindis = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(BinStrategy()->BinDiscret()->lRowElement(0));
      firstbinindis->AddNode(newParticle.get());
      homelessparticles.clear();
    }
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();
  UpdateStates();

  // reset/set of the pointers after updatestestesafterparticletransfer
  dispn = particles_->WriteAccessDispnp();
  mass = particles_->WriteAccessMass();
  radiusn = particles_->WriteAccessRadiusnp();
  specEnthalpyn = particles_->WriteAccessSpecEnthalpynp();
  Teuchos::RCP<Epetra_Vector> densityn = particles_->WriteAccessDensitynp();
  Teuchos::RCP<Epetra_Vector> veln = particles_->WriteAccessVelnp();
  Teuchos::RCP<Epetra_Vector> accn = particles_->WriteAccessAccnp();
  Teuchos::RCP<Epetra_Vector> inertia = particles_->WriteAccessInertia();

  int lidNodeCounter = 0;
  for (std::list<homelessParticleTemp >::const_iterator iNodeList = newParticleList.begin(); iNodeList != newParticleList.end(); ++iNodeList)
  {
    // get node lids
    const int lidNode_new = maxLidNode_old+lidNodeCounter;
    DRT::Node* currParticle_new = BinStrategy()->BinDiscret()->lRowNode(lidNode_new);
    const int lidDof_new = BinStrategy()->BinDiscret()->DofRowMap()->LID(BinStrategy()->BinDiscret()->Dof(currParticle_new, 0));

    const int lidNode_old = iNodeList->lidNode_old;
    const int lidDof_old = iNodeList->lidDof_old;

    // check
    if (lidNode_new == -1)
      dserror("invalid new node lid");
    if (lidDof_new == -1)
      dserror("invalid new dof lid");
    if (lidNode_old == -1)
      dserror("invalid old node lid");
    if (lidDof_old == -1)
      dserror("invalid old dof lid");

    // new masses and densities (to conserve the overall mass)
    MassDensityUpdaterForParticleDismemberer(mass, densityn, radiusn, lidNode_new, lidNode_old, listOrganizer[lidNode_old]);
    (*radiusn)[lidNode_new] = dismemberRadius;
    (*specEnthalpyn)[lidNode_new] = (*specEnthalpyn)[lidNode_old];
    (*inertia)[lidNode_new] = PARTICLE::Utils::ComputeInertia((*radiusn)[lidNode_new], (*mass)[lidNode_new]);

    for(int d=0; d<3; ++d)
    {
      (*dispn)[lidDof_new + d] = currParticle_new->X()[d];
      (*veln)[lidDof_new + d] = (*veln)[lidDof_old + d];
      (*accn)[lidDof_new + d] = (*accn)[lidDof_old + d];
    }

    ++lidNodeCounter;
  }

  // ------------------------------
  // update of the central particle
  // ------------------------------
  for (int lidNode_old = 0; lidNode_old < maxLidNode_old; ++lidNode_old)
  {
    if (Step() > 0 && (*specEnthalpy)[lidNode_old] <= specEnthalpyTL &&
        (*specEnthalpyn)[lidNode_old]>specEnthalpyTL && dismemberRadius <= ((*radiusn)[lidNode_old]/3))
    {
      MassDensityUpdaterForParticleDismemberer(mass, densityn, radiusn, lidNode_old, lidNode_old, listOrganizer[lidNode_old]);
      // radius MUST be updated after MassDensityUpdaterForParticleDismemberer
      (*radiusn)[lidNode_old] = dismemberRadius;
      (*inertia)[lidNode_old] = PARTICLE::Utils::ComputeInertia((*radiusn)[lidNode_old], (*mass)[lidNode_old]);
    }
  }
}


/*----------------------------------------------------------------------*
 | ComputeSemiLengthInParticlesForParticleDismemberer      catta 06/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::MassDensityUpdaterForParticleDismemberer(
    Teuchos::RCP<Epetra_Vector> &mass,
    Teuchos::RCP<Epetra_Vector> &densitynp,
    Teuchos::RCP<Epetra_Vector> &radius,
    const int &lidNode_new,
    const int &lidNode_old,
    const int &nlist)
{
  const double dismemberRadius = extParticleMat_->dismemberRadius_;
  // new masses and densities (to conserve the overall mass)
  (*mass)[lidNode_new] = ((*mass)[lidNode_old])/(nlist+1); // the +1 is due to the central node that is resized
  const double radiusOld = (*radius)[lidNode_old];
  (*densitynp)[lidNode_new] = (*densitynp)[lidNode_old] * radiusOld * radiusOld * radiusOld /((nlist + 1) * dismemberRadius * dismemberRadius * dismemberRadius);
}

/*------------------------------------------------------------------------*
 | compute thermodynamic expansion - new densities and radii catta 06/16  |
 *------------------------------------------------------------------------*/
void PARTICLE::Algorithm::ThermalExpansion()
{
  // extract the interesting state vectors
  Teuchos::RCP<const Epetra_Vector> mass = particles_->Mass();
  Teuchos::RCP<const Epetra_Vector> specEnthalpy = particles_->SpecEnthalpyn();
  Teuchos::RCP<const Epetra_Vector> specEnthalpyn = particles_->SpecEnthalpynp();
  Teuchos::RCP<Epetra_Vector> densityn = particles_->WriteAccessDensitynp();
  Teuchos::RCP<Epetra_Vector> radiusn = particles_->WriteAccessRadiusnp();

  // extract the material parameters for easy access
  const double specEnthalpyST = extParticleMat_->SpecEnthalpyST();
  const double specEnthalpyTL = extParticleMat_->SpecEnthalpyTL();
  const double CPS = extParticleMat_->CPS_;
  const double inv_CPS = 1/CPS;
  const double CPL = extParticleMat_->CPL_;
  const double inv_CPL = 1/CPL;
  const double latentHeat = extParticleMat_->latentHeat_;
  const double thermalExpansionS = extParticleMat_->thermalExpansionS_;
  const double thermalExpansionL = extParticleMat_->thermalExpansionL_;
  const double thermalExpansionT = extParticleMat_->thermalExpansionT_;



  // update the other state vectors (\rho and R)
  for (int lidNode = 0; lidNode < mass->MyLength(); ++lidNode)
  {
    const double oldSpecEnthalpy = (*specEnthalpy)[lidNode];
    const double newSpecEnthalpy = (*specEnthalpyn)[lidNode];
    // skip in case the specEnthalpy did not change
    if (newSpecEnthalpy != oldSpecEnthalpy)
    {
      // compute the current volume
      double volume = PARTICLE::Utils::Radius2Volume((*radiusn)[lidNode]);
      // specEnthalpy difference
      double deltaSpecEnthalpy = newSpecEnthalpy - oldSpecEnthalpy;

      // --- compute the new volume --- //

      // WAS it solid?
      if (oldSpecEnthalpy <= specEnthalpyST)
      {
        // IS it solid?
        if (newSpecEnthalpy <= specEnthalpyST)
        {
          volume *= inv_CPS * thermalExpansionS * deltaSpecEnthalpy + 1;
        }
        // IS it liquid?
        else if (newSpecEnthalpy >= specEnthalpyTL)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyST - oldSpecEnthalpy;

          // expansion in solid state
          volume *= inv_CPS * thermalExpansionS * deltaSpecEnthalpyUpToTransition + 1;
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in transition state
          volume *= thermalExpansionT * latentHeat + 1;
          deltaSpecEnthalpy -= latentHeat;

          // expansion in liquid state
          volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpy + 1;
        }
        // it IS transition state
        else
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyST - oldSpecEnthalpy;

          // expansion in solid state
          volume *= inv_CPS * thermalExpansionS * deltaSpecEnthalpyUpToTransition + 1;
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in transition state
          volume *= thermalExpansionT * deltaSpecEnthalpy + 1;
        }
      }
      // WAS it liquid?
      else if (oldSpecEnthalpy >= specEnthalpyTL)
      {
        // IS it solid?
        if (newSpecEnthalpy <= specEnthalpyST)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyTL - oldSpecEnthalpy;

          // expansion in liquid state
          volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpyUpToTransition + 1;
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in the transition state
          volume *= thermalExpansionT *(- latentHeat) + 1;
          deltaSpecEnthalpy -= -latentHeat;

          // expansion in liquid state
          volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpy + 1;
        }
        // IS it liquid?
        else if (newSpecEnthalpy >= specEnthalpyTL)
        {
          volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpy + 1;
        }
        // it IS transition state
        else
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyTL - oldSpecEnthalpy;

          // expansion in liquid state
          volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpyUpToTransition + 1;
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in transition state
          volume *= thermalExpansionT * deltaSpecEnthalpy + 1;
        }
      }
      // it WAS in transition state
      else
      {
        // IS it solid?
        if (newSpecEnthalpy <= specEnthalpyST)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyST - oldSpecEnthalpy;

          // expansion in transition state
          volume *= thermalExpansionT * deltaSpecEnthalpyUpToTransition + 1;
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in liquid state
          volume *= inv_CPS * thermalExpansionS * deltaSpecEnthalpy + 1;
        }
        // IS it liquid?
        else if (newSpecEnthalpy >= specEnthalpyTL)
        {
          const double deltaSpecEnthalpyUpToTransition = specEnthalpyTL - oldSpecEnthalpy;

          // expansion in transition state
          volume *= thermalExpansionT * deltaSpecEnthalpyUpToTransition + 1;
          deltaSpecEnthalpy -= deltaSpecEnthalpyUpToTransition;

          // expansion in liquid state
          volume *= inv_CPL * thermalExpansionL * deltaSpecEnthalpy + 1;
        }
        // it IS transition state
        else
        {
          volume *= thermalExpansionT * deltaSpecEnthalpy + 1;
        }
      }

      // --- compute the new volume --- //

      // updates
      (*radiusn)[lidNode] = PARTICLE::Utils::Volume2Radius(volume);
      (*densityn)[lidNode] = (*mass)[lidNode]/volume;
    }
  }
}

/*------------------------------------------------------------------------*
 | set up wall discretizations                               catta 06/16  |
 *------------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetUpWallDiscret()
{
  if(particlewalldis_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> walldisn = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> walldisnp = Teuchos::null;
    Teuchos::RCP<Epetra_Vector> wallvelnp = Teuchos::null;

    // solve for structural (wall) problem
    if(moving_walls_)
    {
      // extract displacement and velocity from full structural field to obtain wall states
      walldisn = wallextractor_->ExtractCondVector(structdispn_);
      walldisnp = wallextractor_->ExtractCondVector(structdispnp_);
      wallvelnp = wallextractor_->ExtractCondVector(structvelnp_);
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
  }
}
