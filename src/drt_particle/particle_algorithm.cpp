/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.cpp

\brief Algorithm to control particle simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"
#include "../drt_adapter/ad_str_structure.H"
#include "particle_timint_centrdiff.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_inpar/inpar_meshfree.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_gmsh.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
PARTICLE::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : AlgorithmBase(comm,params),
  particledis_(Teuchos::null),
  particles_(Teuchos::null),
  bincolmap_(Teuchos::null),
  particlewalldis_(Teuchos::null),
  wall_output_(Teuchos::null),
  myrank_(comm.MyPID())
{
  const Teuchos::ParameterList& meshfreeparams = DRT::Problem::Instance()->MeshfreeParams();
  // safety check
  INPAR::MESHFREE::meshfreetype meshfreetype = DRT::INPUT::IntegralValue<INPAR::MESHFREE::meshfreetype>(meshfreeparams,"TYPE");
  if (meshfreetype!=INPAR::MESHFREE::particle)
    dserror("MESHFREE -> TYPE must be Particle in input file.");

  cutoff_radius_ = meshfreeparams.get<double>("CUTOFF_RADIUS");

  XAABB_.PutScalar(1.0e12);
  // get bounding box specified in the input file
  std::istringstream xaabbstream(Teuchos::getNumericStringParameter(meshfreeparams,"BOUNDINGBOX"));
  for(int col=0; col<2; col++)
  {
    for(int row=0; row<3; row++)
    {
      double value = 1.0e12;
      if(xaabbstream >> value)
        XAABB_(row,col) = value;
    }
  }

  for(int dim=0; dim < 3; dim++)
  {
    bin_size_[dim] = 0.0;
  }

  const Teuchos::ParameterList& cavitationparams = DRT::Problem::Instance()->CavitationParams();
  gravity_acc_.PutScalar(0.0);
  // get acceleration vector due to gravity for particles
  std::istringstream accstream(Teuchos::getNumericStringParameter(cavitationparams,"GRAVITY_ACCELERATION"));
  for(int dim=0; dim<3; dim++)
  {
    double value = 0.0;
    if(accstream >> value)
      gravity_acc_(dim) = value;
  }

  // initial setup of particle discretization
  particledis_ = DRT::Problem::Instance()->GetDis("particle");
  // new dofs are numbered from zero, minnodgid is ignored and it does not register in static_dofsets_
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset = Teuchos::rcp(new DRT::IndependentDofSet(true));
  particledis_->ReplaceDofSet(independentdofset);

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
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // transfer particles into their correct bins
    TransferParticles();

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
  CreateBins();

  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcs();

  if(binrowmap->NumGlobalElements() > particlerowmap->NumGlobalElements() / 4.0)
    IO::cout << "\n\n\n WARNING: Reduction of number of bins recommended!! Increase cutoff radius. \n\n\n" << IO::endl;

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,PARTICLE::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins and particles according to the bins --> final FillComplete() call included
  SetupGhosting(binrowmap);

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // add fully redundant discret for particle walls with identical dofs to full structural discret
    Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
    if (not structdis->Filled() or not structdis->HaveDofs())
      structdis->FillComplete();
    SetupParticleWalls(structdis);

    // access structural dynamic params list which will be possibly modified while creating the time integrator
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    // create time integrator based on structural time integration
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams(), const_cast<Teuchos::ParameterList&>(sdyn), particledis_));
    particles_ = particles->StructureFieldrcp();

    // determine consistent initial acceleration for the particles
    CalculateAndApplyForcesToParticles();
    Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->DetermineMassDampConsistAccel();
  }

  // some output
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
  // PrepareTimeStep() does nothing for explicit time integrators
//  particles_->PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Integrate()
{
  CalculateAndApplyForcesToParticles();

  {
    Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("PARTICLE::Algorithm::Integrate");
    Teuchos::TimeMonitor monitor(*t);
    particles_->Solve();
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate forces on particle and apply it               ghamm 02/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::CalculateAndApplyForcesToParticles()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("PARTICLE::Algorithm::CalculateAndApplyForcesToParticles");
  Teuchos::TimeMonitor monitor(*t);

  particledis_->ClearState();

  // particle radius of layout nodal col map --> SetState not possible
  Teuchos::RCP<Epetra_Vector> particleradius = LINALG::CreateVector(*particledis_->NodeColMap(),false);
  LINALG::Export(*Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp(),*particleradius);

  // vector to be filled with forces
  Teuchos::RCP<Epetra_Vector> particleforces = LINALG::CreateVector(*particledis_->DofColMap(),true);

  double particledensity = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ParticleDensity();

  // all particles (incl ghost) are evaluated
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ibin++)
  {
    DRT::Element* actbin = particledis_->lColElement(ibin);
    DRT::Node** particlesinbin = actbin->Nodes();

    // calculate forces for all particles in this bin
    for(int iparticle=0; iparticle<actbin->NumNode(); iparticle++)
    {
      const DRT::Node* currparticle = particlesinbin[iparticle];
      std::vector<int> lm_b = particledis_->Dof(currparticle);

      // get particle radius
      double r_p = (*particleradius)[ particledis_->NodeColMap()->LID(currparticle->Id()) ];

      // variable to sum forces for the current particle under observation
      Epetra_SerialDenseVector forcecurrparticle(3);

      /*------------------------------------------------------------------*/
      //// gravity forces = volume_p * rho_p * g
      double vol_bub = 4.0 / 3.0 * M_PI * r_p * r_p* r_p;
      double mass = vol_bub * particledensity;
      LINALG::Matrix<3,1> gravityforce(true);
      gravityforce.Update(mass, gravity_acc_);
      //assemble
      for(int dim=0; dim<3; dim++)
        forcecurrparticle[dim] += gravityforce(dim);
      /*------------------------------------------------------------------*/
      // assemble of particle forces can be done without further need of LINALG::Export (ghost particles also evaluated)
      std::vector<int> lmowner_b(lm_b.size(), myrank_);
      LINALG::Assemble(*particleforces,forcecurrparticle,lm_b,lmowner_b);
    } // end iparticle
  } // end ibin


  // apply forces to particles
  particledis_->SetState("particleforces", particleforces);

  return;
}


/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Update()
{
  // update of state vectors to the new maps
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
  // TODO: must be combined
  particles_->Update();
  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ReadRestart(int restart)
{
  // particledis_ needs to be cleared to remove initial particles (bins could stay!?!)
  particledis_->ClearDiscret();

  // read in particles for restart
  {
    IO::DiscretizationReader reader(particledis_, restart);
    reader.ReadNodesOnly(restart);
  }

  // Init() is needed to obtain connectivity -> includes FillComplete())
  Init(true);

  // now, correct map layouts are available and states can be read
  particles_->ReadRestart(restart);
  SetTimeStep(particles_->GetTime(),restart);

  // CAREFUL when deformable walls are handled in case of restart (history data)

  return;
}


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::CreateBins()
{
  // if not yet specified, get XAABB_ from underlying discretization
  if( XAABB_(2,1) > 0.9e12  and  XAABB_(2,1) < 1.1e12 )
  {
    IO::cout << "XAABB is computed based on the underlying particle discretization" << IO::endl;
    XAABB_ = GEO::getXAABBofNodes(*particledis_);
    // local bounding box
    double locmin[3] = {XAABB_(0,0), XAABB_(1,0), XAABB_(2,0)};
    double locmax[3] = {XAABB_(0,1), XAABB_(1,1), XAABB_(2,1)};
    // global bounding box
    double globmin[3];
    double globmax[3];
    // do the necessary communication
    Comm().MinAll(&locmin[0], &globmin[0], 3);
    Comm().MaxAll(&locmax[0], &globmax[0], 3);

    for(int dim=0; dim<3; dim++)
    {
      XAABB_(dim,0) = globmin[dim];
      XAABB_(dim,1) = globmax[dim];
    }
  }

  // divide global bounding box into bins
  for (int dim = 0; dim < 3; dim++)
  {
    // std::floor leads to bins that are at least of size cutoff_radius
    bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));
    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

  IO::cout << "Global bounding box size: " << XAABB_;
  IO::cout << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << IO::endl;

  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> PARTICLE::Algorithm::DistributeBinsToProcs()
{
  std::vector<int> rowbins;

  // preliminary distribution for testing
  // TODO: Parmetis has to be used or distribution of bins dependent on background mesh
  int bin_per_proc = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2] / Comm().NumProc();
  // number of bins on last proc is adapted
  if(myrank_ != Comm().NumProc()-1)
  {
    for(int i=0; i<bin_per_proc; i++)
    {
      int gid = i + bin_per_proc*myrank_;
      Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
      particledis_->AddElement(bin);
      rowbins.push_back(gid);
    }
  }
  else
  {
    int last_bins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2] - (Comm().NumProc()-1)*bin_per_proc;
    for(int i=0; i<last_bins; i++)
    {
      int gid = i + bin_per_proc*myrank_;
      Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
      particledis_->AddElement(bin);
      rowbins.push_back(gid);
    }
  }

  // return binrowmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)rowbins.size(),&rowbins[0],0,Comm()));
}


/*----------------------------------------------------------------------*
| fill particles into their correct bin on according proc   ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::FillParticlesIntoBins(std::set<Teuchos::RCP<DRT::Node>,Less>& homelessparticles)
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
    std::cout << " There are " << homelessparticles.size() << " particles which have left the computational domain on rank " << myrank << std::endl;
  // erase everything that is left
  homelessparticles.clear();

  return;
}


/*----------------------------------------------------------------------*
| node is placed into the correct row bin                   ghamm 09/12 |
 *----------------------------------------------------------------------*/
bool PARTICLE::Algorithm::PlaceNodeCorrectly(Teuchos::RCP<DRT::Node> node, const double* currpos, std::set<Teuchos::RCP<DRT::Node>,Less>& homelessparticles)
{
//  cout << "node with ID: " << node->Id() << " and owner: " << node->Owner() << " arrived in PlaceNodeCorrectly" << endl;
  int binId = ConvertPosToGid(currpos);

  // check whether the current node belongs into a bin on this proc
  bool found = particledis_->HaveGlobalElement(binId);

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
//      cout << "for node " << node->Id() << " a row bin was found on proc " << myrank_ << endl;
      // node already exists (either row or ghost)
      if( particledis_->HaveGlobalNode(node->Id()) == true)
      {
        DRT::Node* existingnode = particledis_->gNode(node->Id());
        // existing node is a row node, this means that node is equal existingnode
        if(existingnode->Owner() == myrank_)
        {
//          cout << "existingnode row node " << existingnode->Id() << " (ID from outside node: " << node->Id() << ") is added to element: " << currbin->Id() << " on proc " << myrank_ << endl;

          // assign node to the correct bin
          currbin->AddNode(existingnode);
        }
        else // ghost node becomes row node and node from outside is trashed
        {
//          cout << "existingnode ghost node " << existingnode->Id() << " (ID from outside node: " << node->Id() << ") is added to element: " << currbin->Id() << " on proc " << myrank_ << " after setting ownership" << endl;

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
//        cout << "node " << node->Id() << " is added to the discretization and assigned to element: " << currbin->Id() << " on proc " << myrank_ << endl;
        // assign node to the correct bin
        currbin->AddNode(node.get());
      }

      return true;
    }
    else // ghost bin
    {
      homelessparticles.insert(node);
      return false;
    }
  }
  else // bin not found on this proc
  {
    homelessparticles.insert(node);
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
      int gid = binrowmap->GID(lid);
      int ijk[3] = {-1,-1,-1};
      ConvertGidToijk(gid, ijk);

      // get all neighboring cells, including the element itself: one layer ghosting
      for(int i=-1;i<2;i++)
      {
        for(int j=-1;j<2;j++)
        {
          for(int k=-1;k<2;k++)
          {
            int ijk_neighbor[3] = {ijk[0]+i, ijk[1]+j, ijk[2]+k};

            int neighborgid = ConvertijkToGid(&ijk_neighbor[0]);
            if(neighborgid != -1)
            {
              bins.insert(neighborgid);
            }
          } // end for int k
        } // end for int j
      } // end for int i
    } // end for lid

    // copy bingids to a vector and create bincolmap
    std::vector<int> bincolmap(bins.begin(),bins.end());
    bincolmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)bincolmap.size(),&bincolmap[0],0,Comm()));

    // create ghosting for bins (each knowing its particle ids)
    particledis_->ExportColumnElements(*bincolmap_);
  }

  // 2st step: ghosting of particles according to bin distribution
  {
    // create a set of particle IDs for each proc (row + ghost)
    std::set<int> particles;
    for (int lid=0;lid<bincolmap_->NumMyElements();++lid)
    {
      DRT::Element* bin = particledis_->gElement(bincolmap_->GID(lid));
      const int* particleids = bin->NodeIds();
      for(int iparticle=0;iparticle<bin->NumNode();iparticle++)
        particles.insert(particleids[iparticle]);
    }

    // copy particlegids to a vector and create particlecolmap
    std::vector<int> colparticles(particles.begin(),particles.end());
    Teuchos::RCP<Epetra_Map> particlecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colparticles.size(),&colparticles[0],0,Comm()));

    // create ghosting for particles
    particledis_->ExportColumnNodes(*particlecolmap);

    // Note: IndependentDofSet is used; new dofs are numbered from zero, minnodgid is ignored and it does not register in static_dofsets_
    particledis_->FillComplete(true, false, true);


#ifdef DEBUG
    // check whether each proc has only particles that are within bins on this proc
    for(int k=0; k<particledis_->NumMyColElements(); k++)
    {
      int binid = particledis_->lColElement(k)->Id();
      DRT::Node** particles = particledis_->lColElement(k)->Nodes();

      for(int iparticle=0; iparticle<particledis_->lColElement(k)->NumNode(); iparticle++)
      {
        int ijk[3] = {-1,-1,-1};
        for(int dim=0; dim < 3; dim++)
        {
          ijk[dim] = (int)((particles[iparticle]->X()[dim]-XAABB_(dim,0)) / bin_size_[dim]);
        }

        int gidofbin = ConvertijkToGid(&ijk[0]);
        if(gidofbin != binid)
          dserror("after ghosting: particle which should be in bin no. %i is in %i",gidofbin,binid);
      }
    }
#endif
  }

  return;
}


/*----------------------------------------------------------------------*
 | particles are checked and transferred if necessary       ghamm 10/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::TransferParticles()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("PARTICLE::Algorithm::TransferParticles");
  Teuchos::TimeMonitor monitor(*t);

  // set of homeless particles
  std::set<Teuchos::RCP<DRT::Node>,Less> homelessparticles;

  // current positions of particles
  Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();

  // check in each bin whether particles have moved out
  for(int ibin=0; ibin<particledis_->NumMyRowElements(); ibin++)
  {
    DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->lRowElement(ibin));
#ifdef DEBUG
    if(currbin == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    int binId = currbin->Id();
    DRT::Node** particles = currbin->Nodes();
    std::vector<int> tobemoved(0);
    for(int iparticle=0; iparticle<currbin->NumNode(); iparticle++)
    {
      DRT::Node* currnode = particles[iparticle];
      // get the first gid of a node and convert it into a LID
      int gid = particledis_->Dof(currnode, 0);
      int lid = disnp->Map().LID(gid);

      double currpos[3];
      for(int dim=0; dim<3; dim++)
        currpos[dim] = (*disnp)[lid+dim];

      // update reference configuration of particle for correct output and correct placement via MPI
      {
        std::vector<double> update(3,0.0);
        const double* refposparticle = currnode->X();
        for(int dim=0; dim<3; dim++)
        {
          update[dim] = currpos[dim] - refposparticle[dim];
        }
        // change X() of current particle
        currnode->ChangePos(update);
//        std::cout << "particle (Id: " << currnode->Id() << " ) position is updated to" << currnode->X()[0] << "  "<< currnode->X()[1] << "  "<< currnode->X()[2] << "  " << std::endl;
      }

      int gidofbin = ConvertPosToGid(currpos);
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
      currbin->DeleteNode(tobemoved[iter]);
    }

  } // end for ibin

  if(homelessparticles.size())
    cout << "There are " << homelessparticles.size() << " homeless particles on proc" << myrank_ << endl;

  // homeless particles are sent to their new processors where they are inserted into their correct bin
  FillParticlesIntoBins(homelessparticles);

  // check whether all procs have a filled particledis_,
  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
  particledis_->CheckFilledGlobally();

  // new ghosting is necessary
  {
    particledis_->ExportColumnElements(*bincolmap_);

    // create a set of particle IDs for each proc (row plus ghost)
    std::set<int> particles;
    for (int lid=0;lid<bincolmap_->NumMyElements();++lid)
    {
      DRT::Element* ele = particledis_->gElement(bincolmap_->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0;inode<ele->NumNode();inode++)
      {
//        cout << "ele with ID: " << ele->Id() << " on proc: " << myrank_ << " contains node: " << nodeids[inode] << endl;
        particles.insert(nodeids[inode]);
      }
    }

    // copy nodegids to a vector and create nodecolmap
    std::vector<int> colparticles(particles.begin(),particles.end());
    Teuchos::RCP<Epetra_Map> particlecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colparticles.size(),&colparticles[0],0,Comm()));

    // create ghosting for particles
    particledis_->ExportColumnNodes(*particlecolmap);
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

  // reconstruct element -> bin pointers for particle wall elements
  BuildWallElementToBinPointers();

  return;
}


/*----------------------------------------------------------------------*
| particle walls are added from the structural discret      ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::SetupParticleWalls(Teuchos::RCP<DRT::Discretization> structdis)
{
  //--------------------------------------------------------------------
  // 1st step: build fully redundant discretization with wall elements
  //--------------------------------------------------------------------

  // declare struct objects in wall condition
  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > > structgelements; // col map of structure elements
  std::map<int, std::map<int, DRT::Node*> > dummy2;  // dummy map
  std::map<int, std::map<int, DRT::Node*> > structgnodes; // col map of structure nodes

  //initialize struct objects in wall condition
  DRT::UTILS::FindConditionObjects(*structdis, dummy2, structgnodes, structgelements,"ParticleWall");

  std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit;

  // initialize new particle wall discretizations
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(structdis->Comm().Clone());
  const std::string discret_name = "particlewalls";
  Teuchos::RCP<DRT::Discretization> particlewalldis = Teuchos::rcp(new DRT::Discretization(discret_name,com));

  std::vector<int> nodeids;
  std::vector<int> eleids;
  // loop over all particle wall nodes and elements and fill new discretization
  for(std::map<int, std::map<int, Teuchos::RCP<DRT::Element> > >::iterator meit=structgelements.begin(); meit!=structgelements.end(); ++meit)
  {
    // care about particle wall nodes
    std::map<int, DRT::Node*> wallgnodes = structgnodes[meit->first];
    for (std::map<int, DRT::Node* >::iterator nit=wallgnodes.begin(); nit != wallgnodes.end(); ++nit)
    {
      DRT::Node* currnode = (*nit).second;
      if (currnode->Owner() == myrank_)
      {
        nodeids.push_back(currnode->Id());
        particlewalldis->AddNode(Teuchos::rcp(new DRT::Node(currnode->Id(), currnode->X(), currnode->Owner())));
      }
    }

    // care about particle wall eles
    std::map<int, Teuchos::RCP<DRT::Element> > structelementsinterf = structgelements[meit->first];
    for (std::map<int, Teuchos::RCP<DRT::Element> >::iterator eit=structelementsinterf.begin(); eit != structelementsinterf.end(); ++eit)
    {
      Teuchos::RCP<DRT::Element> currele = eit->second;
      if (currele->Owner() == myrank_)
      {
        eleids.push_back(currele->Id() );
        // structural surface elements cannot be distributed --> Bele3 element is used
        Teuchos::RCP<DRT::Element> wallele = DRT::UTILS::Factory("BELE3","Polynomial", currele->Id(), currele->Owner());
        wallele->SetNodeIds(currele->NumNode(), currele->NodeIds());
        particlewalldis->AddElement( wallele );
      }
    }
  }

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

  // find out if we are in parallel; needed for TransparentDofSet
  bool parallel = (particlewalldis->Comm().NumProc() == 1) ? false : true;

  // dofs of the original discretization are used to set same dofs for the new particle wall discretization
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentDofSet(structdis,parallel));
  particlewalldis->ReplaceDofSet(newdofset);
  newdofset=Teuchos::null;

  // final fill complete to reorganize everything in the discretization
  particlewalldis->FillComplete(true, false, false);
  particlewalldis_ = particlewalldis;

  // some output to screen and initialization of binary output
  IO::cout << "after adding particle walls" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particlewalldis_);
  wall_output_ = Teuchos::rcp(new IO::DiscretizationWriter(particlewalldis_));


  //--------------------------------------------------------------------
  // 2nd step: assign wall elements to bins
  //--------------------------------------------------------------------

  // so far we have only fixed walls
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for (int lid = 0; lid < particlewalldis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = particlewalldis_->lColNode(lid);
    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];
    currentpositions[node->Id()] = currpos;
  }

  // find bins for all wall elements
  for (int lid = 0; lid < particlewalldis_->NumMyColElements(); ++lid)
  {
    DRT::Element* wallele = particlewalldis_->lColElement(lid);
    DRT::Node** wallnodes = wallele->Nodes();
    const int numnode = wallele->NumNode();
    // variable to store bin ids in which this wall element is located
    std::set<int> binIds;

    //// 2.1) do a positive search and get all bins enclosed in the bounding box of each wall element
    {
      // initialize ijk_range with ijk of first node of wall element
      int ijk[3];
      {
        const DRT::Node* node = wallnodes[0];
        const double* coords = node->X();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; j++)
      {
        const DRT::Node* node = wallnodes[j];
        const double* coords = node->X();
        int ijk[3];
        ConvertPosToijk(coords, ijk);

        for(int dim=0; dim<3; dim++)
        {
          if(ijk[dim]<ijk_range[dim*2])
            ijk_range[dim*2]=ijk[dim];
          if(ijk[dim]>ijk_range[dim*2+1])
            ijk_range[dim*2+1]=ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range and fill them into binIds
      GidsInijkRange(&ijk_range[0], binIds);
    }

    //// 2.2) do a first negative search and remove bins that are not on this processor
    {
      std::set<int> binnotfound;
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
      {
        if(not particledis_->HaveGlobalElement(*biniter))
          binnotfound.insert(*biniter);
      }
      for(std::set<int>::const_iterator biniter=binnotfound.begin(); biniter!=binnotfound.end(); ++biniter)
        binIds.erase(*biniter);
    }

    // if already all bins have been removed, next wall element can be processed
    if(binIds.empty())
      break;

    //// 2.3) do a second negative search and remove bins that are too far away from the wall element
    // all corners of the close bin are projected onto the wall element: if at least one projection
    // point is inside the bin, it won't be removed from the list
    {
      Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("PARTICLE::Algorithm::SetupParticleWalls::step2.3_check");
      Teuchos::TimeMonitor monitor(*t);

      std::map<int, Teuchos::RCP<DRT::Element> > elements;
      elements[wallele->Id()] = Teuchos::rcp(wallele, false);
      std::map<int, std::set<int> > elementList;
      elementList[-1].insert(wallele->Id());
      std::set<int> binfaraway;

      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
      {
        std::vector<LINALG::Matrix<3,1> > bincorners;
        GetBinCorners(*biniter, bincorners);

        // in case wall element is axis aligned, it might not be detected as inside because projection points
        // are located on the edges of the bin --> Remedy: bin centroid is tested as well
        bincorners.push_back(GetBinCentroid(*biniter));

        bool projpointinsidebin = false;
        // loop over all corners + centroid of one bin and project them onto the wall element
        for(size_t corner=0; corner<bincorners.size(); corner++)
        {
          //search for the closest object, more exactly it's coordinates
          LINALG::Matrix<3,1> minDistCoords;
          int eleid = GEO::nearest3DObjectInNode(particlewalldis_, elements, currentpositions,
            elementList, bincorners[corner], minDistCoords);
          if(eleid == -1)
            dserror("Closest surface not found. This should be the neighborhood of the element already!?!");

          int gid = ConvertPosToGid(minDistCoords);
          if(gid == *biniter)
          {
            projpointinsidebin = true;
            break;
          }
        }
        if(projpointinsidebin == false)
          binfaraway.insert(*biniter);
      }

      for(std::set<int>::const_iterator biniter=binfaraway.begin(); biniter!=binfaraway.end(); ++biniter)
        binIds.erase(*biniter);
    }

    //// 2.3) assign wall element to remaining bins
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
void PARTICLE::Algorithm::BuildWallElementToBinPointers()
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
    for(int iwall=0; iwall<numwallele; iwall++)
    {
      const int wallid = walleleids[iwall];
      wallelements[iwall] = particlewalldis_->gElement(wallid);
    }
    actbin->BuildWallElePointers(&wallelements[0]);
  }
  return;
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 01/13 |
 *----------------------------------------------------------------------*/
int PARTICLE::Algorithm::ConvertPosToGid(const std::vector<double>& pos)
{
  int ijk[3] = {0,0,0};
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
  }

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 01/13 |
 *----------------------------------------------------------------------*/
int PARTICLE::Algorithm::ConvertPosToGid(const double* pos)
{
  int ijk[3] = {0,0,0};
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
  }

  return ConvertijkToGid(&ijk[0]);
}

/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 02/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ConvertPosToijk(const double* pos, int* ijk)
{
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
  }
  return;
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
int PARTICLE::Algorithm::ConvertPosToGid(const LINALG::Matrix<3,1> pos)
{
  int ijk[3] = {0,0,0};
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)((pos(dim)-XAABB_(dim,0)) / bin_size_[dim]);
  }

  return ConvertijkToGid(&ijk[0]);
}

/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ConvertPosToijk(const LINALG::Matrix<3,1> pos, int* ijk)
{
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)((pos(dim)-XAABB_(dim,0)) / bin_size_[dim]);
  }
  return;
}


/*----------------------------------------------------------------------*
| convert i,j,k into bin id                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
int PARTICLE::Algorithm::ConvertijkToGid(int* ijk)
{
  // given ijk is outside of XAABB
  if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 || ijk[0]>=bin_per_dir_[0] || ijk[1]>=bin_per_dir_[1] || ijk[2]>=bin_per_dir_[2] )
    return -1;

  return ijk[0] + ijk[1]*bin_per_dir_[0] + ijk[2]*bin_per_dir_[0]*bin_per_dir_[1];
}


/*----------------------------------------------------------------------*
| convert bin id into i,j,k                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ConvertGidToijk(int gid, int* ijk)
{
  ijk[2] = gid / (bin_per_dir_[0]*bin_per_dir_[1]);

  int tmp = gid - ijk[2]*bin_per_dir_[0]*bin_per_dir_[1];

  ijk[1] = tmp / bin_per_dir_[0];

  ijk[0] = tmp - ijk[1]*bin_per_dir_[0];

  // found ijk is outside of XAABB
  if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 || ijk[0]>=bin_per_dir_[0] || ijk[1]>=bin_per_dir_[1] || ijk[2]>=bin_per_dir_[2] )
    ijk[0] = -1;

  return;
}


/*----------------------------------------------------------------------*
 | get all bins in ijk range                               ghamm 02/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GidsInijkRange(int* ijk_range, std::set<int>& binIds)
{
  for(int i=ijk_range[0]; i<=ijk_range[1]; i++)
  {
    for(int j=ijk_range[2]; j<=ijk_range[3]; j++)
    {
      for(int k=ijk_range[4]; k<=ijk_range[5]; k++)
      {
        int ijk[3] = {i,j,k};

        int gid = ConvertijkToGid(&ijk[0]);
        if(gid != -1)
        {
          binIds.insert(gid);
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return;
}


/*----------------------------------------------------------------------*
| corner position for given bin id                          ghamm 03/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::GetBinCorners(int binId, std::vector<LINALG::Matrix<3,1> >& bincorners)
{
  bincorners.clear();
  bincorners.reserve(8);
  int ijk_base[3];
  ConvertGidToijk(binId, &ijk_base[0]);

  // order in bincorners is identical to ordering of i,j and k
  for(int k=ijk_base[2]; k<(ijk_base[2]+2); k++)
  {
    for(int j=ijk_base[1]; j<(ijk_base[1]+2); j++)
    {
      for(int i=ijk_base[0]; i<(ijk_base[0]+2); i++)
      {
        int ijk_curr[] = {i,j,k};
        LINALG::Matrix<3,1> curr_corner;
        for(int dim=0; dim<3; dim++)
        {
          curr_corner(dim) = XAABB_(dim,0) + bin_size_[dim]*ijk_curr[dim];
        }
        bincorners.push_back(curr_corner);

      } // end for int k
    } // end for int j
  } // end for int i

  return;
}


/*----------------------------------------------------------------------*
| centroid position for given bin id                        ghamm 04/13 |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> PARTICLE::Algorithm::GetBinCentroid(int binId)
{
  int ijk[3];
  ConvertGidToijk(binId, ijk);
  if(ijk[0] == -1)
    dserror("given bin id is outside of bins; centroid of bin is does not make sense");

  LINALG::Matrix<3,1> centroid(false);
  for(int dim=0; dim<3; dim++)
  {
    centroid(dim) = XAABB_(dim,0) + bin_size_[dim]*ijk[dim] + bin_size_[dim]/2.0;
  }

  return centroid;
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
 | output particle time step                                ghamm 10/12  |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::Output()
{
  // INFO regarding output: Bins are not written to file because they cannot
  // be post-processed anyway (no nodes and connectivity available)
  particles_->Output();

  // so far: fixed walls
  if(Step() == 1 and wall_output_ != Teuchos::null)
  {
    wall_output_->WriteMesh(Step(),Time());
    wall_output_->NewStep(Step(),Time());
    Teuchos::RCP<Epetra_Vector> walldisp = Teuchos::rcp(new Epetra_Vector(*particlewalldis_->DofRowMap(), true));
    wall_output_->WriteVector("displacement", walldisp);
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
//      int lid = particles_->ExtractDispnp()->Map().LID(gid);
//      Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();
//      Teuchos::RCP<Epetra_Vector> velnp = particles_->ExtractVelnp();
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
//      int lid = particles_->ExtractDispnp()->Map().LID(gid);
//      Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();
//      LINALG::Matrix<3,1> posXYZDomain(true);
//      for (int dim=0; dim<3; dim++)
//      {
//        posXYZDomain(dim) = (*disnp)[lid+dim];
//      }
//
//      double density = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ParticleDensity();
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
//      int lid = particles_->ExtractDispnp()->Map().LID(gid);
//      Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();
//      LINALG::Matrix<3,1> posXYZDomain(true);
//      for (int dim=0; dim<3; dim++)
//      {
//        posXYZDomain(dim) = (*disnp)[lid+dim];
//      }
//
//      double radius = (*Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp())[n];
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
 | Class for comparing Teuchos::RCP<DRT::Node> in std::set ghamm 10/12  |
 *----------------------------------------------------------------------*/
bool PARTICLE::Less::operator()(const Teuchos::RCP<const DRT::Node>& first, const Teuchos::RCP<const DRT::Node>& second) const
{
  return first->Id() < second->Id();
}
