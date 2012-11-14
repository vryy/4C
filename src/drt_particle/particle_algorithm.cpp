/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.cpp

\brief Algorithm to control particle simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
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
#include "../drt_meshfree_discret/drt_meshfree_bin.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io_gmsh.H"

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 09/12 |
 *----------------------------------------------------------------------*/
PARTICLE::Algorithm::Algorithm(
  const Epetra_Comm& comm
  ) : AlgorithmBase(comm,DRT::Problem::Instance()->StructuralDynamicParams()),
  particledis_(Teuchos::null),
  particles_(Teuchos::null),
  elecolmap_(Teuchos::null),
  myrank_(comm.MyPID())
{
  const Teuchos::ParameterList& meshfreeparams = DRT::Problem::Instance()->MeshfreeParams();
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

  particledis_ = DRT::Problem::Instance()->GetDis("particle");

  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  CreateBins();

  Teuchos::RCP<Epetra_Map> elerowmap = DistributeBinsToProcs();

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,Less> homelessparticles;

  for (int lid = 0; lid < noderowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(noderowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins and particles according to the bins
  SetupGhosting(elerowmap);

  // some output
  if(myrank_ == 0)
    cout << "after ghosting" << endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);

  // create time integrator based on structural time integration
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> particles =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->StructuralDynamicParams(), particledis_));
  particles_ = particles->StructureFieldrcp();

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
  particles_->Solve();
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
| read restart information for given time step              ghamm 09/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::Algorithm::ReadRestart(int restart)
{
  if (restart)
  {
    dserror("restart is not yet available");
  }

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
    bin_per_dir_[dim] = (int)std::floor( (XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_ );
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
  std::vector<int> roweles;

  // preliminary distribution for testing
  // TODO: Parmetis has to be used or distribution of bins dependent on background mesh
  int bin_per_proc = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2] / Comm().NumProc();
  // number of bins on last proc is adapted
  if(myrank_ != Comm().NumProc()-1)
  {
    for(int i=0; i<bin_per_proc; i++)
    {
      int gid = i + bin_per_proc*myrank_;
      Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEBIN","dummy", gid, myrank_);
      particledis_->AddElement(bin);
      roweles.push_back(gid);
    }
  }
  else
  {
    int last_bins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2] - (Comm().NumProc()-1)*bin_per_proc;
    for(int i=0; i<last_bins; i++)
    {
      int gid = i + bin_per_proc*myrank_;
      Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEBIN","dummy", gid, myrank_);
      particledis_->AddElement(bin);
      roweles.push_back(gid);
    }
  }

  // return elerowmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)roweles.size(),&roweles[0],0,Comm()));
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

  return;
}


/*----------------------------------------------------------------------*
| node is placed into the correct row bin                   ghamm 09/12 |
 *----------------------------------------------------------------------*/
bool PARTICLE::Algorithm::PlaceNodeCorrectly(Teuchos::RCP<DRT::Node> node, const double* currpos, std::set<Teuchos::RCP<DRT::Node>,Less>& homelessparticles)
{
//  cout << "node with ID: " << node->Id() << " and owner: " << node->Owner() << " arrived in PlaceNodeCorrectly" << endl;
  int ijk[3] = {0,0,0};
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)std::floor( (currpos[dim]-XAABB_(dim,0)) / bin_size_[dim] );
  }

  int binId = ConvertijkToGid(&ijk[0]);
  // check whether the current node belongs into a bin on this proc
  bool found = particledis_->HaveGlobalElement(binId);

  // either fill particle into correct bin on this proc or mark it as homeless
  if(found == true)
  {
    DRT::MESHFREE::MeshfreeBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeBin*>( particledis_->gElement(binId) );
#ifdef DEBUG
    if(currbin == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeBin failed");
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

          // received node is no longer needed, careful, X() does not have to be updated for the former ghost node
          // NOTE: here it is done only for output of correct discret with discret->Print(cout)
          {
            std::vector<double> update(3,0.0);
            for(int dim=0; dim < 3; dim++)
            {
              update[dim] = node->X()[dim] - existingnode->X()[dim];
            }
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
void PARTICLE::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> elerowmap)
{
  // gather elements of rowmap and all its neighbors (row + ghost)
  std::set<int> elements;
  for (int lid=0;lid<elerowmap->NumMyElements();++lid)
  {
    int elegid = elerowmap->GID(lid);
    int ijk[3] = {-1,-1,-1};
    ConvertGidToijk(elegid, ijk);

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
            elements.insert(neighborgid);
          }
        } // end for int k
      } // end for int j
    } // end for int i
  } // end for lid

  // copy elegids to a vector and create elecolmap
  std::vector<int> coleles(elements.begin(),elements.end());
  elecolmap_ = Teuchos::rcp(new Epetra_Map(-1,(int)coleles.size(),&coleles[0],0,Comm()));

  // create ghosting for bins (each knowing its particle ids)
  particledis_->ExportColumnElements(*elecolmap_);

  // create a set of node IDs for each proc (row + ghost)
  std::set<int> nodes;
  for (int lid=0;lid<elecolmap_->NumMyElements();++lid)
  {
    DRT::Element* ele = particledis_->gElement(elecolmap_->GID(lid));
    const int* nodeids = ele->NodeIds();
    for(int inode=0;inode<ele->NumNode();inode++)
      nodes.insert(nodeids[inode]);
  }

  // copy nodegids to a vector and create nodecolmap
  std::vector<int> colnodes(nodes.begin(),nodes.end());
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

  // create ghosting for particles
  particledis_->ExportColumnNodes(*nodecolmap);

  //new dofs are numbered from zero, minnodgid is ignored and it does not register in static_dofsets_
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset = Teuchos::rcp(new DRT::IndependentDofSet(true));
  particledis_->ReplaceDofSet(independentdofset);
  particledis_->FillComplete(true, false, true);


#ifdef DEBUG
  // check whether each proc has only particles that are within bins on this proc
  for(int k=0; k<particledis_->NumMyColElements(); k++)
  {
    int binid = particledis_->lColElement(k)->Id();
    DRT::Node** nodes = particledis_->lColElement(k)->Nodes();

    for(int inode=0; inode<particledis_->lColElement(k)->NumNode(); inode++)
    {
      int ijk[3] = {-1,-1,-1};
      for(int dim=0; dim < 3; dim++)
      {
        ijk[dim] = (int)std::floor((nodes[inode]->X()[dim]-XAABB_(dim,0)) / bin_size_[dim]);
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
void PARTICLE::Algorithm::TransferParticles()
{
  // set of homeless particles
  std::set<Teuchos::RCP<DRT::Node>,Less> homelessparticles;

  // current positions of particles
  Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();

  // check in each bin whether particles have moved out
  for(int ibin=0; ibin<particledis_->NumMyRowElements(); ibin++)
  {
    DRT::MESHFREE::MeshfreeBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeBin*>(particledis_->lRowElement(ibin));
#ifdef DEBUG
    if(currbin == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeBin failed");
#endif
    int binId = currbin->Id();
    DRT::Node** nodes = currbin->Nodes();
    std::vector<int> tobemoved(0);
    for(int inode=0; inode<currbin->NumNode(); inode++)
    {
      DRT::Node* currnode = nodes[inode];
      int ijk[3] = {-1,-1,-1};
      // get the first gid of a node and convert it into a LID
      int gid = particledis_->Dof(currnode, 0);
      int lid = particles_->ExtractDispnp()->Map().LID(gid);
      for(int dim=0; dim < 3; dim++)
      {
        ijk[dim] = (int)std::floor(((*disnp)[lid+dim]-XAABB_(dim,0)) / bin_size_[dim]);
      }

      int gidofbin = ConvertijkToGid(&ijk[0]);
      if(gidofbin != binId) // particle has left current bin
      {
        // gather all node Ids that will be removed and remove them afterwards
        // (looping over nodes and deleting at the same time is detrimental)
        tobemoved.push_back(currnode->Id());
        // find new bin for particle
        double currpos[3];
        for(int dim=0; dim < 3; dim++)
        {
          currpos[dim] = (*disnp)[lid+dim];
        }
        /*bool placed = */PlaceNodeCorrectly(Teuchos::rcp(currnode,false), currpos, homelessparticles);


      } // end if(gidofbin != binId)

      double currpos[3];
      for(int dim=0; dim < 3; dim++)
      {
        currpos[dim] = (*disnp)[lid+dim];
      }
      // update reference configuration of homeless particles for more efficient round robin loop
      // TO BE ADDED AGAIN; JUST FOR BETTER TESTING, SHOULD BE IN INNERST LOOP if(gidofbin != binId)
//        if(placed == false)
      {
        std::vector<double> update(3,0.0);
        for(int dim=0; dim < 3; dim++)
        {
          update[dim] = currpos[dim] - currnode->X()[dim];
        }
        // change X() of current node
        currnode->ChangePos(update);
//        std::cout << "particle (Id: " << currnode->Id() << " ) position is updated to" << currnode->X()[0] << "  "<< currnode->X()[1] << "  "<< currnode->X()[2] << "  " << std::endl;
      }


    } // end for inode

    // finally remove nodes from their old bin
    for(size_t iter=0; iter<tobemoved.size(); iter++)
    {
      currbin->DeleteNode(tobemoved[iter]);
    }

  } // end for ibin

  cout << "There are " << homelessparticles.size() << " homeless particles on proc" << myrank_ << endl;

  // homeless particles are sent to their new processors where they are inserted into their correct bin
  FillParticlesIntoBins(homelessparticles);

  // check whether all procs have a filled particledis_,
  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
  particledis_->CheckFilledGlobally();

  // new ghosting is necessary
  {
    particledis_->ExportColumnElements(*elecolmap_);

    // create a set of node IDs for each proc (row plus ghost)
    std::set<int> nodes;
    for (int lid=0;lid<elecolmap_->NumMyElements();++lid)
    {
      DRT::Element* ele = particledis_->gElement(elecolmap_->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0;inode<ele->NumNode();inode++)
      {
//        cout << "ele with ID: " << ele->Id() << " on proc: " << myrank_ << " contains node: " << nodeids[inode] << endl;
        nodes.insert(nodeids[inode]);
      }
    }

    // copy nodegids to a vector and create nodecolmap
    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

    // create ghosting for particles
    particledis_->ExportColumnNodes(*nodecolmap);
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

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
//  particles_->Output();

  const std::string filename = IO::GMSH::GetFileName("particle_data", Step(), true, Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  // velocity
  {
    gmshfilecontent << "View \" " << "velocity" << " \" {\n";
    LINALG::Matrix<3,1> vectorvalue(true);

    for(int n=0; n<particledis_->NumMyRowNodes(); n++)
    {
      DRT::Node* actnode = particledis_->lRowNode(n);
      // get the first gid of a node and convert it into a LID
      int gid = particledis_->Dof(actnode, 0);
      int lid = particles_->ExtractDispnp()->Map().LID(gid);
      Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();
      Teuchos::RCP<Epetra_Vector> velnp = particles_->ExtractVelnp();
      LINALG::Matrix<3,1> posXYZDomain(true);
      for(int dim=0; dim < 3; dim++)
      {
        posXYZDomain(dim) = (*disnp)[lid+dim];
        vectorvalue(dim) = (*velnp)[lid+dim];
      }

      // write data to Gmsh file
      IO::GMSH::VectorToStream(posXYZDomain, vectorvalue, gmshfilecontent);
    }

    gmshfilecontent << "};\n";
  }

  // density
  {
    gmshfilecontent << "View \" " << "density" << " \" {\n";

    for(int n=0; n<particledis_->NumMyRowNodes(); n++)
    {
      DRT::Node* actnode = particledis_->lRowNode(n);
      // get the first gid of a node and convert it into a LID
      int gid = particledis_->Dof(actnode, 0);
      int lid = particles_->ExtractDispnp()->Map().LID(gid);
      Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();
      LINALG::Matrix<3,1> posXYZDomain(true);
      for (int dim=0; dim<3; dim++)
      {
        posXYZDomain(dim) = (*disnp)[lid+dim];
      }

      double density = (*rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractDensitynp())[n];

      // write data to Gmsh file
      IO::GMSH::ScalarToStream(posXYZDomain, density, gmshfilecontent);
    }

    gmshfilecontent << "};\n";
  }

  // radius
  {
    gmshfilecontent << "View \" " << "radius" << " \" {\n";

    for(int n=0; n<particledis_->NumMyRowNodes(); n++)
    {
      DRT::Node* actnode = particledis_->lRowNode(n);
      // get the first gid of a node and convert it into a LID
      int gid = particledis_->Dof(actnode, 0);
      int lid = particles_->ExtractDispnp()->Map().LID(gid);
      Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();
      LINALG::Matrix<3,1> posXYZDomain(true);
      for (int dim=0; dim<3; dim++)
      {
        posXYZDomain(dim) = (*disnp)[lid+dim];
      }

      double radius = (*rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp())[n];

      // write data to Gmsh file
      IO::GMSH::ScalarToStream(posXYZDomain, radius, gmshfilecontent);
    }

    gmshfilecontent << "};\n";
  }

  gmshfilecontent.close();

  return;
}


/*----------------------------------------------------------------------*
 | Class for comparing Teuchos::RCP<DRT::Node> in std::set ghamm 10/12  |
 *----------------------------------------------------------------------*/
bool PARTICLE::Less::operator()(const Teuchos::RCP<const DRT::Node>& first, const Teuchos::RCP<const DRT::Node>& second) const
{
  return first->Id() < second->Id();
}
