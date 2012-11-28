/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
CAVITATION::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : PARTICLE::Algorithm(comm,params),
  fluiddis_(Teuchos::null)
{

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Timeloop()
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
 | setup of the system                                      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Init()
{
  particledis_ = DRT::Problem::Instance()->GetDis("particle");

  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),false));
  fluid_ = fluid->FluidFieldrcp();

  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  CreateBins();

  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcs();

  if(binrowmap->NumGlobalElements() > fluiddis_->NumGlobalElements() / 4.0)
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

  // ghost bins and particles according to the bins
  SetupGhosting(binrowmap);

  // some output
  IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*fluiddis_);

  // create time integrator based on structural time integration
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> particles =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(), particledis_));
  particles_ = particles->StructureFieldrcp();

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareTimeStep()
{
  PARTICLE::Algorithm::PrepareTimeStep();
  fluid_->PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate()
{
  fluid_->NonlinearSolve();
  PARTICLE::Algorithm::Integrate();
  return;
}


/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Update()
{
  fluid_->Update();
  PARTICLE::Algorithm::Update();
  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ReadRestart(int restart)
{
  fluid_->ReadRestart(restart);
  PARTICLE::Algorithm::ReadRestart(restart);
  return;
}


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CreateBins()
{
  // if not yet specified, get XAABB_ from underlying discretization
  if( XAABB_(2,1) > 0.9e12  and  XAABB_(2,1) < 1.1e12 )
  {
    IO::cout << "XAABB is computed based on the underlying fluid discretization" << IO::endl;
    XAABB_ = GEO::getXAABBofDis(*fluiddis_);
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
    bin_per_dir_[dim] = (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_);
    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

  IO::cout << "Global bounding box size: " << XAABB_;
  IO::cout << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << IO::endl;

  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CAVITATION::Algorithm::DistributeBinsToProcs()
{
  std::vector<int> rowbins;

  int numbins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];

  // loop over all fluid nodes and determine their number in each bin (init with 0)
  std::vector<int> mynumnodes_per_bin(numbins,0);
  for(int inode=0; inode<fluiddis_->NumMyRowNodes(); inode++)
  {
    int ijk[3] = {0,0,0};
    const double* currpos = fluiddis_->lRowNode(inode)->X();
    for(int dim=0; dim < 3; dim++)
    {
      ijk[dim] = (int)((currpos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
    }

    int binId = ConvertijkToGid(&ijk[0]);
    ++mynumnodes_per_bin[binId];
  }

  // find maximum number of nodes in each bin over all procs (init with -1)
  std::vector<int> maxnumnodes_per_bin(numbins,-1);
  fluiddis_->Comm().MaxAll(&mynumnodes_per_bin[0], &maxnumnodes_per_bin[0], numbins);

  // it is possible that several procs have the same number of nodes in a bin
  // only proc which has maximum number of nodes in a bin writes its rank
  std::vector<int> myrank_per_bin(numbins,-1);
  for(int i=0; i<numbins; i++)
  {
    if(mynumnodes_per_bin[i] == maxnumnodes_per_bin[i])
      myrank_per_bin[i] = myrank_;
  }

  mynumnodes_per_bin.clear();
  maxnumnodes_per_bin.clear();

  // find maximum myrank for each bin over all procs (init with -1)
  std::vector<int> maxmyrank_per_bin(numbins,-1);
  fluiddis_->Comm().MaxAll(&myrank_per_bin[0], &maxmyrank_per_bin[0], numbins);

  // distribute bins to proc with highest rank
  for(int gid=0; gid<numbins; gid++)
  {
    if(myrank_ == maxmyrank_per_bin[gid])
    {
      Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEBIN","dummy", gid, myrank_);
      particledis_->AddElement(bin);
      rowbins.push_back(gid);
    }
  }

  myrank_per_bin.clear();
  maxmyrank_per_bin.clear();

  // return binrowmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)rowbins.size(),&rowbins[0],0,Comm()));
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap)
{
  // 1st and 2nd step
  PARTICLE::Algorithm::SetupGhosting(binrowmap);

  // 3st step: extend ghosting of underlying fluid discretization according to bin distribution
  {
    // gather all fluid roweles /*and their nodes*/ in each bin
    std::map<int, std::set<int> > fluideles;
    /*std::map<int, std::set<int> > fluidnodes;*/
    for (int lid=0;lid<fluiddis_->NumMyRowElements();++lid)
    {
      DRT::Element* fluidele = fluiddis_->lRowElement(lid);
      DRT::Node** nodes = fluidele->Nodes();

      /*const int numnode = fluidele->NumNode();*/
      for(int inode=0; inode<fluidele->NumNode(); inode++)
      {
        DRT::Node* currnode = nodes[inode];
        const double* pos = currnode->X();
        int ijk[3] = {-1,-1,-1};
        for(int dim=0; dim < 3; dim++)
        {
          ijk[dim] = (int)(((pos)[dim]-XAABB_(dim,0)) / bin_size_[dim]);
        }
        int gidofbin = ConvertijkToGid(&ijk[0]);
        /*
        const int* nodeids = fluidele->NodeIds();
        // insert node ids into set
        fluidnodes[gidofbin].insert(nodeids,nodeids+numnode);
        */
        // insert ele id into set; one element can be part of several bins
        fluideles[gidofbin].insert(fluidele->Id());
      }
    }

    // do communication to gather all elements for extended ghosting
    const int numproc = fluiddis_->Comm().NumProc();
    std::map<int, std::set<int> > extendedfluidghosting;

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many col bins it has
      int numbin = bincolmap_->NumMyElements();
      fluiddis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which col bins it has
      std::vector<int> binid(numbin,0);
      if(iproc == myrank_)
      {
        int* bincolmap = bincolmap_->MyGlobalElements();
        for (int i=0; i<numbin; ++i)
          binid[i] = bincolmap[i];
      }
      fluiddis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int> > sdata;
      std::map<int, std::set<int> > rdata;

      for(int i=0; i<numbin; i++)
      {
        sdata[binid[i]].insert(fluideles[binid[i]].begin(),fluideles[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, fluiddis_->Comm());

      // proc i has to store the received data
      if(iproc == myrank_)
      {
        extendedfluidghosting = rdata;
      }
    }

    //reduce map of sets to one set and copy to a vector to create fluidcolmap
    std::set<int> redufluideleset;
    std::map<int, std::set<int> >::iterator iter;
    for(iter=extendedfluidghosting.begin(); iter!= extendedfluidghosting.end(); ++iter)
    {
      redufluideleset.insert(iter->second.begin(),iter->second.end());
    }
    std::vector<int> fluidcolgids(redufluideleset.begin(),redufluideleset.end());
    Teuchos::RCP<Epetra_Map> fluidcolmap = Teuchos::rcp(new Epetra_Map(-1,(int)fluidcolgids.size(),&fluidcolgids[0],0,Comm()));

    // create ghosting for fluid eles (each knowing its node ids)
    fluiddis_->ExportColumnElements(*fluidcolmap);

    // create a set of node IDs for each proc (row + ghost)
    std::set<int> nodes;
    for (int lid=0;lid<fluidcolmap->NumMyElements();++lid)
    {
      DRT::Element* ele = fluiddis_->gElement(fluidcolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0;inode<ele->NumNode();inode++)
        nodes.insert(nodeids[inode]);
    }

    // copy nodegids to a vector and create nodecolmap
    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

    // create ghosting for nodes
    fluiddis_->ExportColumnNodes(*nodecolmap);

    // do a final fillcomplete to build connectivity
    fluiddis_->FillComplete(true,true,true);

  }

#ifdef DEBUG
  // check whether each particle has an underlying fluid element
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for(int i=0; i<fluiddis_->NumMyColNodes(); i++)
  {
    DRT::Node* node = fluiddis_->lColNode(i);
    LINALG::Matrix<3,1> currpos;

    for (int a=0; a<3; a++)
    {
      currpos(a) = node->X()[a];
    }
    currentpositions.insert(std::pair<int,LINALG::Matrix<3,1> >(node->Id(),currpos));
  }
  // start loop over all particles
  for(int k=0; k<particledis_->NumMyColNodes(); k++)
  {
    DRT::Node* particle = particledis_->lColNode(k);
    const double* pos = particle->X();
    LINALG::Matrix<3,1> projpoint;
    for(int dim=0; dim<3; dim++)
      projpoint(dim) = pos[dim];
    bool foundele = false;
    for(int i=0; i<fluiddis_->NumMyColElements(); i++)
    {
      DRT::Element* fluidele = fluiddis_->lColElement(i);

      LINALG::Matrix<3,1> elecoord(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, elecoord);

      // THIS IS JUST TO CHECK WHETHER GEO::currentToVolumeElementCoordinates delivers correct results
      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::Matrix<3,numnode> xyze_linalg;
      for(int dim=0;dim<3;dim++)
        for(size_t n=0;n<numnode;n++)
          xyze_linalg(dim,n) = xyze(dim,n);
      GEO::CUT::Position<DRT::Element::hex8> pos(xyze_linalg, projpoint);
      bool withinele = pos.ComputeTol(GEO::TOL7);
      LINALG::Matrix<3,1> elecoordCut = pos.LocalCoordinates();

      if(abs(elecoordCut(0) - elecoord(0)) > GEO::TOL12 or abs(elecoordCut(1) - elecoord(1)) > GEO::TOL12
          or abs(elecoordCut(2) - elecoord(2)) > GEO::TOL12 or withinele != foundele)
        dserror("GEO::currentToVolumeElementCoordinates delivers different results compared to GEO::CUT::Position");
      // END: THIS IS JUST TO CHECK WHETHER GEO::currentToVolumeElementCoordinates delivers correct results

      if(foundele == true)
        break;
    }
    if(foundele == false)
      dserror("particle (Id:%d) was found which does not have fluid support", particle->Id());
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  PARTICLE::Algorithm::TestResults(comm);
  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Output()
{
  fluid_->Output();
  PARTICLE::Algorithm::Output();
  return;
}

