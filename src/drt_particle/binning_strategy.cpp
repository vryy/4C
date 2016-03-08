/*----------------------------------------------------------------------*/
/*!
\file binning_strategy.cpp

\brief Binning strategy for neighborhood search

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/13 |
 *----------------------------------------------------------------------*/
//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_node.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_discret.H"

#include "binning_strategy.H"

/*----------------------------------------------------------------------*
 | Binning strategy constructor                             ghamm 11/13 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy(
  const Epetra_Comm& comm,
  double cutoff_radius,
  LINALG::Matrix<3,2> XAABB
  ) :
  particledis_(Teuchos::null),
  cutoff_radius_(cutoff_radius),
  XAABB_(XAABB),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")),
  sparse_binning_(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CavitationParams(),"SPARSE_BIN_DISTRIBUTION")),
  myrank_(comm.MyPID())
{
  if( XAABB_(0,0) >= XAABB_(0,1) or XAABB_(1,0) >= XAABB_(1,1) or XAABB_(2,0) >= XAABB_(2,1))
    dserror("XAABB is not computed correctly");

  if(cutoff_radius_ <= 0.0)
    dserror("cutoff radius cannot be zero or negativ");

  // compute bins
  CreateBins(Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | Binning strategy constructor                             ghamm 11/13 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy(
  const Epetra_Comm& comm
  ) :
  particledis_(Teuchos::null),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")),
  sparse_binning_(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CavitationParams(),"SPARSE_BIN_DISTRIBUTION")),
  myrank_(comm.MyPID())
{
  const Teuchos::ParameterList& meshfreeparams = DRT::Problem::Instance()->MeshfreeParams();

  // get cutoff radius
  cutoff_radius_ = meshfreeparams.get<double>("CUTOFF_RADIUS");
  // get number of bins per direction
  std::istringstream binstream(Teuchos::getNumericStringParameter(meshfreeparams,"BIN_PER_DIR"));
  for(int idim=0; idim<3; idim++)
  {
    int val = -1;
    if (binstream >> val)
      bin_per_dir_[idim] = val;
  }
  // check input: either the cutoff_radius_ or the number of bins per direction have to be set
  if (cutoff_radius_<0.0 and bin_per_dir_[0]<0.0 and bin_per_dir_[1]<0.0 and bin_per_dir_[2]<0.0)
    dserror("Cutoff radius and number of bins per direction have not been set in the input file. Please prescribe the cutoff radius or define the number of bins for each spatial direction.");
  if (cutoff_radius_>0.0 and bin_per_dir_[0]>0.0 and bin_per_dir_[1]>0.0 and bin_per_dir_[2]>0.0)
    dserror("Cutoff radius and number of bins per direction have been set in the input file. Please prescribe only one of the two options");

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

  // initialize bin size
  for(int dim=0; dim<3; ++dim)
  {
    bin_size_[dim] = 0.0;
  }

  return;
}


/*----------------------------------------------------------------------*
 | Repartitioning Binning strategy constructor              ghamm 06/14 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy(
  std::vector<Teuchos::RCP<DRT::Discretization> > dis,
  std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
  std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap
  ) :
  particledis_(Teuchos::null),
  cutoff_radius_(0.0),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")),
  sparse_binning_(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CavitationParams(),"SPARSE_BIN_DISTRIBUTION")),
  myrank_(dis[0]->Comm().MyPID())
{
  WeightedRepartitioning(dis,stdelecolmap,stdnodecolmap);

  return;
}


/*----------------------------------------------------------------------*
| assign elements into bins                                 ghamm 11/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeElesToBins(
  const DRT::Discretization& mortardis,
  std::map<int, std::set<int> >& binelemap,
  bool isslave)
{
  // exploit bounding box idea for elements and bins
  for (int lid = 0; lid<mortardis.NumMyColElements(); ++lid)
  {
    DRT::Element* ele = mortardis.lColElement(lid);
    if(dynamic_cast<MORTAR::MortarElement*>(ele)->IsSlave() == isslave)
    {
      DRT::Node** nodes = ele->Nodes();
      const int numnode = ele->NumNode();

      // initialize ijk_range with ijk of first node of element
      int ijk[3];
      {
        DRT::Node* node = nodes[0];
        const double* coords = dynamic_cast<MORTAR::MortarNode*>(node)->xspatial();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; ++j)
      {
        DRT::Node* node = nodes[j];
        const double* coords = dynamic_cast<MORTAR::MortarNode*>(node)->xspatial();
        int ijk[3];
        ConvertPosToijk(coords, ijk);

        for(int dim=0; dim<3; ++dim)
        {
          if(ijk[dim]<ijk_range[dim*2])
            ijk_range[dim*2]=ijk[dim];
          if(ijk[dim]>ijk_range[dim*2+1])
            ijk_range[dim*2+1]=ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range
      std::vector<int> binIds;
      binIds.reserve((ijk_range[1]-ijk_range[0]+1) * (ijk_range[3]-ijk_range[2]+1) * (ijk_range[5]-ijk_range[4]+1));
      GidsInijkRange(&ijk_range[0], binIds, false);

      // assign element to bins
      for(std::vector<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        binelemap[*biniter].insert(ele->Id());
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| assign elements into bins                                 ghamm 11/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeElesToBins(
  Teuchos::RCP<DRT::Discretization> underlyingdis,
  std::map<int, std::set<int> >& rowelesinbin
  )
{
  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for (int lid = 0; lid < underlyingdis->NumMyRowElements(); ++lid)
  {
    DRT::Element* ele = underlyingdis->lRowElement(lid);
    DRT::Node** nodes = ele->Nodes();
    const int numnode = ele->NumNode();

    // initialize ijk_range with ijk of first node of fluid element
    int ijk[3];
    {
      const DRT::Node* node = nodes[0];
      const double* coords = node->X();
      ConvertPosToijk(coords, ijk);
    }

    // ijk_range contains: i_min i_max j_min j_max k_min k_max
    int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

    // fill in remaining nodes
    for (int j=1; j<numnode; ++j)
    {
      const DRT::Node* node = nodes[j];
      const double* coords = node->X();
      int ijk[3];
      ConvertPosToijk(coords, ijk);

      for(int dim=0; dim<3; ++dim)
      {
        if(ijk[dim]<ijk_range[dim*2])
          ijk_range[dim*2]=ijk[dim];
        if(ijk[dim]>ijk_range[dim*2+1])
          ijk_range[dim*2+1]=ijk[dim];
      }
    }

    // get corresponding bin ids in ijk range
    std::vector<int> binIds;
    binIds.reserve((ijk_range[1]-ijk_range[0]+1) * (ijk_range[3]-ijk_range[2]+1) * (ijk_range[5]-ijk_range[4]+1));
    GidsInijkRange(&ijk_range[0], binIds, false);

    // assign element to bins
    for(std::vector<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
      rowelesinbin[*biniter].insert(ele->Id());
  }

  return;
}

/*----------------------------------------------------------------------*
| assign nodes into bins                                    ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeNodesToBins(
  Teuchos::RCP<DRT::Discretization> underlyingdis,
  std::map<int, std::vector<int> >& nodesinbin
  )
{
  for (int lid = 0; lid < underlyingdis->NumMyRowNodes(); ++lid)
  {
    DRT::Node* node = underlyingdis->lRowNode(lid);

    const double* coords = node->X();
    int ijk[3];
    ConvertPosToijk(coords, ijk);
    int binid = ConvertijkToGid(&ijk[0]);

    // assign node to bin
    nodesinbin[binid].push_back(node->Id());
  }

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::WeightedRepartitioning(
  std::vector<Teuchos::RCP<DRT::Discretization> > dis,
  std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
  std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap
  )
{
  // Careful: At the moment only reference configuration is considered
  {
    // initialize XAABB as rectangle around the first node of dis
    const DRT::Node* node = dis[0]->lRowNode(0);
    for(int dim=0; dim<3; ++dim)
    {
      XAABB_(dim, 0) = node->X()[dim] - GEO::TOL7;
      XAABB_(dim, 1) = node->X()[dim] + GEO::TOL7;
    }
  }

  // build bounding box for all discrets and determine maximal element extension
  for(size_t i=0; i<dis.size(); ++i)
  {
    LINALG::Matrix<3,2> XAABB;
    double cutoff;
    CreateXAABB(dis[i], XAABB, cutoff);

    cutoff_radius_ = std::max(cutoff, cutoff_radius_);

    for(int dim=0; dim < 3; dim++)
    {
      XAABB_(dim, 0) = std::min( XAABB_(dim, 0), XAABB(dim, 0) );
      XAABB_(dim, 1) = std::max( XAABB_(dim, 1), XAABB(dim, 1) );
    }
  }

  // enlarge cutoff a little bit for safety reasons
  cutoff_radius_ += GEO::TOL7;

  // build bins on existing xaabb and cutoff radius
  CreateBins(Teuchos::null);

  // distribution of bins equally on all procs
  int numbin = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];
  if(numbin<dis[0]->Comm().NumProc() && myrank_ == 0)
    dserror("ERROR: Too many processors to distribute your bins properly!!!");
  if(numbin < 8*dis[0]->Comm().NumProc() && myrank_==0)
    std::cout << "\n\nWARNING: partitioning not useful, choose less procs. Owner distribution may be inefficient!\n\n" << std::endl;

  // dummy row bin distribution covering the whole XAABB domain
  Teuchos::RCP<Epetra_Map> rowbins = Teuchos::rcp(new Epetra_Map(numbin,0,dis[0]->Comm()));

  // create nodal graph of existing problem
   Teuchos::RCP<Epetra_CrsGraph> bingraph = Teuchos::rcp( new Epetra_CrsGraph(Copy,*rowbins,108,false));
   // Now we're going to create a Epetra_Vector with vertex weights to be used in the partitioning operation.
   // weights must be at least one for zoltan
   Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*rowbins, true);

   // set weights of bins related to the number of nodes that are contained
   // empty bins have weight of 1
   vweights->PutScalar(1.0);

   // determine which node is in which bin
   std::vector<std::map<int, std::vector<int> > > nodesinbin(dis.size());
  for(size_t i=0; i<dis.size(); ++i)
  {
    DistributeNodesToBins(dis[i], nodesinbin[i]);

    std::map<int, std::vector<int> > mynodesinbin;
    // gather information of bin content from other procs (bin is owned by this proc and there are
    // some nodes on other procs which are located in this bin here)
    CollectInformation(rowbins, nodesinbin[i], mynodesinbin);

    for(std::map<int, std::vector<int> >::const_iterator biniter=mynodesinbin.begin(); biniter!=mynodesinbin.end(); ++biniter)
    {
      int lid = rowbins->LID(biniter->first);
      (*vweights)[lid] += 10.0*(double)biniter->second.size();
    }
  }

   // fill bin connectivity into bin graph
   for (int lid=0; lid<rowbins->NumMyElements(); ++lid)
   {
     int rowbinid = rowbins->GID(lid);
     // insert neighbors
     std::vector<int> neighbors;
     GetBinConnectivity(rowbinid,neighbors);

     int err = bingraph->InsertGlobalIndices(rowbinid,(int)neighbors.size(),&neighbors[0]);
     if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,rowbinid);
   }

   // complete graph
   int err = bingraph->FillComplete();
   if (err) dserror("graph->FillComplete() returned err=%d",err);
   err = bingraph->OptimizeStorage();
   if (err) dserror("graph->OptimizeStorage() returned err=%d",err);

   // call redistribution of bin graph using bin weights
  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs = Teuchos::rcp(new Isorropia::Epetra::CostDescriber);
  costs->setVertexWeights(vweights);

  Teuchos::ParameterList paramlist;
  paramlist.set("PARTITIONING METHOD", "GRAPH");
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  sublist.set("LB_APPROACH", "PARTITION");
  sublist.set("GRAPH_PACKAGE", "PHG");
  sublist.set("OBJ_WEIGHT_DIM","1");
  sublist.set("EDGE_WEIGHT_DIM","0");

  Teuchos::RCP<const Epetra_CrsGraph> constbingraph(bingraph);

  // Now create the partitioner object
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(constbingraph, costs, paramlist));

  Isorropia::Epetra::Redistributor rd(partitioner);

  // redistribute bingraph
  Teuchos::RCP<Epetra_CrsGraph> balanced_bingraph = rd.redistribute(*bingraph);

  // extract repartitioned bin row map
  const Epetra_BlockMap& rbinstmp = balanced_bingraph->RowMap();
  Teuchos::RCP<Epetra_Map> newrowbins = Teuchos::rcp(new Epetra_Map(-1,rbinstmp.NumMyElements(),rbinstmp.MyGlobalElements(),0,dis[0]->Comm()));

  stdelecolmap.resize(dis.size());
  stdnodecolmap.resize(dis.size());

  // rebuild discretizations including extended ghosting
  for(size_t i=0; i<dis.size(); ++i)
  {
    //----------------------------
    // start with standard ghosting
    //----------------------------
    Teuchos::RCP<Epetra_CrsGraph> initgraph = dis[i]->BuildNodeGraph();

    std::map<int, std::vector<int> > mynodesinbin;
    // gather information of bin content from other procs (bin is owned by this proc and there are
    // some nodes on other procs which are located in this bin here)
    CollectInformation(newrowbins, nodesinbin[i], mynodesinbin);

    std::vector<int> mynewrownodes;
    for(std::map<int, std::vector<int> >::const_iterator biniter=mynodesinbin.begin(); biniter!=mynodesinbin.end(); ++biniter)
    {
      for(std::vector<int>::const_iterator nodeiter=biniter->second.begin(); nodeiter!=biniter->second.end(); ++nodeiter)
      {
        mynewrownodes.push_back(*nodeiter);
      }
    }
    mynodesinbin.clear();

    Teuchos::RCP<Epetra_Map> newnoderowmap = Teuchos::rcp(new Epetra_Map(-1,mynewrownodes.size(),&mynewrownodes[0],0,dis[0]->Comm()));

    // create the new graph and export to it
    Teuchos::RCP<Epetra_CrsGraph> newnodegraph;

    newnodegraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*newnoderowmap,108,false));
    Epetra_Export exporter2(initgraph->RowMap(),*newnoderowmap);
    err = newnodegraph->Export(*initgraph,exporter2,Add);
    if (err<0)
      dserror("Graph export returned err=%d",err);
    newnodegraph->FillComplete();
    newnodegraph->OptimizeStorage();

    // the column map will become the new ghosted distribution of nodes (standard ghosting)
    const Epetra_BlockMap cntmp = newnodegraph->ColMap();
    stdnodecolmap[i] = Teuchos::rcp(new Epetra_Map(-1,cntmp.NumMyElements(),cntmp.MyGlobalElements(),0,dis[0]->Comm()));

    // rebuild of the discretizations with new maps for standard ghosting
    Teuchos::RCP<Epetra_Map> roweles;
    dis[i]->BuildElementRowColumn(*newnoderowmap,*stdnodecolmap[i],roweles,stdelecolmap[i]);
    dis[i]->ExportRowNodes(*newnoderowmap);
    dis[i]->ExportRowElements(*roweles);
    dis[i]->ExportColumnNodes(*stdnodecolmap[i]);
    dis[i]->ExportColumnElements(*stdelecolmap[i]);
    dis[i]->FillComplete(false,false,false);
    if(myrank_ == 0)
      std::cout << "parallel distribution with standard ghosting" << std::endl;
    DRT::UTILS::PrintParallelDistribution(*dis[i]);

    //----------------------------
    // start with extended ghosting
    //----------------------------
    // fill elements into bins
    std::map<int, std::set<int> > binelemap;
    DistributeElesToBins(dis[i], binelemap);

    // ghosting is extended
    std::map<int, std::set<int> > dummy;
    Teuchos::RCP<Epetra_Map> extendedelecolmap = ExtendGhosting(dis[i]->ElementColMap(), binelemap, dummy);

    // adapt layout to extended ghosting in discret
    // first export the elements according to the processor local element column maps
    dis[i]->ExportColumnElements(*extendedelecolmap);

    // get the node ids of the elements that are to be ghosted and create a proper node column map for their export
    std::set<int> nodes;
    for (int lid=0; lid<extendedelecolmap->NumMyElements(); ++lid)
    {
      DRT::Element* ele = dis[i]->gElement(extendedelecolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0; inode<ele->NumNode(); ++inode)
        nodes.insert(nodeids[inode]);
    }

    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,dis[i]->Comm()));

    // now ghost the nodes
    dis[i]->ExportColumnNodes(*nodecolmap);

    // fillcomplete discret with extended ghosting
    dis[i]->FillComplete();
    if(myrank_ == 0)
      std::cout << "parallel distribution with extended ghosting" << std::endl;
    DRT::UTILS::PrintParallelDistribution(*dis[i]);
  }

  return;
}


/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 11/13 |
 *-------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
  DRT::Discretization& mortardis,
  Teuchos::RCP<Epetra_Map> initial_elecolmap,
  std::map<int, std::set<int> >& slavebinelemap,
  std::map<int, std::set<int> >& masterbinelemap)
{
  std::map<int, std::set<int> > extendedghosting;

  // do communication to gather all elements for extended ghosting
  const int numproc = mortardis.Comm().NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // get all neighboring bins around bins that contain slave elements
    std::set<int> binset;
    if(iproc == myrank_)
    {
      for(std::map<int, std::set<int> >::const_iterator iter=slavebinelemap.begin(); iter!=slavebinelemap.end(); ++iter)
      {
        int binId = iter->first;
        std::vector<int> bins;
        // get neighboring bins
        GetBinConnectivity(binId, bins);
        binset.insert(bins.begin(), bins.end());
        // insert bin itself
        binset.insert(binId);
      }
    }
    // copy set to vector in order to broadcast data
    std::vector<int> binids(binset.begin(),binset.end());

    // first: proc i tells all procs how many bins it has
    int numbin = binids.size();
    mortardis.Comm().Broadcast(&numbin, 1, iproc);
    // second: proc i tells all procs which bins it has
    binids.resize(numbin);

    mortardis.Comm().Broadcast(&binids[0], numbin, iproc);

    // loop over all own bins and find requested ones, fill in master elements in these bins
    std::map<int, std::set<int> > sdata;
    std::map<int, std::set<int> > rdata;

    for(int i=0; i<numbin; ++i)
    {
      sdata[binids[i]].insert(masterbinelemap[binids[i]].begin(),masterbinelemap[binids[i]].end());
    }

    LINALG::Gather<int>(sdata, rdata, 1, &iproc, mortardis.Comm());

    // proc i has to store the received data
    if(iproc == myrank_)
    {
      extendedghosting = rdata;
    }
  }

  // reduce map of sets to one set and copy to a vector to create extended elecolmap
  std::set<int> mastereleset;
  std::map<int, std::set<int> >::iterator iter;
  for(iter=extendedghosting.begin(); iter!= extendedghosting.end(); ++iter)
  {
    mastereleset.insert(iter->second.begin(),iter->second.end());
  }

  // insert standard ghosting for master and slave side
  for(int lid=0; lid<initial_elecolmap->NumMyElements(); ++lid)
  {
    mastereleset.insert(initial_elecolmap->GID(lid));
  }

  std::vector<int> mastercolgids(mastereleset.begin(),mastereleset.end());

  // return extendedmastercolmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)mastercolgids.size(),&mastercolgids[0],0,mortardis.Comm()));
}


/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 02/14 |
 *-------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ExtendGhosting(
  Teuchos::RCP<DRT::Discretization> scatradis,
  std::map<int, std::set<int> >& escapedpartelemap,
  std::map<int, std::set<int> >& myescapedpartelemap)
{
  // get fully redundant map of escaped particles
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // CAUTION: We chose this way here, since we expect that the map escapedpartelemap only contains
  //          a small faction of elements!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LINALG::GatherAll(escapedpartelemap,scatradis->Comm());

  // filter particles required on this proc
  std::set<int> myescapedparticles;
  for(std::map<int, std::set<int> >::iterator iter=escapedpartelemap.begin(); iter!= escapedpartelemap.end(); ++iter)
  {
    if (scatradis->HaveGlobalElement(iter->first))
    {
      // gather all particles to create col map of particles
      myescapedparticles.insert(iter->second.begin(),iter->second.end());

      // insert data to map of all elements on this proc with escaped particles
      myescapedpartelemap[iter->first].insert(iter->second.begin(),iter->second.end());
    }
  }

  // insert standard row particle distribution
  const Epetra_Map* particlerowmap = particledis_->NodeRowMap();
  for(int lid=0; lid<particlerowmap->NumMyElements(); ++lid)
    myescapedparticles.insert(particlerowmap->GID(lid));

  // copy to a vector and create extended particle colmap
  std::vector<int> myparticlecolgids(myescapedparticles.begin(),myescapedparticles.end());
  Teuchos::RCP<Epetra_Map> particlecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)myparticlecolgids.size(),&myparticlecolgids[0],0,particledis_->Comm()));

  // now ghost the nodes
  particledis_->ExportColumnNodes(*particlecolmap);

  // call fillcomplete
  particledis_->FillComplete();

  return;
}


/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 06/14 |
 *-------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
  const Epetra_Map* initial_elecolmap,
  std::map<int, std::set<int> >& binelemap,
  std::map<int, std::set<int> >& extendedghosting,
  Teuchos::RCP<Epetra_Map> bincolmap)
{
  // do communication to gather all elements for extended ghosting
  const int numproc = initial_elecolmap->Comm().NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // gather set of column bins for each proc
    std::set<int> binset;
    if(iproc == myrank_)
    {
      // either use given column layout of bins ...
      if(bincolmap != Teuchos::null)
      {
        int nummyeles = bincolmap->NumMyElements();
        int* entries = bincolmap->MyGlobalElements();
        binset.insert(entries, entries+nummyeles);
      }
      else // ... or add an extra layer to the given bin distribution
      {
        for(std::map<int, std::set<int> >::const_iterator iter=binelemap.begin(); iter!=binelemap.end(); ++iter)
        {
          int binId = iter->first;
          std::vector<int> bins;
          // get neighboring bins
          GetBinConnectivity(binId, bins);
          binset.insert(bins.begin(), bins.end());
          // insert bin itself
          binset.insert(binId);
        }
      }
    }
    // copy set to vector in order to broadcast data
    std::vector<int> binids(binset.begin(),binset.end());

    // first: proc i tells all procs how many bins it has
    int numbin = binids.size();
    initial_elecolmap->Comm().Broadcast(&numbin, 1, iproc);
    // second: proc i tells all procs which bins it has
    binids.resize(numbin);

    initial_elecolmap->Comm().Broadcast(&binids[0], numbin, iproc);

    // loop over all own bins and find requested ones, fill in master elements in these bins
    std::map<int, std::set<int> > sdata;
    std::map<int, std::set<int> > rdata;

    for(int i=0; i<numbin; ++i)
    {
      if(binelemap.find(binids[i]) != binelemap.end())
        sdata[binids[i]].insert(binelemap[binids[i]].begin(),binelemap[binids[i]].end());
    }

    LINALG::Gather<int>(sdata, rdata, 1, &iproc, initial_elecolmap->Comm());

    // proc i has to store the received data
    if(iproc == myrank_)
    {
      extendedghosting = rdata;
    }
  }

  // reduce map of sets to one set and copy to a vector to create extended elecolmap
  std::set<int> coleleset;
  std::map<int, std::set<int> >::iterator iter;
  for(iter=extendedghosting.begin(); iter!= extendedghosting.end(); ++iter)
  {
    coleleset.insert(iter->second.begin(),iter->second.end());
  }

  // insert standard ghosting
  for(int lid=0; lid<initial_elecolmap->NumMyElements(); ++lid)
  {
    coleleset.insert(initial_elecolmap->GID(lid));
  }

  std::vector<int> colgids(coleleset.begin(),coleleset.end());

  // return extended elecolmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)colgids.size(),&colgids[0],0,initial_elecolmap->Comm()));
}


/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 11/13 |
 *-------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ExtendGhosting(
    std::vector<Teuchos::RCP<DRT::Discretization> > dis)
{
  for(size_t i=0; i<dis.size(); ++i)
  {
    //----------------------------
    // start with extended ghosting
    //----------------------------
    // fill elements into bins
    std::map<int, std::set<int> > binelemap;
    DistributeElesToBins(dis[i], binelemap);

    // ghosting is extended
    std::map<int, std::set<int> > dummy;
    Teuchos::RCP<Epetra_Map> extendedelecolmap = ExtendGhosting(dis[i]->ElementColMap(), binelemap, dummy);

    // adapt layout to extended ghosting in discret
    // first export the elements according to the processor local element column maps
    dis[i]->ExportColumnElements(*extendedelecolmap);

    // get the node ids of the elements that are to be ghosted and create a proper node column map for their export
    std::set<int> nodes;
    for (int lid=0; lid<extendedelecolmap->NumMyElements(); ++lid)
    {
      DRT::Element* ele = dis[i]->gElement(extendedelecolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0; inode<ele->NumNode(); ++inode)
        nodes.insert(nodeids[inode]);
    }

    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,dis[i]->Comm()));

    // now ghost the nodes
    dis[i]->ExportColumnNodes(*nodecolmap);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::RevertExtendedGhosting(
  std::vector<Teuchos::RCP<DRT::Discretization> > dis,
  std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
  std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap
  )
{
  for(size_t i=0; i<dis.size(); ++i)
  {
    //----------------------------
    // revert extended ghosting
    //----------------------------

    // adapt layout to standard ghosting in discret
    // first export the elements according to the processor local element column maps
    dis[i]->ExportColumnElements(*(stdelecolmap[i]));

    // now ghost the nodes
    dis[i]->ExportColumnNodes(*(stdnodecolmap[i]));

    // fillcomplete discret with standard ghosting
    dis[i]->FillComplete();
    if(myrank_ == 0)
      std::cout << "parallel distribution with reverted ghosting" << std::endl;
    DRT::UTILS::PrintParallelDistribution(*dis[i]);
  }
}

/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 11/13 |
 *-------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CollectInformation(
  Teuchos::RCP<Epetra_Map> rowbins,
  std::map<int, std::vector<int> >& nodesinbin,
  std::map<int, std::vector<int> >& mynodesinbin)
{
  // do communication to gather all nodes
  const int numproc = rowbins->Comm().NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // vector with row bins on this proc
    std::vector<int> binids;
    int numbin;
    if(iproc == myrank_)
    {
      int* myrowbinsdata = rowbins->MyGlobalElements();
      numbin = rowbins->NumMyElements();
      binids.insert(binids.begin(), myrowbinsdata, myrowbinsdata+numbin);
    }

    // first: proc i tells all procs how many bins it has
    rowbins->Comm().Broadcast(&numbin, 1, iproc);
    binids.resize(numbin);
    // second: proc i tells all procs which bins it has
    rowbins->Comm().Broadcast(&binids[0], numbin, iproc);

    // loop over all own bins and find requested ones, fill in master elements in these bins
    std::map<int, std::vector<int> > sdata;
    std::map<int, std::vector<int> > rdata;

    for(int i=0; i<numbin; ++i)
    {
      if(nodesinbin.find(binids[i]) != nodesinbin.end())
        sdata[binids[i]].insert(sdata[binids[i]].begin(), nodesinbin[binids[i]].begin(), nodesinbin[binids[i]].end());
    }

    LINALG::Gather<int>(sdata, rdata, 1, &iproc, rowbins->Comm());

    // proc i has to store the received data
    if(iproc == myrank_)
    {
      // clear data and refill
      mynodesinbin.clear();
      mynodesinbin.insert(rdata.begin(), rdata.end());
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateBins(Teuchos::RCP<DRT::Discretization> dis)
{
  if(dis != Teuchos::null)
  {
    // create XAABB for discretization
    CreateXAABB(dis);
  }

  // divide global bounding box into bins
  for (int dim=0; dim<3; ++dim)
  {
    // determine number of bins per direction for prescribed cutoff radius
    // std::floor leads to bins that are at least of size cutoff_radius
    if (cutoff_radius_>0.0)
      bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));

    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

  if(particle_dim_ != INPAR::PARTICLE::particle_3D)
  {
    int entry = -1;
    switch (particle_dim_)
    {
    case INPAR::PARTICLE::particle_2Dx:
      entry = 0;
      break;
    case INPAR::PARTICLE::particle_2Dy:
      entry = 1;
      break;
    case INPAR::PARTICLE::particle_2Dz:
      entry = 2;
      break;
    default:
      dserror("number of particle dimensions not yet implemented");
      break;
    }

    // one bin in pseudo direction is enough
    bin_per_dir_[entry] = 1;
    bin_size_[entry] = (XAABB_(entry,1)-XAABB_(entry,0))/bin_per_dir_[entry];
  }

//  if(myrank_ == 0)
//  {
//    std::cout << "Global bounding box size: " << XAABB_
//        << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << std::endl;
//  }

  return;
}


/*----------------------------------------------------------------------*
| find XAABB and cutoff                                     ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateXAABB(
  Teuchos::RCP<DRT::Discretization> dis,
  LINALG::Matrix<3,2>& XAABB,
  double& cutoff
  )
{

  double locmaxcutoff = 0.0;
  // initialize XAABB as rectangle around the first node of dis
  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim, 0) = dis->lRowNode(0)->X()[dim] - GEO::TOL7;
    XAABB(dim, 1) = dis->lRowNode(0)->X()[dim] + GEO::TOL7;
  }

  // loop
  for(int i=0; i<dis->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = dis->lRowElement(i);

    LINALG::Matrix<3,2> eleXAABB(false);

     {
       // initialize XAABB as rectangle around the first node of ele
       for(int dim=0; dim<3; ++dim)
       {
         eleXAABB(dim, 0) = ele->Nodes()[0]->X()[dim] - GEO::TOL7;
         eleXAABB(dim, 1) = ele->Nodes()[0]->X()[dim] + GEO::TOL7;
       }
     }

     // loop over remaining nodes and merge XAABB with their eXtendedAxisAlignedBoundingBox
     for (int lid = 1; lid < ele->NumNode(); ++lid)
     {
       const DRT::Node* node = ele->Nodes()[lid];

       for(int dim=0; dim < 3; dim++)
       {
         eleXAABB(dim, 0) = std::min( eleXAABB(dim, 0), node->X()[dim] - GEO::TOL7);
         eleXAABB(dim, 1) = std::max( eleXAABB(dim, 1), node->X()[dim] + GEO::TOL7);
       }
     }

     for(int dim=0; dim<3; ++dim)
     {
       locmaxcutoff = std::max( locmaxcutoff, eleXAABB(dim, 1) - eleXAABB(dim, 0));
     }

     // merge XAABB of elements

     for(int dim=0; dim < 3; dim++)
     {
       XAABB(dim, 0) = std::min( XAABB(dim, 0), eleXAABB(dim, 0));
       XAABB(dim, 1) = std::max( XAABB(dim, 1), eleXAABB(dim, 1));
     }
  }

  // local bounding box
  double locmin[3] = {XAABB(0,0), XAABB(1,0), XAABB(2,0)};
  double locmax[3] = {XAABB(0,1), XAABB(1,1), XAABB(2,1)};
  // global bounding box
  double globmin[3];
  double globmax[3];
  // do the necessary communication
  dis->Comm().MinAll(&locmin[0], &globmin[0], 3);
  dis->Comm().MaxAll(&locmax[0], &globmax[0], 3);

  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim,0) = globmin[dim];
    XAABB(dim,1) = globmax[dim];
  }

  // maxall of cutoff
  cutoff = 0.0;
  dis->Comm().MaxAll(&locmaxcutoff, &cutoff, 1);

  return;
}


/*----------------------------------------------------------------------*
| find XAABB                                                ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateXAABB(
  Teuchos::RCP<DRT::Discretization> dis
  )
{
  // if not yet specified, get XAABB_ from underlying discretization
  if( XAABB_(2,1) > 0.9e12  and  XAABB_(2,1) < 1.1e12 )
  {
    if(myrank_ == 0)
      std::cout << "XAABB is computed based on the underlying discretization" << std::endl;
    XAABB_ = GEO::getXAABBofNodes(*dis);
    // local bounding box
    double locmin[3] = {XAABB_(0,0), XAABB_(1,0), XAABB_(2,0)};
    double locmax[3] = {XAABB_(0,1), XAABB_(1,1), XAABB_(2,1)};
    // global bounding box
    double globmin[3];
    double globmax[3];
    // do the necessary communication
    dis->Comm().MinAll(&locmin[0], &globmin[0], 3);
    dis->Comm().MaxAll(&locmax[0], &globmax[0], 3);

    for(int dim=0; dim<3; ++dim)
    {
      XAABB_(dim,0) = globmin[dim];
      XAABB_(dim,1) = globmax[dim];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 01/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid(const std::vector<double>& pos)
{
  int ijk[3] = {0,0,0};
  for(int dim=0; dim < 3; dim++)
  {
    ijk[dim] = (int)(std::floor((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]));
  }

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 01/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid(const double* pos)
{
  int ijk[3];
  for(int dim=0; dim<3; ++dim)
  {
    ijk[dim] = (int)(std::floor((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]));
  }

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 02/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertPosToijk(const double* pos, int* ijk)
{
  for(int dim=0; dim<3; ++dim)
  {
    ijk[dim] = (int)(std::floor((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]));
  }
  return;
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid(const LINALG::Matrix<3,1> pos)
{
  int ijk[3];
  for(int dim=0; dim<3; ++dim)
  {
    ijk[dim] = (int)(std::floor((pos(dim)-XAABB_(dim,0)) / bin_size_[dim]));
  }

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertPosToijk(const LINALG::Matrix<3,1> pos, int* ijk)
{
  for(int dim=0; dim<3; ++dim)
  {
    ijk[dim] = (int)(std::floor((pos(dim)-XAABB_(dim,0)) / bin_size_[dim]));
  }
  return;
}


/*----------------------------------------------------------------------*
| convert i,j,k into bin id                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertijkToGid(int* ijk)
{
  // given ijk is outside of XAABB
  if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 || ijk[0]>=bin_per_dir_[0] || ijk[1]>=bin_per_dir_[1] || ijk[2]>=bin_per_dir_[2] )
    return -1;

  return ijk[0] + ijk[1]*bin_per_dir_[0] + ijk[2]*bin_per_dir_[0]*bin_per_dir_[1];
}


/*----------------------------------------------------------------------*
| convert bin id into i,j,k                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertGidToijk(int gid, int* ijk)
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
void BINSTRATEGY::BinningStrategy::GidsInijkRange(int* ijk_range, std::set<int>& binIds, bool checkexistence)
{
  if(checkexistence == true and particledis_ == Teuchos::null)
    dserror("particle discretization is not set up correctly");

  for(int i=ijk_range[0]; i<=ijk_range[1]; ++i)
  {
    for(int j=ijk_range[2]; j<=ijk_range[3]; ++j)
    {
      for(int k=ijk_range[4]; k<=ijk_range[5]; ++k)
      {
        int ijk[3] = {i,j,k};

        int gid = ConvertijkToGid(&ijk[0]);
        if(gid != -1)
        {
          if(checkexistence)
          {
            if(particledis_->HaveGlobalElement(gid))
              binIds.insert(gid);
          }
          else
          {
            binIds.insert(gid);
          }
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return;
}


/*----------------------------------------------------------------------*
 | get all bins in ijk range                               ghamm 03/16  |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GidsInijkRange(int* ijk_range, std::vector<int>& binIds, bool checkexistence)
{
  if(checkexistence == true and particledis_ == Teuchos::null)
    dserror("particle discretization is not set up correctly");

  for(int i=ijk_range[0]; i<=ijk_range[1]; ++i)
  {
    for(int j=ijk_range[2]; j<=ijk_range[3]; ++j)
    {
      for(int k=ijk_range[4]; k<=ijk_range[5]; ++k)
      {
        int ijk[3] = {i,j,k};

        int gid = ConvertijkToGid(&ijk[0]);
        if(gid != -1)
        {
          if(checkexistence)
          {
            if(particledis_->HaveGlobalElement(gid))
              binIds.push_back(gid);
          }
          else
          {
            binIds.push_back(gid);
          }
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return;
}


/*----------------------------------------------------------------------*
 | get 26 neighboring bin ids to binId (if existing)       ghamm 08/13  |
*-----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetBinConnectivity(int binId, std::vector<int>& binIds)
{
  int ijk_base[3];
  ConvertGidToijk(binId, &ijk_base[0]);
  for(int i=ijk_base[0]-1; i<=ijk_base[0]+1; ++i)
  {
    for(int j=ijk_base[1]-1; j<=ijk_base[1]+1; ++j)
    {
      for(int k=ijk_base[2]-1; k<=ijk_base[2]+1; ++k)
      {
        int ijk[3] = {i,j,k};
        int gid = ConvertijkToGid(&ijk[0]);
        if(gid!=-1 and gid!=binId)
        {
          binIds.push_back(gid);
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return;
}


/*----------------------------------------------------------------------*
| corner position for given bin id                          ghamm 03/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetBinCorners(int binId, std::vector<LINALG::Matrix<3,1> >& bincorners)
{
  bincorners.clear();
  bincorners.reserve(8);
  int ijk_base[3];
  ConvertGidToijk(binId, &ijk_base[0]);

  // order in bincorners is identical to ordering of i,j and k
  for(int k=ijk_base[2]; k<(ijk_base[2]+2); ++k)
  {
    for(int j=ijk_base[1]; j<(ijk_base[1]+2); ++j)
    {
      for(int i=ijk_base[0]; i<(ijk_base[0]+2); ++i)
      {
        int ijk_curr[] = {i,j,k};
        LINALG::Matrix<3,1> curr_corner;
        for(int dim=0; dim<3; ++dim)
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
LINALG::Matrix<3,1> BINSTRATEGY::BinningStrategy::GetBinCentroid(int binId)
{
  int ijk[3];
  ConvertGidToijk(binId, ijk);
  if(ijk[0] == -1)
    dserror("given bin id is outside of bins; centroid of bin is does not make sense");

  LINALG::Matrix<3,1> centroid;
  for(int dim=0; dim<3; ++dim)
    centroid(dim) = XAABB_(dim,0) + bin_size_[dim]*(ijk[dim] + 0.5);

  return centroid;
}


/*----------------------------------------------------------------------*
 | Class for comparing Teuchos::RCP<DRT::Node> in std::set ghamm 10/12  |
 *----------------------------------------------------------------------*/
bool BINSTRATEGY::Less::operator()(const Teuchos::RCP<const DRT::Node>& first, const Teuchos::RCP<const DRT::Node>& second) const
{
  return first->Id() < second->Id();
}
