/*----------------------------------------------------------------------*/
/*!
\file binning_strategy.cpp

\brief Binning strategy for neighborhood search

\level 2

\maintainer Georg Hammerl
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
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"

#include "binning_strategy.H"

/*----------------------------------------------------------------------*
 | Binning strategy constructor                             ghamm 11/13 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy(
  const Epetra_Comm& comm,
  double cutoff_radius,
  LINALG::Matrix<3,2> XAABB
  ) :
  bindis_(Teuchos::null),
  cutoff_radius_(cutoff_radius),
  XAABB_(XAABB),
  havepbc_(false),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")),
  myrank_(comm.MyPID())
{
  if( XAABB_(0,0) >= XAABB_(0,1) or XAABB_(1,0) >= XAABB_(1,1) or XAABB_(2,0) >= XAABB_(2,1))
    dserror("XAABB is not computed correctly");

  if(cutoff_radius_ <= 0.0)
    dserror("Cutoff radius cannot be zero or negative!");

  // initialize arrays
  for(unsigned idim=0; idim<3; ++idim)
  {
    bin_size_[idim] = 0.0;
    inv_bin_size_[idim] = 0.0;
    bin_per_dir_[idim] = 0;
    pbconoff_[idim] = false;
    pbcdeltas_[idim] = 0.0;
  }

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
  bindis_(Teuchos::null),
  havepbc_(false),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")),
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

  // initialize arrays
  for(int idim=0; idim<3; ++idim)
  {
    bin_size_[idim] = 0.0;
    inv_bin_size_[idim] = 0.0;
    pbconoff_[idim] = false;
    pbcdeltas_[idim] = 0.0;
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
  bindis_(Teuchos::null),
  cutoff_radius_(0.0),
  havepbc_(false),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION")),
  myrank_(dis[0]->Comm().MyPID())
{
  // initialize arrays
  for(int idim=0; idim<3; ++idim)
  {
    bin_size_[idim] = 0.0;
    inv_bin_size_[idim] = 0.0;
    bin_per_dir_[idim] = 0;
    pbconoff_[idim] = false;
    pbcdeltas_[idim] = 0.0;
  }

  WeightedPartitioning(dis,stdelecolmap,stdnodecolmap);

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
  Teuchos::RCP<DRT::Discretization> discret,
  std::map<int, std::set<int> >&    rowelesinbin,
  Teuchos::RCP<Epetra_Vector>       disnp
  )
{
  // current node position
  double currpos[3] = {0.0,0.0,0.0};
  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for (int lid = 0; lid < discret->NumMyRowElements(); ++lid)
  {
    DRT::Element* ele = discret->lRowElement(lid);
    DRT::Node** nodes = ele->Nodes();
    const int numnode = ele->NumNode();

    // initialize ijk_range with ijk of first node of element
    int ijk[3];
    {
      const DRT::Node* node = nodes[0];
      GetCurrentNodePos(discret,node,disnp,currpos);
      const double* coords = currpos;
      ConvertPosToijk(coords, ijk);
    }

    // ijk_range contains: i_min i_max j_min j_max k_min k_max
    int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

    // fill in remaining nodes
    for (int j=1; j<numnode; ++j)
    {
      const DRT::Node* node = nodes[j];
      GetCurrentNodePos(discret,node,disnp,currpos);
      const double* coords = currpos;
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

/*-----------------------------------------------------------------------------*
| assign beam elements (taking into account cut elements due to periodic       |
| boundary conditions) into bins                               eichinger 09/16 |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeCutElesToBins(
  Teuchos::RCP<DRT::Discretization> discret,
  std::map<int, std::set<int> >&    rowelesinbin,
  Teuchos::RCP<Epetra_Vector>       disnp
  )
{
  // current node position
  double currpos[3] = {0.0,0.0,0.0};
  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for (int lid = 0; lid < discret->NumMyRowElements(); ++lid)
  {
    DRT::Element* ele = discret->lRowElement(lid);
    DRT::Node** nodes = ele->Nodes();
    const int numnode = ele->NumNode();

    // initialize ijk_range with ijk of first node of element
    int ijk[3];
    {
      const DRT::Node* node = nodes[0];
      GetCurrentNodePos(discret,node,disnp,currpos);
      const double* coords = currpos;
      ConvertPosToijk(coords, ijk);
    }

    // ijk_range contains: i_min i_max j_min j_max k_min k_max
    int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

    // fill in remaining nodes
    for (int j=1; j<numnode; ++j)
    {
      const DRT::Node* node = nodes[j];
      GetCurrentNodePos(discret,node,disnp,currpos);
      const double* coords = currpos;
      int ijk[3];
      ConvertPosToijk(coords, ijk);

      for(int dim=0; dim<3; ++dim)
      {
        if(ijk[dim]<ijk_range[dim*2])
        {
          if((ijk[dim] == 0) && (abs(ijk[dim]-ijk_range[dim*2])>4))
          {
            ijk_range[dim*2+1]=bin_per_dir_[dim];
            continue;
          }
          else
            ijk_range[dim*2]=ijk[dim];
        }
        if(ijk[dim]>ijk_range[dim*2+1])
        {
          if((ijk[dim] == bin_per_dir_[dim] - 1) &&(abs(ijk[dim]-ijk_range[dim*2+1])>4))
            ijk_range[dim*2]= -1;
          else
            ijk_range[dim*2+1]=ijk[dim];
        }
      } // loop over dim
    } // loop over nodes > 1

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
| assign elements into bins                                 ghamm 08/16 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::AssignElesToBins(
    Teuchos::RCP<DRT::Discretization> discret,
    std::map<int, std::set<int> >     extendedfieldghosting,
    INPAR::BINSTRATEGY::BinContent    bincontent
  )
{
  // loop over bins
  std::map<int, std::set<int> >::const_iterator biniter;
  for(biniter = extendedfieldghosting.begin(); biniter!=extendedfieldghosting.end(); ++biniter)
  {
    // get current bin
    DRT::MESHFREE::MeshfreeMultiBin* currbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(bindis_->gElement(biniter->first));

    // loop over ele content of this bin
    std::set<int>::const_iterator eleiter;
    for(eleiter = biniter->second.begin(); eleiter!=biniter->second.end(); ++eleiter)
    {
      int eleid = *eleiter;
      // add eleid and elepointer to current bin
      currbin->AddAssociatedEle(
          bincontent,eleid, discret->gElement(eleid));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
| assign elements into bins                                 ghamm 08/16 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::RemoveElesFromBins(
    INPAR::BINSTRATEGY::BinContent bincontent
  )
{
  // loop over all bins and remove assigned elements
  const int numcolbins = bindis_->ElementColMap()->NumMyElements();
  for(int binlid=0; binlid<numcolbins; ++binlid)
  {
    DRT::Element *currentbin = bindis_->lColElement(binlid);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currentbin)->RemoveAssociatedEles(bincontent);
  }

  return;
}

/*----------------------------------------------------------------------*
| assign nodes into bins                                    ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeNodesToBins(
  Teuchos::RCP<DRT::Discretization> discret,
  std::map<int, std::vector<int> >& nodesinbin,
  Teuchos::RCP<Epetra_Vector> disnp
  )
{
  // current position of nodes
  double currpos[3] = {0.0,0.0,0.0};

  // loop over row nodes
  for (int lid = 0; lid < discret->NumMyRowNodes(); ++lid)
  {
    DRT::Node* node = discret->lRowNode(lid);
    GetCurrentNodePos(discret,node,disnp,currpos);

    const double* coords = currpos;
    int ijk[3];
    ConvertPosToijk(coords, ijk);
    const int binid = ConvertijkToGid(&ijk[0]);

    // assign node to bin
    nodesinbin[binid].push_back(node->Id());
  }

  return;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::WeightedPartitioning(
  std::vector<Teuchos::RCP<DRT::Discretization> > discret,
  std::vector<Teuchos::RCP<Epetra_Map> >&         stdelecolmap,
  std::vector<Teuchos::RCP<Epetra_Map> >&         stdnodecolmap
  )
{
  // initialize dummys
  std::vector<std::map<int, std::set<int> > > dummy1(discret.size());
  std::vector<Teuchos::RCP<Epetra_Vector> >   dummy2(discret.size());

  // create XAABB_ and set cutoff
  CreateXAABB(discret,dummy2,true);

  // create bins
  CreateBins(Teuchos::null);

  // ------------------------------------------------------------------------
  // create bins, weight them according to number of nodes (of discrets) they
  // contain, account for bin connectivity. Then an optimal distribution of
  // bins to procs can be obtained
  // ------------------------------------------------------------------------
  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int> > > nodesinbin(discret.size());
  Teuchos::RCP<Epetra_Map> newrowbins =
      WeightedDistributionOfBinsToProcs(discret,dummy2,nodesinbin,false);

  stdelecolmap.resize(discret.size());
  stdnodecolmap.resize(discret.size());

  // ------------------------------------------------------------------------
  // now we have an optimal distribution of bins (with respect to their content
  // and connectivity). Now we have to apply it by rebuilding the input discret,
  // i.e. we need to change the ownership of the nodes/elements according to
  // the bin they belong to (each proc then owns the nodes/eles laying in its
  // bins.
  // ------------------------------------------------------------------------
  // rebuild discretizations including extended ghosting
  for(size_t i=0; i<discret.size(); ++i)
  {
   // ----------------------------------------------------------------------
   // start with standard ghosting
   // ----------------------------------------------------------------------
    StandardGhosting(discret[i],newrowbins,dummy2[i],
        stdelecolmap[i],stdnodecolmap[i],nodesinbin[i]);

    // some output after standard ghosting
    if(myrank_ == 0)
      std::cout << "parallel distribution with standard ghosting" << std::endl;
    DRT::UTILS::PrintParallelDistribution(*discret[i]);

   // ----------------------------------------------------------------------
   // extended ghosting
   // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    // start with extended ghosting. Here this means the following: Each proc
    // ghosts all elements whose XAABB cuts a bin that is next to a bin that is
    // owned by a proc an not empty. All associated nodes are ghosted as well
    // ----------------------------------------------------------------------
    // here each proc assignes his owned elements in the means of a XAABB to
    // the global binids that do not need be owned by this proc.
    // binelemap on each proc than contains all bins (not neccesarily owned by
    // this proc) that are cut by the procs row elements
    std::map<int, std::set<int> > bintoelemap;
    DistributeElesToBins(discret[i],bintoelemap,dummy2[i]);

    // ghosting is extended to one layer (two layer ghosting is excluded as it
    // is not needed, this case is covered by other procs then) around bins that
    // actually contain elements.
    // extbintoelemap[i] than contains all bins and its corresponding elements
    // that need to be owned or ghosted to ensure correct interaction handling
    // of the elements in the range of one layer
    Teuchos::RCP<Epetra_Map> extendedelecolmap =
        ExtendGhosting(discret[i]->ElementColMap(), bintoelemap, dummy1[i], newrowbins);

    // adapt layout to extended ghosting in discret
    // first export the elements according to the processor local element column maps
    discret[i]->ExportColumnElements(*extendedelecolmap);

    // get the node ids of the elements that are to be ghosted
    // and create a proper node column map for their export
    std::set<int> nodes;
    for (int lid=0; lid<extendedelecolmap->NumMyElements(); ++lid)
    {
      DRT::Element* ele = discret[i]->gElement(extendedelecolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0; inode<ele->NumNode(); ++inode)
        nodes.insert(nodeids[inode]);
    }

    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap =
        Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,discret[i]->Comm()));

    // now ghost the nodes
    discret[i]->ExportColumnNodes(*nodecolmap);

    // fillcomplete discret with extended ghosting
    discret[i]->FillComplete();
    if(myrank_ == 0)
      std::cout << "parallel distribution with extended ghosting" << std::endl;
    DRT::UTILS::PrintParallelDistribution(*discret[i]);
  }

  return newrowbins;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::WeightedDistributionOfBinsToProcs(
  std::vector<Teuchos::RCP<DRT::Discretization> >  discret,
  std::vector<Teuchos::RCP<Epetra_Vector> >        disnp,
  std::vector<std::map<int, std::vector<int> > > & nodesinbin,
  bool repartition
  )
{
  // calculate total number of bins
  const int numbin = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];

  // some safety checks to ensure efficiency
  if(numbin<discret[0]->Comm().NumProc() && myrank_ == 0)
    dserror("ERROR:NumProc > NumBin. Too many processors to "
            "distribute your bins properly!!!");
  if(numbin < 8*discret[0]->Comm().NumProc() && myrank_==0)
    std::cout << "\n\nWARNING: partitioning not useful, choose less procs. "
                 " Owner distribution may be inefficient!\n\n" << std::endl;

  // row bin distribution
  Teuchos::RCP<Epetra_Map> rowbins = Teuchos::null;
  Teuchos::RCP< Epetra_CrsGraph> bingraph;
  if(repartition)
  {
    // use old bin distribution
    rowbins = Teuchos::rcp(const_cast<Epetra_Map*>(bindis_->ElementRowMap()));
    const Epetra_Map* oldrowmap = bindis_->ElementRowMap();

    const int maxband = 26;
    bingraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*oldrowmap,maxband,false));

    // fill all local entries into the graph
    {
      for (int lid=0; lid<oldrowmap->NumMyElements(); ++lid)
      {
        const int binId = oldrowmap->GID(lid);

        std::vector<int> neighbors;
        GetBinConnectivity(binId,neighbors);

        int err = bingraph->InsertGlobalIndices(binId,(int)neighbors.size(),&neighbors[0]);
        if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,binId);
      }
    }
  }
  else
  {
    // dummy row bin distribution (equally distributed over all procs as no
    // weighting done so far)
    rowbins = Teuchos::rcp(new Epetra_Map(numbin,0,discret[0]->Comm()));
    // create nodal graph
     bingraph = Teuchos::rcp( new Epetra_CrsGraph(Copy,*rowbins,108,false));
  }

   // Now we're going to create a Epetra_Vector with vertex/node weights to be
   // used for the partitioning operation (weights must be at least one for zoltan)
   Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*rowbins, true);

   // set weights of bins related to the number of nodes of discrets that are contained
   // empty bins have weight of 1
   vweights->PutScalar(1.0);

   // determine which node is in which bin and weight each bin according to
   // 10 times the number of nodes it contains
   // assign all node gids to their corresponding bin
  for(size_t i=0; i<discret.size(); ++i)
  {
    // distribute nodes, that are owned by a proc, to the bins of this proc
    DistributeNodesToBins(discret[i], nodesinbin[i], disnp[i]);

    std::map<int, std::vector<int> > mynodesinbin;
    // gather information of bin content from other procs (bin is owned by this
    // proc and there are some nodes on other procs which are located in this bin)
    CollectInformation(rowbins, nodesinbin[i], mynodesinbin);

    // weight each bin with 10 times the number of node it contains
    // empty bins remain with weight one
    std::map<int, std::vector<int> >::const_iterator biniter;
    for(biniter=mynodesinbin.begin(); biniter!=mynodesinbin.end(); ++biniter)
    {
      int lid = rowbins->LID(biniter->first);
      if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",discret[i]->Comm().MyPID(),biniter->first);
      (*vweights)[lid] += 10.0*(double)biniter->second.size();
    }
  }

   // fill bin connectivity into bin graph
   for (int lid=0; lid<rowbins->NumMyElements(); ++lid)
   {
     int rowbinid = rowbins->GID(lid);
     // insert 26 (one level) neighboring bins to graph
     // (if active, periodic boundary conditions are considered here)
     std::vector<int> neighbors;
     GetBinConnectivity(rowbinid,neighbors);

     int err = bingraph->InsertGlobalIndices(rowbinid,(int)neighbors.size(),&neighbors[0]);
     if (err<0)
       dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,rowbinid);
   }

   // complete graph
   int err = bingraph->FillComplete();
   if (err) dserror("graph->FillComplete() returned err=%d",err);
   err = bingraph->OptimizeStorage();
   if (err) dserror("graph->OptimizeStorage() returned err=%d",err);

   // call redistribution of bin graph using bin weights
  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
      Teuchos::rcp(new Isorropia::Epetra::CostDescriber);
  costs->setVertexWeights(vweights);

  Teuchos::ParameterList paramlist;
  paramlist.set("PARTITIONING METHOD", "GRAPH");
  Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  if(repartition)
    sublist.set("LB_APPROACH", "REPARTITION");
  else
    sublist.set("LB_APPROACH", "PARTITION");

  Teuchos::RCP<const Epetra_CrsGraph> constbingraph(bingraph);

  // Now create the partitioner object
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(constbingraph, costs, paramlist));

  Isorropia::Epetra::Redistributor rd(partitioner);

  // redistribute bingraph
  Teuchos::RCP<Epetra_CrsGraph> balanced_bingraph = rd.redistribute(*bingraph);

  // extract repartitioned bin row map
  const Epetra_BlockMap& rbinstmp = balanced_bingraph->RowMap();
  Teuchos::RCP<Epetra_Map> newrowbins =
      Teuchos::rcp(new Epetra_Map(-1,rbinstmp.NumMyElements(),rbinstmp.MyGlobalElements(),0,discret[0]->Comm()));

  return newrowbins;
}

/*--------------------------------------------------------------------------*
| standard ghosting according to bin distribution               ghamm 06/14 |
 *--------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::StandardGhosting(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<Epetra_Map>          rowbins,
    Teuchos::RCP<Epetra_Vector>&      disnp,
    Teuchos::RCP<Epetra_Map>&         stdelecolmap,
    Teuchos::RCP<Epetra_Map>&         stdnodecolmap,
    std::map<int, std::vector<int> >  nodesinbin)
{
  // ----------------------------------------------------------------------
  // start with standard ghosting
  // ----------------------------------------------------------------------
  // each owner of a bin gets owner of the nodes this bin contains
  // all other nodes of elements, of which proc is owner of at least one
  // node, are ghosted
  Teuchos::RCP<Epetra_CrsGraph> initgraph = discret->BuildNodeGraph();

  // distribute nodes, that are owned by a proc, to the bins of this proc
  if(nodesinbin.empty())
    DistributeNodesToBins(discret, nodesinbin, disnp);

  std::map<int, std::vector<int> > mynodesinbin;
  // gather information of bin content from other procs (bin is owned by
  // this proc and there are some nodes on other procs which are located
  // in this bin here)
  CollectInformation(rowbins, nodesinbin, mynodesinbin);

  // build new node row map
  std::vector<int> mynewrownodes;
  std::map<int, std::vector<int> >::const_iterator biniter;
  for(biniter=mynodesinbin.begin(); biniter!=mynodesinbin.end(); ++biniter)
  {
    std::vector<int>::const_iterator nodeiter;
    for(nodeiter=biniter->second.begin(); nodeiter!=biniter->second.end(); ++nodeiter)
    {
      mynewrownodes.push_back(*nodeiter);
    }
  }
  mynodesinbin.clear();

  Teuchos::RCP<Epetra_Map> newnoderowmap =
      Teuchos::rcp(new Epetra_Map(-1,mynewrownodes.size(),&mynewrownodes[0],0,discret->Comm()));

  // create the new graph and export to it
  Teuchos::RCP<Epetra_CrsGraph> newnodegraph;

  newnodegraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*newnoderowmap,108,false));
  Epetra_Export exporter2(initgraph->RowMap(),*newnoderowmap);
  int err = newnodegraph->Export(*initgraph,exporter2,Add);
  if (err<0)
    dserror("Graph export returned err=%d",err);
  newnodegraph->FillComplete();
  newnodegraph->OptimizeStorage();

  // the column map will become the new ghosted distribution of nodes (standard ghosting)
  const Epetra_BlockMap cntmp = newnodegraph->ColMap();
  stdnodecolmap =
      Teuchos::rcp(new Epetra_Map(-1,cntmp.NumMyElements(),cntmp.MyGlobalElements(),0,discret->Comm()));

  // rebuild of the discretizations with new maps for standard ghosting
  Teuchos::RCP<Epetra_Map> roweles;
  discret->BuildElementRowColumn(*newnoderowmap,*stdnodecolmap,roweles,stdelecolmap);
  discret->ExportRowNodes(*newnoderowmap);
  discret->ExportRowElements(*roweles);
  discret->ExportColumnNodes(*stdnodecolmap);
  discret->ExportColumnElements(*stdelecolmap);
  // in case we have a state vector, we need to build the dof map to enable its rebuild
  if(disnp==Teuchos::null)
  {
    discret->FillComplete(false,false,false);
  }
  else
  {
    discret->FillComplete(true,false,false);
    Teuchos::RCP<Epetra_Vector> old;
    old = disnp;
    disnp = LINALG::CreateVector(*discret->DofRowMap(),true);
    LINALG::Export(*old, *disnp);
  }

  return;
}

/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 11/13 |
 *-------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
  DRT::Discretization&           mortardis,
  Teuchos::RCP<Epetra_Map>       initial_elecolmap,
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
  const Epetra_Map* particlerowmap = bindis_->NodeRowMap();
  for(int lid=0; lid<particlerowmap->NumMyElements(); ++lid)
    myescapedparticles.insert(particlerowmap->GID(lid));

  // copy to a vector and create extended particle colmap
  std::vector<int> myparticlecolgids(myescapedparticles.begin(),myescapedparticles.end());
  Teuchos::RCP<Epetra_Map> particlecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)myparticlecolgids.size(),&myparticlecolgids[0],0,bindis_->Comm()));

  // now ghost the nodes
  bindis_->ExportColumnNodes(*particlecolmap);

  // call fillcomplete
  bindis_->FillComplete();

  return;
}

/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 06/14 |
 *-------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
  const Epetra_Map*              initial_elecolmap,
  std::map<int, std::set<int> >& binelemap,
  std::map<int, std::set<int> >& ext_bintoele_ghosting,
  Teuchos::RCP<Epetra_Map>       rowbins,
  Teuchos::RCP<Epetra_Map>       bincolmap)
{
  // do communication to gather all elements for extended ghosting
  const int numproc = initial_elecolmap->Comm().NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // gather set of column bins for each proc
    std::set<int> bins;
    if(iproc == myrank_)
    {
      // either use given column layout of bins ...
      if(bincolmap != Teuchos::null)
      {
        int nummyeles = bincolmap->NumMyElements();
        int* entries = bincolmap->MyGlobalElements();
        bins.insert(entries, entries+nummyeles);
      }
      else // ... or add an extra layer to the given bin distribution
      {
        for(std::map<int, std::set<int> >::const_iterator iter=binelemap.begin(); iter!=binelemap.end(); ++iter)
        {
          int binId = iter->first;
          // avoid getting two layer ghosting as this is not needed
          if(rowbins!=Teuchos::null)
          {
            const int lid = rowbins->LID(binId);
            if(lid<0)
              continue;
          }
          std::vector<int> binvec;
          // get neighboring bins
          GetBinConnectivity(binId, binvec);
          bins.insert(binvec.begin(), binvec.end());
          // insert bin itself
          bins.insert(binId);
        }
      }
    }
    // copy set to vector in order to broadcast data
    std::vector<int> binids(bins.begin(),bins.end());

    // first: proc i tells all procs how many bins it has
    int numbin = binids.size();
    initial_elecolmap->Comm().Broadcast(&numbin, 1, iproc);
    // second: proc i tells all procs which bins it has
    binids.resize(numbin);

    initial_elecolmap->Comm().Broadcast(&binids[0], numbin, iproc);

    // loop over all own bins and find requested ones, fill in elements in these bins
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
      ext_bintoele_ghosting = rdata;
    }
  }

  // reduce map of sets to one set and copy to a vector to create extended elecolmap
  std::set<int> coleleset;
  std::map<int, std::set<int> >::iterator iter;
  for(iter=ext_bintoele_ghosting.begin(); iter!= ext_bintoele_ghosting.end(); ++iter)
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

/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 11/16 |
 *-------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ExtendGhosting(
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<Epetra_Map> initial_elecolmap,
    Teuchos::RCP<Epetra_Map> bincolmap,
    bool assigndegreesoffreedom,
    bool initelements,
    bool doboundaryconditions,
    bool checkghosting)
{
  std::map<int, std::set<int> > rowelesinbin;
  DistributeElesToBins(dis, rowelesinbin);

  // get extended column map elements
  std::map<int, std::set<int> > dummy;
   Teuchos::RCP<Epetra_Map> elecolmapextended =
       ExtendGhosting(&*initial_elecolmap, rowelesinbin, dummy, Teuchos::null, bincolmap);

  // extend ghosting (add nodes/elements) according to the new column layout
  dis->ExtendedGhosting(*elecolmapextended, assigndegreesoffreedom, initelements, doboundaryconditions, checkghosting);

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
  // create XAABB for discretization
  if(dis != Teuchos::null)
    CreateXAABB(dis);

  // divide global bounding box into bins
  for (int dim=0; dim<3; ++dim)
  {
    // determine number of bins per direction for prescribed cutoff radius
    // std::floor leads to bins that are at least of size cutoff_radius
    if (cutoff_radius_>0.0)
      bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));
    else
      dserror("cutoff_radius <= zero");

    // for detailed description of the difference between bin_per_dir
    // and id_calc_bin_per_dir_ see BinningStrategy::ConvertGidToijk;
    int n=0;
    do
    {
      id_calc_bin_per_dir_[dim] = std::pow(2,n);
      id_calc_exp_bin_per_dir_[dim] = n;
      ++n;
    } while (id_calc_bin_per_dir_[dim] < bin_per_dir_[dim]);

    // calculate size of bins in each direction
    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0)) / bin_per_dir_[dim];
    // calculate inverse of size of bins in each direction
    inv_bin_size_[dim] = 1.0/bin_size_[dim];
  }

  if(id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1] * id_calc_bin_per_dir_[2] > std::numeric_limits<int>::max())
    dserror("number of bins is larger than an integer can hold! Reduce number of bins by increasing the cutoff radius");

  // determine cutoff radius for prescribed number of bins per direction
  // fixme: this cannot be reached
  if(cutoff_radius_ <= 0.)
    cutoff_radius_ = std::min(bin_size_[0],std::min(bin_size_[1],bin_size_[2]));

  // 2D case
  if(particle_dim_ != INPAR::PARTICLE::particle_3D)
    CreateBins2D();

  return;
}

/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateBins2D()
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
  id_calc_bin_per_dir_[entry] = 1;
  id_calc_exp_bin_per_dir_[entry] = 0;
  bin_size_[entry] = (XAABB_(entry,1)-XAABB_(entry,0));
  inv_bin_size_[entry] = 1.0/bin_size_[entry];

  return;
}

/*----------------------------------------------------------------------*
| find XAABB and cutoff                                     ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateXAABB(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<Epetra_Vector> disnp,
  LINALG::Matrix<3,2>& XAABB,
  bool setcutoff
  )
{
  // cutoff as largest element in discret on each proc
  double locmaxcutoff = 0.0;
  double currpos[3] = {0.0,0.0,0.0};
  // initialize XAABB of discret as rectangle around the first node of
  // discret on each proc
  GetCurrentNodePos(discret,discret->lRowNode(0),disnp,currpos);
  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim, 0) = currpos[dim] - GEO::TOL7;
    XAABB(dim, 1) = currpos[dim] + GEO::TOL7;
  }

  // loop over row elements of each proc
  for(int i=0; i<discret->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = discret->lRowElement(i);

    // eleXAABB for each row element
    LINALG::Matrix<3,2> eleXAABB(false);

    // initialize eleXAABB as rectangle around the first node of ele
    GetCurrentNodePos(discret,ele->Nodes()[0],disnp,currpos);
    for(int dim=0; dim<3; ++dim)
    {
     eleXAABB(dim, 0) = currpos[dim] - GEO::TOL7;
     eleXAABB(dim, 1) = currpos[dim] + GEO::TOL7;
    }

    // loop over remaining nodes of current rowele
    for (int lid = 1; lid < ele->NumNode(); ++lid)
    {
      const DRT::Node* node = ele->Nodes()[lid];
      GetCurrentNodePos(discret,node,disnp,currpos);

      //  merge eleXAABB of all nodes of this element
      for(int dim=0; dim < 3; dim++)
      {
        eleXAABB(dim, 0) = std::min( eleXAABB(dim, 0), currpos[dim] - GEO::TOL7);
        eleXAABB(dim, 1) = std::max( eleXAABB(dim, 1), currpos[dim] + GEO::TOL7);
      }
    }

    // compute cutoff as largest element in discret
    if(setcutoff)
      for(int dim=0; dim<3; ++dim)
        locmaxcutoff = std::max( locmaxcutoff, eleXAABB(dim, 1) - eleXAABB(dim, 0));

     // merge XAABB of all roweles
     for(int dim=0; dim < 3; dim++)
     {
       XAABB(dim, 0) = std::min( XAABB(dim, 0), eleXAABB(dim, 0));
       XAABB(dim, 1) = std::max( XAABB(dim, 1), eleXAABB(dim, 1));
     }
  }

  // local bounding box on each proc
  double locmin[3] = {XAABB(0,0), XAABB(1,0), XAABB(2,0)};
  double locmax[3] = {XAABB(0,1), XAABB(1,1), XAABB(2,1)};
  // global bounding box over all procs
  double globmin[3];
  double globmax[3];
  // do the necessary communication
  discret->Comm().MinAll(&locmin[0], &globmin[0], 3);
  discret->Comm().MaxAll(&locmax[0], &globmax[0], 3);

  // set global XAABB for discret
  for(int dim=0; dim<3; ++dim)
  {
    XAABB(dim,0) = globmin[dim];
    XAABB(dim,1) = globmax[dim];
  }

  // maxall of cutoff
  if(setcutoff)
  {
    double globmaxcutoff = 0.0;
    discret->Comm().MaxAll(&locmaxcutoff, &globmaxcutoff, 1);
    // this is necessary if more than one discret is relevant
    cutoff_radius_ = std::max(globmaxcutoff, cutoff_radius_);
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------*
| compute max cutoff as largest element in discret                      |
| in current configuration                              eichinger 09/16 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ComputeMaxCutoff(
  std::vector<Teuchos::RCP<DRT::Discretization> > discret,
  std::vector<Teuchos::RCP<Epetra_Vector> >disnp
  )
{
  // reset cutoff_radius
  cutoff_radius_= 0.0;

  // loop over all input discrets
  for(size_t ndis=0; ndis<discret.size(); ++ndis)
  {
    // cutoff as largest element in discret
    double locmaxcutoff = 0.0;
    double currpos[3] = {0.0,0.0,0.0};

    // loop over row elements of each proc
    for(int i=0; i<discret[ndis]->NumMyRowElements(); ++i)
    {
      DRT::Element* ele = discret[ndis]->lRowElement(i);

      // eleXAABB for each row element
      LINALG::Matrix<3,2> eleXAABB(false);

      // initialize eleXAABB as rectangle around the first node of ele
      GetCurrentNodePos(discret[ndis],ele->Nodes()[0],disnp[ndis],currpos);
      for(int dim=0; dim<3; ++dim)
      {
       eleXAABB(dim, 0) = currpos[dim] - GEO::TOL7;
       eleXAABB(dim, 1) = currpos[dim] + GEO::TOL7;
      }

      // loop over remaining nodes of current rowele
      for (int lid = 1; lid < ele->NumNode(); ++lid)
      {
        const DRT::Node* node = ele->Nodes()[lid];
        GetCurrentNodePos(discret[ndis],node,disnp[ndis],currpos);

        //  merge eleXAABB of all nodes of this element
        for(int dim=0; dim < 3; dim++)
        {
          eleXAABB(dim, 0) = std::min( eleXAABB(dim, 0), currpos[dim] - GEO::TOL7);
          eleXAABB(dim, 1) = std::max( eleXAABB(dim, 1), currpos[dim] + GEO::TOL7);
        }
      }

      // compute cutoff as largest element in discret
      for(int dim=0; dim<3; ++dim)
        locmaxcutoff = std::max( locmaxcutoff, eleXAABB(dim, 1) - eleXAABB(dim, 0));
    }

    double globmaxcutoff = 0.0;
    discret[ndis]->Comm().MaxAll(&locmaxcutoff, &globmaxcutoff, 1);
    // this is necessary if more than one discret is relevant
    cutoff_radius_ = std::max(globmaxcutoff, cutoff_radius_);
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------*
| find XAABB and cutoff                                     ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateXAABB(
  std::vector<Teuchos::RCP<DRT::Discretization> > discret,
  std::vector<Teuchos::RCP<Epetra_Vector> >disnp,
  bool setcutoff
  )
{
  // reset cutoff
  if(setcutoff)
    cutoff_radius_ = 0.0;

  // initialize XAABB_ as rectangle around the first node of first discret
  const DRT::Node* node = discret[0]->lRowNode(0);
  // calculate current position of this node
  double currpos[3] = {0.0,0.0,0.0};
  GetCurrentNodePos(discret[0],node,disnp[0],currpos);

  for(int dim=0; dim<3; ++dim)
  {
    XAABB_(dim, 0) = currpos[dim] - GEO::TOL7;
    XAABB_(dim, 1) = currpos[dim] + GEO::TOL7;
  }

  // build XAABB_ from XAABB of all discrets and determine maximal element extension
  // to use as new cutoff
  for(size_t i=0; i<discret.size(); ++i)
  {
    LINALG::Matrix<3,2> XAABB;
    CreateXAABB(discret[i],disnp[0],XAABB, setcutoff);

    // set XAABB_ considering all input discrets
    for(int dim=0; dim < 3; dim++)
    {
      XAABB_(dim, 0) = std::min( XAABB_(dim, 0), XAABB(dim, 0) );
      XAABB_(dim, 1) = std::max( XAABB_(dim, 1), XAABB(dim, 1) );
    }
  }

  // enlarge cutoff a little bit for safety reasons
  if(setcutoff) cutoff_radius_ += GEO::TOL7;

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
 | build periodic boundary conditions                       ghamm 04/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::BuildPeriodicBC()
{
  // build periodic boundary condition
  std::vector<DRT::Condition*> conds;
  bindis_->GetCondition("ParticlePeriodic", conds);
  if(conds.size() > 1)
    dserror("only one periodic boundary condition allowed for particles");

  // leave when no pbc available
  if(conds.size() == 0)
    return;
  else
    havepbc_ = true;

  // now read in the available condition
  const std::vector<int>* onoff = conds[0]->Get<std::vector<int> >("ONOFF");

  // loop over all spatial directions
  for(int dim=0; dim<3; ++dim)
  {
    if((*onoff)[dim])
    {
      // output pbc bounds based on XAABB of bins
      if(myrank_ == 0)
        std::cout << "INFO: PBC bounds for particles is computed automatically for direction " << dim
                  << " based on XAABB of bins (left: " <<  XAABB_(dim,0) << " , right: " <<  XAABB_(dim,1) << " )" << std::endl;

      // additional safety check whether at least some bin layers exist in pbc direction
      // --> facilitates neighbor search --> see contact
      if(bin_per_dir_[dim] < 3)
        dserror("There are just very few bins in pbc direction -> maybe nasty for neighborhood search (especially in contact)");

      // set flag
      pbconoff_[dim] = true;

      // offset delta for pbc direction
      pbcdeltas_[dim] = XAABB_(dim,1) - XAABB_(dim,0);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::BuildPeriodicBC(Teuchos::RCP<std::vector<double> > periodlength)
{

  if(periodlength->size()<3)
    dserror("Size of periodlength of periodic bounding box < 3 for 3D Problem");

  // loop over all spatial directions
  for(int dim=0; dim<3; ++dim)
  {
    if(periodlength->at(dim) > 0.0)
    {
      // output pbc bounds based on XAABB of bins
      if(myrank_ == 0)
        std::cout << "INFO: PBC bounds are computed automatically for direction " << dim
                  << " based on XAABB of bins (left: " <<  XAABB_(dim,0) << " , right: " <<  XAABB_(dim,1) << " )" << std::endl;

      // to use the following, bins have to be created before so that bin_per_dir is set. This function
      // is intended to be called before bins are created. Therefore the following has to be ensured somewhere else
//      // additional safety check whether at least some bin layers exist in pbc direction
//      // --> facilitates neighbor search --> see contact
//      if(bin_per_dir_[dim] < 3)
//        dserror("There are just very few bins in pbc direction -> maybe nasty for neighborhood search (especially in contact)");

      // set flag
      pbconoff_[dim] = true;

      // at least one direction has periodic boundary conditions
      havepbc_ = true;

      // offset delta for pbc direction
      pbcdeltas_[dim] = periodlength->at(dim);
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
    ijk[dim] = (int)(std::floor((pos[dim]-XAABB_(dim,0)) * inv_bin_size_[dim]));

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 01/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid(const double* pos)
{
  int ijk[3];
  for(int dim=0; dim<3; ++dim)
    ijk[dim] = (int)(std::floor((pos[dim]-XAABB_(dim,0)) * inv_bin_size_[dim]));

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 02/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertPosToijk(const double* pos, int* ijk)
{
  for(int dim=0; dim<3; ++dim)
    ijk[dim] = (int)(std::floor((pos[dim]-XAABB_(dim,0)) * inv_bin_size_[dim]));

  return;
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid(const LINALG::Matrix<3,1>& pos)
{
  int ijk[3];
  for(int dim=0; dim<3; ++dim)
    ijk[dim] = (int)(std::floor((pos(dim)-XAABB_(dim,0)) * inv_bin_size_[dim]));

  return ConvertijkToGid(&ijk[0]);
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertPosToijk(const LINALG::Matrix<3,1>& pos, int* ijk)
{
  for(int dim=0; dim<3; ++dim)
    ijk[dim] = (int)(std::floor((pos(dim)-XAABB_(dim,0)) * inv_bin_size_[dim]));

  return;
}


/*----------------------------------------------------------------------*
| convert i,j,k into bin id                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertijkToGid(int* ijk)
{
  // might need to modify ijk connectivity in the presence of periodic boundary conditions
  if(havepbc_)
  {
    for(unsigned idim=0; idim<3; ++idim)
    {
      if(pbconoff_[idim])
      {
        if(ijk[idim] == -1)
          ijk[idim] = bin_per_dir_[idim] - 1;
        else if(ijk[idim] == bin_per_dir_[idim])
          ijk[idim] = 0;
      }
    }
  }

  // given ijk is outside of XAABB
  if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 || ijk[0]>=bin_per_dir_[0] || ijk[1]>=bin_per_dir_[1] || ijk[2]>=bin_per_dir_[2] )
    return -1;

  return ijk[0] + ijk[1]*id_calc_bin_per_dir_[0] + ijk[2]*id_calc_bin_per_dir_[0]*id_calc_bin_per_dir_[1];
}


/*----------------------------------------------------------------------*
| convert bin id into i,j,k                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertGidToijk(const int gid, int* ijk)
{
  // in order to efficiently compute the ijk triple from a given bin id,
  // use of the shift operator is made (right shift by one equals division by two)
  // therefore it is necessary that the number of bins per direction is
  // divisible by 2
  // (shift operation costs one cycle vs division (or modulo) costs 20--40 cycles on cpu)
  // Hence, two different number of bins per direction are needed
  // one for the used domain and the other one for converting gid <-> ijk

  // example: 2^n = bin_per_dir
  // example: gid >> n = (int)gid/bin_per_dir

  ijk[2] = gid >> (id_calc_exp_bin_per_dir_[0] + id_calc_exp_bin_per_dir_[1]);

  const int tmp = gid - ijk[2]*id_calc_bin_per_dir_[0]*id_calc_bin_per_dir_[1];

  ijk[1] = tmp >> id_calc_exp_bin_per_dir_[0];

  ijk[0] = tmp - ijk[1]*id_calc_bin_per_dir_[0];

  // alternative method - more expensive but only based on integer operations:
//  {
//    const int tmp1 = gid % (id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1]);
//    // compute i
//    ijk[0] = tmp1 % id_calc_bin_per_dir_[0];
//    // compute j
//    ijk[1] = (tmp1 - ijk[0]) / id_calc_bin_per_dir_[0];
//    // compute k
//    ijk[2] = (gid - ijk[0] - ijk[1]*id_calc_bin_per_dir_[0]) / (id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1]);
//  }

  // found ijk is outside of XAABB
  if( ijk[0]<0 || ijk[1]<0 || ijk[2]<0 || ijk[0]>=bin_per_dir_[0] || ijk[1]>=bin_per_dir_[1] || ijk[2]>=bin_per_dir_[2] )
    dserror("ijk (%d %d %d) for given gid: %d is outside of range (bin per dir: %d %d %d)",
        ijk[0], ijk[1], ijk[2], gid, bin_per_dir_[0], bin_per_dir_[1], bin_per_dir_[2]);

  return;
}


/*----------------------------------------------------------------------*
 | get all bins in ijk range                               ghamm 02/13  |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GidsInijkRange(const int* ijk_range, std::set<int>& binIds, bool checkexistence)
{
  if(checkexistence == true and bindis_ == Teuchos::null)
    dserror("particle discretization is not set up correctly");

  for(int i=ijk_range[0]; i<=ijk_range[1]; ++i)
  {
    for(int j=ijk_range[2]; j<=ijk_range[3]; ++j)
    {
      for(int k=ijk_range[4]; k<=ijk_range[5]; ++k)
      {
        int ijk[3] = {i,j,k};

        const int gid = ConvertijkToGid(&ijk[0]);
        if(gid != -1)
        {
          if(checkexistence)
          {
            if(bindis_->HaveGlobalElement(gid))
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
void BINSTRATEGY::BinningStrategy::GidsInijkRange(const int* ijk_range, std::vector<int>& binIds, bool checkexistence)
{
  if(checkexistence == true and bindis_ == Teuchos::null)
    dserror("particle discretization is not set up correctly");

  for(int i=ijk_range[0]; i<=ijk_range[1]; ++i)
  {
    for(int j=ijk_range[2]; j<=ijk_range[3]; ++j)
    {
      for(int k=ijk_range[4]; k<=ijk_range[5]; ++k)
      {
        int ijk[3] = {i,j,k};

        const int gid = ConvertijkToGid(&ijk[0]);
        if(gid != -1)
        {
          if(checkexistence)
          {
            if(bindis_->HaveGlobalElement(gid))
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
void BINSTRATEGY::BinningStrategy::GetBinConnectivity(const int binId, std::vector<int>& binIds)
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
        const int gid = ConvertijkToGid(&ijk[0]);
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
void BINSTRATEGY::BinningStrategy::GetBinCorners(const int binId, std::vector<LINALG::Matrix<3,1> >& bincorners)
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
        const int ijk_curr[] = {i,j,k};
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
LINALG::Matrix<3,1> BINSTRATEGY::BinningStrategy::GetBinCentroid(const int binId)
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetCurrentNodePos(
  Teuchos::RCP<DRT::Discretization> discret,
  const DRT::Node* node,
  Teuchos::RCP<Epetra_Vector> disnp,
  double* currpos
  )
{
  if(disnp!=Teuchos::null)
  {
    const int gid = discret->Dof(node, 0);
    const int lid = disnp->Map().LID(gid);
    if(lid<0)
      dserror("Your displacement is incomplete (need to be based on a column map"
              " as this function is also called from a loop over elements and "
              "each proc does (usually) not own all nodes of his row elements ");
    for(int dim=0; dim<3; ++dim)
      currpos[dim] = node->X()[dim] + (*disnp)[lid+dim];
  }
  else
  {
    for(int dim=0; dim<3; ++dim)
      currpos[dim] = node->X()[dim];
  }

  return;
}
