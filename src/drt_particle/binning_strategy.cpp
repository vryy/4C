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
#include "binning_strategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_node.H"


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
    if(static_cast<MORTAR::MortarElement*>(ele)->IsSlave() == isslave)
    {
      DRT::Node** nodes = ele->Nodes();
      const int numnode = ele->NumNode();

      // initialize ijk_range with ijk of first node of element
      int ijk[3];
      {
        DRT::Node* node = nodes[0];
        const double* coords = static_cast<MORTAR::MortarNode*>(node)->xspatial();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; ++j)
      {
        DRT::Node* node = nodes[j];
        const double* coords = static_cast<MORTAR::MortarNode*>(node)->xspatial();
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
      std::set<int> binIds;
      GidsInijkRange(&ijk_range[0], binIds, false);

      // assign element to bins
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
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
  std::map<int, std::set<int> >& elesinbin
  )
{
  // exploit bounding box idea for elements in underlying discretization and bins
  {
    // loop over row elements is enough, extended ghosting always includes standard ghosting
    for (int lid = 0; lid < underlyingdis->NumMyRowElements(); ++lid)
    {
      DRT::Element* fluidele = underlyingdis->lRowElement(lid);
      DRT::Node** fluidnodes = fluidele->Nodes();
      const int numnode = fluidele->NumNode();

      // initialize ijk_range with ijk of first node of fluid element
      int ijk[3];
      {
        const DRT::Node* node = fluidnodes[0];
        const double* coords = node->X();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; ++j)
      {
        const DRT::Node* node = fluidnodes[j];
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
      std::set<int> binIds;
      GidsInijkRange(&ijk_range[0], binIds, false);

      // assign fluid element to bins
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        elesinbin[*biniter].insert(fluidele->Id());
    }
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


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateBins(Teuchos::RCP<DRT::Discretization> dis)
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

  // divide global bounding box into bins
  for (int dim=0; dim<3; ++dim)
  {
    // determine number of bins per direction for prescribed cutoff radius
    // std::floor leads to bins that are at least of size cutoff_radius
    if (cutoff_radius_>0.0)
      bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));

    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

//  if(myrank_ == 0)
//  {
//    std::cout << "Global bounding box size: " << XAABB_
//        << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << std::endl;
//  }

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
    ijk[dim] = (int)((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
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
    ijk[dim] = (int)((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
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
    ijk[dim] = (int)((pos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
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
    ijk[dim] = (int)((pos(dim)-XAABB_(dim,0)) / bin_size_[dim]);
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
    ijk[dim] = (int)((pos(dim)-XAABB_(dim,0)) / bin_size_[dim]);
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
            binIds.insert(gid);
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
