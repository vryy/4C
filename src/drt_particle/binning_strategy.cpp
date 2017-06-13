/*----------------------------------------------------------------------*/
/*!
\file binning_strategy.cpp

\brief Binning strategy for neighborhood search

\level 2

\maintainer Jonas Eichinger
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
#include <Teuchos_TimeMonitor.hpp>

#include "binning_strategy.H"
#include "binning_strategy_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_inpar/inpar_meshfree.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_dofset_independent.H"

#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_node.H"

#include "../drt_particle/particle_algorithm.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_beam3/beam3_base.H"


/*----------------------------------------------------------------------*
 | standard constructor                                                 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy() :
    bindis_(Teuchos::null),
    visbindis_(Teuchos::null),
    cutoff_radius_(0.0),
    XAABB_(true),
    deforming_simulation_domain_handler(Teuchos::null),
    writebinstype_(DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(DRT::Problem::Instance()->MeshfreeParams(),("WRITEBINS"))),
    havepbc_(false),
    particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(DRT::Problem::Instance()->ParticleParams(),"DIMENSION"))
{
  // initialize arrays
  for( int idim = 0; idim < 3; ++idim )
  {
    bin_size_[idim] = 0.0;
    inv_bin_size_[idim] = 0.0;
    bin_per_dir_[idim] = 0;
    pbconoff_[idim] = false;
    pbcdeltas_[idim] = 0.0;
  }
  boundaryrowbins_.clear();
  boundarycolbins_.clear();
}

/*----------------------------------------------------------------------------*
 | Init                                                       eichinger 11/16 |
 *----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::Init(
    Teuchos::RCP<DRT::Discretization>& bindis,
    Teuchos::RCP<DRT::Discretization> const discret,
    Teuchos::RCP<Epetra_Vector> const disnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox>const pbb)
{
  // myrank
  myrank_ = bindis->Comm().MyPID();
  // binning discretization
  bindis_ = bindis;
  //
  deforming_simulation_domain_handler = pbb;
  // meshfree params
  const Teuchos::ParameterList& meshfreeparams = DRT::Problem::Instance()->MeshfreeParams();

  // get type of bounding box specification
  INPAR::MESHFREE::xaabbspectype xaabbpectype =
      DRT::INPUT::IntegralValue<INPAR::MESHFREE::xaabbspectype>(meshfreeparams,"DEFINEXAABBPER");

  switch (xaabbpectype)
  {
    case INPAR::MESHFREE::input:
    {
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
          else
            dserror("specify six values for bounding box in three dimensional problem. Fix input file");
        }
      }

      break;
    }
    case INPAR::MESHFREE::dynamic:
    {
      CreateXAABB(discret, disnp, XAABB_);
      break;
    }
    default :
    {
      dserror("You should not be here");
      break;
    }
  }

  // get type for bin specification
  INPAR::MESHFREE::binspectype binspectype =
      DRT::INPUT::IntegralValue<INPAR::MESHFREE::binspectype>(meshfreeparams,"DEFINEBINSPER");

  switch (binspectype)
  {
    case INPAR::MESHFREE::cutoff:
    {
      // get cutoff radius
      cutoff_radius_ = meshfreeparams.get<double>("CUTOFF_RADIUS");
      if(cutoff_radius_<0.0)
        dserror("Negative cutoff radius set in input file for definition of bins. Fix it ...");

      // some check
      std::istringstream binstream(Teuchos::getNumericStringParameter(meshfreeparams,"BIN_PER_DIR"));
      for(int idim=0; idim<3; idim++)
      {
        int val = -1;
        if (binstream >> val)
        {
          if(val > 0 && myrank_ == 0)
            std::cout<<"\n WARNING: specified number of bins per direction not used "
                       " as you choose DEFINEBINSPER cutoff"<<std::endl;
        }
      }
      break;
    }
    case INPAR::MESHFREE::binsperdir:
    {
      // get number of bins per direction
      std::istringstream binstream(Teuchos::getNumericStringParameter(meshfreeparams,"BIN_PER_DIR"));
      for(int idim=0; idim<3; idim++)
      {
        int val = -1;
        if (binstream >> val)
        {
          if(val>0)
            bin_per_dir_[idim] = val;
          else
            dserror("Negative number of bins in direction %i does not make sense", idim);
        }
        else
        {
          dserror("You need to specify three figures for BIN_PER_DIR in input file for three dimensional problem. ");
        }
      }

      // some check
      if(cutoff_radius_ > 0.0)
        std::cout<<"\n WARNING: specified cutoff radius not used "
                   " as you choose DEFINEBINSPER binsperdir"<<std::endl;

      break;
    }
    case INPAR::MESHFREE::largestele:
    {
      // todo:
      dserror("Biopolynet: unshifted configuration is needed (not yet here) for calculation of cutoff.");
      // store structure discretization in vector
      std::vector<Teuchos::RCP<DRT::Discretization> > discret_vec(1);
      discret_vec[0] = discret;
      // displacement vector according to periodic boundary conditions
      std::vector<Teuchos::RCP<Epetra_Vector> > disnp_vec(1);
      disnp_vec[0] = disnp;
      cutoff_radius_ = ComputeMinCutoffAsMaxEdgeLengthOfXAABBOfLargestEle( discret_vec, disnp_vec );

      break;
    }
    default :
    {
      dserror("You should not be here");
      break;
    }
  }
  // done
  return;
}

/*----------------------------------------------------------------------------*
 | Setup                                                      eichinger 11/13 |
 *----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::Setup()
{
  // create bins
  CreateBinsBasedOnCutoffAndXAABB();

  // build periodic boundary condition
  BuildPeriodicBC();

  // done
  return;
}

/*----------------------------------------------------------------------*
 | Binning strategy constructor                             ghamm 11/13 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy(
    const Epetra_Comm& comm,
    double cutoff_radius,
    LINALG::Matrix<3,2> XAABB
    ) :
    bindis_(Teuchos::null),
    visbindis_(Teuchos::null),
    cutoff_radius_(cutoff_radius),
    XAABB_(XAABB),
    deforming_simulation_domain_handler(Teuchos::null),
    writebinstype_(DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(DRT::Problem::Instance()->MeshfreeParams(),("WRITEBINS"))),
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
  CreateBinsBasedOnCutoffAndXAABB();

  return;
}


/*----------------------------------------------------------------------*
 | Binning strategy constructor                             ghamm 11/13 |
 *----------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy::BinningStrategy(
    const Epetra_Comm& comm
    ) :
    bindis_(Teuchos::null),
    visbindis_(Teuchos::null),
    cutoff_radius_(0.0),
    XAABB_(true),
    deforming_simulation_domain_handler(Teuchos::null),
    writebinstype_(DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(DRT::Problem::Instance()->MeshfreeParams(),("WRITEBINS"))),
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
    XAABB_(true),
    deforming_simulation_domain_handler(Teuchos::null),
    writebinstype_(DRT::INPUT::IntegralValue<INPAR::MESHFREE::compltype>(DRT::Problem::Instance()->MeshfreeParams(),("WRITEBINS"))),
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
}

/*----------------------------------------------------------------------*
| fill bins into bin discretization                         ghamm 08/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::FillBinsIntoBinDiscretization(
    Teuchos::RCP<Epetra_Map> const& rowbins)
{
  // fill bins into bindis_
  for( int i = 0; i < rowbins->NumMyElements(); ++i )
  {
    const int gid = rowbins->GID(i);
    Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
    bindis_->AddElement(bin);
  }
}

/*----------------------------------------------------------------------*
| assign elements into bins                                 ghamm 11/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeElesToBins(
    const DRT::Discretization& mortardis,
    std::map<int, std::set<int> >& binelemap,
    bool isslave) const
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
      binIds.reserve( GetNumberOfBinsInijkRange(ijk_range) );
      GidsInijkRange(&ijk_range[0], binIds, false);

      // assign element to bins
      for(std::vector<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        binelemap[*biniter].insert(ele->Id());
    }
  }

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeRowElementsToBinsUsingEleXAABB(
    Teuchos::RCP< DRT::Discretization > const& discret,
    std::map< int, std::set< int > >&   bintorowelemap,
    Teuchos::RCP< Epetra_Vector >       disnp
    ) const
{
  bintorowelemap.clear();

  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for ( int lid = 0; lid < discret->NumMyRowElements(); ++lid )
  {
    DRT::Element* eleptr = discret->lRowElement(lid);
    // get corresponding bin ids in ijk range
    std::vector<int> binIds;
    DistributeElementToBinsUsingEleXAABB( discret, eleptr, binIds, disnp );

   // assign element to bins
    std::vector< int >::const_iterator biniter;
    for( biniter = binIds.begin(); biniter != binIds.end(); ++biniter )
      bintorowelemap[*biniter].insert( eleptr->Id() );
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeElementToBinsUsingEleXAABB(
    Teuchos::RCP< DRT::Discretization > const& discret,
    DRT::Element* eleptr,
    std::vector<int>& binIds,
    Teuchos::RCP< Epetra_Vector > const& disnp
    ) const
{
  binIds.clear();
  DRT::Node** nodes = eleptr->Nodes();

  // initialize ijk_range with ijk of first node of element
  int ijk[3];
  DRT::Node const * const node = nodes[0];
  DistributeNodeToijk( discret, node, disnp, ijk );

  // ijk_range contains: i_min i_max j_min j_max k_min k_max
  int ijk_range[] = { ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2] };

  // bounding box idea for rigid sphere element with just one node needs
  // some special treatment
  if ( eleptr->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
  {
    BuildAxisAlignedijkRangeForRigidSphere( discret, eleptr, disnp, ijk, ijk_range );
  }
  else if ( dynamic_cast<DRT::ELEMENTS::Beam3Base*>(eleptr) != NULL )
  {
    // fill in remaining nodes
    for ( int j = 1; j < eleptr->NumNode(); ++j )
    {
      DRT::Node const * const node = nodes[j];
      DistributeNodeToijk( discret, node, disnp, ijk );
      AddijkToAxisAlignedijkRangeOfBeamElement( ijk, ijk_range );
    }
  }
  else
  {
    // fill in remaining nodes
    for ( int j = 1; j < eleptr->NumNode(); ++j )
    {
      DRT::Node const * const node = nodes[j];
      DistributeNodeToijk( discret, node, disnp, ijk );
      AddijkToAxisAlignedijkRangeOfElement( ijk, ijk_range );
    }
  }

  // get corresponding bin ids in ijk range
  binIds.reserve( GetNumberOfBinsInijkRange(ijk_range) );
  GidsInijkRange( &ijk_range[0], binIds, false );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeNodeToijk(
    Teuchos::RCP< DRT::Discretization> const& discret,
    DRT::Node const * const node,
    Teuchos::RCP< Epetra_Vector > const& disnp,
    int ijk[3]
  ) const
{
  double currpos[3] = { 0.0, 0.0, 0.0 };
  GetCurrentNodePos( discret, node, disnp, currpos );
  double const* coords = currpos;
  ConvertPosToijk( coords, ijk );
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::GetNumberOfBinsInijkRange( int const ijk_range[6] ) const
{
  return ( ( ijk_range[1] - ijk_range[0] + 1 ) * ( ijk_range[3] - ijk_range[2] + 1 )
       * ( ijk_range[5] - ijk_range[4] + 1) ) ;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::AddijkToAxisAlignedijkRangeOfElement(
    int const ijk[3],
    int ijk_range[6]
  ) const
{
  for( int dim = 0; dim < 3; ++dim )
  {
    if( ijk[dim] < ijk_range[dim * 2] )
    {
      ijk_range[dim * 2] = ijk[dim];
    }

    if( ijk[dim] > ijk_range[dim * 2 + 1] )
    {
      ijk_range[dim * 2 + 1] = ijk[dim];
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::AddijkToAxisAlignedijkRangeOfBeamElement(
    int const ijk[3],
    int ijk_range[6]
  ) const
{
  // this should be large enough
  int const cutcheckfac = 5;

  for( int dim = 0; dim < 3; ++dim )
  {
    if( ijk[dim] < ijk_range[dim * 2] )
    {
      // this is needed if your element is cut by a periodic boundary and you don't
      // want a bounding box that is as big as the the whole binning domain
      if( ( ijk[dim] == 0 ) && ( abs(ijk[dim] - ijk_range[dim * 2] ) > cutcheckfac ) )
      {
        ijk_range[dim * 2 + 1] = bin_per_dir_[dim];
        continue;
      }
      else
      {
        ijk_range[dim * 2] = ijk[dim];
      }
    }
    if( ijk[dim] > ijk_range[dim * 2 + 1] )
    {
      // cut check
      if( ( ijk[dim] == bin_per_dir_[dim] - 1) && ( abs(ijk[dim] - ijk_range[dim * 2 + 1]) > cutcheckfac ) )
      {
        ijk_range[dim * 2] = -1;
      }
      else
      {
        ijk_range[dim * 2 + 1] = ijk[dim];
      }
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::BuildAxisAlignedijkRangeForRigidSphere(
    Teuchos::RCP< DRT::Discretization> const& discret,
    DRT::Element const * const eleptr,
    Teuchos::RCP< Epetra_Vector > const& disnp,
    int ijk[3],
    int ijk_range[6]
  ) const
{
  double const& radius = dynamic_cast< const DRT::ELEMENTS::Rigidsphere* >(eleptr)->Radius();
  double currpos[3] = { 0.0, 0.0, 0.0 };
  GetCurrentNodePos( discret, eleptr->Nodes()[0], disnp, currpos );

  for ( int j = 0; j < 3; ++j )
  {
    double* coords = currpos;
    coords[j] += radius;
    ConvertPosToijk( coords, ijk );
    AddijkToAxisAlignedijkRangeOfElement( ijk, ijk_range );
    coords[j] -= (2.0 * radius);
    ConvertPosToijk( coords, ijk );
    AddijkToAxisAlignedijkRangeOfElement( ijk, ijk_range );
    coords[j] += radius;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::AssignElesToBins(
    Teuchos::RCP<DRT::Discretization> discret,
    std::map<int, std::set<int> >     extendedfieldghosting
  ) const
{
  // loop over bins
  std::map<int, std::set<int> >::const_iterator biniter;
  for( biniter = extendedfieldghosting.begin(); biniter != extendedfieldghosting.end(); ++biniter )
  {
    // extract bins from discretization after checking on existence
    const int lid = bindis_->ElementColMap()->LID( biniter->first );
    if( lid < 0 )
      continue;

    // get current bin
    DRT::MESHFREE::MeshfreeMultiBin* currbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>( bindis_->gElement( biniter->first ) );

    // loop over ele content of this bin
    std::set<int>::const_iterator eleiter;
    for( eleiter = biniter->second.begin(); eleiter != biniter->second.end(); ++eleiter )
    {
      // add eleid and elepointer to current bin
      currbin->AddAssociatedEle( BINSTRATEGY::UTILS::ConvertElementToBinContentType( discret->gElement( *eleiter ) ),
          *eleiter, discret->gElement( *eleiter ) );
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetBinContent(
    std::set<DRT::Element*> &eles,
    std::vector< BINSTRATEGY::UTILS::BinContentType > bincontent,
    std::vector<int>& binIds,
    bool roweles
)
{
  // loop over all bins
  std::vector<int>::const_iterator biniter;
  for( biniter = binIds.begin(); biniter != binIds.end(); ++biniter )
  {
    // extract bins from discretization after checking on existence
    const int lid = bindis_->ElementColMap()->LID(*biniter);
    if( lid < 0 )
      continue;

    // get content of current bin
    GetBinContent( bindis_->lColElement(lid), eles, bincontent, roweles );
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetBinContent(
    DRT::Element* binptr,
    std::set<DRT::Element*> &eles,
    std::vector< BINSTRATEGY::UTILS::BinContentType > bincontent,
    bool roweles
)
{
  DRT::MESHFREE::MeshfreeMultiBin* bin =
      static_cast<DRT::MESHFREE::MeshfreeMultiBin*>( binptr );

#ifdef DEBUG
  // safety check
  if( bin == NULL )
    dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif

  // loop over bincontent you want to get
  for( int bc_i = 0; bc_i < static_cast<int>( bincontent.size() ); ++bc_i )
  {
    // gather elements of with specific bincontent type
    DRT::Element** elements = bin->AssociatedEles(bincontent[bc_i]);
    const int numeles = bin->NumAssociatedEle(bincontent[bc_i]);
    for( int iele = 0;iele < numeles; ++iele )
    {
      if( roweles && elements[iele]->Owner() != myrank_ )
        continue;
      eles.insert(elements[iele]);
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::RemoveSpecificElesFromBins(
    BINSTRATEGY::UTILS::BinContentType bincontent )
{
  // loop over all bins and remove assigned elements
  const int numcolbins = bindis_->NumMyColElements();
  for( int binlid = 0; binlid < numcolbins; ++binlid )
  {
    DRT::Element *currentbin = bindis_->lColElement(binlid);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currentbin)->RemoveSpecificAssociatedEles( bincontent );
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::RemoveAllElesFromBins()
{
  // loop over all bins and remove assigned elements
  const int numcolbins = bindis_->NumMyColElements();
  for( int binlid = 0; binlid < numcolbins; ++binlid )
  {
    DRT::Element *currentbin = bindis_->lColElement(binlid);
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currentbin)->RemoveAllAssociatedEles();
  }
}

/*----------------------------------------------------------------------*
| assign nodes into bins                                    ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DistributeNodesToBins(
    Teuchos::RCP<DRT::Discretization> discret,
    std::map<int, std::vector<int> >& nodesinbin,
    Teuchos::RCP<Epetra_Vector> disnp
    ) const
{
  // current position of nodes
  double currpos[3] = {0.0,0.0,0.0};

  // loop over row nodes
  for ( int lid = 0; lid < discret->NumMyRowNodes(); ++lid )
  {
    DRT::Node* node = discret->lRowNode(lid);
    GetCurrentNodePos( discret, node, disnp, currpos );

    const double* coords = currpos;
    int ijk[3];
    ConvertPosToijk( coords, ijk );
    const int binid = ConvertijkToGid(&ijk[0]);

    if ( binid == -1 )
      dserror( "There are nodes in your discretization that reside outside the binning \n"
               "domain, this does not work at this point.");

    // assign node to bin
    nodesinbin[binid].push_back( node->Id() );
  }

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::WeightedPartitioning(
    std::vector<Teuchos::RCP<DRT::Discretization> > discret,
    std::vector<Teuchos::RCP<Epetra_Map> >&         stdelecolmap,
    std::vector<Teuchos::RCP<Epetra_Map> >&         stdnodecolmap
    )
{
  // initialize dummys
  std::vector<std::map<int, std::set<int> > > dummy1(discret.size());
  std::vector<Teuchos::RCP<Epetra_Vector> >   dummy2(discret.size());

  ComputeMinXAABBContainingAllElementsOfInputDiscrets( discret, dummy2, XAABB_, true );

  // create bins
  CreateBinsBasedOnCutoffAndXAABB(Teuchos::null);

  // ------------------------------------------------------------------------
  // create bins, weight them according to number of nodes (of discrets) they
  // contain, account for bin connectivity. Then an optimal distribution of
  // bins to procs can be obtained
  // ------------------------------------------------------------------------
  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int> > > nodesinbin(discret.size());
  // default weight 10.0
  double const weight = 10.0;
  Teuchos::RCP<Epetra_Map> newrowbins =
      WeightedDistributionOfBinsToProcs( discret, dummy2, nodesinbin, weight );

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
  for( size_t i = 0; i < discret.size(); ++i )
  {
   // ----------------------------------------------------------------------
   // start with standard ghosting
   // ----------------------------------------------------------------------
    StandardDiscretizationGhosting( discret[i], newrowbins, dummy2[i],
        stdelecolmap[i], stdnodecolmap[i] );

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
    DistributeRowElementsToBinsUsingEleXAABB( discret[i], bintoelemap, dummy2[i] );

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
    std::vector< Teuchos::RCP<DRT::Discretization> >&  discret,
    std::vector< Teuchos::RCP<Epetra_Vector> >&        disnp,
    std::vector< std::map<int, std::vector<int> > > &  nodesinbin,
    double const& weight,
    bool repartition
    ) const
{
  // calculate total number of bins
  const int numbin = bin_per_dir_[0] * bin_per_dir_[1] * bin_per_dir_[2];

  // some safety checks to ensure efficiency
  {
    if( numbin < discret[0]->Comm().NumProc() && myrank_ == 0 )
      dserror("ERROR:NumProc > NumBin. Too many processors to "
              "distribute your bins properly!!!");

    if( numbin < 8 * discret[0]->Comm().NumProc() && myrank_ == 0 )
      std::cout << "\n\nWARNING: partitioning not useful, choose less procs. "
                   " Owner distribution may be inefficient!\n\n" << std::endl;
  }

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
    for (int lid=0; lid<oldrowmap->NumMyElements(); ++lid)
    {
      const int binId = oldrowmap->GID(lid);

      std::vector<int> neighbors;
      GetNeighborBinIds(binId,neighbors);

      int err = bingraph->InsertGlobalIndices(binId,(int)neighbors.size(),&neighbors[0]);
      if (err<0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d for global row %d",err,binId);
    }
  }
  else
  {
    // dummy row bin distribution (equally distributed over all procs as no
    // weighting done so far)
    rowbins = CreateLinearMapForNumbin(discret[0]->Comm());
    // create nodal graph
    bingraph = Teuchos::rcp( new Epetra_CrsGraph( Copy, *rowbins, 108, false) );
  }

  // Now we're going to create a Epetra_Vector with vertex/node weights to be
  // used for the partitioning operation (weights must be at least one for zoltan)
  Teuchos::RCP<Epetra_Vector> vweights = LINALG::CreateVector(*rowbins, true);

  // set weights of bins related to the number of nodes of discrets that are contained
  // empty bins have weight of 1
  vweights->PutScalar(1.0);

  // determine which node is in which bin and weight each bin according to
  // "weight" times the number of nodes it contains
  // assign all node gids to their corresponding bin
  for( int i = 0; i < static_cast<int>( discret.size() ); ++i )
  {
    // distribute nodes, that are owned by a proc, to the bins of this proc
    DistributeNodesToBins( discret[i], nodesinbin[i], disnp[i] );

    std::map<int, std::vector<int> > nodesinmybins;
    // gather information of bin content from other procs (bin is owned by this
    // proc and there are some nodes on other procs which are located in this bin)
    // mynodesinbin then contains all node gids (vector) that reside in a owned bin (gid is map key)
    CollectInformation( rowbins, nodesinbin[i], nodesinmybins );

    // weight each bin with 10 times the number of node it contains
    // empty bins remain with weight one
    std::map< int, std::vector<int> >::const_iterator biniter;
    for( biniter = nodesinmybins.begin(); biniter != nodesinmybins.end(); ++biniter)
    {
      int lid = rowbins->LID(biniter->first);
      // safety check
      if (lid <0 )
        dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",discret[i]->Comm().MyPID(),biniter->first);

      // weighting
      (*vweights)[lid] += weight * static_cast<double> ( biniter->second.size() );
    }
  }

   // fill bin connectivity into bin graph
   for ( int lid = 0; lid < rowbins->NumMyElements(); ++lid )
   {
     int rowbinid = rowbins->GID(lid);
     // insert 26 (one level) neighboring bins to graph
     // (if active, periodic boundary conditions are considered here)
     std::vector<int> neighbors;
     GetNeighborBinIds( rowbinid, neighbors );

     int err = bingraph->InsertGlobalIndices( rowbinid, static_cast<int>( neighbors.size()), &neighbors[0] );
     if ( err < 0 )
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

/*-----------------------------------------------------------------------------*
| determine boundary row bins                                  eichinger 01/17 |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DetermineBoundaryRowBins()
{
  // clear old content
  boundaryrowbins_.clear();

  // fill discret if necessary to obtain bin row map
  if( bindis_->Filled() == false )
    bindis_->FillComplete( false, false, false );

  // determine maximal possible number of neighbors
  size_t nummaxneighbors = -1;
  if( particle_dim_ == INPAR::PARTICLE::particle_3D )
    nummaxneighbors = 26;
  else
    nummaxneighbors = 8;

  // loop over row bins and decide whether they are located at the boundary
  const Epetra_Map* binrowmap = bindis_->ElementRowMap();
  for ( int lid = 0; lid < binrowmap->NumMyElements(); ++lid )
  {
    DRT::Element* currbin = bindis_->lRowElement(lid);
    std::vector<int> binvec;
    binvec.reserve( nummaxneighbors );
    // get neighboring bins
    GetNeighborBinIds( currbin->Id(), binvec );

    // a bin with less than 26 (or 8 in 2D) neighbors is a boundary bin
    if( binvec.size() < nummaxneighbors )
    {
      boundaryrowbins_.push_back( currbin );
      continue;
    }

    // a bin with less than 26 (or 8 in 2D) row neighbors is a boundary bin between processors
    std::vector<int> rowbinvec;
    rowbinvec.reserve( nummaxneighbors );
    for( std::vector<int>::const_iterator it = binvec.begin(); it != binvec.end(); ++it )
    {
      if( binrowmap->LID( *it ) != -1 )
        rowbinvec.push_back(*it);
    }

    if( rowbinvec.size() < nummaxneighbors )
      boundaryrowbins_.push_back( currbin );
  }
}

/*-----------------------------------------------------------------------------*
 | determine one layer ghosting around boundary row bins      eichinger 01/17  |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::DetermineBoundaryColBinsIds()
{
  boundarycolbins_.clear();

  if( boundaryrowbins_.size() == 0 )
    DetermineBoundaryRowBins();

  // loop over boundary row bins and add neighbors
  std::list<DRT::Element*>::const_iterator it;
  for( it = boundaryrowbins_.begin(); it != boundaryrowbins_.end(); ++it )
  {
      std::vector<int> binvec;
      binvec.reserve(26);
      // get neighboring bins
      GetNeighborBinIds( (*it)->Id(), binvec );
      boundarycolbins_.insert( binvec.begin(), binvec.end() );
  }
}

/*--------------------------------------------------------------------------*
| standard ghosting according to bin distribution               ghamm 06/14 |
 *--------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::StandardDiscretizationGhosting(
    Teuchos::RCP<DRT::Discretization> & discret,
    Teuchos::RCP<Epetra_Map> const&     rowbins,
    Teuchos::RCP<Epetra_Vector>&        disnp,
    Teuchos::RCP<Epetra_Map>&           stdelecolmap,
    Teuchos::RCP<Epetra_Map>&           stdnodecolmap ) const
{
  // each owner of a bin gets owner of the nodes this bin contains
  // all other nodes of elements, of which proc is owner of at least one
  // node, are ghosted
  Teuchos::RCP<Epetra_CrsGraph> initgraph = discret->BuildNodeGraph();

  // Todo introduced this export to column map due to special handling of
  //      beam nodes without own position DoFs in GetCurrentNodePos()
  Teuchos::RCP<Epetra_Vector> disnp_col = Teuchos::null;
  if ( discret->HaveDofs() and disnp != Teuchos::null)
  {
    disnp_col = Teuchos::rcp(new Epetra_Vector( *discret->DofColMap() ) );
    LINALG::Export( *disnp, *disnp_col );
  }

  // distribute nodes, that are owned by a proc, to the bins of this proc
  std::map<int, std::vector<int> > nodesinbin;
  DistributeNodesToBins( discret, nodesinbin, disnp_col );

  std::map<int, std::vector<int> > nodesinmybins;
  // gather information of bin content from other procs (bin is owned by
  // this proc and there are some nodes on other procs which are located
  // in this bin here)
  CollectInformation( rowbins, nodesinbin, nodesinmybins );

  // build new node row map
  std::vector<int> mynewrownodes;
  std::map<int, std::vector<int> >::const_iterator biniter;
  for( biniter = nodesinmybins.begin(); biniter != nodesinmybins.end(); ++biniter )
  {
    std::vector<int>::const_iterator nodeiter;
    for(nodeiter=biniter->second.begin(); nodeiter!=biniter->second.end(); ++nodeiter)
    {
      mynewrownodes.push_back(*nodeiter);
    }
  }
  nodesinmybins.clear();

  Teuchos::RCP<Epetra_Map> newnoderowmap =
      Teuchos::rcp(new Epetra_Map( -1, mynewrownodes.size(), &mynewrownodes[0], 0, discret->Comm() ) );

  // create the new graph and export to it
  Teuchos::RCP<Epetra_CrsGraph> newnodegraph;

  newnodegraph = Teuchos::rcp(new Epetra_CrsGraph( Copy, *newnoderowmap, 108, false ) );
  Epetra_Export exporter( initgraph->RowMap(), *newnoderowmap );
  int err = newnodegraph->Export( *initgraph, exporter, Add );
  if ( err < 0 )
    dserror("Graph export returned err=%d",err);
  newnodegraph->FillComplete();
  newnodegraph->OptimizeStorage();

  // the column map will become the new ghosted distribution of nodes (standard ghosting)
  const Epetra_BlockMap cntmp = newnodegraph->ColMap();
  stdnodecolmap =
      Teuchos::rcp(new Epetra_Map(-1,cntmp.NumMyElements(),cntmp.MyGlobalElements(),0,discret->Comm()));

  // rebuild of the discretizations with new maps for standard ghosting
  Teuchos::RCP<Epetra_Map> roweles;
  discret->BuildElementRowColumn( *newnoderowmap, *stdnodecolmap, roweles, stdelecolmap );
  discret->ExportRowNodes( *newnoderowmap );
  discret->ExportRowElements( *roweles );
  discret->ExportColumnNodes( *stdnodecolmap );
  discret->ExportColumnElements( *stdelecolmap );
  // in case we have a state vector, we need to build the dof map to enable its rebuild
  if( disnp == Teuchos::null )
  {
    discret->FillComplete(false,false,false);
  }
  else
  {
    discret->FillComplete(true,false,false);
    Teuchos::RCP<Epetra_Vector> old;
    old = disnp;
    disnp = LINALG::CreateVector( *discret->DofRowMap(), true );
    LINALG::Export( *old, *disnp );
  }

#ifdef DEBUG
  // print distribution after standard ghosting
  // some output after standard ghosting
  if( myrank_ == 0 )
    std::cout << "parallel distribution with standard ghosting" << std::endl;
  DRT::UTILS::PrintParallelDistribution(*discret);
#endif

  return;
}

/*-----------------------------------------------------------------------------*
| extend ghosting according to bin distribution                eichinger 09/16 |
*-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ExtendBinGhosting(
    Teuchos::RCP<Epetra_Map> rowbins,
    std::set< int > const&   colbins,
    bool assigndegreesoffreedom)
{
  // gather set of bins that need to be ghosted by myrank
  std::set<int> bins( colbins.begin(), colbins.end() );

  // insert row bins
  for( int i = 0; i < rowbins->NumMyElements(); ++i )
    bins.insert( rowbins->GID(i) );

  std::vector<int> bincolmapvec( bins.begin(), bins.end() );
  Teuchos::RCP<Epetra_Map> bincolmap = Teuchos::rcp( new Epetra_Map(
      -1, static_cast<int>( bincolmapvec.size() ), &bincolmapvec[0], 0, bindis_->Comm() ) );

  if( bincolmap->NumGlobalElements() == 1 && bindis_->Comm().NumProc() > 1 )
    dserror("one bin cannot be run in parallel -> reduce CUTOFF_RADIUS");

  BINSTRATEGY::UTILS::ExtendDiscretizationGhosting( bindis_ , bincolmap , assigndegreesoffreedom , false , true );

#ifdef DEBUG
  // check whether each proc has only particles that are within bins on this proc
  for ( int k = 0; k<bindis_->NumMyColElements(); ++k )
  {
    int binid = bindis_->lColElement(k)->Id();
    DRT::Node** particles = bindis_->lColElement(k)->Nodes();

    for ( int iparticle = 0; iparticle < bindis_->lColElement(k)->NumNode(); ++iparticle )
    {
      double const* pos = particles[iparticle]->X();
      int ijk[3] = { -1, -1, -1 };
      ConvertPosToijk( pos, ijk );

      int gidofbin = ConvertijkToGid( &ijk[0] );
      if ( gidofbin != binid )
        dserror("after ghosting: particle which should be in bin no. %i is in %i", gidofbin, binid );
    }
  }
#endif
}

/*-------------------------------------------------------------------*
| extend ghosting according to bin distribution          ghamm 11/13 |
 *-------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
    DRT::Discretization&           mortardis,
    Teuchos::RCP<Epetra_Map>       initial_elecolmap,
    std::map<int, std::set<int> >& slavebinelemap,
    std::map<int, std::set<int> >& masterbinelemap) const
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
        GetNeighborAndOwnBinIds(binId, bins);
        binset.insert(bins.begin(), bins.end());
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
    std::map<int, std::set<int> >& myescapedpartelemap) const
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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
    std::map<int, std::set<int> >& binelemap,
    std::map<int, std::set<int> >& ext_bintoele_ghosting,
    Teuchos::RCP<Epetra_Map>       bincolmap) const
{
  // do communication to gather all elements for extended ghosting
  const int numproc = bindis_->Comm().NumProc();
  for ( int iproc = 0; iproc < numproc; ++iproc )
  {
    // gather set of column bins for each proc
    std::set<int> bins;
    if( iproc == myrank_ )
    {
      int nummyeles = bincolmap->NumMyElements();
      int* entries = bincolmap->MyGlobalElements();
      bins.insert(entries, entries+nummyeles);
    }

    // copy set to vector in order to broadcast data
    std::vector<int> binids( bins.begin(), bins.end() );

    // first: proc i tells all procs how many bins it has
    int numbin = binids.size();
    bindis_->Comm().Broadcast( &numbin, 1, iproc );
    // second: proc i tells all procs which bins it has
    binids.resize(numbin);

    bindis_->Comm().Broadcast( &binids[0], numbin, iproc );

    // loop over all own bins and find requested ones, fill in elements in these bins
    std::map<int, std::set<int> > sdata;
    std::map<int, std::set<int> > rdata;

    for( int i = 0; i < numbin; ++i )
    {
      if( binelemap.find(binids[i]) != binelemap.end() )
        sdata[binids[i]].insert( binelemap[binids[i]].begin(), binelemap[binids[i]].end() );
    }

    LINALG::Gather<int>(sdata, rdata, 1, &iproc, bindis_->Comm());

    // proc i has to store the received data
    if(iproc == myrank_)
      ext_bintoele_ghosting = rdata;
  }

  // reduce map of sets to one set and copy to a vector to create extended elecolmap
  std::set<int> coleleset;
  std::map<int, std::set<int> >::iterator iter;
  for( iter = ext_bintoele_ghosting.begin(); iter != ext_bintoele_ghosting.end(); ++iter )
    coleleset.insert(iter->second.begin(),iter->second.end());

  std::vector<int> colgids(coleleset.begin(),coleleset.end());

  // return extended elecolmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)colgids.size(),&colgids[0],0,bindis_->Comm()));
}

/*----------------------------------------------------------------------------*
| extend ghosting according to bin distribution                   ghamm 06/14 |
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::ExtendGhosting(
    const Epetra_Map*              initial_elecolmap,
    std::map<int, std::set<int> >& binelemap,
    std::map<int, std::set<int> >& ext_bintoele_ghosting,
    Teuchos::RCP<Epetra_Map>       rowbins,
    Teuchos::RCP<Epetra_Map>       bincolmap) const
{
  // do communication to gather all elements for extended ghosting
  const int numproc = initial_elecolmap->Comm().NumProc();
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    // gather set of column bins for each proc
    std::set<int> bins;
    if( iproc == myrank_ )
    {
      // either use given column layout of bins ...
      if( bincolmap != Teuchos::null )
      {
        int nummyeles = bincolmap->NumMyElements();
        int* entries = bincolmap->MyGlobalElements();
        bins.insert(entries, entries+nummyeles);
      }
      else // ... or add an extra layer to the given bin distribution
      {
        std::map<int, std::set<int> >::const_iterator iter;
        for( iter = binelemap.begin(); iter != binelemap.end(); ++iter )
        {
          int binId = iter->first;
          // avoid getting two layer ghosting as this is not needed
          if( rowbins != Teuchos::null )
          {
            const int lid = rowbins->LID(binId);
            if( lid < 0 )
              continue;
          }
          std::vector<int> binvec;
          // get neighboring bins
          GetNeighborAndOwnBinIds( binId, binvec );
          bins.insert( binvec.begin(), binvec.end() );
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

    for( int i = 0; i < numbin; ++i )
    {
      if( binelemap.find(binids[i]) != binelemap.end() )
        sdata[binids[i]].insert( binelemap[binids[i]].begin(), binelemap[binids[i]].end() );
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
    std::vector<Teuchos::RCP<DRT::Discretization> > dis) const
{
  for(size_t i=0; i<dis.size(); ++i)
  {
    //----------------------------
    // start with extended ghosting
    //----------------------------
    // fill elements into bins
    std::map<int, std::set<int> > binelemap;
    DistributeRowElementsToBinsUsingEleXAABB(dis[i], binelemap);

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
 | extend ghosting according to bin distribution         ghamm 11/16 |
 *-------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ExtendEleGhosting(
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<Epetra_Map> initial_elecolmap,
    Teuchos::RCP<Epetra_Map> bincolmap,
    bool assigndegreesoffreedom,
    bool initelements,
    bool doboundaryconditions) const
{
  std::map<int, std::set<int> > rowelesinbin;
  DistributeRowElementsToBinsUsingEleXAABB(dis, rowelesinbin);

  // get extended column map elements
  std::map<int, std::set<int> > dummy;
   Teuchos::RCP<Epetra_Map> elecolmapextended =
       ExtendGhosting(&*initial_elecolmap, rowelesinbin, dummy, Teuchos::null, bincolmap);

  // extend ghosting (add nodes/elements) according to the new column layout
   BINSTRATEGY::UTILS::ExtendDiscretizationGhosting(dis, elecolmapextended,
       assigndegreesoffreedom, initelements, doboundaryconditions);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::RevertExtendedGhosting(
    std::vector<Teuchos::RCP<DRT::Discretization> > dis,
    std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
    std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap
    ) const
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
    std::map<int, std::vector<int> >& nodesinmybins) const
{
  // do communication to gather all nodes
  const int numproc = rowbins->Comm().NumProc();
  for ( int iproc = 0; iproc < numproc; ++iproc )
  {
    // vector with row bins on this proc
    std::vector<int> binids;
    int numbin;
    if( iproc == myrank_ )
    {
      int* myrowbinsdata = rowbins->MyGlobalElements();
      numbin = rowbins->NumMyElements();
      binids.insert( binids.begin(), myrowbinsdata, myrowbinsdata + numbin );
    }

    // first: proc i tells all procs how many bins it has
    rowbins->Comm().Broadcast( &numbin, 1, iproc );
    binids.resize(numbin);
    // second: proc i tells all procs which bins it has, now each proc contains
    // rowbingids of iproc in vector binids
    rowbins->Comm().Broadcast( &binids[0], numbin, iproc );

    // loop over all own bins and find requested ones, fill in master elements in these bins
    // (map key is bin gid owned by iproc, vector contains all node gids of all procs in this bin)
    std::map< int, std::vector<int> > sdata;
    std::map< int, std::vector<int> > rdata;

    for( int i = 0; i < numbin; ++i )
    {
      // now each procs checks if row nodes lie in bins of iproc ...
      if( nodesinbin.find( binids[i] ) != nodesinbin.end() )
        // ... if so, each proc assignes its node gids to iprocs bins
        sdata[binids[i]].insert(sdata[binids[i]].begin(), nodesinbin[binids[i]].begin(),
            nodesinbin[binids[i]].end() );
    }

    // iprocs gathers all this information from other procs
    LINALG::Gather<int>(sdata, rdata, 1, &iproc, rowbins->Comm());

    // iproc has to store the received data
    if( iproc == myrank_ )
    {
      // clear data and refill
      nodesinmybins.clear();
      nodesinmybins.insert( rdata.begin(), rdata.end() );
    }
  }

  return;
}

/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 09/12 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateBinsBasedOnCutoffAndXAABB(Teuchos::RCP<DRT::Discretization> dis)
{
  // create XAABB for discretization
  if(dis != Teuchos::null)
    CreateXAABB(dis);

  // divide global bounding box into bins
  for (int dim=0; dim<3; ++dim)
  {
    // determine number of bins per direction for prescribed cutoff radius
    // std::floor leads to bins that are at least of size cutoff_radius
    if (cutoff_radius_ > 0.0)
    {
      bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));
    }

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

  // determine cutoff radius if number of bins per dir was prescribed
  // was prescribed in input file
  if(cutoff_radius_ <= 0.0 )
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
| find XAABB and divide into bins                       rasthofer 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateBinsScatra(Teuchos::RCP<DRT::Discretization> dis)
{
  double tolbox = 1.0e-7;

  // for level-set problems, we always determine the bounding box based on the discretization
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

  // ensure that bounding box is indeed equal to the domain or slightly smaller but never larger!!!
  for(int dim=0; dim<3; ++dim)
  {
    XAABB_(dim,0) = globmin[dim] + tolbox;
    XAABB_(dim,1) = globmax[dim] - tolbox;
   }

  // divide global bounding box into bins
  for (int dim=0; dim<3; ++dim)
  {
    // determine number of bins per direction for prescribed cutoff radius
    // std::floor leads to bins that are at least of size cutoff_radius
    if (cutoff_radius_>0.0)
      dserror("Do not use the cutoff radius for level-set problems!");
      //bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));

    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
    inv_bin_size_[dim] = 1.0/bin_size_[dim];

    // for detailed description of the difference between bin_per_dir
    // and id_calc_bin_per_dir_ see BinningStrategy::ConvertGidToijk;
    int n=0;
    do
    {
      id_calc_bin_per_dir_[dim] = std::pow(2,n);
      id_calc_exp_bin_per_dir_[dim] = n;
      ++n;
    } while (id_calc_bin_per_dir_[dim] < bin_per_dir_[dim]);
  }

  if(id_calc_bin_per_dir_[0] * id_calc_bin_per_dir_[1] * id_calc_bin_per_dir_[2] > std::numeric_limits<int>::max())
    dserror("number of bins is larger than an integer can hold! Reduce number of bins by increasing the cutoff radius");

  if(myrank_ == 0)
  {
    std::cout << "Global bounding box size: " << XAABB_
        << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << std::endl;
  }

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
double BINSTRATEGY::BinningStrategy::ComputeMinCutoffAsMaxEdgeLengthOfXAABBOfLargestEle(
    std::vector<Teuchos::RCP<DRT::Discretization> > discret,
    std::vector<Teuchos::RCP<Epetra_Vector> > disnp
    )
{
  double cutoff_radius = 0.0;

  // loop over all input discrets
  for( size_t ndis = 0; ndis < discret.size(); ++ndis )
  {
    // cutoff as largest element in discret
    double locmaxcutoff = 0.0;
    double currpos[3] = {0.0,0.0,0.0};

    // loop over row elements of each proc
    for( int i = 0; i < discret[ndis]->NumMyRowElements(); ++i )
    {
      DRT::Element* ele = discret[ndis]->lRowElement(i);

      // eleXAABB for each row element
      LINALG::Matrix<3,2> eleXAABB(false);

      // initialize eleXAABB as rectangle around the first node of ele
      GetCurrentNodePos( discret[ndis], ele->Nodes()[0], disnp[ndis], currpos );
      for( int dim = 0; dim < 3; ++dim )
      {
       eleXAABB(dim, 0) = currpos[dim] - GEO::TOL7;
       eleXAABB(dim, 1) = currpos[dim] + GEO::TOL7;
      }

      // rigid sphere elements needs to consider its radius
      if( ele->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
      {
        double radius = dynamic_cast<DRT::ELEMENTS::Rigidsphere*>(ele)->Radius();
        for( int dim = 0; dim < 3; ++dim )
        {
          eleXAABB(dim, 0) = std::min( eleXAABB(dim, 0), eleXAABB(dim, 0) - radius - GEO::TOL7 );
          eleXAABB(dim, 1) = std::max( eleXAABB(dim, 1), eleXAABB(dim, 0) + radius + GEO::TOL7 );
        }
      }
      else
      {
        // loop over remaining nodes of current rowele
        for ( int lid = 1; lid < ele->NumNode(); ++lid )
        {
          const DRT::Node* node = ele->Nodes()[lid];
          GetCurrentNodePos( discret[ndis], node, disnp[ndis], currpos );

          //  merge eleXAABB of all nodes of this element
          for( int dim = 0; dim < 3; ++dim )
          {
            eleXAABB(dim, 0) = std::min( eleXAABB(dim, 0), currpos[dim] - GEO::TOL7 );
            eleXAABB(dim, 1) = std::max( eleXAABB(dim, 1), currpos[dim] + GEO::TOL7 );
          }
        }
      }
      // compute cutoff as largest element in discret
      for( int dim = 0; dim < 3; ++dim )
        locmaxcutoff = std::max( locmaxcutoff, eleXAABB(dim, 1) - eleXAABB(dim, 0) );
    }

    double globmaxcutoff = 0.0;
    discret[ndis]->Comm().MaxAll( &locmaxcutoff, &globmaxcutoff, 1);
    // this is necessary if more than one discret is relevant
    cutoff_radius = std::max( globmaxcutoff, cutoff_radius );
  }

  return cutoff_radius;
}

/*----------------------------------------------------------------------*
| find XAABB and cutoff                                     ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ComputeMinXAABBContainingAllElementsOfInputDiscrets(
    std::vector<Teuchos::RCP<DRT::Discretization> > discret,
    std::vector<Teuchos::RCP<Epetra_Vector> >disnp,
    LINALG::Matrix<3,2>& XAABB,
    bool setmincutoff)
{
  // reset cutoff
  if(setmincutoff)
    cutoff_radius_ = 0.0;

  // initialize XAABB_ as rectangle around the first node of first discret
  const DRT::Node* node = discret[0]->lRowNode(0);
  // calculate current position of this node
  double currpos[3] = { 0.0, 0.0, 0.0 };
  GetCurrentNodePos( discret[0], node, disnp[0], currpos );

  for( int dim = 0; dim < 3 ; ++dim )
  {
    XAABB(dim, 0) = currpos[dim] - GEO::TOL7;
    XAABB(dim, 1) = currpos[dim] + GEO::TOL7;
  }

  // build XAABB_ from XAABB of all discrets and determine maximal element extension
  // to use as new cutoff
  for( size_t i = 0; i < discret.size(); ++i )
  {
    LINALG::Matrix<3,2> locXAABB;
    CreateXAABB( discret[i], disnp[0], locXAABB, setmincutoff );

    // set XAABB_ considering all input discrets
    for( int dim = 0; dim < 3; ++dim )
    {
      XAABB(dim, 0) = std::min( XAABB(dim, 0), locXAABB(dim, 0) );
      XAABB(dim, 1) = std::max( XAABB(dim, 1), locXAABB(dim, 1) );
    }
  }

  // enlarge cutoff a little bit for safety reasons
  if(setmincutoff) cutoff_radius_ += GEO::TOL7;
}

/*----------------------------------------------------------------------*
| find XAABB                                                ghamm 06/14 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CreateXAABB(
    Teuchos::RCP<DRT::Discretization> dis )
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
  std::istringstream periodicbc(Teuchos::getNumericStringParameter(
      DRT::Problem::Instance()->MeshfreeParams(),"PERIODICONOFF"));

  // loop over all spatial directions
  for(int dim=0; dim<3; ++dim)
  {
    int val = -1;
    if (periodicbc >> val)
    {
      if( val )
      {
        // output pbc bounds based on XAABB of bins
        if(myrank_ == 0)
          std::cout << "INFO: PBC bounds for particles is computed automatically for direction " << dim
                    << " based on XAABB of bins (left: " <<  XAABB_(dim,0) << " , right: " <<  XAABB_(dim,1) << " )" << std::endl;

        // set flag
        pbconoff_[dim] = true;

        // offset delta for pbc direction
        pbcdeltas_[dim] = XAABB_(dim,1) - XAABB_(dim,0);

        // set global flag
        havepbc_ = true;
      }
    }
    else
    {
      dserror("Enter three values to specify each direction as periodic or non periodic. Fix input file ...");
    }
  }

  return;
}

/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 01/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid( const double* pos ) const
{
  int ijk[3];
  double pos_ud[3];
  if ( deforming_simulation_domain_handler != Teuchos::null )
  {
    deforming_simulation_domain_handler->TransformFromGlobalToUndeformedBoundingBoxSystem( pos, pos_ud );
    for ( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos_ud[dim] - XAABB_( dim, 0 ) ) * inv_bin_size_[dim] ) );
  }
  else
  {
    for ( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos[dim] - XAABB_( dim, 0 ) ) * inv_bin_size_[dim] ) );
  }

  return ConvertijkToGid( &ijk[0] );
}


/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 02/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertPosToijk(const double* pos, int* ijk) const
{
  double pos_ud[3];
  if ( deforming_simulation_domain_handler != Teuchos::null )
  {
    deforming_simulation_domain_handler->TransformFromGlobalToUndeformedBoundingBoxSystem( pos, pos_ud );
    for ( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos_ud[dim] - XAABB_( dim, 0 ) ) * inv_bin_size_[dim] ) );
  }
  else
  {
    for ( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos[dim] - XAABB_( dim, 0 ) ) * inv_bin_size_[dim] ) );
  }

  return;
}

/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::ConvertPosToijk( const LINALG::Matrix<3,1>& pos, int* ijk ) const
{
  LINALG::Matrix<3,1> pos_ud;
  if ( deforming_simulation_domain_handler != Teuchos::null )
  {
    deforming_simulation_domain_handler->TransformFromGlobalToUndeformedBoundingBoxSystem( pos, pos_ud );
    for( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos_ud(dim) - XAABB_( dim, 0) ) * inv_bin_size_[dim] ) );
  }
  else
  {
    for( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos(dim) - XAABB_( dim, 0) ) * inv_bin_size_[dim] ) );
  }

  return;
}

/*----------------------------------------------------------------------*
| convert position first to i,j,k, then into bin id         ghamm 03/13 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertPosToGid( const LINALG::Matrix<3,1>& pos ) const
{
  int ijk[3];
  LINALG::Matrix<3,1> pos_ud;
  if ( deforming_simulation_domain_handler != Teuchos::null )
  {
    deforming_simulation_domain_handler->TransformFromGlobalToUndeformedBoundingBoxSystem( pos, pos_ud );
    for ( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos_ud(dim) - XAABB_( dim, 0 ) ) * inv_bin_size_[dim] ) );
  }
  else
  {
    for ( int dim = 0; dim < 3; ++dim )
      ijk[dim] = static_cast<int>( std::floor( ( pos(dim) - XAABB_( dim, 0 ) ) * inv_bin_size_[dim] ) );
  }

  return ConvertijkToGid( &ijk[0] );
}

/*----------------------------------------------------------------------*
| convert i,j,k into bin id                                 ghamm 09/12 |
 *----------------------------------------------------------------------*/
int BINSTRATEGY::BinningStrategy::ConvertijkToGid(int* ijk) const
{
  // might need to modify ijk connectivity in the presence of periodic boundary conditions
  if( havepbc_ )
  {
    for( unsigned idim = 0; idim < 3; ++idim )
    {
      if( pbconoff_[idim] )
      {
        if( ijk[idim] == -1 )
          ijk[idim] = bin_per_dir_[idim] - 1;
        else if( ijk[idim] == bin_per_dir_[idim] )
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
void BINSTRATEGY::BinningStrategy::ConvertGidToijk(const int gid, int* ijk) const
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
void BINSTRATEGY::BinningStrategy::GidsInijkRange(
    const int* ijk_range, std::set<int>& binIds, bool checkexistence) const
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
void BINSTRATEGY::BinningStrategy::GidsInijkRange(
    const int* ijk_range, std::vector<int>& binIds, bool checkexistence) const
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
void BINSTRATEGY::BinningStrategy::GetNeighborBinIds(
    const int binId, std::vector<int>& binIds) const
{
  int ijk_base[3];
  ConvertGidToijk( binId, &ijk_base[0] );

  for( int i = ijk_base[0] - 1; i <= ijk_base[0] + 1; ++i )
  {
    for( int j = ijk_base[1] - 1; j <= ijk_base[1] + 1 ; ++j )
    {
      for( int k = ijk_base[2] - 1; k <= ijk_base[2] + 1; ++k )
      {
        int ijk[3] = { i, j, k };
        const int gid = ConvertijkToGid( &ijk[0] );
        if( gid != -1 and gid != binId )
        {
          binIds.push_back(gid);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | get 26 neighboring bin ids and myself                   ghamm 08/13  |
*-----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetNeighborAndOwnBinIds(
    const int binId, std::vector<int>& binIds) const
{
  // get neighbors
  GetNeighborBinIds(binId, binIds);

  // add myself
  binIds.push_back(binId);

   //in case of less than two bins in pbc direction, this is needed
   //to avoid double contact evaluation
  if( havepbc_ )
  {
    std::sort( binIds.begin(), binIds.end() );
    binIds.erase( std::unique( binIds.begin(), binIds.end() ), binIds.end() );
  }

  return;
}

/*----------------------------------------------------------------------*
| corner position for given bin id                          ghamm 03/13 |
 *----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetBinCorners(
    const int binId, std::vector<LINALG::Matrix<3,1> >& bincorners) const
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

  // change entries to get node numbering according to baci convention
  std::swap( bincorners[2], bincorners[3] );
  std::swap( bincorners[6], bincorners[7] );

  return;
}


/*----------------------------------------------------------------------*
| centroid position for given bin id                        ghamm 04/13 |
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,1> BINSTRATEGY::BinningStrategy::GetBinCentroid(const int binId) const
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

/*-----------------------------------------------------------------------------*
| write bin output                                             eichinger 11/16 |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::WriteBinOutput(int const step, double const time)
{
  // no bin output
  if( writebinstype_ == INPAR::MESHFREE::none )
    return;

  // -------------------------------------------------------------------------
  // note: this is a debug feature only (as very expensive)
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // create new discretization containing the bins as elements
  // -------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp( bindis_->Comm().Clone() );
  visbindis_ = Teuchos::rcp(new DRT::Discretization( "bins" , com ) );
  // create discretization writer
  visbindis_->SetWriter( Teuchos::rcp( new IO::DiscretizationWriter( visbindis_ ) ) );

  // store gids of ghosted elements
  std::map<int, std::vector<LINALG::Matrix<3,1> > > ghostcorners;
  // add elements and nodes
  for ( int i = 0; i < bindis_->NumMyColElements(); ++i )
  {
    DRT::Element* ele = bindis_->lColElement(i);
    if (!ele) dserror("Cannot find element with lid %", i );

    // get corner position as node positions
    const int numcorner = 8;
    int bingid = ele->Id();
    std::vector<LINALG::Matrix<3,1> > bincorners;
    GetBinCorners( bingid , bincorners );

    // if element is a ghost
    if( ele->Owner() != myrank_ )
    {
      ghostcorners[ele->Id()] = bincorners;
      continue;
    }

    // add new node
    std::vector<double> cornerpos(3,0.0);
    std::vector<int> nids(8,0);
    for( int corner_i = 0; corner_i < numcorner; ++corner_i )
    {
      for( int dim = 0; dim < 3; ++dim )
        cornerpos[dim] = bincorners[corner_i](dim);

      nids[corner_i] = (bingid*numcorner) + corner_i;
      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp( new DRT::Node( nids[corner_i], &cornerpos[0], myrank_ ) );
      visbindis_->AddNode(newnode);
    }

    // assign nodes to elements
    Teuchos::RCP<DRT::Element> newele = DRT::UTILS::Factory( "VELE3", "Polynomial", ele->Id(), myrank_ );
    newele->SetNodeIds(nids.size(), &nids[0]);
    visbindis_->AddElement(newele);
  }

  // get max gid before adding elements
  int maxgid = bindis_->ElementRowMap()->MaxAllGID() + 1;
  if( writebinstype_ == INPAR::MESHFREE::cols )
  {
    // gather all numbers of ghosted bins that are going to be row eles
    std::vector<int> nummycol(1);
    nummycol[0] = static_cast<int>( ghostcorners.size() );
    // initialize std::vector for communication
    std::vector<int> numcol(com->NumProc(),0);
    // communicate
    com->GatherAll(&nummycol[0],&numcol[0],nummycol.size());
    com->Barrier();

    // calculate starting index on myrank
    int startnewgid = 0;
    for(int i=0; i<myrank_; ++i)
      startnewgid += numcol[i];

    // loop over all ghosted bins
    std::map<int, std::vector<LINALG::Matrix<3,1> > >::const_iterator iter;
    int counter = 0;
    for( iter = ghostcorners.begin(); iter != ghostcorners.end(); ++iter )
    {
      // new elegid (unique over all procs)
      int newelegid = maxgid + startnewgid + counter;
      counter++;

      // add new node
      // get corner position as node positions
      const int numcorner = 8;
      std::vector<double> cornerpos(3,0.0);
      std::vector<int> nids(8,0);
      for(int corner_i=0; corner_i<numcorner; ++corner_i)
      {
        for( int dim = 0; dim < 3; ++dim )
          cornerpos[dim] = iter->second[corner_i](dim);

        nids[corner_i] = (newelegid*numcorner) + corner_i;
        Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(new DRT::Node(nids[corner_i], &cornerpos[0], myrank_));
        visbindis_->AddNode(newnode);
      }

      // assign nodes to elements
      Teuchos::RCP<DRT::Element> newele = DRT::UTILS::Factory("VELE3","Polynomial", newelegid, myrank_);
      newele->SetNodeIds(nids.size(), &nids[0]);
      visbindis_->AddElement(newele);
    }
  }

  // complete new dis
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset = Teuchos::rcp( new DRT::IndependentDofSet(true) );
  visbindis_->ReplaceDofSet(independentdofset);
  visbindis_->FillComplete( true, true, false );

  // create vector that shows ghosting
  Teuchos::RCP<Epetra_Vector> ownedghostsvec = LINALG::CreateVector(*visbindis_->ElementRowMap(),true);
  for( int i = 0; i < visbindis_->NumMyRowElements(); ++i )
  {
    DRT::Element* ele = visbindis_->lRowElement(i);

    if( ele->Id() <= maxgid )
      (*ownedghostsvec)[i] = 0; // owned
    else
      (*ownedghostsvec)[i] = 1; // ghost
  }

  // write output
  visbindis_->Writer()->WriteMesh( step, time );
  visbindis_->Writer()->NewStep( step, time );
  visbindis_->Writer()->WriteVector( "owner0ghost1", ownedghostsvec, IO::elementvector );
  visbindis_->Writer()->WriteElementData( true );

  visbindis_->ClearDiscret();
}

/*----------------------------------------------------------------------*
| create linear map with bin ids                            ghamm 11/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> BINSTRATEGY::BinningStrategy::CreateLinearMapForNumbin(
    const Epetra_Comm& comm
    ) const
{
  // initial dummy distribution using a linear map
  const int numproc = comm.NumProc();
  const int numbin = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];
  const int start = numbin / numproc * myrank_;
  int end;
  // special treatment for last proc
  if(myrank_ != numproc-1)
    end = (int)(numbin / numproc * (myrank_+1));
  else
    end = numbin;

  std::vector<int> linearmap;
  linearmap.reserve(end-start);
  for(int k=0; k<bin_per_dir_[2]; ++k)
  {
    for(int j=0; j<bin_per_dir_[1]; ++j)
    {
      for(int i=0; i<bin_per_dir_[0]; ++i)
      {
        int curr = i + j*bin_per_dir_[0] + k*bin_per_dir_[0]*bin_per_dir_[1];
        if(start <= curr and curr < end)
          linearmap.push_back(i + j*id_calc_bin_per_dir_[0] + k*id_calc_bin_per_dir_[0]*id_calc_bin_per_dir_[1]);
      }
    }
  }

  return Teuchos::rcp(new Epetra_Map(numbin, linearmap.size(), &linearmap[0], 0, comm));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::GetCurrentNodePos(
    Teuchos::RCP< const DRT::Discretization> const discret,
    DRT::Node const * node,
    Teuchos::RCP<const Epetra_Vector>const  disnp,
    double* currpos ) const
{
  // Todo make this nicer

  // the problem is that we might have nodes without position DoFs
  // (e.g. for beam elements with 'interior' nodes that are only used for
  // triad interpolation)
  // instead of the node position itself, we return the position of the
  // first node of the  element here (for the sake of binning)

  // standard case
  DRT::Node const* node_with_position_Dofs = node;

  const DRT::Element* element = node->Elements()[0];
  const DRT::ELEMENTS::Beam3Base* beamelement =
      dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element);

  // if the node does not have position DoFs, we return the position of the first
  // node of the corresponding element
  if ( beamelement != NULL and not beamelement->IsCenterlineNode(*node) )
  {
    node_with_position_Dofs = beamelement->Nodes()[0];
  }

  if( disnp != Teuchos::null )
  {
    const int gid = discret->Dof(node_with_position_Dofs, 0);
    const int lid = disnp->Map().LID(gid);
    if( lid < 0 )
      dserror("Your displacement is incomplete (need to be based on a column map"
              " as this function is also called from a loop over elements and "
              "each proc does (usually) not own all nodes of his row elements ");
    for( int dim = 0; dim < 3; ++dim )
    {
      currpos[dim] = node_with_position_Dofs->X()[dim] + (*disnp)[lid+dim];
    }
  }
  else
  {
    for( int dim = 0; dim < 3; ++dim )
      currpos[dim] = node_with_position_Dofs->X()[dim];
  }

  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::PeriodicBoundaryShift3D( LINALG::Matrix< 3, 1 >& d,
    LINALG::Matrix< 3, 1 > const X ) const
{
  if ( not havepbc_ )
    return;

  if ( deforming_simulation_domain_handler != Teuchos::null )
  {
    deforming_simulation_domain_handler->Shift3D( d, X );
  }
  else
  {
    for( int dim = 0; dim < 3 ; ++dim )
    {
      if ( not pbconoff_[dim] )
        continue;

      while( d(dim) + X(dim) < XAABB_( dim, 0 ) )
        d(dim) += pbcdeltas_[dim];

      while( d(dim) + X(dim) > XAABB_( dim, 1 ) )
        d(dim) -= pbcdeltas_[dim];
    }
  }
  return;
}

/*-----------------------------------------------------------------------------*
| transfer nodes and elements of input discret                 eichinger 11/16 |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::TransferNodesAndElements(
    Teuchos::RCP<DRT::Discretization> & discret,
    Teuchos::RCP<Epetra_Vector> disnp,
    std::map<int, std::set<int> >& bintorowelemap)
{
  TEUCHOS_FUNC_TIME_MONITOR("BINSTRATEGY::BinningStrategy::TransferNodesAndElements");

  // store elements that need treatment
  std::set< DRT::Element* > elestoupdate;

  double currpos[3] = { 0.0, 0.0, 0.0 };
  // loop over all column nodes and check ownership
  for (int i = 0; i < discret->NodeColMap()->NumMyElements(); ++i )
  {
    // get current node and position
    DRT::Node* currnode = discret->lColNode(i);
    GetCurrentNodePos( discret, currnode, disnp, currpos );

    int const gidofbin = ConvertPosToGid(currpos);

    if( bindis_->HaveGlobalElement(gidofbin) )
    {
      int const hostbinowner = bindis_->gElement(gidofbin)->Owner();
      if( currnode->Owner() != hostbinowner )
      {
        // set new owner of node
        currnode->SetOwner( hostbinowner );
        // in case myrank is owner of associated element, add it to set
        DRT::Element** curreles = currnode->Elements();
        for( int j = 0; j < currnode->NumElement(); ++j )
          if( curreles[j]->Owner() == myrank_ )
            elestoupdate.insert(curreles[j]);
      }
    }
    /*else: in this case myrank was not owner of node and a corresponding element had and
    will only have ghost nodes on myrank, therefore we can leave the old owner because
    during the built up of the node col map all ghost nodes get deleted anyway  */
  }

  // store elements that need to be communicated
  std::map< int, std::vector<DRT::Element*> > toranktosendeles;
  std::map< int, std::vector<std::pair< int, std::vector<int> > > > toranktosendbinids;

  // loop over row elements whose ownership may need an update
  std::set< DRT::Element* >::const_iterator eleiter;
  for( eleiter = elestoupdate.begin(); eleiter != elestoupdate.end(); ++eleiter )
  {
    DRT::Element* currele = *eleiter;
    DRT::Node** nodes = currele->Nodes();
    std::map< int, int > owner;
    for( int inode = 0; inode < currele->NumNode(); ++inode )
      owner[nodes[inode]->Owner()] += 1;

    // check if any proc owns more nodes than myrank (for same number myrank_ stays owner)
    int newowner = myrank_;
    int numowned = ( owner.find( myrank_ ) != owner.end() ) ? owner[myrank_] : -1;
    std::map< int, int >::const_iterator i;
    for( i = owner.begin(); i != owner.end(); ++i )
      if( i->second > numowned )
        newowner = i->first;

    if ( newowner != myrank_ )
    {
      currele->SetOwner(newowner);
      toranktosendeles[newowner].push_back(currele);
      std::vector<int> binids;
      DistributeElementToBinsUsingEleXAABB( discret, currele, binids, disnp );
      std::pair<int, std::vector<int> > dummy(currele->Id(), binids);
      toranktosendbinids[newowner].push_back(dummy);
    }
  }

  // exploit bounding box idea for elements in underlying discretization and bins
  // loop over all row elements
  for ( int lid = 0; lid < discret->NumMyRowElements(); ++lid )
  {
    DRT::Element* eleptr = discret->lRowElement(lid);

    // check if owner did change
    if( eleptr->Owner() != myrank_ )
      continue;

    // get corresponding bin ids in ijk range
    std::vector<int> binIds;
    DistributeElementToBinsUsingEleXAABB( discret, eleptr, binIds, disnp );

   // assign element to bins
    std::vector< int >::const_iterator biniter;
    for( biniter = binIds.begin(); biniter != binIds.end(); ++biniter )
      bintorowelemap[*biniter].insert( eleptr->Id() );
  }

  // todo: send in one package
  // send and receive elements
  CommunicateElements( discret, toranktosendeles );
  // send and receive new elements to bin relation, like this no fillcomplete call necessary here
  CommunicateDistributionOfTransferredElementsToBins( discret, toranktosendbinids, bintorowelemap );
}

/*-----------------------------------------------------------------------------*
 | communicte elements to transfer                             eichinger 11/16 |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CommunicateElements(
    Teuchos::RCP<DRT::Discretization> & discret,
    std::map< int, std::vector<DRT::Element*> > const& toranktosendeles) const
{
  // build exporter
  DRT::Exporter exporter( discret->Comm() );
  int const numproc = discret->Comm().NumProc();

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // ---- pack data for sending -----
  std::map<int, std::vector<char> > sdata;
  std::vector<int> targetprocs( numproc, 0 );
  std::map< int, std::vector<DRT::Element*> >::const_iterator p;
  for( p = toranktosendeles.begin(); p != toranktosendeles.end(); ++p )
  {
    std::vector<DRT::Element*>::const_iterator iter;
    for( iter = p->second.begin(); iter != p->second.end(); ++iter )
    {
     DRT::PackBuffer data;
     (*iter)->Pack(data);
     data.StartPacking();
     (*iter)->Pack(data);
     sdata[p->first].insert( sdata[p->first].end(), data().begin(), data().end() );
    }
    targetprocs[p->first] = 1;
  }

  // ---- send ----
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for(std::map<int, std::vector<char> >::const_iterator p = sdata.begin(); p != sdata.end(); ++p )
  {
    exporter.ISend( myrank_, p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets( numproc, 0) ;
  discret->Comm().SumAll( targetprocs.data(), summedtargets.data(), numproc );

  // ---- receive ----
  for( int rec = 0; rec < summedtargets[myrank_]; ++rec )
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
        DRT::ParObject::ExtractfromPack( index, rdata, data );
        // this Teuchos::rcp holds the memory of the node
        Teuchos::RCP< DRT::ParObject > object = Teuchos::rcp( DRT::UTILS::Factory(data), true );
        Teuchos::RCP<DRT::Element> element = Teuchos::rcp_dynamic_cast< DRT::Element >(object);
        if (element == Teuchos::null) dserror("Received object is not a element");

        // safety check
        if( discret->HaveGlobalElement( element->Id() ) != true)
          dserror("%i is getting owner of element %i without having it ghosted before, "
                  "this is not intended.", myrank_, element->Id() );

        // delete already existing element (as it has wrong internal variables)
        discret->DeleteElement( element->Id() );
        // add node (ownership already adapted on sending proc)
        discret->AddElement(element);
      }
      if ( index != rdata.size() )
        dserror("Mismatch in size of data %d <-> %d",static_cast<int>( rdata.size() ), index );
    }
  }

  // wait for all communications to finish
  for ( int i = 0; i < length; ++i )
    exporter.Wait( request[i] );
  // safety, should be a no time operation if everything works fine before
  discret->Comm().Barrier();
}

/*-----------------------------------------------------------------------------*
 | communicate element to bin distribution                     eichinger 11/16 |
 *-----------------------------------------------------------------------------*/
void BINSTRATEGY::BinningStrategy::CommunicateDistributionOfTransferredElementsToBins(
    Teuchos::RCP<DRT::Discretization> & discret,
    std::map< int, std::vector< std::pair< int, std::vector<int> > > > const& toranktosendbinids,
    std::map<int, std::set<int> >& bintorowelemap) const
{
  // build exporter
  DRT::Exporter exporter( discret->Comm() );
  int const numproc = discret->Comm().NumProc();

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // ---- pack data for sending -----
  std::map<int, std::vector<char> > sdata;
  std::vector<int> targetprocs( numproc, 0 );
  std::map< int, std::vector<std::pair< int, std::vector<int> > > >::const_iterator p;
  for( p = toranktosendbinids.begin(); p != toranktosendbinids.end(); ++p )
  {
    std::vector<std::pair< int, std::vector<int> > >::const_iterator iter;
    for( iter = p->second.begin(); iter != p->second.end(); ++iter )
    {
      DRT::PackBuffer data;
      DRT::ParObject::AddtoPack( data, *iter );
      data.StartPacking();
      DRT::ParObject::AddtoPack( data, *iter );
      sdata[p->first].insert( sdata[p->first].end(), data().begin(), data().end() );
    }
   targetprocs[p->first] = 1;
  }

  // ---- send ----
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for(std::map<int, std::vector<char> >::const_iterator p = sdata.begin(); p != sdata.end(); ++p )
  {
    exporter.ISend( myrank_, p->first, &((p->second)[0]), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets( numproc, 0) ;
  discret->Comm().SumAll( targetprocs.data(), summedtargets.data(), numproc );

  // ---- receive ----
  for( int rec = 0; rec < summedtargets[myrank_]; ++rec )
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
        std::pair< int, std::vector<int> > pair;
        DRT::ParObject::ExtractfromPack( index, rdata, pair );
        std::vector<int>::const_iterator j;
        for( j = pair.second.begin(); j != pair.second.end(); ++j )
          bintorowelemap[*j].insert(pair.first);
      }
      if ( index != rdata.size() )
        dserror("Mismatch in size of data %d <-> %d",static_cast<int>( rdata.size() ), index );
    }
  }

  // wait for all communications to finish
  for ( int i = 0; i < length; ++i )
    exporter.Wait( request[i] );
  // safety, should be a no time operation if everything works fine before
  discret->Comm().Barrier();
}

