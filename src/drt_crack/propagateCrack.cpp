/*----------------------------------------------------------------------*/
/*!
\file PropagateCrack.cpp

\brief After each time step, check the structure field and propagate
crack in the structure if necessary.

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "propagateCrack.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_factory.H"

#include "../drt_mat/elasthyper.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_matelast/elast_coupneohooke.H"

#include "../linalg/linalg_utils.H"

/*------------------------------------------------------------------------------------*
 * Constructor                                                                     sudhakar 12/13
 * Read crack tip nodes and nodes falling on the crack surface from input data
 *------------------------------------------------------------------------------------*/
DRT::CRACK::PropagateCrack::PropagateCrack( Teuchos::RCP<DRT::Discretization>& discret )
:discret_(discret),
 comm_( discret_->Comm() ),
 myrank_( comm_.MyPID() )
{

  PI_ = 22.0 / 7.0;

  // get the initial crack tip nodes when analyzing propagation of an existing crack
  DRT::Condition* crackpts = discret_->GetCondition( "CrackInitiationPoints" );

  const std::vector<int>* tipnodes = const_cast<std::vector<int>* >(crackpts->Nodes());
  tipnodes_.resize( tipnodes->size() );

  for( unsigned i=0; i< tipnodes->size(); i++ )
  {
    int val = (*tipnodes)[i];
    tipnodes_[i] = val;
  }

  DRT::Condition* maspts = discret_->GetCondition( "masterCrackSurface" );
  DRT::Condition* slapts = discret_->GetCondition( "slaveCrackSurface" );

  const std::vector<int>* masternodes = const_cast<std::vector<int>* >(maspts->Nodes());
  const std::vector<int>* slavenodes = const_cast<std::vector<int>* >(slapts->Nodes());

  if( masternodes->size() != slavenodes->size() )
    dserror("There should be equal number of master and slave nodes\n");

  if( masternodes->size() == 0 )
    dserror("No master nodes defined. Are you dreaming of simulating crack initiation?\n");

  for( unsigned i=0; i<masternodes->size(); i++ )
  {
    cracknodes_.insert( (*masternodes)[i] );
    cracknodes_.insert( (*slavenodes)[i] );

    // since old tip nodes are not available in the initial time step, we
    // copy all the crack surface nodes here so that appropriate neighbors can be found
    oldTipnodes_.insert( (*masternodes)[i] );
    oldTipnodes_.insert( (*slavenodes)[i] );
  }

  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  DRT::Element * tempele = discret_->lRowElement(0);

  Teuchos::RCP<MAT::Material> actmat = tempele->Material();

  switch(actmat->MaterialType())
  {
    case INPAR::MAT::m_elasthyper:
    {
      MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
      if ( not params ) dserror("Cannot cast material parameters");
      if ( not params->nummat_ == 1 ) dserror("At the moment, not possible");
      int matid = (*(params->matids_))[0];

      const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(matid)->second;

      switch (actelastmat->Type())
      {
        case INPAR::MAT::mes_coupneohooke:
        {
          MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                        dynamic_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
          young_ = params2->youngs_;
          poisson_ = params2->nue_;
          break;
        }
        default:
          dserror("material model not supported");
          break;
      }
      break;
    }
    default:
      dserror("material type  not supported for crack simulations");
      break;
  }

  const Teuchos::ParameterList& crackparam = DRT::Problem::Instance()->CrackParams();

  critical_K_I_ = crackparam.get<double>("CRITICAL_K1");
  critical_K_II_ = crackparam.get<double>("CRITICAL_K2");

  std::string thick(Teuchos::getNumericStringParameter(crackparam,"THICKNESS_ASSUMPTION"));

  // Calculate kappa_ based on thickness assumption
  if( thick == "plane_stress" )
    kappa_ = ( 3.0 - poisson_ ) / ( 1.0 + poisson_ );
  else if ( thick == "plane_strain" )
    kappa_ = 3.0 - poisson_;

  fac_ = young_ / ( 1.0 + poisson_ ) / ( 1.0 + kappa_ );
}

/*------------------------------------------------------------------------------------*
 * Perform all the operations related to crack propagation                  sudhakar 11/13
 * This calculates stress intensity factor (K), and if it is  higher than
 * the critical value, crack is introduced into discretization
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::propagateOperations( Teuchos::RCP<const Epetra_Vector>& displace )
{
  oldnew_.clear();
  tip_phi_.clear();
  tip_mphi_.clear();

  std::cout<<"tip nodes = "<<tipnodes_[0]<<"\t"<<tipnodes_[1]<<"\n";

  //-------------
  // STEP 1 : Compute stress-intensity factors at crack tip
  //-------------
  findStressIntensityFactor( displace );

  // this arises at the very beginning of time integration
  // physically uncorrect
  // experienced this when using gen-alpha method  (sudhakar 12/13)
  if( K_I_ < 0.0 )
    return;

  //-------------
  // STEP 2 : Check whether propagation criterion is satisfied
  //-------------
  bool isProp = DoCrackPropagate();
  if( not isProp )
    return;

  //-------------
  // STEP 3 : Decide crack propagation angle from stress-intensity factors
  //-------------
  decidePropagationAngle();

  //-------------
  // STEP 4 : Update crack information
  //-------------
  std::vector<int> newTip = findNewCrackTip();

  //-------------
  // STEP 5 : Update crack information
  //-------------
  updateCrack( newTip );
}

/*------------------------------------------------------------------------------------*
 * In order to calculate stress intensity factors at the crack tip, we        sudhakar 12/13
 * need appropriate connectivity of few nodes at the crack tip
 * VERY IMPORTANT : at the moment it works only for brick (hex) elements
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::buildNeighborInfoTipNodes( DRT::Node * node1, DRT::Node * node2, DRT::Node * tipnode )
{

#ifdef DEBUG
  if( node1->Id() == node2->Id() )
    dserror( "We should have got two different node ids\n" );

  const double * cord1 = node1->X();
  const double * cord2 = node2->X();

  if( fabs( cord1[0] - cord2[0] ) > 1e-15 or
      fabs( cord1[1] - cord2[1] ) > 1e-15 or
      fabs( cord1[2] - cord2[2] ) > 1e-15 )
  {
    dserror( "two neighboring nodes should be at the same position\n" );
  }
#endif

  // we are considering only the first node to construct normal
  // because we assume through-thickness crack with crack properties same at each thickness plane
  //TODO: while extending to 3D crack, neighbors has to be calculated at each crack tip node

  DRT::Node * attach1 = findAttachedNode( node1, tipnode );
  DRT::Node * attach2 = findAttachedNode( node2, tipnode );

  const double * tipcord = tipnode->X();
  const double * atcord1 = attach1->X();
  const double * atcord2 = attach2->X();

  double tangDist[3];
  for( int dim = 0; dim < 3; dim++ )
    tangDist[dim] = tipcord[dim] + tangent_(dim,0);

  double dist1 = sqrt( pow((tangDist[0]-atcord1[0]),2) + pow((tangDist[1]-atcord1[1]),2) + pow((tangDist[2]-atcord1[2]),2) );
  double dist2 = sqrt( pow((tangDist[0]-atcord2[0]),2) + pow((tangDist[1]-atcord2[1]),2) + pow((tangDist[2]-atcord2[2]),2) );

  std::pair<DRT::Node*, DRT::Node*> surNodes;

  if( dist1 < dist2 )
  {
    tip_phi_[tipnode->Id()] = node1->Id();
    tip_mphi_[tipnode->Id()] = node2->Id();
  }
  else
  {
    tip_phi_[tipnode->Id()] = node2->Id();
    tip_mphi_[tipnode->Id()] = node1->Id();
  }
}

/*------------------------------------------------------------------------------------*
 * Build normal coordinate system at the crack tip                            sudhakar 12/13
 *------------------------------------------------------------------------------------*/
DRT::Node * DRT::CRACK::PropagateCrack::findNormal( const DRT::Node * tipnode )
{
  // we are considering only the first node to construct normal
  // because we assume through-thickness crack with crack properties same at each thickness plane
  //TODO: while extending to 3D crack, normal has to be calculated at each crack tip segment

  // In order to construct local coodinate system at crack tip, we consider one crack tip node
  // We need to find another node that is connected to the tip node and is located in the crack
  // surface (but not in the crack tip)
  // REMEMBER : This procedure works only for through-thickness crack

  const double * tipcoord = tipnode->X();

  DRT::Node * surnode = findNeighboringCrackNode( tipnode, false, NULL );
  const double * surcoord = surnode->X();

  normal_(0,0) = tipcoord[0] - surcoord[0];
  normal_(1,0) = tipcoord[1] - surcoord[1];
  normal_(2,0) = 0.0;

  double fac = sqrt( normal_(0,0)*normal_(0,0) + normal_(1,0)*normal_(1,0) + normal_(2,0)*normal_(2,0) );
  normal_(0,0) = normal_(0,0) / fac;
  normal_(1,0) = normal_(1,0) / fac;
  normal_(2,0) = normal_(2,0) / fac;

  // ---------------------
  // calculating the proper direction of tangent vector is not straight forward (it can be in either direction)
  // We do the following to decide correct tangent direction
  // If the crack is located in horizontal direction, the normal and tangent are (1,0) and (0,1)
  // we calculated the proper normal direction (indirectly) using crack surface orientation
  //
  //
  //
  //                                                        ^ (nx,ny)
  //                         ^ (0,1)                       .
  //                         |  tangent                   .\ theta
  //                         |                           .----------->(1,0)
  //                         |                          //
  //                         |         (1,0)           //
  //          ===============----------> normal       //
  //                                                 //
  //
  //
  // we decide the angle between (nx,ny) and (1,0) --> theta
  // then we rotate (0,1) to angle theta, and we get the required tangent unit vector
  // ---------------------

  double theta = atan2( normal_(1,0), normal_(0,0) ) - atan2( 0.0, 1.0 );

  // apply linear transformation to rotate (0,1) to angle theta
  tangent_(0,0) = -1.0 * sin(theta);
  tangent_(1,0) = cos(theta);
  tangent_(2,0) = 0.0;

  return surnode;
}

/*------------------------------------------------------------------------------------*
 * Calculate stress intensity factors at crack tip                             sudhakar 12/13
 *
 *
  // For a crack tip node, we need neighboring nodes ("o") to calculate the stress intensity factor
  // In order to do so, we get the elements that are attached with "o" which should be 4
  // For one node "o" it should be ele1, ele2 and two other elements in z-direction
  // For another node "o" it should be ele3, ele4 and two other elements in z-direction
  // we choose the elements that are close to crack tip, they are ele1 and ele3 for the two nodes
  //
  //
  //                                                                       =====  crack surface
  //                                                                           *  tip node
  //                                               ^ tangent                   o  nodes on crack surface (two nodes at same position)
  //                                               !                          "o" neighboring nodes of crack tip node
  //                     ..........................!...                        #  attached nodes
  //                    /             /            ! /|
  //                   /_____________#_____________!/ |
  //                   |             |             !  |
  //                   |             |             !  |
  //                   |  (ele2)     |   (ele1)    !  |
  //                   |             |             ! /
  //                ===o============"o"============*--------------> normal
  //                   |             |             !  |
  //                   |             |             !  |
  //                   |  (ele4)     |    (ele3)   !  |
  //                   |             |             ! /
  //                   ..............#..............
  //
  //
  // Nodes marked as "o" are used to determine stress intensity factors
  // However in order to determine the correct sign for K_II, we need to decide which node is at theta=pi and which is at -pi
  // In order for this, we decide the nodes (#) attached with "o"
  // we mark a point whose coordinates are newpt = X_tip + X_tangent_
  // calculate distance between newpt and #1, #2
  // Which node is closer to newpt is at theta=pi and another is at -pi
   *
   *
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::findStressIntensityFactor( Teuchos::RCP<const Epetra_Vector>& displace )
{
  DRT::Node * tipnode = NULL;
  int lmaster = 0, gmaster = 0;

  if( discret_->HaveGlobalNode( tipnodes_[0] ) )
  {
    tipnode = discret_->gNode( tipnodes_[0] );
    if( tipnode->Owner() == myrank_ )
    {
      lmaster = myrank_;
      DRT::Node * node1 = findNormal( tipnode );
      DRT::Node * node2 = findNeighboringCrackNode( tipnode, true, node1 );

      // find neighbor points to calculate stress-intensity factor
      buildNeighborInfoTipNodes( node1, node2, tipnode );
    }
  }

  // making sure that master processor id is available at all processors
  comm_.SumAll( &lmaster, &gmaster, 1 );

  // normal and tangent vectors are computed only on master processor
  // broadcast them to all processors
  comm_.Broadcast( &normal_(0,0), 3, gmaster );
  comm_.Broadcast( &tangent_(0,0), 3, gmaster );

  LINALG::GatherAll( tip_phi_, discret_->Comm() );
  LINALG::GatherAll( tip_mphi_, discret_->Comm() );

  // export "displacement" to column map
  Teuchos::RCP<Epetra_Vector> disp_col = LINALG::CreateVector( *discret_->DofColMap(), true );
  LINALG::Export(*displace,*disp_col);

  K_I_ = 0.0; K_II_ = 0.0;
  double local_k1 = 0.0, local_k2 = 0.0;
  for( std::map<int,int>::iterator itip = tip_phi_.begin(); itip != tip_phi_.end(); itip++ )
  {
    int tip_id = itip->first;
    int phi_id = itip->second;
    int mphi_id = tip_mphi_.find( tip_id )->second;

    // we do calculation on the owner of a tipnode, as only this processor
    // has the information about the neighboring nodes
    if( discret_->HaveGlobalNode( tip_id ) )
    {
      tipnode = discret_->gNode( tip_id );
      if( tipnode->Owner() == myrank_ )
      {
        if( not discret_->HaveGlobalNode( phi_id ) or
            not discret_->HaveGlobalNode( mphi_id ) )
        {
          dserror( "owner of tip node must know the neighboring nodes\n" );
        }

        DRT::Node * node_phi = discret_->gNode( phi_id );
        DRT::Node * node_mphi = discret_->gNode( mphi_id );

        std::vector<double> disp_phi = getDisplacementNode( node_phi, disp_col );
        std::vector<double> disp_mphi = getDisplacementNode( node_mphi, disp_col );

        const double * tipcord = tipnode->X();
        const double * cord_phi = node_phi->X();
        const double * cord_mphi = node_mphi->X();

        double dist1 = sqrt( pow((tipcord[0]-cord_phi[0]),2) +
                             pow((tipcord[1]-cord_phi[1]),2) +
                             pow((tipcord[2]-cord_phi[2]),2) );

        double dist2 = sqrt( pow((tipcord[0]-cord_mphi[0]),2) +
                             pow((tipcord[1]-cord_mphi[1]),2) +
                             pow((tipcord[2]-cord_mphi[2]),2) );

        double avg_dist = 0.5*( dist1 + dist2 );

        double tang_phi = disp_phi[0] * tangent_(0,0) +
                          disp_phi[1] * tangent_(1,0) +
                          disp_phi[2] * tangent_(2,0);

        double tang_mphi = disp_mphi[0] * tangent_(0,0) +
                           disp_mphi[1] * tangent_(1,0) +
                           disp_mphi[2] * tangent_(2,0);

        double norm_phi = disp_phi[0] * normal_(0,0) +
                          disp_phi[1] * normal_(1,0) +
                          disp_phi[2] * normal_(2,0);

        double norm_mphi = disp_mphi[0] * normal_(0,0) +
                           disp_mphi[1] * normal_(1,0) +
                           disp_mphi[2] * normal_(2,0);

        local_k1 += fac_*sqrt( 0.5*PI_/avg_dist ) * ( tang_phi - tang_mphi );
        local_k2 += fac_*sqrt( 0.5*PI_/avg_dist ) * ( norm_phi - norm_mphi );
      }
    }
  }

  comm_.SumAll( &local_k1, &K_I_, 1 );
  comm_.SumAll( &local_k2, &K_II_, 1 );

  K_I_ = K_I_ / tip_phi_.size();
  K_II_ = K_II_ / tip_phi_.size();

  std::cout<<"stress intensity factors = "<<K_I_<<"\t"<<K_II_<<"\n";
}

/*------------------------------------------------------------------------------------*
 * Extract displacement at the given node                                    sudhakar 12/13
 *------------------------------------------------------------------------------------*/
std::vector<double> DRT::CRACK::PropagateCrack::getDisplacementNode( DRT::Node * node,
                                                                     Teuchos::RCP<Epetra_Vector>& disp )
{
  std::vector<int> lm;
  std::vector<double> mydisp;
  LINALG::Matrix<3,1> displ;

  discret_->Dof( node, lm );

  DRT::UTILS::ExtractMyValues( *disp, mydisp, lm );

  return mydisp;
}

/*------------------------------------------------------------------------------------*
 * Find neighboring node associated with the first crack tip node
 * The neighboring node is located on the crack surface, not on tip             sudhakar 12/13
 * if bool second = true, it means we already have one node available and
 * searching for another node located at the same position as the first node
 *------------------------------------------------------------------------------------*/
DRT::Node * DRT::CRACK::PropagateCrack::findNeighboringCrackNode( const DRT::Node * tipnode,
                                                                  bool second,
                                                                  const DRT::Node * first_node )
{
  if( second and first_node == NULL )
    dserror( "While finding second neighbor, one should supply the first node\n" );

  int tipnodeid = tipnode->Id();
  const double * tipcoord = tipnode->X();

  DRT::Node * surnode = NULL;
  int surnodeid;
  const double * surcoord = NULL;

  // get all elements that has this tipnode
  const DRT::Element* const * ElementsPtr = tipnode->Elements();
  int NumElement = tipnode->NumElement();

  bool foundnode = false;

  for (int jEle = 0; jEle < NumElement; jEle++)
  {
    const DRT::Element* Element = ElementsPtr[jEle];

    const int* nodes = Element->NodeIds();
    const int numNodes = Element->NumNode();

    for( int iNode = 0; iNode < numNodes; iNode++ )
    {
      surnodeid = nodes[iNode];

      // check to make sure not the same node
      if( surnodeid == tipnodeid )
        continue;

      // make sure we are not getting the same node again
      // here we break the loop, because first and second nodes are always in different elements
      if( second and surnodeid == first_node->Id() )
        break;

      // neighboring node should not be in crack tip
      if( std::find( tipnodes_.begin(), tipnodes_.end(), surnodeid ) != tipnodes_.end() )
        continue;

      // the node is not located
      if( std::find( oldTipnodes_.begin(), oldTipnodes_.end(), surnodeid ) != oldTipnodes_.end() )
      {
        if(not  discret_->HaveGlobalNode( surnodeid ) )
          dserror("surrounding node not found on this processor\n");

        surnode = discret_->gNode( surnodeid );
        surcoord = surnode->X();

        if( fabs( tipcoord[2] - surcoord[2] ) > 1e-12 )
          continue;

        foundnode = true;
        break;
      }
    }

    if( foundnode )
      break;
  }

  if( not foundnode )
    dserror("For the tip node, not able to find connected surface node\n");

  return surnode;
}

/*--------------------------------------------------------------------------------------------------*
 * For each neighboring node, we find an attached node that satisfied certain criterion      sudhakar 12/13
 * These attached nodes are used to decide the proper sign of K_II_.
 * Refer to the explanatory figure in this file for clear understanding
 *--------------------------------------------------------------------------------------------------*/
DRT::Node * DRT::CRACK::PropagateCrack::findAttachedNode( DRT::Node * neigh, DRT::Node * tipnode )
{
  DRT::Node * attach = NULL;
  const double * neighcord = neigh->X();

  // find all the elements attached with this neighboring node
  DRT::Element** ElementsPtr = neigh->Elements();
  int NumElement = neigh->NumElement();

  if( NumElement > 4 )
    dserror( "For Hex8 elements, the crack point should have max four elements attached to it\n" );

  bool foundele = false;

  DRT::Element * atele = NULL;
  if( NumElement == 1 )
  {
    atele = ElementsPtr[0];
    foundele = true;
  }

  // if there are more than one elements, we choose the one which is shares the crack tip node
  // there can be more than one element and it does not matter which one we choose
  else
  {
    // calculation of centre coordinate by checking whether tipnode is present
    for( int num = 0; num < NumElement; num++ )
    {
      DRT::Element * ele = ElementsPtr[num];

      if( discret_->HaveGlobalElement( ele->Id() ) )
      {
        const int * nodeids = ele->NodeIds();

        for( int no = 0; no < ele->NumNode(); no++ )
        {
          const int id = nodeids[no];
          if( std::find( tipnodes_.begin(), tipnodes_.end(), id ) != tipnodes_.end() )
          {
            foundele = true;
            atele = ele;
            break;
          }
        }

        if( foundele )
          break;
      }
    }
  }

  if( not foundele )
    dserror( "atleast one neighboring element should contain the tipnode\n" );

  // Take all surfaces of "atele"
  // we already confirmed this element is in the present processor
  std::vector< Teuchos::RCP< DRT::Element > > Surfaces = atele->Surfaces();

  Teuchos::RCP<DRT::Element> surf = getSurfaceSameZplane( Surfaces, neighcord );

  // There are two nodes attached with this neighboring node in this elements
  // One node is the tip node, and other is the attached node
  // once we check which is the tipnode, we can find the attached node
  DRT::Node ** searchnodes = surf->Nodes();
  int neigh_index = 0;
  bool found_neigh = false;
  int numnodes = surf->NumNode();
  for (int num = 0; num < numnodes; num++)
  {
    DRT::Node * nod = searchnodes[num];
    if( nod == neigh )
    {
      found_neigh = true;
      neigh_index = num;
    }
  }

  if( not found_neigh )
    dserror("The surface we got does not contain the neighboring node that is given as input\n");

  DRT::Node * temp = searchnodes[(neigh_index+1)%numnodes];

  if( std::find( tipnodes_.begin(), tipnodes_.end(), temp->Id() ) == tipnodes_.end() )
    attach = temp;

  else
  {
    if( neigh_index == 0 )
      attach = searchnodes[numnodes-1];
    else
      attach = searchnodes[neigh_index-1];
  }

  return attach;
}

/*------------------------------------------------------------------------------------*
 * For the given element, get the surface that lies in the z-plane of "coord      sudhakar 12/13
 *------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::CRACK::PropagateCrack::getSurfaceSameZplane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces,
                                                                             const double * coord )
{
  Teuchos::RCP<DRT::Element> surf = Teuchos::null;
  bool foundsurf = true;

  // which surfaces has all the nodes in considered z-plane
  for (unsigned int surfele = 0; surfele < Surfaces.size(); surfele++)
  {
    foundsurf = true;
    surf = Surfaces[surfele];

    DRT::Node ** searchnodes = surf->Nodes();
#ifdef DEBUG
    if (searchnodes==NULL) dserror("No Nodes in Surface Element");
#endif
    for (int num = 0; num < surf->NumNode(); num++)
    {
      DRT::Node * nod = searchnodes[num];
      const double * atcord = nod->X();

      if( fabs( atcord[2] - coord[2] ) > 1e-12 )
      {
        foundsurf = false;
        break;
      }
    }

    if( foundsurf )
      break;
  }

  if( not foundsurf )
    dserror( "the required surface is not found\n" );

  if( surf == Teuchos::null )
    dserror( "the required surface is not found\n" );

  return surf;
}

/*------------------------------------------------------------------------------------*
 * Decide crack propagation angle from stress-intensity factors               sudhakar 01/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::decidePropagationAngle()
{
  double deno = K_I_ + sqrt( K_I_ * K_I_ + 8.0 * K_II_ * K_II_ );
  propAngle_ = 2.0*atan( -2.0 * K_II_ / deno );

  // atan2() returns values in the range (-pi,pi)
  // in order to make meaningful comparisons we need it in range (0,2*pi)
  if( fabs(propAngle_) < 1e-10 )
    propAngle_ = 0.0;
  else if( propAngle_ < 0.0 )
    propAngle_ = 2 * PI_ + propAngle_;

  std::cout<<"propagation angle = "<<propAngle_ * 180.0 / PI_ <<"deg\n";
}

/*------------------------------------------------------------------------------------*
 * Returns true if the crack propagation criterion is satisfied               sudhakar 12/13
 *------------------------------------------------------------------------------------*/
bool DRT::CRACK::PropagateCrack::DoCrackPropagate()
{
  double check = pow( (K_I_ / critical_K_I_), 2.0 ) + pow( (K_II_ / critical_K_II_), 2.0 );
  if( check < 0.99999999999999999 )
    return false;
  return true;
}

/*------------------------------------------------------------------------------------*
 * Make all the modifications in the discretization related to crack propagation
 * Duplicate crack tip nodes, and modify element connectivity               sudhakar 12/13
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::updateCrack( std::vector<int>& newTip )
{
  //
  // Now we need to duplicate the crack tip nodes, and modify the connectivity
  // In 2D, 4 elements share crack tip nodes, out of which two elements (above the crack tip)
  // retains the same node, and other 2 elements below crack tip get new duplicated node
  // To decide which nodes retain, and which nodes get new nodes. we do the following
  // 1. construct a vector from crack tip to in the direction of normal (vec1) -> normal vector
  // 2. Construct another vector from crack tip to centre point of each element (vec2)
  // If angle between vec1 and vec2 is more than PI, then the element gets a new duplicate node
  //
  //                               _________________________
  //                              |            |   ^        |
  //                              |            |  /vec2     |
  //                              |            | /          |
  //                ==============o============*.......>....o
  //                              |            |\    vec1   |
  //                              |            | \          |
  //                              |            |            |
  //                              o------------o------------o
  //
  //
  // The procedure is same in 3D, except that there are 8 elements, and 4 gets new duplicate node
  //

  int totalNodes = discret_->NumGlobalNodes(); // to decide the ids of new nodes

  std::vector<double> lmtAngle = getLimitAngles();

  // map of element ids to be modified with the new node
  // "key" contains the element id, and "value" is dummy here
  // map is used to make sure that one element is stored only once
  std::map<int, int> delEle;

  for( unsigned num = 0; num < tipnodes_.size(); num++ )
  {
    if( discret_->HaveGlobalNode( tipnodes_[num] ) )
    {
      DRT::Node * tipnode = discret_->gNode( tipnodes_[num] );

      Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( totalNodes + num, tipnode->X(), tipnode->Owner() ) );

      oldnew_[tipnode->Id()] = dupnode->Id();

      DRT::Element** ElementsPtr = tipnode->Elements();
      int NumElement = tipnode->NumElement();

      for (int jEle = 0; jEle < NumElement; jEle++)
      {
        DRT::Element* Element = ElementsPtr[jEle];
        if( toReplaceNode( Element, tipnode, lmtAngle ) )
          delEle[ Element->Id() ] = 0;
      }

      discret_->AddNode( dupnode );
    }
  }

  LINALG::GatherAll( oldnew_, discret_->Comm() );
  LINALG::GatherAll( delEle, discret_->Comm() );

  for( std::map<int,int>::iterator num = delEle.begin(); num != delEle.end(); num++ )
  {
    int eleid = num->first;
    if( discret_->HaveGlobalElement( eleid ) )
    {
      bool del = false;
      DRT::Element * ele = discret_->gElement( eleid );
      const int * oldnodes = ele->NodeIds();
      std::vector<int> newnodes( ele->NumNode() );

      for( int i = 0; i < ele->NumNode(); i++ )
      {
        std::map<int, int>::iterator delnod = oldnew_.find( oldnodes[i] );
        if( delnod != oldnew_.end() )
        {
          del = true;
          newnodes[i] = delnod->second;
        }
        else
          newnodes[i] = oldnodes[i];
      }

      if( not del )
        dserror("This element should have atleast one replaceable node\n");

      if( not newnodes.size() == ele->NumNode() )
        dserror("Check the number of new nodes\n");

#if 1
      //-----------------
      // modifying the nodes of an element is easy
      // just modify the node ids of the element
      // when fillcomplete is called the corresponding nodes will be
      //  set through DRT::Element::BuildNodalPointers()
      //-----------------
      ele->SetNodeIds(newnodes.size(), &newnodes[0]);
#endif

#if 0 // old way of deleting an element and adding a new element with correct nodes
      int own = ele->Owner();
      DRT::Element::DiscretizationType distype = ele->Shape();
      int material_id = ele->Material()->Parameter()->Id();

      discret_->DeleteElement( eleid );

      Teuchos::RCP<DRT::Element> newele = Teuchos::null;
      if( distype == DRT::Element::hex8 )
      {
        newele = DRT::UTILS::Factory( "SOLIDH8", "HEX8", eleid, own );
        // since we need "Evaluate" functionality on this element, we should
        // initialize the element with proper material
        newele->SetMaterial( material_id );
      }
      else
        dserror("This shape is not yet supported for crack simulations\n");

      newele->SetNodeIds(newnodes.size(), &newnodes[0]);

      if( myrank_ == own )
        discret_->AddElement( newele );
#endif
    }
  }

  modifyConditions();

  discret_->FillComplete();

  // update crack tip nodes and add new crack tip nodes to cracknodes_
  tipnodes_.clear();
  tipnodes_.resize( newTip.size() );
  std::copy( newTip.begin(), newTip.end(), tipnodes_.begin() );
  std::copy( newTip.begin(), newTip.end(), std::inserter( cracknodes_, cracknodes_.end() ) );


  oldTipnodes_.clear();
  for( std::map<int,int>::iterator it = oldnew_.begin(); it != oldnew_.end(); it++ )
  {
    oldTipnodes_.insert( it->first );
    oldTipnodes_.insert( it->second );
  }
}

/*------------------------------------------------------------------------------------*
 * Get limiting angles for this crack geometry and propagation angle         sudhakar 01/14
 * these angles are used to determine which elements get new nodes and which
 * of them keep the old nodes
 *------------------------------------------------------------------------------------*/
std::vector<double> DRT::CRACK::PropagateCrack::getLimitAngles()
{
  std::vector<double> ang;

  //the two angles are
  //1. negative of angle formed by normal
  //2. crack propagation angle
  double norm_ang = atan2( normal_(1,0), normal_(0,0) ) + PI_;
  if( norm_ang > 2*PI_ )
    norm_ang = norm_ang - 2*PI_;

  ang.push_back( norm_ang );
  ang.push_back( propAngle_ );

  std::sort( ang.begin(), ang.end() );

  return ang;
}

/*------------------------------------------------------------------------------------*
 * condition maps are for the discretization is built in the initial setup         sudhakar 01/14
 * After introducing new crack tip nodes, we modify the condition maps accordingly
 * (This is just a work-around. The vector copying operations can be avoided if
 * we can modify the way conditions are generated in Baci)
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::modifyConditions()
{
  std::multimap<std::string,Teuchos::RCP<Condition> > allcondn;
  discret_->GetAllConditions( allcondn );

  for( std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit = allcondn.begin();
                                                                     conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<Condition> cond = conit->second;

    for( std::map<int,int>::iterator it = oldnew_.begin(); it != oldnew_.end(); it++ )
    {
      int oldid = it->first;
      if( cond->ContainsNode( oldid ) )
      {
        const std::vector<int>* conNodes = cond->Nodes();

        Teuchos::RCP<std::vector<int> > storage = Teuchos::rcp(new std::vector<int>(conNodes->size()+1));
        std::vector<int>& access = *storage;

        for( unsigned ii = 0; ii < conNodes->size(); ii++ )
          access[ii] = (*conNodes)[ii];
        access[ conNodes->size() ] = it->second;

        cond->Add( "Node Ids", storage );
      }
    }
  }
}

/*------------------------------------------------------------------------------------*
 * Print all conditions of the discretization (just for debugging)           sudhakar 01/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::printConditions(std::multimap<std::string,Teuchos::RCP<Condition> > allcondn)
{
  std::cout<<"number of conditions = "<<allcondn.size()<<"\n";
  for( std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit = allcondn.begin();
                                                                       conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<Condition> cond = conit->second;
    std::cout<<"condn type = "<<cond->GType()<<" geom disc = "<<cond->GeometryDescription()<<"\n";
    const std::vector<int>* nodeids = cond->Nodes();

    std::cout<<"condition nodes are = ";
    for( unsigned i=0;i<nodeids->size();i++ )
    {
      std::cout<<(*nodeids)[i]<<"  ";
    }
    std::cout<<"\n";
  }
}

/*------------------------------------------------------------------------------------*
 * Returns true if the criteron to replace the tipnode with a new            sudhakar 01/14
 * duplicate node is satisfied
 *------------------------------------------------------------------------------------*/
bool DRT::CRACK::PropagateCrack::toReplaceNode( DRT::Element * ele, DRT::Node * tip, const std::vector<double>& lmtAngle )
{
  std::vector<double> cen = ele->ElementCenterRefeCoords();
  const double * tipcord = tip->X();

  // construct vector from tip node to middle of this element
  // For the procedure ref. to the figure in DRT::CRACK::PropagateCrack::updateCrack()
  double vec2[3];
  for( unsigned i=0; i<3; i++ )
    vec2[i] = cen[i] - tipcord[i];

  // find the angle between this vector and the normal vector
  double theta = atan2( vec2[1], vec2[0] );

  // atan2() returns values in the range (-pi,pi)
  // in order to make comparisons we need it in range (0,2*pi)
  if( fabs(theta) < 1e-12 )
    theta = 0.0;
  else if( theta < 0.0 )
    theta = 2 * PI_ + theta;

  // If this angle is within the limiting angles defined, then new node is allocated
  // this is already a sorted vector
  // so this comparison is fine
  if( theta > lmtAngle[0] and theta < lmtAngle[1] )
    return true;

  return false;
}

/*------------------------------------------------------------------------------------*
 *  find the new crack tip nodes.                                             sudhakar 12/13
 *  We do not update the crack tip nodes here, because we still need the
 *  old values to propagate crack from old position
 *------------------------------------------------------------------------------------*/
std::vector<int> DRT::CRACK::PropagateCrack::findNewCrackTip()
{
  std::map<int,int> oldnew_tip;

  for( unsigned tip = 0; tip < tipnodes_.size(); tip++ )
  {
    int tipid = tipnodes_[tip];
    int gnewtip = 0, lnewtip = 0;

    if( discret_->HaveGlobalNode( tipid ) )
    {
      DRT::Node* tipnode = discret_->gNode( tipid );
      if( tipnode->Owner() == myrank_ )
      {
        // nodes that are already used to calculate the crack angle
        std::set<int> procNodes;

        const double * tipcord = tipnode->X();
        double diff = 1000.0; //just a large number to compare the angle

        DRT::Element** ElementsPtr = tipnode->Elements();
        int NumElement = tipnode->NumElement();

        for( int eleno = 0; eleno < NumElement; eleno++ )
        {
          DRT::Element* Element = ElementsPtr[eleno];
          std::vector< Teuchos::RCP< DRT::Element > > Surfaces = Element->Surfaces();

          Teuchos::RCP<DRT::Element> surf = getSurfaceSameZplane( Surfaces, tipcord );

          if( surf == Teuchos::null )
            dserror("Did not found the surface on same z-plane\n");

          DRT::Node ** searchnodes = surf->Nodes();
          int tip_index = 0;
          bool found_tip = false;
          int numnodes = surf->NumNode();
          for (int num = 0; num < numnodes; num++)
          {
            DRT::Node * nod = searchnodes[num];
            if( nod == tipnode )
            {
              found_tip = true;
              tip_index = num;
            }
          }
          if( not found_tip )
            dserror("The surface we got does not contain the neighboring node that is given as input\n");

          DRT::Node * temp1 = searchnodes[(tip_index+1)%numnodes];
          int temp1id = temp1->Id();
          if( std::find( procNodes.begin(), procNodes.end(), temp1id ) == procNodes.end() )
          {
            const double * tempcord = temp1->X();
            double angl = atan2( (tempcord[1]-tipcord[1]), (tempcord[0]-tipcord[0]) );
            // atan2() returns values in the range (-pi,pi)
            // in order to make comparisons we need it in range (0,2*pi)
            if( fabs(angl) < 1e-12 )
              angl = 0.0;
            else if( angl < 0.0 )
              angl = 2 * PI_ + angl;

            if( fabs( angl - propAngle_ ) < diff )
            {
              diff = fabs( angl - propAngle_ );
              lnewtip = temp1id;
            }
            procNodes.insert( temp1id );
          }

          int temp_index = 0;
          if( tip_index == 0 )
            temp_index = numnodes-1;
          else
            temp_index = tip_index-1;
          DRT::Node * temp2 = searchnodes[temp_index];
          int temp2id = temp2->Id();
          if( std::find( procNodes.begin(), procNodes.end(), temp2id ) == procNodes.end() )
          {
            const double * tempcord = temp2->X();
            double angl = atan2( (tempcord[1]-tipcord[1]), (tempcord[0]-tipcord[0]) );

            // atan2() returns values in the range (-pi,pi)
            // in order to make comparisons we need it in range (0,2*pi)
            if( fabs(angl) < 1e-12 )
              angl = 0.0;
            else if( angl < 0.0 )
              angl = 2 * PI_ + angl;

            if( fabs( angl - propAngle_ ) < diff )
            {
              diff = fabs( angl - propAngle_ );
              lnewtip = temp2id;
            }
            procNodes.insert( temp2id );
          }
        }
      }
    }

    comm_.SumAll( &lnewtip, &gnewtip, 1 );
    oldnew_tip[tipid] = gnewtip;
  }

  if( not tipnodes_.size() == oldnew_tip.size() )
    dserror("for each node, we should have a new tip node\n");

  std::vector<int> newTip;
  for( std::map<int,int>::iterator it = oldnew_tip.begin(); it != oldnew_tip.end(); it++ )
    newTip.push_back( it->second );

  return newTip;
}
