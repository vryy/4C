/*----------------------------------------------------------------------*/
/*!
\file PropagateTip.cpp

\brief After each time step, check the structure field and propagate
this crack tip segment in the structure if necessary.

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "propagateTip.H"
#include "crack_tolerance.H"
#include "crackUtils.H"
#include "SplitHexIntoTwoWedges.H"
#include "j_Integral.H"

#include "../drt_mat/elasthyper.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_matelast/elast_coupneohooke.H"
#include "../drt_mat/stvenantkirchhoff.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils.H"

/*-----------------------------------------------------------------------------------------*
 *        Constructor                                                              sudhakar 09/14
 *-----------------------------------------------------------------------------------------*/
DRT::CRACK::PropagateTip::PropagateTip( Teuchos::RCP<DRT::Discretization> discret,
                                        std::vector<int> tipnodes,
                                        std::set<int> cracknodes,
                                        int segment_id )
:discret_(discret),
 tipnodes_(tipnodes),
 cracknodes_(cracknodes),
 segment_id_(segment_id),
 disp_col_(Teuchos::null),
 comm_( discret_->Comm() ),
 myrank_( comm_.MyPID() )
{
  PI_ = 22.0 / 7.0;
  min_angle_tol_ = MIN_PROP_ANGLE * PI_ / 180.0;

  alreadySplit_ = false;

  //---
  // Store ALE Dirichlet nodes as boundary nodes
  //---
  std::vector<DRT::Condition*> linecond;
  discret_->GetCondition( "ALEDirichlet", linecond );

  if( linecond.size() == 0 )
    dserror("No Dirichlet condition found for this discretization\n");

  for(unsigned i=0;i<linecond.size();i++)
  {
    DRT::Condition* con = linecond[i];
    if( con->Type() == DRT::Condition::LineDirichlet )
    {
      const std::vector<int>* nodes = con->Nodes();
      for( unsigned nodno = 0; nodno < nodes->size(); nodno++ )
        boun_nodes_.insert( (*nodes)[nodno] );
    }
  }

  //---
  // copy all surface nodes as oldtipnodes --> this will be used to construct normal and tangent vectors
  //---
  for( std::set<int>::iterator it = cracknodes.begin(); it != cracknodes.end(); it++ )
    oldTipnodes_.insert(*it);

  //---
  // get material id from element materials
  // This will be used to set materials approp when splitting HEX into WEDGEs
  //---
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  DRT::Element * tempele = discret_->lRowElement(0);

  Teuchos::RCP<MAT::Material> actmat = tempele->Material();

  bool mat_found = false;
  for( std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator itm = mats.begin(); itm != mats.end(); itm++ )
  {
    if( itm->second->Type() == actmat->MaterialType() )
    {
      material_id_ = itm->first;
      mat_found = true;
      break;
    }
  }

  if( not mat_found )
    dserror("given material is not found in material maps?\n");

  switch(actmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:
    {
      MAT::PAR::StVenantKirchhoff* params = dynamic_cast<MAT::PAR::StVenantKirchhoff*>(actmat->Parameter());
      if ( not params ) dserror("Cannot cast material parameters");
      young_ = params->youngs_;
      poisson_ = params->poissonratio_;
      break;
    }
    case INPAR::MAT::m_elasthyper:
    {
      MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
      if ( not params ) dserror("Cannot cast material parameters");
      if ( params->nummat_ != 1 ) dserror("At the moment, not possible");
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
    {
      std::cout<<"material type = "<<actmat->MaterialType()<<"\n";
      dserror("This material type not supported for crack simulations");
      break;
    }
  }

  //---
  // Calculate kappa_ based on thickness assumption
  //---
  const Teuchos::ParameterList& crackparam = DRT::Problem::Instance()->CrackParams();
  std::string thick(Teuchos::getNumericStringParameter(crackparam,"THICKNESS_ASSUMPTION"));
  if( thick == "plane_stress" )
  {
    kappa_ = ( 3.0 - poisson_ ) / ( 1.0 + poisson_ );
  }
  else if ( thick == "plane_strain" )
  {
    kappa_ = 3.0 - poisson_;
  }
  fac_ = young_ / ( 1.0 + poisson_ ) / ( 1.0 + kappa_ );

  // get propagation criterion
  propcrit_ = DRT::INPUT::IntegralValue<INPAR::CRACK::propagationCriterion>(crackparam, "CRACK_PROPAGATION_CRITERION");

  //---
  // Read critical stress intensity factor and critical J-integral from input file
  //---
  critical_K_I_ = crackparam.get<double>("CRITICAL_K1");
  critical_K_II_ = crackparam.get<double>("CRITICAL_K2");
  critical_J_ = crackparam.get<double>("CRITICAL_J");

  //---
  // Preparation for writing crack tip location into a file
  // makes sure already existing files are overwritten
  //---
  {
    const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());
    std::ostringstream pid_stream;
    pid_stream << ".crackTip" << segment_id_ <<".txt";

    filename_<< filebase<< pid_stream.str();

    std::ofstream f;
    if( myrank_ == 0 )
      f.open( filename_.str().c_str(), std::fstream::trunc );

    f.close();
  }

  //---
  // Preparation for writing J-integral
  // makes sure already existing files are overwritten
  //---
  {
    const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());
    std::ostringstream pid_stream;
    pid_stream << ".Jintegral" << segment_id_ <<".txt";

    filename_jint_<< filebase<< pid_stream.str();

    std::ofstream f;
    if( myrank_ == 0 )
      f.open( filename_jint_.str().c_str(), std::fstream::trunc );
    f.close();
  }
}

/*----------------------------------------------------------------------------------------------*
 * Perform operations related to propagation of this crack tip segment                  sudhakar 09/14
 *----------------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::propagateThisTip( Teuchos::RCP<const Epetra_Vector>& displace )
{
  // general geometric paramters
  oldnew_.clear();
  tip_phi_.clear();
  tip_mphi_.clear();
  newTip_.clear();
  tip_bc_disp_.clear();

  all_split_ele_.clear();

  isProp_ = false;

  // Crack has already propagated that the structure is completely split into two
  if( alreadySplit_ )
    return;

  newTip_.clear();

  // export "displacement" to column map
  disp_col_ = LINALG::CreateVector( *discret_->DofColMap(), true );
  LINALG::Export(*displace,*disp_col_);

  //WriteCrackSurfacePoints();

  //-------------
  // STEP 1 : Compute stress-intensity factors at crack tip
  //-------------

  if( propcrit_ == INPAR::CRACK::displacementCorrelation )
  {
    findStressIntensityFactor();
    // this arises at the very beginning of time integration
    // physically uncorrect
    // experienced this when using gen-alpha method  (sudhakar 12/13)
    if( K_I_ < 0.0 )
      return;
  }
  else if( propcrit_ == INPAR::CRACK::J_Integral )
  {
    Compute_J_integral();
  }
  else
    dserror("unknown crack propagation method\n");

  //-------------
  // STEP 2 : Check whether propagation criterion is satisfied
  //-------------
  isProp_ = DoCrackPropagate();
  if( not isProp_ )
    return;

  //-------------
  // STEP 3 : Decide crack propagation angle from stress-intensity factors
  //-------------
  decidePropagationAngle();

  //-------------
  // STEP 4 : Update crack information
  //-------------
  newTip_ = findNewCrackTip();

  //-------------
  // STEP 5 : Write the location of crack tip to a file
  //-------------
  WriteCrackTipLocation( newTip_ );
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
void DRT::CRACK::PropagateTip::findStressIntensityFactor()
{
  DRT::Node * tipnode = NULL;
  int lmaster = 0, gmaster = 0;

  if( discret_->HaveGlobalNode( tipnodes_[0] ) )
  {
    tipnode = discret_->gNode( tipnodes_[0] );
    if( tipnode->Owner() == myrank_ )
    {
      lmaster = myrank_;

      int attach1_id = 0, attach2_id = 0;

      DRT::Node * node1 = findNeighboringCrackNode( tipnode, false, NULL, attach1_id );
      DRT::Node * node2 = findNeighboringCrackNode( tipnode, true, node1, attach2_id );

      findNormal( tipnode, node1, node2 );

      // find neighbor points to calculate stress-intensity factor
      buildNeighborInfoTipNodes( node1, node2, tipnode, attach1_id, attach2_id );
    }
  }

  // making sure that master processor id is available at all processors
  comm_.SumAll( &lmaster, &gmaster, 1 );

  // normal and tangent vectors are computed only on master processor
  // broadcast them to all processors
  comm_.Broadcast( &normal_(0,0), 3, gmaster );
  comm_.Broadcast( &tangent_(0,0), 3, gmaster );

  LINALG::GatherAll( tip_phi_, comm_ );
  LINALG::GatherAll( tip_mphi_, comm_ );

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

        std::vector<double> disp_phi = DRT::CRACK::UTILS::getDisplacementNode( discret_, node_phi, disp_col_ );
        std::vector<double> disp_mphi = DRT::CRACK::UTILS::getDisplacementNode( discret_, node_mphi, disp_col_ );

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

  //TODO: can we simulate mixed mode and couple it to FSI-crack interaction
  //K_II_ = 0.0; //at the moment we simulate only K_1 mode crack??

  if( myrank_ == 0 )
  {
    std::cout<<"stress intensity factors = "<<K_I_<<"\t"<<K_II_<<"\n";
  }
}

/*-----------------------------------------------------------------------------------------*
 * Compute domain form of vector J-integral                                        sudhakar 11/14
 *-----------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::Compute_J_integral()
{
  DRT::Node * tipnode = NULL;
  int lmaster = 0, gmaster = 0;

  if( discret_->HaveGlobalNode( tipnodes_[0] ) )
  {
    tipnode = discret_->gNode( tipnodes_[0] );
    if( tipnode->Owner() == myrank_ )
    {
      lmaster = myrank_;

      int attach1_id = 0, attach2_id = 0;

      DRT::Node * node1 = findNeighboringCrackNode( tipnode, false, NULL, attach1_id );
      DRT::Node * node2 = findNeighboringCrackNode( tipnode, true, node1, attach2_id );

      findNormal( tipnode, node1, node2 );
    }
  }

  // making sure that master processor id is available at all processors
  comm_.SumAll( &lmaster, &gmaster, 1 );

  // normal and tangent vectors are computed only on master processor
  // broadcast them to all processors
  comm_.Broadcast( &normal_(0,0), 3, gmaster );
  comm_.Broadcast( &tangent_(0,0), 3, gmaster );
  comm_.Broadcast( &normal_ref_(0,0), 3, gmaster );
  comm_.Broadcast( &tangent_ref_(0,0), 3, gmaster );

  Jvector_.clear();
  J_Integral j_int( discret_, gausspts_stress_, gausspts_strain_, disp_col_, iforce_col_, cracknodes_, boun_nodes_, tipnodes_,
                    normal_, tangent_, segment_id_ );
  Jvector_ = j_int.compute_j_integral();

  if( myrank_ == 0 )
    std::cout<<"jvector = "<<Jvector_[0]<<" "<<Jvector_[1]<<"\n";

  // Write J-integral into a file
  std::ofstream f;
  if( myrank_ == 0 )
    f.open( filename_jint_.str().c_str(), std::fstream::ate | std::fstream::app );

  static int itno = 0;
  itno++;
  if( myrank_ == 0 )
  {
    std::ostringstream s;
    s<<itno<<"\t"<<Jvector_[0]<<"\t"<<Jvector_[1];
    f<<s.str()<<"\n";
  }

  if( myrank_ == 0 )
    f.close();
}

/*------------------------------------------------------------------------------------*
 * In order to calculate stress intensity factors at the crack tip, we        sudhakar 12/13
 * need appropriate connectivity of few nodes at the crack tip
 * VERY IMPORTANT : at the moment it works only for brick (hex) elements
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::buildNeighborInfoTipNodes( DRT::Node * node1, DRT::Node * node2, DRT::Node * tipnode, int& attach1_id, int& attach2_id )
{
#ifdef DEBUG
  if( node1->Id() == node2->Id() )
    dserror( "We should have got two different node ids\n" );
#endif

  // we are considering only the first node to construct normal
  // because we assume through-thickness crack with crack properties same at each thickness plane
  //TODO: while extending to 3D crack, neighbors has to be calculated at each crack tip node

  //DRT::Node * attach1 = findAttachedNode( node1, tipnode );
  //DRT::Node * attach2 = findAttachedNode( node2, tipnode );

  DRT::Node * attach1 = discret_->gNode( attach1_id );
  DRT::Node * attach2 = discret_->gNode( attach2_id );

  if( attach1 == NULL or attach2 == NULL)
    dserror("either one of the attached node is not found\n");

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
void DRT::CRACK::PropagateTip::findNormal( const DRT::Node * tipnode,
                                           const DRT::Node * surnode1,
                                           const DRT::Node * surnode2 )
{
  // we are considering only the first node to construct normal
  // because we assume through-thickness crack with crack properties same at each thickness plane
  //TODO: while extending to 3D crack, normal has to be calculated at each crack tip segment

  // In order to construct local coodinate system at crack tip, we consider one crack tip node
  // and find two neighboring nodes that are falling on the surfaces of the crack (surnode1 and surnode2 here)
  // REMEMBER : This procedure works only for through-thickness crack

  const double * tipcoord = tipnode->X();

  const double * surcoord1 = surnode1->X();
  const double * surcoord2 = surnode2->X();

  std::vector<double> tip_disp = DRT::CRACK::UTILS::getDisplacementNode( discret_, tipnode, disp_col_ );
  std::vector<double> sur1_disp = DRT::CRACK::UTILS::getDisplacementNode( discret_, surnode1, disp_col_ );
  std::vector<double> sur2_disp = DRT::CRACK::UTILS::getDisplacementNode( discret_, surnode2, disp_col_ );

  for( unsigned i=0; i<3; i++ )
  {
    tip_disp[i] += tipcoord[i];
    sur1_disp[i] += surcoord1[i];
    sur2_disp[i] += surcoord2[i];
  }

  double surcoord[3];
  for( int s = 0; s < 3; s++ )
    surcoord[s] = 0.5*( sur1_disp[s] + sur2_disp[s] );

  normal_(0,0) = tip_disp[0] - surcoord[0];
  normal_(1,0) = tip_disp[1] - surcoord[1];
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

  if( myrank_ == 0 )
  {
    std::cout<<"normal vector = "<<normal_(0,0)<<"\t"<<normal_(1,0)<<"\t"<<normal_(2,0)<<"\n";
    std::cout<<"tangent vector = "<<tangent_(0,0)<<"\t"<<tangent_(1,0)<<"\t"<<tangent_(2,0)<<"\n";
  }

  //---
  // Computation of normal and tangents in reference frame
  //---
  std::vector<double> surcoord_ref(3);
  for( int dim=0; dim<3; dim++ )
    surcoord_ref[dim] = 0.5*( surcoord1[dim] + surcoord2[dim] );

  normal_ref_(0,0) = tipcoord[0] - surcoord_ref[0];
  normal_ref_(1,0) = tipcoord[1] - surcoord_ref[1];
  normal_ref_(2,0) = 0.0;

  double fac_ref = sqrt( normal_ref_(0,0)*normal_ref_(0,0) + normal_ref_(1,0)*normal_ref_(1,0) + normal_ref_(2,0)*normal_ref_(2,0) );
  normal_ref_(0,0) = normal_ref_(0,0) / fac_ref;
  normal_ref_(1,0) = normal_ref_(1,0) / fac_ref;
  normal_ref_(2,0) = normal_ref_(2,0) / fac_ref;

  double theta_ref = atan2( normal_ref_(1,0), normal_ref_(0,0) ) - atan2( 0.0, 1.0 );

  // apply linear transformation to rotate (0,1) to angle theta
  tangent_ref_(0,0) = -1.0 * sin(theta_ref);
  tangent_ref_(1,0) = cos(theta_ref);
  tangent_ref_(2,0) = 0.0;
}

/*------------------------------------------------------------------------------------*
 * Find neighboring node associated with the first crack tip node
 * The neighboring node is located on the crack surface, not on tip             sudhakar 12/13
 * if bool second = true, it means we already have one node available and
 * searching for another node located at the same position as the first node
 *------------------------------------------------------------------------------------*/
DRT::Node * DRT::CRACK::PropagateTip::findNeighboringCrackNode( DRT::Node * tipnode,
                                                                bool second,
                                                                const DRT::Node * first_node,
                                                                int & attach_id )
{
  if( second and first_node == NULL )
    dserror( "While finding second neighbor, one should supply the first node\n" );

  int tipnodeid = tipnode->Id();
  DRT::Node * surnode = NULL;

  // get all elements that has this tipnode
  DRT::Element** ElementsPtr = tipnode->Elements();
  int NumElement = tipnode->NumElement();

  bool foundnode = false;

#if 0
  int surnodeid=0;
  const double * tipcoord = tipnode->X();
  const double * surcoord = NULL;
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
#endif

  for (int jEle = 0; jEle < NumElement; jEle++)
  {
    DRT::Element* Element = ElementsPtr[jEle];

    std::vector<Teuchos::RCP<DRT::Element > >allsurf = Element->Surfaces();
    for( std::vector<Teuchos::RCP<DRT::Element> >::iterator it = allsurf.begin(); it != allsurf.end(); it++ )
    {
      Teuchos::RCP<DRT::Element> surfele = *it;

      // if surface does not have this element, do not consider it
      if( not DRT::CRACK::UTILS::ElementHasThisNodeId( surfele, tipnodeid ) )
        continue;

      const int* surnodes = surfele->NodeIds();

      int conn1_index = 0, conn2_index = 0;
      DRT::CRACK::UTILS::getConnectedNodeIdIndex( surfele, tipnodeid, conn1_index, conn2_index );

      int conn1_id = surnodes[conn1_index];
      int conn2_id = surnodes[conn2_index];

      // if already found first node, no  connecting node of this, should be a first node
      if( second )
      {
        int firstid = first_node->Id();
        if( conn1_id == firstid or conn2_id == firstid )
          continue;
      }

      // we have got wrong side if one of the connecting nodes falls in tip node
      if( std::find( tipnodes_.begin(), tipnodes_.end(), conn1_id ) != tipnodes_.end() or
          std::find( tipnodes_.begin(), tipnodes_.end(), conn2_id ) != tipnodes_.end() )
        continue;

      int surnodeid = 0;

      if( std::find( oldTipnodes_.begin(), oldTipnodes_.end(), conn1_id ) != oldTipnodes_.end() )
      {
        surnodeid = conn1_id;
        foundnode = true;
      }
      else if( std::find( oldTipnodes_.begin(), oldTipnodes_.end(), conn2_id ) != oldTipnodes_.end() )
      {
        surnodeid = conn2_id;
        foundnode = true;
      }

      if( foundnode )
      {
        if( surnodeid == conn1_id )
        {
          surnode = discret_->gNode( conn1_id );
          attach_id = conn2_id;
          break;
        }
        else if( surnodeid == conn2_id )
        {
          surnode = discret_->gNode( conn2_id );
          attach_id = conn1_id;
          break;
        }
        else
          dserror("Surnodeid is wrong\n");
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
DRT::Node * DRT::CRACK::PropagateTip::findAttachedNode( DRT::Node * neigh, DRT::Node * tipnode )
{
  DRT::Node * attach = NULL;
  //const double * neighcord = neigh->X();

  int tipid = tipnode->Id();

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
          if( id == tipid )
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

  Teuchos::RCP<DRT::Element> surf = getSurfaceThisPlane( Surfaces, tipid );

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

/*---------------------------------------------------------------------------------------------------------*
 * Given as input of all surfaces of an element, we find the surface one among these,              sudhakar 09/14
 * that contains only one crack tip node which is given as input "tipid
 *---------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::CRACK::PropagateTip::getSurfaceThisPlane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces,
                                                                          int tipid )
{
  bool foundsur = false;
  Teuchos::RCP< DRT::Element > reqd_surf = Teuchos::null;

  for( unsigned surno = 0; surno<Surfaces.size(); surno++ )
  {
    Teuchos::RCP<DRT::Element> surele = Surfaces[surno];
    const int * nodeids = surele->NodeIds();

    // make sure given tip node is available on this surface
    bool found_this_tip = false;
    for( int nodno = 0; nodno < surele->NumNode(); nodno++ )
    {
      if( nodeids[nodno] == tipid )
      {
        found_this_tip = true;
        break;
      }
    }

    if( not found_this_tip )
      continue;

    // we need a surface that contains given crack tip --> we got them
    // surface that is on this plane should not contain any other tip nodes
    bool found_other_tip = false;
    for( int nodno = 0; nodno < surele->NumNode(); nodno++ )
    {
      int thisid = nodeids[nodno];
      if( thisid == tipid )
        continue;

      if( std::find( tipnodes_.begin(), tipnodes_.end(), thisid ) != tipnodes_.end() )
      {
        found_other_tip = true;
        break;
      }
    }

    if( not found_other_tip )
    {
      reqd_surf = surele;
      foundsur = true;
      break;
    }

  }

  if( not foundsur )
    dserror( "Required surface is not found\n" );

  return reqd_surf;
}

/*------------------------------------------------------------------------------------*
 * Returns true if the crack propagation criterion is satisfied               sudhakar 12/13
 *------------------------------------------------------------------------------------*/
bool DRT::CRACK::PropagateTip::DoCrackPropagate()
{
  if ( DRT::Problem::Instance()->ProblemType() == prb_fsi_crack )
  {
    if( propcrit_ == INPAR::CRACK::displacementCorrelation )
    {
      if( K_I_ > critical_K_I_ )
        return true;
      return false;
    }
    else if( propcrit_ == INPAR::CRACK::J_Integral )
    {
      if( Jvector_[0] > critical_J_ )
        return true;
      return false;
    }
    else
      dserror("unknown propagation criterion\n");
  }

  else
  {
    if( propcrit_ == INPAR::CRACK::displacementCorrelation )
    {
      double check = pow( (K_I_ / critical_K_I_), 2.0 ) + pow( (K_II_ / critical_K_II_), 2.0 );
      if( check < 0.99999999999999999 )
        return false;
      return true;
    }
    else if( propcrit_ == INPAR::CRACK::J_Integral )
    {
      double jtot = sqrt( pow(Jvector_[0],2) + pow(Jvector_[1],2) );
      if( jtot > critical_J_ )
        return true;
      return false;
    }
    else
      dserror("unknown propagation criterion\n");
  }

  dserror("should not have reached here");
  return true; // just for warning-free compilation
}

/*------------------------------------------------------------------------------------*
 * Decide crack propagation angle from stress-intensity factors               sudhakar 01/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::decidePropagationAngle()
{
  if( propcrit_ == INPAR::CRACK::displacementCorrelation )
  {
    double deno = K_I_ + sqrt( K_I_ * K_I_ + 8.0 * K_II_ * K_II_ );
    propAngle_ = 2.0*atan( -2.0 * K_II_ / deno );
  }
  else if( propcrit_ == INPAR::CRACK::J_Integral )
  {
    propAngle_ = atan2( Jvector_[1], Jvector_[0] );
  }
  else
    dserror("unknown crack propagation method\n");

  //propAngle_ = 0.0;
  /*********************************************************/
  /*static int k = 1;
  if( k==1 )
  {
    std::cout<<"------WARNING: crack propagation angle set----------\n";
    //propAngle_ = atan2( normal_(1,0), normal_(0,0) );
    propAngle_ = 0.0;
  }
  k++;*/
  /*********************************************************/

  // Need to confirm whether this is generalized for all propagation criterion
  // read the paper "Numerical modelling of crack propagation: automatic remeshing and comparison of
  // different criteria" by Bouchard, Chastel in CMAME 2003 -->Maximum circumferential stress criterion
  if( propAngle_ * 180 / PI_ > 70.54 )
    dserror("Predicted crack propagation angle = %f but the limiting prop angle = 70.54",propAngle_ * 180 / PI_);

  // the obtained propagation angle is w.r.t normal direction
  // in order to get absolute propagation angle, add the normal angle to this
  double norm_ang = atan2( normal_(1,0), normal_(0,0) );
  DRT::CRACK::UTILS::convertAngleTo_02PI_range( norm_ang );

  propAngle_ += norm_ang;

  if( (fabs(propAngle_) < ANGLE_TOL_ZERO) or (fabs(propAngle_-2*PI_) < ANGLE_TOL_ZERO) )
    propAngle_ = 0.0;
  else if( propAngle_ < 0.0 )
    propAngle_ = 2 * PI_ + propAngle_;
  else if( propAngle_ > 2*PI_ )
    propAngle_ = propAngle_ - 2*PI_;

  if( myrank_ == 0 )
    std::cout<<"propagation angle = "<<propAngle_ * 180.0 / PI_ <<"deg\n";
}

/*---------------------------------------------------------------------------------------------------*
 * Procedure to find new crack tip nodes. In addition, this also decides
 * whether to split the HEX into two, and whether ALE step is                               sudhakar 08/14
 * mandatory
 *---------------------------------------------------------------------------------------------------*/
std::vector<int> DRT::CRACK::PropagateTip::findNewCrackTip()
{
  std::map<int,int> oldnew_tip;

  //TODO: In the first step of crack propagation, we should not reach the surface
  // that is located on the bundary of the domain
  // Implement a procedure to take care of this

  for( unsigned tip = 0; tip < tipnodes_.size(); tip++ )
  {
    int tipid = tipnodes_[tip];
    int gnewtip = 0, lnewtip = 0;

    int lis_split=0, gis_split=0;       // is any hex is split at this crack tip
    int lspl_node_id=0, gspl_node_id=0; // node other than the tip at which it splits
    int lspl_ele_id=0, gspl_ele_id=0;   // Id of HEX element that is split

    if( discret_->HaveGlobalNode( tipid ) )
    {
      DRT::Node* tipnode = discret_->gNode( tipid );
      if( tipnode->Owner() == myrank_ )
      {
        const double * tipcord = tipnode->X();

        std::vector<double> disp_tip = DRT::CRACK::UTILS::getDisplacementNode( discret_, tipnode, disp_col_ );
        for( unsigned dim=0; dim<3; dim++ )
          disp_tip[dim] += tipcord[dim];

        // construct a unit vector in the direction of propagation angle from crack tip
        std::vector<double> prop_vec(3);
        prop_vec[0] = disp_tip[0] + LENGTH_PROP_VECTOR * cos( propAngle_ );
        prop_vec[1] = disp_tip[1] + LENGTH_PROP_VECTOR * sin( propAngle_ );
        prop_vec[2] = disp_tip[2];

        // This is same as the propagation angle but in the range [-pi,pi]
        double projangle = propAngle_;
        DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( projangle );

        DRT::Element** ElementsPtr = tipnode->Elements();
        int NumElement = tipnode->NumElement();

        bool found_edge = false;

        for( int eleno = 0; eleno < NumElement; eleno++ )
        {
          DRT::Element* Element = ElementsPtr[eleno];
          std::vector< Teuchos::RCP< DRT::Element > > Surfaces = Element->Surfaces();

          Teuchos::RCP<DRT::Element> surf = getSurfaceThisPlane( Surfaces, tipid );

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

          /* FInd possible node ids. This means through which ids the crack can propagate
             For Tri surface there is no problem
             But if the surface is a Quad, then we must avoid the diagonal node as possible id

                                                                              === crack surface
                        o---------------                    o                     * crack tip
                        |               |                   | \                   o possible id nodes
                        |               |                   |  \
                        |               |                   |   \
                        |               |                   |    \
                        |               |                   |     \
                   =====*---------------o                ===*------o

          */

          DRT::Node * possible1 = searchnodes[(tip_index+1)%numnodes];
          int temp_index = 0;
          if( tip_index == 0 )
            temp_index = numnodes-1;
          else
            temp_index = tip_index-1;
          DRT::Node * possible2 = searchnodes[temp_index];

          std::vector<Teuchos::RCP<DRT::Element> > ele_lines = surf->Lines();

          if( ele_lines.size() == 0 )
            dserror("Surface element does not contain any lines");

          for( std::vector<Teuchos::RCP<DRT::Element> >::iterator linit = ele_lines.begin();
                                                                  linit != ele_lines.end(); linit++ )
          {
            Teuchos::RCP<DRT::Element> linele = *linit;

            DRT::Node ** searchnodes = linele->Nodes();
            DRT::Node* node1 = searchnodes[0];
            DRT::Node* node2 = searchnodes[1];

            // crack cannot propagate through the edge that has the crack node
            if( node1 == tipnode or node2 == tipnode )
              continue;

            /*if( ( std::find( cracknodes_.begin(), cracknodes_.end(), node1->Id() ) != cracknodes_.end() ) or
                ( std::find( cracknodes_.begin(), cracknodes_.end(), node2->Id() ) != cracknodes_.end() )  )
              continue;*/

            const double * node1_cord = node1->X();
            const double * node2_cord = node2->X();

            std::vector<double> disp1 = DRT::CRACK::UTILS::getDisplacementNode( discret_, node1, disp_col_ );
            std::vector<double> disp2 = DRT::CRACK::UTILS::getDisplacementNode( discret_, node2, disp_col_ );

            for( unsigned dim=0; dim<3; dim++ )
            {
              disp1[dim] += node1_cord[dim];
              disp2[dim] += node2_cord[dim];
            }

            int num = 0;

            double angle1 = atan2( disp1[1] - disp_tip[1], disp1[0] - disp_tip[0] ) - projangle;
            double angle2 = atan2( disp2[1] - disp_tip[1], disp2[0] - disp_tip[0] ) - projangle;

            DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( angle1 );
            DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( angle2 );

            //TODO: should this condition exist?
            if( fabs(angle1) > 0.5*PI_ and fabs(angle2) > 0.5*PI_ )
              continue;

            // it may be possible that the crack can pass through the edge itself
            if( ( fabs( angle1 ) < min_angle_tol_ ) and ( node1->Id() == possible1->Id() or node1->Id() == possible2->Id() ) )
            {
              found_edge = true;
              lnewtip = node1->Id();

              // Now we are checking this only once because we simulate pseudo-3D crack propagation
              // For real 3D this has to  be set for all crack nodes
              break;
            }

            if( ( fabs( angle2 ) < min_angle_tol_ ) and ( node2->Id() == possible1->Id() or node2->Id() == possible2->Id() ) )
            {
              found_edge = true;
              lnewtip = node2->Id();

              // Now we are checking this only once because we simulate pseudo-3D crack propagation
              // For real 3D this has to  be set for all crack nodes
              break;
            }

            if( fabs( angle1 ) < TOL_DIAG_NODE and ( node1->Id() != possible1->Id() and node1->Id() != possible2->Id() ))
            {
              if( (angle1 < 0.0 and angle2 > 0.0 ) or
                (angle1 > 0.0 and angle2 < 0.0 ) )
              {
                found_edge = true;
                lnewtip = node1->Id();

                lis_split = 1;
                lspl_node_id = node1->Id();
                lspl_ele_id = Element->Id();
                //SplitEleData( Element->Id(), node1->Id() );

                //TODO: Check this
                //if( is_A_BoundaryNode( lnewtip ) )
                //  break;

                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                double ratio = fabs(angle1 / tot_ang);
                for( unsigned dim=0; dim<3; dim++ )
                  disp_bc[dim] = ( ratio * disp2[dim] + (1.0-ratio) * disp1[dim]) - disp1[dim];//node1_cord[dim];
                tip_bc_disp_[node1->Id()] = disp_bc;

                //all_split_ele_ids_.insert( Element->Id() );

                break;
              }
            }

            if( fabs( angle2 ) < TOL_DIAG_NODE and ( node2->Id() != possible1->Id() and node2->Id() != possible2->Id() ))
            {
              if( (angle1 < 0.0 and angle2 > 0.0 ) or
                (angle1 > 0.0 and angle2 < 0.0 ) )
              {
                found_edge = true;
                lnewtip = node2->Id();

                lis_split = 1;
                lspl_node_id = node2->Id();
                lspl_ele_id = Element->Id();
                //SplitEleData( Element->Id(), node2->Id() );

                //TODO: Check this
                //if( is_A_BoundaryNode( lnewtip ) )
                //  break;

                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                double ratio = fabs(angle2 / tot_ang);
                for( unsigned dim=0; dim<3; dim++ )
                  disp_bc[dim] = ( ratio * disp1[dim] + (1.0-ratio) * disp2[dim] ) - disp2[dim];//node2_cord[dim];
                tip_bc_disp_[node2->Id()] = disp_bc;

                //all_split_ele_ids_.insert( Element->Id() );

                break;
              }
            }

            if( (angle1 < 0.0 and angle2 > 0.0 ) or
                (angle1 > 0.0 and angle2 < 0.0 ) )
            {
              std::string whichnode = "";

              if( (node1->Id() == possible1->Id() and node2->Id() == possible2->Id()) or
                  (node1->Id() == possible2->Id() and node2->Id() == possible1->Id()) )
              {
                if( fabs(angle1) < fabs(angle2) )
                  whichnode = "node1";
                else
                  whichnode = "node2";
              }

              else if( node1->Id() == possible1->Id() or node1->Id() == possible2->Id() )
                whichnode = "node1";
              else if( node2->Id() == possible1->Id() or node2->Id() == possible2->Id() )
                whichnode = "node2";
              else
                dserror(" one of the nodes must be a possible node id ");

              if( whichnode == "node1" )
              {
                if( angle1 * 180 / PI_ > 90.0 )
                  continue;
                num = 1;
                lnewtip = node1->Id();
              }
              else if( whichnode == "node2" )
              {
                if( angle2 * 180 / PI_ > 90.0 )
                  continue;
                num = 2;
                lnewtip = node2->Id();
              }

              found_edge = true;

              //TODO: Check this
              //if( is_A_BoundaryNode( lnewtip ) )
              //  break;

              if( found_edge )
              {
                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                if( num == 1 )
                {
                  double ratio = fabs(angle1 / tot_ang);
                  for( unsigned dim=0; dim<3; dim++ )
                    disp_bc[dim] = ( ratio * disp2[dim] + (1.0-ratio) * disp1[dim]) - disp1[dim];//node1_cord[dim];
                  tip_bc_disp_[node1->Id()] = disp_bc;
                }
                else if( num == 2 )
                {
                  double ratio = fabs(angle2 / tot_ang);
                  for( unsigned dim=0; dim<3; dim++ )
                    disp_bc[dim] = ( ratio * disp1[dim] + (1.0-ratio) * disp2[dim] ) - disp2[dim];//node2_cord[dim];
                  tip_bc_disp_[node2->Id()] = disp_bc;
                }
                else
                  dserror("the correct index not found");
              }

              break;
            }
          }

          if( found_edge )
          {
            break;
          }
        }

        if( not found_edge )
          dserror( "not found the new crack tip for nodeid = %d", tipid );
      }
    }

    comm_.SumAll( &lnewtip, &gnewtip, 1 );
    oldnew_tip[tipid] = gnewtip;

    // Generate data for HEX elements that has to be split
    // It is mandatory to have this data redundant on all processors to make
    // sure that new elements get identical Ids on all proc.
    // Also make sure all proc that has this element perform the splitting
    comm_.SumAll( &lis_split, &gis_split, 1 );
    if( gis_split > 0 )
    {
      comm_.SumAll( &lspl_ele_id, &gspl_ele_id, 1 );
      comm_.SumAll( &lspl_node_id, &gspl_node_id, 1 );
      SplitEleData( gspl_ele_id, gspl_node_id );
    }
  }


  LINALG::GatherAll( tip_bc_disp_, comm_ );

  if(myrank_==0)
  {
    std::cout<<"----------------printing ALE displacement boundary condition---------------\n";
    for(std::map<int,std::vector<double> >::iterator it = tip_bc_disp_.begin(); it != tip_bc_disp_.end(); it++)
      std::cout<<"node id = "<<it->first<<" ALE disp = "<<(it->second)[0]<<" "<<(it->second)[1]<<" "<<(it->second)[2]<<"\n";
    std::cout<<"---------------------------------------------------------------------------\n";
  }


  if( tipnodes_.size() != oldnew_tip.size() )
    dserror("for each node, we should have a new tip node\n");

  std::vector<int> newTip;
  for( std::map<int,int>::iterator it = oldnew_tip.begin(); it != oldnew_tip.end(); it++ )
  {
    newTip.push_back( it->second );
  }

  return newTip;
}

/*---------------------------------------------------------------------------------------------------------*
 * Store all data that is necessary to split a HEX into WEDGE elements                           sudhakar 07/14
 *---------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::SplitEleData( std::map<int,Teuchos::RCP<DRT::Element> >& zele,
                                               int add_id )
{
  bool hasnode = false;
  int ori_ele_id = 0;

  // Get original element id that contains this node
  // Remember : zele --> contains surface elements that has particular z-coordinate
  for( std::map<int,Teuchos::RCP<DRT::Element> >::iterator it = zele.begin(); it != zele.end(); it++ )
  {
    Teuchos::RCP<DRT::Element> ele = it->second;
    if( DRT::CRACK::UTILS::ElementHasThisNodeId( ele, add_id ) )
    {
      hasnode = true;
      ori_ele_id = it->first;
    }
  }

  if( not hasnode )
    dserror("The z-plane elements that are passed do not contain the given node id\n");

  SplitEleData( ori_ele_id, add_id );
}

/*---------------------------------------------------------------------------------------------------------*
 * Store all data that is necessary to split a HEX into WEDGE elements                           sudhakar 07/14
 *---------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::SplitEleData( int ele_id, int node_id )
{
  bool add_data = false;

  // check if already a structure created for this element
  // if so, just add this node id to it
  if( all_split_ele_.find( ele_id ) != all_split_ele_.end() )
  {
    splitThisEle_& spl = all_split_ele_[ele_id];
    if( spl.element_id_ == ele_id )
    {
      (spl.node_ids_ele_).push_back( node_id );
      add_data = true;
    }
  }

  // This element has not dealt with before, hence create a new entry
  if( not add_data )
  {
    splitThisEle_ spl;
    spl.element_id_ = ele_id;
    (spl.node_ids_ele_).push_back( node_id );
    all_split_ele_[ele_id] = spl;
  }
}

/*------------------------------------------------------------------------------------*
 * Write crack tip location to an output file to trace crack path              sudhakar 07/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::WriteCrackTipLocation( const std::vector<int>& newTip )
{
  std::ofstream f;
  if( myrank_ == 0 )
    f.open( filename_.str().c_str(), std::fstream::ate | std::fstream::app );

  int lmaster = 0, gmaster = 0;
  std::vector<double> disp(3,0.0);
  int nodeid = newTip[0];

  if( not (tip_bc_disp_.size() == 0) )
  {
    std::map<int, std::vector<double> >::iterator it = tip_bc_disp_.begin();
    nodeid = it->first;
    disp=it->second;
  }

  if( discret_->HaveGlobalNode( nodeid ) )
  {
    DRT::Node * node = discret_->gNode( nodeid );
    if( node->Owner() == myrank_ )
    {
      const double* x = discret_->gNode( nodeid )->X();
      for( unsigned dim = 0; dim < 3; dim++ )
        disp[dim] += x[dim];

      lmaster = myrank_;
    }
  }

  // making sure that master processor id is available at all processors
  comm_.SumAll( &lmaster, &gmaster, 1 );

  // normal and tangent vectors are computed only on master processor
  // broadcast them to all processors
  comm_.Broadcast( &disp[0], 3, gmaster );

  if( myrank_ == 0 )
  {
    std::ostringstream s;
    s<<disp[0]<<"\t"<<disp[1];
    f<<s.str()<<"\n";
  }

  if( myrank_ == 0 )
    f.close();
}

/*---------------------------------------------------------------------------------------------------------*
 * Loop over all HEX elements to be split, perform splitting operations                          sudhakar 08/14
 * and add them to discretization
 *
 *                   ---------------*
 *                 /|              / |                     o tip nodes
 *                / |             /  |                     * split_nodes
 *               |---------------*   |                   === crack surface
 *               |  |            |   |
 *               |  |            |   |
 *               |  |            |   |
 *               |  o------------|--/
 *               | /             | /
 *               |/              |/
 *      =========o---------------
 *
 *  only split_nodes are stored into all_split_ele_ data structure
 *  we can get tip nodes from the element details itself
 *
 *  HEX is split into two WEDGEs by cutting along the diagonal given by these nodes
 *---------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::split_All_HEX_Elements( int &new_ele_id )
{
  if( all_split_ele_.size() == 0 )
    return;

  SplitHexIntoTwoWedges split( discret_, material_id_ );

  for( std::map<int,splitThisEle_>::iterator it = all_split_ele_.begin(); it != all_split_ele_.end(); it++ )
  {
    splitThisEle_ spl = it->second;
    int eleid = spl.element_id_;
    std::vector<int> split_nodes = spl.node_ids_ele_;

    // There should be exactly 2 nodes along which to split this ele
    if( split_nodes.size() != 2 )
      dserror( "There should only be two split nodes for HEX elements" );

    if( discret_->HaveGlobalElement( eleid ) )
    {
      DRT::Element* ele = discret_->gElement( eleid );
      const int* elenodes = ele->NodeIds();

      // get tip nodes from element details
      std::vector<int> tip_ids;
      for( int i=0; i<ele->NumNode(); i++ )
      {
        if( std::find( tipnodes_.begin(), tipnodes_.end(), elenodes[i] ) != tipnodes_.end() )
          tip_ids.push_back( elenodes[i] );
      }

      if( tip_ids.size() != 2 )
        dserror( "Assuming HEX element, there should be only 2 tip nodes\n" );

      split.DoAllSplittingOperations( ele, tip_ids, split_nodes, ele->Id(), new_ele_id );
    }

    // it is necessary to have same value on all proc to have idential Ids of same element
    // on all proc
    new_ele_id++;

    discret_->FillComplete( false, true, false );
  }


  // elerowmap, elecolmap, noderowmap and nodecolmap are in some kind of relation in BACI
  // After adding a new element, this is disturbed
  // To restore them, the following operations are performed.
  // Othewise, we get problems when creating ALE discretization, where this
  // compatible check is performed
  const Epetra_Map * noderowmap = discret_->NodeRowMap();
  const Epetra_Map * nodecolmap = discret_->NodeColMap();

  Teuchos::RCP< Epetra_Map > elerowmap;
  Teuchos::RCP< Epetra_Map > elecolmap;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  discret_->BuildElementRowColumn(*noderowmap, *nodecolmap, elerowmap, elecolmap);
  // we can now export elements to resonable row element distribution
  discret_->ExportRowElements(*elerowmap);
  // export to the column map / create ghosting of elements
  discret_->ExportColumnElements(*elecolmap);
  discret_->FillComplete();
}

/*------------------------------------------------------------------------------------*
 * Make all the modifications in the discretization related to crack propagation
 * Duplicate crack tip nodes, and modify element connectivity               sudhakar 12/13
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::updateCrack( int& start_new_ele_id, int & start_new_node_id )
{
  // split a hex element into two wedge elements when necessary
  split_All_HEX_Elements( start_new_ele_id );

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

  // int totalNodes = discret_->NumGlobalNodes();

  std::vector<double> lmtAngle = getLimitAngles( newTip_ );

  //std::cout<<"limit angle = "<<lmtAngle[0]*180.0/PI_<<"\t"<<lmtAngle[1]*180.0/PI_<<"\n";

  // map of element ids to be modified with the new node
  // "key" contains the element id, and "value" is dummy here
  // map is used to make sure that one element is stored only once
  std::map<int, int> delEle;

  for( unsigned num = 0; num < tipnodes_.size(); num++ )
  {
    if( discret_->HaveGlobalNode( tipnodes_[num] ) )
    {
      DRT::Node * tipnode = discret_->gNode( tipnodes_[num] );

      Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( start_new_node_id, tipnode->X(), tipnode->Owner() ) );
      oldnew_[tipnode->Id()] = dupnode->Id();

      DRT::Element** ElementsPtr = tipnode->Elements();
      int NumElement = tipnode->NumElement();

      for (int jEle = 0; jEle < NumElement; jEle++)
      {
        DRT::Element* Element = ElementsPtr[jEle];
        if( toReplaceNode( Element, tipnode, lmtAngle ) )
        {
          delEle[ Element->Id() ] = 0;
        }
      }
      discret_->AddNode( dupnode );
    }

    start_new_node_id++;
  }

  LINALG::GatherAll( oldnew_, comm_ );
  LINALG::GatherAll( delEle, comm_ );

  if( delEle.size() == 0 )
    dserror("propagated crack. but no elements found to attach new nodes. This leads to incorrect discretization\n");

  DRT::CRACK::UTILS::ModifyElementConnectivity( discret_, delEle, oldnew_ );

  // add newly generated nodes to appropriate conditions
  DRT::CRACK::UTILS::AddConditions( discret_, oldnew_ );

  discret_->FillComplete();

  // update crack tip nodes and add new crack tip nodes to cracknodes_
  tipnodes_.clear();
  tipnodes_.resize( newTip_.size() );
  std::copy( newTip_.begin(), newTip_.end(), tipnodes_.begin() );
  std::copy( newTip_.begin(), newTip_.end(), std::inserter( cracknodes_, cracknodes_.end() ) );


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
std::vector<double> DRT::CRACK::PropagateTip::getLimitAngles( const std::vector<int>& newTip )
{
  std::vector<double> ang;

  //the two angles are
  //1. negative of angle formed by normal
  //2. angle between new crack tip nodes, and old tip nodes (shuld be crack propagation angle if we move our nodes
  // to accommodate crack propagation)
  double norm_ang = atan2( normal_(1,0), normal_(0,0) ) + PI_;
  if( norm_ang > 2*PI_ )
    norm_ang = norm_ang - 2*PI_;

  ang.push_back( norm_ang );

  double tempAng = 0.0;

  if( discret_->HaveGlobalNode( tipnodes_[0] ) )
  {
    DRT::Node * tipnode = discret_->gNode( tipnodes_[0] );
    if( tipnode->Owner() == myrank_ )
    {
      const double * tipco = tipnode->X();

      DRT::Node * newnode = discret_->gNode( newTip[0] );
      const double * newco = newnode->X();

      tempAng = atan2( (newco[1]-tipco[1]), (newco[0]-tipco[0]) );

      if( fabs(tempAng) < 1e-12 )
        tempAng = 0.0;
      else if( tempAng < 0.0 )
        tempAng = 2 * PI_ + tempAng;
    }
  }

  double secAng = 0.0;
  comm_.SumAll( &tempAng, &secAng, 1 );

#if 1 // this shud be the case if we move crack nodes
  ang.push_back( propAngle_ );
#else
  ang.push_back( secAng );
#endif

  std::sort( ang.begin(), ang.end() );

  return ang;
}

/*------------------------------------------------------------------------------------*
 * Returns true if the criteron to replace the tipnode with a new            sudhakar 01/14
 * duplicate node is satisfied
 *------------------------------------------------------------------------------------*/
bool DRT::CRACK::PropagateTip::toReplaceNode( DRT::Element * ele, DRT::Node * tip, const std::vector<double>& lmtAngle )
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
  if( fabs(theta) < ANGLE_TOL_ZERO )
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
 *  The tip of crack has reached the boundary of the domain                      sudhakar 12/13
 *  This means two parts of structure are joined only along a line. In this case
 *  we remove the connection between two portions, and allow the structure to
 *  separate into two distinct parts
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::CheckCompleteSplit( int & start_new_node_id )
{

  // Count no of tip nodes that are fallling on the boundary
  unsigned nodes_on_boundary = 0;
  for( std::vector<int>::iterator itip = tipnodes_.begin(); itip != tipnodes_.end(); itip++ )
  {
    int tip_id = *itip;
    if( boun_nodes_.find( tip_id ) != boun_nodes_.end() )
      nodes_on_boundary++;
  }

  // no nodes are falling on the boundary. Nothing to do
  if( nodes_on_boundary == 0 )
    return;

  // This corresponds to 3D crack, or multiple cracks in 2D or 3D
  else if( nodes_on_boundary != tipnodes_.size() )
    dserror("Only few tip nodes are falling on the boundary. This means it is not through-thickness crack."
                  "Current crack implementation does not support this\n");

  std::cout<<"================splitting the body into two================\n";

  // All tip nodes on the boundary. Split the structure into two parts by creating new nodes
  // and modifying the connectivity

  /*                 |\                                 |\
                     | \                                | \
                     |  \                               |  \            * Tip node
                     |   \                              |   \           o new node generated to split the body
                     \   |                               \   |
                      \  |                                \  |
                       \ |                                 \ |
                        \|                                  \|
                         *       ===============>            *
                        /|
                       / |                                   o
                      /  |                                  /|
                     /   |                                 / |
                     |   /                                /  |
                     |  /                                /   |
                     | /                                 |   /
                     |/                                  |  /
                                                         | /
                                                         |/

                                                         */

  // calculation of the direction of crack propagation vector
  // This is different than the crack propagation angle at the moment because
  // we propaagate crack only along the crack edges, and not on with the correct propagation angle

  // map of element ids to be modified with the new node
  // "key" contains the element id, and "value" is dummy here
  // map is used to make sure that one element is stored only once
  std::map<int, int> delEle;

  // stores the old new node map only for split elements
  // in the end they are copied to oldnew_
  std::map<int,int> oldnew_split;

  for( unsigned num = 0; num < tipnodes_.size(); num++ )
  {
    int tip_id = tipnodes_[num];
    std::map<int,int>::iterator iter_map = oldnew_.begin();
    std::advance( iter_map, num );
    int old_tipid1 = iter_map->first;
    int old_tipid2 = iter_map->second;

    if( discret_->HaveGlobalNode( tip_id ) )
    {
      DRT::Node * tipnode = discret_->gNode( tip_id );
      const double * tipco = tipnode->X();

      if( not discret_->HaveGlobalNode( old_tipid1 ) )
      {
        dserror("old tip id not found on this processor");
        //continue;
      }

      DRT::Node* oldnode = discret_->gNode( old_tipid1 );
      const double * oldco = oldnode->X();
      double main_angle = DRT::CRACK::UTILS::FindAngle( tipco, oldco );

      // Get the elements that share the present and old tip node
      // Only these element's connectivity needs to be modified
      DRT::Element** ElementsPtr = tipnode->Elements();
      int NumElement = tipnode->NumElement();
      for (int jEle = 0; jEle < NumElement; jEle++)
      {
        DRT::Element* Element = ElementsPtr[jEle];
        if( DRT::CRACK::UTILS::ElementHasThisNodeId( Element, old_tipid1 ) or
            DRT::CRACK::UTILS::ElementHasThisNodeId( Element, old_tipid2 ) )
        {
          std::vector<double> cenco = Element->ElementCenterRefeCoords();
          double cen[3];
          for( unsigned ice=0; ice<3; ice++ )
            cen[ice] = cenco[ice];

          double angdiff = main_angle - DRT::CRACK::UTILS::FindAngle( tipco, cen );
          DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( angdiff );

          if( fabs(angdiff) < ANGLE_TOL_ZERO )
            dserror("vector with element centre coordinates coincides with vector with old tip node? Impossible!");

          if( angdiff > 0.0 )
          {
            delEle[ Element->Id() ] = 0;
          }
        }
      }
      Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( start_new_node_id, tipnode->X(), tipnode->Owner() ) );
      oldnew_split[tipnode->Id()] = dupnode->Id();
      discret_->AddNode( dupnode );
    }
    start_new_node_id++;
  }

  LINALG::GatherAll( oldnew_split, comm_ );
  LINALG::GatherAll( delEle, comm_ );

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
        std::map<int, int>::iterator delnod = oldnew_split.find( oldnodes[i] );
        if( delnod != oldnew_split.end() )
        {
          del = true;
          newnodes[i] = delnod->second;
        }
        else
          newnodes[i] = oldnodes[i];
      }

      if( not del )
        dserror("This element should have atleast one replaceable node\n");

      if( newnodes.size() != static_cast<std::size_t>(ele->NumNode()) )
        dserror("Check the number of new nodes\n");

      //-----------------
      // modifying the nodes of an element is easy
      // just modify the node ids of the element
      // when fillcomplete is called the corresponding nodes will be
      //  set through DRT::Element::BuildNodalPointers()
      //-----------------
      ele->SetNodeIds(newnodes.size(), &newnodes[0]);
    }
  }

  // copy oldnew for split elements into oldnew_
  oldnew_.insert(oldnew_split.begin(), oldnew_split.end());

  DRT::CRACK::UTILS::AddConditions( discret_, oldnew_ );

  discret_->FillComplete();


  /************************************************************/ //-------------
  // This is to add a prescribed displacement when the body separates completely into two
  //TODO: check to remove the lines marked as //-----------------
  /*std::map<int,std::vector<double> > disppl;//-----------------
  std::vector<double> disp(3,0.0);//-----------------
  disp[0] = -0.05;//------------------

  for( unsigned num = 0; num < tipnodes_.size(); num++ )
  {
    int tip_id = tipnodes_[num];
    if( discret_->HaveGlobalNode( tip_id ) )
    {
      DRT::Node * tipnode = discret_->gNode( tip_id );
      std::vector<double> disp_temp = DRT::CRACK::UTILS::getDisplacementNode( discret_, tipnode, disp_col_ );//-------------

      for( unsigned dispp = 0; dispp < 3; dispp++ )//-------------
        disp_temp[dispp] += disp[dispp];//-------------

      int dupid = 0;
      //std::map<int,int>::iterator itmap = std::find( oldnew_.begin(), oldnew_.end(), tip_id );
      std::map<int,int>::iterator itmap = oldnew_.find( tip_id );
      if( itmap == oldnew_.end() )
        dserror( "tip node must be present in oldnew_ data \n" );
      dupid = itmap->second;
      disppl[dupid] = disp_temp;
    }
  }

  LINALG::GatherAll( disppl, comm_ ); //-------------

  clearCondns_ = true;

  DRT::CRACK::UTILS::AddConditions( discret_, disppl );

  discret_->FillComplete();*/
  /************************************************************/ //-------------

  //TODO: how to handle tip nodes?
  //tipnodes_.clear();
  for( std::map<int,int>::iterator it = oldnew_.begin(); it!= oldnew_.end(); it++ )
  {
    int idcheck = it->first;
    if( std::find( tipnodes_.begin(), tipnodes_.end(), idcheck ) != tipnodes_.end() )
      tipnodes_.push_back( it->second );
  }

  //TODO: push the newly formed nodes into tip nodes --> for fsi-crack problem

  alreadySplit_ = true;
}

/*-----------------------------------------------------------------------*
 * Returns true if the given node id is on the boundary          sudhakar 01/15
 *-----------------------------------------------------------------------*/
bool DRT::CRACK::PropagateTip::is_A_BoundaryNode( int & id )
{
  if( boun_nodes_.find( id ) == boun_nodes_.end() )
    return false;
  return true;
}

/*-----------------------------------------------------------------------*
 * Output the location all crack surface nodes in current conf   sudhakar 01/15
 *-----------------------------------------------------------------------*/
void DRT::CRACK::PropagateTip::WriteCrackSurfacePoints()
{
  const std::string filebase = "cracklocation";
  std::ostringstream pid_stream;
  pid_stream << ".crackTip" << segment_id_ <<".txt";

  std::ostringstream filename;
  filename<< filebase<< pid_stream.str();

  std::ofstream f;
  f.open( filename.str().c_str() );

  for( std::set<int>::iterator it = cracknodes_.begin(); it != cracknodes_.end(); it++ )
  {
    int nid = *it;
    DRT::Node * node = discret_->gNode( nid );
    const double * coo = node->X();

    std::vector<double> disp_node = DRT::CRACK::UTILS::getDisplacementNode( discret_, node, disp_col_ );

    f << coo[1] + disp_node[1] <<" "<<coo[0]<<"\n";
  }
  f.close();

}
