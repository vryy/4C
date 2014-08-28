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
#include "crack_tolerance.H"
#include "crackUtils.H"
#include "aleCrack.H"
#include "SplitHexIntoTwoWedges.H"

#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_ale/ale.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_createdis.H"

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
 disp_col_(Teuchos::null),
 comm_( discret_->Comm() ),
 clearCondns_( false ),
 myrank_( comm_.MyPID() )
{
  PI_ = 22.0 / 7.0;
  min_angle_tol_ = MIN_PROP_ANGLE * PI_ / 180.0;

  strIsSplit_ = false;
  crackInitiated_ = false;

  // get the initial crack tip nodes when analyzing propagation of an existing crack
  DRT::Condition* boundpts = discret_->GetCondition( "CrackBoundaryPoints" );

  if( boundpts == NULL )
    dserror( "Crack boundary nodes unspecified\n" );

  const std::vector<int>* bounnodes = const_cast<std::vector<int>* >(boundpts->Nodes());
  boun_nodes_.resize( bounnodes->size() );

  for( unsigned i=0; i< bounnodes->size(); i++ )
  {
    int val = (*bounnodes)[i];
    boun_nodes_[i] = val;
  }

  // get the initial crack tip nodes when analyzing propagation of an existing crack
  DRT::Condition* crackpts = discret_->GetCondition( "CrackInitiationPoints" );

  if( crackpts != NULL )
  {
    const std::vector<int>* tipnodes = const_cast<std::vector<int>* >(crackpts->Nodes());
    tipnodes_.resize( tipnodes->size() );

    for( unsigned i=0; i< tipnodes->size(); i++ )
    {
      int val = (*tipnodes)[i];
      tipnodes_[i] = val;
    }

    crackInitiated_ = true;
  }

  if( crackInitiated_  )
  {
    DRT::Condition* maspts = discret_->GetCondition( "masterCrackSurface" );
    DRT::Condition* slapts = discret_->GetCondition( "slaveCrackSurface" );

    const std::vector<int>* masternodes = const_cast<std::vector<int>* >(maspts->Nodes());
    const std::vector<int>* slavenodes = const_cast<std::vector<int>* >(slapts->Nodes());

    //if( masternodes->size() != slavenodes->size() )
    //  dserror("There should be equal number of master and slave nodes\n");

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
  }

  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  DRT::Element * tempele = discret_->lRowElement(0);

  Teuchos::RCP<MAT::Material> actmat = tempele->Material();

  // get material id from element materials
  // This will be used to set materials approp when splitting HEX into WEDGEs
  bool mat_found = false;
  for( std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator itm = mats.begin(); itm != mats.end(); itm++ )
  {
    if( itm->second->Type() == actmat->MaterialType() )
    {
      material_id_ = itm->first;
      mat_found = true;
    }
  }

  if( not mat_found )
    dserror("given material is not found in material maps?\n");

  switch(actmat->MaterialType())
  {
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
      dserror("material type  not supported for crack simulations");
      break;
  }

  step_ = 1;

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

  startNewNodeId_ = crackparam.get<int>("START_NEW_NODE_ID");
  startNewEleId_ = crackparam.get<int>("START_NEW_ELE_ID");


  //-------------------------------------------------------------------------
  //---Preparation for writing crack tip location into a file----------------
  tip_string_ = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                    + ".crackTip.txt";
  std::ofstream f;
  f.open(tip_string_.c_str(),std::fstream::trunc);
  //--------------------------------------------------------------------------

  //-------------------------------------------------------------------------
  //---Store ALE line condition nodes-----------------------------------------
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
        ALE_line_nodes_.insert( (*nodes)[nodno] );
    }
  }
  //-------------------------------------------------------------------------

  //this->WriteNodesDebug();
}


/*------------------------------------------------------------------------------------*
 * Perform all the operations related to crack propagation                  sudhakar 11/13
 * This calculates stress intensity factor (K), and if it is  higher than
 * the critical value, crack is introduced into discretization
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::propagateOperations( Teuchos::RCP<const Epetra_Vector>& displace,
                                                      Teuchos::RCP<std::vector<char> >& strdata )
{
  // general geometric paramters
  oldnew_.clear();
  tip_phi_.clear();
  tip_mphi_.clear();

  // ALE boundary conditions on the tip
  tip_bc_disp_.clear();

  //parameters related to splitting hex into wedges
  DoHexSplit_ = false;
  all_split_ele_.clear();

  // parameters related to clearing all conditions
  justClearedCondns_ = false;
  if( clearCondns_ )
  {
    DeleteConditions();
    justClearedCondns_ = true;
    return;
  }

  // Crack has already propagated that the structure is completely split into two
  if( strIsSplit_ )
    return;

  // export "displacement" to column map
  disp_col_ = LINALG::CreateVector( *discret_->DofColMap(), true );
  LINALG::Export(*displace,*disp_col_);

  //-------------
  // STEP 1 : Initiate crack if not already existing
  //-------------
  InitiateCrack( strdata );

  // crack has not yet initiated. nothing more to do here
  if( not crackInitiated_ )
    return;

  std::cout<<"tip nodes = "<<tipnodes_[0]<<"\t"<<tipnodes_[1]<<"\n";

  //-------------
  // STEP 2 : Compute stress-intensity factors at crack tip
  //-------------
  findStressIntensityFactor();

  // this arises at the very beginning of time integration
  // physically uncorrect
  // experienced this when using gen-alpha method  (sudhakar 12/13)
  if( K_I_ < 0.0 )
    return;

  //-------------
  // STEP 3 : Check whether propagation criterion is satisfied
  //-------------
  bool isProp = DoCrackPropagate();
  if( not isProp )
    return;

  //-------------
  // STEP 4 : Decide crack propagation angle from stress-intensity factors
  //-------------
  decidePropagationAngle();

  //-------------
  // STEP 5 : Update crack information
  //-------------
  //std::vector<int> newTip = findNewCrackTip1(); // old method --> working
  //std::vector<int> newTip = findNewCrackTip2(); // second new one --> works only for hex --> search this two correct node etc
  std::vector<int> newTip = findNewCrackTip3(); // combines method 1 and hex element splitting

  //-------------
  // STEP 6 : Write the location of crack tip to a file
  //-------------
  WriteCrackTipLocation( newTip );

  //-------------
  // STEP 7 : ALE step : move new crack tip nodes to correct location
  //-------------
  Perform_ALE_Step();

  //-------------
  // STEP 8 : Introduced new nodes into structural discretization
  //-------------
  updateCrack( newTip );

  //-------------
  // STEP 9 : Check whether the crack has split the body completely into two
  //-------------
  CheckCompleteSplit();

  //-------------
  // STEP 10 : Analyze the discretization, and delete unnecessary nodes in each processor
  //-------------
  analyzeAndCleanDiscretization();
}

/*------------------------------------------------------------------------------------*
 * If not already present, initiate crack into the structure                sudhakar 05/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::InitiateCrack( Teuchos::RCP<std::vector<char> >& strdata )
{

  if( crackInitiated_ )
    return;

  if( strdata == Teuchos::null )
  {
    dserror( "stress data is empty!" );
  }

  //TODO: Is this row or column map? which is right? check!
  const Epetra_Map* elemap = discret_->ElementRowMap();
  Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > mapdata = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  std::vector<char>::size_type position=0;

  for (int i=0;i<elemap->NumMyElements();++i)
  {
    Teuchos::RCP<Epetra_SerialDenseMatrix> gpstress = Teuchos::rcp(new Epetra_SerialDenseMatrix);
    DRT::ParObject::ExtractfromPack(position, *strdata, *gpstress);
    (*mapdata)[elemap->GID(i)]=gpstress;
  }

  const Epetra_Map& elecolmap = *(discret_->ElementColMap());
  DRT::Exporter ex( *elemap, elecolmap, comm_ );
  ex.Export(*mapdata);

  const Epetra_Map* noderowmap = discret_->NodeRowMap();

  Teuchos::ParameterList p;
  p.set("action","postprocess_stress");
  p.set("stresstype","ndxyz");
  p.set("gpstressmap", mapdata);
  Teuchos::RCP<Epetra_MultiVector> nodal_stress = Teuchos::rcp(new Epetra_MultiVector(*noderowmap,6,true));
  p.set("poststress",nodal_stress);
  discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  if (nodal_stress==Teuchos::null)
  {
    dserror("vector containing nodal stresses/strains not available");
  }

  double max_str = 0.0;
  std::vector<int >node_no;

  const int numnodes = discret_->NumMyRowNodes();

  for (int i=0;i<numnodes;++i)
  {
    Epetra_SerialDenseMatrix eigenvec(3,3);
    Epetra_SerialDenseVector eigenval(3);

    int this_node = noderowmap->GID(i);
    if( std::find( boun_nodes_.begin(), boun_nodes_.end(), this_node ) == boun_nodes_.end() )
      continue;

    eigenvec(0,0) = (*((*nodal_stress)(0)))[i];
    eigenvec(0,1) = (*((*nodal_stress)(3)))[i];
    eigenvec(0,2) = (*((*nodal_stress)(5)))[i];
    eigenvec(1,0) = eigenvec(0,1);
    eigenvec(1,1) = (*((*nodal_stress)(1)))[i];
    eigenvec(1,2) = (*((*nodal_stress)(4)))[i];
    eigenvec(2,0) = eigenvec(0,2);
    eigenvec(2,1) = eigenvec(1,2);
    eigenvec(2,2) = (*((*nodal_stress)(2)))[i];

    LINALG::SymmetricEigenProblem(eigenvec, eigenval, true);

    //calculate von-mises stress
    /*double von_mis = pow( (eigenval(0)-eigenval(1)), 2 ) +
                     pow( (eigenval(1)-eigenval(2)), 2 ) +
                     pow( (eigenval(2)-eigenval(0)), 2 ) ;

    von_mis = von_mis * 0.5;*/

    double von_mis = fabs( eigenval(0) );

    if( von_mis > max_str )
    {
      max_str = von_mis;
      node_no.clear();
      node_no.push_back( this_node );
    }
    else if( fabs(von_mis - max_str) < 1e-12 )
      node_no.push_back( this_node );
  }
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
void DRT::CRACK::PropagateCrack::findNormal( const DRT::Node * tipnode,
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

  std::vector<double> tip_disp = getDisplacementNode( tipnode, disp_col_ );
  std::vector<double> sur1_disp = getDisplacementNode( surnode1, disp_col_ );
  std::vector<double> sur2_disp = getDisplacementNode( surnode2, disp_col_ );

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
void DRT::CRACK::PropagateCrack::findStressIntensityFactor()
{
  DRT::Node * tipnode = NULL;
  int lmaster = 0, gmaster = 0;

  if( discret_->HaveGlobalNode( tipnodes_[0] ) )
  {
    tipnode = discret_->gNode( tipnodes_[0] );
    if( tipnode->Owner() == myrank_ )
    {
      lmaster = myrank_;

      DRT::Node * node1 = findNeighboringCrackNode( tipnode, false, NULL );
      DRT::Node * node2 = findNeighboringCrackNode( tipnode, true, node1 );

      findNormal( tipnode, node1, node2 );

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

        std::vector<double> disp_phi = getDisplacementNode( node_phi, disp_col_ );
        std::vector<double> disp_mphi = getDisplacementNode( node_mphi, disp_col_ );

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
    std::cout<<"stress intensity factors = "<<K_I_<<"\t"<<K_II_<<"\n";
}

/*------------------------------------------------------------------------------------*
 * Extract displacement at the given node                                    sudhakar 12/13
 *------------------------------------------------------------------------------------*/
std::vector<double> DRT::CRACK::PropagateCrack::getDisplacementNode( const DRT::Node * node,
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
  int surnodeid=0;
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

  //propAngle_ = 44.5 * PI_ / 180.0;

  if( myrank_ == 0 )
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
  // split a hex element into two wedge elements when necessary
  split_All_HEX_Elements();

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

  std::vector<double> lmtAngle = getLimitAngles( newTip );

  //std::cout<<"limit angle = "<<lmtAngle[0]*180.0/PI_<<"\t"<<lmtAngle[1]*180.0/PI_<<"\n";

  // map of element ids to be modified with the new node
  // "key" contains the element id, and "value" is dummy here
  // map is used to make sure that one element is stored only once
  std::map<int, int> delEle;

  for( unsigned num = 0; num < tipnodes_.size(); num++ )
  {
    // map that stores <proc_id, 1 + id of new node created>
    std::map<int,int> startnewid_map;

    if( discret_->HaveGlobalNode( tipnodes_[num] ) )
    {
      DRT::Node * tipnode = discret_->gNode( tipnodes_[num] );

      Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( startNewNodeId_++, tipnode->X(), tipnode->Owner() ) );
      //Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( totalNodes + num, tipnode->X(), tipnode->Owner() ) );

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

      startnewid_map[myrank_] = startNewNodeId_;
    }

    // the following procedure make sure that a node created on different processor gets same id
    // We get new node id created on each processor, and the maximum value among all
    // processors are stored as the new startNewNodeId_
    startNewNodeId_ = 0;
    LINALG::GatherAll( startnewid_map, comm_ );
    for( std::map<int,int>::iterator startit = startnewid_map.begin(); startit != startnewid_map.end(); startit++ )
    {
      int val = startit->second;
      if( val > startNewNodeId_ )
        startNewNodeId_ = val;
    }
  }

  LINALG::GatherAll( oldnew_, comm_ );
  LINALG::GatherAll( delEle, comm_ );

  if( delEle.size() == 0 )
    dserror("propagated crack. but no elements found to attach new nodes. This leads to incorrect discretization\n");

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

      if( newnodes.size() != static_cast<std::size_t>(ele->NumNode()) )
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

  AddConditions();

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
    cracknodes_.insert( it->second );
    new_ale_bc_nodes_.insert( it->second );
    ALE_line_nodes_.insert( it->second );
  }

  for( unsigned i=0; i< newTip.size(); i++ )
  {
    new_ale_bc_nodes_.insert( newTip[i] );
    ALE_line_nodes_.insert( newTip[i] );
  }
}

/*------------------------------------------------------------------------------------*
 * Get limiting angles for this crack geometry and propagation angle         sudhakar 01/14
 * these angles are used to determine which elements get new nodes and which
 * of them keep the old nodes
 *------------------------------------------------------------------------------------*/
std::vector<double> DRT::CRACK::PropagateCrack::getLimitAngles( const std::vector<int>& newTip )
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

#if 1 // this shud be the case if we move crack nodes (I believe) -- sudhakar 02/14
  ang.push_back( propAngle_ );
#else
  ang.push_back( secAng );
#endif


  std::sort( ang.begin(), ang.end() );

  return ang;
}

/*------------------------------------------------------------------------------------*
 * condition maps are for the discretization is built in the initial setup         sudhakar 01/14
 * After introducing new crack tip nodes, we modify the condition maps accordingly
 * (This is just a work-around. The vector copying operations can be avoided if
 * we can modify the way conditions are generated in Baci)
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::AddConditions()
{
  std::multimap<std::string,Teuchos::RCP<Condition> >& allcondn = discret_->GetAllConditions();

  for( std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit = allcondn.begin();
                                                                     conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<Condition> cond = conit->second;

    //TODO: check how to control adding new nodes to FSIcondition

    // Do not include the new nodes into FSICoupling and XFEMcoupling conditions
    // This is because the interface is built based on these conditions, and we want
    // to control what nodes are added during FSI-crack problem
    // Appropriate nodes are added when building new interface after crack propagation

     /*DRT::Condition::ConditionType ct = cond->Type();
     if( ct == DRT::Condition::FSICoupling or ct == DRT::Condition::XFEMCoupling )
      continue;*/

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

        // sorting is mandatory, as Baci assumes condition nodes are always sorted
        std::sort( storage->begin(), storage->end() );
        cond->Add( "Node Ids", storage );
      }
    }
  }

  // Add the new Dirichlet conditions at the tip node to make sure that the crack
  // propagation angle is maintained
  /*clearCondns_ = true;
  std::map<int, std::vector<double> >::iterator iter;
  for( iter = tip_bc_disp_.begin(); iter != tip_bc_disp_.end(); iter++ )
  {
    Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp( new DRT::Condition() );
    std::vector<int> onoff(3,1);

    cond->SetConditionType(  DRT::Condition::PointDirichlet );
    cond->Add( "Node Ids", iter->first );
    cond->Add("onoff",onoff);
    cond->Add("val",iter->second);

    discret_->SetCondition( "Dirichlet", cond );
  }*/
}

/*------------------------------------------------------------------------------------*
 * Delete the Dirichlet conditions existing at the previous crack tip nodes      sudhakar 01/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::DeleteConditions()
{
  if( not clearCondns_ )
    return;

  DRT::CRACK::UTILS::deleteConditions( discret_ );

  /*std::vector<std::multimap<std::string,Teuchos::RCP<Condition> >::iterator> del;

  std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit;
  std::multimap<std::string,Teuchos::RCP<Condition> >& allcondn = discret_->GetAllConditions();
  for( conit = allcondn.begin(); conit != allcondn.end(); conit++ )
  {
    Teuchos::RCP<Condition> cond = conit->second;
    if( cond->Nodes()->size() == 1 )
      del.push_back( conit );
  }

  for( unsigned i=0; i< del.size(); i++ )
  {
    conit = del[i];
    allcondn.erase( conit );
  }*/

  // already we have cleared it now. unless new condns are set, no need to delete anyting
  clearCondns_ = false;
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
    std::cout<<"Id = "<<cond->Id()<<" condn type = "<<cond->Type()<<" geom disc = "<<cond->GeometryDescription()<<"\n";
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

          DRT::Node * temp1 = searchnodes[(tip_index+1)%numnodes];
          int temp1id = temp1->Id();
          if( std::find( procNodes.begin(), procNodes.end(), temp1id ) == procNodes.end() )
          {
            const double * tempcord = temp1->X();
            double angl = atan2( (tempcord[1]-tipcord[1]), (tempcord[0]-tipcord[0]) );
            DRT::CRACK::UTILS::convertAngleTo_02PI_range( angl );

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
            DRT::CRACK::UTILS::convertAngleTo_02PI_range( angl );

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

  if( tipnodes_.size() != oldnew_tip.size() )
    dserror("for each node, we should have a new tip node\n");

  std::vector<int> newTip;
  for( std::map<int,int>::iterator it = oldnew_tip.begin(); it != oldnew_tip.end(); it++ )
    newTip.push_back( it->second );

  return newTip;
}

/*------------------------------------------------------------------------------------*
 *  find the new crack tip nodes.                                             sudhakar 12/13
 *  We do not update the crack tip nodes here, because we still need the
 *  old values to propagate crack from old position
 *------------------------------------------------------------------------------------*/
std::vector<int> DRT::CRACK::PropagateCrack::findNewCrackTip2()
{
  moveNodes_ = true;

  // distance between present crack tip and the tip in next step, in material coordinates
  double scale_loc = 0.0, scale_glo=0.0;

  std::map<int,int> oldnew_tip;

  // for each possible node, store its corresponding neighbor nodes
  std::map<int,std::vector<int> > poss_neigh_nodes;//TODO: Is this necessary?
  std::set<int> not_ALE_nodes; //nodes that should not be added to ALE boundary conditions

  //TODO: In the first step of crack propagation, we should not reach the surface
  // that is located on the bundary of the domain
  // Implement a procedure to take care of this
  for( unsigned tip = 0; tip < tipnodes_.size(); tip++ )
  {
    int tipid = tipnodes_[tip];
    int gnewtip = 0, lnewtip = 0;

    if( discret_->HaveGlobalNode( tipid ) )
    {
      DRT::Node* tipnode = discret_->gNode( tipid );
      if( tipnode->Owner() == myrank_ )
      {
        const double * tipcord = tipnode->X();

        std::vector<double> disp_tip = getDisplacementNode( tipnode, disp_col_ );
        for( unsigned dim=0; dim<3; dim++ )
          disp_tip[dim] += tipcord[dim];

        // This is same as the propagation angle but in the range [-pi,pi]
        double projangle = propAngle_;
        DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( projangle );

        DRT::Element** ElementsPtr = tipnode->Elements();
        int NumElement = tipnode->NumElement();

        bool found_edge = false;

        // Get all the surface elements that are in same z-value as this node
        // We eliminate all the nodes that has atleast one ALE_line_dirich condition
        // Actually this element also has the crack tip node which has not yet added to ALE_line_dirich
        // If we have the new crack tip node in this element, this results in ALE conditions at
        // 3 nodes in a single element --> may result in concave Quad element
        std::map<int,Teuchos::RCP<DRT::Element> > eleNewTip; //possible surface elements that can have new crack tip elements

        // nodes that are impossible to be new crack tip
        // These are nothing but all nodes of element that already contains ALE line Dirich conditions
        std::set<int> impossiTipNodes;

        for( int eleno = 0; eleno < NumElement; eleno++ )
        {
          DRT::Element* Element = ElementsPtr[eleno];
          std::vector< Teuchos::RCP< DRT::Element > > Surfaces = Element->Surfaces();

          Teuchos::RCP<DRT::Element> surf = getSurfaceThisPlane( Surfaces, tipid );
          if( surf == Teuchos::null )
            dserror("Did not found the surface on same z-plane\n");

          const int* nodes_surf = surf->NodeIds();
          int numnodes = surf->NumNode();

          bool node_line_condn = false;

          for (int num = 0; num < numnodes; num++)
          {

            int this_id = nodes_surf[num];
            if( this_id == tipid )
              continue;

            if( std::find( cracknodes_.begin(), cracknodes_.end(), this_id ) != cracknodes_.end() )
            {
              node_line_condn = true;
              break;
            }
          }

          if( node_line_condn )
          {
            for (int num = 0; num < numnodes; num++)
              impossiTipNodes.insert( nodes_surf[num] );
            continue;
          }

          eleNewTip[Element->Id()] = surf;
        }

        not_ALE_nodes.insert( impossiTipNodes.begin(), impossiTipNodes.end() );

        GetTwoCorrectNewTipEle( eleNewTip, tipid );

        if( eleNewTip.size() != 2 )
          dserror("Assuming hex elements, there should only be 2 elements but we got only %u elements\n",eleNewTip.size());

        DRT::Node* poss = NULL;
        DRT::Node* neigh1 = NULL;
        DRT::Node* neigh2 = NULL;

        int surfeleno = 0;
        for( std::map<int,Teuchos::RCP<DRT::Element> >::iterator itsurf = eleNewTip.begin(); itsurf != eleNewTip.end(); itsurf++ )
        {
          surfeleno++;
          Teuchos::RCP<DRT::Element> surf = itsurf->second;

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

          int poss_id=0;
          if( std::find( impossiTipNodes.begin(), impossiTipNodes.end(), possible1->Id() ) == impossiTipNodes.end() )
          {
            poss_id = possible1->Id();
            poss = possible1;
          }
          else if( std::find( impossiTipNodes.begin(), impossiTipNodes.end(), possible2->Id() ) == impossiTipNodes.end() )
          {
            poss_id = possible2->Id();
            poss = possible2;
          }
          else
            dserror("There should only be one possible node, considering hex elements\n");

          lnewtip = poss_id;

          //TODO: What to do for boundary nodes?
          // Find distance between tip node and possible node in material configuration
          const double * pos_cord = poss->X();
          scale_loc = pow( (pos_cord[0]-tipcord[0]), 2 ) + pow( (pos_cord[1]-tipcord[1]), 2 ) + pow( (pos_cord[2]-tipcord[2]), 2 );
          scale_loc = sqrt(scale_loc);

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

            if( node1 != poss and node2 != poss )
              continue;

            if( node1 == poss )
            {
              if( surfeleno == 1 )
                neigh1 = node2;
              else if( surfeleno == 2 )
                neigh2 = node2;
              else
                dserror( "There should be only two elements\n" );
            }
            else
            {
              if( surfeleno == 1 )
                neigh1 = node1;
              else if( surfeleno == 2 )
                neigh2 = node1;
              else
                dserror( "There should be only two elements\n" );
            }
          }
        }

        if( poss == NULL )
          dserror( "new crack tip node not found" );
        if( neigh1 == NULL or neigh2 == NULL )
          dserror( "Neighboring node for the new tip node not found" );

        const double * newtip_cord = poss->X();
        const double * neigh1_cord = neigh1->X();
        const double * neigh2_cord = neigh2->X();

        std::vector<double> disp1 = getDisplacementNode( neigh1, disp_col_ );
        std::vector<double> disp2 = getDisplacementNode( neigh2, disp_col_ );
        std::vector<double> disp_new = getDisplacementNode( poss, disp_col_ );

        for( unsigned dim=0; dim<3; dim++ )
        {
          disp1[dim] += neigh1_cord[dim];
          disp2[dim] += neigh2_cord[dim];
          disp_new[dim] += newtip_cord[dim];
        }

        double angle1 = atan2( disp1[1] - disp_tip[1], disp1[0] - disp_tip[0] ) - projangle;
        double angle2 = atan2( disp2[1] - disp_tip[1], disp2[0] - disp_tip[0] ) - projangle;
        double angle_new = atan2( disp_new[1] - disp_tip[1], disp_new[0] - disp_tip[0] ) - projangle;

        DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( angle1 );
        DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( angle2 );

        // crack is naturally propagating through an existing element surface
        if( fabs( angle_new ) < min_angle_tol_ )
        {
          found_edge = true;
          // Now we are checking this only once because we simulate pseudo-3D crack propagation
          // For real 3D this has to  be set for all crack nodes
          moveNodes_ = false;
          //break;
        }

        std::cout<<"angle1 = "<<angle1<<" angle2 = "<<angle2<<" newangle = "<<angle_new<<"\n";

        found_edge = true;

        std::vector<double> disp_bc(3,0.0);
        disp_bc[2] = 0.0;

        // Crack propagating along neigh1 --> which is a diagonal
        // Corresponding HEX should be split into two WEDGE elements
        if( fabs(angle1) < TOL_DIAG_NODE )
        {
          if( angle1 * angle_new < 0.0 )
          {
            double tot_ang = fabs( angle1 - angle_new );
            double ratio = fabs(angle_new / tot_ang);
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( ratio * disp1[dim] + (1.0-ratio) * disp_new[dim]) - disp1[dim];
          }
          else
          {
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( disp1[dim] * angle_new - disp_new[dim] * angle1 )/( angle_new - angle1 ) - disp1[dim];
          }
          tip_bc_disp_[neigh1->Id()] = disp_bc;
          lnewtip = neigh1->Id();
          DoHexSplit_ = true;


          SplitEleData( eleNewTip, neigh1->Id() );
        }

        // Crack propagating along neigh2 --> which is a diagonal
        // Corresponding HEX should be split into two WEDGE elements
        else if( fabs(angle2) < TOL_DIAG_NODE )
        {
          if( angle2 * angle_new < 0.0 )
          {
            double tot_ang = fabs( angle2 - angle_new );
            double ratio = fabs(angle_new / tot_ang);
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( ratio * disp2[dim] + (1.0-ratio) * disp_new[dim]) - disp2[dim];
          }
          else
          {
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( disp2[dim] * angle_new - disp_new[dim] * angle2 )/( angle_new - angle2 ) - disp2[dim];
          }
          tip_bc_disp_[neigh2->Id()] = disp_bc;
          lnewtip = neigh2->Id();
          DoHexSplit_ = true;

          SplitEleData( eleNewTip, neigh2->Id() );
        }

        // Crack does not propagate through diagonal, but through a normal node
        // No need of splitting any HEX elements
        else
        {
          if( angle1 * angle_new < 0.0 )
          {
            double tot_ang = fabs( angle1 - angle_new );
            double ratio = fabs(angle1 / tot_ang);
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( ratio * disp_new[dim] + (1.0-ratio) * disp1[dim]) - newtip_cord[dim];
          }
          else if( angle2 * angle_new < 0.0 )
          {
            double tot_ang = fabs( angle2 - angle_new );
            double ratio = fabs(angle2 / tot_ang);
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( ratio * disp_new[dim] + (1.0-ratio) * disp2[dim] ) - newtip_cord[dim];
          }
          else
          {
            for( unsigned dim=0; dim<2; dim++ )
              disp_bc[dim] = ( disp2[dim] * angle_new - disp_new[dim] * angle2 )/( angle_new - angle2 ) - newtip_cord[dim];
          }

          tip_bc_disp_[poss->Id()] = disp_bc;
        }


        /*****************/
        /*if( all_split_ele_.size() == 0 )
        {
          std::vector<double> half_disp_bc(3,0.0);
          for(unsigned dim=0;dim<2;dim++)
            half_disp_bc[dim] = 0.7*disp_bc[dim];
          std::vector<int> vec_neigh;
          if( ALE_line_nodes_.find( neigh1->Id() ) == ALE_line_nodes_.end() )
          {
            tip_bc_disp_[neigh1->Id()] = half_disp_bc;
            vec_neigh.push_back( neigh1->Id() );
          }
          if( ALE_line_nodes_.find( neigh2->Id() ) == ALE_line_nodes_.end() )
          {
            tip_bc_disp_[neigh2->Id()] = half_disp_bc;
            vec_neigh.push_back( neigh2->Id() );
          }
          poss_neigh_nodes[poss->Id()] = vec_neigh;
        }*/
        /*****************/
        if( not found_edge )
          dserror("New crack tip not found\n");
      }
    }
    comm_.SumAll( &lnewtip, &gnewtip, 1 );
    oldnew_tip[tipid] = gnewtip;
  }

  comm_.SumAll( &scale_loc, &scale_glo, 1 );
  LINALG::GatherAll( tip_bc_disp_, comm_ );
  LINALG::GatherAll( not_ALE_nodes, comm_ );

  /************/
  /*for( std::map<int, std::vector<int> >::iterator itm = poss_neigh_nodes.begin(); itm != poss_neigh_nodes.end(); itm++ )
  {
    std::vector<int> vec_neigh = itm->second;
    for( unsigned neno=0; neno< vec_neigh.size(); neno++ )
      new_ale_bc_nodes_.insert( vec_neigh[neno] );
  }*/
  /************/

  /*******************************************************************/
  // Get displacement of nodes that are located around the crack tip
  // construct a unit vector in the direction of propagation angle from crack tip
  //if( all_split_ele_.size() == 0 )
  /*{
    if( tip_bc_disp_.size() > 0 )
    {
      std::vector<double> newTipMatCord(3,0.0);
      int lmaster=0, gmaster=0;
      {
        std::map<int, std::vector<double> >::iterator it = tip_bc_disp_.begin();
        int id = it->first;
        if( discret_->HaveGlobalNode( id ) )
        {
          DRT::Node* tipnode = discret_->gNode( id );
          if( tipnode->Owner() == myrank_ )
          {
            const double * cord = tipnode->X();
            for( unsigned i=0; i < 3; i++ )
              newTipMatCord[i] = cord[i];

            lmaster = myrank_;
          }
        }
      }

      // making sure that master processor id is available at all processors
      comm_.SumAll( &lmaster, &gmaster, 1 );
      comm_.Broadcast( &newTipMatCord[0], 3, gmaster );

      std::vector<double> tip_ale_dist = tip_bc_disp_.begin()->second;

      //std::vector<double> new_ale_pos(3,0.0);
      //new_ale_pos[2] = 0.0;
      //for( unsigned dim = 0; dim < 2; dim++ )
      //  new_ale_pos[dim] = newTipMatCord[dim] + tip_ale_dist[dim];
      std::vector<double> prop_vec(3);
      prop_vec[0] = cos( propAngle_ );
      prop_vec[1] = sin( propAngle_ );
      prop_vec[2] = 0.0;

      // Choose all the nodes that are located within 5*scale_loc distance from new tip nodes
      for( int inode=0; inode<discret_->NodeRowMap()->NumMyElements(); ++inode )
      {
        DRT::Node* currnode = discret_->lRowNode(inode);

        // already processed (TODO: Is this necessary?)
        if( tip_bc_disp_.find( currnode->Id() ) != tip_bc_disp_.end() )
          continue;

        // there are some nodes that are attached to old crack tip, cant be assigned to ALE condition
        if( not_ALE_nodes.find( currnode->Id() ) != not_ALE_nodes.end() )
          continue;

        // to make sure we do not move boundary nodes
        if( std::find( ALE_line_nodes_.begin(), ALE_line_nodes_.end(), currnode->Id() ) != ALE_line_nodes_.end() )
          continue;

        const double* currloc = currnode->X();
        std::vector<double> disp_curr = getDisplacementNode( currnode, disp_col_ );
        for( unsigned dim = 0; dim < 2; dim++ )
          disp_curr[dim] += currloc[dim];
        disp_curr[2] = 0.0;

        double dx = disp_curr[0]-newTipMatCord[0];
        double dy = disp_curr[1]-newTipMatCord[1];

        double fwd = dx * prop_vec[0] + dy * prop_vec[1];
        if( fabs(fwd) > APPROX_ZERO and fwd < 0.0 )
          continue;

        double currdis = pow( dx, 2 ) + pow( dy, 2 );
        currdis = sqrt( currdis ) / scale_glo;

        double n = 10.0;

        if( currdis > ( n + APPROX_ZERO ) )
          continue;

        //double dirac = DRT::CRACK::UTILS::ComputeDiracDelta( dx, dy, scale_glo );

        //double dirac = DRT::CRACK::UTILS::HatFunction( currdis );

        // hat function defined as hat = 0.5*(1+cos(r*pi/n)) from Immersed boundary paper of Peskin
        //double dirac = DRT::CRACK::UTILS::SineHatFunction( currdis, n );

        //double dirac = DRT::CRACK::UTILS::QuadraticHatFunction( currdis, n );

        //TODO: For all boundary nodes --> may be again call ALE displacement condition

        //double dirac = 1.0 - currdis / 5.0 ;

        double dirac_x = DRT::CRACK::UTILS::QuadraticHatFunction( dx/scale_glo, n );
        double dirac_y = DRT::CRACK::UTILS::QuadraticHatFunction( dx/scale_glo, n );

        //if( dirac < 0.0 or dirac > 1.0 )
        //  dserror("dirac = %lf; but it can never be less than zero or more than 1\n");

        std::vector<double> add_disp(3,0.0);
        add_disp[0] = dirac_x * tip_ale_dist[0];
        add_disp[1] = dirac_y * tip_ale_dist[1];
        add_disp[2] = 0.0;

        tip_bc_disp_[currnode->Id()] = add_disp;
      }


      LINALG::GatherAll( tip_bc_disp_, comm_ );
    }
  }*/
  /*******************************************************************/


  if( tipnodes_.size() != oldnew_tip.size() )
    dserror("for each node, we should have a new tip node\n");

  std::vector<int> newTip;
  for( std::map<int,int>::iterator it = oldnew_tip.begin(); it != oldnew_tip.end(); it++ )
  {
    newTip.push_back( it->second );
  }

  /*if( myrank_ == 0 )
  {
    std::cout<<"----------all tip ale displacement boundar conditions------------\n";
    for( std::map<int,std::vector<double> >::iterator it = tip_bc_disp_.begin(); it != tip_bc_disp_.end(); it++ )
      std::cout<<"node = "<<it->first<<" disp = "<<(it->second)[0]<<" "<<(it->second)[1]<<"\n";
  }*/

  if( myrank_ == 0 )
  {
    std::cout<<"Printing the splitting hex elements\n";
    for( unsigned i=0;i< all_split_ele_.size(); i++ )
      all_split_ele_[i].print();
  }

  //LINALG::GatherAll( all_split_ele_, comm_ );

  return newTip;
}

/*------------------------------------------------------------------------------------*
 *  find the new crack tip nodes.                                             sudhakar 07/13
 *  We do not update the crack tip nodes here, because we still need the
 *  old values to propagate crack from old position
 *------------------------------------------------------------------------------------*/
std::vector<int> DRT::CRACK::PropagateCrack::findNewCrackTip1()
{
  moveNodes_ = true;

  std::map<int,int> oldnew_tip;

  //TODO: In the first step of crack propagation, we should not reach the surface
  // that is located on the bundary of the domain
  // Implement a procedure to take care of this

  for( unsigned tip = 0; tip < tipnodes_.size(); tip++ )
  {
    int tipid = tipnodes_[tip];
    int gnewtip = 0, lnewtip = 0;

    if( discret_->HaveGlobalNode( tipid ) )
    {
      DRT::Node* tipnode = discret_->gNode( tipid );
      if( tipnode->Owner() == myrank_ )
      {
        const double * tipcord = tipnode->X();

        std::vector<double> disp_tip = getDisplacementNode( tipnode, disp_col_ );
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

            std::vector<double> disp1 = getDisplacementNode( node1, disp_col_ );
            std::vector<double> disp2 = getDisplacementNode( node2, disp_col_ );

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
              moveNodes_ = false;
              break;
            }

            if( ( fabs( angle2 ) < min_angle_tol_ ) and ( node2->Id() == possible1->Id() or node2->Id() == possible2->Id() ) )
            {
              found_edge = true;
              lnewtip = node2->Id();

              // Now we are checking this only once because we simulate pseudo-3D crack propagation
              // For real 3D this has to  be set for all crack nodes
              moveNodes_ = false;
              break;
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

              if( found_edge )
              {
                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                if( num == 1 )
                {
                  double ratio = fabs(angle1 / tot_ang);
                  for( unsigned dim=0; dim<3; dim++ )
                    disp_bc[dim] = ( ratio * disp2[dim] + (1.0-ratio) * disp1[dim]) - node1_cord[dim];
                  tip_bc_disp_[node1->Id()] = disp_bc;
                }
                else if( num == 2 )
                {
                  double ratio = fabs(angle2 / tot_ang);
                  for( unsigned dim=0; dim<3; dim++ )
                    disp_bc[dim] = ( ratio * disp1[dim] + (1.0-ratio) * disp2[dim] ) - node2_cord[dim];
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

/*---------------------------------------------------------------------------------------------------*
 * Procedure to find new crack tip nodes. In addition, this also decides
 * whether to split the HEX into two, and whether ALE step is                               sudhakar 08/14
 * mandatory
 *---------------------------------------------------------------------------------------------------*/
std::vector<int> DRT::CRACK::PropagateCrack::findNewCrackTip3()
{
  moveNodes_ = true;

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

        std::vector<double> disp_tip = getDisplacementNode( tipnode, disp_col_ );
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

            std::vector<double> disp1 = getDisplacementNode( node1, disp_col_ );
            std::vector<double> disp2 = getDisplacementNode( node2, disp_col_ );

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
              moveNodes_ = false;
              break;
            }

            if( ( fabs( angle2 ) < min_angle_tol_ ) and ( node2->Id() == possible1->Id() or node2->Id() == possible2->Id() ) )
            {
              found_edge = true;
              lnewtip = node2->Id();

              // Now we are checking this only once because we simulate pseudo-3D crack propagation
              // For real 3D this has to  be set for all crack nodes
              moveNodes_ = false;
              break;
            }

            if( fabs( angle1 ) < TOL_DIAG_NODE and ( node1->Id() != possible1->Id() and node1->Id() != possible2->Id() ))
            {
              if( (angle1 < 0.0 and angle2 > 0.0 ) or
                (angle1 > 0.0 and angle2 < 0.0 ) )
              {
                found_edge = true;
                lnewtip = node1->Id();

                DoHexSplit_ = true;
                lis_split = 1;
                lspl_node_id = node1->Id();
                lspl_ele_id = Element->Id();
                //SplitEleData( Element->Id(), node1->Id() );

                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                double ratio = fabs(angle1 / tot_ang);
                for( unsigned dim=0; dim<3; dim++ )
                  disp_bc[dim] = ( ratio * disp2[dim] + (1.0-ratio) * disp1[dim]) - node1_cord[dim];
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

                DoHexSplit_ = true;
                lis_split = 1;
                lspl_node_id = node2->Id();
                lspl_ele_id = Element->Id();
                //SplitEleData( Element->Id(), node2->Id() );

                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                double ratio = fabs(angle2 / tot_ang);
                for( unsigned dim=0; dim<3; dim++ )
                  disp_bc[dim] = ( ratio * disp1[dim] + (1.0-ratio) * disp2[dim] ) - node2_cord[dim];
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

              if( found_edge )
              {
                double tot_ang = fabs( angle1 - angle2 );
                std::vector<double> disp_bc(3);
                if( num == 1 )
                {
                  double ratio = fabs(angle1 / tot_ang);
                  for( unsigned dim=0; dim<3; dim++ )
                    disp_bc[dim] = ( ratio * disp2[dim] + (1.0-ratio) * disp1[dim]) - node1_cord[dim];
                  tip_bc_disp_[node1->Id()] = disp_bc;
                }
                else if( num == 2 )
                {
                  double ratio = fabs(angle2 / tot_ang);
                  for( unsigned dim=0; dim<3; dim++ )
                    disp_bc[dim] = ( ratio * disp1[dim] + (1.0-ratio) * disp2[dim] ) - node2_cord[dim];
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


void DRT::CRACK::PropagateCrack::GetTwoCorrectNewTipEle( std::map<int,Teuchos::RCP<DRT::Element> >& eles,
                                                         int tipid )
{
  if( eles.size() == 2 )
    return;

  if( eles.size() != 3 )
    dserror( "There should maximum be 3 elements\n" );

  // Get the current position of crack tip node = material coordinate + displacement
  DRT::Node * tipnode = discret_->gNode( tipid );
  std::vector<double> tipcoord = getDisplacementNode( tipnode, disp_col_ );
  const double * tipx = tipnode->X();
  for( unsigned dim = 0; dim < 3; dim++ )
    tipcoord[dim] += tipx[dim];

  std::map<int,Teuchos::RCP<DRT::Element> >::iterator delit;

  double min = 0.0;
  for( std::map<int,Teuchos::RCP<DRT::Element> >::iterator it=eles.begin(); it!=eles.end(); it++ )
  {
    double ang = GetAbsoluteMinAngle( it->second, tipid, tipcoord );
    if( ang > min )
    {
      delit = it;
      min = ang;
    }
  }

  eles.erase( delit );

}

double DRT::CRACK::PropagateCrack::GetAbsoluteMinAngle( Teuchos::RCP<DRT::Element>& ele, int tipid, std::vector<double>& tipcoord )
{
  double minAngle = 10000.0;

  double projangle = propAngle_;
  DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( projangle );

  DRT::Node ** nodes = ele->Nodes();
  for( int inod = 0; inod < ele->NumNode(); inod++ )
  {
    const DRT::Node* nod = nodes[inod];
    if( nod->Id() == tipid )
      continue;

   std::vector<double> nodcoord = getDisplacementNode( nod, disp_col_ );
   const double * nodx = nod->X();
   for( unsigned dim = 0; dim < 3; dim++ )
     nodcoord[dim] += nodx[dim];

   double angle = atan2( nodcoord[1] - tipcoord[1], nodcoord[0] - tipcoord[0] ) - projangle;
   if( fabs(angle) < minAngle )
     minAngle = angle;
  }

  return fabs(minAngle);
}

/*------------------------------------------------------------------------------------*
 *  The tip of crack has reached the boundary of the domain                      sudhakar 12/13
 *  This means two parts of structure are joined only along a line. In this case
 *  we remove the connection between two portions, and allow the structure to
 *  separate into two distinct parts
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::CheckCompleteSplit()
{
  // Count no of tip nodes that are fallling on the boundary
  unsigned nodes_on_boundary = 0;
  for( std::vector<int>::iterator itip = tipnodes_.begin(); itip != tipnodes_.end(); itip++ )
  {
    int tip_id = *itip;
    if( std::find( boun_nodes_.begin(), boun_nodes_.end(), tip_id ) != boun_nodes_.end() )
      nodes_on_boundary++;
  }

  // no nodes are falling on the boundary. Nothing to do
  if( nodes_on_boundary == 0 )
    return;

  // This corresponds to 3D crack, or multiple cracks in 2D or 3D
  else if( nodes_on_boundary != tipnodes_.size() )
    dserror("Only few tip nodes are falling on the boundary. Current crack implementation does not support this\n");

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
        dserror("old tip id not found on this processor");

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
            delEle[ Element->Id() ] = 0;
        }
      }
      Teuchos::RCP<DRT::Node> dupnode = Teuchos::rcp( new DRT::Node( startNewNodeId_++, tipnode->X(), tipnode->Owner() ) );
      oldnew_[tipnode->Id()] = dupnode->Id();
      discret_->AddNode( dupnode );
    }
  }

  LINALG::GatherAll( oldnew_, comm_ );
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

  AddConditions();

  discret_->FillComplete();


  /************************************************************/ //-------------
  // This is to add a prescribed displacement when the body separates completely into two
  //TODO: check to remove the lines marked as //-----------------
  std::map<int,std::vector<double> > disppl;//-----------------
  std::vector<double> disp(3,0.0);//-----------------
  disp[1] = 0.05;//------------------

  for( unsigned num = 0; num < tipnodes_.size(); num++ )
  {
    int tip_id = tipnodes_[num];
    if( discret_->HaveGlobalNode( tip_id ) )
    {
      DRT::Node * tipnode = discret_->gNode( tip_id );
      std::vector<double> disp_temp = this->getDisplacementNode( tipnode, disp_col_ );//-------------
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

  //DRT::CRACK::UTILS::AddConditions( discret_, disppl );

  discret_->FillComplete();
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

  strIsSplit_ = true;
}

/*------------------------------------------------------------------------------------*
 * Perform all operations related to ALE step                                 sudhakar 05/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::Perform_ALE_Step()
{
  if( tip_bc_disp_.size() == 0 )
    return;

  if( 0 )
  {
    for( std::map<int,std::vector<double> >::iterator it = tip_bc_disp_.begin(); it != tip_bc_disp_.end(); it++ )
    {
      int tipid = it->first;
      if( discret_->HaveGlobalNode( tipid ) )
      {
        DRT::Node * tipnode = discret_->gNode( tipid );
        tipnode->ChangePos( it->second );
      }
    }
  }
  else
  {
    aleCrack ale( discret_ );
    ale.ALE_step( tip_bc_disp_, new_ale_bc_nodes_ );


    ALE_line_nodes_ = ale.getLineDirichNodes();
    if( ( step_ == 1 ) and ALE_line_nodes_.size() == 0 )
      dserror( "No Line Dirichlet found on ALE discretization\n" );

    ale.clearALE_discret();
  }
}

/*------------------------------------------------------------------------------------*
 * Write crack tip location to an output file to trace crack path              sudhakar 07/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::WriteCrackTipLocation( const std::vector<int>& newTip )
{
  std::ofstream f;
  if( myrank_ == 0 )
    f.open( tip_string_.c_str(), std::fstream::ate | std::fstream::app );

  int lmaster = 0, gmaster = 0;
  std::vector<double> disp(3,0.0);
  int nodeid = newTip[0];

  if( not tip_bc_disp_.size() == 0 )
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

/*------------------------------------------------------------------------------------*
 * Analyze --> Make sure all nodes are available on each proc for considered
 *             element distribution                                             sudhakar 06/14
 * Clean   --> Delete all extra nodes on each proc that are not needed for
 *             given element distribution
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::analyzeAndCleanDiscretization()
{
  std::set<int> colnodesReqd;  // all nodes that are required in this processor for available element distribution
  std::set<int> colnodesAvail; // all nodes that are availble in this processor for available element distribution

  // get all required nodes on this proc
  int numelements = discret_->NumMyColElements();
  for (int i=0; i<numelements; i++ )
  {
    DRT::Element* actele = discret_->lColElement(i);
    const int* nodes = actele->NodeIds();
    for( int nno = 0; nno < actele->NumNode(); nno++ )
      colnodesReqd.insert( nodes[nno] );
  }

  // get all available nodes on this proc
  int numnodes = discret_->NumMyColNodes();
  const Epetra_Map * colnodemap = discret_->NodeColMap();
  for( int i=0; i<numnodes; i++ )
  {
    int nid = colnodemap->GID( i );
    colnodesAvail.insert( nid );
  }

  std::vector<int> delNodeIds;
  if( colnodesReqd == colnodesAvail )
  {
  }
  else
  {
    std::set<int>::iterator itavail, itreqd;

    // check whether all node ids are available on this proc
    std::vector<int> notavail;
    for( itreqd = colnodesReqd.begin(); itreqd != colnodesReqd.end(); itreqd++ )
    {
      int nid = *itreqd;
      itavail = colnodesAvail.find( nid );
      if( itavail == colnodesAvail.end() )
        notavail.push_back( nid );
    }

    if( not notavail.size() == 0 )
    {
      std::cout<<"The following nodeids are not available in proc "<<myrank_<<" : ";
      for( unsigned i=0; i< notavail.size(); i++ )
        std::cout<<notavail[i]<<" ";
      std::cout<<"\n";
      dserror("failed!\n");
    }

    // node ids that is not required but available in this proc
    // Hence these will  be erased

    for( itavail = colnodesAvail.begin(); itavail != colnodesAvail.end(); itavail++ )
    {
      int nid = *itavail;
      itreqd = colnodesReqd.find( nid );
      if( itreqd == colnodesReqd.end() )
        delNodeIds.push_back( nid );
    }
  }

  for( std::vector<int>::iterator it = delNodeIds.begin(); it != delNodeIds.end(); it++ )
  {
    discret_->DeleteNode( *it );
  }

  discret_->FillComplete();

  step_++;
}

/*-------------------------------------------------------------------------------------------------*
 * Write discretization with nodal Ids for debugging reasons                               sudhakar 06/14
 *-------------------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::WriteNodesDebug()
{
  std::ofstream f;
  std::string nodfil = "nodeids";

  f.open( nodfil.c_str());
  f<<"View \"Element geo\" {\n";

  for( int iele=0; iele<discret_->ElementRowMap()->NumMyElements(); ++iele )
  {
    DRT::Element * curele = discret_->lRowElement( iele );
    f<<"SH(";
    for(int j=0; j<curele->NumNode(); ++j)
    {
      DRT::Node* node = curele->Nodes()[j];
      const double * co = node->X();
      f<<co[0]<<","<<co[1]<<","<<co[2];
      if( j != curele->NumNode()-1 )
        f<<",";
    }
    f<<"){";
    for(int j=0; j<curele->NumNode(); ++j)
    {
      f<<curele->NodeIds()[j];
      if( j != curele->NumNode()-1 )
        f<<",";
    }
    f<<"};\n";
  }
  f<<"};\n";



  f<<"View \"Node Ids\" {\n";
  for( int inod = 0; inod < discret_->NodeRowMap()->NumMyElements(); inod++ )
  {
    DRT::Node * curnod = discret_->lRowNode( inod );
    f<<"SP(";
    const double * co=curnod->X();
    f<<co[0]<<","<<co[1]<<","<<co[2]<<"){"<<curnod->Id()<<"};\n";
  }
  f<<"};\n";

  f.close();

  dserror("done");
}

/*---------------------------------------------------------------------------------------------------------*
 * Store all data that is necessary to split a HEX into WEDGE elements                           sudhakar 07/14
 *---------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::PropagateCrack::SplitEleData( std::map<int,Teuchos::RCP<DRT::Element> >& zele,
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
void DRT::CRACK::PropagateCrack::SplitEleData( int ele_id, int node_id )
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
void DRT::CRACK::PropagateCrack::split_All_HEX_Elements()
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

      split.DoAllSplittingOperations( ele, tip_ids, split_nodes, ele->Id(), startNewEleId_ );
    }

    // it is necessary to have same value on all proc to have idential Ids of same element
    // on all proc
    startNewEleId_++;

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

Teuchos::RCP<DRT::Element> DRT::CRACK::PropagateCrack::getSurfaceThisPlane( std::vector< Teuchos::RCP< DRT::Element > >& Surfaces,
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
