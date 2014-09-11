/*----------------------------------------------------------------------*/
/*!
\file crackUtils.cpp

\brief Utility functions for crack propagation problem

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "crackUtils.H"
#include "crack_tolerance.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils.H"

/*----------------------------------------------------------------------------------------*
 * When a new node is introduced in the discret, DOFs corresponding to the new node are zeros.
 * We copy the field values for this node from the old node that is already         sudhakar 03/14
 * existing at the same position
 *----------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( const Teuchos::RCP<DRT::Discretization>& discret,
                                                     Teuchos::RCP<Epetra_Vector>& vec,
                                                     const std::map<int,int>& oldnewIds )
{
  Teuchos::RCP<Epetra_Vector> old = vec;
  vec = LINALG::CreateVector(*discret->DofRowMap(),true);

  LINALG::Export( *old, *vec );

  for( std::map<int,int>::const_iterator it = oldnewIds.begin(); it != oldnewIds.end(); it++ )
  {
    int oldid = it->first;
    int newid = it->second;

    DRT::UTILS::EquateValuesAtTheseNodes( *vec, discret, oldid, newid );
  }
}

/*--------------------------------------------------------------------------------------*
 * Find angle of a vector pointing from point a to b and horizontal plane       sudhakar 04/14
 * The calculated angle is in the range [0,2pi]
 *--------------------------------------------------------------------------------------*/
double DRT::CRACK::UTILS::FindAngle( const double a[3], const double b[3] )
{
  double normal[3];
  normal[0] = b[0] - a[0];
  normal[1] = b[1] - a[1];
  normal[2] = 0.0;

  double angle = atan2( normal[1], normal[0] );

  convertAngleTo_02PI_range( angle );

  return angle;
}

/*-----------------------------------------------------------------------------------------*
 * convert the angle in range [-pi,pi] to [0,2pi]                                 sudhakar 04/14
 * This function is usually called after using atan2() to find angle
 *-----------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::convertAngleTo_02PI_range( double & ang )
{
  double PI = 22.0/7.0;
  if( (fabs(ang) < ANGLE_TOL_ZERO) or (fabs(ang - 2.0*PI) < ANGLE_TOL_ZERO) )
    ang = 0.0;
  else if( ang < 0.0 )
    ang = 2.0 * PI + ang;
}

/*-----------------------------------------------------------------------------------------*
 * convert the angle in range [0,2pi] to [-pi,pi]                                  sudhakar 07/14
 *-----------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::convertAngleTo_PI_mPI_range( double & ang )
{
  double PI = 22.0/7.0;

  // angle is zero
  if( (fabs(ang) < ANGLE_TOL_ZERO) or (fabs(ang - 2.0*PI) < ANGLE_TOL_ZERO) )
  {
    ang = 0.0;
    return;
  }

  // already angle is in the correct form
  if( fabs(ang) < PI )
    return;

  if( ang > 0.0 and ang > PI )
  {
    ang = -1.0* ( 2.0 * PI - ang );
    return;
  }

  if( ang < 0.0 and fabs(ang) > PI )
  {
    ang = 2.0 * PI + ang;
    return;
  }

  else
    dserror( "Angle should satisfy one of the above cases\n" );

  return;
}

/*-----------------------------------------------------------------------------------------*
 * Returns true if the given element has given nodeid                              sudhakar 04/14
 *-----------------------------------------------------------------------------------------*/
bool DRT::CRACK::UTILS::ElementHasThisNodeId( const DRT::Element* ele, int nodeid )
{
  bool found = false;

  const int* nodes = ele->NodeIds();
  const int numNodes = ele->NumNode();

  for( int iNode = 0; iNode < numNodes; iNode++ )
  {
    int surnodeid = nodes[iNode];

    if( surnodeid == nodeid )
    {
      found = true;
      break;
    }
  }
  return found;
}

/*-----------------------------------------------------------------------------------------*
 * Returns true if the given element has given nodeid                              sudhakar 07/14
 *-----------------------------------------------------------------------------------------*/
bool DRT::CRACK::UTILS::ElementHasThisNodeId( const Teuchos::RCP<DRT::Element>& ele, int nodeid )
{
  bool found = false;

  const int* nodes = ele->NodeIds();
  const int numNodes = ele->NumNode();

  for( int iNode = 0; iNode < numNodes; iNode++ )
  {
    int surnodeid = nodes[iNode];

    if( surnodeid == nodeid )
    {
      found = true;
      break;
    }
  }
  return found;
}

/*-----------------------------------------------------------------------------------------*
 * Delete all Dirichlet conditions that has only one associated node                sudhakar 06/14
 * These are the conditions that we added to fix crack tip to desire location
 *-----------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::deleteConditions( Teuchos::RCP<DRT::Discretization> discret )
{
  std::vector<std::multimap<std::string,Teuchos::RCP<Condition> >::iterator> del;

  std::multimap<std::string,Teuchos::RCP<Condition> >::iterator conit;
  std::multimap<std::string,Teuchos::RCP<Condition> >& allcondn = discret->GetAllConditions();
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
  }
}

/*------------------------------------------------------------------------------------*
 * Add the given node ids and Displacement conditions to discretization        sudhakar 06/14
 * Type of the condition added is PointDirichlet
 * For each node, a separate condition is added. This is helpful to identify
 * these conditions later and delete them at appropriate time
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::AddConditions( Teuchos::RCP<DRT::Discretization> discret,
                                       const std::map<int, std::vector<double> >& ale_bc_nodes )
{
  std::map<int, std::vector<double> >::const_iterator iter;
  int id = 101;
  for( iter = ale_bc_nodes.begin(); iter != ale_bc_nodes.end(); iter++ )
  {
    Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp( new DRT::Condition( id++, DRT::Condition::PointDirichlet, false, DRT::Condition::Point ) );
    std::vector<int> onoff(3,1);

    cond->Add( "Node Ids", iter->first );
    cond->Add("onoff",onoff);
    cond->Add("val",iter->second);

    discret->SetCondition( "Dirichlet", cond );
  }
}

/*------------------------------------------------------------------------------------*
 * Add some new nodes to the given condition                                  sudhakar 03/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::addNodesToConditions( DRT::Condition * cond,
                                              std::map<int,int> oldnew )
{
  std::vector<int> add;

  for( std::map<int,int>::iterator itm = oldnew.begin(); itm != oldnew.end(); itm++ )
  {
    add.push_back(itm->first);
    add.push_back(itm->second);
  }
  addNodesToConditions( cond, add );
}

/*------------------------------------------------------------------------------------*
 * Add some new nodes to the given condition                                  sudhakar 03/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::addNodesToConditions( DRT::Condition * cond,
                                              std::vector<int> add )
{
  if( add.size() == 0 )
    dserror("No nodes passed for adding to condition\n");

  const std::vector<int>* conNodes = cond->Nodes();
  for( unsigned ii = 0; ii < conNodes->size(); ii++ )
      add.push_back( (*conNodes)[ii] );

  Teuchos::RCP<std::vector<int> > storage = Teuchos::rcp( new std::vector<int>(add) );

  std::sort( storage->begin(), storage->end() );
  cond->Add( "Node Ids", storage );
}

/*------------------------------------------------------------------------------------*
 * Add some new nodes to the given condition                                  sudhakar 05/14
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::addNodesToConditions( DRT::Condition * cond,
                                              std::set<int> add )
{
  std::vector<int> vec (add.size());

  int count = 0;
  for( std::set<int>::iterator it = add.begin(); it != add.end(); it++ )
  {
    vec[count] = *it;
    count++;
  }

  addNodesToConditions( cond, vec );
}

double DRT::CRACK::UTILS::ComputeDiracDelta( const double & dx, const double & dy, const double & h )
{
  double rx = fabs(dx)/h;
  double ry = fabs(dy)/h;

  if( rx > ( 2.0 + APPROX_ZERO ) or ry > ( 2.0 + APPROX_ZERO ) )
    return 0.0;

  double dirac = HatFunction( rx ) * HatFunction( ry );
  return dirac;
}

double DRT::CRACK::UTILS::HatFunction( const double & r )
{
  double hat = 0.0;

  if( r < 1.0 )
    hat = 3.0 - 2.0*r + sqrt( 1.0 + 4.0*r-4.0*r*r );
  else if( r < 2.0 + APPROX_ZERO )
    hat = 5.0 - 2.0*r + sqrt( -7.0 + 12.0*r-4.0*r*r );
  else
    dserror( "Already this condition is dealt with in DRT::CRACK::UTILS::ComputeDiracDelta()\n" );

  hat = hat * 0.125;
  return hat;
}

double DRT::CRACK::UTILS::SineHatFunction( const double &r, double n )
{
  double PI = 22.0/7.0;

  if( r > n )
    return 0.0;

  double hat = 0.5 * (1.0+cos(PI*fabs(r)/n));
  return hat;
}

double DRT::CRACK::UTILS::QuadraticHatFunction( const double &r, double n )
{
  if( r > n )
    return 0.0;

  double hat = (n*n-r*r)/(n*n);
  return hat;
}

/*----------------------------------------------------------------------------------------------------------------------------------*
 * Modify element connectivity of given discretization
 * discret ---> given discretization
 * delEle  ---> elements for which connectivity is modified (key -- element id, value -- dummy)                             sudhakar 09/14
 * oldnew  ---> key -- contains old nodes , val -- contains new nodes
 * This discretization, for delEle elements gets new nodes in place of old nodes
 *----------------------------------------------------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::ModifyElementConnectivity( Teuchos::RCP<DRT::Discretization>& discret, std::map<int,int> delEle, std::map<int,int> oldnew )
{
  for( std::map<int,int>::iterator num = delEle.begin(); num != delEle.end(); num++ )
  {
    int eleid = num->first;
    if( discret->HaveGlobalElement( eleid ) )
    {
      bool del = false;
      DRT::Element * ele = discret->gElement( eleid );
      const int * oldnodes = ele->NodeIds();
      std::vector<int> newnodes( ele->NumNode() );

      for( int i = 0; i < ele->NumNode(); i++ )
      {
        std::map<int, int>::iterator delnod = oldnew.find( oldnodes[i] );
        if( delnod != oldnew.end() )
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
}

/*------------------------------------------------------------------------------------*
 * condition maps are for the discretization is built in the initial setup         sudhakar 01/14
 * After introducing new crack tip nodes, we modify the condition maps accordingly
 * (This is just a work-around. The vector copying operations can be avoided if
 * we can modify the way conditions are generated in Baci)
 *------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::AddConditions( Teuchos::RCP<DRT::Discretization>& discret, std::map<int,int> oldnew )
{
  std::multimap<std::string,Teuchos::RCP<Condition> >& allcondn = discret->GetAllConditions();

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

    for( std::map<int,int>::iterator it = oldnew.begin(); it != oldnew.end(); it++ )
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

/*--------------------------------------------------------------------------------------------------*
 * Get two node ids that are in direct connection with th given tip_id                      sudhakar 09/14
 *
 *
 *        3 o---------  2                     * tip node
 *          |         |                       o connected nodes
 *          |         |
 *          |         |                       conn1_index -- 1
 *          *---------o                       conn2_index -- 3
 *         0           1
 *--------------------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::getConnectedNodeIdIndex( const Teuchos::RCP<DRT::Element> & surfele,
                                                 int & tip_id,
                                                 int & conn1_index,
                                                 int & conn2_index )
{
  const int * surfnodeids = surfele->NodeIds();
  int numnode = surfele->NumNode();

  int spl_index = 0;
  for (int num = 0; num < numnode; num++)
  {
    if( tip_id == surfnodeids[num] )
    {
      spl_index = num;
      break;
    }
  }

  conn1_index = (spl_index+1)%numnode;
  if( spl_index == 0 )
    conn2_index = numnode-1;
  else
    conn2_index = spl_index-1;
}
