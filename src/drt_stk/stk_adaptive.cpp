#ifdef STKADAPTIVE

#include <iostream>
#include <fstream>
#include <map>

#include <boost/bind.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_dirichletextractor.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"

#include "../drt_inpar/inpar_material.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"

#include "stk_adaptive.H"
#include "../stk_refine/stk_mesh.H"
#include "../stk_refine/stk_refine.H"
#include "../stk_refine/stk_unrefine.H"
#include "../stk_refine/stk_utils.H"
#include "../stk_lib/stk_field_comm.H"
#include "stk_fixedsparsematrix.H"

STK::Adaptive::Adaptive( DRT::Discretization & dis )
  : dis_( dis ),
    refine_log_( "refine.log" )
{
}


STK::Adaptive::~Adaptive()
{
}


void STK::Adaptive::SetupSTKMesh()
{
  meta_ = Teuchos::rcp( new MetaMesh );

  stk::mesh::MetaData & meta_data = meta_->MetaData();

  std::vector<std::string> condition_names;
  dis_.GetConditionNames( condition_names );

  std::vector<std::string> condition_parts;

  // declare a part for each type of condition

  for ( std::vector<std::string>::iterator i=condition_names.begin();
        i!=condition_names.end();
        ++i )
  {
    std::string & name = *i;

    std::vector<DRT::Condition*> conditions;
    dis_.GetCondition( name, conditions );
    for ( std::vector<DRT::Condition*>::iterator j=conditions.begin();
          j!=conditions.end();
          ++j )
    {
      DRT::Condition & cond = **j;
      DRT::Condition::GeometryType gtype = DRT::Condition::Point;
      if ( cond.GeometryDescription() )
      {
        gtype = cond.GType();
      }

      stk::mesh::EntityRank condition_rank;

      switch ( gtype )
      {
      case DRT::Condition::Point:
        condition_rank = stk::mesh::Node;
        break;
      case DRT::Condition::Line:
        condition_rank = stk::mesh::Edge;
        break;
      case DRT::Condition::Surface:
        condition_rank = stk::mesh::Face;
        break;
      case DRT::Condition::Volume:
        condition_rank = stk::mesh::Element;
        break;
      default:
        dserror( "no such geometry type: %d", gtype );
      }

      std::stringstream condname;
      condname << name << " " << cond.Id();
      meta_data.declare_part( condname.str(), condition_rank );
    }
  }

  DefineFields( meta_data );

  // done with meta data

  meta_->Commit();

  mesh_ = Teuchos::rcp( new Mesh( *meta_, MPI_COMM_WORLD ) );

  stk::mesh::BulkData & bulk_data = mesh_->BulkData();

  //const Epetra_Map & noderowmap = *dis_.NodeRowMap();
  //const Epetra_Map & nodecolmap = *dis_.NodeColMap();
  const Epetra_Map & elementrowmap = *dis_.ElementRowMap();
  //const Epetra_Map & elementcolmap = *dis_.ElementColMap();

  int numele = elementrowmap.NumMyElements();
  int * elements = elementrowmap.MyGlobalElements();

  stk::mesh::Part & active = mesh_->ActivePart();
  stk::mesh::PartVector empty ;
  stk::mesh::PartVector add ;
  add.push_back( &active );

  bulk_data.modification_begin();

  for ( int elid = 0; elid < numele; ++elid )
  {
    int egid = elements[elid];
    const DRT::Element* e = dis_.gElement( egid );

    int numnode = e->NumNode();
    const DRT::Node * const * nodes = e->Nodes();

    std::vector<int> nids;

    for ( int nlid = 0; nlid < numnode; ++nlid )
    {
      const DRT::Node * n = nodes[nlid];
      int ngid = n->Id()+1;
      nids.push_back( ngid );
      stk::mesh::Entity * stk_n = bulk_data.get_entity( stk::mesh::Node, ngid );
      if ( stk_n==NULL )
      {
        stk_n = & bulk_data.declare_entity( stk::mesh::Node, ngid, empty );
        bulk_data.change_entity_parts( *stk_n, add, empty );
        double * const coord = stk::mesh::field_data( mesh_->Coordinates() , *stk_n );
        const double* x = n->X();
        coord[0] = x[0];
        coord[1] = x[1];
        coord[2] = x[2];
      }
    }

    stk::mesh::Part * topology_part = NULL;

    switch ( e->Shape() )
    {
    case DRT::Element::quad4:    topology_part = & mesh_->Quad4(); break;
    case DRT::Element::quad8:    topology_part = & mesh_->Quad8(); break;
    case DRT::Element::quad9:    topology_part = & mesh_->Quad9(); break;
    case DRT::Element::tri3:     topology_part = & mesh_->Tri3(); break;
    case DRT::Element::tri6:     topology_part = & mesh_->Tri6(); break;
    case DRT::Element::hex8:     topology_part = & mesh_->Hex8(); break;
    case DRT::Element::hex20:    topology_part = & mesh_->Hex20(); break;
    case DRT::Element::hex27:    topology_part = & mesh_->Hex27(); break;
    case DRT::Element::tet4:     topology_part = & mesh_->Tet4(); break;
    case DRT::Element::tet10:    topology_part = & mesh_->Tet10(); break;
    case DRT::Element::wedge6:   topology_part = & mesh_->Wedge6(); break;
    case DRT::Element::wedge15:  topology_part = & mesh_->Wedge15(); break;
    case DRT::Element::pyramid5: topology_part = & mesh_->Pyramid5(); break;
    case DRT::Element::line2:    topology_part = & mesh_->Line2(); break;
    case DRT::Element::line3:    topology_part = & mesh_->Line3(); break;
    default:
      dserror( "unsupported element shape %d", e->Shape() );
    }

    stk::mesh::Entity * stk_e = & stk::mesh::declare_element( bulk_data, *topology_part, egid+1, &nids[0] );

    bulk_data.change_entity_parts( *stk_e, add, empty );
  }

  // fill condition parts

  for ( std::vector<std::string>::iterator i=condition_names.begin();
        i!=condition_names.end();
        ++i )
  {
    std::string & name = *i;

    std::vector<DRT::Condition*> conditions;
    dis_.GetCondition( name, conditions );
    for ( std::vector<DRT::Condition*>::iterator j=conditions.begin();
          j!=conditions.end();
          ++j )
    {
      DRT::Condition & cond = **j;

      std::stringstream condname;
      condname << name << " " << cond.Id();
      stk::mesh::Part * condition_part = meta_data.get_part( condname.str() );
      if ( condition_part==NULL )
      {
        dserror( "no part for condition '%s'", condname.str().c_str() );
      }

      part_condition_[condition_part] = &cond;

      switch ( condition_part->primary_entity_rank() )
      {
      case stk::mesh::Node:
      {
        add.clear();
        add.push_back( condition_part );
        const std::vector<int> & nodes = *cond.Nodes();
        for ( std::vector<int>::const_iterator k=nodes.begin();
              k!=nodes.end();
              ++k )
        {
          stk::mesh::Entity * n = bulk_data.get_entity( stk::mesh::Node, *k + 1 );
          if ( n!=NULL and n->owner_rank()==bulk_data.parallel_rank() )
          {
            bulk_data.change_entity_parts( *n, add, empty );
          }
        }
        break;
      }
      case stk::mesh::Edge:
        dserror( "not yet" );
        break;
      case stk::mesh::Face:
        dserror( "not yet" );
        break;
      case stk::mesh::Element:
      {
        add.clear();
        add.push_back( condition_part );

        if ( not cond.GeometryDescription() or
             cond.GType() != DRT::Condition::Volume )
        {
          dserror( "expect condition on volume" );
        }

        std::map<int,Teuchos::RCP<DRT::Element> > & elements = cond.Geometry();
        for ( std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator k=elements.begin();
              k!=elements.end();
              ++k )
        {
          int egid = k->first+1;
          stk::mesh::Entity * e = bulk_data.get_entity( stk::mesh::Element, egid );
          if ( e!=NULL and e->owner_rank()==bulk_data.parallel_rank() )
          {
            bulk_data.change_entity_parts( *e, add, empty );
          }
        }
        break;
      }
      default:
        dserror( "unexpected primary entity rank: %d", condition_part->primary_entity_rank() );
      }
    }
  }

  // Dirichlet conditions require special treatment

  SetupDirichletNodeMap();

  // done with data

  bulk_data.modification_end();

  //mesh_->Output( "converted_mesh.exo", "baci mesh" );

  //mesh_->Print( std::cout );
}


void STK::Adaptive::RefineAll()
{
  // select all local elements

  const Epetra_Map & elementrowmap = *dis_.ElementRowMap();

  int numele = elementrowmap.NumMyElements();
  int * elements = elementrowmap.MyGlobalElements();

  std::vector<stk::mesh::EntityKey> eids;
  eids.reserve( numele );
  for ( int i=0; i<numele; ++i )
  {
    stk::mesh::EntityKey key( stk::mesh::Element, elements[i]+1 );
    eids.push_back( key );
  }

  Refine( eids );
}


void STK::Adaptive::RefineHalf()
{
  // select all local elements

  const Epetra_Map & elementrowmap = *dis_.ElementRowMap();

  int numele = elementrowmap.NumMyElements();
  int * elements = elementrowmap.MyGlobalElements();

  std::vector<stk::mesh::EntityKey> eids;
  eids.reserve( numele );
  for ( int i=numele/4; i<3*numele/4; ++i )
  {
    stk::mesh::EntityKey key( stk::mesh::Element, elements[i]+1 );
    eids.push_back( key );
  }

  Refine( eids );
}


void STK::Adaptive::RefineFirst()
{
  const Epetra_Map & elementrowmap = *dis_.ElementRowMap();

  //int numele = elementrowmap.NumMyElements();
  int * elements = elementrowmap.MyGlobalElements();

  std::vector<stk::mesh::EntityKey> eids;

  stk::mesh::EntityKey key( stk::mesh::Element, elements[0]+1 );
  eids.push_back( key );

  Refine( eids );
}


void STK::Adaptive::Refine( const std::vector<stk::mesh::EntityKey> & eids )
{
#if 1
  refine_log_ << "Refine: ";
  for ( std::vector<stk::mesh::EntityKey>::const_iterator i=eids.begin();
        i!=eids.end(); ++i )
  {
    stk::mesh::EntityKey key = *i;
    stk::mesh::print_entity_key( refine_log_ , GetMesh().MetaData() , key );
    refine_log_ << " ";
  }
  refine_log_ << std::endl;
#endif

  RefineSet rs( &*mesh_ );
  mesh_->Modify();
  rs.Refine( eids );
  mesh_->Commit();

  SyncDRT( true );
}


void STK::Adaptive::Unrefine( const std::vector<stk::mesh::EntityKey> & eids )
{
#if 1
  refine_log_ << "Unrefine: ";
  for ( std::vector<stk::mesh::EntityKey>::const_iterator i=eids.begin();
        i!=eids.end(); ++i )
  {
    stk::mesh::EntityKey key = *i;
    stk::mesh::print_entity_key( refine_log_ , GetMesh().MetaData() , key );
    refine_log_ << " ";
  }
  refine_log_ << std::endl;
#endif

  UnrefineSet rs( &*mesh_ );
  mesh_->Modify();
  rs.Unrefine( eids );
  mesh_->Commit();

  SyncDRT( false );
}


void STK::Adaptive::SyncDRT( bool refine )
{
  const Epetra_Map & elementrowmap = *dis_.ElementRowMap();

  int numele = elementrowmap.NumMyElements();
  int * elements = elementrowmap.MyGlobalElements();

  // copy changes back to DRT Discretization

  // collect new and deleted nodes and elements

  std::vector<stk::mesh::EntityKey> add;
  std::vector<stk::mesh::EntityKey> remove;

  stk::mesh::Selector sel = mesh_->OwnedPart() & mesh_->ActivePart();

//   std::map<stk::mesh::EntityRank, const Epetra_Map*> rank_map;
//   rank_map[stk::mesh::Node]    = & noderowmap;
//   rank_map[stk::mesh::Element] = & elementrowmap;

  {
    EntitySet nodes;
    std::set<stk::mesh::EntityKey> element_found;

    const std::vector<stk::mesh::Bucket*> & ebuckets = mesh_->BulkData().buckets( stk::mesh::Element );
    for ( std::vector<stk::mesh::Bucket*>::const_iterator i=ebuckets.begin();
          i!=ebuckets.end();
          ++i )
    {
      stk::mesh::Bucket & bucket = **i;
      if ( sel( bucket ) )
      {
        for ( stk::mesh::Bucket::iterator j=bucket.begin();
              j!=bucket.end();
              ++j )
        {
          stk::mesh::Entity & e = *j;
          int egid = e.key().id()-1;
          if ( elementrowmap.LID( egid )==-1 )
          {
            add.push_back( e.key() );
          }
          else
          {
            element_found.insert( e.key() );
          }

          // keep nodes of owned elements
          stk::mesh::PairIterRelation rn = e.relations( stk::mesh::Node );
          std::transform( rn.begin(), rn.end(),
                          std::inserter( nodes, nodes.begin() ),
                          boost::bind( &stk::mesh::Relation::entity, _1 ) );
        }
      }
    }

    for ( int i=0; i<numele; ++i )
    {
      stk::mesh::EntityKey key( stk::mesh::Element, elements[i]+1 );
      if ( element_found.find( key )==element_found.end() )
      {
        stk::mesh::Entity * e = mesh_->BulkData().get_entity( key );
        if ( e==NULL or not sel( e->bucket() ) )
        {
          remove.push_back( key );
        }
      }
    }

    const Epetra_Map & nodecolmap = *dis_.NodeColMap();
    std::set<stk::mesh::EntityKey> node_found;

    for ( EntitySet::iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      stk::mesh::Entity & n = **i;
      int ngid = n.key().id()-1;

      if ( nodecolmap.LID( ngid )==-1 )
      {
        add.push_back( n.key() );
      }
      else
      {
        node_found.insert( n.key() );
      }
    }

    int numnodes = nodecolmap.NumMyElements();
    int * node_ids = nodecolmap.MyGlobalElements();

    for ( int i=0; i<numnodes; ++i )
    {
      stk::mesh::EntityKey key( stk::mesh::Node, node_ids[i]+1 );
      if ( node_found.find( key )==node_found.end() )
      {
        stk::mesh::Entity * n = mesh_->BulkData().get_entity( key );
        //if ( n==NULL or not sel( n->bucket() ) )
        if ( n==NULL or n->bucket().capacity()==0 )
        {
          remove.push_back( key );
        }
      }
    }

    std::sort( add.begin(), add.end() );
    std::sort( remove.begin(), remove.end() );
  }

#if 0
  std::cout << "   add: ";
  for ( std::vector<stk::mesh::EntityKey>::iterator i=add.begin();
        i!=add.end(); ++i )
  {
    stk::mesh::EntityKey key = *i;
    stk::mesh::print_entity_key( std::cout , GetMesh().MetaData() , key );
    std::cout << " ";
  }
  std::cout << "\n";

  std::cout << "remove: ";
  for ( std::vector<stk::mesh::EntityKey>::iterator i=remove.begin();
        i!=remove.end(); ++i )
  {
    stk::mesh::EntityKey key = *i;
    stk::mesh::print_entity_key( std::cout , GetMesh().MetaData() , key );
    std::cout << " ";
  }
  std::cout << "\n";
#endif

  // add new nodes and elements
  // needs access to parents

  for ( std::vector<stk::mesh::EntityKey>::iterator i=add.begin();
        i!=add.end(); ++i )
  {
    stk::mesh::EntityKey key = *i;
    switch ( key.rank() )
    {
    case stk::mesh::Node:
    {
      stk::mesh::Entity * n = mesh_->BulkData().get_entity( key );

      stk::mesh::VectorField & coordinates = mesh_->Coordinates();
      double * const c = stk::mesh::field_data( coordinates , *n );

      Teuchos::RCP<DRT::Node> node = Teuchos::rcp( new DRT::Node( key.id()-1, c, n->owner_rank() ) );
      dis_.AddNode( node );

      // change conditions

      stk::mesh::PartVector ps;
      n->bucket().supersets( ps );
      for ( stk::mesh::PartVector::iterator i=ps.begin(); i!=ps.end(); ++i )
      {
        stk::mesh::Part * p = *i;

        std::map<stk::mesh::Part*, DRT::Condition*>::iterator pc = part_condition_.find( p );
        if ( pc!=part_condition_.end() )
        {
          DRT::Condition * cond = pc->second;
          std::vector<int> * nodes = cond->GetMutable<std::vector<int> >( "Node Ids" );
          nodes->push_back( node->Id() );
          std::sort( nodes->begin(), nodes->end() );
        }
      }

      break;
    }
    case stk::mesh::Element:
    {
      stk::mesh::Entity * e = mesh_->BulkData().get_entity( key );
      std::vector<int> nids;
      for ( stk::mesh::PairIterRelation rn = e->relations( stk::mesh::Node );
            not rn.empty();
            ++rn )
      {
        nids.push_back( rn->entity()->key().id()-1 );
      }

      // take the nodes of one child element and find its parent

      EntitySet nodes;
      stk::mesh::PairIterRelation rn = e->relations( stk::mesh::Node );
      std::transform( rn.begin(), rn.end(),
                      std::inserter( nodes, nodes.begin() ),
                      boost::bind( &stk::mesh::Relation::entity, _1 ) );

      if ( refine )
      {
        stk::mesh::Entity * parent_e = InactiveParent( mesh_->ActivePart(), nodes, stk::mesh::Element );

        int pgid = parent_e->key().id()-1;
        DRT::Element * parent = dis_.gElement( pgid );

        // create new element from parent's element type

        Teuchos::RCP<DRT::Element> ele = Teuchos::rcp( parent->Clone() );

        // owner ?!

        ele->SetId( key.id()-1 );
        ele->SetNodeIds( nids.size(), &nids[0] );


        dis_.AddElement( ele );
      }
      else
      {
        // find (removed) child element to clone
        //
        // get DRT Node
        // get DRT Element that is meant to be removed and has the same number
        // of nodes

        bool done = false;
        for ( EntitySet::iterator i=nodes.begin(); not done and i!=nodes.end(); ++i )
        {
          stk::mesh::Entity & n = **i;
          DRT::Node * node = dis_.gNode( n.key().id()-1 );
          if ( node==NULL )
            dserror( "Discretization mismatch" );

          int numele = node->NumElement();
          DRT::Element** ele = node->Elements();
          for ( int j=0; j<numele; ++j )
          {
            DRT::Element* element = ele[j];
            if ( element->NumNode()==nodes.size() )
            {
              stk::mesh::EntityKey elementkey( stk::mesh::Element, element->Id()+1 );
              std::vector<stk::mesh::EntityKey>::iterator e =
                std::lower_bound( remove.begin(), remove.end(), elementkey );
              if ( e != remove.end() and *e == elementkey )
              {
                Teuchos::RCP<DRT::Element> ele = Teuchos::rcp( element->Clone() );

                // owner ?!

                ele->SetId( key.id()-1 );
                ele->SetNodeIds( nids.size(), &nids[0] );

                dis_.AddElement( ele );

                done = true;
                break;
              }
            }
          }
        }
        if ( not done )
          dserror( "failed to clone element" );
      }
      break;
    }
    default:
      dserror( "unexpected rank %d", key.rank() );
    }
  }

  // remove unused nodes and elements

  for ( std::vector<stk::mesh::EntityKey>::iterator i=remove.end();
        i!=remove.begin(); )
  {
    --i;

    stk::mesh::EntityKey key = *i;
    switch ( key.rank() )
    {
    case stk::mesh::Node:
    {
      int gid = key.id()-1;
      dis_.DeleteNode( gid );

      // change conditions

      for ( std::map<stk::mesh::Part*, DRT::Condition*>::iterator pc=part_condition_.begin();
            pc!=part_condition_.end();
            ++pc )
      {
        DRT::Condition * cond = pc->second;
        std::vector<int> * nodes = cond->GetMutable<std::vector<int> >( "Node Ids" );

        std::vector<int>::iterator j = std::lower_bound( nodes->begin(), nodes->end(), gid );
        if ( j!=nodes->end() and *j==gid )
        {
          nodes->erase( j );
        }
      }

      break;
    }
    case stk::mesh::Element:
      dis_.DeleteElement( key.id()-1 );
      break;
    default:
      dserror( "unexpected rank %d", key.rank() );
    }
  }

  dis_.CheckFilledGlobally();
  dis_.FillComplete();
}


Teuchos::RCP<Epetra_Vector> STK::Adaptive::GatherFieldData( const std::vector<stk::mesh::FieldBase*> & fields )
{
  Epetra_Map rowmap = *dis_.DofRowMap();
  Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp( new Epetra_Vector( rowmap ) );

  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();
  stk::mesh::Part & owned  = GetMesh().OwnedPart();

  const std::vector<stk::mesh::Bucket*> & nodes = bulk_data.buckets( stk::mesh::Node );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=nodes.begin();
        i!=nodes.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;
    if ( has_superset( bucket, owned ) )
    {
      for ( stk::mesh::Bucket::iterator j=bucket.begin();
            j!=bucket.end();
            ++j )
      {
        stk::mesh::Entity & n = *j;
        DRT::Node * node = Discretization().gNode( n.key().id()-1 );

        if ( node==NULL )
        {
          dserror( "no node %d found", n.key().id() );
        }

        std::vector<int> dofs = Discretization().Dof( node );

        unsigned dofcount = 0;
        for ( std::vector<stk::mesh::FieldBase*>::const_iterator k=fields.begin();
              k!=fields.end();
              ++k )
        {
          const stk::mesh::FieldBase & f = **k;
          unsigned field_data_size = stk::mesh::field_data_size( f, n.bucket() ) / sizeof( double );
          if ( field_data_size == 0 )
          {
            dserror( "no field data on node %d", n.key().id() );
          }

          if ( dofcount+field_data_size > dofs.size() )
          {
            dserror( "buffer overflow" );
          }

          double * data = reinterpret_cast<double*>( stk::mesh::field_data( f , n ) );
          for ( unsigned l=0; l<field_data_size; ++l )
          {
            int lid = rowmap.LID( dofs[dofcount+l] );
            double value = data[l];
            ( *v )[lid] = value;
          }
          dofcount += field_data_size;
        }
      }
    }
  }

  return v;
}

void STK::Adaptive::ScatterFieldData( Teuchos::RCP<Epetra_Vector> v, const std::vector<stk::mesh::FieldBase*> & fields )
{
  Epetra_Map rowmap = *dis_.DofRowMap();

  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();
  stk::mesh::Part & owned  = GetMesh().OwnedPart();

  // loop all owned nodes and copy dof values from Epetra_Vector to given fields

  const std::vector<stk::mesh::Bucket*> & nodes = bulk_data.buckets( stk::mesh::Node );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=nodes.begin();
        i!=nodes.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;
    if ( has_superset( bucket, owned ) )
    {
      for ( stk::mesh::Bucket::iterator j=bucket.begin();
            j!=bucket.end();
            ++j )
      {
        stk::mesh::Entity & n = *j;

        DRT::Node * node = Discretization().gNode( n.key().id()-1 );
        if ( node==NULL )
        {
          dserror( "no node %d found", n.key().id() );
        }

        std::vector<int> dofs = Discretization().Dof( node );

        unsigned dofcount = 0;
        for ( std::vector<stk::mesh::FieldBase*>::const_iterator k=fields.begin();
              k!=fields.end();
              ++k )
        {
          const stk::mesh::FieldBase & f = **k;
          unsigned field_data_size = stk::mesh::field_data_size( f, n.bucket() ) / sizeof( double );
          if ( field_data_size == 0 )
          {
            dserror( "no field data on node %d", n.key().id() );
          }

          if ( dofcount+field_data_size > dofs.size() )
          {
            dserror( "buffer overflow" );
          }

          double * data = reinterpret_cast<double*>( stk::mesh::field_data( f , n ) );
          for ( unsigned l=0; l<field_data_size; ++l )
          {
            int lid = rowmap.LID( dofs[dofcount+l] );
            data[l] = ( *v )[lid];
          }
          dofcount += field_data_size;
        }
      }
    }
  }

  // communicate field data from owner to ghosts

  CommunicateFieldGhosting( bulk_data, fields );

#if 0
  // now that all real nodal values are available get the hanging node values

  const std::vector<stk::mesh::Bucket*> & constraints = bulk_data.buckets( stk::mesh::Constraint );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=constraints.begin();
        i!=constraints.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;

    // If I know the constraint I know the hanging node. And vice versa?

    //if ( has_superset( bucket, owned ) )
    {
      for ( stk::mesh::Bucket::iterator j=bucket.begin();
            j!=bucket.end();
            ++j )
      {
        stk::mesh::Entity & c = *j;

        stk::mesh::PairIterRelation rel = c.relations( stk::mesh::Node );
        if ( not rel.empty() )
        {
          stk::mesh::Entity & hn = *rel[0].entity();
          ++rel;

          int numrealnodes = rel.size();
          double fact = 1.0/numrealnodes;

          for ( std::vector<stk::mesh::FieldBase*>::const_iterator k=fields.begin();
                k!=fields.end();
                ++k )
          {
            const stk::mesh::FieldBase & f = **k;

            unsigned hn_field_data_size = stk::mesh::field_data_size( f, hn.bucket() ) / sizeof( double );
            if ( hn_field_data_size == 0 )
            {
              dserror( "no field data on node %d", hn.key().id() );
            }
            double * hn_data = reinterpret_cast<double*>( stk::mesh::field_data( f , hn ) );
            std::fill( hn_data, hn_data+hn_field_data_size, 0 );

            for ( stk::mesh::PairIterRelation::iterator r=rel.begin(); r!=rel.end(); ++r )
            {
              stk::mesh::Entity & n = *r->entity();
              unsigned n_field_data_size = stk::mesh::field_data_size( f, n.bucket() ) / sizeof( double );
              if ( hn_field_data_size!=n_field_data_size )
              {
                dserror( "field sizes are different" );
              }
              double * n_data = reinterpret_cast<double*>( stk::mesh::field_data( f , n ) );
              for ( unsigned l=0; l<n_field_data_size; ++l )
              {
                hn_data[l] += fact*n_data[l];
              }
            }
          }
        }
      }
    }
  }
#endif
}

#if 0
void STK::Adaptive::Evaluate( Teuchos::ParameterList&              params,
                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                              Teuchos::RCP<Epetra_Vector>          systemvector1,
                              Teuchos::RCP<Epetra_Vector>          systemvector2,
                              Teuchos::RCP<Epetra_Vector>          systemvector3 )

{
  if (!dis_.Filled()) dserror("FillComplete() was not called");
  if (!dis_.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // see what we have for input
  bool assemblemat1 = systemmatrix1!=Teuchos::null;
  bool assemblemat2 = systemmatrix2!=Teuchos::null;
  bool assemblevec1 = systemvector1!=Teuchos::null;
  bool assemblevec2 = systemvector2!=Teuchos::null;
  bool assemblevec3 = systemvector3!=Teuchos::null;

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // call the element's register class preevaluation method
  // for each type of element
  // for most element types, just the base class dummy is called
  // that does nothing
  DRT::ParObjectFactory::Instance().PreEvaluate(dis_,params,systemmatrix1,systemmatrix2,
                                                systemvector1,systemvector2,systemvector3);

  //Element::LocationArray la(dofsets_.size());
  DRT::Element::LocationArray la( 1 );

  //stk::mesh::BulkData & bulk_data = GetMesh().BulkData();
  stk::mesh::Part & owned  = GetMesh().OwnedPart();
  stk::mesh::Part & active = GetMesh().ActivePart();
  stk::mesh::Selector sel = owned & active;

//   const std::vector<stk::mesh::Bucket*> & elements = bulk_data.buckets( stk::mesh::Element );

  // loop over column elements
  const int numcolele = dis_.NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = dis_.lColElement(i);

    actele->LocationVector(dis_,la,false);

    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = la[0].Size();
    if (assemblemat1)
    {
      if (elematrix1.M()!=eledim or elematrix1.N()!=eledim)
        elematrix1.Shape(eledim,eledim);
      else
        memset(elematrix1.A(),0,eledim*eledim*sizeof(double));
    }
    if (assemblemat2)
    {
      if (elematrix2.M()!=eledim or elematrix2.N()!=eledim)
        elematrix2.Shape(eledim,eledim);
      else
        memset(elematrix2.A(),0,eledim*eledim*sizeof(double));
    }
    if (assemblevec1)
    {
      if (elevector1.Length()!=eledim)
        elevector1.Size(eledim);
      else
        memset(elevector1.Values(),0,eledim*sizeof(double));
    }
    if (assemblevec2)
    {
      if (elevector2.Length()!=eledim)
        elevector2.Size(eledim);
      else
        memset(elevector2.Values(),0,eledim*sizeof(double));
    }
    if (assemblevec3)
    {
      if (elevector3.Length()!=eledim)
        elevector3.Size(eledim);
      else
        memset(elevector3.Values(),0,eledim*sizeof(double));
    }

    int err = actele->Evaluate( params,dis_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3 );
    if (err)
      dserror("Proc %d: Element %d returned err=%d",dis_.Comm().MyPID(),actele->Id(),err);

    int eid = actele->Id();

#if 0
    int numnodes = actele->NumNode();
    DRT::Node** nodes = actele->Nodes();
    for ( int i=0; i<numnodes; ++i )
    {
      DRT::Node * node = nodes[i];
      stk::mesh::Entity * n = bulk_data.get_entity( stk::mesh::Node, node->Id()+1 );
      if ( n!=NULL )
      {
        // we demand that all relevant nodes are available in STK as well
        stk::mesh::PairIterRelation rel = n->relations( stk::mesh::Constraint );
        if ( not rel.empty() and rel[0]==n )
        {
          // hanging node

          int numrealnodes = rel.size()-1;
          double fact = 1.0/numrealnodes;
          std::vector<int> dofs = dis_.Dof( node );

          for ( ++rel; not rel.empty(); ++rel )
          {
            stk::mesh::Entity * nn = rel->entity();
            DRT::Node * nnode = dis_.gNode( nn->key().id()-1 );
            std::vector<int> ndofs = dis_.Dof( nnode );

            for ( unsigned j=0; j<dofs.size(); ++j )
            {
              hanging[dofs[j]].push_back( ndofs[j] );
            }
          }
        }
      }
    }

    if ( hangingnodes )
    {
    }
    else
#endif

    {
      if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,la[0].lm_,la[0].lmowner_);
      if (assemblemat2) systemmatrix2->Assemble(eid,elematrix2,la[0].lm_,la[0].lmowner_);
      if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,la[0].lm_,la[0].lmowner_);
      if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,la[0].lm_,la[0].lmowner_);
      if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,la[0].lm_,la[0].lmowner_);
    }
  }
}
#endif

#if 0
void STK::Adaptive::EvaluateDirichlet( double time,
                                       const std::vector<stk::mesh::FieldBase*> & v,
                                       const std::vector<stk::mesh::FieldBase*> * dv,
                                       const std::vector<stk::mesh::FieldBase*> * ddv )
{
  bool usetime = time >= 0.0;

  // We have the list of nodes with Dirichlet constraints. Thus we do not need
  // to search the conditions.

  for ( std::map<stk::mesh::Entity*, stk::mesh::Part*>::iterator i=dirichlet_condition_.begin();
        i!=dirichlet_condition_.end();
        ++i )
  {
    stk::mesh::Entity * n = i->first;
    DRT::Node* actnode    = dis_.gNode( n->key().id()-1 );
    if (!actnode) dserror( "Cannot find global node %d", n->key().id()-1 );

    stk::mesh::Part * bc_part = i->second;
    DRT::Condition & cond     = * part_condition_[bc_part];

    const vector<int>*    curve  = cond.Get<vector<int> >("curve");
    const vector<int>*    funct  = cond.Get<vector<int> >("funct");
    //const vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
    const vector<double>* val    = cond.Get<vector<double> >("val");

    unsigned deg = 0;  // highest degree of requested time derivative
    if ( dv!=NULL )
    {
      if ( v.size() != dv->size() )
        dserror( "nonmatching derivative field" );
      deg += 1;
      if ( ddv!=NULL )
      {
        if ( v.size() != ddv->size() )
          dserror( "nonmatching derivative field" );
        deg += 1;
      }
    }

    // call explicitly the main dofset, i.e. the first column
    std::vector<int> dofs = dis_.Dof(0,actnode);
    const unsigned numdf = dofs.size();

    unsigned dofcount = 0;
    for ( unsigned j=0; j<v.size(); ++j )
    {
      const stk::mesh::FieldBase & f = *v[j];
      unsigned field_data_size = stk::mesh::field_data_size( f, n->bucket() ) / sizeof( double );
      if ( field_data_size == 0 )
      {
        dserror( "no field data on node %d", n->key().id() );
      }

      if ( dofcount+field_data_size > numdf )
      {
        dserror( "buffer overflow" );
      }

      double * data = reinterpret_cast<double*>( stk::mesh::field_data( f , *n ) );
      for ( unsigned l=0; l<field_data_size; ++l )
      {
        int dofnum = dofcount+l;

        // factor given by time curve
        std::vector<double> curvefac( deg+1, 0.0 );
        curvefac[0] = 1.0;
        if (curve and usetime)
        {
          int curvenum = (*curve)[dofnum];
          if (curvenum>=0)
            curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer( time, deg );
        }

        // factor given by spatial function
        double functfac = 1.0;
        if (funct)
        {
          int funct_num = (*funct)[dofnum];
          if (funct_num>0)
            functfac = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(dofnum,
                                                                             actnode->X(),
                                                                             time,
                                                                             &dis_ );
        }

        data[l] = (*val)[dofnum] * functfac * curvefac[0];
        if ( dv!=NULL )
        {
          // assume these vectors contain the same type of fields

          const stk::mesh::FieldBase & df = *( *dv )[j];
          if ( field_data_size != stk::mesh::field_data_size( df, n->bucket() ) / sizeof( double ) )
            dserror( "nonmatching derivative field" );
          double * ddata = reinterpret_cast<double*>( stk::mesh::field_data( df , *n ) );
          ddata[l] = (*val)[dofnum] * functfac * curvefac[1];

          if ( ddv!=NULL )
          {
            const stk::mesh::FieldBase & ddf = *( *ddv )[j];
            if ( field_data_size != stk::mesh::field_data_size( ddf, n->bucket() ) / sizeof( double ) )
              dserror( "nonmatching derivative field" );
            double * dddata = reinterpret_cast<double*>( stk::mesh::field_data( ddf , *n ) );
            dddata[l] = (*val)[dofnum] * functfac * curvefac[2];
          }
        }
      }
      dofcount += field_data_size;
    }
  }
}
#endif

#if 0
void STK::Adaptive::Assemble( stk::mesh::Entity & e,
                              DRT::Element * actele,
                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                              const Epetra_SerialDenseMatrix & Aele,
                              const std::vector<int> & lm,
                              const std::vector<int> & lmowner )
{
  Teuchos::RCP<LINALG::SparseMatrix> sm = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>( systemmatrix );

  stk::mesh::PairIterRelation rel = e.relations( stk::mesh::Node );
  std::map<stk::mesh::Entity*, stk::mesh::Part*>::iterator conditer;

  //const std::vector<int>*    curve = NULL;
  //const std::vector<int>*    funct = NULL;
  const std::vector<int>*    onoff = NULL;
  //const std::vector<double>* val   = NULL;

  int row = 0;

  for ( ; not rel.empty(); ++rel )
  {
    stk::mesh::Entity * n = rel->entity();
    conditer = dirichlet_condition_.find( n );
    if ( conditer != dirichlet_condition_.end() )
    {
      stk::mesh::Part & part = *conditer->second;
      DRT::Condition & cond = *part_condition_[&part];

      //curve  = cond.Get<std::vector<int> >("curve");
      //funct  = cond.Get<std::vector<int> >("funct");
      onoff  = cond.Get<std::vector<int> >("onoff");
      //val    = cond.Get<std::vector<double> >("val");
      break;
    }
    row += dis_.NumDof( dis_.gNode( n->key().id()-1 ) );
  }

  if ( rel.empty() )
  {
    sm->FEAssemble( actele->Id(), Aele, lm, lmowner, lm );
  }
  else
  {
    std::vector<int> dirichlet( lm.size(), 0 );

    stk::mesh::Entity * n = rel->entity();
    std::vector<int> rdofs = dis_.Dof( dis_.gNode( n->key().id()-1 ) );
    for ( unsigned rj=0; rj<rdofs.size(); ++rj, ++row )
    {
      if ( ( *onoff )[rj]!=0 )
      {
        dirichlet[row] = 1;
      }
    }

    for ( ; not rel.empty(); ++rel )
    {
      stk::mesh::Entity * n = rel->entity();
      conditer = dirichlet_condition_.find( n );
      if ( conditer != dirichlet_condition_.end() )
      {
        stk::mesh::Part & part = *conditer->second;
        DRT::Condition & cond = *part_condition_[&part];

        //curve  = cond.Get<std::vector<int> >("curve");
        //funct  = cond.Get<std::vector<int> >("funct");
        onoff  = cond.Get<std::vector<int> >("onoff");
        //val    = cond.Get<std::vector<double> >("val");

        std::vector<int> rdofs = dis_.Dof( dis_.gNode( n->key().id()-1 ) );
        for ( unsigned rj=0; rj<rdofs.size(); ++rj, ++row )
        {
          if ( ( *onoff )[rj]!=0 )
          {
            dirichlet[row] = 1;
          }
        }
      }
    }

    // now call FEAssemble with dirichlet flags
  }
}
#endif


#if 0
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<STK::FixedSparseMatrix> STK::Adaptive::SetupMatrix()
{
  // initialize standard (stabilized) system matrix

  stk::mesh::BulkData & bulk_data = GetMesh().BulkData();
  DRT::Discretization & discret = Discretization();

  const Epetra_Map & dofrowmap  = * discret.DofRowMap();
  const Epetra_Map & noderowmap = * discret.NodeRowMap();

  int numnode = noderowmap.NumMyElements();
  int * nodeids = noderowmap.MyGlobalElements();

  std::map<int, std::set<int> > graph;

  int row = 0;
  for ( int i=0; i<numnode; ++i )
  {
    int rngid = nodeids[i];
    DRT::Node * rn = discret.gNode( rngid );

    stk::mesh::Entity * n = bulk_data.get_entity( stk::mesh::Node, rngid+1 );
    if ( n==NULL )
      dserror( "stk mesh mismatch" );

    const std::vector<int>*    curve = NULL;
    const std::vector<int>*    funct = NULL;
    const std::vector<int>*    onoff = NULL;
    const std::vector<double>* val   = NULL;

    std::map<stk::mesh::Entity*, stk::mesh::Part*>::iterator conditer = dirichlet_condition_.find( n );
    if ( conditer!=dirichlet_condition_.end() )
    {
      stk::mesh::Part & part = *conditer->second;
      DRT::Condition & cond = *part_condition_[&part];

      curve  = cond.Get<std::vector<int> >("curve");
      funct  = cond.Get<std::vector<int> >("funct");
      onoff  = cond.Get<std::vector<int> >("onoff");
      val    = cond.Get<std::vector<double> >("val");
    }

    bool rhanging = false;

    stk::mesh::PairIterRelation constraints = n->relations( stk::mesh::Constraint );
    for ( ; not constraints.empty(); ++constraints )
    {
      stk::mesh::Entity & constr = * constraints->entity();
      stk::mesh::PairIterRelation rrel = constr.relations( stk::mesh::Node );
      if ( not rrel.empty() and rrel[0].entity()==n )
      {
        // a hanging node

        // All lines that are connected to the supporting nodes need to be
        // extended. Unless those lines are Dirichlet constraint.

        rhanging = true;
        //break;
      }
    }

    int numelements = rn->NumElement();
    DRT::Element** elements = rn->Elements();
    if ( elements==NULL )
      dserror( "no elements at node %d", rn->Id() );

    std::vector<int> rdofs = discret.Dof( rn );
    for ( unsigned rj=0; rj<rdofs.size(); ++rj, ++row )
    {
      //std::set<int> & rowset = graph[row];
      std::set<int> & rowset = graph[rdofs[rj]];
      if ( not rhanging and
           ( onoff==NULL or ( *onoff )[rj]==0 ) )
      {
        // non-Dirichlet row
        for ( int k=0; k<numelements; ++k )
        {
          DRT::Element * e = elements[k];
          int numnodes = e->NumNode();
          DRT::Node ** nodes = e->Nodes();
          for ( int l=0; l<numnodes; ++l )
          {
            DRT::Node * cn = nodes[l];

            stk::mesh::Entity * n = bulk_data.get_entity( stk::mesh::Node, cn->Id()+1 );
            if ( n==NULL )
              dserror( "stk mesh mismatch" );

            stk::mesh::PairIterRelation constraints = n->relations( stk::mesh::Constraint );
            for ( ; not constraints.empty(); ++constraints )
            {
              stk::mesh::Entity & constr = * constraints->entity();
              stk::mesh::PairIterRelation crel = constr.relations( stk::mesh::Node );

              if ( not crel.empty() and crel[0].entity()==n )
              {
                // a hanging node
                for ( ++crel; not crel.empty(); ++crel )
                {
                  n = crel->entity();
                  DRT::Node * node = discret.gNode( n->key().id()-1 );
                  std::vector<int> dofs = discret.Dof( node );
                  rowset.insert( dofs.begin(), dofs.end() );
                }
              }
              else
              {
                std::vector<int> cdofs = discret.Dof( cn );
                rowset.insert( cdofs.begin(), cdofs.end() );
              }
            }
          }
        }
      }
      else
      {
        // just diagonal entry on Dirichlet rows and hanging node rows
        rowset.insert( rdofs[rj] );
      }
    }
  }

  // setup graph with row length and column indices per row

  std::vector<int> sizes;
  sizes.reserve( graph.size() );
  for ( std::map<int, std::set<int> >::iterator i=graph.begin(); i!=graph.end(); ++i )
  {
    std::set<int> & rowset = i->second;
    unsigned s = rowset.size();
    sizes.push_back( s );
  }

  Teuchos::RCP<Epetra_CrsGraph> crsgraph =
    Teuchos::rcp( new Epetra_CrsGraph( Copy, dofrowmap, &sizes[0], true ) );

#if 0
  std::ofstream out( "graph.txt" );
#endif

  for ( std::map<int, std::set<int> >::iterator i=graph.begin(); i!=graph.end(); ++i )
  {
    int gid = i->first;
    std::set<int> & rowset = i->second;
    unsigned s = rowset.size();
    std::vector<int> row;
    row.reserve( s );
    row.assign( rowset.begin(), rowset.end() );

    int err = crsgraph->InsertGlobalIndices( gid, row.size(), &row[0] );
    if ( err )
      dserror( "InsertGlobalIndices failed: %d", err );

#if 0
    for ( std::vector<int>::iterator j=row.begin(); j!=row.end(); ++j )
    {
      out << gid << " " << ( *j ) << "\n";
    }
#endif
  }

  crsgraph->FillComplete();

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new DRT::DirichletExtractor());
  dbcmaps_->Setup( dis_ );

  return Teuchos::rcp( new LINALG::SparseMatrix( crsgraph, dbcmaps_ ) );

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new DRT::DirichletExtractor());
  dbcmaps_->Setup( dis_ );

  return Teuchos::rcp( new STK::FixedSparseMatrix( dis_, *mesh_, dbcmaps_->DirichletMap() ) );
}
#endif


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::Adaptive::SetupDirichletNodeMap()
{
  stk::mesh::BulkData & bulk_data = mesh_->BulkData();

  // we store a node-to-condition-part map for fast Dirichlet condition lookup

  dirichlet_condition_.clear();

  // Dirichlet conditions on lesser entities have preference over Dirichlet
  // conditions on higher entities.

  DRT::Condition::ConditionType dirichlet_types[] = {
    DRT::Condition::VolumeDirichlet,
    DRT::Condition::SurfaceDirichlet,
    DRT::Condition::LineDirichlet,
    DRT::Condition::PointDirichlet,
    DRT::Condition::none
  };

  for ( DRT::Condition::ConditionType * condtype = dirichlet_types;
        *condtype!=DRT::Condition::none;
        ++condtype )
  {
    // loop all registered conditions
    for ( std::map<stk::mesh::Part*, DRT::Condition*>::iterator i=part_condition_.begin();
          i!=part_condition_.end();
          ++i )
    {
      stk::mesh::Part & part = *i->first;
      DRT::Condition & cond = *i->second;

      if ( cond.Type() == *condtype )
      {
        // loop all nodal buckets
        const std::vector<stk::mesh::Bucket*> & nodes = bulk_data.buckets( stk::mesh::Node );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator j=nodes.begin();
              j!=nodes.end();
              ++j )
        {
          stk::mesh::Bucket & bucket = **j;
          if ( has_superset( bucket, part ) )
          {
            // insert all nodes of a matching bucket to the
            // dirichlet_condition map
            for ( stk::mesh::Bucket::iterator k=bucket.begin(); k!=bucket.end(); ++k )
            {
              stk::mesh::Entity & n = *k;
              dirichlet_condition_[ &n ] = & part;
            }
          }
        }
      }
    }
  }
}

#endif
