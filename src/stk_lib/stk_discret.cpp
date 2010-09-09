
#ifdef STKADAPTIVE

#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>

#include "stk_algorithm.H"
#include "stk_discret.H"
#include "stk_field_comm.H"
#include "../stk_refine/stk_mesh.H"
#include "../stk_refine/stk_refine.H"
#include "../stk_refine/stk_unrefine.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_dirichletextractor.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"

namespace STK
{

template <class T>
void CopyCondition( std::pair<typename std::map<std::string,T>::const_iterator,
                              typename std::map<std::string,T>::const_iterator> range,
                    Teuchos::ParameterList & condflags )

{
  for ( typename std::map<std::string,T>::const_iterator i=range.first; i!=range.second; ++i )
  {
    std::string name = i->first;
    if ( name!="Node Ids" )
    {
      T value = i->second;      // must be fine to copy value anyway
      condflags.set( name, value );
    }
  }
}

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::Discretization::Discretization( const Epetra_Comm & comm )
  : comm_( comm ),
    algo_( NULL )
{
#if 0
  std::stringstream str;
  str << "refine-" << comm.MyPID() << ".log";
  refine_log_.open( str.str().c_str() );
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::Setup( DRT::Discretization & dis, STK::Algorithm & algo )
{
  algo_ = & algo;
  meta_ = Teuchos::rcp( new MetaMesh );

  stk::mesh::MetaData & meta_data = meta_->MetaData();

  std::vector<std::string> condition_names;
  dis.GetConditionNames( condition_names );

  std::vector<std::string> condition_parts;

  // declare a part for each type of condition

  for ( std::vector<std::string>::iterator i=condition_names.begin();
        i!=condition_names.end();
        ++i )
  {
    std::string & name = *i;

    std::vector<DRT::Condition*> conditions;
    dis.GetCondition( name, conditions );
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
      condname << name << "-" << cond.Type() << "-" << cond.Id();
      std::string cname = condname.str();
      meta_data.declare_part( cname, condition_rank );

      Teuchos::ParameterList & condflags = conditions_[cond.Type()][cond.Id()];
      condflags.set( "Id", cond.Id() );
      condflags.set( "Name", cname );
      condflags.set( "GeometryType", cond.GType() );
      condflags.set( "GeometryDescription", cond.GeometryDescription() );
      condflags.set( "Type", cond.Type() );

      //std::cout << "Dirichlet condition:\n" << condflags << "\n";

      STK::CopyCondition<Teuchos::RCP<std::vector<int> > >       ( cond.IntRange(),    condflags );
      STK::CopyCondition<Teuchos::RCP<std::vector<double> > >    ( cond.DoubleRange(), condflags );
      STK::CopyCondition<std::string>                            ( cond.StringRange(), condflags );
      STK::CopyCondition<Teuchos::RCP<Epetra_SerialDenseMatrix> >( cond.MatRange(),    condflags );
      STK::CopyCondition<Teuchos::RCP<Epetra_MultiVector> >      ( cond.EvecRange(),   condflags );
    }
  }

  // declare a part for each material

  materials_ = DRT::Problem::Instance()->Materials();
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> > & matmap = * materials_->Map();

  for ( std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator i=matmap.begin();
        i!=matmap.end();
        ++i )
  {
    MAT::PAR::Material & mat = *i->second;
    MAT::PAR::Parameter & par = *mat.Parameter();
    std::stringstream str;
    str << par.Name() << par.Id();
    meta_data.declare_part( str.str(), stk::mesh::Element );
  }

  // declare fields

  algo_->declare_fields( meta_data );

  // done with meta data

  meta_->Commit();

  mesh_ = Teuchos::rcp( new Mesh( *meta_, MPI_COMM_WORLD ) );

  stk::mesh::BulkData & bulk_data = mesh_->BulkData();

  //const Epetra_Map & noderowmap = *dis.NodeRowMap();
  //const Epetra_Map & nodecolmap = *dis.NodeColMap();
  const Epetra_Map & elementrowmap = *dis.ElementRowMap();
  //const Epetra_Map & elementcolmap = *dis.ElementColMap();

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
    const DRT::Element* e = dis.gElement( egid );

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

    // consider material parts

    {
      MAT::Material & mat = * e->Material();
      MAT::PAR::Parameter & par = * mat.Parameter();
      std::stringstream str;
      str << par.Name() << par.Id();

      stk::mesh::Part * part = meta_data.get_part( str.str() );
      if ( part==NULL )
        dserror( "part '%s' not found", str.str().c_str() );

      matpar_[part] = & par;

      stk::mesh::PartVector matpart( 1, part );
      bulk_data.change_entity_parts( *stk_e, matpart, empty );
    }
  }

  // fill condition parts

  for ( std::vector<std::string>::iterator i=condition_names.begin();
        i!=condition_names.end();
        ++i )
  {
    std::string & name = *i;

    std::vector<DRT::Condition*> conditions;
    dis.GetCondition( name, conditions );
    for ( std::vector<DRT::Condition*>::iterator j=conditions.begin();
          j!=conditions.end();
          ++j )
    {
      DRT::Condition & cond = **j;

      std::stringstream condname;
      condname << name << "-" << cond.Type() << "-" << cond.Id();
      stk::mesh::Part * condition_part = meta_data.get_part( condname.str() );
      if ( condition_part==NULL )
      {
        dserror( "no part for condition '%s'", condname.str().c_str() );
      }

      //part_condition_[condition_part] = &cond;

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

  //SetupDirichletNodeMap();

  // done with data

  bulk_data.modification_end();

#if 0

  bulk_data.modification_begin();

  std::vector<stk::mesh::EntityProc> ep;

  if ( bulk_data.parallel_rank()==0 )
  {
    ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 271 ), 1 ) );

//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 382 ), 1 ) );
//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 383 ), 1 ) );
  }
  if ( bulk_data.parallel_rank()==3 )
  {
//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 281 ), 1 ) );
    ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 314 ), 1 ) );
//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Node, 342 ), 1 ) );

//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 321 ), 1 ) );
//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 322 ), 1 ) );
//     ep.push_back( stk::mesh::EntityProc( bulk_data.get_entity( stk::mesh::Element, 323 ), 1 ) );
  }

  bulk_data.change_entity_owner(ep);

  bulk_data.modification_end();

#endif

  // create new state and propagate information
  CreateState();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const std::map<int, Teuchos::ParameterList> * STK::Discretization::Condition( DRT::Condition::ConditionType condtype ) const
{
  std::map< DRT::Condition::ConditionType, std::map<int, Teuchos::ParameterList> >::const_iterator j = conditions_.find( condtype );
  if ( j!=conditions_.end() )
  {
    return & j->second;
  }
  return NULL;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::AdaptMesh( const std::vector<stk::mesh::EntityKey> & refine,
                                     const std::vector<stk::mesh::EntityKey> & unrefine )
{
  //GetMesh().Dump( "refine" );
  Refine( refine );
  //GetMesh().Dump( "unrefine" );
  Unrefine( unrefine );
  //GetMesh().Dump( "done" );
  Rebalance();
  CreateState();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::Refine( const std::vector<stk::mesh::EntityKey> & eids )
{
#if 0
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::Unrefine( const std::vector<stk::mesh::EntityKey> & eids )
{
#if 0
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::Rebalance()
{
  // get current graph
  Teuchos::RCP<Epetra_CrsGraph> graph = ElementGraph();
  const Epetra_BlockMap& src_map = graph->Map();

  Teuchos::ParameterList params;
  Teuchos::ParameterList& sublist = params.sublist("Zoltan");

  params.set("PARTITIONING METHOD", "GRAPH");
  sublist.set("LB_APPROACH", "PARTITION");
  sublist.set("GRAPH_PACKAGE", "PARMETIS");
  //params.set("PARTITIONING METHOD", "HYPERGRAPH");

  //sublist.set("DEBUG_LEVEL", "1"); // Zoltan will print out parameters
  //sublist.set("DEBUG_LEVEL", "5");   // proc 0 will trace Zoltan calls
  //sublist.set("DEBUG_MEMORY", "2");  // Zoltan will trace alloc & free

  // do partitioning
  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner( Teuchos::rcp_implicit_cast<const Epetra_CrsGraph>( graph ), params ) );

  // get epetra communication pattern
  Teuchos::RCP<Epetra_Map> target_map = partitioner->createNewMap();
  Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*target_map, src_map));

  // extract communication information
  std::vector<stk::mesh::EntityProc> ep;

  int numsend = importer->NumSend();
  int* exportlids = importer->ExportLIDs();
  int* exportpids = importer->ExportPIDs();

  ep.reserve(numsend);

  Epetra_Map nodemap = NodeRowMap();
  Epetra_IntVector nodeproc(nodemap);
  nodeproc.PutValue(-1);

  STK::Mesh & mesh = GetMesh();
  stk::mesh::BulkData & bulk = mesh.BulkData();
  unsigned rank = bulk.parallel_rank();

  for (int i=0; i<numsend; ++i)
  {
    int lid = exportlids[i];
    int pid = exportpids[i];
    int gid = src_map.GID(lid);

    //std::cout << lid << " " << pid << " " << gid << "\n";

    stk::mesh::Entity * e = bulk.get_entity( stk::mesh::Element , gid );
    if ( has_superset( e->bucket(), mesh.OwnedPart() ) )
    {
      ep.push_back(stk::mesh::EntityProc(e,pid));

      // send edge/face elements along with the main element

      for ( stk::mesh::PairIterRelation rel = e->relations( stk::mesh::Edge );
            not rel.empty();
            ++rel )
      {
        ep.push_back( stk::mesh::EntityProc( rel->entity(), pid ) );
      }

      for ( stk::mesh::PairIterRelation rel = e->relations( stk::mesh::Face );
            not rel.empty();
            ++rel )
      {
        ep.push_back( stk::mesh::EntityProc( rel->entity(), pid ) );
      }
    }

    // look at the nodes of this element
    stk::mesh::PairIterRelation nrel = e->relations( stk::mesh::Node );
    for (stk::mesh::PairIterRelation::iterator inr=nrel.begin(); inr!=nrel.end(); ++inr)
    {
      stk::mesh::Entity * n = inr->entity();
      if (rank==n->owner_rank())
      {
        int id = static_cast<int>( n->identifier() );
        int nlid = nodemap.LID( id );
        if (nlid<0)
          throw std::logic_error("not a lid");
        nodeproc[nlid] = std::max(nodeproc[nlid],pid);
      }
    }
  }

  // add node communications
  for (int i=0; i<nodeproc.MyLength(); ++i)
  {
    int pid = nodeproc[i];
    if (pid>-1)
    {
      int gid = nodemap.GID(i);

      stk::mesh::Entity * n = bulk.get_entity( stk::mesh::Node , gid );
      if ( stk::mesh::has_superset( n->bucket(), mesh.OwnedPart() ) )
      {
        ep.push_back( stk::mesh::EntityProc( n, pid ) );

        // move hanging node constraints along with the node
        stk::mesh::PairIterRelation pi = n->relations( stk::mesh::Constraint );
        for ( ; not pi.empty(); ++pi )
        {
          // a hanging node always is the first node in its constraint
          if ( pi->entity()->relations( stk::mesh::Node )->entity() == n )
          {
            ep.push_back( stk::mesh::EntityProc( pi->entity(), pid ) );
            break;
          }
        }
      }
    }
  }

  //stk::mesh::sort_unique(ep);

  std::sort( ep.begin() , ep.end() , stk::mesh::EntityLess() );
  std::vector<stk::mesh::EntityProc>::iterator i = std::unique( ep.begin() , ep.end() );
  ep.erase( i , ep.end() );

  // do communication
  mesh.Modify();

  bulk.change_entity_owner(ep);

  mesh.Dump( "rebalance" );

  mesh.Commit();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::CreateState()
{
  state_ = Teuchos::rcp( new STK::DiscretizationState( *this ) );

  std::vector<stk::mesh::FieldBase*> fields;

  algo_->collect_unknowns( fields );

  state_->Setup( fields );

  algo_->notify_state_changed();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> STK::Discretization::NodeGraph()
{
  STK::Mesh & mesh = GetMesh();
  std::map<int,std::set<int> > graph;

  // Get all owned nodes
  stk::mesh::Selector ownedselector = mesh.OwnedPart() & mesh.ActivePart();
  const std::vector<stk::mesh::Bucket*> & buckets = mesh.BulkData().buckets( stk::mesh::Node );

  // loop buckets and entries
  for (std::vector<stk::mesh::Bucket*>::const_iterator i=buckets.begin();
       i!=buckets.end();
       ++i)
  {
    stk::mesh::Bucket & bucket = **i;
    if ( ownedselector( bucket ) )
    {
      for (stk::mesh::Bucket::iterator ib=bucket.begin();
           ib!=bucket.end();
           ++ib)
      {
        stk::mesh::Entity & e = *ib;
        int id = static_cast<int>( e.identifier() );
        std::set<int> & row = graph[id];

        // all elements connected
        // we expect proper ghosting here!
        stk::mesh::PairIterRelation erel = e.relations( stk::mesh::Element );
        for (stk::mesh::PairIterRelation::iterator ier=erel.begin(); ier!=erel.end(); ++ier)
        {
          // all nodes connected to the element
          stk::mesh::PairIterRelation nrel = ier->entity()->relations( stk::mesh::Node );
          for (stk::mesh::PairIterRelation::iterator inr=nrel.begin(); inr!=nrel.end(); ++inr)
          {
            row.insert( static_cast<int>( inr->entity()->identifier() ) );
          }
        }
      }
    }
  }

  const Epetra_Map & nodemap = NodeRowMap();
  std::vector<int> numIndicesPerRow( graph.size() );
  for (std::map<int,std::set<int> >::iterator i=graph.begin();
       i!=graph.end();
       ++i)
  {
    int gid = i->first;
    int lid = nodemap.LID(gid);
    if (lid<0)
      dserror( "illegal lid" );
    numIndicesPerRow[lid] = i->second.size();
  }

  Teuchos::RCP<Epetra_CrsGraph> g = Teuchos::rcp(new Epetra_CrsGraph(Copy,nodemap,&numIndicesPerRow[0],true));
  for (std::map<int,std::set<int> >::iterator i=graph.begin();
       i!=graph.end();
       ++i)
  {
    int gid = i->first;
    int lid = nodemap.LID(gid);
    if (lid<0)
      dserror( "illegal lid" );
    std::vector<int> row;
    row.reserve(i->second.size());
    row.assign(i->second.begin(),i->second.end());
    g->InsertGlobalIndices(gid,row.size(),&row[0]);
  }
  g->FillComplete();
  return g;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> STK::Discretization::ElementGraph()
{
  STK::Mesh & mesh = GetMesh();
  std::map<int,std::set<int> > graph;

  std::set<int> eids;

  stk::mesh::Part & active = mesh.ActivePart();

  // Get all owned nodes
  stk::mesh::Selector ownedselector = mesh.OwnedPart() & mesh.ActivePart();
  const std::vector<stk::mesh::Bucket*> & buckets = mesh.BulkData().buckets( stk::mesh::Element );

  // loop buckets and entries
  for (std::vector<stk::mesh::Bucket*>::const_iterator i=buckets.begin();
       i!=buckets.end();
       ++i)
  {
    stk::mesh::Bucket & bucket = **i;
    if ( ownedselector( bucket ) )
    {
      for (stk::mesh::Bucket::iterator ib=bucket.begin();
           ib!=bucket.end();
           ++ib)
      {
        stk::mesh::Entity & e = *ib;
        int id = static_cast<int>( e.identifier() );
        std::set<int> & row = graph[id];

        eids.insert( id );

        // all elements connected
        // we expect proper ghosting here!
        stk::mesh::PairIterRelation erel = e.relations( stk::mesh::Node );
        for (stk::mesh::PairIterRelation::iterator ier=erel.begin(); ier!=erel.end(); ++ier)
        {
          // all nodes connected to the element
          stk::mesh::PairIterRelation nrel = ier->entity()->relations( stk::mesh::Element );
          for (stk::mesh::PairIterRelation::iterator inr=nrel.begin(); inr!=nrel.end(); ++inr)
          {
            stk::mesh::Entity & e = *inr->entity();
            if ( has_superset( e.bucket(), active ) )
            {
              row.insert( static_cast<int>( e.identifier() ) );
            }
          }
        }
      }
    }
  }

  std::vector<int> veids;
  veids.reserve( eids.size() );
  veids.assign( eids.begin(), eids.end() );

  Epetra_Map elementmap( -1,veids.size(),&veids[0],0,Comm() );

  std::vector<int> numIndicesPerRow( graph.size() );
  for (std::map<int,std::set<int> >::iterator i=graph.begin();
       i!=graph.end();
       ++i)
  {
    int gid = i->first;
    int lid = elementmap.LID(gid);
    if (lid<0)
      dserror( "illegal lid" );
    numIndicesPerRow[lid] = i->second.size();
  }

  Teuchos::RCP<Epetra_CrsGraph> g = Teuchos::rcp(new Epetra_CrsGraph(Copy,elementmap,&numIndicesPerRow[0],true));
  for (std::map<int,std::set<int> >::iterator i=graph.begin();
       i!=graph.end();
       ++i)
  {
    int gid = i->first;
    int lid = elementmap.LID(gid);
    if (lid<0)
      dserror( "illegal lid" );
    std::vector<int> row;
    row.reserve(i->second.size());
    row.assign(i->second.begin(),i->second.end());
    g->InsertGlobalIndices(gid,row.size(),&row[0]);
  }
  g->FillComplete();
  return g;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::ScatterFieldData( const Epetra_Vector & v, const std::vector<stk::mesh::FieldBase*> & fields )
{
  const Epetra_Map & rowmap = DofRowMap();

  stk::mesh::BulkData & bulk = GetMesh().BulkData();
  stk::mesh::Part & owned    = GetMesh().OwnedPart();

  // loop all owned nodes and copy dof values from Epetra_Vector to given fields

  const std::vector<stk::mesh::Bucket*> & nodes = bulk.buckets( stk::mesh::Node );
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

        std::vector<int> dofs;
        Dof( n.key(), dofs );

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
            if ( lid < 0 )
              dserror( "illegal lid" );
            data[l] = v[lid];
          }
          dofcount += field_data_size;
        }
      }
    }
  }

  // communicate field data from owner to ghosts

  CommunicateFieldGhosting( bulk, fields );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::Discretization::GatherFieldData( const std::vector<stk::mesh::FieldBase*> & fields, Epetra_Vector & v )
{
  const Epetra_Map & rowmap = DofRowMap();

  stk::mesh::BulkData & bulk = GetMesh().BulkData();
  stk::mesh::Part & owned    = GetMesh().OwnedPart();

  // loop all owned nodes and copy dof values from Epetra_Vector to given fields

  const std::vector<stk::mesh::Bucket*> & nodes = bulk.buckets( stk::mesh::Node );
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

        std::vector<int> dofs;
        Dof( n.key(), dofs );

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
            if ( lid < 0 )
              dserror( "illegal lid" );
            v[lid] = data[l];
          }
          dofcount += field_data_size;
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
stk::mesh::Selector STK::Discretization::DirichletSelector()
{
  STK::Mesh & mesh = GetMesh();
  stk::mesh::MetaData & meta = mesh.MetaData();

  stk::mesh::Selector select;

  static DRT::Condition::ConditionType dirichlet_types[] = {
    DRT::Condition::PointDirichlet,
    DRT::Condition::LineDirichlet,
    DRT::Condition::SurfaceDirichlet,
    DRT::Condition::VolumeDirichlet,
    DRT::Condition::none
  };

  for ( DRT::Condition::ConditionType * condtype = dirichlet_types;
        *condtype!=DRT::Condition::none;
        ++condtype )
  {
    const std::map<int, Teuchos::ParameterList> * dirichlet = Condition( *condtype );
    if ( dirichlet!=NULL )
    {
      for ( std::map<int, Teuchos::ParameterList>::const_iterator i=dirichlet->begin();
            i!=dirichlet->end();
            ++i )
      {
        const Teuchos::ParameterList & condflags = i->second;
        std::string name = condflags.get<std::string>( "Name" );

        stk::mesh::Part * condpart = meta.get_part( name );
        if ( condpart==NULL )
          dserror( "No condition part '%s'", name.c_str() );

        select |= *condpart;
      }
    }
  }

  return select;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Parameter* STK::Discretization::MaterialParameter( const stk::mesh::Bucket & bucket )
{
  for ( std::map<stk::mesh::Part*, MAT::PAR::Parameter*>::iterator i=matpar_.begin();
        i!=matpar_.end();
        ++i )
  {
    stk::mesh::Part & part = * i->first;
    if ( has_superset( bucket, part ) )
    {
      return i->second;
    }
  }
  return NULL;
}

#endif
