
#ifdef STKADAPTIVE

#include "stk_algorithm.H"
#include "stk_discret.H"
#include "stk_fei.H"
#include "stk_field_comm.H"
#include "stk_mesh.H"
#include "stk_refine.H"
#include "stk_unrefine.H"

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
  Refine( refine );
  Unrefine( unrefine );
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
void STK::Discretization::CreateState()
{
  state_ = Teuchos::rcp( new STK::FEI::DiscretizationState( *this ) );

  std::vector<stk::mesh::FieldBase*> fields;

  algo_->collect_unknowns( fields );

  state_->Setup( fields );

  algo_->notify_state_changed();
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
MAT::PAR::Parameter* STK::Discretization::MaterialParameter( stk::mesh::Bucket & bucket )
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
