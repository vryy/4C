#ifdef STKADAPTIVE

#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>

#include "stk_mesh.H"


STK::MetaMesh::MetaMesh()
  : mesh_meta_data_( stk::mesh::fem_entity_rank_names() )
{
  using namespace stk;

  active_           = & mesh_meta_data_.declare_part( "active" );
  hanging_          = & mesh_meta_data_.declare_part( "hanging"  , mesh::Constraint );

  // element types

  quad4_    = & mesh_meta_data_.declare_part("quad4"       , mesh::Element );
  quad8_    = & mesh_meta_data_.declare_part("quad8"       , mesh::Element );
  quad9_    = & mesh_meta_data_.declare_part("quad9"       , mesh::Element );
  tri3_     = & mesh_meta_data_.declare_part("tri3"        , mesh::Element );
  tri6_     = & mesh_meta_data_.declare_part("tri6"        , mesh::Element );
  hex8_     = & mesh_meta_data_.declare_part("hex8"        , mesh::Element );
  hex20_    = & mesh_meta_data_.declare_part("hex20"       , mesh::Element );
  hex27_    = & mesh_meta_data_.declare_part("hex27"       , mesh::Element );
  tet4_     = & mesh_meta_data_.declare_part("tet4"        , mesh::Element );
  tet10_    = & mesh_meta_data_.declare_part("tet10"       , mesh::Element );
  wedge6_   = & mesh_meta_data_.declare_part("wedge6"      , mesh::Element );
  wedge15_  = & mesh_meta_data_.declare_part("wedge15"     , mesh::Element );
  pyramid5_ = & mesh_meta_data_.declare_part("pyramid5"    , mesh::Element );
  line2_    = & mesh_meta_data_.declare_part("line2"       , mesh::Element );
  line3_    = & mesh_meta_data_.declare_part("line3"       , mesh::Element );

  //face_quad_        = & mesh_meta_data_.declare_part("face_quad" , mesh::Face );
  //edge_line_        = & mesh_meta_data_.declare_part("edge_line" , mesh::Edge );

  mesh::set_cell_topology< shards::Quadrilateral<4> >( *quad4_ );
  mesh::set_cell_topology< shards::Quadrilateral<8> >( *quad8_ );
  mesh::set_cell_topology< shards::Quadrilateral<9> >( *quad9_ );
  mesh::set_cell_topology< shards::Triangle<3>      >( *tri3_ );
  mesh::set_cell_topology< shards::Triangle<6>      >( *tri6_ );
  mesh::set_cell_topology< shards::Hexahedron<8>    >( *hex8_ );
  mesh::set_cell_topology< shards::Hexahedron<20>   >( *hex20_ );
  mesh::set_cell_topology< shards::Hexahedron<27>   >( *hex27_ );
  mesh::set_cell_topology< shards::Tetrahedron<4>   >( *tet4_ );
  mesh::set_cell_topology< shards::Tetrahedron<10>  >( *tet10_ );
  mesh::set_cell_topology< shards::Wedge<6>         >( *wedge6_ );
  mesh::set_cell_topology< shards::Wedge<15>        >( *wedge15_ );
  mesh::set_cell_topology< shards::Pyramid<5>       >( *pyramid5_ );
  mesh::set_cell_topology< shards::Line<2>          >( *line2_ );
  mesh::set_cell_topology< shards::Line<3>          >( *line3_ );

  //mesh::set_cell_topology< shards::Quadrilateral<4>       >( *face_quad_ );
  //mesh::set_cell_topology< shards::Line<2>                >( *edge_line_ );

  coordinates_field_ = & mesh::declare_vector_field_on_all_nodes( mesh_meta_data_, "coordinates", 3 );

  element_attr_ = & mesh_meta_data_.declare_field<AttributeField>( std::string("element_attribute") );

  file_schema_ = new phdmesh::exodus::FileSchema( mesh_meta_data_ , *coordinates_field_ , *element_attr_ , *active_ );

  file_schema_->declare_part( *quad4_    , 1 ); //
                                                // quad4_->mesh_meta_data_ordinal() ???
  file_schema_->declare_part( *quad8_    , 2 );
  file_schema_->declare_part( *quad9_    , 3 );
  file_schema_->declare_part( *tri3_     , 4 );
  file_schema_->declare_part( *tri6_     , 5 );
  file_schema_->declare_part( *hex8_     , 6 );
  file_schema_->declare_part( *hex20_    , 7 );
  file_schema_->declare_part( *hex27_    , 8 );
  file_schema_->declare_part( *tet4_     , 9 );
  file_schema_->declare_part( *tet10_    , 10 );
  file_schema_->declare_part( *wedge6_   , 11 );
  file_schema_->declare_part( *wedge15_  , 12 );
  file_schema_->declare_part( *pyramid5_ , 13 );
  file_schema_->declare_part( *line2_    , 14 );
  file_schema_->declare_part( *line3_    , 15 );
}


void STK::MetaMesh::Commit()
{
  mesh_meta_data_.commit();
}


STK::Mesh::Mesh( MetaMesh & meta, stk::ParallelMachine pm )
  : meta_( meta ),
    mesh_bulk_data_( meta.MetaData(), pm, 200 )
{
}


void STK::Mesh::Output( std::string filename, std::string title )
{
  std::vector< const stk::mesh::FieldBase * > out_fields ;
  Output( filename, title, out_fields );
}



void STK::Mesh::Output( std::string filename, std::string title, std::vector< const stk::mesh::FieldBase * > & out_fields )
{
  Teuchos::RCP<phdmesh::exodus::FileOutput> exo = OutputContext( filename, title, out_fields );
  exo->write( 0, 0 );
}


Teuchos::RCP<phdmesh::exodus::FileOutput> STK::Mesh::OutputContext( std::string filename, std::string title, std::vector< const stk::mesh::FieldBase * > & out_fields )
{
  meta_.FileSchema().assign_indices( mesh_bulk_data_ );

  int flags[ stk::mesh::EntityRankEnd ] = { 1 , 0 , 0 , 1 , 0 , 0 };
  Teuchos::RCP<phdmesh::exodus::FileOutput> exo =
    Teuchos::rcp( new phdmesh::exodus::FileOutput( meta_.FileSchema() , mesh_bulk_data_ ,
                                                   meta_.Active() ,
                                                   filename , title ,
                                                   false , out_fields , flags ) );
  return exo;
}


void STK::Mesh::Modify()
{
  mesh_bulk_data_.modification_begin();
}


void STK::Mesh::Commit()
{
  mesh_bulk_data_.modification_end();
}

void STK::Mesh::Print(std::ostream& stream)
{
  {
    const std::vector<stk::mesh::Bucket*> & nodes = mesh_bulk_data_.buckets( stk::mesh::Node );
    for (std::vector<stk::mesh::Bucket*>::const_iterator i=nodes.begin();
         i!=nodes.end();
         ++i)
    {
      //stream << (**i) << " (" << (*i)->size() << ")" <<"\n";
      for (stk::mesh::Bucket::iterator j=(*i)->begin();
           j!=(*i)->end();
           ++j)
      {
        double * const coord = stk::mesh::field_data( Coordinates() , *j );
        stream << coord[0] << " " << coord[1] << " " << coord[2]
               << " # ";
        stk::mesh::print_entity_key( stream , mesh_bulk_data_.mesh_meta_data() , j->key() );
        stream << " (" << j->owner_rank() << ")\n";
      }
    }
  }

  for ( int entity_type=0; entity_type<stk::mesh::EntityRankEnd; ++entity_type )
  {
    const std::vector<stk::mesh::Bucket*> & nodes = mesh_bulk_data_.buckets( entity_type );
    for (std::vector<stk::mesh::Bucket*>::const_iterator i=nodes.begin();
         i!=nodes.end();
         ++i)
    {
      stream << (**i) << " (" << (*i)->size() << ")" <<"\n";
      for (stk::mesh::Bucket::iterator j=(*i)->begin();
           j!=(*i)->end();
           ++j)
      {
        stk::mesh::print_entity_key( stream , mesh_bulk_data_.mesh_meta_data() , j->key() );
        stream << " (" << j->owner_rank() << ") : ";
        for (stk::mesh::PairIterRelation r=j->relations(); not r.empty(); ++r)
        {
          stream << (*r) << ", ";
        }
        for ( stk::mesh::PairIterEntityComm comm = j->comm(); not comm.empty(); ++comm )
        {
          stream << "(" << comm->ghost_id << "," << comm->proc << "), ";
        }
        stream << "\n";
      }
      stream << "\n";
    }
  }
  stream << std::flush;
}

void STK::Mesh::Dump(std::string filename)
{
  std::ostringstream str;
  str << filename << "-" << parallel_rank() << ".dat";
  std::ofstream f(str.str().c_str());
  Print(f);
}


void STK::Mesh::OutputCounters( std::string name, std::vector<stk::mesh::EntityRank> & count )
{
  unsigned myrank = parallel_rank();
  unsigned nprocs = parallel_size();
  unsigned size = count.size();

  std::vector<stk::mesh::EntityRank> localcount( size*( nprocs+1 ), 0 );
  std::vector<stk::mesh::EntityRank> globalcount( size*( nprocs+1 ), 0 );
  std::copy( count.begin(), count.end(), &localcount[size*myrank] );

  stk::all_reduce_sum( parallel(), &localcount[0], &globalcount[0], localcount.size() );

  for ( unsigned p = 0 ; p < nprocs ; ++p )
  {
    for ( unsigned i=0; i<size; ++i )
    {
      globalcount[size*nprocs+i] += globalcount[size*p+i];
    }
  }

  if ( myrank == 0 )
  {
    std::cout << std::setw( 15 ) << name << " : ";
    for ( unsigned p = 0 ; p < nprocs+1 ; ++p )
    {
      for ( stk::mesh::EntityRank* ptr=&globalcount[p*size];
            ptr!=&globalcount[( p+1 )*size];
            ++ptr )
      {
        std::cout << std::setw( 3 ) << ( *ptr ) << " ";
      }
      if ( p<nprocs )
        std::cout << "| ";
    }
    std::cout << "\n";
  }
}


void STK::Mesh::Statistics()
{
  std::vector<stk::mesh::EntityRank> count;

  if ( parallel_rank() == 0 )
    std::cout << "\n";

  stk::mesh::Selector ownedselector = OwnedPart();
  count_entities(ownedselector,mesh_bulk_data_,count);
  OutputCounters( "owns", count );

  stk::mesh::Selector sharedselector = SharedPart();
  count_entities(sharedselector,mesh_bulk_data_,count);
  OutputCounters( "shares", count );

  stk::mesh::Selector activeownedselector = OwnedPart() & ActivePart();
  count_entities(activeownedselector,mesh_bulk_data_,count);
  OutputCounters( "active owns", count );

  stk::mesh::Selector activesharedselector = SharedPart() & ActivePart();
  count_entities(activesharedselector,mesh_bulk_data_,count);
  OutputCounters( "active shares", count );

  const std::vector<stk::mesh::Entity*> & entity_comm = mesh_bulk_data_.entity_comm();

  std::fill(count.begin(),count.end(),0);
  for ( std::vector<stk::mesh::Entity*>::const_iterator i = entity_comm.begin();
        i != entity_comm.end();
        ++i )
  {
    if ( stk::mesh::in_shared( **i ) )
    {
      count[( *i )->entity_rank()] += 1;
    }
  }
  OutputCounters( "shared", count );

  std::fill(count.begin(),count.end(),0);
  for ( std::vector<stk::mesh::Entity*>::const_iterator i = entity_comm.begin();
        i != entity_comm.end();
        ++i )
  {
    if ( stk::mesh::in_send_ghost( **i ) )
    {
      count[( *i )->entity_rank()] += 1;
    }
  }
  OutputCounters( "send", count );

  std::fill(count.begin(),count.end(),0);
  for ( std::vector<stk::mesh::Entity*>::const_iterator i = entity_comm.begin();
        i != entity_comm.end();
        ++i )
  {
    if ( stk::mesh::in_receive_ghost( **i ) )
    {
      count[( *i )->entity_rank()] += 1;
    }
  }
  OutputCounters( "recv", count );
}

#endif
