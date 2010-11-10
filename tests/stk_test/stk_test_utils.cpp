
#include "stk_test_utils.H"

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

void Print( stk::mesh::BulkData & bulk_data, std::ostream & stream )
{
  {
    const std::vector<stk::mesh::Bucket*> & nodes = bulk_data.buckets( stk::mesh::Node );
    for (std::vector<stk::mesh::Bucket*>::const_iterator i=nodes.begin();
         i!=nodes.end();
         ++i)
    {
      for (stk::mesh::Bucket::iterator j=(*i)->begin();
           j!=(*i)->end();
           ++j)
      {
        double * const coord = stk::mesh::field_data( *bulk_data.mesh_meta_data().get_field<stk::mesh::VectorField>( "coordinates" ) , *j );
        stream << coord[0] << " " << coord[1] << " " << coord[2]
               << " # ";
        stk::mesh::print_entity_key( stream , bulk_data.mesh_meta_data() , j->key() );
        stream << " (" << j->owner_rank() << ")\n";
      }
    }
  }

  stream << "\n";

  for ( int entity_type=0; entity_type<stk::mesh::EntityRankEnd; ++entity_type )
  {
    const std::vector<stk::mesh::Bucket*> & nodes = bulk_data.buckets( entity_type );
    for (std::vector<stk::mesh::Bucket*>::const_iterator i=nodes.begin();
         i!=nodes.end();
         ++i)
    {
      stream << (**i) << " (" << (*i)->size() << ")" <<"\n";
      for (stk::mesh::Bucket::iterator j=(*i)->begin();
           j!=(*i)->end();
           ++j)
      {
        stk::mesh::print_entity_key( stream , bulk_data.mesh_meta_data() , j->key() );
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
}


void Dump( stk::mesh::BulkData & bulk_data, std::string filename )
{
  std::ostringstream str;
  str << filename << "-" << bulk_data.parallel_rank() << ".dat";
  std::ofstream f( str.str().c_str() );
  Print( bulk_data, f );
}


void ShowSync( stk::mesh::BulkData & bulk_data, unsigned nx, unsigned ny )
{
  for ( unsigned p=0; p<bulk_data.parallel_size(); ++p )
  {
    if ( p==bulk_data.parallel_rank() )
    {
      Show( bulk_data, nx, ny );
      std::cout.flush();
    }
    MPI_Barrier( bulk_data.parallel() );
  }
}


// far from perfect but good enough for now
void Show( stk::mesh::BulkData & bulk_data, unsigned nx, unsigned ny )
{
  unsigned p_rank = bulk_data.parallel_rank();
  std::cout << "P" << p_rank << ":\n";
  int iy = ny-1;
  for ( unsigned ix = 0; ix < nx; ++ix )
  {
    stk::mesh::EntityId elem = 1 + ix + iy * nx ;
    stk::mesh::Entity * e = bulk_data.get_entity( stk::mesh::Element, elem );
    bool firstcol = ix==0 or bulk_data.get_entity( stk::mesh::Element, elem-1 )==NULL or bulk_data.get_entity( stk::mesh::Element, elem-1 )->bucket().capacity()==0;
    if ( e!=NULL and e->bucket().capacity() > 0 )
    {
      stk::mesh::PairIterRelation rel = e->relations( stk::mesh::Node );
      if ( firstcol )
      {
        stk::mesh::Entity * n = rel[3].entity();
        if ( n->owner_rank()==p_rank )
        {
          std::cout << "-" << std::setw( 2 ) << n->identifier() << "-";
        }
        else
        {
          std::cout << "(" << std::setw( 2 ) << n->identifier() << ")";
        }
      }
      std::cout << "---";
      stk::mesh::Entity * n = rel[2].entity();
      if ( n->owner_rank()==p_rank )
      {
        std::cout << "-" << std::setw( 2 ) << n->identifier() << "-";
      }
      else
      {
        std::cout << "(" << std::setw( 2 ) << n->identifier() << ")";
      }
    }
    else
    {
      if ( firstcol )
        std::cout << "    ";
      std::cout << "   ";
    }
  }
  std::cout << "\n";
  for ( ; iy >= 0; --iy )
  {
    std::cout << "  ";
    for ( unsigned ix = 0; ix < nx; ++ix )
    {
      stk::mesh::EntityId elem = 1 + ix + iy * nx ;
      stk::mesh::Entity * e = bulk_data.get_entity( stk::mesh::Element, elem );
      bool firstcol = ix==0 or bulk_data.get_entity( stk::mesh::Element, elem-1 )==NULL or bulk_data.get_entity( stk::mesh::Element, elem-1 )->bucket().capacity()==0;
      if ( e!=NULL and e->bucket().capacity() > 0 )
      {
        if ( firstcol )
        {
          std::cout << "|";
        }
        std::cout << "      |";
      }
      else
      {
        if ( firstcol ) std::cout << " ";
        std::cout << "      ";
      }
    }
    std::cout << "\n";
    std::cout << "  ";
    for ( unsigned ix = 0; ix < nx; ++ix )
    {
      stk::mesh::EntityId elem = 1 + ix + iy * nx ;
      stk::mesh::Entity * e = bulk_data.get_entity( stk::mesh::Element, elem );
      bool firstcol = ix==0 or bulk_data.get_entity( stk::mesh::Element, elem-1 )==NULL or bulk_data.get_entity( stk::mesh::Element, elem-1 )->bucket().capacity()==0;
      if ( e!=NULL and e->bucket().capacity() > 0 )
      {
        if ( firstcol )
        {
          std::cout << "|";
        }
        if ( e->owner_rank()==p_rank )
        {
          std::cout << "  " << std::setw( 2 ) << elem << "  |";
        }
        else
        {
          std::cout << " (" << std::setw( 2 ) << elem << ") |";
        }
      }
      else
      {
        if ( firstcol ) std::cout << " ";
        std::cout << "      ";
      }
    }
    std::cout << "\n";
    std::cout << "  ";
    for ( unsigned ix = 0; ix < nx; ++ix )
    {
      stk::mesh::EntityId elem = 1 + ix + iy * nx ;
      stk::mesh::Entity * e = bulk_data.get_entity( stk::mesh::Element, elem );
      bool firstcol = ix==0 or bulk_data.get_entity( stk::mesh::Element, elem-1 )==NULL or bulk_data.get_entity( stk::mesh::Element, elem-1 )->bucket().capacity()==0;
      if ( e!=NULL and e->bucket().capacity() > 0 )
      {
        if ( firstcol )
        {
          std::cout << "|";
        }
        std::cout << "      |";
      }
      else
      {
        if ( firstcol ) std::cout << " ";
        std::cout << "      ";
      }
    }
    std::cout << "\n";
    for ( unsigned ix = 0; ix < nx; ++ix )
    {
      stk::mesh::EntityId elem = 1 + ix + iy * nx ;
      stk::mesh::Entity * e = bulk_data.get_entity( stk::mesh::Element, elem );
      bool firstcol = ix==0 or bulk_data.get_entity( stk::mesh::Element, elem-1 )==NULL or bulk_data.get_entity( stk::mesh::Element, elem-1 )->bucket().capacity()==0;
      if ( e!=NULL and e->bucket().capacity() > 0 )
      {
        stk::mesh::PairIterRelation rel = e->relations( stk::mesh::Node );
        if ( firstcol )
        {
          stk::mesh::Entity * n = rel[0].entity();
          if ( n->owner_rank()==p_rank )
          {
            std::cout << "-" << std::setw( 2 ) << n->identifier() << "-";
          }
          else
          {
            std::cout << "(" << std::setw( 2 ) << n->identifier() << ")";
          }
        }
        std::cout << "---";
        stk::mesh::Entity * n = rel[1].entity();
        if ( n->owner_rank()==p_rank )
        {
          std::cout << "-" << std::setw( 2 ) << n->identifier() << "-";
        }
        else
        {
          std::cout << "(" << std::setw( 2 ) << n->identifier() << ")";
        }
      }
      else
      {
        if ( firstcol )
          std::cout << "    ";
        std::cout << "   ";
      }
    }
    std::cout << "\n";
  }
}
