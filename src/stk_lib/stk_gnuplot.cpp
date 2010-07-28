#ifdef STKADAPTIVE

#include <fstream>

#include "stk_gnuplot.H"
#include "stk_mesh.H"

void STK::GnuplotDump( STK::Mesh & mesh, std::string filename, stk::mesh::FieldBase & field )
{
  std::ofstream out( filename.c_str() );

  stk::mesh::BulkData & bulk = mesh.BulkData();
  stk::mesh::Part & active = mesh.ActivePart();

  const std::vector<stk::mesh::Bucket*> & element_buckets = bulk.buckets( stk::mesh::Element );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=element_buckets.begin();
        i!=element_buckets.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;

    if ( has_superset( bucket, active ) )
    {

      for ( stk::mesh::Bucket::iterator j=bucket.begin();
            j!=bucket.end();
            ++j )
      {
        stk::mesh::Entity & e = *j;

        stk::mesh::PairIterRelation rel = e.relations( stk::mesh::Node );
        for ( ; not rel.empty(); ++rel )
        {
          stk::mesh::Entity & n = *rel->entity();
          double * const coord = stk::mesh::field_data( mesh.Coordinates() , n );
          double * const f = reinterpret_cast<double*>( stk::mesh::field_data( field , n ) );

          out << coord[0] << " "
              << coord[1] << " "
              << f[0] << " "
              << f[1] << " "
              << "\n";
        }

        rel = e.relations( stk::mesh::Node );
        if ( not rel.empty() )
        {
          stk::mesh::Entity & n = *rel->entity();
          double * const coord = stk::mesh::field_data( mesh.Coordinates() , n );
          double * const f = reinterpret_cast<double*>( stk::mesh::field_data( field , n ) );

          out << coord[0] << " "
              << coord[1] << " "
              << f[0] << " "
              << f[1] << " "
              << "\n";
        }

        out << "\n\n";
      }
    }
  }
}

#endif
