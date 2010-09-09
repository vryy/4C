
#ifdef STKADAPTIVE

#include "stk_dofset.H"
#include "stk_discret.H"
#include "../stk_refine/stk_mesh.H"

#include <Isorropia_EpetraOrderer.hpp>
#include <Epetra_Import.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_condition.H"

#include "../linalg/linalg_sparseoperator.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::DofSet::Setup( STK::Discretization& dis, const std::vector<stk::mesh::FieldBase*>& fields )
{
  stk::mesh::BulkData & bulk = dis.GetMesh().BulkData();

  const std::vector<stk::mesh::Bucket*>& buckets = bulk.buckets( stk::mesh::Node );

  stk::mesh::Selector activeselector = dis.GetMesh().ActivePart();
  stk::mesh::Selector ownedselector = dis.GetMesh().OwnedPart() & dis.GetMesh().ActivePart();

  Teuchos::RCP<Epetra_CrsGraph> graph = dis.NodeGraph();

  Epetra_BlockMap noderowmap = graph->RowMap();
  Epetra_BlockMap nodecolmap = graph->ColMap();

  std::vector<int> stknodecolvector;

  for ( std::size_t i=0; i<buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *buckets[i];
    if ( activeselector( bucket ) )
    {
      for (stk::mesh::Bucket::iterator ib=bucket.begin();
           ib!=bucket.end();
           ++ib)
      {
        stk::mesh::Entity & n = *ib;

        int gid = static_cast<int>( n.identifier() );
        stknodecolvector.push_back( gid );
      }
    }
  }

  // the stk ghost distribution might be different from the graph column map
  Epetra_Map stknodecolmap( -1, stknodecolvector.size(), &stknodecolvector[0], 0, graph->Comm() );

  numdfcolnodes_ = Teuchos::rcp( new Epetra_IntVector( stknodecolmap ) );
  idxcolnodes_   = Teuchos::rcp( new Epetra_IntVector( stknodecolmap ) );

  for ( std::size_t i=0; i<buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *buckets[i];
    if ( activeselector( bucket ) )
    {
      for ( std::vector<stk::mesh::FieldBase*>::const_iterator j=fields.begin();
            j!=fields.end();
            ++j )
      {
        stk::mesh::FieldBase & field = **j;

        int field_data_size = stk::mesh::field_data_size( field, bucket ) / sizeof( double );
        if ( field_data_size == 0 ) continue;

        for (stk::mesh::Bucket::iterator ib=bucket.begin();
             ib!=bucket.end();
             ++ib)
        {
          stk::mesh::Entity & n = *ib;

          int gid = static_cast<int>( n.identifier() );
          int lid = stknodecolmap.LID( gid );
          if ( lid<0 )
            dserror( "map failure" );

          ( *numdfcolnodes_ )[lid] += field_data_size;
        }
      }
    }
  }

  int max_count = numdfcolnodes_->MaxValue();

  Epetra_IntVector colpermutation( stknodecolmap );

#if 1

  Teuchos::ParameterList params;
  Teuchos::ParameterList& sublist = params.sublist("ZOLTAN");

  //params.set("PARTITIONING METHOD", "GRAPH");
  //sublist.set("LB_APPROACH", "ORDER");
  //sublist.set("GRAPH_PACKAGE", "PARMETIS");
  //
  // metis does not work too well for reordering ... use scotch instead
  //sublist.set( "ORDER_METHOD", "PARMETIS" );
  //sublist.set( "ORDER_METHOD", "METIS" );
  sublist.set( "ORDER_METHOD", "PTSCOTCH" );
  //sublist.set( "ORDER_METHOD", "SCOTCH" );
  //sublist.set( "USE_ORDER_INFO", "1" );
  //sublist.set( "DEBUG_LEVEL", "11" );

  Isorropia::Epetra::Orderer orderer( Teuchos::rcp_implicit_cast<const Epetra_CrsGraph>( graph ), params );

  int size;
  const int *array;
  orderer.extractPermutationView( size, array );

  if ( size != noderowmap.NumMyElements() )
    dserror( "confusion!" );

  Epetra_IntVector rowpermutation( noderowmap );
  std::copy( array, array+size, &rowpermutation[0] );

  Epetra_Import nodemapimporter( stknodecolmap, noderowmap );

  int err = colpermutation.Import( rowpermutation, nodemapimporter, Insert );
  if ( err )
    dserror( "Import failed: err=%d", err );

#else

  // No reordering. For debugging.

  for ( int i=0; i<stknodecolmap.NumMyElements(); ++i )
  {
    int gid = stknodecolmap.GID( i );
    colpermutation[i] = gid;
  }

#endif

#if 0
  {
    std::stringstream fs;
    fs << "orig" << dis.GetMesh().parallel_rank() << ".graph";
    std::ofstream f( fs.str().c_str() );

    for ( int i=0; i<noderowmap.NumMyElements(); ++i )
    {
      int rgid = noderowmap.GID( i );
      int NumIndices;
      int *Indices;
      int err = graph->ExtractMyRowView(i, NumIndices, Indices);
      if ( err!=0 )
        dserror( "nope! err=%d", err );
      for ( int j=0; j<NumIndices; ++j )
      {
        int cgid = nodecolmap.GID( Indices[j] );
        f << cgid << " " << rgid << "\n";
      }
    }

    std::stringstream ns;
    ns << "new" << dis.GetMesh().parallel_rank() << ".graph";
    std::ofstream n( ns.str().c_str() );
    for ( int i=0; i<noderowmap.NumMyElements(); ++i )
    {
      int rgid = noderowmap.GID( i );
      int NumIndices;
      int *Indices;
      int err = graph->ExtractMyRowView(i, NumIndices, Indices);
      if ( err!=0 )
        dserror( "nope! err=%d", err );

      int stk_rlid = stknodecolmap.LID( rgid );

      for ( int j=0; j<NumIndices; ++j )
      {
        int cgid = nodecolmap.GID( Indices[j] );
        int stk_clid = stknodecolmap.LID( cgid );
        n << colpermutation[stk_clid] << " " << colpermutation[stk_rlid] << "\n";
      }
    }
  }
#endif

  std::map<int, stk::mesh::Entity *> nodes;

  for ( std::size_t i=0; i<buckets.size(); ++i )
  {
    stk::mesh::Bucket & bucket = *buckets[i];
    if ( activeselector( bucket ) )
    {
      for (stk::mesh::Bucket::iterator ib=bucket.begin();
           ib!=bucket.end();
           ++ib)
      {
        stk::mesh::Entity & n = *ib;

        int gid = static_cast<int>( n.identifier() );
        int lid = stknodecolmap.LID( gid );
        if ( lid<0 )
          dserror( "map failure" );

        int permid = colpermutation[lid];
        nodes[permid] = &n;
      }
    }
  }

  std::vector<int> dofrowmap;
  std::vector<int> dofcolmap;

  dofrowmap.reserve( noderowmap.NumMyElements()*max_count );
  dofcolmap.reserve( stknodecolmap.NumMyElements()*max_count );

  for ( std::map<int, stk::mesh::Entity *>::iterator i=nodes.begin();
        i!=nodes.end();
        ++i )
  {
    int permid = i->first;
    stk::mesh::Entity & n = * i->second;
    int gid = static_cast<int>( n.identifier() );
    int lid = stknodecolmap.LID( gid );
    if ( lid<0 )
      dserror( "map failure" );

    ( *idxcolnodes_ )[lid] = dofcolmap.size();

    bool isown = noderowmap.MyGID( gid );

    // scale new permutation to get independent dof vales
    int dofbase = permid*max_count;
    int base = 0;

    for ( std::vector<stk::mesh::FieldBase*>::const_iterator j=fields.begin();
          j!=fields.end();
          ++j )
    {
      stk::mesh::FieldBase & field = **j;

      int field_data_size = stk::mesh::field_data_size( field, n.bucket() );
      field_data_size /= sizeof( double );

      for ( int k=0; k<field_data_size; ++k )
      {
        // Datenstrukturen um zwischen beiden zu vermitteln

        int dofgid = dofbase + base + k;

        DofName & dn = reverse_[dofgid];
        dn.key = n.key();
        dn.field = &field;
        dn.pos = k;

        if ( isown )
        {
          dofrowmap.push_back( dofgid );
        }
        dofcolmap.push_back( dofgid );
      }
      base += field_data_size;
    }
  }

  dofrowmap_ = Teuchos::rcp( new Epetra_Map( -1, dofrowmap.size(), &dofrowmap[0], 0, graph->Comm() ) );
  dofcolmap_ = Teuchos::rcp( new Epetra_Map( -1, dofcolmap.size(), &dofcolmap[0], 0, graph->Comm() ) );
}

#endif
