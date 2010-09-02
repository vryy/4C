
#ifdef STKADAPTIVE

#include "stk_fei.H"
#include "stk_discret.H"
#include "stk_mesh.H"

#include <Isorropia_EpetraOrderer.hpp>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_Permutation.h>
#include <Epetra_Import.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_condition.H"

#include "../linalg/linalg_sparseoperator.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::DofSet::Setup( STK::Discretization& dis, const std::vector<stk::mesh::FieldBase*>& fields )
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::FEI::AssembleStrategy::AssembleStrategy( STK::Discretization & dis,
                                              Teuchos::RCP<const Epetra_Map> dbcmap,
                                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                              Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                              Teuchos::RCP<Epetra_Vector> systemvector1,
                                              Teuchos::RCP<Epetra_Vector> systemvector2,
                                              Teuchos::RCP<Epetra_Vector> systemvector3 )
  : DRT::AssembleStrategy( systemmatrix1,
                           systemmatrix2,
                           systemvector1,
                           systemvector2,
                           systemvector3 ),
    dbcmap_( dbcmap )
{
  // hanging node dof maps
  // extract dof information to be able to assemble later on

  Mesh & mesh = dis.GetMesh();
  stk::mesh::BulkData & bulk_data = mesh.BulkData();
  //stk::mesh::Part & owned = mesh.OwnedPart();

  const std::vector<stk::mesh::Bucket*> & constraints = bulk_data.buckets( stk::mesh::Constraint );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=constraints.begin();
        i!=constraints.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;

    // All hanging nodes are active. Shared and ghosted nodes need to be
    // handled as well.

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

          std::vector<int> hndofs;
          dis.Dof( hn.key(), hndofs );
          for ( unsigned k=0; k!=hndofs.size(); ++k )
          {
            int d = hndofs[k];
            DofData & dd = hndofmap_[d];
            dd.fact = fact;
            std::set<int> & realdofmap = dd.realdofmap;
            for ( stk::mesh::PairIterRelation::iterator r=rel.begin(); r!=rel.end(); ++r )
            {
              stk::mesh::Entity & n = *r->entity();

              std::vector<int> dofs;
              dis.Dof( n.key(), dofs );
              realdofmap.insert( dofs[k] );
            }
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> STK::FEI::AssembleStrategy::MatrixGraph( DRT::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap )
{
  dserror( "no matrix graph based on DRT::Discretization here" );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> STK::FEI::AssembleStrategy::MatrixGraph( STK::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap )
{
  STK::Mesh & mesh = dis.GetMesh();
  const Epetra_Map & dofrowmap = dis.DofRowMap();

  // build graph

  std::map<int, std::set<int> > graph;

  stk::mesh::Selector selector = mesh.OwnedPart() & mesh.ActivePart();
  const std::vector<stk::mesh::Bucket*> & buckets = mesh.BulkData().buckets( stk::mesh::Node );

  // loop buckets and entries
  for (std::vector<stk::mesh::Bucket*>::const_iterator i=buckets.begin();
       i!=buckets.end();
       ++i)
  {
    stk::mesh::Bucket & bucket = **i;
    if ( selector( bucket ) )
    {
      for ( stk::mesh::Bucket::iterator ib=bucket.begin();
            ib!=bucket.end();
            ++ib )
      {
        stk::mesh::Entity & rn = *ib;

        std::vector<int> rdofs;
        dis.Dof( rn.key(), rdofs );
        for ( unsigned rj=0; rj<rdofs.size(); ++rj )
        {
          std::set<int> & rowset = graph[rdofs[rj]];
          if ( not dbcmap_->MyGID( rdofs[rj] ) )
          {
            // non-Dirichlet row

            std::map<int, DofData>::iterator rhndofmap = hndofmap_.find( rdofs[rj] );
            if ( rhndofmap==hndofmap_.end() )
            {
              for ( stk::mesh::PairIterRelation k=rn.relations( stk::mesh::Element );
                    not k.empty();
                    ++k )
              {
                stk::mesh::Entity & e = *k->entity();
                for ( stk::mesh::PairIterRelation l=e.relations( stk::mesh::Node );
                      not l.empty();
                      ++l )
                {
                  stk::mesh::Entity & cn = *l->entity();

                  std::vector<int> cdofs;
                  dis.Dof( cn.key(), cdofs );
                  //rowset.insert( cdofs.begin(), cdofs.end() );
                  for ( unsigned cj=0; cj<cdofs.size(); ++cj )
                  {
                    std::map<int, DofData>::iterator chndofmap = hndofmap_.find( cdofs[cj] );
                    if ( chndofmap==hndofmap_.end() )
                    {
                      rowset.insert( cdofs[cj] );
                    }
                    else
                    {
                      // Instead of the hanging node dof we are connected with all
                      // those real node dofs.
                      DofData & cdd = chndofmap->second;
                      rowset.insert( cdd.realdofmap.begin(), cdd.realdofmap.end() );
                    }
                  }
                }
              }
            }
            else
            {
              // have an innocent (unconnected) diagonal entry in a hanging node row
              rowset.insert( rdofs[rj] );

              DofData & rdd = rhndofmap->second;
              std::set<int> & realdofmap = rdd.realdofmap;

              rowset.insert( realdofmap.begin(), realdofmap.end() );

              // Add my line to all those connected lines, if those lines are no
              // Dirichlet lines. Here we do not have to watch for hanging nodes
              // any more.
              for ( std::set<int>::iterator r=realdofmap.begin(); r!=realdofmap.end(); ++r )
              {
                int rd = *r;
                if ( not dbcmap_->MyGID( rd ) and dofrowmap.MyGID( rd ) )
                {
                  std::set<int> & hnrowset = graph[rd];

                  for ( stk::mesh::PairIterRelation k=rn.relations( stk::mesh::Element );
                        not k.empty();
                        ++k )
                  {
                    stk::mesh::Entity & e = *k->entity();
                    for ( stk::mesh::PairIterRelation l=e.relations( stk::mesh::Node );
                          not l.empty();
                          ++l )
                    {
                      stk::mesh::Entity & cn = *l->entity();

                      std::vector<int> cdofs;
                      dis.Dof( cn.key(), cdofs );
                      for ( unsigned cj=0; cj<cdofs.size(); ++cj )
                      {
                        std::map<int, DofData>::iterator chndofmap = hndofmap_.find( cdofs[cj] );
                        if ( chndofmap==hndofmap_.end() )
                        {
                          hnrowset.insert( cdofs[cj] );
                        }
                        else
                        {
                          // Instead of the hanging node dof we are connected with all
                          // those real node dofs.
                          DofData & cdd = chndofmap->second;
                          hnrowset.insert( cdd.realdofmap.begin(), cdd.realdofmap.end() );
                        }
                      }
                    }
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
    }
  }

  // setup graph with row length and column indices per row

  std::vector<int> sizes( graph.size() );
  for ( std::map<int, std::set<int> >::iterator i=graph.begin(); i!=graph.end(); ++i )
  {
    int gid = i->first;
    int lid = dofrowmap.LID( gid );
    std::set<int> & rowset = i->second;
    unsigned s = rowset.size();
    sizes[lid] = s;
  }

  Teuchos::RCP<Epetra_CrsGraph> crsgraph =
    Teuchos::rcp( new Epetra_CrsGraph( Copy, dofrowmap, &sizes[0], true ) );

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
  }

  crsgraph->FillComplete();

  return crsgraph;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat,
                                          int eid,
                                          const Epetra_SerialDenseMatrix& Aele,
                                          const std::vector<int>& lm,
                                          const std::vector<int>& lmowner)
{
  Assemble( sysmat, eid, Aele, lm, lmowner, lm );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat,
                                          int eid,
                                          const Epetra_SerialDenseMatrix& Aele,
                                          const std::vector<int>& lmrow,
                                          const std::vector<int>& lmrowowner,
                                          const std::vector<int>& lmcol)
{
  Epetra_CrsMatrix * mat = dynamic_cast<Epetra_CrsMatrix*>( &*sysmat.EpetraOperator() );
  if ( mat==NULL )
    dserror( "unsupported SparseOperator" );

  if (not mat->Filled())
    dserror( "not filled" );

  const unsigned lrowdim = lmrow.size();
  const unsigned lcoldim = lmcol.size();

  const int myrank = mat->Comm().MyPID();
  const Epetra_Map& rowmap = mat->RowMap();
  const Epetra_Map& colmap = mat->ColMap();

  // loop rows of local matrix
  for (unsigned lrow=0; lrow<lrowdim; ++lrow)
  {
    const int rgid = lmrow[lrow];

    // check if this row is a Dirichlet row
    if ( dbcmap_->MyGID( rgid ) ) continue;

    // check for hanging node row
    std::map<int, DofData>::iterator rhndofmap = hndofmap_.find( rgid );

    if ( rhndofmap==hndofmap_.end() )
    {
      // check ownership of row
      if (lmrowowner[lrow] == myrank)
      {
        const int rlid = rowmap.LID(rgid);

        for (unsigned lcol=0; lcol<lcoldim; ++lcol)
        {
          double value = Aele(lrow,lcol);
          int cgid = lmcol[lcol];

          // check for hanging node row
          std::map<int, DofData>::iterator chndofmap = hndofmap_.find( cgid );

          if ( chndofmap==hndofmap_.end() )
          {
            int idx = colmap.LID( cgid );
            const int errone = mat->SumIntoMyValues(rlid,1,&value,&idx);
            if (errone)
              dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);
          }
          else
          {
            DofData & cdd = chndofmap->second;
            double fvalue = value*cdd.fact; // linear elements only for now
            for ( std::set<int>::iterator dof=cdd.realdofmap.begin();
                  dof!=cdd.realdofmap.end();
                  ++dof )
            {
              int gid = *dof;
              int idx = colmap.LID( gid );
              const int errone = mat->SumIntoMyValues(rlid,1,&fvalue,&idx);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);
            }
          }
        }
      }
    }
    else
    {
      // do innocent decouled value first

      int rlid = rowmap.LID(rgid);
      double value = 1;

      const int errone = mat->ReplaceMyValues(rlid,1,&value,&rlid);
      if (errone)
        dserror("Epetra_CrsMatrix::ReplaceMyValues returned error code %d",errone);

      // put hanging node row to real node rows

      DofData & rdd = rhndofmap->second;
      std::set<int> & realdofmap = rdd.realdofmap;

      double fvalue = -rdd.fact;

      for ( std::set<int>::iterator c=realdofmap.begin(); c!=realdofmap.end(); ++c )
      {
        int cgid = *c;
        int clid = rowmap.LID(cgid);
        const int errone = mat->ReplaceMyValues(rlid,1,&fvalue,&clid);
        if (errone)
          dserror("Epetra_CrsMatrix::ReplaceMyValues returned error code %d",errone);
      }

      for ( std::set<int>::iterator r=realdofmap.begin(); r!=realdofmap.end(); ++r )
      {
        int rgid = *r;
        if ( not dbcmap_->MyGID( rgid ) and rowmap.MyGID( rgid ) )
        {
          const int rlid = rowmap.LID(rgid);

          for (unsigned lcol=0; lcol<lcoldim; ++lcol)
          {
            double value = Aele(lrow,lcol)*rdd.fact;
            int cgid = lmcol[lcol];

            // check for hanging node row
            std::map<int, DofData>::iterator chndofmap = hndofmap_.find( cgid );

            if ( chndofmap==hndofmap_.end() )
            {
              int idx = colmap.LID( cgid );
              const int errone = mat->SumIntoMyValues(rlid,1,&value,&idx);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);
            }
            else
            {
              DofData & cdd = chndofmap->second;
              double fvalue = value*cdd.fact; // linear elements only for now
              for ( std::set<int>::iterator dof=cdd.realdofmap.begin();
                    dof!=cdd.realdofmap.end();
                    ++dof )
              {
                int gid = *dof;
                int idx = colmap.LID( gid );
                const int errone = mat->SumIntoMyValues(rlid,1,&fvalue,&idx);
                if (errone)
                  dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);
              }
            }
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid)
{
  Epetra_CrsMatrix * mat = dynamic_cast<Epetra_CrsMatrix*>( &*sysmat.EpetraOperator() );
  if ( mat==NULL )
    dserror( "unsupported SparseOperator" );

  dserror( "not supported" );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::AssembleStrategy::Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                                          const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  const unsigned ldim = lm.size();
  if ( ldim!=lmowner.size() or ldim!=static_cast<unsigned>( Vele.Length() ) )
    dserror("Mismatch in dimensions");

  const int myrank = V.Comm().MyPID();

  const Epetra_BlockMap& rowmap = V.Map();

  for (unsigned lrow=0; lrow<ldim; ++lrow)
  {
    if (lmowner[lrow] == myrank)
    {
      int rgid = lm[lrow];
      if ( not rowmap.MyGID( rgid ) )
        dserror("Sparse vector V does not have global row %d", rgid);
      int rlid = rowmap.LID( rgid );

      // check for hanging node row
      std::map<int, DofData>::iterator rhndofmap = hndofmap_.find( rgid );

      if ( rhndofmap==hndofmap_.end() )
      {
        V[rlid] += Vele[lrow];
      }
      else
      {
        V[rlid] = 0.;

        DofData & rdd = rhndofmap->second;
        std::set<int> & realdofmap = rdd.realdofmap;
        double fvalue = Vele[lrow] * rdd.fact;

        for ( std::set<int>::iterator r=realdofmap.begin(); r!=realdofmap.end(); ++r )
        {
          int rgid = *r;
          if ( not dbcmap_->MyGID( rgid ) and rowmap.MyGID( rgid ) )
          {
            const int rlid = rowmap.LID( rgid );
            V[rlid] += fvalue;
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::AssembleStrategy::Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                                          const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  dserror( "not supported" );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::FEI::DiscretizationState::DiscretizationState( STK::Discretization& dis )
  : dis_( dis )
{
  CreateNodeMaps();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::DiscretizationState::Setup( const std::vector<stk::mesh::FieldBase*>& fields )
{
  dofset_.Setup( dis_, fields );

  // build dirichlet map

  Teuchos::RCP<const Epetra_Map> dirichletmap;

  STK::Mesh & mesh = dis_.GetMesh();
  stk::mesh::MetaData & meta = mesh.MetaData();
  stk::mesh::BulkData & bulk = mesh.BulkData();

  const Epetra_Map & dofrowmap = DofRowMap();
  Epetra_IntVector flag( dofrowmap );

  // Loop dirichlet conditions from higher to lower dimension since the
  // lower dimensions overwrite higher ones.

  static DRT::Condition::ConditionType dirichlet_types[] = {
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
    const std::map<int, Teuchos::ParameterList> * dirichlet = dis_.Condition( *condtype );
    if ( dirichlet!=NULL )
    {
      for ( std::map<int, Teuchos::ParameterList>::const_iterator i=dirichlet->begin();
            i!=dirichlet->end();
            ++i )
      {
        //int id = i->first;
        const Teuchos::ParameterList & condflags = i->second;

        std::string name = condflags.get<std::string>( "Name" );
        const std::vector<int> & onoff = * condflags.get<Teuchos::RCP<std::vector<int> > >( "onoff" );

        stk::mesh::Part * condpart = meta.get_part( name );
        if ( condpart==NULL )
          dserror( "No condition part '%s'", name.c_str() );

        stk::mesh::Selector selector = *condpart & mesh.OwnedPart();

        const std::vector<stk::mesh::Bucket*> & nodes = bulk.buckets( stk::mesh::Node );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator j=nodes.begin();
              j!=nodes.end();
              ++j )
        {
          stk::mesh::Bucket & bucket = **j;
          if ( selector( bucket ) )
          {
            for ( stk::mesh::Bucket::iterator k=bucket.begin(); k!=bucket.end(); ++k )
            {
              stk::mesh::Entity & n = *k;
              std::vector<int> dofids;
              Dof( n.key(), dofids );
              if ( dofids.size() > onoff.size() )
                dserror( "too few Dirichlet flags" );
#if 0
              if ( condflags.get<DRT::Condition::ConditionType>( "Type" )==DRT::Condition::PointDirichlet )
              {
                std::copy( dofids.begin(), dofids.end(), std::ostream_iterator<int>( std::cout, " " ) );
                std::cout << " - ";
                std::copy( onoff.begin(), onoff.end(), std::ostream_iterator<int>( std::cout, " " ) );
                std::cout << "\n";
              }
#endif
              for ( unsigned k=0; k<dofids.size(); ++k )
              {
                int lid = dofrowmap.LID( dofids[k] );
                if ( lid < 0 )
                  dserror( "illegal lid" );
                flag[lid] = onoff[k];
              }
            }
          }
        }
      }
    }
  }

  std::vector<int> dbc;
  int num = dofrowmap.NumMyElements();
  for ( int i=0; i<num; ++i )
  {
    if ( flag[i]!=0 )
    {
      dbc.push_back( dofrowmap.GID( i ) );
    }
  }

  dirichletmap = Teuchos::rcp( new Epetra_Map( -1, dbc.size(), &dbc[0], 0, dis_.Comm() ) );

  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back( Teuchos::null ); // non-dirichlet map is not needed?!
  maps.push_back( dirichletmap );
  dirichlet_.Setup( DofRowMap(), maps );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::DiscretizationState::EvaluateDirichlet( double time,
                                                       const std::vector<stk::mesh::FieldBase*> * v,
                                                       const std::vector<stk::mesh::FieldBase*> * dv,
                                                       const std::vector<stk::mesh::FieldBase*> * ddv )
{
  STK::Mesh & mesh = dis_.GetMesh();
  stk::mesh::MetaData & meta = mesh.MetaData();
  stk::mesh::BulkData & bulk = mesh.BulkData();

  //const Epetra_Map & dofrowmap = DofRowMap();

  unsigned deg = 0;  // highest degree of requested time derivative
  if ( dv!=NULL )
  {
    if ( v->size()!=dv->size() )
      dserror( "existing fields must match" );
    if ( ddv!=NULL )
    {
      if ( v->size()!=ddv->size() )
        dserror( "existing fields must match" );
      deg = 2;
    }
    else
    {
      deg = 1;
    }
  }

  static DRT::Condition::ConditionType dirichlet_types[] = {
    DRT::Condition::PointDirichlet,
    DRT::Condition::LineDirichlet,
    DRT::Condition::SurfaceDirichlet,
    DRT::Condition::VolumeDirichlet,
    DRT::Condition::none
  };

  stk::mesh::PartVector done;

  for ( DRT::Condition::ConditionType * condtype = dirichlet_types;
        *condtype!=DRT::Condition::none;
        ++condtype )
  {
    const std::map<int, Teuchos::ParameterList> * dirichlet = dis_.Condition( *condtype );
    if ( dirichlet!=NULL )
    {
      for ( std::map<int, Teuchos::ParameterList>::const_iterator i=dirichlet->begin();
            i!=dirichlet->end();
            ++i )
      {
        //int id = i->first;
        const Teuchos::ParameterList & condflags = i->second;

        if ( condflags.get<DRT::Condition::ConditionType>( "Type" )==*condtype )
        {
          // find the part, find which node is covered by which part

          std::string name = condflags.get<std::string>( "Name" );
          const std::vector<int> & onoff = * condflags.get<Teuchos::RCP<std::vector<int> > >( "onoff" );
          const std::vector<int> & curve = * condflags.get<Teuchos::RCP<std::vector<int> > >( "curve" );
          const std::vector<int> & funct = * condflags.get<Teuchos::RCP<std::vector<int> > >( "funct" );
          const std::vector<double> & val = * condflags.get<Teuchos::RCP<std::vector<double> > >( "val" );

          stk::mesh::Part * condpart = meta.get_part( name );
          if ( condpart==NULL )
            dserror( "No condition part '%s'", name.c_str() );

          stk::mesh::Selector selector = *condpart & mesh.OwnedPart();

          const std::vector<stk::mesh::Bucket*> & nodes = bulk.buckets( stk::mesh::Node );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator j=nodes.begin();
                j!=nodes.end();
                ++j )
          {
            stk::mesh::Bucket & bucket = **j;
            if ( selector( bucket ) )
            {
              bool isdone = false;
              for ( stk::mesh::PartVector::iterator k=done.begin(); k!=done.end(); ++k )
              {
                if ( has_superset( bucket, **k ) )
                {
                  isdone = true;
                  break;
                }
              }

              if ( not isdone )
              {
                for ( stk::mesh::Bucket::iterator k=bucket.begin(); k!=bucket.end(); ++k )
                {
                  stk::mesh::Entity & n = *k;

                  unsigned ifield = 0;
                  int field_base = 0;

                  int vfield_data_size = stk::mesh::field_data_size( *( *v )[ifield], bucket ) / sizeof( double );
                  double * vdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );

                  // apply Dirichlet condition to node dofs

                  std::vector<int> dofs;
                  dis_.Dof( n.key(), dofs );
                  const int numdf = dofs.size();
                  for ( int j=0; j<numdf; ++j )
                  {
                    if ( j >= field_base+vfield_data_size )
                    {
                      field_base += vfield_data_size;
                      ifield += 1;
                      if ( ifield >= v->size() )
                        dserror( "to few fields for Dirichlet evaluation" );
                      vfield_data_size = stk::mesh::field_data_size( *( *v )[ifield], bucket ) / sizeof( double );
                      vdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );
                    }

                    if ( onoff[j]!=0 )
                    {
                      //const int gid = dofs[j];
                      std::vector<double> value( deg+1, val[j] );

                      // factor given by time curve
                      std::vector<double> curvefac( deg+1, 1.0 );
                      int curvenum = curve[j];
                      if ( curvenum>=0 and time>=0.0 )
                      {
                        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer( time, deg );
                      }
                      else
                      {
                        for (unsigned i=1; i<(deg+1); ++i)
                          curvefac[i] = 0.0;
                      }

                      // factor given by spatial function
                      double functfac = 1.0;
                      int funct_num = funct[j];
                      if (funct_num>0)
                      {
                        stk::mesh::VectorField & coords = dis_.GetMesh().Coordinates();
                        double * x = stk::mesh::field_data( coords , n );
                        functfac = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(j,x,time,NULL);
                      }

                      // apply factors to Dirichlet value
                      for (unsigned i=0; i<deg+1; ++i)
                      {
                        value[i] *= functfac * curvefac[i];
                      }

                      vdata[j - field_base] = value[0];
                      if ( deg > 0 )
                      {
                        double * dvdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );
                        dvdata[j - field_base] = value[1];
                        if ( deg > 1 )
                        {
                          double * ddvdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );
                          ddvdata[j - field_base] = value[2];
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          done.push_back( condpart );
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> STK::FEI::DiscretizationState::NodeGraph()
{
  STK::Mesh & mesh = dis_.GetMesh();
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
      throw "illegal lid";
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
      throw "illegal lid";
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
void STK::FEI::DiscretizationState::CreateNodeMaps()
{
  {
    std::vector<int> vnids;
    stk::mesh::Selector selector = dis_.GetMesh().OwnedPart() & dis_.GetMesh().ActivePart();
    CreateNodeMap( selector, vnids );
    rownodemap_ = Teuchos::rcp( new Epetra_Map( -1,vnids.size(),&vnids[0],0,dis_.Comm() ) );
  }

  {
    std::vector<int> vnids;
    stk::mesh::Selector selector = dis_.GetMesh().ActivePart();
    CreateNodeMap( selector, vnids );
    colnodemap_ = Teuchos::rcp( new Epetra_Map( -1,vnids.size(),&vnids[0],0,dis_.Comm() ) );
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::DiscretizationState::CreateNodeMap( const stk::mesh::Selector & selector, std::vector<int> & vnids )
{
  // Use a set to gather nodal ids. This changes the order.
  std::set<int> nids;

  const std::vector<stk::mesh::Bucket*> & buckets = dis_.GetMesh().BulkData().buckets( stk::mesh::Node );

  // loop buckets and entries
  for (std::vector<stk::mesh::Bucket*>::const_iterator i=buckets.begin();
       i!=buckets.end();
       ++i)
  {
    stk::mesh::Bucket & bucket = **i;
    if ( selector( bucket ) )
    {
      for (stk::mesh::Bucket::iterator ib=bucket.begin();
           ib!=bucket.end();
           ++ib)
      {
        stk::mesh::Entity & e = *ib;
        int id = static_cast<int>( e.identifier() );
        nids.insert(id);
      }
    }
  }

  vnids.reserve(nids.size());
  vnids.assign(nids.begin(),nids.end());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::FEI::DiscretizationState::LocationVector( stk::mesh::Entity & e,
                                                    std::vector<int> & lm,
                                                    std::vector<int> & lmowner )
{
  //stk::mesh::Bucket & bucket = e.bucket();
  //std::vector<stk::mesh::FieldBase*> fields;
  //algo.collect_unknowns( fields );

  lm.clear();
  lmowner.clear();

  for ( stk::mesh::PairIterRelation nodes = e.relations( stk::mesh::Node );
        not nodes.empty();
        ++nodes )
  {
    stk::mesh::Entity & n = * nodes->entity();

    std::vector<int> dof;
    Dof( n.key(), dof );
    std::copy( dof.begin(), dof.end(), std::back_inserter( lm ) );
    int owner = n.owner_rank();
    lmowner.resize( lm.size(), owner );
  }
}

#endif

