
#ifdef STKADAPTIVE

#include "stk_linearsystem.H"
#include "stk_discret.H"
#include "../stk_refine/stk_mesh.H"

#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::LinearSystem::LinearSystem( STK::Discretization & dis )
  : rowdbcmap_( dis.DofRowDirichletMap() ),
    coldbcmap_( dis.DofColDirichletMap() ),
    importer_( dis.DofColMap(), dis.DofRowMap() )
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
              if ( dofs.size() <= k )
                dserror( "too few dofs in supporting node. expected %d but got %d", k+1, dofs.size() );
              realdofmap.insert( dofs[k] );
            }
          }
        }
        else
        {
          dserror( "no nodes in constraint" );
        }
      }
    }
  }

  sysmat_ = Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *MatrixGraph( dis ) ) );
  rowrhs_ = Teuchos::rcp( new Epetra_Vector( dis.DofRowMap() ) );
  rhs_    = Teuchos::rcp( new Epetra_Vector( dis.DofColMap() ) );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsGraph> STK::LinearSystem::MatrixGraph( STK::Discretization & dis )
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
          if ( not rowdbcmap_.MyGID( rdofs[rj] ) )
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
                if ( not rowdbcmap_.MyGID( rd ) and dofrowmap.MyGID( rd ) )
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
void STK::LinearSystem::Assemble(int eid,
                                 const Epetra_SerialDenseMatrix& Aele,
                                 const std::vector<int>& lmrow,
                                 const std::vector<int>& lmrowowner,
                                 const std::vector<int>& lmcol)
{
  Epetra_FECrsMatrix * mat = &*sysmat_;

  if (not mat->Filled())
    dserror( "not filled" );

#if 0
  if ( eid==382 or eid==438 or eid==383 )
  {
    std::cout << "P" << dbcmap_->Comm().MyPID() << ": assemble element " << eid << "\n";
  }
#endif

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
    if ( coldbcmap_.MyGID( rgid ) ) continue;

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
      else
      {
        // off-row assemble
        // here we use the global id versions of the code above
#if 1
        for (unsigned lcol=0; lcol<lcoldim; ++lcol)
        {
          double value = Aele(lrow,lcol);
          int cgid = lmcol[lcol];

          // check for hanging node row
          std::map<int, DofData>::iterator chndofmap = hndofmap_.find( cgid );

          if ( chndofmap==hndofmap_.end() )
          {
            const int errone = mat->SumIntoGlobalValues(rgid,1,&value,&cgid);
            if (errone)
              dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
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
              const int errone = mat->SumIntoGlobalValues(rgid,1,&fvalue,&gid);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
            }
          }
        }
#endif
      }
    }
    else
    {
      // put hanging node row to real node rows

      DofData & rdd = rhndofmap->second;
      std::set<int> & realdofmap = rdd.realdofmap;

      double fvalue = -rdd.fact;

      // check ownership of row
      if (lmrowowner[lrow] == myrank)
      {
        // do innocent decoupled value first

        int rlid = rowmap.LID(rgid);
        double value = 1;

        const int errone = mat->ReplaceMyValues(rlid,1,&value,&rlid);
        if (errone)
          dserror("Epetra_CrsMatrix::ReplaceMyValues returned error code %d",errone);

        for ( std::set<int>::iterator c=realdofmap.begin(); c!=realdofmap.end(); ++c )
        {
          int cgid = *c;
          int clid = colmap.LID(cgid);
          const int errone = mat->ReplaceMyValues(rlid,1,&fvalue,&clid);
          if (errone)
            dserror("Epetra_CrsMatrix::ReplaceMyValues(%d,1,%f,%d) returned error code %d",rlid,fvalue,clid,errone);
        }
      }
      else
      {
        // off-row assemble
#if 1
        double value = 1;

        const int errone = mat->ReplaceGlobalValues(rgid,1,&value,&rgid);
        if (errone)
          dserror("Epetra_CrsMatrix::ReplaceGlobalValues returned error code %d",errone);

        for ( std::set<int>::iterator c=realdofmap.begin(); c!=realdofmap.end(); ++c )
        {
          int cgid = *c;
          const int errone = mat->ReplaceGlobalValues(rgid,1,&fvalue,&cgid);
          if (errone)
            dserror("Epetra_CrsMatrix::ReplaceGlobalValues returned error code %d",errone);
        }
#endif
      }

      for ( std::set<int>::iterator r=realdofmap.begin(); r!=realdofmap.end(); ++r )
      {
        int rgid = *r;
        if ( not coldbcmap_.MyGID( rgid ) )
        {
          if ( rowmap.MyGID( rgid ) )
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
          else
          {
            // off-row assemble
#if 1
            for (unsigned lcol=0; lcol<lcoldim; ++lcol)
            {
              double value = Aele(lrow,lcol)*rdd.fact;
              int cgid = lmcol[lcol];

              // check for hanging node row
              std::map<int, DofData>::iterator chndofmap = hndofmap_.find( cgid );

              if ( chndofmap==hndofmap_.end() )
              {
                const int errone = mat->SumIntoGlobalValues(rgid,1,&value,&cgid);
                if (errone)
                  dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
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
                  const int errone = mat->SumIntoGlobalValues(rgid,1,&fvalue,&gid);
                  if (errone)
                    dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
                }
              }
            }
#endif
          }
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::LinearSystem::Assemble(int eid,
                                 const Epetra_SerialDenseVector& Vele,
                                 const std::vector<int>& lm)
{
  Epetra_Vector & V = *rhs_;

  const unsigned ldim = lm.size();
  if ( ldim!=static_cast<unsigned>( Vele.Length() ) )
    dserror("Mismatch in dimensions");

  const Epetra_BlockMap& colmap = V.Map();

  for (unsigned lrow=0; lrow<ldim; ++lrow)
  {
    int rgid = lm[lrow];

    int rlid = colmap.LID( rgid );
    if ( rlid < 0 )
      dserror("Sparse vector V does not have global row %d", rgid);

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
        int rlid = colmap.LID( rgid );
        if ( rlid < 0 )
          dserror("Sparse vector V does not have global row %d", rgid);
        V[rlid] += fvalue;
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::LinearSystem::Complete()
{
  int err = sysmat_->GlobalAssemble();
  if ( err )
    dserror( "Epetra_FECrsMatrix::GlobalAssemble failed: err=%d", err );

  err = rowrhs_->Export( *rhs_, importer_, Add );
  if (err)
    dserror("Export using importer returned err=%d", err);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::LinearSystem::Zero()
{
  sysmat_->PutScalar( 0.0 );
  rowrhs_->PutScalar( 0.0 );    // just to make sure nobody gets confused
  rhs_   ->PutScalar( 0.0 );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::LinearSystem::ApplyDirichlet( Epetra_Vector & dirichlet )
{
  Epetra_Vector & rhs = *rowrhs_;
  Epetra_FECrsMatrix & mat = *sysmat_;

  const Epetra_BlockMap & rowmap = rhs.Map();

  if ( not rowmap.PointSameAs( dirichlet.Map() ) )
    dserror( "expect same maps here" );

  int length = rowdbcmap_.NumMyElements();
  int * gids = rowdbcmap_.MyGlobalElements();

  double v = 1.0;

  for (int i=0; i<length; ++i)
  {
    int gid = gids[i];
    int lid = rowmap.LID( gid );
    if ( lid<0 )
      dserror( "illegal dirichlet map" );
    rhs[lid] = dirichlet[lid];

    int err = mat.ReplaceMyValues(lid,1,&v,&lid);
    if (err<0)
      dserror("Epetra_CrsMatrix::ReplaceGlobalValues returned err=%d",err);
  }
}

#endif

