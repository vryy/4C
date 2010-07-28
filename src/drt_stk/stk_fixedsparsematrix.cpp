#ifdef STKADAPTIVE

#include "stk_fixedsparsematrix.H"
#include "../stk_lib/stk_mesh.H"

#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FixedSparseMatrix::FixedSparseMatrix( Teuchos::RCP<const Epetra_Map> dbcmap )
{
  dbcmap_ = dbcmap;

#if 0
  const Epetra_Map & noderowmap = * dis.NodeRowMap();
  const Epetra_Map & dofrowmap  = * dis.DofRowMap();

  int numnode = noderowmap.NumMyElements();

  // hanging node dof maps
  // extract dof information to be able to assemble later on

  stk::mesh::BulkData & bulk_data = mesh.BulkData();
  //stk::mesh::Part & owned = mesh.OwnedPart();

  const std::vector<stk::mesh::Bucket*> & constraints = bulk_data.buckets( stk::mesh::Constraint );
  for ( std::vector<stk::mesh::Bucket*>::const_iterator i=constraints.begin();
        i!=constraints.end();
        ++i )
  {
    stk::mesh::Bucket & bucket = **i;

    // if I own the constraint I own the hanging node and vice versa
    //if ( has_superset( bucket, owned ) )

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

          DRT::Node * hnode = dis.gNode( hn.key().id()-1 );
          if ( hnode==NULL )
            dserror( "discretization mismatch" );

          std::vector<int> hndofs = dis.Dof( hnode );
          for ( std::vector<int>::iterator k=hndofs.begin();
                k!=hndofs.end(); ++k )
          {
            int d = *k;
            DofData & dd = hndofmap_[d];
            dd.fact = fact;
            std::set<int> & realdofmap = dd.realdofmap;
            for ( stk::mesh::PairIterRelation::iterator r=rel.begin(); r!=rel.end(); ++r )
            {
              stk::mesh::Entity & n = *r->entity();

              DRT::Node * node = dis.gNode( n.key().id()-1 );
              if ( node==NULL )
                dserror( "discretization mismatch" );

              std::vector<int> dofs = dis.Dof( node );
              realdofmap.insert( dofs.begin(), dofs.end() );
            }
          }
        }
      }
    }
  }

  // build graph

  std::map<int, std::set<int> > graph;

  int row = 0;
  for ( int i=0; i<numnode; ++i )
  {
    DRT::Node * rn = dis.lRowNode( i );

    int numelements = rn->NumElement();
    DRT::Element** elements = rn->Elements();
    if ( elements==NULL )
      dserror( "no elements at node %d", rn->Id() );

    std::vector<int> rdofs = dis.Dof( rn );
    for ( unsigned rj=0; rj<rdofs.size(); ++rj, ++row )
    {
      std::set<int> & rowset = graph[rdofs[rj]];
      if ( not dbcmap->MyGID( rdofs[rj] ) )
      {
        // non-Dirichlet row

        std::map<int, DofData>::iterator rhndofmap = hndofmap_.find( rdofs[rj] );
        if ( rhndofmap==hndofmap_.end() )
        {
          for ( int k=0; k<numelements; ++k )
          {
            DRT::Element * e = elements[k];
            int numnodes = e->NumNode();
            DRT::Node ** nodes = e->Nodes();
            for ( int l=0; l<numnodes; ++l )
            {
              DRT::Node * cn = nodes[l];

              std::vector<int> cdofs = dis.Dof( cn );
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

          // Add my line to all those connected lines, if those lines are no
          // Dirichlet lines. Here we do not have to watch for hanging nodes
          // any more.
          for ( std::set<int>::iterator r=realdofmap.begin(); r!=realdofmap.end(); ++r )
          {
            int rd = *r;
            if ( not dbcmap->MyGID( rd ) and dofrowmap.MyGID( rd ) )
            {
              std::set<int> & hnrowset = graph[rd];

              for ( int k=0; k<numelements; ++k )
              {
                DRT::Element * e = elements[k];
                int numnodes = e->NumNode();
                DRT::Node ** nodes = e->Nodes();
                for ( int l=0; l<numnodes; ++l )
                {
                  DRT::Node * cn = nodes[l];

                  std::vector<int> cdofs = dis.Dof( cn );
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

  sysmat_ = Teuchos::rcp( new Epetra_CrsMatrix( Copy, *crsgraph ) );
#endif
}


#if 0
void STK::FixedSparseMatrix::Assemble(int eid,
                                      const Epetra_SerialDenseMatrix& Aele,
                                      const std::vector<int>& lmrow,
                                      const std::vector<int>& lmrowowner,
                                      const std::vector<int>& lmcol)
{
  if (not sysmat_->Filled())
    dserror( "not filled" );

  const unsigned lrowdim = lmrow.size();
  const unsigned lcoldim = lmcol.size();

  const int myrank = sysmat_->Comm().MyPID();
  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();

  // loop rows of local matrix
  for (unsigned lrow=0; lrow<lrowdim; ++lrow)
  {
    const int rgid = lmrow[lrow];

    // if we have a Dirichlet map check if this row is a Dirichlet row
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
            const int errone = sysmat_->SumIntoMyValues(rlid,1,&value,&idx);
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
              const int errone = sysmat_->SumIntoMyValues(rlid,1,&fvalue,&idx);
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

      const int errone = sysmat_->SumIntoMyValues(rlid,1,&value,&rlid);
      if (errone)
        dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);

      // put hanging node row to real node rows

      DofData & rdd = rhndofmap->second;
      std::set<int> & realdofmap = rdd.realdofmap;

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
              const int errone = sysmat_->SumIntoMyValues(rlid,1,&value,&idx);
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
                const int errone = sysmat_->SumIntoMyValues(rlid,1,&fvalue,&idx);
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
#endif
#endif
