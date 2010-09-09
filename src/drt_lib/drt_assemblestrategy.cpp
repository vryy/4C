
#include "drt_assemblestrategy.H"
#include "drt_discret.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_utils.H"

DRT::AssembleStrategy::AssembleStrategy( Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                                         Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                                         Teuchos::RCP<Epetra_Vector> systemvector1,
                                         Teuchos::RCP<Epetra_Vector> systemvector2,
                                         Teuchos::RCP<Epetra_Vector> systemvector3 )
  : systemmatrix1_( systemmatrix1 ),
    systemmatrix2_( systemmatrix2 ),
    systemvector1_( systemvector1 ),
    systemvector2_( systemvector2 ),
    systemvector3_( systemvector3 )
{
}

DRT::AssembleStrategy::~AssembleStrategy()
{
}

Teuchos::RCP<Epetra_CrsGraph> DRT::AssembleStrategy::MatrixGraph( DRT::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap )
{
  const Epetra_Map & noderowmap = * dis.NodeRowMap();
  const Epetra_Map & dofrowmap  = * dis.DofRowMap();

  int numnode = noderowmap.NumMyElements();

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
        for ( int k=0; k<numelements; ++k )
        {
          DRT::Element * e = elements[k];
          int numnodes = e->NumNode();
          DRT::Node ** nodes = e->Nodes();
          for ( int l=0; l<numnodes; ++l )
          {
            DRT::Node * cn = nodes[l];

            std::vector<int> cdofs = dis.Dof( cn );
            rowset.insert( cdofs.begin(), cdofs.end() );
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

  return crsgraph;
}

void DRT::AssembleStrategy::Zero()
{
  if ( Assemblemat1() )
  {
    systemmatrix1_->Zero();
  }
  if ( Assemblemat2() )
  {
    systemmatrix1_->Zero();
  }
  if ( Assemblevec1() )
  {
    systemvector1_->PutScalar( 0.0 );
  }
  if ( Assemblevec2() )
  {
    systemvector2_->PutScalar( 0.0 );
  }
  if ( Assemblevec3() )
  {
    systemvector3_->PutScalar( 0.0 );
  }
}

void DRT::AssembleStrategy::Complete()
{
  if ( Assemblemat1() )
  {
    systemmatrix1_->Complete();
  }
  if ( Assemblemat2() )
  {
    systemmatrix1_->Complete();
  }
}

void DRT::AssembleStrategy::ClearElementStorage( int eledim )
{
  if (Assemblemat1())
  {
    if (elematrix1_.M()!=eledim or elematrix1_.N()!=eledim)
      elematrix1_.Shape(eledim,eledim);
    else
      memset(elematrix1_.A(),0,eledim*eledim*sizeof(double));
  }
  if (Assemblemat2())
  {
    if (elematrix2_.M()!=eledim or elematrix2_.N()!=eledim)
      elematrix2_.Shape(eledim,eledim);
    else
      memset(elematrix2_.A(),0,eledim*eledim*sizeof(double));
  }
  if (Assemblevec1())
  {
    if (elevector1_.Length()!=eledim)
      elevector1_.Size(eledim);
    else
      memset(elevector1_.Values(),0,eledim*sizeof(double));
  }
  if (Assemblevec2())
  {
    if (elevector2_.Length()!=eledim)
      elevector2_.Size(eledim);
    else
      memset(elevector2_.Values(),0,eledim*sizeof(double));
  }
  if (Assemblevec3())
  {
    if (elevector3_.Length()!=eledim)
      elevector3_.Size(eledim);
    else
      memset(elevector3_.Values(),0,eledim*sizeof(double));
  }
}

void DRT::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat,
                                     int eid,
                                     const Epetra_SerialDenseMatrix& Aele,
                                     const std::vector<int>& lm,
                                     const std::vector<int>& lmowner)
{
  sysmat.Assemble( eid, Aele, lm, lmowner );
}

void DRT::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat,
                                     int eid,
                                     const Epetra_SerialDenseMatrix& Aele,
                                     const std::vector<int>& lmrow,
                                     const std::vector<int>& lmrowowner,
                                     const std::vector<int>& lmcol)
{
  sysmat.Assemble( eid, Aele, lmrow, lmrowowner, lmcol );
}

void DRT::AssembleStrategy::Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid)
{
  sysmat.Assemble( val, rgid, cgid );
}

void DRT::AssembleStrategy::Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                                     const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  LINALG::Assemble( V, Vele, lm, lmowner );
}

void DRT::AssembleStrategy::Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                                     const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  LINALG::Assemble( V, n, Vele, lm, lmowner );
}
