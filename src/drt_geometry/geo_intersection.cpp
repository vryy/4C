
#ifdef QHULL

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "../drt_fluid/xfluid_defines.H"

#include "geo_intersection.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io_control.H"

// #include "searchtree.H"
// #include "searchtree_geometry_service.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_facet.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_mesh.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_side.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_parallel.H"


// int GEO::CutDofSet::NumDofPerNode( const DRT::Node & node, unsigned dspos ) const
// {
//   GEO::CUT::Node * n = mesh_->GetNode( node.Id() );
//   if ( n!=NULL )
//   {
//     int numdofpernode = DRT::DofSet::NumDofPerNode( node, dspos );
//     return numdofpernode * n->NumDofSets( include_inner_ );
//   }
//   return DRT::DofSet::NumDofPerNode( node, dspos );
// }

// int GEO::CutDofSet::PlainNumDofPerNode( const DRT::Node & node ) const
// {
//   return DRT::DofSet::NumDofPerNode( node, dspos_ );
// }

// void GEO::CutDofSet::Dof( DRT::Node & node, int nodaldofset, std::vector<int> & dofs ) const
// {
//   const int lid = node.LID();
//   if (lid==-1)
//     return;
//   int numdf = PlainNumDofPerNode( node );
//   const int idx = (*idxcolnodes_)[lid] + nodaldofset*numdf;
//   dofs.reserve( numdf );
//   for ( int i=0; i<numdf; ++i )
//   {
//     dofs.push_back( idx+i );
//   }
// }

GEO::CutWizard::CutWizard( DRT::Discretization & dis, bool include_inner, int numcutmesh )
  : dis_( dis ),
    include_inner_( include_inner )
{
  mesh_ = Teuchos::rcp( new GEO::CUT::MeshIntersection( numcutmesh ) );
}

void GEO::CutWizard::SetFindPositions( bool positions )
{
  mesh_->SetFindPositions( positions );
}

void GEO::CutWizard::AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze )
{
  const int numnode = ele->NumNode();
  const int * nodeids = ele->NodeIds();

  std::vector<int> nids( nodeids, nodeids+numnode );
  mesh_->AddCutSide( ele->Id(), nids, xyze, ele->Shape(), mi );
}

void GEO::CutWizard::AddElement( DRT::Element * ele )
{
  const int numnode = ele->NumNode();
  const DRT::Node * const * nodes = ele->Nodes();
  const int * nodeids = ele->NodeIds();

  Epetra_SerialDenseMatrix xyze( 3, numnode );

  for ( int i=0; i < numnode; ++i )
  {
    const DRT::Node & node = *nodes[i];
    std::copy( node.X(), node.X()+3, &xyze( 0, i ) );
  }

  std::vector<int> nids( nodeids, nodeids+numnode );
  mesh_->AddElement( ele->Id(), nids, xyze, ele->Shape() );
}

GEO::CUT::ElementHandle * GEO::CutWizard::GetElement( DRT::Element * ele )
{
  return mesh_->GetElement( ele->Id() );
}

GEO::CUT::Node * GEO::CutWizard::GetNode( int nid )
{
  return mesh_->GetNode( nid );
}




/*------------------------------------------------------------------------------------------------*
 * cut routine for parallel framework in XFSI and XFLUIDFLUID                        schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::CutParallel( bool include_inner, std::string VCellgausstype, std::string BCellgausstype )
{
  // for XFSI and XFLUIDFLUID we have communicate node positions and dofset data
  bool parallel = true;

  mesh_->Status();


  // FIRST step: cut the mesh
  mesh_->Cut_Mesh( include_inner );

  // SECOND step: find node positions and create dofset in PARALLEL
  CutParallel_FindPositionDofSets( include_inner, parallel );

  // THIRD step: perform tesselation or moment fitting on the mesh
  mesh_->Cut_Finalize( include_inner, VCellgausstype, BCellgausstype );


  mesh_->Status(VCellgausstype);
}


/*------------------------------------------------------------------------------------------------*
 * routine for finding node positions and computing vc dofsets in a parallel way     schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::CutParallel_FindPositionDofSets(bool include_inner, bool parallel)
{
  GEO::CUT::Options options;
  mesh_->GetOptions(options);

  if ( options.FindPositions() )
  {

    GEO::CUT::Mesh & m = mesh_->NormalMesh();

    // find inside and outside positions of nodes
    m.FindNodePositions();

    // find undecided nodes
    // * for serial simulations all node positions should be set
    // * for parallel simulations there can be some undecided nodes

    // create a parallel Cut object for the current background mesh to communicate missing data
    Teuchos::RCP<GEO::CUT::Parallel> cut_parallel = Teuchos::rcp( new GEO::CUT::Parallel( dis_, m, *mesh_ ) );

    if(parallel) cut_parallel->CommunicateNodePositions();

    m.FindFacetPositions();

    // find number and connection of dofsets at nodes from cut volumes
    #ifdef DOFSETS_NEW
      mesh_->CreateNodalDofSetNEW( include_inner, dis_);
    #else
      m.FindNodalDOFSets( include_inner );
    #endif

    if(parallel) cut_parallel->CommunicateNodeDofSetNumbers();

  }

}


/*------------------------------------------------------------------------------------------------*
 *cut routine for standard non-parallel framework (only for cuttest)                 schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::Cut( bool include_inner, std::string VCellgausstype, std::string BCellgausstype )
{
  mesh_->Cut( include_inner, VCellgausstype, BCellgausstype );
}


/*------------------------------------------------------------------------------------------------*
 * cut routine for standard combustion and two phase flow framework                  schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void GEO::CutWizard::Cut( std::map< int, DomainIntCells >& domainintcells,
                          std::map< int, BoundaryIntCells >& boundaryintcells,
                          std::string VCellgausstype,
                          std::string BCellgausstype)
{
  mesh_->Cut( true, VCellgausstype, BCellgausstype );

  std::vector<int> cutelements;

  for ( int k = 0; k < dis_.NumMyColElements(); ++k )
  {
    DRT::Element* ele = dis_.lColElement( k );
    int eid = ele->Id();
    GEO::CUT::ElementHandle * e = mesh_->GetElement( eid );

    if ( e!=NULL and e->IsCut() )
    {
      cutelements.push_back( eid );

      GEO::CUT::plain_integrationcell_set cells;
      e->GetIntegrationCells( cells );

      GEO::CUT::plain_boundarycell_set bcells;
      e->GetBoundaryCells( bcells );

      BoundaryIntCells & bics = boundaryintcells[eid];
      DomainIntCells   & dics = domainintcells  [eid];

      LINALG::Matrix<3,1> physCoordCorner;
      LINALG::Matrix<3,1> eleCoordDomainCorner;
      LINALG::Matrix<3,1> eleCoordBoundaryCorner;

      for ( GEO::CUT::plain_boundarycell_set::iterator i=bcells.begin();
            i!=bcells.end();
            ++i )
      {
        GEO::CUT::BoundaryCell * bc = *i;
        DRT::Element::DiscretizationType distype = bc->Shape();
        int numnodes = DRT::UTILS::getNumberOfElementNodes( distype );

        GEO::CUT::Facet * f = bc->GetFacet();
        GEO::CUT::Side * s = f->ParentSide();

        LINALG::SerialDenseMatrix physDomainCoord = bc->Coordinates();

        LINALG::SerialDenseMatrix eleDomainCoord(3, numnodes); // in xfem parent element domain
        LINALG::SerialDenseMatrix eleBoundaryCoord(3, numnodes);

        for ( int j=0; j<numnodes; ++j )
        {
          std::copy( &physDomainCoord( 0, j ), &physDomainCoord( 0, j ) + 3, physCoordCorner.A() );

          e->LocalCoordinates( physCoordCorner, eleCoordDomainCorner );
          s->LocalCoordinates( physCoordCorner, eleCoordBoundaryCorner );

          std::copy( eleCoordDomainCorner  .A(), eleCoordDomainCorner  .A()+3, &eleDomainCoord  ( 0, j ) );
          std::copy( eleCoordBoundaryCorner.A(), eleCoordBoundaryCorner.A()+3, &eleBoundaryCoord( 0, j ) );
        }
        bics.push_back( BoundaryIntCell( distype, s->Id(), eleDomainCoord, eleBoundaryCoord, physDomainCoord ) );
      }

      for ( GEO::CUT::plain_integrationcell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::IntegrationCell * ic = *i;
        DRT::Element::DiscretizationType distype = ic->Shape();
        int numnodes = DRT::UTILS::getNumberOfElementNodes( distype );

        LINALG::SerialDenseMatrix physCoord = ic->Coordinates();
        LINALG::SerialDenseMatrix coord( 3, numnodes );

        for ( int j=0; j<numnodes; ++j )
        {
          std::copy( &physCoord( 0, j ), &physCoord( 0, j ) + 3, physCoordCorner.A() );

          e->LocalCoordinates( physCoordCorner, eleCoordDomainCorner );

          std::copy( eleCoordDomainCorner.A(), eleCoordDomainCorner.A()+3, &coord( 0, j ) );
        }

        if ( distype==DRT::Element::tet4 )
        {
          // create planes consisting of 3 nodes each
          LINALG::Matrix<3,1> p0( coord.A()  , true );
          LINALG::Matrix<3,1> p1( coord.A()+3, true );
          LINALG::Matrix<3,1> p2( coord.A()+6, true );
          LINALG::Matrix<3,1> p3( coord.A()+9, true );

          LINALG::Matrix<3,1> v01;
          LINALG::Matrix<3,1> v02;
          LINALG::Matrix<3,1> v03;

          v01.Update( 1, p1, -1, p0, 0 );
          v02.Update( 1, p2, -1, p0, 0 );
          v03.Update( 1, p3, -1, p0, 0 );

          // create 4 normal vectors to each tet surface plane
          LINALG::Matrix<3,1> nplane012;

          // cross product
          nplane012(0) = v01(1)*v02(2) - v01(2)*v02(1);
          nplane012(1) = v01(2)*v02(0) - v01(0)*v02(2);
          nplane012(2) = v01(0)*v02(1) - v01(1)*v02(0);

          // compute norm (area) of plane
          double norm012 = nplane012.Norm2();

          // compute normal distance of point to plane of the three remaining points
          double distance = nplane012.Dot( v03 );

          double vol_tet = distance / 6.0;

          // smallest volume
          if ( fabs( vol_tet ) < 1e-10 )
            continue;

          if ( fabs( distance / norm012 ) < 1e-7 )
            continue;

          // tet numbering wrong exchange 1 with 3
          if ( distance < 0 )
          {
            for ( int i = 0; i < 3; ++i )
            {
              std::swap( coord    ( i, 1 ), coord    ( i, 3 ) );
              std::swap( physCoord( i, 1 ), physCoord( i, 3 ) );
            }
          }
        }

        dics.push_back( DomainIntCell( distype, coord, physCoord, ic->Position()==GEO::CUT::Point::outside ) );
      }
    }
  }

  cutelementmap_ = Teuchos::rcp( new Epetra_Map( -1, cutelements.size(), &cutelements[0], 0, dis_.Comm() ) );
}

void GEO::CutWizard::PrintCellStats()
{
  mesh_->PrintCellStats();
}

void GEO::CutWizard::DumpGmshNumDOFSets( bool include_inner)
{
  std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << filename;

  mesh_->DumpGmshNumDOFSets( str.str(), include_inner, dis_ );
}

void GEO::CutWizard::DumpGmshVolumeCells( bool include_inner )
{
  std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << name
      << ".volumecells."
      << dis_.Comm().MyPID()
      << ".pos";
  mesh_->DumpGmshVolumeCells( str.str(), include_inner );
}

void GEO::CutWizard::DumpGmshIntegrationCells()
{
  std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
  std::stringstream str;
  str << name
      << ".integrationcells."
      << dis_.Comm().MyPID()
      << ".pos";
  mesh_->DumpGmshIntegrationCells( str.str() );
}

#if 0
Teuchos::RCP<Epetra_CrsGraph> GEO::CutWizard::MatrixGraph( const CutDofSet & dofset, const Epetra_Map & dbcmap )
{
  const Epetra_Map & noderowmap = * dis_.NodeRowMap();
  const Epetra_Map & dofrowmap  = * dofset.DofRowMap();

  int numnode = noderowmap.NumMyElements();

  // build graph

  std::map<int, std::set<int> > graph;

  int row = 0;
  for ( int i=0; i<numnode; ++i )
  {
    DRT::Node * rn = dis_.lRowNode( i );

    int numelements = rn->NumElement();
    DRT::Element** elements = rn->Elements();
    if ( elements==NULL )
      dserror( "no elements at node %d", rn->Id() );

    GEO::CUT::Node * crn = mesh_->GetNode( rn->Id() );
    if ( crn!=NULL )
    {
      std::vector<int> rdofs = dis_.Dof( rn );
      int rplain_numdf = dofset.PlainNumDofPerNode( *rn );
      for ( int k=0; k<numelements; ++k )
      {
        DRT::Element * e = elements[k];
        int numnodes = e->NumNode();
        DRT::Node ** nodes = e->Nodes();

        GEO::CUT::ElementHandle * ce = mesh_->GetElement( e->Id() );
        if ( ce!=NULL )
        {
          std::set<GEO::CUT::VolumeCell*> cells;
          ce->GetVolumeCells( cells );
          for ( std::set<GEO::CUT::VolumeCell*>::iterator i=cells.begin();
                i!=cells.end();
                ++i )
          {
            GEO::CUT::VolumeCell * c = *i;
            const std::vector<int> & nodaldofset = c->NodalDofSet();
            if ( nodaldofset.size()!=numnodes )
            {
              throw std::runtime_error( "number of nodes mismatch" );
            }
            for ( int l=0; l<numnodes; ++l )
            {
              DRT::Node * cn = nodes[l];
              std::vector<int> cdofs;
              dofset.Dof( *cn, nodaldofset[l], cdofs );
            }
          }
        }
        else
        {
          for ( int l=0; l<numnodes; ++l )
          {
            DRT::Node * cn = nodes[l];

            std::vector<int> cdofs = dis_.Dof( cn );
          }
        }
      }
    }
    else
    {
      std::vector<int> rdofs = dis_.Dof( rn );
      for ( unsigned rj=0; rj<rdofs.size(); ++rj, ++row )
      {
        std::set<int> & rowset = graph[rdofs[rj]];
        if ( not dbcmap.MyGID( rdofs[rj] ) )
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

              std::vector<int> cdofs = dis_.Dof( cn );
              rowset.insert( cdofs.begin(), cdofs.end() );
            }
          }
        }
        else
        {
          // just diagonal entry on Dirichlet rows
          rowset.insert( rdofs[rj] );
        }
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
  return Teuchos::null;
}
#endif


void GEO::computeIntersection( const Teuchos::RCP<DRT::Discretization> xfemdis,
                               const Teuchos::RCP<DRT::Discretization> cutterdis,
                               const std::map<int,LINALG::Matrix<3,1> >& currentcutterpositions,
                               const std::map<int,LINALG::Matrix<3,2> >& currentXAABBs,
                               std::map< int, DomainIntCells >& domainintcells,
                               std::map< int, BoundaryIntCells >& boundaryintcells,
                               const std::map<int,int>& labelPerElementId,
                               const std::vector<int>& MovingFluideleGIDs,
                               std::string VCellgausstype,
                               std::string BCellgausstype)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::computeIntersection" );

  if ( xfemdis->Comm().MyPID() == 0 )
    std::cout << "\nGEO::Intersection:" << std::flush;

  const double t_start = Teuchos::Time::wallTime();

  CutWizard wizard( *xfemdis, false );

  for ( int k = 0; k < cutterdis->NumMyColElements(); ++k )
  {
    DRT::Element* ele = cutterdis->lColElement( k );

    std::map<int,int>::const_iterator k = labelPerElementId.find( ele->Id() );
    if ( k==labelPerElementId.end() )
    {
      dserror( "no label for cutter element %d", ele->Id() );
    }
    if ( k->second < 1 )
      continue;

    const int numnode = ele->NumNode();
    const DRT::Node * const * nodes = ele->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      const DRT::Node & node = *nodes[i];
      std::map<int,LINALG::Matrix<3,1> >::const_iterator j = currentcutterpositions.find( node.Id() );
      if ( j==currentcutterpositions.end() )
      {
        dserror( "no positions to node %d", node.Id() );
      }
      const LINALG::Matrix<3,1> & x = j->second;
      std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
    }

    wizard.AddCutSide( 0, ele, xyze );
  }

  for ( int k = 0; k < xfemdis->NumMyColElements(); ++k )
  {
    DRT::Element* xfemElement = xfemdis->lColElement( k );

    // for fluid-fluid-coupling consider just the elements of background fluid
    if ( cutterdis->Name() == "FluidFluidboundary" or
         cutterdis->Name() == "ALEFluidboundary" )
    {
      if ( std::find( MovingFluideleGIDs.begin(),
                      MovingFluideleGIDs.end(),
                      xfemElement->Id() ) != MovingFluideleGIDs.end() )
      {
        continue;
      }
    }

    wizard.AddElement( xfemElement );
  }

  // Call tetgen on all cut elements.

  wizard.Cut( domainintcells, boundaryintcells, VCellgausstype, BCellgausstype );

  // cleanup

  int localcells = domainintcells.size();
  int globalcells;
  xfemdis->Comm().SumAll( &localcells, &globalcells, 1 );

  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( xfemdis->Comm().MyPID() == 0 )
  {
    std::cout << " Success (" << t_end  <<  " secs), intersected elements: " << globalcells
              << "\n";
  }

  wizard.PrintCellStats();
  wizard.DumpGmshIntegrationCells();
}

#endif
