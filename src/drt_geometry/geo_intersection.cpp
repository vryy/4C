
#ifdef QHULL

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "geo_intersection.H"
#include "../drt_lib/drt_discret.H"

// #include "searchtree.H"
// #include "searchtree_geometry_service.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_element.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_side.H"
#include "../drt_cut/cut_facet.H"

void GEO::computeIntersection( const Teuchos::RCP<DRT::Discretization> xfemdis,
                               const Teuchos::RCP<DRT::Discretization> cutterdis,
                               const std::map<int,LINALG::Matrix<3,1> >& currentcutterpositions,
                               const std::map<int,LINALG::Matrix<3,2> >& currentXAABBs,
                               std::map< int, DomainIntCells >& domainintcells,
                               std::map< int, BoundaryIntCells >& boundaryintcells,
                               const std::map<int,int>& labelPerElementId,
                               const std::vector<int>& MovingFluideleGIDs )
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::computeIntersection" );

  if ( xfemdis->Comm().MyPID() == 0 )
    std::cout << "\nGEO::Intersection:" << std::flush;

  const double t_start = Teuchos::Time::wallTime();

  GEO::CUT::MeshIntersection intersection;

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
      //std::copy( node.X(), node.X()+3, &xyze( 0, i ) );
      std::map<int,LINALG::Matrix<3,1> >::const_iterator j = currentcutterpositions.find( node.Id() );
      if ( j==currentcutterpositions.end() )
      {
        dserror( "no positions to node %d", node.Id() );
      }
      const LINALG::Matrix<3,1> & x = j->second;
      std::copy( x.A(), x.A()+3, &xyze( 0, i ) );
    }

    std::vector<int> nids( ele->NodeIds(), ele->NodeIds()+numnode );

    intersection.AddCutSide( ele->Id(), nids, xyze, ele->Shape() );
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

    const int numnode = xfemElement->NumNode();
    const DRT::Node * const * nodes = xfemElement->Nodes();

    Epetra_SerialDenseMatrix xyze( 3, numnode );

    for ( int i=0; i < numnode; ++i )
    {
      const DRT::Node & node = *nodes[i];
      std::copy( node.X(), node.X()+3, &xyze( 0, i ) );
    }

    std::vector<int> nids( xfemElement->NodeIds(), xfemElement->NodeIds()+numnode );

    intersection.AddElement( xfemElement->Id(), nids, xyze, xfemElement->Shape() );
  }

  // Call tetgen on all cut elements. The integration cells are gathered via
  // the callback generator object

  //CellGenerator generator( *xfemdis, intersection, domainintcells, boundaryintcells );
  intersection.Cut();

  for ( int k = 0; k < xfemdis->NumMyColElements(); ++k )
  {
    DRT::Element* xfemElement = xfemdis->lColElement( k );
    int eid = xfemElement->Id();
    GEO::CUT::Element * e = intersection.GetElement( eid );

    if ( e!=NULL and e->IsCut() )
    {
      std::set<GEO::CUT::IntegrationCell*> cells;
      e->GetIntegrationCells( cells );

      std::set<GEO::CUT::BoundaryCell*> bcells;
      e->GetBoundaryCells( bcells );

      BoundaryIntCells & bics = boundaryintcells[eid];
      DomainIntCells   & dics = domainintcells  [eid];

      LINALG::Matrix<3,1> physCoordCorner;
      LINALG::Matrix<3,1> eleCoordDomainCorner;
      LINALG::Matrix<3,1> eleCoordBoundaryCorner;

      for ( std::set<GEO::CUT::BoundaryCell*>::iterator i=bcells.begin();
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

      for ( std::set<GEO::CUT::IntegrationCell*>::iterator i=cells.begin(); i!=cells.end(); ++i )
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

  // cleanup

  int localcells = domainintcells.size();
  int globalcells;
  xfemdis->Comm().SumAll( &localcells, &globalcells, 1 );

  const double t_end = Teuchos::Time::wallTime()-t_start;
  if ( xfemdis->Comm().MyPID() == 0 )
  {
    std::cout << " Success (" << t_end  <<  " secs), intersected elements: " << globalcells;
    std::cout << endl;
  }
}

void GEO::CellGenerator::Generate( GEO::CUT::Element* element, const tetgenio & out )
{
  if ( out.numberoftetrahedra==0 )
    dserror( "cut element without cut cells" );

  BoundaryIntCells & bics = boundaryintcells_[element->Id()];
  DomainIntCells   & dics = domainintcells_  [element->Id()];

  LINALG::SerialDenseMatrix eleDomainCoord(3, 3); // in xfem parent element domain
  LINALG::SerialDenseMatrix eleBoundaryCoord(3, 3);
  LINALG::SerialDenseMatrix physDomainCoord(3, 3); // in physical domain

  LINALG::Matrix<3,1> physCoordCorner;
  LINALG::Matrix<3,1> eleCoordDomainCorner;
  LINALG::Matrix<3,1> eleCoordBoundaryCorner;

  for ( int i=0; i<out.numberoftrifaces; ++i )
  {
    if ( out.trifacemarkerlist[i] > -1 )
    {
      GEO::CUT::Side * cut_side = intersection_.GetCutSide( out.trifacemarkerlist[i] );
      for ( int j=0; j<3; ++j )
      {
        int pointidx = out.trifacelist[i*3+j] * 3;
        std::copy( &out.pointlist[pointidx], &out.pointlist[pointidx+3], physCoordCorner.A() );
        physCoordCorner.Scale( 1./TETGENPOINTSCALE );

        element ->LocalCoordinates( physCoordCorner, eleCoordDomainCorner );
        cut_side->LocalCoordinates( physCoordCorner, eleCoordBoundaryCorner );

        std::copy( physCoordCorner       .A(), physCoordCorner       .A()+3, &physDomainCoord ( 0, j ) );
        std::copy( eleCoordDomainCorner  .A(), eleCoordDomainCorner  .A()+3, &eleDomainCoord  ( 0, j ) );
        std::copy( eleCoordBoundaryCorner.A(), eleCoordBoundaryCorner.A()+3, &eleBoundaryCoord( 0, j ) );
      }
      bics.push_back( BoundaryIntCell( DRT::Element::tri3,
                                       cut_side->Id(),
                                       eleDomainCoord,
                                       eleBoundaryCoord,
                                       physDomainCoord ) );
    }
  }

//   if ( bics.size()==0 )
//     dserror( "cut element without cut surface" );

  DRT::Element::DiscretizationType distype = DRT::Element::tet4;

  const int numTetNodes = DRT::UTILS::getNumberOfElementNodes(distype);
  if ( out.numberofcorners < numTetNodes )
  {
    dserror( "quadratic tets?" );
  }

  LINALG::SerialDenseMatrix tetrahedronCoord( 3, numTetNodes );
  LINALG::SerialDenseMatrix physTetrahedronCoord( 3, numTetNodes );

  for ( int i=0; i<out.numberoftetrahedra; ++i )
  {
    for ( int j=0; j<numTetNodes; ++j )
    {
      int pointidx = out.tetrahedronlist[i*out.numberofcorners+j] * 3;
      std::copy( &out.pointlist[pointidx], &out.pointlist[pointidx+3], physCoordCorner.A() );
      physCoordCorner.Scale( 1./TETGENPOINTSCALE );

      element ->LocalCoordinates( physCoordCorner, eleCoordDomainCorner );

      std::copy( physCoordCorner     .A(), physCoordCorner     .A()+3, &physTetrahedronCoord( 0, j ) );
      std::copy( eleCoordDomainCorner.A(), eleCoordDomainCorner.A()+3, &tetrahedronCoord    ( 0, j ) );
    }

    // create planes consisting of 3 nodes each
    LINALG::Matrix<3,1> p0( tetrahedronCoord.A()  , true );
    LINALG::Matrix<3,1> p1( tetrahedronCoord.A()+3, true );
    LINALG::Matrix<3,1> p2( tetrahedronCoord.A()+6, true );
    LINALG::Matrix<3,1> p3( tetrahedronCoord.A()+9, true );

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
        std::swap( tetrahedronCoord    ( i, 1 ), tetrahedronCoord    ( i, 3 ) );
        std::swap( physTetrahedronCoord( i, 1 ), physTetrahedronCoord( i, 3 ) );
      }
    }

    dics.push_back( DomainIntCell( distype, tetrahedronCoord, physTetrahedronCoord ) );
  }
}

#endif
