#include "cut_point.H"
#include "cut_volumecell.H"
#include "cut_boundarycell.H"
#include "cut_integrationcell.H"
#include "cut_tetmesh.H"
#include "cut_mesh.H"
#include "cut_options.H"
#include "cut_kernel.H"
#include "cut_triangulateFacet.H"
#include "volume_integration.H"
#include "boundarycell_integration.H"
#include "direct_divergence.H"

#include<algorithm>


int GEO::CUT::VolumeCell::hex8totet4[5][4] = {
  {0, 1, 3, 4},
  {1, 2, 3, 6},
  {4, 5, 1, 6},
  {6, 7, 3, 4},
  {1, 6, 3, 4}
};

int GEO::CUT::VolumeCell::wedge6totet4[3][4] = {
  {0, 1, 2, 3},
  {3, 4, 1, 5},
  {1, 5, 2, 3}
};


int GEO::CUT::VolumeCell::pyramid5totet4[2][4] = {
  {0, 1, 3, 4},
  {1, 2, 3, 4}
};


GEO::CUT::VolumeCell::VolumeCell( const plain_facet_set & facets,
                                  const std::map<std::pair<Point*, Point*>, plain_facet_set > & volume_lines,
                                  Element * element )
  : element_( element ),
    position_( Point::undecided ),
    facets_( facets )
{
  for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Register( this );
  }
}

void GEO::CUT::VolumeCell::Neighbors( Point * p,
                                      const plain_volumecell_set & cells,
                                      const plain_volumecell_set & done,
                                      plain_volumecell_set & connected,
                                      plain_element_set & elements )
{
  if ( done.count( this )==0 )
  {
    // this volume is included
    connected.insert( this );
    elements.insert( element_ );

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      if ( p==NULL or f->Contains( p ) )
      {
        f->Neighbors( p, cells, done, connected, elements );
      }
    }

    if ( p!=NULL )
    {
      for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
      {
        Facet * f = *i;
        if ( not f->Contains( p ) )
        {
          f->Neighbors( p, cells, done, connected, elements );
        }
      }
    }
  }
}

// without check for elements
void GEO::CUT::VolumeCell::Neighbors( Point * p,
                                      const plain_volumecell_set & cells,
                                      const plain_volumecell_set & done,
                                      plain_volumecell_set & connected)
{
  if ( done.count( this )==0 )
  {
    // this volume is included
    connected.insert( this );

    // Do the facets that include the point first. This ensures we choose the
    // right volumes (the ones attached to the point), if there are multiple
    // connections possible (we are faced with a thin structure cut.)

    for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      if ( p==NULL or f->Contains( p ) )
      {
        f->Neighbors( p, cells, done, connected);
      }
    }

    if ( p!=NULL )
    {
      for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
      {
        Facet * f = *i;
        if ( not f->Contains( p ) )
        {
          f->Neighbors( p, cells, done, connected);
        }
      }
    }
  }
}


void GEO::CUT::VolumeCell::GetAllPoints( Mesh & mesh, PointSet & cut_points )
{
  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->GetAllPoints( mesh, cut_points );
  }
}

bool GEO::CUT::VolumeCell::Contains( Point * p )
{
  for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->Contains( p ) )
    {
      return true;
    }
  }
  return false;
}

bool GEO::CUT::VolumeCell::Contains( LINALG::Matrix<3,1>& x)
{

  if(integrationcells_.size() == 0) dserror("no integrationcells for volumecell stored, implement Contains check without integrationcells");

  for(GEO::CUT::plain_integrationcell_set::iterator it=integrationcells_.begin(); it!=integrationcells_.end(); it++)
  {
    GEO::CUT::IntegrationCell* intcell = *it;

    if(intcell->Contains(x)) return true;
  }

  return false;
}


void GEO::CUT::VolumeCell::CreateTet4IntegrationCells( Mesh & mesh,
                                                       const std::vector<std::vector<Point*> > & tets,
                                                       const std::map<Facet*, std::vector<Point*> > & sides_xyz )
{
  for ( std::vector<std::vector<Point*> >::const_iterator i=tets.begin();
        i!=tets.end();
        ++i )
  {
    const std::vector<Point*> & tet = *i;
    if ( tet.size()!=4 )
    {
      throw std::runtime_error( "tet expected" );
    }
    NewTet4Cell( mesh, tet );
  }

  for ( std::map<Facet*, std::vector<Point*> >::const_iterator i=sides_xyz.begin();
        i!=sides_xyz.end();
        ++i )
  {
    Facet * f = i->first;
    const std::vector<Point*> & points = i->second;

    std::size_t length = points.size();
    if ( length % 3 != 0 )
      throw std::runtime_error( "expect list of triangles" );

    length /= 3;
    std::vector<Point*> p( 3 );
    for ( std::size_t i=0; i<length; ++i ) // loop the list of triangles
    {
      std::copy( &points[3*i], &points[3*( i+1 )], &p[0] );
      //Tri3BoundaryCell::CreateCell( mesh, this, f, p );
      NewTri3Cell( mesh, f, p ); // create tri3 cell
    }
  }
}

void GEO::CUT::VolumeCell::GetIntegrationCells( plain_integrationcell_set & cells )
{
  std::copy( integrationcells_.begin(), integrationcells_.end(), std::inserter( cells, cells.begin() ) );
}

void GEO::CUT::VolumeCell::GetBoundaryCells( std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells )
{
  for ( plain_boundarycell_set::iterator i=bcells_.begin(); i!=bcells_.end(); ++i )
  {
    BoundaryCell * bc = *i;
    Facet * f = bc->GetFacet();
    int sid = f->SideId();
    if ( sid > -1 )
    {
      bcells[sid].push_back( bc );
    }
  }
}

void GEO::CUT::VolumeCell::ConnectNodalDOFSets( bool include_inner )
{
//   if ( Empty() )
//     return;
  if ( not include_inner and Position()!=Point::outside )
    return;

  const std::vector<Node*> & nodes = element_->Nodes();
  nodaldofset_.reserve( nodes.size() );

  for ( std::vector<Node*>::const_iterator i=nodes.begin();
        i!=nodes.end();
        ++i )
  {
    Node * n = *i;
    nodaldofset_.push_back( n->DofSetNumber( this ) );
  }
}

void GEO::CUT::VolumeCell::Position( Point::PointPosition position )
{
  if ( position_ != position )
  {
    position_ = position;

    for ( plain_facet_set::const_iterator i=facets_.begin(); i!=facets_.end(); ++i )
    {
      Facet * f = *i;
      Point::PointPosition fp = f->Position();
      if ( fp==Point::undecided )
      {
        f->Position( position );
      }
    }
  }
}

void GEO::CUT::VolumeCell::Print( std::ostream & stream )
{
  stream << "# VolumeCell: "
         << " pos: "      << position_ << " "
         << "#facets: "   << facets_.size() << " "
         << "#intcells: " << integrationcells_.size() << " "
         << "#bcells: "   << bcells_.size()
         << "\n";
  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->Print( stream );
  }
}

void GEO::CUT::VolumeCell::NewBoundaryCell( Mesh & mesh, DRT::Element::DiscretizationType shape, Facet * f, const std::vector<Point*> & x )
{
  if ( facets_.count( f )==0 )
  {
    throw std::runtime_error( "facet does not belong to volume cell" );
  }
  switch ( shape )
  {
  case DRT::Element::tri3:
    NewTri3Cell( mesh, f, x );
    break;
  case DRT::Element::quad4:
    NewQuad4Cell( mesh, f, x );
    break;
  default:
    throw std::runtime_error( "unknown shape" );
  }
}

void GEO::CUT::VolumeCell::NewTri3Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x )
{
  f->NewTri3Cell( mesh, this, x, bcells_ );
}

void GEO::CUT::VolumeCell::NewQuad4Cell( Mesh & mesh, Facet * f, const std::vector<Point*> & x )
{
  f->NewQuad4Cell( mesh, this, x, bcells_ );
}

void GEO::CUT::VolumeCell::NewArbitraryCell( Mesh & mesh, Facet * f, const std::vector<Point*> & x,
    const DRT::UTILS::GaussIntegration& gp, const LINALG::Matrix<3,1>& normal )
{
  f->NewArbitraryCell( mesh, this, x, bcells_, gp, normal );
}

/*double GEO::CUT::VolumeCell::Volume()
{
  double volume = 0;
  for ( plain_integrationcell_set::iterator i=integrationcells_.begin(); i!=integrationcells_.end(); ++i )
  {
    IntegrationCell * ic = *i;
    volume += ic->Volume();
  }
  return volume;
}*/

int GEO::CUT::VolumeCell::NumGaussPoints( DRT::Element::DiscretizationType shape )
{
  int numgp = 0;

  for ( plain_integrationcell_set::const_iterator i=integrationcells_.begin(); i!=integrationcells_.end(); ++i )
  {
    IntegrationCell * ic = *i;

    // Create (unmodified) gauss points for integration cell with requested
    // polynomial order. This is supposed to be fast, since there is a cache.
    DRT::UTILS::GaussIntegration gi( ic->Shape(), ic->CubatureDegree( shape ) );

    // we just need the number of points per cell
    numgp += gi.NumPoints();
  }

  return numgp;
}

void GEO::CUT::VolumeCell::Disconnect()
{
  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    f->DisconnectVolume( this );
  }
}

void GEO::CUT::VolumeCell::NewIntegrationCell( Mesh & mesh, DRT::Element::DiscretizationType shape, const std::vector<Point*> & x )
{
  switch ( shape )
  {
  case DRT::Element::hex8:
    NewHex8Cell( mesh, x );
    break;
  case DRT::Element::tet4:
    NewTet4Cell( mesh, x );
    break;
  case DRT::Element::wedge6:
    NewWedge6Cell( mesh, x );
    break;
  case DRT::Element::pyramid5:
    NewPyramid5Cell( mesh, x );
    break;
  default:
    throw std::runtime_error( "unknown shape" );
  }
}

void GEO::CUT::VolumeCell::NewHex8Cell( Mesh & mesh, const std::vector<Point*> & points )
{
  Point::PointPosition position = Position();
  if ( mesh.CreateOptions().GenHex8() )
  {
    integrationcells_.insert( mesh.NewHex8Cell( position, points, this ) );
  }
  else
  {
    std::vector<Point*> tet4_points( 4 );
    for ( int i=0; i<5; ++i )
    {
      SetTetPoints( hex8totet4[i], points, tet4_points );
      integrationcells_.insert( mesh.NewTet4Cell( position, tet4_points, this ) );
    }
  }
}

GEO::CUT::IntegrationCell * GEO::CUT::VolumeCell::NewTet4Cell( Mesh & mesh, const std::vector<Point*> & points )
{
  Point::PointPosition position = Position();
  IntegrationCell * ic = mesh.NewTet4Cell( position, points, this );
  integrationcells_.insert( ic );
  return ic;
}

void GEO::CUT::VolumeCell::NewWedge6Cell( Mesh & mesh, const std::vector<Point*> & points )
{
  Point::PointPosition position = Position();
  if ( mesh.CreateOptions().GenWedge6() )
  {
    integrationcells_.insert( mesh.NewWedge6Cell( position, points, this ) );
  }
  else
  {
    std::vector<Point*> tet4_points( 4 );
    for ( int i=0; i<3; ++i )
    {
      SetTetPoints( wedge6totet4[i], points, tet4_points );
      integrationcells_.insert( mesh.NewTet4Cell( position, tet4_points, this ) );
    }
  }
}

void GEO::CUT::VolumeCell::NewPyramid5Cell( Mesh & mesh, const std::vector<Point*> & points )
{
  Point::PointPosition position = Position();
  if ( mesh.CreateOptions().GenPyramid5() )
  {
    integrationcells_.insert( mesh.NewPyramid5Cell( position, points, this ) );
  }
  else
  {
    std::vector<Point*> tet4_points( 4 );
    for ( int i=0; i<2; ++i )
    {
      SetTetPoints( pyramid5totet4[i], points, tet4_points );
      integrationcells_.insert( mesh.NewTet4Cell( position, tet4_points, this ) );
    }
  }
}

void GEO::CUT::VolumeCell::SimplifyIntegrationCells( Mesh & mesh )
{
  // do whatever can be done to get simpler cells
  //

  std::map<int, std::vector<Facet*> > side_facets;

  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;
    if ( f->OnCutSide() )
    {
      side_facets[f->SideId()].push_back( f );
    }
  }

  for ( std::map<int, std::vector<Facet*> >::iterator i=side_facets.begin();
        i!=side_facets.end();
        ++i )
  {
    int sideid = i->first;
    std::vector<Facet*> & facets = i->second;
    std::vector<BoundaryCell*> bcs;
    sorted_vector<std::pair<Point*, Point*> > lines;
    for ( plain_boundarycell_set::iterator i=bcells_.begin(); i!=bcells_.end(); ++i )
    {
      BoundaryCell * bc = *i;
      if ( bc->GetFacet()->SideId()==sideid )
      {
        const Cycle & cycle = bc->PointCycle();
        cycle.Add( lines );
        bcs.push_back( bc );
      }
    }
    if ( bcs.size() > 1 )
    {
      Cycle cycle;
      if ( Cycle::MakeCycle( lines, cycle ) )
      {
        std::vector<Point*> corner_points;
        DRT::Element::DiscretizationType shape = KERNEL::CalculateShape( cycle(), corner_points );

        if ( shape!=DRT::Element::dis_none )
        {
          for ( std::vector<BoundaryCell*>::iterator i=bcs.begin(); i!=bcs.end(); ++i )
          {
            BoundaryCell * bc = *i;
            bcells_.erase( bc );
            bc->Clear();
          }
          switch ( shape )
          {
          case DRT::Element::quad4:
            // the facet is too small, but it knows the right side
            if ( mesh.CreateOptions().GenQuad4() )
            {
              mesh.NewQuad4Cell( this, facets[0], corner_points );
            }
            else
            {
              std::vector<Point*> tri3_points = corner_points;
              tri3_points.pop_back();
              mesh.NewTri3Cell( this, facets[0], tri3_points );
              tri3_points.erase( tri3_points.begin()+1 );
              tri3_points.push_back( corner_points.back() );
              mesh.NewTri3Cell( this, facets[0], tri3_points );
            }
            break;
          case DRT::Element::tri3:
            // the facet is too small, but it knows the right side
            mesh.NewTri3Cell( this, facets[0], corner_points );
            break;
          default:
            throw std::runtime_error( "unsupported boundary cell type" );
          }
        }
#if 0
        std::cout << "found cycle with " << cycle.size()
                  << " points on cut side " << sideid
                  << " out of " << numbc
                  << " boundary cells: shape=" << shape
                  << " with " << line_points.size()
                  << " points\n";
#endif
      }
    }
  }
}


/*--------------------------------------------------------------------*
 * Check wheter the point is inside, outside or on the boundary
 * of this volumecelll                                    sudhakar 07/12
 *--------------------------------------------------------------------*/
std::string GEO::CUT::VolumeCell::IsThisPointInside( Point *pt )
{
  LINALG::Matrix<3,1> xglo;
  pt->Coordinates(xglo.A());
  std::string inside = IsThisPointInside( xglo );
  return inside;
}

/*-----------------------------------------------------------------------------------------------*
 * Check whether the point with this global coordinates is inside, outside or on the boundary
 * of this volumecell                                                               sudhakar 07/12
 *-----------------------------------------------------------------------------------------------*/
std::string GEO::CUT::VolumeCell::IsThisPointInside( LINALG::Matrix<3,1>& xglo )
{
  LINALG::Matrix<3,1> xloc;
  element_->LocalCoordinates( xglo, xloc );

  const GEO::CUT::Point::PointPosition posi = Position();
  if( posi==0 )
    dserror( "undefined position for the volumecell" );

  VolumeIntegration vc(this,element_,posi,0);
  std::string inside = vc.IsPointInside( xloc );
  return inside;
}

void GEO::CUT::VolumeCell::TestSurface()
{
  if ( Empty() )
  {
    // This is an artificial cell with zero volume. It should not exist in the
    // first place.
    return;
  }

  // see if all lines are closed
  //
  // This finds all the degenerated cases that where dropped before. Thus the
  // test complains a lot.

  for ( plain_facet_set::iterator i=facets_.begin(); i!=facets_.end(); ++i )
  {
    Facet * f = *i;

    if ( f->OnCutSide() )
    {
      if ( f->IsTriangulated() )
      {
        //
      }
      if ( f->HasHoles() )
      {
        //
      }

      point_line_set lines;

      const std::vector<Point*> & points = f->Points();
      Cycle cycle( points );
      cycle.Add( lines );

      for ( plain_boundarycell_set::iterator i=bcells_.begin(); i!=bcells_.end(); ++i )
      {
        BoundaryCell * bc = *i;
        if ( bc->GetFacet() == f )
        {
          const std::vector<Point*> & points = bc->Points();
          Cycle cycle( points );
          cycle.Add( lines );
        }
      }

      if ( lines.size()!=0 )
      {
        throw std::runtime_error( "volume cut facets not closed" );
      }
    }
  }
}

/*-------------------------------------------------------------------------------------*
                Write the volumecell details for visualization
                Gausspoints of moment fitting are not included
*--------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DumpGmsh( std::ofstream& file )
{
    const plain_facet_set & facete = Facets();

    file<<"View \"Volume Cell \" {\n";
    for( unsigned j=0;j<facete.size();j++ )
    {
      Facet * ref = facete[j];
      const std::vector<std::vector<double> > corners = ref->CornerPointsLocal(ParentElement());
      for( unsigned i=0;i<corners.size();i++ )
      {
        const std::vector<double> coords1 = corners[i];
        const std::vector<double> coords2 = corners[(i+1)%corners.size()];
        file<<"SL("<<coords1[0]<<","<<coords1[1]<<","<<coords1[2]<<","<<
            coords2[0]<<","<<coords2[1]<<","<<coords2[2]<<")"<<"{0,0};\n";
      }
    }
    file<<"};\n";
    file<<"View[PostProcessing.NbViews-1].ColorTable = { {0,0,255} };\n"; // Changing color to red
    file<<"View[PostProcessing.NbViews-1].Light=0;\n";    // Disable the lighting
    file<<"View[PostProcessing.NbViews-1].ShowScale=0;\n";  // Disable legend
    file<<"View[PostProcessing.NbViews-1].LineWidth = 3.0;"; // increase line width
}

/*--------------------------------------------------------------------------------------------------------*
        write the boundaries of volumecell and the positions of Gauss points for visualization
        a separate file with "side" prefix is generated for every volumecell as the gausspoint
        distribution can be clearly seen
*---------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DumpGmshGaussPoints(const std::vector<std::vector<double> >&gauspts)
{

  static int sideno = 0;
  sideno++;

  std::stringstream str;
  str << "side" << sideno << ".pos";
  std::ofstream file( str.str().c_str() );

  DumpGmsh(file);

  file<<"Geometry.PointSize=6.0;\n";      // Increase the point size
  int pointno=1;
  for(unsigned i=0;i<gauspts.size();i++)
  {
     file<<"Point("<<pointno<<")={"<<gauspts[i][0]<<","<<gauspts[i][1]<<","<<gauspts[i][2]<<","<<"1"<<"};"<<std::endl;
     pointno++;
  }
  file.close();
}

void GEO::CUT::VolumeCell::integrateSpecificFunctionsTessellation()
{

  Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
                  Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

  const plain_integrationcell_set & cells = IntegrationCells();
  for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
  {
    GEO::CUT::IntegrationCell * ic = *i;
    switch ( ic->Shape() )
    {
    case DRT::Element::hex8:
    {
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>( ic );
      gpc->Append( gp );
      break;
    }
    case DRT::Element::tet4:
    {
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>( ic );
      gpc->Append( gp );
      break;
    }
    default:
    {
      dserror("Include this element here");
      break;
    }
    }
  }

  DRT::UTILS::GaussIntegration gpv( gpc );

  double intVal = 0.0;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=gpv.begin(); iquad!=gpv.end(); ++iquad )
  {
    double weight = iquad.Weight();

    const LINALG::Matrix<3,1> eta( iquad.Point() );
    double xx = eta(0,0);
    double yy = eta(1,0);
    double zz = eta(2,0);

    intVal += (pow(xx,6)+xx*pow(yy,4)*zz+xx*xx*yy*yy*zz*zz+pow(zz,6))*weight;
  }
  std::cout<<setprecision(20)<<"TESSELLATION Integration = "<<intVal<<"\n";

}

template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::VolumeCell::CreateProjected( GEO::CUT::IntegrationCell * ic )
{
  const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<3, nen> xie;

  const std::vector<GEO::CUT::Point*> & cpoints = ic->Points();
  if ( cpoints.size() != nen )
    throw std::runtime_error( "non-matching number of points" );

  for ( unsigned i=0; i<nen; ++i )
  {
    GEO::CUT::Point * p = cpoints[i];
    LINALG::Matrix<3,1> xg,xi;
    p->Coordinates(xg.A());
    element_->LocalCoordinates(xg,xi);
    std::copy( xi.A(), xi.A()+3, &xie( 0, i ) );
  }

  Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
    DRT::UTILS::GaussIntegration::CreateProjected<distype>( xie, ic->CubatureDegree( element_->Shape() ) );
  return gp;
}

/*------------------------------------------------------------------------------------------------------*
    convert the Gaussian points and weights into appropriate Gauss rule as per BACI implementation
*-------------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::VolumeCell::GaussPointsFitting()
{
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new
                      DRT::UTILS::CollectedGaussPoints( 0 ) );

  for(unsigned i=0;i<gausPts_.size();i++)
  {
    LINALG::Matrix<3,1> xe,xei;
    xe(0,0) = gausPts_[i][0];
    xe(1,0) = gausPts_[i][1];
    xe(2,0) = gausPts_[i][2];

    cgp->Append( xe, weights_(i) );
  }

  return cgp;
}

/*--------------------------------------------------------------------------------------------*
                 Generate boundary cells for the cut facets of the volumecell
*---------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::GenerateBoundaryCells( Mesh &mesh,
                                                  const GEO::CUT::Point::PointPosition posi,
                                                  Element *elem,
                                                  int BaseNos,
                                                  std::string BCellgausstype )
{
  const plain_facet_set & facete = Facets();
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
    Facet *fac = *i;

    if( fac->OnCutSide() == false )             //we need boundary cells only for the cut facets
      continue;

    //--------------------------------------------------------------------
    // Normal vector from parent side is used to identify whether normals
    // from facet is in appropriate direction or not
    //--------------------------------------------------------------------
    const Side* parside = fac->ParentSide();
    const std::vector<Node*> &par_nodes = parside->Nodes();
    std::vector<Point*> parpts(3);

    parpts[0] = par_nodes[0]->point();
    parpts[1] = par_nodes[1]->point();
    parpts[2] = par_nodes[2]->point();

    std::vector<double> eqnpar(4),eqnfac(4);
    // equation of plane denotes normal direction
    eqnpar = KERNEL::EqnPlane( parpts[0], parpts[1], parpts[2] );

    std::vector<Point*> corners = fac->CornerPoints();
    std::vector<Point*> cornersTemp (corners);

    // when finding eqn of plane for the facet, inline points should not be taken
   CUT::KERNEL::DeleteInlinePts( cornersTemp );

   bool rever = false;
   if( cornersTemp.size()!=0 )
   {
     eqnfac = KERNEL::EqnPlanePolygon( cornersTemp );
     rever = ToReverse( posi, eqnpar, eqnfac );
   }

   if(rever)                                       // normal from facet is in wrong direction
   {
     std::reverse(corners.begin(),corners.end());  // change ordering to correct this
     std::reverse(cornersTemp.begin(),cornersTemp.end());
   }

   //if no of corners are 3 or 4, just add them as boundary integrationcells directly
   if(corners.size()==3)
   {
     NewTri3Cell( mesh, fac, corners );
   }
   else if(corners.size()==4)
   {
     NewQuad4Cell(mesh,fac,corners);
   }
   else
   {
      if(BCellgausstype=="Tessellation")//generate boundarycell gausspoints by triangulation
      {
#if 1 // create only triangles - result in more number of Gauss points
        // Use "corners" for triangulation - no points are deleted
        if(!fac->IsTriangulated())
          fac->DoTriangulation( mesh, corners );
        const std::vector<std::vector<Point*> > & triangulation = fac->Triangulation();
#else  // creates both tri and quad. less no of Gauss points - but deleted some points leads to error

        if( cornersTemp.size()==0 )
              continue;
        // Use "cornersTemp" for triangulation - deleted some points results in error
        if( !fac->IsFacetSplit() )
          fac->SplitFacet(  cornersTemp );

        const std::vector<std::vector<Point*> > triangulation = fac->GetSplitCells();
#endif

        for ( std::vector<std::vector<Point*> >::const_iterator j=triangulation.begin();
                          j!=triangulation.end(); ++j )
        {
          std::vector<Point*> tri = *j;

          if(tri.size()==3)
            NewTri3Cell(mesh,fac,tri);
          else if(tri.size()==4)
            NewQuad4Cell(mesh,fac,tri);
          else
            dserror("Triangulation created neither tri3 or quad4");
        }
     }

     else if(BCellgausstype=="MomentFitting")//generate boundarycell gausspoints by solving moment fitting equations
     {
        BoundarycellIntegration bcell_inte(elem,fac,posi,BaseNos);
        Bcellweights_ = bcell_inte.GenerateBoundaryCellIntegrationRule();
        BcellgausPts_ = bcell_inte.getBcellGaussPointLocation();

        //the boundarycell integration is carriedout in the local coord of the element
        //to project the coordinates of Gauss points, shape functions of element can be used
        //
        //                                            area of facet in global coordinates
        //but to transform the weight, the jacobian = -----------------------------------
        //                                            area of facet in local coordinates
        FacetIntegration bcellLocal(fac,elem,posi,true,false);
        bcellLocal.set_integ_number(1);
        double areaLocal = bcellLocal.integrate_facet();

        FacetIntegration bcellGlobal(fac,elem,posi,true,true);
        bcellGlobal.set_integ_number(1);
        double areaGlobal = bcellGlobal.integrate_facet();
        double jaco = areaGlobal/areaLocal;

        int numBcellpts = BcellgausPts_.size();
        Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints( numBcellpts ) );

        LINALG::Matrix<3,1> xeLocal,xeGlobal;
        for(unsigned i=0;i<BcellgausPts_.size();i++)
        {
          xeLocal(0,0) = BcellgausPts_[i][0];
          xeLocal(1,0) = BcellgausPts_[i][1];
          xeLocal(2,0) = BcellgausPts_[i][2];

          elem->GlobalCoordinates( xeLocal, xeGlobal );

          cgp->Append( xeGlobal, Bcellweights_(i)*jaco );
        }

        LINALG::Matrix<3,1> normal;
        double normalFac;
        if(rever)
        {
          //std::reverse(corners.begin(),corners.end());
          normalFac = -1.0;
        }
        else
          normalFac = 1.0;

        normalFac = normalFac*sqrt(eqnfac[0]*eqnfac[0]+eqnfac[1]*eqnfac[1]+eqnfac[2]*eqnfac[2]);
        for(unsigned i=0;i<3;i++)
          normal(i,0) = eqnfac[i]/normalFac;

        DRT::UTILS::GaussIntegration gi(cgp);
        NewArbitraryCell(mesh, fac, corners, gi, normal);
      }
    }
  }
}

/*--------------------------------------------------------------------------------------------------------*
    This is to check whether the corner points of the cut side facet is aligned to give outward normal
*---------------------------------------------------------------------------------------------------------*/
bool GEO::CUT::VolumeCell::ToReverse(const GEO::CUT::Point::PointPosition posi,
                                     std::vector<double> parEqn,
                                     std::vector<double> facetEqn)
{
	bool rever = false;

	// position is inside
	if(posi==-3)
	{
    if(fabs(parEqn[0])>TOL_EQN_PLANE && parEqn[0]*facetEqn[0]>0.0)
      rever = true;
    else if(fabs(parEqn[1])>TOL_EQN_PLANE && parEqn[1]*facetEqn[1]>0.0)
      rever = true;
    else if(fabs(parEqn[2])>TOL_EQN_PLANE && parEqn[2]*facetEqn[2]>0.0)
      rever = true;
    else
      rever = false;
	}

	// position is outside
	else if(posi==-2)
	{
		if(fabs(parEqn[0])>TOL_EQN_PLANE && parEqn[0]*facetEqn[0]<0.0)
      rever = true;
		else if(fabs(parEqn[1])>TOL_EQN_PLANE && parEqn[1]*facetEqn[1]<0.0)
      rever = true;
		else if(fabs(parEqn[2])>TOL_EQN_PLANE && parEqn[2]*facetEqn[2]<0.0)
      rever = true;
		else
			rever = false;
	}
	return rever;
}

/*------------------------------------------------------------------------------------------*
   When DirectDivergence method is used for gauss point generation, for every gauss point
   on the facet, an internal gauss rule is to be generated to find the modified integrand
*-------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::GenerateInternalGaussRule()
{
  DRT::UTILS::GaussIntegration grule(gp_);

  intGP_.resize( grule.NumPoints(), grule );

  int num = 0;
  for ( DRT::UTILS::GaussIntegration::iterator quadint=grule.begin(); quadint!=grule.end(); ++quadint )
  {
    const LINALG::Matrix<3,1> etaFacet( quadint.Point() );  //coordinates and weight of main gauss point
    LINALG::Matrix<3,1> intpt( etaFacet );

    DRT::UTILS::GaussIntegration gi( DRT::Element::line2, 6 ); //internal gauss rule for interval (-1,1)

    Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new
                         DRT::UTILS::CollectedGaussPoints( 0 ) );

    //x-coordinate of main Gauss point is projected in the reference plane
    double xbegin = (RefEqnPlane_[3]-RefEqnPlane_[1]*etaFacet(1,0)-
                    RefEqnPlane_[2]*etaFacet(2,0))/RefEqnPlane_[0];

    double jac = fabs(xbegin-etaFacet(0,0))*0.5; // jacobian for 1D transformation rule

    // -----------------------------------------------------------------------------
    // project internal gauss point from interval (-1,1) to the actual interval
    // -----------------------------------------------------------------------------
    for ( DRT::UTILS::GaussIntegration::iterator iqu=gi.begin(); iqu!=gi.end(); ++iqu )
    {
      const LINALG::Matrix<1,1> eta( iqu.Point() );
      double weight = iqu.Weight();

      double xmid = 0.5*(xbegin+etaFacet(0,0));
      intpt(0,0) = (xmid-xbegin)*eta(0,0)+xmid;    // location of internal gauss points

      weight = weight*jac;                         // weight of internal gauss points
      if( xbegin>etaFacet(0,0) )
        weight = -1.0*weight;

      cgp->Append( intpt, weight );
    }

    DRT::UTILS::GaussIntegration gint(cgp);

    intGP_[num] = gint;
    num++;
  }

  if( grule.NumPoints() != num )
    dserror( "some facet points missed?" );
}

/*------------------------------------------------------------------------------------------*
   Moment fitting equations are solved at each volume cell to construct integration rules
*-------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::MomentFitGaussWeights(Element *elem,
                                                 Mesh & mesh,
                                                 bool include_inner,
                                                 std::string BCellgausstype)
{

	//position is used to decide whether the ordering of points are in clockwise or not
	const GEO::CUT::Point::PointPosition posi = Position();

	if( posi==0 )
	  dserror( "undefined position for the volumecell" );

	//if the volumecell is inside and includeinner is false, no need to compute the Gaussian points
	//as this vc will never be computed in xfem algorithm
	if(posi==-2 && include_inner==false)
		return;

	int BaseNos=84;                                     // number of base functions to be used in the integration
  VolumeIntegration vc_inte(this,elem,posi,BaseNos);

  weights_ = vc_inte.compute_weights();              // obtain the integration weight at all points
  gausPts_ = vc_inte.getGaussPointLocation();        // get the coordinates of all the Gauss points

  gp_ = GaussPointsFitting();                        // convert the weight and the location to Gauss rule

  // generate boundary cells -- when using tessellation this is automatically done
  GenerateBoundaryCells( mesh, posi, elem, BaseNos, BCellgausstype );

  //std::cout<<"MOMENT FITTING ::: Number of points = "<<weights_.Length()<<"\n";
}

/*---------------------------------------------------------------------------------------------------------------*
                     The facets that have non-zero x-component normal is triangulated.              sudhakar 03/12
                  The gauss integration rules are generated by applying divergence theorem
      The reference facet is identified which will be used to find the modified integral in fluid integration
*----------------------------------------------------------------------------------------------------------------*/
void GEO::CUT::VolumeCell::DirectDivergenceGaussRule( Element *elem,
                                                      Mesh & mesh,
                                                      bool include_inner,
                                                      std::string BCellgausstype )
{

  //position is used to decide whether the ordering of points are in clockwise or not
  const GEO::CUT::Point::PointPosition posi = Position();

  if( posi==0 )
    dserror( "undefined position for the volumecell" );

  //if the volumecell is inside and includeinner is false, no need to compute the Gaussian points
  //as this vc will never be computed in xfem algorithm
  if(posi == Point::inside && include_inner==false)
    return;

  DirectDivergence dd(this,elem,posi,mesh);

  RefEqnPlane_.reserve(4);                   //it has to store a,b,c,d in ax+by+cz=d
  gp_ = dd.VCIntegrationRule( RefEqnPlane_ );// compute main gauss points

  GenerateInternalGaussRule();               // compute internal gauss points for every main gauss point

  // compute volume of this cell
  // also check whether generated gauss rule predicts volume accurately
  DRT::UTILS::GaussIntegration gpi(gp_);
  dd.DebugVolume( gpi, RefEqnPlane_, intGP_ );

#if 0 // integrate a predefined function
  dd.IntegrateSpecificFuntions( gpi, RefEqnPlane_, intGP_ );
#endif

#ifdef DEBUGCUTLIBRARY  // write volumecell, main and internal Gauss points
  dd.DivengenceCellsGMSH( gpi, intGP_ );
#endif

  // generate boundary cells -- when using tessellation this is automatically done
  GenerateBoundaryCells( mesh, posi, elem, 0, "Tessellation" );
}

/*-------------------------------------------------------------------------------------*
| Return Ids of all the points associated with this volumecell           shahmiri 06/12
*--------------------------------------------------------------------------------------*/
std::set<int> GEO::CUT::VolumeCell::VolumeCellPointIds()
{
  if ( vcpoints_ids_.size() != 0)
  {
    return vcpoints_ids_;
  }
  else
  {
    const plain_facet_set & facete = Facets();

    // loop over facets
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
      Facet *fe = *i;
      std::vector<Point*> corners = fe->CornerPoints();

      for(std::vector<Point*>::const_iterator c=corners.begin(); c!=corners.end(); c++)
      {
        Point* pt = *c;
        vcpoints_ids_.insert(pt->Id());
      }
    }
  }

  if ( vcpoints_ids_.size() == 0)
    dserror("The size of volumecell points is zero!!");

  return vcpoints_ids_;
}
