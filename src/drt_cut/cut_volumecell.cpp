
#include "cut_volumecell.H"
#include "cut_boundarycell.H"
#include "cut_integrationcell.H"
#include "cut_facet.H"
#include "cut_tetmesh.H"
#include "cut_mesh.H"
#include "cut_options.H"
#include "cut_kernel.H"
#include "volume_integration.H"

#include<algorithm>
#include "../../src/drt_fem_general/drt_utils_gausspoints.H"


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
    for ( std::size_t i=0; i<length; ++i )
    {
      std::copy( &points[3*i], &points[3*( i+1 )], &p[0] );
      //Tri3BoundaryCell::CreateCell( mesh, this, f, p );
      NewTri3Cell( mesh, f, p );
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
         << position_ << " "
         << facets_.size() << " "
         << integrationcells_.size() << " "
         << bcells_.size()
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

double GEO::CUT::VolumeCell::Volume()
{
  double volume = 0;
  for ( plain_integrationcell_set::iterator i=integrationcells_.begin(); i!=integrationcells_.end(); ++i )
  {
    IntegrationCell * ic = *i;
    volume += ic->Volume();
  }
  return volume;
}

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

//write the volumecell details for visualization
//gausspoints are not included
void GEO::CUT::VolumeCell::DumpGmsh(std::ofstream& file)
{
    const plain_facet_set & facete = Facets();

    static int pointno=1,point_begin,point_end,lineno=1,line_begin,line_end,surf=1,surf_begin, surf_end;
    surf_begin = surf;
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
         point_begin = pointno;
         Facet *fe = *i;
         const std::vector<std::vector<double> > corners = fe->CornerPointsLocal(ParentElement(),0);
         for(std::vector<std::vector<double> >::const_iterator k=corners.begin();k!=corners.end();k++)
         {
             const std::vector<double> coords = *k;
             file<<"Point("<<pointno<<")={"<<coords[0]<<","<<coords[1]<<","<<coords[2]<<","<<"1"<<"};"<<std::endl;
             pointno++;
         }
         point_end = pointno;
         line_begin = lineno;
         for(int i=point_begin;i!=point_end;i++)
         {
        	 if(i!=point_end-1)
            	 file<<"Line("<<lineno<<")={"<<i<<","<<i+1<<"};"<<std::endl;
             else
                 file<<"Line("<<lineno<<")={"<<i<<","<<point_begin<<"};"<<std::endl;
             lineno++;
         }
         line_end = lineno;
         file<<"Line Loop("<<surf<<")={";
         for(int i=line_begin;i!=line_end;i++)
		 {
			 file<<i;
			 if(i!=line_end-1)
				 file<<",";
		 }
         file<<"};\n";
         file<<"Plane Surface("<<surf<<")={"<<surf<<"};"<<std::endl;

         surf++;
    }
    surf_end = surf;
    file<<"Surface Loop("<<surf<<")={";
    for(int i=surf_begin;i<surf;i++)
    {
		 file<<i;
		 if(i!=(surf-1))
			 file<<",";
	 }
	 file<<"};\n";
	 file<<"Volume("<<surf<<")={"<<surf<<"};"<<std::endl;
}

//write the boundaries of volumecell and the positions of Gauss points for visualization
//a separate file with "side" prefix is generated for every volumecell as the gausspoint
//distribution can be clearly seen
void GEO::CUT::VolumeCell::DumpGmshGaussPoints(const std::vector<std::vector<double> >&gauspts)
{
    const plain_facet_set & facete = Facets();
    std::string filename="side";
    std::ofstream file;

    int pointno=1,point_begin,point_end,lineno=1;
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
         point_begin = pointno;
         Facet *fe = *i;
//the side id of the first facet is used as the file name for every volumecell
         if(i==facete.begin())
         {
                 static int sideno = 0;
                 sideno++;
                 std::stringstream out;
                 out <<"side"<<sideno<<".pos";
                 filename = out.str();
                 file.open(filename.c_str());
         }
         const std::vector<std::vector<double> > corners = fe->CornerPointsLocal(ParentElement(),0);
         for(std::vector<std::vector<double> >::const_iterator k=corners.begin();k!=corners.end();k++)
         {
             const std::vector<double> coords = *k;
             file<<"Point("<<pointno<<")={"<<coords[0]<<","<<coords[1]<<","<<coords[2]<<","<<"1"<<"};"<<std::endl;
             pointno++;
         }
         point_end = pointno;
         for(int i=point_begin;i!=point_end;i++)
         {
                 if(i!=point_end-1)
                         file<<"Line("<<lineno<<")={"<<i<<","<<i+1<<"};"<<std::endl;
                 else
                         file<<"Line("<<lineno<<")={"<<i<<","<<point_begin<<"};"<<std::endl;
                 lineno++;
         }
    }

    for(unsigned i=0;i<gauspts.size();i++)
    {
             file<<"Point("<<pointno<<")={"<<gauspts[i][0]<<","<<gauspts[i][1]<<","<<gauspts[i][2]<<","<<"1"<<"};"<<std::endl;
             pointno++;
    }
    file.close();
}

Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::VolumeCell::GaussPointsFitting()
//void GEO::CUT::VolumeCell::GaussPointsFitting()
{
    Element *ele1 =ParentElement();
//    const DRT::Element::DiscretizationType distyp = ele1->Shape();
//    std::cout<<"distype"<<distyp<<std::endl;
//    const int nsd = DRT::UTILS::DisTypeToDim<DRT::Element::hex8>::dim;
    const int nen = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    LINALG::Matrix<3,nen> xyze;
    ele1->Coordinates(xyze.A());

    Teuchos::RCP<DRT::UTILS::GaussPoints> gp = DRT::UTILS::GaussPointCache::Instance().Create( DRT::Element::hex8, 
                    5/*DRT::UTILS::CubatureDegree(DRT::Element::hex8) */);
    Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new 
                    DRT::UTILS::CollectedGaussPoints( gp->NumPoints() ) );

    for(unsigned i=0;i<gausPts_.size();i++)
    {
        LINALG::Matrix<3,1> xe,xei;
        xe(0,0) = gausPts_[i][0];
        xe(1,0) = gausPts_[i][1];
        xe(2,0) = gausPts_[i][2];

 /*       ele1->LocalCoordinates( xe, xei );

//transformation from global to element local coordinate system
        LINALG::Matrix<nsd,nen> deriv;
        LINALG::Matrix<nsd,nsd> xjm;
        DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(xei,deriv);

        xjm.MultiplyNT( deriv, xyze );

        double det = xjm.Determinant();
 //       std::cout<<weights_[i]<<"\t"<<weights_[i]/det<<"\n";
        double wei = weights_(i)/det; 
//      std::cout<<xei(0,0)<<"\t"<<xei(1,0)<<"\t"<<xei(2,0)<<"\n";*/
//
        cgp->Append( xe, weights_(i) );
    }
/*    double sum=0.0;
    for(unsigned i=0;i<weights_.size();i++)
            sum += weights_[i];
    std::cout<<"volume"<<sum<<"\n";*/
    return cgp;
}

//generate boundary cells for the volumecell
void GEO::CUT::VolumeCell::GenerateBoundaryCells(Mesh &mesh, const GEO::CUT::Point::PointPosition posi)
{
//    std::cout<<"boundary cell size in this volume = "<<bcells_.size()<<"\n";
	const plain_facet_set & facete = Facets();
	for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
	{
		Facet *fac = *i;
		const Side* parside = fac->ParentSide();//order checking
		const std::vector<Node*> &par_nodes = parside->Nodes();
		std::vector<Point*> parpts(3);
		int ptscount=0;
		for(std::vector<Node*>::const_iterator j=par_nodes.begin();
				j!=par_nodes.end();j++)
		{
			Node *parnode = *j;
			parpts[ptscount] = parnode->point();
			ptscount++;
			if(ptscount==3)
				break;
		}
		double parOri[3],sideOri[3];
		OrientationFacet(parpts,parOri);

		if(fac->OnCutSide())
		{
			std::vector<Point*> corners = fac->CornerPoints();
//            std::cout<<"no of corners = "<<corners.size()<<"\n";
//if no of corners are 3 or 4, strightforward
			if(corners.size()==3)
			{
				OrientationFacet(corners,sideOri);
				//the corner points of the boundary cell must be ordered in anti-clockwise manner.
				//because normal is computed in xfem calculations
				bool rever = ToReverse(posi,parOri,sideOri);
				if(rever)
					std::reverse(corners.begin(),corners.end());
				NewTri3Cell( mesh, fac, corners );
			}
			else if(corners.size()==4)
			{
				OrientationFacet(corners,sideOri);
				bool rever = ToReverse(posi,parOri,sideOri);
				if(rever)
					std::reverse(corners.begin(),corners.end());
				NewQuad4Cell(mesh,fac,corners);
			}
//for more than four corners, the facet is triangulated to generate boundary cells
			else
			{
				if(!fac->IsTriangulated())
					fac->DoTriangulation( mesh, corners );
				const std::vector<std::vector<Point*> > & triangulation = fac->Triangulation();
				for ( std::vector<std::vector<Point*> >::const_iterator j=triangulation.begin();
													j!=triangulation.end(); ++j )
				{
					std::vector<Point*> tri = *j;
					OrientationFacet(tri,sideOri);
					bool rever = ToReverse(posi,parOri,sideOri);
					if(rever)
						std::reverse(tri.begin(),tri.end());
					if(tri.size()==3)
						NewTri3Cell(mesh,fac,tri);
					else if(tri.size()==4)
						NewQuad4Cell(mesh,fac,tri);
					else
						dserror("Triangulation created neither tri3 or quad4");
				}
			}
		}
	}
}

bool GEO::CUT::VolumeCell::ToReverse(const GEO::CUT::Point::PointPosition posi,double* parOri,double* sideOri)
{
	bool rever = false;
	if(posi==-3)
	{
		if(fabs(sideOri[0])>0.0000001 && sideOri[0]*parOri[0]>0.0)
				rever = true;
		else if(fabs(sideOri[1])>0.0000001 && sideOri[1]*parOri[1]>0.0)
				rever = true;
		else if(fabs(sideOri[2])>0.0000001 && sideOri[2]*parOri[2]>0.0)
				rever = true;
		else
			//dserror("Point is given instead of side");
			rever = false;
	}
	if(posi==-2)
	{
		if(fabs(sideOri[0])>0.0000001 && sideOri[0]*parOri[0]<0.0)
				rever = true;
		else if(fabs(sideOri[1])>0.0000001 && sideOri[1]*parOri[1]<0.0)
				rever = true;
		else if(fabs(sideOri[2])>0.0000001 && sideOri[2]*parOri[2]<0.0)
				rever = true;
		else
			//dserror("Point is given instead of side");
			rever = false;
	}
	return rever;
}

void GEO::CUT::VolumeCell::OrientationFacet(const std::vector<Point*>pts, double *coef)
{
	int count=0;
	double x1[3]={0.0,0.0,0.0},y1[3]={0.0,0.0,0.0},z1[3]={0.0,0.0,0.0};
	for(std::vector<Point*>::const_iterator i=pts.begin();i!=pts.end();i++)
	{
		Point* pt1 = *i;
		double xm[3];
		pt1->Coordinates(xm);
		x1[count] = xm[0];
		y1[count] = xm[1];
		z1[count] = xm[2];
		count++;
		if(count==3)
			break;
	}

	coef[0] = y1[0]*(z1[1]-z1[2])+y1[1]*(z1[2]-z1[0])+y1[2]*(z1[0]-z1[1]);
	coef[1] = z1[0]*(x1[1]-x1[2])+z1[1]*(x1[2]-x1[0])+z1[2]*(x1[0]-x1[1]);
	coef[2] = x1[0]*(y1[1]-y1[2])+x1[1]*(y1[2]-y1[0])+x1[2]*(y1[0]-y1[1]);
}

//Find Gaussian points and weights by moment fitting equations
void GEO::CUT::VolumeCell::MomentFitGaussWeights(Element *elem, Mesh & mesh, bool include_inner)
{
//        std::cout<<"volume"<<std::endl;
//
//        static int k=0; //blockkk or remove
//if(k>0)        //blockkk or remove
//	return;//blockkk or remove
//const GEO::CUT::Point::PointPosition posu = Position();//blockkk or remove
//if(posu==-3)//blokkk or remove
//	k++;//blockkk or remove*/

	//position is used to decide whether the ordering of points are in clockwise or not
	const GEO::CUT::Point::PointPosition posi = Position();
	std::cout<<"position = "<<posi<<"\n";
//if the volumecell is inside and includeinner is false, no need to compute the Gaussian points
//as this vc will never be computed in xfem algorithm
	if(posi==-2 && include_inner==false)//unblockkk
		return;
    VolumeIntegration vc_inte(this,elem,posi,84); //change the number of equations

    weights_ = vc_inte.compute_weights();
    gausPts_ = vc_inte.getGaussPointLocation(); 

//generate boundary cells. if Tessellation option is used instead of MomentFitting,
//this happens inside "createintegrationcells"
    GenerateBoundaryCells(mesh,posi);

    /*static int kk=1;
    std::ofstream file; //blockkk this group
    std::stringstream out;
    out <<"volume"<<kk<<".pos";
    std::string filename = out.str();
    file.open(filename.c_str());
    DumpGmsh(file);
    file.close();
    kk++;*/

//    k++;               //blockkk or remove
//}
    
    /*std::cout<<"volume"<<"\n";
    const plain_facet_set & facete = Facets();
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
        Facet *fac = *i;
        const std::vector<vector<double> > corners = fac->CornerPointsLocal(elem,0);
        std::cout<<"facet"<<std::endl;
        for(std::vector<vector<double> >::const_iterator k=corners.begin();k!=corners.end();k++)
        {
            vector<double> coords = *k;
            std::cout<<coords[0]<<"\t"<<coords[1]<<"\t"<<coords[2]<<std::endl;
        }
    }*/

//}                      //blockkk or remove 


/*    std::cout<<"volume"<<std::endl;
    for(unsigned i=0;i<weights_.size();i++)
    {
            std::cout<<gausPts_[i][0]<<"\t"<<gausPts_[i][1]<<"\t"<<gausPts_[i][2]<<"\t"<<weights_[i]<<std::endl;
    }*/

/*  const plain_facet_set & facete = Facets();
    for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
        Facet *fe = *i;
        FacetIntegration faee1(fe,elem,posi);
        double dee = faee1.integrate_facet();
        std::cout<<dee<<std::endl;
    //  break;
    }*/


	  /*static int sideno = 0;
          sideno++;
	  std::string filename="wrong";
   	  std::ofstream file;

          std::stringstream out;
          out <<"parent"<<sideno<<".dat";
          filename = out.str();
          file.open(filename.c_str());
	  for (unsigned i=0;i<gausPts_.size();i++)
	  { 
		  file<<gausPts_[i][0]<<"\t"<<gausPts_[i][1]<<"\t"<<gausPts_[i][2]<<"\t"<<weights_[i]<<std::endl;
	  }
	  file.close();*/

}
