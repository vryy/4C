/*!----------------------------------------------------------------------
\file cut_elementhandle.cpp
\brief Outside world interface to element. Converts quadratic to linear element. This provides the
  Gaussian rules generated from the cut

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>

*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "cut_integrationcell.H"
#include "cut_meshintersection.H"
#include "cut_position.H"
#include "cut_point.H"
#include "cut_volumecell.H"

#include "../drt_inpar/inpar_xfem.H"

#include<fstream>


/*----------------------------------------------------------------------*/
// Project the integration rule available in the local coordinates of the
// integation-cells to the local coordinates of background element
/*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::ElementHandle::CreateProjected(
    const std::vector<GEO::CUT::Point*> & cpoints,
    Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic
)
{
  const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<3, nen> xie;
  if ( cpoints.size() != nen )
    throw std::runtime_error( "non-matching number of points" );

  // Find the local coordinates of given corner points w.r to background ElementHandle
  for ( unsigned i=0; i<nen; ++i )
  {
    GEO::CUT::Point * p = cpoints[i];
    const LINALG::Matrix<3,1> & xi = LocalCoordinates( p );
    std::copy( xi.A(), xi.A()+3, &xie( 0, i ) );
  }

  DRT::UTILS::GaussIntegration intpoints( gp_ic );
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp = Teuchos::rcp( new DRT::UTILS::CollectedGaussPoints( gp_ic->NumPoints() ) );

  // Perform actual mapping to correct local coordinates
  DRT::UTILS::GaussIntegration::ProjectGaussPoints<distype> ( xie, intpoints, cgp );
  return cgp;
}



/*----------------------------------------------------------------------*/
// Collect the Gaussian points of all volume-cells belonging to this element in such a way
// that Gaussian rule for every volume-cell can be separated
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::VolumeCellGaussPoints(
    plain_volumecell_set &                      cells,
    std::vector<DRT::UTILS::GaussIntegration> & intpoints,
    INPAR::CUT::VCellGaussPts                   gausstype
)
{
  /*cells.clear(); //check whether any problems with gmsh output
  GetVolumeCells( cells );*/

  intpoints.clear();
  intpoints.reserve( cells.size() );

  //---------------
  // For tessellation, we have Gauss points calculated at local coordinates of each integrationcells
  // we transform this to local coordinates of background ElementHandle
  //----------------
  if(gausstype == INPAR::CUT::VCellGaussPts_Tessellation)
  {
    for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::VolumeCell * vc = *i;

      Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
          Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

      const plain_integrationcell_set & cells = vc->IntegrationCells();
      for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::IntegrationCell * ic = *i;

        Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic = DRT::UTILS::GaussPointCache::Instance().
                                                                    Create( ic->Shape(), ic->CubatureDegree( ic->Shape() ) );
        const std::vector<GEO::CUT::Point*> & cpoints = ic->Points();

        switch ( ic->Shape() )
        {
        case DRT::Element::hex8:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        case DRT::Element::tet4:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        case DRT::Element::wedge6:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::wedge6>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        case DRT::Element::pyramid5:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::pyramid5>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        default:
          throw std::runtime_error( "unsupported integration cell type" );
        }
      }

      intpoints.push_back( DRT::UTILS::GaussIntegration( gpc ) );
    }
  }

  //-------------------
  // For MomentFitting, we have Gauss points that are calculated w.r to local coordinates of linear shadow element
  // If background ElementHandle is linear, then no need for any mapping
  // Else, we map these points to local coordinates of corresponding Quad element
  //-------------------
  else if( gausstype == INPAR::CUT::VCellGaussPts_MomentFitting)
  {
    for(plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i)
    {
      GEO::CUT::VolumeCell *vc = *i;

      Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic = vc->GetGaussRule();
      Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
                      Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

      const std::vector<GEO::CUT::Point*> & cpoints = vc->ParentElement()->Points();

      switch( Shape() )
      {
      case DRT::Element::hex8:
      case DRT::Element::tet4:
      case DRT::Element::wedge6:
      case DRT::Element::pyramid5:
      {
        gpc->Append(gp_ic);
        break;
      }

      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>( cpoints, gp_ic );
        gpc->Append( gp );
        break;
      }
      case DRT::Element::tet10:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>( cpoints, gp_ic );
        gpc->Append( gp );
        break;
      }
      default:
      {
        dserror("element handle for this element is not available\n");
        break;
      }
      }

      intpoints.push_back( DRT::UTILS::GaussIntegration( gpc ) );
    }
  }

  //-------------------
  // For DirectDivergence, we calculate Gauss points at the correct local coord. during construction itself
  // This method is handled separately because
  // 1. main Gauss pts should be mapped w.r to each facet of vcell
  //         --> element volume mapping as done for tessellation and moment fitting do not work
  // 2. Internal Gauss pts can be obtained only if we have correctly mapped main Gauss points
  //-------------------
  else if( gausstype==INPAR::CUT::VCellGaussPts_DirectDivergence )
  {
    for(plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i)
    {
      GEO::CUT::VolumeCell *vc = *i;
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp = vc->GetGaussRule();
      Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
               Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );
      gpc->Append(gp);
      intpoints.push_back( DRT::UTILS::GaussIntegration( gpc ) );
    }
  }

  /*if(IsCut())
  {
    static int eeno=0;
    eeno++;
    if(1)//eeno==1 || eeno==2 || eeno==3)
    {
  for(std::vector<DRT::UTILS::GaussIntegration>::iterator i=intpoints.begin();i!=intpoints.end();i++)
  {
    DRT::UTILS::GaussIntegration ga = *i;
    static int sideno = 0;
          sideno++;
    std::string filename="wrong";
      std::ofstream file;

          std::stringstream out;
          out <<"gauss"<<sideno<<".dat";
          filename = out.str();
          file.open(filename.c_str());
    for ( DRT::UTILS::GaussIntegration::const_iterator iquad=ga.begin(); iquad!=ga.end(); ++iquad )
    {
      const double* gpp = iquad.Point();
      file<<gpp[0]<<"\t"<<gpp[1]<<"\t"<<gpp[2]<<"\t";
      file<<iquad.Weight()<<std::endl;
    }
    file.close();
  }
  }
  }*/

#if 0

  // assume normal integration of element
  DRT::UTILS::GaussIntegration element_intpoints( Shape() );

  // see if we benefit from a "negative" integartion
  int numgps = element_intpoints.NumPoints();
  for ( std::vector<DRT::UTILS::GaussIntegration>::iterator i=intpoints.begin(); i!=intpoints.end(); ++i )
  {
    const DRT::UTILS::GaussIntegration & gi = *i;
    numgps += gi.NumPoints();
  }

  for ( std::vector<DRT::UTILS::GaussIntegration>::iterator i=intpoints.begin(); i!=intpoints.end(); ++i )
  {
    DRT::UTILS::GaussIntegration & gi = *i;
    if ( gi.NumPoints() > numgps / 2 )
    {
      // Here we have the expensive bastard.
      // Exchange gauss points. Use all other points instead.

      Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
        Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( intpoints.size() ) );

      gpc->Append( element_intpoints.Points() );

      for ( std::vector<DRT::UTILS::GaussIntegration>::iterator j=intpoints.begin(); j!=intpoints.end(); ++j )
      {
        const DRT::UTILS::GaussIntegration & gi = *j;
        if ( i!=j )
        {
          gpc->Append( gi.Points() );
        }
      }

      gi.SetPoints( gpc );

      break;
    }
  }
#endif

//   if ( not include_inner )
//   {
//     int count = 0;
//     for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
//     {
//       GEO::CUT::VolumeCell * vc = *i;
//       if ( vc->Position()==Point::inside )
//       {
//         intpoints[count].Clear();
//       }
//       count += 1;
//     }
//   }
}


/*----------------------------------------------------------------------*/
// Collect the Gaussian points of all the volume-cells belonging to this element.
// The integration rules over all the volume-cells are connected.
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPointsComposite> GEO::CUT::ElementHandle::GaussPointsConnected(
    plain_volumecell_set & cells,
    INPAR::CUT::VCellGaussPts gausstype )
{

  Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
      Teuchos::rcp( new DRT::UTILS::GaussPointsComposite( 0 ) );

  for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
  {
    GEO::CUT::VolumeCell * vc = *i;

    const plain_integrationcell_set & cells = vc->IntegrationCells();


    if(gausstype == INPAR::CUT::VCellGaussPts_Tessellation)
    {
      for ( plain_integrationcell_set::const_iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::IntegrationCell * ic = *i;

        Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic = DRT::UTILS::GaussPointCache::Instance().
            Create( ic->Shape(), ic->CubatureDegree( ic->Shape() ) );
        const std::vector<GEO::CUT::Point*> & cpoints = ic->Points();

        switch ( ic->Shape() )
        {
        case DRT::Element::hex8:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::hex8>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        case DRT::Element::tet4:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::tet4>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        case DRT::Element::wedge6:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::wedge6>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        case DRT::Element::pyramid5:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = CreateProjected<DRT::Element::pyramid5>( cpoints, gp_ic );
          gpc->Append( gp );
          break;
        }
        default:
          throw std::runtime_error( "unsupported integration cell type" );
        }
      }
    }
    else if( gausstype == INPAR::CUT::VCellGaussPts_MomentFitting || gausstype==INPAR::CUT::VCellGaussPts_DirectDivergence )
    {
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp = vc->GetGaussRule();
      gpc->Append(gp);
    }
  }

  return gpc;
}

#if(0)
/*----------------------------------------------------------------------*/
// Unused
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::BoundaryCellGaussPoints(
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
    std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints
)
{
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::tri3:
      {
        cell_points.push_back( DRT::UTILS::GaussIntegration( DRT::Element::tri3,
                                                             bc->CubatureDegree( DRT::Element::quad9 ) ) );
        break;
      }
      case DRT::Element::quad4:
      {
        cell_points.push_back( DRT::UTILS::GaussIntegration( DRT::Element::quad4,
                                                             bc->CubatureDegree( DRT::Element::quad9 ) ) );
        break;
      }
      default:
        throw std::runtime_error( "unsupported integration cell type" );
      }
    }
  }
}


/*----------------------------------------------------------------------*/
// unused: Old implementation of boundarycell Gauss points collection
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::BoundaryCellGaussPoints(
    MeshIntersection & mesh,
    int mi,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
    std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints
)
{
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    SideHandle * side = mesh.GetCutSide( sid, mi );
    if ( side==NULL )
    {
      throw std::runtime_error( "no such side" );
    }

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::tri3:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::tri3>( bc );
        cell_points.push_back( DRT::UTILS::GaussIntegration( gp ) );
        break;
      }
      case DRT::Element::quad4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::quad4>( bc );
        cell_points.push_back( DRT::UTILS::GaussIntegration( gp ) );
        break;
      }
      default:
        throw std::runtime_error( "unsupported integration cell type" );
      }
    }
  }
}
#endif

/*----------------------------------------------------------------------*/
// Collect the Gauss points of all the boundary-cells belong to this element.
// This is the method used now in the new implementation
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::BoundaryCellGaussPointsLin(
    MeshIntersection & mesh,
    int mi,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
    std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints
    )
{
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    SideHandle * side = mesh.GetCutSide( sid, mi );
    if ( side==NULL )
    {
      throw std::runtime_error( "no such side" );
    }

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      // Create (unmodified) gauss points for integration cell with requested
      // polynomial order. This is supposed to be fast, since there is a cache.
      cell_points.push_back(bc->gaussRule());
    }
  }
}

/*----------------------------------------------------------------------*/
// Collect the Gauss points of all the boundary-cells belong to this element.
// Used for Level set applications as cutside does not exist here.
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::BoundaryCellGaussPointsLin( const std::map<int, std::vector<GEO::CUT::BoundaryCell*> > & bcells,
                                                               std::map<int, std::vector<DRT::UTILS::GaussIntegration> > & intpoints )
{

  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::const_iterator i=bcells.begin(); i!=bcells.end(); ++i )
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*> & cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration> & cell_points = intpoints[sid];

    cell_points.clear();
    cell_points.reserve( cells.size() );

    for ( std::vector<GEO::CUT::BoundaryCell*>::const_iterator i=cells.begin(); i!=cells.end(); ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      // Create (unmodified) gauss points for integration cell with requested
      // polynomial order. This is supposed to be fast, since there is a cache.
      cell_points.push_back(bc->gaussRule());
    }
  }
}



/*----------------------------------------------------------------------*/
// get all the element' sets of volume-cells and nds-vectors
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::GetVolumeCellsDofSets (
    std::vector<plain_volumecell_set > & cellsets ,
    std::vector< std::vector< int > >  & nds_sets,
    bool include_inner
)
{

  const std::vector<plain_volumecell_set> & ele_vc_sets_inside = GetVcSetsInside();
  const std::vector<plain_volumecell_set> & ele_vc_sets_outside = GetVcSetsOutside();

  std::vector<std::vector<int> > & nodaldofset_vc_sets_inside = GetNodalDofSet_VcSets_Inside();
  std::vector<std::vector<int> > & nodaldofset_vc_sets_outside = GetNodalDofSet_VcSets_Outside();

  if(include_inner)
  {
    std::copy( ele_vc_sets_inside.begin(), ele_vc_sets_inside.end(), std::inserter(cellsets,cellsets.end()) );
    std::copy( nodaldofset_vc_sets_inside.begin(), nodaldofset_vc_sets_inside.end(), std::inserter(nds_sets,nds_sets.end()) );
  }

  std::copy( ele_vc_sets_outside.begin(), ele_vc_sets_outside.end(), std::inserter(cellsets,cellsets.end()) );
  std::copy( nodaldofset_vc_sets_outside.begin(), nodaldofset_vc_sets_outside.end(), std::inserter(nds_sets,nds_sets.end()) );

}





/*----------------------------------------------------------------------*/
//! Collect all volume-cells belonging to this elements
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::GetVolumeCells( plain_volumecell_set & cells )
{
  const plain_volumecell_set & cs = element_->VolumeCells();
  std::copy( cs.begin(), cs.end(), std::inserter( cells, cells.begin() ) );
}


/*----------------------------------------------------------------------*/
//! Collect all volume-cells belonging to this element ordered by position
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::CollectVolumeCells (
    bool include_inner,
    plain_volumecell_set & cells_inside,
    plain_volumecell_set & cells_outside
)
{

  const plain_volumecell_set & ecells = element_->VolumeCells();

  // sort for inside and outside volume cells
  for (plain_volumecell_set::const_iterator i = ecells.begin(); i!=ecells.end(); i++)
  {
    if( (*i)->Position() == GEO::CUT::Point::outside )
    {
      cells_outside.insert( *i );
    }
    else // inside vc
    {
      if( include_inner)
      {
        cells_inside.insert( *i );
      }
    }
  }

}


/*----------------------------------------------------------------------*/
// get all the element' sets of volume-cells, nds-vectors and integration points
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::GetCellSets_DofSets_GaussPoints (
    std::vector< plain_volumecell_set >                        & cell_sets ,
    std::vector< std::vector< int > >                          & nds_sets,
    std::vector< std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
    INPAR::CUT::VCellGaussPts                                    gausstype,
    bool                                                         include_inner
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Get_CellSets_nds_GaussPoints" );

  GetVolumeCellsDofSets( cell_sets, nds_sets, include_inner);

  intpoints_sets.clear();

  for(std::vector<plain_volumecell_set>::iterator i= cell_sets.begin(); i!=cell_sets.end(); i++)
  {
    plain_volumecell_set & cells = *i;

    if(cells.size() != 1) dserror("number of volumecells in set not equal 1, this should not be for linear elements!");

    std::vector<DRT::UTILS::GaussIntegration> gaussCellsets;
    VolumeCellGaussPoints( cells, gaussCellsets, gausstype );

    intpoints_sets.push_back( gaussCellsets );
  }
}


/*----------------------------------------------------------------------*/
// get the element's sets of volume-cells ordered by inside/outside position
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::VolumeCellSets (
    bool include_inner,
    std::vector<plain_volumecell_set> & ele_vc_sets_inside,
    std::vector<plain_volumecell_set> & ele_vc_sets_outside
)
{

  if(!cells_set_)
  {
    const plain_volumecell_set & ecells = element_->VolumeCells();

    // sort for inside and outside volume cells
    for (plain_volumecell_set::const_iterator i = ecells.begin(); i!=ecells.end(); i++)
    {
      if( (*i)->Position() == GEO::CUT::Point::outside )
      {
        plain_volumecell_set s;
        s.insert(*i);
        vc_sets_outside_.push_back(s);
      }
      else // inside vc
      {
        if( include_inner)
        {
          plain_volumecell_set s;
          s.insert(*i);
          vc_sets_inside_.push_back(s);
        }
      }
    }

    cells_set_ = true;
  }



  if(include_inner)
  {
    std::copy( vc_sets_inside_.begin(), vc_sets_inside_.end(), std::inserter(ele_vc_sets_inside, ele_vc_sets_inside.begin()));
  }

  std::copy( vc_sets_outside_.begin(), vc_sets_outside_.end(), std::inserter(ele_vc_sets_outside, ele_vc_sets_outside.begin()));

}



/*----------------------------------------------------------------------*/
// returns true in case that any cut-side cut with the element produces cut points,
// i.e. also for touched cases (at points, edges or sides),
// or when an element side has more than one facet or is touched by fully/partially by the cut side
/*----------------------------------------------------------------------*/
bool GEO::CUT::QuadraticElementHandle::IsCut()
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    if ( e->IsCut() )
    {
      return true;
    }
  }
  return false;
}


/*----------------------------------------------------------------------*/
// return true if one of the sub-elements is intersected or the sub-elements
// have different positions and therefore the global element is intersected by a cut-side
/*----------------------------------------------------------------------*/
bool GEO::CUT::QuadraticElementHandle::IsIntersected()
{
  GEO::CUT::Point::PointPosition unique_pos = GEO::CUT::Point::undecided;

  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    if ( e->IsIntersected() )
    {
      return true;
    }

    // do the sub-elements have different positions ? (e.g. when the cut side directly cuts between sub-elements)
    if(     (unique_pos != GEO::CUT::Point::undecided)
        and (e->VolumeCells()[0]->Position() != unique_pos) )
      return true;
  }
  return false;
}


/*----------------------------------------------------------------------*/
// Collect all volume-cells belonging to this elements
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetVolumeCells( plain_volumecell_set & cells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    const plain_volumecell_set & cs = e->VolumeCells();
    std::copy( cs.begin(), cs.end(), std::inserter( cells, cells.begin() ) );
  }
}


/*----------------------------------------------------------------------*/
// Collect all volume-cells belonging to this element ordered by position
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::CollectVolumeCells (
    bool include_inner,
    plain_volumecell_set & cells_inside,
    plain_volumecell_set & cells_outside
)
{
  for ( std::vector<Element*>::const_iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    const plain_volumecell_set & ecells = e->VolumeCells();

    // sort for inside and outside volume-cells
    for (plain_volumecell_set::const_iterator i = ecells.begin(); i!=ecells.end(); i++)
    {
      if( (*i)->Position() == GEO::CUT::Point::outside )
      {
        cells_outside.insert( *i );
      }
      else // inside vc
      {
        if( include_inner)
        {
          cells_inside.insert( *i );
        }
      }
    } // volume-cells
  } // sub-elements
}


/*----------------------------------------------------------------------*/
//  get the quadratic element's volumetric integration cells (just for Tessellation)
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetIntegrationCells( plain_integrationcell_set & cells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->GetIntegrationCells( cells );
  }
}


/*----------------------------------------------------------------------*/
//  get all the quadratic element's boundary integration cells
// TODO: this has to be corrected such that just the bcs which belong the outside vcs will be returned
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetBoundaryCells( plain_boundarycell_set & bcells )
{
  for ( std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * e = *i;
    e->GetBoundaryCells( bcells );
  }
}


/*----------------------------------------------------------------------*/
// get all the element' sets of volume-cells, nds-vectors and integration points
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetCellSets_DofSets_GaussPoints (
    std::vector< plain_volumecell_set >                        & cell_sets ,
    std::vector< std::vector< int > >                          & nds_sets,
    std::vector< std::vector< DRT::UTILS::GaussIntegration > > & intpoints_sets,
    INPAR::CUT::VCellGaussPts                                    gausstype,
    bool                                                         include_inner
    )
{
  GetVolumeCellsDofSets( cell_sets, nds_sets, include_inner );

  intpoints_sets.clear();
  intpoints_sets.reserve( cell_sets.size() );

  for(std::vector<plain_volumecell_set>::iterator i= cell_sets.begin(); i!=cell_sets.end(); i++)
  {
    plain_volumecell_set & cells = *i;

    std::vector<DRT::UTILS::GaussIntegration> gaussCellsets;
    VolumeCellGaussPoints( cells, gaussCellsets, gausstype );

    intpoints_sets.push_back( gaussCellsets );
  }

}


/*----------------------------------------------------------------------*/
//! get the element's sets of volume-cells ordered by inside/outside position
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::VolumeCellSets (
    bool include_inner,
    std::vector<plain_volumecell_set> & ele_vc_sets_inside,
    std::vector<plain_volumecell_set> & ele_vc_set_ouside
)
{

  // connect volumecells of subelements
  ConnectVolumeCells( include_inner );

  // get inside and outside connected volumecellsets for the current element
  GetConnectedVolumeCellSets( include_inner, ele_vc_sets_inside, ele_vc_set_ouside );

}


/*----------------------------------------------------------------------*/
//! connect volume-cells to sets of volume-cells
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::ConnectVolumeCells ( bool include_inner )
{
  // find the connection between volumecells of all subelements for the current element (hex8, hex20, tet4, tet10 etc.)
  // remark: this function determines not the connection outside this element

  // get the volumecells of all subelements stored in two plainvolume_sets (inside and outside)
  if( !cells_connected_ )
  {

    plain_volumecell_set e_vcs_inside;
    plain_volumecell_set e_vcs_outside;

    CollectVolumeCells(include_inner, e_vcs_inside, e_vcs_outside);

    if( include_inner )
    {
      BuildCellSets(e_vcs_inside , connected_vc_sets_inside_);
    }
    BuildCellSets(e_vcs_outside , connected_vc_sets_outside_);

    cells_connected_ = true;

  }

}


/*----------------------------------------------------------------------*/
//! get connections/sets of volume-cells between sub-elements
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetConnectedVolumeCellSets ( bool include_inner, std::vector<plain_volumecell_set> & connected_vc_sets_inside, std::vector<plain_volumecell_set> & connected_vc_sets_outside)
{

  if(include_inner)
  {
    std::copy(connected_vc_sets_inside_.begin(), connected_vc_sets_inside_.end(), std::inserter(connected_vc_sets_inside, connected_vc_sets_inside.end()));
  }

  std::copy(connected_vc_sets_outside_.begin(), connected_vc_sets_outside_.end(), std::inserter(connected_vc_sets_outside, connected_vc_sets_outside.end()));
}


/*----------------------------------------------------------------------*/
//! build sets
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::BuildCellSets (
    plain_volumecell_set & cells_to_connect,
    std::vector<plain_volumecell_set> & connected_sets
)
{
  plain_volumecell_set done;

  for ( plain_volumecell_set::const_iterator i=cells_to_connect.begin();
      i!=cells_to_connect.end();
      ++i )
  {
    VolumeCell * cell = *i;
    if ( done.count( cell )==0 ) // cell currently not-done
    {
      plain_volumecell_set connected;
      // REMARK: here use the version without! elements check:
      // here we build cell sets within one global element with vcs of subelements
      // maybe the vcs of one subelement are not connected within one subelement, but within one global element,
      // therefore more than one vc of one subelements may be connected.
      //        cell->Neighbors( NULL, cells_to_connect, done, connected, elements );
      cell->Neighbors( NULL, cells_to_connect, done, connected);

      if ( connected.size()>0 )
      {
        connected_sets.push_back( connected );
        std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );
      }
    }
  }
}


/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Hex20ElementHandle::Hex20ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nodes )
: QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  // create middle nodes

  LINALG::Matrix<3,1> xyz;
  LINALG::Matrix<1,1> lsv;

  plain_int_set node_nids;

  LINALG::Matrix<3,8> side_xyze;
  LINALG::Matrix<1,8> side_lsvs;
  LINALG::Matrix<8,1> side_funct;
  std::vector<Node*> side_nodes( 8 );

  std::vector<Node*> center_nodes( 6 );

  // get the real nodes of the hex20 element and get the center nodes (shadow nodes) of its quadratic sides
  // which stored in the shadow_nodes list in cut_mesh.H using the eight side nodes of the side as key
  // special handling for the inner center node of hex20 element required, see below,
  // the inner center node is stored in the shadow_nodes map using the 20 nodes of the hex20 element as key

  // loop the six sides of the hex 20 element
  for ( int localsideid = 0; localsideid < 6; ++localsideid )
  {
    node_nids.clear();
    // loop the eight nodes of each quad8 side of the hex20 element
    for ( int i=0; i<8; ++i )
    {
      int localnodeid = DRT::UTILS::eleNodeNumbering_hex27_surfaces[localsideid][i];
      Node * n = mesh.GetNode( nodes[localnodeid], static_cast<double*>( NULL ) );
      side_nodes[i] = n;
      node_nids.insert( nodes[localnodeid] );
      n->Coordinates( &side_xyze( 0, i ) );
      side_lsvs( i ) = n->LSV();
    }

    DRT::UTILS::shape_function_2D( side_funct, 0.0, 0.0, DRT::Element::quad8 );
    xyz.Multiply( side_xyze, side_funct );
    lsv.Multiply( side_lsvs, side_funct );

    // find the unique center node of the quadratic quad8 sides
    center_nodes[localsideid] = mesh.GetNode( node_nids, xyz.A(), lsv( 0 ) );
  }

  Node* node20 = center_nodes[0];
  int node20_id = node20->Id();

  Node* node21 = center_nodes[1];
  int node21_id = node21->Id();

  Node* node22 = center_nodes[2];
  int node22_id = node22->Id();

  Node* node23 = center_nodes[3];
  int node23_id = node23->Id();

  Node* node24 = center_nodes[4];
  int node24_id = node24->Id();

  Node* node25 = center_nodes[5];
  int node25_id = node25->Id();

  LINALG::Matrix<3,20> xyze;
  LINALG::Matrix<1,20> lsvs;
  nodes_.reserve( 20 );
  for ( int i=0; i<20; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], static_cast<double*>( NULL ) );
    nodes_.push_back( n );
    n->Coordinates( &xyze( 0,i ) );
    lsvs( i ) = n->LSV();
  }

  // special handling for the inner center node of hex20 element
  // Remark: this node is also stored as shadow node in cut_mesh, however, the key for this shadow node
  // is the set of all 20 nodes of the hex20 element
  // in contrast to the shadow nodes of sides, for that the key are the eight nodes of the quad7 side
  LINALG::Matrix<20,1> funct;
  DRT::UTILS::shape_function_3D( funct, 0.0, 0.0, 0.0, DRT::Element::hex20 );

  xyz.Multiply( xyze, funct );
  lsv.Multiply( lsvs, funct );
  node_nids.clear();
  std::copy( nodes.begin(), nodes.end(), std::inserter( node_nids, node_nids.begin() ) );
  Node* node26 = mesh.GetNode( node_nids, xyz.A(), lsv( 0 ) );
  int node26_id = node26->Id();


  std::vector<int> nids( 8 );

  nids[0] = nodes[ 0];
  nids[1] = nodes[ 8];
  nids[2] = node20_id;
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = node21_id;
  nids[6] = node26_id;
  nids[7] = node24_id;
  Element* sub1 = mesh.GetElement( -1, nids, *top_data );
  sub1->setAsShadowElem();
  sub1->setQuadCorners( mesh, nodes );
  sub1->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub1 );

  nids[0] = nodes[ 8];
  nids[1] = nodes[ 1];
  nids[2] = nodes[ 9];
  nids[3] = node20_id;
  nids[4] = node21_id;
  nids[5] = nodes[13];
  nids[6] = node22_id;
  nids[7] = node26_id;
  Element* sub2 = mesh.GetElement( -1, nids, *top_data );
  sub2->setAsShadowElem();
  sub2->setQuadCorners( mesh, nodes );
  sub2->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub2 );

  nids[0] = node20_id;
  nids[1] = nodes[ 9];
  nids[2] = nodes[ 2];
  nids[3] = nodes[10];
  nids[4] = node26_id;
  nids[5] = node22_id;
  nids[6] = nodes[14];
  nids[7] = node23_id;
  Element* sub3 = mesh.GetElement( -1, nids, *top_data );
  sub3->setAsShadowElem();
  sub3->setQuadCorners( mesh, nodes );
  sub3->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub3 );

  nids[0] = nodes[11];
  nids[1] = node20_id;
  nids[2] = nodes[10];
  nids[3] = nodes[ 3];
  nids[4] = node24_id;
  nids[5] = node26_id;
  nids[6] = node23_id;
  nids[7] = nodes[15];
  Element* sub4 = mesh.GetElement( -1, nids, *top_data );
  sub4->setAsShadowElem();
  sub4->setQuadCorners( mesh, nodes );
  sub4->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub4 );

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = node21_id;
  nids[2] = node26_id;
  nids[3] = node24_id;
  nids[4] = nodes[ 4];
  nids[5] = nodes[16];
  nids[6] = node25_id;
  nids[7] = nodes[19];
  Element* sub5 = mesh.GetElement( -1, nids, *top_data );
  sub5->setAsShadowElem();
  sub5->setQuadCorners( mesh, nodes );
  sub5->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub5 );

  nids[0] = node21_id;
  nids[1] = nodes[13];
  nids[2] = node22_id;
  nids[3] = node26_id;
  nids[4] = nodes[16];
  nids[5] = nodes[ 5];
  nids[6] = nodes[17];
  nids[7] = node25_id;
  Element* sub6 = mesh.GetElement( -1, nids, *top_data );
  sub6->setAsShadowElem();
  sub6->setQuadCorners( mesh, nodes );
  sub6->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub6 );

  nids[0] = node26_id;
  nids[1] = node22_id;
  nids[2] = nodes[14];
  nids[3] = node23_id;
  nids[4] = node25_id;
  nids[5] = nodes[17];
  nids[6] = nodes[ 6];
  nids[7] = nodes[18];
  Element* sub7 = mesh.GetElement( -1, nids, *top_data );
  sub7->setAsShadowElem();
  sub7->setQuadCorners( mesh, nodes );
  sub7->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub7 );

  nids[0] = node24_id;
  nids[1] = node26_id;
  nids[2] = node23_id;
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = node25_id;
  nids[6] = nodes[18];
  nids[7] = nodes[ 7];
  Element* sub8 = mesh.GetElement( -1, nids, *top_data );
  sub8->setAsShadowElem();
  sub8->setQuadCorners( mesh, nodes );
  sub8->setQuadShape( DRT::Element::hex20 );
  subelements_.push_back( sub8 );

  // each subelement should know its parents id
  for(std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * subelement = *i;
    subelement->ParentId(eid);
  }

}


/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Hex27ElementHandle::Hex27ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nodes )
: QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  nodes_.reserve( 27 );
  for ( int i=0; i<27; ++i )
  {
    Node * n = mesh.GetNode( nodes[i], static_cast<double*>( NULL ) );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Hexahedron<8> >();

  std::vector<int> nids( 8 );

  nids[0] = nodes[ 0];
  nids[1] = nodes[ 8];
  nids[2] = nodes[20];
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = nodes[21];
  nids[6] = nodes[26];
  nids[7] = nodes[24];
  Element* sub1 = mesh.GetElement( -1, nids, *top_data );
  sub1->setAsShadowElem();
  sub1->setQuadCorners( mesh, nodes );
  sub1->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub1 );

  nids[0] = nodes[ 8];
  nids[1] = nodes[ 1];
  nids[2] = nodes[ 9];
  nids[3] = nodes[20];
  nids[4] = nodes[21];
  nids[5] = nodes[13];
  nids[6] = nodes[22];
  nids[7] = nodes[26];
  Element* sub2 = mesh.GetElement( -1, nids, *top_data );
  sub2->setAsShadowElem();
  sub2->setQuadCorners( mesh, nodes );
  sub2->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub2 );

  nids[0] = nodes[20];
  nids[1] = nodes[ 9];
  nids[2] = nodes[ 2];
  nids[3] = nodes[10];
  nids[4] = nodes[26];
  nids[5] = nodes[22];
  nids[6] = nodes[14];
  nids[7] = nodes[23];
  Element* sub3 = mesh.GetElement( -1, nids, *top_data );
  sub3->setAsShadowElem();
  sub3->setQuadCorners( mesh, nodes );
  sub3->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub3 );

  nids[0] = nodes[11];
  nids[1] = nodes[20];
  nids[2] = nodes[10];
  nids[3] = nodes[ 3];
  nids[4] = nodes[24];
  nids[5] = nodes[26];
  nids[6] = nodes[23];
  nids[7] = nodes[15];
  Element* sub4 = mesh.GetElement( -1, nids, *top_data );
  sub4->setAsShadowElem();
  sub4->setQuadCorners( mesh, nodes );
  sub4->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub4 );

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = nodes[21];
  nids[2] = nodes[26];
  nids[3] = nodes[24];
  nids[4] = nodes[ 4];
  nids[5] = nodes[16];
  nids[6] = nodes[25];
  nids[7] = nodes[19];
  Element* sub5 = mesh.GetElement( -1, nids, *top_data );
  sub5->setAsShadowElem();
  sub5->setQuadCorners( mesh, nodes );
  sub5->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub5 );

  nids[0] = nodes[21];
  nids[1] = nodes[13];
  nids[2] = nodes[22];
  nids[3] = nodes[26];
  nids[4] = nodes[16];
  nids[5] = nodes[ 5];
  nids[6] = nodes[17];
  nids[7] = nodes[25];
  Element* sub6 = mesh.GetElement( -1, nids, *top_data );
  sub6->setAsShadowElem();
  sub6->setQuadCorners( mesh, nodes );
  sub6->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub6 );

  nids[0] = nodes[26];
  nids[1] = nodes[22];
  nids[2] = nodes[14];
  nids[3] = nodes[23];
  nids[4] = nodes[25];
  nids[5] = nodes[17];
  nids[6] = nodes[ 6];
  nids[7] = nodes[18];
  Element* sub7 = mesh.GetElement( -1, nids, *top_data );
  sub7->setAsShadowElem();
  sub7->setQuadCorners( mesh, nodes );
  sub7->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub7 );

  nids[0] = nodes[24];
  nids[1] = nodes[26];
  nids[2] = nodes[23];
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = nodes[25];
  nids[6] = nodes[18];
  nids[7] = nodes[ 7];
  Element* sub8 = mesh.GetElement( -1, nids, *top_data );
  sub8->setAsShadowElem();
  sub8->setQuadCorners( mesh, nodes );
  sub8->setQuadShape( DRT::Element::hex27 );
  subelements_.push_back( sub8 );

  // each subelement should know its parents id
  for(std::vector<Element*>::iterator i=subelements_.begin(); i!=subelements_.end(); ++i )
  {
    Element * subelement = *i;
    subelement->ParentId(eid);
  }

}


/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Tet10ElementHandle::Tet10ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids )
  : QuadraticElementHandle()
{
  subelements_.reserve( 8 );

  nodes_.reserve( 10 );
  for ( int i=0; i<10; ++i )
  {
    Node * n = mesh.GetNode( nids[i], static_cast<double*>( NULL ) );
    nodes_.push_back( n );
  }

  const CellTopologyData * top_data = shards::getCellTopologyData< shards::Tetrahedron<4> >();

  std::vector<int> subnids( 4 );

  subnids[0] = nids[ 0];
  subnids[1] = nids[ 4];
  subnids[2] = nids[ 6];
  subnids[3] = nids[ 7];
  Element* sub1 = mesh.GetElement( -1, subnids, *top_data );
  sub1->setAsShadowElem();
  sub1->setQuadCorners( mesh, nids );
  sub1->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub1 );

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 1];
  subnids[2] = nids[ 5];
  subnids[3] = nids[ 8];
  Element* sub2 = mesh.GetElement( -1, subnids, *top_data );
  sub2->setAsShadowElem();
  sub2->setQuadCorners( mesh, nids );
  sub2->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub2 );

  subnids[0] = nids[ 6];
  subnids[1] = nids[ 5];
  subnids[2] = nids[ 2];
  subnids[3] = nids[ 9];
  Element* sub3 = mesh.GetElement( -1, subnids, *top_data );
  sub3->setAsShadowElem();
  sub3->setQuadCorners( mesh, nids );
  sub3->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub3 );

  subnids[0] = nids[ 7];
  subnids[1] = nids[ 8];
  subnids[2] = nids[ 9];
  subnids[3] = nids[ 3];
  Element* sub4 = mesh.GetElement( -1, subnids, *top_data );
  sub4->setAsShadowElem();
  sub4->setQuadCorners( mesh, nids );
  sub4->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub4 );

  /////////////////////////////////////////////////////////////////

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 5];
  subnids[2] = nids[ 6];
  subnids[3] = nids[ 8];
  Element* sub5 = mesh.GetElement( -1, subnids, *top_data );
  sub5->setAsShadowElem();
  sub5->setQuadCorners( mesh, nids );
  sub5->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub5 );

  subnids[0] = nids[ 6];
  subnids[1] = nids[ 9];
  subnids[2] = nids[ 7];
  subnids[3] = nids[ 8];
  Element* sub6 = mesh.GetElement( -1, subnids, *top_data );
  sub6->setAsShadowElem();
  sub6->setQuadCorners( mesh, nids );
  sub6->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub6 );

  subnids[0] = nids[ 4];
  subnids[1] = nids[ 6];
  subnids[2] = nids[ 7];
  subnids[3] = nids[ 8];
  Element* sub7 = mesh.GetElement( -1, subnids, *top_data );
  sub7->setAsShadowElem();
  sub7->setQuadCorners( mesh, nids );
  sub7->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub7 );

  subnids[0] = nids[ 9];
  subnids[1] = nids[ 6];
  subnids[2] = nids[ 5];
  subnids[3] = nids[ 8];
  Element* sub8 = mesh.GetElement( -1, subnids, *top_data );
  sub8->setAsShadowElem();
  sub8->setQuadCorners( mesh, nids );
  sub8->setQuadShape( DRT::Element::tet10 );
  subelements_.push_back( sub8 );
}


/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Hex20ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex20> pos( nodes_, xyz );

  bool success = pos.ComputeTol(1e-10);;
  if ( not success )
  {
    std::cout << "local coordinates for hex20 element could not be determined" << std::endl;
    for(int i=0; i< (int)(nodes_.size()); i++)
    {
      std::cout << " node " << i << std::endl;
      nodes_[i]->Plot(std::cout);
    }

    std::cout << "point in xyz: " << xyz << std::endl;

    dserror("local coordinates for hex20 element could not be determined");
  }
  rst = pos.LocalCoordinates();
}


/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Hex27ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::hex27> pos( nodes_, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}


/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Tet10ElementHandle::LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
{
  Position<DRT::Element::tet10> pos( nodes_, xyz );
  bool success = pos.Compute();
  if ( not success )
  {
  }
  rst = pos.LocalCoordinates();
}

