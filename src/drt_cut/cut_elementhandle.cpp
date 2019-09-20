/*----------------------------------------------------------------------*/
/*! \file
\brief Outside world interface to element. Converts quadratic to linear element. This provides the
  Gaussian rules generated from the cut

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

\level 2
*/
/*----------------------------------------------------------------------*/


#include <Teuchos_TimeMonitor.hpp>

#include "cut_boundarycell.H"
#include "cut_integrationcell.H"
#include "cut_position.H"
#include "cut_point.H"
#include "cut_node.H"
#include "cut_volumecell.H"
#include "cut_mesh.H"
#include "cut_tolerance.H"
#include "quadrature_compression.H"


#include "../drt_inpar/inpar_xfem.H"

#include "cut_elementhandle.H"

#include <fstream>


/*----------------------------------------------------------------------*/
// Project the integration rule available in the local coordinates of the
// integation-cells to the local coordinates of background element
/*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<DRT::UTILS::GaussPoints> GEO::CUT::ElementHandle::CreateProjected(
    const std::vector<GEO::CUT::Point*>& cpoints, Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic)
{
  const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  const unsigned dim = DRT::UTILS::DisTypeToDim<distype>::dim;
  LINALG::Matrix<dim, nen> xie;
  if (cpoints.size() != nen) run_time_error("non-matching number of points");

  // Find the local coordinates of given corner points w.r. to background ElementHandle
  for (unsigned i = 0; i < nen; ++i)
  {
    GEO::CUT::Point* p = cpoints[i];
    const LINALG::Matrix<3, 1>& xi = LocalCoordinates(p);

    // copy first dim entries into xie
    std::copy(xi.A(), xi.A() + dim, &xie(0, i));
  }

  DRT::UTILS::GaussIntegration intpoints(gp_ic);
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> cgp =
      Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints(gp_ic->NumPoints()));

  // Perform actual mapping to correct local coordinates
  DRT::UTILS::GaussIntegration::ProjectGaussPointsLocalToGlobal<distype>(xie, intpoints, cgp);
  return cgp;
}


/*----------------------------------------------------------------------*/
// Collect the Gaussian points of all volume-cells belonging to this element in such a way
// that Gaussian rule for every volume-cell can be separated
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::VolumeCellGaussPoints(
    plain_volumecell_set& cells, std::vector<DRT::UTILS::GaussIntegration>& intpoints)
{
  intpoints.clear();
  intpoints.reserve(cells.size());

  for (plain_volumecell_set::iterator i = cells.begin(); i != cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = *i;

    Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
        Teuchos::rcp(new DRT::UTILS::GaussPointsComposite(0));

    switch (vc->ParentElement()->GetElementIntegrationType())
    {
      case INPAR::CUT::EleIntType_Tessellation:
      {
        AppendVolumeCellGaussPoints_Tessellation(gpc, vc);
        break;
      }
      case INPAR::CUT::EleIntType_MomentFitting:
      {
        AppendVolumeCellGaussPoints_MomentFitting(gpc, vc);
        break;
      }
      case INPAR::CUT::EleIntType_DirectDivergence:
      {
        AppendVolumeCellGaussPoints_DirectDivergence(gpc, vc);
        break;
      }
      default:
      {
        dserror("non supported element integration type for given volume-cell %i",
            vc->ParentElement()->GetElementIntegrationType());
        exit(EXIT_FAILURE);
      }
    }

#ifdef QUADCOMP
    if (gpc->NumPoints() > 56)
    {
      QuadratureCompression qc;
      bool quad_comp_success = qc.PerformCompressionOfQuadrature(*gpc, vc);
      if (quad_comp_success)
      {
        intpoints.push_back(DRT::UTILS::GaussIntegration(qc.GetCompressedQuadrature()));

        // reset the Gauss points for the volumecell so that the compression need not be performed
        // for each iteration within the Newton loop
        vc->SetGaussRule(qc.GetCompressedQuadrature());
      }
      else
        intpoints.push_back(DRT::UTILS::GaussIntegration(gpc));
    }
    else
    {
      intpoints.push_back(DRT::UTILS::GaussIntegration(gpc));
    }
#else
    intpoints.push_back(DRT::UTILS::GaussIntegration(gpc));
#endif
  }
}


void GEO::CUT::ElementHandle::AppendVolumeCellGaussPoints_Tessellation(
    Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc, GEO::CUT::VolumeCell* vc)
{
  //---------------
  // For tessellation, we have Gauss points calculated at local coordinates of each integrationcells
  // we transform this to local coordinates of background ElementHandle
  //----------------
  const GEO::CUT::plain_integrationcell_set& cells = vc->IntegrationCells();
  for (GEO::CUT::plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
  {
    GEO::CUT::IntegrationCell* ic = *i;

    Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic = DRT::UTILS::GaussPointCache::Instance().Create(
        ic->Shape(), ic->CubatureDegree(ic->Shape()));
    const std::vector<GEO::CUT::Point*>& cpoints = ic->Points();

    switch (ic->Shape())
    {
      case DRT::Element::tri3:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            CreateProjected<DRT::Element::tri3>(cpoints, gp_ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::quad4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            CreateProjected<DRT::Element::quad4>(cpoints, gp_ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::hex8:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            CreateProjected<DRT::Element::hex8>(cpoints, gp_ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::tet4:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            CreateProjected<DRT::Element::tet4>(cpoints, gp_ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::wedge6:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            CreateProjected<DRT::Element::wedge6>(cpoints, gp_ic);
        gpc->Append(gp);
        break;
      }
      case DRT::Element::pyramid5:
      {
        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            CreateProjected<DRT::Element::pyramid5>(cpoints, gp_ic);
        gpc->Append(gp);
        break;
      }
      default:
        dserror("unsupported integration cell type ( cell type = %s )",
            DRT::DistypeToString(ic->Shape()).c_str());
        exit(EXIT_FAILURE);
    }
  }
}

void GEO::CUT::ElementHandle::AppendVolumeCellGaussPoints_MomentFitting(
    Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc, GEO::CUT::VolumeCell* vc)
{
  //-------------------
  // For MomentFitting, we have Gauss points that are calculated w.r to local coordinates of linear
  // shadow element If background ElementHandle is linear, then no need for any mapping Else, we map
  // these points to local coordinates of corresponding Quad element
  //-------------------

  //---------------------------------------------
  const std::vector<GEO::CUT::Point*>& cpoints = vc->ParentElement()->Points();
  Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic = vc->GetGaussRule();


  switch (Shape())
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
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
          CreateProjected<DRT::Element::hex8>(cpoints, gp_ic);
      gpc->Append(gp);
      break;
    }
    case DRT::Element::tet10:
    {
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
          CreateProjected<DRT::Element::tet4>(cpoints, gp_ic);
      gpc->Append(gp);
      break;
    }
    default:
    {
      dserror("element handle for this element is not available\n");
      break;
    }
  }
}


void GEO::CUT::ElementHandle::AppendVolumeCellGaussPoints_DirectDivergence(
    Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc, GEO::CUT::VolumeCell* vc)
{
  //-------------------
  // For DirectDivergence, we calculate Gauss points at the correct local coord. during construction
  // itself This method is handled separately because
  // 1. main Gauss pts should be mapped w.r to each facet of vcell
  //         --> element volume mapping as done for tessellation and moment fitting do not work
  // 2. Internal Gauss pts can be obtained only if we have correctly mapped main Gauss points
  //-------------------
  Teuchos::RCP<DRT::UTILS::GaussPoints> gp = vc->GetGaussRule();

  // volume cell gausspoints are identified to be negligible in
  // GEO::CUT::VolumeCell::DirectDivergenceGaussRule
  if (gp == Teuchos::null) return;

  gpc->Append(gp);
}

/*----------------------------------------------------------------------*/
// Collect the Gaussian points of all the volume-cells belonging to this element.
// The integration rules over all the volume-cells are connected.
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::GaussPointsComposite> GEO::CUT::ElementHandle::GaussPointsConnected(
    plain_volumecell_set& cells, INPAR::CUT::VCellGaussPts gausstype)
{
  Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc =
      Teuchos::rcp(new DRT::UTILS::GaussPointsComposite(0));

  for (plain_volumecell_set::iterator i = cells.begin(); i != cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = *i;

    const plain_integrationcell_set& cells = vc->IntegrationCells();


    if (gausstype == INPAR::CUT::VCellGaussPts_Tessellation)
    {
      for (plain_integrationcell_set::const_iterator i = cells.begin(); i != cells.end(); ++i)
      {
        GEO::CUT::IntegrationCell* ic = *i;

        Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic =
            DRT::UTILS::GaussPointCache::Instance().Create(
                ic->Shape(), ic->CubatureDegree(ic->Shape()));
        const std::vector<GEO::CUT::Point*>& cpoints = ic->Points();

        switch (ic->Shape())
        {
          case DRT::Element::hex8:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
                CreateProjected<DRT::Element::hex8>(cpoints, gp_ic);
            gpc->Append(gp);
            break;
          }
          case DRT::Element::tet4:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
                CreateProjected<DRT::Element::tet4>(cpoints, gp_ic);
            gpc->Append(gp);
            break;
          }
          case DRT::Element::wedge6:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
                CreateProjected<DRT::Element::wedge6>(cpoints, gp_ic);
            gpc->Append(gp);
            break;
          }
          case DRT::Element::pyramid5:
          {
            Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
                CreateProjected<DRT::Element::pyramid5>(cpoints, gp_ic);
            gpc->Append(gp);
            break;
          }
          default:
            dserror("unsupported integration cell type ( cell type = %s )",
                DRT::DistypeToString(ic->Shape()).c_str());
            exit(EXIT_FAILURE);
        }
      }
    }
    else if (gausstype == INPAR::CUT::VCellGaussPts_MomentFitting ||
             gausstype == INPAR::CUT::VCellGaussPts_DirectDivergence)
    {
      Teuchos::RCP<DRT::UTILS::GaussPoints> gp = vc->GetGaussRule();
      gpc->Append(gp);
    }
  }

  return gpc;
}

#if (0)
/*----------------------------------------------------------------------*/
// Unused
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::BoundaryCellGaussPoints(
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,
    std::map<int, std::vector<DRT::UTILS::GaussIntegration>>& intpoints)
{
  for (std::map<int, std::vector<GEO::CUT::BoundaryCell*>>::const_iterator i = bcells.begin();
       i != bcells.end(); ++i)
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*>& cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration>& cell_points = intpoints[sid];

    cell_points.clear();
    cell_points.reserve(cells.size());

    for (std::vector<GEO::CUT::BoundaryCell*>::const_iterator i = cells.begin(); i != cells.end();
         ++i)
    {
      GEO::CUT::BoundaryCell* bc = *i;

      switch (bc->Shape())
      {
        case DRT::Element::tri3:
        {
          cell_points.push_back(
              DRT::UTILS::GaussIntegration(DRT::Element::tri3, bc->CubatureDegree()));
          break;
        }
        case DRT::Element::quad4:
        {
          cell_points.push_back(
              DRT::UTILS::GaussIntegration(DRT::Element::quad4, bc->CubatureDegree()));
          break;
        }
        default:
          throw std::runtime_error("unsupported integration cell type");
      }
    }
  }
}


/*----------------------------------------------------------------------*/
// unused: Old implementation of boundarycell Gauss points collection
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::BoundaryCellGaussPoints(MeshIntersection& mesh, int mi,
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,
    std::map<int, std::vector<DRT::UTILS::GaussIntegration>>& intpoints)
{
  for (std::map<int, std::vector<GEO::CUT::BoundaryCell*>>::const_iterator i = bcells.begin();
       i != bcells.end(); ++i)
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*>& cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration>& cell_points = intpoints[sid];

    SideHandle* side = mesh.GetCutSide(sid, mi);
    if (side == NULL)
    {
      throw std::runtime_error("no such side");
    }

    cell_points.clear();
    cell_points.reserve(cells.size());

    for (std::vector<GEO::CUT::BoundaryCell*>::const_iterator i = cells.begin(); i != cells.end();
         ++i)
    {
      GEO::CUT::BoundaryCell* bc = *i;

      switch (bc->Shape())
      {
        case DRT::Element::tri3:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::tri3>(bc);
          cell_points.push_back(DRT::UTILS::GaussIntegration(gp));
          break;
        }
        case DRT::Element::quad4:
        {
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp = side->CreateProjected<DRT::Element::quad4>(bc);
          cell_points.push_back(DRT::UTILS::GaussIntegration(gp));
          break;
        }
        default:
          dserror("unsupported integration cell type ( cell type = %s )",
              DRT::DistypeToString(bc->Shape()).c_str());
          exit(EXIT_FAILURE);
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
    const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,
    std::map<int, std::vector<DRT::UTILS::GaussIntegration>>& intpoints)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT::ElementHandle::BoundaryCellGaussPointsLin" );

  for (std::map<int, std::vector<GEO::CUT::BoundaryCell*>>::const_iterator i = bcells.begin();
       i != bcells.end(); ++i)
  {
    int sid = i->first;
    const std::vector<GEO::CUT::BoundaryCell*>& cells = i->second;
    std::vector<DRT::UTILS::GaussIntegration>& cell_points = intpoints[sid];

    //    // safety check
    //    if(sid < 0)
    //    {
    //      dserror("invalid sid: %i", sid);
    //    }
    //    else
    //    {
    //      // ask for a mesh cutting side
    //      SideHandle * side = wizard->GetMeshCuttingSide( sid, 0 );
    //
    //      // ask for level-set cutting side
    //      bool is_ls_side = wizard->HasLSCuttingSide( sid );
    //      if ( side==NULL and !is_ls_side )
    //      {
    //        throw std::runtime_error( "no side with given id available fo combined mesh and
    //        level-set cut" );
    //      }
    //    }

    cell_points.clear();
    cell_points.reserve(cells.size());

    for (std::vector<GEO::CUT::BoundaryCell*>::const_iterator i = cells.begin(); i != cells.end();
         ++i)
    {
      GEO::CUT::BoundaryCell* bc = *i;

      // Create (unmodified) gauss points for integration cell with requested
      // polynomial order. This is supposed to be fast, since there is a cache.
      cell_points.push_back(bc->gaussRule());
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::GetBoundaryCellSets(
    const std::vector<GEO::CUT::Point::PointPosition>& desired_positions,
    std::vector<plain_boundarycell_set>& bcellsets)
{
  for (std::vector<GEO::CUT::Point::PointPosition>::const_iterator ip = desired_positions.begin();
       ip != desired_positions.end(); ++ip)
  {
    const std::vector<plain_boundarycell_set>& ele_bcellsets = GetBoundaryCellSet(*ip);

    std::copy(
        ele_bcellsets.begin(), ele_bcellsets.end(), std::inserter(bcellsets, bcellsets.end()));
  }
}

/*----------------------------------------------------------------------*/
// get all the element' sets of volume-cells and nds-vectors
/*----------------------------------------------------------------------*/
void GEO::CUT::ElementHandle::GetVolumeCellsDofSets(std::vector<plain_volumecell_set>& cellsets,
    std::vector<std::vector<int>>& nds_sets, bool include_inner)
{
  const std::vector<plain_volumecell_set>& ele_vc_sets_inside = GetVcSetsInside();
  const std::vector<plain_volumecell_set>& ele_vc_sets_outside = GetVcSetsOutside();

  std::vector<std::vector<int>>& nodaldofset_vc_sets_inside = GetNodalDofSet_VcSets_Inside();
  std::vector<std::vector<int>>& nodaldofset_vc_sets_outside = GetNodalDofSet_VcSets_Outside();

  if (include_inner)
  {
    std::copy(ele_vc_sets_inside.begin(), ele_vc_sets_inside.end(),
        std::inserter(cellsets, cellsets.end()));
    std::copy(nodaldofset_vc_sets_inside.begin(), nodaldofset_vc_sets_inside.end(),
        std::inserter(nds_sets, nds_sets.end()));
  }

  std::copy(ele_vc_sets_outside.begin(), ele_vc_sets_outside.end(),
      std::inserter(cellsets, cellsets.end()));
  std::copy(nodaldofset_vc_sets_outside.begin(), nodaldofset_vc_sets_outside.end(),
      std::inserter(nds_sets, nds_sets.end()));
}



/*----------------------------------------------------------------------*/
//! Collect all volume-cells belonging to this elements
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::GetVolumeCells(plain_volumecell_set& cells)
{
  const plain_volumecell_set& cs = element_->VolumeCells();
  std::copy(cs.begin(), cs.end(), std::inserter(cells, cells.begin()));
}


/*----------------------------------------------------------------------*/
//! Collect all volume-cells belonging to this element ordered by position
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::CollectVolumeCells(
    plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside)
{
  const plain_volumecell_set& ecells = element_->VolumeCells();

  // sort for inside and outside volume cells
  for (plain_volumecell_set::const_iterator i = ecells.begin(); i != ecells.end(); i++)
  {
    if ((*i)->Position() == GEO::CUT::Point::outside)
    {
      cells_outside.insert(*i);
    }
    else  // inside vc
    {
      cells_inside.insert(*i);
    }
  }
}


/*----------------------------------------------------------------------------*
  get all the element sets of volume-cells, nds-vectors and integration points
  return true if a specific XFEM-Gaussrule is available and necessary
 *----------------------------------------------------------------------------*/
bool GEO::CUT::ElementHandle::GetCellSets_DofSets_GaussPoints(
    std::vector<plain_volumecell_set>& cell_sets, std::vector<std::vector<int>>& nds_sets,
    std::vector<std::vector<DRT::UTILS::GaussIntegration>>& intpoints_sets, bool include_inner)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT::ElementHandle::GetCellSets_DofSets_GaussPoints");

  GetVolumeCellsDofSets(cell_sets, nds_sets, include_inner);

  intpoints_sets.clear();

  /* switch this on to only integrate cut elements
   * only one cell_sets for current element and element is not intersected, then
   * use a non-XFEM Gaussrule */
  if (!IsIntersected())
  {
    // perform  standard integration on an uncut element
    if (cell_sets.size() == 1)
      return false;
    else if (cell_sets.size() == 0)
    {
      /* the element does not have a physical but only a ghost set, therefore
       * no integration is necessary and intpoints remains empty */
    }
    else
      dserror("number of cell_sets for a non-intersected element is invalid: %i", cell_sets.size());
  }


  intpoints_sets.reserve(cell_sets.size());

  for (std::vector<plain_volumecell_set>::iterator i = cell_sets.begin(); i != cell_sets.end(); i++)
  {
    plain_volumecell_set& cells = *i;

    std::vector<DRT::UTILS::GaussIntegration> gaussCellsets;
    VolumeCellGaussPoints(cells, gaussCellsets);

    intpoints_sets.push_back(gaussCellsets);
  }

  // return true if specific XFEM-integration rule available
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::BoundaryCellSet(Point::PointPosition position)
{
  // boundary cell sets were already added
  if (bcell_sets_.find(position) != bcell_sets_.end()) return;

  // increase the size of the paired vector and extract the std::vector
  // corresponding to the given position
  const unsigned curr_size = bcell_sets_.size();
  bcell_sets_.resize(curr_size + 1);
  std::vector<plain_boundarycell_set>& bcell_sets = bcell_sets_[position];
  bcell_sets.resize(1, plain_boundarycell_set());

  // get the volume cells of this linear element
  const plain_volumecell_set& evolcells = element_->VolumeCells();

  plain_boundarycell_set& bcells = bcell_sets[0];
  for (plain_volumecell_set::const_iterator citvol = evolcells.begin(); citvol != evolcells.end();
       ++citvol)
  {
    const VolumeCell& evolcell = **citvol;
    if (evolcell.Position() == position)
    {
      const plain_boundarycell_set& ebcells = evolcell.BoundaryCells();
      for (plain_boundarycell_set::const_iterator citbc = ebcells.begin(); citbc != ebcells.end();
           ++citbc)
      {
        // avoid to add boundary cells twice
        if (bcells.find(*citbc) == bcells.end()) bcells.insert(*citbc);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
// get the element's sets of volume-cells ordered by inside/outside position
/*----------------------------------------------------------------------*/
void GEO::CUT::LinearElementHandle::VolumeCellSets()
{
  if (!cells_set_)
  {
    const plain_volumecell_set& ecells = element_->VolumeCells();

    // sort for inside and outside volume cells
    for (plain_volumecell_set::const_iterator i = ecells.begin(); i != ecells.end(); i++)
    {
      if ((*i)->Position() == GEO::CUT::Point::outside)
      {
        plain_volumecell_set s;  // plain volume cell set with only one entry
        s.insert(*i);
        vc_sets_outside_.push_back(s);
      }
      else  // inside vc
      {
        plain_volumecell_set s;  // plain volume cell set with only one entry
        s.insert(*i);
        vc_sets_inside_.push_back(s);
      }
    }

    cells_set_ = true;
  }
}



/*----------------------------------------------------------------------*/
// returns true in case that any cut-side cut with the element produces cut points,
// i.e. also for touched cases (at points, edges or sides),
// or when an element side has more than one facet or is touched by fully/partially by the cut side
/*----------------------------------------------------------------------*/
bool GEO::CUT::QuadraticElementHandle::IsCut()
{
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    if (e->IsCut())
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

  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    if (e->IsIntersected())
    {
      return true;
    }


    // do the sub-elements have different positions ? (e.g. when the cut side directly cuts between
    // sub-elements)
    if (unique_pos == GEO::CUT::Point::undecided)
    {
      // assume a new unique position for all sub elements
      unique_pos = (*e->VolumeCells().begin())->Position();
    }
    else if ((*e->VolumeCells().begin())->Position() != unique_pos)
    {
      return true;
    }
  }
  return false;
}


/*----------------------------------------------------------------------*/
// Collect all volume-cells belonging to this elements
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetVolumeCells(plain_volumecell_set& cells)
{
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    const plain_volumecell_set& cs = e->VolumeCells();
    std::copy(cs.begin(), cs.end(), std::inserter(cells, cells.begin()));
  }
}


/*----------------------------------------------------------------------*/
// Collect all volume-cells belonging to this element ordered by position
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::CollectVolumeCells(
    plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside)
{
  for (std::vector<Element*>::const_iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    const plain_volumecell_set& ecells = e->VolumeCells();

    // sort for inside and outside volume-cells
    for (plain_volumecell_set::const_iterator i = ecells.begin(); i != ecells.end(); i++)
    {
      if ((*i)->Position() == GEO::CUT::Point::outside)
      {
        cells_outside.insert(*i);
      }
      else  // inside vc
      {
        cells_inside.insert(*i);
      }
    }  // volume-cells
  }    // sub-elements
}


/*----------------------------------------------------------------------*/
//  get the quadratic element's volumetric integration cells (just for Tessellation)
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetIntegrationCells(plain_integrationcell_set& cells)
{
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    e->GetIntegrationCells(cells);
  }
}


/*----------------------------------------------------------------------*/
//  get all the quadratic element's boundary integration cells
//  TODO: this has to be corrected such that just the bcs which belong the
//  outside vcs will be returned
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::GetBoundaryCells(plain_boundarycell_set& bcells)
{
  dserror("Deprecated version!");
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    e->GetBoundaryCells(bcells);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::BoundaryCellSet(Point::PointPosition position)
{
  if (connected_bcell_sets_.find(position) != connected_bcell_sets_.end()) return;

  ConnectBoundaryCells(position);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::ConnectBoundaryCells(Point::PointPosition position)
{
  plain_volumecell_set evolcells_position;
  CollectVolumeCells(position, evolcells_position);

  std::vector<plain_volumecell_set> connected_evolcells_position;
  BuildCellSets(evolcells_position, connected_evolcells_position);

  BuildBoundaryCellSets(connected_evolcells_position, connected_bcell_sets_[position]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::BuildBoundaryCellSets(
    const std::vector<plain_volumecell_set>& connected_vcell_set,
    std::vector<plain_boundarycell_set>& connected_bcell_set) const
{
  connected_bcell_set.resize(connected_vcell_set.size(), plain_boundarycell_set());

  unsigned vcell_set_count = 0;
  for (std::vector<plain_volumecell_set>::const_iterator citvset = connected_vcell_set.begin();
       citvset != connected_vcell_set.end(); ++citvset)
  {
    const plain_volumecell_set& vcell_set = *citvset;
    plain_boundarycell_set& bcell_set = connected_bcell_set[vcell_set_count++];
    for (plain_volumecell_set::const_iterator citvc = vcell_set.begin(); citvc != vcell_set.end();
         ++citvc)
    {
      const VolumeCell& vcell = **citvc;
      const plain_boundarycell_set& bcells = vcell.BoundaryCells();
      for (plain_boundarycell_set::const_iterator citbc = bcells.begin(); citbc != bcells.end();
           ++citbc)
      {
        // avoid to add bcells twice
        if (bcell_set.find(*citbc) == bcell_set.end()) bcell_set.insert(*citbc);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::CollectVolumeCells(
    Point::PointPosition position, plain_volumecell_set& evolcells_position) const
{
  for (std::vector<Element*>::const_iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* e = *i;
    const plain_volumecell_set& ecells = e->VolumeCells();

    // sort for inside and outside volume-cells
    for (plain_volumecell_set::const_iterator i = ecells.begin(); i != ecells.end(); i++)
    {
      VolumeCell* evolcell = *i;
      if (evolcell->Position() == position)
      {
        evolcells_position.insert(evolcell);
      }
    }
  }
}



/*----------------------------------------------------------------------*/
//! get the element's sets of volume-cells ordered by inside/outside position
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::VolumeCellSets()
{
  // connect volumecells of subelements
  ConnectVolumeCells();
}


/*----------------------------------------------------------------------*/
//! connect volume-cells to sets of volume-cells
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::ConnectVolumeCells()
{
  // find the connection between volumecells of all subelements for the current element (hex8,
  // hex20, tet4, tet10 etc.) remark: this function determines not the connection outside this
  // element

  // get the volumecells of all subelements stored in two plainvolume_sets (inside and outside)
  if (!cells_connected_)
  {
    plain_volumecell_set e_vcs_inside;
    plain_volumecell_set e_vcs_outside;

    CollectVolumeCells(e_vcs_inside, e_vcs_outside);

    BuildCellSets(e_vcs_inside, connected_vc_sets_inside_);
    BuildCellSets(e_vcs_outside, connected_vc_sets_outside_);

    cells_connected_ = true;
  }
}

/*----------------------------------------------------------------------*/
//! build sets
/*----------------------------------------------------------------------*/
void GEO::CUT::QuadraticElementHandle::BuildCellSets(
    plain_volumecell_set& cells_to_connect, std::vector<plain_volumecell_set>& connected_sets)
{
  plain_volumecell_set done;

  for (plain_volumecell_set::const_iterator i = cells_to_connect.begin();
       i != cells_to_connect.end(); ++i)
  {
    VolumeCell* cell = *i;
    if (done.count(cell) == 0)  // cell currently not-done
    {
      plain_volumecell_set connected;
      // REMARK: here use the version without! elements check:
      // here we build cell sets within one global element with vcs of subelements
      // maybe the vcs of one subelement are not connected within one subelement,
      // but within one global element, therefore more than one vc of one subelements
      // may be connected.
      //        cell->Neighbors( NULL, cells_to_connect, done, connected, elements );
      cell->Neighbors(NULL, cells_to_connect, done, connected);

      if (connected.size() > 0)
      {
        connected_sets.push_back(connected);
        std::copy(connected.begin(), connected.end(), std::inserter(done, done.begin()));
      }
    }
  }
}


/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Hex20ElementHandle::Hex20ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nodes)
    : QuadraticElementHandle()
{
  subelements_.reserve(8);

  const CellTopologyData* top_data = shards::getCellTopologyData<shards::Hexahedron<8>>();

  // create middle nodes

  LINALG::Matrix<3, 1> xyz;
  LINALG::Matrix<1, 1> lsv;

  plain_int_set node_nids;

  LINALG::Matrix<3, 8> side_xyze;
  LINALG::Matrix<1, 8> side_lsvs;
  LINALG::Matrix<8, 1> side_funct;
  std::vector<Node*> side_nodes(8);

  std::vector<Node*> center_nodes(6);

  // get the real nodes of the hex20 element and get the center nodes (shadow nodes) of its
  // quadratic sides which stored in the shadow_nodes list in cut_mesh.H using the eight side nodes
  // of the side as key special handling for the inner center node of hex20 element required, see
  // below, the inner center node is stored in the shadow_nodes map using the 20 nodes of the hex20
  // element as key

  // loop the six sides of the hex 20 element
  for (int localsideid = 0; localsideid < 6; ++localsideid)
  {
    node_nids.clear();
    // loop the eight nodes of each quad8 side of the hex20 element
    for (int i = 0; i < 8; ++i)
    {
      int localnodeid = DRT::UTILS::eleNodeNumbering_hex27_surfaces[localsideid][i];
      Node* n = mesh.GetNode(nodes[localnodeid], static_cast<double*>(NULL));
      side_nodes[i] = n;
      node_nids.insert(nodes[localnodeid]);
      n->Coordinates(&side_xyze(0, i));
      side_lsvs(i) = n->LSV();
    }

    DRT::UTILS::shape_function_2D(side_funct, 0.0, 0.0, DRT::Element::quad8);
    xyz.Multiply(side_xyze, side_funct);
    lsv.Multiply(side_lsvs, side_funct);

    // find the unique center node of the quadratic quad8 sides
    center_nodes[localsideid] = mesh.GetNode(node_nids, xyz.A(), lsv(0));
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

  LINALG::Matrix<3, 20> xyze;
  LINALG::Matrix<1, 20> lsvs;
  nodes_.reserve(20);
  for (int i = 0; i < 20; ++i)
  {
    Node* n = mesh.GetNode(nodes[i], static_cast<double*>(NULL));
    nodes_.push_back(n);
    n->Coordinates(&xyze(0, i));
    lsvs(i) = n->LSV();
  }

  // special handling for the inner center node of hex20 element
  // Remark: this node is also stored as shadow node in cut_mesh, however, the key for this shadow
  // node is the set of all 20 nodes of the hex20 element in contrast to the shadow nodes of sides,
  // for that the key are the eight nodes of the quad7 side
  LINALG::Matrix<20, 1> funct;
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, DRT::Element::hex20);

  xyz.Multiply(xyze, funct);
  lsv.Multiply(lsvs, funct);
  node_nids.clear();
  std::copy(nodes.begin(), nodes.end(), std::inserter(node_nids, node_nids.begin()));
  Node* node26 = mesh.GetNode(node_nids, xyz.A(), lsv(0));
  int node26_id = node26->Id();


  std::vector<int> nids(8);

  nids[0] = nodes[0];
  nids[1] = nodes[8];
  nids[2] = node20_id;
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = node21_id;
  nids[6] = node26_id;
  nids[7] = node24_id;
  Element* sub1 = mesh.GetElement(-1, nids, *top_data);
  sub1->setAsShadowElem();
  sub1->setQuadCorners(mesh, nodes);
  sub1->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub1);

  nids[0] = nodes[8];
  nids[1] = nodes[1];
  nids[2] = nodes[9];
  nids[3] = node20_id;
  nids[4] = node21_id;
  nids[5] = nodes[13];
  nids[6] = node22_id;
  nids[7] = node26_id;
  Element* sub2 = mesh.GetElement(-1, nids, *top_data);
  sub2->setAsShadowElem();
  sub2->setQuadCorners(mesh, nodes);
  sub2->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub2);

  nids[0] = node20_id;
  nids[1] = nodes[9];
  nids[2] = nodes[2];
  nids[3] = nodes[10];
  nids[4] = node26_id;
  nids[5] = node22_id;
  nids[6] = nodes[14];
  nids[7] = node23_id;
  Element* sub3 = mesh.GetElement(-1, nids, *top_data);
  sub3->setAsShadowElem();
  sub3->setQuadCorners(mesh, nodes);
  sub3->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub3);

  nids[0] = nodes[11];
  nids[1] = node20_id;
  nids[2] = nodes[10];
  nids[3] = nodes[3];
  nids[4] = node24_id;
  nids[5] = node26_id;
  nids[6] = node23_id;
  nids[7] = nodes[15];
  Element* sub4 = mesh.GetElement(-1, nids, *top_data);
  sub4->setAsShadowElem();
  sub4->setQuadCorners(mesh, nodes);
  sub4->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub4);

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = node21_id;
  nids[2] = node26_id;
  nids[3] = node24_id;
  nids[4] = nodes[4];
  nids[5] = nodes[16];
  nids[6] = node25_id;
  nids[7] = nodes[19];
  Element* sub5 = mesh.GetElement(-1, nids, *top_data);
  sub5->setAsShadowElem();
  sub5->setQuadCorners(mesh, nodes);
  sub5->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub5);

  nids[0] = node21_id;
  nids[1] = nodes[13];
  nids[2] = node22_id;
  nids[3] = node26_id;
  nids[4] = nodes[16];
  nids[5] = nodes[5];
  nids[6] = nodes[17];
  nids[7] = node25_id;
  Element* sub6 = mesh.GetElement(-1, nids, *top_data);
  sub6->setAsShadowElem();
  sub6->setQuadCorners(mesh, nodes);
  sub6->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub6);

  nids[0] = node26_id;
  nids[1] = node22_id;
  nids[2] = nodes[14];
  nids[3] = node23_id;
  nids[4] = node25_id;
  nids[5] = nodes[17];
  nids[6] = nodes[6];
  nids[7] = nodes[18];
  Element* sub7 = mesh.GetElement(-1, nids, *top_data);
  sub7->setAsShadowElem();
  sub7->setQuadCorners(mesh, nodes);
  sub7->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub7);

  nids[0] = node24_id;
  nids[1] = node26_id;
  nids[2] = node23_id;
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = node25_id;
  nids[6] = nodes[18];
  nids[7] = nodes[7];
  Element* sub8 = mesh.GetElement(-1, nids, *top_data);
  sub8->setAsShadowElem();
  sub8->setQuadCorners(mesh, nodes);
  sub8->setQuadShape(DRT::Element::hex20);
  subelements_.push_back(sub8);

  // each subelement should know its parents id
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* subelement = *i;
    subelement->ParentId(eid);
  }
}


/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Hex27ElementHandle::Hex27ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nodes)
    : QuadraticElementHandle()
{
  subelements_.reserve(8);

  nodes_.reserve(27);
  for (int i = 0; i < 27; ++i)
  {
    Node* n = mesh.GetNode(nodes[i], static_cast<double*>(NULL));
    nodes_.push_back(n);
  }

  const CellTopologyData* top_data = shards::getCellTopologyData<shards::Hexahedron<8>>();

  std::vector<int> nids(8);

  nids[0] = nodes[0];
  nids[1] = nodes[8];
  nids[2] = nodes[20];
  nids[3] = nodes[11];
  nids[4] = nodes[12];
  nids[5] = nodes[21];
  nids[6] = nodes[26];
  nids[7] = nodes[24];
  Element* sub1 = mesh.GetElement(-1, nids, *top_data);
  sub1->setAsShadowElem();
  sub1->setQuadCorners(mesh, nodes);
  sub1->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub1);

  nids[0] = nodes[8];
  nids[1] = nodes[1];
  nids[2] = nodes[9];
  nids[3] = nodes[20];
  nids[4] = nodes[21];
  nids[5] = nodes[13];
  nids[6] = nodes[22];
  nids[7] = nodes[26];
  Element* sub2 = mesh.GetElement(-1, nids, *top_data);
  sub2->setAsShadowElem();
  sub2->setQuadCorners(mesh, nodes);
  sub2->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub2);

  nids[0] = nodes[20];
  nids[1] = nodes[9];
  nids[2] = nodes[2];
  nids[3] = nodes[10];
  nids[4] = nodes[26];
  nids[5] = nodes[22];
  nids[6] = nodes[14];
  nids[7] = nodes[23];
  Element* sub3 = mesh.GetElement(-1, nids, *top_data);
  sub3->setAsShadowElem();
  sub3->setQuadCorners(mesh, nodes);
  sub3->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub3);

  nids[0] = nodes[11];
  nids[1] = nodes[20];
  nids[2] = nodes[10];
  nids[3] = nodes[3];
  nids[4] = nodes[24];
  nids[5] = nodes[26];
  nids[6] = nodes[23];
  nids[7] = nodes[15];
  Element* sub4 = mesh.GetElement(-1, nids, *top_data);
  sub4->setAsShadowElem();
  sub4->setQuadCorners(mesh, nodes);
  sub4->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub4);

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[12];
  nids[1] = nodes[21];
  nids[2] = nodes[26];
  nids[3] = nodes[24];
  nids[4] = nodes[4];
  nids[5] = nodes[16];
  nids[6] = nodes[25];
  nids[7] = nodes[19];
  Element* sub5 = mesh.GetElement(-1, nids, *top_data);
  sub5->setAsShadowElem();
  sub5->setQuadCorners(mesh, nodes);
  sub5->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub5);

  nids[0] = nodes[21];
  nids[1] = nodes[13];
  nids[2] = nodes[22];
  nids[3] = nodes[26];
  nids[4] = nodes[16];
  nids[5] = nodes[5];
  nids[6] = nodes[17];
  nids[7] = nodes[25];
  Element* sub6 = mesh.GetElement(-1, nids, *top_data);
  sub6->setAsShadowElem();
  sub6->setQuadCorners(mesh, nodes);
  sub6->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub6);

  nids[0] = nodes[26];
  nids[1] = nodes[22];
  nids[2] = nodes[14];
  nids[3] = nodes[23];
  nids[4] = nodes[25];
  nids[5] = nodes[17];
  nids[6] = nodes[6];
  nids[7] = nodes[18];
  Element* sub7 = mesh.GetElement(-1, nids, *top_data);
  sub7->setAsShadowElem();
  sub7->setQuadCorners(mesh, nodes);
  sub7->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub7);

  nids[0] = nodes[24];
  nids[1] = nodes[26];
  nids[2] = nodes[23];
  nids[3] = nodes[15];
  nids[4] = nodes[19];
  nids[5] = nodes[25];
  nids[6] = nodes[18];
  nids[7] = nodes[7];
  Element* sub8 = mesh.GetElement(-1, nids, *top_data);
  sub8->setAsShadowElem();
  sub8->setQuadCorners(mesh, nodes);
  sub8->setQuadShape(DRT::Element::hex27);
  subelements_.push_back(sub8);

  // each subelement should know its parents id
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* subelement = *i;
    subelement->ParentId(eid);
  }
}


/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Tet10ElementHandle::Tet10ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nids)
    : QuadraticElementHandle()
{
  subelements_.reserve(8);

  nodes_.reserve(10);
  for (int i = 0; i < 10; ++i)
  {
    Node* n = mesh.GetNode(nids[i], static_cast<double*>(NULL));
    nodes_.push_back(n);
  }

  const CellTopologyData* top_data = shards::getCellTopologyData<shards::Tetrahedron<4>>();

  std::vector<int> subnids(4);

  subnids[0] = nids[0];
  subnids[1] = nids[4];
  subnids[2] = nids[6];
  subnids[3] = nids[7];
  Element* sub1 = mesh.GetElement(-1, subnids, *top_data);
  sub1->setAsShadowElem();
  sub1->setQuadCorners(mesh, nids);
  sub1->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub1);

  subnids[0] = nids[4];
  subnids[1] = nids[1];
  subnids[2] = nids[5];
  subnids[3] = nids[8];
  Element* sub2 = mesh.GetElement(-1, subnids, *top_data);
  sub2->setAsShadowElem();
  sub2->setQuadCorners(mesh, nids);
  sub2->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub2);

  subnids[0] = nids[6];
  subnids[1] = nids[5];
  subnids[2] = nids[2];
  subnids[3] = nids[9];
  Element* sub3 = mesh.GetElement(-1, subnids, *top_data);
  sub3->setAsShadowElem();
  sub3->setQuadCorners(mesh, nids);
  sub3->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub3);

  subnids[0] = nids[7];
  subnids[1] = nids[8];
  subnids[2] = nids[9];
  subnids[3] = nids[3];
  Element* sub4 = mesh.GetElement(-1, subnids, *top_data);
  sub4->setAsShadowElem();
  sub4->setQuadCorners(mesh, nids);
  sub4->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub4);

  /////////////////////////////////////////////////////////////////

  subnids[0] = nids[4];
  subnids[1] = nids[5];
  subnids[2] = nids[6];
  subnids[3] = nids[8];
  Element* sub5 = mesh.GetElement(-1, subnids, *top_data);
  sub5->setAsShadowElem();
  sub5->setQuadCorners(mesh, nids);
  sub5->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub5);

  subnids[0] = nids[6];
  subnids[1] = nids[9];
  subnids[2] = nids[7];
  subnids[3] = nids[8];
  Element* sub6 = mesh.GetElement(-1, subnids, *top_data);
  sub6->setAsShadowElem();
  sub6->setQuadCorners(mesh, nids);
  sub6->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub6);

  subnids[0] = nids[4];
  subnids[1] = nids[6];
  subnids[2] = nids[7];
  subnids[3] = nids[8];
  Element* sub7 = mesh.GetElement(-1, subnids, *top_data);
  sub7->setAsShadowElem();
  sub7->setQuadCorners(mesh, nids);
  sub7->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub7);

  subnids[0] = nids[9];
  subnids[1] = nids[6];
  subnids[2] = nids[5];
  subnids[3] = nids[8];
  Element* sub8 = mesh.GetElement(-1, subnids, *top_data);
  sub8->setAsShadowElem();
  sub8->setQuadCorners(mesh, nids);
  sub8->setQuadShape(DRT::Element::tet10);
  subelements_.push_back(sub8);

  // each subelement should know its parents id
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* subelement = *i;
    subelement->ParentId(eid);
  }
}

/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
GEO::CUT::Wedge15ElementHandle::Wedge15ElementHandle(
    Mesh& mesh, int eid, const std::vector<int>& nodes)
    : QuadraticElementHandle()
{
  subelements_.reserve(8);  // subdivide into 8 wedge 6 elements

  const CellTopologyData* top_data = shards::getCellTopologyData<shards::Wedge<6>>();

  // create middle nodes

  LINALG::Matrix<3, 1> xyz;
  LINALG::Matrix<1, 1> lsv;

  plain_int_set node_nids;

  LINALG::Matrix<3, 8> side_xyze;
  LINALG::Matrix<1, 8> side_lsvs;
  LINALG::Matrix<8, 1> side_funct;
  std::vector<Node*> side_nodes(8);

  std::vector<Node*> center_nodes(3);

  // get the real nodes of the wedge15 element and get the center nodes (shadow nodes) of its
  // quadratic quad8 sides which stored in the shadow_nodes list in cut_mesh.H using the eight side
  // nodes of the side as key special handling for the inner center node of wedge15 element
  // required, see below,

  // loop the 5 sides of the wedge15 element

  // three quadratic sides
  for (int localsideid = 0; localsideid < 3; ++localsideid)
  {
    node_nids.clear();
    // loop the eight nodes of each quad8 side of the wedge15 element
    for (int i = 0; i < 8; ++i)
    {
      int localnodeid = DRT::UTILS::eleNodeNumbering_wedge18_quadsurfaces[localsideid][i];
      Node* n = mesh.GetNode(nodes[localnodeid], static_cast<double*>(NULL));
      side_nodes[i] = n;
      node_nids.insert(nodes[localnodeid]);
      n->Coordinates(&side_xyze(0, i));
      side_lsvs(i) = n->LSV();
    }

    DRT::UTILS::shape_function_2D(side_funct, 0.0, 0.0, DRT::Element::quad8);
    xyz.Multiply(side_xyze, side_funct);
    lsv.Multiply(side_lsvs, side_funct);

    // find the unique center node of the quadratic quad8 sides
    center_nodes[localsideid] = mesh.GetNode(node_nids, xyz.A(), lsv(0));
  }

  LINALG::Matrix<3, 6> tb_side_xyze;  // top_bottom_sides
  LINALG::Matrix<1, 6> tb_side_lsvs;
  LINALG::Matrix<6, 1> tb_side_funct;
  std::vector<Node*> tb_side_nodes(6);

  // two quadratic sides on top and bottom
  for (int localsideid = 0; localsideid < 2; ++localsideid)
  {
    node_nids.clear();
    // loop the 6 nodes of each tri6 side of the wedge15 element
    for (int i = 0; i < 6; ++i)
    {
      int localnodeid = DRT::UTILS::eleNodeNumbering_wedge18_trisurfaces[localsideid][i];
      Node* n = mesh.GetNode(nodes[localnodeid], static_cast<double*>(NULL));
      tb_side_nodes[i] = n;
      node_nids.insert(nodes[localnodeid]);
      n->Coordinates(&tb_side_xyze(0, i));
      tb_side_lsvs(i) = n->LSV();
    }
  }

  Node* node15 = center_nodes[0];
  int node15_id = node15->Id();

  Node* node16 = center_nodes[1];
  int node16_id = node16->Id();

  Node* node17 = center_nodes[2];
  int node17_id = node17->Id();


  LINALG::Matrix<3, 15> xyze;
  LINALG::Matrix<1, 15> lsvs;
  nodes_.reserve(15);
  for (int i = 0; i < 15; ++i)
  {
    Node* n = mesh.GetNode(nodes[i], static_cast<double*>(NULL));
    nodes_.push_back(n);
    n->Coordinates(&xyze(0, i));
    lsvs(i) = n->LSV();
  }


  std::vector<int> nids(6);

  nids[0] = nodes[0];
  nids[1] = nodes[6];
  nids[2] = nodes[8];
  nids[3] = nodes[9];
  nids[4] = node15_id;
  nids[5] = node17_id;
  Element* sub1 = mesh.GetElement(-1, nids, *top_data);
  sub1->setAsShadowElem();
  sub1->setQuadCorners(mesh, nodes);
  sub1->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub1);

  nids[0] = nodes[6];
  nids[1] = nodes[1];
  nids[2] = nodes[7];
  nids[3] = node15_id;
  nids[4] = nodes[10];
  nids[5] = node16_id;
  Element* sub2 = mesh.GetElement(-1, nids, *top_data);
  sub2->setAsShadowElem();
  sub2->setQuadCorners(mesh, nodes);
  sub2->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub2);

  nids[0] = nodes[6];
  nids[1] = nodes[7];
  nids[2] = nodes[8];
  nids[3] = node15_id;
  nids[4] = node16_id;
  nids[5] = node17_id;
  Element* sub3 = mesh.GetElement(-1, nids, *top_data);
  sub3->setAsShadowElem();
  sub3->setQuadCorners(mesh, nodes);
  sub3->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub3);

  nids[0] = nodes[8];
  nids[1] = nodes[7];
  nids[2] = nodes[2];
  nids[3] = node17_id;
  nids[4] = node16_id;
  nids[5] = nodes[11];
  Element* sub4 = mesh.GetElement(-1, nids, *top_data);
  sub4->setAsShadowElem();
  sub4->setQuadCorners(mesh, nodes);
  sub4->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub4);

  /////////////////////////////////////////////////////////////////

  nids[0] = nodes[9];
  nids[1] = node15_id;
  nids[2] = node17_id;
  nids[3] = nodes[3];
  nids[4] = nodes[12];
  nids[5] = nodes[14];
  Element* sub5 = mesh.GetElement(-1, nids, *top_data);
  sub5->setAsShadowElem();
  sub5->setQuadCorners(mesh, nodes);
  sub5->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub5);

  nids[0] = node15_id;
  nids[1] = nodes[10];
  nids[2] = node16_id;
  nids[3] = nodes[12];
  nids[4] = nodes[4];
  nids[5] = nodes[13];
  Element* sub6 = mesh.GetElement(-1, nids, *top_data);
  sub6->setAsShadowElem();
  sub6->setQuadCorners(mesh, nodes);
  sub6->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub6);

  nids[0] = node15_id;
  nids[1] = node16_id;
  nids[2] = node17_id;
  nids[3] = nodes[12];
  nids[4] = nodes[13];
  nids[5] = nodes[14];
  Element* sub7 = mesh.GetElement(-1, nids, *top_data);
  sub7->setAsShadowElem();
  sub7->setQuadCorners(mesh, nodes);
  sub7->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub7);

  nids[0] = node17_id;
  nids[1] = node16_id;
  nids[2] = nodes[11];
  nids[3] = nodes[14];
  nids[4] = nodes[13];
  nids[5] = nodes[5];
  Element* sub8 = mesh.GetElement(-1, nids, *top_data);
  sub8->setAsShadowElem();
  sub8->setQuadCorners(mesh, nodes);
  sub8->setQuadShape(DRT::Element::wedge6);
  subelements_.push_back(sub8);

  // each subelement should know its parents id
  for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
  {
    Element* subelement = *i;
    subelement->ParentId(eid);
  }
}


/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Hex20ElementHandle::LocalCoordinates(
    const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst)
{
  Teuchos::RCP<GEO::CUT::Position> pos =
      GEO::CUT::PositionFactory::BuildPosition<3, DRT::Element::hex20>(nodes_, xyz);

  bool success = pos->Compute(1e-10);
  ;
  if (not success)
  {
    std::cout << "local coordinates for hex20 element could not be determined" << std::endl;
    for (int i = 0; i < (int)(nodes_.size()); i++)
    {
      std::cout << " node " << i << std::endl;
      nodes_[i]->Plot(std::cout);
    }

    std::cout << "point in xyz: " << xyz << std::endl;

    dserror("local coordinates for hex20 element could not be determined");
  }
  pos->LocalCoordinates(rst);
}


/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Hex27ElementHandle::LocalCoordinates(
    const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst)
{
  Teuchos::RCP<GEO::CUT::Position> pos =
      GEO::CUT::PositionFactory::BuildPosition<3, DRT::Element::hex27>(nodes_, xyz);

  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}


/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Tet10ElementHandle::LocalCoordinates(
    const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst)
{
  Teuchos::RCP<GEO::CUT::Position> pos =
      GEO::CUT::PositionFactory::BuildPosition<3, DRT::Element::tet10>(nodes_, xyz);

  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}

/*----------------------------------------------------------------------*/
// compute local coordinates of the element for given global coordinates
/*----------------------------------------------------------------------*/
void GEO::CUT::Wedge15ElementHandle::LocalCoordinates(
    const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst)
{
  Teuchos::RCP<GEO::CUT::Position> pos =
      GEO::CUT::PositionFactory::BuildPosition<3, DRT::Element::wedge15>(nodes_, xyz);

  bool success = pos->Compute();
  if (not success)
  {
  }
  pos->LocalCoordinates(rst);
}
