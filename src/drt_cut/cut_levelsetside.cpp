/*---------------------------------------------------------------------*/
/*!

\brief for intersection with an levelset, levelsetside represents the surface described by the
levelset

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include "cut_levelsetside.H"
#include "cut_mesh.H"
#include "cut_element.H"
#include "cut_pointgraph.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int probdim>
bool GEO::CUT::LevelSetSide<probdim>::Cut(Mesh& mesh, Edge& edge, PointSet& cut_points)
{
  return edge.LevelSetCut(mesh, *this, cut_points);
}

template <int probdim>
bool GEO::CUT::LevelSetSide<probdim>::FindCutPointsDispatch(
    Mesh& mesh, Element* element, Side& side, Edge& e)
{
  return e.FindCutPointsLevelSet(mesh, element, side, *this);
}


#if 0
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template < int probdim >
void GEO::CUT::LevelSetSide<probdim>::MakeSideCutFacets(
    Mesh & mesh, Element * element, plain_facet_set & facets )
{
  Side::MakeSideCutFacets( mesh, element, facets );
}
#endif

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int probdim>
void GEO::CUT::LevelSetSide<probdim>::MakeInternalFacets(
    Mesh& mesh, Element* element, plain_facet_set& facets)
{
  Teuchos::RCP<IMPL::PointGraph> pg = Teuchos::rcp(IMPL::PointGraph::Create(
      mesh, element, this, IMPL::PointGraph::cut_side, IMPL::PointGraph::own_lines));

  for (IMPL::PointGraph::facet_iterator i = pg->fbegin(); i != pg->fend(); ++i)
  {
    const Cycle& points = *i;
    Side::MakeInternalFacets(mesh, element, points, facets);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int probdim>
bool GEO::CUT::LevelSetSide<probdim>::FindAmbiguousCutLines(
    Mesh& mesh, Element* element, Side& side, const PointSet& cut)
{
  // More than two cut points show a touch.
  //
  // (1) If all nodes are caught and nothing else, the cut surface has hit this
  //     surface exactly. No need to cut anything. However, the surface might be
  //     required for integration.
  if (GEO::CUT::Side::FindTouchingCutLines(mesh, element, side, cut)) return true;

  switch (side.Shape())
  {
    case DRT::Element::line2:
    case DRT::Element::tri3:
      return false;
    case DRT::Element::quad4:
    {
      switch (cut.size())
      {
        case 3:
        {
          std::vector<Point*> edge_points;
          edge_points.reserve(2);
          for (PointSet::const_iterator i = cut.begin(); i != cut.end(); ++i)
          {
            Point* p = *i;
            if (not p->NodalPoint(element->Nodes()))
            {
              edge_points.push_back(p);
            }
          }
          if (edge_points.size() == 2)
          {
            mesh.NewLine(edge_points[0], edge_points[1], &side, this, element);
            std::cout << "WARNING: levelset cut on node not defined\n";
            return true;
          }
          else if (edge_points.size() == 0)
          {
            const std::vector<Edge*>& edges = side.Edges();
            for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end(); ++i)
            {
              Edge* e = *i;
              Point* p1 = e->BeginNode()->point();
              Point* p2 = e->EndNode()->point();

              if (cut.count(p1) > 0 and cut.count(p2) == 0)
              {
                edge_points.push_back(p1);
              }
              else if (cut.count(p1) == 0 and cut.count(p2) > 0)
              {
                edge_points.push_back(p2);
              }
            }
            if (edge_points.size() == 2)
            {
              mesh.NewLine(edge_points[0], edge_points[1], &side, this, element);
              return true;
            }
          }
          dserror("expect two edge points");
          break;
        }
        case 4:
        {
          // First, check if all cutpoints cut the edges of the side-element.
          std::vector<Point*> edge_points;
          edge_points.reserve(4);
          const std::vector<Edge*>& edges = side.Edges();
          PointSet cut_points(cut);
          for (std::vector<Edge*>::const_iterator i = edges.begin(); i != edges.end(); ++i)
          {
            Edge* e = *i;
            for (PointSet::iterator i = cut_points.begin(); i != cut_points.end();)
            {
              Point* p = *i;
              if (p->IsCut(e))
              {
                edge_points.push_back(p);
                // cut_points.erase( i++ );
                set_erase(cut_points, i);
                break;
              }
              else
              {
                ++i;
              }
            }
          }
          if (edge_points.size() != 4 or cut_points.size() != 0)
          {
            throw std::runtime_error("failed to associate cut points with edges");
          }

          // find levelset value at side center
          LINALG::Matrix<4, 1> lsv;
          LINALG::Matrix<4, 1> funct;
          DRT::UTILS::shape_function_2D(funct, 0., 0., DRT::Element::quad4);
          const std::vector<Node*>& nodes = side.Nodes();
          std::vector<int> zero_positions;
          for (unsigned i = 0; i < 4; ++i)
          {
            lsv(i) = nodes[i]->LSV();
            if (lsv(i) == 0) zero_positions.push_back(i);
          }

          LINALG::Matrix<1, 1> midlsv;
          midlsv.MultiplyTN(lsv, funct);

          for (unsigned i = 1; i < zero_positions.size(); ++i)
          {
            int diff = zero_positions[i] - zero_positions[i - 1];
            if (diff == 1 or diff == 3)
            {
              throw std::runtime_error("cannot have adjacent zeros");
            }
          }

          bool negativemiddle;

          if (midlsv(0) < 0)
          {
            negativemiddle = true;
          }
          else if (midlsv(0) > 0)
          {
            negativemiddle = false;
          }
          else
          {
            // throw std::runtime_error( "side center at interface on multiple cuts: undefined" );
            return false;
          }

#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
          LINALG::Matrix<3, 1> coords0;
          bool connect01and23;

          edge_points[0]->Coordinates(&coords0(0, 0));
          std::vector<double> grad_phi0 = element->GetLevelSetGradient(coords0);

          LINALG::Matrix<3, 1> coords1;
          edge_points[1]->Coordinates(&coords1(0, 0));
          std::vector<double> grad_phi1 = element->GetLevelSetGradient(coords1);

          double dotProduct01 = grad_phi0[0] * grad_phi1[0] + grad_phi0[1] * grad_phi1[1] +
                                grad_phi0[2] * grad_phi1[2];

          LINALG::Matrix<3, 1> coords2;
          edge_points[2]->Coordinates(&coords2(0, 0));
          std::vector<double> grad_phi2 = element->GetLevelSetGradient(coords2);

          LINALG::Matrix<3, 1> coords3;
          edge_points[3]->Coordinates(&coords3(0, 0));
          std::vector<double> grad_phi3 = element->GetLevelSetGradient(coords3);

          double dotProduct23 = grad_phi2[0] * grad_phi3[0] + grad_phi2[1] * grad_phi3[1] +
                                grad_phi2[2] * grad_phi3[2];

          double dotProduct;

          if (fabs(dotProduct01) > fabs(dotProduct23))
            dotProduct = dotProduct01;
          else
            dotProduct = dotProduct23;

#ifdef DEBUGCUTLIBRARY
          std::cout << "dotProduct01: " << dotProduct01 << ", dotProduct23 " << dotProduct23
                    << std::endl;
          if (dotProduct01 * dotProduct23 < 0.0)
            std::cout << "WARNING: dotProduct not unique!!!" << std::endl;
#endif

          if (dotProduct > 0.0)
            connect01and23 = true;
          else
            connect01and23 = false;
#endif

          if (lsv(0) <= 0 and lsv(1) >= 0 and lsv(2) <= 0 and lsv(3) >= 0)
          {
#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
#ifdef DEBUGCUTLIBRARY
            if (negativemiddle != connect01and23)
              std::cout << "Changed from previous configuration of midpoint evaluation!"
                        << std::endl;
#endif

            negativemiddle = connect01and23;
#endif

            if (negativemiddle)
            {
              mesh.NewLine(edge_points[0], edge_points[1], &side, this, element);
              mesh.NewLine(edge_points[2], edge_points[3], &side, this, element);
            }
            else
            {
              mesh.NewLine(edge_points[0], edge_points[3], &side, this, element);
              mesh.NewLine(edge_points[2], edge_points[1], &side, this, element);
            }
            return true;
          }
          else if (lsv(0) >= 0 and lsv(1) <= 0 and lsv(2) >= 0 and lsv(3) <= 0)
          {
#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
#ifdef DEBUGCUTLIBRARY
            if (negativemiddle == connect01and23)
              std::cout << "Changed from previous configuration of midpoint evaluation!"
                        << std::endl;
#endif

            negativemiddle = !connect01and23;
#endif
            if (negativemiddle)
            {
              mesh.NewLine(edge_points[0], edge_points[3], &side, this, element);
              mesh.NewLine(edge_points[2], edge_points[1], &side, this, element);
            }
            else
            {
              mesh.NewLine(edge_points[0], edge_points[1], &side, this, element);
              mesh.NewLine(edge_points[2], edge_points[3], &side, this, element);
            }
            return true;
          }
          else
          {
            throw std::runtime_error("illegal levelset pattern");
          }

          return false;
          break;
        }
        default:
          return false;
      }
      break;
    }  // case DRT::Element::quad4:
    default:
      dserror("Unsupported side shape! (shape = %d | %s )", side.Shape(),
          DRT::DistypeToString(side.Shape()).c_str());
      break;
  }
  return false;
}

template class GEO::CUT::LevelSetSide<2>;
template class GEO::CUT::LevelSetSide<3>;
