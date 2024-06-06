/*---------------------------------------------------------------------*/
/*! \file

\brief for intersection with an levelset, levelsetside represents the surface described by the
levelset

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_levelsetside.hpp"

#include "4C_cut_element.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_pointgraph.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int probdim>
bool Core::Geo::Cut::LevelSetSide<probdim>::Cut(Mesh& mesh, Edge& edge, PointSet& cut_points)
{
  return edge.LevelSetCut(mesh, *this, cut_points);
}

template <int probdim>
bool Core::Geo::Cut::LevelSetSide<probdim>::find_cut_points_dispatch(
    Mesh& mesh, Element* element, Side& side, Edge& e)
{
  return e.find_cut_points_level_set(mesh, element, side, *this);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int probdim>
void Core::Geo::Cut::LevelSetSide<probdim>::MakeInternalFacets(
    Mesh& mesh, Element* element, plain_facet_set& facets)
{
  Teuchos::RCP<Impl::PointGraph> pg = Teuchos::rcp(Impl::PointGraph::Create(
      mesh, element, this, Impl::PointGraph::cut_side, Impl::PointGraph::own_lines));

  for (Impl::PointGraph::facet_iterator i = pg->fbegin(); i != pg->fend(); ++i)
  {
    const Cycle& points = *i;
    Side::MakeInternalFacets(mesh, element, points, facets);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <int probdim>
bool Core::Geo::Cut::LevelSetSide<probdim>::find_ambiguous_cut_lines(
    Mesh& mesh, Element* element, Side& side, const PointSet& cut)
{
  // More than two cut points show a touch.
  //
  // (1) If all nodes are caught and nothing else, the cut surface has hit this
  //     surface exactly. No need to cut anything. However, the surface might be
  //     required for integration.
  if (Core::Geo::Cut::Side::find_touching_cut_lines(mesh, element, side, cut)) return true;

  switch (side.Shape())
  {
    case Core::FE::CellType::line2:
    case Core::FE::CellType::tri3:
      return false;
    case Core::FE::CellType::quad4:
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
          FOUR_C_THROW("expect two edge points");
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
            FOUR_C_THROW("failed to associate cut points with edges");
          }

          // find levelset value at side center
          Core::LinAlg::Matrix<4, 1> lsv;
          Core::LinAlg::Matrix<4, 1> funct;
          Core::FE::shape_function_2D(funct, 0., 0., Core::FE::CellType::quad4);
          const std::vector<Node*>& nodes = side.Nodes();
          std::vector<int> zero_positions;
          for (unsigned i = 0; i < 4; ++i)
          {
            lsv(i) = nodes[i]->LSV();
            if (lsv(i) == 0) zero_positions.push_back(i);
          }

          Core::LinAlg::Matrix<1, 1> midlsv;
          midlsv.MultiplyTN(lsv, funct);

          for (unsigned i = 1; i < zero_positions.size(); ++i)
          {
            int diff = zero_positions[i] - zero_positions[i - 1];
            if (diff == 1 or diff == 3)
            {
              FOUR_C_THROW("cannot have adjacent zeros");
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
            // FOUR_C_THROW( "side center at interface on multiple cuts: undefined" );
            return false;
          }

#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
          Core::LinAlg::Matrix<3, 1> coords0;
          bool connect01and23;

          edge_points[0]->Coordinates(&coords0(0, 0));
          std::vector<double> grad_phi0 = element->get_level_set_gradient(coords0);

          Core::LinAlg::Matrix<3, 1> coords1;
          edge_points[1]->Coordinates(&coords1(0, 0));
          std::vector<double> grad_phi1 = element->get_level_set_gradient(coords1);

          double dotProduct01 = grad_phi0[0] * grad_phi1[0] + grad_phi0[1] * grad_phi1[1] +
                                grad_phi0[2] * grad_phi1[2];

          Core::LinAlg::Matrix<3, 1> coords2;
          edge_points[2]->Coordinates(&coords2(0, 0));
          std::vector<double> grad_phi2 = element->get_level_set_gradient(coords2);

          Core::LinAlg::Matrix<3, 1> coords3;
          edge_points[3]->Coordinates(&coords3(0, 0));
          std::vector<double> grad_phi3 = element->get_level_set_gradient(coords3);

          double dotProduct23 = grad_phi2[0] * grad_phi3[0] + grad_phi2[1] * grad_phi3[1] +
                                grad_phi2[2] * grad_phi3[2];

          double dotProduct;

          if (fabs(dotProduct01) > fabs(dotProduct23))
            dotProduct = dotProduct01;
          else
            dotProduct = dotProduct23;

          if (dotProduct > 0.0)
            connect01and23 = true;
          else
            connect01and23 = false;
#endif

          if (lsv(0) <= 0 and lsv(1) >= 0 and lsv(2) <= 0 and lsv(3) >= 0)
          {
#ifdef USE_PHIDERIV_FOR_CUT_DETERMINATION
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
            FOUR_C_THROW("illegal levelset pattern");
          }

          return false;
          break;
        }
        default:
          return false;
      }
      break;
    }  // case Core::FE::CellType::quad4:
    default:
      FOUR_C_THROW("Unsupported side shape! (shape = %d | %s )", side.Shape(),
          Core::FE::CellTypeToString(side.Shape()).c_str());
      break;
  }
  return false;
}

template class Core::Geo::Cut::LevelSetSide<2>;
template class Core::Geo::Cut::LevelSetSide<3>;

FOUR_C_NAMESPACE_CLOSE
