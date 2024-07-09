/*---------------------------------------------------------------------*/
/*! \file

\brief bounding box for cut

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_boundingbox.hpp"

#include "4C_cut_element.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create()
{
  static const int probdim = Global::Problem::instance()->n_dim();
  switch (probdim)
  {
    case 2:
      return new ConcreteBoundingBox<2>();
    case 3:
      return new ConcreteBoundingBox<3>();
    default:
      FOUR_C_THROW("Unsupported problem dimension!");
  }
}


Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create(Node& node)
{
  BoundingBox* box = create();
  box->init(node);
  return box;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create(Edge& edge)
{
  BoundingBox* box = create();
  box->init(edge);
  return box;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create(Side& side)
{
  BoundingBox* box = create();
  box->init(side);
  return box;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create(VolumeCell& volcell)
{
  BoundingBox* box = create();
  box->init(volcell);
  return box;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create(
    VolumeCell& volcell, Element* elem1)
{
  BoundingBox* box = create();
  box->init(volcell, elem1);
  return box;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::Cut::BoundingBox* Core::Geo::Cut::BoundingBox::create(Element& element)
{
  BoundingBox* box = create();
  box->init(element);
  return box;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::init(Node& node)
{
  double x[3];
  node.coordinates(x);
  add_point(x);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::init(Edge& edge)
{
  const std::vector<Node*>& nodes = edge.nodes();
  add_points(nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::init(Side& side)
{
  const std::vector<Node*>& nodes = side.nodes();
  add_points(nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::init(VolumeCell& volcell)
{
  const plain_facet_set& facete = volcell.facets();
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fac = *i;
    const std::vector<Point*>& corners = fac->corner_points();
    for (std::vector<Point*>::const_iterator k = corners.begin(); k != corners.end(); k++)
    {
      const Point* po = *k;
      const double* coords = po->x();
      add_point(coords);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::init(VolumeCell& volcell, Element* elem1)
{
  const plain_facet_set& facete = volcell.facets();
  double x[3];
  for (plain_facet_set::const_iterator i = facete.begin(); i != facete.end(); i++)
  {
    Facet* fac = *i;
    std::vector<std::vector<double>> corLocal;
    fac->corner_points_local(elem1, corLocal);
    for (std::vector<std::vector<double>>::const_iterator m = corLocal.begin(); m != corLocal.end();
         m++)
    {
      std::vector<double> loc = *m;
      for (int j = 0; j < 3; j++) x[j] = loc[j];
      add_point(x);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::init(Element& element)
{
  const std::vector<Node*>& nodes = element.nodes();
  add_points(nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::assign(Side& side)
{
  empty_ = true;
  const std::vector<Node*>& nodes = side.nodes();
  add_points(nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::assign(Edge& edge)
{
  empty_ = true;
  const std::vector<Node*>& nodes = edge.nodes();
  add_points(nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::assign(Element& element)
{
  empty_ = true;
  const std::vector<Node*>& nodes = element.nodes();
  add_points(nodes);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::add_points(const std::vector<Node*>& nodes)
{
  for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
  {
    Node* n = *i;
    double x[3];
    n->coordinates(x);
    add_point(x);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::BoundingBox::within(double norm, const BoundingBox& b) const
{
  if (empty_) return true;
  return (in_between(norm, minx(), maxx(), b.minx(), b.maxx()) and
          in_between(norm, miny(), maxy(), b.miny(), b.maxy()) and
          in_between(norm, minz(), maxz(), b.minz(), b.maxz()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::BoundingBox::within(double norm, const double* x) const
{
  if (empty_) return true;
  return (in_between(norm, minx(), maxx(), x[0], x[0]) and
          in_between(norm, miny(), maxy(), x[1], x[1]) and
          in_between(norm, minz(), maxz(), x[2], x[2]));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::BoundingBox::within(
    double norm, const Core::LinAlg::SerialDenseMatrix& xyz) const
{
  Teuchos::RCP<BoundingBox> bb = Teuchos::rcp(create());
  int numnode = xyz.numCols();
  for (int i = 0; i < numnode; ++i)
  {
    bb->add_point(&xyz(0, i));
  }
  return within(norm, *bb);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::BoundingBox::within(double norm, Element& element) const
{
  Teuchos::RCP<BoundingBox> bb = Teuchos::rcp(create(element));
  return within(norm, *bb);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool Core::Geo::Cut::BoundingBox::in_between(const double& norm, const double& smin,
    const double& smax, const double& omin, const double& omax) const
{
  double tol = BOXOVERLAP * norm;
  return ((omax > smin - tol) and (smax > omin - tol));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::print()
{
  if (empty_)
  {
    std::cout << "  BB: {}\n";
  }
  else
  {
    std::cout << "  BB: {(" << box_(0, 0) << "," << box_(1, 0) << "," << box_(2, 0) << ")-("
              << box_(0, 1) << "," << box_(1, 1) << "," << box_(2, 1) << ")}\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::BoundingBox::corner_point(int i, double* x)
{
  x[0] = ((i & 1) == 1) ? maxx() : minx();
  x[1] = ((i & 2) == 2) ? maxy() : miny();
  x[2] = ((i & 4) == 4) ? maxz() : minz();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim>
void Core::Geo::Cut::ConcreteBoundingBox<probdim>::add_point(const double* x)
{
  if (empty_)
  {
    empty_ = false;
    box_(0, 0) = box_(0, 1) = x[0];
    box_(1, 0) = box_(1, 1) = x[1];
    if (probdim > 2) box_(2, 0) = box_(2, 1) = x[2];
  }
  else
  {
    box_(0, 0) = std::min(box_(0, 0), x[0]);
    box_(1, 0) = std::min(box_(1, 0), x[1]);
    // special treatment of z-coordinate
    if (probdim > 2) box_(2, 0) = std::min(box_(2, 0), x[2]);

    box_(0, 1) = std::max(box_(0, 1), x[0]);
    box_(1, 1) = std::max(box_(1, 1), x[1]);
    // special treatment of z-coordinate
    if (probdim > 2) box_(2, 1) = std::max(box_(2, 1), x[2]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim>
bool Core::Geo::Cut::ConcreteBoundingBox<probdim>::within(double norm, const double* x) const
{
  if (probdim > 2) return BoundingBox::within(norm, x);

  if (empty_) return true;
  return (in_between(norm, minx(), maxx(), x[0], x[0]) and
          in_between(norm, miny(), maxy(), x[1], x[1]));
}


template class Core::Geo::Cut::ConcreteBoundingBox<2>;
template class Core::Geo::Cut::ConcreteBoundingBox<3>;

FOUR_C_NAMESPACE_CLOSE
