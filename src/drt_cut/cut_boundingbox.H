/*---------------------------------------------------------------------*/
/*!
\file cut_boundingbox.H

\brief bounding box for cut

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_BOUNDINGBOX_H
#define CUT_BOUNDINGBOX_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "cut_tolerance.H"

// only forward declarations are not working when
// building a release version of BACI with the new
// Trilinos. The cut module in BACI is a very weird
// piece of code and needs some major review of the
// code organization...
#include "cut_node.H"  // compilation fix
#include "cut_edge.H"  // compilation fix
namespace GEO
{
namespace CUT
{
class Node;
class Edge;
class Side;
class Facet;
class VolumeCell;
class Element;

/*!
 \brief Construction of boundingbox over the considered geometrical element for fast overlap detection
 */
class BoundingBox
{
public:

  BoundingBox() :
      empty_(true)
  {
  }

  explicit BoundingBox(Edge & edge);

  explicit BoundingBox(Side & side);

  explicit BoundingBox(VolumeCell & volcell);

  explicit BoundingBox(VolumeCell & volcell, Element* elem1);

  explicit BoundingBox(Element & element);

  /*!
   \brief Make the boundingbox larger to include this edge within the box
   */
  void Assign(Edge & edge);

  /*!
   \brief Make the boundingbox larger to include this side within the box
   */
  void Assign(Side & side);

  /*!
   \brief Make the boundingbox larger to include this element within the box
   */
  void Assign(Element & element);

  double operator()(int i, int j) const
  {
    return box_(i, j);
  }

  /*!
   \brief If necessary make the boundingbox larger to include this point as one of the corners of the box
   */
  void AddPoint(const double * x);

  /*!
   \brief If necessary make the boundingbox larger to include this point as one of the corners of the box
   */
  void AddPoint(const LINALG::Matrix<3, 1> & p)
  {
    AddPoint(p.A());
  }

  /*!
   \brief If necessary make the boundingbox larger to include all these nodes as one of the corners of the box
   */
  void AddPoints(const std::vector<Node*> & nodes);

  /*!
   \brief Check whether "b" is within this boundingbox
   */
  bool Within(double norm, const BoundingBox & b) const;

  /*!
   \brief Check the point is within this boundingbox
   */
  bool Within(double norm, const double * x) const;

  /*!
   \brief Check these points are within this boundingbox
   */
  bool Within(double norm, const Epetra_SerialDenseMatrix & xyz) const;

  /*!
   \brief Check the element is within this boundingbox
   */
  bool Within(double norm, Element & element) const;

  /*!
   \brief Print the corner points of boundingbox on the screen
   */
  void Print();

  double minx() const
  {
    return box_(0, 0);
  }
  double miny() const
  {
    return box_(1, 0);
  }
  double minz() const
  {
    return box_(2, 0);
  }

  double maxx() const
  {
    return box_(0, 1);
  }
  double maxy() const
  {
    return box_(1, 1);
  }
  double maxz() const
  {
    return box_(2, 1);
  }

  /*!
   \brief Get the outmost point of the boundingbox
   */
  void CornerPoint(int i, double * x);

  /*!
   \brief Get the boundingbox
   */
  LINALG::Matrix<3, 2> GetBoundingVolume()
  {
    return box_;
  }

private:

  bool InBetween(double norm, double smin, double smax, double omin,
      double omax) const
  {
    double tol = BOXOVERLAP * norm;
    return ((omax > smin - tol) and (smax > omin - tol));
  }

  bool empty_;
  LINALG::Matrix<3, 2> box_;
};

}
}

#endif
