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

  /*! \brief Construction of boundingbox over the considered geometrical element
   *  for fast overlap detection */
  class BoundingBox
  {
  public:
    /// @name Different Create methods
    /// @{

    static BoundingBox * Create( );
    static BoundingBox * Create( Edge & edge );
    static BoundingBox * Create( Side & side );
    static BoundingBox * Create( VolumeCell & volcell );
    static BoundingBox * Create( VolumeCell & volcell, Element* elem1 );
    static BoundingBox * Create( Element & element );

    /// @}
  protected:
    /// empty constructor
    BoundingBox() : empty_(true) { };

    /// Initialize the BoundingBox by using information of an edge
    void Init(Edge & edge);

    /// Initialize the BoundingBox by using information of a side
    void Init(Side & side);

    /// Initialize the BoundingBox by using information of a volume cell
    void Init(VolumeCell & volcell);

    /** \brief initializes bounding box over the volumecell in local
     *  coordinates (rst) with respect to elem1 */
    void Init(VolumeCell & volcell, Element* elem1);

    /// Initialize the BoundingBox by using the information of an element
    void Init(Element & element);

  public:
    /// destructor
    virtual ~BoundingBox() { };

    /*!
     \brief Make the boundingbox larger to include this edge within the box
     */
    void Assign(Edge & edge);

    /*! \brief Make the boundingbox larger to include this side within the
     *  box */
    void Assign(Side & side);

    /*! \brief Make the boundingbox larger to include this element
     * within the box */
    void Assign(Element& element);

    double operator()(int i, int j) const
    {
      return box_(i, j);
    }

    /*!  \brief If necessary make the boundingbox larger to include
     * this point as one of the corners of the box */
    virtual void AddPoint(const double * x) = 0;

    /*! \brief If necessary make the boundingbox larger to include
     *  this point as one of the corners of the box */
    void AddPoint(const LINALG::Matrix<3, 1> & p)
    {
      AddPoint(p.A());
    }

    /*! \brief If necessary make the boundingbox larger to include
     *  all these nodes as one of the corners of the box */
    void AddPoints(const std::vector<Node*> & nodes);

    /*! \brief Check whether "b" is within this boundingbox */
    bool Within(double norm, const BoundingBox & b) const;

    /*! \brief Check the point is within this boundingbox */
    virtual bool Within(double norm, const double * x) const;

    /*! \brief Check these points are within this boundingbox */
    bool Within(double norm, const Epetra_SerialDenseMatrix & xyz) const;

    /*! \brief Check the element is within this boundingbox */
    bool Within(double norm, Element & element) const;

    /*! \brief Print the corner points of boundingbox on the screen */
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

    /*! \brief Get the out-most point of the boundingbox */
    void CornerPoint(int i, double * x);

    /*! \brief Get the boundingbox */
    const LINALG::Matrix<3, 2>& GetBoundingVolume() const
    {
      return box_;
    }

  protected:

    bool InBetween(const double& norm,
        const double& smin,
        const double& smax,
        const double& omin,
        const double& omax) const;

    bool empty_;
    LINALG::Matrix<3, 2> box_;
  };  // class BoundingBox

  /*--------------------------------------------------------------------------*/
  /** Simple modification of the BoundingBox implementation, such that
   *  a correct and deterministic evaluation is ensured for all problem
   *  dimensions.
   *
   *  \author hiermeier \date 11/16 */
  template <unsigned probdim >
  class ConcreteBoundingBox : public BoundingBox
  {
  public:
    /// empty constructor
    ConcreteBoundingBox() : BoundingBox() { };

    /// destructor
    virtual ~ConcreteBoundingBox() {};

    /** add a new point under consideration of the problem dimension. If a
     *  2-dimensional problem is observed, we ignore the z-coordinate and keep
     *  a value of zero. */
    virtual void AddPoint( const double * x );

    /** \brief [derived]
     *
     *  special treatment of the 2-D case if necessary. */
    virtual bool Within( double norm, const double * x ) const;
  };  // class ConcreteBoundingBox

}
}

#endif
