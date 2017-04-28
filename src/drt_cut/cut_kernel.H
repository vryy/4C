/*---------------------------------------------------------------------*/
/*!
\file cut_kernel.H

\brief The cut kernel computes basic geometric operation, implemented are
    - Intersection of Surface and line or line and line
    - Calculate local coordinates inside an element
    - Compute Distance from a point to an embedded geometrical object
      ( surface or line )

\level 1

\maintainer Michael Hiermeier

*/
/*----------------------------------------------------------------------*/

#ifndef CUT_KERNEL_H
#define CUT_KERNEL_H

#include "cut_tolerance.H"
#include "cut_utils.H"
#include "cut_enum.H"

#include "../linalg/linalg_gauss.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

//activate debug intersection output in cut_kernel
//#define DEBUG_CUTKERNEL_OUTPUT

#ifdef DEBUG_CUTKERNEL_OUTPUT
  #include "cut_output.H"
#endif

namespace GEO {
namespace CUT {
/*! \brief Core geometrical operations like finding the cut points,
 *   distance calculations takes place here */
namespace KERNEL {

  void FindCornerPoints(
      const std::vector<Point*> & points,
      std::vector<Point*> & corner_points );

  unsigned FindNextCornerPoint(
      const std::vector<Point*> & points,
      LINALG::Matrix<3,1> & x1,
      LINALG::Matrix<3,1> & x2,
      LINALG::Matrix<3,1> & x3,
      LINALG::Matrix<3,1> & b1,
      LINALG::Matrix<3,1> & b2,
      LINALG::Matrix<3,1> & b3,
      unsigned i );

  /// check if the given corner points can belong to a point1 "element"
  bool IsValidPoint1( const std::vector<Point*> & corner_points );

  /// check if the given corner points can belong to a line2 element
  bool IsValidLine2( const std::vector<Point*> & corner_points );

  /// check if the given corner points can belong to a tri3 element
  inline bool IsValidTri3( const std::vector<Point*> & points )
  { return points.size()==3; }

  /// check if the given corner points can belong to a quad4 element
  bool IsValidQuad4( const std::vector<Point*> & points );

  DRT::Element::DiscretizationType CalculateShape(
      const std::vector<Point*> & points,
      std::vector<Point*> & line_points );

  template <int points>
  double FindL2Scaling( LINALG::Matrix<3,points> & xyze )
  {
    double scale = 0;
    LINALG::Matrix<3,1> d;
    for ( int i=0; i<points; ++i )
    {
      LINALG::Matrix<3,1> x1( &xyze( 0, i ), true );
      LINALG::Matrix<3,1> x2( &xyze( 0, ( i+1 )%points ), true );
      d.Update( 1, x2, -1, x1, 0 );
      scale += d.Norm2();
    }
    scale /= points;
    return scale;
  }

  /** \brief check if the given local point \c xsi of the side is on an edge
   *
   *  The function returns \TRUE if the given point lies on one of the edges
   *  of the side.
   *
   *  \param xsi (in) : local parameter space coordinates
   *  \param tol (in) : tolerance for the check (default: \c tol=REFERENCETOL) */
  template < DRT::Element::DiscretizationType elementType, class T >
  bool AtEdge( const T & xsi, const double & tol )
  {
    // sanity check
    if (xsi.Rows() < DRT::UTILS::DisTypeToDim<elementType>::dim)
      dserror("The given local coordinate has the wrong dimension!\n"
          "xsi.Rows() < Dim <===> %d<%d",xsi.Rows(),
          DRT::UTILS::DisTypeToDim<elementType>::dim);

    switch ( elementType )
    {
      case DRT::Element::line2:
      {
        return ( std::abs(xsi( 0 )+1)<tol or std::abs(xsi( 0 )-1)<tol );
        break;
      }
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        return ( std::abs( xsi( 0 )+1 )<tol or std::abs( xsi( 1 )+1 )<tol or
                 std::abs( xsi( 0 )-1 )<tol or std::abs( xsi( 1 )-1 )<tol );
        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        return ( std::abs( xsi( 0 )+0 )<tol or std::abs( xsi( 1 )+0 )<tol or
                 std::abs( xsi( 1 )+xsi( 0 )-1 )<tol );
        break;
      }
      default:
        dserror( "unsupported element type: %i | %s",elementType,
            DRT::DistypeToString(elementType).c_str() );
        break;
    }
  }
  template < DRT::Element::DiscretizationType elementType, class T >
  bool AtEdge( const T & xsi )
  {
    return AtEdge<elementType>(xsi,REFERENCETOL);
  }

  /** \brief check if the local coordinate is within the specified limits
   *
   *  Returns \TRUE if the given point lies inside the element.
   *
   *  \param xsi (in) : local parameter space coordinates
   *  \param tol (in) : tolerance for the check (default: \c tol=REFERENCETOL) */
  template < DRT::Element::DiscretizationType elementType, class T >
  bool WithinLimits( const T & xsi, const double & tol)
  {
    // sanity check
    if (xsi.M() < DRT::UTILS::DisTypeToDim<elementType>::dim)
      dserror("The given local coordinate has the wrong dimension!\n"
          "xsi.Rows() < eleDim <===> %d<%d",xsi.M(),
          DRT::UTILS::DisTypeToDim<elementType>::dim);

    switch ( elementType )
    {
      // ---------------------------------------------------------------
      // 1-D elements
      // ---------------------------------------------------------------
      case DRT::Element::line2:
      {
        return (xsi( 0 ) >= -1-tol and xsi( 0 ) <= 1+tol);
        break;
      }
      // ---------------------------------------------------------------
      // 2-D elements
      // ---------------------------------------------------------------
      case DRT::Element::quad4:
      case DRT::Element::quad8:
      case DRT::Element::quad9:
      {
        return xsi( 0 ) >= -1-tol and xsi( 1 ) >= -1-tol and
               xsi( 0 ) <=  1+tol and xsi( 1 ) <=  1+tol ;
        break;
      }
      case DRT::Element::tri3:
      case DRT::Element::tri6:
      {
        return xsi( 0 ) >= 0-tol and xsi( 1 ) >= 0-tol and
               xsi( 0 ) <= 1+tol and xsi( 1 ) <= 1+tol and
               xsi( 1 ) <= ( 1 - xsi( 0 ) )+tol ;
        break;
      }
      // ---------------------------------------------------------------
      // 3-D elements
      // ---------------------------------------------------------------
      case DRT::Element::hex8:
      case DRT::Element::hex16:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        return ( xsi( 0 ) >= -1-tol and xsi( 1 ) >= -1-tol and xsi( 2 ) >= -1-tol and
                 xsi( 0 ) <=  1+tol and xsi( 1 ) <=  1+tol and xsi( 2 ) <=  1+tol );
        break;
      }
      case DRT::Element::tet4:
      case DRT::Element::tet10:
      {
        return ( xsi( 0 ) >=   -tol and xsi( 1 ) >=   -tol and xsi( 2 ) >=   -tol and
                 xsi( 0 )+xsi( 1 )+xsi( 2 ) <= 1+tol );
        break;
      }
      case DRT::Element::pyramid5:
      {
        return ( xsi( 0 ) >= -1-tol and xsi( 1 ) >= -1-tol and xsi( 2 ) >=   -tol and
                 xsi( 0 ) <=  1+tol and xsi( 1 ) <=  1+tol and
                 ( std::abs( xsi( 0 ) ) > std::abs( xsi( 1 ) ) ?
                   ( xsi( 2 )+std::abs( xsi( 0 ) ) <= 1+tol ) :
                   ( xsi( 2 )+std::abs( xsi( 1 ) ) <= 1+tol ) ));
        break;
      }
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
      {
        return ( xsi( 0 ) >=   -tol and xsi( 1 ) >=   -tol and xsi( 2 ) >= -1-tol and
                 xsi( 2 ) <=  1+tol and
                 xsi( 0 )+xsi( 1 ) <= 1+tol );
        break;
      }
      default:
        dserror( "unsupported element type: %i | %s",elementType,
            DRT::DistypeToString(elementType).c_str() );
        break;
    } // switch ( elementType )
    return false;
  }
  template < DRT::Element::DiscretizationType elementType, class T >
  bool WithinLimits( const T & xsi)
  {
    return WithinLimits<elementType>(xsi,REFERENCETOL);
  }

  /** \brief check if the local augmented coordinates are within the specified limits
   *
   *  \author hiermeier */
#if __cplusplus < 201103L

  template <unsigned probDim,
              DRT::Element::DiscretizationType sideType>
    bool WithinLimitsEmbeddedManifold( const Epetra_SerialDenseMatrix & xsi_aug,
        double tol, bool allow_dist, double tol2)
    {
      const LINALG::Matrix<probDim,1> xsi_aug_linalg( xsi_aug.A(), true );
      return WithinLimitsEmbeddedManifold<probDim,sideType>( xsi_aug_linalg, tol, allow_dist, tol2 );
    }

  template <unsigned probDim,
            DRT::Element::DiscretizationType sideType>
  bool WithinLimitsEmbeddedManifold( const LINALG::Matrix<probDim,1> & xsi_aug,
      double tol, bool allow_dist, double tol2)
  {
    const unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim;
#else
  template <unsigned probDim,
            DRT::Element::DiscretizationType sideType,
            unsigned dimSide = DRT::UTILS::DisTypeToDim<sideType>::dim,
            class T>
  bool WithinLimitsEmbeddedManifold( const T & xsi_aug,
      double tol, bool allow_dist, double tol2)
  {
#endif
    if (probDim==1)
      dserror("probDim is not allowed to be equal 1!");
    if (xsi_aug.M()!=probDim)
      dserror("Wrong dimension of the augmented xsi input variable. \n"
          "The dimension is supposed to be equal to the problem dimension. \n"
          "The last coordinates are holding the distances from the global \n"
          "point to the side element. ");

    bool check = true;
    // check the distances if desired
    if (not allow_dist)
      for ( unsigned r = dimSide; r<probDim; ++r )
        check = ( check and (xsi_aug(r,0) >= -tol2 and xsi_aug(r,0) <= tol2) ? true : false );

    /* Since the first dimSide vector entries are the important ones,
     * we do not have to perform any extraction. */
    return ( check and WithinLimits<sideType>( xsi_aug, tol ) );
  }
  template <unsigned probDim,
            DRT::Element::DiscretizationType sideType,
            class T>
  bool WithinLimitsEmbeddedManifold( const T & xsi_aug,
      const bool& allow_dist, const double & tol2)
  {
    return KERNEL::WithinLimitsEmbeddedManifold<probDim,sideType>(xsi_aug,
        REFERENCETOL,allow_dist,tol2);
  }

  /*! \brief Check whether three points are on the same line
   *  by checking whether the cross product is zero
   *
   *  \author Sudhakar
   *  \date 04/12 */
  bool IsOnLine( Point*& pt1, Point*& pt2, Point*& pt3, bool DeleteInlinePts = false );

  /*! \brief Check whether the polygon defined by the set of points is convex */
  std::vector<int> CheckConvexity( const std::vector<Point*>& ptlist,
                                   GEO::CUT::FacetShape& geomType,
                                   bool InSplit = true,
                                   bool DeleteInlinePts = false );

  //std::vector<double> EqnPlanePolygon( const std::vector<Point*>& ptlist, bool DeleteInlinePts = false );


  /*! \brief Compute the equation of plane of this polygon using Newell's method */
  std::vector<double> EqnPlaneOfPolygon( const std::vector<Point*>& ptlist );

  /*! \brief Compute the equation of plane of this polygon using Newell's method */
  std::vector<double> EqnPlaneOfPolygon( const std::vector<std::vector<double> >& vertices );

  /*! \brief Find the equation of plane that contains these non-collinear points
   *
   *  It must be noted while using this function to find equation of facets,
   *  none of these 3 points must be a reflex (concave) point */
  std::vector<double> EqnPlane( Point* & pt1, Point* & pt2, Point* & pt3 );

  /*! \brief Check whether the point named \c check is inside the triangle
   *  formed by tri */
  bool PtInsideTriangle( std::vector<Point*> tri, Point* check, bool DeleteInlinePts = false );

  /*! \brief Check whether the point \c check is inside the Quad */
  bool PtInsideQuad( std::vector<Point*> quad, Point* check );

  /*! \brief Return \TRUE if the points of the polygon are ordered clockwise
   *
   *  Polygon in 3D space is first projected into 2D plane, and the plane
   *  of projection is returned in projType.
   *
   *  \author sudhakar
   *  \date 05/12 */
  bool IsClockwiseOrderedPolygon( std::vector<Point*>polyPoints, std::string& projType );

  /*! \brief Delete unnecessary in-line points.
   *
   *  If more than two points are on a line, all points except the end points
   *  are deleted. This is checked for all the lines for a facet. So once this is
   *  called the facet is free of more than 2 inline points
   *
   *  \author sudhakar
   *  \date 06/12 */
  void DeleteInlinePts( std::vector<Point*>& poly );

  /*! \brief Returns true if at least 3 points are collinear
   *
   * \author wirtz
   * \date 05/13 */
  bool HaveInlinePts( std::vector<Point*>& poly );

  /*! \brief Finds tree points of the polygon which are not collinear
   *
   *  \author wirtz
   *  \date 05/13 */
  std::vector<Point*> Get3NoncollinearPts( std::vector<Point*> & polyPoints );

  /*! \brief Find appropriate projection plane
   *
   *  In several cases, it is appropriate to project the surface in 3D space into appropriate
   *  coordinate plane. It is better to project over the plane which is max normal component
   *  because this will reduce the round-off error in further calculations.
   *
   *  \author sudhakar
   *  \date 06/12 */
  void FindProjectionPlane( std::string& projPlane, const std::vector<double>& eqn );

  /*! \brief Split a QUAD4 element into two TRI3 elements
   *
   *  Here it is important that the triangle is created in the same rotation
   *  as the QUAD4 is, to get normal in the same direction and therefore the
   *  same signed distance:
   *
   *  tri3_id = 0 ---> QUAD4 nodes = {0 1 2}
   *  tri3_id = 1 ---> QUAD4 nodes = {2 3 0}
   *
   *  \param xyze_quad4 (in) : nodal coordinates of a QUAD4 element
   *  \param tri3_id    (in) : id of the desired TRI3 element \f& id \in \{0,\; 1\}\f$
   *  \param xyze_tri3  (out): created/filled TRI3 nodal coordinate matrix
   *
   *  \author hiermeier \date 11/16 */
  template < class T1, class T2 >
  void SplitQuad4IntoTri3(const T1 & xyze_quad4, const unsigned & tri3_id, T2 & xyze_tri3)
  {
     if (tri3_id > 1)
       dserror("A QUAD4 is supposed to be split into 2 TRI3 elements, \n"
           "therefore you have the choice between the tri3_id's 0 and 1. \n"
           "Nevertheless, you tried to access the id %d.",tri3_id);

     unsigned n0 = 2 * tri3_id;
     for (unsigned r = 0; r < static_cast<unsigned>( xyze_quad4.M() ); ++r)
       for (unsigned c = 0; c < 3 ; ++c)
         xyze_tri3( r, c ) = xyze_quad4( r,( (c+n0) % 4 ) );
  };

  /*! \brief Get area of triangle in 3D space
   *
   *  \author sudhakar
   *  \date 11/14 */
  double getAreaTri(const std::vector<Point*> & poly,
                    LINALG::Matrix<3,1>* normalvec = NULL);
  double getAreaTri( const double * a, const double * b, const double * c,
      LINALG::Matrix<3,1>* normalvec = NULL );
  template < class T >
  double getAreaTri( const T & xyze,
      LINALG::Matrix<3,1>* normalvec = NULL )
  {
    if ( xyze.M() != 3 )
      run_time_error( "Currently unsupported element dimension!" );
    if ( xyze.N() != 3 )
      run_time_error( "Wrong node number!" );

    return getAreaTri( & xyze(0,0), & xyze(0,1), & xyze(0,2), normalvec );
  }

  /*! \brief Get area of convex Quad in 3D space
   *
   *  Quad is split into two triangles and area of tri are summed up
   *
   *  \param poly (in): vector containting four pointers to the quad points
   *
   *  \author sudhakar
   *  \date 11/14 */
  double getAreaConvexQuad( std::vector<Point*> & poly );

  enum NewtonStatus
  {
    converged = 1,
    unconverged = 0,
    failed = -1
  };

  /** \brief Build an adaptive combined Newton tolerance
   *
   *  First a scaled absolute tolerance is calculated, based on the
   *  inf-norm of the absolute coordinate values of the considered side
   *  element and point.
   *
   *  In a second step a relative tolerance is added, which accounts for the
   *  fact, that the point px_ is maybe far outside of the considered side.
   *  If the point is inside or near the side, the additional term will stay
   *  very low and will not distort the used tolerance too much.
   *
   *  \author hiermeier \date 02/17 */
  template < class T1, class T2, class T3 >
  double AdaptiveCombinedNewtonTolerance(
      const T1 & xyze,
      const T2 & px,
      const T3 & initial_rhs )
  {
    /* --- Build the absolute tolerance */
    double tol = xyze.NormInf();
    double linescale = px.NormInf();
    if (linescale > tol )
      tol = linescale;
    tol *= LINSOLVETOL;

    /* --- Add the relative tolerance */
    tol += LINSOLVETOL * initial_rhs.NormInf();

    return tol;
  }

  /// Convert Newton Status enumerator to string
  static inline std::string NewtonStatus2String( const enum NewtonStatus & status )
  {
    switch ( status )
    {
      case converged:
        return "CONVERGED";
      case unconverged:
        return "UNCONVERGED";
      case failed:
        return "FAILED";
      default:
        dserror("Unknown Newton status!");
        exit(EXIT_FAILURE);
    }
    exit(EXIT_FAILURE);
  }

  /*--------------------------------------------------------------------------*/
  /// generic Newton algorithm
  template <class Strategy,
            unsigned dim=3,
            unsigned maxiter=20>
  class NewtonSolve : public Strategy
  {
  public:
    /// constructor
    NewtonSolve( LINALG::Matrix<dim,1> & xsi,
                 bool checklimits )
        : Strategy( xsi, checklimits )
    {
    }

    bool Solve()
    {
      this->SetupSolve();

      for ( unsigned iter = 0; iter < maxiter; ++iter )
      {
        this->SetupStep( iter );

        switch ( this->TestConverged( iter ) )
        {
          // Newton iteration was successful
          case converged:
          {
            return true;
          }
          // Newton did not yet converge
          case unconverged:
          {
            break;
          }
          // Newton failed, thus we can stop here
          case failed:
          {
            return this->NewtonFailed();
          }
        }

        if ( not this->LinearSolve( iter ) )
        {
          return false;
        }

        if ( not this->Update( iter ) )
        {
          return false;
        }
      }

      return this->NewtonFailed();
    }

  };  // class class NewtonSolve

  /*--------------------------------------------------------------------------*/
  /// empty strategy for generic Newton algorithm
  class EmptyNewtonStrategy
  {
  public:

    void SetupSolve() {}

    void SetupStep( int iter ) {}

    enum NewtonStatus TestConverged( int iter )
    {
      return unconverged;
    }

    bool LinearSolve( int iter )
    {
      return true;
    }

    bool Update( int iter )
    {
      return true;
    }

    bool NewtonFailed()
    {
      return false;
    }

    double GetTolerance()
    {
      dserror("Try to get tolerance from EmptyNewtonStrategy!");
      return 0.0; //just to make compiler happy!
    }

    bool ZeroArea()
    {
      dserror("Try to get ZeroArea Information from EmptyNewtonStrategy!");
      return false; //just to make compiler happy!
    }

    void WritetoGmsh(std::ofstream& file)
    {
      dserror("Try WritetoGmsh() from EmptyNewtonStrategy!");
    }
  }; // class EmptyNewtonStrategy

  /*--------------------------------------------------------------------------*/
  /// Debug helper strategy for generic Newton algorithm
  template <class Strategy,
            unsigned dim=3>
  class DebugNewtonStrategy : public Strategy
  {
  public:
    /// constructor
    DebugNewtonStrategy( LINALG::Matrix<dim,1> & xsi,
                         bool checklimits )
        : Strategy( xsi, checklimits )
    {
    }

    void SetupSolve()
    {
      std::cout << "SetupSolve()\n";
      Strategy::SetupSolve();
    }

    void SetupStep( int iter )
    {
      std::cout << "SetupStep( iter = " << std::setw(2) << iter << " )\n";
      Strategy::SetupStep( iter );
    }

    enum NewtonStatus TestConverged( int iter )
    {
      enum NewtonStatus status = Strategy::TestConverged( iter );
      std::cout << "TestConverged( iter = " << std::setw(2) << iter << " ) = "
                << NewtonStatus2String( status ) << "\n" << std::flush;
      return status;
    }

    bool LinearSolve( int iter )
    {
      bool res = Strategy::LinearSolve( iter );
      std::cout << "LinearSolve( iter = " << std::setw(2) << iter << " )   = "
                << (res ? "SUCCESS" : "FAILED") << "\n" << std::flush;
      return res;
    }

    bool Update( int iter )
    {
      bool res = Strategy::Update( iter );
      std::cout << "Update( iter = " << std::setw(2) << iter << " )        = "
                << (res ? "SUCCESS" : "FAILED") << "\n" << std::flush;
      return res;
    }

    bool NewtonFailed()
    {
      bool res = Strategy::NewtonFailed();
      std::cout << "NewtonFailed() = "
                << (res ? "SUCCESS" : "FAILED") << "\n" << std::flush;
      return res;
    }
  }; // class DebugNewtonStrategy

  /*--------------------------------------------------------------------------*/
  /** \brief strategy for position of point within element
   *
   *  inheritance diagram:
   *
   *  ComputePosition --> GenericComputePosition --> NewtonSolve
   *  --> ComputePositionStrategy --> EmptyNewtonStrategy
   *      ^^^^^^^^^^^^^^^^^^^^^^^
   */
  template <bool debug,
            unsigned probDim,
            DRT::Element::DiscretizationType elementType,
            unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
            unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
  class ComputePositionStrategy : public EmptyNewtonStrategy
  {
  public:
    /// constructor
    ComputePositionStrategy(
        LINALG::Matrix<dim,1> & xsi,
        bool checklimits )
        : xsi_( xsi ),
          xyze_( NULL ),
          px_( NULL ),
          tol_( 0.0 )
    {
    }

    /** setup the computation
     *
     *  \param xyze (in) : global nodal coordinates
     *  \param px   (in) : global position vector of the searched parameter space point */
    void Setup( const LINALG::Matrix<probDim, numNodesElement> & xyze,
                const LINALG::Matrix<probDim,1> & px )
    {
      xyze_ = &xyze;
      px_ = &px;

      /* use tol_ first just to store the scaling for the tolerance ..
       * set it to the maximal absolute coordinate value of the involved element
       * or point!
       *
       * --> this defines accuracy of the intersection (for good conditioned
       *     linear systems) */
      tol_ = xyze_->NormInf();
      double linescale = px_->NormInf();

      if (linescale > tol_)
        tol_ = linescale;

      tol_ *= LINSOLVETOL;

      if ( debug )
      {
        std::cout << "\n=== ComputePosition ===\n";
        std::cout << "Element    = " << DRT::DistypeToString( elementType ) << "\n";
        std::cout << "xyze_      = " << *xyze_;
        std::cout << "px_        = " << *px_;
        std::cout << "\n";
      }
    }

    /// get the local solution coordinates
    const LINALG::Matrix<dim,1> & LocalCoordinates() { return xsi_; }

    /// initialize the solution variable to zero
    void SetupSolve()
    {
      xsi_ = 0;

      // calculate initial right hand side
      PositionRHS( *xyze_, *px_, b_ );

      tol_ = AdaptiveCombinedNewtonTolerance(*xyze_, *px_, b_ );
    }

    /** \brief evaluate the current right hand side
     *
     *  \f[
     *      b_ = - (px_ - x(xsi_))
     *  \f] */
    void SetupStep( int iter )
    {
      if ( iter > 0)
      {
        PositionRHS( *xyze_, *px_, b_ );
      }
    }

    /// calculate the position rhs
    void PositionRHS(
        const LINALG::Matrix<probDim,numNodesElement> & xyze,
        const LINALG::Matrix<probDim,1> & px,
        LINALG::Matrix<probDim, 1> & b
        )
    {
      DRT::UTILS::shape_function<elementType>( xsi_, funct_ );
      b = px;
      b.Multiply( -1, xyze, funct_, 1 );
    }

    /// check for convergence
    enum NewtonStatus TestConverged( int iter )
    {
      double residual = b_.Norm2();

      if ( debug )
      {
        std::cout << "WithinLimits             = "
                  << (WithinLimits<elementType>( xsi_ ) ? "TRUE" : "FALSE") << "\n"
                  << "rsd_.Norm2()             = " << xsi_.Norm2() << "\n"
                  << "dx_.Norm2() / rsd.Norm2()= "
                  << ( iter != 0 ? ( dx_.Norm2() / xsi_.Norm2() ) : -1.0 ) << "\n"
                  << "b_.Norm2() ( RHS )       = " << b_.Norm2() << "\n"
                  << "b_.Norm2() / rsd.Norm2() = "
                  << (iter != 0 ? ( b_.Norm2() / xsi_.Norm2() ) : -1.0 ) << "\n"
                  << "tol_                     = " << tol_ << "\n"
                  << "rsd_ ( aka xsi_ )        = " << xsi_
                  << "b_ (RHS)                 = " << b_ << "\n" << std::flush;
      }

      return ( std::abs( residual ) < tol_ ? converged : unconverged );
    }

    /** \brief Build and solve the linear system
     *
     *  \param iter (in): Currently unused. Maximal number of linear iterations,
     *                    but we do a direct solve here.
     *
     * \f[
     * \underline{\underline{J}} \Delta \underline{\xi} = - (\underline{x}_{p} - \underline{x(\underline{\xi})},
     * \f]
     * where
     * \f[
     * \underline{\underline{J}} =
     * \begin{pmatrix}
     * \frac{\partial x(\underline{\xi})}{\partial \xi} & \frac{\partial x(\underline{\xi})}{\partial \eta} & \frac{\partial x(\underline{\xi})}{\partial \zeta} \\
     * \frac{\partial x(\underline{\xi})}{\partial \xi} & \frac{\partial y(\underline{\xi})}{\partial \eta} & \frac{\partial y(\underline{\xi})}{\partial \zeta} \\
     * \frac{\partial x(\underline{\xi})}{\partial \xi} & \frac{\partial z(\underline{\xi})}{\partial \eta} & \frac{\partial z(\underline{\xi})}{\partial \zeta}
     * \end{pmatrix}.
     * \f]
     *
     * For the special case that the element dimension is smaller than the
     * problem dimension (i.e. manifold), see the ComputeDistance classes. */
    bool LinearSolve( int iter )
    {
      double det = 0.0;
      // ToDo check if it is possible to switch completely to the manifold case,
      // without too much loss in performance.
      DRT::UTILS::shape_function_deriv1<elementType>( xsi_ , deriv1_ );
      A_.MultiplyNT( *xyze_, deriv1_ );

      dx_ = 0;
      det = LINALG::gaussElimination<true, probDim>( A_, b_, dx_ );

      if ( debug )
      {
        std::cout << "det                        = " << det << "\n";
        std::cout << "dx_ ( solution increment ) = " << dx_;
        std::cout << "dx_.Norm2()                = " << dx_.Norm2() << "\n";
      }

      return std::abs( det ) >= LINSOLVETOL;
    }

    /// Update the solution parameter space coordinates
    bool Update( int iter )
    {
      // update the solution variables (element dimension)
      for ( unsigned i=0; i<dim; ++i )
        xsi_(i,0) += dx_(i,0);

      return true;
    }

    bool NewtonFailed()
    {
#if 1
      return false;
#else
      std::stringstream str;
      str << "ComputePosition: Newton scheme did not converge:\n"
          << ( *xyze_ )
          << ( *px_ )
          << xsi_
        ;
      throw std::runtime_error( str.str() );
#endif
    }

    /// get the Newton tolerance
    double GetTolerance() const
    {
      return tol_;
    }

  private:

    /** parameter space coordinates (final result of the calculation)
     * corresponding to the global position px_. */
    LINALG::Matrix<dim,1> & xsi_;

    /// pointer to the global nodal positions
    const LINALG::Matrix<probDim, numNodesElement> * xyze_;
    /// global position vector (we are looking for the corresponding parameter coordinates)
    const LINALG::Matrix<probDim,1> * px_;

    /// nodal shape function values at the position xsi_
    static LINALG::Matrix<numNodesElement,1> funct_;
    /** nodal first derivative shape function values at the position xsi_ */
    static LINALG::Matrix<probDim,numNodesElement> deriv1_;

    /** \brief (extended) jacobian matrix
     *
     *  Keep in mind, that the jacobian has to be extended if dim < probDim! */
    static LINALG::Matrix<probDim,probDim> A_;
    /// right hand side vector b_ = -(px_ - x_(xsi_))
    static LINALG::Matrix<probDim,1> b_;
    /// newton increment (extended) parameter space coordinates
    static LINALG::Matrix<probDim,1> dx_;

    /// adapted Newton tolerance
    double tol_;
  }; // class ComputePositionStrategy

  /*--------------------------------------------------------------------------*/
  /** \brief generic strategy for position of point within element
   *
   *  inheritance diagram:
   *
   *  ComputePosition --> GenericComputePosition --> NewtonSolve
   *                      ^^^^^^^^^^^^^^^^^^^^^^
   *  --> ComputePositionStrategy --> EmptyNewtonStrategy
   */
  template <class Strategy,
            unsigned probDim,
            DRT::Element::DiscretizationType elementType,
            unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
            unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
  class GenericComputePosition : Strategy
  {
  public:

    GenericComputePosition( LINALG::Matrix<dim,1> & xsi )
        : Strategy( xsi, true )
    {
    }

    bool operator()( const LINALG::Matrix<probDim, numNodesElement> & xyze,
                     const LINALG::Matrix<probDim,1> & px )
    {
      this->Setup( xyze, px );
      if ( this->Solve() )
      {
        return true;
      }
      return false;
    }

    const LINALG::Matrix<dim,1> & LocalCoordinates()
    {
      return Strategy::LocalCoordinates();
    }

    double GetTolerance() const
    {
      return Strategy::GetTolerance();
    }

  }; // class GenericComputePosition

  /*--------------------------------------------------------------------------*/
  /** \brief most derived strategy class for position of point within element
   *
   *  inheritance diagram:
   *
   *  ComputePosition --> GenericComputePosition --> NewtonSolve
   *  ^^^^^^^^^^^^^^^
   *  --> ComputePositionStrategy --> EmptyNewtonStrategy
   */
  template <unsigned probDim,
            DRT::Element::DiscretizationType elementType,
            unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
            unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
  class ComputePosition
      : public GenericComputePosition<NewtonSolve<ComputePositionStrategy<false,probDim,elementType>,
                                                  dim >,
                                      probDim,
                                      elementType>
  {
  public:
    ComputePosition( LINALG::Matrix<dim,1> & xsi )
      : GenericComputePosition<NewtonSolve<ComputePositionStrategy<false,probDim,elementType,numNodesElement>,dim >,
                               probDim,
                               elementType>( xsi )
    {
      if (dim!=probDim)
        dserror("You called the wrong object! Use the ComputeDistance class "
            "for the embedded case.");
    }
  }; // class ComputePosition

  /*--------------------------------------------------------------------------*/
  template <unsigned probDim,
            DRT::Element::DiscretizationType elementType,
            unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementType>::numNodePerElement,
            unsigned dim = DRT::UTILS::DisTypeToDim<elementType>::dim>
  class DebugComputePosition
      : public GenericComputePosition<NewtonSolve<DebugNewtonStrategy<ComputePositionStrategy<true,probDim,elementType,numNodesElement>,dim >,dim >,
                                      probDim,
                                      elementType>
  {
  public:

    DebugComputePosition( LINALG::Matrix<dim,1> & xsi )
        : GenericComputePosition<NewtonSolve<DebugNewtonStrategy<ComputePositionStrategy<true,probDim,elementType,numNodesElement>,dim >,dim >,
                                 probDim,
                                 elementType>( xsi )
    {
    }
  }; // class DebugComputePosition

  /*--------------------------------------------------------------------------*/
  /** \brief strategy for distance between side and point
   *
   *  inheritance diagram:
   *
   *  ComputeDistance --> GenericComputeDistance --> NewtonSolve
   *  --> ComputeDistanceStrategy --> EmptyNewtonStrategy
   *      ^^^^^^^^^^^^^^^^^^^^^^^
   */
  template <bool debug,
            unsigned probDim,
            DRT::Element::DiscretizationType sideType,
            unsigned dimSide      = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class ComputeDistanceStrategy : public EmptyNewtonStrategy
  {
  public:
    /// constructor
    ComputeDistanceStrategy( LINALG::Matrix<probDim,1> & xsi, bool checklimits )
      : xsi_( xsi ),
        distance_( xsi.A()+dimSide ),
        xyze_side_( NULL ),
        px_( NULL ),
        tol_(0.0),
        zeroarea_(false)
    {
    }

    /** \brief Setup the distance calculation problem
     *
     *  \param xyze_side (in) : global coordinates of the considered side
     *  \param px        (in) : global coordinates of the considered point */
    void Setup( const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
                const LINALG::Matrix<probDim,1> & px )
    {
      xyze_side_ = &xyze_side;
      px_ = &px;

      if ( debug )
      {
        std::cout << "\n=== ComputeDistance ===\n";
        std::cout << "Side       = " << DRT::DistypeToString( sideType ) << "\n";
        std::cout << "xyze_side_ = " << *xyze_side_;
        std::cout << "px_        = " << *px_;
        std::cout << "\n";
      }
    }

    /// return the absolute value of the distance of the point \c px_ to the given side
    double Distance() const
    {
      switch ( probDim-dimSide )
      {
        case 1:
          return std::abs( distance_[ 0 ] );
        case 2:
          return std::sqrt( distance_[ 0 ] * distance_[ 0 ] + distance_[ 1 ] * distance_[ 1 ] );
        default:
          dserror( "Unsupported probDim and dimSide combination!" );
          exit(EXIT_FAILURE);
      }
      exit(EXIT_FAILURE);
    }

    /// return the signed distance of the point px_ to the given side
    const double * SignedDistance() const
    {
      return distance_;
    }

    /// return the local solution vector (parameter space coordinates + distance)
    const LINALG::Matrix<probDim,1> & LocalCoordinates() const { return xsi_; }

    /// set the initial solution vector to zero
    void SetupSolve()
    {
      xsi_ = 0;

      // evaluate initial rhs value (w/o distance contribution)
      DistanceRHS( *xyze_side_, *px_, b_ );

      tol_ = AdaptiveCombinedNewtonTolerance(*xyze_side_, *px_, b_ );
    }

    /** \brief Setup routine for a new Newton step
     *
     *  Setup the system of equations. */
    void SetupStep( int iter )
    {
      if ( iter > 0 )
      {
        DistanceRHS( *xyze_side_, *px_, b_ );
      }

      // build the linear system of equations
      DistanceSystem( *xyze_side_,*px_, distance_,
          A_, B_, C_, N_, b_ );
    }

    /// compute the rhs value ( w/o distance contributions )
    void DistanceRHS(
        const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
        const LINALG::Matrix<probDim,1> & px,
        LINALG::Matrix<probDim,1> & b )
    {
      const LINALG::Matrix<dimSide,1> xsi_side( xsi_.A(), true );
      DRT::UTILS::shape_function<sideType>( xsi_side, sideFunct_ );
      b = px;
      b.Multiply( -1.0, xyze_side, sideFunct_, 1.0 );
    }

    /// Test the stop / convergence criterion of the Newton scheme
    enum NewtonStatus TestConverged( int iter )
    {
      if ( debug )
      {
        std::cout << "WithinLimits             = "
                  << (WithinLimits<sideType>( xsi_ ) ? "TRUE" : "FALSE") << "\n"
                  << "rsd_.Norm2()             = "
                  << xsi_.Norm2() << "\n"
                  << "dx_.Norm2() / rsd.Norm2()= "
                  << ( iter != 0 ? ( dx_.Norm2() / xsi_.Norm2() ) : -1.0 ) << "\n"
                  << "b_.Norm2() ( RHS )       = "
                  << b_.Norm2() << "\n"
                  << "b_.Norm2() / rsd.Norm2() = "
                  << (iter != 0 ? ( b_.Norm2() / xsi_.Norm2() ) : -1.0 ) << "\n"
                  << "tol_                     = "
                  << tol_ << "\n"
                  << "rsd_ ( aka xsi_ )        = "
                  << xsi_
                  << "b_ (RHS)                 = "
                  << b_ << "\n"
                  << std::flush;
      }

      if ( zeroarea_ )
        return failed;

      double residual = b_.Norm2();
      if ( residual < tol_ )
      {
        // actually I'm not sure why we are doing this here... ( hiermeier 08/16 )
        tol_ *= sqrt( probDim );
        return converged;
      }
      else
        return unconverged;
    }

    /** solve the linear system:
     *  \f[
     *      A \; \Delta \xi = b
     *  \f] */
    bool LinearSolve( int iter )
    {
      if ( zeroarea_ ) return false;
      dx_ = 0;
      //std::cout << "GAUSSELEMIN FOR COMPUTE DISTANCE ... " << std::flush;
      double det = LINALG::gaussElimination<true, probDim>( A_, b_, dx_ );
      if ( debug )
      {
        std::cout << "det                        = " << det << "\n";
        std::cout << "dx_ ( solution increment ) = " << dx_;
        std::cout << "dx_.Norm2()                = " << dx_.Norm2() << "\n";
      }

      return ( det != 0 ); //here det = 0 will just happen if the surface element is distorted!
    }

    /// update the local coordinates solution vector
    bool Update( int iter )
    {
      xsi_ += dx_;
      return true;
    }

    /// fall back routine, if the Newton failed
    bool NewtonFailed()
    {
      tol_ = b_.Norm2();
      #ifdef DEBUG_CUTKERNEL_OUTPUT
       {
       std::stringstream str;
       str << "ComputeDistance: Newton scheme did not converge:\n"
          << std::setprecision( 16 )
          << ( *xyze_side_ )
          << ( *px_ )
          << xsi_
          ;

      std::string filename = OUTPUT::GenerateGmshOutputFilename(".NewtonFailed_distance.pos");
      std::ofstream file(filename.c_str());
      WritetoGmsh(file);
      file.close();
      #endif
      return false;
    }

    /// Gmsh debug output
    void WritetoGmsh(std::ofstream& file)
    {
      file.precision(32); //higher precision!
      char elementType;
      switch (numNodesSide )
      {
      case 3:
        elementType = 'T';
        break;
      case 4:
        elementType = 'Q';
        break;
      default:
        dserror( "unsupported element type in GmshSideDump."
            " Please feel free to extend the functionality if necessary." );
      }

      file << "View \"" << "Side" << "\" {\n";
      {
        file << "S" << elementType
             << "(";
        for ( uint i = 0; i < numNodesSide; ++i )
        {
          if ( i != 0 )
            file << ",";
          file << ( *xyze_side_ )(0,i) << "," << ( *xyze_side_ )(1,i) << "," << ( *xyze_side_ )(2,i);
        }
        file << "){";
        for ( uint i = 0; i < numNodesSide; ++i )
        {
          if ( i != 0 )
            file << ",";
          file << "1";
        }
        file << "};\n";
      }
      file<<"};\n";

      //calculate the points
      LINALG::Matrix<numNodesSide,1> surface;

      DRT::UTILS::shape_function_2D( surface, xsi_( 0 ), xsi_( 1 ), sideType );

      LINALG::Matrix<3,1> x1;;

      x1.Multiply( *xyze_side_, surface );

      file << "View \"" << "Point" << "\" {\n";
      {
        file << "SP (";
        file << ( *px_ )(0,0) << "," << ( *px_ )(1,0)
            << "," << ( *px_ )(2,0);
        file << "){";
        file << 1;
        file << "};\n";

        file << "SP (";
        file << x1(0,0) << "," << x1(1,0)
            << "," << x1(2,0);
        file << "){";
        file << 2;
        file << "};\n";
      }
      file<<"};";
    }

    /// get the adapted Newton tolerance
    double GetTolerance() const
    {
      return tol_;
    }

    /// get actual reached residual
    double GetResidualL2Norm() const
    {
      return b_.Norm2();
    }

    /// return \TRUE if the side metric vectors is zero
    bool ZeroArea() const
    {
      return zeroarea_;
    }

  private:
    enum NormalPlane
    {
      normal_in_xy_plane,
      normal_in_yz_plane,
      normal_plane_undefined
    };

    enum NormalPlane DetectNormalPlane(
        const LINALG::Matrix<probDim,1> & n ) const
    {
      if ( probDim < 3 )
        dserror( "This function makes only sense for the 3-D case!" );

      if ( n( 2 ) == 0.0 )
        return normal_in_xy_plane;
      else if ( n( 0 ) == 0.0 )
        return normal_in_yz_plane;

      dserror( "Couldn't detect a feasible plane for the given normal vector!" );
      exit( EXIT_FAILURE );
    }

    /** \build Build the linear system of equations for the Newton scheme
     *
     *  Right-hand-side:
     *  \f[
     *    b = \hat{x} - x(\xi) - \frac{\tilde{n}}{\| \tilde{n} \|} d,
     *  \f]
     *  where \f$ \hat{x} \in \mathbb{R}^{d}\f$ is the given point,
     *        \f$ \xi \in \mathbb{R}^{d-1}\f$ is the projected in-plane
     *        parameter space coordinate vector,
     *        \f$ \tilde{n} \f$ is the non-unit normal vector on the side
     *        in the point \f$ \xi \f$,
     *        \f$ d \f$ is the distance in normal direction between the point
     *        \f$\hat{x}\f$ and the side.
     *
     *  The matrix \f$ A \f$ holds the consistent linearization.
     *
     *  \author hiermeier
     *  \date 08/16    */
    bool DistanceSystem( const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
        const LINALG::Matrix<probDim,1> & px,
        const double * distance,
        LINALG::Matrix<probDim,probDim> & A,
        LINALG::Matrix<probDim,probDim> & B,
        LINALG::Matrix<probDim,2*dimSide-1> & C,
        LINALG::Matrix<probDim,2> & N,
        LINALG::Matrix<probDim,1> & b )
    {
      LINALG::Matrix<probDim,1> n1( & N( 0, 0 ), true );
      LINALG::Matrix<probDim,1> n2( & N( 0, 1 ), true );

      /* get only the local side coordinates and use them to evaluate
       * the 1-st and 2-nd derivatives */
      const LINALG::Matrix<dimSide,1> xsi_side( xsi_.A(), true );

      // pre-evaluate some important variables
      B = 0;
      double det = EvalDerivsInParameterSpace<probDim,sideType>(xyze_side,
          xsi_side, sideDeriv1_, B, NULL, &n1, &n2, false);
      A.UpdateT( B );

      if ( debug )
      {
        std::cout << "det-metric               = " << det << std::endl;
        if ( det < 1.0e-16)
          std::cout << "!!! Determinant of jacobian is smaller than 1.0e-16 !!!\n";
      }

      if( std::abs( det ) == 0.0 )
      {
        /* then calculation of a normal to a line like this makes no sense at
         * all! */
        zeroarea_ = true;
        return false;
      }

      // inverse of the normal direction norm
      double n1norm_inv = 1.0 / n1.Norm2();

      // --- add the distance part to the right-hand-side
      b.Update(- distance[ 0 ] * n1norm_inv, n1, 1.0);

      // --- evaluate the linearization
      // scale the 1-st normal in A to unit length
      for ( unsigned r = 0; r < probDim; ++r )
        A( r, dimSide ) *= n1norm_inv;

      // --- evaluate 2-nd derivatives at xsi w.r.t. the parameter space coordinates
      C = 0;
      DRT::UTILS::shape_function_deriv2<sideType>( xsi_side, sideDeriv2_ );
      C.MultiplyNT( xyze_side, sideDeriv2_ );

      // reset B-matrix
      B = 0;

      // 1-dimensional side ( a.k.a. edge or line ) embedded in 3-dimensional space
      if ( probDim==3 and dimSide==1 )
      {
        double n2norm_inv = 1.0 / n2.Norm2();

        // scale the 2-nd normal in the last column of A to unit length
        for ( unsigned r = 0; r < probDim; ++r )
          A( r, 2 ) *= n2norm_inv;

        // the normal linearization can be skipped for this case
        if ( Distance() == 0.0 )
          return true;

        // add extra term to the right hand side
        b.Update(- distance[ 1 ] * n2norm_inv, n2, 1.0 );

        // linearization of the 1-st unscaled normal vector
        switch ( DetectNormalPlane( n1 ) )
        {
          case normal_in_xy_plane:
          {
            // d(n1_0) / dr =   d^2(y) / dr^2
            B( 0, 0 ) =   C( 1, 0 );
            // d(n1_1) / dr = - d^2(x) / dr^2
            B( 1, 0 ) = - C( 0, 0 );
            // d(n1_2) / dr =   0.0
            break;
          }
          case normal_in_yz_plane:
          {
            // d(n1_0) / dr =    0.0
            // d(n1_1) / dr =  - d^2(z) / dr^2
            B( 1, 0 ) = - C( 2, 0 );
            // d(n1_2) / dr =    d^2(y) / dr^2
            B( 2, 0 ) =   C( 1, 0 );
            break;
          }
          default:
          {
            dserror( "Shouldn't happen!" );
            exit( EXIT_FAILURE );
          }
        }

        // linearization of the 2-nd unscaled normal vector
        // d(n2_0) / dr
        B( 0, 1 ) = ( C( 1, 0 ) * n1( 2 ) + A( 1, 0 ) * B( 2, 0 )
                  - ( C( 2, 0 ) * n1( 1 ) + A( 2, 0 ) * B( 1, 0 ) ) );
        // d(n2_1) / dr
        B( 1, 1 ) = ( C( 2, 0 ) * n1( 0 ) + A( 2, 0 ) * B( 0, 0 )
                  - ( C( 0, 0 ) * n1( 2 ) + A( 0, 0 ) * B( 2, 0 ) ) );
        // d(n2_2) / dr
        B( 2, 1 ) = ( C( 0, 0 ) * n1( 1 ) + A( 0, 0 ) * B( 1, 0 )
                  - ( C( 1, 0 ) * n1( 0 ) + A( 1, 0 ) * B( 0, 0 ) ) );

        // linearization of the scaling factor ( 1-st normal )
        double lin_nnorm = 0.0;
        lin_nnorm = ( n1( 0 ) * B( 0, 0 ) + n1( 1 ) * B( 1, 0 )
                    + n1( 2 ) * B( 2, 0 ) ) * n1norm_inv * n1norm_inv;
        for ( unsigned r=0; r<probDim; ++r )
          B( r, 0 ) -= lin_nnorm * A( r, 1 );

        // linearization of the scaling factor ( 2-nd normal )
        lin_nnorm = ( n2( 0 ) * B( 0, 1 ) + n2( 1 ) * B( 1, 1 )
                    + n2( 2 ) * B( 2, 1 ) ) * n2norm_inv * n2norm_inv;
        for ( unsigned r=0; r<probDim; ++r )
          B( r, 1 ) -= lin_nnorm * A( r, 2 );

        // final scaling and add
        for ( unsigned r=0; r<probDim; ++r )
        {
          A( r, 0 ) += distance[ 0 ] * B( r, 0 );
          A( r, 0 ) += distance[ 1 ] * B( r, 1 );
        }
      }
      // 2-dimensional side embedded in 3-dimensional space
      else if ( probDim==3 and dimSide==2 )
      {
        // the normal linearization can be skipped for this case
        if ( Distance() == 0.0 )
          return true;

        // linearization of the unscaled normal vector
        // d(n_0) / dr
        B( 0, 0 ) = ( C( 1, 0 ) * A( 2, 1 ) + A( 1, 0 ) * C( 2, 2 )
                  - ( C( 2, 0 ) * A( 1, 1 ) + A( 2, 0 ) * C( 1, 2 ) ) );
        // d(n_0) / ds
        B( 0, 1 ) = ( C( 1, 2 ) * A( 2, 1 ) + A( 1, 0 ) * C( 2, 1 )
                  - ( C( 2, 2 ) * A( 1, 1 ) + A( 2, 0 ) * C( 1, 1 ) ) );

        // d(n_1) / dr
        B( 1, 0 ) = ( C( 2, 0 ) * A( 0, 1 ) + A( 2, 0 ) * C( 0, 2 )
                  - ( C( 0, 0 ) * A( 2, 1 ) + A( 0, 0 ) * C( 2, 2 ) ) );
        // d(n_1) / ds
        B( 1, 1 ) = ( C( 2, 2 ) * A( 0, 1 ) + A( 2, 0 ) * C( 0, 1 )
                  - ( C( 0, 2 ) * A( 2, 1 ) + A( 0, 0 ) * C( 2, 1 ) ) );

        // d(n_2) / dr
        B( 2, 0 ) = ( C( 0, 0 ) * A( 1, 1 ) + A( 0, 0 ) * C( 1, 2 )
                  - ( C( 1, 0 ) * A( 0, 1 ) + A( 1, 0 ) * C( 0, 2 ) ) );
        // d(n_2) / ds
        B( 2, 1 ) = ( C( 0, 2 ) * A( 1, 1 ) + A( 0, 0 ) * C( 1, 1 )
                  - ( C( 1, 2 ) * A( 0, 1 ) + A( 1, 0 ) * C( 0, 1 ) ) );

        // linearization of the scaling factor fact
        double lin_nnorm = 0.0;
        for ( unsigned c=0; c<dimSide; ++c )
        {
          lin_nnorm = ( n1( 0 ) * B( 0, c ) + n1( 1 ) * B( 1, c )
                      + n1( 2 ) * B( 2, c ) ) * n1norm_inv * n1norm_inv;
          for ( unsigned r=0; r<probDim; ++r )
            B( r, c ) -= lin_nnorm * A( r, 2 );
        }

        // final scaling and add
        A.Update( distance[ 0 ] , B , 1.0 );
      }
      // 1-dimensional side/line element embedded in 2-dimensional space
      else if ( probDim==2 and dimSide==1 )
      {
        // the normal linearization can be skipped for this case
        if ( Distance() == 0.0 )
          return true;

        // linearization of the unscaled normal vector
        B(0,0) =  C(1,0);
        B(1,0) = -C(0,0);

        // linearization of the scaling factor fact
        double lin_nnorm = ( C( 1, 0 )+C( 0, 0 ) ) * n1norm_inv * n1norm_inv;
        for ( unsigned r=0; r<probDim; ++r )
          B(r,0) -= lin_nnorm * A( r, 1 );

        // final scaling and add
        A.Update( distance[ 0 ] , B , 1.0 );
      }
      else
        dserror("Unsupported dim <--> probDim relation!");

      return true;
    }

    /// solution vector (parameter space coordinates + distance values)
    LINALG::Matrix<probDim,1> & xsi_;

    /// distance values
    double* distance_;

    /// nodal global positions of the side
    const LINALG::Matrix<probDim,numNodesSide> * xyze_side_;
    /// given global point coordinates
    const LINALG::Matrix<probDim,1> * px_;

    /// nodal shape function values at \c xsi
    static LINALG::Matrix<numNodesSide,1> sideFunct_;

    /** nodal 1-st derivative values at \c xsi
     *
     *  \remark This is an augmented derivatives matrix. The actual
     *  dimension is \f$ dim \times numNodesSide \f$. The remaining
     *  entries are filled with zeros. This is due to the use of an
     *  UTILS function, which expects this input. */
    static LINALG::Matrix<probDim,numNodesSide> sideDeriv1_;

    /** \brief nodal 2-nd derivatives at \c xsi
     *
     *  1-D case: \f$ \x_{,\xi \xi}
     *                \;\in \mathbb{R}^{1 \times N}\f$,
     *
     *  2-D case: \f$ \left( x_{,\xi \xi}, \;x_{,\eta \eta}, \;
     *                       x_{,\xi\eta} \right)
     *                \;\in \mathbb{R}^{3 \times N} \f$ */
    static LINALG::Matrix<2*dimSide-1,numNodesSide> sideDeriv2_;

    /// complete linearization matrix
    static LINALG::Matrix<probDim,probDim> A_;
    /// auxiliary matrix
    static LINALG::Matrix<probDim,probDim> B_;
    /** auxiliary matrix holding 2-nd derivatives at \c xsi_
     *  w.r.t. the parameter space coordinates */
    static LINALG::Matrix<probDim,2*dimSide-1> C_;
    /// right hand side vector
    static LINALG::Matrix<probDim,1> b_;
    /** \brief solution increment
     *
     * (parameter space coordinates of the side + distance increment) */
    static LINALG::Matrix<probDim,1> dx_;
    /// (unscaled) normal vectors
    static LINALG::Matrix<probDim,2> N_;

    /// adapted Newton tolerance
    double tol_;

    /** \brief  indicates that the surface element has a area value of zero
     *
     *  --> cross product for normal is zero!!! */
    bool zeroarea_;
  }; // class ComputeDistanceStrategy


  /*--------------------------------------------------------------------------*/
  /** \brief generic distance between side and point
   *
   *  inheritance diagram:
   *
   *  ComputeDistance --> GenericComputeDistance --> NewtonSolve
   *                      ^^^^^^^^^^^^^^^^^^^^^^
   *  --> ComputeDistanceStrategy --> EmptyNewtonStrategy  */
  template <class Strategy,
            unsigned probDim,
            DRT::Element::DiscretizationType sideType,
            unsigned dimSide      = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class GenericComputeDistance : Strategy
  {
  public:
    /// constructor
    GenericComputeDistance( LINALG::Matrix<probDim,1> & xsi, bool checklimits )
      : Strategy( xsi, checklimits )
    {
    }

    /// start and solve the distance calculation
    bool operator()( const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
                     const LINALG::Matrix<probDim,1> & px,
                     double & distance,
                     bool signeddistance = false)
    {
      this->Setup( xyze_side, px );

      bool conv = this->Solve();

      if ( not signeddistance )
        distance = this->Distance();
      else
      {
        switch ( probDim-dimSide )
        {
          case 1:
            distance = this->SignedDistance()[0];
            break;
          default:
            dserror( "A scalar signed distance value is not available!" );
            exit( EXIT_FAILURE );
        }
      }

      return conv;
    }

    // get the local coordinates
    const LINALG::Matrix<probDim,1> & LocalCoordinates() const
    {
      return Strategy::LocalCoordinates();
    }

    bool SurfaceWithinLimits() const
    {
      return GEO::CUT::KERNEL::WithinLimits<sideType>( LocalCoordinates() );
    }

    /// access the Newton tolerance
    double GetTolerance() const
    {
      return Strategy::GetTolerance();
    }

    /// access the l2 norm of the reached residual
    double GetResidualL2Norm() const
    {
      return Strategy::GetResidualL2Norm();
    }

    bool ZeroArea() const
    {
      return Strategy::ZeroArea();
    }
  };  // class GenericComputeDistance

  /*--------------------------------------------------------------------------*/
  /** \brief most derived strategy class for distance between side and point
   *
   *  inheritance diagram:
   *
   *  ComputeDistance --> GenericComputeDistance --> NewtonSolve
   *  ^^^^^^^^^^^^^^^
   *  --> ComputeDistanceStrategy --> EmptyNewtonStrategy  */
  template <unsigned probDim,
            DRT::Element::DiscretizationType sideType,
            unsigned dimSide       = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesSide  = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class ComputeDistance
      : public GenericComputeDistance<NewtonSolve<ComputeDistanceStrategy<false,
                                                                          probDim,
                                                                          sideType>,
                                                  probDim>,
                                      probDim,sideType>
  {
  public:
    /// constructor
    ComputeDistance( LINALG::Matrix<probDim,1> & xsi, bool checklimits =true )
        : GenericComputeDistance<NewtonSolve<ComputeDistanceStrategy<false,
                                                                   probDim,
                                                                   sideType>,
                                           probDim>,
                               probDim,sideType>( xsi, checklimits )
    {
    }
  };  //class ComputeDistance

  /*--------------------------------------------------------------------------*/
  template <unsigned probDim,
            DRT::Element::DiscretizationType sideType,
            unsigned dimSide      = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class DebugComputeDistance
      : public GenericComputeDistance<NewtonSolve<DebugNewtonStrategy<ComputeDistanceStrategy<true,
                                                                                              probDim,
                                                                                              sideType>,
                                                                      probDim>,
                                                  probDim>,
                                      probDim,
                                      sideType>
  {
  public:
    /// constructor
    DebugComputeDistance( LINALG::Matrix<probDim,1> & xsi, bool checklimits = true )
      : GenericComputeDistance<NewtonSolve<DebugNewtonStrategy<ComputeDistanceStrategy<true,
                                                                                       probDim,
                                                                                       sideType>,
                                                               probDim>,
                                           probDim>,
                               probDim,
                               sideType>( xsi, checklimits )
    {
    }
  };  // class DebugComputeDistance


  /*--------------------------------------------------------------------------*/
  /** \brief strategy class for intersection of side and edge/line
   *
   *  inheritance diagram:
   *
   *  ComputeIntersection --> GenericComputeIntersection --> NewtonSolve
   *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy
   *      ^^^^^^^^^^^^^^^^^^^^^^^^^^^  */
  template <bool debug,
            unsigned probDim,
            DRT::Element::DiscretizationType edgeType,
            DRT::Element::DiscretizationType sideType,
            unsigned dimEdge       = DRT::UTILS::DisTypeToDim<edgeType>::dim,
            unsigned dimSide       = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesEdge  = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
            unsigned numNodesSide  = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class ComputeIntersectionStrategy : public EmptyNewtonStrategy
  {
  public:
    /// constructor
    ComputeIntersectionStrategy( LINALG::Matrix<dimEdge+dimSide,1> & xsi, bool checklimits )
        : xsi_( xsi ),
          xyze_side_( NULL ),
          xyze_edge_( NULL ),
          residual_(0.0),
          off_count_(0),
          check_limits_( checklimits ),
          tol_(0.0)
    {
    }

    void Setup( const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
                const LINALG::Matrix<probDim,numNodesEdge> & xyze_edge )
    {
      xyze_side_ = &xyze_side;
      xyze_edge_ = &xyze_edge;

      if ( debug )
      {
        std::cout << "\n\n === ComputeIntersection ===\n";
        std::cout << "--- Setup()\n";
        std::cout << "  Edge = " << DRT::DistypeToString( edgeType ) << "\n";
        std::cout << "  Side = " << DRT::DistypeToString( sideType ) << "\n";
        std::cout << "  xyze_side = " << std::setprecision( 15 ) << xyze_side;
        std::cout << "  xyze_edge = " << std::setprecision( 15 ) << xyze_edge;
      }
    }

    const LINALG::Matrix<dimEdge+dimSide,1> & LocalCoordinates() const
    {
      return xsi_;
    }

    void SetupSolve()
    {
      xsi_ = 0;
      dx_ = 0;
      off_count_ = 0;
      residual_ = 0;

      // evaluate initial rhs value
      const LINALG::Matrix<dimSide,1> xsi_side(xsi_.A(),true);
      const LINALG::Matrix<dimEdge,1> xsi_edge(xsi_.A()+dimSide,true);
      IntersectionRHS( xsi_edge, xsi_side, *xyze_edge_, *xyze_side_, c_ );

      tol_ = AdaptiveCombinedNewtonTolerance(*xyze_side_, *xyze_edge_, c_ );
    }

    void SetupStep( int iter )
    {
      // build linear system of equations
      const LINALG::Matrix<dimSide,1> xsi_side(xsi_.A(),true);
      const LINALG::Matrix<dimEdge,1> xsi_edge(xsi_.A()+dimSide,true);

      if ( iter > 0 )
      {
        IntersectionRHS( xsi_edge, xsi_side, *xyze_edge_, *xyze_side_, c_ );
      }

      IntersectionSystem(xsi_edge,xsi_side,*xyze_edge_,*xyze_side_,c_,A_,B_,b_);
    }

    /// compute the right-hand-side
    void IntersectionRHS(
        const LINALG::Matrix< dimEdge,1 >& xsi_edge,
        const LINALG::Matrix< dimSide,1 >& xsi_side,
        const LINALG::Matrix<probDim,numNodesEdge> & xyze_edge,
        const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
        LINALG::Matrix< probDim,1 >& c )
    {
      DRT::UTILS::shape_function<sideType>( xsi_side, sideFunct_ );
      DRT::UTILS::shape_function<edgeType>( xsi_edge, edgeFunct_ );

      c.Multiply( xyze_edge, edgeFunct_ );
      c.Multiply( -1, xyze_side, sideFunct_, 1 );
    }


    enum NewtonStatus TestConverged( int iter )
    {
      residual_ = b_.Norm2();

      if ( debug )
      {
        std::cout << "--- TestConverged\n";
        double incr = dx_.Norm2();
        std::cout << "  residual = " << std::setprecision( 15 ) << residual_
                  << "  incr = " << std::setprecision( 15 ) << incr
                  << "\n";
      }

      if ( c_.Norm2() < tol_ && iter > 0)
      {
        tol_ *= sqrt( probDim );
        return converged;
      }

      /* If the norm of c_ does not but the norm of b_ does fulfill the criterion,
       * the least square scheme failed! */
      if ( residual_ < tol_ && iter > 0)
      {
        return failed;
      }

      return unconverged;
    }

    bool LinearSolve( int iter )
    {
      dx_ = 0;

      double det = LINALG::gaussElimination<true, dimSide+dimEdge>( A_, b_, dx_ );

      if ( debug )
      {
        std::cout << "--- LinearSolve\n";
        std::cout << "  A_ = " << A_;
        std::cout << "  b_ = " << b_;
        std::cout << "  det = " << std::setprecision( 15 ) << det << "\n";
        std::cout << "  dxi_ = " << dx_;
      }

      return ( det != 0.0 );
    }

    bool Update( int iter )
    {
      xsi_.Update( 1.0, dx_, 1.0 );

//      std::cout << "xsi: " << xsi_ << "   dx: " << dx_ << std::endl;

      if ( debug )
      {
        std::cout << "--- Update\n";
        const LINALG::Matrix<dimSide,1> xsi_side( xsi_.A(), true );
        const LINALG::Matrix<dimEdge,1> xsi_edge( xsi_.A()+dimSide, true );
        std::cout << "  off_count=" << off_count_
                  << "\n  Side within limits = " << ( WithinLimits<sideType>( xsi_side ) ? "TRUE" : "FALSE" )
                  << "\n  Edge within limits = " << ( WithinLimits<edgeType>( xsi_edge ) ? "TRUE" : "FALSE" )
                  << "\n";
        std::cout << "  xsi = " << xsi_;
        std::cout << "  xsi_side = " << xsi_side;
        std::cout << "  xsi_edge = " << xsi_edge;

        LINALG::Matrix<numNodesSide,1> sideFunct;
        LINALG::Matrix<numNodesEdge,1> edgeFunct;

        DRT::UTILS::shape_function<sideType>( xsi_side, sideFunct );
        DRT::UTILS::shape_function<edgeType>( xsi_edge, edgeFunct );

        LINALG::Matrix<probDim,1> x1;
        LINALG::Matrix<probDim,1> x2;

        x1.Multiply( *xyze_side_, sideFunct );
        x2.Multiply( *xyze_edge_, edgeFunct );

        std::cout << "  x_side =" << x1;
        std::cout << "  x_edge =" << x2;

        x2.Update( 1, x1, -1 );
        std::cout << "  (x_side - x_edge) =" << x2;
      }

      return true;
    }

    bool NewtonFailed()
    {
      tol_ = c_.Norm2()*sqrt(3);
      //if (debug)
      #ifdef DEBUG_CUTKERNEL_OUTPUT
        {
        std::stringstream str;
        str << "Newton scheme did not converge:\n  "
            << ( *xyze_side_ )
            << ( *xyze_edge_ )
          ;
        DumpDoubles( str, xyze_side_->A(), xyze_side_->M()*xyze_side_->N() );
        str << "\n";
        DumpDoubles( str, xyze_edge_   ->A(), xyze_edge_   ->M()*xyze_edge_   ->N() );


        std::stringstream str;
        str  <<".NewtonFailed_intersection.pos";
        std::string filename(GEO::CUT::OUTPUT::GenerateGmshOutputFilename(str.str()));
        std::ofstream file(filename.c_str());
        WritetoGmsh(file);
        }
      #endif
      return false;
    }

    double GetTolerance() const
    {
      return tol_;
    }

    void WritetoGmsh(std::ofstream& file)
    {
      file.precision(32); //higher precision!
      char elementType = '\0';
      switch (numNodesSide )
      {
      case 3:
        elementType = 'T';
        break;
      case 4:
        elementType = 'Q';
        break;
      default:
        dserror( "unsupported element type in WritetoGmsh."
            " Please feel free to extend the functionality if necessary." );
      }

      file << "View \"" << "Side" << "\" {\n";
      {
        file << "S" << elementType
             << "(";
        for ( uint i = 0; i < numNodesSide; ++i )
        {
          if ( i != 0 )
            file << ",";
          file << ( *xyze_side_ )(0,i) << "," << ( *xyze_side_ )(1,i) << "," << ( *xyze_side_ )(2,i);
        }
        file << "){";
        for ( uint i = 0; i < numNodesSide; ++i )
        {
          if ( i != 0 )
            file << ",";
          file << "1";
        }
        file << "};\n";
      }
      file<<"};\n";

      file << "View \"" << "Line" << "\" {\n";
      for (unsigned iline = 0; iline < numNodesEdge -1; ++iline)
      {
        file << "SL (";
        file << ( *xyze_edge_ )(0,iline) << "," << ( *xyze_edge_ )(1,iline)
            << "," << ( *xyze_edge_ )(2,iline) << ",";
        file << ( *xyze_edge_ )(0,iline+1) << "," << ( *xyze_edge_ )(1,iline+1)
            << "," << ( *xyze_edge_ )(2,iline+1);
        file << "){";
        file << 1<< ",";
        file << 1;
        file << "};\n";
      }
      file<<"};\n";


      //calculate the points
      LINALG::Matrix<numNodesSide,1> surface;
      LINALG::Matrix<numNodesEdge,1> line;

      DRT::UTILS::shape_function_2D( surface, xsi_( 0 ), xsi_( 1 ), sideType );
      DRT::UTILS::shape_function_1D( line, xsi_( 2 ), edgeType );

      LINALG::Matrix<probDim,1> x1;
      LINALG::Matrix<probDim,1> x2;

      x1.Multiply( *xyze_side_, surface );
      x2.Multiply( *xyze_edge_, line );
      //endof


      file << "View \"" << "Point" << "\" {\n";
      {
        file << "SP (";
        file << x1(0,0) << "," << x1(1,0)
            << "," << x1(2,0);
        file << "){";
        file << 1;
        file << "};\n";

        file << "SP (";
        file << x2(0,0) << "," << x2(1,0)
            << "," << x2(2,0);
        file << "){";
        file << 2;
        file << "};\n";
      }
      file<<"};";
    }

  private:
    /** \brief Build the linear system of equations
     *
     *  Right-hand-side:
     *  \f[
     *    b = x(\xi) - x(\tilde{\xi}),
     *  \f]
     *  where \f$\xi \in \mathbb{R}^{dimEgde}\f$ is the parameter space
     *  coordinate of the edge-element and \f$ \tilde{\xi} \in \mathbb{R}^{dimSide} \f$
     *  are the parameter space coordinates of the side element.
     *
     *  The matrix \f$ A \f$ holds the consistent linearization. If the problem
     *  dimension is larger than the accumulated dimensions of the side and edge
     *  elements, a least square approach will be used.
     *
     *  \author hiermeier
     *  \date 08/16 */
    void IntersectionSystem(
        const LINALG::Matrix< dimEdge,1 >& xsi_edge,
        const LINALG::Matrix< dimSide,1 >& xsi_side,
        const LINALG::Matrix< probDim,numNodesEdge >& xyze_edge,
        const LINALG::Matrix< probDim,numNodesSide >& xyze_side,
        LINALG::Matrix< probDim,1 >& c,
        LINALG::Matrix< dimEdge+dimSide,dimEdge+dimSide >& A,
        LINALG::Matrix< probDim,dimEdge+dimSide >& B,
        LINALG::Matrix< dimEdge+dimSide,1 >& b
        )
    {
      // --- compute the linearization
      DRT::UTILS::shape_function_deriv1<sideType>( xsi_side, sideDeriv1_ );
      DRT::UTILS::shape_function_deriv1<edgeType>( xsi_edge, edgeDeriv1_ );

      A = 0.0;
      B = 0.0;
      // linearization w.r.t. the parameter space coordinate(s) of the side
      for ( unsigned inode=0; inode<numNodesSide; ++inode )
        for( unsigned row=0; row<probDim; ++row )
          for ( unsigned col=0; col<dimSide; ++col )
            B(row,col) += xyze_side(row,inode) * sideDeriv1_(col,inode);

      // linearization w.r.t. the parameter space coordinate(s) of the edge
      for ( unsigned inode=0; inode<numNodesEdge; ++inode )
        for ( unsigned row=0; row<probDim; ++row )
          for ( unsigned col=0; col<dimEdge; ++col)
            B(row,dimSide+col) -= xyze_edge(row,inode) * edgeDeriv1_(col,inode);

      // default case
      if (probDim == dimSide+dimEdge)
      {
        A.SetView(B.A());
        b.SetView(c.A());
      }
      /* Actually there is only one relevant case, which comes into my mind
       * and that is the intersection of two lines in 3 dimensions. Anyway
       * we use a more general if-statement here.
       *
       * If we have to handle such a configuration we use a Least-Squares
       * approach. */
      else if (probDim > dimSide+dimEdge)
      {
        A.MultiplyTN( B, B );
        b.MultiplyTN( B, c );
      }
      else
        dserror("The problem dimension is smaller than the combination of"
            "the side and edge element dimensions. This case is currently "
            "unsupported!");
    }

  private:
    /// local coordinates (xi1_side, xi2_side, xi_line)
    LINALG::Matrix<dimEdge+dimSide,1> & xsi_;

    const LINALG::Matrix<probDim,numNodesSide> * xyze_side_;
    const LINALG::Matrix<probDim,numNodesEdge> * xyze_edge_;

    static LINALG::Matrix<numNodesSide,1> sideFunct_;
    static LINALG::Matrix<numNodesEdge,1> edgeFunct_;
    static LINALG::Matrix<dimSide,numNodesSide> sideDeriv1_;
    static LINALG::Matrix<dimEdge,numNodesEdge> edgeDeriv1_;

    static LINALG::Matrix<dimEdge+dimSide,dimEdge+dimSide> A_;
    static LINALG::Matrix<probDim,dimEdge+dimSide> B_;
    static LINALG::Matrix<dimEdge+dimSide,1> b_;
    static LINALG::Matrix<probDim,1> c_;

    /// increment in local coordinates d(xi1_side, xi2_side, xi_line)
    static LINALG::Matrix<dimEdge+dimSide,1> dx_;

    /// residual norm
    double residual_;

    /** count how many iterations during the Newton the local coordinates
     *  are not within limits! Stop the newton scheme if the iteration is
     *  too often outside the limits */
    int off_count_;

    /** shall we check the limits? ( seems out-dated, since we donot check
     *  limits by default ) */
    bool check_limits_;

    //basic tolerance!!!
    double tol_;
  }; // class ComputeIntersectionStrategy

  /*--------------------------------------------------------------------------*/
  /** \brief generic strategy class for intersection of side and line/edge
   *
   *  inheritance diagram:
   *
   *  ComputeIntersection --> GenericComputeIntersection --> NewtonSolve
   *                          ^^^^^^^^^^^^^^^^^^^^^^^^^^
   *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy  */
  template <class Strategy,
            unsigned probDim,
            DRT::Element::DiscretizationType edgeType,
            DRT::Element::DiscretizationType sideType,
            unsigned dimEdge      = DRT::UTILS::DisTypeToDim<edgeType>::dim,
            unsigned dimSide      = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
            unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class GenericComputeIntersection : Strategy
  {
  public:
    //! constructor
    GenericComputeIntersection( LINALG::Matrix<dimEdge+dimSide,1> & xsi, bool checklimits)
      : Strategy( xsi, checklimits)
    {
    }

    bool operator()( const LINALG::Matrix<probDim,numNodesSide> & xyze_side,
                     const LINALG::Matrix<probDim,numNodesEdge> & xyze_edge )
    {
      this->Setup( xyze_side, xyze_edge );
      if ( this->Solve() )
      {
        return true;
      }
      return false;
    }

    const LINALG::Matrix<dimEdge+dimSide,1> & LocalCoordinates() const
    {
      return Strategy::LocalCoordinates();
    }

    bool SurfaceWithinLimits() const
    {
      return WithinLimits<sideType>( LocalCoordinates() );
    }

    bool LineWithinLimits() const
    {
      const LINALG::Matrix<dimEdge,1> xsi_line( LocalCoordinates().A()+dimSide, true );
      return WithinLimits<edgeType>( xsi_line );
    }

    double GetTolerance() const
    {
      return Strategy::GetTolerance();
    }

    void WritetoGmsh(std::ofstream& file)
    {
      Strategy::WritetoGmsh(file);
    }
  }; // class GenericComputeIntersection


  /*--------------------------------------------------------------------------*/
  /** \brief most derived strategy class for intersection of side and line
   *
   *  inheritance diagram:
   *
   *  ComputeIntersection --> GenericComputeIntersection --> NewtonSolve
   *  ^^^^^^^^^^^^^^^^^^^
   *  --> ComputeIntersectionStrategy --> EmptyNewtonStrategy  */
  template <unsigned probDim,
            DRT::Element::DiscretizationType edgeType,
            DRT::Element::DiscretizationType sideType,
            unsigned dimEdge      = DRT::UTILS::DisTypeToDim<edgeType>::dim,
            unsigned dimSide      = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
            unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class ComputeIntersection :
      public GenericComputeIntersection<NewtonSolve<ComputeIntersectionStrategy<false,
                                                                                probDim,
                                                                                edgeType,
                                                                                sideType>,
                                                    dimEdge+dimSide>,
                                        probDim,edgeType,sideType>
  {
  public:
    //! constructor
    ComputeIntersection( LINALG::Matrix<dimEdge+dimSide,1> & xsi, bool checklimits = true )
        : GenericComputeIntersection<NewtonSolve<ComputeIntersectionStrategy<false,
                                                                             probDim,
                                                                             edgeType,
                                                                             sideType>,
                                                 dimEdge+dimSide>,
                                     probDim,edgeType,sideType>( xsi, checklimits )
    {
    }
  };  // class ComputeIntersection

  /*--------------------------------------------------------------------------*/
  /** \brief most derived strategy class for debugging intersection of side
   *  and line */
  template <unsigned probDim,
            DRT::Element::DiscretizationType edgeType,
            DRT::Element::DiscretizationType sideType,
            unsigned dimEdge      = DRT::UTILS::DisTypeToDim<edgeType>::dim,
            unsigned dimSide      = DRT::UTILS::DisTypeToDim<sideType>::dim,
            unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgeType>::numNodePerElement,
            unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sideType>::numNodePerElement>
  class DebugComputeIntersection
      : public GenericComputeIntersection<NewtonSolve<DebugNewtonStrategy<ComputeIntersectionStrategy<true,
                                                                                                      probDim,
                                                                                                      edgeType,
                                                                                                      sideType>,
                                                                          dimEdge+dimSide>,
                                                      dimEdge+dimSide >,
                                          probDim,edgeType,sideType>
  {
  public:
    //! constructor
    DebugComputeIntersection( LINALG::Matrix<dimEdge+dimSide,1> & xsi, bool checklimits = true )
      : GenericComputeIntersection<NewtonSolve<DebugNewtonStrategy<ComputeIntersectionStrategy<true,
                                                                                               probDim,
                                                                                               edgeType,
                                                                                               sideType>,
                                                                   dimEdge+dimSide>,
                                               dimEdge+dimSide >,
                                   probDim,edgeType,sideType>( xsi, checklimits )
    {
    }

  }; // class DebugComputeIntersection

} // namespace KERNEL
} // namespace CUT
} // namespace GEO


//static members in KERNEL
//in compute position strategy
template <bool debug,unsigned probDim,DRT::Element::DiscretizationType elementType, unsigned numNodesElement, unsigned dim > LINALG::Matrix<numNodesElement, 1> GEO::CUT::KERNEL::ComputePositionStrategy<debug, probDim, elementType, numNodesElement, dim >::funct_;
template <bool debug,unsigned probDim,DRT::Element::DiscretizationType elementType, unsigned numNodesElement, unsigned dim > LINALG::Matrix<probDim, numNodesElement> GEO::CUT::KERNEL::ComputePositionStrategy<debug, probDim, elementType, numNodesElement, dim >::deriv1_;
template <bool debug,unsigned probDim,DRT::Element::DiscretizationType elementType, unsigned numNodesElement, unsigned dim > LINALG::Matrix<probDim, probDim> GEO::CUT::KERNEL::ComputePositionStrategy<debug, probDim, elementType, numNodesElement, dim >::A_;
template <bool debug,unsigned probDim,DRT::Element::DiscretizationType elementType, unsigned numNodesElement, unsigned dim > LINALG::Matrix<probDim, 1> GEO::CUT::KERNEL::ComputePositionStrategy<debug, probDim, elementType, numNodesElement, dim>::b_;
template <bool debug,unsigned probDim,DRT::Element::DiscretizationType elementType, unsigned numNodesElement, unsigned dim > LINALG::Matrix<probDim, 1> GEO::CUT::KERNEL::ComputePositionStrategy<debug, probDim, elementType, numNodesElement, dim >::dx_;

//in compute distance strategy
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<numNodesSide, 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::sideFunct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, numNodesSide> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::sideDeriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<2*dimSide-1, numNodesSide> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::sideDeriv2_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, probDim> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::A_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, probDim> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::B_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, 2*dimSide-1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::C_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::b_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::dx_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType sideType, unsigned dimSide, unsigned numNodesSide > LINALG::Matrix<probDim, 2> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, probDim, sideType, dimSide, numNodesSide >::N_;

//in compute unsignedersection strategy
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<numNodesSide, 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::sideFunct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<numNodesEdge, 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::edgeFunct_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<dimSide, numNodesSide> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::sideDeriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<dimEdge, numNodesEdge> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::edgeDeriv1_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<dimEdge+dimSide, dimEdge+dimSide> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::A_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<probDim, dimEdge+dimSide> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::B_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<dimEdge+dimSide, 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::b_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<probDim, 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::c_;
template <bool debug, unsigned probDim, DRT::Element::DiscretizationType edgeType, DRT::Element::DiscretizationType sideType, unsigned dimEdge, unsigned dimSide, unsigned numNodesEdge , unsigned numNodesSide > LINALG::Matrix<dimEdge+dimSide, 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, probDim, edgeType, sideType, dimEdge, dimSide, numNodesEdge, numNodesSide >::dx_;

#endif
