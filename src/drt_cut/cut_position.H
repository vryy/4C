/*----------------------------------------------------------------------------*/
/*!
\file cut_position.H

\brief Class to check whether a point lies inside an element, and also to
       compute local coordinates of a point with respect to the element.
       The standard as well as the embedded cases are treated here. Feel
       free to extend the content for your needs.

\level 2

\maintainer Michael Hiermeier
*/
/*----------------------------------------------------------------------------*/

#ifndef CUT_POSITION_H
#define CUT_POSITION_H

#include "cut_element.H"
#include "cut_kernel.H"

namespace GEO {
namespace CUT {
  class BoundingBox;
  class PositionFactory;

  /*----------------------------------------------------------------------------*/
  /** \brief Base class of all \e Position objects
   *
   *  \author hiermeier \date 08/16 */
  class Position
  {
  public:
    /** \brief Status of the position computation.
     *
     *  Please not that the ordering is important, since the assigned integers are
     *  used for comparison operations. As soon as all position cases are working
     *  flawlessly this enum list can and should be adapted.
     *
     *  \author hiermeier \date 01/17 */
    enum Status
    {
      position_invalid          = -2,   ///< The calculated position is invalid.
      position_zero_area        = -1,   ///< A zero element was detected during the position calculation.
      position_unevaluated      =  0,   ///< The Compute() routine has not been called, yet.
      position_distance_valid   =  1,   ///< The distance but not the position is valid.
      position_valid            =  2    ///< The distance and the position are valid.
    };

    //! Map Status enum to std::string
    static inline std::string Status2String( enum Status pstatus )
    {
      switch ( pstatus )
      {
        case position_unevaluated :
          return "position_unevaluated";
        case position_zero_area :
          return "position_zero_area";
        case position_invalid :
          return "position_invalid";
        case position_valid :
          return "position_valid";
        case position_distance_valid:
          return "position_distance_valid";
        default :
          return "Unknown PositionStatus";
      }
      exit(EXIT_FAILURE);
    };

  public:
    /// @ name create methods
    /// @{

    /** \brief build variant #1
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param point   (in) : Given global point object. */
    static Teuchos::RCP<Position> Create(
        const Element & element,
        const Point & point );

    /** \brief build variant #2
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param xyz     (in) : Global coordinates of the given point.
     *
     *  \author hiermeier \date 08/16 */
    template < unsigned rdim >
    static Teuchos::RCP<Position> Create(
        const Element & element,
        const LINALG::Matrix<rdim, 1> & xyz );

    /** \brief build variant #3
     *
     *  \param xyze    (in) : Global nodal positions of the element.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type.
     *
     *  \author hiermeier \date 08/16 */
    template < unsigned rdim, unsigned cdim, unsigned rdim_2 >
    static Teuchos::RCP<Position> Create(
        const LINALG::Matrix< rdim, cdim > & xyze,
        const LINALG::Matrix< rdim_2, 1 > & xyz,
        const DRT::Element::DiscretizationType & distype );

    /// \brief build variant #3-1 (Epetra_SerialDenseMatrix)
    template < unsigned rdim >
    static Teuchos::RCP<Position> Create(
        const Epetra_SerialDenseMatrix & xyze,
        const LINALG::Matrix< rdim, 1 > & xyz,
        const DRT::Element::DiscretizationType & distype );

    /** \brief build variant #4
     *
     *  \param nodes   (in) : Nodes of the element we want to check.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type (optional).
     *
     *  \author hiermeier \date 08/16 */
    template < unsigned rdim >
    static Teuchos::RCP<Position> Create(
        const std::vector<Node*> nodes,
        const LINALG::Matrix<rdim, 1> & xyz,
        DRT::Element::DiscretizationType distype = DRT::Element::dis_none);

    /// @}
  public:
    /// constructor
    Position() { /* intentionally left blank */ };

    /// destructor
    virtual ~Position() {};

    virtual unsigned Dim() const = 0;

    virtual unsigned ProbDim() const = 0;

    virtual enum Status Status() const = 0;

    bool Compute()
    { return Compute(POSITIONTOL); }
    virtual bool Compute( const double& tol ) = 0;
    bool Compute(const bool& allow_dist)
    { return Compute(POSITIONTOL,allow_dist); }
    virtual bool Compute( const double& tol,
        const bool& allow_dist)
    {
      dserror("Unsupported for the standard case!");
      exit(EXIT_FAILURE);
    }

    virtual bool IsGivenPointWithinElement() = 0;

    virtual const double & NewtonTolerance() const = 0;

    /*! \brief Return the local coordinates of the given point \c px_ */
    template < class T >
    void LocalCoordinates(T & rst)
    {
      if (rst.M() < Dim())
        dserror("rst has the wrong row number! ( DIM = %d ( rst.M() = %d < DIM ) )",Dim(),
            rst.M());
      LocalCoordinates(rst.A());

      std::fill( rst.A() + Dim(), rst.A() + rst.M(), 0.0 );
    }

    /*! \brief Return the scalar signed perpendicular distance between
     *         given point and embedded element */
    virtual double Distance() const
    {
      dserror("Unsupported for the standard case!");
      exit(EXIT_FAILURE);
    }

    /*! \brief Return the signed perpendicular distance vector between
     *         given point and embedded element.
     *
     *  If you are interested in the distance between a line and a point
     *  in 3-D, you will get 2 distance values by calling this function.
     *  The scalar variant won't work, because it is designed to give you
     *  a signed scalar distance value, what is not possible in this
     *  special case.
     *
     *  \param distance (in) : access calculated distance value(s)
     *
     *  \author hiermeier \date 01/17 */
    template < class T >
    void Distance( T & distance ) const
    {
      if (distance.M() < ( ProbDim() - Dim() ) )
        dserror("The distance vector has the wrong row number! "
            "( DIM = %d ( rst.M() = %d < DIM ) )",Dim(), distance.M());

      Distance( distance.A() );
    }

    bool WithinLimits() const
    { return WithinLimitsTol(POSITIONTOL); }
    virtual bool WithinLimitsTol(const double& Tol) const
    { return WithinLimitsTol(Tol,false); }
    bool WithinLimits(const bool& allow_dist) const
    { return WithinLimitsTol(POSITIONTOL, allow_dist); }
    virtual bool WithinLimitsTol(const double& Tol, const bool& allow_dist) const
    {
      dserror("Unsupported for the standard case!");
      exit(EXIT_FAILURE);
    }

  protected:
    /*! \brief Return the local coordinates of the given point \c px_ */
    virtual void LocalCoordinates(double * rst) = 0;

    virtual void Distance( double * distance ) const = 0;

    /*! \brief Scaling and shifting of the input nodal positions and the location
     *         of the point \c px_
     *
     *  The corner points of the side and the location of the point
     *  are converted into the corresponding range of the natural coordinate
     *  system of the element. This also gives the proper initial choice for
     *  the Newton scheme to find the local coordinate system. */
    virtual void ScaleAndShift() = 0;
  };  // class Position

  /*----------------------------------------------------------------------------*/
  /*! \brief Class to check whether a point lies inside an element, and to
   *  compute local coordinates of a point with respect to the element
   *
   *  This class is a template and capable of different kinds of embedded
   *  elements in a higher dimensional space (manifolds). Considered are 1-D elements
   *  embedded in a 2-D space, as well 2-D elements embedded in a 3-D space. Not
   *  considered is the case of 1-D elements embedded in a 3-D space. Nevertheless,
   *  the extension is possible.
   *
   *  The class is a template on the problem dimension \c probdim,
   *                             the element/side type \c eletype,
   *                             the number of nodes per element/side \c numNodesElement,
   *                             the element dimension \c dim.
   *
   *
   *  \author hiermeier
   *  \date 08/16 */
  template <unsigned probdim,
      DRT::Element::DiscretizationType eletype,
      unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
  class PositionGeneric : public Position
  {
  public:
    /** \brief constructor
     *
     *  \param xyze (in) : Global nodal positions of the element.
     *  \param xyz  (in) : Global coordinates of the given point. */
    PositionGeneric( const LINALG::Matrix<probdim,numNodesElement> & xyze,
              const LINALG::Matrix<probdim,1> & xyz )
      : px_( xyz ),
        xyze_( xyze ),
        scale_( 1.0 ),
        pos_status_( position_unevaluated ),
        compute_tolerance_( -1.0 )
    {
      ScaleAndShift();
    }

    /// destructor
    virtual ~PositionGeneric() {};

    /* Don't call any of these methods directly! Use the base class public
     * methods, instead. */
  protected:
    /// return the element dimension
    virtual unsigned Dim() const
    {
      return dim;
    }

    /// return the problem dimension
    virtual unsigned ProbDim() const
    {
      return probdim;
    }

    virtual enum Status Status() const
    {
      return pos_status_;
    }

    /*! \brief Return the local coordinates of the given point \c px_
     *
     *  The xsi_is_valid flag indicates, whether the Compute() function has
     *  been called successfully or not. In the embedded case, it is possible, that
     *  the fall-back method of the Compute function succeeded. Anyway, the
     *  local coordinates are not updated consistently in this case. If you
     *  run into this, you will have to take a closer look.
     *
     *  \param rst: local coordinates.
     *
     *  \author hiermeier */
    virtual void LocalCoordinates( double * rst )
    {
      if ( pos_status_ != position_valid )
      {
        std::ostringstream msg;
        msg << "The local coordinates are not valid. "
            "( Position::Status = " << Status2String( pos_status_ ) << " )";
        run_time_error( msg.str() );
      }

     std::copy( xsi_.A(), xsi_.A()+xsi_.M(), rst );
    }

    virtual void Distance ( double * distance ) const = 0;

    virtual bool Compute( const double& tol ) = 0;

    virtual bool IsGivenPointWithinElement() = 0;

    const double & NewtonTolerance() const
    {
      if ( compute_tolerance_ < 0.0 )
        dserror( "Call the Compute() routine first!" );
      if ( compute_tolerance_ > 1.0 )
        dserror( "The compute_tolerance seems not trustworthy, "
            "properly the Compute() routine failed in some way. "
            "( tol = %e )", compute_tolerance_ );
      return compute_tolerance_;
    }

    /*! \brief Scaling and shifting of the input nodal positions and the location
     *         of the point \c px_
     *
     *  The corner points of the side and the location of the point
     *  are converted into the corresponding range of the natural coordinate
     *  system of the element. This leads to a  proper initial choice for
     *  the Newton scheme to find the local coordinate system. */
    virtual void ScaleAndShift()
    {
      GetElementScale<probdim>( xyze_, scale_ );

      px_.Scale( 1.0 / scale_ );
      xyze_.Scale( 1.0 / scale_ );

      GetElementShift<probdim>( xyze_, shift_ );

      for ( unsigned i=0; i<numNodesElement; ++i )
      {
        LINALG::Matrix<probdim,1> x1( &xyze_( 0, i ),  true );
        x1.Update( -1, shift_, 1 );
      }
      px_.Update( -1, shift_, 1 );
    }

  protected:
    /// @name class internal variables
    /// @{

    /// given point location vector (scaled and shifted after constructor call)
    LINALG::Matrix<probdim,1> px_;

    /// spatial nodal positions (scaled and shifted after constructor call)
    LINALG::Matrix<probdim, numNodesElement> xyze_;

    /** initial scaling of the spatial nodal coordinates \c xyze_
     *  and the location vector \c px_ */
    double scale_;

    /** initial shifting of the spatial nodal coordinates \c xyze_
     *  and the location vector \c px_ */
    LINALG::Matrix<probdim,1> shift_;

    /// local coordinates corresponding to \c px_
    LINALG::Matrix<dim,1> xsi_;

    /// computation status of the position calculation
    enum Status pos_status_;

    /** contains the tolerance used for the internal Newton method
     *  ( see the Compute() routines ) */
    double compute_tolerance_;

    /// @}

  }; // class Position

  /*----------------------------------------------------------------------------*/
  /** \brief class for the position computation in the standard/non-embedded case
   *
   *  \author hiermeier \date 08/16 */
  template <unsigned probdim,
      DRT::Element::DiscretizationType eletype,
      unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
  class ComputePosition : public PositionGeneric<probdim,eletype>
  {
  public:
    /// constructor
    ComputePosition (
        const LINALG::Matrix<probdim,numNodesElement> & xyze,
        const LINALG::Matrix<probdim,1> & xyz )
        : PositionGeneric<probdim,eletype>(xyze,xyz)
    {};

    /// destructor
    virtual ~ComputePosition() {};

    /* Don't call any of these methods directly! Use the base class public
     * methods, instead. */
  protected:

    virtual void Distance( double * distance ) const
    {
      dserror("Unsupported for the standard case!");
    }

    /** \brief Compute method for the standard case with user defined tolerance
     *
     *  Check whether the n-dimensional point lies inside the n-dimensional
     *  element.
     *
     *  \param Tol (in) : User defined tolerance for the WithinLimits check. */
    virtual bool Compute(const double& Tol)
    {
//      KERNEL::DebugComputePosition<probdim,eletype> cp( this->xsi_ );
      KERNEL::ComputePosition<probdim,eletype> cp( this->xsi_ );
      this->pos_status_ =
          ( cp( this->xyze_, this->px_ ) ? Position::position_valid : Position::position_invalid );
      this->compute_tolerance_ = cp.GetTolerance();

      if ( this->pos_status_ == Position::position_valid )
        return WithinLimitsTol(Tol);

      return false;
    }

    virtual bool IsGivenPointWithinElement()
    {
      return Position::Compute();
    }

    virtual bool WithinLimitsTol(const double& Tol) const
    {
      return KERNEL::WithinLimits<eletype>(this->xsi_,Tol);
    }
  }; // class ComputePosition

  /*----------------------------------------------------------------------------*/
  /** \brief class for the position computation in the embedded case
   *
   *  \author hiermeier \date 08/16 */
  template <unsigned probdim,
      DRT::Element::DiscretizationType eletype,
      unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim>
  class ComputeEmbeddedPosition : public PositionGeneric<probdim,eletype>
  {
  public:
    /// constructor
    ComputeEmbeddedPosition (
        const LINALG::Matrix<probdim,numNodesElement> & xyze,
        const LINALG::Matrix<probdim,1> & xyz )
        : PositionGeneric<probdim,eletype>(xyze,xyz)
    {
      this->xsi_.SetView(xsi_aug_.A());
      ConstructBoundingBox();
    };

    /// destructor
    virtual ~ComputeEmbeddedPosition() {};

    /* Don't call any of these methods directly! Use the base class public
     * methods, instead. */
  protected:
    /** \brief Default Compute method for the embedded case
     *
     *  \param Tol (in)        : User defined tolerance for the WithinLimits check.
     *  \param allow_dist (in) : If TRUE, the method allows an offset in normal direction,
     *                           otherwise the higher-dimensional point \c px_ has to lie
     *                           on the lower dimensional element. */
    virtual bool Compute( const double& tol )
    { return Compute(tol,false); }
    virtual bool Compute(const double& Tol, const bool& allow_dist);

    virtual bool IsGivenPointWithinElement();

    /*! \brief Return the perpendicular distance between given point to the side */
    virtual double Distance() const
    {
      if ( this->pos_status_ < Position::position_distance_valid )
      {
        std::ostringstream msg;
        msg << "Neither the position nor the distance value is valid! "
            "( Position::Status = "
            << Position::Status2String( this->pos_status_ ) << " )";
        run_time_error( msg.str() );
      }

      switch ( probdim - dim )
      {
        case 1:
          return xsi_aug_( dim, 0 );
        default:
          dserror( "A scalar signed distance value is not available!" );
          exit( EXIT_FAILURE );
      }
    }

    virtual void Distance( double * distance ) const
    {
      if ( this->pos_status_ < Position::position_distance_valid )
      {
        std::ostringstream msg;
        msg << "Neither the position nor the distance value is valid! "
            "( Position::Status = "
            << Position::Status2String( this->pos_status_ ) << " )";
        run_time_error( msg.str() );
      }

      std::copy( xsi_aug_.A() + dim, xsi_aug_.A() + probdim, distance );
    }

    virtual bool WithinLimitsTol(const double& Tol, const bool& allow_dist) const
    {
      double tol2 = Tol;
      double tol_xyze = this->xyze_.NormInf();
      double tol_px = this->px_.NormInf();
      // choose the weaker tolerance
      if (tol_xyze > tol_px)
        tol2 *= tol_xyze;
      else
        tol2 *= tol_px;

      return KERNEL::WithinLimitsEmbeddedManifold<probdim,eletype>(xsi_aug_,
          Tol,allow_dist,tol2);
    }

  private:
    /** \brief Construct bounding box over the given embedded element
     *  (after scaling) */
    void ConstructBoundingBox();

    //not used in the actual version 02/17
#if (0)

    /*! \brief Compute the perpendicular distance between point p3 and the line joining p1 and p2
     *
     *                  \f$(p_3)\f$ x_____
     *                             /:    ^
     *                            / :    |  distance d
     *                           /  :    v
     *               f$(p_1)\f$ o---+-------------o \f$(p_2)\f$
     *                           (\f$\hat{p_3}\f$) */
    bool LinePointDistance( const LINALG::Matrix<probdim,1> & p1,
                            const LINALG::Matrix<probdim,1> & p2,
                            const LINALG::Matrix<probdim,1> & p3,
                            double & distance )
    {
      // vector between p1 and p2
      LINALG::Matrix<probdim,1> n;
      n.Update( 1, p2, -1, p1, 0 );

      // vector between p1 and p3
      LINALG::Matrix<probdim,1> m;
      m.Update( 1, p3, -1, p1, 0 );

      // calculate the scaled distance factor from point p1 to the projection
      // point \hat{p3}.
      double a = m.Dot( n );
      double b = n.Dot(n);
      double c = a/b;

      /* If c is smaller than 0.0, the projection will be on the left side
       * of p1 (see sketch). If it is larger than 1.0 the projection lies
       * on the right side of p2. */
      if ( c < 0 or c > 1)
        return false;

      // calculate the distance from p3 to the projection point \hat{p3}
      LINALG::Matrix<probdim,1> d;
      d.Update( 1.0, p1, c, n, 0.0 );
      d.Update( -1.0, p3, 1.0 );

      distance = d.Norm2();
      return true;
    }
#endif

  private:
    /** local coordinates corresponding to \c px_, the last n = ( probdim - dim )
     *  coordinates are the distance values (embedded case only) */
    LINALG::Matrix<probdim,1> xsi_aug_;

    //! Bounding box over the given embedded element (after scaling is performed)
    Teuchos::RCP<BoundingBox> bbside_;
  }; // class ComputeEmbeddedPosition

  /*----------------------------------------------------------------------------*/
  /** \class Position factory
   *
   *  \author hiermeier \date 11/16 */
  class PositionFactory
  {
  public:
    /// constructor
    PositionFactory();

    /// get the current problem dimension
    unsigned ProbDim() const
    {
      return probdim_;
    }

    /** \brief build variant #1
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param point   (in) : Given global point */
    Teuchos::RCP<Position> CreatePosition(
        const Element & element,
        const Point & point ) const;

    /** \brief build variant #2
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param xyz     (in) : Global coordinates of the given point. */
    Teuchos::RCP<Position> CreatePosition(
        const Element & element,
        const double * xyz ) const;

    /** \brief build variant #3
     *
     *  \param xyze    (in) : Global nodal positions of the element.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type. */
    Teuchos::RCP<Position> CreatePosition(
        const double * xyze,
        const double * xyz,
        const DRT::Element::DiscretizationType & distype ) const;

    /** \brief build variant #4
     *
     *  \param nodes   (in) : Nodes of the element we want to check.
     *  \param xyz     (in) : Global coordinates of the given point.
     *  \param distype (in) : element discretization type. (optional)*/
    Teuchos::RCP<Position> CreatePosition(
        const std::vector<Node*> nodes,
        const double * xyz,
        DRT::Element::DiscretizationType distype = DRT::Element::dis_none) const;

  private:
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #1                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{

    /** build variant #1
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param point   (in) : Given global point
     *
     *  \author hiermeier \date 08/16 */
#if __cplusplus < 201103L
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype
        >
    Position * BuildPosition(
        const Element & element,
        const Point & point ) const
    {
      const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      const unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
      const bool isembedded = ( probdim > dim );
#else
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype,
        unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
        bool isembedded = (probdim>dim)
        >
    Position * BuildPosition(
        const Element & element,
        const Point & point ) const
    {
#endif
      LINALG::Matrix<probdim,numNodesElement> xyze;
      element.Coordinates( xyze );
      LINALG::Matrix<probdim,1> px;
      point.Coordinates( px.A() );

      Position * position_ptr = NULL;
      if ( isembedded )
        position_ptr = new ComputeEmbeddedPosition<probdim,eletype>( xyze, px );
      else
        position_ptr = new ComputePosition<probdim,eletype>( xyze, px );

      return position_ptr;
    }

    /// concrete create variant #1
    template< DRT::Element::DiscretizationType eletype>
    Position * CreateConcretePosition(
        const Element & element,
        const Point & point) const
    {
      const int dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      if ( dim > probdim_ )
        dserror( "The element dimension is larger than the problem dimension! \n"
            "dim = %d, probdim = %d", dim, probdim_ );

      switch (probdim_)
      {
        case 2:
          return BuildPosition<2,eletype>( element, point );
        case 3:
          return BuildPosition<3,eletype>( element, point );
        default:
          dserror("Unsupported problem dimension! (probdim = %d)", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #2                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{
  public:
    /** \brief build variant #2
     *
     *  \param element (in) : Cut element. We check whether the given
     *                        point lies inside / on it.
     *  \param xyz     (in) : Global coordinates of the given point.
     *
     *  \author hiermeier \date 08/16 */
#if __cplusplus < 201103L
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype
        >
    static Teuchos::RCP<Position> BuildPosition(
        const Element & element,
        const LINALG::Matrix<probdim,1> & xyz )
    {
      const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      const unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
      const bool isembedded = ( probdim > dim );
#else
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype,
        unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
        bool isembedded = (probdim>dim)
        >
    static Teuchos::RCP<Position> BuildPosition(
        const Element & element,
        const LINALG::Matrix<probdim,1> & xyz )
    {
#endif
      LINALG::Matrix<probdim,numNodesElement> xyze;
      element.Coordinates( xyze );

      Teuchos::RCP<Position> position_ptr = Teuchos::null;
      if ( isembedded )
        position_ptr = Teuchos::rcp(new ComputeEmbeddedPosition<probdim,eletype>(xyze,xyz));
      else
        position_ptr = Teuchos::rcp(new ComputePosition<probdim,eletype>(xyze,xyz));

      return position_ptr;
    }

  private:
    /// concrete create variant #2
    template< DRT::Element::DiscretizationType eletype>
    Teuchos::RCP<Position> CreateConcretePosition(
        const Element & element,
        const double * xyz) const
    {
      const int dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      if ( dim > probdim_ )
        dserror( "The element dimension is larger than the problem dimension! \n"
            "dim = %d, probdim = %d", dim, probdim_ );

      switch (probdim_)
      {
        case 2:
        {
          LINALG::Matrix<2,1> xyz_mat( xyz, true );
          return BuildPosition<2,eletype>( element, xyz_mat );
        }
        case 3:
        {
          LINALG::Matrix<3,1> xyz_mat( xyz, true );
          return BuildPosition<3,eletype>( element, xyz_mat );
        }
        default:
          dserror("Unsupported problem dimension! (probdim = %d)", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #3                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{
  public:
    /** \brief build variant #3
     *
     *  \param xyze    (in) : Global nodal positions of the element.
     *  \param xyz     (in) : Global coordinates of the given point.
     *
     *  \author hiermeier \date 08/16 */
#if __cplusplus < 201103L
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype
        >
    static Teuchos::RCP<Position> BuildPosition(
        const LINALG::Matrix<probdim, DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement> & xyze,
        const LINALG::Matrix<probdim, 1> & xyz )
    {
      const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      const bool isembedded = ( probdim > dim );
#else
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype,
        unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
        bool isembedded = (probdim>dim)
        >
    static Teuchos::RCP<Position> BuildPosition(
        const LINALG::Matrix<probdim, numNodesElement> & xyze,
        const LINALG::Matrix<probdim, 1> & xyz )
    {
#endif
      Teuchos::RCP<Position> position_ptr = Teuchos::null;
      if ( isembedded )
        position_ptr = Teuchos::rcp(new ComputeEmbeddedPosition<probdim,eletype>(xyze,xyz));
      else
        position_ptr = Teuchos::rcp(new ComputePosition<probdim,eletype>(xyze,xyz));

      return position_ptr;
    }

  private:
    /// concrete create variant #3
    template< DRT::Element::DiscretizationType eletype>
    Teuchos::RCP<Position> CreateConcretePosition(
        const double * xyze,
        const double * xyz ) const
    {
      const int dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      if ( dim > probdim_ )
        dserror( "The element dimension is larger than the problem dimension! \n"
            "dim = %d, probdim = %d", dim, probdim_ );

      const unsigned num_nodes_ele = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
      switch (probdim_)
      {
        case 2:
        {
          LINALG::Matrix<2, num_nodes_ele> xyze_mat( xyze, true );
          LINALG::Matrix<2, 1> xyz_mat( xyz, true );
          return BuildPosition<2,eletype>( xyze_mat, xyz_mat );
        }
        case 3:
        {
          LINALG::Matrix<3, num_nodes_ele> xyze_mat( xyze, true );
          LINALG::Matrix<3, 1> xyz_mat( xyz, true );
          return BuildPosition<3,eletype>( xyze_mat, xyz_mat );
        }
        default:
          dserror("Unsupported problem dimension! (probdim = %d)", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}
    /*----------------------------------------------------------------------------*/
    /** @name VARIANT #4                                                          */
    /*----------------------------------------------------------------------------*/
    /// @{
  public:
    /** \brief build variant #4
     *
     *  \param nodes (in) : Nodes of the element we want to check.
     *  \param xyz   (in) : Global coordinates of the given point.
     *
     *  \author hiermeier \date 08/16 */
#if __cplusplus < 201103L
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype
        >
    static Teuchos::RCP<Position> BuildPosition(
        const std::vector<Node*> nodes,
        const LINALG::Matrix<probdim,1> & xyz )
    {
      const unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      const unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement;
      const bool isembedded = ( probdim > dim );
#else
    template <unsigned probdim,
        DRT::Element::DiscretizationType eletype,
        unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
        unsigned numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement,
        bool isembedded = (probdim>dim)
        >
    static Teuchos::RCP<Position> BuildPosition(
        const std::vector<Node*> nodes,
        const LINALG::Matrix<probdim,1> & xyz )
    {
#endif
      if( not (nodes.size() == numNodesElement) )
        dserror( "node number for this element is not correct\n"
            "numNodesElement = %d, nodes.size() = %d", numNodesElement, nodes.size() );

      LINALG::Matrix<probdim,numNodesElement> xyze;
      for ( unsigned i=0; i<numNodesElement; ++i )
      {
        Node * n = nodes[i];
        n->Coordinates( &xyze( 0, i ) );
      }

      Teuchos::RCP<Position> position_ptr = Teuchos::null;
      if ( isembedded )
        position_ptr = Teuchos::rcp( new ComputeEmbeddedPosition<probdim,eletype>(xyze,xyz) );
      else
        position_ptr = Teuchos::rcp( new ComputePosition<probdim,eletype>(xyze,xyz) );

      return position_ptr;
    }

  private:
    /// concrete create variant #4
    template< DRT::Element::DiscretizationType eletype>
    Teuchos::RCP<Position> CreateConcretePosition(
        const std::vector<Node*> nodes,
        const double * xyz) const
    {
      const int dim = DRT::UTILS::DisTypeToDim<eletype>::dim;
      if ( dim > probdim_ )
        dserror( "The element dimension is larger than the problem dimension! \n"
            "dim = %d, probdim = %d", dim, probdim_ );

      switch ( probdim_ )
      {
        case 2:
        {
          LINALG::Matrix<2,1> xyz_mat( xyz, true );
          return BuildPosition<2,eletype>( nodes, xyz_mat );
        }
        case 3:
        {
          LINALG::Matrix<3,1> xyz_mat( xyz, true );
          return BuildPosition<3,eletype>( nodes, xyz_mat );
        }
        default:
          dserror("Unsupported problem dimension! (probdim = %d)", probdim_);
          exit(EXIT_FAILURE);
      }

      exit(EXIT_FAILURE);
    }

    /// @}

    /// problem dimension
    unsigned probdim_;
  }; // class PositionFactory


} // namespace CUT
} // namespace GEO


#endif
