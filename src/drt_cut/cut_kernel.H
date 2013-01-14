#ifndef CUT_KERNEL_H
#define CUT_KERNEL_H

#include "cut_tolerance.H"
#include "cut_utils.H"

#include "../linalg/linalg_gauss.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

namespace GEO
{
  namespace CUT
  {
/*!
\brief Core geometrical operations like finding the cut points, distance calculations takes place here
 */
    namespace KERNEL
    {

      void FindCornerPoints( const std::vector<Point*> & points,
                             std::vector<Point*> & corner_points );

      unsigned FindNextCornerPoint( const std::vector<Point*> & points,
                                    LINALG::Matrix<3,1> & x1,
                                    LINALG::Matrix<3,1> & x2,
                                    LINALG::Matrix<3,1> & x3,
                                    LINALG::Matrix<3,1> & b1,
                                    LINALG::Matrix<3,1> & b2,
                                    LINALG::Matrix<3,1> & b3,
                                    unsigned i );

      inline bool IsValidTri3( const std::vector<Point*> & points )
      { return points.size()==3; }

      bool IsValidQuad4( const std::vector<Point*> & points );

      DRT::Element::DiscretizationType CalculateShape( const std::vector<Point*> & points,
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


      template <DRT::Element::DiscretizationType sidetype>
      bool AtEdge( const LINALG::Matrix<3,1> & xsi )
      {
        switch ( sidetype )
        {
        case DRT::Element::quad4:
        case DRT::Element::quad8:
        case DRT::Element::quad9:
        {
          return ( fabs( xsi( 0 )+1 )<TOLERANCE or fabs( xsi( 1 )+1 )<TOLERANCE or
                   fabs( xsi( 0 )-1 )<TOLERANCE or fabs( xsi( 1 )-1 )<TOLERANCE );
        }
        case DRT::Element::tri3:
        case DRT::Element::tri6:
        {
          return ( fabs( xsi( 0 )+0 )<TOLERANCE or fabs( xsi( 1 )+0 )<TOLERANCE or
                   fabs( xsi( 1 )+xsi( 0 )-1 )<TOLERANCE );
        }
        default:
          throw std::runtime_error( "unsupported side type" );
        }
      }


      template <DRT::Element::DiscretizationType sidetype>
      bool SurfaceWithinLimits( const LINALG::Matrix<3,1> & xsi )
      {
        bool intersects = false;
        switch ( sidetype )
        {
        case DRT::Element::quad4:
        case DRT::Element::quad8:
        case DRT::Element::quad9:
        {
          intersects = ( xsi( 0 ) >= -1-TOLERANCE and xsi( 1 ) >= -1-TOLERANCE and
                         xsi( 0 ) <=  1+TOLERANCE and xsi( 1 ) <=  1+TOLERANCE );
          break;
        }
        case DRT::Element::tri3:
        case DRT::Element::tri6:
        {
          intersects = ( xsi( 0 ) >= 0-TOLERANCE and xsi( 1 ) >= 0-TOLERANCE and
                         xsi( 0 ) <= 1+TOLERANCE and xsi( 1 ) <= 1+TOLERANCE and
                         xsi( 1 ) <= ( 1 - xsi( 0 ) )+TOLERANCE );
          break;
        }
        default:
          throw std::runtime_error( "unsupported side type" );
        }
        return intersects;
      }


      inline bool LineWithinLimits( const LINALG::Matrix<3,1> & xsi )
      {
        return xsi( 2 ) >= -1-TOLERANCE and xsi( 2 ) <=  1+TOLERANCE;
      }


      template <DRT::Element::DiscretizationType elementtype>
      bool WithinLimits( const LINALG::Matrix<3,1> & xsi )
      {
        switch ( elementtype )
        {
        case DRT::Element::hex8:
        case DRT::Element::hex20:
        case DRT::Element::hex27:
          return ( xsi( 0 ) >= -1-TOLERANCE and xsi( 1 ) >= -1-TOLERANCE and xsi( 2 ) >= -1-TOLERANCE and
                   xsi( 0 ) <=  1+TOLERANCE and xsi( 1 ) <=  1+TOLERANCE and xsi( 2 ) <=  1+TOLERANCE );
        case DRT::Element::tet4:
        case DRT::Element::tet10:
          return ( xsi( 0 ) >=   -TOLERANCE and xsi( 1 ) >=   -TOLERANCE and xsi( 2 ) >=   -TOLERANCE and
                   xsi( 0 )+xsi( 1 )+xsi( 2 ) <= 1+TOLERANCE );
        case DRT::Element::pyramid5:
          return ( xsi( 0 ) >= -1-TOLERANCE and xsi( 1 ) >= -1-TOLERANCE and xsi( 2 ) >=   -TOLERANCE and
                   xsi( 0 ) <=  1+TOLERANCE and xsi( 1 ) <=  1+TOLERANCE and
                   ( fabs( xsi( 0 ) ) > fabs( xsi( 1 ) ) ?
                     ( xsi( 2 )+fabs( xsi( 0 ) ) <= 1+TOLERANCE ) :
                     ( xsi( 2 )+fabs( xsi( 1 ) ) <= 1+TOLERANCE ) )
            );
        case DRT::Element::wedge6:
        case DRT::Element::wedge15:
          return ( xsi( 0 ) >=   -TOLERANCE and xsi( 1 ) >=   -TOLERANCE and xsi( 2 ) >= -1-TOLERANCE and
                   xsi( 2 ) <=  1+TOLERANCE and
                   xsi( 0 )+xsi( 1 ) <= 1+TOLERANCE );
        default:
          throw std::runtime_error( "unsupported element type" );
        }
        return false;
      }

      /*!
      \brief Check whether three points are on the same line
       */
      bool IsOnLine( Point*& pt1, Point*& pt2, Point*& pt3 );

      /*!
      \brief Check whether the polygon defined by the set of points is convex
       */
      std::vector<int> CheckConvexity( const std::vector<Point*>& ptlist,
                                       std::string& geoType,
                                       bool InSplit = true );

      std::vector<double> EqnPlanePolygon( const std::vector<Point*>& ptlist );

      /*!
      \brief Find the equation of plane that contains these non-collinear points
      It must be noted while using this function to find equation of facets,
             none of these 3 points must be a reflex (concave) point
       */
      std::vector<double> EqnPlane( Point* & pt1, Point* & pt2, Point* & pt3 );

      /*!
      \brief Check whether the point "check" is inside the triangle formed by tri
      */
      bool PtInsideTriangle( std::vector<Point*> tri, Point* check );

      /*!
      \brief Check whether the point "check" is inside the Quad
      */
      bool PtInsideQuad( std::vector<Point*> quad, Point* check );

      /*!
      \brief Return true if the points of the polygon are ordered clockwise
       */
      bool IsClockwiseOrderedPolygon( std::vector<Point*>polyPoints, std::string& projType );

      /*!
      \brief If more than two points are on a line, all points except the end points are deleted
      */
      void DeleteInlinePts( std::vector<Point*>& poly );


      /// generic Newton algorithm
      template <class Strategy, int maxiter=20>
      class NewtonSolve : public Strategy
      {
      public:

        NewtonSolve( LINALG::Matrix<3,1> & xsi )
          : Strategy( xsi )
        {
        }

        bool Solve()
        {
          this->SetupSolve();

          for ( int iter = 0; iter < maxiter; ++iter )
          {
            this->SetupStep( iter );

            if ( this->TestConverged( iter ) )
            {
              return true;
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

      };

      /// empty strategy for generic Newton algorithm
      class EmptyNewtonStrategy
      {
      public:

        void SetupSolve() {}

        void SetupStep( int iter ) {}

        bool TestConverged( int iter )
        {
          return false;
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
      };

      /// Debug helper strategy for generic Newton algorithm
      template <class Strategy>
      class DebugNewtonStrategy : public Strategy
      {
      public:

        DebugNewtonStrategy( LINALG::Matrix<3,1> & xsi )
          : Strategy( xsi )
        {
        }

        void SetupSolve()
        {
          std::cout << "SetupSolve()\n";
          Strategy::SetupSolve();
        }

        void SetupStep( int iter )
        {
          std::cout << "SetupStep( " << iter << " )\n";
          Strategy::SetupStep( iter );
        }

        bool TestConverged( int iter )
        {
          bool res = Strategy::TestConverged( iter );
          std::cout << "TestConverged( " << iter << " ) = " << res << "\n";
          return res;
        }

        bool LinearSolve( int iter )
        {
          bool res = Strategy::LinearSolve( iter );
          std::cout << "LinearSolve( " << iter << " ) = " << res << "\n";
          return res;
        }

        bool Update( int iter )
        {
          bool res = Strategy::Update( iter );
          std::cout << "Update( " << iter << " ) = " << res << "\n";
          return res;
        }

        bool NewtonFailed()
        {
          bool res = Strategy::NewtonFailed();
          std::cout << "NewtonFailed() = " << res << "\n";
          return res;
        }
      };

      /// strategy for position of point within element
      template <bool debug,
                DRT::Element::DiscretizationType elementtype,
                int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
      >
      class ComputePositionStrategy : public EmptyNewtonStrategy
      {
      public:

        ComputePositionStrategy( LINALG::Matrix<3,1> & xsi )
          : xsi_( xsi ),
            xyze_( NULL ),
            px_( NULL )
        {
        }

        void Setup( const LINALG::Matrix<3, numNodesElement> & xyze,
                    const LINALG::Matrix<3,1> & px )
        {
          xyze_ = &xyze;
          px_ = &px;
        }

        const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

        void SetupSolve()
        {
          xsi_ = 0;
        }

        void SetupStep( int iter )
        {
          DRT::UTILS::shape_function_3D( funct_, xsi_( 0 ), xsi_( 1 ), xsi_( 2 ), elementtype );

          b_ = *px_;
          b_.Multiply( -1, *xyze_, funct_, 1 );
        }

        bool TestConverged( int iter )
        {
          double residual = b_.Norm2();
          return fabs( residual ) < TOLERANCE;
        }

        bool LinearSolve( int iter )
        {
          DRT::UTILS::shape_function_3D_deriv1( deriv1_, xsi_( 0 ), xsi_( 1 ), xsi_( 2 ), elementtype );

          A_.MultiplyNT( *xyze_, deriv1_ );

          dx_ = 0;
          double det = LINALG::gaussElimination<true, 3>( A_, b_, dx_ );
          return fabs( det ) >= LINSOLVETOL;
        }

        bool Update( int iter )
        {
          xsi_ += dx_;
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

      private:

        LINALG::Matrix<3,1> & xsi_;

        const LINALG::Matrix<3, numNodesElement> * xyze_;
        const LINALG::Matrix<3,1> * px_;

        LINALG::Matrix<numNodesElement,1> funct_;
        LINALG::Matrix<3,numNodesElement> deriv1_;

        LINALG::Matrix<3,3> A_;
        LINALG::Matrix<3,1> b_;
        LINALG::Matrix<3,1> dx_;
      };


      /// distance between side and point
      template <class Strategy,
                DRT::Element::DiscretizationType elementtype,
                int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
      >
      class GenericComputePosition : Strategy
      {
      public:

        GenericComputePosition( LINALG::Matrix<3,1> & xsi )
          : Strategy( xsi )
        {
        }

        bool operator()( const LINALG::Matrix<3, numNodesElement> & xyze,
                         const LINALG::Matrix<3,1> & px )
        {
          this->Setup( xyze, px );
          if ( this->Solve() )
          {
            return true;
          }
          return false;
        }

        const LINALG::Matrix<3,1> & LocalCoordinates()
        {
          return Strategy::LocalCoordinates();
        }
      };

      template <DRT::Element::DiscretizationType elementtype,
                int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
      >
      class ComputePosition : public GenericComputePosition<NewtonSolve<ComputePositionStrategy<false,elementtype,numNodesElement> >,
                                                            elementtype>
      {
      public:

        ComputePosition( LINALG::Matrix<3,1> & xsi )
          : GenericComputePosition<NewtonSolve<ComputePositionStrategy<false,elementtype,numNodesElement> >,
                                   elementtype>( xsi )
        {
        }
      };

      template <DRT::Element::DiscretizationType elementtype,
                int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement
      >
      class DebugComputePosition : public GenericComputePosition<NewtonSolve<DebugNewtonStrategy<ComputePositionStrategy<true,elementtype,numNodesElement> > >,
                                                                 elementtype>
      {
      public:

        DebugComputePosition( LINALG::Matrix<3,1> & xsi )
          : GenericComputePosition<NewtonSolve<DebugNewtonStrategy<ComputePositionStrategy<true,elementtype,numNodesElement> > >,
                                   elementtype>( xsi )
        {
        }
      };


      /// strategy for distance between side and point
      template <bool debug,
                DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class ComputeDistanceStrategy : public EmptyNewtonStrategy
      {
      public:

        ComputeDistanceStrategy( LINALG::Matrix<3,1> & xsi )
          : xsi_( xsi ),
            xyze_surface_( NULL ),
            px_( NULL )
        {
        }

        void Setup( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                    const LINALG::Matrix<3,1> & px )
        {
          xyze_surface_ = &xyze_surface;
          px_ = &px;
        }

        double Distance() { return fabs( xsi_( 2 ) ); }

        const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

        void SetupSolve()
        {
          xsi_ = 0;
        }

        void SetupStep( int iter )
        {
          DRT::UTILS::shape_function_2D( surfaceFunct_, xsi_( 0 ), xsi_( 1 ), sidetype );
          DRT::UTILS::shape_function_2D_deriv1( surfaceDeriv1_, xsi_( 0 ), xsi_( 1 ), sidetype );
          DRT::UTILS::shape_function_2D_deriv2( surfaceDeriv2_, xsi_( 0 ), xsi_( 1 ), sidetype );

          DistanceSystem( surfaceDeriv2_, surfaceDeriv1_, surfaceFunct_,
                          *xyze_surface_, *px_,
                          xsi_( 2 ),
                          A_, B_, C_, cross_, b_ );
        }

        bool TestConverged( int iter )
        {
          if ( debug )
          {
            std::cout << "  " << SurfaceWithinLimits<sidetype>( xsi_ )
                      << "  " << xsi_.Norm2()
                      << "  " << b_.Norm2()
                      << "  " << ( b_.Norm2() / xsi_.Norm2() )
                      << "\n"
                      << "  " << xsi_
                      << "  " << b_;
          }

          if ( iter > 0 and not SurfaceWithinLimits<sidetype>( xsi_ ) )
          {
            // if we are outside our side, just go
            return true;
          }

          double solres = xsi_.Norm2();
          double residual = b_.Norm2();

          const double tol = LINSOLVETOL; // / scale_;

          return ( ( solres >= std::numeric_limits<double>::min() and
                     residual / solres < tol ) or
                   residual < tol );
        }

        bool LinearSolve( int iter )
        {
          dx_ = 0;
          double det = LINALG::gaussElimination<true, 3>( A_, b_, dx_ );
          return det != 0;
        }

        bool Update( int iter )
        {
          xsi_ += dx_;
          return true;
        }

        bool NewtonFailed()
        {
#if 1
          return true;
#else
          std::stringstream str;
          str << "ComputeDistance: Newton scheme did not converge:\n"
              << std::setprecision( 16 )
              << ( *xyze_surface_ )
              << ( *px_ )
              << xsi_
            ;
          throw std::runtime_error( str.str() );
#endif
        }

      private:

        void DistanceSystem( const LINALG::Matrix<3,numNodesSurface> & surfaceDeriv2,
                             const LINALG::Matrix<2,numNodesSurface> & surfaceDeriv1,
                             const LINALG::Matrix<numNodesSurface,1> & surface,
                             const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                             const LINALG::Matrix<3,1> & px,
                             double distance,
                             LINALG::Matrix<3,3> & A,
                             LINALG::Matrix<3,3> & B,
                             LINALG::Matrix<3,3> & C,
                             LINALG::Matrix<3,1> & cross,
                             LINALG::Matrix<3,1> & b )
        {
          b = px;
          b.Multiply( -1, xyze_surface, surface, 1 );

          A = 0;

          for ( int inode=0; inode<numNodesSurface; ++inode )
            for( int isd=0; isd<3; ++isd )
            {
              A(isd,0) += xyze_surface(isd,inode) * surfaceDeriv1(0,inode);
              A(isd,1) += xyze_surface(isd,inode) * surfaceDeriv1(1,inode);
            }

          // cross product to get the normal at the point
          cross( 0 ) = A( 1, 0 )*A( 2, 1 ) - A( 2, 0 )*A( 1, 1 );
          cross( 1 ) = A( 2, 0 )*A( 0, 1 ) - A( 0, 0 )*A( 2, 1 );
          cross( 2 ) = A( 0, 0 )*A( 1, 1 ) - A( 1, 0 )*A( 0, 1 );

          // norm normal direction
          double cnorm = cross.Norm2();
          double fact1 = 1.0;

          if(fabs(cnorm) < 1e-14)
          {
            fact1=1.0;
          }
          else
          {
            fact1 = 1./cnorm;
          }

          //double fact2 = -0.5/( cnorm*sqrt( cnorm ) );
          cross.Scale( fact1 );

          C = 0;

          for ( int inode=0; inode<numNodesSurface; ++inode )
            for( int isd=0; isd<3; ++isd )
            {
              C(isd,0) += xyze_surface(isd,inode) * surfaceDeriv2(0,inode);
              C(isd,1) += xyze_surface(isd,inode) * surfaceDeriv2(1,inode);
              C(isd,2) += xyze_surface(isd,inode) * surfaceDeriv2(2,inode);
            }

          b( 0 ) -= cross( 0 )*distance;
          b( 1 ) -= cross( 1 )*distance;
          b( 2 ) -= cross( 2 )*distance;

          A( 0, 2 ) = cross( 0 );
          A( 1, 2 ) = cross( 1 );
          A( 2, 2 ) = cross( 2 );

          B = 0;

          B( 0, 0 ) = ( C( 1, 0 )*A( 2, 1 )+A( 1, 0 )*C( 2, 2 ) - ( C( 2, 0 )*A( 1, 1 )+A( 2, 0 )*C( 1, 2 ) ) )*distance;
          B( 0, 1 ) = ( C( 1, 2 )*A( 2, 1 )+A( 1, 0 )*C( 2, 1 ) - ( C( 2, 2 )*A( 1, 1 )+A( 2, 0 )*C( 1, 1 ) ) )*distance;

          B( 1, 0 ) = ( C( 2, 0 )*A( 0, 1 )+A( 2, 0 )*C( 0, 2 ) - ( C( 0, 0 )*A( 2, 1 )+A( 0, 0 )*C( 2, 2 ) ) )*distance;
          B( 1, 1 ) = ( C( 2, 2 )*A( 0, 1 )+A( 2, 0 )*C( 0, 1 ) - ( C( 0, 2 )*A( 2, 1 )+A( 0, 0 )*C( 2, 1 ) ) )*distance;

          B( 2, 0 ) = ( C( 0, 0 )*A( 1, 1 )+A( 0, 0 )*C( 1, 2 ) - ( C( 1, 0 )*A( 0, 1 )+A( 1, 0 )*C( 0, 2 ) ) )*distance;
          B( 2, 1 ) = ( C( 0, 2 )*A( 1, 1 )+A( 0, 0 )*C( 1, 1 ) - ( C( 1, 2 )*A( 0, 1 )+A( 1, 0 )*C( 0, 1 ) ) )*distance;

          // product rule for normed cross product

          B( 0, 0 ) *= fact1; //+ fact2 * A( 0, 2 );
          B( 0, 1 ) *= fact1; //+ fact2 * A( 0, 2 );

          B( 1, 0 ) *= fact1; //+ fact2 * A( 1, 2 );
          B( 1, 1 ) *= fact1; //+ fact2 * A( 1, 2 );

          B( 2, 0 ) *= fact1; //+ fact2 * A( 2, 2 );
          B( 2, 1 ) *= fact1; //+ fact2 * A( 2, 2 );

          A.Update( 1, B, 1 );
        }

        LINALG::Matrix<3,1> & xsi_;

        const LINALG::Matrix<3,numNodesSurface> * xyze_surface_;
        const LINALG::Matrix<3,1> * px_;

        LINALG::Matrix<numNodesSurface,1> surfaceFunct_;
        LINALG::Matrix<2,numNodesSurface> surfaceDeriv1_;
        LINALG::Matrix<3,numNodesSurface> surfaceDeriv2_;

        LINALG::Matrix<3,3> A_;
        LINALG::Matrix<3,3> B_;
        LINALG::Matrix<3,3> C_;
        LINALG::Matrix<3,1> b_;
        LINALG::Matrix<3,1> dx_;
        LINALG::Matrix<3,1> cross_;
      };


      /// distance between side and point
      template <class Strategy,
                DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class GenericComputeDistance : Strategy
      {
      public:

        GenericComputeDistance( LINALG::Matrix<3,1> & xsi )
          : Strategy( xsi )
        {
        }

        bool operator()( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                         const LINALG::Matrix<3,1> & px,
                         double & distance )
        {
          this->Setup( xyze_surface, px );
          if ( this->Solve() )
          {
            distance = this->Distance();
            return true;
          }
          return false;
        }

        const LINALG::Matrix<3,1> & LocalCoordinates()
        {
          return Strategy::LocalCoordinates();
        }

        bool SurfaceWithinLimits()
        {
          return SurfaceWithinLimits<sidetype>( LocalCoordinates() );
        }
      };

      template <DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class ComputeDistance : public GenericComputeDistance<NewtonSolve<ComputeDistanceStrategy<false,edgetype,sidetype> >,
                                                            edgetype,sidetype>
      {
      public:

        ComputeDistance( LINALG::Matrix<3,1> & xsi )
          : GenericComputeDistance<NewtonSolve<ComputeDistanceStrategy<false,edgetype,sidetype> >,
                                   edgetype,sidetype>( xsi )
        {
        }
      };

      template <DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class DebugComputeDistance : public GenericComputeDistance<NewtonSolve<DebugNewtonStrategy<ComputeDistanceStrategy<true,edgetype,sidetype> > >,
                                                                 edgetype,sidetype>
      {
      public:

        DebugComputeDistance( LINALG::Matrix<3,1> & xsi )
          : GenericComputeDistance<NewtonSolve<DebugNewtonStrategy<ComputeDistanceStrategy<true,edgetype,sidetype> > >,
                                   edgetype,sidetype>( xsi )
        {
        }
      };


      /// strategy for intersection of side and line
      template <bool debug,
                DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class ComputeIntersectionStrategy : public EmptyNewtonStrategy
      {
      public:

        ComputeIntersectionStrategy( LINALG::Matrix<3,1> & xsi )
          :  xsi_( xsi ),
             xyze_surface_( NULL ),
             xyze_line_( NULL )
        {
        }

        void Setup( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                    const LINALG::Matrix<3,numNodesLine> & xyze_line )
        {
          xyze_surface_ = &xyze_surface;
          xyze_line_    = &xyze_line;

          if ( debug )
          {
            std::cout << std::setprecision( 15 ) << xyze_surface << xyze_line;
          }
        }

        const LINALG::Matrix<3,1> & LocalCoordinates() { return xsi_; }

        void SetupSolve()
        {
          xsi_ = 0;                // local coordinates (xi1_side, xi2_side, xi_line)
          dx_ = 0;                 // increment in local coordinates d(xi1_side, xi2_side, xi_line)
          old_b_ = 0;              // rhs of linear system
          off_count_ = 0;

          residual_ = 0;           // residual
          incr_ = 1;               // increment 2-norm || dx_ = dxi_ ||

          // Newton steps are supposed to converge within the side and within the line. Provide a
          // maximum legal step size. Larger steps are supposed to diverge.
          // REMARK: local coordinates for lines and sides are supposed to be in [-1,1], a radius of 5 should be okay
          maxincr_ = 10.0; //   Uli's version: 10 * xyze_surface_->Norm2(); (unreasonable, update is done in local coordinates, xyz is global!!!)
        }

        void SetupStep( int iter )
        {
          DRT::UTILS::shape_function_2D( surfaceFunct_, xsi_( 0 ), xsi_( 1 ), sidetype );
          DRT::UTILS::shape_function_1D( lineFunct_, xsi_( 2 ), edgetype );

          b_.Multiply( *xyze_line_, lineFunct_ );
          b_.Multiply( -1, *xyze_surface_, surfaceFunct_, 1 );
        }

        bool TestConverged( int iter )
        {
          const double tol = LINSOLVETOL; // / scale_;

          residual_ = b_.Norm2();

          if ( debug )
            std::cout << "  residual=" << std::setprecision( 15 ) << residual_
                      << "  incr=" << std::setprecision( 15 ) << incr_
                      << "\n";

          if ( residual_ < tol and incr_ < tol )
          {
            return true;
          }
          if ( iter > 0 and off_count_ == 0 )
          {
            // test for stagnation
            old_b_.Update( 1, b_, -1 );
            if ( old_b_.Norm2() < tol )
            {
              return true;
            }
          }
          old_b_ = b_;
          return false;
        }

        bool LinearSolve( int iter )
        {
          DRT::UTILS::shape_function_2D_deriv1( surfaceDeriv1_, xsi_( 0 ), xsi_( 1 ), sidetype );
          DRT::UTILS::shape_function_1D_deriv1( lineDeriv1_, xsi_( 2 ), edgetype );

          A_ = 0;

          for ( int inode=0; inode<numNodesSurface; ++inode )
            for( int isd=0; isd<3; ++isd )
            {
              A_(isd,0) += ( *xyze_surface_ )(isd,inode) * surfaceDeriv1_(0,inode);
              A_(isd,1) += ( *xyze_surface_ )(isd,inode) * surfaceDeriv1_(1,inode);
            }

          for ( int inode=0; inode<numNodesLine; ++inode )
            for ( int isd=0; isd<3; ++isd )
            {
              A_(isd,2) -= ( *xyze_line_ )(isd,inode) * lineDeriv1_(0,inode);
            }


          dx_ = 0;
          double det = LINALG::gaussElimination<true, 3>( A_, b_, dx_ );

          if ( debug )
            std::cout << "  det=" << std::setprecision( 15 ) << det << "\n";

          if ( fabs( det ) < PARALLEL_DET_TOL )
          {
            if ( iter < 3 )
            {
              // try a different start point
              dx_( ( 0+iter ) % 3 ) = .1;
              dx_( ( 1+iter ) % 3 ) = .12;
              dx_( ( 2+iter ) % 3 ) = .09;
            }
            else
            {
              // There is no solution. The line is parallel to the surface.
              return false;
            }
          }

          // beware of large steps that diverge
          incr_ = dx_.Norm2();
          if ( incr_ > maxincr_ )
          {
            return false;
          }

          return true;
        }

        bool Update( int iter )
        {
          xsi_ += dx_;
          if ( debug )
          {
            std::cout << "  off_count=" << off_count_
                      << "  " << SurfaceWithinLimits<sidetype>( xsi_ )
                      << "  " << LineWithinLimits( xsi_ )
                      << "\n";
            std::cout << "  xsi=" << xsi_;

            LINALG::Matrix<numNodesSurface,1> surface;
            LINALG::Matrix<numNodesLine,1> line;

            DRT::UTILS::shape_function_2D( surface, xsi_( 0 ), xsi_( 1 ), sidetype );
            DRT::UTILS::shape_function_1D( line, xsi_( 2 ), edgetype );

            LINALG::Matrix<3,1> x1;
            LINALG::Matrix<3,1> x2;

            x1.Multiply( *xyze_surface_, surface );
            x2.Multiply( *xyze_line_, line );

            std::cout << "  x1 =" << x1;
            std::cout << "  x2 =" << x2;

            x2.Update( 1, x1, -1 );
            std::cout << "  dx =" << x2;
          }
          if ( not SurfaceWithinLimits<sidetype>( xsi_ ) or not LineWithinLimits( xsi_ ) )
          {
            off_count_ += 1;
            if ( off_count_==3 )
              return false;
          }
          else
          {
            off_count_ = 0;
          }
          return true;
        }

        bool NewtonFailed()
        {
          std::stringstream str;
          str << "Newton scheme did not converge:\n  "
              << ( *xyze_surface_ )
              << ( *xyze_line_ )
            ;
          DumpDoubles( str, xyze_surface_->A(), xyze_surface_->M()*xyze_surface_->N() );
          str << "\n";
          DumpDoubles( str, xyze_line_   ->A(), xyze_line_   ->M()*xyze_line_   ->N() );
          str << "\n";
          throw std::runtime_error( str.str() );
        }

      private:

        LINALG::Matrix<3,1> & xsi_;

        const LINALG::Matrix<3,numNodesSurface> * xyze_surface_;
        const LINALG::Matrix<3,numNodesLine> * xyze_line_;

        LINALG::Matrix<numNodesSurface,1> surfaceFunct_;
        LINALG::Matrix<numNodesLine,1> lineFunct_;
        LINALG::Matrix<2,numNodesSurface> surfaceDeriv1_;
        LINALG::Matrix<1,numNodesLine> lineDeriv1_;

        LINALG::Matrix<3,3> A_;
        LINALG::Matrix<3,1> b_;
        LINALG::Matrix<3,1> dx_;

        LINALG::Matrix<3,1> old_b_;

        double residual_;
        double incr_;
        double maxincr_;
        int off_count_;
      };

      /// intersection of side and line
      template <class Strategy,
                DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class GenericComputeIntersection : Strategy
      {
      public:

        GenericComputeIntersection( LINALG::Matrix<3,1> & xsi )
          : Strategy( xsi )
        {
        }

        bool operator()( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                         const LINALG::Matrix<3,numNodesLine> & xyze_line )
        {
          this->Setup( xyze_surface, xyze_line );
          if ( this->Solve() )
          {
            return true;
          }
          return false;
        }

        const LINALG::Matrix<3,1> & LocalCoordinates()
        {
          return Strategy::LocalCoordinates();
        }

        bool SurfaceWithinLimits()
        {
          return SurfaceWithinLimits<sidetype>( LocalCoordinates() );
        }

        bool LineWithinLimits()
        {
          return LineWithinLimits( LocalCoordinates() );
        }
      };


      template <DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class ComputeIntersection : public GenericComputeIntersection<NewtonSolve<ComputeIntersectionStrategy<false,edgetype,sidetype> >,
                                                                    edgetype,sidetype>
      {
      public:

        ComputeIntersection( LINALG::Matrix<3,1> & xsi )
          : GenericComputeIntersection<NewtonSolve<ComputeIntersectionStrategy<false,edgetype,sidetype> >,
                                       edgetype,sidetype>( xsi )
        {
        }
      };

      template <DRT::Element::DiscretizationType edgetype,
                DRT::Element::DiscretizationType sidetype,
                int numNodesLine    = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
                int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement
      >
      class DebugComputeIntersection : public GenericComputeIntersection<NewtonSolve<DebugNewtonStrategy<ComputeIntersectionStrategy<true,edgetype,sidetype> > >,
                                                                         edgetype,sidetype>
      {
      public:

        DebugComputeIntersection( LINALG::Matrix<3,1> & xsi )
          : GenericComputeIntersection<NewtonSolve<DebugNewtonStrategy<ComputeIntersectionStrategy<true,edgetype,sidetype> > >,
                                       edgetype,sidetype>( xsi )
        {
        }
      };

    }
  }
}

#endif
