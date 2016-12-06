/*---------------------------------------------------------------------*/
/*!
\file cut_kernel.H

\brief the cut kernel computes basic geometric operation, implemeted are
    - Intersection of Surface and line
    - Calculate local coordinates inside an element
    - Compute Distance from a point to a surface

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

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
          return ( fabs( xsi( 0 )+1 )<REFERENCETOL or fabs( xsi( 1 )+1 )<REFERENCETOL or
                   fabs( xsi( 0 )-1 )<REFERENCETOL or fabs( xsi( 1 )-1 )<REFERENCETOL );
        }
        case DRT::Element::tri3:
        case DRT::Element::tri6:
        {
          return ( fabs( xsi( 0 )+0 )<REFERENCETOL or fabs( xsi( 1 )+0 )<REFERENCETOL or
                   fabs( xsi( 1 )+xsi( 0 )-1 )<REFERENCETOL );
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
          intersects = ( xsi( 0 ) >= -1-REFERENCETOL and xsi( 1 ) >= -1-REFERENCETOL and
                         xsi( 0 ) <=  1+REFERENCETOL and xsi( 1 ) <=  1+REFERENCETOL );
          break;
        }
        case DRT::Element::tri3:
        case DRT::Element::tri6:
        {
          intersects = ( xsi( 0 ) >= 0-REFERENCETOL and xsi( 1 ) >= 0-REFERENCETOL and
                         xsi( 0 ) <= 1+REFERENCETOL and xsi( 1 ) <= 1+REFERENCETOL and
                         xsi( 1 ) <= ( 1 - xsi( 0 ) )+REFERENCETOL );
          break;
        }
        default:
          throw std::runtime_error( "unsupported side type" );
        }
        return intersects;
      }


      inline bool LineWithinLimits( const LINALG::Matrix<3,1> & xsi )
      {
        return xsi( 2 ) >= -1-REFERENCETOL and xsi( 2 ) <=  1+REFERENCETOL;
      }


      template <DRT::Element::DiscretizationType elementtype>
      bool WithinLimits( const LINALG::Matrix<3,1> & xsi )
      {
        switch ( elementtype )
        {
        case DRT::Element::hex8:
        case DRT::Element::hex20:
        case DRT::Element::hex27:
          return ( xsi( 0 ) >= -1-REFERENCETOL and xsi( 1 ) >= -1-REFERENCETOL and xsi( 2 ) >= -1-REFERENCETOL and
                   xsi( 0 ) <=  1+REFERENCETOL and xsi( 1 ) <=  1+REFERENCETOL and xsi( 2 ) <=  1+REFERENCETOL );
        case DRT::Element::tet4:
        case DRT::Element::tet10:
          return ( xsi( 0 ) >=   -REFERENCETOL and xsi( 1 ) >=   -REFERENCETOL and xsi( 2 ) >=   -REFERENCETOL and
                   xsi( 0 )+xsi( 1 )+xsi( 2 ) <= 1+REFERENCETOL );
        case DRT::Element::pyramid5:
          return ( xsi( 0 ) >= -1-REFERENCETOL and xsi( 1 ) >= -1-REFERENCETOL and xsi( 2 ) >=   -REFERENCETOL and
                   xsi( 0 ) <=  1+REFERENCETOL and xsi( 1 ) <=  1+REFERENCETOL and
                   ( fabs( xsi( 0 ) ) > fabs( xsi( 1 ) ) ?
                     ( xsi( 2 )+fabs( xsi( 0 ) ) <= 1+REFERENCETOL ) :
                     ( xsi( 2 )+fabs( xsi( 1 ) ) <= 1+REFERENCETOL ) )
            );
        case DRT::Element::wedge6:
        case DRT::Element::wedge15:
          return ( xsi( 0 ) >=   -REFERENCETOL and xsi( 1 ) >=   -REFERENCETOL and xsi( 2 ) >= -1-REFERENCETOL and
                   xsi( 2 ) <=  1+REFERENCETOL and
                   xsi( 0 )+xsi( 1 ) <= 1+REFERENCETOL );
        default:
          throw std::runtime_error( "unsupported element type" );
        }
        return false;
      }

      /*!
      \brief Check whether three points are on the same line
       */
      bool IsOnLine( Point*& pt1, Point*& pt2, Point*& pt3, bool DeleteInlinePts = false );

      /*!
      \brief Check whether the polygon defined by the set of points is convex
       */
      std::vector<int> CheckConvexity( const std::vector<Point*>& ptlist,
                                       GEO::CUT::FacetShape& geomType,
                                       bool InSplit = true,
                                       bool DeleteInlinePts = false );

      //std::vector<double> EqnPlanePolygon( const std::vector<Point*>& ptlist, bool DeleteInlinePts = false );


      /*!
      \brief Compute the equation of plane of this polygon using Newell's method
       */
      std::vector<double> EqnPlaneOfPolygon( const std::vector<Point*>& ptlist );

      /*!
      \brief Compute the equation of plane of this polygon using Newell's method
       */
      std::vector<double> EqnPlaneOfPolygon( const std::vector<std::vector<double> >& vertices );

      /*!
      \brief Find the equation of plane that contains these non-collinear points
      It must be noted while using this function to find equation of facets,
             none of these 3 points must be a reflex (concave) point
       */
      std::vector<double> EqnPlane( Point* & pt1, Point* & pt2, Point* & pt3 );

      /*!
      \brief Check whether the point "check" is inside the triangle formed by tri
      */
      bool PtInsideTriangle( std::vector<Point*> tri, Point* check, bool DeleteInlinePts = false );

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

      /*!
      \brief Returns true if at least 3 points are collinear
      */
      bool HaveInlinePts( std::vector<Point*>& poly );

      /*!
      \brief Finds tree points of the polygon which are not collinear
      */
      std::vector<Point*> Get3NoncollinearPts( std::vector<Point*> & polyPoints );

      /*!
       \brief Find the appropriate coordinate plane to which the given geometry can be
       projected given the equation of plane containing this geometry
       */
      void FindProjectionPlane( std::string& projPlane, const std::vector<double>& eqn );

      /*!
      \brief Get area of triangle in 3D space
       */
      double getAreaTri( std::vector<Point*> & poly );

      /*!
      \brief Get area of convex Quad in 3D space
       */
      double getAreaConvexQuad( std::vector<Point*> & poly );

      /// generic Newton algorithm
      template <class Strategy, int maxiter=20>
      class NewtonSolve : public Strategy
      {
      public:

        NewtonSolve( LINALG::Matrix<3,1> & xsi, bool checklimits )
          : Strategy( xsi, checklimits )
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
      };

      /// Debug helper strategy for generic Newton algorithm
      template <class Strategy>
      class DebugNewtonStrategy : public Strategy
      {
      public:

        DebugNewtonStrategy( LINALG::Matrix<3,1> & xsi, bool checklimits )
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

        ComputePositionStrategy( LINALG::Matrix<3,1> & xsi, bool checklimits )
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

          if ( debug )
          {
            std::cout << " residual " << residual
                      << "\n"
                      << " xsi " << xsi_
                      << " increment " << dx_
                      << " rhs b " << b_;
          }

          return fabs( residual ) < REFERENCETOL;
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

        static LINALG::Matrix<numNodesElement,1> funct_;
        static LINALG::Matrix<3,numNodesElement> deriv1_;

        static LINALG::Matrix<3,3> A_;
        static LINALG::Matrix<3,1> b_;
        static LINALG::Matrix<3,1> dx_;
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
          : Strategy( xsi, true )
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

        ComputeDistanceStrategy( LINALG::Matrix<3,1> & xsi, bool checklimits )
          : xsi_( xsi ),
            xyze_surface_( NULL ),
            px_( NULL ),
            tol_(0.0),
            zeroarea_(false)
        {
        }

        void Setup( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                    const LINALG::Matrix<3,1> & px )
        {
          xyze_surface_ = &xyze_surface;
          px_ = &px;

          //use tol_ first just to store the scaling for the tolerance ..
          //set it to the maximal coordinate of the involved surface or point! -- this defines accuracy of the intersection (for good conditioned linear systems)
          tol_ = xyze_surface.NormInf();
          double linescale = px.NormInf();
          if (linescale > tol_)
            tol_ = linescale;

          tol_ *= LINSOLVETOL;
        }

        double Distance() { return fabs( xsi_( 2 ) ); }

        double SignedDistance() { return  xsi_( 2 ); }

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

          if (zeroarea_) return false;

          double residual = b_.Norm2();
          if (residual < tol_)
          {
            tol_ *= sqrt(3);
            return true;
          }
          else
            return false;
        }

        bool LinearSolve( int iter )
        {
          if (zeroarea_) return false;
          dx_ = 0;
          //std::cout << "GAUSSELEMIN FOR COMPUTE DISTANCE ... " << std::flush;
          double det = LINALG::gaussElimination<true, 3>( A_, b_, dx_ );
          return det != 0; //here det = 0 will just happen if the surface element is distorted!
        }

        bool Update( int iter )
        {
          xsi_ += dx_;
          return true;
        }

        bool NewtonFailed()
        {
          tol_ = b_.Norm2();
          #ifdef DEBUG_CUTKERNEL_OUTPUT
           {
           std::stringstream str;
           str << "ComputeDistance: Newton scheme did not converge:\n"
              << std::setprecision( 16 )
              << ( *xyze_surface_ )
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

        void WritetoGmsh(std::ofstream& file)
        {
          file.precision(32); //higher precicion!
          char elementtype;
          switch (numNodesSurface )
          {
          case 3:
            elementtype = 'T';
            break;
          case 4:
            elementtype = 'Q';
            break;
          default:
            throw std::runtime_error( "unknown element type in GmshSideDump" );
          }

          file << "View \"" << "Side" << "\" {\n";
          {
            file << "S" << elementtype
                 << "(";
            for ( uint i = 0; i < numNodesSurface; ++i )
            {
              if ( i != 0 )
                file << ",";
              file << ( *xyze_surface_ )(0,i) << "," << ( *xyze_surface_ )(1,i) << "," << ( *xyze_surface_ )(2,i);
            }
            file << "){";
            for ( uint i = 0; i < numNodesSurface; ++i )
            {
              if ( i != 0 )
                file << ",";
              file << "1";
            }
            file << "};\n";
          }
          file<<"};\n";

          //calculate the points
          LINALG::Matrix<numNodesSurface,1> surface;

          DRT::UTILS::shape_function_2D( surface, xsi_( 0 ), xsi_( 1 ), sidetype );

          LINALG::Matrix<3,1> x1;;

          x1.Multiply( *xyze_surface_, surface );

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

        double GetTolerance()
        {
          return tol_;
        }

        bool ZeroArea()
        {
          return zeroarea_;
        }

      private:

        bool DistanceSystem( const LINALG::Matrix<3,numNodesSurface> & surfaceDeriv2,
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

          if(fabs(cnorm) == 0.0)
          {
            zeroarea_ = true; //then calculation of a normal to a line like this make no sense at all!
            return false;
          }

          fact1 = 1./cnorm;
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

          return true;
        }

        LINALG::Matrix<3,1> & xsi_;

        const LINALG::Matrix<3,numNodesSurface> * xyze_surface_;
        const LINALG::Matrix<3,1> * px_;

        static LINALG::Matrix<numNodesSurface,1> surfaceFunct_;
        static LINALG::Matrix<2,numNodesSurface> surfaceDeriv1_;
        static LINALG::Matrix<3,numNodesSurface> surfaceDeriv2_;

        static LINALG::Matrix<3,3> A_;
        static LINALG::Matrix<3,3> B_;
        static LINALG::Matrix<3,3> C_;
        static LINALG::Matrix<3,1> b_;
        static LINALG::Matrix<3,1> dx_;
        static LINALG::Matrix<3,1> cross_;

        double tol_; //tolerance!!!

        bool zeroarea_; //surface element has zero area --> cross product for normal is zero!!!
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

        GenericComputeDistance( LINALG::Matrix<3,1> & xsi, bool checklimits )
          : Strategy( xsi, checklimits )
        {
        }

        bool operator()( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                         const LINALG::Matrix<3,1> & px,
                         double & distance,
                         bool signeddistance = false)
        {
          this->Setup( xyze_surface, px );

          bool conv = this->Solve();

          if (!signeddistance)
            distance = this->Distance();
          else
            distance = this->SignedDistance();

          return conv;
        }

        const LINALG::Matrix<3,1> & LocalCoordinates()
        {
          return Strategy::LocalCoordinates();
        }

        bool SurfaceWithinLimits()
        {
          return GEO::CUT::KERNEL::SurfaceWithinLimits<sidetype>( LocalCoordinates() );
        }

        double GetTolerance()
        {
          return Strategy::GetTolerance();
        }

        bool ZeroArea()
        {
          return Strategy::ZeroArea();
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

        ComputeDistance( LINALG::Matrix<3,1> & xsi, bool checklimits =true )
          : GenericComputeDistance<NewtonSolve<ComputeDistanceStrategy<false,edgetype,sidetype> >,
                                   edgetype,sidetype>( xsi, checklimits )
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

        DebugComputeDistance( LINALG::Matrix<3,1> & xsi, bool checklimits = true )
          : GenericComputeDistance<NewtonSolve<DebugNewtonStrategy<ComputeDistanceStrategy<true,edgetype,sidetype> > >,
                                   edgetype,sidetype>( xsi, checklimits )
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

        ComputeIntersectionStrategy( LINALG::Matrix<3,1> & xsi, bool checklimits )
          :  xsi_( xsi ),
             xyze_surface_( NULL ),
             xyze_line_( NULL ),
             check_limits_( checklimits ),
             tol_(0.0)
        {
        }

        void Setup( const LINALG::Matrix<3,numNodesSurface> & xyze_surface,
                    const LINALG::Matrix<3,numNodesLine> & xyze_line )
        {
          xyze_surface_ = &xyze_surface;
          xyze_line_    = &xyze_line;

          //use tol_ first just to store the scaling for the tolerance ..
          //set it to the maximal coordinate of the involved surface or line! -- this defines accuracy of the intersection (for good conditioned linear systems)
          tol_ = xyze_surface.NormInf();
          double linescale = xyze_line.NormInf();
          if (linescale > tol_)
            tol_ = linescale;

          tol_ *= LINSOLVETOL;

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
          off_count_ = 0;          // count how many iterations during the newton, the local coordinates are not within limits!
                                   // stop the newton scheme if the iteration is too often outside the limits

          residual_ = 0;           // residual

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
          residual_ = b_.Norm2();

          if ( debug )
          {
            double incr = dx_.Norm2();
            std::cout << "  residual=" << std::setprecision( 15 ) << residual_
                      << "  incr=" << std::setprecision( 15 ) << incr
                      << "\n";
          }

          if ( residual_ < tol_ && iter > 1)
          {
            tol_ *= sqrt(3);
            return true;
          }
          return false;
        }

        bool LinearSolve( int iter )
        {
          DRT::UTILS::shape_function_2D_deriv1( surfaceDeriv1_, xsi_( 0 ), xsi_( 1 ), sidetype );
          DRT::UTILS::shape_function_1D_deriv1( lineDeriv1_, xsi_( 2 ), edgetype );

          A_ = 0.0;  //Fill with double, same as PutScalar(other).

          //I guess some improvement is possible here... (matrix mult.!!!)
          for ( int inode=0; inode<numNodesSurface; ++inode )
          {
            for( int isd=0; isd<3; ++isd )
            {
              A_(isd,0) += ( *xyze_surface_ )(isd,inode) * surfaceDeriv1_(0,inode);
              A_(isd,1) += ( *xyze_surface_ )(isd,inode) * surfaceDeriv1_(1,inode);
            }
          }

          for ( int inode=0; inode<numNodesLine; ++inode )
          {
            for ( int isd=0; isd<3; ++isd )
            {
              A_(isd,2) -= ( *xyze_line_ )(isd,inode) * lineDeriv1_(0,inode);
            }
          }


          dx_ = 0;

          double det = LINALG::gaussElimination<true, 3>( A_, b_, dx_ );

          if ( debug )
            std::cout << "  det=" << std::setprecision( 15 ) << det << "\n";

          return true;
        }

        bool Update( int iter )
        {
          xsi_ += dx_;

//          std::cout << "xsi: " << xsi_ << "   dx: " << dx_ << std::endl;

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

          return true;
        }

        bool NewtonFailed()
        {
          tol_ = b_.Norm2()*sqrt(3);
          //if (debug)
          #ifdef DEBUG_CUTKERNEL_OUTPUT
            {
            std::stringstream str;
            str << "Newton scheme did not converge:\n  "
                << ( *xyze_surface_ )
                << ( *xyze_line_ )
              ;
            DumpDoubles( str, xyze_surface_->A(), xyze_surface_->M()*xyze_surface_->N() );
            str << "\n";
            DumpDoubles( str, xyze_line_   ->A(), xyze_line_   ->M()*xyze_line_   ->N() );


            std::stringstream str;
            str  <<".NewtonFailed_intersection.pos";
            std::string filename(GEO::CUT::OUTPUT::GenerateGmshOutputFilename(str.str()));
            std::ofstream file(filename.c_str());
            WritetoGmsh(file);
            }
          #endif
          return false;
        }

        double GetTolerance()
        {
          return tol_;
        }

        void WritetoGmsh(std::ofstream& file)
        {
          file.precision(32); //higher precicion!
          char elementtype;
          switch (numNodesSurface )
          {
          case 3:
            elementtype = 'T';
            break;
          case 4:
            elementtype = 'Q';
            break;
          default:
            throw std::runtime_error( "unknown element type in GmshSideDump" );
          }

          file << "View \"" << "Side" << "\" {\n";
          {
            file << "S" << elementtype
                 << "(";
            for ( uint i = 0; i < numNodesSurface; ++i )
            {
              if ( i != 0 )
                file << ",";
              file << ( *xyze_surface_ )(0,i) << "," << ( *xyze_surface_ )(1,i) << "," << ( *xyze_surface_ )(2,i);
            }
            file << "){";
            for ( uint i = 0; i < numNodesSurface; ++i )
            {
              if ( i != 0 )
                file << ",";
              file << "1";
            }
            file << "};\n";
          }
          file<<"};\n";

          file << "View \"" << "Line" << "\" {\n";
          for (int iline = 0; iline < numNodesLine -1; ++iline)
          {
            file << "SL (";
            file << ( *xyze_line_ )(0,iline) << "," << ( *xyze_line_ )(1,iline)
                << "," << ( *xyze_line_ )(2,iline) << ",";
            file << ( *xyze_line_ )(0,iline+1) << "," << ( *xyze_line_ )(1,iline+1)
                << "," << ( *xyze_line_ )(2,iline+1);
            file << "){";
            file << 1<< ",";
            file << 1;
            file << "};\n";
          }
          file<<"};\n";


          //calculate the points
          LINALG::Matrix<numNodesSurface,1> surface;
          LINALG::Matrix<numNodesLine,1> line;

          DRT::UTILS::shape_function_2D( surface, xsi_( 0 ), xsi_( 1 ), sidetype );
          DRT::UTILS::shape_function_1D( line, xsi_( 2 ), edgetype );

          LINALG::Matrix<3,1> x1;
          LINALG::Matrix<3,1> x2;

          x1.Multiply( *xyze_surface_, surface );
          x2.Multiply( *xyze_line_, line );
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

        LINALG::Matrix<3,1> & xsi_;

        const LINALG::Matrix<3,numNodesSurface> * xyze_surface_;
        const LINALG::Matrix<3,numNodesLine> * xyze_line_;

         static LINALG::Matrix<numNodesSurface,1> surfaceFunct_;
         static LINALG::Matrix<numNodesLine,1> lineFunct_;
         static LINALG::Matrix<2,numNodesSurface> surfaceDeriv1_;
         static LINALG::Matrix<1,numNodesLine> lineDeriv1_;

         static LINALG::Matrix<3,3> A_;
         static LINALG::Matrix<3,1> b_;
         static LINALG::Matrix<3,1> dx_;

        double residual_;
        int off_count_;
        bool check_limits_;
        double tol_; //basic tolerance!!!
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

        GenericComputeIntersection( LINALG::Matrix<3,1> & xsi, bool checklimits)
          : Strategy( xsi, checklimits)
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
          return GEO::CUT::KERNEL::SurfaceWithinLimits<sidetype>( LocalCoordinates() );
        }

        bool LineWithinLimits()
        {
          return LineWithinLimits( LocalCoordinates() );
        }

        double GetTolerance()
        {
          return Strategy::GetTolerance();
        }

        void WritetoGmsh(std::ofstream& file)
        {
          return Strategy::WritetoGmsh(file);
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

        ComputeIntersection( LINALG::Matrix<3,1> & xsi, bool checklimits = true )
          : GenericComputeIntersection<NewtonSolve<ComputeIntersectionStrategy<false,edgetype,sidetype> >,
                                       edgetype,sidetype>( xsi, checklimits )
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

        DebugComputeIntersection( LINALG::Matrix<3,1> & xsi, bool checklimits = true )
          : GenericComputeIntersection<NewtonSolve<DebugNewtonStrategy<ComputeIntersectionStrategy<true,edgetype,sidetype> > >,
                                       edgetype,sidetype>( xsi, checklimits )
        {
        }

      };

    }
  }
}


//static members in KERNEL
//in compute position strategy
template <bool debug,DRT::Element::DiscretizationType elementtype, int numNodesElement > LINALG::Matrix<numNodesElement , 1> GEO::CUT::KERNEL::ComputePositionStrategy<debug, elementtype, numNodesElement >::funct_;
template <bool debug,DRT::Element::DiscretizationType elementtype, int numNodesElement > LINALG::Matrix<3 , numNodesElement> GEO::CUT::KERNEL::ComputePositionStrategy<debug, elementtype, numNodesElement >::deriv1_;
template <bool debug,DRT::Element::DiscretizationType elementtype, int numNodesElement > LINALG::Matrix<3 , 3> GEO::CUT::KERNEL::ComputePositionStrategy<debug, elementtype, numNodesElement >::A_;
template <bool debug,DRT::Element::DiscretizationType elementtype, int numNodesElement > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputePositionStrategy<debug, elementtype, numNodesElement >::b_;
template <bool debug,DRT::Element::DiscretizationType elementtype, int numNodesElement > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputePositionStrategy<debug, elementtype, numNodesElement >::dx_;

//in compute distance strategy
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<numNodesSurface , 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::surfaceFunct_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<2 , numNodesSurface> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::surfaceDeriv1_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , numNodesSurface> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::surfaceDeriv2_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , 3> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::A_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , 3> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::B_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , 3> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::C_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::b_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::dx_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine, int numNodesSurface > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputeDistanceStrategy<debug, edgetype, sidetype, numNodesLine, numNodesSurface >::cross_;

  //in compute intersection strategy
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<numNodesSurface , 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::surfaceFunct_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<numNodesLine , 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::lineFunct_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<2 , numNodesSurface> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::surfaceDeriv1_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<1 , numNodesLine> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::lineDeriv1_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<3 , 3> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::A_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::b_;
template <bool debug, DRT::Element::DiscretizationType edgetype, DRT::Element::DiscretizationType sidetype, int numNodesLine , int numNodesSurface > LINALG::Matrix<3 , 1> GEO::CUT::KERNEL::ComputeIntersectionStrategy<debug, edgetype, sidetype,numNodesLine, numNodesSurface>::dx_;

#endif
