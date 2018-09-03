/*---------------------------------------------------------------------*/
/*!
\file volume_integration.H

\brief Integrates base functions over volume, distribute Gauss points and solve moment fitting
equations

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/


#ifndef VOLUME_INTEGRATION_H
#define VOLUME_INTEGRATION_H

#include <vector>
#include "facet_integration.H"
#include "cut_volumecell.H"
#include "cut_enum.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Facet;
    class VolumeCell;
    // While performing the volume integration, the points of the facet should be arranged in
    // anti-clockwise manner when looking the surface away from the body this ensures outward normal
    // vector
    /*!
    \brief Class to construct integration rule over the volumecell in element local coordinates.
    It performs integration of base functions over the volume, distribution of Gauss points and
    solution of moment fitting matrix to arrive at the Gauss weights
     */
    class VolumeIntegration
    {
     public:
      VolumeIntegration(VolumeCell* volcell, Element* elem,
          const GEO::CUT::Point::PointPosition posi, int num_func)
          : volcell_(volcell), elem1_(elem), position_(posi), num_func_(num_func)
      {
      }

      ~VolumeIntegration(){};

      /*!
      \brief Compute Gauss point weights by solving the moment fitting equations and returns the
      coordinates of Gauss points and their corresponding weights
      */
      Epetra_SerialDenseVector compute_weights();

      /*!
      \brief Computes the RHS of the moment fitting matrix (performs integration of base functions
      over the volumecell)
      */
      Epetra_SerialDenseVector compute_rhs_moment();

      /*!
      \brief Returns the location of Gauss points distributed over the volumecell
      */
      std::vector<std::vector<double>> getGaussPointLocation() { return gaus_pts_; }

      /*!
      \brief Check whether the point with this element Local coordinates is inside, outisde or
      boundary of this volumecell. The output std::string will be either "outside", "inside" or
      "onBoundary"
      */
      std::string IsPointInside(LINALG::Matrix<3, 1>& rst);

     private:
      /*!
      \brief Distribute Gaussian points over the volumecell with "numeach" points in each direction
      */
      bool compute_Gaussian_points(int numeach);

      /*!
      \brief Computes the moment fitting matrix
      */
      void moment_fitting_matrix(
          std::vector<std::vector<double>>& mom, std::vector<std::vector<double>> gauspts);

      /*!
      \brief Check whether the generated ray intersect the facets of the volumecell, if so
      distribute Gauss points along this ray
      */
      bool IsIntersect(double* pt, double* mini, double* maxi,
          std::vector<std::vector<double>>& linePts, std::vector<std::vector<double>> zcoord,
          std::vector<std::vector<double>> ycoord, double toler, int numeach);

      /*!
      \brief Check whether the particular z-plane of the volumecell contains significant area so as
      to distribute the Gauss points in that plane
      */
      bool IsContainArea(double minn[3], double maxx[3], double& zmin,
          std::vector<std::vector<double>>& pts, std::vector<std::vector<double>> zcoord,
          std::vector<std::vector<double>> ycoord, double toler, int numeach);

      /*!
      \brief Writes the Geometry of volumecell and location of Gauss points in GMSH format output
      file
      */
      void GaussPointGmsh();

      /*!
      \brief Generates equally spaced "num" number of points on the line whose end points are
      specified by inter1 and inter2
      */
      void OnLine(std::vector<double> inter1, std::vector<double> inter2,
          std::vector<std::vector<double>>& linePts, int num);

      /*!
      \brief Store the z- and y-coordinates of the all corner points which will be used to find
      whether the intersection point lies inside the volume or not
      */
      void get_zcoordinates(
          std::vector<std::vector<double>>& zcoord, std::vector<std::vector<double>>& ycoord);

      /*
      \brief Check whether the intersection point, which is in the plane containing the facet,
      actually lies with in the facet area
      */
      int pnpoly(int npol, std::vector<double> xp, std::vector<double> yp, double x, double y);

      /*
      \brief Check whether the intersection point, which is in the plane containing the facet,
      actually lies with in the facet area
      */
      int pnpoly(const std::vector<std::vector<double>>& xp, const LINALG::Matrix<3, 1>& pt,
          GEO::CUT::ProjectionDirection projType);

      /*!
      \brief Adds linear combination of first order base function to the moment fitting system of
      equations
      */
      void FirstOrderAdditionalTerms(
          std::vector<std::vector<double>>& mat, Epetra_SerialDenseVector& rhs);

      /*!
      \brief Adds linear combination of second order base function to the moment fitting system of
      equations
      */
      void SecondOrderAdditionalTerms(
          std::vector<std::vector<double>>& mat, Epetra_SerialDenseVector& rhs);

      /*!
      \brief Adds linear combination of third order base function to the moment fitting system of
      equations
      */
      void ThirdOrderAdditionalTerms(
          std::vector<std::vector<double>>& mat, Epetra_SerialDenseVector& rhs);

      /*!
      \brief Adds linear combination of fourth order base function to the moment fitting system of
      equations
      */
      void FourthOrderAdditionalTerms(
          std::vector<std::vector<double>>& mat, Epetra_SerialDenseVector& rhs);

      /*!
      \brief Adds linear combination of fifth order base function to the moment fitting system of
      equations
      */
      void FifthOrderAdditionalTerms(
          std::vector<std::vector<double>>& mat, Epetra_SerialDenseVector& rhs);

      /*!
      \brief Adds linear combination of sixth order base function to the moment fitting system of
      equations
      */
      void SixthOrderAdditionalTerms(
          std::vector<std::vector<double>>& mat, Epetra_SerialDenseVector& rhs);

      /*!
      \brief Computes the error introduced by the generated integration rule for integrating some
      specific functions
      */
      void ErrorForSpecificFunction(
          Epetra_SerialDenseVector rhs_moment, Epetra_SerialDenseVector weights, int numeach);

      //! considered volumecell
      VolumeCell* volcell_;

      //! background element that contains this volumecell
      Element* elem1_;

      //! position (inside or outside) of volumecell
      const GEO::CUT::Point::PointPosition position_;

      //! defines the base function to be integrated
      int num_func_;

      //! position of Gauss points
      std::vector<std::vector<double>> gaus_pts_;

      //! equations of plane in which facets of volumecell are contained
      std::vector<std::vector<double>> eqn_facets_;
    };
  }  // namespace CUT
}  // namespace GEO

#endif
