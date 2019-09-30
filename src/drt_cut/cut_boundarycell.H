/*---------------------------------------------------------------------*/
/*! \file

\brief cut boundary cell header file

\level 3

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/
#ifndef CUT_BOUNDARYCELL_H
#define CUT_BOUNDARYCELL_H

#include "cut_tolerance.H"
#include "cut_kernel.H"

#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include <Teuchos_TimeMonitor.hpp>


namespace GEO
{
  namespace CUT
  {
    class Facet;
    class Element;
    class VolumeCell;
    class Mesh;
    class Point;
    class Cycle;

    /*----------------------------------------------------------------------------*/
    /*! \brief Base class for boundary cells. Boundary cells are used to represent
     *  the cut surface. Each volume cell has its own boundary cells
     *  at any cut surface with outward normals. */
    class BoundaryCell
    {
     public:
      /// constructor
      BoundaryCell(
          const Epetra_SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points);

      /// destructor
      virtual ~BoundaryCell() {}

      /*!
      \brief Returns the shape of the boundarycell
       */

      virtual DRT::Element::DiscretizationType Shape() const = 0;

      /*! \brief Returns the cubature degree to generate quadrature rule for the cell
       *
       *  This is the maximal polynomial degree which is integrate exactly by the used
       *  gaussion quadrature rule. */
      virtual int CubatureDegree() const = 0;

      /*!
      \brief Writes the geometry of boundarycell, and the constant scalar "value" in GMSH format
       */
      void DumpGmsh(std::ofstream& file, int* value = NULL);

      /*!
      \brief Writes the geometry of boundarycell, and the constant scalar "value" in GMSH format
       */
      virtual void DumpGmshNormal(std::ofstream& file) = 0;

      /*!
      \brief Returns the area of tri3 and quad4 boundarycell
       */
      virtual double Area() = 0;

      /*!
      \brief Returns the center of tri3 and quad4 boundarycell
       */
      virtual void ElementCenter(LINALG::Matrix<3, 1>& midpoint) = 0;

      /*!
      \brief Get the outward normal vector for the tri3 and quad4 boundarycell
       */
      virtual void Normal(const LINALG::Matrix<2, 1>& xsi, LINALG::Matrix<3, 1>& normal) const = 0;

      /*!
      \brief Get the corner points of the boundarycell in Cylce for geometrical operations
       */
      const Cycle& PointCycle() const { return *points_; }

      /*!
      \brief Get the corner points of the boundarycell as vector of Points
       */
      const std::vector<Point*>& Points() const;

      /*!
      \brief Get the coodinates of the corner points of boundarycell
       */
      const Epetra_SerialDenseMatrix& Coordinates() const { return xyz_; }

      /*!
      \brief Move the corner points of the boundary cell by a specific offset
       */
      template <DRT::Element::DiscretizationType celldistype>
      void AssignOffset(int idx, double offset)
      {
        LINALG::Matrix<3, DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement>
            xyz_mat(xyz_, true);
        for (uint n = 0; n < DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement;
             ++n)
          xyz_mat(idx, n) += offset;
      }

      /*! \brief Get the coodinates of the corner points of boundarycell
       *
       *  To make it easier with DD data structures (should be removed at some point)
       *
       *  \author winter
       *  \date 07/15 */
      std::vector<std::vector<double>> CoordinatesV();

      /*!
      \brief Get the Facet which represent the boundarycell
       */
      Facet* GetFacet() { return facet_; }
      const Facet* GetFacet() const { return facet_; }

      /*!
      \brief Delete all the points of this boundarycell
       */
      void Clear();

      bool IsValid() const;

      /*!
      \brief Function to test if the distance between points is within point tolerance.
       */
      virtual bool IsValidBoundaryCell() = 0;

      /*!
      \brief Get the Gaussian integration rule for the boundarycell
       */
      virtual DRT::UTILS::GaussIntegration gaussRule(int cubaturedegree) = 0;

      /*!
      \brief Get the Gaussian integration rule for the boundarycell
       */
      virtual DRT::UTILS::GaussIntegration gaussRule() { return gaussRule(CubatureDegree()); }

      /*!
      \brief Get the normal vector for the arbitrary boundarycell alone
       */
      virtual LINALG::Matrix<3, 1> GetNormalVector() = 0;

      /*! \brief Print the corner points on screen
       *
       *  just for debugging
       *
       *  \author sudhakar
       *  \date 10/14
       */
      void Print(std::ostream& stream);
      void Print() { Print(std::cout); }

      /*!
      \brief Computes the location of Gauss points on the boundarycell (x_gp_lin) from the standard
      Gauss point location (eta) corresponding to the shape of the boundarycell. Also computes the
      factor to be multiplied with integration weight and normal vector of the boundarycell
       */
      template <DRT::Element::DiscretizationType celldistype>
      void Transform(const LINALG::Matrix<2, 1>& eta, LINALG::Matrix<3, 1>& x_gp_lin,
          LINALG::Matrix<3, 1>& normal, double& drs, bool referencepos = false)
      {
        const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement;
        LINALG::Matrix<3, numnodes> xyze(xyz_, true);
        if (referencepos) xyze = LINALG::Matrix<3, numnodes>(xyz_ref_, true);

        LINALG::Matrix<numnodes, 1> funct(false);
        LINALG::Matrix<2, numnodes> deriv(false);
        LINALG::Matrix<2, 2> metrictensor(false);

        DRT::UTILS::shape_function_2D(funct, eta(0), eta(1), celldistype);

        if (celldistype != DRT::Element::tri3)
        {
          DRT::UTILS::shape_function_2D_deriv1(deriv, eta(0), eta(1), celldistype);
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<celldistype>(
              xyze, deriv, metrictensor, drs, &normal);
        }
        else
        {
          // For Tri's this method of determining the area and thus the gp-weights is more robust.
          //  It is needed for TRI's which are small/ill-conditioned but large enough to affect the
          //  simulation.
          static LINALG::Matrix<3, 1> p0(true);
          static LINALG::Matrix<3, 1> p1(true);
          static LINALG::Matrix<3, 1> p2(true);
          for (uint dim = 0; dim < 3; ++dim)
          {
            p0(dim) = xyze(dim, 0);
            p1(dim) = xyze(dim, 1);
            p2(dim) = xyze(dim, 2);
          }
          drs = 2.0 * (CUT::KERNEL::getAreaTri(p0.A(), p1.A(), p2.A(), &normal));
        }

        x_gp_lin.Multiply(xyze, funct);

        return;
      }

      /// Reset the point with local index lid
      void ResetPos(int lid, LINALG::Matrix<3, 1> newpos)
      {
        if (lid > xyz_.N()) dserror("Index out of range! %d > %d", lid, xyz_.N());

        xyz_(0, lid) = newpos(0, 0);
        xyz_(1, lid) = newpos(1, 0);
        xyz_(2, lid) = newpos(2, 0);
      }

      /*!
      \brief Compute the location of Gauss points on the background element's local coordinate
      system setting shadow = true means the mapping is done w.r. to the parent quad element from
      which elem1 is derived
       */
      template <DRT::Element::DiscretizationType celldistype>
      void TransformLocalCoords(Element* elem1, const LINALG::Matrix<2, 1>& eta,
          LINALG::Matrix<3, 1>& x_gp_lin, LINALG::Matrix<3, 1>& normal, double& drs,
          bool shadow = false);

     protected:
      virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() = 0;

      template <DRT::Element::DiscretizationType distype>
      double MyArea()
      {
        const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        LINALG::Matrix<3, numnodes> xyze(this->xyz_, true);
        LINALG::Matrix<numnodes, 1> funct;
        LINALG::Matrix<2, numnodes> deriv;
        LINALG::Matrix<2, 2> metrictensor;

        DRT::UTILS::GaussRule2D gaussrule = this->MySimpleGaussRule();
        DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

        double area = 0;
        double drs;
        for (int i = 0; i < intpoints.nquad; ++i)
        {
          double* eta = intpoints.qxg[i];
          DRT::UTILS::shape_function_2D(funct, eta[0], eta[1], distype);
          DRT::UTILS::shape_function_2D_deriv1(deriv, eta[0], eta[1], distype);
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(
              xyze, deriv, metrictensor, drs, NULL);
          if (not std::isnan(drs)) area += intpoints.qwgt[i] * drs;
        }
        return area;
      }

      /** \brief To calculate the element center of a boundary cell
       *
       *  \author shahmiri
       *  \date 07/12 */
      template <DRT::Element::DiscretizationType distype>
      void MyElementCenter(LINALG::Matrix<3, 1>& center, LINALG::Matrix<3, 1>& midpoint)
      {
        const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        LINALG::Matrix<3, numnodes> xyze(this->xyz_, true);
        LINALG::Matrix<numnodes, 1> funct;
        DRT::UTILS::shape_function<distype>(center, funct);
        midpoint.Multiply(xyze, funct);
      }

      /// Current position of the boundary cell
      Epetra_SerialDenseMatrix xyz_;

      /// Reference position of the boundary cell
      Epetra_SerialDenseMatrix xyz_ref_;
      Facet* facet_;
      Teuchos::RCP<Cycle> points_;
    };

    /*----------------------------------------------------------------------------*/
    /// point1 boundary cell
    class Point1BoundaryCell : public BoundaryCell
    {
     public:
      Point1BoundaryCell(
          const Epetra_SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
        /* intentionally left blank */
      }

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::point1; }

      virtual int CubatureDegree() const { return 0; }

      virtual void DumpGmshNormal(std::ofstream& file);

      virtual double Area() { return 0.0; };

      virtual void ElementCenter(LINALG::Matrix<3, 1>& midpoint);

      virtual void Normal(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<3, 1>& normal) const;

      virtual DRT::UTILS::GaussIntegration gaussRule(int cubaturedegree);

      virtual LINALG::Matrix<3, 1> GetNormalVector();

      virtual bool IsValidBoundaryCell() { return true; };

     protected:
      virtual DRT::UTILS::GaussRule2D MySimpleGaussRule()
      {
        return DRT::UTILS::intrule2D_undefined;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// point1 boundary cell
    class Line2BoundaryCell : public BoundaryCell
    {
     public:
      Line2BoundaryCell(
          const Epetra_SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
        /* intentionally left blank */
      }

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::line2; }

      virtual int CubatureDegree() const { return 4; }

      virtual void DumpGmshNormal(std::ofstream& file);

      virtual double Area();

      virtual void ElementCenter(LINALG::Matrix<3, 1>& midpoint);

      virtual void Normal(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<3, 1>& normal) const;

      virtual DRT::UTILS::GaussIntegration gaussRule(int cubaturedegree);

      virtual LINALG::Matrix<3, 1> GetNormalVector();

      virtual bool IsValidBoundaryCell() { return (Area() > REF_AREA_BCELL); };

     protected:
      virtual DRT::UTILS::GaussRule2D MySimpleGaussRule()
      {
        return DRT::UTILS::intrule2D_undefined;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// tri3 boundary cell
    class Tri3BoundaryCell : public BoundaryCell
    {
     public:
      Tri3BoundaryCell(
          const Epetra_SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
      }

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tri3; }

      virtual int CubatureDegree() const { return 20; }

      virtual void DumpGmshNormal(std::ofstream& file);

      virtual double Area();

      virtual void ElementCenter(LINALG::Matrix<3, 1>& midpoint);

      virtual void Normal(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<3, 1>& normal) const;

      virtual DRT::UTILS::GaussIntegration gaussRule(int cubaturedegree);

      virtual LINALG::Matrix<3, 1> GetNormalVector();

      /** \brief  A first step to validate if a boundary cell is valid.
       *
       *  \author winter
       *  \date 11/15 */
      virtual bool IsValidBoundaryCell();

     protected:
      virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() { return DRT::UTILS::intrule_tri_3point; }
    };

    /*----------------------------------------------------------------------------*/
    /// quad4 boundary cell
    class Quad4BoundaryCell : public BoundaryCell
    {
     public:
      Quad4BoundaryCell(
          const Epetra_SerialDenseMatrix& xyz, Facet* facet, const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
      }

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::quad4; }

      virtual int CubatureDegree() const { return 20; }

      virtual void DumpGmshNormal(std::ofstream& file);

      // Maybe shoelace theorem can be used here?
      virtual double Area() { return MyArea<DRT::Element::quad4>(); }

      virtual void ElementCenter(LINALG::Matrix<3, 1>& midpoint);

      virtual void Normal(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<3, 1>& normal) const;

      virtual DRT::UTILS::GaussIntegration gaussRule(int cubaturedegree);

      virtual LINALG::Matrix<3, 1> GetNormalVector();

      // Probably not the best way...
      virtual bool IsValidBoundaryCell() { return (Area() > REF_AREA_BCELL); }

     protected:
      virtual DRT::UTILS::GaussRule2D MySimpleGaussRule()
      {
        return DRT::UTILS::intrule_quad_4point;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// Irregular boundary cell generated during the cut
    class ArbitraryBoundaryCell : public BoundaryCell
    {
     public:
      ArbitraryBoundaryCell(const Epetra_SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points, const DRT::UTILS::GaussIntegration& gaussRule,
          const LINALG::Matrix<3, 1>& normal)
          : BoundaryCell(xyz, facet, points), gaussRule_(gaussRule), normal_(normal)
      {
      }

      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::dis_none; }

      virtual int CubatureDegree() const { return 0; }

      virtual void DumpGmshNormal(std::ofstream& file);

      virtual double Area() { return 0.0; }

      virtual void ElementCenter(LINALG::Matrix<3, 1>& midpoint);

      virtual void Normal(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<3, 1>& normal) const;

      virtual DRT::UTILS::GaussIntegration gaussRule(int cubaturedegree);

      virtual LINALG::Matrix<3, 1> GetNormalVector();

      virtual bool IsValidBoundaryCell() { return (Area() > REF_AREA_BCELL); }

     protected:
      virtual DRT::UTILS::GaussRule2D MySimpleGaussRule()
      {
        return DRT::UTILS::intrule_quad_4point;
      }

     private:
      DRT::UTILS::GaussIntegration gaussRule_;
      LINALG::Matrix<3, 1> normal_;
    };

  }  // namespace CUT
}  // namespace GEO

#endif
