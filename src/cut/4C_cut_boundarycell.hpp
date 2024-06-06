/*---------------------------------------------------------------------*/
/*! \file

\brief cut boundary cell header file

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_CUT_BOUNDARYCELL_HPP
#define FOUR_C_CUT_BOUNDARYCELL_HPP

#include "4C_config.hpp"

#include "4C_cut_kernel.hpp"
#include "4C_cut_tolerance.hpp"
#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
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
      BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points);

      /// destructor
      virtual ~BoundaryCell() = default;
      /*!
      \brief Returns the shape of the boundarycell
       */

      virtual Core::FE::CellType Shape() const = 0;

      /*! \brief Returns the cubature degree to generate quadrature rule for the cell
       *
       *  This is the maximal polynomial degree which is integrate exactly by the used
       *  gaussion quadrature rule. */
      virtual int CubatureDegree() const = 0;

      /*!
      \brief Writes the geometry of boundarycell, and the constant scalar "value" in GMSH format
       */
      void DumpGmsh(std::ofstream& file, int* value = nullptr);

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
      virtual void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) = 0;

      /*!
      \brief Get the outward normal vector for the tri3 and quad4 boundarycell
       */
      virtual void Normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const = 0;

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
      const Core::LinAlg::SerialDenseMatrix& Coordinates() const { return xyz_; }

      /*!
      \brief Move the corner points of the boundary cell by a specific offset
       */
      template <Core::FE::CellType celldistype>
      void AssignOffset(int idx, double offset)
      {
        Core::LinAlg::Matrix<3, Core::FE::num_nodes<celldistype>> xyz_mat(xyz_, true);
        for (unsigned int n = 0; n < Core::FE::num_nodes<celldistype>; ++n)
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
      virtual Core::FE::GaussIntegration gaussRule(int cubaturedegree) = 0;

      /*!
      \brief Get the Gaussian integration rule for the boundarycell
       */
      virtual Core::FE::GaussIntegration gaussRule() { return gaussRule(CubatureDegree()); }

      /*!
      \brief Get the normal vector for the arbitrary boundarycell alone
       */
      virtual Core::LinAlg::Matrix<3, 1> GetNormalVector() = 0;

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
      template <Core::FE::CellType celldistype>
      void Transform(const Core::LinAlg::Matrix<2, 1>& eta, Core::LinAlg::Matrix<3, 1>& x_gp_lin,
          Core::LinAlg::Matrix<3, 1>& normal, double& drs, bool referencepos = false)
      {
        const int numnodes = Core::FE::num_nodes<celldistype>;
        Core::LinAlg::Matrix<3, numnodes> xyze(xyz_, true);
        if (referencepos) xyze = Core::LinAlg::Matrix<3, numnodes>(xyz_ref_, true);

        Core::LinAlg::Matrix<numnodes, 1> funct(false);
        Core::LinAlg::Matrix<2, numnodes> deriv(false);
        Core::LinAlg::Matrix<2, 2> metrictensor(false);

        Core::FE::shape_function_2D(funct, eta(0), eta(1), celldistype);

        if (celldistype != Core::FE::CellType::tri3)
        {
          Core::FE::shape_function_2D_deriv1(deriv, eta(0), eta(1), celldistype);
          Core::FE::ComputeMetricTensorForBoundaryEle<celldistype>(
              xyze, deriv, metrictensor, drs, &normal);
        }
        else
        {
          // For Tri's this method of determining the area and thus the gp-weights is more robust.
          //  It is needed for TRI's which are small/ill-conditioned but large enough to affect the
          //  simulation.
          static Core::LinAlg::Matrix<3, 1> p0(true);
          static Core::LinAlg::Matrix<3, 1> p1(true);
          static Core::LinAlg::Matrix<3, 1> p2(true);
          for (unsigned dim = 0; dim < 3; ++dim)
          {
            p0(dim) = xyze(dim, 0);
            p1(dim) = xyze(dim, 1);
            p2(dim) = xyze(dim, 2);
          }
          drs = 2.0 * (Cut::Kernel::getAreaTri(p0.A(), p1.A(), p2.A(), &normal));
        }

        x_gp_lin.Multiply(xyze, funct);

        return;
      }

      /// Reset the point with local index lid
      void ResetPos(int lid, Core::LinAlg::Matrix<3, 1> newpos)
      {
        if (lid > xyz_.numCols()) FOUR_C_THROW("Index out of range! %d > %d", lid, xyz_.numCols());

        xyz_(0, lid) = newpos(0, 0);
        xyz_(1, lid) = newpos(1, 0);
        xyz_(2, lid) = newpos(2, 0);
      }

      /*!
      \brief Compute the location of Gauss points on the background element's local coordinate
      system setting shadow = true means the mapping is done w.r. to the parent quad element from
      which elem1 is derived
       */
      template <Core::FE::CellType celldistype>
      void transform_local_coords(Element* elem1, const Core::LinAlg::Matrix<2, 1>& eta,
          Core::LinAlg::Matrix<3, 1>& x_gp_lin, Core::LinAlg::Matrix<3, 1>& normal, double& drs,
          bool shadow = false);

     protected:
      virtual Core::FE::GaussRule2D my_simple_gauss_rule() = 0;

      template <Core::FE::CellType distype>
      double my_area()
      {
        const int numnodes = Core::FE::num_nodes<distype>;
        Core::LinAlg::Matrix<3, numnodes> xyze(this->xyz_, true);
        Core::LinAlg::Matrix<numnodes, 1> funct;
        Core::LinAlg::Matrix<2, numnodes> deriv;
        Core::LinAlg::Matrix<2, 2> metrictensor;

        Core::FE::GaussRule2D gaussrule = this->my_simple_gauss_rule();
        Core::FE::IntegrationPoints2D intpoints(gaussrule);

        double area = 0;
        double drs;
        for (int i = 0; i < intpoints.nquad; ++i)
        {
          double* eta = intpoints.qxg[i];
          Core::FE::shape_function_2D(funct, eta[0], eta[1], distype);
          Core::FE::shape_function_2D_deriv1(deriv, eta[0], eta[1], distype);
          Core::FE::ComputeMetricTensorForBoundaryEle<distype>(
              xyze, deriv, metrictensor, drs, nullptr);
          if (not std::isnan(drs)) area += intpoints.qwgt[i] * drs;
        }
        return area;
      }

      /** \brief To calculate the element center of a boundary cell
       *
       *  \author shahmiri
       *  \date 07/12 */
      template <Core::FE::CellType distype>
      void my_element_center(
          Core::LinAlg::Matrix<3, 1>& center, Core::LinAlg::Matrix<3, 1>& midpoint)
      {
        const int numnodes = Core::FE::num_nodes<distype>;
        Core::LinAlg::Matrix<3, numnodes> xyze(this->xyz_, true);
        Core::LinAlg::Matrix<numnodes, 1> funct;
        Core::FE::shape_function<distype>(center, funct);
        midpoint.Multiply(xyze, funct);
      }

      /// Current position of the boundary cell
      Core::LinAlg::SerialDenseMatrix xyz_;

      /// Reference position of the boundary cell
      Core::LinAlg::SerialDenseMatrix xyz_ref_;
      Facet* facet_;
      Teuchos::RCP<Cycle> points_;
    };

    /*----------------------------------------------------------------------------*/
    /// point1 boundary cell
    class Point1BoundaryCell : public BoundaryCell
    {
     public:
      Point1BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
        /* intentionally left blank */
      }

      Core::FE::CellType Shape() const override { return Core::FE::CellType::point1; }

      int CubatureDegree() const override { return 0; }

      void DumpGmshNormal(std::ofstream& file) override;

      double Area() override { return 0.0; };

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void Normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gaussRule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> GetNormalVector() override;

      bool IsValidBoundaryCell() override { return true; };

     protected:
      Core::FE::GaussRule2D my_simple_gauss_rule() override
      {
        return Core::FE::GaussRule2D::undefined;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// point1 boundary cell
    class Line2BoundaryCell : public BoundaryCell
    {
     public:
      Line2BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
        /* intentionally left blank */
      }

      Core::FE::CellType Shape() const override { return Core::FE::CellType::line2; }

      int CubatureDegree() const override { return 4; }

      void DumpGmshNormal(std::ofstream& file) override;

      double Area() override;

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void Normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gaussRule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> GetNormalVector() override;

      bool IsValidBoundaryCell() override { return (Area() > REF_AREA_BCELL); };

     protected:
      Core::FE::GaussRule2D my_simple_gauss_rule() override
      {
        return Core::FE::GaussRule2D::undefined;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// tri3 boundary cell
    class Tri3BoundaryCell : public BoundaryCell
    {
     public:
      Tri3BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
      }

      Core::FE::CellType Shape() const override { return Core::FE::CellType::tri3; }

      int CubatureDegree() const override { return 20; }

      void DumpGmshNormal(std::ofstream& file) override;

      double Area() override;

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void Normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gaussRule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> GetNormalVector() override;

      /** \brief  A first step to validate if a boundary cell is valid.
       *
       *  \author winter
       *  \date 11/15 */
      bool IsValidBoundaryCell() override;

     protected:
      Core::FE::GaussRule2D my_simple_gauss_rule() override
      {
        return Core::FE::GaussRule2D::tri_3point;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// quad4 boundary cell
    class Quad4BoundaryCell : public BoundaryCell
    {
     public:
      Quad4BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points)
      {
      }

      Core::FE::CellType Shape() const override { return Core::FE::CellType::quad4; }

      int CubatureDegree() const override { return 20; }

      void DumpGmshNormal(std::ofstream& file) override;

      // Maybe shoelace theorem can be used here?
      double Area() override { return my_area<Core::FE::CellType::quad4>(); }

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void Normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gaussRule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> GetNormalVector() override;

      // Probably not the best way...
      bool IsValidBoundaryCell() override { return (Area() > REF_AREA_BCELL); }

     protected:
      Core::FE::GaussRule2D my_simple_gauss_rule() override
      {
        return Core::FE::GaussRule2D::quad_4point;
      }
    };

    /*----------------------------------------------------------------------------*/
    /// Irregular boundary cell generated during the cut
    class ArbitraryBoundaryCell : public BoundaryCell
    {
     public:
      ArbitraryBoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points, const Core::FE::GaussIntegration& gaussRule,
          const Core::LinAlg::Matrix<3, 1>& normal)
          : BoundaryCell(xyz, facet, points), gauss_rule_(gaussRule), normal_(normal)
      {
      }

      Core::FE::CellType Shape() const override { return Core::FE::CellType::dis_none; }

      int CubatureDegree() const override { return 0; }

      void DumpGmshNormal(std::ofstream& file) override;

      double Area() override { return 0.0; }

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void Normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gaussRule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> GetNormalVector() override;

      bool IsValidBoundaryCell() override { return (Area() > REF_AREA_BCELL); }

     protected:
      Core::FE::GaussRule2D my_simple_gauss_rule() override
      {
        return Core::FE::GaussRule2D::quad_4point;
      }

     private:
      Core::FE::GaussIntegration gauss_rule_;
      Core::LinAlg::Matrix<3, 1> normal_;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
