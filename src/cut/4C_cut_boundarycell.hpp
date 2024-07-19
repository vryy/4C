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
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_integration.hpp"

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
          const std::vector<Point*>& points, int cubature_degree);

      /// destructor
      virtual ~BoundaryCell() = default;
      /*!
      \brief Returns the shape of the boundarycell
       */

      virtual Core::FE::CellType shape() const = 0;

      /*!
      \brief Writes the geometry of boundarycell, and the constant scalar "value" in GMSH format
       */
      void dump_gmsh(std::ofstream& file, int* value = nullptr);

      /*!
      \brief Writes the geometry of boundarycell, and the constant scalar "value" in GMSH format
       */
      virtual void dump_gmsh_normal(std::ofstream& file) = 0;

      /*!
      \brief Returns the area of tri3 and quad4 boundarycell
       */
      virtual double area() = 0;

      /*!
      \brief Returns the center of tri3 and quad4 boundarycell
       */
      virtual void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) = 0;

      /*!
      \brief Get the outward normal vector for the tri3 and quad4 boundarycell
       */
      virtual void normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const = 0;

      /*!
      \brief Get the corner points of the boundarycell in Cylce for geometrical operations
       */
      const Cycle& point_cycle() const { return *points_; }

      /*!
      \brief Get the corner points of the boundarycell as vector of Points
       */
      const std::vector<Point*>& points() const;

      /*!
      \brief Get the coodinates of the corner points of boundarycell
       */
      const Core::LinAlg::SerialDenseMatrix& coordinates() const { return xyz_; }

      /*!
      \brief Move the corner points of the boundary cell by a specific offset
       */
      template <Core::FE::CellType celldistype>
      void assign_offset(int idx, double offset)
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
      std::vector<std::vector<double>> coordinates_v();

      /*!
      \brief Get the Facet which represent the boundarycell
       */
      Facet* get_facet() { return facet_; }
      const Facet* get_facet() const { return facet_; }

      /*!
      \brief Delete all the points of this boundarycell
       */
      void clear();

      bool is_valid() const;

      /*!
      \brief Function to test if the distance between points is within point tolerance.
       */
      virtual bool is_valid_boundary_cell() = 0;

      /*!
      \brief Get the Gaussian integration rule for the boundarycell
       */
      virtual Core::FE::GaussIntegration gauss_rule(int cubaturedegree) = 0;

      /*!
      \brief Get the Gaussian integration rule for the boundarycell
       */
      virtual Core::FE::GaussIntegration gauss_rule() { return gauss_rule(get_cubature_degree()); }
      /*!
      \brief Get the normal vector for the arbitrary boundarycell alone
       */
      virtual Core::LinAlg::Matrix<3, 1> get_normal_vector() = 0;

      /*! \brief Print the corner points on screen
       *
       *  just for debugging
       *
       *  \author sudhakar
       *  \date 10/14
       */
      void print(std::ostream& stream);
      void print() { print(std::cout); }

      /*!
      \brief Set a different cubature degree
       */
      void set_new_cubature_degree(int cubature_degree) { cubature_degree_ = cubature_degree; }

      /*!
      \brief Returns the cubature degree to generate quadrature rule for the cell
       */
      int get_cubature_degree() { return cubature_degree_; }

      /*!
      \brief Computes the location of Gauss points on the boundarycell (x_gp_lin) from the standard
      Gauss point location (eta) corresponding to the shape of the boundarycell. Also computes the
      factor to be multiplied with integration weight and normal vector of the boundarycell
       */
      template <Core::FE::CellType celldistype>
      void transform(const Core::LinAlg::Matrix<2, 1>& eta, Core::LinAlg::Matrix<3, 1>& x_gp_lin,
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
          drs = 2.0 * (Cut::Kernel::getAreaTri(p0.data(), p1.data(), p2.data(), &normal));
        }

        x_gp_lin.multiply(xyze, funct);

        return;
      }

      /// Reset the point with local index lid
      void reset_pos(int lid, Core::LinAlg::Matrix<3, 1> newpos)
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

      /*!
      \brief Get the global id of this boundary cell
       */
      int get_global_boundary_cell_id();

      /*!
      \brief Set a pointer to the background element of this boundary cell
       */
      void set_background_ele_ptr(Teuchos::RCP<Geo::Cut::Element> background_ele_ptr)
      {
        background_ele_ptr_ = background_ele_ptr;
      }

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
        midpoint.multiply(xyze, funct);
      }

      /// Current position of the boundary cell
      Core::LinAlg::SerialDenseMatrix xyz_;

      /// Reference position of the boundary cell
      Core::LinAlg::SerialDenseMatrix xyz_ref_;
      Facet* facet_;
      Teuchos::RCP<Cycle> points_;

      /// Pointer to the background element of this boundary cell
      Teuchos::RCP<Geo::Cut::Element> background_ele_ptr_;

      /// Cubature degree
      int cubature_degree_;
    };

    /*----------------------------------------------------------------------------*/
    /// point1 boundary cell
    class Point1BoundaryCell : public BoundaryCell
    {
     public:
      /// Default value for the cubature degree of a Point1BoundaryCell
      static constexpr int cubature_degree_default = 0;

      Point1BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points, cubature_degree_default)
      {
      }

      Core::FE::CellType shape() const override { return Core::FE::CellType::point1; }

      void dump_gmsh_normal(std::ofstream& file) override;

      double area() override { return 0.0; };

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gauss_rule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> get_normal_vector() override;

      bool is_valid_boundary_cell() override { return true; };

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
      /// Default value for the cubature degree of a Line2BoundaryCell
      static constexpr int cubature_degree_default = 4;

      Line2BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points, cubature_degree_default)
      {
      }

      Core::FE::CellType shape() const override { return Core::FE::CellType::line2; }

      void dump_gmsh_normal(std::ofstream& file) override;

      double area() override;

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gauss_rule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> get_normal_vector() override;

      bool is_valid_boundary_cell() override { return (area() > REF_AREA_BCELL); };

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
      /// Default value for the cubature degree of a Tri3BoundaryCell
      static constexpr int cubature_degree_default = 20;

      Tri3BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points, cubature_degree_default)
      {
      }

      Core::FE::CellType shape() const override { return Core::FE::CellType::tri3; }

      void dump_gmsh_normal(std::ofstream& file) override;

      double area() override;

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gauss_rule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> get_normal_vector() override;

      /** \brief  A first step to validate if a boundary cell is valid.
       *
       *  \author winter
       *  \date 11/15 */
      bool is_valid_boundary_cell() override;

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
      /// Default value for the cubature degree of a Quad4BoundaryCell
      static constexpr int cubature_degree_default = 20;

      Quad4BoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points)
          : BoundaryCell(xyz, facet, points, cubature_degree_default)
      {
      }

      Core::FE::CellType shape() const override { return Core::FE::CellType::quad4; }

      void dump_gmsh_normal(std::ofstream& file) override;

      // Maybe shoelace theorem can be used here?
      double area() override { return my_area<Core::FE::CellType::quad4>(); }

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gauss_rule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> get_normal_vector() override;

      // Probably not the best way...
      bool is_valid_boundary_cell() override { return (area() > REF_AREA_BCELL); }

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
      /// Default value for the cubature degree of an ArbitraryBoundaryCell
      static constexpr int cubature_degree_default = 0;

      ArbitraryBoundaryCell(const Core::LinAlg::SerialDenseMatrix& xyz, Facet* facet,
          const std::vector<Point*>& points, const Core::FE::GaussIntegration& gaussRule,
          const Core::LinAlg::Matrix<3, 1>& normal)
          : BoundaryCell(xyz, facet, points, cubature_degree_default),
            gauss_rule_(gaussRule),
            normal_(normal)
      {
      }

      Core::FE::CellType shape() const override { return Core::FE::CellType::dis_none; }

      void dump_gmsh_normal(std::ofstream& file) override;

      double area() override { return 0.0; }

      void element_center(Core::LinAlg::Matrix<3, 1>& midpoint) override;

      void normal(
          const Core::LinAlg::Matrix<2, 1>& xsi, Core::LinAlg::Matrix<3, 1>& normal) const override;

      Core::FE::GaussIntegration gauss_rule(int cubaturedegree) override;

      Core::LinAlg::Matrix<3, 1> get_normal_vector() override;

      bool is_valid_boundary_cell() override { return (area() > REF_AREA_BCELL); }

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
