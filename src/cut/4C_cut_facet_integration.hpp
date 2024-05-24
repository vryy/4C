/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over the facet for both volumecell facets and for boundarycells
equations

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_CUT_FACET_INTEGRATION_HPP
#define FOUR_C_CUT_FACET_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_cut_element.hpp"
#include "4C_cut_enum.hpp"
#include "4C_cut_mesh.hpp"

FOUR_C_NAMESPACE_OPEN

// #define DIRECTDIV_EXTENDED_DEBUG_OUTPUT
// #define TRIANGULATE_ALL_FACETS_FOR_DIVERGENCECELLS

namespace CORE::GEO
{
  namespace CUT
  {
    class Element;
    class Facet;

    /*!
    \brief This class performs the integration of base functions over the facet. The points of the
    facet should be arranged in anti-clockwise manner when looking the facet away from the volume
    this ensures outward normal vector when divergence theorem is used
    */
    class FacetIntegration
    {
     public:
      FacetIntegration(Facet *face1, Element *element1,
          const CORE::GEO::CUT::Point::PointPosition posi, bool bcellInt, bool global)
          : face1_(face1),             // facet under consideration
            elem1_(element1),          // the element for which the facet is a part of
            position_(posi),           // position
            bcell_int_(bcellInt),      //"true" if it is boundarycell integration
            ordering_computed_(false)  // whether cw or acw ordering of vertices computed
      {
      }

      /*!
      \brief Select the base function to be integrated
      */
      void set_integ_number(int inte_num) { inte_num_ = inte_num; }

      /*!
      \brief Performs the integration of a function over the facet
      */
      double integrate_facet();

      /*!
      \brief Computes the equation of the plane that contains this facet
      */
      std::vector<double> equation_plane(const std::vector<std::vector<double>> &cornersLocal);

      /*!
      \brief Returns the equation of plane that contains this facet
      */
      std::vector<double> get_equation() { return eqn_plane_; }

      /*!
      \brief Return whether the vertices numbering of the facet is clockwise
      */
      bool IsClockwiseOrdering();

      /*!
      \brief Generate Gaussian points over the considered facet by triangulating it. This is used
      when DirectDivergence option is used for Gauss point generation
      */
      void divergence_integration_rule(
          Mesh &mesh, Teuchos::RCP<CORE::FE::CollectedGaussPoints> &cgp);

      /*!
      \brief Generate Gaussian points over the considered facet by triangulating it. This is used
      when DirectDivergence option is used for Gauss point generation
      */
      void divergence_integration_rule_new(
          Mesh &mesh, Teuchos::RCP<CORE::FE::CollectedGaussPoints> &cgp);

     private:
      /*!
      \brief Check whether the vertices numbering of the facet is clockwise
      */
      void is_clockwise(const std::vector<double> &eqn_plane,
          const std::vector<std::vector<double>> &cornersLocal);

      /*
      \brief Compute the function which replaces "x" when projecting the facet into coordinate plane
      */
      std::vector<double> compute_alpha(
          std::vector<double> &eqn_plane, CORE::GEO::CUT::ProjectionDirection intType);

      /*!
      \brief Get normal of the considered facet in a particular coordinate direction defined by
      intType
      */
      double get_normal(CORE::GEO::CUT::ProjectionDirection intType);

      /*!
      \brief Perform integration of base functions over boundarycells
      */
      void boundary_facet_integration(const std::vector<std::vector<double>> &cornersLocal,
          double &facet_integ, CORE::GEO::CUT::ProjectionDirection intType);

      /*!
      \brief Generate boundary cells for the considered facet. May need to perform triangulatio
      */
      void generate_divergence_cells(
          bool divergenceRule, Mesh &mesh, std::list<Teuchos::RCP<BoundaryCell>> &divCells);

      /*!
      \brief Generate boundary cells for the considered facet. May need to perform triangulatio
      */
      void generate_divergence_cells_new(bool divergenceRule, Mesh &mesh,
          std::list<Teuchos::RCP<BoundaryCell>> &divCells,
          const std::vector<Point *> &cornersGlobal);


      /*!
      \brief Temporarily create Tri3 cell. This is not stored in Mesh
      */
      void temporary_tri3(
          const std::vector<Point *> &corners, std::list<Teuchos::RCP<BoundaryCell>> &divCells);

      /*!
      \brief Temporarily create Quad4 cell. This is not stored in Mesh
      */
      void temporary_quad4(
          const std::vector<Point *> &corners, std::list<Teuchos::RCP<BoundaryCell>> &divCells);

      //! considered facet
      Facet *face1_;

      //! background element which was cut to produce this facet
      Element *elem1_;

      //! position of the facet
      const CORE::GEO::CUT::Point::PointPosition position_;

      //! True for boundarycell integration
      bool bcell_int_;

      //! set the base function to be integrated
      int inte_num_;

      //! True if nodes of facet are arranged to give inward normal
      bool clockwise_;

      //! whether the clockwise or ACW ordering of facet computed already
      bool ordering_computed_;

      //! equation of plane that contains the facet
      std::vector<double> eqn_plane_;

      std::list<Teuchos::RCP<BoundaryCell>> boundarycells_;
    };
  }  // namespace CUT
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif
