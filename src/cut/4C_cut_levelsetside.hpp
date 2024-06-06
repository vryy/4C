/*---------------------------------------------------------------------*/
/*! \file

\brief for intersection with an levelset, levelsetside represents the surface described by the
levelset

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_LEVELSETSIDE_HPP
#define FOUR_C_CUT_LEVELSETSIDE_HPP

#include "4C_config.hpp"

#include "4C_cut_side.hpp"

// Use derivatives of LevelSet field to determine the Cut configuration
#define USE_PHIDERIV_FOR_CUT_DETERMINATION

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    /*! \brief Class to handle level-set cut side which does not have a regular
     *  geometric shape
     *
     *  The class is a template on the problem dimension \c probdim. */
    template <int probdim>
    class LevelSetSide : public Side
    {
     public:
      /// constructor
      LevelSetSide(int sid) : Side(sid, std::vector<Node*>(), std::vector<Edge*>())
      {
        if (sid < 0) FOUR_C_THROW("The level-set side must have a positive side id!");
      }

      /// Returns the geometric shape of this side
      Core::FE::CellType Shape() const override { return Core::FE::CellType::dis_none; }

      /// element dimension
      unsigned Dim() const override
      {
        FOUR_C_THROW(
            "No dimension information for level set sides. It's "
            "likely that you can't call the calling function for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /// problem dimension
      unsigned ProbDim() const override { return probdim; }

      /// number of nodes
      unsigned NumNodes() const override
      {
        FOUR_C_THROW(
            "No number of nodes information for level set sides. It's "
            "likely that you can't call the calling function for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /// \brief Returns the topology data for the side from Shards library
      const CellTopologyData* Topology() const override
      {
        FOUR_C_THROW("No topology data for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /** Get the cut points between the levelset side and the specified edge */
      bool Cut(Mesh& mesh, Edge& edge, PointSet& cut_points) override;

      /** In the level-set case, the level-set side is the cut side and will
       *  be divided into facets. Among other places, this becomes important for
       *  the boundary integration cell creation. */
      void MakeInternalFacets(Mesh& mesh, Element* element, plain_facet_set& facets) override;

      //   virtual bool DoTriangulation() { return true; }

      bool find_ambiguous_cut_lines(
          Mesh& mesh, Element* element, Side& side, const PointSet& cut) override;

      // a levelset-side returns true
      bool IsLevelSetSide() override { return true; };


      bool find_cut_points_dispatch(Mesh& mesh, Element* element, Side& side, Edge& e) override;

     protected:
      /*! \brief is this side closer to the start-point as the other side? */
      bool is_closer_side(
          const double* startpoint_xyz, Core::Geo::Cut::Side* other, bool& is_closer) override
      {
        FOUR_C_THROW("no is_closer_side routine for level set cut side");
        exit(EXIT_FAILURE);
      }

      /*! \brief Returns the global coordinates of the nodes of this side */
      void Coordinates(double* xyze) const override
      {
        FOUR_C_THROW("no coordinates on level set cut side");
      }

      /*! \brief get all edges adjacent to given local coordinates */
      void edge_at(const double* rs, std::vector<Edge*>& edges) override
      {
        FOUR_C_THROW("no edges on level set cut side");
      }

      /*! \brief get the global coordinates on side at given local coordinates */
      void point_at(const double* rs, double* xyz) override
      {
        FOUR_C_THROW("no PointAt on level set cut side defined");
      }

      /*! \brief get global coordinates of the center of the side */
      void side_center(double* midpoint) override
      {
        FOUR_C_THROW("no SideCenter on level set cut side defined");
      }

      bool within_side(const double* xyz, double* rs, double& dist) override
      {
        FOUR_C_THROW("no WithinSide check implemented");
        exit(EXIT_FAILURE);
      }

      bool ray_cut(const double* p1_xyz, const double* p2_xyz, double* rs, double& line_xi) override
      {
        FOUR_C_THROW("no RayCut with level set cut side implemented");
        exit(EXIT_FAILURE);
      }

      /*! \brief Calculates the local coordinates (rsd) with respect to the element shape
       *  from its global coordinates (xyz), return TRUE if successful. The last coordinate
       *  of \c rsd is the distance of the n-dimensional point \c xyz to the embedded
       *  side. */
      bool local_coordinates(const double* xyz, double* rst, bool allow_dist = false,
          double tol = POSITIONTOL) override
      {
        FOUR_C_THROW("no local coordinates on level set cut side");
        exit(EXIT_FAILURE);
      }

      /*! \brief get local coordinates (rst) with respect to the element shape
       * for all the corner points */
      void local_corner_coordinates(double* rst_corners) override
      {
        FOUR_C_THROW("no local coordinates of corner points on level set cut side");
      }

      /*! \brief Calculates the normal vector with respect to the element shape
       *  at local coordinates \c rs */
      void normal(const double* xsi, double* normal, bool unitnormal = true) override
      {
        FOUR_C_THROW("no normal vector on level set cut side implemented");
      }

      /* \brief Calculates a Basis of two tangential vectors (non-orthogonal!) and
       * the normal vector with respect to the element shape at local coordinates rs,
       * basis vectors have norm=1 */
      void basis_at_center(double* t1, double* t2, double* n) override
      {
        FOUR_C_THROW("no basis_at_center on level set cut side implemented");
      }

      /* \brief Calculates a Basis of two tangential vectors (non-orthogonal!) and
       * the normal vector with respect to the element shape at local coordinates rs,
       * basis vectors have norm=1. */
      void basis(const double* xsi, double* t1, double* t2, double* n) override
      {
        FOUR_C_THROW("no Basis on level set cut side implemented");
      }
    };  // class LevelSetSide

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
