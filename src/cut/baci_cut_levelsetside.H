/*---------------------------------------------------------------------*/
/*! \file

\brief for intersection with an levelset, levelsetside represents the surface described by the
levelset

\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_CUT_LEVELSETSIDE_H
#define BACI_CUT_LEVELSETSIDE_H

#include "baci_config.H"

#include "baci_cut_side.H"

// Use derivatives of LevelSet field to determine the Cut configuration
#define USE_PHIDERIV_FOR_CUT_DETERMINATION

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
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
        if (sid < 0) dserror("The level-set side must have a positive side id!");
      }

      /// Returns the geometric shape of this side
      CORE::FE::CellType Shape() const override { return CORE::FE::CellType::dis_none; }

      /// element dimension
      unsigned Dim() const override
      {
        dserror(
            "No dimension information for level set sides. It's "
            "likely that you can't call the calling function for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /// problem dimension
      unsigned ProbDim() const override { return probdim; }

      /// number of nodes
      unsigned NumNodes() const override
      {
        dserror(
            "No number of nodes information for level set sides. It's "
            "likely that you can't call the calling function for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /// \brief Returns the topology data for the side from Shards library
      const CellTopologyData* Topology() const override
      {
        dserror("No topology data for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /** Get the cut points between the levelset side and the specified edge */
      bool Cut(Mesh& mesh, Edge& edge, PointSet& cut_points) override;

      /** In the level-set case, the level-set side is the cut side and will
       *  be divided into facets. Among other places, this becomes important for
       *  the boundary integration cell creation. */
      void MakeInternalFacets(Mesh& mesh, Element* element, plain_facet_set& facets) override;

      //   virtual bool DoTriangulation() { return true; }

      bool FindAmbiguousCutLines(
          Mesh& mesh, Element* element, Side& side, const PointSet& cut) override;

      // a levelset-side returns true
      bool IsLevelSetSide() override { return true; };


      bool FindCutPointsDispatch(Mesh& mesh, Element* element, Side& side, Edge& e) override;

     protected:
      /*! \brief is this side closer to the start-point as the other side? */
      bool IsCloserSide(
          const double* startpoint_xyz, CORE::GEO::CUT::Side* other, bool& is_closer) override
      {
        dserror("no IsCloserSide routine for level set cut side");
        exit(EXIT_FAILURE);
      }

      /*! \brief Returns the global coordinates of the nodes of this side */
      void Coordinates(double* xyze) const override
      {
        dserror("no coordinates on level set cut side");
      }

      /*! \brief get all edges adjacent to given local coordinates */
      void EdgeAt(const double* rs, std::vector<Edge*>& edges) override
      {
        dserror("no edges on level set cut side");
      }

      /*! \brief get the global coordinates on side at given local coordinates */
      void PointAt(const double* rs, double* xyz) override
      {
        dserror("no PointAt on level set cut side defined");
      }

      /*! \brief get global coordinates of the center of the side */
      void SideCenter(double* midpoint) override
      {
        dserror("no SideCenter on level set cut side defined");
      }

      bool WithinSide(const double* xyz, double* rs, double& dist) override
      {
        dserror("no WithinSide check implemented");
        exit(EXIT_FAILURE);
      }

      bool RayCut(const double* p1_xyz, const double* p2_xyz, double* rs, double& line_xi) override
      {
        dserror("no RayCut with level set cut side implemented");
        exit(EXIT_FAILURE);
      }

      /*! \brief Calculates the local coordinates (rsd) with respect to the element shape
       *  from its global coordinates (xyz), return TRUE if successful. The last coordinate
       *  of \c rsd is the distance of the n-dimensional point \c xyz to the embedded
       *  side. */
      bool LocalCoordinates(const double* xyz, double* rst, bool allow_dist = false,
          double tol = POSITIONTOL) override
      {
        dserror("no local coordinates on level set cut side");
        exit(EXIT_FAILURE);
      }

      /*! \brief get local coordinates (rst) with respect to the element shape
       * for all the corner points */
      void LocalCornerCoordinates(double* rst_corners) override
      {
        dserror("no local coordinates of corner points on level set cut side");
      }

      /*! \brief Calculates the normal vector with respect to the element shape
       *  at local coordinates \c rs */
      void Normal(const double* xsi, double* normal, bool unitnormal = true) override
      {
        dserror("no normal vector on level set cut side implemented");
      }

      /* \brief Calculates a Basis of two tangential vectors (non-orthogonal!) and
       * the normal vector with respect to the element shape at local coordinates rs,
       * basis vectors have norm=1 */
      void BasisAtCenter(double* t1, double* t2, double* n) override
      {
        dserror("no BasisAtCenter on level set cut side implemented");
      }

      /* \brief Calculates a Basis of two tangential vectors (non-orthogonal!) and
       * the normal vector with respect to the element shape at local coordinates rs,
       * basis vectors have norm=1. */
      void Basis(const double* xsi, double* t1, double* t2, double* n) override
      {
        dserror("no Basis on level set cut side implemented");
      }
    };  // class LevelSetSide

  }  // namespace CUT
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif
