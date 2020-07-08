/*---------------------------------------------------------------------*/
/*! \file

\brief for intersection with an levelset, levelsetside represents the surface described by the
levelset

\level 2


*----------------------------------------------------------------------*/

#ifndef CUT_LEVELSETSIDE_H
#define CUT_LEVELSETSIDE_H

#include "cut_side.H"

// Use derivatives of LevelSet field to determine the Cut configuration
#define USE_PHIDERIV_FOR_CUT_DETERMINATION

namespace GEO
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

      /// destructor
      virtual ~LevelSetSide(){};

      /// Returns the geometric shape of this side
      virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::dis_none; }

      /// element dimension
      virtual unsigned Dim() const
      {
        dserror(
            "No dimension information for level set sides. It's "
            "likely that you can't call the calling function for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /// problem dimension
      virtual unsigned ProbDim() const { return probdim; }

      /// number of nodes
      virtual unsigned NumNodes() const
      {
        dserror(
            "No number of nodes information for level set sides. It's "
            "likely that you can't call the calling function for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /// \brief Returns the topology data for the side from Shards library
      virtual const CellTopologyData* Topology() const
      {
        dserror("No topology data for level-set sides!");
        exit(EXIT_FAILURE);
      }

      /** Get the cut points between the levelset side and the specified edge */
      virtual bool Cut(Mesh& mesh, Edge& edge, PointSet& cut_points);

      //   virtual void MakeSideCutFacets( Mesh & mesh, Element * element, plain_facet_set & facets
      //   );

      /** In the level-set case, the level-set side is the cut side and will
       *  be divided into facets. Among other places, this becomes important for
       *  the boundary integration cell creation. */
      virtual void MakeInternalFacets(Mesh& mesh, Element* element, plain_facet_set& facets);

      //   virtual bool DoTriangulation() { return true; }

      virtual bool FindAmbiguousCutLines(
          Mesh& mesh, Element* element, Side& side, const PointSet& cut);

      // a levelset-side returns true
      virtual bool IsLevelSetSide() { return true; };


      virtual bool FindCutPointsDispatch(Mesh& mesh, Element* element, Side& side, Edge& e);

     protected:
      /*! \brief is this side closer to the start-point as the other side? */
      virtual bool IsCloserSide(
          const double* startpoint_xyz, GEO::CUT::Side* other, bool& is_closer)
      {
        dserror("no IsCloserSide routine for level set cut side");
        exit(EXIT_FAILURE);
      }

      /*! \brief Returns the global coordinates of the nodes of this side */
      virtual void Coordinates(double* xyze) const
      {
        dserror("no coordinates on level set cut side");
      }

      /*! \brief get all edges adjacent to given local coordinates */
      virtual void EdgeAt(const double* rs, std::vector<Edge*>& edges)
      {
        dserror("no edges on level set cut side");
      }

      /*! \brief get the global coordinates on side at given local coordinates */
      virtual void PointAt(const double* rs, double* xyz)
      {
        dserror("no PointAt on level set cut side defined");
      }

      /*! \brief get global coordinates of the center of the side */
      virtual void SideCenter(double* midpoint)
      {
        dserror("no SideCenter on level set cut side defined");
      }

      virtual bool WithinSide(const double* xyz, double* rs, double& dist)
      {
        dserror("no WithinSide check implemented");
        exit(EXIT_FAILURE);
      }

      virtual bool RayCut(const double* p1_xyz, const double* p2_xyz, double* rs, double& line_xi)
      {
        dserror("no RayCut with level set cut side implemented");
        exit(EXIT_FAILURE);
      }

      /*! \brief Calculates the local coordinates (rsd) with respect to the element shape
       *  from its global coordinates (xyz), return TRUE if successful. The last coordinate
       *  of \c rsd is the distance of the n-dimensional point \c xyz to the embedded
       *  side. */
      virtual bool LocalCoordinates(
          const double* xyz, double* rst, bool allow_dist = false, double tol = POSITIONTOL)
      {
        dserror("no local coordinates on level set cut side");
        exit(EXIT_FAILURE);
      }

      /*! \brief get local coordinates (rst) with respect to the element shape
       * for all the corner points */
      virtual void LocalCornerCoordinates(double* rst_corners)
      {
        dserror("no local coordinates of corner points on level set cut side");
      }

      /*! \brief Calculates the normal vector with respect to the element shape
       *  at local coordinates \c rs */
      virtual void Normal(const double* xsi, double* normal, bool unitnormal = true)
      {
        dserror("no normal vector on level set cut side implemented");
      }

      /* \brief Calculates a Basis of two tangential vectors (non-orthogonal!) and
       * the normal vector with respect to the element shape at local coordinates rs,
       * basis vectors have norm=1 */
      virtual void BasisAtCenter(double* t1, double* t2, double* n)
      {
        dserror("no BasisAtCenter on level set cut side implemented");
      }

      /* \brief Calculates a Basis of two tangential vectors (non-orthogonal!) and
       * the normal vector with respect to the element shape at local coordinates rs,
       * basis vectors have norm=1. */
      virtual void Basis(const double* xsi, double* t1, double* t2, double* n)
      {
        dserror("no Basis on level set cut side implemented");
      }
    };  // class LevelSetSide

  }  // namespace CUT
}  // namespace GEO

#endif
