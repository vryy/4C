/*---------------------------------------------------------------------*/
/*! \file

\brief Sidehandle represents a side original loaded into the cut, internal it can be split into
subsides

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#ifndef CUT_SIDEHANDLE_H
#define CUT_SIDEHANDLE_H

#include "cut_side.H"
#include "cut_boundarycell.H"

namespace GEO
{
  namespace CUT
  {
    class Node;
    class Side;
    class Mesh;

    /*!
    \brief Outside world interface to the side. This breaks the quadratic side into linear sides
     */
    class SideHandle
    {
     public:
      virtual ~SideHandle() {}

      /*!
      \brief Get the shape of this sides
       */
      virtual DRT::Element::DiscretizationType Shape() = 0;

      /*!
      \brief Get the coordinates of the nodes of this side
       */
      virtual void Coordinates(Epetra_SerialDenseMatrix& xyze) = 0;

      /*!
      \brief Get the local coordinates "rst"from the global coordinates "xyz" with respect to this
      side. Since side is 2D, the local coordinates will have only two coordinates
      */
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<2, 1>& rst) = 0;

      /*!
      \brief For a quadratic side, get the resulting linear sides
       */
      virtual void CollectSides(plain_side_set& sides) = 0;


      /*!
      \brief Gets all facets of a side
       */
      virtual void Facets(std::vector<Facet*>& facets) = 0;

      /*!
      \brief Get the Gaussian rule projected on the side. Not used now
       */
      template <DRT::Element::DiscretizationType distype>
      Teuchos::RCP<DRT::UTILS::GaussPoints> CreateProjected(BoundaryCell* bc)
      {
        const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

        LINALG::Matrix<2, nen> xie;

        const std::vector<GEO::CUT::Point*>& cpoints = bc->Points();
        if (cpoints.size() != nen) throw std::runtime_error("non-matching number of points");

        for (unsigned i = 0; i < nen; ++i)
        {
          GEO::CUT::Point* p = cpoints[i];
          const LINALG::Matrix<2, 1>& xi = LocalCoordinates(p);
          std::copy(xi.A(), xi.A() + 2, &xie(0, i));
        }

        Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
            DRT::UTILS::GaussIntegration::CreateProjected<distype>(
                xie, bc->CubatureDegree(Shape()));
        return gp;
      }

      /*!
      \brief Get the local coordinates of point p with respect to this side. Since side is 2D, the
      local coordinates will have only two coordinates
       */
      const LINALG::Matrix<2, 1>& LocalCoordinates(Point* p)
      {
        std::map<Point*, LINALG::Matrix<2, 1>>::iterator i = local_coordinates_.find(p);
        if (i != local_coordinates_.end())
        {
          return i->second;
        }
        LINALG::Matrix<2, 1>& rst = local_coordinates_[p];
        LINALG::Matrix<3, 1> xyz;
        p->Coordinates(xyz.A());
        LocalCoordinates(xyz, rst);
        return rst;
      }

     private:
      std::map<Point*, LINALG::Matrix<2, 1>> local_coordinates_;
    };

    /// linear side handle
    class LinearSideHandle : public SideHandle
    {
     public:
      LinearSideHandle() : side_(NULL) {}

      explicit LinearSideHandle(Side* s) : side_(s) {}

      virtual DRT::Element::DiscretizationType Shape() { return side_->Shape(); }

      virtual void Coordinates(Epetra_SerialDenseMatrix& xyze)
      {
        xyze.Reshape(3, side_->Nodes().size());
        side_->Coordinates(xyze.A());
      }

      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<2, 1>& rs)
      {
        LINALG::Matrix<3, 1> rst;
        side_->LocalCoordinates(xyz, rst);
        rs(0) = rst(0);
        rs(1) = rst(1);
      }

      virtual void CollectSides(plain_side_set& sides) { sides.insert(side_); }

      /// gets all facets of a linear side
      virtual void Facets(std::vector<Facet*>& facets)
      {
        std::vector<Facet*> sidefacets = side_->Facets();
        for (std::vector<Facet*>::iterator i = sidefacets.begin(); i != sidefacets.end(); ++i)
        {
          Facet* sidefacet = *i;
          facets.push_back(sidefacet);
        }
      }


     private:
      Side* side_;
    };

    /// quadratic side handle
    class QuadraticSideHandle : public SideHandle
    {
     public:
      virtual void Coordinates(Epetra_SerialDenseMatrix& xyze)
      {
        xyze.Reshape(3, nodes_.size());
        for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
        {
          Node* n = *i;
          n->Coordinates(&xyze(0, i - nodes_.begin()));
        }
      }

      virtual void CollectSides(plain_side_set& sides)
      {
        std::copy(subsides_.begin(), subsides_.end(), std::inserter(sides, sides.begin()));
      }

      /// gets all facets of a quadratic side
      virtual void Facets(std::vector<Facet*>& facets)
      {
        for (std::vector<Side*>::iterator i = subsides_.begin(); i != subsides_.end(); ++i)
        {
          Side* subside = *i;
          std::vector<Facet*> sidefacets = subside->Facets();
          for (std::vector<Facet*>::iterator i = sidefacets.begin(); i != sidefacets.end(); ++i)
          {
            Facet* sidefacet = *i;
            facets.push_back(sidefacet);
          }
        }
      }

     protected:
      std::vector<Side*> subsides_;
      std::vector<Node*> nodes_;
    };

    /// tri6 side handle
    class Tri6SideHandle : public QuadraticSideHandle
    {
     public:
      Tri6SideHandle(Mesh& mesh, int sid, const std::vector<int>& nids);

      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::tri6; }

      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<2, 1>& rst);
    };

    /// quad4 side handle
    /*!
     * quad4 need to be split in 3 tri3 in order to avoid subtle ambiguities.
     */
    class Quad4SideHandle : public QuadraticSideHandle
    {
     public:
      Quad4SideHandle(Mesh& mesh, int sid, const std::vector<int>& nids);

      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::quad4; }

      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<2, 1>& rst);
    };

    /// quad8 side handle
    class Quad8SideHandle : public QuadraticSideHandle
    {
     public:
      Quad8SideHandle(Mesh& mesh, int sid, const std::vector<int>& nids, bool iscutside = true);

      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::quad8; }

      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<2, 1>& rst);
    };

    /// quad9 side handle
    class Quad9SideHandle : public QuadraticSideHandle
    {
     public:
      Quad9SideHandle(Mesh& mesh, int sid, const std::vector<int>& nids, bool iscutside = true);

      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::quad9; }

      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<2, 1>& rst);
    };

  }  // namespace CUT
}  // namespace GEO

#endif
