/*---------------------------------------------------------------------*/
/*! \file

\brief Sidehandle represents a side original loaded into the cut, internal it can be split into
subsides

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_SIDEHANDLE_HPP
#define FOUR_C_CUT_SIDEHANDLE_HPP

#include "4C_config.hpp"

#include "4C_cut_boundarycell.hpp"
#include "4C_cut_side.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
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
      virtual ~SideHandle() = default;
      /*!
      \brief Get the shape of this sides
       */
      virtual Core::FE::CellType shape() = 0;

      /*!
      \brief Get the coordinates of the nodes of this side
       */
      virtual void coordinates(Core::LinAlg::SerialDenseMatrix& xyze) = 0;

      /*!
      \brief Get the local coordinates "rst"from the global coordinates "xyz" with respect to this
      side. Since side is 2D, the local coordinates will have only two coordinates
      */
      virtual void local_coordinates(
          const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<2, 1>& rst) = 0;

      /*!
      \brief For a quadratic side, get the resulting linear sides
       */
      virtual void collect_sides(plain_side_set& sides) = 0;


      /*!
      \brief Gets all facets of a side
       */
      virtual void facets(std::vector<Facet*>& facets) = 0;

      /*!
      \brief Get the Gaussian rule projected on the side. Not used now
       */
      template <Core::FE::CellType distype>
      Teuchos::RCP<Core::FE::GaussPoints> create_projected(BoundaryCell* bc)
      {
        const unsigned nen = Core::FE::num_nodes<distype>;

        Core::LinAlg::Matrix<2, nen> xie;

        const std::vector<Core::Geo::Cut::Point*>& cpoints = bc->points();
        if (cpoints.size() != nen) FOUR_C_THROW("non-matching number of points");

        for (unsigned i = 0; i < nen; ++i)
        {
          Core::Geo::Cut::Point* p = cpoints[i];
          const Core::LinAlg::Matrix<2, 1>& xi = local_coordinates(p);
          std::copy(xi.data(), xi.data() + 2, &xie(0, i));
        }

        Teuchos::RCP<Core::FE::GaussPoints> gp =
            Core::FE::GaussIntegration::create_projected<distype>(xie, bc->get_cubature_degree());
        return gp;
      }

      /*!
      \brief Get the local coordinates of point p with respect to this side. Since side is 2D, the
      local coordinates will have only two coordinates
       */
      const Core::LinAlg::Matrix<2, 1>& local_coordinates(Point* p)
      {
        std::map<Point*, Core::LinAlg::Matrix<2, 1>>::iterator i = local_coordinates_.find(p);
        if (i != local_coordinates_.end())
        {
          return i->second;
        }
        Core::LinAlg::Matrix<2, 1>& rst = local_coordinates_[p];
        Core::LinAlg::Matrix<3, 1> xyz;
        p->coordinates(xyz.data());
        local_coordinates(xyz, rst);
        return rst;
      }

      /// Remove the SubSidePointer of given side from this Sidehandle
      virtual void remove_sub_side_pointer(const Side* side)
      {
        FOUR_C_THROW("remove_sub_side_pointer: Not available in base class!");
      }

      /// Add the SubSidePointer of given side to this Sidehandle
      virtual void add_sub_side_pointer(Side* side)
      {
        FOUR_C_THROW("AddSubSidePointer: Not available in base class!");
      }

      /// Add the SubSide in to the unphysical list
      virtual void mark_sub_sideunphysical(Side* side)
      {
        FOUR_C_THROW("SetSubSidePointertounphysical: Not available in base class!");
      }

      /// Is this side and unphysical subside
      virtual bool isunphysical_sub_side(Side* side)
      {
        FOUR_C_THROW("IsunphysicalSubSide: Not available in base class!");
        return false;  // dummy
      }

      /// Does this sidehandle have unphysical subsides
      virtual bool hasunphysical_sub_side()
      {
        FOUR_C_THROW("hasunphysical_sub_side: Not available in base class!");
        return false;  // dummy
      }

      virtual const std::vector<Node*>& get_nodes() const
      {
        FOUR_C_THROW("GetNodes: Not available in base class!");
        static const std::vector<Node*> dummy;
        return dummy;  // dummy
      }

     private:
      std::map<Point*, Core::LinAlg::Matrix<2, 1>> local_coordinates_;
    };

    /// linear side handle
    class LinearSideHandle : public SideHandle
    {
     public:
      LinearSideHandle() : side_(nullptr) {}

      explicit LinearSideHandle(Side* s) : side_(s) {}

      Core::FE::CellType shape() override { return side_->shape(); }

      void coordinates(Core::LinAlg::SerialDenseMatrix& xyze) override
      {
        xyze.reshape(3, side_->nodes().size());
        side_->coordinates(xyze.values());
      }

      void local_coordinates(
          const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<2, 1>& rs) override
      {
        Core::LinAlg::Matrix<3, 1> rst;
        side_->local_coordinates(xyz, rst);
        rs(0) = rst(0);
        rs(1) = rst(1);
      }

      void collect_sides(plain_side_set& sides) override { sides.insert(side_); }

      /// gets all facets of a linear side
      void facets(std::vector<Facet*>& facets) override
      {
        std::vector<Facet*> sidefacets = side_->facets();
        for (std::vector<Facet*>::iterator i = sidefacets.begin(); i != sidefacets.end(); ++i)
        {
          Facet* sidefacet = *i;
          facets.push_back(sidefacet);
        }
      }

      /// Get the nodes of the Sidehandle
      const std::vector<Node*>& get_nodes() const override { return side_->nodes(); }


     private:
      Side* side_;
    };

    /// quadratic side handle
    class QuadraticSideHandle : public SideHandle
    {
     public:
      void coordinates(Core::LinAlg::SerialDenseMatrix& xyze) override
      {
        xyze.reshape(3, nodes_.size());
        for (std::vector<Node*>::iterator i = nodes_.begin(); i != nodes_.end(); ++i)
        {
          Node* n = *i;
          n->coordinates(&xyze(0, i - nodes_.begin()));
        }
      }

      void collect_sides(plain_side_set& sides) override
      {
        std::copy(subsides_.begin(), subsides_.end(), std::inserter(sides, sides.begin()));
      }

      /// gets all facets of a quadratic side
      void facets(std::vector<Facet*>& facets) override
      {
        for (std::vector<Side*>::iterator i = subsides_.begin(); i != subsides_.end(); ++i)
        {
          Side* subside = *i;
          std::vector<Facet*> sidefacets = subside->facets();
          for (std::vector<Facet*>::iterator i = sidefacets.begin(); i != sidefacets.end(); ++i)
          {
            Facet* sidefacet = *i;
            facets.push_back(sidefacet);
          }
        }
      }

      /// Remove the SubSidePointer of given side from this Sidehandle
      void remove_sub_side_pointer(const Side* side) override
      {
        std::vector<Side*>::iterator tmpssit = subsides_.end();
        for (std::vector<Side*>::iterator ssit = subsides_.begin(); ssit != subsides_.end(); ++ssit)
        {
          if ((*ssit) == side)
          {
            tmpssit = ssit;
            break;
          }
        }
        if (tmpssit != subsides_.end())
          subsides_.erase(tmpssit);
        else
          FOUR_C_THROW("remove_sub_side_pointer: Couldn't identify subside!");
      }

      /// Add the SubSidePointer of given side to this Sidehandle
      void add_sub_side_pointer(Side* side) override
      {
        std::vector<Side*>::iterator tmpssit = subsides_.end();
        for (std::vector<Side*>::iterator ssit = subsides_.begin(); ssit != subsides_.end(); ++ssit)
        {
          if ((*ssit) == side)
          {
            tmpssit = ssit;
            break;
          }
        }
        if (tmpssit == subsides_.end())
          subsides_.push_back(side);
        else
          return;
      }

      /// Add the SubSide in to the unphysical list
      void mark_sub_sideunphysical(Side* side) override
      {
        std::vector<Side*>::iterator tmpssit = subsides_.end();
        for (std::vector<Side*>::iterator ssit = subsides_.begin(); ssit != subsides_.end(); ++ssit)
        {
          if ((*ssit) == side)
          {
            tmpssit = ssit;
            break;
          }
        }
        if (tmpssit == subsides_.end())
          FOUR_C_THROW(
              "mark_sub_sideunphysical failed, your side is not a Subside of the "
              "QuadraticSideHandle!");
        else
        {
          unphysical_subsides_.push_back(side);
        }
        return;
      }

      /// Is this side and unphysical subside
      bool isunphysical_sub_side(Side* side) override
      {
        for (std::vector<Side*>::iterator ssit = unphysical_subsides_.begin();
             ssit != unphysical_subsides_.end(); ++ssit)
        {
          if ((*ssit) == side) return true;
        }
        return false;
      }

      /// Does this sidehandle have unphysical subsides
      bool hasunphysical_sub_side() override { return unphysical_subsides_.size(); }

      /// Get the nodes of the Sidehandle
      const std::vector<Node*>& get_nodes() const override { return nodes_; }

     protected:
      std::vector<Side*> subsides_;
      std::vector<Node*> nodes_;

      /// all unphysical subsides of the side handle
      std::vector<Side*> unphysical_subsides_;
    };

    /// tri6 side handle
    class Tri6SideHandle : public QuadraticSideHandle
    {
     public:
      Tri6SideHandle(Mesh& mesh, int sid, const std::vector<int>& node_ids);

      Core::FE::CellType shape() override { return Core::FE::CellType::tri6; }

      void local_coordinates(
          const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<2, 1>& rst) override;
    };

    /// quad4 side handle
    /*!
     * quad4 need to be split in 3 tri3 in order to avoid subtle ambiguities.
     */
    class Quad4SideHandle : public QuadraticSideHandle
    {
     public:
      Quad4SideHandle(Mesh& mesh, int sid, const std::vector<int>& node_ids);

      Core::FE::CellType shape() override { return Core::FE::CellType::quad4; }

      void local_coordinates(
          const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<2, 1>& rst) override;
    };

    /// quad8 side handle
    class Quad8SideHandle : public QuadraticSideHandle
    {
     public:
      Quad8SideHandle(Mesh& mesh, int sid, const std::vector<int>& node_ids, bool iscutside = true);

      Core::FE::CellType shape() override { return Core::FE::CellType::quad8; }

      void local_coordinates(
          const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<2, 1>& rst) override;
    };

    /// quad9 side handle
    class Quad9SideHandle : public QuadraticSideHandle
    {
     public:
      Quad9SideHandle(Mesh& mesh, int sid, const std::vector<int>& node_ids, bool iscutside = true);

      Core::FE::CellType shape() override { return Core::FE::CellType::quad9; }

      void local_coordinates(
          const Core::LinAlg::Matrix<3, 1>& xyz, Core::LinAlg::Matrix<2, 1>& rst) override;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
