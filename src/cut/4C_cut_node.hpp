/*----------------------------------------------------------------------*/
/*! \file

\brief class representing a geometrical node


\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_NODE_HPP
#define FOUR_C_CUT_NODE_HPP

#include "4C_config.hpp"

#include "4C_cut_point.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class VolumeCell;

    /*----------------------------------------------------------------------------*/
    /*!
      \brief internal comparator class used for sorting sets of plain_volumecell_sets
             and plain_volumecell_sets. This sorting is necessary to get the same order
             of sets, important for restarts, restarts of Newton schemes in monolithic
             solves but also for debugging
     */
    class Cmp
    {
     public:
      /// compare operator for plain_volumecell_sets
      bool operator()(const plain_volumecell_set& s1, const plain_volumecell_set& s2) const;

      /// compare operator for sets of plain_volumecell_sets
      bool operator()(const std::set<plain_volumecell_set, Cmp>& set1,
          const std::set<plain_volumecell_set, Cmp>& set2) const;

      /// compare routine for sets of plain_volumecell_set
      bool compare(const std::set<plain_volumecell_set, Cmp>& set1,
          const std::set<plain_volumecell_set, Cmp>& set2) const;

     private:
      /// compare routine for plain_volumecell_sets
      bool compare(const plain_volumecell_set& s1, const plain_volumecell_set& s2) const;

      /// compare routine for volumecells based on point ids
      bool compare(VolumeCell* vc1, VolumeCell* vc2) const;
    };

    /*----------------------------------------------------------------------------*/
    class NodalDofSet
    {
     public:
      /// constructor for a nodal dofset
      explicit NodalDofSet(
          std::set<plain_volumecell_set, Cmp>& connected_volumecells, bool is_std_dofset);

      /// constructor for CompositeNodalDofSet
      explicit NodalDofSet(bool is_std_dofset, Core::Geo::Cut::Point::PointPosition pos);

      virtual ~NodalDofSet() = default;

      /// return the composite of all connected volumecells between elements and the element's
      /// sub-elements
      const std::set<plain_volumecell_set, Cmp>& volume_cell_composite() const
      {
        return volumecell_composite_;
      }

      /// collect cut sides of the entire volumecell composite
      void collect_cut_sides(plain_int_set& cutside_ids) const;

      /// is the nodal dofset a standard dofset?
      bool is_standard_dof_set() const { return is_std_dofset_; }

      /// get the position of the nodal dofset which is the same for all its contained volumecells
      Core::Geo::Cut::Point::PointPosition position() const { return position_; };

      /// does the nodal dofset's composite of volumecells contain the point?
      virtual bool contains(Core::Geo::Cut::Point* p);

      virtual void print();

     protected:
      /** stores for each element (and its subelements) in a 1-ring around the node
       *  all its volumecells which has been connected between the surrounding elements
       *  for linear elements the plain_volumecell_sets contain only one volumecell,
       *  for quadratic elements usually more vc's */
      std::set<plain_volumecell_set, Cmp> volumecell_composite_;

      /// is the dofset a standard physical dofset or a ghost dofset?
      bool is_std_dofset_;

      Core::Geo::Cut::Point::PointPosition position_;
    };

    /*----------------------------------------------------------------------------*/
    class CompositeNodalDofSet : public NodalDofSet
    {
     public:
      /// constructor
      CompositeNodalDofSet(bool is_std_dofset, Core::Geo::Cut::Point::PointPosition pos)
          : NodalDofSet(is_std_dofset, pos)
      {
        nodal_dofsets_.clear();
      };

      /// add a nodal dofset to the composite of nodal dofsets
      void add(Teuchos::RCP<NodalDofSet> nds, bool allow_connect_std_and_ghost_sets)
      {
        if (allow_connect_std_and_ghost_sets)
        {
          if (position_ != nds->position())
            FOUR_C_THROW(
                "NodalDofSet you want to combine to a CompositeNodalDofSet do not have the same "
                "position! Invalid!");
        }
        else  // require same type of dofsets: std/ghost
        {
          if (is_std_dofset_ != nds->is_standard_dof_set() or position_ != nds->position())
            FOUR_C_THROW(
                "NodalDofSet you want to combine to a CompositeNodalDofSet do not have the same "
                "properties! Invalid!");
        }

        nodal_dofsets_.push_back(nds);

        std::copy(nds->volume_cell_composite().begin(), nds->volume_cell_composite().end(),
            std::inserter(volumecell_composite_, volumecell_composite_.end()));
      }

      /// does one of the nodal dofset's contain the point?
      bool contains(Core::Geo::Cut::Point* p) override
      {
        for (unsigned int i = 0; i < nodal_dofsets_.size(); i++)
          if (nodal_dofsets_[i]->contains(p)) return true;

        return false;
      }

      void print() override;

     private:
      std::vector<Teuchos::RCP<NodalDofSet>>
          nodal_dofsets_;  ///< a collection of NodalDofSet objects combined to one composite
    };

    /*----------------------------------------------------------------------------*/
    /// One node in a cut mesh
    /*!
     \brief Class to deal with nodes in the cut mesh. Nodes have a mesh-unique id, a
            levelset-value and a point. The point contains the major part of the
            information.

      Furthermore, nodes have some idea how many dofsets are required. This is
      optional information that can be used by the xfem code.
     */
    class Node
    {
     public:
      /// constructor
      Node(int nid, Point* point, double lsv)
          : nid_(nid), point_(point), lsv_(lsv), selfcutposition_(Point::undecided)
      {
      }

      /// get node's node Id
      int id() const { return nid_; }

      /// register an edge adjacent to this node
      void register_entity(Edge* edge)
      {
        edges_.insert(edge);
        point_->add_edge(edge);
      }

      /// register a side adjacent to this node
      void register_entity(Side* side) { point_->add_side(side); }

      /// register an element adjacent to this node
      void register_entity(Element* element) { point_->add_element(element); }

      /// register cuts
      void register_cuts();


      /*========================================================================*/
      //! @name get routines for adjacent objects
      /*========================================================================*/

      /// get edges adjacent to this node
      const plain_edge_set& edges() { return edges_; }

      /// Get the coordinates of the node from its point's information
      void coordinates(double* x) const { point_->coordinates(x); }

      /// Get the position of the node whether it is in fluid, structure or on the cut face
      Point::PointPosition position() const { return point_->position(); }

      /// Returns the point that defines the node
      Point* point() const { return point_; }

      /// Returns the levelset value at this node (if it is a levelset node)
      double lsv() const { return lsv_; }

      /// Returns sides that are connected at this node
      const plain_side_set& sides() const { return point_->cut_sides(); }

      /// get all cut elements adjacent to this node's point
      const plain_element_set& elements() const { return point_->elements(); }


      /*========================================================================*/
      //! @name print and plot routines
      /*========================================================================*/

      /// print node's or its point's information to the stream
      void print(std::ostream& f = std::cout) { point_->print(); }

      /// plot node's or its point's information to the stream
      void plot(std::ostream& f) { point_->plot(f, this->id()); };


      /*========================================================================*/
      //! @name DofSet-management
      /*========================================================================*/

      /// Assign the vc_sets to the node if possible
      void assign_nodal_cell_set(const std::vector<plain_volumecell_set>& ele_vc_sets,
          std::map<Node*, std::vector<plain_volumecell_set>>& nodal_cell_sets);

      /// Find the dofsets required at this node. (old unused version)
      void find_dof_sets(bool include_inner);

      /// Find the dofsets required at this node.
      void find_dof_sets_new(std::map<Node*, std::vector<plain_volumecell_set>>& nodal_cell_sets,
          std::vector<plain_volumecell_set>& cell_sets);

      /// get the dofset number of the Volumecell w.r.t this node (old unused version)
      int dof_set_number(VolumeCell* cell);

      /// get the dofset number of the Volumecell w.r.t this node
      int dof_set_number_new(const plain_volumecell_set& cells);

      /// return the sets of volumecells (old unsed version)
      const std::vector<plain_volumecell_set>& dof_sets() const { return dofsets_; }


      /*========================================================================*/
      //! @name Selfcut related routines
      /*========================================================================*/

      /// Returns the selfcutposition of this node
      Point::PointPosition self_cut_position() { return selfcutposition_; }

      /// Gives this node a selfcutposition and spreads the positional information
      void self_cut_position(Point::PointPosition p);

      /// Changes the selfcutposition of this node and spreads the positional information
      void change_self_cut_position(Point::PointPosition p);

      /// Erase the cutsideedge from this node because it is deleted in the selfcut
      void erase_cut_side_edge(Edge* cutsideedge) { edges_.erase(cutsideedge); }

      /// Return true is the given node has same position coordinates as this node
      bool is_at_same_location(const Node* nod1) const;


      /*========================================================================*/
      //! @name nodal DofSets
      /*========================================================================*/

      /// get the number of dofsets at this node
      int num_dof_sets() const;

      /// return a vector of all nodal dofsets
      const std::vector<Teuchos::RCP<NodalDofSet>>& nodal_dof_sets() const { return nodaldofsets_; }

      /// return a vector of all nodal dofsets
      NodalDofSet* get_nodal_dof_set(const int nds) const { return &*nodaldofsets_[nds]; }

      /// remove non-standard nodal dofsets
      void remove_non_standard_nodal_dof_sets();

      /// get the unique standard NodalDofSet for a given nodal dofset position
      int get_standard_nodal_dof_set(Point::PointPosition pos);


      /*========================================================================*/
      //! @name sorting and collecting of nodal DofSets
      /*========================================================================*/

      /// sort all nodal dofsets via xyz point coordinates (use compare functions in cut_node.H)
      void sort_nodal_dof_sets();

      /// collect the (ghost) dofsets for this node w.r.t each phase to avoid multiple ghost nodal
      /// dofsets for a certain phase
      void collect_nodal_dof_sets(bool connect_ghost_with_standard_nds);

      /// compare operator for sorting NodalDofSets
      bool operator()(NodalDofSet* nodaldofset1, NodalDofSet* nodaldofset2) const;



     private:
      /// build sets of connected volumecells in a 1-ring around the node
      void build_dof_cell_sets(Point* p, const std::vector<plain_volumecell_set>& cell_sets,
          const plain_volumecell_set& cells,
          const std::vector<plain_volumecell_set>& nodal_cell_sets, plain_volumecell_set& done,
          bool isnodalcellset = false);

      /// build sets of connected volumecells in a 1-ring around the node (old unused version)
      void build_dof_cell_sets(Point* p, const plain_volumecell_set& cells,
          const plain_volumecell_set& nodal_cells, plain_volumecell_set& done);


      /*========================================================================*/
      //! @name private class variables
      /*========================================================================*/

      int nid_;       ///< node's id (the same as in  Core::FE::Discretization)
      Point* point_;  ///< pointer to the Core::Geo::Cut::Point which defines the node
      double lsv_;    ///< levelset value if it is a levelset cut

      plain_edge_set edges_;  ///< all adjacent edges

      std::vector<plain_volumecell_set> dofsets_;  ///< set of Volumecells (unused old version)

      std::vector<Teuchos::RCP<NodalDofSet>> nodaldofsets_;  ///< nodal dofsets

      Point::PointPosition selfcutposition_;  ///< every cutsidenode knows its selfcutposition
    };

    /*----------------------------------------------------------------------------*/
    /*!
      \brief internal comparator class used for sorting nodal dofsets.
             This sorting is necessary to get the
             same order of sets, important for restarts, restarts of Newton schemes in
             monolithic solves but also for debugging
     */
    class NodalDofSetCmp
    {
     public:
      bool operator()(
          Teuchos::RCP<NodalDofSet> nodaldofset1, Teuchos::RCP<NodalDofSet> nodaldofset2);
    };

    /** \brief Find if the nodes in nelement share a common element,
     *  (i.e. do all nodes lie on the same element?)
     *
     *  A node knows by the help of its underlying point which cut_elements
     *  it is associated to.
     *
     *  \author hiermeier \date 11/16 */
    void FindCommonElements(const std::vector<Node*>& nelement, plain_element_set& elements);


  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
