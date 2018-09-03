/*!-----------------------------------------------------------------------------------------------*
 \file cut_node.H

 \brief class representing a geometrical node

 <pre>
\level 3
\maintainer Andy Wirtz
 wirtz@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15270
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_NODE_H
#define CUT_NODE_H

#include "cut_point.H"

namespace GEO
{
  namespace CUT
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
      bool operator()(const plain_volumecell_set& s1, const plain_volumecell_set& s2);

      /// compare operator for sets of plain_volumecell_sets
      bool operator()(const std::set<plain_volumecell_set, Cmp>& set1,
          const std::set<plain_volumecell_set, Cmp>& set2);

      /// compare routine for sets of plain_volumecell_set
      bool Compare(const std::set<plain_volumecell_set, Cmp>& set1,
          const std::set<plain_volumecell_set, Cmp>& set2);

     private:
      /// compare routine for plain_volumecell_sets
      bool Compare(const plain_volumecell_set& s1, const plain_volumecell_set& s2);

      /// compare routine for volumecells based on point ids
      bool Compare(VolumeCell* vc1, VolumeCell* vc2);
    };

    /*----------------------------------------------------------------------------*/
    class NodalDofSet
    {
     public:
      /// constructor for a nodal dofset
      explicit NodalDofSet(
          std::set<plain_volumecell_set, Cmp>& connected_volumecells, bool is_std_dofset);

      /// constructor for CompositeNodalDofSet
      explicit NodalDofSet(bool is_std_dofset, GEO::CUT::Point::PointPosition pos);

      virtual ~NodalDofSet(){};

      /// return the composite of all connected volumecells between elements and the element's
      /// sub-elements
      const std::set<plain_volumecell_set, Cmp>& VolumeCellComposite() const
      {
        return volumecell_composite_;
      }

      /// collect cut sides of the entire volumecell composite
      void CollectCutSides(plain_int_set& cutside_ids) const;

      /// is the nodal dofset a standard dofset?
      bool Is_Standard_DofSet() const { return is_std_dofset_; }

      /// get the position of the nodal dofset which is the same for all its contained volumecells
      GEO::CUT::Point::PointPosition Position() const { return position_; };

      /// does the nodal dofset's composite of volumecells contain the point?
      virtual bool Contains(GEO::CUT::Point* p);

      virtual void Print();

     protected:
      /** stores for each element (and its subelements) in a 1-ring around the node
       *  all its volumecells which has been connected between the surrounding elements
       *  for linear elements the plain_volumecell_sets contain only one volumecell,
       *  for quadratic elements usually more vc's */
      std::set<plain_volumecell_set, Cmp> volumecell_composite_;

      /// is the dofset a standard physical dofset or a ghost dofset?
      bool is_std_dofset_;

      GEO::CUT::Point::PointPosition position_;
    };

    /*----------------------------------------------------------------------------*/
    class CompositeNodalDofSet : public NodalDofSet
    {
     public:
      /// constructor
      CompositeNodalDofSet(bool is_std_dofset, GEO::CUT::Point::PointPosition pos)
          : NodalDofSet(is_std_dofset, pos)
      {
        nodal_dofsets_.clear();
      };

      /// add a nodal dofset to the composite of nodal dofsets
      void add(Teuchos::RCP<NodalDofSet> nds, bool allow_connect_std_and_ghost_sets)
      {
        if (allow_connect_std_and_ghost_sets)
        {
          if (position_ != nds->Position())
            dserror(
                "NodalDofSet you want to combine to a CompositeNodalDofSet do not have the same "
                "position! Invalid!");
        }
        else  // require same type of dofsets: std/ghost
        {
          if (is_std_dofset_ != nds->Is_Standard_DofSet() or position_ != nds->Position())
            dserror(
                "NodalDofSet you want to combine to a CompositeNodalDofSet do not have the same "
                "properties! Invalid!");
        }

        nodal_dofsets_.push_back(nds);

        std::copy(nds->VolumeCellComposite().begin(), nds->VolumeCellComposite().end(),
            std::inserter(volumecell_composite_, volumecell_composite_.end()));
      }

      /// does one of the nodal dofset's contain the point?
      bool Contains(GEO::CUT::Point* p)
      {
        for (unsigned int i = 0; i < nodal_dofsets_.size(); i++)
          if (nodal_dofsets_[i]->Contains(p)) return true;

        return false;
      }

      void Print();

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
      int Id() const { return nid_; }

      /// register an edge adjacent to this node
      void Register(Edge* edge)
      {
        edges_.insert(edge);
        point_->AddEdge(edge);
      }

      /// register a side adjacent to this node
      void Register(Side* side) { point_->AddSide(side); }

      /// register an element adjacent to this node
      void Register(Element* element) { point_->AddElement(element); }

      /// register cuts
      void RegisterCuts();


      /*========================================================================*/
      //! @name get routines for adjacent objects
      /*========================================================================*/

      /// get edges adjacent to this node
      const plain_edge_set& Edges() { return edges_; }

      /// Get the coordinates of the node from its point's information
      void Coordinates(double* x) const { point_->Coordinates(x); }

      /// Get the position of the node whether it is in fluid, structure or on the cut face
      Point::PointPosition Position() const { return point_->Position(); }

      /// Returns the point that defines the node
      Point* point() const { return point_; }

      /// Returns the levelset value at this node (if it is a levelset node)
      double LSV() const { return lsv_; }

      /// Returns sides that are connected at this node
      const plain_side_set& Sides() const { return point_->CutSides(); }

      /// get all cut elements adjacent to this node's point
      const plain_element_set& Elements() const { return point_->Elements(); }


      /*========================================================================*/
      //! @name print and plot routines
      /*========================================================================*/

      /// print node's or its point's information to the stream
      void Print(std::ostream& f = std::cout) { point_->Print(); }

      /// plot node's or its point's information to the stream
      void Plot(std::ostream& f) { point_->Plot(f, this->Id()); };


      /*========================================================================*/
      //! @name DofSet-management
      /*========================================================================*/

      /// Assign the vc_sets to the node if possible
      void AssignNodalCellSet(const std::vector<plain_volumecell_set>& ele_vc_sets,
          std::map<Node*, std::vector<plain_volumecell_set>>& nodal_cell_sets);

      /// Find the dofsets required at this node. (old unused version)
      void FindDOFSets(bool include_inner);

      /// Find the dofsets required at this node.
      void FindDOFSetsNEW(std::map<Node*, std::vector<plain_volumecell_set>>& nodal_cell_sets,
          std::vector<plain_volumecell_set>& cell_sets);

      /// get the dofset number of the Volumecell w.r.t this node (old unused version)
      int DofSetNumber(VolumeCell* cell);

      /// get the dofset number of the Volumecell w.r.t this node
      int DofSetNumberNEW(const plain_volumecell_set& cells);

      /// return the sets of volumecells (old unsed version)
      const std::vector<plain_volumecell_set>& DofSets() const { return dofsets_; }


      /*========================================================================*/
      //! @name Selfcut related routines
      /*========================================================================*/

      /// Returns the selfcutposition of this node
      Point::PointPosition SelfCutPosition() { return selfcutposition_; }

      /// Gives this node a selfcutposition and spreads the positional information
      void SelfCutPosition(Point::PointPosition p);

      /// Changes the selfcutposition of this node and spreads the positional information
      void ChangeSelfCutPosition(Point::PointPosition p);

      /// Erase the cutsideedge from this node because it is deleted in the selfcut
      void EraseCutSideEdge(Edge* cutsideedge) { edges_.erase(cutsideedge); }

      /// Return true is the given node has same position coordinates as this node
      bool isAtSameLocation(const Node* nod1) const;


      /*========================================================================*/
      //! @name nodal DofSets
      /*========================================================================*/

      /// get the number of dofsets at this node
      int NumDofSets() const;

      /// return a vector of all nodal dofsets
      const std::vector<Teuchos::RCP<NodalDofSet>>& NodalDofSets() const { return nodaldofsets_; }

      /// return a vector of all nodal dofsets
      NodalDofSet* GetNodalDofSet(const int nds) const { return &*nodaldofsets_[nds]; }

      /// remove non-standard nodal dofsets
      void RemoveNonStandardNodalDofSets();

      /// get the unique standard NodalDofSet for a given nodal dofset position
      int GetStandardNodalDofSet(Point::PointPosition pos);


      /*========================================================================*/
      //! @name sorting and collecting of nodal DofSets
      /*========================================================================*/

      /// sort all nodal dofsets via xyz point coordinates (use compare functions in cut_node.H)
      void SortNodalDofSets();

      /// collect the (ghost) dofsets for this node w.r.t each phase to avoid multiple ghost nodal
      /// dofsets for a certain phase
      void CollectNodalDofSets(bool connect_ghost_with_standard_nds);

      /// compare operator for sorting NodalDofSets
      bool operator()(NodalDofSet* nodaldofset1, NodalDofSet* nodaldofset2);



     private:
      /// build sets of connected volumecells in a 1-ring around the node
      void BuildDOFCellSets(Point* p, const std::vector<plain_volumecell_set>& cell_sets,
          const plain_volumecell_set& cells, const std::vector<plain_volumecell_set>& nodal_cells,
          plain_volumecell_set& done, bool isnodalcellset = false);

      /// build sets of connected volumecells in a 1-ring around the node (old unused version)
      void BuildDOFCellSets(Point* p, const plain_volumecell_set& cells,
          const plain_volumecell_set& nodal_cells, plain_volumecell_set& done);


      /*========================================================================*/
      //! @name private class variables
      /*========================================================================*/

      int nid_;       ///< node's id (the same as in DRT::Discretization)
      Point* point_;  ///< pointer to the GEO::CUT::Point which defines the node
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


  }  // namespace CUT
}  // namespace GEO

#endif
