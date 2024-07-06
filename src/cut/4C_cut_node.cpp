/*----------------------------------------------------------------------*/
/*! \file

\brief class representing a geometrical node


\level 3

 *------------------------------------------------------------------------------------------------*/

#include "4C_cut_edge.hpp"
#include "4C_cut_element.hpp"
#include "4C_cut_volumecell.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*
Functions to catch implementation erros in debug mode
*/
namespace
{
  [[maybe_unused]] bool IsCutPositionUnchanged(
      Core::Geo::Cut::Point::PointPosition cutposition, Core::Geo::Cut::Point::PointPosition pos)
  {
    if ((cutposition == Core::Geo::Cut::Point::inside and pos == Core::Geo::Cut::Point::outside) or
        (cutposition == Core::Geo::Cut::Point::outside and pos == Core::Geo::Cut::Point::inside))
      return false;
    else
      return true;
  }
}  // namespace


/*------------------------------------------------------------------------------*
  | Operator () compare operator for plain_volumecell_sets
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Cmp::operator()(
    const plain_volumecell_set& s1, const plain_volumecell_set& s2) const
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets ---
  // Core::Geo::Cut::Cmp::operator()" );

  // call Compare function for two plain_volumecell_sets
  return compare(s1, s2);
}



/*------------------------------------------------------------------------------*
  | Compare() to compare two sets of plain_volumecell_set via comparing their first
 plain_volumecell_sets |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Cmp::compare(const std::set<plain_volumecell_set, Cmp>& set1,
    const std::set<plain_volumecell_set, Cmp>& set2) const
{
  // compare two sets of plain_volumecell_set
  // take the first plain_volumecell_set of each set and compare them

  std::set<plain_volumecell_set, Cmp>::iterator it1 = set1.begin();
  const plain_volumecell_set& s1 = *it1;

  std::set<plain_volumecell_set, Cmp>::iterator it2 = set2.begin();
  const plain_volumecell_set& s2 = *it2;

  return compare(s1, s2);
}

/*------------------------------------------------------------------------------*
  | Compare() to compare two plain_volumecell_set via the ids of their first volumecell's points
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Cmp::compare(
    const plain_volumecell_set& s1, const plain_volumecell_set& s2) const
{
  // take the first vc in plain_volumecell_set. In case of linear elements
  // this is the only volumecell. For quadratic elements this is a set of vcs
  // but we just compare the first vc of each plain_volumecell_set. So for
  // quadratic elements the  plain_volumecell_set is not sorted!

  plain_volumecell_set::const_iterator it1 = s1.begin();
  VolumeCell* vc1 = *it1;

  plain_volumecell_set::const_iterator it2 = s2.begin();
  VolumeCell* vc2 = *it2;

  return compare(vc1, vc2);
}


/*------------------------------------------------------------------------------*
  | Operator () to compare two volume cells via the ids of their points
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool Core::Geo::Cut::Cmp::compare(VolumeCell* vc1, VolumeCell* vc2) const
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets ---
  // Core::Geo::Cut::Cmp::Compare(vc,vc)" );


  if (vc1 == vc2) return false;

  // compare the point ids of the two volumecells

  const std::set<int>& vc1points = vc1->volume_cell_point_ids();
  const std::set<int>& vc2points = vc2->volume_cell_point_ids();


  // during reducing the two sets to two minimized/disjoint sets which don't have any points in
  // common the first non-common point ids can be used to sort the sets

  std::set<int>::iterator it1 = vc1points.begin();
  std::set<int>::iterator it2 = vc2points.begin();

  while (true)  // sets are already sorted!
  {
    if (it1 == vc1points.end() or it2 == vc2points.end()) break;

    if (*it1 == *it2)
    {
      // identical points, go to the next one
      it1++;
      it2++;
    }
    else  // compare the points now via their point ids
    {
      if (*it1 < *it2)
        return true;
      else
        return false;  // (*it1 > *it2)
    }
  }

  FOUR_C_THROW(
      "sorting failed: one volume-cell is completely contained in the other volume-cell or both "
      "vcs are equal!");

  return false;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
Core::Geo::Cut::NodalDofSet::NodalDofSet(
    std::set<plain_volumecell_set, Cmp>& connected_volumecells, bool is_std_dofset)
    : is_std_dofset_(is_std_dofset)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets ---
  // Core::Geo::Cut::NodalDofSet::NodalDofSet" );

  std::copy(connected_volumecells.begin(), connected_volumecells.end(),
      std::inserter(volumecell_composite_, volumecell_composite_.end()));

  // set the position of the nodal dofset
  const plain_volumecell_set& set = *(volumecell_composite_.begin());
  position_ = (*set.begin())->position();
}


Core::Geo::Cut::NodalDofSet::NodalDofSet(
    bool is_std_dofset, Core::Geo::Cut::Point::PointPosition pos)
    : is_std_dofset_(is_std_dofset), position_(pos)
{
  volumecell_composite_.clear();
}

/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
bool Core::Geo::Cut::NodalDofSet::contains(Core::Geo::Cut::Point* p)
{
  // check if any volume-cell of the volumecell_composite contains this point
  for (std::set<plain_volumecell_set, Cmp>::iterator it = volumecell_composite_.begin();
       it != volumecell_composite_.end(); it++)
  {
    const plain_volumecell_set& vc_set = *it;
    for (plain_volumecell_set::const_iterator vcs = vc_set.begin(); vcs != vc_set.end(); vcs++)
    {
      if ((*vcs)->contains(p)) return true;
    }
  }

  return false;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::NodalDofSet::collect_cut_sides(
    Core::Geo::Cut::plain_int_set& cutside_ids) const
{
  // collect all cut sides
  for (std::set<plain_volumecell_set, Cmp>::iterator it = volumecell_composite_.begin();
       it != volumecell_composite_.end(); it++)
  {
    const plain_volumecell_set& vc_set = *it;
    for (plain_volumecell_set::const_iterator vcs = vc_set.begin(); vcs != vc_set.end(); vcs++)
    {
      (*vcs)->collect_cut_sides(cutside_ids);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::NodalDofSet::print()
{
  std::cout << "Core::Geo::Cut::NodalDofSet:\n"
            << "STD dofset = " << (this->is_standard_dof_set() ? "TRUE\n" : "FALSE\n")
            << "Position   = " << Point::point_position2_string(this->position()) << std::endl;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::CompositeNodalDofSet::print()
{
  std::cout << "Core::Geo::Cut::CompositeNodalDofSet which contains " << nodal_dofsets_.size()
            << " combined Core::Geo::Cut::NodalDofSet:\n "
            << "STD dofset = " << (this->is_standard_dof_set() ? "TRUE\n" : "FALSE\n")
            << "Position   = " << Point::point_position2_string(this->position()) << std::endl;
}


/*-----------------------------------------------------------------------------------------*
 * register cuts
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::register_cuts()
{
  if (position() == Point::oncutsurface)
  {
    for (plain_edge_set::iterator i = edges_.begin(); i != edges_.end(); ++i)
    {
      Edge* e = *i;
      point_->add_edge(e);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * Assign the vc_sets to the node if possible
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::assign_nodal_cell_set(
    const std::vector<plain_volumecell_set>& ele_vc_sets,
    std::map<Node*, std::vector<plain_volumecell_set>>& nodal_cell_sets)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets --- AssignNodalCellSet");

  std::vector<plain_volumecell_set>& nodal_cell_set = nodal_cell_sets[this];

  for (std::vector<plain_volumecell_set>::const_iterator s = ele_vc_sets.begin();
       s != ele_vc_sets.end(); s++)
  {
    const plain_volumecell_set& cell_set = *s;

    for (plain_volumecell_set::const_iterator i = cell_set.begin(); i != cell_set.end(); ++i)
    {
      VolumeCell* cell = *i;

      // if at least one cell of this cell_set contains the point, then the whole cell_set contains
      // the point
      bool contains = false;

      {
        TEUCHOS_FUNC_TIME_MONITOR(
            "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets --- cell->Contains( point() )");

        contains = cell->contains(point());
      }

      if (contains)
      {
        nodal_cell_set.push_back(cell_set);

        // the rest of cells in this set has not to be checked for this node
        break;  // finish the cell_set loop, breaks the inner for loop!
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * Find the dofsets required at this node. (old unused version)
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::find_dof_sets(bool include_inner)
{
  const plain_element_set& elements = Node::elements();

  std::map<Node*, plain_volumecell_set> nodal_cells;
  plain_volumecell_set cells;

  for (plain_element_set::const_iterator i = elements.begin(); i != elements.end(); ++i)
  {
    Element* e = *i;
    {
      const plain_volumecell_set& element_cells = e->volume_cells();
      if (include_inner)
      {
        std::copy(element_cells.begin(), element_cells.end(), std::inserter(cells, cells.begin()));
      }
      else
      {
        for (plain_volumecell_set::const_iterator i = element_cells.begin();
             i != element_cells.end(); ++i)
        {
          VolumeCell* vc = *i;
          if (vc->position() == Point::outside)
          {
            cells.insert(vc);
          }
        }
      }

      const std::vector<Node*>& nodes = e->nodes();
      for (std::vector<Node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i)
      {
        Node* n = *i;
        for (plain_volumecell_set::const_iterator i = element_cells.begin();
             i != element_cells.end(); ++i)
        {
          VolumeCell* cell = *i;

          if (not include_inner and cell->position() != Point::outside) continue;

          if (cell->contains(n->point()))
          {
            nodal_cells[n].insert(cell);
          }
        }
      }
    }
  }

  //   std::cout << "id=" << Id()
  //             << " #ele=" << elements.size()
  //             << " #cells=" << cells.size()
  //             << " #nodal=" << nodal_cells.size()
  //     ;

  // First, get the nodal cells that make up the first dofset. In most cases
  // this loop has one pass only. But if the node is cut, there will be more
  // than one set of cells that are attached to this node.

  plain_volumecell_set done;

  build_dof_cell_sets(point(), cells, nodal_cells[this], done);

  nodal_cells.erase(this);

  for (std::map<Node*, plain_volumecell_set>::iterator i = nodal_cells.begin();
       i != nodal_cells.end(); ++i)
  {
    Node* n = i->first;
    plain_volumecell_set& cellset = i->second;
    build_dof_cell_sets(n->point(), cells, cellset, done);
  }

  // do any remaining internal volumes that are not connected to any node
  build_dof_cell_sets(nullptr, cells, cells, done);

  //   std::cout << " #dofsets: " << dofsets_.size() << " [";
  //   for ( unsigned i=0; i<dofsets_.size(); ++i )
  //   {
  //     std::cout << " " << dofsets_[i].size();
  //   }
  //   std::cout << " ]\n";
}


/*-----------------------------------------------------------------------------------------*
 * Find the dofsets required at this node.
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::find_dof_sets_new(
    std::map<Node*, std::vector<plain_volumecell_set>>& nodal_cell_sets,
    std::vector<plain_volumecell_set>& cell_sets)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets ---
  // FindDOFSetsNEW" );


  // finally: fill dof_cellsets_

  // do the connection between elements
  plain_volumecell_set done;
  plain_volumecell_set cells;

  // get cell_sets as a plain_volume_set
  for (std::vector<plain_volumecell_set>::iterator i = cell_sets.begin(); i != cell_sets.end(); i++)
  {
    std::copy((*i).begin(), (*i).end(), std::inserter(cells, cells.end()));
  }


  // First, get the nodal cells that make up the first dofset. In most cases
  // this loop has one pass only. But if the node is cut, there will be more
  // than one set of cells that are attached to this node.

  // call this function with isnodalcellset=true flag to identify the first std set
  build_dof_cell_sets(point(), cell_sets, cells, nodal_cell_sets[this], done, true);

  nodal_cell_sets.erase(this);


  for (std::map<Node*, std::vector<plain_volumecell_set>>::iterator i = nodal_cell_sets.begin();
       i != nodal_cell_sets.end(); ++i)
  {
    Node* n = i->first;

    std::vector<plain_volumecell_set>& cellset = i->second;
    build_dof_cell_sets(n->point(), cell_sets, cells, cellset, done);
  }

  // do any remaining internal volumes that are not connected to any node
  build_dof_cell_sets(nullptr, cell_sets, cells, cell_sets, done);
}


/*-----------------------------------------------------------------------------------------*
 * get the dofset number of the Volumecell w.r.t this node (old unused version)
 *-----------------------------------------------------------------------------------------*/
int Core::Geo::Cut::Node::dof_set_number(VolumeCell* cell)
{
  int dofset = -1;
  for (unsigned i = 0; i < dofsets_.size(); ++i)
  {
    plain_volumecell_set& cells = dofsets_[i];
    if (cells.count(cell) > 0)
    {
      if (dofset == -1)
      {
        dofset = i;
      }
      else
      {
        FOUR_C_THROW("volume dofset not unique");
      }
    }
  }
  if (dofset == -1)
  {
    std::cout << "dofset not found for node " << this->id() << std::endl;
    FOUR_C_THROW("volume dofset not found");
  }
  return dofset;
}


/*-----------------------------------------------------------------------------------------*
 * get number of dofsets
 *-----------------------------------------------------------------------------------------*/
int Core::Geo::Cut::Node::num_dof_sets() const { return nodal_dof_sets().size(); }



/*-----------------------------------------------------------------------------------------*
 * get the dofset number of the Volumecell w.r.t this node
 *-----------------------------------------------------------------------------------------*/
int Core::Geo::Cut::Node::dof_set_number_new(const plain_volumecell_set& cells)
{
  int dofset = -1;

  // find the first cell of cells, this is only a volume cell of a subelement
  if (cells.size() == 0) FOUR_C_THROW("cells is empty");

  //  VolumeCell* cell = cells[0];
  VolumeCell* cell = *(cells.begin());

  for (unsigned int i = 0; i < nodaldofsets_.size(); ++i)  // loop over sets
  {
    const std::set<plain_volumecell_set, Core::Geo::Cut::Cmp>& cellsets =
        nodaldofsets_[i]->volume_cell_composite();

    for (std::set<plain_volumecell_set, Core::Geo::Cut::Cmp>::const_iterator j = cellsets.begin();
         j != cellsets.end(); j++)
    {
      if (j->count(cell) > 0)
      {
        if (dofset == -1)
        {
          dofset = i;
        }
        else
        {
          std::cout << "node: " << id() << std::endl;
          std::cout << "first dofset id: " << dofset << std::endl;
          std::cout << "new dofset id: " << i << std::endl;
          cell->print(std::cout);
          FOUR_C_THROW("volume dofset not unique");
        }
      }
    }
  }
  if (dofset == -1)
  {
    std::cout << "dofset not found for node " << this->id() << std::endl;
    //    FOUR_C_THROW( "volume dofset not found" );
  }
  return dofset;
}



/*-----------------------------------------------------------------------------------------*
 * sort all dofsets via xyz point coordinates (use compare functions in cut_node.H)
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::sort_nodal_dof_sets()
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets ---
  // SortNodalDofSets" );


  if (nodaldofsets_.size() > 1)
  {
    std::vector<Teuchos::RCP<Core::Geo::Cut::NodalDofSet>>::iterator it_start =
        nodaldofsets_.begin();

    if (nodaldofsets_[0]->is_standard_dof_set())
    {
      // REMARK:
      // if the first nodal dofset is a standard dofset, then the first set must not be changed
      // during sorting because elements without eh assume the first set as std set!

      it_start++;  // exclude the standard set from sorting
    }

    // sort the cellsets w.r.t point ids in first vc in first set of sorted sets of plain volume
    // cells sets REMARK: do not sort the first dofset, it has to be kept the standard dofset
    sort(it_start, nodaldofsets_.end(), NodalDofSetCmp());
  }

  return;
}


/*-----------------------------------------------------------------------------------------*
 * collect the (ghost) dofsets for this node w.r.t each phase to avoid multiple ghost nodal dofsets
 *for a certain phase
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::collect_nodal_dof_sets(bool connect_ghost_with_standard_nds)
{
  // TEUCHOS_FUNC_TIME_MONITOR( "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets ---
  // CollectNodalDofSets"
  // );

  // assume that the nodal dofsets have been sorted in a step before,
  // such that ghost sets to be combined are stored consecutively in the vector of sorted nodal
  // dofsets

  std::vector<Teuchos::RCP<CompositeNodalDofSet>> collected_nodaldofsets;

  for (std::vector<Teuchos::RCP<NodalDofSet>>::iterator it = nodaldofsets_.begin();
       it != nodaldofsets_.end(); it++)
  {
    Teuchos::RCP<NodalDofSet> nds = *it;

    bool is_std_dofset = nds->is_standard_dof_set();
    Core::Geo::Cut::Point::PointPosition pos = nds->position();

    // already an appropriate composite of nodal dofsets found, the current nodal dofset can be
    // combined with?
    Teuchos::RCP<CompositeNodalDofSet> cnds = Teuchos::null;

    if (is_std_dofset)  // do not combine standard dofsets as they are unique for each phase
    {
      cnds = Teuchos::rcp(new Core::Geo::Cut::CompositeNodalDofSet(is_std_dofset, pos));
      collected_nodaldofsets.push_back(cnds);
    }
    else  // ghost set -> create new collected set or append to an already existing one
    {
      if (collected_nodaldofsets.size() == 0)  // no composite added yet
      {
        cnds = Teuchos::rcp(new Core::Geo::Cut::CompositeNodalDofSet(
            is_std_dofset, pos));  // if first, then create a new composite
        collected_nodaldofsets.push_back(cnds);
      }
      else
      {
        // assume that the nodal dofsets have been sorted in a step before
        // then we potentially combine the current nodal dofset with the last CompositeNodalDofSet
        // at most
        Teuchos::RCP<CompositeNodalDofSet> cnds_last = collected_nodaldofsets.back();
        if (cnds_last == Teuchos::null)
          FOUR_C_THROW("there should be a valid CompositeNodalDofSet");

        if (connect_ghost_with_standard_nds)  // classical std-FEM based cut approximation --
                                              // combine ghost and std sets with same position
        {
          if (cnds_last->position() == pos)  // FIX!!!! same position (phase) and combine also
                                             // standard and ghost sets! Might change results!
            cnds = cnds_last;
          else  // different position, then create a new one!
          {
            cnds = Teuchos::rcp(new Core::Geo::Cut::CompositeNodalDofSet(
                is_std_dofset, pos));  // if first, then create a new composite
            collected_nodaldofsets.push_back(cnds);
          }
        }
        else  // combine only ghost dofset with each other, do not combine std with ghost sets as
              // usual in standard FEM
        {
          // is_std_dofset=false in this case
          if (cnds_last->is_standard_dof_set() == is_std_dofset and
              cnds_last->position() == pos)  // same position (phase) and also ghost dofset
            cnds = cnds_last;
          else
          {
            cnds = Teuchos::rcp(new Core::Geo::Cut::CompositeNodalDofSet(
                is_std_dofset, pos));  // if first, then create a new composite
            collected_nodaldofsets.push_back(cnds);
          }
        }
      }
    }

    cnds->add(nds, connect_ghost_with_standard_nds);
  }

  // set the composite of nodal dofsets for the node
  nodaldofsets_.clear();

  std::copy(collected_nodaldofsets.begin(), collected_nodaldofsets.end(),
      std::inserter(nodaldofsets_, nodaldofsets_.begin()));
}



/*-----------------------------------------------------------------------------------------*
 * build sets of connected volumecells in a 1-ring around the node
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::build_dof_cell_sets(Point* p,
    const std::vector<plain_volumecell_set>& cell_sets, const plain_volumecell_set& cells,
    const std::vector<plain_volumecell_set>& nodal_cell_sets, plain_volumecell_set& done,
    bool isnodalcellset)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "Core::Geo::CUT --- 5/6 --- cut_positions_dofsets --- build_dof_cell_sets");


  for (std::vector<plain_volumecell_set>::const_iterator s = nodal_cell_sets.begin();
       s != nodal_cell_sets.end(); s++)
  {
    const plain_volumecell_set& nodal_cells = *s;

    for (plain_volumecell_set::const_iterator i = nodal_cells.begin(); i != nodal_cells.end(); ++i)
    {
      VolumeCell* cell = *i;
      if (done.count(cell) == 0)
      {
        plain_volumecell_set connected;
        /* REMARK: here use the version without! elements check:
         * here we build cell sets within one global element with vcs of sub-elements
         * maybe the vcs of one sub-element are not connected within one sub-element,
         * but within one global element, therefore more than one vc of one
         * sub-element may be connected. */
        cell->neighbors(p, cells, done, connected);

        if (connected.size() > 0)
        {
          std::set<plain_volumecell_set, Cmp> connected_sets;

          // find all cells of connected in cell_sets and add the corresponding cell_sets
          for (plain_volumecell_set::iterator c = connected.begin(); c != connected.end(); c++)
          {
            VolumeCell* connected_cell = *c;

            for (std::vector<plain_volumecell_set>::const_iterator i = cell_sets.begin();
                 i != cell_sets.end(); i++)
            {
              // contains the current cell_it
              if ((*i).count(connected_cell) > 0)
              {
                connected_sets.insert(*i);
              }
            }
          }

          // set if this set is a std set
          nodaldofsets_.push_back(Teuchos::rcp(new NodalDofSet(connected_sets, isnodalcellset)));

          std::copy(connected.begin(), connected.end(), std::inserter(done, done.end()));
        }  // connected.size() > 0
      }    // done.count( cell )==0
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * build sets of connected volumecells in a 1-ring around the node (old unused version)
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::build_dof_cell_sets(Point* p, const plain_volumecell_set& cells,
    const plain_volumecell_set& nodal_cells, plain_volumecell_set& done)
{
  for (plain_volumecell_set::const_iterator i = nodal_cells.begin(); i != nodal_cells.end(); ++i)
  {
    VolumeCell* cell = *i;
    if (done.count(cell) == 0)
    {
      plain_volumecell_set connected;
      plain_element_set elements;
      cell->neighbors(p, cells, done, connected, elements);

      if (connected.size() > 0)
      {
        dofsets_.push_back(connected);

        std::copy(connected.begin(), connected.end(), std::inserter(done, done.begin()));
      }
    }
  }
}



/*-----------------------------------------------------------------------------------------*
 *  Gives this node a selfcutposition and spreads the positional information    wirtz 05/13
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::self_cut_position(Point::PointPosition pos)
{
  if (selfcutposition_ != pos)
  {
    FOUR_C_ASSERT(IsCutPositionUnchanged(selfcutposition_, pos),
        "Are you sure that you want to change the selfcut-node-position from inside to outside "
        "or vice versa?");

    // do not overwrite oncutsurface nodes
    if (selfcutposition_ == Point::oncutsurface) return;

    // change position for points just in case of undecided node and do not change oncutsurface
    // nodes
    if (selfcutposition_ == Point::undecided)
    {
      selfcutposition_ = pos;
      if (pos == Point::outside or pos == Point::inside)
      {
        for (plain_edge_set::iterator i = edges_.begin(); i != edges_.end(); ++i)
        {
          Edge* e = *i;
          e->self_cut_position(pos);
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------*
 *  Changes the selfcutposition of this node and spreads the positional information
 *                                                                               wirtz 07/16
 *-----------------------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::change_self_cut_position(Point::PointPosition pos)
{
  if (selfcutposition_ != pos)
  {
    selfcutposition_ = pos;
    for (plain_edge_set::iterator i = edges_.begin(); i != edges_.end(); ++i)
    {
      Edge* e = *i;
      e->change_self_cut_position(pos);
    }
  }
}

/*-----------------------------------------------------------------------------------------*
 *  Returns true if the given node is exactly at the same position as this node     sudhakar 09/13
 *-----------------------------------------------------------------------------------------*/
bool Core::Geo::Cut::Node::is_at_same_location(const Node* nod) const
{
  Core::LinAlg::Matrix<3, 1> nx1, nx2;

  point_->coordinates(nx1.data());
  coordinates(nx1.data());
  nod->coordinates(nx2.data());

  nx1.update(-1, nx2, 1);

  if (nx1.norm2() < (this->point()->tolerance() + nod->point()->tolerance())) return true;

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::Node::remove_non_standard_nodal_dof_sets()
{
  std::vector<Teuchos::RCP<NodalDofSet>> std_nodaldofset;
  std_nodaldofset.reserve(1);
  for (unsigned i = 0; i < static_cast<unsigned>(num_dof_sets()); ++i)
  {
    Teuchos::RCP<NodalDofSet>& nodaldofset = nodaldofsets_[i];
    if (nodaldofset->is_standard_dof_set())
    {
      std_nodaldofset.push_back(nodaldofset);
    }
    else
    {
      if (nodaldofset.strong_count() > 1)
        FOUR_C_THROW(
            "nodaldofset cannot be destroyed! (strong_count = %d)", nodaldofset.strong_count());
      nodaldofset = Teuchos::null;
    }
  }
  nodaldofsets_.swap(std_nodaldofset);
}



/*-----------------------------------------------------------------------------------------*
 * get the unique standard NodalDofSet for a given nodal dofset position
 *-----------------------------------------------------------------------------------------*/
int Core::Geo::Cut::Node::get_standard_nodal_dof_set(Point::PointPosition pos)
{
  for (int i = 0; i < num_dof_sets(); i++)
  {
    Teuchos::RCP<NodalDofSet> nodaldofset = nodaldofsets_[i];
    if (nodaldofset->is_standard_dof_set())
    {
      if (nodaldofset->position() == pos) return i;
    }
  }

  return -1;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
bool Core::Geo::Cut::NodalDofSetCmp::operator()(
    Teuchos::RCP<NodalDofSet> nodaldofset1, Teuchos::RCP<NodalDofSet> nodaldofset2)
{
  //==============================================
  // classical sorting of dofsets with:
  // - possible positions: outside/inside
  // in combination with
  // - at most one standard dofset per position (phase)
  // - arbitrary number of ghost dofsets per position (phase)
  //==============================================
  // nds-order starting from 0: STD(outside)/STD(inside) // GHOST(outside)_1 / GHOST(outside)_2 ...
  // GHOST(outside)_nout // GHOST(inside)_1 / GHOST(inside)_2 ... GHOST(inside)_nin where GHOST(*)_1
  // ... GHOST(*)_n are sorted by PointIds of contained volumecells
  //==============================================


  // FIRST: sort by std vs ghost nodal dofset if possible
  if (nodaldofset1->is_standard_dof_set() !=
      nodaldofset2->is_standard_dof_set())  // one set is standard, the other is ghost
    return nodaldofset1->is_standard_dof_set() <
           nodaldofset2->is_standard_dof_set();  // std=0, ghost=1 => STD before GHOST

  // now: both nodal dofset are standard dofsets or both nodal dofsets are ghost dofsets

  // SECOND: sort nodal dofset by position
  if (nodaldofset1->position() !=
      nodaldofset2->position())  // the sets belong to two different phases
    return nodaldofset1->position() <
           nodaldofset2->position();  // compare the enum: outside=-3 , inside=-2

  // now: both nodal dofset are standard dofsets or both nodal dofsets are ghost dofsets
  //  AND both nodal dofsets have the same position

  // THIRD: sort by point ids of the contained volume-cell's points
  const std::set<plain_volumecell_set, Cmp>& composite1 = nodaldofset1->volume_cell_composite();
  const std::set<plain_volumecell_set, Cmp>& composite2 = nodaldofset2->volume_cell_composite();

  Core::Geo::Cut::Cmp comp;

  return comp.compare(composite1, composite2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Core::Geo::Cut::FindCommonElements(
    const std::vector<Node*>& nelement, plain_element_set& elements)
{
  // find the element defined by the given nodes
  std::vector<Core::Geo::Cut::Point*> pelement;
  pelement.reserve(nelement.size());

  for (std::vector<Core::Geo::Cut::Node*>::const_iterator cit = nelement.begin();
       cit != nelement.end(); ++cit)
    pelement.push_back((*cit)->point());

  FindCommonElements(pelement, elements);
}

FOUR_C_NAMESPACE_CLOSE
