/*!-----------------------------------------------------------------------------------------------*
 \file cut_node.cpp

 \brief class representing a geometrical node

 <pre>
\level 3
\maintainer Andy Wirtz
 wirtz@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15270
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "cut_element.H"
#include "cut_volumecell.H"

#include <Teuchos_TimeMonitor.hpp>



/*------------------------------------------------------------------------------*
  | Operator () compare operator for plain_volumecell_sets
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool GEO::CUT::Cmp::operator()(
    const plain_volumecell_set& s1,
    const plain_volumecell_set& s2
)
{
  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- GEO::CUT::Cmp::operator()" );

  // call Compare function for two plain_volumecell_sets
  return Compare(s1, s2);
}



/*------------------------------------------------------------------------------*
  | Compare() to compare two sets of plain_volumecell_set via comparing their first plain_volumecell_sets
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool GEO::CUT::Cmp::Compare(
    const std::set<plain_volumecell_set, Cmp>& set1,
    const std::set<plain_volumecell_set, Cmp>& set2 )
{

  // compare two sets of plain_volumecell_set
  // take the first plain_volumecell_set of each set and compare them

  std::set<plain_volumecell_set, Cmp>::iterator it1 = set1.begin();
  const plain_volumecell_set & s1 = *it1;

  std::set<plain_volumecell_set, Cmp>::iterator it2 = set2.begin();
  const plain_volumecell_set & s2 = *it2;

  return Compare(s1, s2);
}

/*------------------------------------------------------------------------------*
  | Compare() to compare two plain_volumecell_set via the ids of their first volumecell's points
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool GEO::CUT::Cmp::Compare(
    const plain_volumecell_set& s1,
    const plain_volumecell_set& s2 )
{

  // take the first vc in plain_volumecell_set. In case of linear elements
  // this is the only volumecell. For quadratic elements this is a set of vcs
  // but we just compare the first vc of each plain_volumecell_set. So for
  // quadratic elements the  plain_volumecell_set is not sorted!

  plain_volumecell_set::const_iterator it1 = s1.begin();
  VolumeCell * vc1 = *it1;

  plain_volumecell_set::const_iterator it2 = s2.begin();
  VolumeCell * vc2 = *it2;

  return Compare(vc1, vc2);
}


/*------------------------------------------------------------------------------*
  | Operator () to compare two volume cells via the ids of their points
  |                                                             shahmiri 06/12
 *-----------------------------------------------------------------------------*/
bool GEO::CUT::Cmp::Compare(
    VolumeCell* vc1,
    VolumeCell* vc2 )
{
  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- GEO::CUT::Cmp::Compare(vc,vc)" );


  if (vc1 == vc2)
    return false;

  // compare the point ids of the two volumecells

  const std::set<int> & vc1points = vc1->VolumeCellPointIds();
  const std::set<int> & vc2points = vc2->VolumeCellPointIds();


  // during reducing the two sets to two minimized/disjoint sets which don't have any points in common
  // the first non-common point ids can be used to sort the sets

  std::set<int>::iterator it1 = vc1points.begin();
  std::set<int>::iterator it2 = vc2points.begin();

  while(true) // sets are already sorted!
  {
    if(it1 == vc1points.end() or it2 == vc2points.end()) break;

    if(*it1 == *it2)
    {
      // identical points, go to the next one
      it1++;
      it2++;
    }
    else // compare the points now via their point ids
    {
      if( *it1 < *it2 ) return true;
      else return false; // (*it1 > *it2)
    }
  }

  dserror("sorting failed: one volume-cell is completely contained in the other volume-cell or both vcs are equal!");

  return false;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::NodalDofSet::NodalDofSet(
    std::set<plain_volumecell_set,Cmp>& connected_volumecells,
    bool is_std_dofset
    ) : is_std_dofset_(is_std_dofset)
{
  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- GEO::CUT::NodalDofSet::NodalDofSet" );

  std::copy(
      connected_volumecells.begin(),
      connected_volumecells.end(),
      std::inserter(volumecell_composite_,volumecell_composite_.end())
  );

  // set the position of the nodal dofset
  const plain_volumecell_set & set = *(volumecell_composite_.begin());
  position_ = (*set.begin())->Position();

}


GEO::CUT::NodalDofSet::NodalDofSet(
    bool is_std_dofset,
    GEO::CUT::Point::PointPosition pos
) : is_std_dofset_(is_std_dofset), position_(pos)
{
  volumecell_composite_.clear();
}

/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
bool GEO::CUT::NodalDofSet::Contains( GEO::CUT::Point* p)
{
  // check if any volume-cell of the volumecell_composite contains this point
  for(std::set<plain_volumecell_set,Cmp>::iterator it = volumecell_composite_.begin();
      it!=volumecell_composite_.end();
      it++)
  {
    const plain_volumecell_set & vc_set = *it;
    for(plain_volumecell_set::const_iterator vcs = vc_set.begin();
        vcs!=vc_set.end();
        vcs++)
    {
      if((*vcs)->Contains(p)) return true;
    }
  }

  return false;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::NodalDofSet::CollectCutSides( GEO::CUT::plain_int_set & cutside_ids ) const
{
  // collect all cut sides
  for(std::set<plain_volumecell_set,Cmp>::iterator it = volumecell_composite_.begin();
      it!=volumecell_composite_.end();
      it++)
  {
    const plain_volumecell_set & vc_set = *it;
    for(plain_volumecell_set::const_iterator vcs = vc_set.begin();
        vcs!=vc_set.end();
        vcs++)
    {
      (*vcs)->CollectCutSides(cutside_ids);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::NodalDofSet::Print()
{
  std::cout << "GEO::CUT::NodalDofSet: "
      << "STD dofset? " << this->Is_Standard_DofSet()
      << " Pos: " << this->Position()
      << std::endl;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::CompositeNodalDofSet::Print()
{
  std::cout << "GEO::CUT::CompositeNodalDofSet which contains " << nodal_dofsets_.size() << " combined GEO::CUT::NodalDofSet: "
      << "STD dofset? " << this->Is_Standard_DofSet()
      << " Pos: " << this->Position()
      << std::endl;
}


/*-----------------------------------------------------------------------------------------*
 * register cuts
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::RegisterCuts()
{
  if ( Position()==Point::oncutsurface )
  {
    for ( plain_edge_set::iterator i=edges_.begin(); i!=edges_.end(); ++i )
    {
      Edge * e = *i;
      point_->AddEdge( e );
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * Assign the vc_sets to the node if possible
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::AssignNodalCellSet(
    const std::vector<plain_volumecell_set> & ele_vc_sets,
    std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- AssignNodalCellSet" );

  std::vector<plain_volumecell_set> & nodal_cell_set = nodal_cell_sets[this];

  for( std::vector<plain_volumecell_set>::const_iterator s=ele_vc_sets.begin();
      s!=ele_vc_sets.end();
      s++)
  {
    const plain_volumecell_set & cell_set = *s;

    for ( plain_volumecell_set::const_iterator i=cell_set.begin(); i!=cell_set.end(); ++i )
    {
      VolumeCell * cell = *i;

      // if at least one cell of this cell_set contains the point, then the whole cell_set contains the point
      bool contains = false;

      {
        TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- cell->Contains( point() )" );

        contains = cell->Contains( point() );
      }

      if ( contains )
      {
        nodal_cell_set.push_back( cell_set );

        // the rest of cells in this set has not to be checked for this node
        break; // finish the cell_set loop, breaks the inner for loop!
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * Find the dofsets required at this node. (old unused version)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::FindDOFSets( bool include_inner )
{
  const plain_element_set & elements = Elements();

  std::map<Node*, plain_volumecell_set > nodal_cells;
  plain_volumecell_set cells;

  for ( plain_element_set::const_iterator i=elements.begin(); i!=elements.end(); ++i )
  {
    Element * e = *i;
    {
      const plain_volumecell_set & element_cells = e->VolumeCells();
      if ( include_inner )
      {
        std::copy( element_cells.begin(), element_cells.end(), std::inserter( cells, cells.begin() ) );
      }
      else
      {
        for ( plain_volumecell_set::const_iterator i=element_cells.begin(); i!=element_cells.end(); ++i )
        {
          VolumeCell * vc = *i;
          if ( vc->Position()==Point::outside )
          {
            cells.insert( vc );
          }
        }
      }

      const std::vector<Node*> & nodes = e->Nodes();
      for ( std::vector<Node*>::const_iterator i=nodes.begin();
          i!=nodes.end();
          ++i )
      {
        Node * n = *i;
        for ( plain_volumecell_set::const_iterator i=element_cells.begin(); i!=element_cells.end(); ++i )
        {
          VolumeCell * cell = *i;

          if ( not include_inner and cell->Position()!=Point::outside )
            continue;

          if ( cell->Contains( n->point() ) )
          {
            nodal_cells[n].insert( cell );
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

  BuildDOFCellSets( point(), cells, nodal_cells[this], done );

  nodal_cells.erase( this );

  for ( std::map<Node*, plain_volumecell_set >::iterator i=nodal_cells.begin();
      i!=nodal_cells.end();
      ++i )
  {
    Node * n = i->first;
    plain_volumecell_set & cellset = i->second;
    BuildDOFCellSets( n->point(), cells, cellset, done );
  }

  // do any remaining internal volumes that are not connected to any node
  BuildDOFCellSets( NULL, cells, cells, done );

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
void GEO::CUT::Node::FindDOFSetsNEW( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets,
    std::vector<plain_volumecell_set> & cell_sets)
{

  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- FindDOFSetsNEW" );


  // finally: fill dof_cellsets_

  // do the connection between elements
  plain_volumecell_set done;
  plain_volumecell_set cells;

  // get cell_sets as a plain_volume_set
  for(std::vector<plain_volumecell_set>::iterator i=cell_sets.begin(); i!=cell_sets.end(); i++)
  {
    std::copy((*i).begin(), (*i).end(), std::inserter( cells, cells.end() ) );
  }


  // First, get the nodal cells that make up the first dofset. In most cases
  // this loop has one pass only. But if the node is cut, there will be more
  // than one set of cells that are attached to this node.

  // call this function with isnodalcellset=true flag to identify the first std set
  BuildDOFCellSets( point(), cell_sets, cells, nodal_cell_sets[this], done, true );

  nodal_cell_sets.erase( this );


  for ( std::map<Node*, std::vector<plain_volumecell_set> >::iterator i=nodal_cell_sets.begin();
      i!=nodal_cell_sets.end();
      ++i )
  {
    Node * n = i->first;

    std::vector<plain_volumecell_set> & cellset = i->second;
    BuildDOFCellSets( n->point(), cell_sets, cells, cellset, done );

  }

  // do any remaining internal volumes that are not connected to any node
  BuildDOFCellSets( NULL, cell_sets, cells, cell_sets, done );


}


/*-----------------------------------------------------------------------------------------*
 * get the dofset number of the Volumecell w.r.t this node (old unused version)
 *-----------------------------------------------------------------------------------------*/
int GEO::CUT::Node::DofSetNumber( VolumeCell * cell )
{
  int dofset = -1;
  for ( unsigned i=0; i<dofsets_.size(); ++i )
  {
    plain_volumecell_set & cells = dofsets_[i];
    if ( cells.count( cell ) > 0 )
    {
      if ( dofset==-1 )
      {
        dofset = i;
      }
      else
      {
        throw std::runtime_error( "volume dofset not unique" );
      }
    }
  }
  if ( dofset==-1 )
  {
    std::cout << "dofset not found for node " << this->Id() << std::endl;
    throw std::runtime_error( "volume dofset not found" );
  }
  return dofset;
}


/*-----------------------------------------------------------------------------------------*
 * get number of dofsets
 *-----------------------------------------------------------------------------------------*/
int GEO::CUT::Node::NumDofSets() const{
  return NodalDofSets().size();
}




/*-----------------------------------------------------------------------------------------*
 * get the dofset number of the Volumecell w.r.t this node
 *-----------------------------------------------------------------------------------------*/
int GEO::CUT::Node::DofSetNumberNEW( const plain_volumecell_set & cells )
{
  int dofset = -1;

  // find the first cell of cells, this is only a volume cell of a subelement
  if(cells.size() == 0) dserror( "cells is empty");

  //  VolumeCell* cell = cells[0];
  VolumeCell* cell = *(cells.begin());

  for ( unsigned int i=0; i<nodaldofsets_.size(); ++i ) // loop over sets
  {
    const std::set< plain_volumecell_set, GEO::CUT::Cmp >& cellsets = nodaldofsets_[i]->VolumeCellComposite();

    for(std::set<plain_volumecell_set, GEO::CUT::Cmp >::const_iterator j = cellsets.begin();
        j!=cellsets.end();
        j++)
    {

      if ( j->count( cell ) > 0 )
      {
        if ( dofset==-1 )
        {
          dofset = i;
        }
        else
        {
          std::cout << "node: " << Id() << std::endl;
          std::cout << "first dofset id: " << dofset << std::endl;
          std::cout << "new dofset id: " << i << std::endl;
          cell->Print(std::cout);
          throw std::runtime_error( "volume dofset not unique" );
        }
      }
    }
  }
  if ( dofset==-1 )
  {
    std::cout << "dofset not found for node " << this->Id() << std::endl;
    //    throw std::runtime_error( "volume dofset not found" );
  }
  return dofset;
}



/*-----------------------------------------------------------------------------------------*
 * sort all dofsets via xyz point coordinates (use compare functions in cut_node.H)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::SortNodalDofSets()
{
  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- SortNodalDofSets" );


  if(nodaldofsets_.size() > 1)
  {
    std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> >::iterator it_start = nodaldofsets_.begin();

    if(nodaldofsets_[0]->Is_Standard_DofSet())
    {
      // REMARK:
      // if the first nodal dofset is a standard dofset, then the first set must not be changed during sorting
      // because elements without eh assume the first set as std set!

      it_start ++;  // exclude the standard set from sorting
    }

    // sort the cellsets w.r.t point ids in first vc in first set of sorted sets of plain volume cells sets
    // REMARK: do not sort the first dofset, it has to be kept the standard dofset
    sort(it_start, nodaldofsets_.end(), NodalDofSetCmp());
  }

#if(0)
  // print the sorted dofsets:
  std::cout << "Sorted DOFSETs for node: " << Id() << std::endl;
  for(int i=0; i< NumDofSets(); i++)
  {
    NodalDofSet * dofset = GetNodalDofSet(i);
    dofset->Print();
  }
#endif

  return;
}


/*-----------------------------------------------------------------------------------------*
 * collect the (ghost) dofsets for this node w.r.t each phase to avoid multiple ghost nodal dofsets for a certain phase
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::CollectNodalDofSets()
{
  //TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- CollectNodalDofSets" );

  // assume that the nodal dofsets have been sorted in a step before,
  // such that ghost sets to be combined are stored consecutively in the vector of sorted nodal dofsets

  std::vector<Teuchos::RCP<CompositeNodalDofSet> > collected_nodaldofsets;

  for(std::vector<Teuchos::RCP<NodalDofSet> >::iterator it=nodaldofsets_.begin();
      it!=nodaldofsets_.end();
      it++)
  {
    Teuchos::RCP<NodalDofSet> nds = *it;
    bool is_std_dofset = nds->Is_Standard_DofSet();
    GEO::CUT::Point::PointPosition pos = nds->Position();

    // already an appropriate composite of nodal dofsets found, the current nodal dofset can be combined with?
    Teuchos::RCP<CompositeNodalDofSet> cnds = Teuchos::null;

    if(is_std_dofset) // do not combine standard dofsets as they are unique for each phase
    {
      cnds = Teuchos::rcp(new GEO::CUT::CompositeNodalDofSet(is_std_dofset, pos));
      collected_nodaldofsets.push_back(cnds);
    }
    else // ghost set -> create new collected set or append to an already existing one
    {

      if( collected_nodaldofsets.size() == 0 ) // no composite added yet
      {
        cnds = Teuchos::rcp(new GEO::CUT::CompositeNodalDofSet(is_std_dofset, pos)); // if first, then create a new composite
        collected_nodaldofsets.push_back(cnds);
      }
      else
      {
        // assume that the nodal dofsets have been sorted in a step before
        // then we potentially combine the current nodal dofset with the last CompositeNodalDofSet at most
        Teuchos::RCP<CompositeNodalDofSet> cnds_last = collected_nodaldofsets.back();

        if(cnds_last->Is_Standard_DofSet() == is_std_dofset and
            cnds_last->Position() == pos) // same position (phase) and also ghost
          cnds = cnds_last;
        else
        {
          cnds = Teuchos::rcp(new GEO::CUT::CompositeNodalDofSet(is_std_dofset, pos)); // if first, then create a new composite
          collected_nodaldofsets.push_back(cnds);
        }
      }
    }

    cnds->add(nds);
  }

  // set the composite of nodal dofsets for the node
  nodaldofsets_.clear();

  std::copy(
      collected_nodaldofsets.begin(),
      collected_nodaldofsets.end(),
      std::inserter(nodaldofsets_,nodaldofsets_.begin())
  );

  // safety check for number of allowed sets (one std and one ghost per position)


#if(0)
  // print the sorted dofsets:
  std::cout << "Collected DOFSETs for node: " << Id() << std::endl;
  for(int i=0; i< NumDofSets(); i++)
  {
    NodalDofSet * dofset = GetNodalDofSet(i);
    dofset->Print();
  }
#endif
}



/*-----------------------------------------------------------------------------------------*
 * build sets of connected volumecells in a 1-ring around the node
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::BuildDOFCellSets(
    Point * p,
    const std::vector<plain_volumecell_set> & cell_sets,
    const plain_volumecell_set & cells,
    const std::vector<plain_volumecell_set> & nodal_cell_sets,
    plain_volumecell_set & done,
    bool isnodalcellset
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "GEO::CUT --- 5/6 --- Cut_Positions_Dofsets --- BuildDOFCellSets" );


  for( std::vector<plain_volumecell_set>::const_iterator s=nodal_cell_sets.begin(); s!=nodal_cell_sets.end(); s++)
  {
    const plain_volumecell_set & nodal_cells = *s;

    for ( plain_volumecell_set::const_iterator i=nodal_cells.begin();
        i!=nodal_cells.end();
        ++i )
    {
      VolumeCell * cell = *i;
      if ( done.count( cell )==0 )
      {
        plain_volumecell_set connected;
        // REMARK: here use the version without! elements check:
        // here we build cell sets within one global element with vcs of subelements
        // maybe the vcs of one subelement are not connected within one subelement, but within one global element,
        // therefore more than one vc of one subelements may be connected.
        cell->Neighbors( p, cells, done, connected);

        if ( connected.size()>0 )
        {

          std::set<plain_volumecell_set,Cmp> connected_sets;

          int count=0;
          // find all cells of connected in cell_sets and add the corresponding cell_sets
          for(plain_volumecell_set::iterator c=connected.begin(); c!= connected.end(); c++)
          {
            VolumeCell* connected_cell=*c;
            count++;

            int cell_it = 0;
            for(std::vector<plain_volumecell_set>::const_iterator i=cell_sets.begin();
                i!=cell_sets.end();
                i++ )
            {
              cell_it++;

              // contains the current cell_it
              if((*i).count( connected_cell ) > 0)
              {
                connected_sets.insert(*i);
              }
            }
          }

          // set if this set is a std set
          nodaldofsets_.push_back(Teuchos::rcp(new NodalDofSet(connected_sets, isnodalcellset)));

          std::copy( connected.begin(), connected.end(), std::inserter( done, done.end() ) );
        } // connected.size() > 0
      } // done.count( cell )==0
    }
  }

}


/*-----------------------------------------------------------------------------------------*
 * build sets of connected volumecells in a 1-ring around the node (old unused version)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::BuildDOFCellSets(
    Point * p,
    const plain_volumecell_set & cells,
    const plain_volumecell_set & nodal_cells,
    plain_volumecell_set & done )
{
  for ( plain_volumecell_set::const_iterator i=nodal_cells.begin();
      i!=nodal_cells.end();
      ++i )
  {
    VolumeCell * cell = *i;
    if ( done.count( cell )==0 )
    {
      plain_volumecell_set connected;
      plain_element_set elements;
      cell->Neighbors( p, cells, done, connected, elements );

      if ( connected.size()>0 )
      {
        dofsets_.push_back( connected );

        std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );
      }
    }
  }
}



/*-----------------------------------------------------------------------------------------*
 *  Gives this node a selfcutposition and spreads the positional information    wirtz 05/13
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::SelfCutPosition( Point::PointPosition pos )
{
  if ( selfcutposition_ != pos )
  {

    //#ifdef DEBUGCUTLIBRARY
    if( (selfcutposition_ == Point::inside and pos == Point::outside) or
        (selfcutposition_ == Point::outside and pos == Point::inside) )
    {

      std::cout << "selfcutnode with changing position inside->outside or vice versa " << nid_ << std::endl;
      throw std::runtime_error("Are you sure that you want to change the selfcut-node-position from inside to outside or vice versa?");
    }
    //#endif

    // do not overwrite oncutsurface nodes
    if(selfcutposition_ == Point::oncutsurface) return;

    // change position for points just in case of undecided node and do not change oncutsurface nodes
    if( selfcutposition_ == Point::undecided )
    {
      selfcutposition_ = pos;
      if ( pos==Point::outside or pos==Point::inside )
      {
        for ( plain_edge_set::iterator i=edges_.begin(); i!=edges_.end(); ++i )
        {
          Edge * e = *i;
          e->SelfCutPosition( pos );
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------*
 *  Changes the selfcutposition of this node and spreads the positional information
 *                                                                               wirtz 07/16
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::ChangeSelfCutPosition( Point::PointPosition pos )
{
  if ( selfcutposition_ != pos )
  {
    selfcutposition_ = pos;
    for (plain_edge_set::iterator i = edges_.begin(); i != edges_.end(); ++i)
    {
      Edge * e = *i;
      e->ChangeSelfCutPosition(pos);
    }
  }
}

/*-----------------------------------------------------------------------------------------*
 *  Returns true if the given node is exactly at the same position as this node     sudhakar 09/13
 *-----------------------------------------------------------------------------------------*/
bool GEO::CUT::Node::isAtSameLocation( const Node * nod ) const
{
  LINALG::Matrix<3,1> nx1, nx2;

  point_->Coordinates( nx1.A() );
  Coordinates( nx1.A() );
  nod->Coordinates( nx2.A() );

  nx1.Update( -1, nx2, 1 );

  if ( nx1.Norm2() < (this->point()->Tolerance() + nod->point()->Tolerance() ))
    return true;

  return false;
}

/*-----------------------------------------------------------------------------------------*
 * get the unique standard NodalDofSet for a given nodal dofset position
 *-----------------------------------------------------------------------------------------*/
int GEO::CUT::Node::GetStandardNodalDofSet( Point::PointPosition pos)
{
  for(int i=0; i< NumDofSets(); i++)
  {
    Teuchos::RCP<NodalDofSet> nodaldofset = nodaldofsets_[i];
    if( nodaldofset->Is_Standard_DofSet() )
    {
      if( nodaldofset->Position() == pos)
        return i;
    }
  }

  return -1;
}


/*-----------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------*/
bool GEO::CUT::NodalDofSetCmp::operator()(
    Teuchos::RCP<NodalDofSet> nodaldofset1,
    Teuchos::RCP<NodalDofSet> nodaldofset2
)
{
  //==============================================
  // classical sorting of dofsets with:
  // - possible positions: outside/inside
  // in combination with
  // - at most one standard dofset per position (phase)
  // - arbitrary number of ghost dofsets per position (phase)
  //==============================================
  // nds-order starting from 0: STD(outside)/STD(inside) // GHOST(outside)_1 / GHOST(outside)_2 ... GHOST(outside)_nout // GHOST(inside)_1 / GHOST(inside)_2 ... GHOST(inside)_nin
  // where GHOST(*)_1 ... GHOST(*)_n are sorted by PointIds of contained volumecells
  //==============================================


  // FIRST: sort by std vs ghost nodal dofset if possible
  if( nodaldofset1->Is_Standard_DofSet() != nodaldofset2->Is_Standard_DofSet() )  // one set is standard, the other is ghost
    return nodaldofset1->Is_Standard_DofSet() < nodaldofset2->Is_Standard_DofSet(); // std=0, ghost=1 => STD before GHOST

  // now: both nodal dofset are standard dofsets or both nodal dofsets are ghost dofsets

  // SECOND: sort nodal dofset by position
  if( nodaldofset1->Position() != nodaldofset2->Position() ) // the sets belong to two different phases
    return nodaldofset1->Position() < nodaldofset2->Position(); // compare the enum: outside=-3 , inside=-2

  // now: both nodal dofset are standard dofsets or both nodal dofsets are ghost dofsets
  //  AND both nodal dofsets have the same position

  // THIRD: sort by point ids of the contained volume-cell's points
  const std::set<plain_volumecell_set, Cmp>&  composite1 = nodaldofset1->VolumeCellComposite();
  const std::set<plain_volumecell_set, Cmp>&  composite2 = nodaldofset2->VolumeCellComposite();

  GEO::CUT::Cmp comp;

  return comp.Compare(composite1, composite2);
}
