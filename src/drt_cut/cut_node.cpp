
#include "cut_element.H"
#include "cut_volumecell.H"




/*------------------------------------------------------------------------------*
  | Operator () compare operator for plain_volumecell_sets
  |                                                             shahmiri 06/12
  *-----------------------------------------------------------------------------*/
bool GEO::CUT::Cmp::operator()(
    plain_volumecell_set s1,
    plain_volumecell_set s2
)
{
  // call Compare function for two plain_volumecell_sets
  return Compare(s1, s2);
}


/*------------------------------------------------------------------------------*
  | Operator () compare operator for sets of plain_volumecell_sets
  |                                                             shahmiri 06/12
  *-----------------------------------------------------------------------------*/
bool GEO::CUT::Cmp::operator()(
    const std::set<plain_volumecell_set, Cmp>& set1,
    const std::set<plain_volumecell_set, Cmp>& set2
)
{
  // call Compare function for two sets of plain_volumecell_sets
  return Compare(set1, set2);
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
  plain_volumecell_set s1 = *it1;

  std::set<plain_volumecell_set, Cmp>::iterator it2 = set2.begin();
  plain_volumecell_set s2 = *it2;

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
  // compare the point ids of the two volumecells

  std::set<int> vc1points = vc1->VolumeCellPointIds();
  std::set<int> vc2points = vc2->VolumeCellPointIds();


  if (vc1 == vc2)
     return false;

  // first build two minimized sets which don't have any points in common
  std::set<int> vc1pointsmin;
  std::set<int> vc2pointsmin;
  for( std::set<int>::iterator iter=vc1points.begin(); iter!=vc1points.end(); iter++ )
  {
    std::set<int>::iterator iter2 = vc2points.find(*iter);
    if (iter2 == vc2points.end())
    {
      vc1pointsmin.insert(*iter);
    }
  }

  for( std::set<int>::iterator iter=vc2points.begin(); iter!=vc2points.end(); iter++ )
  {
    std::set<int>::iterator iter2 = vc1points.find(*iter);
    if (iter2 == vc1points.end())
    {
      vc2pointsmin.insert(*iter);
    }
  }

  std::set<int>::iterator p1=vc1pointsmin.begin();
  std::set<int>::iterator p2=vc2pointsmin.begin();

  while ( p1!=vc1pointsmin.end() and p2!=vc2pointsmin.end() )
  {
    int id1 = *p1;
    int id2 = *p2;

    if( id1 < id2)
    {
      return true;
      break;
    }
    if( id1 > id2)
    {
      return false;
      break;
    }

    p1++;
    p2++;
  }

   return false;
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
void GEO::CUT::Node::AssignNodalCellSet( std::vector<plain_volumecell_set> & ele_vc_sets,
                                         std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets)
{
    for( std::vector<plain_volumecell_set>::iterator s=ele_vc_sets.begin();
         s!=ele_vc_sets.end();
         s++)
    {
        plain_volumecell_set  cell_set = *s;

        for ( plain_volumecell_set::const_iterator i=cell_set.begin(); i!=cell_set.end(); ++i )
        {
            VolumeCell * cell = *i;

            // if at least one cell of this cell_set contains the point, then the whole cell_set contains the point
            if ( cell->Contains( point() ) )
            {
                nodal_cell_sets[this].push_back( cell_set );

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

    // finally: fill dof_cellsets_

    // do the connection between elements
    plain_volumecell_set done;
    plain_volumecell_set cells;

    // get cell_sets as a plain_volume_set
    for(std::vector<plain_volumecell_set>::iterator i=cell_sets.begin(); i!=cell_sets.end(); i++)
    {
       std::copy((*i).begin(), (*i).end(), std::inserter( cells, cells.begin() ) );
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
 * get the dofset number of the Volumecell w.r.t this node
 *-----------------------------------------------------------------------------------------*/
int GEO::CUT::Node::DofSetNumberNEW( plain_volumecell_set & cells )
{
  int dofset = -1;

  // find the first cell of cells, this is only a volume cell of a subelement
  if(cells.size() == 0) dserror( "cells is empty");

//  VolumeCell* cell = cells[0];
  VolumeCell* cell = *(cells.begin());

  for ( unsigned int i=0; i<dof_cellsets_.size(); ++i ) // loop over sets
  {
    std::set< plain_volumecell_set, GEO::CUT::Cmp >& cellsets = dof_cellsets_[i];

    for(std::set<plain_volumecell_set, GEO::CUT::Cmp >::iterator j = cellsets.begin();
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
void GEO::CUT::Node::SortDOFCellSets()
{
  // check if the first set is std set

  // REMARK:
  // first_set_is_std_set_ = boolian if the first set in dof_cellsets_ is a std set or not
  // if it is then the first set must not be changed during sorting
  // because elements without eh assume the first set as std set!

  if(dof_cellsets_.size() > 1)
  {

    // set the start iterator for sorting
    if(first_set_is_std_set_)
    {
      std::vector<std::set<plain_volumecell_set, Cmp> >::iterator it_start = (dof_cellsets_.begin())+1;

      // sort the cellsets w.r.t point ids in first vc in first set of sorted sets of plain volume cells sets
      // REMARK: do not sort the first dofset, it has to be kept the standard dofset
      sort(it_start, dof_cellsets_.end(), Cmp());
    }
    else
    {

      std::vector<std::set<plain_volumecell_set, Cmp> >::iterator it_start = dof_cellsets_.begin();

      // sort the cellsets w.r.t point ids in first vc in first set of sorted sets of plain volume cells sets
      sort(it_start, dof_cellsets_.end(), Cmp());

    }
  }

  return;
}


/*-----------------------------------------------------------------------------------------*
 * build sets of connected volumecells in a 1-ring around the node
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::BuildDOFCellSets( Point * p,
                                       const std::vector<plain_volumecell_set> & cell_sets,
                                       const plain_volumecell_set & cells,
                                       const std::vector<plain_volumecell_set> & nodal_cell_sets,
                                       plain_volumecell_set & done,
                                       bool isnodalcellset)
{
    for( std::vector<plain_volumecell_set>::const_iterator s=nodal_cell_sets.begin(); s!=nodal_cell_sets.end(); s++)
    {
        plain_volumecell_set nodal_cells = *s;

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

                  dof_cellsets_.push_back(connected_sets);

                  // set if this set is a std set, but only if this function was called with a nodalcellset
                  if(isnodalcellset) first_set_is_std_set_=true;

                  std::copy( connected.begin(), connected.end(), std::inserter( done, done.begin() ) );

                }
            }
        }
    }

}


/*-----------------------------------------------------------------------------------------*
 * build sets of connected volumecells in a 1-ring around the node (old unused version)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::Node::BuildDOFCellSets( Point * p,
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



#if 0
int GEO::CUT::Node::NumDofSets( bool include_inner )
{
  if ( include_inner )
  {
    return DofSets().size();
  }
  else
  {
    int numdofsets = 0;
    for ( std::vector<plain_volumecell_set >::iterator i=dofsets_.begin();
          i!=dofsets_.end();
          ++i )
    {
      plain_volumecell_set & cells = *i;
      GEO::CUT::Point::PointPosition position = GEO::CUT::Point::undecided;
      for ( plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        VolumeCell * c = *i;
        GEO::CUT::Point::PointPosition cp = c->Position();
        if ( cp == GEO::CUT::Point::undecided )
        {
          throw std::runtime_error( "undecided volume cell position" );
        }
        if ( position!=GEO::CUT::Point::undecided and position!=cp )
        {
          throw std::runtime_error( "mixed volume cell set" );
        }
        position = cp;
      }
      if ( position==GEO::CUT::Point::outside )
      {
        numdofsets += 1;
      }
    }
    return numdofsets;
  }
}
#endif

/*-----------------------------------------------------------------------------------------*
 *  Gives this node a selfcutposition and spread the positional information     wirtz 05/13
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
    if( (selfcutposition_ == Point::undecided) )
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
 *  Returns true if the given node is exactly at the same position as this node     sudhakar 09/13
 *-----------------------------------------------------------------------------------------*/
bool GEO::CUT::Node::isAtSameLocation( const Node * nod ) const
{
  LINALG::Matrix<3,1> nx1, nx2;

  point_->Coordinates( nx1.A() );
  Coordinates( nx1.A() );
  nod->Coordinates( nx2.A() );

  nx1.Update( -1, nx2, 1 );

  if ( nx1.Norm2() < MINIMALTOL )
    return true;

  return false;
}
