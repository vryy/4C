
#include "cut_element.H"
#include "cut_volumecell.H"
#include "cut_node.H"


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
    throw std::runtime_error( "volume dofset not found" );
  }
  return dofset;
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
