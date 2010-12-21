
#include "cut_element.H"
#include "cut_volumecell.H"
#include "cut_node.H"

void GEO::CUT::Node::FindDOFSets()
{
  const std::set<Element*> & elements = Elements();

  std::map<Node*, std::set<VolumeCell*> > nodal_cells;
  std::set<VolumeCell*> cells;

  for ( std::set<Element*>::const_iterator i=elements.begin(); i!=elements.end(); ++i )
  {
    Element * e = *i;
    {
      const std::set<VolumeCell*> & element_cells = e->VolumeCells();
      std::copy( element_cells.begin(), element_cells.end(), std::inserter( cells, cells.begin() ) );

      const std::vector<Node*> & nodes = e->Nodes();
      for ( std::vector<Node*>::const_iterator i=nodes.begin();
            i!=nodes.end();
            ++i )
      {
        Node * n = *i;
        for ( std::set<VolumeCell*>::const_iterator i=element_cells.begin(); i!=element_cells.end(); ++i )
        {
          VolumeCell * cell = *i;
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

  std::set<VolumeCell*> done;

  BuildDOFCellSets( point(), cells, nodal_cells[this], done );

  nodal_cells.erase( this );

  for ( std::map<Node*, std::set<VolumeCell*> >::iterator i=nodal_cells.begin();
        i!=nodal_cells.end();
        ++i )
  {
    Node * n = i->first;
    std::set<VolumeCell*> & cellset = i->second;
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
                                       const std::set<VolumeCell*> & cells,
                                       const std::set<VolumeCell*> & nodal_cells,
                                       std::set<VolumeCell*> & done )
{
  for ( std::set<VolumeCell*>::const_iterator i=nodal_cells.begin();
        i!=nodal_cells.end();
        ++i )
  {
    VolumeCell * cell = *i;
    if ( done.count( cell )==0 )
    {
      std::set<VolumeCell*> connected;
      std::set<Element*> elements;
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
    std::set<VolumeCell*> & cells = dofsets_[i];
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
