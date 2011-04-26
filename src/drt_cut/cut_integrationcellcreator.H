#ifndef CUT_INTEGRATIONCELLCREATOR_H
#define CUT_INTEGRATIONCELLCREATOR_H

#include "cut_volumecell.H"

namespace GEO
{
namespace CUT
{
  class Mesh;

class IntegrationCellCreator
{
public:

  static bool CreateCells( Mesh & mesh, const std::set<VolumeCell*> & cells );

  static bool CreateCell( Mesh & mesh, DRT::Element::DiscretizationType shape, VolumeCell * cell );

private:

  void Execute( Mesh & mesh );

  bool CreateTet4Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets );
  bool CreateHex8Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets );
  bool CreateWedge6Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets );
  bool CreatePyramid5Cell( Mesh & mesh, VolumeCell * cell, const std::set<Facet*> & facets );

  void Add( VolumeCell * vc,
            DRT::Element::DiscretizationType shape,
            const std::vector<Point*> & points )
  {
    ic & cell = cells_[vc];
    cell.shape_ = shape;
    cell.points_.reserve( points.size() );
    cell.points_.assign( points.begin(), points.end() );
  }

  void AddSide( VolumeCell * vc,
                Facet * facet,
                DRT::Element::DiscretizationType shape,
                const std::vector<Point*> & side )
  {
    ic & cell = cells_[vc];
    std::vector<bc> & bcells = cell.boundary_;
    bcells.push_back( bc() );
    bc & bcell = bcells.back();
    bcell.shape_ = shape;
    bcell.facet_ = facet;
    bcell.side_.reserve( side.size() );
    bcell.side_.assign( side.begin(), side.end() );
  }

  struct bc
  {
    DRT::Element::DiscretizationType shape_;
    std::vector<Point*> side_;
    Facet * facet_;

    void Execute( Mesh & mesh, VolumeCell * vc )
    {
      vc->NewBoundaryCell( mesh, shape_, facet_, side_ );
    }
  };

  struct ic
  {
    DRT::Element::DiscretizationType shape_;
    std::vector<Point*> points_;
    //VolumeCell * cell_;
    std::vector<bc> boundary_;

    void Execute( Mesh & mesh, VolumeCell * vc )
    {
      vc->NewIntegrationCell( mesh, shape_, points_ );
    }
  };

  std::map<VolumeCell*, ic> cells_;
};

}
}

#endif
