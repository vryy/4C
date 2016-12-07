/*---------------------------------------------------------------------*/
/*!
\file cut_integrationcellcreator.H

\brief Create and handle integrationcells for the Tessellation routine

\level 3

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_INTEGRATIONCELLCREATOR_H
#define CUT_INTEGRATIONCELLCREATOR_H

#include "cut_volumecell.H"

namespace GEO
{
namespace CUT
{
  class Mesh;

  /*!
  \brief Special cases for integration cell creation: a library of cuts.

  Some cuts lead to known shapes those meshes are known. Here we collect
  these cases. Feel free to add your own.
 */
class IntegrationCellCreator
{
public:

  static bool CreateCells( Mesh & mesh, Element * element, const plain_volumecell_set & cells );

  static bool CreateCell( Mesh & mesh, DRT::Element::DiscretizationType shape, VolumeCell * cell );

private:

  void Execute( Mesh & mesh );

  bool CreateTet4Cell( Mesh & mesh, VolumeCell * cell, const plain_facet_set & facets );
  bool CreateHex8Cell( Mesh & mesh, VolumeCell * cell, const plain_facet_set & facets );
  bool CreateWedge6Cell( Mesh & mesh, VolumeCell * cell, const plain_facet_set & facets );
  bool CreatePyramid5Cell( Mesh & mesh, VolumeCell * cell, const plain_facet_set & facets );
  bool CreateSpecialCases( Mesh & mesh, VolumeCell * cell, const plain_facet_set & facets );

  bool Hex8HorizontalCut( Mesh & mesh,
                          Element * element,
                          VolumeCell * cell,
                          const plain_facet_set & facets,
                          int axis,
                          double r );

  void Add( VolumeCell * vc,
            DRT::Element::DiscretizationType shape,
            const std::vector<Point*> & points )
  {
    volume & v = cells_[vc];
    std::vector<ic> & cells = v.domain_;
    cells.push_back( ic() );
    ic & cell = cells.back();
    cell.shape_ = shape;
    cell.points_.reserve( points.size() );
    cell.points_.assign( points.begin(), points.end() );
  }

  void AddSide( VolumeCell * vc,
                Facet * facet,
                DRT::Element::DiscretizationType shape,
                const std::vector<Point*> & side )
  {
    volume & cell = cells_[vc];
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

    void Execute( Mesh & mesh, VolumeCell * vc )
    {
      vc->NewIntegrationCell( mesh, shape_, points_ );
    }
  };

  struct volume
  {
    std::vector<ic> domain_;
    std::vector<bc> boundary_;

    void Execute( Mesh & mesh, VolumeCell * vc )
    {
      for ( std::vector<ic>::iterator i=domain_.begin(); i!=domain_.end(); ++i )
      {
        ic & cell = *i;
        cell.Execute( mesh, vc );
      }
      for ( std::vector<bc>::iterator i=boundary_.begin(); i!=boundary_.end(); ++i )
      {
        bc & bcell = *i;
        bcell.Execute( mesh, vc );
      }
    }
  };

  std::map<VolumeCell*, volume> cells_;
};

}
}

#endif
