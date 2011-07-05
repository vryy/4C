#ifndef CUT_VOLUMECELLGENERATOR_H
#define CUT_VOLUMECELLGENERATOR_H

#include <algorithm>
#include <vector>
#include <map>
#include <list>

#include "cut_utils.H"
#include "cut_point.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Facet;
  class Side;
  class Mesh;
  class Element;
  class VolumeCell;

class VolumeCellGenerator
{
public:

  VolumeCellGenerator( const std::vector<Side*> & sides, const plain_facet_set & facets );

  void Print()
  {
  }

  void CreateVolumeCells( Mesh & mesh, Element * element, plain_volumecell_set & cells );

private:

  bool Done();

  struct GenPoint;
  struct GenLine;
  struct GenFacet;
  struct GenVolume;

#ifdef CUT_USE_SORTED_VECTOR
  typedef sorted_vector<GenPoint*> GenPointSet;
  typedef sorted_vector<GenLine*> GenLineSet;
  typedef sorted_vector<GenFacet*> GenFacetSet;
  typedef sorted_vector<GenVolume*> GenVolumeSet;
#else
  typedef std::set<GenPoint*> GenPointSet;
  typedef std::set<GenLine*> GenLineSet;
  typedef std::set<GenFacet*> GenFacetSet;
  typedef std::set<GenVolume*> GenVolumeSet;
#endif

  struct GenPoint
  {
    GenPoint()
      : p_( NULL )
    {
    }

    Point * p_;

    void AddLine( GenLine * l )
    {
      lines_.insert( l );
    }

    void RemoveLine( GenLine * l )
    {
      lines_.erase( l );
    }

    bool TestCorner( std::list<GenVolume> & volumes );

  private:
    GenLineSet lines_;
  };

  struct GenLine
  {
    GenLine()
      : p1_( NULL ),
        p2_( NULL )
    {
    }

    void SetPoints( GenPoint * p1, GenPoint * p2 )
    {
      if ( p1_ != NULL )
        p1_->RemoveLine( this );
      if ( p2_ != NULL )
        p2_->RemoveLine( this );
      p1_ = p1;
      p2_ = p2;
      p1_->AddLine( this );
      p2_->AddLine( this );
    }

    void AddFacet( GenFacet * facet )
    {
      facets_.insert( facet );
    }

    void RemoveFacet( GenFacet * facet )
    {
      facets_.erase( facet );
      if ( facets_.size()==0 )
      {
        p1_->RemoveLine( this );
        p2_->RemoveLine( this );
        //p1_ = NULL;
        //p2_ = NULL;
      }
    }

    void CommonFacets( GenLine * other, GenFacetSet & facets );

    GenFacet * CommonFacet( GenLine * other );

    GenPoint * P1() const { return p1_; }
    GenPoint * P2() const { return p2_; }

    const GenFacetSet & Facets() const { return facets_; }

  private:
    GenPoint * p1_;
    GenPoint * p2_;
    GenFacetSet facets_;
  };

  struct GenFacet
  {
    GenFacet()
      : facet_( NULL ),
        openvolumes_( 0 )
    {
    }

    void AddLine( GenLine * l )
    {
      lines_.insert( l );
      l->AddFacet( this );
    }

    void AddVolume( GenVolume * v )
    {
      if ( Done() )
        throw std::runtime_error( "cannot add any more volume" );
      volumes_.insert( v );
      openvolumes_ -= 1;
      if ( Done() )
      {
        for ( GenLineSet::iterator i=lines_.begin(); i!=lines_.end(); ++i )
        {
          GenLine * gl = *i;
          gl->RemoveFacet( this );
        }
      }
    }

    bool Done() { return openvolumes_==0; }

    const GenLineSet & Lines() const { return lines_; }

    Facet * facet_;
    int openvolumes_;

  private:
    GenLineSet lines_;
    GenVolumeSet volumes_;
  };

  struct GenVolume
  {
    void AddFacet( GenFacet * f )
    {
      if ( facets_.count( f )==0 )
      {
        facets_.insert( f );
        f->AddVolume( this );
      }
    }

    bool HasFacet( GenFacet * f )
    {
      return facets_.count( f ) > 0;
    }

    void Finalize();

    const GenFacetSet & Facets() const { return facets_; }

  private:
    GenFacetSet facets_;
  };

  struct GenVolumeFinalizer
  {
    explicit GenVolumeFinalizer( GenVolume & volume )
      : volume_( volume )
    {
    }

    void AddFacet( GenFacet * f );

    void CloseVolume();

    GenFacet * CommonFacet( GenPoint * p );

  private:
    GenVolume & volume_;
    std::map<GenPoint*, GenLineSet> points_;
    std::map<GenLine*, GenFacetSet> lines_;
    GenLineSet openlines_;
    std::map<GenPoint*, GenLineSet> openpoints_;
  };

  std::map<Point*, GenPoint> points_;
  std::vector<GenLine> lines_;
  std::map<Facet*, GenFacet> facets_;
  std::list<GenVolume> volumes_;
};

}
}


#endif
