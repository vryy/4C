#ifndef DRT_UTILS_GAUSSPOINTS_H
#define DRT_UTILS_GAUSSPOINTS_H

#include <iostream>

#include <vector>

#include "../../src/drt_lib/drt_element.H"

namespace DRT
{
namespace UTILS
{

/// base class of gauss point collection
class GaussPoints
{
public:

  virtual ~GaussPoints() {}

  /// number of gauss points
  virtual int NumPoints() const = 0;

  /// spacial dimension
  virtual int NumDimension() const = 0;

  /// gauss point coordinates
  virtual const double * Point( int point ) const = 0;

  /// gauss weight
  virtual double Weight( int point ) const = 0;

  /// debug print
  virtual void Print() const = 0;
};


/// specific collected gauss points for xfem usage
class CollectedGaussPoints : public GaussPoints
{
public:

  CollectedGaussPoints( int size=0 )
  {
    gp_.reserve( size );
  }

  void Append( double x, double y, double z, double w )
  {
    gp_.push_back( Entry( x,y,z,w ) );
  }

  virtual int NumPoints() const { return gp_.size(); }

  virtual int NumDimension() const { return 3; }

  virtual const double * Point( int point ) const { return &gp_[point].data[0]; }

  virtual double Weight( int point ) const { return gp_[point].data[3]; }

  virtual void Print() const
  {
    std::cout << " collected gauss points:\n";
    for ( int i=0; i<NumPoints(); ++i )
    {
      std::cout << "    ";
      for ( int j=0; j<NumDimension(); ++j )
        std::cout << gp_[i].data[j] << " ";
      std::cout << gp_[i].data[3] << "\n";
    }
  }

private:

  struct Entry
  {
    Entry( double x, double y, double z, double w )
    {
      data[0] = x;
      data[1] = y;
      data[2] = z;
      data[3] = w;
    }
    double data[4];
  };
  std::vector<Entry> gp_;
};


/// gauss integration interface
class GaussIntegration
{
public:

  /// internal iterator
  class GaussPointIterator
  {
  public:

    GaussPointIterator( GaussPoints * gp, int point )
      : gp_( gp ),
        point_( point )
    {}

    void operator++()
    {
      point_ += 1;
    }

    std::pair<const double*, double> operator*()
    {
      return std::make_pair( gp_->Point( point_ ), gp_->Weight( point_ ) );
    }

    bool operator==( const GaussPointIterator & other ) const
    {
      return gp_==other.gp_ and point_==other.point_;
    }

    bool operator!=( const GaussPointIterator & other ) const
    {
      return gp_!=other.gp_ or point_!=other.point_;
    }

  private:
    GaussPoints * gp_;
    int point_;
  };

  typedef GaussPointIterator iterator;
  typedef const GaussPointIterator const_iterator;

  GaussIntegration( DRT::Element::DiscretizationType distype );

  explicit GaussIntegration( Teuchos::RCP<GaussPoints> gp )
    : gp_( gp )
  {
  }

  GaussPointIterator begin() { return GaussPointIterator( &*gp_, 0 ); }

  GaussPointIterator end() { return GaussPointIterator( &*gp_, gp_->NumPoints() ); }

  /// number of gauss points
  int NumPoints() const { return gp_->NumPoints(); }

  /// spacial dimension
  int NumDimension() const { return gp_->NumDimension(); }

  /// gauss point coordinates
  const double * Point( int point ) const { return gp_->Point( point ); }

  /// gauss weight
  double Weight( int point ) const { return gp_->Weight( point ); }

  /// debug print
  void Print() const { gp_->Print(); }

private:

  /// internal collection
  Teuchos::RCP<GaussPoints> gp_;
};

}
}

#endif
