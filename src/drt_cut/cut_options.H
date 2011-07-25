#ifndef CUT_OPTIONS_H
#define CUT_OPTIONS_H

namespace GEO
{
namespace CUT
{

class Options
{
public:

  Options( bool positions=true,
           bool simpleshapes=true,
           bool genquad4=true,
           bool genhex8=true,
           bool genwedge6=false,
           bool genpyramid5=false )
    : positions_(positions),
      simpleshapes_(simpleshapes),
      genquad4_(genquad4),
      genhex8_(genhex8),
      genwedge6_(genwedge6),
      genpyramid5_(genpyramid5)
  {
  }

  void SetFindPositions( bool positions ) { positions_ = positions; }

  bool FindPositions() const { return positions_; }

  bool SimpleShapes() const { return simpleshapes_; }

  void SetSimpleShapes( bool simpleshapes ) { simpleshapes_ = simpleshapes; }

  bool GenQuad4() const { return genquad4_; }
  bool GenHex8() const { return genhex8_; }
  bool GenWedge6() const { return genwedge6_; }
  bool GenPyramid5() const { return genpyramid5_; }

private:

  bool positions_;
  bool simpleshapes_;
  bool genquad4_;
  bool genhex8_;
  bool genwedge6_;
  bool genpyramid5_;
};

}
}

#endif
