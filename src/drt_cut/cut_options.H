/*---------------------------------------------------------------------*/
/*!
\file cut_options.H

\brief options to set up the cut

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_OPTIONS_H
#define CUT_OPTIONS_H

#include "../drt_inpar/inpar_cut.H"

namespace GEO
{
namespace CUT
{

/*!
 \brief Options defined to configure the cutting behavior
 */
class Options
{
public:

  Options( INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy = INPAR::CUT::NDS_Strategy_full,
           bool positions=true,
           bool simpleshapes=true,
           bool genquad4=true,
           bool genhex8=true,
           bool genwedge6=false,
           bool genpyramid5=false )
    : nodal_dofset_strategy_(nodal_dofset_strategy),
      positions_(positions),
      simpleshapes_(simpleshapes),
      genquad4_(genquad4),
      genhex8_(genhex8),
      genwedge6_(genwedge6),
      genpyramid5_(genpyramid5)
  {
  }

  void SetFindPositions( bool positions ) { positions_ = positions; }

  void SetNodalDofSetStrategy( INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy) { nodal_dofset_strategy_ = nodal_dofset_strategy; }

  INPAR::CUT::NodalDofSetStrategy GetNodalDofSetStrategy(){ return nodal_dofset_strategy_; }

  bool FindPositions() const { return positions_; }

  bool SimpleShapes() const { return simpleshapes_; }

  void SetSimpleShapes( bool simpleshapes ) { simpleshapes_ = simpleshapes; }

  bool GenQuad4() const { return genquad4_; }
  bool GenHex8() const { return genhex8_; }
  bool GenWedge6() const { return genwedge6_; }
  bool GenPyramid5() const { return genpyramid5_; }

private:

  INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy_;

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
