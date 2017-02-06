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
           bool genline2=true,
           bool genquad4=true,
           bool genhex8=true,
           bool genwedge6=false,
           bool genpyramid5=false )
    : nodal_dofset_strategy_( nodal_dofset_strategy ),
      positions_( positions ),
      simpleshapes_( simpleshapes ),
      genquad4_( genquad4 ),
      genhex8_( genhex8 ),
      genwedge6_( genwedge6 ),
      genpyramid5_( genpyramid5 ),
      gen_bcell_position_( INPAR::CUT::bcells_on_cut_side )
  {
  }

  void SetFindPositions( bool positions ) { positions_ = positions; }

  void SetNodalDofSetStrategy( INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy) { nodal_dofset_strategy_ = nodal_dofset_strategy; }

  INPAR::CUT::NodalDofSetStrategy GetNodalDofSetStrategy(){ return nodal_dofset_strategy_; }

  bool FindPositions() const { return positions_; }

  bool SimpleShapes() const { return simpleshapes_; }

  void SetSimpleShapes( bool simpleshapes ) { simpleshapes_ = simpleshapes; }

  /** \brief Set the position for the boundary cell creation
   *
   *  \author hiermeier \date 01/17 */
  void SetGenBoundaryCellPosition( INPAR::CUT::BoundaryCellPosition gen_bcell_position )
  { gen_bcell_position_ = gen_bcell_position; };

  enum INPAR::CUT::BoundaryCellPosition GenBoundaryCellPosition() const
  {
    return gen_bcell_position_;
  }

  /** \brief get the quad4 integration cell generation indicator
   *  (TRUE: QUAD4, FALSE: TRI3's are used) */
  bool GenQuad4() const { return genquad4_; }

  /** \brief get the hex8 integration cell generation indicator
   *  (TRUE: HEX8, FALSE: TRI4's are used) */
  bool GenHex8() const { return genhex8_; }

  /** \brief get the wedge6 integration cell generation indicator
   *  (TRUE: WEDGE6, FALSE: TRI4's are used) */
  bool GenWedge6() const { return genwedge6_; }

  /** \brief get the pyramid5 integration cell generation indicator
   *  (TRUE: PYRAMID5, FALSE: TRI4's are used) */
  bool GenPyramid5() const { return genpyramid5_; }

private:

  INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy_;

  bool positions_;

  bool simpleshapes_;

  /** \brief Indicator if QUAD4 integrations cells are going to be used or not
   *
   *  If TRUE a QUAD4 element is used as integration cell or otherwise the
   *  QUAD4 is splitted into 2 TRI3 elements. */
  bool genquad4_;

  /** \brief Indicator if HEX8 integrations cells are going to be used or not
   *
   *  If TRUE a HEX8 element is used as integration cell or otherwise the
   *  HEX8 is splitted into 5 TRI4 elements. */
  bool genhex8_;

  /** \brief Indicator if WEDGE6 integrations cells are going to be used or not
   *
   *  If TRUE a WEDGE6 element is used as integration cell or otherwise the
   *  WEDGE6 is splitted into 3 TRI4 elements. */
  bool genwedge6_;

  /** \brief Indicator if PYRAMID5 integrations cells are going to be used or not
   *
   *  If TRUE a PYRAMID5 element is used as integration cell or otherwise the
   *  PYRAMID5 is splitted into 2 TRI4 elements. */
  bool genpyramid5_;

  enum INPAR::CUT::BoundaryCellPosition gen_bcell_position_;
};

}
}

#endif
