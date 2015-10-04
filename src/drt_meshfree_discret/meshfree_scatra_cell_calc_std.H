/*--------------------------------------------------------------------------*/
/*!
\file meshfree_scatra_cell_calc_std.H

\brief Internal standard implementation of meshfree scalar transport cells

  <pre>
  Maintainer: Keijo Nissen
  nissen@lnm.mw.tum.de
  http://www.lnm.mw.tum.de
  089 - 289-15253
  </pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef MESHFREE_SCATRA_Cell_CALC_STD_H
#define MESHFREE_SCATRA_Cell_CALC_STD_H

#include "meshfree_scatra_cell_calc.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class MeshfreeScaTraCellCalcStd : public MeshfreeScaTraCellCalc<distype>
{
  /// private constructor, since we are a Singleton.
  MeshfreeScaTraCellCalcStd(const int numdofpernode,const int numscal,const std::string& disname);

  typedef MeshfreeScaTraCellCalc<distype> my;

public:
  virtual ~MeshfreeScaTraCellCalcStd() {}

  /// Singleton access method
  static MeshfreeScaTraCellCalcStd<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const MeshfreeScaTraCellCalcStd* delete_me = NULL );

  /// called upon destruction
  virtual void Done();

private:

};
}
}

#endif
