/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_graddesc.H
 *
<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANA_GRADDESC_H
#define STAT_INV_ANA_GRADDESC_H

#include "stat_inv_analysis.H"

namespace STR
{
namespace INVANA
{

class StatInvAnaGradDesc : public StatInvAnalysis
{

private:

  //! do the update on the parameter vector
  void UpdateOptStep();

  //! step size for the gradient scheme
  double stepsize_;

  //! maximum number of iterations
  int maxiter_;


public:

  //! standard constructor
  StatInvAnaGradDesc(Teuchos::RCP<DRT::Discretization> dis);

  ~StatInvAnaGradDesc() {};

  //! start the optimization routine
  void Optimize();

}; // class StatInvAnaGradDesc
}  // namespace INVANA
}  // namespace STR

#endif /*INV_ANALYSIS_GRADDESC_H*/
