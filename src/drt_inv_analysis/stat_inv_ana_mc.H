/*----------------------------------------------------------------------*/
/*!
 * \file stat_inv_ana_mc.H
 *
<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef STAT_INV_ANA_MC_H
#define STAT_INV_ANA_MC_H

#include "stat_inv_analysis.H"

namespace STR
{
namespace INVANA
{

class StatInvAnaMC : public StatInvAnalysis
{

private:

  //! do the update on the parameter vector
  void UpdateOptStep(Epetra_MultiVector* objgrad, int nstep);

public:

  //! standard constructor
  StatInvAnaMC(Teuchos::RCP<DRT::Discretization> dis);

  //! start the optimization routine
  void Optimize();

}; // class StatInvAnaGradDesc
}  // namespace INVANA
}  // namespace STR

#endif /*INV_ANALYSIS_MC_H*/
