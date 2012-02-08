/*----------------------------------------------------------------------*/
/*!
\file inpar_statmech.H

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_MLMC_H
#define INPAR_MLMC_H

#include "inpar_parameterlist_utils.H"
/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in multi level monte carlo
  namespace MLMC
  {
  enum CorrStruct
      {
        corr_gaussian, ///< gaussian correlation structure
        // more to follow
      };
  enum MarginalPdf
        {
          pdf_gaussian,// normal distribution
          pdf_beta,    // beta distribution
          pdf_lognormal,// lognormal distribution
          // more to follow
        };
    /// To be filled with input params
  }  // namespace MLMC

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_MLMC_H
#endif  // #ifdef CCADISCRET
