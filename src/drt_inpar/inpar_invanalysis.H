/*----------------------------------------------------------------------*/
/*!
\file inpar_invanalysis.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_INVANALYSIS_H
#define INPAR_INVANALYSIS_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for Inverse Analysis
  namespace STR
  {

    /// Type of Inverse Analysis
    enum InvAnalysisType
    {
      inv_none,          // no inverse analysis
      inv_lung,          // sophies lung inverse analysis
      inv_coll,          // sophie's collagenase inverse analysis
      inv_airways,       // sophies airways inverse analysis
      inv_generalized    // generalized inverse analysis
    };



    /// Type of experiment
    enum InvAnalysisExpType
    {
      inv_exp_none,      // no treatment
      inv_exp_col,       // collagenase treated
      inv_exp_ela,       // elastase trested
      inv_exp_colfull,   // first trested with collagenase and
                         // afterwards with elastase
      inv_exp_elafull,   // first trested with elastase and
                         // afterwards with collagenase
      inv_exp_full       // treated with both right away
    };

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_INVANALYSIS_H
#endif  // #ifdef CCADISCRET
