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

  /// Input parameters for FSI
  namespace STR
  {

    /// Type of partitioned coupling for FSI problems
    enum InvAnalysisType
    {
      inv_none,          // no inverse analysis
      inv_lung,          // sophies lung inverse analysis
      inv_generalized    // generalized inverse analysis
    };

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_INVANALYSIS_H
#endif  // #ifdef CCADISCRET
