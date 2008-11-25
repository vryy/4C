/*----------------------------------------------------------------------*/
/*!
\file inpar_combust.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_COMBUST_H_
#define INPAR_COMBUST_H_

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for combustion
  namespace COMBUST
  {
    /// Type of reinitialization for level set function (G-function)
    enum ReInitialActionGfunc
    {
      reinitialize_by_function,
      compute_signeddistancefunction
    };
  }  // namespace COMBUST

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_COMBUST_H
#endif  // #ifdef CCADISCRET

