/*----------------------------------------------------------------------*/
/*!
\file inpar_combust.H

\brief Input parameters for combustion

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef INPAR_COMBUST_H
#define INPAR_COMBUST_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  namespace COMBUST
  {
    /// type of combustion problem
    enum CombustionType
    {
      combusttype_premixedcombustion,
      combusttype_twophaseflow
    };

    /// type of reinitialization for level set function (G-function)
    enum ReInitialActionGfunc
    {
      reinitaction_none,
      reinitaction_byfunction,
      reinitaction_signeddistancefunction
    };
  }  // namespace COMBUST

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_COMBUST_H
#endif  // #ifdef CCADISCRET

