/*!----------------------------------------------------------------------
\file xfluid_dyn_nln_drt.H
\brief Control routine for fluid (time) integration. Calls

     o Single step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step

     o Genalpha time integartion

     o Stationary solver


<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef XFLUID_DYN_NLN_DRT_H
#define XFLUID_DYN_NLN_DRT_H


/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus
extern "C"
{
#endif
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


/*!
 * Main control routine for fluid including various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha
 *        o stationary
 *
 */
void xdyn_fluid_drt();




/* BEGIN   C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/
#ifdef __cplusplus
}
#endif
/* END     C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++C++*/


#endif  /* #ifndef XFLUID_DYN_NLN_DRT_H */
#endif  /* #ifdef CCADISCRET */
