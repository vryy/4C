/*!
 * \file xfluid_dyn_nln_drt.H
 *
 * contains control routine for fluid (time) integration
 * 
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 * 
 */
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
 *     o Single step one-step-theta time integration
 *     
 *     o Two step BDF2 Gear's methode with one-step-theta start step
 *     
 *     o Genalpha time integartion
 *     
 *     o Stationary solver
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
