/*!----------------------------------------------------------------------
\file tsi_dyn.H
\brief Control routine for thermo-structure-interaction problems.


<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TSI_DYN_H
#define TSI_DYN_H

#ifdef __cplusplus
extern "C"
{
#endif
/* headers and prototypes of functions also accessible from c */

/* entry point for the solution of thermo-structure-interaction problems */
void tsi_dyn_drt();

/*! prints the BACI tsi-logo on the screen */
/*void printlogo();*/

#ifdef __cplusplus
}
#endif

/*----------------------------------------------------------------------*/
#endif /* #ifndef TSI_DYN_H */
#endif /* #ifdef CCADISCRET */
