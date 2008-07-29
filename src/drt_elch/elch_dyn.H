/*!----------------------------------------------------------------------
\file elch_dyn.H
\brief Control routine for Electrochemistry module.


<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ELCH_DYN_H
#define ELCH_DYN_H

#ifdef __cplusplus
extern "C"
{
#endif
/* headers and prototypes of functions also accessible from c */

/*! entry point for the solution of electrochemistry problems */
void elch_dyn(
    int disnumff,      /* fluid discretization number, typically 0 */
    int disnumscatra,  /* scatra discretization number, typically 1 */
    int restart        /* do we have to perform a restart?  */
    );

/*! prints the BACI electrochemistry-module logo on the screen */
void printlogo();

#ifdef __cplusplus
}
#endif

#endif /* #ifndef ELCH_DYN_H */

#endif /* #ifdef CCADISCRET */
