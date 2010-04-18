/*!----------------------------------------------------------------------
\file loma_dyn.H
\brief Control routine for low-Mach-number flow module.


<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef LOMA_DYN_H
#define LOMA_DYN_H

#ifdef __cplusplus
extern "C"
{
#endif
/* headers and prototypes of functions also accessible from c */

/*! entry point for the solution of low-Mach-number flow problems */
void loma_dyn(
    int disnumff,      /* fluid discretization number, typically 0 */
    int disnumscatra,  /* scatra discretization number, typically 1 */
    int restart        /* do we have to perform a restart?  */
    );

#ifdef __cplusplus
}
#endif

#endif /* #ifndef LOMA_DYN_H */

#endif /* #ifdef CCADISCRET */
