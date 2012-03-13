/*!----------------------------------------------------------------------
\file scatra_dyn.H
\brief

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_DYN_H
#define SCATRA_DYN_H

/*! entry point for the solution of scalar transport problems */
void scatra_dyn(
    int disnumff,      /* fluid discretization number, typically 0 */
    int disnumscatra,  /* scatra discretization number, typically 1 */
    int restart        /* do we have to perform a restart?  */
    );

#endif  /* #ifndef SCATRA_DYN_H */
