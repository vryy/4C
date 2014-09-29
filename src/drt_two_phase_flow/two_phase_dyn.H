/*!----------------------------------------------------------------------
\file two_phase_dyn.H
\brief Control routine for two phase flow module.


<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915236
</pre>

*----------------------------------------------------------------------*/

#ifndef TWO_PHASE_DYN_H
#define TWO_PHASE_DYN_H

/*! entry point for the solution of two phase flow problems */
void two_phase_dyn(
    int restart        /* do we have to perform a restart?  */
    );

/*! entry point for the solution of two phase flow problems utilizing XFEM */
void fluid_xfem_ls_drt(
    );

#endif /* #ifndef TWO_PHASE_DYN_H */

