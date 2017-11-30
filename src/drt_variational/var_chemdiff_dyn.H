/*!----------------------------------------------------------------------
\file var_chemdiff_dyn.H
\brief Control routine for variational formulations

\level 2

<pre>
\maintainer Jorge De Anda Salazar
            DeAnda@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#ifndef VAR_CHEMDIFF_DYN_H
#define VAR_CHEMDIFF_DYN_H

/*! entry point for the solution of Variational formulations problems */
void var_chemdiff_dyn(
    int restart        /* do we have to perform a restart?  */
    );

/*! prints the BACI Variational formulations-module logo on the screen */
void print_var_logo();


#endif /* #ifndef VAR_CHEMDIFF_DYN_H */
