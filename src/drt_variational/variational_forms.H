/*!----------------------------------------------------------------------
\file variational_forms.H
\brief Control routine for variational formulations

\level 2

<pre>
\maintainer Jorge De Anda Salazar
            DeAnda@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*----------------------------------------------------------------------*/

#ifndef VARIATIONAL_FORMS_H
#define VARIATIONAL_FORMS_H

/*! entry point for the solution of Variational formulations problems */
void variational_forms(
    int restart        /* do we have to perform a restart?  */
    );

/*! prints the BACI Variational formulations-module logo on the screen */
void printlogo();


#endif /* #ifndef VARIATIONAL_FORMS_H */
