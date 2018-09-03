/*!----------------------------------------------------------------------
\file scatra_cardiac_monodomain_dyn.H

\brief entry point for cardiac monodomain scalar transport problems

\level 2

<pre>
\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15264
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_CARDIAC_MONODOMAIN_DYN_H
#define SCATRA_CARDIAC_MONODOMAIN_DYN_H

/*! entry point for the solution of electrochemistry problems */
void scatra_cardiac_monodomain_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI cardiac monodomain-module logo on the screen */
void printheartlogo();


#endif /* #ifndef SCATRA_CARDIAC_MONODOMAIN_DYN_H */
