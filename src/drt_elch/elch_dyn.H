/*!----------------------------------------------------------------------
\file elch_dyn.H
\brief Control routine for Electrochemistry module.

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#ifndef ELCH_DYN_H
#define ELCH_DYN_H

/*! entry point for the solution of electrochemistry problems */
void elch_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI electrochemistry-module logo on the screen */
void printlogo();


#endif /* #ifndef ELCH_DYN_H */
