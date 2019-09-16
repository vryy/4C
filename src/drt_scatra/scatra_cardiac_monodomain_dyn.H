/*----------------------------------------------------------------------*/
/*! \file

\brief entry point for cardiac monodomain scalar transport problems

\level 2

\maintainer Anh-Tu Vuong

*----------------------------------------------------------------------*/

#ifndef SCATRA_CARDIAC_MONODOMAIN_DYN_H
#define SCATRA_CARDIAC_MONODOMAIN_DYN_H

/*! entry point for the solution of electrochemistry problems */
void scatra_cardiac_monodomain_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI cardiac monodomain-module logo on the screen */
void printheartlogo();


#endif /* #ifndef SCATRA_CARDIAC_MONODOMAIN_DYN_H */
