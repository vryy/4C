/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for Electrochemistry module.

\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_ELCH_DYN_H
#define BACI_ELCH_DYN_H

/*! entry point for the solution of electrochemistry problems */
void elch_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI electrochemistry-module logo on the screen */
void printlogo();


#endif  // ELCH_DYN_H
