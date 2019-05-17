/*---------------------------------------------------------------------*/
/*!

\brief Control routine for variational formulations

\level 2

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef VARIATIONAL_FORMS_H
#define VARIATIONAL_FORMS_H

/*! entry point for the solution of Variational formulations problems */
void variational_forms(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI Variational formulations-module logo on the screen */
void printlogo();


#endif /* #ifndef VARIATIONAL_FORMS_H */
