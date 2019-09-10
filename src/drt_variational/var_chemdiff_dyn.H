/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for variational formulations

\level 2

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef VAR_CHEMDIFF_DYN_H
#define VAR_CHEMDIFF_DYN_H

/*! entry point for the solution of Variational formulations problems */
void var_chemdiff_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the BACI Variational formulations-module logo on the screen */
void print_var_logo();


#endif /* #ifndef VAR_CHEMDIFF_DYN_H */
