/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for low-Mach-number flow module.

\level 2


*/
/*---------------------------------------------------------------------*/

#ifndef LOMA_DYN_H
#define LOMA_DYN_H

/*! entry point for the solution of low-Mach-number flow problems */
void loma_dyn(int restart /* do we have to perform a restart?  */
);


#endif /* #ifndef LOMA_DYN_H */
