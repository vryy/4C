/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for poroelasticity multiphase flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROMULTIPHASE_DYN_H
#define BACI_POROMULTIPHASE_DYN_H

BACI_NAMESPACE_OPEN

/*! entry point for the solution of poro multiphase problems problems */
void poromultiphase_dyn(int restart /* do we have to perform a restart?  */
);



BACI_NAMESPACE_CLOSE

#endif  // POROMULTIPHASE_DYN_H
