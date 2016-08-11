/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_dyn.H

 \brief entry point (global control routine) for poroelasticity multiphase flow

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_DYN_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_DYN_H_

/*! entry point for the solution of Lubrication problems */
void poromultiphase_dyn(
    int restart        /* do we have to perform a restart?  */
    );



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_DYN_H_ */
