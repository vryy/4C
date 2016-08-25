/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_dyn.H

 \brief global access method for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_DYN_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_DYN_H_


/*! entry point for the solution of poro multiphase problems problems */
void poromultiphasescatra_dyn(
    int restart        /* do we have to perform a restart?  */
    );


#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_DYN_H_ */
