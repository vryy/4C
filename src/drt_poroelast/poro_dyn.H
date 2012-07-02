/*!------------------------------------------------------------------------------------------------*
 \file poroelast.H

 \brief control routine of poroelasticity problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef POROELAST_H
#define POROELAST_H

#ifdef __cplusplus
extern "C"
{
#endif
  /* headers and prototypes of functions also accessible from c */

  void poroelast_drt();  // global control routine for standard porous media problem

  void poro_scatra_drt(); // global control routine for Poro + ScaTra problem.

  #ifdef __cplusplus
}
#endif

/*----------------------------------------------------------------------*/
#endif /* #ifndef POROELAST_H */
