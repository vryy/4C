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

  //! global control routine for standard porous media problem
  void poroelast_drt();

  //! global control routine for Poro + ScaTra problem.
  void poro_scatra_drt();


/*----------------------------------------------------------------------*/
#endif /* #ifndef POROELAST_H */
