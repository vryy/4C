/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport control algorithm
\level 1


*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_DYN_H
#define BACI_SCATRA_DYN_H

BACI_NAMESPACE_OPEN

/*! entry point for the solution of scalar transport problems */
void scatra_dyn(int restart /* do we have to perform a restart?  */
);

BACI_NAMESPACE_CLOSE

#endif  // SCATRA_DYN_H
