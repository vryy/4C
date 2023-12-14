/*----------------------------------------------------------------------*/
/*! \file
\brief entry point for level-set transport problems
\level 2

*----------------------------------------------------------------------*/

#ifndef BACI_LEVELSET_DYN_H
#define BACI_LEVELSET_DYN_H

BACI_NAMESPACE_OPEN

/*! entry point for the solution of level-set problems */
void levelset_dyn(int restart /* do we have to perform a restart?  */
);

BACI_NAMESPACE_CLOSE

#endif  // LEVELSET_DYN_H
