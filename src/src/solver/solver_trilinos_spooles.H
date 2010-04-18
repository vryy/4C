#ifndef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef SOLVER_TRILINOS_SPOOLES_H
#define SOLVER_TRILINOS_SPOOLES_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "../headers/standardtypes.h"
#include "../solver/solver.h"


#ifdef __cplusplus
}
#endif

void solver_trilinos_spooles(
                             struct _SOLVAR*         actsolv,
                             struct _INTRA*          actintra,
                             struct _TRILINOSMATRIX* tri,
                             struct _SPOOLMAT*       spo,
                             struct _DIST_VECTOR*    sol,
                             struct _DIST_VECTOR*    rhs,
                             int                     option
                             );


#endif /* SOLVER_TRILINOS_SPOOLES_H */
#endif /* TRILINOS_PACKAGE */
#endif
