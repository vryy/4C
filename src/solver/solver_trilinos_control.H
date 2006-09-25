#ifdef TRILINOS_PACKAGE
#ifndef SOLVER_TRILINOS_CONTROL_H
#define SOLVER_TRILINOS_CONTROL_H

#ifdef __cplusplus
extern "C" 
{
#endif

#include "../headers/standardtypes.h"
#include "../solver/solver.h"

void solver_trilinos_control(struct _FIELD          *actfield,
                             int                     disnum,
                             struct _SOLVAR         *actsolv,
                             struct _INTRA          *actintra,
                             enum   _SPARSE_TYP     *sysarray_typ,
                             union  _SPARSE_ARRAY   *sysarray,
                             struct _DIST_VECTOR    *sol,
                             struct _DIST_VECTOR    *rhs,
                             INT                     option);

#ifdef __cplusplus
}
#endif



#endif /* SOLVER_TRILINOS_CONTROL_H */
#endif /* TRILINOS_PACKAGE */
