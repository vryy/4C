#ifdef TRILINOS_PACKAGE
#ifndef SOLVER_TRILINOS_SERVICE_H
#define SOLVER_TRILINOS_SERVICE_H

#ifdef __cplusplus
extern "C" 
{
#endif

#include "../headers/standardtypes.h"
#include "../solver/solver.h"

void init_trilinos_matrix(FIELD *actfield,PARTITION *actpart,SOLVAR *actsolv,
                          INTRA *actintra,TRILINOSMATRIX *trimatrix,int disnum);
  
void trilinos_cp_matrixmask(TRILINOSMATRIX  *from, TRILINOSMATRIX* to);

void trilinos_zero_matrix(TRILINOSMATRIX *tri);

void add_trilinos(
    struct _PARTITION       *actpart,
    struct _SOLVAR          *actsolv,
    struct _INTRA           *actintra,
    struct _ELEMENT         *actele,
    struct _TRILINOSMATRIX  *tri1,
    struct _TRILINOSMATRIX  *tri2,
    struct _ARRAY           *elearray1,
    struct _ARRAY           *elearray2
    );

void close_trilinos_matrix(struct _TRILINOSMATRIX *tri);
void exchange_coup_trilinos(
    PARTITION*      actpart,
    SOLVAR*         actsolv,
    INTRA*          actintra,
    TRILINOSMATRIX* tri
    );
void add_trilinos_matrix(TRILINOSMATRIX* from, TRILINOSMATRIX* to, double factor);
void matvec_trilinos(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A);
void scale_trilinos_matrix(TRILINOSMATRIX* A, double factor);

#ifdef __cplusplus
}
#endif


#endif /* SOLVER_TRILINOS_SERVICE_H */
#endif /* TRILINOS_PACKAGE */
