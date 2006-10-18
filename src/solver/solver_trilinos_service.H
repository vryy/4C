#ifdef TRILINOS_PACKAGE
#ifndef SOLVER_TRILINOS_SERVICE_H
#define SOLVER_TRILINOS_SERVICE_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "../headers/standardtypes.h"
#include "../solver/solver.h"


#ifdef __cplusplus

/*----------------------------------------------------------------------*/
/*!
  \brief Simplify dstrc_enter and dstrc_exit calls in c++
 */
/*----------------------------------------------------------------------*/
class DSTraceHelper
{
public:
#ifdef DEBUG
  explicit DSTraceHelper(char* functname)
    { dstrc_enter(functname); }
  ~DSTraceHelper()
    { dstrc_exit(); }
#else
  explicit DSTraceHelper(char* functname)
    {}
#endif
};

#endif


void construct_trilinos_matrix(INTRA *actintra, TRILINOSMATRIX *trimatrix);
void construct_trilinos_diagonal_matrix(INTRA *actintra,
                                        TRILINOSMATRIX *trimatrix);
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
void close_nonquad_trilinos_matrix(TRILINOSMATRIX *A, TRILINOSMATRIX *cmat);
void exchange_coup_trilinos(
    PARTITION*      actpart,
    SOLVAR*         actsolv,
    INTRA*          actintra,
    TRILINOSMATRIX* tri
    );
void add_trilinos_value(struct _TRILINOSMATRIX *tri, DOUBLE v, INT row, INT col);
void add_trilinos_matrix(TRILINOSMATRIX* from, TRILINOSMATRIX* to, double factor);
void matvec_trilinos(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A);
  void matvec_trilinos_trans(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A);
void scale_trilinos_matrix(TRILINOSMATRIX* A, double factor);
void invert_trilinos_diagonal_matrix(TRILINOSMATRIX* A);
void mult_trilinos_mmm(TRILINOSMATRIX* dest,
                       TRILINOSMATRIX* A, INT transA,
                       TRILINOSMATRIX* B, INT transB,
                       TRILINOSMATRIX* C, INT transC);
void mult_trilinos_mmm_cont(TRILINOSMATRIX* dest,
                            TRILINOSMATRIX* A, INT transA,
                            TRILINOSMATRIX* B, INT transB,
                            TRILINOSMATRIX* C, INT transC);

#ifdef DEBUG
void extractGlobalRow(FILE* out,TRILINOSMATRIX* mat,INT row);
#endif

#ifdef __cplusplus
}
#endif


#endif /* SOLVER_TRILINOS_SERVICE_H */
#endif /* TRILINOS_PACKAGE */
