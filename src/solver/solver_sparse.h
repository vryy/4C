/*!
\file
\brief a simple a general parallel sparse matrix

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This attempt stems from the need to have efficient solver independent
global sparse matrices for Gresho's version of the projection
method. It's just a sparse matrix, there is no dependency on any
fields or other variables.

*/

#ifndef SOLVER_SPARSE_H
#define SOLVER_SPARSE_H

#include "../headers/standardtypes.h"
#include "../pss_full/pss_set.h"


/*----------------------------------------------------------------------*/
/*!
  \brief support structures used to build up a sparse mask

  To create a sparse mask we need to gather the information. We do so
  here using sorted fields.

  The current implementation uses dynamic arrays. Therefore we need to
  copy a lot. Lists, however, are annoying to search. A better way
  would be some tree structure. But right now this seems to be too
  much work. However, it could be done and the change would be local
  to this structure. :)

  \author u.kue
  \date 12/05
*/
/*----------------------------------------------------------------------*/
typedef struct _SPARSE_MASK_LIST
{
  INT bandwidth;
  INTSET* columns;

  INT rows;
  INT cols;			/* columns in this matrix */
} SPARSE_MASK_LIST;


/*----------------------------------------------------------------------*/
/*!
  \brief Customized sparse matrix for the pressure gradient

  To be independent of any solvers we use your own structure
  here. It's a compressed column based sparse matrix with sorted row
  entries, the umfpack style. Here, however, we use local dofs numbers
  in parallel execution.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
typedef struct _SPARSE
{
  INT rows;
  INT cols;			/* columns in this matrix */
  INT* Ap;			/* column indices */
  INT* Ai;			/* row numbers */
  DOUBLE* Ax;			/* matrix entries */
} SPARSE;


/*----------------------------------------------------------------------*/
/*!
  \brief Slicing in parallel execution

  When executed parallel each processor owns a slice of the matrix. We
  can slice horizontally or vertically.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
typedef enum _SLICING_DIRECTION
{
  sd_slice_horizontal,
  sd_slice_vertical
} SLICING_DIRECTION;


/*----------------------------------------------------------------------*/
/*!
  \brief Parallel version of the compressed column matrix.

  Each processor holds just a slice (some colunms or some rows)

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
typedef struct _PARALLEL_SPARSE
{
  INT total_cols;

#ifdef PARALLEL

  SLICING_DIRECTION slicing;

  /* list of dofs updated on this proc */
  /* global (column or row) ids of local columns */
  /* The user has to fill this array. This way we keep the sparse
   * matrix independent of the ccarat mesh structures. */
  INT* update;

#endif

  SPARSE slice;
} PARALLEL_SPARSE;


/* Prototypes */

void sparse_mask_list_init(SPARSE_MASK_LIST* s, INT rows, INT cols, INT bandwidth);
void sparse_mask_list_destroy(SPARSE_MASK_LIST* s);
void sparse_mask_list_mark(SPARSE_MASK_LIST* s, INT row, INT col);
#ifdef DEBUG
void sparse_mask_list_dump(SPARSE_MASK_LIST* s, char* filename);
#endif
void sparse_mask_list_mask_solver(SPARSE_MASK_LIST* s, SOLVAR *solv);

void parallel_sparse_init(PARALLEL_SPARSE* ps, INT rows, INT cols,
                          SLICING_DIRECTION slicing, INT local_cols);
void parallel_sparse_destroy(PARALLEL_SPARSE* ps);
void sparse_init(SPARSE* s, INT rows, INT cols);
void sparse_destroy(SPARSE* s);
void sparse_fix_mask(SPARSE* s, SPARSE_MASK_LIST* sml);
void sparse_zero(SPARSE* s);
DOUBLE* sparse_entry(SPARSE* s, INT row, INT col);
INT sparse_has_entry(SPARSE* s, INT row, INT col);

#ifdef DEBUG
void sparse_debugdump(SPARSE* s, INT* update, INT update_length, FILE* f, CHAR* name);
#endif


/* Special support for the projection method */
#ifdef D_FLUID_PM
void parallel_sparse_pm_matmat(PARALLEL_SPARSE* pmat, PARALLEL_SPARSE* grad,
                               DOUBLE* lmass, INTRA* actintra);
#ifdef AZTEC_PACKAGE
void parallel_sparse_convert_aztec(PARALLEL_SPARSE* ps, SOLVAR* solvar);
void parallel_sparse_copy_aztec(PARALLEL_SPARSE* ps, SOLVAR* solvar);
#endif
#endif

#endif
