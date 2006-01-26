/*!
\file
\brief a simple a general parallel sparse matrix

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
*/

#include "solver_sparse.h"


/*----------------------------------------------------------------------*/
/*!
  \brief set up an empty parallel sparse matrix

  In the parallel case we store just some rows or some columns on each
  processor, depending on the slicing direction. The user has to fill
  the update array in that case so that the global dofs can be mapped
  to the local ones. Filling the update array requires knowledge about
  the environment where the matrix is used, that is why it cannot be
  done by a matrix implementation that is not connected to any
  environment whatsoever.

  Note that on vertical slicing the parallel matrix does not know the
  total numbers of rows.

  \param ps           (o) mask of parallel sparse matrix
  \param rows         (i) global number of rows
  \param cols         (i) global number of columns
  \param slicing      (i) slicing direction
  \param local_cols   (i) local number of columns

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void parallel_sparse_init(PARALLEL_SPARSE* ps, INT rows, INT cols,
                          SLICING_DIRECTION slicing, INT local_cols)
{
#ifdef DEBUG
  dstrc_enter("parallel_sparse_init");
#endif

  ps->total_cols = cols;
  sparse_init(&(ps->slice), rows, local_cols);

#ifdef PARALLEL
  ps->slicing = slicing;
  if (slicing == sd_slice_horizontal)
  {
    ps->update = (INT*)CCACALLOC(rows, sizeof(INT));
  }
  else /* slice_vertical */
  {
    ps->update = (INT*)CCACALLOC(cols, sizeof(INT));
  }
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief free the memory of a parallel sparse matrix

  \param ps           (i) mask of parallel sparse matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void parallel_sparse_destroy(PARALLEL_SPARSE* ps)
{
#ifdef DEBUG
  dstrc_enter("parallel_sparse_destroy");
#endif

#ifdef PARALLEL
  CCAFREE(ps->update);
#endif

  sparse_destroy(&(ps->slice));

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief set up an empty sparse matrix

  The sparse mask needs to be build before it can be used. Umfpack
  style.

  This one is always local. It does not know anything about any other
  processor. In a parallel execution it holds just one slice.

  \param s            (o) mask of local sparse matrix
  \param rows         (i) global number of rows
  \param cols         (i) global number of columns

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_init(SPARSE* s, INT rows, INT cols)
{
#ifdef DEBUG
  dstrc_enter("sparse_init");
#endif

  s->rows = rows;
  s->cols = cols;

  s->Ai = NULL;
  s->Ap = NULL;
  s->Ax = NULL;

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief free the memory of a sparse matrix

  \param s            (i) mask of local sparse matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_destroy(SPARSE* s)
{
#ifdef DEBUG
  dstrc_enter("sparse_destroy");
#endif

  if (s->Ax) CCAFREE(s->Ax);
  if (s->Ap) CCAFREE(s->Ap);
  if (s->Ai) CCAFREE(s->Ai);

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief start mask phase.

  \param s            (o) mask list of sparse matrix
  \param rows         (i) number of rows
  \param cols         (i) number of columns
  \param bandwidth    (i) estimated bandwidth. a guess.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_mask_list_init(SPARSE_MASK_LIST* s,
                           INT rows, INT cols, INT bandwidth)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("sparse_mask_list_init");
#endif

  s->bandwidth = bandwidth;
  s->rows = rows;
  s->cols = cols;

  s->columns = (INTSET*)CCAMALLOC(cols*sizeof(INTSET));
  for (i=0; i<cols; ++i)
  {
    intset_init(&(s->columns[i]), bandwidth);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief finish mask phase

  \param s            (i) mask list of sparse matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_mask_list_destroy(SPARSE_MASK_LIST* s)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("sparse_mask_list_destroy");
#endif

  for (i=0; i<s->cols; ++i)
  {
    intset_destroy(&(s->columns[i]));
  }
  CCAFREE(s->columns);
  s->columns = NULL;

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief mark one entry

  \param s            (o) mask list of sparse matrix
  \param row          (i) row position
  \param col          (i) column position

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_mask_list_mark(SPARSE_MASK_LIST* s, INT row, INT col)
{
#ifdef DEBUG
  dstrc_enter("sparse_mark");
#endif

  dsassert(row >= 0 && row < s->rows, "illegal row");
  dsassert(col >= 0 && col < s->cols, "illegal col");

  if (!intset_contains(&(s->columns[col]),row))
  {
    intset_add(&(s->columns[col]), row);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief create the mask out of the mask list.

  \param s            (o) sparse matrix
  \param sml          (i) mask list of sparse matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_fix_mask(SPARSE* s, SPARSE_MASK_LIST* sml)
{
  INT i;
  INT accumulate = 0;

#ifdef DEBUG
  dstrc_enter("sparse_fix_mask");
#endif

  s->rows = sml->rows;
  s->cols = sml->cols;

  s->Ap = (INT*)CCACALLOC(s->cols+1, sizeof(INT));
  for (i=0; i<s->cols; ++i)
  {
    s->Ap[i] = accumulate;
    accumulate += sml->columns[i].count;
  }
  s->Ap[s->cols] = accumulate;

  s->Ai = (INT*)CCACALLOC(accumulate, sizeof(INT));
  for (i=0; i<s->cols; ++i)
  {
    INT j;
    INT* p;
    p = &(s->Ai[s->Ap[i]]);
    for (j=0; j<sml->columns[i].count; ++j)
    {
      *p++ = sml->columns[i].value[j];
    }
  }

  /* If there were memory issues we could allocate this one after the
   * mask list has been freed. */
  s->Ax = (DOUBLE*)CCAMALLOC(accumulate*sizeof(DOUBLE));

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief zero out the matrix

  \param s            (o) sparse matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_zero(SPARSE* s)
{
#ifdef DEBUG
  dstrc_enter("sparse_zero");
#endif

  memset(s->Ax, 0, s->Ap[s->cols]*sizeof(DOUBLE));

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief get matrix entry

  \param s            (i) sparse matrix
  \param row          (i) row position
  \param col          (i) column position

  \return pointer to matrix entry

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
DOUBLE* sparse_entry(SPARSE* s, INT row, INT col)
{
  INT* p;
#ifdef DEBUG
  dstrc_enter("sparse_entry");
#endif

  dsassert(row >= 0 && row < s->rows, "illegal row");
  dsassert(col >= 0 && col < s->cols, "illegal col");

  p = bsearch(&row, &(s->Ai[s->Ap[col]]), s->Ap[col+1] - s->Ap[col], sizeof(INT), cmp_int);
  dsassert(p != NULL, "row not in matrix");

#ifdef DEBUG
  dstrc_exit();
#endif

  return &(s->Ax[p - s->Ai]);
}


/*----------------------------------------------------------------------*/
/*!
  \brief test if matrix contains entry (row,col)

  \param s            (i) sparse matrix
  \param row          (i) row position
  \param col          (i) column position

  \return True if entry is nonzero.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
INT sparse_has_entry(SPARSE* s, INT row, INT col)
{
  INT* p;
#ifdef DEBUG
  dstrc_enter("sparse_has_entry");
#endif

  dsassert(row >= 0 && row < s->rows, "illegal row");
  dsassert(col >= 0 && col < s->cols, "illegal col");

  p = bsearch(&row, &(s->Ai[s->Ap[col]]), s->Ap[col+1] - s->Ap[col], sizeof(INT), cmp_int);

#ifdef DEBUG
  dstrc_exit();
#endif
  return p != NULL;
}


#ifdef DEBUG

/*----------------------------------------------------------------------*/
/*!
  \brief dump the matrix to a python file for debugging

  For debugging only!

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void sparse_debugdump(SPARSE* s,
                      INT* update,
                      INT update_length,
                      FILE* f,
                      CHAR* name)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("sparse_debugdump");
#endif

  fprintf(f, "%s = {\n", name);

  fprintf(f, "    'rows' : %d,\n", s->rows);
  fprintf(f, "    'cols' : %d,\n", s->cols);

  fprintf(f, "    'update' : [");
  for (i=0; i<update_length; ++i)
  {
    fprintf(f, "%d, ", update[i]);
  }
  fprintf(f, "],\n");

  fprintf(f, "    'Ap' : [");
  for (i=0; i<s->cols+1; ++i)
  {
    fprintf(f, "%d, ", s->Ap[i]);
  }
  fprintf(f, "],\n");

  fprintf(f, "    'Ai' : [");
  for (i=0; i<s->Ap[s->cols]; ++i)
  {
    fprintf(f, "%d, ", s->Ai[i]);
  }
  fprintf(f, "],\n");

  fprintf(f, "    'Ax' : [");
  for (i=0; i<s->Ap[s->cols]; ++i)
  {
    fprintf(f, "%e, ", s->Ax[i]);
  }
  fprintf(f, "],\n");

  fprintf(f, "}\n");

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


#ifdef D_FLUID_PM

/* The parallel version is more difficult. We need to do an extra
 * indirection because the columns can have arbitrary ids. */

#ifdef PARALLEL

/*----------------------------------------------------------------------*/
/*!
  \brief store C^T*ML^-1*C in pmat.

  lmass is a simple array of all mass dofs, velocity dirichlet dofs
  removed.

  C^T is the discrete divergence operator.
  C is the discrete gradient operator.

  parallel version

  \param pmat         (o) pressure sparse matrix slice
  \param grad_col_ids (i) update array of gradient slice
  \param grad         (i) gradient slice
  \param div_col_ids  (i) update array of divergence slice
  \param div          (i) divergence slice
  \param lmass        (i) global inverted lumped mass matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
static void sparse_pm_matmat(SPARSE* pmat,
                             INT* grad_col_ids,
                             SPARSE* grad,
                             INT* div_col_ids,
                             SPARSE* div,
                             DOUBLE* lmass)
{
  INT col;

  /* Loop all nonempty columns and rows of the destination matrix. */
  /* For the parallel version we only have those columns that are
   * local in this slice of the gradient matrix. */

  for (col=0; col<grad->cols; ++col)
  {
    INT Ae;
    INT Ac;

    Ac = grad_col_ids[col];

    for (Ae=pmat->Ap[Ac]; Ae<pmat->Ap[Ac+1]; ++Ae)
    {
      INT e1, e2;

      /*
       * the row numbers are always local to the matrix, independent
       * of the slicing direction. */
      INT Ar = pmat->Ai[Ae];

#if 0
      // Wir sind in C.
      // Die Zeile Ar von C^T (Spalte Ar von C) mit der Spalte Ac
      // von C multiplizieren.

      // Also zwei Sparse-Spalten miteinander skalar
      // multiplizieren.
#endif

      e1 = div->Ap[Ar];
      e2 = grad->Ap[col];
      while ((e1 < div->Ap[Ar+1]) && (e2 < grad->Ap[col+1]))
      {
        INT r1, r2;
        r1 = div->Ai[e1];
        r2 = grad->Ai[e2];
        if (r1 < r2)
        {
          ++e1;
        }
        else if (r1 > r2)
        {
          ++e2;
        }
        else /* if (r1==r2) */
        {
          pmat->Ax[Ae] += div->Ax[e1]*lmass[r1]*grad->Ax[e2];
          ++e1;
          ++e2;
        }
      }
    }
  }
}

#else

/*----------------------------------------------------------------------*/
/*!
  \brief store C^T*ML^-1*C in pmat.

  lmass is a simple array, velocity dirichlet dofs removed.

  sequential version.

  We visit all entries in the new matrix that we know are not
  zero. For these we take the two columns from the grad matrix and
  scalar multiply them.

  \param pmat         (o) pressure sparse matrix slice
  \param grad         (i) gradient matrix
  \param lmass        (i) global inverted lumped mass matrix

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
static void sparse_pm_matmat(SPARSE* pmat, SPARSE* grad, DOUBLE* lmass)
{
  INT Ac;

  /* Loop all nonempty columns and rows of the destination matrix. */

  for (Ac=0; Ac<pmat->cols; ++Ac)
  {
    INT Ae;

    for (Ae=pmat->Ap[Ac]; Ae<pmat->Ap[Ac+1]; ++Ae)
    {
      INT e1, e2;
      INT Ar = pmat->Ai[Ae];

#if 0
      // Wir sind in C.
      // Die Zeile Ar von C^T (Spalte Ar von C) mit der Spalte Ac
      // von C multiplizieren.

      // Also zwei Sparse-Spalten miteinander skalar
      // multiplizieren.
#endif

      e1 = grad->Ap[Ar];
      e2 = grad->Ap[Ac];
      while ((e1 < grad->Ap[Ar+1]) && (e2 < grad->Ap[Ac+1]))
      {
        INT r1, r2;
        r1 = grad->Ai[e1];
        r2 = grad->Ai[e2];
        if (r1 < r2)
        {
          ++e1;
        }
        else if (r1 > r2)
        {
          ++e2;
        }
        else /* if (r1==r2) */
        {
          pmat->Ax[Ae] += grad->Ax[e1]*lmass[r1]*grad->Ax[e2];
          ++e1;
          ++e2;
        }
      }
    }
  }
}

#endif

/*----------------------------------------------------------------------*/
/*!
  \brief calculate C^T*ML^-1*C for the projection method

  This is a very specific support function for the projection
  method. Yet it does not depend on any external structures, thus it
  is completely fine to have it here.

  In the parallel case we have to communicate all the slices of the
  gradient operator matrix to calculate the matrix product.

  \param pmat         (o) parallel pressure sparse matrix
  \param grad         (i) parallel gradient operator matrix
  \param lmass        (i) global inverted lumped mass matrix
  \param actintra     (i) communicator

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void parallel_sparse_pm_matmat(PARALLEL_SPARSE* pmat, PARALLEL_SPARSE* grad,
                               DOUBLE* lmass, INTRA* actintra)
{
#ifdef PARALLEL
  INT local_size[2];
  INT global_size[2];

  INT* grad_col_ids;
  SPARSE grad_slice;

  INT i;
  INT err;
#endif

#ifdef DEBUG
  dstrc_enter("parallel_sparse_pm_matmat");
#endif

  dsassert(pmat->slice.rows == grad->slice.cols, "dimension mismatch");
  dsassert(pmat->slice.cols == grad->total_cols, "dimension mismatch");

  sparse_zero(&(pmat->slice));

#ifdef PARALLEL

  dsassert(grad->slicing==sd_slice_vertical, "odd grad slicing");
  dsassert(pmat->slicing==sd_slice_horizontal, "odd pmat slicing");

  /* For projections on a fixed grid (and these are the only ones
   * supported) we need this product just once. Thus we don't do the
   * setup beforehand. */

  /* We need to communicate a whole sparse matrix. So lets allocate
   * space big enough for the biggest slice. */

  local_size[0] = grad->slice.cols;
  local_size[1] = grad->slice.Ap[grad->slice.cols];
  err = MPI_Allreduce(local_size, global_size, 2, MPI_INT, MPI_MAX, actintra->MPI_INTRA_COMM);
  if (err != 0)
  {
    dserror("mpi sendrecv error %d", err);
  }

  grad_col_ids = (INT*)CCAMALLOC(global_size[0]*sizeof(INT));
  grad_slice.Ap = (INT*)CCAMALLOC((global_size[0]+1)*sizeof(INT));
  grad_slice.Ai = (INT*)CCAMALLOC(global_size[1]*sizeof(INT));
  grad_slice.Ax = (DOUBLE*)CCAMALLOC(global_size[1]*sizeof(DOUBLE));
  grad_slice.rows = grad->slice.rows;

  /* the curious communication pattern without local communication */
  /* See BINIO documentation. */
  for (i=0; i<actintra->intra_nprocs-1; ++i)
  {
    MPI_Status status;
    INT dst;
    INT src;
    INT recv_size[2];
    dst = (actintra->intra_rank + i + 1) % actintra->intra_nprocs;
    src = (actintra->intra_nprocs + actintra->intra_rank - i - 1) % actintra->intra_nprocs;

    /* We need to receive the whole matrix. */
    /* If we would be very clever we could know in advance which
     * columns we need and which ones only produce zeros... */

    err = MPI_Sendrecv(local_size, 2, MPI_INT, dst, 123+i,
                       recv_size,  2, MPI_INT, src, 123+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0)
    {
      dserror("mpi sendrecv error %d", err);
    }

    err = MPI_Sendrecv(grad->update, local_size[0], MPI_INT, dst, 2345+i,
                       grad_col_ids, recv_size[0], MPI_INT, src, 2345+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0)
    {
      dserror("mpi sendrecv error %d", err);
    }

    err = MPI_Sendrecv(grad->slice.Ap, local_size[0]+1, MPI_INT, dst, 234+i,
                       grad_slice.Ap, recv_size[0]+1, MPI_INT, src, 234+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0)
    {
      dserror("mpi sendrecv error %d", err);
    }

    err = MPI_Sendrecv(grad->slice.Ai, local_size[1], MPI_INT, dst, 345+i,
                       grad_slice.Ai, recv_size[1], MPI_INT, src, 345+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0)
    {
      dserror("mpi sendrecv error %d", err);
    }

    err = MPI_Sendrecv(grad->slice.Ax, local_size[1], MPI_DOUBLE, dst, 3456+i,
                       grad_slice.Ax, recv_size[1], MPI_DOUBLE, src, 3456+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0)
    {
      dserror("mpi sendrecv error %d", err);
    }

    grad_slice.cols = recv_size[0];

    /* Do the multiplication with the currently available slice. */
    sparse_pm_matmat(&(pmat->slice),
                     grad_col_ids, &grad_slice,
                     grad->update, &(grad->slice), lmass);
  }

  sparse_destroy(&grad_slice);
  CCAFREE(grad_col_ids);

  /* finally do the diagonal blocks */
  sparse_pm_matmat(&(pmat->slice),
                   grad->update, &(grad->slice),
                   grad->update, &(grad->slice), lmass);

#ifdef DEBUG
  /* Check diagonal */
  for (i=0; i<pmat->slice.rows; ++i)
  {
    INT global_row;
    global_row = pmat->update[i];
    if (*sparse_entry(&(pmat->slice), i, global_row) == 0.0)
    {
      dserror("zero diagonal at pmat(%d,%d)", i, global_row);
    }
  }
#endif

#else
  sparse_pm_matmat(&(pmat->slice), &(grad->slice), lmass);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif

#ifdef D_FLUID_PM
#ifdef AZTEC_PACKAGE

/*----------------------------------------------------------------------*/
/*!
  \brief convert the given sparse matrix to an aztec solver object

  The actual values are not considered here, those might not be set
  yet. This function is concerned with the sparse mask only.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void parallel_sparse_convert_aztec(PARALLEL_SPARSE* ps, SOLVAR* solvar)
{
#ifdef FAST_ASS
#error No way! We want to generate a sparse matrix without connection to the element/node structure
#endif

  AZ_ARRAY_MSR* array_msr;
  INT* bindx;
  DOUBLE* val;
  INT i;
  INT Ac;
  INT accumulate;
  INT* fill;

#ifdef DEBUG
  dstrc_enter("parallel_sparse_convert_aztec");
#endif

#ifdef PARALLEL

  dsassert(ps->slicing==sd_slice_horizontal, "odd ps slicing");

#endif

  /* Not very nice... */
  solvar->fieldtyp = fluid;

  solvar->solvertyp = aztec_msr;
  solvar->parttyp = cut_elements;
  solvar->matrixtyp = matrix_none;
  solvar->azvar = (AZVAR*)CCACALLOC(1,sizeof(AZVAR));

  /* Default values. There is no input we could read... */

  /*solvar->azvar->azsolvertyp = azsolv_CG;*/
  solvar->azvar->azsolvertyp = azsolv_GMRES;
  solvar->azvar->azprectyp = azprec_SymmGaussSeidel;

  solvar->azvar->azreuse   = 0;
  solvar->azvar->azgfill   = 0;
  solvar->azvar->aziter    = 2500;
  solvar->azvar->azsub     = 50;
  solvar->azvar->azgraph   = 1;
  solvar->azvar->azpoly    = 5;
  solvar->azvar->blockdiag = 0;

  solvar->azvar->azdrop  = 0.;
  solvar->azvar->azfill  = 2.;
  solvar->azvar->aztol   = 1.e-12;
  solvar->azvar->azomega = 1.;

  /* We create one array */

  solvar->nsysarray = 1;
  solvar->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(solvar->nsysarray,sizeof(SPARSE_TYP));
  solvar->sysarray     = (SPARSE_ARRAY*)CCACALLOC(solvar->nsysarray,sizeof(SPARSE_ARRAY));

  solvar->sysarray_typ[0] = msr;
  solvar->sysarray[0].msr = (AZ_ARRAY_MSR*)CCACALLOC(1,sizeof(AZ_ARRAY_MSR));
  array_msr = solvar->sysarray[0].msr;

  /* Build up the sparse mask. */

  array_msr->numeq_total = ps->total_cols;
  array_msr->numeq = ps->slice.rows;

  array_msr->nnz = ps->slice.Ap[ps->slice.cols];

  amdef("update",&(array_msr->update),array_msr->numeq,1,"IV");
#ifdef PARALLEL
  memcpy(array_msr->update.a.iv, ps->update, array_msr->numeq*sizeof(INT));
#else
  for (i=0; i<array_msr->numeq; ++i)
  {
    array_msr->update.a.iv[i] = i;
  }
#endif

  bindx =  (INT*)amdef("bindx",&(array_msr->bindx),(array_msr->nnz+1),1,"IV");
  val = (DOUBLE*)amdef("val"  ,&(array_msr->val)  ,(array_msr->nnz+1),1,"DV");
  array_msr->bindx_backup.Typ = cca_XX;

  bindx[0] = ps->slice.rows + 1;
  memset(bindx+1, 0, ps->slice.rows*sizeof(INT));

  /* The number of entries in a row */

  for (Ac=0; Ac<ps->slice.cols; ++Ac)
  {
    INT Ae;
    for (Ae=ps->slice.Ap[Ac]; Ae<ps->slice.Ap[Ac+1]; ++Ae)
    {
      INT Ar = ps->slice.Ai[Ae];

      /* Exclude diagonal entries. */
#ifdef PARALLEL
      if (Ac != ps->update[Ar])
#else
      if (Ac != Ar)
#endif
      {
        bindx[Ar+1] += 1;
      }
    }
  }

  /* aztec stores the number per row incrementally, excluding the
   * diagonal entries. */

  accumulate = 0;
  for (i=0; i<ps->slice.rows; ++i)
  {
    bindx[i] += accumulate;
    accumulate = bindx[i];
  }
  bindx[ps->slice.rows] += accumulate;

  /* fill in all column positions */
  /* Copy values */

  fill = (INT*)CCACALLOC(ps->slice.rows, sizeof(INT));

  for (Ac=0; Ac<ps->slice.cols; ++Ac)
  {
    INT Ae;
    for (Ae=ps->slice.Ap[Ac]; Ae<ps->slice.Ap[Ac+1]; ++Ae)
    {
      INT Ar = ps->slice.Ai[Ae];

#ifdef PARALLEL
      if (Ac != ps->update[Ar])
#else
      if (Ac != Ar)
#endif
      {
        bindx[bindx[Ar]+fill[Ar]] = Ac;
        /*val  [bindx[Ar]+fill[Ar]] = ps->slice.Ax[Ae];*/
        fill[Ar] += 1;
      }
      else
      {
        /*val[Ar] = ps->slice.Ax[Ae];*/
      }
    }
  }

#ifdef DEBUG

  /* self check */
  for (i=0; i<ps->slice.rows; ++i)
  {
    if (fill[i] != bindx[i+1]-bindx[i])
    {
      dserror("filled %d entries in line %d, expected %d", fill[i], i, bindx[i+1]-bindx[i]);
    }
  }

#endif

  CCAFREE(fill);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief copy all values from a sparse matrix to an aztec solver object

  Both matrices must have the same structure. Most often the aztec
  object should be generated by \ref parallel_sparse_convert_aztec.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void parallel_sparse_copy_aztec(PARALLEL_SPARSE* ps, SOLVAR* solvar)
{
#ifdef FAST_ASS
#error No way! We want to generate a sparse matrix without connection to the element/node structure
#endif

  AZ_ARRAY_MSR* array_msr;
  INT* bindx;
  DOUBLE* val;
  INT* fill;
  INT i;
  INT Ac;

#ifdef DEBUG
  dstrc_enter("parallel_sparse_convert_aztec");
#endif

  /* Just a few tests */

  dsassert(solvar->solvertyp == aztec_msr, "wrong solver type");
  dsassert(solvar->sysarray_typ[0] == msr, "wrong sparse type");

  array_msr = solvar->sysarray[0].msr;
  dsassert(array_msr->numeq_total == ps->total_cols, "dimension mismatch");
  dsassert(array_msr->numeq == ps->slice.rows, "dimension mismatch");

#ifdef PARALLEL

  dsassert(ps->slicing==sd_slice_horizontal, "odd ps slicing");

#endif

  /* We reset the values. So we reinitialize the solver. No reuse
   * whatsoever. */
  array_msr->is_factored=0;
  array_msr->is_transformed=0;

  bindx = array_msr->bindx.a.iv;
  val   = array_msr->val.a.dv;

  /* Copy values */

  /* We still use an auxiliary fill array here. This way we avoid
   * searching the column entry. We can do this because the columns in
   * the original matrix are orderned. */

  fill = (INT*)CCACALLOC(ps->slice.rows, sizeof(INT));

  for (Ac=0; Ac<ps->slice.cols; ++Ac)
  {
    INT Ae;
    for (Ae=ps->slice.Ap[Ac]; Ae<ps->slice.Ap[Ac+1]; ++Ae)
    {
      INT Ar = ps->slice.Ai[Ae];

#ifdef PARALLEL
      if (Ac != ps->update[Ar])
#else
      if (Ac != Ar)
#endif
      {
        val[bindx[Ar]+fill[Ar]] = ps->slice.Ax[Ae];
        fill[Ar] += 1;
        dsassert(fill[Ar] <= bindx[Ar+1]-bindx[Ar], "fill overflow");
      }
      else
      {
        val[Ar] = ps->slice.Ax[Ae];
      }
    }
  }

#ifdef DEBUG

  /* self check */
  for (i=0; i<ps->slice.rows; ++i)
  {
    if (fill[i] != bindx[i+1]-bindx[i])
    {
      dserror("filled %d entries in line %d, expected %d", fill[i], i, bindx[i+1]-bindx[i]);
    }
  }

  for (Ac=0; Ac<ps->slice.rows; ++Ac)
  {
    if (val[Ac] == 0.0)
    {
      dserror("zero diagonal at row %d of %d", Ac, ps->slice.rows);
    }
  }

#endif

  CCAFREE(fill);

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
#endif
