/*!----------------------------------------------------------------------
\file
\brief contains functions to handle oll matrices

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
/* #include "../headers/prototypes_sol.h" */
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
/*!
\addtogroup OLL
*//*! @{ (documentation module open)*/


/*!---------------------------------------------------------------------
\brief opens a matrix in oll format

<pre>                                                         mn 02/03
This function creates an oll matrix for the given number of equations,
building the update vector and allocating the spare vector.
</pre>
\param *matrix       OLL  (i/o) the matrix object to be created
\param  numeq        INT  (i)   local number of equations
\param  numeq_total  INT  (i)   global number of equations
\param *actfield     FIELD  (i)   the active field
\param *actpart      PARTITION  (i)   the active partition
\param *actintra     INTRA  (i)   the active communicator

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_open(OLL       *matrix,
              INT        numeq,
              INT        numeq_total,
              FIELD     *actfield,
              PARTITION *actpart,
              INTRA     *actintra,
              INT        dis)
{
  INT       *update;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_open");
#endif
  /*----------------------------------------------------------------------*/
  matrix->numeq       = numeq;
  matrix->numeq_total = numeq_total;
  matrix->nnz         =0;
  matrix->is_copied   =0;
  matrix->is_masked   =0;
  /*------------------------------------- if the ll storage doesn't exist */
  if (matrix->row==NULL)
  {
    matrix->rdim = matrix->numeq;
    matrix->cdim = matrix->numeq_total;
    matrix->row  = (MATENTRY**)CCACALLOC(matrix->rdim,sizeof(MATENTRY*));
    matrix->col  = (MATENTRY**)CCACALLOC(matrix->cdim,sizeof(MATENTRY*));
    update = amdef("update", &(matrix->update), numeq, 1, "IV");

    oll_update(actfield, actpart, actintra, dis, matrix);
    oll_nnz_topology(actfield, actpart, actintra, matrix, dis);

    matrix->total = matrix->nnz;
    matrix->used  = 0;
    matrix->spare = (MATENTRY*)CCACALLOC(matrix->total,sizeof(MATENTRY));

  }
  else
  {
    dserror("matrix already exists");
  }
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_open */



/*!---------------------------------------------------------------------
\brief fetches a matrix entry out of the spare vector

<pre>                                                        mn 02/03
This function fetches a matrix entry from the spare vector if there are
still unused entries. It allocates an new entry otherwise.
</pre>
\param *matrix       OLL  (i)   the oll matrix
\param *ret          MATENTRY  (o)   the matrix entry

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_getentry(
    OLL*  matrix,
    MATENTRY** ret)
{
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_getentry");
#endif
  /*----------------------------------------------------------------------*/
  if((matrix->total - matrix->used)>0)
  {
    /* there are still unused MATENTRIES in spare */
    *ret = (&(matrix->spare[matrix->used]));
    matrix->used++;
  }
  else
  {
    *ret  = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));
  }

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}



/*!---------------------------------------------------------------------
\brief zeros an oll matrix

<pre>                                                        mn 03/03
This functions sets all entries in the given oll matrix to zero, the
poiter structure is preserved.
</pre>
\param *matrix       OLL  (i)   the oll matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_zero(OLL*  oll)
{
  MATENTRY  **row;
  MATENTRY   *actentry;
  INT         i;
  DOUBLE      null = 0.0;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_zero");
#endif
  /*----------------------------------------------------------------------*/
  row = oll->row;

  for (i=0; i<oll->rdim; i++)
  {
    actentry = row[i];
    while(actentry != NULL )
    {
      actentry->val = null;
      actentry = actentry->rnext;
    } /* end while row i */
  } /* end for all rows */

  oll->is_copied = 0;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_zero */




/*!---------------------------------------------------------------------
\brief adds two oll matrices

<pre>                                                        mn 03/03
This function performs
oll1 = oll1 + oll2 * factor
for two oll matrices.  These must have the same mask.
</pre>
\param *oll1         OLL    (i/o) the first oll matrix
\param *oll2         OLL    (i)   the second oll matrix
\param  factor       DOUBLE (i)   the multiplication factor

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_add(
    OLL    *oll1,
    OLL    *oll2,
    DOUBLE  factor)
{
  INT      i;
  MATENTRY  **row1;
  MATENTRY   *actentry1;
  MATENTRY  **row2;
  MATENTRY   *actentry2;

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_add");
#endif
  /*----------------------------------------------------------------------*/
  row1 = oll1->row;
  row2 = oll2->row;

  for (i=0; i<oll1->rdim; i++)
  {
    actentry1 = row1[i];
    actentry2 = row2[i];
    while(actentry1 != NULL )
    {
      if (actentry1->c != actentry2->c) dserror("Incompatible matrices");
      actentry1->val += actentry2->val*factor;
      actentry1 = actentry1->rnext;
      actentry2 = actentry2->rnext;
    } /* end while row i */
  } /* end for all rows */
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_add */



/*!---------------------------------------------------------------------
\brief scales an oll matrix

<pre>                                                        mn 03/03
This function multiplies all entries of the oll matrix with the given
factor.
</pre>
\param *oll          OLL    (i)   the second oll matrix
\param  factor       DOUBLE (i)   the multiplication factor

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_scal(
    OLL    *oll,
    DOUBLE  factor)
{
  INT      i;
  MATENTRY  **row;
  MATENTRY   *actentry;

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_scal");
#endif
  /*----------------------------------------------------------------------*/
  row = oll->row;

  for (i=0; i<oll->rdim; i++)
  {
    actentry = row[i];
    while(actentry != NULL )
    {
      actentry->val = actentry->val*factor;
      actentry = actentry->rnext;
    } /* end while row i */
  } /* end for all rows */
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_scal */




/*!---------------------------------------------------------------------
\brief copies the structure of an oll matrix

<pre>                                                        mn 03/03
This function copies the structure of an oll matrix, not the values or
the pointers. Only the update vector is copied and the row- and column
pointer vectors and the spare vector are allocated.
</pre>
\param *from        OLL    (i)   the oll matrix copied from
\param *to          OLL    (i)   the oll matrix copied to

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_cp_mask(
    OLL *from,
    OLL *to)
{

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_cp_mask");
#endif
  /*----------------------------------------------------------------------*/
  /* copy all information, which is directly included in the structure */
  *to = *from;

  /* alloccopy update */
  am_alloc_copy(&(from->update),&(to->update));
  to->row   = (MATENTRY**)CCACALLOC(to->rdim,sizeof(MATENTRY*));
  to->col   = (MATENTRY**)CCACALLOC(to->cdim,sizeof(MATENTRY*));
  to->spare = (MATENTRY*)CCACALLOC(to->total,sizeof(MATENTRY));
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_cp_mask */



/*!---------------------------------------------------------------------
\brief deletes an oll matrix

<pre>                                                        mn 03/03
This function deletes the contents of the structure OLL
</pre>
\param *matrix        OLL    (i)   the oll matrix to be deleted

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_delete(OLL*  matrix)
{
  INT        i,dim;
  MATENTRY  *rowptr,*nextptr;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_delete");
#endif
  /*----------------------------------------------------------------------*/
  /*-------------------------------------------- destroy the general part */
  if (matrix->update.Typ != cca_XX)
    amdel(&(matrix->update));
  if (matrix->gdofrecv.Typ != cca_XX)
    amdel(&(matrix->gdofrecv));
  if (matrix->recvbuff.Typ != cca_XX)
    amdel(&(matrix->recvbuff));
  if (matrix->computebuff.Typ != cca_XX)
    amdel(&(matrix->computebuff));
  if (matrix->gdofsend.Typ != cca_XX)
    amdel(&(matrix->gdofsend));
  if (matrix->sendbuff.Typ != cca_XX)
    amdel(&(matrix->sendbuff));
#ifdef PARALLEL
  if (matrix->status != NULL)
    CCAFREE(matrix->status);
  if (matrix->request != NULL)
    CCAFREE(matrix->request);
#endif
  /*------------------------------------------- destroy the matrix itself */
  if (matrix->col != NULL) matrix->col = CCAFREE(matrix->col);
  if (matrix->row != NULL)
  {
    dim = matrix->rdim;
    for (i=0; i<dim; i++)
    {
      rowptr = matrix->row[i];
      if (rowptr == NULL) continue;
      /* move pointer forwards to the end of row */
      while (rowptr != NULL)
      {
        nextptr = rowptr->rnext;
        CCAFREE(rowptr);
        rowptr = nextptr;
      }
      /* we are now positioned on the last entry in row */
    }
  }
  matrix->row = CCAFREE(matrix->row);
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_delete */



/*!---------------------------------------------------------------------
\brief sets a value in an oll matrix

<pre>                                                        mn 02/03
This function sets a value in an oll matrix. If the matrix entry already
exists, the existing value is overwritten, otherwise a new entry is
fetched.
</pre>
\param *matrix    OLL    (i/o)   the matrix
\param  actrow    INT      (i)   the true row index
\param  actcol    INT      (i)   the true col index
\param  val       DOUBLE   (i)   the value to set

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_setval(
    OLL*  matrix,
    INT actrow,
    INT actcol,
    DOUBLE val)
{
  MATENTRY  **row;
  MATENTRY  **col;
  INT         rindex;
  MATENTRY   *actptr;
  MATENTRY   *tmpptr;
  MATENTRY   *cptr;
  INT         success;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_setval");
#endif
  /*----------------------------------------------------------------------*/
  row = matrix->row;
  col = matrix->col;
  /*------------------------------------------------- get the correct row */
  rindex = find_index(actrow,matrix->update.a.iv,matrix->numeq);
  if (rindex==-1) dserror("row not in update");
  /*------------------------------------------ check whether row is empty */
  success=0;
  if (row[rindex]==NULL)
  {
    row[rindex]      = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));
    matrix->nnz++;
    row[rindex]->r   = actrow;
    row[rindex]->c   = actcol;
    row[rindex]->val = val;
    actptr           = row[rindex];
    success++;
  }
  else /*----------------------------------------------- row is not empty */
  {
    actptr           = row[rindex];
startrow:
    if (actptr->c > actcol) /* first entry larger then actcol */
    {
      dsassert(actptr==row[rindex],"Matrix got severly mixed up");
      tmpptr        = actptr;
      row[rindex]   = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));
      matrix->nnz++;
      actptr        = row[rindex];
      actptr->r     = actrow;
      actptr->c     = actcol;
      actptr->val   = val;
      actptr->rnext = tmpptr;
      success++;
      goto endrow;
    }
    if (actptr->c == actcol) /* this is correct place */
    {
      dsassert(actptr->r==actrow,"Matrix got mixed up");
      actptr->val = val;
      success++;
      goto endrow;
    }
    if (actptr->rnext == NULL) /* next element not exist */
    {
      actptr->rnext = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));
      matrix->nnz++;
      actptr        = actptr->rnext;
      actptr->r     = actrow;
      actptr->c     = actcol;
      actptr->val   = val;
      success++;
      goto endrow;
    }
    if (actptr->rnext->c == actcol) /* next element correct element */
    {
      actptr        = actptr->rnext;
      dsassert(actptr->c==actcol && actptr->r==actrow,"Matrix got mixed up");
      actptr->val   = val;
      success++;
      goto endrow;
    }
    if (actptr->rnext->c > actcol &&
        actptr->c        < actcol ) /* next element larger and active smaller */
    {
      tmpptr        = actptr->rnext;
      actptr->rnext = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));
      matrix->nnz++;
      actptr        = actptr->rnext;
      actptr->rnext = tmpptr;
      actptr->r     = actrow;
      actptr->c     = actcol;
      actptr->val   = val;
      success++;
      goto endrow;
    }
    /* did not find the proper place */
    actptr = actptr->rnext;
    goto startrow;
endrow:;
  }
  if (!success) dserror("Failed to set entry");
  success=0;
  /*---------------------------------------------- check dimension of col */
  if (actcol >= matrix->cdim) dserror("Overvlow in column dimension");
  /*--------------------------------------- check whether column is empty */
  if (col[actcol]==NULL)
  {
    col[actcol] = actptr;
    success++;
  }
  else /*-------------------------------------------- column is not empty */
  {
    cptr = col[actcol];
startcol:
    if (cptr->r >  actrow)
    {
      dsassert(cptr==col[actcol],"Matrix got severly mixed up");
      dsassert(actptr->cnext==NULL,"Matrix got severly mixed up");
      tmpptr        = cptr;
      col[actcol]   = actptr;
      actptr->cnext = tmpptr;
      success++;
      goto endcol;
    }
    if (cptr->r == actptr->r) /* this entry is already linked in */
    {
      success++;
      goto endcol;
    }
    if (cptr->cnext==NULL) /* next entry is not assigned */
    {
      cptr->cnext   = actptr;
      success++;
      goto endcol;
    }
    if (cptr->r < actrow && cptr->cnext->r > actrow) /* active smaller, next larger */
    {
      tmpptr        = cptr->cnext;
      cptr->cnext   = actptr;
      actptr->cnext = tmpptr;
      success++;
      goto endcol;
    }
    /* dis not find the proper place */
    cptr = cptr->cnext;
    goto startcol;
endcol:;
  }
  if (!success) dserror("Failed to set entry");
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_setval */




/*!---------------------------------------------------------------------
\brief adds a value in an oll matrix

<pre>                                                        mn 02/03
This function adds a value in an oll matrix. If the matrix entry already
exists, the values are added together, otherwise a new entry is
fetched.
</pre>
\param *matrix    OLL    (i/o)   the matrix
\param  actrow    INT      (i)   the true row index
\param  actcol    INT      (i)   the true col index
\param  val       DOUBLE   (i)   the value to set

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_addval(
                OLL*  matrix,
                INT actrow,
                INT actcol,
                DOUBLE val)
{
  MATENTRY  **row;
  MATENTRY  **col;
  INT         rindex;
  MATENTRY   *actptr;
  MATENTRY   *tmpptr;
  MATENTRY   *cptr;
  INT         success;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_addval");
#endif
  /*----------------------------------------------------------------------*/
  row = matrix->row;
  col = matrix->col;
  /*------------------------------------------------- get the correct row */
  rindex = find_index(actrow,matrix->update.a.iv,matrix->numeq);
  if (rindex==-1) dserror("row not in update");
  /*------------------------------------------ check whether row is empty */
  success=0;
  if (row[rindex]==NULL)
  {
    /*row[rindex]      = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
    oll_getentry(matrix,&row[rindex]);
    matrix->nnz++;
    row[rindex]->r   = actrow;
    row[rindex]->c   = actcol;
    row[rindex]->val = val;
    actptr           = row[rindex];
    success++;
  }
  else /*----------------------------------------------- row is not empty */
  {
    actptr           = row[rindex];
startrow:
    if (actptr->c > actcol) /* first entry larger then actcol */
    {
      dsassert(actptr==row[rindex],"Matrix got severly mixed up");
      tmpptr        = actptr;
      /*   row[rindex]   = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
      oll_getentry(matrix,&row[rindex]);
      matrix->nnz++;
      actptr        = row[rindex];
      actptr->r     = actrow;
      actptr->c     = actcol;
      actptr->val   = val;
      actptr->rnext = tmpptr;
      success++;
      goto endrow;
    }
    if (actptr->c == actcol) /* this is correct place */
    {
      dsassert(actptr->r==actrow,"Matrix got mixed up");
      actptr->val += val;
      success++;
      goto endrow;
    }
    if (actptr->rnext == NULL) /* next element not exist */
    {
      /*actptr->rnext = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
      oll_getentry(matrix,&actptr->rnext);
      matrix->nnz++;
      actptr        = actptr->rnext;
      actptr->r     = actrow;
      actptr->c     = actcol;
      actptr->val   = val;
      success++;
      goto endrow;
    }
    if (actptr->rnext->c == actcol) /* next element correct element */
    {
      actptr        = actptr->rnext;
      dsassert(actptr->c==actcol && actptr->r==actrow,"Matrix got mixed up");
      actptr->val  += val;
      success++;
      goto endrow;
    }
    if (actptr->rnext->c > actcol &&
        actptr->c        < actcol ) /* next element larger and active smaller */
    {
      tmpptr        = actptr->rnext;
      /*actptr->rnext = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
      oll_getentry(matrix,&actptr->rnext);
      matrix->nnz++;
      actptr        = actptr->rnext;
      actptr->rnext = tmpptr;
      actptr->r     = actrow;
      actptr->c     = actcol;
      actptr->val   = val;
      success++;
      goto endrow;
    }
    /* did not find the proper place */
    actptr = actptr->rnext;
    goto startrow;
endrow:;
  }
  if (!success) dserror("Failed to set entry");
  success=0;
  /*---------------------------------------------- check dimension of col */
  if (actcol >= matrix->cdim) dserror("Overvlow in column dimension");
  /*--------------------------------------- check whether column is empty */
  if (col[actcol]==NULL)
  {
    col[actcol] = actptr;
    success++;
  }
  else /*-------------------------------------------- column is not empty */
  {
    cptr = col[actcol];
startcol:
    if (cptr->r >  actrow)
    {
      dsassert(cptr==col[actcol],"Matrix got severly mixed up");
      dsassert(actptr->cnext==NULL,"Matrix got severly mixed up");
      tmpptr        = cptr;
      col[actcol]   = actptr;
      actptr->cnext = tmpptr;
      success++;
      goto endcol;
    }
    if (cptr->r == actptr->r) /* this entry is already linked in */
    {
      success++;
      goto endcol;
    }
    if (cptr->cnext==NULL) /* next entry is not assigned */
    {
      cptr->cnext   = actptr;
      success++;
      goto endcol;
    }
    if (cptr->r < actrow && cptr->cnext->r > actrow) /* active smaller, next larger */
    {
      tmpptr        = cptr->cnext;
      cptr->cnext   = actptr;
      actptr->cnext = tmpptr;
      success++;
      goto endcol;
    }
    /* dis not find the proper place */
    cptr = cptr->cnext;
    goto startcol;
endcol:;
  }
  if (!success) dserror("Failed to set entry");
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_addval */




/*!---------------------------------------------------------------------
\brief adds a row in an oll matrix

<pre>                                                        mn 02/03
This function adds a row of values in an oll matrix. If one matrix entry
already exists, the values are added together, otherwise a new entry is
fetched.
</pre>
\param *matrix    OLL    (i/o)   the matrix
\param  actrow    INT      (i)   the true row index
\param  lm        INT[]    (i)   the location vactor
\param  val       DOUBLE[] (i)   the values to set
\param  nd        INT      (i)   the length of the row

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_addrow(
    OLL*  matrix,
    INT actrow,
    INT lm[],
    DOUBLE val[],
    INT nd)
{
  INT         i;
  MATENTRY  **row;
  MATENTRY  **col;
  INT         rindex;
  MATENTRY   *actptr;
  MATENTRY   *tmpptr;
  MATENTRY   *cptr;
  INT         success=0;
  INT         lastcol;
  INT         actcol;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_addrow");
#endif
  /*----------------------------------------------------------------------*/
  row = matrix->row;
  col = matrix->col;
  /*------------------------------------------------- get the correct row */
  rindex = find_index(actrow,matrix->update.a.iv,matrix->numeq);
  if (rindex==-1) dserror("row not in update");
  /* --------------------------check wether the mask of the matrix exists */
  if (matrix->is_masked)
    /* ---------------------------------------------- mask already exists */
  {
    /* ----------------------------------------------for all element dofs */
    lastcol=-1;
    actptr           = row[rindex];
    for(i=0; i<nd; i++)
    {
      actcol  = lm[i];
      if( actcol >= matrix->cdim) continue;
      if(actcol < lastcol)
        actptr           = row[rindex];
startrow1:
      if (actptr->c == actcol) /* this is correct place */
      {
        dsassert(actptr->r==actrow,"Matrix got mixed up");
        actptr->val += val[i];
        lastcol      = actcol;
        goto endrow1;
      }
      /* did not find the proper place */
      if(actptr->rnext == NULL)
        dserror("Failed to find correct place");
      else
        actptr = actptr->rnext;
      goto startrow1;
endrow1:;
    } /* -- END OF for all element dofs */
  } /* -- END OF mask already exists */
  else
    /* ------------------------------------------------ mask does not exist */
  {
    /* ----------------------------------------------for all element dofs */
    for(i=0; i<nd; i++)
    {
      actcol  = lm[i];
      if( actcol >= matrix->cdim) continue;
      /*------------------------------------------ check whether row is empty */
      if (row[rindex]==NULL)
      {
        oll_getentry(matrix,&row[rindex]);
        actptr           = row[rindex];
        actptr->r     = actrow;
        actptr->c     = actcol;
        actptr->val   = val[i];
      } /* END OF row is empty */

      else /*----------------------------------------------- row is not empty */
      {
        actptr           = row[rindex];
startrow:
        if (actptr->c > actcol) /* first entry larger then actcol */
        {
          dsassert(actptr==row[rindex],"Matrix got severly mixed up");
          tmpptr        = actptr;
          /*   row[rindex]   = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
          oll_getentry(matrix,&row[rindex]);
          actptr        = row[rindex];
          actptr->r     = actrow;
          actptr->c     = actcol;
          actptr->val   = val[i];
          actptr->rnext = tmpptr;
          goto endrow;
        }
        if (actptr->c == actcol) /* this is correct place */
        {
          dsassert(actptr->r==actrow,"Matrix got mixed up");
          actptr->val += val[i];
          goto endrow;
        }
        if (actptr->rnext == NULL) /* next element not exist */
        {
          /*actptr->rnext = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
          oll_getentry(matrix,&actptr->rnext);
          actptr        = actptr->rnext;
          actptr->r     = actrow;
          actptr->c     = actcol;
          actptr->val   = val[i];
          goto endrow;
        }
        if (actptr->rnext->c == actcol) /* next element correct element */
        {
          actptr        = actptr->rnext;
          dsassert(actptr->c==actcol && actptr->r==actrow,"Matrix got mixed up");
          actptr->val  += val[i];
          goto endrow;
        }
        if (actptr->rnext->c > actcol &&
            actptr->c        < actcol ) /* next element larger and active smaller */
        {
          tmpptr        = actptr->rnext;
          /*actptr->rnext = (MATENTRY*)CCACALLOC(1,sizeof(MATENTRY));*/
          oll_getentry(matrix,&actptr->rnext);
          actptr        = actptr->rnext;
          actptr->rnext = tmpptr;
          actptr->r     = actrow;
          actptr->c     = actcol;
          actptr->val   = val[i];
          goto endrow;
        }
        /* did not find the proper place */
        actptr = actptr->rnext;
        goto startrow;
endrow:;
      } /* -- END OF row is not empty */
      /*---------------------------------------------- check dimension of col */
      if (actcol >= matrix->cdim) dserror("Overvlow in column dimension");
      /*--------------------------------------- check whether column is empty */
      if (col[actcol]==NULL)
      {
        col[actcol] = actptr;
        success++;
      }
      else /*-------------------------------------------- column is not empty */
      {
        cptr = col[actcol];
startcol:
        if (cptr->r >  actrow)
        {
          dsassert(cptr==col[actcol],"Matrix got severly mixed up");
          dsassert(actptr->cnext==NULL,"Matrix got severly mixed up");
          tmpptr        = cptr;
          col[actcol]   = actptr;
          actptr->cnext = tmpptr;
          success++;
          goto endcol;
        }
        if (cptr->r == actptr->r) /* this entry is already linked in */
        {
          success++;
          goto endcol;
        }
        if (cptr->cnext==NULL) /* next entry is not assigned */
        {
          cptr->cnext   = actptr;
          success++;
          goto endcol;
        }
        if (cptr->r < actrow && cptr->cnext->r > actrow) /* active smaller, next larger */
        {
          tmpptr        = cptr->cnext;
          cptr->cnext   = actptr;
          actptr->cnext = tmpptr;
          success++;
          goto endcol;
        }
        /* dis not find the proper place */
        cptr = cptr->cnext;
        goto startcol;
endcol:;
      }
      if (!success) dserror("Failed to set entry");
    } /* -- END OF for all element dofs */
  } /* -- END OF mask does not exist */
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_addrow */



/*!---------------------------------------------------------------------
\brief copies an oll matrix to the skyline format

<pre>                                                        mn 02/03
This function copies an oll matrix to a matrix in skyline format. The maxa
vector is already in fortran convention.
</pre>
\param *oll       OLL          (i)   the oll matrix
\param *sysarray  SPARSE_ARRAY (o)   the sysarry matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_to_sky(
    OLL*  oll,
    SPARSE_ARRAY*  sysarray)
{

  SKYMATRIX      *sky;

  INT            *maxaf;
  DOUBLE         *a;

  INT             adim;

  MATENTRY     **row;
  MATENTRY     **col;
  MATENTRY      *actentry;

  INT            i,j;
  INT            counter,height,lastrow,gap;

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_to_sky");
#endif
  /*----------------------------------------------------------------------*/
  sky = sysarray->sky;

  sky->nnz = oll->nnz;
  sky->numeq_total = oll->numeq_total;
  sky->numeq = oll->numeq;
  sky->is_init = 1;
  sky->is_factored = 0;
  sky->ncall = 0;

  am_alloc_copy(&(oll->update),&(sky->update));

  /* --------------------------------------------- create skyline mtrix */
  row = oll->row;
  col = oll->col;
  /*----------------------------------------------------- allocate maxaf */
  maxaf = amdef("maxaf",&(sky->maxaf),sky->numeq_total+1,1,"IV");
  /*---------------------------------------------------------allocate A */
  /* calculate size of a */
  adim = 0;
  for(i=0; i<oll->cdim; i++)
  {
    adim += (i - col[i]->r + 1 );
  }
  a = amdef("A",&(sky->A),adim,1,"DV");

  /* ------------------------------ copy oll matrix into skyline format */
  counter = 0;
  for (i=0; i<oll->cdim; i++)
  {
    /* -------------------------------- get the first entry of column i */
    actentry = col[i];
    /* -------------------------------- calculate height of this column */
    height = i - actentry->r + 1;
    /* ------------------ counter is the postion of the diagonal entry  */
    maxaf[i]= counter+1;
    /* ------------------------------------find position of first entry */
    counter += (height-1);
    /* -------------------- now copy all entries in this column into A, *
     * ---------------------------------until diagonal entry is reached */
    lastrow = actentry->r-1;
    while (actentry!=NULL && actentry->r<=i)
    {
      gap = actentry->r - lastrow;
      if(gap>1)
      {
        /* there were (gap-1) zero entries before this entry */
        for (j=0;j<gap-1;j++)
        {
          a[counter] = 0.0;
          counter--;
        }
      }
      a[counter] = actentry->val;
      lastrow = actentry->r;
      counter--;
      actentry = actentry->cnext;
    } /* end of column i */
    counter+= height+1;
  } /* end of the matrix */
  maxaf[i] = counter + 1;

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_to_sky */




/*!---------------------------------------------------------------------
\brief copies an oll matrix to the spooles format

<pre>                                                        mn 02/03
This function copies an oll matrix to a matrix in spooles format. The
values are not copied!!
</pre>
\param *oll       OLL          (i)   the oll matrix
\param *sysarray  SPARSE_ARRAY (o)   the sysarry matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_to_spo(
    OLL*  oll,
    SPARSE_ARRAY*  sysarray)
{

  SPOOLMAT      *spo;

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_to_spo");
#endif
  /*----------------------------------------------------------------------*/
  spo = sysarray->spo;

  spo->nnz = oll->nnz;
  spo->numeq_total = oll->numeq_total;
  spo->numeq = oll->numeq;
  spo->is_init = 1;
  spo->is_factored = 0;
  spo->ncall = 0;

  am_alloc_copy(&(oll->update),&(spo->update));

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_to_spo */



/*!---------------------------------------------------------------------
\brief copies an oll matrix to the msr format

<pre>                                                        mn 02/03
This function copies an oll matrix to a matrix in msr format.
</pre>
\param *oll       OLL          (i)   the oll matrix
\param *sysarray  SPARSE_ARRAY (o)   the sysarry matrix

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ---

------------------------------------------------------------------------*/
void oll_to_msr(
    OLL*  oll,
    SPARSE_ARRAY*  sysarray)
{

  AZ_ARRAY_MSR   *msr;
  INT            *bindx;
  DOUBLE         *val;
  INT            *bindx_b;
  DOUBLE         *val_b;


  MATENTRY     **row;
  MATENTRY     **col;
  MATENTRY      *actentry;

  INT            i;
  INT            rindex;
  INT            nnz, numeq, counter;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("oll_to_msr");
#endif
  /*----------------------------------------------------------------------*/
  msr     = sysarray->msr;

  msr->nnz = oll->nnz;
  msr->numeq_total = oll->numeq_total;
  msr->numeq = oll->numeq;
  msr->is_init = 1;
  msr->is_factored = 0;
  msr->ncall = 0;

  am_alloc_copy(&(oll->update),&(msr->update));

  /* ------------------------------------------------- create msr matrix */
  row   = oll->row;
  col   = oll->col;
  nnz   = oll->nnz;
  numeq = oll->numeq;
  /*---------------------------------------------- allocate bindx and val */
  bindx     = amdef("bindx",&(msr->bindx)       ,(msr->nnz+1),1,"IV");
  val       = amdef("val"  ,&(msr->val)         ,(msr->nnz+1),1,"DV");
  bindx_b   = amdef("bindx",&(msr->bindx_backup),(msr->nnz+1),1,"IV");
  val_b     = amdef("val"  ,&(msr->val_backup)  ,(msr->nnz+1),1,"DV");

  /* ----------------------------------- copy oll matrix into msr format */
  counter = numeq+1;
  for (i=0; i<oll->rdim; i++)
  {
    /* -------------------- write index of beginning of row in bindx */
    bindx[i] = counter;
    /* -------------------------------- get the first entry of row i */
    actentry = row[i];
    while (actentry!=NULL)
    {
      /*if(oll->update.a.iv[actentry->r]==actentry->c)*/
      if(actentry->r==actentry->c)
      {
        /* --------------------------- this is a diagonal entry */
        /* -------------------------------- get the correct row */
        rindex = find_index(actentry->r,oll->update.a.iv,oll->numeq);
        if (rindex==-1) dserror("row not in update");
        val[rindex] = actentry->val;
      }
      else
      {
        /* --------------------------- this is no diagonal entry */
        val[counter]   = actentry->val;
        bindx[counter]  = actentry->c;
        counter++;
      }
      actentry = actentry->rnext;
    } /* END OF while row */
  } /* END OF for all rows */
  bindx[i] = counter;
  val[i]   = 0.0;

  if(counter!=msr->nnz+1)
    dserror("Error in copying oll to msr");

  /*--------- make backup copy of bindx, as it is permuted in solution */
  am_alloc_copy(&(msr->bindx),&(msr->bindx_backup));
  am_alloc_copy(&(msr->val),&(msr->val_backup));

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of oll_to_msr */


/*!---------------------------------------------------------------------
\brief copies an oll matrix to the ccf format

<pre>                                                       genk 10/03

This function copies an oll matrix to a matrix in ccf format.

</pre>
\param *oll       OLL          (i)   the oll matrix
\param *sysarray  SPARSE_ARRAY (o)   the sysarry matrix
\return void

------------------------------------------------------------------------*/
void oll_to_ccf(
                OLL  *oll,
                SPARSE_ARRAY *sysarray
	       )
{
CCF	   *ccf;         /* a sparse matrix in compressed column format */
INT	   *Ap;          /* column pointer vector			*/
INT	   *Ai;          /* row pointer vector			        */
DOUBLE     *Ax;          /* values of the matrix			*/
MATENTRY **col;          /* matrix column                               */
MATENTRY  *actentry;     /* actual matrix entry                         */
INT	   i;            /* simply a counter                            */
INT	   nnz;          /* number of nonzero matrix entries            */
INT        numeq;        /* number of equations on this proc            */
INT        numeq_total;  /* total number of equations                   */
INT        counter;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("oll_to_ccf");
#endif
#ifdef PARALLEL
dserror("No UMFPACK for parallel OLL!\n");
#endif
/*----------------------------------------------------------------------*/
ccf     = sysarray->ccf;
/* -------------------------------------------------- create ccf matrix */
col   = oll->col;
nnz   = oll->nnz;
numeq = oll->numeq;
numeq_total = oll->numeq_total;

/*--------------------------------------------- redefine matrix vectors */
Ap = amredef(&(oll->sysarray[0].ccf->Ap),numeq_total+1,1,"IV");
Ai = amredef(&(oll->sysarray[0].ccf->Ai),nnz 	 ,1,"IV");
Ax = amredef(&(oll->sysarray[0].ccf->Ax),nnz 	 ,1,"DV");
amredef(&(oll->sysarray[0].ccf->update),numeq,1,"IV");

/*-------------------------------------------------- copy update vector */
amcopy(&(oll->update),&(oll->sysarray[0].ccf->update));

/* ------------------------------------ copy oll matrix into ccf format */
counter = 0;
Ap[0]=0;
for (i=0; i<oll->cdim; i++)
{
   /* --------------------------------- get the first entry of column i */
   actentry = col[i];
   while (actentry!=NULL)
   {
      Ax[counter] = actentry->val;
      Ai[counter] = actentry->r;
      counter++;
      actentry = actentry->cnext;
   } /* END OF while col */
   Ap[i+1]=counter;
} /* END OF for all cols */

if (counter!=oll->nnz)
dserror("Error in copying oll to ccf!\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of oll_to_msr */

/*!---------------------------------------------------------------------
\brief copy a matrix in oll format

<pre>                                                        mn 03/03
copy a matrix in oll format, the matrix to has to be correcly opened
using oll_open
</pre>
\param from    OLL*    (i/o)   the matrix copied from
\param to      OLL*    (i/o)   the matrix copied to
\return void
\sa  mlpcg_ll_open
------------------------------------------------------------------------*/
void oll_copy(
              OLL   *from,
              OLL   *to)
{
INT        i;
MATENTRY  *actptr;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("oll_copy");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------------------------- check dimensions */
if (from->rdim != to->rdim || from->cdim != to->cdim) dserror("Mismatch of dimensions");
/*------------------------------------------------------ make copy loop */
for (i=0; i<from->rdim; i++)
{
   actptr = from->row[i];
   while (actptr != NULL)
   {
      oll_setval(to,actptr->r,actptr->c,actptr->val);
      actptr = actptr->rnext;
   }
}

to->is_masked = 1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of oll_copy */

/*! @} (documentation module close)*/
