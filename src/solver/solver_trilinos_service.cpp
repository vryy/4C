/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef TRILINOS_PACKAGE

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include "mrtr_utils.H"


#include "../solver/solver_trilinos_service.H"

INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
static void make_update(
    FIELD           *actfield,
    PARTITION       *actpart,
    SOLVAR          *actsolv,
    INTRA           *actintra,
    TRILINOSMATRIX  *tri,
    int              disnum
    );
#ifdef PARALLEL
static void add_tri_checkcouple(
    const int   ii,
    const int** cdofs,
    const int   ncdofs,
    int*        iscouple,
    int*        isowner,
    const int   nprocs
    );
#endif
static void add_tri_sendbuff(
    const int ii,
    const int jj,
    const int i,
    const int j,
    const int ii_owner,
    int**     isend,
    double**  dsend,
    double**  estif,
    const int numsend
    );


using namespace std;


/*----------------------------------------------------------------------*/
/*!
  \brief construct the map of an Epetra_CrsMatrix from known update vector

  \param actintra       (i)  the intra-communicator we do not need
  \param trimatrix    (i/o)  the matrix object to set up

  \author kuettler
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void construct_trilinos_matrix(
                            INTRA           *actintra,
                            TRILINOSMATRIX  *trimatrix
                         )
{
  DSTraceHelper dst("construct_trilinos_matrix");

  // create an epetra comm
#ifdef PARALLEL
  Epetra_MpiComm*    comm = new Epetra_MpiComm(actintra->MPI_INTRA_COMM);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif
  trimatrix->epetracomm = (void*)comm;

  // create an Epetra_Map
  int nglobal = trimatrix->numeq_total;
  int nlocal  = trimatrix->numeq;
  int* update = trimatrix->update.a.iv;
  Epetra_Map* map = new Epetra_Map(nglobal,nlocal,update,0,*comm);
  trimatrix->rowmap = (void*)map;

  // Allocate the Epetra_CrsMatrix
  trimatrix->NumEntriesPerRow = 81;
  trimatrix->StaticProfile = false;
  Epetra_CrsMatrix* matrix = new Epetra_CrsMatrix(Copy,*map,81,false);
  trimatrix->matrix = (void*)matrix;

  // set flags indicating status of matrix
  trimatrix->is_init=1;

  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief construct the map of a diagonal Epetra_CrsMatrix from known update vector

  \param actintra       (i)  the intra-communicator we do not need
  \param trimatrix    (i/o)  the matrix object to set up

  \author kuettler
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void construct_trilinos_diagonal_matrix(INTRA *actintra,
                                        TRILINOSMATRIX *trimatrix)
{
  DSTraceHelper dst("construct_trilinos_diagonal_matrix");

  // create an epetra comm
#ifdef PARALLEL
  Epetra_MpiComm*    comm = new Epetra_MpiComm(actintra->MPI_INTRA_COMM);
#else
  Epetra_SerialComm* comm = new Epetra_SerialComm();
#endif
  trimatrix->epetracomm = (void*)comm;

  // create an Epetra_Map
  int nglobal = trimatrix->numeq_total;
  int nlocal  = trimatrix->numeq;
  int* update = trimatrix->update.a.iv;
  Epetra_Map* map = new Epetra_Map(nglobal,nlocal,update,0,*comm);
  trimatrix->rowmap = (void*)map;

  // Allocate the Epetra_CrsMatrix
  trimatrix->NumEntriesPerRow = 1;
  trimatrix->StaticProfile = true;
  Epetra_CrsMatrix* matrix = new Epetra_CrsMatrix(Copy,*map,1,true);
  trimatrix->matrix = (void*)matrix;

  // set flags indicating status of matrix
  trimatrix->is_init=1;

  return;
}


/*----------------------------------------------------------------------*
  |  calculate the map of an Epetra_CrsMatrix from a discretization     |
  |                                                        m.gee 9/06   |
 *----------------------------------------------------------------------*/
void init_trilinos_matrix(
                            FIELD           *actfield,
                            PARTITION       *actpart,
                            SOLVAR          *actsolv,
                            INTRA           *actintra,
                            TRILINOSMATRIX  *trimatrix,
                            int              disnum
                         )
{
  DSTraceHelper dst("init_trilinos_matrix");

  // global number of dofs
  trimatrix->numeq_total = actfield->dis[disnum].numeq;
  // local number of dofs
  mask_numeq(actfield,actpart,actsolv,actintra,&(trimatrix->numeq),disnum);

  // allocate an update vector and fill it
  amdef("update",&(trimatrix->update),trimatrix->numeq,1,"IV");
  amzero(&(trimatrix->update));
  make_update(actfield,actpart,actsolv,actintra,trimatrix,disnum);

  construct_trilinos_matrix(actintra,trimatrix);

  return;
}

/*----------------------------------------------------------------------*
  |  allocate update put dofs in update in ascending order  m.gee 9/06 |
 *----------------------------------------------------------------------*/
void make_update(
    FIELD           *actfield,
    PARTITION       *actpart,
    SOLVAR          *actsolv,
    INTRA           *actintra,
    TRILINOSMATRIX  *tri,
    int              disnum
    )

{
  INT       i,k,l;
  INT       counter;
  INT      *update;
  INT       dof;
  INT       foundit;
  INT       imyrank;
  INT       inprocs;
  NODE     *actnode;
  ARRAY     coupledofs;

  DSTraceHelper dst("make_update");

  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  memset(&coupledofs, 0, sizeof(ARRAY));
  if (actpart->pdis[disnum].coupledofs.Typ != ARRAY::cca_XX)
    am_alloc_copy(&(actpart->pdis[disnum].coupledofs),&coupledofs);
  /*----------------------------------------------------------------------*/
  update = tri->update.a.iv;
  counter=0;
  /*------------------------------------- loop the nodes on the partition */
  for (i=0; i<actpart->pdis[disnum].numnp; i++)
  {
    actnode = actpart->pdis[disnum].node[i];
    for (l=0; l<actnode->numdf; l++)
    {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->dis[disnum].numeq) continue;
      /* no coupling on dof */
      if (actnode->gnode->couple==NULL)
      {
        update[counter] = dof;
        counter++;
        continue;
      }
      else /* coupling on node */
      {
        foundit=0;
        /* find dof in coupledofs */
        for (k=0; k<coupledofs.fdim; k++)
        {
          if (dof == coupledofs.a.ia[k][0])
          {
            /* am I owner of this dof or not */
            if (coupledofs.a.ia[k][imyrank+1]==2)
              foundit=2;
            else if (coupledofs.a.ia[k][imyrank+1]==1)
              foundit=1;
            break;
          }
        }
        /* dof found in coupledofs */
        if (foundit==2)/* I am master owner of this coupled dof */
        {
          update[counter] = dof;
          counter++;
          coupledofs.a.ia[k][imyrank+1]=1;
          continue;
        }
        else if (foundit==1)/* I am slave owner of this coupled dof */
        {
          /* do nothing, this dof doesn't exist for me (no more)*/
        }
        else /* this dof is not a coupled one */
        {
          update[counter] = dof;
          counter++;
          continue;
        }
      }

    }
  }
  /*----------- check whether the correct number of dofs has been counted */
  if (counter != tri->numeq) dserror("Number of dofs in vector update wrong");
  /*---------------------------- sort the vector update just to make sure */
  qsort((INT*) update, counter, sizeof(INT), cmp_int);
  /*----------------------------------------------------------------------*/
  if (coupledofs.fdim > 0)
    amdel(&coupledofs);
  /*----------------------------------------------------------------------*/

  return;
} /* end of make_update */


/*----------------------------------------------------------------------*
  |  copy a trilinos matrix  mask                            m.gee 9/06 |
 *----------------------------------------------------------------------*/
void trilinos_cp_matrixmask(TRILINOSMATRIX  *from, TRILINOSMATRIX* to)
{
  DSTraceHelper dst("trilinos_cp_matrixmask");

  if (!(from->epetracomm)  || !(from->rowmap) || !(from->matrix))
    dserror("TRILINOSMATRIX from not properly initialized");

  am_alloc_copy(&(from->update),&(to->update));

  // Do Epetra_Comm
#ifdef PARALLEL
  Epetra_MpiComm*    comm = (Epetra_MpiComm*)from->epetracomm;
  Epetra_MpiComm*    newcomm = new Epetra_MpiComm(*comm);
#else
  Epetra_SerialComm* comm = (Epetra_SerialComm*)from->epetracomm;
  Epetra_SerialComm* newcomm = new Epetra_SerialComm(*comm);
#endif
  to->epetracomm = (void*)newcomm;

  // Do Epetra_Map
  int nglobal = from->numeq_total;
  int nlocal  = from->numeq;
  int* update = to->update.a.iv;
  Epetra_Map *rowmap = new Epetra_Map(nglobal,nlocal,update,0,*newcomm);
  to->rowmap = (void*)rowmap;

  // Do Epetra_CrsMatrix
  Epetra_CrsMatrix* matrix = new Epetra_CrsMatrix(Copy,*rowmap,
                                                  from->NumEntriesPerRow,
                                                  from->StaticProfile);
  to->matrix = (void*)matrix;

  // do dimensions and flags
  to->is_init     = 1;
  to->numeq_total = from->numeq_total;
  to->numeq       = from->numeq;

  return;
} /* end of trilinos_cp_matrixmask */


/*----------------------------------------------------------------------*
  |  zero a trilinos matrix                                  m.gee 9/06 |
 *----------------------------------------------------------------------*/
void trilinos_zero_matrix(TRILINOSMATRIX *tri)
{
  DSTraceHelper dst("trilinos_zero_matrix");

  if (!(tri->epetracomm)  || !(tri->rowmap) || !(tri->matrix))
    dserror("TRILINOSMATRIX tri not properly initialized");

  // destroy the old matrix
  if (tri->matrix)
  {
    Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)(tri->matrix);
    delete matrix;
  }

  // get the Epetra_Map
  Epetra_Map *rowmap = (Epetra_Map*)(tri->rowmap);

  // Do Epetra_CrsMatrix
  Epetra_CrsMatrix* matrix = new Epetra_CrsMatrix(Copy,*rowmap,
                                                  tri->NumEntriesPerRow,
                                                  tri->StaticProfile);
  tri->matrix = (void*)matrix;

  // do flag
  tri->is_factored = 0;

  // See if matrix was used with spooles
#ifdef SPOOLES_PACKAGE
  if (tri->sysarray_typ)
  if (tri->sysarray_typ[0]==spoolmatrix)
  if (tri->sysarray)
  if (tri->sysarray->spo)
  {
    SPOOLMAT* spo = tri->sysarray[0].spo;
    spo->is_factored=0;
    if (spo->ncall > 0)
    {
      FrontMtx_free(spo->frontmtx);
      InpMtx_free(spo->newA);
      DenseMtx_free(spo->newY);
      ETree_free(spo->frontETree);
      SubMtxManager_free(spo->mtxmanager);
      IV_free(spo->newToOldIV);
      IV_free(spo->oldToNewIV);
      IV_free(spo->ownersIV);
      IV_free(spo->vtxmapIV);
      IV_free(spo->ownedColumnsIV);
      SolveMap_free(spo->solvemap);
      IVL_free(spo->symbfacIVL);
    }
  }
#endif

  return;
} /* end of trilinos_zero_matrix */


/*----------------------------------------------------------------------*/
/*!
  \brief add one value to a trilinos matrix

  The matrix must not be completed.

  \param tri  (i) matrix
  \param v    (i) value to add
  \param row  (i) global row number
  \param col  (i) global col number

  \author kuettler
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void add_trilinos_value(struct _TRILINOSMATRIX *tri, DOUBLE v, INT row, INT col)
{
  DSTraceHelper dst("add_trilinos_value");

  Epetra_CrsMatrix* mat = (Epetra_CrsMatrix*)tri->matrix;
  if (mat->Filled())
    dserror("Epetra_CrsMatrix::FillComplete() has been called, cannot assemble anymore");
  if (mat->IndicesAreLocal())
    dserror("IndicesAreLocal()");

  int err = mat->SumIntoGlobalValues(row,1,&v,&col);
  if (err)
    err = mat->InsertGlobalValues(row,1,&v,&col);
  if (err<0)
    dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code %d",err);

}


/*----------------------------------------------------------------------*/
/*!
  \brief assemble into a Epetra_CrsMatrix matrix

  This routine assembles one or two element matrices (elearray1 and
  elearray2) into the global matrices in the Epetra format.

  \param actpart   *PARTITION    (i)  the partition we are working on
  \param actsolv   *SOLVAR       (i)  the solver we are using
  \param actintra  *INTRA        (i)  the intra-communicator we do not need
  \param actele    *ELEMENT      (i)  the element we would like to work with
  \param msr1      *TRILINOSMATRIX (i)  one sparse matrix we will assemble into
  \param msr2      *TRILINOSMATRIX (i)  the other sparse matrix we will assemble into

  \author gee
  \date 09/06

*/
/*----------------------------------------------------------------------*/
void  add_trilinos(
    struct _PARTITION       *actpart,
    struct _SOLVAR          *actsolv,
    struct _INTRA           *actintra,
    struct _ELEMENT         *actele,
    struct _TRILINOSMATRIX  *tri1,
    struct _TRILINOSMATRIX  *tri2,
    struct _ARRAY           *elearray1,
    struct _ARRAY           *elearray2
    )

{
  DSTraceHelper dst("add_trilinos");

  /* set some pointers and variables */
  const int myrank      = actintra->intra_rank;
  double** estif        = elearray1->a.da;
  double** emass        = NULL;
  if (tri2) emass       = elearray2->a.da;
  const int numeq_total = tri1->numeq_total;
  int lm[MAXDOFPERELE]; /* location vector for this element */
#ifdef PARALLEL
  int         owner[MAXDOFPERELE];  /* the owner of every dof */
  const int   nprocs = actintra->intra_nprocs;
  const int** cdofs  = (const int**)actpart->pdis[0].coupledofs.a.ia;
  const int   ncdofs = (const int)actpart->pdis[0].coupledofs.fdim;
#endif
  int    **isend1 = NULL;        /* p to sendbuffer to communicate coupling cond */
  double **dsend1 = NULL;        /* p to sendbuffer to communicate coupling cond */
  int    **isend2 = NULL;        /* p to sendbuffer to communicate coupling cond */
  double **dsend2 = NULL;        /* p to sendbuffer to communicate coupling cond */
  int      nsend  = 0;
  Epetra_CrsMatrix* mat1 = (Epetra_CrsMatrix*)tri1->matrix;
  if (mat1->Filled()) dserror("Epetra_CrsMatrix::FillComplete() has been called, cannot assemble anymore");

  Epetra_CrsMatrix* mat2 = NULL;
  if (tri2)
  {
    mat2 = (Epetra_CrsMatrix*)tri2->matrix;
    if (mat2->Filled()) dserror("Epetra_CrsMatrix::FillComplete() has been called, cannot assemble anymore");
  }

  /* put pointers to sendbuffers if any */
#ifdef PARALLEL
  if (tri1->couple_i_send)
  {
    isend1 = tri1->couple_i_send->a.ia;
    dsend1 = tri1->couple_d_send->a.da;
    nsend  = tri1->couple_i_send->fdim;
    if (tri2)
    {
      isend2 = tri2->couple_i_send->a.ia;
      dsend2 = tri2->couple_d_send->a.da;
    }
  }
#endif

  /* make location vector lm*/
  int counter=0;
  const int numnp = actele->numnp;
  for (int i=0; i<numnp; ++i)
  {
    for (int j=0; j<actele->node[i]->numdf; ++j)
    {
      lm[counter]    = actele->node[i]->dof[j];
#ifdef PARALLEL
      owner[counter] = actele->node[i]->proc;
#endif
      ++counter;
    }
  }
  const int nd = counter; // no. degrees of freedom on this element <-> size of element matrices

  /* now start looping the dofs */
  /* loop over i (the element row) */
  int ii_iscouple = 0;
  int ii_owner    = myrank;
  for (int i=0; i<nd; ++i)
  {
    const int ii = lm[i];
    /* loop only my own rows */
#ifdef PARALLEL
    if (owner[i]!=myrank) continue;
#endif

    /* check for boundary condition */
    if (ii>=numeq_total) continue;

    /* check for coupling condition */
#ifdef PARALLEL
    if (ncdofs)
    {
      ii_iscouple = 0;
      ii_owner    = -1;
      add_tri_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
    }
#endif

    /* loop over j (the element column) */
    /* This is the full unsymmetric version ! */
    for (int j=0; j<nd; ++j)
    {
      int jj = lm[j];

      /* check for boundary condition */
      if (jj>=numeq_total) continue;

      /* (either not a coupled dof or I am master owner) */
      if (!ii_iscouple || ii_owner==myrank)
      {
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
        // do stiffness matrix (nonzero values only)
        int err = 0;
        if (abs(estif[i][j])>EPS10 || i==j)
        {
          err = mat1->SumIntoGlobalValues(ii,1,&(estif[i][j]),&jj);
          if (err)
            err = mat1->InsertGlobalValues(ii,1,&(estif[i][j]),&jj);
          if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code");
        }
        // do mass matrix if present (nonzero values only)
        if (mat2)
        if (abs(emass[i][j])>EPS10 || i==j)
        {
          err = mat2->SumIntoGlobalValues(ii,1,&(emass[i][j]),&jj);
          if (err)
            err = mat2->InsertGlobalValues(ii,1,&(emass[i][j]),&jj);
          if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code");
        }
#else
        // do stiffness matrix (all values)
        int err = mat1->SumIntoGlobalValues(ii,1,&(estif[i][j]),&jj);
        if (err)
          err = mat1->InsertGlobalValues(ii,1,&(estif[i][j]),&jj);
        if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code");

        // do mass matrix if present (all values)
        if (mat2)
        {
          err = mat2->SumIntoGlobalValues(ii,1,&(emass[i][j]),&jj);
          if (err)
            err = mat2->InsertGlobalValues(ii,1,&(emass[i][j]),&jj);
          if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code");
        }
#endif
      }
      /* (a coupled dof and I am slave owner) */
      else
      {
        add_tri_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
        if (mat2)
          add_tri_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }
    } /* end loop over j */
  }/* end loop over i */


  return;
} /* end of add_trilinos */


#ifdef PARALLEL
/*----------------------------------------------------------------------*
  |  checks coupling for the add_trilinos routine            m.gee 9/06 |
 *----------------------------------------------------------------------*/
void add_tri_checkcouple(
    const int   ii,
    const int** cdofs,
    const int   ncdofs,
    int*        iscouple,
    int*        isowner,
    const int   nprocs
    )

{
  int         i,k;
  DSTraceHelper dst("add_tri_checkcouple");

  /*----------------------------------------------------------------------*/
  for (k=0; k<ncdofs; k++)
  {
    if (ii==cdofs[k][0])
    {
      *iscouple=1;
      for (i=1; i<=nprocs; i++)
      {
        if (cdofs[k][i]==2)
        {
          *isowner=i-1;
          break;
        }
      }
    }
  }
  /*----------------------------------------------------------------------*/

  return;
} /* end of add_tri_checkcouple */
#endif // PARALLEL


/*----------------------------------------------------------------------*
  |  fill sendbuffer isend and dsend                         m.gee 9/06|
 *----------------------------------------------------------------------*/
void add_tri_sendbuff(
    const int      ii,
    const int      jj,
    const int      i,
    const int      j,
    const int      ii_owner,
    int**          isend,
    double**       dsend,
    double**       estif,
    const int      numsend
    )

{
  DSTraceHelper dst("add_tri_sendbuff");
  /*----------------------------------------------------------------------*/
  int k;
  for (k=0; k<numsend; ++k)
  {
    if (isend[k][0]==ii) break;
  }
  isend[k][1]  = ii_owner;
  dsend[k][jj]+= estif[i][j];
  /*----------------------------------------------------------------------*/

  return;
} /* end of add_tri_sendbuff */


/*----------------------------------------------------------------------*
  |  exchange coupled dofs and add to epetra matrix           m.gee 9/06|
 *----------------------------------------------------------------------*/
void exchange_coup_trilinos(
    PARTITION*      actpart,
    SOLVAR*         actsolv,
    INTRA*          actintra,
    TRILINOSMATRIX* tri
    )

{
  DSTraceHelper dst("exchange_coup_trilinos");

#ifdef PARALLEL
  int            tag;
  int            source;
  int            numeq,numeq_total;
  int            numsend;
  int            numrecv;
  int           *update;
  int          **isend = NULL;
  double       **dsend = NULL;
  int          **irecv = NULL;
  double       **drecv = NULL;
  int            imyrank;
  int            inprocs;

  MPI_Status    *irecv_status = NULL;
  MPI_Status    *drecv_status = NULL;

  MPI_Request   *isendrequest = NULL;
  MPI_Request   *dsendrequest = NULL;

  MPI_Comm      *ACTCOMM;
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  ACTCOMM = &(actintra->MPI_INTRA_COMM);
  /*---------------------------------------- set some pointers and values */
  numsend     = tri->numcoupsend;
  numrecv     = tri->numcouprecv;
  update      = tri->update.a.iv;
  numeq_total = tri->numeq_total;
  numeq       = tri->numeq;
  if (tri->couple_i_send) isend   = tri->couple_i_send->a.ia;
  if (tri->couple_d_send) dsend   = tri->couple_d_send->a.da;
  if (tri->couple_i_recv) irecv   = tri->couple_i_recv->a.ia;
  if (tri->couple_d_recv) drecv   = tri->couple_d_recv->a.da;
  /*--------------------------------------------- allocate some envelopes */
  if (numrecv)
  {
    irecv_status = (MPI_Status*)CCACALLOC(numrecv,sizeof(MPI_Status));
    drecv_status = (MPI_Status*)CCACALLOC(numrecv,sizeof(MPI_Status));
    if (!irecv_status || !drecv_status) dserror("Allocation of memory failed");
  }
  if (numsend)
  {
    isendrequest = (MPI_Request*)CCACALLOC(numsend,sizeof(MPI_Request));
    dsendrequest = (MPI_Request*)CCACALLOC(numsend,sizeof(MPI_Request));
    if ( !isendrequest || !dsendrequest) dserror("Allocation of memory failed");
  }
  /*-------------------------------------------- loop the dofs to be send */
  /* do all non-blocking sends and don't care about arrival (wird scho' klappe)*/
  /*     the only thing to care for is the order in which things are send */
  for (int i=0; i<numsend; ++i)
  {
    /*            sendbuffer       lenght    typ        dest.        tag          comm      request-handle */
    MPI_Isend(&(isend[i][0]),          2,MPI_INT   ,isend[i][1],isend[i][0],(*ACTCOMM),&(isendrequest[i]));
    MPI_Isend(&(dsend[i][0]),numeq_total,MPI_DOUBLE,isend[i][1],isend[i][0],(*ACTCOMM),&(dsendrequest[i]));
  }/*------------------------------------------------ end of sending loop */
  /*------------------------------- now loop over the dofs to be received */
  /*
     do blocking receives, 'cause one can't add something to the system
     matrix, which has not yet arrived, easy, isn't it?
     */
  if (!tri->matrix) dserror("tri->matrix is NULL");
  Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;
  if (matrix->Filled()) dserror("Epetra_CrsMatrix::FilComplete() was called on matrix, cannot assemble anymore");

  for (int i=0; i<numrecv; ++i)
  {
    /*--------------------------- use wildcards to receive first to come */
    /*          recv-buf  lenght typ     source           tag       comm              status */
    MPI_Recv(&(irecv[i][0]),2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,(*ACTCOMM),&(irecv_status[i]));
    if (irecv_status[i].MPI_ERROR) dserror("An error in MPI - communication occured !");

    /*---------------------- the dof number was sent as tag and as entry */
    tag    = irecv_status[i].MPI_TAG;
    if (tag != irecv[i][0]) dserror("MPI messages somehow got mixed up");
    source = irecv_status[i].MPI_SOURCE;

    /* do not use wildcards for second recv, we know now where it should come from */
    MPI_Recv(&(drecv[i][0]),numeq_total,MPI_DOUBLE,source,tag,(*ACTCOMM),&(drecv_status[i]));
    if (drecv_status[i].MPI_ERROR) dserror("An error in MPI - communication occured !");

    /* now add the received data properly to my own piece of sparse matrix */
    int ii = tag;
    //ii_index = AZ_quick_find(ii,update,numeq,shift,bins);
    //if (ii_index==-1) dserror("dof ii not found on this proc");
    // my rowmap
    const Epetra_Map& map = matrix->RowMap();
    // my local receive length
    const int recvlength = map.NumMyElements();
    for (int j=0; j<recvlength; ++j)
    {
      int jj = map.GID(j);
      if (jj<0) dserror("Cannot find global for local dof number on this proc");
      if (abs(drecv[i][jj])<EPS10) continue;
      int err = matrix->SumIntoGlobalValues(ii,1,&drecv[i][jj],&jj);
      if (err)
        err = matrix->InsertGlobalValues(ii,1,&drecv[i][jj],&jj);
      if (err<0)
        dserror("Epetra_CrsMatrix::InsertGlobalValues(...) returned an error");
    }
  }/*---------------------------------------------- end of receiving loop */
  /*-------------------------------------------- free allocated MPI-stuff */
  if (numrecv){CCAFREE(irecv_status);CCAFREE(drecv_status);}
  if (numsend){CCAFREE(isendrequest);CCAFREE(dsendrequest);}
  /*----------------------------------------------------------------------
    do a barrier, because this is the end of the assembly, the msr matrix
    is now ready for solve
    */
  MPI_Barrier(*ACTCOMM);
#endif /*---------------------------------------------- end of PARALLEL */
  /*----------------------------------------------------------------------*/
  return;
} /* end of exchange_coup_trilinos */





/*----------------------------------------------------------------------*
  |  finalize the assembly of a Epetra_CrsMatrix              m.gee 9/06|
 *----------------------------------------------------------------------*/
void close_trilinos_matrix(struct _TRILINOSMATRIX *tri)

{
  DSTraceHelper dsh("close_trilinos_matrix");
  /*----------------------------------------------------------------------*/
  if (!tri) return;
  Epetra_CrsMatrix* matrix = (Epetra_CrsMatrix*)tri->matrix;
  if (matrix->Filled()) return;

  if (!tri->rowmap) dserror("Matrix has no row map");
  Epetra_Map* map = (Epetra_Map*)tri->rowmap;

  int err = matrix->FillComplete(*map,*map);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domainmap,rowmap) returned an error");

  err = matrix->OptimizeStorage();
  if (err) dserror("Epetra_CrsMatrix::OptimizeStorage() returned an error");
  /*----------------------------------------------------------------------*/
  return;
} /* end of close_trilinos_matrix */


/*----------------------------------------------------------------------*/
/*!
  \brief call FillComplete with explicit arguments

  Complete a non-square matrix. This requires an explicit
  DomainMap. This is assumed to be the RowMap of another matrix. The
  RowMap of the matrix to be closed is assumed to be the right
  one. Thus one additional matrix has to be given.

  \param A    (i/o) matrix to be completed
  \param cmat   (i) matrix that brings additional rowmap

  \author u.kue
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void close_nonquad_trilinos_matrix(TRILINOSMATRIX *A, TRILINOSMATRIX *cmat)
{
  DSTraceHelper dsh("close_nonquad_trilinos_matrix");

  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;
  if (!Amat) dserror("trilinos matrix not set");
  if (Amat->Filled()) return;

  Epetra_CrsMatrix* colmat = (Epetra_CrsMatrix*)cmat->matrix;
  if (!colmat) dserror("trilinos matrix not set");

  int err = Amat->FillComplete(colmat->RowMap(),Amat->RowMap());
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domainmap,rowmap) returned an error");

  err = Amat->OptimizeStorage();
  if (err) dserror("Epetra_CrsMatrix::OptimizeStorage() returned an error");
}

/*----------------------------------------------------------------------*
  |                                                           m.gee 9/06|
  | add to = to + from*factor                                           |
  | do not FillComplete to upon exit if not Filled before               |
  | do call FillComplete to upon exit if Filled before                  |
 *----------------------------------------------------------------------*/
void add_trilinos_matrix(TRILINOSMATRIX* from, TRILINOSMATRIX* to, double factor)

{
  DSTraceHelper dst("add_trilinos_matrix");
  /*----------------------------------------------------------------------*/
  if (!from->matrix || !to->matrix) dserror("Either from or to matrix is NULL");


  Epetra_CrsMatrix* mfrom = (Epetra_CrsMatrix*)from->matrix;
  Epetra_CrsMatrix* mto   = (Epetra_CrsMatrix*)to->matrix;

  // Matrix from has to be filled
  if (!mfrom->Filled()) dserror("FillComplete() was not called on matrix mfrom");

  // Matrix to must NOT be filled
  if (mto->Filled())
  {
    // target has been called FillComplete(), we can't add to it anymore
    // Create a new one and add both old ones
    Epetra_CrsMatrix* newmatrix = new Epetra_CrsMatrix(Copy,mto->RowMap(),
                                                       to->NumEntriesPerRow,
                                                       to->StaticProfile);
    MOERTEL::MatrixMatrixAdd(*mto,false,1.0,*newmatrix,0.0);
    MOERTEL::MatrixMatrixAdd(*mfrom,false,factor,*newmatrix,1.0);
    newmatrix->FillComplete(mto->OperatorDomainMap(),mto->OperatorRangeMap());
    newmatrix->OptimizeStorage();
    delete mto;
    to->matrix = (void*)newmatrix;
  }
  else
    MOERTEL::MatrixMatrixAdd(*mfrom,false,factor,*mto,1.0);


  /*----------------------------------------------------------------------*/
  return;
} /* end of add_trilinos_matrix */


/*----------------------------------------------------------------------*
 |                                                            m.gee 9/06|
 | computes y = A*x
 *----------------------------------------------------------------------*/
void matvec_trilinos(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A)
{
  DSTraceHelper dst("matvec_trilinos");
  /*----------------------------------------------------------------------*/
  // get Epetra_CrsMatrix
  if (!A->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;

  // test Amat
  if (!Amat->Filled()) dserror("FillComplete() was not called on Amat");

  // wrap y and x in Epetra_Vector classes
  Epetra_Vector ex(View,Amat->OperatorDomainMap(),x->vec.a.dv);
  Epetra_Vector ey(View,Amat->OperatorRangeMap(),y->vec.a.dv);

  // do multiply
  int err = Amat->Multiply(false,ex,ey);
  if (err) dserror("Epetra_CrsMatrix::Multiply returned an error");
  /*----------------------------------------------------------------------*/
  return;
} /* end of matvec_trilinos */


/*----------------------------------------------------------------------*
 |                                                            m.gee 9/06|
 | computes y = A^T*x
 *----------------------------------------------------------------------*/
void matvec_trilinos_trans(DIST_VECTOR* y, DIST_VECTOR* x, TRILINOSMATRIX* A)
{
  DSTraceHelper dst("matvec_trilinos_trans");
  /*----------------------------------------------------------------------*/
  // get Epetra_CrsMatrix
  if (!A->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;

  // test Amat
  if (!Amat->Filled()) dserror("FillComplete() was not called on Amat");

  // wrap y and x in Epetra_Vector classes
  Epetra_Vector ex(View,Amat->OperatorRangeMap(),x->vec.a.dv);
  Epetra_Vector ey(View,Amat->OperatorDomainMap(),y->vec.a.dv);

  // do multiply
  int err = Amat->Multiply(true,ex,ey);
  if (err) dserror("Epetra_CrsMatrix::Multiply returned an error");
  /*----------------------------------------------------------------------*/
  return;
} /* end of matvec_trilinos_trans */


/*----------------------------------------------------------------------*
 |                                                            m.gee 9/06|
 | computes A = A*factor
 *----------------------------------------------------------------------*/
void scale_trilinos_matrix(TRILINOSMATRIX* A, double factor)
{
  DSTraceHelper dst("scale_trilinos_matrix");
  /*----------------------------------------------------------------------*/
  // get Epetra_CrsMatrix
  if (!A->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;

  Amat->Scale(factor);
  /*----------------------------------------------------------------------*/
  return;
} /* end of scale_trilinos_matrix */


/*----------------------------------------------------------------------*/
/*!
  \brief invert a diagonal matrix

  \param tri       (i/o)  the diagonal matrix

  \warning It is not tested if the matrix is really diagonal.

  \author kuettler
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void invert_trilinos_diagonal_matrix(TRILINOSMATRIX* A)
{
  DSTraceHelper dst("invert_trilinos_diagonal_matrix");

  if (!A->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;

  Epetra_Vector Diagonal(Amat->RowMap(),false);
  if (Amat->ExtractDiagonalCopy(Diagonal)!=0)
    dserror("error in ExtractDiagonalCopy(Diagonal)");
  if (Diagonal.Reciprocal(Diagonal)!=0)
    dserror("error in Reciprocal(Diagonal)");
  if (Amat->ReplaceDiagonalValues(Diagonal)!=0)
    dserror("error in ReplaceDiagonalValues(Diagonal)");
}


/*----------------------------------------------------------------------*/
/*!
  \brief matrix-matrix-matrix multiplication

  Needed for a particular implementation of the projection method. The
  middle matrix is diagonal, the first is the transposed of the
  last within that method. But this routine works with any set of
  matrices with matching dimensions.

  \param dest   (o) destination matrix
  \param A      (i) input matrix
  \param transA (i) transposed flag of A
  \param B      (i) input matrix
  \param transB (i) transposed flag of B
  \param C      (i) input matrix
  \param transC (i) transposed flag of C

  \warning the ccarat wrapper dest is not filled completely, comm and
  rowmap are missing. You probably do not need them.

  \author u.kue
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void mult_trilinos_mmm(TRILINOSMATRIX* dest,
                       TRILINOSMATRIX* A,
                       INT transA,
                       TRILINOSMATRIX* B,
                       INT transB,
                       TRILINOSMATRIX* C,
                       INT transC)
{
  DSTraceHelper dst("mult_trilinos_mmm");

  if (!A->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;
  if (!B->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Bmat = (Epetra_CrsMatrix*)B->matrix;
  if (!C->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Cmat = (Epetra_CrsMatrix*)C->matrix;

  if (dest->matrix) dserror("destination matrix already set");

  RefCountPtr<Epetra_CrsMatrix> AB = rcp(MOERTEL::MatMatMult(*Amat,transA,*Bmat,transB,0));

  Epetra_CrsMatrix* ABC = MOERTEL::MatMatMult(*AB,false,*Cmat,transC,0);
  dest->matrix = (void*)ABC;

  // if these are ever used, we can set them
  dest->epetracomm = NULL;
  dest->rowmap = NULL;

  dest->is_init=1;
}


/*----------------------------------------------------------------------*/
/*!
  \brief matrix-matrix-matrix multiplication

  Needed for a particular implementation of the projection method. The
  middle matrix is diagonal, the first is the transposed of the
  last within that method. But this routine works with any set of
  matrices with matching dimensions.

  \param dest   (o) destination matrix
  \param A      (i) input matrix
  \param transA (i) transposed flag of A
  \param B      (i) input matrix
  \param transB (i) transposed flag of B
  \param C      (i) input matrix
  \param transC (i) transposed flag of C

  \warning the ccarat wrapper dest is not filled completely, comm and
  rowmap are missing. You probably do not need them.

  \author u.kue
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void mult_trilinos_mmm_cont(TRILINOSMATRIX* dest,
                            TRILINOSMATRIX* A,
                            INT transA,
                            TRILINOSMATRIX* B,
                            INT transB,
                            TRILINOSMATRIX* C,
                            INT transC)
{
  DSTraceHelper dst("mult_trilinos_mmm_cont");

  if (!A->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Amat = (Epetra_CrsMatrix*)A->matrix;
  if (!B->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Bmat = (Epetra_CrsMatrix*)B->matrix;
  if (!C->matrix) dserror("Matrix is NULL");
  Epetra_CrsMatrix* Cmat = (Epetra_CrsMatrix*)C->matrix;

  if (!dest->matrix) dserror("destination matrix must be set");

  RefCountPtr<Epetra_CrsMatrix> AB = rcp(MOERTEL::MatMatMult(*Amat,transA,*Bmat,transB,0));

  Epetra_CrsMatrix* ABC = (Epetra_CrsMatrix*)dest->matrix;

  int err = EpetraExt::MatrixMatrix::Multiply(*AB,false,*Cmat,transC,*ABC);
  if (err) dserror("error %d in MatrixMatrix()",err);

  // We expect to obtain the same RowMap, but we do not check that.

  dest->is_init=1;
}


#ifdef DEBUG

void extractGlobalRow(FILE* out,TRILINOSMATRIX* mat,INT row)
{
  DSTraceHelper dst("extractGlobalRow");

  int length=1000;

  vector<double> values(length);
  vector<int> indices(length);
  int NumEntries=0;

  Epetra_CrsMatrix* A = (Epetra_CrsMatrix*)mat->matrix;
  int err = A->ExtractGlobalRowCopy(row,length,NumEntries,&values[0],&indices[0]);
  if (err)
    dserror("ExtractGlobalRowCopy failed err=%d",err);

  fprintf(out,"NumGlobalRows=%d\n",A->NumGlobalRows());
  fprintf(out,"NumGlobalCols=%d\n",A->NumGlobalCols());
  fprintf(out,"indices = ");
  for (int i=0; i<NumEntries; ++i)
  {
    if (i>0)
      fprintf(out,",");
    fprintf(out,"%d",indices[i]);
  }
  fprintf(out,"\n");

  fprintf(out,"values = ");
  for (int i=0; i<NumEntries; ++i)
  {
    if (i>0)
      fprintf(out,",");
    fprintf(out,"%f",values[i]);
  }
  fprintf(out,"\n");
}

#endif

#endif // TRILINOS_PACKAGE
