#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  create number of distributed vectors - collective call ! m.gee 10/01|
 *----------------------------------------------------------------------*/
void solserv_create_vec(
                           DIST_VECTOR         **vector,
                           int                   numvectors,
                           int                   numeq_total,
                           int                   numeq,
                           char                  typstr[])
{
int                  i;
DIST_VECTOR *actvector;
#ifdef DEBUG 
dstrc_enter("solserv_create_vec");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ allocate the vectors */
*vector = (DIST_VECTOR*)calloc(numvectors,sizeof(DIST_VECTOR));
if (*vector==NULL) dserror("Allocation of DIST_VECTOR failed");
/*--------------------------- loop the created vectors and perfrom init */
for (i=0; i<numvectors; i++)
{
   actvector = &((*vector)[i]);
   actvector->numeq_total = numeq_total;
   actvector->numeq       = numeq;
   amdef("dist_vec",&(actvector->vec),numeq,1,typstr);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_create_vec */




/*----------------------------------------------------------------------*
 |  init a distributed vector to zero - collective call !    m.gee 10/01|
 *----------------------------------------------------------------------*/
void solserv_zero_vec(DIST_VECTOR *disvector)
{
int                  i;
#ifdef DEBUG 
dstrc_enter("solserv_zero_vec");
#endif
/*----------------------------------------------------------------------*/
amzero(&(disvector->vec));
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_zero_vec */



/*----------------------------------------------------------------------*
 |  init a distributed matrix to zero - collective call !    m.gee 10/01|
 *----------------------------------------------------------------------*/
void solserv_zero_mat(INTRA *actintra, SPARSE_ARRAY *mat,SPARSE_TYP *mattyp)
{
int                  imyrank;
int                  inprocs;
int                  ilower,iupper,jlower,jupper;
int                  err=1;
#ifdef DEBUG 
dstrc_enter("solserv_zero_mat");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
switch (*mattyp)
{
case mds:
break;
case msr:
   amzero(&(mat->msr->val));
   mat->msr->is_factored=0;
break;
case parcsr:/*---- this stupid package does not have a zero function!!!!*/
#ifdef HYPRE_PACKAGE
   err=HYPRE_IJMatrixDestroy(mat->parcsr->ij_matrix);
   if (err) dserror("Cannot destroy PARCSR matrix");
   ilower = jlower = mat->parcsr->perm.a.ia[imyrank][0];
   iupper = jupper = mat->parcsr->perm.a.ia[imyrank][mat->parcsr->perm_sizes.a.iv[imyrank]-1];
   err=HYPRE_IJMatrixCreate(actintra->MPI_INTRA_COMM,ilower,iupper,jlower,jupper,&(mat->parcsr->ij_matrix));
   if (err) dserror("Cannot Create PARCSR matrix");
   err=HYPRE_IJMatrixSetObjectType(mat->parcsr->ij_matrix,HYPRE_PARCSR);
   if (err) dserror("Cannot Set PARCSR matrix object type");
   err=HYPRE_IJMatrixInitialize(mat->parcsr->ij_matrix);
   if (err) dserror("Cannot init PARCSR matrix");
#endif
   mat->parcsr->is_factored=0;
break;
case ucchb:
   amzero(&(mat->ucchb->a));
   mat->ucchb->is_factored=0;
break;
case dense:
   amzero(&(mat->dense->A));
   mat->dense->is_factored=0;
break;
case rc_ptr:
   amzero(&(mat->rc_ptr->A_loc));
   mat->rc_ptr->is_factored=0;
break;
case skymatrix:
   amzero(&(mat->sky->A));
   mat->sky->is_factored=0;
break;
case sparse_none:
   dserror("Unknown typ of sparse distributed system matrix");
break;
default:
   dserror("Unknown typ of sparse distributed system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_zero_mat */



/*----------------------------------------------------------------------*
 |  add contents of the vector vec_from to vec_to            m.gee 10/01|
 *----------------------------------------------------------------------*/
void solserv_add_vec(DIST_VECTOR *vec_from, 
                        DIST_VECTOR *vec_to)
{
int                  i,dim;
double              *dfrom;
double              *dto;
#ifdef DEBUG 
dstrc_enter("solserv_add_vec");
#endif
/*----------------------------------------------------------------------*/
if (vec_from->vec.fdim != vec_to->vec.fdim)
dserror("Cannot copy distributed vectors, not same dimension");
/*----------------------------------------------------------------------*/
dim   = vec_from->vec.fdim;
dfrom = vec_from->vec.a.dv;
dto   = vec_to->vec.a.dv;
for (i=0; i<dim; i++) *(dto++) += *(dfrom++);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_add_vec */


/*----------------------------------------------------------------------*
 |  copy contents of the vector vec_from to vec_to           m.gee 11/01|
 *----------------------------------------------------------------------*/
void solserv_copy_vec(DIST_VECTOR *vec_from, 
                        DIST_VECTOR *vec_to)
{
int                  i,dim;
double              *dfrom;
double              *dto;
#ifdef DEBUG 
dstrc_enter("solserv_copy_vec");
#endif
/*----------------------------------------------------------------------*/
if (vec_from->vec.fdim != vec_to->vec.fdim)
dserror("Cannot copy distributed vectors, not same dimension");
/*----------------------------------------------------------------------*/
dim   = vec_from->vec.fdim;
dfrom = vec_from->vec.a.dv;
dto   = vec_to->vec.a.dv;
for (i=0; i<dim; i++) *(dto++) = *(dfrom++);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_copy_vec */




/*----------------------------------------------------------------------*
 |  make euclidian norm of a distributed vector              m.gee 11/01|
 *----------------------------------------------------------------------*/
void solserv_vecnorm_euclid(INTRA       *actintra,
                            DIST_VECTOR *dist_vec, 
                            double      *result)
{
int                  i;
double              *vec;
double               sendbuff;
int                  numeq;
#ifdef DEBUG 
dstrc_enter("solserv_vecnorm_euclid");
#endif
/*----------------------------------------------------------------------*/
vec      = dist_vec->vec.a.dv;
numeq    = dist_vec->numeq;
sendbuff = 0.0;
for (i=0; i<dist_vec->numeq; i++)
{
   sendbuff += vec[i]*vec[i];
}
#ifdef PARALLEL 
MPI_Allreduce(&sendbuff,result,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
*result = sendbuff;
#endif
*result = sqrt(*result);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_vecnorm_euclid */

/*----------------------------------------------------------------------*
 |  get a certain entry from a distr. vector to all procs    m.gee 11/01|
 *----------------------------------------------------------------------*/
void solserv_getele_vec(INTRA       *actintra,
                        SPARSE_TYP   *sysarray_typ,
                        SPARSE_ARRAY *sysarray,
                        DIST_VECTOR *dist_vec,
                        int          indiz, 
                        double      *result)
{
int                  i;
int                  imyrank;
int                  inprocs;
int                  bcaster;
int                  recvbuff;
int                  index;
int                 *update;
double              *vec;
int                  numeq;
#ifdef DEBUG 
dstrc_enter("solserv_getele_vec");
#endif
/*----------------------------------------------------------------------*/
imyrank  = actintra->intra_rank;
inprocs  = actintra->intra_nprocs;
vec      = dist_vec->vec.a.dv;
numeq    = dist_vec->numeq;
switch(*sysarray_typ)
{
case msr:
   update = sysarray->msr->update.a.iv;
break;
case parcsr:
   update = sysarray->parcsr->update.a.ia[imyrank];
break;
case ucchb:
   update = sysarray->ucchb->update.a.iv;
break;
case dense:
   update = sysarray->dense->update.a.iv;
break;
case rc_ptr:
   update = sysarray->rc_ptr->update.a.iv;
break;
case skymatrix:
   update = sysarray->sky->update.a.iv;
break;
default:
   dserror("Unknown typ of system matrix given");
break;
}
index = find_index(indiz,update,numeq);
#ifndef PARALLEL /* this is sequentiell */
if (index==-1) dserror("Cannot find indize in distributed vector");
*result = dist_vec->vec.a.dv[index]; 
#else            /* this is parallel */
bcaster=-1;
if (index != -1) 
{
   bcaster = imyrank;
   *result = dist_vec->vec.a.dv[index]; 
}
MPI_Allreduce(&bcaster,&recvbuff,1,MPI_INT,MPI_MAX,actintra->MPI_INTRA_COMM);
bcaster = recvbuff;
if (bcaster==-1) dserror("Cannot find indize in distributed vector");

MPI_Bcast(result,1,MPI_DOUBLE,bcaster,actintra->MPI_INTRA_COMM);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_getele_vec */

/*----------------------------------------------------------------------*
 |  make dot product between 2 distr. vectors                m.gee 11/01|
 *----------------------------------------------------------------------*/
void solserv_dot_vec(INTRA       *actintra,
                     DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,
                     double      *dot)
{
int                  i;
double               localsum;
double               globalsum;
double              *vec1;
double              *vec2;

#ifdef DEBUG 
dstrc_enter("solserv_dot_vec");
#endif
/*----------------------------------------------------------------------*/
if (dist_vec1->numeq != dist_vec2->numeq)
   dserror("Mismatch in dimensions, cannot do dot-product");
/*----------------------------------------------------------------------*/
vec1 = dist_vec1->vec.a.dv;
vec2 = dist_vec2->vec.a.dv;
/*----------------------------------------------------------------------*/
localsum = 0.0;   

for (i=0; i<dist_vec1->numeq; i++) localsum += vec1[i]*vec2[i];

#ifdef PARALLEL 
MPI_Allreduce(&localsum,&globalsum,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
*dot = globalsum;
#else
*dot = localsum;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_dot_vec */


/*----------------------------------------------------------------------*
 |  make product between scalar and distr. vector            m.gee 11/01|
 *----------------------------------------------------------------------*/
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,
                            double       scalar)
{
int                  i;
int                  dim;
double              *dptr;

#ifdef DEBUG 
dstrc_enter("solserv_dot_vec");
#endif
/*----------------------------------------------------------------------*/
dptr = dist_vec->vec.a.dv;
dim  = dist_vec->numeq;
for (i=0; i<dim; i++) *(dptr++) *= scalar;
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_dot_vec */



/*----------------------------------------------------------------------*
 |  Allreduce a distributed vector in an INTRACOMM           m.gee 10/01|
 |  This is a collective call!                                          |
 *----------------------------------------------------------------------*/
void solserv_reddistvec(DIST_VECTOR  *distvec,
                          SPARSE_ARRAY *sysarray,
                          SPARSE_TYP   *sysarray_typ,
                          double       *fullvec,
                          int           dim,
                          INTRA        *actintra)
{
int             i;
int             dof,dofperm;
double         *dfrom;
int             imyrank;
int             inprocs;
#ifdef PARALLEL 
static double  *recvbuff;
static ARRAY    recv;    
#endif

#ifdef DEBUG 
dstrc_enter("solserv_reddistvec");
#endif
/*----------------------------------------------------------------------*/
if (dim != distvec->numeq_total) dserror("Dimension mismatch");
/*------------------------- allocate communication buffer, if necessary */
#ifdef PARALLEL 
if (recv.Typ != DV) 
{
   recvbuff = amdef("recvbuff",&recv,dim,1,"DV");
}
if (dim > recv.fdim)
{
   amdel(&recv);
   recvbuff = amdef("recvbuff",&recv,dim,1,"DV");
}
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
switch(*sysarray_typ)
{
case msr:
   for (i=0; i<sysarray->msr->numeq; i++)
   {
      dof = sysarray->msr->update.a.iv[i];
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
#endif
break;
case parcsr:
/*
   for (i=0; i<sysarray->parcsr->numeq; i++)
   {
      dof     = sysarray->parcsr->update.a.ia[imyrank][i];
      dofperm = sysarray->parcsr->perm.a.ia[imyrank][i];
      fullvec[dof] = distvec->vec.a.dv[dofperm];
   }
*/
   for (i=0; i<sysarray->parcsr->numeq; i++)
   {
      dof     = sysarray->parcsr->update.a.ia[imyrank][i];
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
#endif
break;
case ucchb:
   for (i=0; i<sysarray->ucchb->numeq; i++)
   {
      dof = sysarray->ucchb->update.a.iv[i];
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
#endif
break;
case dense:
   for (i=0; i<sysarray->dense->numeq; i++)
   {
      dof = sysarray->dense->update.a.iv[i];
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
#endif
break;
case mds:
   for (i=0; i<sysarray->mds->numeq; i++)
   {
      fullvec[i] = distvec->vec.a.dv[i];
   }
break;
case rc_ptr:
   for (i=0; i<sysarray->rc_ptr->numeq; i++)
   {
      dof = sysarray->rc_ptr->update.a.iv[i];
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
#endif
break;
case skymatrix:
   for (i=0; i<sysarray->sky->numeq; i++)
   {
      dof = sysarray->sky->update.a.iv[i];
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
#endif
break;
default:
   dserror("Unknown typ of system matrix given");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_reddistvec */




/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 10/01|
 |  certain place  in ARRAY sol                                         |
 |  Result has to be allreduced and are put to the whole                |
 |  field on each proc                                                  |
 *----------------------------------------------------------------------*/
void solserv_result_total(
                          FIELD          *actfield,
                          INTRA          *actintra,
                          DIST_VECTOR    *sol,
                          int             place,
                          SPARSE_ARRAY   *sysarray,
                          SPARSE_TYP     *sysarray_typ
                         )
{
int      i,j;
int      max;
int      diff;
int      dof;

int      numeq_total;
NODE    *actnode;
ARRAY    result_a;
double  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_total");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->numnp; i++)
{
   actnode = &(actfield->node[i]);
   /*---------------------------------------- enlarge sol, if necessary */
   if (place >= actnode->sol.fdim)
   {
      diff = place - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max,actnode->sol.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
      actnode->sol.a.da[place][j] = result[dof];
   }   
   
}
/*----------------------------------------------------------------------*/
amdel(&result_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_result_total */


/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_increment                                |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 *----------------------------------------------------------------------*/
void solserv_result_incre(
                          FIELD          *actfield,
                          INTRA          *actintra,
                          DIST_VECTOR    *sol,
                          int             place,
                          SPARSE_ARRAY   *sysarray,
                          SPARSE_TYP     *sysarray_typ
                         )
{
int      i,j;
int      max;
int      diff;
int      dof;

int      numeq_total;
NODE    *actnode;
ARRAY    result_a;
double  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_incre");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->numnp; i++)
{
   actnode = &(actfield->node[i]);
   /*------------------------------ enlarge sol_increment, if necessary */
   if (place >= actnode->sol.fdim)
   {
      diff = place - actnode->sol_increment.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
      actnode->sol_increment.a.da[place][j] = result[dof];
   }   
   
}
/*----------------------------------------------------------------------*/
amdel(&result_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_result_incre */



/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 11/01|
 |  certain place in ARRAY sol_residual                                 |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 *----------------------------------------------------------------------*/
void solserv_result_resid(
                          FIELD          *actfield,
                          INTRA          *actintra,
                          DIST_VECTOR    *sol,
                          int             place,
                          SPARSE_ARRAY   *sysarray,
                          SPARSE_TYP     *sysarray_typ
                         )
{
int      i,j;
int      max;
int      diff;
int      dof;

int      numeq_total;
NODE    *actnode;
ARRAY    result_a;
double  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_resid");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->numnp; i++)
{
   actnode = &(actfield->node[i]);
   /*------------------------------- enlarge sol_residual, if necessary */
   if (place >= actnode->sol.fdim)
   {
      diff = place - actnode->sol_residual.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol_residual),actnode->sol_residual.fdim+max,actnode->sol_residual.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
      actnode->sol_residual.a.da[place][j] = result[dof];
   }   
   
}
/*----------------------------------------------------------------------*/
amdel(&result_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_result_resid */
