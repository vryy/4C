#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  create number of distributed vectors - collective call ! m.gee 10/01|
 |  DIST_VECTOR **vector (i/o) adress of pointer a vector of            |
 |                             DIST_VECTORs will be allocated to        |
 |  INT numvectors       (i)   number of DIST_VECTORs to allocate       |
 |  INT numeq_total      (i)   proc-global dimension of the DIST_VECTORs|
 |  INT numeq            (i)   proc_local  dimension of the DIST_VECTORs|
 |  char typstr[]        (i)   ="DV" for DOUBLE-DIST_VECTORs            |
 |  the values in the DIST_VECTORs is NOT initialized                   |
 *----------------------------------------------------------------------*/
void solserv_create_vec(DIST_VECTOR **vector,INT numvectors,INT numeq_total,
                        INT numeq,char typstr[])
{
INT                  i;
DIST_VECTOR *actvector;
#ifdef DEBUG 
dstrc_enter("solserv_create_vec");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ allocate the vectors */
*vector = (DIST_VECTOR*)CCACALLOC(numvectors,sizeof(DIST_VECTOR));
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
 |   delete number of distributed vectors - collective call ! m.gee 2/02|
 |  DIST_VECTOR **vector (i/o) adress of pointer a vector of            |
 |                             DIST_VECTORs is allocated to             |
 |  INT numvectors       (i)   number of DIST_VECTORs to free           |
 |  the routine frees all DIST_VECTORs in vector and sets vector=NULL   |
 *----------------------------------------------------------------------*/
void solserv_del_vec(DIST_VECTOR **vector,INT numvectors)
{
INT                  i;
DIST_VECTOR *actvector;
#ifdef DEBUG 
dstrc_enter("solserv_del_vec");
#endif
/*--------------------------- loop the created vectors and delete them */
for (i=0; i<numvectors; i++)
{
   actvector = &((*vector)[i]);
   amdel(&(actvector->vec));
}
CCAFREE(*vector);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_del_vec */




/*----------------------------------------------------------------------*
 |  init a distributed vector to zero - collective call !    m.gee 10/01|
 |  DIST_VECTOR *disvector (i/o) adress of a DIST_VECTOR to be set to 0.0|
 *----------------------------------------------------------------------*/
void solserv_zero_vec(DIST_VECTOR *disvector)
{
INT                  i;
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
 |  add contents of the vector vec_from to vec_to            m.gee 10/01|
 |  vec_to->vec.a.dv[i] += vec_from->vec.a.dv[i]*factor                 |
 |  DIST_VECTOR *vec_from (i)   vector to be added to another vector    |
 |  DIST_VECTOR *vec_to   (i/o) vector to be added to                   |
 |  DOUBLE factor         (i)   scaling factor                          |
 *----------------------------------------------------------------------*/
void solserv_add_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to,DOUBLE factor)
{
INT                  i,dim;
DOUBLE              *dfrom;
DOUBLE              *dto;
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
for (i=0; i<dim; i++) *(dto++) += (*(dfrom++) * factor);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_add_vec */




/*----------------------------------------------------------------------*
 |  copy contents of the vector vec_from to vec_to           m.gee 11/01|
 |  vec_to->vec.a.dv[i] = vec_from->vec.a.dv[i]                         |
 |  DIST_VECTOR *vec_from (i)   vector to be copied to another vector   |
 |  DIST_VECTOR *vec_to   (i/o) vector to be copied to                  |
 |  user must assure matching dimensions and types                      |
 *----------------------------------------------------------------------*/
void solserv_copy_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to)
{
INT                  i,dim;
DOUBLE              *dfrom;
DOUBLE              *dto;
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
 |  *result = sqrt( sumof(vec[i]*vec[i]) )                              |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  DOUBLE *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_euclid(INTRA *actintra,DIST_VECTOR *dist_vec,DOUBLE *result)
{
INT                  i;
DOUBLE              *vec;
DOUBLE               sendbuff;
INT                  numeq;
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
 |  find  absolute maximum value in a vector (Linf-Norm)     m.gee 02/02|
 |  *result = MAX( ABS(vec[i]) )                                        |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  DOUBLE *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_Linf(INTRA *actintra,DIST_VECTOR *dist_vec,DOUBLE *result)
{
INT                  i;
DOUBLE              *vec;
DOUBLE               sendbuff=0.0;
INT                  numeq;
#ifdef DEBUG 
dstrc_enter("solserv_vecnorm_Linf");
#endif
/*----------------------------------------------------------------------*/
vec      = dist_vec->vec.a.dv;
numeq    = dist_vec->numeq;
sendbuff = 0.0;
for (i=0; i<dist_vec->numeq; i++)
{
   if (FABS(vec[i])>sendbuff) sendbuff = FABS(vec[i]);
}
#ifdef PARALLEL 
MPI_Allreduce(&sendbuff,result,1,MPI_DOUBLE,MPI_MAX,actintra->MPI_INTRA_COMM);
#else
*result = sendbuff;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_vecnorm_Linf */




/*----------------------------------------------------------------------*
 |  get a certain entry from a distr. vector to all procs    m.gee 11/01|
 |  returns the value of dof indiz in the vector dist_vec on all procs  |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  SPARSE_TYP *sysarray_typ (i) sparsity typ of vector-matching matrix |
 |  SPARSE_ARRAY *sysarray   (i) sparse matrix the vector matches in    |
 |                               distribution                           | 
 |  DIST_VECTOR  *dist_vec   (i) vector the value shall be taken from   |
 |  INT           indiz      (i) field-local (unsupported) dof number   |
 |  DOUBLE       *result     (o) value in vector at the given dof       |
 |                               returned redundant on all procs        |
 *----------------------------------------------------------------------*/
void solserv_getele_vec(INTRA*actintra,SPARSE_TYP *sysarray_typ,
                        SPARSE_ARRAY *sysarray,DIST_VECTOR *dist_vec,
                        INT indiz,DOUBLE *result)
{
INT                  i;
INT                  imyrank;
INT                  inprocs;
INT                  bcaster;
INT                  recvbuff;
INT                  index;
INT                 *update;
DOUBLE              *vec;
INT                  numeq;
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
case spoolmatrix:
   update = sysarray->spo->update.a.iv;
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
case ccf:
   update = sysarray->ccf->update.a.iv;
break;
case skymatrix:
   update = sysarray->sky->update.a.iv;
break;
case bdcsr:
   update = sysarray->bdcsr->update.a.iv;
break;
case mds:
   index = indiz;
break;
case oll:
   update = sysarray->oll->update.a.iv;
break;
default:
   dserror("Unknown typ of system matrix given");
break;
}
if(*sysarray_typ!=mds)
{
   index = find_index(indiz,update,numeq);
}
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
 |  *dot = sumover_i( vec1[i]*vec2[i] )                                 |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTORs live on     |
 |  DIST_VECTOR *dist_vec1 (i) first vector to be multiplied            |
 |  DIST_VECTOR *dist_vec2 (i) scnd  vector to be multiplied            |
 |  DOUBLE      *dot       (o) result of vector-vector multiplication   |
 |                             returned redundant on all procs          |
 *----------------------------------------------------------------------*/
void solserv_dot_vec(INTRA *actintra,DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,DOUBLE *dot)
{
INT                  i;
DOUBLE               localsum;
DOUBLE               globalsum;
DOUBLE              *vec1;
DOUBLE              *vec2;

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
 |  vec[i] = vec[i] * scalar                                            |
 |  DIST_VECTOR *dist_vec (i/o) vector to be multiplied by scalar       |
 |  DOUBLE       scalar   (o)   scalar value                            |
 *----------------------------------------------------------------------*/
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,DOUBLE scalar)
{
INT                  i;
INT                  dim;
DOUBLE              *dptr;

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
 |  distributed vector to full redundant vector                         |
 |                                                                      |
 |  note that the disributed vectors match a certain type of sparse     |
 |  matrix in the layout of distribution and values. This means, that   |
 |  the value of a certain dof are NOT in distvec->vec.a.dv[dof]!!!!!   |
 |                                                                      |
 |  the redundant vector fullvec holds values of a certain dof in       |
 |  fullvec[dof]                                                        |
 |                                                                      |
 |  the values in the given DIST_VECTOR are copied to a vector of       |
 |  size numeq_total, which is redundant on all procs                   |
 |  DIST_VECTOR *distvec (i) DIST_VECTORto be 'allreduced'              |
 |  SPARSE_ARRAY *sysarray (i) sparse matrix matching the distribution  |
 |                             of distvec                               |
 |  SPARSE_TYP *sysarray_typ (i) type of sparse matrix                  |
 |  DOUBLE *fullvec (o) vector of lenght numeq_total will be holding    |
 |                      the values from distvec in correct dof-ordering:|
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_reddistvec(DIST_VECTOR *distvec,SPARSE_ARRAY *sysarray,
                        SPARSE_TYP *sysarray_typ,DOUBLE *fullvec,
                        INT dim,INTRA *actintra)
{
INT             i;
INT             dof,dofperm;
DOUBLE         *dfrom;
INT             imyrank;
INT             inprocs;
#ifdef PARALLEL 
static DOUBLE  *recvbuff;
static ARRAY    recv;    
#endif

#ifdef DEBUG 
dstrc_enter("solserv_reddistvec");
#endif
/*----------------------------------------------------------------------*/
if (dim != distvec->numeq_total) dserror("Dimension mismatch");
/*------------------------- allocate communication buffer, if necessary */
#ifdef PARALLEL 
if (recv.Typ != cca_DV) 
{
   recvbuff = amdef("recvbuff",&recv,dim,1,"DV");
}
if (dim > recv.fdim)
{
   amdel(&recv);
   recvbuff = amdef("recvbuff",&recv,dim,1,"DV");
}
amzero(&recv);
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<dim; i++) fullvec[i] = 0.0;
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
#ifdef PARALLEL 
      recvbuff[dof] = distvec->vec.a.dv[i];
#else
      fullvec[dof] = distvec->vec.a.dv[i];
#endif
   }
#ifdef PARALLEL 
   MPI_Allreduce(recvbuff,fullvec,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
break;



case parcsr:
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

case spoolmatrix:
   for (i=0; i<sysarray->spo->numeq; i++)
   {
      dof = sysarray->spo->update.a.iv[i];
#ifdef PARALLEL 
      recvbuff[dof] = distvec->vec.a.dv[i];
#else
      fullvec[dof] = distvec->vec.a.dv[i];
#endif
   }
#ifdef PARALLEL 
   MPI_Allreduce(recvbuff,fullvec,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
break;


case ccf:
   for (i=0; i<sysarray->ccf->numeq; i++)
   {
      dof = sysarray->ccf->update.a.iv[i];
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


case bdcsr:
   for (i=0; i<sysarray->bdcsr->numeq; i++)
   {
      dof = sysarray->bdcsr->update.a.iv[i];
#ifdef PARALLEL 
      recvbuff[dof] = distvec->vec.a.dv[i];
#else
      fullvec[dof] = distvec->vec.a.dv[i];
#endif
   }
#ifdef PARALLEL 
   MPI_Allreduce(recvbuff,fullvec,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
break;


case oll:
   for (i=0; i<sysarray->oll->numeq; i++)
   {
      dof = sysarray->oll->update.a.iv[i];
#ifdef PARALLEL 
      recvbuff[dof] = distvec->vec.a.dv[i];
#else
      fullvec[dof] = distvec->vec.a.dv[i];
#endif
   }
#ifdef PARALLEL 
   MPI_Allreduce(recvbuff,fullvec,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
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
 |  distribute a full redundant vector                       m.gee 02/02|
 |  This is a collective call!                                          |
 |  full redundant vector to distributed vector                         |
 |  this routine is the inverse of solserv_reddistvec                   |
 |  It copies the values in a vector fullvec of size numeq_total, that  |
 |  is ordered such that fullvec[dof] = value at dof                    |
 |  to a distributed vector matching a certain sparse matrix in         |
 |  distribution of values. Note that in the distvec the values         |
 |  value_at_dof are NOT in distvec->vec.a.dv[dof] !!!!                 |
 |                                                                      |
 |  DIST_VECTOR *distvec (o) DIST_VECTOR to be copied to                |
 |  SPARSE_ARRAY *sysarray (i) sparse matrix matching the distribution  |
 |                             of distvec                               |
 |  SPARSE_TYP *sysarray_typ (i) type of sparse matrix                  |
 |  DOUBLE *fullvec (o) vector of lenght numeq_total  holding           |
 |                      the values  correct dof-ordering:               |
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_distribdistvec(DIST_VECTOR  *distvec,SPARSE_ARRAY *sysarray,
                            SPARSE_TYP *sysarray_typ,DOUBLE *fullvec,
                            INT dim,INTRA *actintra)
{
INT             i;
INT             dof,dofperm;
DOUBLE         *dfrom;
INT             imyrank;
INT             inprocs;

#ifdef DEBUG 
dstrc_enter("solserv_distribdistvec");
#endif
/*----------------------------------------------------------------------*/
if (dim != distvec->numeq_total) dserror("Dimension mismatch");
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
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case parcsr:
   for (i=0; i<sysarray->parcsr->numeq; i++)
   {
      dof = sysarray->parcsr->update.a.ia[imyrank][i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case ucchb:
   for (i=0; i<sysarray->ucchb->numeq; i++)
   {
      dof = sysarray->ucchb->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case dense:
   for (i=0; i<sysarray->dense->numeq; i++)
   {
      dof = sysarray->dense->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case mds:
   for (i=0; i<sysarray->mds->numeq; i++)
   {
      distvec->vec.a.dv[i] = fullvec[i];
   }
break;

case rc_ptr:
   for (i=0; i<sysarray->rc_ptr->numeq; i++)
   {
      dof = sysarray->rc_ptr->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case ccf:
   for (i=0; i<sysarray->ccf->numeq; i++)
   {
      dof = sysarray->ccf->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case spoolmatrix:
   for (i=0; i<sysarray->spo->numeq; i++)
   {
      dof = sysarray->spo->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case skymatrix:
   for (i=0; i<sysarray->sky->numeq; i++)
   {
      dof = sysarray->sky->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
break;

case oll:
   for (i=0; i<sysarray->oll->numeq; i++)
   {
      dof = sysarray->oll->update.a.iv[i];
      distvec->vec.a.dv[i] = fullvec[dof];
   }
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
} /* end of solserv_distribdistvec */




/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       m.gee 10/01|
 |  certain place  in ARRAY sol                                         |
 |  Result has to be allreduced and are put to the whole                |
 |  field on each procs                                                 |
 |  FIELD *actfield (i) the active field                                |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *sol (i) vector of values to be put to the nodes        |
 |  INT place        (i) place in the ARRAY node->sol where to put the  |
 |                       values. Every structure NODE has an ARRAY sol  |
 |                       of type sol.a.da[place][0..numdf-1]            |
 |                       if place >= actual dimensions of the ARRAY sol |
 |                       sol is enlarged                                |
 |  SPARSE_ARRAY *sysarray (i) sparse matrix matching the distribution  |
 |                             of DIST_VECTOR *sol                      |
 |  SPARSE_TYP *sysarray_typ (i) type of sparse matrix                  |
 |                                                                      |
 *----------------------------------------------------------------------*/
void solserv_result_total(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ)
{
INT      i,j;
INT      max;
INT      diff;
INT      dof;

INT      numeq_total;
NODE    *actnode;
ARRAY    result_a;
DOUBLE  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_total");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(sol,sysarray,sysarray_typ,result,sol->numeq_total,actintra);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*---------------------------------------- enlarge sol, if necessary */
   if (place >= actnode->sol.fdim)
   {
      diff = place - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max+1,actnode->sol.sdim,"DA");
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
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_incre(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ)
{
INT      i,j;
INT      max;
INT      diff;
INT      dof;

INT      numeq_total;
NODE    *actnode;
ARRAY    result_a;
DOUBLE  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_incre");
#endif

/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(sol,sysarray,sysarray_typ,result,sol->numeq_total,actintra);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*------------------------------ enlarge sol_increment, if necessary */
   if (place >= actnode->sol_increment.fdim)
   {
      diff = place - actnode->sol_increment.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max+1,actnode->sol_increment.sdim,"DA");
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
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_resid(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ)
{
INT      i,j;
INT      max;
INT      diff;
INT      dof;

INT      numeq_total;
NODE    *actnode;
ARRAY    result_a;
DOUBLE  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_resid");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(sol,sysarray,sysarray_typ,result,sol->numeq_total,actintra);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*------------------------------- enlarge sol_residual, if necessary */
   if (place >= actnode->sol_residual.fdim)
   {
      diff = place - actnode->sol_residual.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol_residual),actnode->sol_residual.fdim+max+1,actnode->sol_residual.sdim,"DA");
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

/*----------------------------------------------------------------------*
 |  Put the results of a DIST_VECTOR to the nodes in a       genk 01/03 |
 |  certain place in ARRAY sol_mf                                       |
 |  Result have to bee allreduced and are put to the whole              |
 |  field on each proc                                                  |
 |  Functionality is the same as in solserv_result_total                |
 *----------------------------------------------------------------------*/
void solserv_result_mf(FIELD *actfield,INTRA *actintra,DIST_VECTOR *sol,
                          INT place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ)
{
INT      i,j;
INT      max;
INT      diff;
INT      dof;

INT      numeq_total;
NODE    *actnode;
ARRAY    result_a;
DOUBLE  *result;

#ifdef DEBUG 
dstrc_enter("solserv_result_mf");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);
solserv_reddistvec(sol,sysarray,sysarray_typ,result,sol->numeq_total,actintra);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actfield->dis[0].numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);
   /*------------------------------- enlarge sol_residual, if necessary */
   if (place >= actnode->sol_mf.fdim)
   {
      diff = place - actnode->sol_mf.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol_mf),actnode->sol_mf.fdim+max+1,actnode->sol_mf.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      dof = actnode->dof[j];
      if (dof>=numeq_total) continue;
      actnode->sol_mf.a.da[place][j] = result[dof];
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
