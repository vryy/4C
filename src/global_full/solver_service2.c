#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  create number of distributed vectors - collective call ! m.gee 10/01|
 |  DIST_VECTOR **vector (i/o) adress of pointer a vector of            |
 |                             DIST_VECTORs will be allocated to        |
 |  int numvectors       (i)   number of DIST_VECTORs to allocate       |
 |  int numeq_total      (i)   proc-global dimension of the DIST_VECTORs|
 |  int numeq            (i)   proc_local  dimension of the DIST_VECTORs|
 |  char typstr[]        (i)   ="DV" for double-DIST_VECTORs            |
 |  the values in the DIST_VECTORs is NOT initialized                   |
 *----------------------------------------------------------------------*/
void solserv_create_vec(DIST_VECTOR **vector,int numvectors,int numeq_total,
                        int numeq,char typstr[])
{
int                  i;
DIST_VECTOR *actvector;
#ifdef DEBUG 
dstrc_enter("solserv_create_vec");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ allocate the vectors */
*vector = (DIST_VECTOR*)CALLOC(numvectors,sizeof(DIST_VECTOR));
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
 |  int numvectors       (i)   number of DIST_VECTORs to free           |
 |  the routine frees all DIST_VECTORs in vector and sets vector=NULL   |
 *----------------------------------------------------------------------*/
void solserv_del_vec(DIST_VECTOR **vector,int numvectors)
{
int                  i;
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
FREE(*vector);
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
 |  add contents of the vector vec_from to vec_to            m.gee 10/01|
 |  vec_to->vec.a.dv[i] += vec_from->vec.a.dv[i]*factor                 |
 |  DIST_VECTOR *vec_from (i)   vector to be added to another vector    |
 |  DIST_VECTOR *vec_to   (i/o) vector to be added to                   |
 |  double factor         (i)   scaling factor                          |
 *----------------------------------------------------------------------*/
void solserv_add_vec(DIST_VECTOR *vec_from,DIST_VECTOR *vec_to,double factor)
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
 |  *result = sqrt( sumof(vec[i]*vec[i]) )                              |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  double *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_euclid(INTRA *actintra,DIST_VECTOR *dist_vec,double *result)
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
 |  find  absolute maximum value in a vector (Linf-Norm)     m.gee 02/02|
 |  *result = MAX( ABS(vec[i]) )                                        |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 |  DIST_VECTOR *dist_vec (i) vector to make norm of                    |
 |  double *result        (o) norm of the vector                        |
 *----------------------------------------------------------------------*/
void solserv_vecnorm_Linf(INTRA *actintra,DIST_VECTOR *dist_vec,double *result)
{
int                  i;
double              *vec;
double               sendbuff=0.0;
int                  numeq;
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
 |  int           indiz      (i) field-local (unsupported) dof number   |
 |  double       *result     (o) value in vector at the given dof       |
 |                               returned redundant on all procs        |
 *----------------------------------------------------------------------*/
void solserv_getele_vec(INTRA*actintra,SPARSE_TYP *sysarray_typ,
                        SPARSE_ARRAY *sysarray,DIST_VECTOR *dist_vec,
                        int indiz,double *result)
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
case mds:
   index = indiz;
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
 |  double      *dot       (o) result of vector-vector multiplication   |
 |                             returned redundant on all procs          |
 *----------------------------------------------------------------------*/
void solserv_dot_vec(INTRA *actintra,DIST_VECTOR *dist_vec1,
                     DIST_VECTOR *dist_vec2,double *dot)
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
 |  vec[i] = vec[i] * scalar                                            |
 |  DIST_VECTOR *dist_vec (i/o) vector to be multiplied by scalar       |
 |  double       scalar   (o)   scalar value                            |
 *----------------------------------------------------------------------*/
void solserv_scalarprod_vec(DIST_VECTOR *dist_vec,double scalar)
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
 |  double *fullvec (o) vector of lenght numeq_total will be holding    |
 |                      the values from distvec in correct dof-ordering:|
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_reddistvec(DIST_VECTOR *distvec,SPARSE_ARRAY *sysarray,
                        SPARSE_TYP *sysarray_typ,double *fullvec,
                        int dim,INTRA *actintra)
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
if (recv.Typ != cca_DV) 
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
      fullvec[dof] = distvec->vec.a.dv[i];
   }
#ifdef PARALLEL 
   MPI_Allreduce(fullvec,recvbuff,dim,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   for (i=0; i<dim; i++) fullvec[i] = recvbuff[i];
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
 |  double *fullvec (o) vector of lenght numeq_total  holding           |
 |                      the values  correct dof-ordering:               |
 |                      fullvec[dof] = value of a certain dof           |
 |  INTRA *actintra (i) intra-communicator the DIST_VECTOR lives on     |
 *----------------------------------------------------------------------*/
void solserv_distribdistvec(DIST_VECTOR  *distvec,SPARSE_ARRAY *sysarray,
                            SPARSE_TYP *sysarray_typ,double *fullvec,
                            int dim,INTRA *actintra)
{
int             i;
int             dof,dofperm;
double         *dfrom;
int             imyrank;
int             inprocs;

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
 |  int place        (i) place in the ARRAY node->sol where to put the  |
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
                          int place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ)
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
                          int place,SPARSE_ARRAY *sysarray,SPARSE_TYP *sysarray_typ)
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
                          int place,SPARSE_ARRAY *sysarray,
                          SPARSE_TYP *sysarray_typ)
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

/*---------------------------------------------------------------------*
 | routine to find the maximum value of a distributed vector           |
 | ab =  0 absolut maximum value                                       |
 | ab =  1 maxium value                                                |
 | ab = -1 minimum value                                               |
 |                                                         genk 03/02  |
 *---------------------------------------------------------------------*/ 
/*void solserv_dmax_distvec(
			  DIST_VECTOR  *distvec,
			  double *res,   /* result */
/*			  int ab        /* flag */
/*			  )
{
int i;
int numeq;
double dm;

#ifdef DEBUG 
dstrc_enter("solserv_dmax_distvec");
#endif

numeq=distvec->numeq;

switch(ab)
{
case 1:      /* maximum value */
/*dm=distvec->vec.a.dv[0];
for(i=1;i<numeq;i++)
{
   dm=DMAX(dm,distvec->vec.a.dv[i]);
}
/*------------- get maximum value from all procs to proc 0 */
/*#ifdef PARALLEL
MPI_Reduce(*dm,res,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
#else
*res=dm;
#endif
break;

case 0:    /* maximum absolut value */
/*dm=FABS(distvec->vec.a.dv[0]);
for(i=1;i<numeq;i++)
{
   dm=DMAX(dm,FABS(distvec->vec.a.dv[i]));
}
/*------------- get maximum value from all procs to proc 0 */
/*#ifdef PARALLEL
MPI_Reduce(*dm,res,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
#else
*res=dm;
#endif
break;

case -1:   /* minimum value */
/*dm=distvec->vec.a.dv[0];
for(i=1;i<numeq;i++)
{
   dm=DMIN(dm,distvec->vec.a.dv[i]);
}
/*------------- get minimum value from all procs to proc 0 */
/*#ifdef PARALLEL
MPI_Reduce(*dm,res,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
#else
*res=dm;
#endif
break;
default:
   dserror("flag ab wrong!");
}
/*--------------------- distribute maximum / minumum value to all procs */
/*#ifdef PARALLEL
MPI_Bcast(res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
/*----------------------------------------------------------------------*/
/*#ifdef DEBUG 
dstrc_exit();
#endif
return;
}
/* end of solverv_dmax_distvec */ 
 
 



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | dirichlet conditions are scaled by scale and written to sol in the   |
 | given place place                                                    |
 | This routine takes values from the structure actnode->gnode->dirich  |
 | and scales them by a given factor scale. Then it writes these values |
 | the ARRAY node->sol in the given place                               |
 | actnode->sol.a.da[place][j] = actnode->gnode->dirich_val.a.dv[j]*scale|
 | Nothing is done for dofs or nodes which do not have a dirichlet      |
 | condition                                                            |
 | FIELD *actfield (i) active field                                     |
 | int    disnum   (i) indize of the discretization in actfield to be used|
 | double scale    (i) scaling factor for dirichlet condition           |
 | int place       (i) row to put values in the ARRAY sol               |
 *----------------------------------------------------------------------*/
void solserv_putdirich_to_dof(FIELD *actfield, int disnum, double scale, 
                              int place)
{
int               i,j;
int               diff,max;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG 
dstrc_enter("solserv_putdirich_to_dof");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*--------- do nothing if there is no dirichlet condition on actnode */
   if (actnode->gnode->dirich==NULL) continue;
   /*------------------------------------------ get dirichlet condition */
   dirich = actnode->gnode->dirich;
   /* of the given place is outside the dimensions of ARRAY sol enlarge it */
   if (place >= actnode->sol.fdim)
   {
      diff = place - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max+1,actnode->sol.sdim,"DA");
   }
   /* put values dirich->dirich_val.a.dv[j] * scale to actnode->sol.a.dv[place] */
   for (j=0; j<actnode->numdf; j++)
   {
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
      actnode->sol.a.da[place][j] = dirich->dirich_val.a.dv[j] * scale;
/* this is special for the ortiz example */
/*      actnode->sol.a.da[place][j] *= actnode->x[0];*/
   }   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_putdirich_to_dof */



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | entries in sol in the place  placefrom1 and placefrom2 are added     |
 | node->sol[to][..] = node->sol[from1][..] * facfrom1 +                |
 |                     node->sol[from2][..] * facfrom2                  |
 |                                                                      |
 | This is ONLY performed for dofs which have a dirichlet condition     |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | int    disnum   (i) indize of the discretization in actfield to be used|
 | int    from1    (i) place in ARRAY sol to take the values from       |
 | int    from2    (i) place in ARRAY sol to take the values from       |
 | int    to       (i) place in ARRAY sol to write values to            |
 | double facfrom1 (i) scaling factor vor values from from1             |
 | double facfrom2 (i) scaling factor vor values from from2             |
 *----------------------------------------------------------------------*/
void solserv_adddirich(FIELD *actfield, int disnum,
                              int from1,int from2,int to,
                              double facfrom1, double facfrom2)
{
int               i,j;
int               diff,max;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG 
dstrc_enter("solserv_adddirich");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   dirich = actnode->gnode->dirich;
   max = IMAX(from1,from2);
   max = IMAX(max,to);
   if (max >= actnode->sol.fdim)
   {
      diff = max - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max+1,actnode->sol.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
      actnode->sol.a.da[to][j] = 
      actnode->sol.a.da[from1][j]*facfrom1 + 
      actnode->sol.a.da[from2][j]*facfrom2;
   }   
   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_adddirich */


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | entries in sol in the place  placefrom1 and placefrom2 are added     |
 | node->sol[to][..] += node->sol[from1][..] * facfrom1 +               |
 |                      node->sol[from2][..] * facfrom2                 |
 |                                                                      |
 | same functionality as solserv_adddirich but adds to sol in the place to|
 |                                                                      |
 *----------------------------------------------------------------------*/
void solserv_assdirich_fac(FIELD *actfield, int disnum,
                           int from1,int from2,int to, 
                           double facfrom1, double facfrom2)
{
int               i,j;
int               diff,max;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG 
dstrc_enter("solserv_assdirich_fac");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
/*------------ loop nodes and put the result back to the node structure */
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   dirich = actnode->gnode->dirich;
   max = IMAX(from1,from2);
   max = IMAX(max,to);
   if (max >= actnode->sol.fdim)
   {
      diff = max - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max+1,actnode->sol.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
      actnode->sol.a.da[to][j] += 
      actnode->sol.a.da[from1][j]*facfrom1 + 
      actnode->sol.a.da[from2][j]*facfrom2;
   }   
   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_assdirich_fac */



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | entries in sol in the place  placefrom1  are copied to  place to     |
 | node->sol[to][..] = node->sol[from][..]                              |
 |                                                                      |
 | This is only performed for dofs which have a dirichlet condition on them |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | int    disnum   (i) indize of the discretization in actfield to be used|
 | int    to       (i) place in ARRAY sol to write values to            |
 | int    from     (i) place in ARRAY sol to take the values from       |
 *----------------------------------------------------------------------*/
void solserv_cpdirich(FIELD *actfield, int disnum,
                      int from,int to)
{
int               i,j;
int               diff,max;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG 
dstrc_enter("solserv_cpdirich");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   dirich = actnode->gnode->dirich;
   max = IMAX(from,to);
   if (max >= actnode->sol.fdim)
   {
      diff = max - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max+1,actnode->sol.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
      actnode->sol.a.da[to][j] = actnode->sol.a.da[from][j];
   }   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_cpdirich */



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/02|
 | init sol[place] to zero                                              |
 |                                                                      |
 | This is only performed for dofs which have a dirichlet condition on them |
 |                                                                      |
 | FIELD *actfield (i) active field                                     |
 | int    disnum   (i) indize of the discretization in actfield to be used|
 | int    place    (i) row in ARRAY sol to be set to zero               | 
 *----------------------------------------------------------------------*/
void solserv_zerodirich(FIELD *actfield, int disnum, int place)
{
int               i,j;
int               diff,max;
NODE             *actnode;
DISCRET          *actdis;
DIRICH_CONDITION *dirich;
#ifdef DEBUG 
dstrc_enter("solserv_zerodirich");
#endif
/*----------------------------------------------------------------------*/
actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->gnode->dirich==NULL) continue;
   dirich = actnode->gnode->dirich;
   if (place >= actnode->sol.fdim)
   {
      diff = place - actnode->sol.fdim;
      max  = IMAX(diff,5);
      amredef(&(actnode->sol),actnode->sol.fdim+max+1,actnode->sol.sdim,"DA");
   }
   for (j=0; j<actnode->numdf; j++)
   {
      if (dirich->dirich_onoff.a.iv[j]==0) continue;
      actnode->sol.a.da[place][j] = 0.0;
   }   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_zerodirich */


