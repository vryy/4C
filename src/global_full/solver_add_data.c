#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;

/*----------------------------------------------------------------------*
 |  routine to assemble element arrays to global sparse arrays m.gee 9/01|
 *----------------------------------------------------------------------*/
void assemble(
                 int                    sysarray1, /* number of first sparse system matrix */
                 struct _ARRAY         *elearray1, /* pointer to first dense element matrix */
                 int                    sysarray2, /* number of first sparse system matrix or -1 if not given */
                 struct _ARRAY         *elearray2, /* pointer to second dense element matrix or NULL is not present*/
                 struct _PARTITION     *actpart,   /* my partition of theactive field */
                 struct _SOLVAR        *actsolv,   /* the active SOLVAR */
                 struct _INTRA         *actintra,  /* the active intracommunicator */
                 struct _ELEMENT       *actele,    /* the element to assemble */
                 enum _ASSEMBLE_ACTION  assemble_action,  /* the assembly option */
                 CONTAINER             *container  /* contains variables defined in container.h */
                )
/*----------------------------------------------------------------------*/
{
int         i,j,k;
enum  _SPARSE_TYP    sysa1_typ;
union _SPARSE_ARRAY *sysa1;
enum  _SPARSE_TYP    sysa2_typ;
union _SPARSE_ARRAY *sysa2;
#ifdef DEBUG 
dstrc_enter("assemble");
#endif
/*----------------------------------------------------------------------*/
if (assemble_action==assemble_do_nothing) goto end;
/*----------------------- check for presence and typ of system matrices */
if (sysarray1>=0) 
{
   sysa1       = &(actsolv->sysarray[sysarray1]);
   sysa1_typ   =   actsolv->sysarray_typ[sysarray1]; 
}
else              
{
   sysa1     = NULL;
   sysa1_typ = sparse_none; 
}
if (sysarray2>=0) 
{
   sysa2     = &(actsolv->sysarray[sysarray2]);
   sysa2_typ =   actsolv->sysarray_typ[sysarray2];
}
else              
{
   sysa2     = NULL;
   sysa2_typ = sparse_none;
}
/*----------------------------------------------------------------------*/
if (assemble_action==assemble_two_matrix)
{
   if (sysa1_typ != sysa2_typ)
   dserror("Assembly of element matrices in different types of sparse mat. not impl.");
/*------------------------------------------------ switch typ of matrix */
/*-------------------------------------------- add to 2 system matrices */
   switch(sysa1_typ)
   {
   case mds:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case msr:
      add_msr(actpart,actsolv,actintra,actele,sysa1->msr,sysa2->msr);
   break;
   case parcsr:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case ucchb:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case dense:
      add_dense(actpart,actsolv,actintra,actele,sysa1->dense,sysa2->dense);
   break;
   case rc_ptr:
      add_rc_ptr(actpart,actsolv,actintra,actele,sysa1->rc_ptr,sysa2->rc_ptr);
   break;
   case ccf:
      add_ccf(actpart,actsolv,actintra,actele,sysa1->ccf,sysa2->ccf);
   break;
   case skymatrix:
      add_skyline(actpart,actsolv,actintra,actele,sysa1->sky,sysa2->sky);
   break;
   case spoolmatrix:
      add_spo(actpart,actsolv,actintra,actele,sysa1->spo,sysa2->spo);
   break;
   case bdcsr:
      add_bdcsr(actpart,actsolv,actintra,actele,sysa1->bdcsr,sysa2->bdcsr);
   break;
   case oll:
      add_oll(actpart,actintra,actele,sysa1->oll,sysa2->oll);
   break;
   case sparse_none:
      dserror("Unspecified type of system matrix");
   break;
   default:
      dserror("Unspecified type of system matrix");
   break;
   }
}
/*--------------------------------------------- add to 1 system matrix */
if (assemble_action==assemble_one_matrix)
{
   switch(sysa1_typ)
   {
   case mds:
      add_mds(actpart,actsolv,actele,sysa1->mds);
   break;
   case msr:
      add_msr(actpart,actsolv,actintra,actele,sysa1->msr,NULL);
   break;
   case parcsr:
      add_parcsr(actpart,actsolv,actintra,actele,sysa1->parcsr);
   break;
   case ucchb:
      add_ucchb(actpart,actsolv,actintra,actele,sysa1->ucchb);
   break;
   case dense:
      add_dense(actpart,actsolv,actintra,actele,sysa1->dense,NULL);
   break;
   case rc_ptr:
      add_rc_ptr(actpart,actsolv,actintra,actele,sysa1->rc_ptr,NULL);
   break;
   case ccf:
      add_ccf(actpart,actsolv,actintra,actele,sysa1->ccf,NULL);
   break;
   case skymatrix:
      add_skyline(actpart,actsolv,actintra,actele,sysa1->sky,NULL);
   break;
   case spoolmatrix:
      add_spo(actpart,actsolv,actintra,actele,sysa1->spo,NULL);
   break;
   case bdcsr:
      add_bdcsr(actpart,actsolv,actintra,actele,sysa1->bdcsr,NULL);
   break;
   case oll:
      add_oll(actpart,actintra,actele,sysa1->oll,NULL);
   break;
   case sparse_none:
      dserror("Unspecified typ of system matrix");
   break;
   default:
      dserror("Unspecified typ of system matrix");
   break;
   }
}
/*----------------- close the system matrix, or close two system matices */
if (assemble_action==assemble_close_1matrix)
{
   switch(sysa1_typ)
   {
   case mds:
   break;
   case msr:
   break;
   case parcsr:
   break;
   case ucchb:
   break;
   case dense:
   break;
   case rc_ptr:
   break;
   case ccf:
   break;
   case skymatrix:
   break;
   case spoolmatrix:
      close_spooles_matrix(sysa1->spo,actintra);
   break;
   case oll:
   break;
   case sparse_none:
      dserror("Unspecified typ of system matrix");
   break;
   default:
      dserror("Unspecified typ of system matrix");
   break;
   }
}
if (assemble_action==assemble_close_2matrix)
{
   switch(sysa1_typ)
   {
   case mds:
   break;
   case msr:
   break;
   case parcsr:
   break;
   case ucchb:
   break;
   case dense:
   break;
   case rc_ptr:
   break;
   case ccf:
   break;
   case skymatrix:
   break;
   case spoolmatrix:
      close_spooles_matrix(sysa1->spo,actintra);
      close_spooles_matrix(sysa2->spo,actintra);
   break;
   case oll:
   break;
   case sparse_none:
      dserror("Unspecified typ of system matrix");
   break;
   default:
      dserror("Unspecified typ of system matrix");
   break;
   }
}
/*-------------- option==1 is exchange of coupled dofs among processors */
/*                    (which, of course, only occures in parallel case) */
#ifdef PARALLEL 
/*----------------------------------------exchange of 2 system matrices */
if (assemble_action==assemble_two_exchange)
{
/*------------------------------------------------ switch typ of matrix */
      switch(sysa1_typ)
      {
      case msr:
         exchange_coup_msr(actpart,actsolv,actintra,sysa1->msr);
         exchange_coup_msr(actpart,actsolv,actintra,sysa2->msr);
      break;
      case parcsr:
         dserror("Simultanous assembly of 2 system matrices not yet impl.");
      break;
      case ucchb:
         dserror("Simultanous assembly of 2 system matrices not yet impl.");
      break;
      case dense:
         redundant_dense(actpart,actsolv,actintra,sysa1->dense,sysa2->dense);
      break;
      case rc_ptr:
         exchange_coup_rc_ptr(actpart,actsolv,actintra,sysa1->rc_ptr);
         exchange_coup_rc_ptr(actpart,actsolv,actintra,sysa2->rc_ptr);
      break;
      case spoolmatrix:
         exchange_coup_spo(actpart,actsolv,actintra,sysa1->spo);
         exchange_coup_spo(actpart,actsolv,actintra,sysa2->spo);
      break;
      case ccf:
         redundant_ccf(actpart,actsolv,actintra,sysa1->ccf,sysa2->ccf);
      break;
      case skymatrix:
         redundant_skyline(actpart,actsolv,actintra,sysa1->sky,sysa2->sky);
      break;
      case bdcsr:;
      break;
      case oll:
         exchange_coup_oll(actpart,actintra,sysa1->oll);
         exchange_coup_oll(actpart,actintra,sysa2->oll);
      break;
      case sparse_none:
         dserror("Unspecified type of system matrix");
      break;
      default:
         dserror("Unspecified type of system matrix");
      break;
      }
}
/*-------------------------------------------exchange of 1 system matrix */
if (assemble_action==assemble_one_exchange)
{
      switch(sysa1_typ)
      {
      case msr:
         exchange_coup_msr(actpart,actsolv,actintra,sysa1->msr);
      break;
      case parcsr:
         exchange_coup_parcsr(actpart,actsolv,actintra,sysa1->parcsr);
      break;
      case ucchb:
         redundant_ucchb(actpart,actsolv,actintra,sysa1->ucchb);
      break;
      case dense:
         redundant_dense(actpart,actsolv,actintra,sysa1->dense,NULL);
      break;
      case skymatrix:
         redundant_skyline(actpart,actsolv,actintra,sysa1->sky,NULL);
      break;
      case rc_ptr:
         exchange_coup_rc_ptr(actpart,actsolv,actintra,sysa1->rc_ptr);
      break;
      case spoolmatrix:
         exchange_coup_spo(actpart,actsolv,actintra,sysa1->spo);
      break;
      case ccf:
         redundant_ccf(actpart,actsolv,actintra,sysa1->ccf,NULL);
      break;
      case bdcsr:;
      break;
      case oll:
         exchange_coup_oll(actpart,actintra,sysa1->oll);
      break;
      case sparse_none:
         dserror("Unspecified type of system matrix");
      break;
      default:
         dserror("Unspecified type of system matrix");
      break;
      }
}
#endif
/*------------------------------------- close the dynamic system matrix */
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble */




/*----------------------------------------------------------------------*
 |  routine to                                           m.gee 9/01     |
 |  allocate the send and recv buffers for                              |
 |  coupling conditions                                                 |
 |  and to perform other inits which may become necessary for assembly  |
 | #################################################################### |
 |  to the paremeter list of this function I added                      |
 |  int actndis!!!  - number of the actual discretisation               |
 |  this has to be done for all other calls of init_assembly            |
 |                                                           genk 08/02 |
 | #################################################################### |
 *----------------------------------------------------------------------*/
void init_assembly(
                       struct _PARTITION      *actpart,
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       struct _FIELD          *actfield,
                       int                     actsysarray,
		       int                     actndis
                     )
{
int         i,j,k,counter;
int         numeq;
int         numsend;
int         numrecv;
int         minusone=-1;
int         imyrank;
int         inprocs;
SPARSE_TYP  sysarraytyp;
ARRAY      *coupledofs;
ELEMENT    *actele;

int        *numcoupsend;
int        *numcouprecv;
ARRAY     **couple_d_send_ptr;
ARRAY     **couple_i_send_ptr;
ARRAY     **couple_d_recv_ptr;
ARRAY     **couple_i_recv_ptr;
ARRAY      *dummyarray;

#ifdef DEBUG 
dstrc_enter("init_assembly");
#endif
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------- check typ of sysarray */
sysarraytyp = actsolv->sysarray_typ[actsysarray];
switch(sysarraytyp)
{
case msr:
   numcoupsend       = &(actsolv->sysarray[actsysarray].msr->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].msr->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].msr->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].msr->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].msr->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].msr->couple_i_recv);
break;
case parcsr:
   numcoupsend       = &(actsolv->sysarray[actsysarray].parcsr->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].parcsr->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].parcsr->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].parcsr->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].parcsr->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].parcsr->couple_i_recv);
break;
case ucchb:
   numcoupsend       = &(actsolv->sysarray[actsysarray].ucchb->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].ucchb->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].ucchb->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].ucchb->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].ucchb->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].ucchb->couple_i_recv);
break;
case dense:
   numcoupsend       = &(actsolv->sysarray[actsysarray].dense->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].dense->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].dense->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].dense->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].dense->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].dense->couple_i_recv);
break;
case rc_ptr:
   numcoupsend       = &(actsolv->sysarray[actsysarray].rc_ptr->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].rc_ptr->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].rc_ptr->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].rc_ptr->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].rc_ptr->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].rc_ptr->couple_i_recv);
break;
case ccf:
   numcoupsend       = &(actsolv->sysarray[actsysarray].ccf->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].ccf->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].ccf->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].ccf->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].ccf->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].ccf->couple_i_recv);
break;
case skymatrix:
   numcoupsend       = &(actsolv->sysarray[actsysarray].sky->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].sky->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].sky->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].sky->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].sky->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].sky->couple_i_recv);
break;
case spoolmatrix:
   numcoupsend       = &(actsolv->sysarray[actsysarray].spo->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].spo->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].spo->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].spo->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].spo->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].spo->couple_i_recv);
break;
case oll:
   numcoupsend       = &(actsolv->sysarray[actsysarray].oll->numcoupsend);
   numcouprecv       = &(actsolv->sysarray[actsysarray].oll->numcouprecv);
   couple_d_send_ptr = &(actsolv->sysarray[actsysarray].oll->couple_d_send);
   couple_i_send_ptr = &(actsolv->sysarray[actsysarray].oll->couple_i_send);
   couple_d_recv_ptr = &(actsolv->sysarray[actsysarray].oll->couple_d_recv);
   couple_i_recv_ptr = &(actsolv->sysarray[actsysarray].oll->couple_i_recv);
break;
case bdcsr:
   goto end; /* coupled dofs are not supported in bdcsr */
default:
   dserror("Unknown typ of sparse array");
break;
}
/*---------------- now check for coupling dofs and interdomain coupling */
coupledofs = &(actpart->pdis[actndis].coupledofs);
numsend = 0;
numrecv = 0;
numeq   = actfield->dis[actndis].numeq;
/* 
   An inter-proc coupled equation produces communications calculating the 
   sparsity mask of the matrix
   An inter-proc coupled equation produces communications adding element
   matrices to the system matrix
   An inter-proc coupled equation ruins the bandwith locally
   ->
   Now one processor has to be owner of the coupled equation. 
   Try to distribute the coupled equations equally over the processors

   The matrix has the following style (after allreduce on all procs the same):
   
               ----------------------
               | 12 | 2 | 0 | 1 | 0 |
               | 40 | 2 | 0 | 0 | 0 |
               | 41 | 1 | 2 | 1 | 1 |
               | 76 | 0 | 1 | 2 | 0 |
               ----------------------
               
               column 0                : number of the coupled equation
               column 1 - inprocs+1 : proc has coupled equation or not
                                         2 indicates owner of equation
*/
/* calculate the number of sends and receives to expect during assemblage */
for (i=0; i<coupledofs->fdim; i++)
{
   /*--------------------------- check whether I am master owner of dof */
   if (coupledofs->a.ia[i][imyrank+1]==2)
   {
      /*-------------------------- check whether other procs are slaves */
      for (j=1; j<coupledofs->sdim; j++)
      {
         if (coupledofs->a.ia[i][j]==1) numrecv++;
      }
   }
   /*---------------------------- check whether I am slave owner of dof */
   if (coupledofs->a.ia[i][imyrank+1]==1) numsend++;
}
*numcoupsend=numsend;
*numcouprecv=numrecv;
/*-------------------------- allocate the necessary send and recv buffs */
/* note:
   Es waere sinnvoll, die sends und recv in einem Matrizenkompressionsformat
   durchzufuehren. Damit waere aber die Art dieses send + recv von dem
   verwendeten Loeser abhaengig, und man muesste das fuer jeden Loeser
   gesondert schreiben. Weil es fuer alle Loeser funktioniert wird hier als
   send und recv buffer eine komplette Zeile der Systemmatrix pro coupled dof 
   verwendet 
*/
if (numsend) /*-------- I have to send couple dof entries to other proc */
{
   *couple_d_send_ptr = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
   *couple_i_send_ptr = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));

   if (!(*couple_d_send_ptr) || !(*couple_i_send_ptr))
   dserror("Allocation of send/recv buffers for coupled dofs failed");

   amdef("c_d_send",(*couple_d_send_ptr),numsend,numeq,"DA");
   amdef("c_i_send",(*couple_i_send_ptr),numsend,2,"IA");

   aminit(*couple_i_send_ptr,&minusone);

   /*----------------------- put the dof number to couple_i_send[0] */
   counter=0;
   for (i=0; i<coupledofs->fdim; i++)
   {
       if (coupledofs->a.ia[i][imyrank+1]==1)
       {
          (*couple_i_send_ptr)->a.ia[counter][0] = coupledofs->a.ia[i][0];
          counter++;
       }
   }
} 
else /*----------------------------------------- I have nothing to send */
{
   *couple_d_send_ptr = NULL;
   *couple_i_send_ptr = NULL;
}


if (numrecv) /* I am master of a coupled dof and expect entries from other procs */
{
   *couple_d_recv_ptr = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
   *couple_i_recv_ptr = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));

   if (!(*couple_d_recv_ptr) || !(*couple_i_recv_ptr))
   dserror("Allocation of send/recv buffers for coupled dofs failed");

   amdef("c_d_recv",(*couple_d_recv_ptr),numrecv,numeq,"DA");
   amdef("c_i_recv",(*couple_i_recv_ptr),numrecv,2,"IA");
}
else /*----------------------- I do not expect entries from other procs */
{
   *couple_d_recv_ptr = NULL;
   *couple_i_recv_ptr = NULL;
}
/*----------------------------------------------------------------------*/
end:
#endif /* end of PARALLEL */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of init_assembly */





/*----------------------------------------------------------------------*
 |  routine to assemble a global vector to a dist. vector   m.gee 10/01 |
 |  irhs & drhs are vectors of global lenght, rhs is a DIST_VECTOR      |
 |  and is filled in a style, which is very much dependent on the       |
 |  typ of system matrix it belongs to                                  |
 *----------------------------------------------------------------------*/
void assemble_vec(INTRA        *actintra,
                  SPARSE_TYP   *sysarraytyp,
                  SPARSE_ARRAY *sysarray,
                  DIST_VECTOR  *rhs,
                  double       *drhs,
                  double        factor)
{
int                   i;
int                   dof;
int                   imyrank;
ML_ARRAY_MDS         *mds_array;
AZ_ARRAY_MSR         *msr_array;
H_PARCSR             *parcsr_array;
UCCHB                *ucchb_array;
DENSE                *dense_array;
RC_PTR               *rcptr_array;
CCF                  *ccf_array;
SKYMATRIX            *sky_array;
SPOOLMAT             *spo;
DBCSR                *bdcsr_array;
OLL                  *oll_array;
#ifdef DEBUG 
dstrc_enter("assemble_vec");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
/*----------------------------------------------------------------------*/
switch(*sysarraytyp)
{
case mds:
    mds_array = sysarray->mds;
    for (i=0; i<rhs->numeq; i++)
    {
       rhs->vec.a.dv[i] += drhs[i]*factor;
    }
break;
case msr:
    msr_array = sysarray->msr;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = msr_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case parcsr:
    parcsr_array = sysarray->parcsr;
    for (i=0; i<rhs->numeq; i++)
    {
       dof     = parcsr_array->update.a.ia[imyrank][i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case ucchb:
    ucchb_array = sysarray->ucchb;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = ucchb_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case dense:
    dense_array = sysarray->dense;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = dense_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case skymatrix:
    sky_array = sysarray->sky;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = sky_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case rc_ptr:
    rcptr_array = sysarray->rc_ptr;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = rcptr_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case ccf:
    ccf_array = sysarray->ccf;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = ccf_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case spoolmatrix:
    spo = sysarray->spo;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = spo->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case bdcsr:
    bdcsr_array = sysarray->bdcsr;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = bdcsr_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
case oll:
    oll_array = sysarray->oll;
    for (i=0; i<rhs->numeq; i++)
    {
       dof = oll_array->update.a.iv[i];
       rhs->vec.a.dv[i] += drhs[dof]*factor;
    }
break;
default:
   dserror("Unknown typ of system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble_vec */


/*----------------------------------------------------------------------*
 |  routine to sum a global vector                          m.gee 6/02 |
 *----------------------------------------------------------------------*/
void sum_vec(INTRA        *actintra,
             SPARSE_TYP   *sysarraytyp,
             SPARSE_ARRAY *sysarray,
             double       *drhs,
             int           numeq,
             double       *sum)
{
int                   i;
int                   dof;
int                   imyrank;
double                recv;
ML_ARRAY_MDS         *mds_array;
AZ_ARRAY_MSR         *msr_array;
H_PARCSR             *parcsr_array;
UCCHB                *ucchb_array;
DENSE                *dense_array;
RC_PTR               *rcptr_array;
CCF                  *ccf_array;
SKYMATRIX            *sky_array;
SPOOLMAT             *spo;
OLL                  *oll_array;
#ifdef DEBUG 
dstrc_enter("sum_vec");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
/*----------------------------------------------------------------------*/
*sum=0.0;
/*----------------------------------------------------------------------*/
switch(*sysarraytyp)
{
case mds:
    mds_array = sysarray->mds;
    for (i=0; i<numeq; i++)
    {
       *sum += drhs[i];
    }
break;
case msr:
    msr_array = sysarray->msr;
    for (i=0; i<numeq; i++)
    {
       dof = msr_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case parcsr:
    parcsr_array = sysarray->parcsr;
    for (i=0; i<numeq; i++)
    {
       dof     = parcsr_array->update.a.ia[imyrank][i];
       *sum += drhs[dof];
    }
break;
case ucchb:
    ucchb_array = sysarray->ucchb;
    for (i=0; i<numeq; i++)
    {
       dof = ucchb_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case dense:
    dense_array = sysarray->dense;
    for (i=0; i<numeq; i++)
    {
       dof = dense_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case skymatrix:
    sky_array = sysarray->sky;
    for (i=0; i<numeq; i++)
    {
       dof = sky_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case rc_ptr:
    rcptr_array = sysarray->rc_ptr;
    for (i=0; i<numeq; i++)
    {
       dof = rcptr_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case ccf:
    ccf_array = sysarray->ccf;
    for (i=0; i<numeq; i++)
    {
       dof = ccf_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case spoolmatrix:
    spo = sysarray->spo;
    for (i=0; i<numeq; i++)
    {
       dof = spo->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
case oll:
    oll_array = sysarray->oll;
    for (i=0; i<numeq; i++)
    {
       dof = oll_array->update.a.iv[i];
       *sum += drhs[dof];
    }
break;
default:
   dserror("Unknown typ of system matrix");
break;
}
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
MPI_Allreduce(sum,&recv,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
*sum = recv;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of sum_vec */




/*----------------------------------------------------------------------*
 |  assembles an element vector to a redundant global vector m.gee 3/02 |
 *----------------------------------------------------------------------*/
void assemble_intforce(ELEMENT *actele,ARRAY *elevec_a,CONTAINER *container,
                       INTRA *actintra)
{
int                   i,j;
int                   dof;
int                   numdf;
int                   imyrank;
int                   irow;
double               *elevec;
#ifdef DEBUG 
dstrc_enter("assemble_intforce");
#endif
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
imyrank = actintra->intra_rank;
#endif
/*----------------------------------------------------------------------*/
elevec = elevec_a->a.dv;
irow=-1;
/*----------------------------------------------------------------------*/
for (i=0; i<actele->numnp; i++)
{
   numdf = actele->node[i]->numdf;
#ifdef PARALLEL 
   if(actele->node[i]->proc!= imyrank)
   {
     irow+=numdf;
     continue;
   }
#endif
   for (j=0; j<numdf; j++)
   {
      irow++;
      dof = actele->node[i]->dof[j];
      if (dof >= container->global_numeq) continue;
       container->dvec[dof] += elevec[irow];
/*      container->dvec[dof] += elevec[i*numdf+j]; */
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble_intforce */


/*---------------------------------------------------------------------------*
 |  dirichlet conditions to an element vector elevec_a       m.gee 3/02      |
 |  and then assembles this element vector of cond. dirich.conditions to the |
 |  global vector fullvec                                                    |
 *---------------------------------------------------------------------------*/
void assemble_dirich(ELEMENT *actele, ARRAY *estif_global, CONTAINER *container)
{
int                   i,j;
int                   dof;
int                   numdf;
int                   iel;
int                   nd=0;
double              **estif;
double                dirich[MAXDOFPERELE];
double                dforces[MAXDOFPERELE];
int                   dirich_onoff[MAXDOFPERELE];
int                   lm[MAXDOFPERELE];
GNODE                *actgnode;
#ifdef DEBUG 
dstrc_enter("assemble_dirich");
#endif
/*----------------------------------------------------------------------*/
estif  = estif_global->a.da;
/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;
/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dforces[i] = 0.0;
   dirich_onoff[i] = 0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
for (i=0; i<actele->numnp; i++)
{
   numdf    = actele->node[i]->numdf;
   actgnode = actele->node[i]->gnode;
   for (j=0; j<numdf; j++)
   {
      lm[i*numdf+j] = actele->node[i]->dof[j];
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[i*numdf+j] = actgnode->dirich->dirich_val.a.dv[j];
   }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */
/*-------- now assemble the vector dforces to the global vector fullvec */
for (i=0; i<nd; i++)
{
   if (lm[i] >= container->global_numeq) continue;
   container->dirich[lm[i]] += dforces[i];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble_dirich */
