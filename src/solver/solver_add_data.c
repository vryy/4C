#include "../headers/standardtypes.h"
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
                 int                sysarray1,
                 struct _ARRAY     *elearray1,
                 int                sysarray2,
                 struct _ARRAY     *elearray2,
                 struct _PARTITION *actpart,
                 struct _SOLVAR    *actsolv,
                 struct _INTRA     *actintra,
                 struct _ELEMENT   *actele,
                 int                option
                )
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
/*----------------------- option==0 is normal assembly of system matrix */
if (!option)/*------------------------- option=0 is the normal assembly */
{
/*------------------------------------------------ switch typ of matrix */
/*-------------------------------------------- add to 2 system matrices */
   if (sysarray1>=0 && sysarray2>=0) 
   switch(sysa1_typ)
   {
   case msr:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case parcsr:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case ucchb:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case dense:
      dserror("Simultanous assembly of 2 system matrices not yet impl.");
   break;
   case sparse_none:
      dserror("Unspecified typ of system matrix");
   break;
   default:
      dserror("Unspecified typ of system matrix");
   break;
   }
/*--------------------------------------------- add to 1 system matrix */
   if (sysarray1>=0) 
   switch(sysa1_typ)
   {
   case msr:
      add_msr(actpart,actsolv,actintra,actele,sysa1->msr);
   break;
   case parcsr:
      add_parcsr(actpart,actsolv,actintra,actele,sysa1->parcsr);
   break;
   case ucchb:
      add_ucchb(actpart,actsolv,actintra,actele,sysa1->ucchb);
   break;
   case dense:
      add_dense(actpart,actsolv,actintra,actele,sysa1->dense);
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
else
{
/*------------------------------------------------ switch typ of matrix */
/*----------------------------------------exchange of 2 system matrices */
      if (sysarray1>=0 && sysarray2>=0) 
      switch(sysa1_typ)
      {
      case msr:
         dserror("Simultanous assembly of 2 system matrices not yet impl.");
      break;
      case parcsr:
         dserror("Simultanous assembly of 2 system matrices not yet impl.");
      break;
      case ucchb:
         dserror("Simultanous assembly of 2 system matrices not yet impl.");
      break;
      case sparse_none:
         dserror("Unspecified typ of system matrix");
      break;
      default:
         dserror("Unspecified typ of system matrix");
      break;
      }
/*-------------------------------------------exchange of 1 system matrix */
      if (sysarray1>=0) 
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
         redundant_dense(actpart,actsolv,actintra,sysa1->dense);
      break;
      case sparse_none:
         dserror("Unspecified typ of system matrix");
      break;
      default:
         dserror("Unspecified typ of system matrix");
      break;
      }
}
#endif
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
 *----------------------------------------------------------------------*/
void init_assembly(
                       struct _PARTITION      *actpart,
                       struct _SOLVAR         *actsolv,
                       struct _INTRA          *actintra,
                       struct _FIELD          *actfield,
                       int                     actsysarray
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
default:
   dserror("Unknown typ of sparse array");
break;
}
/*---------------- now check for coupling dofs and interdomain coupling */
coupledofs = &(actpart->coupledofs);
numsend = 0;
numrecv = 0;
numeq   = actfield->numeq;
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
   /*------------------------------------ check for master owner of dof */
   if (coupledofs->a.ia[i][imyrank+1]==2)
   {
      /*------------------------- check whether other procs have slaves */
      for (j=1; j<coupledofs->sdim; j++)
      {
         if (coupledofs->a.ia[i][j]==1) numrecv++;
      }
   }
   /*------------------------------------ check for slave owners of dof */
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
   *couple_d_send_ptr = (ARRAY*)calloc(1,sizeof(ARRAY));
   *couple_i_send_ptr = (ARRAY*)calloc(1,sizeof(ARRAY));

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
   *couple_d_recv_ptr = (ARRAY*)calloc(1,sizeof(ARRAY));
   *couple_i_recv_ptr = (ARRAY*)calloc(1,sizeof(ARRAY));

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
#endif /* end of PARALLEL */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of init_assembly */





/*----------------------------------------------------------------------*
 |  routine to assemble nodal neumann conditions         m.gee 10/01    |
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
AZ_ARRAY_MSR         *msr_array;
H_PARCSR             *parcsr_array;
UCCHB                *ucchb_array;
DENSE                *dense_array;
#ifdef DEBUG 
dstrc_enter("assemble_vec");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
/*----------------------------------------------------------------------*/
switch(*sysarraytyp)
{
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


