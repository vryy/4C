#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  routine to call rhs-routines                         m.gee 10/01    |
 *----------------------------------------------------------------------*/
void calrhs(FIELD        *actfield,     /* the active field */
            SOLVAR       *actsolv,      /* the active SOLVAR */
            PARTITION    *actpart,      /* my partition of this field */
            INTRA        *actintra,     /* the field's intra-communicator */
            int           actsysarray,  /* the active sparse array */
            DIST_VECTOR  *rhs1,         /* 2 dist. vectors for rhs */
            DIST_VECTOR  *rhs2,
            int           kstep,        /* actual time or load incremental step */
            CALC_ACTION  *action)       /* action to be passed to element routines */
{
int i;
SPARSE_TYP   *sysarraytyp;
SPARSE_ARRAY *sysarray;
 #ifdef DEBUG
dstrc_enter("calrhs");
#endif
/*----------------------------------------------------------------------*/
sysarraytyp = &(actsolv->sysarray_typ[actsysarray]);
sysarray    = &(actsolv->sysarray[actsysarray]);            
              
/*---------------------------------make rhs of nodal neumann conditions */
calrhs_nodal_neumann(actpart,
                     actintra,
                     sysarraytyp,
                     sysarray,
                     rhs1);
/*------------------------------ make rhs of element neumann conditions */
/*               (yes I know now, strictly there is no such thing.....) */
calrhs_ele_neumann(actfield,
                   actsolv,
                   actpart,
                   actintra,
                   sysarraytyp,
                   sysarray,
                   actsysarray,
                   rhs2,
                   kstep,
                   action);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calrhs */


 


/*----------------------------------------------------------------------*
 |  routine to assemble nodal neumann conditions         m.gee 10/01    |
 *----------------------------------------------------------------------*/
void calrhs_nodal_neumann(
                           PARTITION    *actpart,    /* my active partition */
                           INTRA        *actintra,   /* the field's communicator */
                           SPARSE_TYP   *sysarraytyp,/* type of sparse matrix */
                           SPARSE_ARRAY *sysarray,   /* sparse matrix number */
                           DIST_VECTOR  *rhs         /* dist. vector */
                          )
{
int                   i,j;               /* some counters */
int                   dof;               
NODE                 *actnode;
COND_NODE            *actcond;

ARRAY                 drhs_send_a;       /* send and receive buffers */
double               *drhs_send;
ARRAY                 drhs_recv_a;
double               *drhs_recv;


#ifdef DEBUG 
dstrc_enter("calrhs");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------ allocate temporary vectors */
/*------------------------------ these vectors are of global dimensions */
drhs_send = amdef("tmp",&drhs_send_a,rhs->numeq_total,1,"DV");
            amzero(&drhs_send_a);
#ifdef PARALLEL 
drhs_recv = amdef("tmp",&drhs_recv_a,rhs->numeq_total,1,"DV");
            amzero(&drhs_recv_a);
#endif
/*---------------------------------------- do assembly to global vector */
assemble_nn(actpart,rhs,drhs_send);
/*----------------------------------------------- allreduce vector drhs */
#ifdef PARALLEL 
MPI_Allreduce(drhs_send,
              drhs_recv,
              rhs->numeq_total,
              MPI_DOUBLE,
              MPI_SUM,
              actintra->MPI_INTRA_COMM);
/*--------------------------------------add data to the DIST_VECTOR rhs */
assemble_vec(actintra,
                sysarraytyp,
                sysarray,
                rhs,
                drhs_recv,
                1.0);
#else
assemble_vec(actintra,
                sysarraytyp,
                sysarray,
                rhs,
                drhs_send,
                1.0);
#endif
/*------------------------------------------------------ delete buffers */
amdel(&drhs_send_a);
#ifdef PARALLEL 
amdel(&drhs_recv_a);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calrhs */




/*----------------------------------------------------------------------*
 |  routine to assemble element neumann conditions       m.gee 10/01    |
 *----------------------------------------------------------------------*/
void calrhs_ele_neumann(FIELD        *actfield,    /* the active field */
                        SOLVAR       *actsolv,     
                        PARTITION    *actpart,
                        INTRA        *actintra,
                        SPARSE_TYP   *sysarraytyp,
                        SPARSE_ARRAY *sysarray,
                        int           actsysarray,
                        DIST_VECTOR  *rhs,
                        int           kstep,
                        CALC_ACTION  *action) /* action to be passsed to element routines */
{
ARRAY                 drhs_send_a;
double               *drhs_send;
ARRAY                 drhs_recv_a;
double               *drhs_recv;

#ifdef DEBUG 
dstrc_enter("calrhs_ele_neumann");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------ allocate temporary vectors */
/*------------------------------ these vectors are of global dimensions */
drhs_send = amdef("tmp",&drhs_send_a,rhs->numeq_total,1,"DV");
            amzero(&drhs_send_a);
#ifdef PARALLEL 
drhs_recv = amdef("tmp",&drhs_recv_a,rhs->numeq_total,1,"DV");
            amzero(&drhs_recv_a);
#endif
/*------------------------------------call element integration routines */
calelm(
         actfield,
         actsolv,
         actpart,
         actintra,
         actsysarray,
         -1,
         drhs_send,
         rhs->numeq_total,
         kstep,
         action
       );
/*-----------------------------------------------Allreduce vectors drhs */
#ifdef PARALLEL 
MPI_Allreduce(drhs_send,
              drhs_recv,
              rhs->numeq_total,
              MPI_DOUBLE,
              MPI_SUM,
              actintra->MPI_INTRA_COMM);
/*--------------------------------------add data to the DIST_VECTOR rhs */
assemble_vec(actintra,
                sysarraytyp,
                sysarray,
                rhs,
                drhs_recv,
                1.0);
#else
assemble_vec(actintra,
                sysarraytyp,
                sysarray,
                rhs,
                drhs_send,
                1.0);
#endif
/*------------------------------------------------------ delete buffers */
amdel(&drhs_send_a);
#ifdef PARALLEL 
amdel(&drhs_recv_a);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calrhs_ele_neumann */
