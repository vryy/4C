/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );
/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/
/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLPRECOND mlprecond;



/*!---------------------------------------------------------------------
\brief number the dofs of a coarse level                                             

<pre>                                                        m.gee 9/02 
number the dofs of a coarse level, but check,
whether this has been done before
</pre>
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\param numdf      int          (i)   number of dofs per supernode
\param actlev     MLLEVEL*     (i/o) the active level of the multilevel hierarchy
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_aggsetdofs(MLLEVEL *actlev,int numdf, INTRA *actintra)
{
int           i,j;
int           myrank,nproc;
int           sendbuff[MAXPROC],recvbuff[MAXPROC];
int           firstdof=0;
int           foundit =0;
int           counter;
int         **blocks;
AGG          *actagg;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_aggsetdofs");
#endif
/*----------------------------------------------------------------------*/
myrank  = actintra->intra_rank;
nproc   = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
exit:;
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_aggsetdofs */





/*!---------------------------------------------------------------------
\brief make the aggregates on the finest grid                                              

<pre>                                                        m.gee 9/02 
IMPORTANT:
This routine makes disconnected patches of nodes (blocks) in parallel.
It uses the information from the DBCSR matrix only, expecially the values
of dbcsr.firstcoupledof and dbcsr.blocks.
The goal is to create aggregates, that have a user chosen number of 
dofs and will later operate as a virtual coarse grid to the multilevel
preconditioner and will form a coarse grid stiffness matrix in DBCSR format. 
The aggregates and there dof numbering have to be chosen 
in a way, that local aggregates have lower dof numbers and aggregates
on domain boundaries have the high dof numbers. This is very important to
fullfill the requirements of the DBCSR matrix, which expects no interproc-
coupled equation below the dof number dbcsr.firstcoupledof on the coarse
level as well. This is also important to the prolongator smoother, which
will have to take extra care of interproc smoothing.

Aggregation is therefore done the following way:
1.) All blocks are split into purely local and interproc blocks.
2.) All purely local blocks are aggregated 
3.) All interproc blocks are aggregated
4.) coarse grid dof numbers are assigned starting with the low aggregates
    and taking care of the parallelity
    
</pre>
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\param actlev     MLLEVEL*     (i/o) the active level of the multilevel hierarchy
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_agg(MLLEVEL *actlev,INTRA *actintra)
{
int           i,j,k,counter;
int           myrank,nproc;
PARTDISCRET  *actpdis;
DISCRET      *actdis;
DBCSR        *A;
AGG          *agg;
AGG          *actagg;
AGG          *neighagg;
int           aggcounter;
int         **blocks;
int           nblock;
int           numdf;
int           dof;
int           fcd;
int          *ia,*ja,*update;
int           numeq,numeq_total;

int           nlblock,niblock,nb;
int         **lblock,**iblock;
int           icounter=0, lcounter=0;

int           nfreeblock;
int         **freeblock;

int          *actblock,*neighblock;
int          *bpatch[100];
int           nbpatch;
int           max;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_agg");
#endif
/*----------------------------------------------------------------------*/
myrank  = actintra->intra_rank;
nproc   = actintra->intra_nprocs;
A       = actlev->csr;
actpdis = mlprecond.partdis;
actdis  = mlprecond.fielddis;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_agg */
/*!---------------------------------------------------------------------
\brief make a neighbourhood to a patch                                              

<pre>                                                        m.gee 9/02 
make a neighbourhood to a patch. Find a neighbour to this neighbourhood 
and return it. The returned neighbour-neighbour has a graph-distance of 2
to the patch and will serve as new starting point for an aggregate
</pre>
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_getneightoagg(int **neighblock, 
                                 int  *bpatch[],
                                 int   nbpatch,
                                 int **freeblock,
                                 int   nfreeblock,
                                 int   numeq,
                                 int  *update,
                                 int  *ia,
                                 int  *ja)
{
int        i,j,k,l,m,n;
int        dof;
int        index;
int        column;
int        colstart,colend;
int       *actblock;
int        counter=0;
int       *neighpatch[200];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_getneightoagg");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
exit:;
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_getneightoagg */
/*!---------------------------------------------------------------------
\brief make a patch of nodal blocks from a list                                              

<pre>                                                        m.gee 9/02 

</pre>
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_getfreenblocks(int  *actblock,
                                  int **freeblock,
                                  int   nfreeblock,
                                  int  *bpatch[],
                                  int  *nbpatch,
                                  int   numeq,
                                  int  *update,
                                  int  *ia,
                                  int  *ja)
{
int        i,j,k,l,m,n;
int        dof;
int        index;
int        column;
int        colstart,colend;
int        foundneigh;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_getfreenblocks");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_getfreenblocks */



/*! @} (documentation module close)*/
