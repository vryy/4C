#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
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
\param numdf      INT          (i)   number of dofs per supernode
\param actlev     MLLEVEL*     (i/o) the active level of the multilevel hierarchy
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_aggsetdofs(MLLEVEL *actlev,INT numdf, INTRA *actintra)
{
INT           i,j;
INT           myrank,nproc;
INT           sendbuff[MAXPROC],recvbuff[MAXPROC];
INT           firstdof=0;
INT           foundit =0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_aggsetdofs");
#endif
/*----------------------------------------------------------------------*/
myrank  = actintra->intra_rank;
nproc   = actintra->intra_nprocs;
/*------------------------ check whether dofs have been numbered before */
for (i=0; i<actlev->nagg; i++)
{
   if (actlev->agg[i].numdf != 0) goto exit;
   /*---------------------------------- set number of dofs in aggregate */
   actlev->agg[i].numdf = numdf;
   /*--------------------------------- allocate vector for dofs numbers */
   actlev->agg[i].dof = (INT*)CCAMALLOC(numdf*sizeof(INT));
}
/*----------------------------- make total number of procs on each proc */
/* This is very easy, 'cause there are no supported dofs on higher levels */
for (i=0; i<nproc; i++) 
{
   if (i==myrank)
      sendbuff[i] = numdf * actlev->nagg;
   else
      sendbuff[i] = 0;
}
#ifdef PARALLEL
MPI_Allreduce(sendbuff,recvbuff,nproc,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
/*------------------------------------- count where my dofnumbers start */
for (i=0; i<myrank; i++)
{
   firstdof += recvbuff[i];
}
/*------------------------------- put the dof numbers to the aggregates */
/*-------------------- check for first appreance of interproc aggregate */
for (i=0; i<actlev->nagg; i++)
{
   for (j=0; j<actlev->agg[i].numdf; j++)
   {
      actlev->agg[i].dof[j] = firstdof;
      firstdof++;
   }
}
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
INT           i,j,k,counter;
INT           myrank,nproc;
PARTDISCRET  *actpdis;
DISCRET      *actdis;
DBCSR        *A;
AGG          *agg;
AGG          *actagg;
AGG          *neighagg;
INT           aggcounter;
INT         **blocks;
INT           nblock;
INT           numdf;
INT           dof;
INT           fcd;
INT          *ia,*ja,*update;
INT           numeq,numeq_total;

INT           nlblock,niblock,nb;
INT         **lblock,**iblock;
INT           icounter=0, lcounter=0;

INT           nfreeblock;
INT         **freeblock;

INT          *actblock,*neighblock;
INT          *bpatch[100];
INT           nbpatch;
INT           max;
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
blocks      = A->blocks.a.ia;
nblock      = A->blocks.fdim;
numdf       = A->blocks.sdim-1;
ia          = A->ia.a.iv;
ja          = A->ja.a.iv;
update      = A->update.a.iv;
fcd         = A->firstcoupledof;
numeq       = A->numeq;
numeq_total = A->numeq_total;
/*--------- calculate an approx. number of aggregates and allocate them */
if (actlev->nagg==0)
{
   actlev->nagg = (INT)((DOUBLE)nblock/4.0+(DOUBLE)nblock/8.0);
   if (actlev->nagg<1) actlev->nagg = 1;
   actlev->agg  = (AGG*)CCACALLOC(actlev->nagg,sizeof(AGG));
}
/*- start with the equations, which are NOT coupled to other processors */
/*- throw interproc blocks out of the list, make separate list for them */
/*---------------- a guess for number of local blocks nlblock is nblock */
nlblock = nblock;
lblock = (INT**)CCACALLOC(nlblock,sizeof(INT*));
/*----------------------------------- a guess for niblock is nblock/2.0 */
niblock = (INT)(nblock/2);
iblock = (INT**)CCACALLOC(niblock,sizeof(INT*));
/*---------------------- sort the blocks to the iblock or lblock vector */
for (i=0; i<nblock; i++)
{
   /* get first dof in the block */
   dof = blocks[i][1];
   /* if there is no unsupported dof in the block, it's an internal block */
   if (blocks[i][0]==0)
   {
      lblock[lcounter] = blocks[i];
      lcounter++;
      continue;
   }
   /* check whether the dof is lower or higherequal to firstcoupledof */
   /* it's an internal block */
   if (dof < fcd)
   {
      lblock[lcounter] = blocks[i];
      lcounter++;
   }
   /* it's an interproc dof */
   else
   {
      iblock[icounter] = blocks[i];
      icounter++;
      /* check size of ptr-vector iblock */
      if (icounter==niblock)
      {
         niblock += (INT)(nblock/2);
         iblock = (INT**)CCAREALLOC(iblock,niblock*sizeof(INT*));
      }
   }
}
if (icounter+lcounter != nblock)
   dserror("Sum of internal and interproc blocks wrong in aggregation");
niblock = icounter;
nlblock = lcounter;
/*------------------------------- we now make the internal blocks first */
aggcounter  = 0;
nfreeblock  = nlblock;
nb          = nlblock;
actblock    = NULL;
neighblock  = NULL;
max         = IMAX(nlblock,niblock);
freeblock   = (INT**)CCACALLOC(max,sizeof(INT*));
for (i=0; i<nb; i++) freeblock[i] = lblock[i];
/*------------------------ loop until all blocks are part of aggregate  */
while (nfreeblock != 0)
{
   /* find a new centernode for an aggregate */
   actblock=NULL;
   /* try to use a graph-distance two neighbour (if exists) to another aggregate */
   if (neighblock==NULL)
   {
      for (i=0; i<nb; i++)
      {
         if (freeblock[i]!=NULL)
         {
            actblock = freeblock[i];
            break;
         }
      }
   }
   else 
   actblock = neighblock;
   if (actblock==NULL && nfreeblock != 0) 
   dserror("Severe error in aggregation");
   /* make aggregate */
   /* get neighbours of this block */
   mlpcg_precond_getfreenblocks(actblock,freeblock,nb,
                                bpatch,&nbpatch,numeq,update,ia,ja);
   /* limit the size of the patch to 4 */
   nbpatch = IMIN(nbpatch,4);
   if (nbpatch==0)
      dserror("Patch of size 0 appeared in aggregation");
   /* add this patch to an aggregate */
   actlev->agg[aggcounter].coupling = mlpcg_agglocal;
   actlev->agg[aggcounter].nblock   = nbpatch;
   actlev->agg[aggcounter].block    = (INT**)CCACALLOC(nbpatch,sizeof(INT*));
   for (i=0; i<nbpatch; i++) 
   actlev->agg[aggcounter].block[i] = bpatch[i];
   aggcounter++;
   if (aggcounter==actlev->nagg)
   {
      actlev->nagg *= 2;
      actlev->agg = (AGG*)CCAREALLOC(actlev->agg,actlev->nagg*sizeof(AGG));
      if (!actlev->agg) dserror("Reallocation of Aggregates failed");
   }
   /* delete this patch from the list of freeblocks */
   for (j=0; j<nbpatch; j++)
   {
      for (i=0; i<nb; i++)
      {
         if (freeblock[i]==NULL) continue;
         if (bpatch[j]==freeblock[i])
         {
            freeblock[i]=NULL;
            break;
         }
      }
   }
   nfreeblock -= nbpatch;
   /* find one free neighbour to the patch */
   if (nfreeblock>0)
   mlpcg_precond_getneightoagg(&neighblock,bpatch,nbpatch,freeblock,nb,
                               numeq,update,ia,ja);
} 
/*---------------------------------- we now make the interproc blocks  */
nfreeblock  = niblock;
nb          = niblock;
actblock    = NULL;
neighblock  = NULL;
for (i=0; i<nb; i++) freeblock[i] = iblock[i];
/*------------------------ loop until all blocks are part of aggregate  */
while (nfreeblock != 0)
{
   /* find a new centernode for an aggregate */
   actblock=NULL;
   /* try to use a direct neighbour (if exists) to another aggregate */
   if (neighblock==NULL)
   {
      for (i=0; i<nb; i++)
      {
         if (freeblock[i]!=NULL)
         {
            actblock = freeblock[i];
            break;
         }
      }
   }
   else 
   actblock = neighblock;
   if (actblock==NULL && nfreeblock != 0) 
   dserror("Severe error in aggregation");
   /* make aggregate */
   /* get neighbours of this block */
   mlpcg_precond_getfreenblocks(actblock,freeblock,nb,
                                bpatch,&nbpatch,numeq,update,ia,ja);
   /* limit the size of the patch to 4 */
   nbpatch = IMIN(nbpatch,4);
   /* add this patch to an aggregate */
   actlev->agg[aggcounter].coupling = mlpcg_agginterproc;
   actlev->agg[aggcounter].nblock = nbpatch;
   actlev->agg[aggcounter].block  = (INT**)CCACALLOC(nbpatch,sizeof(INT*));
   for (i=0; i<nbpatch; i++) 
   actlev->agg[aggcounter].block[i] = bpatch[i];
   aggcounter++;
   if (aggcounter==actlev->nagg && nfreeblock-nbpatch != 0)
   {
      actlev->nagg *= 2;
      actlev->agg = (AGG*)CCAREALLOC(actlev->agg,actlev->nagg*sizeof(AGG));
      if (!actlev->agg) dserror("Reallocation of Aggregates failed");
   }
   /* delete this patch from the list of freeblocks */
   for (j=0; j<nbpatch; j++)
   {
      for (i=0; i<nb; i++)
      {
         if (freeblock[i]==NULL) continue;
         if (bpatch[j]==freeblock[i])
         {
            freeblock[i]=NULL;
            break;
         }
      }
   }
   nfreeblock -= nbpatch;
   /* find one free neighbour to the patch */
   if (nfreeblock>0)
   mlpcg_precond_getneightoagg(&neighblock,bpatch,nbpatch,freeblock,nb,
                               numeq,update,ia,ja);
} 
/* do not allow aggregates with a single block, merge with other aggregate */
for (i=0; i<aggcounter; i++)
{
   actagg = &(actlev->agg[i]);
   if (actagg->nblock>1) continue;
   /*--------------------------- found an aggregate with only one block */
   /*- find the direct graph-distance 1 neighbourhood of this aggregate */
   /* the has to be of same typ (mlpcg_agglocal or mlpcg_agginterproc) */
   if (actagg->coupling == mlpcg_agglocal)
   {
       dsassert(nlblock>1,"Only one local block in aggregation");
       nfreeblock  = nlblock;
       nb          = nlblock;
      for (k=0; k<nb; k++) freeblock[k] = lblock[k];
   }
   else if (actagg->coupling == mlpcg_agginterproc)
   {
       dsassert(niblock>1,"Only one inter block in aggregation");
       nfreeblock  = niblock;
       nb          = niblock;
       for (k=0; k<nb; k++) freeblock[k] = iblock[k];
   }
   else dserror("Unknown coupling type in aggregate");
   actblock = actagg->block[0];
   mlpcg_precond_getfreenblocks(actblock,freeblock,nb,
                                bpatch,&nbpatch,numeq,update,ia,ja);
   /*-- this patch includes myself, choose a neighbour, but not myself  */
   neighblock=NULL;
   for (j=0; j<nbpatch; j++)
      if (bpatch[j]!=actblock)
      {
         neighblock = bpatch[j];
         break;
      } 
   if (neighblock==NULL)
      continue;
   /*------------------- find the aggregate that is owner of this block */   
   neighagg = NULL;
   for (j=0; j<aggcounter; j++)
   {
      if (actagg == &(actlev->agg[j])) continue;
      for (k=0; k<actlev->agg[j].nblock; k++)
         if (actlev->agg[j].block[k] == neighblock)
         {
            neighagg = &(actlev->agg[j]);
            goto foundmergeagg;
         }
   }
   foundmergeagg:;
   if (neighagg==NULL) dserror("Cannot find aggregate to merge with");
   /*--------------------------------------------- enlarge the neighagg */
   neighagg->nblock++; /*----------- there will be one more block in it */
   neighagg->block = (INT**)CCAREALLOC(neighagg->block,neighagg->nblock*sizeof(INT*));
   /*------------------- put the lonely block to the neighbouraggregate */
   neighagg->block[neighagg->nblock-1] = actagg->block[0];
   /*------------------------------- delete the lonely aggregate actagg */
   actagg->nblock   = 0;
   actagg->coupling = mlpcg_aggnone;
   CCAFREE(actagg->block);
}/* end of for (i=0; i<aggcounter; i++) */
/*------------------------ realloc the aggregates to the correct amount */
/*---------- allocate a new vector of aggregates to a temporary pointer */
agg = (AGG*)CCACALLOC(aggcounter,sizeof(AGG));
/*-------------------- loop the aggregates and create them again in agg */
counter=0;
for (i=0; i<aggcounter; i++)
{
   actagg = &(actlev->agg[i]);
   if (actagg->coupling != mlpcg_aggnone)
   {
      agg[counter].coupling = actagg->coupling;
      agg[counter].nblock   = actagg->nblock;
      agg[counter].block = (INT**)CCACALLOC(agg[counter].nblock,sizeof(INT*));
      for (j=0; j<agg[counter].nblock; j++)
         agg[counter].block[j] = actagg->block[j];
      actagg->block = CCAFREE(actagg->block);
      counter++;
   }
}
/*--------------------------------------- get rid of the old agg vector */
CCAFREE(actlev->agg);
/*----------------------------- set the pointer to the actlev structure */
actlev->agg = agg;
/*---------------------------------------------- resize it if necessary */
if (actlev->nagg != counter)
{
   actlev->nagg  = counter;
   actlev->agg = (AGG*)CCAREALLOC(actlev->agg,actlev->nagg*sizeof(AGG));
}
/*----------------------------------------------------------------------*/
CCAFREE(freeblock);
CCAFREE(lblock);
CCAFREE(iblock);
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
void mlpcg_precond_getneightoagg(INT **neighblock, 
                                 INT  *bpatch[],
                                 INT   nbpatch,
                                 INT **freeblock,
                                 INT   nfreeblock,
                                 INT   numeq,
                                 INT  *update,
                                 INT  *ia,
                                 INT  *ja)
{
INT        i,j,k,l,m;
INT        dof;
INT        index;
INT        column;
INT        colstart,colend;
INT       *actblock;
INT        counter=0;
INT       *neighpatch[200];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_getneightoagg");
#endif
/*----------------------------------------------------------------------*/
*neighblock=NULL;
/*----------------------------------------------------------------------*/
/* loop the patch and try to find a free neighbour block which is not   */
/* guaranteed to exist                                                  */
/*----------------------------------------------------------------------*/
for (i=0; i<nbpatch; i++)
{
   actblock = bpatch[i];
   /* check size of block */
   if (actblock[0]==0) continue;
   dof   = actblock[1];
   index = mlpcg_getindex(dof,update,numeq);
   if (index==-1) dserror("Cannot find local dof in update");
   /* loop row of dof */
   colstart = ia[index];
   colend   = ia[index+1];
   for (j=colstart; j<colend; j++)
   {
      column = ja[j];
      /* check whether it belongs to a freeblock */
      /* if not, it is already member of this patch or */
      /* some other aggregate */
      for (k=0; k<nfreeblock; k++)
      {
         if (freeblock[k]==NULL) continue;
         for (l=0; l<freeblock[k][0]; l++)
         {
            if (column==freeblock[k][l+1])
            {
                /* check whether this freeblock[k] already exists in neighpatch */
                for (m=0; m<counter; m++)
                {
                   if (neighpatch[m]==freeblock[k])
                   goto nextj;
                }
                /* it's a new neighbour, put to list */
                neighpatch[counter]=freeblock[k];
                counter++;
                /* check size of neighpatch */
                if (counter==200) dserror("neighpatch too small in mlpcg_precond_getneightoagg");
                goto nextj;
            }
         }
      }      
      nextj:;
   }
}
/*----------------------------------------------------------------------*/
/* 
now loop the neighbourhood and find first block which is
1.) part of freeblock (therefore doesn not belong to the patch 
2.) not part of the neighbourhood
3.) neighbour to the neighbourhood
*/
/*----------------------------------------------------------------------*/
for (i=0; i<counter; i++)
{
   actblock = neighpatch[i];
   dof = actblock[1];
   index = mlpcg_getindex(dof,update,numeq);
   if (index==-1) dserror("Cannot find local dof in update");
   /* loop row of dof */
   colstart = ia[index];
   colend   = ia[index+1];
   for (j=colstart; j<colend; j++)
   {
      column = ja[j];
      /* check whether it belong to the neighbourpatch */
      for (k=0; k<counter; k++)
      {
          for (l=0; l<neighpatch[k][0]; l++)
          {
             if (column==neighpatch[k][l+1])
             goto nextcolumn;
          }
      }
      /* check whether it belong to freeblock */
      for (k=0; k<nfreeblock; k++)
      {
         if (freeblock[k]==NULL) continue;
         for (l=0; l<freeblock[k][0]; l++)
         {
            if (column==freeblock[k][l+1])
            {
               (*neighblock)=freeblock[k];
               goto nexttest;
            }
            
         }
      }
      nextcolumn:;
   }
}
nexttest:;
/*----------------------------------------------------------------------*/
/* 
check whether neiblock is zero or not. if no distance-2 neighbour could be
found, then find a distance-1 neighbour. If this also is not possible, then
return neighblock==NULL
*/
/*----------------------------------------------------------------------*/
if ((*neighblock)==NULL)
{
   for (i=0; i<nbpatch; i++)
   {
      actblock = bpatch[i];
      /* check size of block */
      if (actblock[0]==0) continue;
      dof   = actblock[1];
      index = mlpcg_getindex(dof,update,numeq);
      if (index==-1) dserror("Cannot find local dof in update");
      /* loop row of dof */
      colstart = ia[index];
      colend   = ia[index+1];
      for (j=colstart; j<colend; j++)
      {
         column = ja[j];
         /* check whether it belongs to a freeblock */
         /* if not, it is already member of this patch or */
         /* some other aggregate */
         for (k=0; k<nfreeblock; k++)
         {
            if (freeblock[k]==NULL) continue;
            for (l=0; l<freeblock[k][0]; l++)
            {
               if (column==freeblock[k][l+1])
               {
                   (*neighblock)=freeblock[k];
                   goto exit;
               }
            }
         }      
      
      }
   }
}
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
void mlpcg_precond_getfreenblocks(INT  *actblock,
                                  INT **freeblock,
                                  INT   nfreeblock,
                                  INT  *bpatch[],
                                  INT  *nbpatch,
                                  INT   numeq,
                                  INT  *update,
                                  INT  *ia,
                                  INT  *ja)
{
INT        i,j,k,l,m;
INT        dof;
INT        index;
INT        column;
INT        colstart,colend;
INT        foundneigh;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_getfreenblocks");
#endif
/*----------------------------------------------------------------------*/
*nbpatch=0;
/*------------------------------------- get a dof from the active block */
dof = actblock[1];
/*------------------------------------------ get the indize of this dof */
index = mlpcg_getindex(dof,update,numeq);
if (index==-1) dserror("Cannot find local dof on proc");
/*--------------------------------- get start and end index of this row */
colstart = ia[index];
colend   = ia[index+1];
if (colstart==colend)
   dserror("Found empty row in aggregation");
/*---------------------------------------- loop the columns of this row */
for (j=colstart; j<colend; j++)
{
   column     = ja[j];
   foundneigh = 0;
   /* check whether column is on a block in the freeblock list */
   for (i=0; i<nfreeblock; i++)
   {
      if (freeblock[i]==NULL) continue; /* this entry has already been deleted */
      for (k=0; k<freeblock[i][0]; k++)
      {
         if (column==freeblock[i][k+1]) /* found a block with this column */
         {
            foundneigh = 1;
            /* check whether the block with this column has been found before */
            for (l=0; l<(*nbpatch); l++)
            {
               for (m=0; m<bpatch[l][0]; m++)
               {
                  if (column==bpatch[l][m+1]) /* yes, has been found before */
                  {
                     foundneigh = 0;
                     goto nextj;
                  } 
               }
            }
            if (foundneigh==1)
            {
               bpatch[(*nbpatch)] = freeblock[i];
               (*nbpatch)++;
               /* limit the size to 4 */
               if (*nbpatch == 4) 
                  goto exit;
               if ((*nbpatch) == 100) dserror("Fixed size Block patch too small");
            }
            goto nextj;
         }
      }
   }
   nextj:;
}
exit:;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_getfreenblocks */









#if 0
INT mlpcg_precond_mergeagg(MLLEVEL *actlev, INTRA *actintra)
{
INT           i,j,k,l,n,m,max1,max2,max,counter;
INT           myrank,nproc;
DBCSR        *A;
AGG          *agg,*iagg,*actagg,*aggptr;
INT           nagg;
INT           niaggs[MAXPROC],niagg[MAXPROC];
INT          *update,*ia,*ja,numeq,numeq_total;
INT           owner[MAXPROC][2],owners[MAXPROC][2];
ARRAY4D       blocks_a,block_a;
INT       ****blocks,****block;
ARRAY         bonaggs_a,bonagg_a;
INT         **bonaggs,**bonagg;
INT           nblock1,**bptr1;
INT           nblock2,**bptr2;
INT           needagg,nmatch;
INT           rdofs[1000];
INT           nrdofs;
INT           dof1;
INT           ncouple[MAXPROC];
INT           index,own,colstart,colend;
INT           aggcounter;
#ifdef PARALLEL
MPI_Status    status;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_mergeagg");
#endif
#ifdef PARALLEL
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
A           = actlev->csr;
numeq       = A->numeq;
numeq_total = A->numeq_total;
update      = A->update.a.iv;
ia          = A->ia.a.iv;
ja          = A->ja.a.iv;
nagg        = actlev->nagg;
agg         = actlev->agg;
/*------------------------------------------------------ make ownership */
for (i=0; i<nproc; i++)
{
   owners[i][0] = 0;
   owners[i][1] = 0;
}
owners[myrank][0] = update[0];
owners[myrank][1] = update[numeq-1];
MPI_Allreduce(owners,owner,MAXPROC*2,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*----------------------------------------- find my first interproc agg */
for (n=0; n<nproc; n++) niaggs[n] = 0;
for (i=0; i<nagg; i++)
{
   if (agg[i].coupling==mlpcg_agginterproc) 
   {
      iagg           = &(agg[i]);
      niaggs[myrank] = nagg-i;
      break;
   }
}
/*--------------------------------------------- allreduce number of iaggs */
MPI_Allreduce(niaggs,niagg,nproc,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*------------------------------------------ make maximum number of iaggs */
max1 = 0;
for (n=0; n<nproc; n++)
   if (max1 < niagg[n]) max1 = niagg[n];
/*------------- make an array holding number of blocks for each aggregate */   
bonaggs = amdef("tmp",&bonaggs_a,nproc,max1,"IA");
bonagg  = amdef("tmp",&bonagg_a ,nproc,max1,"IA");
amzero(&bonaggs_a);
/*--------------------------------------------- fill my nblock to bonaggs */
for (i=0; i<niagg[myrank]; i++)
   bonaggs[myrank][i] = iagg[i].nblock;
/*-------------------------------------------------- allreduce this array */   
MPI_Allreduce(&(bonaggs[0][0]),&(bonagg[0][0]),nproc*max1,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
amdel(&bonaggs_a);
/*-------------------------------- find the maximum number of blocks max2 */
max2 = 0;
for (n=0; n<nproc; n++)
{
   for (m=0; m<niagg[n]; m++)
   {
      if (max2 < bonagg[n][m]) max2 = bonagg[n][m];
   }
}
/*- one aggregate consists of max. max2 blocks, each block maximal 6 long */
blocks = am4def("tmp",&blocks_a,nproc,max1,max2,7,"I4");
block  = am4def("tmp",&block_a ,nproc,max1,max2,7,"I4");
am4zero(&blocks_a);
/*-------------------------------------------- fill my own iagg to blocks */
for (i=0; i<niagg[myrank]; i++)
{
   actagg = &(iagg[i]);
   for (j=0; j<actagg->nblock; j++)
   {
      for (k=0; k<actagg->block[j][0]+1; k++)
      {
         blocks[myrank][i][j][k] = actagg->block[j][k];
      }
   }
}
/*-------------------------------------------- allreduce this block array */
MPI_Allreduce(blocks[0][0][0],block[0][0][0],nproc*max1*max2*7,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
am4del(&blocks_a);
/*------ loop processors and aggregates and find aggregates to merge with */
for (n=0; n<nproc; n++)/* loop over processors */
{
   for (i=0; i<niagg[n]; i++)/* loop over aggregates of processor */
   {
      if (myrank==n)
         actagg = &(iagg[i]);
      /* number of blocks on this aggregate */
      nblock1 = bonagg[n][i];
      if (nblock1==0) 
         continue;
      /* set pointer to these blocks */
      bptr1 = block[n][i];
      /* construct a list of remote dofs the aggregate is connected to */
      /* count whith which other proc this is mostly coupled */
      if (myrank==n)
      {
         for (k=0; k<nproc; k++) ncouple[k] = 0;
         for (k=0; k<nblock1; k++)
         {
            for (l=0; l<bptr1[k][0]; l++)
            {
               dof1 = bptr1[k][l+1];
               index = mlpcg_getindex(dof1,update,numeq);
               if (index==-1) dserror("Cannot find local dof1");
               colstart = ia[index];
               colend   = ia[index+1];
               for (j=colstart; j<colend; j++)
               {
                  own = mlpcg_getowner(ja[j],owner,nproc);
                  if (own==myrank) continue;
                  ncouple[own]++;
               }
            }
         }
         /* make list of rdofs from the proc which has most couplings */
         max = 0;
         m   = myrank;
         for (k=0; k<nproc; k++)
         if (ncouple[k]>max)
         {
            max = ncouple[k];
            m   = k;
         }
         if (m==myrank) 
            goto found_none;
         /* make list of dofs on m coupled to this aggregate */
         nrdofs = 0;
         for (k=0; k<nblock1; k++)
         {
            for (l=0; l<bptr1[k][0]; l++)
            {
               dof1 = bptr1[k][l+1];
               index = mlpcg_getindex(dof1,update,numeq);
               if (index==-1) dserror("Cannot find local dof1");
               colstart = ia[index];
               colend   = ia[index+1];
               for (j=colstart; j<colend; j++)
               {
                  own = mlpcg_getowner(ja[j],owner,nproc);
                  if (own != m) continue;
                  rdofs[nrdofs] = ja[j];
                  nrdofs++;
                  if (nrdofs==1000) goto out1;
               }
            }
         }
         out1:;
         mg_sort(rdofs,nrdofs,NULL,NULL);
         /*------------ processor n chooses an aggregate from processor m */
         needagg   = 0;
         nmatch    = 0;
         for (j=0; j<niagg[m]; j++)
         {
            nblock2 = bonagg[m][j];
            if (nblock2==0)
               continue;
            bptr2   = block[m][j];
            counter = 0;
            for (k=0; k<nblock2; k++)
            {
               for (l=0; l<bptr2[k][0]; l++)
               {
                  dof1 = bptr2[k][l+1];
                  index = find_index(dof1,rdofs,nrdofs);
                  if (index != -1) 
                     counter++;
               }
            }
            if (counter > nmatch)
            {
                nmatch  = counter;
                needagg = j;
            }
         }/*j loop over aggregates processor m */
      } /* end of myrank==n */
      /*------------------------------------------------ broadcast this m */
      found_none:
      MPI_Bcast(&m,1,MPI_INT,n,actintra->MPI_INTRA_COMM);
      if (m==n) 
         continue;
      /*---------------------------------- broadcast which aggregate on m */
      MPI_Bcast(&needagg,1,MPI_INT,n,actintra->MPI_INTRA_COMM);
      /*------------------------------------- broadcast number of matches */
      MPI_Bcast(&nmatch,1,MPI_INT,n,actintra->MPI_INTRA_COMM);
      /*------------------------------- if there was no match, do nothing */
      if (myrank==0)
         printf("n %d m %d nmatch %d\n",n,m,nmatch);
      if (nmatch==0) 
         continue;
      /*-------------------------- m sets pointer to approbiate aggregate */
      if (myrank==m)
         actagg = &(iagg[needagg]);
      /*---------------------------------------------- enter only n and m */
      /*----------------------- the common aggregate is always taken by n */
      nblock2 = bonagg[m][needagg];
      bptr2   = block[m][needagg];
      if (myrank==n)/* n gets the blocks */
      {
         actagg->nrblock = nblock2;
         actagg->rblock  = (INT**)CCAMALLOC((actagg->nrblock)*sizeof(INT*));
         for (j=0; j<nblock2; j++)
         {
            actagg->rblock[j] = (INT*)CCACALLOC(bptr2[j][0]+1,sizeof(INT));
            for (k=0; k<bptr2[j][0]+1; k++)
            {
               actagg->rblock[j][k] = bptr2[j][k];
            }
         }
      }
      if (myrank==m)
      {
         actagg->block    = CCAFREE(actagg->block);
         actagg->nblock   = 0;
         actagg->coupling = mlpcg_aggnone;
      }
      /*-------------------- all delete these aggregates from their lists */
      bonagg[m][needagg] = 0;
      bonagg[n][i]       = 0;
   }/*i loop over aggregates of processor n */
}/*n loop over processors */
/*-------------------------------realloc the correct amount of aggregates */
/*------------------------------------------------------ count aggregates */
aggcounter=0;
for (i=0; i<nagg; i++)
{
   if (actlev->agg[i].coupling != mlpcg_aggnone)
      aggcounter++;
}
/*------------------------- do nothing, if number of aggregates unchanged */
if (aggcounter==nagg) 
   goto exit;
/*---------------- allocate a vector of aggregates to a temporary pointer */
aggptr = (AGG*)CCACALLOC(aggcounter,sizeof(AGG));
/*------------------- loop the aggregates and create them again in aggptr */
aggcounter=0;
for (i=0; i<nagg; i++)
{
   actagg = &(actlev->agg[i]);
   if (actagg->coupling == mlpcg_aggnone)
      continue;
   aggptr[aggcounter].coupling = actagg->coupling;
   aggptr[aggcounter].nblock   = actagg->nblock;
   aggptr[aggcounter].block    = (INT**)CCACALLOC(aggptr[aggcounter].nblock,sizeof(INT*));
   if (actagg->nrblock != 0)
   {
      aggptr[aggcounter].nrblock  = actagg->nrblock;
      aggptr[aggcounter].rblock   = (INT**)CCACALLOC(aggptr[aggcounter].nrblock,sizeof(INT*));
      for (j=0; j<aggptr[aggcounter].nrblock; j++)
         aggptr[aggcounter].rblock[j] = actagg->rblock[j];
      actagg->rblock = CCAFREE(actagg->rblock);
   }
   for (j=0; j<aggptr[aggcounter].nblock; j++)
   {
      aggptr[aggcounter].block[j] = actagg->block[j];
   }
   actagg->block = CCAFREE(actagg->block);
   aggcounter++;
}
/*-------------------- free the old agg vector and set pointer to new one */
actlev->agg  = CCAFREE(actlev->agg);
actlev->agg  = aggptr;
actlev->nagg = aggcounter;
exit:
/*--------------------------------------------------------------- tidy up */
amdel(&bonagg_a);
am4del(&block_a);
/*------------------------------------------------------------------------*/
#endif /* end of PARALLEL */
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_mergeagg */
#endif




/*! @} (documentation module close)*/
#endif
