/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |  calculate the mask of an rc_ptr matrix               m.gee 1/02     |
 *----------------------------------------------------------------------*/
void mask_skyline(FIELD         *actfield, 
                  PARTITION     *actpart, 
                  SOLVAR        *actsolv,
                  INTRA         *actintra, 
                  SKYMATRIX     *sky)
{
INT       i;
INT       numeq;
INT     **dof_connect;
ARRAY     red_dof_connect;

ELEMENT  *actele;

#ifdef DEBUG 
dstrc_enter("mask_skyline");
#endif
/*----------------------------------------------------------------------*/
/* remember some facts:
   PARTITION is different on every proc.
   AZ_ARRAY_MSR will be different on every proc
   FIELD is the same everywhere
   In this routine, the vectors update and bindx and val are determined
   in size and allocated, the contents of the vectors update and bindx 
   are calculated
*/
/*------------------------------------------- put total size of problem */
sky->numeq_total = actfield->dis[0].numeq;
/* count number of eqns on proc and build processor-global couplingdof 
                                                                 matrix */
mask_numeq(actfield,actpart,actsolv,actintra,&numeq,0);
sky->numeq = numeq;
/*---------------------------------------------- allocate vector update */
amdef("update",&(sky->update),numeq,1,"IV");
amzero(&(sky->update));
/*--------------------------------put dofs in update in ascending order */
skyline_update(actfield,actpart,actsolv,actintra,sky);
/*------------------------ count number of nonzero entries on partition 
                                    and calculate dof connectivity list */
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
dof_connect = (INT**)CCACALLOC(sky->numeq_total,sizeof(INT*));
if (!dof_connect) dserror("Allocation of dof_connect failed");
skyline_nnz_topology(actfield,actpart,actsolv,actintra,sky,dof_connect);
/*------------------------------------------------------ make nnz_total */
sky->nnz_total=sky->nnz;
/*------------------------------make dof_connect redundant on all procs */
skyline_make_red_dof_connect(actfield,
                             actpart,
                             actsolv,
                             actintra,
                             sky,
                             dof_connect,
                             &red_dof_connect);
/*---------------------------------------- make arrays from dof_connect */
skyline_make_sparsity(sky,&red_dof_connect);
/*---------------------------------------- delete the array dof_connect */
for (i=0; i<sky->numeq_total; i++)
{
   if (dof_connect[i]) CCAFREE(dof_connect[i]);
}
CCAFREE(dof_connect);
amdel(&red_dof_connect);

/* make the index vector for faster assembling */
  for (i=0; i<actpart->pdis[0].numele; i++)
  {
    actele = actpart->pdis[0].element[i];
    sky_make_index(actfield,actpart,actintra,actele,sky);
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_skyline */



/*----------------------------------------------------------------------*
 |  allocate update put dofs in update in ascending order   m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  skyline_update(FIELD         *actfield, 
                    PARTITION     *actpart, 
                    SOLVAR        *actsolv,
                    INTRA         *actintra,
                    SKYMATRIX     *sky)
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
#ifdef DEBUG 
dstrc_enter("skyline_update");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------ make a local copy of the array actpart->coupledofs */
am_alloc_copy(&(actpart->pdis[0].coupledofs),&coupledofs);
/*------------------------------------- loop the nodes on the partition */
update = sky->update.a.iv;
counter=0;
for (i=0; i<actpart->pdis[0].numnp; i++)
{
   actnode = actpart->pdis[0].node[i];
   for (l=0; l<actnode->numdf; l++)
   {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->dis[0].numeq) continue;
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
/*---------- check whether the correct number of dofs have been counted */
if (counter != sky->numeq) dserror("Number of dofs in update wrong");
/*---------------------------- sort the vector update just to make sure */
qsort((INT*) update, counter, sizeof(INT), cmp_int);
/*----------------------------------------------------------------------*/
amdel(&coupledofs);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_update */



/*----------------------------------------------------------------------*
 |  calculate number of nonzero entries and dof topology    m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  skyline_nnz_topology(FIELD      *actfield, 
                         PARTITION    *actpart, 
                         SOLVAR       *actsolv,
                         INTRA        *actintra,
                         SKYMATRIX    *sky,
                         INT         **dof_connect)
{
INT        i,j,k,l,m;
INT        counter,counter2;
INT        dof;
INT        nnz;
INT        iscoupled;
INT       *update;
INT        numeq;
INT        actdof;
NODE      *centernode;
NODE      *actnode;
ELEMENT   *actele;
ARRAY      dofpatch;
ARRAY     *coupledofs;
INT        imyrank;
INT        inprocs;
INT        actndis;

#ifdef PARALLEL 
MPI_Status status;
#endif

#ifdef DEBUG 
dstrc_enter("skyline_nnz_topology");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*------------------------------------------- set actual discretisation */
actndis=0;
/*----------------------------------------------------------- shortcuts */
sky->nnz=0;
/*numeq  = sky->numeq;*/
numeq  = sky->numeq_total;
update = sky->update.a.iv;
for (i=0; i<sky->numeq_total; i++) dof_connect[i]=NULL;
amdef("tmp",&dofpatch,1000,1,"IV");
amzero(&dofpatch);
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   dof = i;
   /*dof = update[i];*/
   /*------------------------------ check whether this is a coupled dof */
   iscoupled=0;
   dof_in_coupledofs(dof,actpart,&iscoupled);
   if (iscoupled==1) continue;
   /*--------------------------------- find the centernode for this dof */
   centernode=NULL;
   for (j=0; j<actfield->dis[actndis].numnp; j++)
   {
      for (k=0; k<actfield->dis[actndis].node[j].numdf; k++)
      {
         if (actfield->dis[actndis].node[j].dof[k] == dof)
         {
            centernode = &actfield->dis[actndis].node[j];
          goto nodefound1;
         }
      }
   }
   nodefound1:
   dsassert(centernode!=NULL,"Cannot find centernode of a patch");
   /*--------------------------------- make dof patch around centernode */
   counter=0;
   for (j=0; j<centernode->numele; j++)
   {
      actele = centernode->element[j];
      for (k=0; k<actele->numnp; k++)
      {
         actnode = actele->node[k];
         for (l=0; l<actnode->numdf; l++)
         {
            if (actnode->dof[l] < actfield->dis[0].numeq)
            {
               if (counter>=dofpatch.fdim) amredef(&dofpatch,dofpatch.fdim+500,1,"IV");
               dofpatch.a.iv[counter] = actnode->dof[l];
               counter++;
            }
         }
      }
   }
   /*----------------------------------------- delete doubles on patch */
   /*------------------------------- also delete dof itself from patch */
   for (j=0; j<counter; j++)
   {
      actdof = dofpatch.a.iv[j];
      if (dofpatch.a.iv[j]==dof) dofpatch.a.iv[j]=-1;
      if (actdof==-1) continue;
      for (k=j+1; k<counter; k++)
      {
         if (dofpatch.a.iv[k] == actdof ||
             dofpatch.a.iv[k] == dof      ) dofpatch.a.iv[k]=-1;
      }
   }
   /*----------------------------------- count number of dofs on patch */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) counter2++;
   }
   /*-------------- allocate the dof_connect vector and put dofs in it */
   dof_connect[dof] = (INT*)CCACALLOC(counter2+3,sizeof(INT));
   if (!dof_connect[dof]) dserror("Allocation of dof connect list failed");
   dof_connect[dof][0] = counter2+3;
   dof_connect[dof][1] = 0; 
   dof_connect[dof][2] = dof;
   /*
      dof_connect[i][0] = lenght of dof_connect[i]
      dof_connect[i][1] = iscoupled ( 1 or 2 ) done later on 
      dof_connect[i][2] = dof
      dof_connect[i][ 2..dof_connect[i][0]-1 ] = connected dofs exluding itself 
   */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) 
      {
         dof_connect[dof][counter2+3] = dofpatch.a.iv[j];
         counter2++;
      }
   }
}  /* end of loop over numeq */ 
/*--------------------------------------------- now do the coupled dofs */
coupledofs = &(actpart->pdis[0].coupledofs);
for (i=0; i<coupledofs->fdim; i++)
{
   dof = coupledofs->a.ia[i][0];
   /*------------------------------------- find all patches to this dof */
   counter=0;
   for (j=0; j<actfield->dis[actndis].numnp; j++)
   {
      centernode=NULL;
      /* */
      for (l=0; l<actfield->dis[actndis].node[j].numdf; l++)
      {
         if (dof == actfield->dis[actndis].node[j].dof[l])
         {
            centernode = &actfield->dis[actndis].node[j];
            break;
         }
      }
      if (centernode !=NULL)
      {
         /*--------------------------- make dof patch around centernode */
         for (k=0; k<centernode->numele; k++)
         {
            actele = centernode->element[k];
            for (m=0; m<actele->numnp; m++)
            {
               actnode = actele->node[m];
               for (l=0; l<actnode->numdf; l++)
               {
                  if (actnode->dof[l] < actfield->dis[0].numeq)
                  {
                     if (counter>=dofpatch.fdim) amredef(&dofpatch,dofpatch.fdim+500,1,"IV");
                     dofpatch.a.iv[counter] = actnode->dof[l];
                     counter++;
                  }
               }
            }
         }
      }
   }/* end of making dofpatch */
   /*----------------------------------------- delete doubles on patch */
   for (j=0; j<counter; j++)
   {
      actdof = dofpatch.a.iv[j];
      if (actdof==-1) continue;
      if (actdof==dof) dofpatch.a.iv[j]=-1;
      for (k=j+1; k<counter; k++)
      {
         if (dofpatch.a.iv[k] == actdof ||
             dofpatch.a.iv[k] == dof      ) dofpatch.a.iv[k]=-1;
      }
   }
   /*----------------------------------- count number of dofs on patch */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) counter2++;
   }
   /*-------------- allocate the dof_connect vector and put dofs in it */
   dof_connect[dof] = (INT*)CCACALLOC(counter2+3,sizeof(INT));
   if (!dof_connect[dof]) dserror("Allocation of dof connect list failed");
   dof_connect[dof][0] = counter2+3;
   dof_connect[dof][1] = 0;
   dof_connect[dof][2] = dof;
   /*-------------------------- put the patch to the dof_connect array */
   counter2=0;
   for (j=0; j<counter; j++)
   {
      if (dofpatch.a.iv[j] != -1) 
      {
         dof_connect[dof][counter2+3] = dofpatch.a.iv[j];
         counter2++;
      }
   }
} /* end of loop over coupled dofs */
/*---------------------------- now go through dof_connect and count nnz */
nnz=0;
for (i=0; i<numeq; i++)
{
   nnz += (dof_connect[i][0]-2);
}
sky->nnz=nnz;
/*--------- last thing to do is to order dof_connect in ascending order */
for (i=0; i<numeq; i++)
{
   qsort((INT*)(&(dof_connect[i][3])), dof_connect[i][0]-3, sizeof(INT), cmp_int);
}
/*----------------------------------------------------------------------*/
amdel(&dofpatch);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_nnz_topology */


/*----------------------------------------------------------------------*
 |  make the dof_connect list redundant                    m.gee 01/02  |
 *----------------------------------------------------------------------*/
void   skyline_make_red_dof_connect(FIELD         *actfield, 
                                   PARTITION     *actpart, 
                                   SOLVAR        *actsolv,
                                   INTRA         *actintra,
                                   SKYMATRIX     *sky,
                                   INT          **dof_connect,
                                   ARRAY         *red_dof_connect)
{
INT        i,j;

INT      **reddof;

INT        max_dof_connect;
/*----------------------------------------------------------------------*/
/* communicate coupled dofs */
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("skyline_make_red_dof_connect");
#endif
/*----------------------------- check for largest row in my dof_connect */
max_dof_connect=0;
for (i=0; i<sky->numeq_total; i++)
{
   if (dof_connect[i])
   if (dof_connect[i][0]>max_dof_connect)
   max_dof_connect=dof_connect[i][0];
}
/*---------------- allocate temporary array to hold global connectivity */
reddof = amdef("tmp",red_dof_connect,sky->numeq_total,max_dof_connect,"IA");
/*----------------------------- put my own dof_connect values to reddof */
for (i=0; i<sky->numeq_total; i++)
{
   if (dof_connect[i])
   {
      for (j=0; j<dof_connect[i][0]; j++) 
         reddof[i][j] = dof_connect[i][j];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_make_red_dof_connect */







/*----------------------------------------------------------------------*
 |  make sparsity mask for skyline matrix                   m.gee 1/02  |
 *----------------------------------------------------------------------*/
void  skyline_make_sparsity(SKYMATRIX  *sky, ARRAY *red_dof_connect)
{
INT        i,j;
INT      **reddof;
INT       *maxa;
INT        actdof;
INT        lenght;
INT        counter=0;
INT        mindof;
#ifdef DEBUG 
dstrc_enter("skyline_make_sparsity");
#endif
/*----------------------------------------------------------------------*/
   /*
      reddof[i][0] = lenght of reddof[i]
      reddof[i][1] = iscoupled ( 1 or 2 ) done later on 
      reddof[i][2] = dof
      reddof[i][ 2..reddof[i][0]-1 ] = connected dofs exluding itself 
   */
reddof = red_dof_connect->a.ia;
/*------------------------------------------------------- allocate maxa */
maxa = amdef("maxa",&(sky->maxa),sky->numeq_total+1,1,"IV");
/*------------------------------------------------------- loop the dofs */       
for (i=0; i<sky->numeq_total; i++)
{
   actdof = reddof[i][2];
   if (actdof != i) printf("Warning: skyline format mixed up!\n");
   
   /*-------------- search reddof[i][2..reddof[i][0] ] for smallest dof */
   mindof=10000000;
   for (j=2; j<reddof[i][0]; j++)
   {
      if (mindof>reddof[i][j]) 
      mindof = reddof[i][j];
   }
   /*-------------------------- check distance between actdof and mindof */
   lenght = actdof-mindof+1;
   /*--------------------------- maxa[i] holds start of column of actdof */
   maxa[i]=counter;
   counter+=lenght;
}
maxa[i]=counter;
/*-----------------------------------------------------------allocate A */
amdef("A",&(sky->A),maxa[i],1,"DV");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of skyline_make_sparsity */




/*----------------------------------------------------------------------*/
/*!
 \brief 

 This routine determines the location vactor for the actele and stores it
 in the element structure.  Furthermore for each component [i][j] in the 
 element stiffness matrix the position in the 1d sparse matrix is 
 calculated and stored in actele->index[i][j]. These can be used later on 
 for the assembling procedure.
  
 \param actfield  *FIELD        (i)  the field we are working on
 \param actpart   *PARTITION    (i)  the partition we are working on
 \param actintra  *INTRA        (i)  the intra-communicator we do not need
 \param actele    *ELEMENT      (i)  the element we would like to work with
 \param sky1      *SKYMATRIX    (i)  the sparse matrix we will assemble into

 \author mn
 \date 07/04

 */
/*----------------------------------------------------------------------*/
void sky_make_index(
    FIELD                 *actfield, 
    PARTITION             *actpart,
    INTRA                 *actintra,
    ELEMENT               *actele,
    struct _SKYMATRIX     *sky1
    )
{

  INT               i,j,counter;           /* some counter variables */
  INT               ii,jj;                 /* counter variables for system matrix */
  INT               nd;                    /* size of estif */
  INT               numeq_total;           /* total number of equations */
  INT               numeq;                 /* number of equations on this proc */
  INT               myrank;                /* my intra-proc number */
  INT               nprocs;                /* my intra- number of processes */
  DOUBLE           *A;                      /* the skyline matrix 1 */
  INT              *maxa;

  INT               startindex;
  INT               height;
  INT               distance;
  INT               index;

  struct _ARRAY ele_index;
  struct _ARRAY ele_locm;
#ifdef PARALLEL
  struct _ARRAY ele_owner;
#endif


#ifdef DEBUG 
  dstrc_enter("sky_make_index");
#endif


  /* set some pointers and variables */
  myrank     = actintra->intra_rank;
  nprocs     = actintra->intra_nprocs;
  numeq_total= sky1->numeq_total;
  numeq      = sky1->numeq;
  A          = sky1->A.a.dv;
  maxa       = sky1->maxa.a.iv;


  /* determine the size of estiff */
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      counter++;
    }
  }
  /* end of loop over element nodes */
  nd = counter;
  actele->nd = counter;


  /* allocate locm, index and owner */
  actele->locm  = amdef("locm" ,&ele_locm ,nd, 1,"IV");
  actele->index = amdef("index",&ele_index,nd,nd,"IA");
#ifdef PARALLEL
  actele->owner = amdef("owner",&ele_owner,nd, 1,"IV");
#endif


  /* make location vector locm */
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      actele->locm[counter]    = actele->node[i]->dof[j];
#ifdef PARALLEL 
      actele->owner[counter]   = actele->node[i]->proc;
#endif
      counter++;
    }
  }
  /* end of loop over element nodes */




  /* loop over i (the element row) */
  for (i=0; i<nd; i++)
  {
    ii = actele->locm[i];

    /* check for boundary condition */
    if (ii>=numeq_total)
    {
      for (j=0; j<nd; j++) actele->index[i][j] = -1;
      continue;
    }

    /* check for ownership of row ii */
#ifdef PARALLEL 
    if (actele->owner[i]!=myrank)
    {
      for (j=0; j<nd; j++) actele->index[i][j] = -1;
      continue;
    }
#endif

    /* start of the skyline of ii is maxa[ii] */
    startindex = maxa[ii];

    /* height of the skyline of ii */
    height     = maxa[ii+1]-maxa[ii];

    /* loop over j (the element column) */
    for (j=0; j<nd; j++)
    {
      jj = actele->locm[j];

      /* check for boundary condition */
      if (jj>=numeq_total)
      {
        actele->index[i][j] = -1;
        continue;
      }

      /* find position [ii][jj] in A */
      distance  = ii-jj;
      if (distance < 0) 
      {
        actele->index[i][j] = -1;
        continue;
      }
      index     = startindex+distance;
      actele->index[i][j] = index;
    } /* end loop over j */
  }/* end loop over i */


#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of sky_make_index */


