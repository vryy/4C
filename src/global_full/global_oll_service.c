/*!----------------------------------------------------------------------
\file
\brief contains functions to handle oll matrices

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
/*! 
\addtogroup OLL 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief numeq for oll matrices 

<pre>                                                              mn 02/03
This function counts the number of equations on this processor
</pre>
\param *actfield       FIELD  (i)   the active field
\param *actpart        PARTITION  (i)   the active partition
\param *actINTra       INTRA  (i)   the active communicator
\param  dis            INT    (i)   number of the active discretization
\param *numeq          INT    (o)   number of equations

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void oll_numeq(
    FIELD         *actfield, 
    PARTITION    *actpart, 
    INTRA        *actintra,
    INT           dis,
    INT          *numeq)
{
  INT       i,j,k,l;
  INT       counter;
  INT       dof;
  INT       iscoupled;
  INT      *sendbuff,*recvbuff, sendsize;
  INT      *tmp;
  INT       inter_proc;
  long int  min;
  INT       proc;
  INT       inprocs;
  INT       imyrank;
  NODE     *actnode;

INT       no_coupling = 0;

#ifdef DEBUG 
  dstrc_enter("oll_numeq");
#endif
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------------- first make a list of dofs which are coupled */
  /*----------------------------------- estimate size of coupdofs to 5000 */
  amdef("coupledofs",&(actpart->pdis[dis].coupledofs),5000,1,"IV");
  amzero(&(actpart->pdis[dis].coupledofs));
  counter=0;
  /*-------------------------------- loop all nodes and find coupled dofs */
  for (i=0; i<actfield->dis[dis].numnp; i++)
  {
    actnode = &(actfield->dis[dis].node[i]);
    if (actnode->gnode->couple==NULL) continue;
    for (l=0; l<actnode->numdf; l++)
    {
      if (actnode->dof[l]>=actfield->dis[dis].numeq) continue;
      /* there is coupling on this dof */
      if (actnode->gnode->couple->couple.a.ia[l][0] != 0 ||
          actnode->gnode->couple->couple.a.ia[l][1] != 0 )
      {
        if (counter>=actpart->pdis[dis].coupledofs.fdim) 
          amredef(&(actpart->pdis[dis].coupledofs),(actpart->pdis[dis].coupledofs.fdim+5000),1,"IV");
        /* the coupled dof could be dirichlet conditioned */
        if (actnode->dof[l]<actfield->dis[dis].numeq)
        {
          actpart->pdis[dis].coupledofs.a.iv[counter] = actnode->dof[l];
          counter++;
        }
      }
    }
  }

if (counter ==0)
  no_coupling = 1;
else
  no_coupling = 0;

  amredef(&(actpart->pdis[dis].coupledofs),counter,1,"IV");
  /*---------------------------------- delete the doubles in coupledofs */

if (!no_coupling)
{

  for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)
  {
    if (actpart->pdis[dis].coupledofs.a.iv[i]==-1) continue;
    dof = actpart->pdis[dis].coupledofs.a.iv[i];
    for (j=i+1; j<actpart->pdis[dis].coupledofs.fdim; j++)
    {
      if (actpart->pdis[dis].coupledofs.a.iv[j]==dof) actpart->pdis[dis].coupledofs.a.iv[j]=-1;
    }
  }
  /*--------- move all remaining coupdofs to the front and redefine again */
  counter=0;
  for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)
  {
    if (actpart->pdis[dis].coupledofs.a.iv[i]!=-1)
    {
      actpart->pdis[dis].coupledofs.a.iv[counter] = actpart->pdis[dis].coupledofs.a.iv[i];
      counter++;
    }
  }
  amredef(&(actpart->pdis[dis].coupledofs),counter,inprocs+1,"IA");
  /*------------------- the newly allocated columns have to be initialized */
  for (i=1; i<actpart->pdis[dis].coupledofs.sdim; i++)
    for (j=0; j<actpart->pdis[dis].coupledofs.fdim; j++) 
      actpart->pdis[dis].coupledofs.a.ia[j][i]=0;

} /* end of if(!no_coupling) */

  /* processor looks on his own domain whether he has some of these coupdofs, 
     puts this information in the array coupledofs in the column myrank+1, so it 
     can be allreduced 

     The matrix has the following style (after allreduce on all procs the same):

     ----------------------
     | 12 | 1 | 0 | 1 | 0 |
     | 40 | 1 | 0 | 0 | 0 |
     | 41 | 1 | 1 | 1 | 1 |
     | 76 | 0 | 1 | 1 | 0 |
     ----------------------

     column 0                : number of the coupled equation
     column 1 - inprocs+1 : proc has coupled equation or not

*/

if (!no_coupling)
{

  if (inprocs==1) /*--------------------------------- sequentiell version */
  {
    for (k=0; k<actpart->pdis[dis].coupledofs.fdim; k++)
    {
      actpart->pdis[dis].coupledofs.a.ia[k][imyrank+1]=2;
    }
  }
  else /*----------------------------------------------- parallel version */
  {
    /*
       actpart->node[i] really loops only nodes with dofs updated on this proc
       */
    for (i=0; i<actpart->pdis[dis].numnp; i++) /* now loop only my nodes */
    {
      for (l=0; l<actpart->pdis[dis].node[i]->numdf; l++)
      {
        dof = actpart->pdis[dis].node[i]->dof[l];
        for (k=0; k<actpart->pdis[dis].coupledofs.fdim; k++)
        {
          if (actpart->pdis[dis].coupledofs.a.ia[k][0]==dof)
          {
            actpart->pdis[dis].coupledofs.a.ia[k][imyrank+1]=1;
            break;
          }
        }
      }
    }
  }
  /* ----- Allreduce the whole array, so every proc knows about where all 
     coupledofs are */
#ifdef PARALLEL
  sendsize = (actpart->pdis[dis].coupledofs.fdim)*(inprocs);
  sendbuff = (INT*)CCACALLOC(sendsize,sizeof(INT));
  recvbuff = (INT*)CCACALLOC(sendsize,sizeof(INT));
  if (sendbuff==NULL || recvbuff==NULL) dserror("Allocation of temporary memory failed");
  counter=0;
  for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)
  {
    for (j=0; j<inprocs; j++)
    {
      sendbuff[counter] = actpart->pdis[dis].coupledofs.a.ia[i][j+1];
      counter++;
    }
  }
  MPI_Allreduce(sendbuff,
      recvbuff,
      sendsize,
      MPI_INT,
      MPI_SUM,
      actintra->MPI_INTRA_COMM);
  counter=0;
  for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)
  {
    for (j=0; j<inprocs; j++)
    {
      actpart->pdis[dis].coupledofs.a.ia[i][j+1] = recvbuff[counter];
      counter++;
    }
  }
  CCAFREE(sendbuff);CCAFREE(recvbuff);
#endif

} /* end of if(!no_coupling) */

  /*------- count number of equations on partition including coupled dofs */
  /*---------------------------------------- count the coupled ones first */
  counter=0;
  for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)
  {
    if (actpart->pdis[dis].coupledofs.a.ia[i][imyrank+1]!=0) counter++;
  }
  /*-------------------------------- count all dofs which are not coupled */
  for (i=0; i<actpart->pdis[dis].numnp; i++)
  {
    actnode = actpart->pdis[dis].node[i];
    for (l=0; l<actnode->numdf; l++)
    {
      dof = actnode->dof[l];
      iscoupled=0;
      for (k=0; k<actpart->pdis[dis].coupledofs.fdim; k++)
      {
        if (dof == actpart->pdis[dis].coupledofs.a.ia[k][0]) 
        {
          iscoupled=1;
          break;
        }
      }
      if (iscoupled==0) 
      {
        if (dof < actfield->dis[dis].numeq)
          counter++;
      }
    }
  }
  /*--- number of equations on this partition including the coupled ones */
  *numeq = counter;
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

if (!no_coupling)
{

  if (inprocs > 1)
  {
    tmp = (INT*)CCACALLOC(inprocs,sizeof(INT));
    if (!tmp) dserror("Allocation of temporary memory failed");
    for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)/*  loop coupled eqns */
    {
      /*--------------------------------- check whether its inter-proc eqn */
      inter_proc=0;
      for (j=0; j<inprocs; j++) inter_proc += actpart->pdis[dis].coupledofs.a.ia[i][j+1];
      if (inter_proc==1)/*----------------- no inter-processor coupling */
      {
        for (j=0; j<inprocs; j++)
        {
          if (actpart->pdis[dis].coupledofs.a.ia[i][j+1]==1) 
          {
            actpart->pdis[dis].coupledofs.a.ia[i][j+1]=2;
            break;
          }
        }
      }
      else/*----------------------------- eqn is an inter-proc equation */
      {
        /* there won't be more than a million procs in the near future....*/
        min=1000000;
        proc=-1;
        for (j=0; j<inprocs; j++)
        {
          if (actpart->pdis[dis].coupledofs.a.ia[i][j+1]==1) 
          {
            if (tmp[j]<=min)
            {
              min = tmp[j];
              proc = j;
            }
          }
        }
        actpart->pdis[dis].coupledofs.a.ia[i][proc+1]=2;
        tmp[proc] += 1;
      }
    }/* end loop over coupling eqns */
    CCAFREE(tmp);
  }
  /* procs who have not become owner of a coupling equation have to reduce there
     number of equations */
  if (inprocs > 1)
  {
    for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)/* loop coupled eqns */
    {
      /* ------Yes, I am slave owner of an inter_proc coupling equation */
      if (actpart->pdis[dis].coupledofs.a.ia[i][imyrank+1]==1)
      {
        (*numeq) = (*numeq)-1;
      }
      /* master owner of equation do nothing, 'cause the equation has been
         counted already */
    }
  }

} /* end of if(!no_coupling) */

  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of oll_numeq */



/*!----------------------------------------------------------------------
\brief check whether dof is in coupledofs 

<pre>                                                              mn 02/03
This function checks whether this dof is in coupled dofs
</pre>
\param  dof            INT        (i)   the dof in question
\param *actpart        PARTITION  (i)   the active partition
\param *iscoupled      INT        (o)   the answer 
\param  dis            INT        (i)   actual disnumber

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void oll_dof_in_coupledofs(
    INT dof,
    PARTITION *actpart,
    INT *iscoupled,
    INT  dis)
{
  INT       i;
#ifdef DEBUG 
  dstrc_enter("oll_dof_in_coupledofs");
#endif
  /*----------------------------------------------------------------------*/
  for (i=0; i<actpart->pdis[dis].coupledofs.fdim; i++)
  {
    if (dof==actpart->pdis[dis].coupledofs.a.ia[i][0])
    {
      *iscoupled = 1;
      break;
    }
  }
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of oll_dof_in_coupledofs */




/*!----------------------------------------------------------------------
\brief nnz and topology for oll 

<pre>                                                              mn 02/03
This function calculates the number of non-zero entrys and the topology
for oll matrices. Only the nnz is needed.
</pre>
\param *actfield       FIELD  (i)   the active field
\param *actpart        PARTITION  (i)   the active partition
\param *actintra       INTRA  (i)   the active communicator
\param *oll            OLL    (i)   the oll matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void oll_nnz_topology(
    FIELD         *actfield, 
    PARTITION     *actpart, 
    INTRA         *actintra,
    OLL           *oll,
    INT            dis)
{
  INT        i,j,k,l,m;
  INT        counter,counter2;
  INT        dof;
  INT        nnz;
  INT        iscoupled;
  INT       *update;
  INT        numeq;
  INT        actdof;
  INT        dofflag;
  INT        dofmaster;
  INT        dofslave;
  INT        recvlength;
  NODE      *centernode;
  NODE      *actnode;
  ELEMENT   *actele;
  ARRAY      dofpatch;
  ARRAY     *coupledofs;
  INT        imyrank;
  INT        inprocs;

  ARRAY      dof_nnz;

  NODE     **node_dof;
  INT        max_dof,dof_id;


#ifdef PARALLEL 
  MPI_Status status;
#endif

#ifdef DEBUG 
  dstrc_enter("oll_nnz_topology");
#endif
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*----------------------------------------------------------- shortcuts */
  oll->nnz=0;
  numeq  = oll->numeq;
  update = oll->update.a.iv;
  amdef("tmp2",&dof_nnz,oll->numeq_total,1,"IV");
  amzero(&dof_nnz);
  amdef("tmp",&dofpatch,MAX_NNZPERROW,1,"IV");
  amzero(&dofpatch);

  /* create node pointers for all dofs in this partition */
  max_dof = 0;
  for (j=0; j<actpart->pdis[dis].numnp; j++)
  {
    for (k=0; k<actpart->pdis[dis].node[j]->numdf; k++)
    {
      if (actpart->pdis[dis].node[j]->dof[k] >= max_dof)
      {
        max_dof = actpart->pdis[dis].node[j]->dof[k];
      }
    }
  }
  /* allocate pointer vector to the nodes */
  node_dof = (NODE**)CCACALLOC(max_dof+1,sizeof(NODE*));
  /* store pointers to nodes in node_dof at position accord. to dof */
  for (j=0; j<actpart->pdis[dis].numnp; j++)
  {
    for (k=0; k<actpart->pdis[dis].node[j]->numdf; k++)
    {
      dof_id = actpart->pdis[dis].node[j]->dof[k];
      dsassert(dof_id <= max_dof,"zu kleiner node_dof Vector");
      node_dof[dof_id] = actpart->pdis[dis].node[j];
    }
  }

  /*----------------------------------------------------------------------*/
  for (i=0; i<numeq; i++)
  {
    dof = update[i];
    /*------------------------------ check whether this is a coupled dof */
    iscoupled=0;
    oll_dof_in_coupledofs(dof,actpart,&iscoupled,dis);
    if (iscoupled==1) continue;
    /*--------------------------------- find the centernode for this dof */
    centernode=NULL;
    centernode = node_dof[dof];
    dsassert(centernode!=NULL,"Cannot make sparsity pattern for Aztec");
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
          if (actnode->dof[l] < actfield->dis[dis].numeq)
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
    /*-------------- allocate the dof_nnz vector and put dofs in it */
    dof_nnz.a.iv[dof] = counter2+1;
  }  /* end of loop over numeq */ 
  /*--------------------------------------------- now do the coupled dofs */
  coupledofs = &(actpart->pdis[dis].coupledofs);
  for (i=0; i<coupledofs->fdim; i++)
  {
    dof = coupledofs->a.ia[i][0];
    /*--------------------------- check for my own ownership of this dof */
    dofflag = coupledofs->a.ia[i][imyrank+1];
    /*----------- if dofflag is zero this dof has nothing to do with me */
    if (dofflag==0) continue;
    /*------------------------------------- find all patches to this dof */
    counter=0;
    for (j=0; j<actpart->pdis[dis].numnp; j++)
    {
      centernode=NULL;
      for (l=0; l<actpart->pdis[dis].node[j]->numdf; l++)
      {
        if (dof == actpart->pdis[dis].node[j]->dof[l])
        {
          centernode = actpart->pdis[dis].node[j];
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
              if (actnode->dof[l] < actfield->dis[dis].numeq)
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
    dof_nnz.a.iv[dof] = counter2+1;
  } /* end of loop over coupled dofs */
  /* make the who-has-to-send-whom-how-much-and-what-arrays and communicate */
#ifdef PARALLEL 
  counter=0;
  for (i=0; i<coupledofs->fdim; i++)
  {
    dof = coupledofs->a.ia[i][0];
    /*-------------------------------------- find the master of this dof */
    for (j=1; j<coupledofs->sdim; j++)
    {
      if (coupledofs->a.ia[i][j]==2) 
      {
        dofmaster = j-1;
        break;
      }
    }
    /*-------------------------------------- find the slaves of this dof */
    for (j=1; j<coupledofs->sdim; j++)
    {
      if (coupledofs->a.ia[i][j]==1)
      {
        dofslave = j-1;
        /*----------------------------------- if I am master I receive */
        if (imyrank==dofmaster)
        {
          /* note:
             This is a nice example to do individual communication
             between two procs without communicating the size
             of the message in advance
             */
          /*--------------------------------- get envelope of message */
          MPI_Probe(dofslave,counter,actintra->MPI_INTRA_COMM,&status);
          /*----------------------------------------- receive message */
          MPI_Recv(&(recvlength),1,MPI_INT,
              dofslave,counter,actintra->MPI_INTRA_COMM,&status);
          /*--------------------------------- put new lenght to array */
          dof_nnz.a.iv[dof] += recvlength;
        }
        if (imyrank==dofslave)
        {
          MPI_Send(
              &(dof_nnz.a.iv[dof]),
              1,
              MPI_INT,
              dofmaster,
              counter,
              actintra->MPI_INTRA_COMM
              );
        }
        counter++;
      }
    }
  }
#endif
  /*--------------------------------- now go through update and count nnz */
  nnz=0;
  for (i=0; i<oll->update.fdim; i++)
  {
    dof = oll->update.a.iv[i];
    nnz += (dof_nnz.a.iv[dof]);
  }
  oll->nnz=nnz;
  /*----------------------------------------------------------------------*/
  amdel(&dofpatch);
  amdel(&dof_nnz);
  CCAFREE(node_dof);
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of oll_nnz_topology */




/*!----------------------------------------------------------------------
\brief create update for oll matrices

<pre>                                                              mn 02/03
This function creates the update vector for oll matrices. 
</pre>
\param *actfield       FIELD  (i)   the active field
\param *actpart        PARTITION  (i)   the active partition
\param *actintra       INTRA  (i)   the active communicator
\param  dis            INT    (i)   number of the active discretization
\param *oll            OLL    (o)   the oll matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void oll_update(
    FIELD         *actfield, 
    PARTITION     *actpart, 
    INTRA         *actintra,
    INT            dis,
    OLL           *oll)
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
  dstrc_enter("oll_update");
#endif
  /*----------------------------------------------------------------------*/
  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;
  /*------------------ make a local copy of the array actpart->coupledofs */
  am_alloc_copy(&(actpart->pdis[dis].coupledofs),&coupledofs);
  /*----------------------------------------------------------------------*/
  update = oll->update.a.iv;
  counter=0;
  /*------------------------------------- loop the nodes on the partition */
  for (i=0; i<actpart->pdis[dis].numnp; i++)
  {
    actnode = actpart->pdis[dis].node[i];
    for (l=0; l<actnode->numdf; l++)
    {
      dof = actnode->dof[l];
      /* dirichlet condition on dof */
      if (dof >= actfield->dis[dis].numeq) continue;
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
  if (counter != oll->numeq) dserror("Number of dofs in OLL-vector update wrong");
  /*---------------------------- sort the vector update just to make sure */
  qsort((INT*) update, counter, sizeof(INT), cmp_int);
  /*----------------------------------------------------------------------*/
  amdel(&coupledofs);
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of oll_update */




/*!---------------------------------------------------------------------
\brief return index of a given dof from local dof list 

<pre>                                                        m.gee 9/02 
return index of a given dof from local dof list. The given vector update
of length length has to be sorted and continous
</pre>
\param dof      INT    (i)           the dof the index is needed for 
\param update   INT*   (i)           the sorted and continous vector update 
\param length   INT    (i)           length of update 
\return the index of dof in update (INT) or -1 if index not in range
\sa 

------------------------------------------------------------------------*/
INT oll_getindex(
    INT dof,
    INT *update,
    INT length)
{
  INT index;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_enter("oll_getindex");
#endif
  /*----------------------------------------------------------------------*/
  index = dof-(*update);
  if (index<0 || index>=length) 
  {
#ifdef DEBUG 
    dstrc_exit();
#endif
    return(-1);
  }
  else
  {
#ifdef DEBUG 
    dstrc_exit();
#endif
    return(index);
  }
  /*----------------------------------------------------------------------*/
} /* end of oll_getindex */




/*!---------------------------------------------------------------------
\brief prints an oll matrix !!

<pre>                                                        mn 02/03 
This function prints the given oll matrix to the screen. Only for numeq
smaller than n_max.
</pre>
\param *oll      OLL    (i)    the oll matrix to print
\param  n_max    INT    (i)    max number of columns to print
\return void
\sa                                        

------------------------------------------------------------------------*/
void oll_print(
    OLL *oll,
    INT n_max)
{
  MATENTRY  **row;
  MATENTRY   *actentry;
  INT         lastcol;
  INT         gap;
  INT         i,j;
  DOUBLE      null = 0.0;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_enter("oll_print");
#endif
  /*----------------------------------------------------------------------*/
  row = oll->row;

  printf("dense System matrix written from oll format: \n");
  printf("============================================ \n");
  printf("rdim =  %4d \n", oll->rdim);
  printf("cdim =  %4d \n", oll->cdim);
  printf("nnz  =  %4d \n", oll->nnz);
  if(oll->rdim <= n_max)
  {
    for (i=0; i<oll->rdim; i++)
    {
      lastcol = -1;
      actentry = row[i];
      while(actentry != NULL )
      {
        gap = actentry->c - lastcol;
        if(gap>1)
        {
          /* there were (gap-1) zero entries before this entry */
          for (j=0;j<gap-1;j++)
          {
            printf(" %10.3E ", null);
          }
        } /* end if gap>1 */
        printf(" %10.3E ", actentry->val);
        lastcol = actentry->c;
        actentry = actentry->rnext;
      } /* end while row i */
      for (j=0; j<oll->numeq_total-lastcol-1; j++)
      {
        printf(" %10.3E ", null);
      }
      printf(" \n");
    } /* end for all rows */
  }
  printf("END of matrix \n");
  printf("============= \n");
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
/*----------------------------------------------------------------------*/
}  /* END of oll_print */




/*!---------------------------------------------------------------------
\brief  prints the pattern of an oll matrix

<pre>                                                        mn 02/03 
This function prints the sparsety pattern of an oll matrix to the screen.
</pre>
\param *oll      OLL    (i)    the oll matrix to print
\param  n_max    INT    (i)    max number of columns to print
\return void
\sa                                        

------------------------------------------------------------------------*/
void oll_pattern(
    OLL *oll,
    INT n_max)
{
  MATENTRY  **row;
  MATENTRY   *actentry;
  INT         lastcol;
  INT         gap;
  INT         i,j;
  INT         lehr,null,nn;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_enter("oll_pattern");
#endif
  /*----------------------------------------------------------------------*/
  lehr = 0;
  null =0;
  nn =0;
  row = oll->row;

  printf("dense System matrix written from oll format: \n");
  printf("============================================ \n");
  printf("rdim =  %4d \n", oll->rdim);
  printf("cdim =  %4d \n", oll->cdim);
  printf("nnz  =  %4d \n", oll->nnz);
  if(oll->rdim <= 2000)
  {
    for (i=0; i<oll->rdim; i++)
    {
      lastcol = -1;
      actentry = row[i];
      while(actentry != NULL )
      {
        gap = actentry->c - lastcol;
        if(gap>1)
        {
          /* there were (gap-1) zero entries before this entry */
          for (j=0;j<gap-1;j++)
          {
            printf(" ");
            lehr++;
          }
        } /* end if gap>1 */
        if( FABS(actentry->val) <= EPS12 ) 
        {
          printf("0");
          null++;
        }
        else
        {
          printf("X");
          nn++;
        }
        lastcol = actentry->c;
        actentry = actentry->rnext;
      } /* end while row i */
      for (j=0; j<oll->numeq_total-lastcol-1; j++)
      {
        printf(" ");
        lehr++;
      }
      printf(" \n");
    } /* end for all rows */
  }
  else
  {
    for (i=0; i<oll->rdim; i++)
    {
      lastcol = -1;
      actentry = row[i];
      while(actentry != NULL )
      {
        gap = actentry->c - lastcol;
        if(gap>1)
        {
          /* there were (gap-1) zero entries before this entry */
          for (j=0;j<gap-1;j++)
          {
            lehr++;
          }
        } /* end if gap>1 */
        if( FABS(actentry->val) <= EPS12 ) 
        {
          null++;
        }
        else
        {
          nn++;
        }
        lastcol = actentry->c;
        actentry = actentry->rnext;
      } /* end while row i */
      for (j=0; j<oll->numeq_total-lastcol-1; j++)
      {
        lehr++;
      }
    } /* end for all rows */
  }
  printf("lehr: %4d \n",lehr);
  printf("null: %4d \n",null);
  printf("nn:   %4d \n",nn);
  printf("END of matrix \n");
  printf("============= \n");
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
/*----------------------------------------------------------------------*/
}  /* END of oll_pattern */


/*! @} (documentation module close)*/
