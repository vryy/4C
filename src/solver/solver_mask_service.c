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

static INT disnum;

/*----------------------------------------------------------------------*
 |  count processor local and global number of equations    m.gee 5/01  |
 *----------------------------------------------------------------------*/
void mask_numeq(
    FIELD         *actfield,
    PARTITION    *actpart,
    SOLVAR       *actsolv,
    INTRA        *actintra,
    INT          *numeq,
    INT           dis_num)
{
  INT       i,j,k,l;
  INT       counter;
  INT       dof;
  INT       iscoupled;
#ifdef PARALLEL
  INT      *sendbuff,*recvbuff, sendsize;
#endif
  INT      *tmp;
  INT       inter_proc;
  long int  min;
  INT       proc;
  INT       inprocs;
  INT       imyrank;
  NODE     *actnode;

  INT       no_coupling = 0;

#ifdef DEBUG
  dstrc_enter("mask_numeq");
#endif

  /* set actual discretisation */
  disnum=dis_num;


  imyrank = actintra->intra_rank;
  inprocs = actintra->intra_nprocs;

  /* first make a list of dofs which are coupled */
  /* estimate size of coupdofs to 5000 */
  amdef("coupledofs",&(actpart->pdis[disnum].coupledofs),5000,1,"IV");
  amzero(&(actpart->pdis[disnum].coupledofs));

  counter=0;
  /* loop all nodes and find coupled dofs */
  for (i=0; i<actfield->dis[disnum].numnp; i++)
  {
    actnode = &(actfield->dis[disnum].node[i]);
    if (actnode->gnode->couple==NULL && actnode->gnode->dirich==NULL) continue;
    if (actnode->gnode->couple==NULL) continue;
    for (l=0; l<actnode->numdf; l++)
    {
      if (actnode->dof[l]>=actfield->dis[disnum].numeq) continue;
      /* there is coupling on this dof */
      if (actnode->gnode->couple->couple.a.ia[l][0] != 0 ||
          actnode->gnode->couple->couple.a.ia[l][1] != 0 )
      {
        if (counter>=actpart->pdis[disnum].coupledofs.fdim)
          amredef(&(actpart->pdis[disnum].coupledofs),(actpart->pdis[disnum].coupledofs.fdim+5000),1,"IV");
        /* the coupled dof could be dirichlet conditioned */
        if (actnode->dof[l]<actfield->dis[disnum].numeq)
        {
          actpart->pdis[disnum].coupledofs.a.iv[counter] = actnode->dof[l];
          counter++;
        }
      }
    }
  }

  if (counter ==0)
  {
    no_coupling = 1;
    amdel(&(actpart->pdis[disnum].coupledofs));
  }
  else
  {
    no_coupling = 0;
    amredef(&(actpart->pdis[disnum].coupledofs),counter,1,"IV");
  }
  
  /* delete the doubles in coupledofs */
  if (!no_coupling)
  {
    for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
    {
      if (actpart->pdis[disnum].coupledofs.a.iv[i]==-1) continue;
      dof = actpart->pdis[disnum].coupledofs.a.iv[i];
      for (j=i+1; j<actpart->pdis[disnum].coupledofs.fdim; j++)
      {
        if (actpart->pdis[disnum].coupledofs.a.iv[j]==dof)
          actpart->pdis[disnum].coupledofs.a.iv[j]=-1;
      }
    }

    /* move all remaining coupdofs to the front and redefine again */
    counter=0;
    for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
    {
      if (actpart->pdis[disnum].coupledofs.a.iv[i]!=-1)
      {
        actpart->pdis[disnum].coupledofs.a.iv[counter] = actpart->pdis[disnum].coupledofs.a.iv[i];
        counter++;
      }
    }
    amredef(&(actpart->pdis[disnum].coupledofs),counter,inprocs+1,"IA");

    /* the newly allocated columns have to be initialized */
    for (i=1; i<actpart->pdis[disnum].coupledofs.sdim; i++)
      for (j=0; j<actpart->pdis[disnum].coupledofs.fdim; j++)
        actpart->pdis[disnum].coupledofs.a.ia[j][i]=0;

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

               column 0             : number of the coupled equation
               column 1 - inprocs+1 : proc has coupled equation or not

*/

  if (!no_coupling)
  {
    if (inprocs==1) /* sequentiell version */
    {
      for (k=0; k<actpart->pdis[disnum].coupledofs.fdim; k++)
      {
        actpart->pdis[disnum].coupledofs.a.ia[k][imyrank+1]=2;
      }
    }
    else /* parallel version */
    {
      /* actpart->node[i] really loops only nodes with dofs updated on this proc */

      /* now loop only my nodes */
      for (i=0; i<actpart->pdis[disnum].numnp; i++)
      {
        for (l=0; l<actpart->pdis[disnum].node[i]->numdf; l++)
        {
          dof = actpart->pdis[disnum].node[i]->dof[l];
          for (k=0; k<actpart->pdis[disnum].coupledofs.fdim; k++)
          {
            if (actpart->pdis[disnum].coupledofs.a.ia[k][0]==dof)
            {
              actpart->pdis[disnum].coupledofs.a.ia[k][imyrank+1]=1;
              break;
            }
          }
        }
      }
    }

    /* Allreduce the whole array, so every proc knows about where all
       coupledofs are */
#ifdef PARALLEL
    sendsize = (actpart->pdis[disnum].coupledofs.fdim)*(inprocs);
    sendbuff = (INT*)CCACALLOC(sendsize,sizeof(INT));
    recvbuff = (INT*)CCACALLOC(sendsize,sizeof(INT));
    if (sendbuff==NULL || recvbuff==NULL) dserror("Allocation of temporary memory failed");
    counter=0;
    for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
    {
      for (j=0; j<inprocs; j++)
      {
        sendbuff[counter] = actpart->pdis[disnum].coupledofs.a.ia[i][j+1];
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
    for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
    {
      for (j=0; j<inprocs; j++)
      {
        actpart->pdis[disnum].coupledofs.a.ia[i][j+1] = recvbuff[counter];
        counter++;
      }
    }
    CCAFREE(sendbuff);CCAFREE(recvbuff);
#endif
  } /* end of if(!no_coupling) */


  /* count number of equations on partition including coupled dofs */
  /* count the coupled ones first */
  counter=0;
  for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
  {
    if (actpart->pdis[disnum].coupledofs.a.ia[i][imyrank+1]!=0) counter++;
  }

  /* count all dofs which are not coupled */
  for (i=0; i<actpart->pdis[disnum].numnp; i++)
  {
    actnode = actpart->pdis[disnum].node[i];
    for (l=0; l<actnode->numdf; l++)
    {
      dof = actnode->dof[l];
      iscoupled=0;
      for (k=0; k<actpart->pdis[disnum].coupledofs.fdim; k++)
      {
        if (dof == actpart->pdis[disnum].coupledofs.a.ia[k][0])
        {
          iscoupled=1;
          break;
        }
      }
      if (iscoupled==0)
      {
        if (dof < actfield->dis[disnum].numeq)
          counter++;
      }
    }
  }

  /* number of equations on this partition including the coupled ones */
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

      /* loop coupled eqns */
      for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
      {
        /* check whether its inter-proc eqn */
        inter_proc=0;
        for (j=0; j<inprocs; j++)
          inter_proc += actpart->pdis[disnum].coupledofs.a.ia[i][j+1];

        if (inter_proc==1)/* no inter-processor coupling */
        {
          for (j=0; j<inprocs; j++)
          {
            if (actpart->pdis[disnum].coupledofs.a.ia[i][j+1]==1)
            {
              actpart->pdis[disnum].coupledofs.a.ia[i][j+1]=2;
              break;
            }
          }
        }
        else/* eqn is an inter-proc equation */
        {
          /* there won't be more than a million procs in the near future....*/
          min=1000000;
          proc=-1;
          for (j=0; j<inprocs; j++)
          {
            if (actpart->pdis[disnum].coupledofs.a.ia[i][j+1]==1)
            {
              if (tmp[j]<=min)
              {
                min = tmp[j];
                proc = j;
              }
            }
          }
          actpart->pdis[disnum].coupledofs.a.ia[i][proc+1]=2;
          tmp[proc] += 1;
        }
      }/* end loop over coupling eqns */
      CCAFREE(tmp);
    }

    /* procs who have not become owner of a coupling equation have to reduce there
       number of equations */
    if (inprocs > 1)
    {
      for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)/* loop coupled eqns */
      {
        /* Yes, I am slave owner of an inter_proc coupling equation */
        if (actpart->pdis[disnum].coupledofs.a.ia[i][imyrank+1]==1)
        {
          (*numeq) = (*numeq)-1;
        }
        /* master owner of equation do nothing, 'cause the equation has been
           counted already */
      }
    }
  } /* end of if(!no_coupling) */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of mask_numeq */





/*----------------------------------------------------------------------*
 |  check whether this dof is in coupledofs              m.gee 6/01     |
 *----------------------------------------------------------------------*/
void dof_in_coupledofs(
    INT           dof,
    PARTITION    *actpart,
    INT          *iscoupled)
{

  INT       i;

#ifdef DEBUG
  dstrc_enter("dof_in_coupledofs");
#endif

  for (i=0; i<actpart->pdis[disnum].coupledofs.fdim; i++)
  {
    if (dof==actpart->pdis[disnum].coupledofs.a.ia[i][0])
    {
      *iscoupled = 1;
      break;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of dof_in_coupledofs */




/*----------------------------------------------------------------------*
 |  find the node to this dof in partition               m.gee 6/01     |
 *----------------------------------------------------------------------*/
void dof_find_centernode(
    INT          dof,
    PARTITION   *actpart,
    NODE       **centernode)
{

  INT       j,k;

#ifdef DEBUG
  dstrc_enter("dof_find_centernode");
#endif

  for (j=0; j<actpart->pdis[disnum].numnp; j++)
  {
    for (k=0; k<actpart->pdis[disnum].node[j]->numdf; k++)
    {
      if (actpart->pdis[disnum].node[j]->dof[k] == dof)
      {
        *centernode = actpart->pdis[disnum].node[j];
        goto nodefound1;
      }
    }
  }
nodefound1:

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of dof_find_centernode */


