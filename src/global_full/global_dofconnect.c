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
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      


/*----------------------------------------------------------------------*
 | prototypes of functions callable only in this file                   |
 *----------------------------------------------------------------------*/
void dofdof_realloc(INT dof1,INT **dof_dof);
INT cmp(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |  put dofs in update in ascending order                 a.lipka 5/01  |
 *----------------------------------------------------------------------*/
/* compare the integers - qsort routine*/
INT cmp(const void *a, const void *b )
{
    return *(INT *)a - * (INT *)b;
}
/*----------------------------------------------------------------------*
 |  realloc memory in dof-dof connectivity list           a.lipka 5/01  |
 *----------------------------------------------------------------------*/
void dofdof_realloc(INT dof1,INT **dof_dof)
{
/*----------------------------------------------------------------------*/
  INT i=0, j=0;
  INT ac=0;
  INT *hivec;
  INT numdofconected    = 100;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("dofdof_realloc");
  #endif
/*----------------------------------------------------------------------*/
  ac = dof_dof[dof1][1];                /* number of allocated integers */
  hivec = (INT*)CCACALLOC(ac ,sizeof(INT));
  for (i=0; i<ac; i++) hivec[i] = dof_dof[dof1][i];/*copy to 2nd vector */
  CCAFREE(dof_dof[dof1]);                                       /* realloc */
  dof_dof[dof1] = (INT*)CCACALLOC(ac+numdofconected ,sizeof(INT));
  for (j=0; j<ac; j++) dof_dof[dof1][j] = hivec[j];         /*copy back */
  dof_dof[dof1][1]  += numdofconected;  /* number of allocated integers */
  CCAFREE(hivec);
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of dofdof_realloc */
/*----------------------------------------------------------------------*
 |  eleminate redundant dof's in connectivity list        a.lipka 5/01  |
 *----------------------------------------------------------------------*/
void delete_redundant_dof(INT dof1,INT **dof_dof)
{
/*----------------------------------------------------------------------*/
  INT m=0, k=0;
  INT dof11, dof22;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("delete_redundant_dof");
  #endif
/*----------------------------------------------------------------------*/
  dof11 = dof_dof[dof1][2];
  for (m=3; m<dof_dof[dof1][0]+2; m++)
  {
    dof22 = dof_dof[dof1][m];
    if(dof11==dof22)
    {
      m--;
      for (k=m; k<dof_dof[dof1][0]+1; k++)
      {
        dof_dof[dof1][k]=dof_dof[dof1][k+1];
      }
      dof_dof[dof1][0]--;
    }
    else
    {
      dof11 = dof22;
    }
  }
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of delete_redundant_dof */

/*----------------------------------------------------------------------*
 |  eleminate one dof in connectivity list                a.lipka 5/01  |
 *----------------------------------------------------------------------*/
void delete_single_dof(INT dof1,INT **dof_dof)
{
/*----------------------------------------------------------------------*/
  INT m=0, k=0;
  INT dof22;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("delete_single_dof");
  #endif
/*----------------------------------------------------------------------*/
  for (m=2; m<dof_dof[dof1][0]+2; m++)
  {
    dof22 = dof_dof[dof1][m];
    if(dof1==dof22)
    {
      for (k=m; k<dof_dof[dof1][0]+1; k++)
      {
        dof_dof[dof1][k]=dof_dof[dof1][k+1];
      }
      dof_dof[dof1][0]--;
    }
  }
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
/*----------------------------------------------------------------------*/
return;
} /* end of delete_single_dof */
/*----------------------------------------------------------------------*
 |  calculate dof topology                                   al  10/01  |
 *----------------------------------------------------------------------*/
void dofconnectivity(FIELD        *actfield, 
                        INT      **dof_dof,
                        INT       *nnz)
{
/*----------------------------------------------------------------------*/
  INT j=0, k=0, l=0, n=0;
  INT dof1       = 0;
  INT dof2       = 0;
  INT nel        = 0;
  INT nod        = 0;
  INT hnod       = 0;
  INT cc         = 0;
  INT numdofconected    = 100;

  NODE    *actnode, *partnernode;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("dofconnectivity");
  #endif
/*----------------------------------------------------------------------*/
  for (j=0; j<actfield->dis[0].numnp; j++)
  {
    actnode = &(actfield->dis[0].node[j]);
    if(actnode->proc!=par.myrank) continue;
    /* loop dofs on this node */
    for (l=0; l<actnode->numdf; l++)
    {/*dofloop01*/
      dof1 = actnode->dof[l];
      if(dof1>=actfield->dis[0].numeq) continue;
      if(dof_dof[dof1]==NULL)
      {
        dof_dof[dof1]    = (INT*)CCACALLOC(numdofconected ,sizeof(INT));
        dof_dof[dof1][0] = 0;
        dof_dof[dof1][1] = numdofconected;
      }
      /* get number of dofs on this node */
      for (k=0; k<actnode->numdf; k++)
      {
        dof2 = actnode->dof[k];
       /* if(dof1!=dof2 && dof2<actfield->numeq)*/
        if(dof2<actfield->dis[0].numeq)
        {
          cc = dof_dof[dof1][0];
          if(cc>=dof_dof[dof1][1]-2)  dofdof_realloc(dof1,dof_dof); /*realloc*/
          dof_dof[dof1][cc+2] = dof2;
          dof_dof[dof1][0]++;
        }
      }
      /* loop neighbour elements */
      for (nel=0; nel<actnode->numele; nel++)
      {
        /* loop neighbour nodes */
        for (nod=0; nod<actnode->element[nel]->numnp; nod++)
        {
 
          partnernode = actnode->element[nel]->node[nod];
          hnod = partnernode->Id;
          if(hnod==actnode->Id) continue;
          for (n=0; n<partnernode->numdf; n++)
          {
            dof2 = partnernode->dof[n];
            if(dof2>=actfield->dis[0].numeq) continue;
            cc = dof_dof[dof1][0];
            if(cc>=dof_dof[dof1][1]-2) dofdof_realloc(dof1, dof_dof); /*realloc*/
            dof_dof[dof1][cc+2] = dof2;
            dof_dof[dof1][0]++;
          }
        }
      }
      /* sort dof numbers */
      qsort((INT*) &dof_dof[dof1][2], dof_dof[dof1][0], sizeof(INT), cmp);     
      /* eliminate redundant numbers */
      delete_redundant_dof(dof1, dof_dof);
      /* eliminate actual dof number from list */
      /*delete_single_dof(dof1,dof_dof);*/
    }/*dofloop01*/
  }
/*---------------------------------------------------------- count nonzeroes---*/

  *nnz = 0;
  for (dof1=0; dof1<actfield->dis[0].numeq; dof1++)
  {
    *nnz += dof_dof[dof1][0];
  }
/*-----------------------------------------------------------------------------*/
  /* matrix: dof - dof connectivity */
/*  sprintf(filename,"msr_proc%d_1.txt",par.myrank);
  
 
  if ((filep = fopen(filename, "w")) != NULL)
  {
    for(k=0;k<actfield->numeq;k++)
    {
      dof1 = k;         
      if(dof_dof[dof1][0]==0)continue;
    
      fprintf(filep, "dof %4d - ndof %3d - dofs  ",dof1,dof_dof[dof1][0]);
      for(m=2;m<dof_dof[dof1][0]+2;m++)
        {
        fprintf(filep, " %3d",dof_dof[dof1][m]);
        }
      fprintf(filep, "\n");
    }
    fclose(filep);
  }
  */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dofconnectivity */
/*----------------------------------------------------------------------*
 |                                                            al  10/01 |
 |  calculate the mask of a column pointer, row index sparse  matrix    |
 |  -symmetric- lower triangular part                                   |
 *----------------------------------------------------------------------*/
void mds_make_colstr_rowind(SOLVAR       *actsolv,
                            ML_ARRAY_MDS  *mds,
                            INT       **dof_dof,
                            INT         numeq)
{
  INT        i,j;
  INT        count;
  INT        dof1, dof2;
  INT        ncdof;
  INT        symm;
  MLVAR        *mlvar;
 
  #ifdef DEBUG
  dstrc_enter("mds_make_colstr_rowind");
  #endif
/*----------------------------------------------------------------------*/
  mlvar = actsolv->mlvar;
  symm = mlvar->symm;
/*----------------------------------------------------------------------*/
    count  = 0;
  for (i=0; i<numeq; i++) /* loop columns */
  {
    dof1 = i+1;
    mds->colstr.a.iv[i]=count+1; 
    
    ncdof = dof_dof[i][0];/* loop connected dofs */
    for (j=0; j<ncdof; j++) 
    {
    
      dof2 = dof_dof[i][j+2]+1;
      if(symm && dof2>=dof1)/*-symmetric- lower triangular part*/
      {
        mds->rowind.a.iv[count++] = dof2;
      }
      else if(!symm&&dof1!=dof2)
      {
        mds->rowind.a.iv[count++] = dof2;
      }

    }
  }
  mds->colstr.a.iv[numeq]=count+1;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
  return;
} /* end of mds_make_colstr_rowind */
/*----------------------------------------------------------------------*/

















