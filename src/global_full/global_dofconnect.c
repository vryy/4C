#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |  put dofs in update in ascending order                 a.lipka 5/01  |
 *----------------------------------------------------------------------*/
/* compare the integers - qsort routine*/
int cmp(const void *a, const void *b )
{
    return *(int *)a - * (int *)b;
}
/*----------------------------------------------------------------------*
 |  realloc memory in dof-dof connectivity list           a.lipka 5/01  |
 *----------------------------------------------------------------------*/
void dofdof_realloc(int dof1,int **dof_dof)
{
/*----------------------------------------------------------------------*/
  int i=0, j=0;
  int ac=0;
  int *hivec;
  int numdofconected    = 100;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("dofdof_realloc");
  #endif
/*----------------------------------------------------------------------*/
  ac = dof_dof[dof1][1];                /* number of allocated integers */
  hivec = (int*)CALLOC(ac ,sizeof(int));
  for (i=0; i<ac; i++) hivec[i] = dof_dof[dof1][i];/*copy to 2nd vector */
  FREE(dof_dof[dof1]);                                       /* realloc */
  dof_dof[dof1] = (int*)CALLOC(ac+numdofconected ,sizeof(int));
  for (j=0; j<ac; j++) dof_dof[dof1][j] = hivec[j];         /*copy back */
  dof_dof[dof1][1]  += numdofconected;  /* number of allocated integers */
  FREE(hivec);
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
void delete_redundant_dof(int dof1,int **dof_dof)
{
/*----------------------------------------------------------------------*/
  int m=0, k=0;
  int dof11, dof22;
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
void delete_single_dof(int dof1,int **dof_dof)
{
/*----------------------------------------------------------------------*/
  int m=0, k=0;
  int dof11, dof22;
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
                        int      **dof_dof,
                        int       *nnz)
{
/*----------------------------------------------------------------------*/
  int i=0, j=0, k=0, l=0, m=0, n=0, o=0;
  int dof1       = 0;
  int dof2       = 0;
  int nel        = 0;
  int nod        = 0;
  int hc1        = 0;
  int hc2        = 0;
  int hnod       = 0;
  int cc         = 0;
  int numdofconected    = 100;

  NODE    *actnode, *partnernode;

  FILE *filep      ;
  char filename[50];
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("dofconnectivity");
  #endif
/*----------------------------------------------------------------------*/
  for (j=0; j<actfield->numnp; j++)
  {
    actnode = &(actfield->node[j]);
    if(actnode->proc!=par.myrank) continue;
    /* loop dofs on this node */
    for (l=0; l<actnode->numdf; l++)
    {/*dofloop01*/
      dof1 = actnode->dof[l];
      if(dof1>=actfield->numeq) continue;
      if(dof_dof[dof1]==NULL)
      {
        dof_dof[dof1]    = (int*)CALLOC(numdofconected ,sizeof(int));
        dof_dof[dof1][0] = 0;
        dof_dof[dof1][1] = numdofconected;
      }
      /* get number of dofs on this node */
      for (k=0; k<actnode->numdf; k++)
      {
        dof2 = actnode->dof[k];
       /* if(dof1!=dof2 && dof2<actfield->numeq)*/
        if(dof2<actfield->numeq)
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
            if(dof2>=actfield->numeq) continue;
            cc = dof_dof[dof1][0];
            if(cc>=dof_dof[dof1][1]-2) dofdof_realloc(dof1, dof_dof); /*realloc*/
            dof_dof[dof1][cc+2] = dof2;
            dof_dof[dof1][0]++;
          }
        }
      }
      /* sort dof numbers */
      qsort((int*) &dof_dof[dof1][2], dof_dof[dof1][0], sizeof(int), cmp);     
      /* eliminate redundant numbers */
      delete_redundant_dof(dof1, dof_dof);
      /* eliminate actual dof number from list */
      /*delete_single_dof(dof1,dof_dof);*/
    }/*dofloop01*/
  }
/*---------------------------------------------------------- count nonzeroes---*/

  *nnz = 0;
  for (dof1=0; dof1<actfield->numeq; dof1++)
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
                            int       **dof_dof,
                            int         numeq)
{
  int        i,j,k,l;
  int        count1,count2;
  int        count;
  int        dof1, dof2;
  int        ncdof;
  int        symm;
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

















