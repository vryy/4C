/*!----------------------------------------------------------------------
\file
\brief contains the routine 'optsmo',
       smoothing routine for gradients, element densities...

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../headers/optimization.h"
/*! 
\addtogroup OPTIMIZATION 
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01   
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;

/*----------------------------------------------------------------------*
 | variational sensitivity analysis                         al 05/01    |
 *----------------------------------------------------------------------*/
void optsmo(double *vvar, int init)
{
/*----------------------------------------------------------------------*/
  int    i, j, k;               /* some counters */
/*----------------------------------------------------------------------*/
  int iloc;
  int    maxneighbour = 1000;
  int    ngbevec[1000];
  double ngbeval[1000];
  double dist, sstab, rstb, rzl, sumo;
  double cpsum[4];
  double **cpele;
  
  static  int     **nbele;
  static  double  **nbelv;
  static  double  *stbobj;
/*----------------------------------------------------------------------*/
  int           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */
  SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
  PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
  FIELD        *actfield;         /* pointer to the structural FIELD */
  INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
  CONTAINER     container;        /* contains variables defined in container.h */
  ELEMENT *actele;                /* active element                            */
  ELEMENT *nctele;                /* another active element                    */
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("optsmo");
  #endif
/*------------ the distributed system matrix, which is used for solving */
  actsysarray=0;
/*--------------------------------------------------- set some pointers */
  actfield    = &(field[0]);
  actsolv     = &(solv[0]);
  actpart     = &(partition[0]);
  action      = &(calc_action[0]);
  #ifdef PARALLEL 
  actintra    = &(par.intra[0]);
  #else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  if (!actintra) dserror("Allocation of INTRA failed");
  actintra->intra_fieldtyp = structure;
  actintra->intra_rank   = 0;
  actintra->intra_nprocs   = 1;
  #endif
  container.fieldtyp  = actfield->fieldtyp;
  container.isdyn = 0;            /* static calculation */
/*----------------------------------------------------------------------*/
  if(init==1)
  {
    stbobj  = (double*)CCACALLOC(actfield->dis[0].numele,sizeof(double));
    /*------------------------ evaluate coordinates of center points ---*/
    cpele  = (double**)CCACALLOC(actfield->dis[0].numele,sizeof(double*));
    for (i=0; i<actfield->dis[0].numele; i++)
    {
      cpele[i]  = (double*)CCACALLOC(3,sizeof(double));
    }
    /*----------*/
    if(opt->smoothtype==sm_grad)
    {
      for (i=0; i<actfield->dis[0].numele; i++)
      { /* loop elements */
        cpsum[0] = 0.;
        cpsum[1] = 0.;
        cpsum[2] = 0.;
        cpsum[3] = 0.;
        actele = &(actfield->dis[0].element[i]);
        for (j=0; j<actele->numnp; j++)
        { /* loop nodes of elements */
          cpsum[0] += actele->node[j]->x[0];
          cpsum[1] += actele->node[j]->x[1];
          cpsum[2] += actele->node[j]->x[2];
          cpsum[3] += 1.;
        }
        cpele[i][0] = cpsum[0]/cpsum[3];
        cpele[i][1] = cpsum[1]/cpsum[3];
        cpele[i][2] = cpsum[2]/cpsum[3];
      }
    }
    /*---- determine all elements in the neigbourhood of one element ---*/
    nbele = (int**   )CCACALLOC(actfield->dis[0].numele,sizeof(int*   ));
    nbelv = (double**)CCACALLOC(actfield->dis[0].numele,sizeof(double*));
    
    
    for (i=0; i<actfield->dis[0].numele; i++)
    { /* loop elements */
      /* */
      actele = &(actfield->dis[0].element[i]);
      if(actele->optdata==NULL) continue; /* element does not take part in opt. */
      if(actele->optdata[0]==0) continue; /* position in variable vector        */
      /* */
      ngbevec[0]=0;
      for (j=0; j<actfield->dis[0].numele; j++)
      { /* loop elements */
        /* */
        nctele = &(actfield->dis[0].element[j]);
        if(nctele->optdata==NULL) continue; /* element does not take part in opt. */
        if(nctele->optdata[0]==0) continue; /* position in variable vector        */
        /* */
        dist=0.;
        for (k=0; k<3; k++) dist+=(cpele[i][k]-cpele[j][k])*(cpele[i][k]-cpele[j][k]);
        dist=sqrt(dist);
        if(dist<opt->smoothrad)
        {
          if(ngbevec[0]>maxneighbour) break;
          ngbevec[ngbevec[0]+1] = j+1;
          sstab = fabs(2.0-dist/opt->smoothrad);
          sstab = pow(sstab,opt->smoothexp);
          
          ngbeval[ngbevec[0]+1] = sstab;
          ngbevec[0]++;
        }
      }
      nbele[i]  = (int*   )CCACALLOC(ngbevec[0]+1,sizeof(int   ));
      nbelv[i]  = (double*)CCACALLOC(ngbevec[0]+1,sizeof(double));
      nbele[i][0]  = ngbevec[0]; /* store number of neighbours */
      for (j=1; j<=ngbevec[0]; j++)
      {
        nbele[i][j] = ngbevec[j];
        nbelv[i][j] = ngbeval[j];
      }
    }




      goto end;
  }
/*----------------------------------------------------------------------*/
/*---------------------------------- smooth objecitves and gradients ---*/
  if(opt->smoothtype==sm_grad)
  {
    for (i=0; i<actfield->dis[0].numele; i++)
    {
      /* */
      actele = &(actfield->dis[0].element[i]);
      if(actele->optdata==NULL) continue; /* element does not take part in opt. */
      if(actele->optdata[0]==0) continue; /* position in variable vector        */
      /* */
      rzl  = 0.0;
      sumo = 0.0;
      for (j=1; j<=nbele[i][0]; j++)
      {
        iloc = nbele[i][j]; /* element Id */
        rstb = nbelv[i][j];
        
        nctele = &(actfield->dis[0].element[iloc-1]);
        iloc   = nctele->optdata[0];
        
        rzl  += rstb;
        sumo += rstb * vvar[iloc-1];
      }
            if(rzl!=0.0) stbobj[i] = sumo / rzl; 
            else         stbobj[i] = 0.0; 
    }
  
    for (i=0; i<actfield->dis[0].numele; i++)
    {
      actele = &(actfield->dis[0].element[i]);
      if(actele->optdata==NULL) continue; /* element does not take part in opt. */
      if(actele->optdata[0]==0) continue; /* position in variable vector        */
      vvar[actele->optdata[0]-1] = stbobj[i];
    }
  }
/*----------------------------------------------------------------------*/
  end: ;
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_exit();
  #endif
return; 
} /* end of optsmo */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
