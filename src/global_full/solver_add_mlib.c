#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include <veclib.h> 
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in globcalelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to assemble element array to global                         |
 |  column pointer, row index sparse  matrix                            |
 |  sequentiell for hp's mlib solver                                    |
 |                                                                      |
 |                                                                      |
 |                                                         al   10/01   |
 *----------------------------------------------------------------------*/
int  add_mds(struct _PARTITION     *actpart,
             struct _SOLVAR        *actsolv,
             struct _ELEMENT       *actele,
             struct _ML_ARRAY_MDS  *mds)
{
  int         i,j,k,l,counter; /* some counter variables */
  int         ii,jj;           /* counter variables for system matrix */
  int         nd;              /* size of estif */
  int         numeq;           /* number of equations on this proc */
  int         lm[MAXDOFPERELE];/* location vector for this element */
  double      value;
  double    **estif;     /* element matrix to be added to system matrix */
  MLVAR       *mlvar;
/*----------------------------------------------------------------------*/
  mlvar = actsolv->mlvar;
/*----------------------------------------------------------------------*/
  #ifdef DEBUG 
  dstrc_enter("add_mds");
  #endif
/*------------------------------------- set some pointers and variables */
  estif      = estif_global.a.da;
  nd         = actele->numnp * actele->node[0]->numdf;
  numeq      = mds->numeq;
/*---------------------------------------------- make location vector lm*/
  counter=0;
  for (i=0; i<actele->numnp; i++)
  {
    for (j=0; j<actele->node[i]->numdf; j++)
    {
      lm[counter]    = actele->node[i]->dof[j];
      #ifdef PARALLEL 
      owner[counter] = actele->node[i]->proc;
      #endif
      counter++;
    }
  }/* end of loop over element nodes */
  if (counter != nd) dserror("assemblage failed due to wrong dof numbering");
/*========================================== now start looping the dofs */

/*======================================= loop over i (the element row) */
  for (i=0; i<nd; i++)
  {
   ii = lm[i]; 
   /*------------------------------------- check for boundary condition */
   if (ii>=numeq) continue;
   /*================================= loop over j (the element column) */
   for (j=0; j<nd; j++)
   {
      jj = lm[j];
      /*---------------------------------- check for boundary condition */
      if (jj>=numeq) continue;
      
      if(ii>=jj && mlvar->symm || !mlvar->symm) /* lower triangular of full*/
      {
        value = estif[i][j];
      
        ii++;
        jj++;
        dslev1 (&ii, &jj, &value, mds->global, &mds->ierr);
        ii--;
        jj--;
      }
      if(mds->ierr!=0) exit; 
    } /* end loop over j */
  }/* end loop over i */
/*-- print additional information, depends on the stage of execution ---*/
  if(mlvar->msglvl==4) dsleps (mds->global); 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return(0);
} /* end of add_mds */


