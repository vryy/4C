#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h" 
/*----------------------------------------------------------------------*
 |  routine to calculate the element dirichlet load vector              |
 |                                                  genk 04/02          |
 *----------------------------------------------------------------------*/
void f2_caldirich(
                  ELEMENT   *actele, 
		  double    *dforces,
                  double   **estif, 
		  int       *hasdirich
		 )     
{

int                   i,j;
int                   dof;
int                   numdf;
int                   nd=0;
double                dirich[MAXDOFPERELE];
int                   dirich_onoff[MAXDOFPERELE];
GNODE                *actgnode;
NODE                 *actnode;

#ifdef DEBUG 
dstrc_enter("f2_caldirich");
#endif  

/*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;   
   if (actgnode->dirich==NULL) 
      continue;
   else
      *hasdirich=1;
      break;
}					  

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already *
/*                           written to the nodes (sol_increment[3][j]) */
for (i=0; i<actele->numnp; i++)
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];   
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++)
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[i*numdf+j] = actnode->sol_increment.a.da[3][j];
   }
}
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f2_caldirich*/ 
#endif
