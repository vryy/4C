#include "../headers/standardtypes.h"
#include "wall1.h"
/*----------------------------------------------------------------------*
 | initialize the element                                    al 6/01    |
 *----------------------------------------------------------------------*/
void w1init(PARTITION *actpart,MATERIAL    *mat )
{
int          i,j,k;
int          size_i, size_j;
ELEMENT     *actele;
NODE        *actnode;
W1_DATA      data;

#ifdef DEBUG 
dstrc_enter("w1init");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->numele; i++)
{
  actele = actpart->element[i];
  if (actele->eltyp != el_wall1) continue;
  /*---------------------------------------- init integration points ---*/
  w1intg(actele,&data,0);
  
  /*--------------------------------------------- init working array ---*/
  if(mat[actele->mat-1].mattyp == m_pl_mises || 
     mat[actele->mat-1].mattyp == m_pl_dp )
  {
    size_i = 1;
    actele->e.w1->elewa = (W1_ELE_WA*)calloc(size_i,sizeof(W1_ELE_WA));
    if (actele->e.w1->elewa==NULL)
    {
      dserror("Allocation of elewa in ELEMENT failed");
      break;
    } 
    
    size_j = actele->e.w1->nGP[0] + actele->e.w1->nGP[1];
    actele->e.w1->elewa[0].ipwa = 
                               (W1_IP_WA*)calloc(size_j,sizeof(W1_IP_WA));
    if (actele->e.w1->elewa[0].ipwa==NULL)
    {
      dserror("Allocation of ipwa in ELEMENT failed");
      break;
    } 
    for (k=0; k<size_j; k++)
    {
      actele->e.w1->elewa[0].ipwa[k].epstn = 0.0;
      actele->e.w1->elewa[0].ipwa[k].yip   = -1;
      for (j=0; j<4; j++)
      {
        actele->e.w1->elewa[0].ipwa[k].sig[j] = 0.0;
        actele->e.w1->elewa[0].ipwa[k].eps[j] = 0.0;
        actele->e.w1->elewa[0].ipwa[k].qn[ j] = 0.0;
      }
    }
    
  
  
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1init */
 
