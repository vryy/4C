/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_init' which initializes the element

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief initializes the element

<pre>                                                              fh 09/02
This routine initializes the actual beam element

</pre>
\param *actpart  PARTITION (i)  actual partition
\param *mat      MATERIAL  (i)  actual material
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_intg()
    called by: beam3() 

*----------------------------------------------------------------------*/
void b3_init(PARTITION *actpart,
             MATERIAL  *mat )
{
INT          i,j,k;
ELEMENT     *actele;
B3_DATA	     data;

#ifdef DEBUG 
dstrc_enter("b3_init");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->pdis[0].numele; i++) 
{ 
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_beam3) continue;
/* check number of hinges and write it to the first entry of hc = hc[0] */
  actele->e.b3->hc[0]=0;  
  for (j=0; j<12; j++)
  {
     if (actele->e.b3->hc[j+1] != 0) actele->e.b3->hc[0]+=1;
  } 
/*------ modifiy hinge code (hinge code only valid at design nodes) ----*/
  for (j=0; j<actele->numnp; j++)
  {
     if (actele->node[j]->gnode->d.dnode->ndline==0)
     {
        actele->e.b3->hc[0]=0;
	for (k=j*6+1; k<(j+1)*6+1; k++) actele->e.b3->hc[k]=0;
     }
  }  
/*---------------------------------------- init integration points -----*/
  b3_intg(actele,&data,0);
/*----------------------- allocate the space for internal forces -------*/
  am4def("force_GP",&(actele->e.b3->force_GP),1,MAXDOFPERNODE,MAXGAUSS,0,"D3");
  am4def("force_ND",&(actele->e.b3->force_ND),1,MAXDOFPERNODE,MAXNOD_BEAM3,0,"D3");
  am4zero(&(actele->e.b3->force_GP));
  am4zero(&(actele->e.b3->force_ND));
   
  /* for plasticity */
  if(mat[actele->mat-1].mattyp == m_pl_mises || 
     mat[actele->mat-1].mattyp == m_pl_dp || 
     mat[actele->mat-1].mattyp == m_pl_epc )
  {
    amdef("elewa",&(actele->e.b3->elewa),actele->e.b3->nGP[0]*5*5,40,"DA");
    amzero(&(actele->e.b3->elewa));
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b3_init */
#endif
/*! @} (documentation module close)*/
