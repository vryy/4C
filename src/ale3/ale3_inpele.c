/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3inp' which reads a 3d ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "ale3.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief reads a 3d ale element from the input file

<pre>                                                              mn 06/02 
This routine reads a 3d ale element from the input file

</pre>
\param *ele  ELEMENT  (o)   the element

\warning There is nothing special to this routine
\return void                                               
\sa caling: ---; called by: inp_ale_field()

*----------------------------------------------------------------------*/
void ale3inp(ELEMENT *ele)
{
INT  i;
INT  ierr=0;
INT  quad;
INT  dum[8];
INT  counter;
INT  lmtmp;
/* wird nicht gebraucht!!  */
/* long int  topology[100];*/
CHAR *colpointer;
CHAR buffer[50];
#ifdef DEBUG 
dstrc_enter("ale3inp");
#endif
/*------------------------------------------------ allocate the element */      
ele->e.ale3 = (ALE3*)CALLOC(1,sizeof(ALE3));
if (ele->e.ale3==NULL) dserror("Allocation of element ALE failed");
/*----------------------------------- read stuff needed for ALE element */
/*---------------------------------------------- read the element nodes */
frchk("HEX8",&ierr);
if (ierr==1)
{
   ele->distyp = hex8;
   ele->numnp=8;
   ele->lm = (INT*)CALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX8",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("HEX20",&ierr);
if (ierr==1)
{
   ele->distyp = hex20;
   ele->numnp=20;
   ele->lm = (INT*)CALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("HEX20",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
frchk("TET4",&ierr);
if (ierr==1)
{
   ele->distyp = tet4;
   ele->numnp=4;
   ele->lm = (INT*)CALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET4",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
   /*-------------------------- rearrange element node numbers for tet4 */
   lmtmp=ele->lm[0];
   ele->lm[0]=ele->lm[1];
   ele->lm[1]=lmtmp;
}
frchk("TET10",&ierr);
if (ierr==1)
{
   ele->distyp = tet10;
   ele->numnp=10;
   ele->lm = (INT*)CALLOC(ele->numnp,sizeof(INT));
   if (ele->lm==NULL) dserror("Allocation of lm in ELEMENT failed");
   frint_n("TET10",&(ele->lm[0]),ele->numnp,&ierr);
   if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
}
/*------------------------------------------ reduce node numbers by one */
for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
/*-------------------------------------------- read the material number */
frint("MAT",&(ele->mat),&ierr);
if (ierr!=1) dserror("Reading of ALE element failed");
/*-------------------------------------------- read the gaussian points */
if (ele->numnp==8 || ele->numnp==20 || ele->numnp==27)
{
   frint_n("GP",&(ele->e.ale3->nGP[0]),3,&ierr);
   if (ierr!=1) dserror("Reading of ALE3 element failed\n");
}   
/*----------------------- read gaussian points for tetrahedral elements */
if (ele->numnp==4 || ele->numnp==10)
{
   frint("GP_TET",&(ele->e.ale3->nGP[0]),&ierr);
   if (ierr!=1) dserror("Reading of ALE3 element failed\n");
}
/*-------------------------- read gaussian points for triangle elements */
frint("JAC",&(ele->e.ale3->jacobi),&ierr);
if (ierr!=1) dserror("Reading of ALE element failed");

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale3inp */




/*!----------------------------------------------------------------------
\brief  copy element info from fluidfield to ale field 

<pre>                                                              mg 04/01 

</pre>
\param *fluidfield  FIELD  (i/o)   
\param *alefield    FIELD  (i/o)   

\warning There is nothing special to this routine
\return void                                               
\sa caling: ---; called by: inpfield()

*----------------------------------------------------------------------*/
void fluid_to_ale(const FIELD *fluidfield, const FIELD *alefield)
{
INT  ierr=0;
INT  i,j,k;
ELEMENT *fluid_ele;
ELEMENT *ale_ele;
#ifdef DEBUG 
dstrc_enter("fluid_to_ale");
#endif
for (i=0; i<fluidfield->dis[0].numele; i++)
{
    fluid_ele = &(fluidfield->dis[0].element[i]);
#ifdef D_FLUID3 
          if (fluid_ele->eltyp==el_fluid3)
          {
             fluid_ele->e.f3->my_ale  = NULL;
          }
#endif
          if (fluid_ele->eltyp==el_fluid2)
          {
             dserror("Fluid2 not yet implemented");
          }
}
for (i=0; i<alefield->dis[0].numele; i++)
{
   ale_ele = &(alefield->dis[0].element[i]);
   ale_ele->e.ale3->my_fluid = NULL;
}
/*----------------------------------------------------------------------*/
for (i=0; i<fluidfield->dis[0].numele; i++)/*------ loop fluid elements */
{
    fluid_ele = &(fluidfield->dis[0].element[i]);
    if (fluid_ele->e.f3->is_ale!=1) continue;
    for (j=0; j<alefield->dis[0].numele; j++)/*------ loop ale elements */
    {
       ale_ele = &(alefield->dis[0].element[j]);
       if (ale_ele->e.ale3->my_fluid!=NULL) continue;
       /*----------------------- check the geometry of the two elements */
       find_compatible_ele(fluid_ele,ale_ele,&ierr);
       if (ierr==0) continue;
/*--------------------------------------- set ale flag in fluid element */
#ifdef D_FLUID3 
       if (fluid_ele->eltyp==el_fluid3)
       {
          fluid_ele->e.f3->my_ale  = ale_ele;
          ale_ele->e.ale3->my_fluid = fluid_ele;
       }
#endif
#ifdef D_FLUID1 
       if (fluid_ele->eltyp==el_fluid2)
       {
          dserror("Fluid2 not yet implemented");
       }
#endif
          break;
    } /* end of loop over alefield */
} /* end of loop over fluidfield */

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_to_ale */




/*!----------------------------------------------------------------------
\brief finds the compatible ale element 

<pre>                                                              mg 04/01 

</pre>
\param *ele1    ELEMENT  (i/o)   
\param *ele2    ELEMENT  (i/o)   
\param *ierr    INT      (i/o)   

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: fluid_to_ale()

*----------------------------------------------------------------------*/
void find_compatible_ele(const ELEMENT *ele1, const ELEMENT *ele2, INT *ierr)
{
INT  i,j,k;
DOUBLE tol=1.0E-08;
DOUBLE x1,y1,z1, x2,y2,z2;
DOUBLE x,y,z;
#ifdef DEBUG 
dstrc_enter("find_compatible_ele");
#endif
/*----------------------------------------------------------------------*/
x1=0.0;y1=0.0;z1=0.0;x2=0.0;y2=0.0;z2=0.0;
for (i=0; i<ele1->numnp; i++)
{
   x1 += ele1->node[i]->x[0];
   y1 += ele1->node[i]->x[1];
   z1 += ele1->node[i]->x[2];
}
x1 /= (ele1->numnp);
y1 /= (ele1->numnp);
z1 /= (ele1->numnp);
for (i=0; i<ele1->numnp; i++)
{
   x2 += ele2->node[i]->x[0];
   y2 += ele2->node[i]->x[1];
   z2 += ele2->node[i]->x[2];
}
x2 /= (ele1->numnp);
y2 /= (ele1->numnp);
z2 /= (ele1->numnp);
x = FABS(x1-x2);
y = FABS(x1-x2);
z = FABS(x1-x2);
if (x<=tol && y<=tol && z<=tol) *ierr=1;
else *ierr=0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of find_compatible_ele */
#endif
/*! @} (documentation module close)*/
