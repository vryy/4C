/*!----------------------------------------------------------------------
\file
\brief ls2_inpele.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*!----------------------------------------------------------------------
\brief input of one ls2 element

<pre>                                                            irhan 05/04
input of one ls2 element
</pre>

*----------------------------------------------------------------------*/
void ls2inp(
  ELEMENT *ele
  )
{
  INT     i;
  INT     ierr=0;
#ifdef DEBUG 
  dstrc_enter("ls2inp");
#endif
/*----------------------------------------------------------------------*/      

  /* allocate the element */      
  ele->e.ls2 = (LS2*)CCACALLOC(1,sizeof(LS2));
  /* READ element parameters */
  /* connectivity */
  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    ele->numnp=4;
    ele->distyp = quad4;
    ele->e.ls2->ntyp = 1;    
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("QUAD4",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    ele->numnp=3;
    ele->distyp = tri3;
    ele->e.ls2->ntyp = 2;    
    ele->lm = (INT*)CCACALLOC(ele->numnp,sizeof(INT));
    frint_n("TRI3",&(ele->lm[0]),ele->numnp,&ierr);
    if (ierr!=1) dserror("Reading of ELEMENT Topology failed");
  }  
  /* reduce node numbers */
  for (i=0; i<ele->numnp; i++) (ele->lm[i])--;
  /* read the material number */
  frint("MAT",&(ele->mat),&ierr);
  if (ierr!=1) dserror("Reading of LS2 element failed");
  /* read number of Gauss points */
  if (ele->numnp==4)
  {
    frint_n("GP",&(ele->e.ls2->nGP[0]),2,&ierr);
    if (ierr!=1) dserror("Reading of LS2 element failed: integration\n");  
  }
  if (ele->numnp==3)
  {
    frint("GP_TRI",&(ele->e.ls2->nGP[0]),&ierr);
    if (ierr!=1) dserror("Reading of LS2 element failed: integration\n");  
  }  
  /* initialize some parameters */
  ele->e.ls2->is_elcut = 0;
  ele->e.ls2->is_elsearch = 0;
  ele->e.ls2->nlayer = 0;
  ele->e.ls2->intcnt   = 0;
  ele->e.ls2->prncnt   = 0;
  ele->e.ls2->rstcnt   = 0;
  ele->e.ls2->is_sedge_set = 0;
  /*ele->e.ls2->intdata = (LS_INT_DATA*)CCACALLOC(2,sizeof(LS_INT_DATA));
  ele->e.ls2->polydata = (LS_POLY_DATA*)CCACALLOC(2,sizeof(LS_POLY_DATA));
  */
/*----------------------------------------------------------------------*/        
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of ls2inp */



/*!----------------------------------------------------------------------
\brief construct connectivity betweeen level set and fluid problems

<pre>                                                            irhan 05/04
construct connectivity betweeen level set and fluid problems
</pre>

*----------------------------------------------------------------------*/
void fluid_to_ls(
  const FIELD *fluidfield,
  const FIELD *lsfield
  )
{
  INT          ii,i,j;
  INT          ierr=0;
  ELEMENT     *fluid_ele;
  ELEMENT     *ls_ele;
  
#ifdef DEBUG 
  dstrc_enter("fluid_to_ls");
#endif
/*----------------------------------------------------------------------*/          
  
  /* initialize */
  for (ii=0; ii<fluidfield->ndis; ii++)
  {
    for (i=0; i<fluidfield->dis[ii].numele; i++)
    {
      fluid_ele = &(fluidfield->dis[0].element[i]);
      if (fluid_ele->eltyp==el_fluid2_xfem)
      {
        fluid_ele->e.f2->my_ls  = NULL;  
      }	
    }
    for (i=0; i<lsfield->dis[0].numele; i++)
    {
      ls_ele = &(lsfield->dis[0].element[i]);
      if (ls_ele->eltyp==el_ls2)
      {
        ls_ele->e.ls2->my_fluid = NULL;
      }   
    }    
  }
  /* fill */
  for (ii=0; ii<fluidfield->ndis; ii++)
  {
    for (i=0; i<fluidfield->dis[ii].numele; i++)
    {
      fluid_ele = &(fluidfield->dis[0].element[i]);
      if (fluid_ele->eltyp==el_fluid2_xfem)
      {
        if (fluid_ele->e.f2->is_ls!=1) continue;
        for (j=0; j<lsfield->dis[0].numele; j++)
        {
          ls_ele = &(lsfield->dis[0].element[j]);
          if (ls_ele->e.ls2->my_fluid!=NULL) continue;
          /* check geometry */
          ls_find_compatible_ele(fluid_ele,ls_ele,&ierr);
          if (ierr==0) continue;
          /* set */       
          fluid_ele->e.f2->my_ls  = ls_ele;
          ls_ele->e.ls2->my_fluid = fluid_ele;
          break;          
        }
      }
    }    
  }
  
/*----------------------------------------------------------------------*/            
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of fluid_to_ls */



/*!----------------------------------------------------------------------
\brief find compatible element 

<pre>                                                            irhan 05/04
find compatible element 
</pre>

*----------------------------------------------------------------------*/
void ls_find_compatible_ele(
  const ELEMENT *ele1,
  const ELEMENT *ele2,
  INT *ierr
  )
{
  INT  i;
  DOUBLE tol = EPS8;
  DOUBLE x1,y1,z1,x2,y2,z2;
  DOUBLE x,y,z;
  
#ifdef DEBUG 
  dstrc_enter("ls_find_compatible_ele");
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
  for (i=0; i<ele2->numnp; i++)
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
} /* end of ls_find_compatible_ele */
/*! @} (documentation module close)*/
#endif
