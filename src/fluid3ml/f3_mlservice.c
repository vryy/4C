/*!----------------------------------------------------------------------
\file
\brief service routines for multilevel fluid3 element 

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3ml_prototypes.h"
#include "../fluid3/fluid3.h"
static INT PREDOF = 3;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!--------------------------------------------------------------------- 
\brief set all arrays for large-scale element for multi-level fluid3

<pre>                                                       gravem 07/03

In this routine, the element velocities, pressure and external loads 
at time steps (n) and (n+1) are set. 

</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param  **eveln    DOUBLE	   (o)    ele vels at time step n
\param  **evel     DOUBLE	   (o)    ele vels at time step n+1
\param   *epren    DOUBLE	   (o)    ele pres at time step n
\param   *epre     DOUBLE	   (o)    ele pres at time step n+1
\param   *edeadn   DOUBLE          (o)    ele dead load at time step n 
\param   *edead    DOUBLE          (o)    ele dead load at time step n+1 
\param   *hasext   INT             (o)    flag for external loads
\return void                                                                       

------------------------------------------------------------------------*/
void f3_lsset(FLUID_DYN_CALC  *dynvar, 
	      ELEMENT	      *ele,	
              DOUBLE	     **eveln,	 
	      DOUBLE	     **evel, 
	      DOUBLE	      *epren,
	      DOUBLE	      *epre,
	      DOUBLE	      *edeadn,
	      DOUBLE	      *edead,
	      INT	      *hasext)
{
INT i,j,irow;       /* simply some counters                             */
INT    actmat  ;    /* material number of the element                   */
DOUBLE dens;        /* density                                          */
DOUBLE visc;        /* viscosity                                        */
NODE  *actnode;     /* actual node                                      */
GVOL  *actgvol;     /* actual g-volume                                  */

#ifdef DEBUG 
dstrc_enter("f3_lsset");
#endif


/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	 sol_increment[1][i]: solution at (n)                          |
 |	 sol_increment[2][i]: solution at (n+g)                        |
 |	 sol_increment[3][i]: solution at (n+1)                        |
 *---------------------------------------------------------------------*/


for(i=0;i<ele->numnp;i++) /* loop nodes of large-scale element */
{
  actnode=ele->node[i];
/*------------------------------------ set element velocities at (n+1) */      
  evel[0][i]=actnode->sol_increment.a.da[3][0];
  evel[1][i]=actnode->sol_increment.a.da[3][1]; 
  evel[2][i]=actnode->sol_increment.a.da[3][2]; 
/*------------------------------------- set element pressures at (n+1) */   
  epre[i]   =actnode->sol_increment.a.da[3][PREDOF];      
} /* end of loop over nodes of large-scale element */

if(dynvar->nif!=0) /* -> computation if time forces "on" --------------
                      -> velocities and pressure at (n) are needed ----*/
{
  for(i=0;i<ele->numnp;i++) /* loop nodes of large-scale element */
  {
    actnode=ele->node[i];
/*-------------------------------------- set element velocities at (n) */	 
    eveln[0][i]=actnode->sol_increment.a.da[1][0];
    eveln[1][i]=actnode->sol_increment.a.da[1][1];    
    eveln[2][i]=actnode->sol_increment.a.da[1][2];    
/*--------------------------------------- set element pressures at (n) */   
    epren[i]=actnode->sol_increment.a.da[1][PREDOF];      
  } /* end of loop over nodes of large-scale element */  
} /* endif (dynvar->nif!=0) */		       

/*------------------------------------------------ check for dead load */
actgvol = ele->g.gvol;
if (actgvol->neum!=NULL)
{
   actmat=ele->mat-1;
   dens = mat[actmat].m.fluid->density;
   for (i=0;i<3;i++)     
   {
      if (actgvol->neum->neum_onoff.a.iv[i]==0)
      {
         edeadn[i] = ZERO;
	 edead[i]  = ZERO;
      }
      if (actgvol->neum->neum_type==neum_dead  &&
          actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadn[i] = actgvol->neum->neum_val.a.dv[i]*dens;
	 edead[i]  = actgvol->neum->neum_val.a.dv[i]*dens;
	 (*hasext)++;
      }
   }
}

/*------------ check for prescribed pressure gradient for channel flow */
if(dynvar->pregrad!=0)
{
/*--------------------------------------------- mean pressure gradient */
  if (dynvar->pregrad==1)
  {
    actmat=ele->mat-1;
    dens = mat[actmat].m.fluid->density;
    visc = mat[actmat].m.fluid->viscosity*dens;/* here we need dyn. visc.! */
    edeadn[0] = TWO*visc;
    edead[0]  = TWO*visc;
  }  
/*------------------------------------------ dynamic pressure gradient */
  else if (dynvar->pregrad==2)
  {
    edeadn[0] = dynvar->washstr;
    edead[0]  = dynvar->washstr;
  }  
  edeadn[1] = ZERO;
  edead[1]  = ZERO;
  edeadn[2] = ZERO;
  edead[2]  = ZERO;
  (*hasext)++;
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_lsset */

/*!--------------------------------------------------------------------- 
\brief set all arrays for (sub-)submesh element for multi-level fluid3

<pre>                                                       gravem 07/03

In this routine, the (sub-)submesh element coordinates, topology array
and location matrix are set. 

</pre>
\param   *smesh    FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *smlme    INT   	   (o)    (sub-)submesh location matrix
\param   *smitope  INT   	   (o)    (sub-)submesh topology array
\param   *smxyze   DOUBLE   	   (o)    (sub-)submesh coordinates
\param   *smxyzep  DOUBLE   	   (o)    (s)sm parent domain coordinates
\param    iele     INT             (i)    actual submesh element number
\param    flag     INT             (i)    flag: submesh or sub-submesh?
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smset(FLUID_ML_SMESH  *smesh, 
	      ELEMENT	      *ele,	
              INT 	      *smlme,	 
	      INT	      *smitope, 
	      DOUBLE	     **smxyze,
	      DOUBLE	     **smxyzep,
	      INT	       iele,
	      INT	       flag)
{
INT i;       /* simply a counter                                       */

#ifdef DEBUG 
dstrc_enter("f3_smset");
#endif

for (i=0; i<smesh->numen; i++)/* loop nodes of (sub-)submesh element */
{
  smitope[i]    = smesh->ien.a.ia[iele][i];
  smlme[i]      = smesh->id.a.iv[smitope[i]];
  smxyzep[0][i] = smesh->xyzpd.a.da[0][smitope[i]];
  smxyzep[1][i] = smesh->xyzpd.a.da[1][smitope[i]];
  smxyzep[2][i] = smesh->xyzpd.a.da[2][smitope[i]];
  switch (flag)
  {
/*-------------------------------------------------------- sub-submesh */
  case 1:
    smxyze[0][i] = ele->e.f3->xyzssm.a.da[0][smitope[i]];
    smxyze[1][i] = ele->e.f3->xyzssm.a.da[1][smitope[i]];
    smxyze[2][i] = ele->e.f3->xyzssm.a.da[2][smitope[i]];
  break;
/*------------------------------------------------------------ submesh */
  default:
    smxyze[0][i] = ele->e.f3->xyzsm.a.da[0][smitope[i]];
    smxyze[1][i] = ele->e.f3->xyzsm.a.da[1][smitope[i]];
    smxyze[2][i] = ele->e.f3->xyzsm.a.da[2][smitope[i]];
  break;
  }
}  
/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_smset */

/*!--------------------------------------------------------------------- 
\brief set all bubble functions for submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh element velocity, pressure and rhs bubble 
functions are set. 

</pre>
\param   *mlvar    FLUID_DYN_ML    (i)
\param   *submesh  FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *smlme    INT   	   (i)    submesh location matrix
\param  **evbub    DOUBLE   	   (o)    submesh element velocity bubble
\param  **epbub    DOUBLE   	   (o)    submesh element pressure bubble
\param  **efbub    DOUBLE   	   (o)    submesh element rhs bubble
\param    flag     INT             (i)    flag: time step (n) or (n+1)?
\return void                                                                       

------------------------------------------------------------------------*/
void f3_bubset(FLUID_DYN_ML    *mlvar,  
               FLUID_ML_SMESH  *submesh, 
	       ELEMENT	       *ele,	
               INT	       *smlme,
	       DOUBLE         **evbub,
	       DOUBLE         **epbub,
	       DOUBLE         **efbub,
	       INT	        flag)
{
INT i,j,ipbub,ifbub;/* simply some counters                            */
INT km;             /* value in location matrix                        */

#ifdef DEBUG 
dstrc_enter("f3_bubset");
#endif

/*---- calculate velocity bubble functions at nodes of submesh element */ 	  
for (i=0;i<mlvar->nvbub;i++)
{
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) evbub[j][i]=ZERO;
    else if (km>=submesh->numeq) continue;
    else 
    {
/*---------------------------------------------------- time step (n+1) */ 	        
      if (flag==0) evbub[j][i]= ele->e.f3->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */ 	        
      else         evbub[j][i]= ele->e.f3->solsmn.a.da[km][i];
    }  
  }      
}      
      
/*---- calculate pressure bubble functions at nodes of submesh element */ 	        
for (i=mlvar->nvbub;i<mlvar->nelbub-3;i++)
{
  ipbub=i-mlvar->nvbub; 
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) epbub[j][ipbub]=ZERO;
    else if (km>submesh->numeq) continue;
    else 
    {
/*---------------------------------------------------- time step (n+1) */ 	        
      if (flag==0) epbub[j][ipbub]= ele->e.f3->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */ 	        
      else         epbub[j][ipbub]= ele->e.f3->solsmn.a.da[km][i];
    }  
  }      
}      

/*--------- calculate rhs bubble functions at nodes of submesh element */ 	        
for (i=mlvar->nelbub-3;i<mlvar->nelbub;i++)
{
  ifbub=i-(mlvar->nvbub+mlvar->npbub); 
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) efbub[j][ifbub]=ZERO;
    else if (km>submesh->numeq) continue;
    else 
    {
/*---------------------------------------------------- time step (n+1) */ 	        
      if (flag==0) efbub[j][ifbub]= ele->e.f3->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */ 	        
      else         efbub[j][ifbub]= ele->e.f3->solsmn.a.da[km][i];
    }  
  }      
}      

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_bubset */

/*!--------------------------------------------------------------------- 
\brief set all arrays for integration of sub-submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the sub-submesh element arrays for the elementwise
integration of a normalized bubble function are set. 

</pre>
\param   *ssmesh   FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *sslme    INT   	   (o)    sub-submesh location matrix
\param   *ssitope  INT   	   (o)    sub-submesh topology array
\param   *ssxyze   DOUBLE   	   (o)    sub-submesh coordinates
\param   *ebub     DOUBLE   	   (o)    sub-submesh element bubble
\param    iele     INT             (i)    actual sub-submesh ele. number
\return void                                                                       

------------------------------------------------------------------------*/
void f3_ssset(FLUID_ML_SMESH  *ssmesh, 
	      ELEMENT	      *ele,	
              INT 	      *sslme,	 
	      INT	      *ssitope, 
	      DOUBLE	     **ssxyze,
	      DOUBLE          *ebub,
	      INT	       iele)
{
INT i;              /* simply a counter                                */
INT km;             /* value in location matrix                        */

#ifdef DEBUG 
dstrc_enter("f3_ssset");
#endif

for (i=0; i<ssmesh->numen; i++)/* loop nodes of sub-submesh element */
{
/*----------------------------------------- sub-submesh element arrays */ 	        
  ssitope[i]   = ssmesh->ien.a.ia[iele][i];
  sslme[i]     = ssmesh->id.a.iv[ssitope[i]];
  ssxyze[0][i] = ele->e.f3->xyzssm.a.da[0][ssitope[i]];
  ssxyze[1][i] = ele->e.f3->xyzssm.a.da[1][ssitope[i]];
  ssxyze[2][i] = ele->e.f3->xyzssm.a.da[2][ssitope[i]];
  
/*- calculate normalized bubble funct. at nodes of sub-submesh element */ 	        
  km=sslme[i];
  if (km==-1) ebub[i]=ZERO;
  else  ebub[i]= ssmesh->rhs.a.dv[km];      
}      

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_ssset */

/*!--------------------------------------------------------------------- 
\brief copy submesh solution to element array for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh solution is copied to the respective
element array. 

</pre>
\param   *smrhs    DOUBLE          (i)    submesh solution array
\param   *ele      ELEMENT	   (o)    actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smcopy(DOUBLE  **smrhs,
               ELEMENT  *ele,
               INT       numeq,
               INT       numrhs)      
{
INT  i,j;

#ifdef DEBUG 
dstrc_enter("f3_smcopy");
#endif

for (i=0;i<numrhs;i++)
{
  for (j=0;j<numeq;j++)
  {
    ele->e.f3->solsm.a.da[j][i] = smrhs[j][i];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_smcopy */ 
   
/*!--------------------------------------------------------------------- 
\brief copy submesh element solution at (n+1) to (n) for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh element solution array at time step (n+1)
is copied to the respective array at time step (n). 

</pre>
\param   *ele      ELEMENT	   (i/o)  actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smcopy2(ELEMENT  *ele,
                INT       numeq,
                INT       numrhs)      
{
INT  i,j;

#ifdef DEBUG 
dstrc_enter("f3_smcopy2");
#endif

for (i=0;i<numrhs;i++)
{
  for (j=0;j<numeq;j++)
  {
    ele->e.f3->solsmn.a.da[j][i] = ele->e.f3->solsm.a.da[j][i];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f3_smcopy2 */ 
   
/*!--------------------------------------------------------------------- 
\brief routine to calculate small-scale pressure at int. p. for fluid3

<pre>                                                       gravem 07/03
				      
</pre>
\param   *smpreint  DOUBLE     (o)   small-scale pressure at int. point
\param  **pbubint   DOUBLE     (i)   pressure bubble functions at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smprei(DOUBLE  *smpreint,     
               DOUBLE **pbubint,    
	       DOUBLE  *epre,     
	       INT      iel) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_smprei");
#endif

for (i=0;i<3;i++) /* loop spatial directions i */
{
   smpreint[i]=ZERO; 
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpreint[i] += pbubint[i][j]*epre[j];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_smprei */

/*!--------------------------------------------------------------------- 
\brief routine to calculate s-s pressure deriv. at int. p. for fluid3

<pre>                                                       gravem 07/03
				      
In this routine, the derivatives of the small-scale pressure w.r.t x/y/z 
are calculated.
				      
</pre>
\param   *smpderxy  DOUBLE     (o)   s-s pressure deriv. at int. point
\param  **pbubderxy DOUBLE     (i)   pre. bubble fun. deriv. at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smpder(DOUBLE  **smpderxy,     
               DOUBLE ***pbubderxy,    
	       DOUBLE   *epre,     
	       INT       iel) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_smpder");
#endif

for (i=0;i<3;i++) /* loop spatial directions i */
{
   smpderxy[0][i]=ZERO;
   smpderxy[1][i]=ZERO;
   smpderxy[2][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpderxy[0][i] += pbubderxy[i][0][j]*epre[j];
      smpderxy[1][i] += pbubderxy[i][1][j]*epre[j];
      smpderxy[2][i] += pbubderxy[i][2][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_smpder */

/*!--------------------------------------------------------------------- 
\brief routine to calc. 2nd s-s pressure deriv. at int. p. for fluid3

<pre>                                                       gravem 07/03
				      
In this routine, the 2nd derivatives of the small-scale pressure w.r.t 
x/y/z are calculated.
				      
</pre>
\param  **smpderxy2  DOUBLE     (o)   2nd s-s pre. deriv. at int. point
\param ***pbubderxy2 DOUBLE     (i)   2nd pre. bub. fun. deriv. at i. p.
\param   *epre       DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	     INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smpder2(DOUBLE  **smpderxy2,    
                DOUBLE ***pbubderxy2,    
	        DOUBLE   *epre,      
	        INT       iel) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_smpder2");
#endif

for (i=0;i<6;i++) /* loop spatial directions i */
{
   smpderxy2[0][i]=ZERO;
   smpderxy2[1][i]=ZERO;
   smpderxy2[2][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpderxy2[0][i] += pbubderxy2[i][0][j]*epre[j];
      smpderxy2[1][i] += pbubderxy2[i][1][j]*epre[j];
      smpderxy2[2][i] += pbubderxy2[i][2][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_smpder2 */


#endif
