/*!----------------------------------------------------------------------
\file
\brief service routines for multi-level fluid2 element 

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
static INT PREDOF = 2;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!--------------------------------------------------------------------- 
\brief set all arrays for large-scale element for multi-level fluid2

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
void f2_lsset(FLUID_DYN_CALC  *dynvar, 
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
GSURF *actgsurf;    /* actual g-surface                                 */

#ifdef DEBUG 
dstrc_enter("f2_lsset");
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
/*--------------------------------------- set element pressures at (n) */   
    epren[i]   =actnode->sol_increment.a.da[1][PREDOF];      
  } /* end of loop over nodes of large-scale element */  
} /* endif (dynvar->nif!=0) */		       

/*------------------------------------------------ check for dead load */
actgsurf = ele->g.gsurf;
if (actgsurf->neum!=NULL)
{
   actmat=ele->mat-1;
   dens = mat[actmat].m.fluid->density;
   for (i=0;i<2;i++)     
   {
      if (actgsurf->neum->neum_onoff.a.iv[i]==0)
      {
         edeadn[i]  = ZERO;
	 edead[i] = ZERO;
      }
      if (actgsurf->neum->neum_type==neum_dead  &&
          actgsurf->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadn[i]  = actgsurf->neum->neum_val.a.dv[i]*dens;
	 edead[i] = actgsurf->neum->neum_val.a.dv[i]*dens;
	 (*hasext)++;
      }
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_lsset */

/*!--------------------------------------------------------------------- 
\brief set all arrays for (sub-)submesh element for multi-level fluid2

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
\param    iele     INT             (i)    actual (sub-)submesh ele. num.
\param    flag     INT             (i)    flag: submesh or sub-submesh?
\return void                                                                       

------------------------------------------------------------------------*/
void f2_smset(FLUID_ML_SMESH  *smesh, 
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
dstrc_enter("f2_smset");
#endif

for (i=0; i<smesh->numen; i++)/* loop nodes of (sub-)submesh element */
{
  smitope[i]    = smesh->ien.a.ia[iele][i];
  smlme[i]      = smesh->id.a.iv[smitope[i]];
  smxyzep[0][i] = smesh->xyzpd.a.da[0][smitope[i]];
  smxyzep[1][i] = smesh->xyzpd.a.da[1][smitope[i]];
  switch (flag)
  {
/*-------------------------------------------------------- sub-submesh */
  case 1:
    smxyze[0][i] = ele->e.f2->xyzssm.a.da[0][smitope[i]];
    smxyze[1][i] = ele->e.f2->xyzssm.a.da[1][smitope[i]];
  break;
/*------------------------------------------------------------ submesh */
  default:
    smxyze[0][i] = ele->e.f2->xyzsm.a.da[0][smitope[i]];
    smxyze[1][i] = ele->e.f2->xyzsm.a.da[1][smitope[i]];
  break;
  }
}  
/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_smset */

/*!--------------------------------------------------------------------- 
\brief set all bubble functions for submesh element for fluid2

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
void f2_bubset(FLUID_DYN_ML    *mlvar,  
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
dstrc_enter("f2_bubset");
#endif

/*---- calculate velocity bubble functions at nodes of submesh element */ 	  
for (i=0;i<mlvar->nvbub;i++)
{
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) evbub[j][i]=ZERO;
    else 
    {
/*---------------------------------------------------- time step (n+1) */ 	        
      if (flag==0) evbub[j][i]= ele->e.f2->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */ 	        
      else         evbub[j][i]= ele->e.f2->solsmn.a.da[km][i];
    }  
  }      
}      

/*---- calculate pressure bubble functions at nodes of submesh element */ 	        
for (i=mlvar->nvbub;i<mlvar->nelbub-2;i++)
{
  ipbub=i-mlvar->nvbub; 
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) epbub[j][ipbub]=ZERO;
    else 
    {
/*---------------------------------------------------- time step (n+1) */ 	        
      if (flag==0) epbub[j][ipbub]= ele->e.f2->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */ 	        
      else         epbub[j][ipbub]= ele->e.f2->solsmn.a.da[km][i];
    }  
  }      
}      

/*--------- calculate rhs bubble functions at nodes of submesh element */ 	        
for (i=mlvar->nelbub-2;i<mlvar->nelbub;i++)
{
  ifbub=i-(mlvar->nvbub+mlvar->npbub); 
  for (j=0;j<submesh->numen;j++)
  {
    km=smlme[j];
    if (km==-1) efbub[j][ifbub]=ZERO;
    else 
    {
/*---------------------------------------------------- time step (n+1) */ 	        
      if (flag==0) efbub[j][ifbub]= ele->e.f2->solsm.a.da[km][i];
/*------------------------------------------------------ time step (n) */ 	        
      else         efbub[j][ifbub]= ele->e.f2->solsmn.a.da[km][i];
    }  
  }      
}      

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_bubset */

/*!--------------------------------------------------------------------- 
\brief set all arrays for integration of sub-submesh element for fluid2

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
void f2_ssset(FLUID_ML_SMESH  *ssmesh, 
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
dstrc_enter("f2_ssset");
#endif

for (i=0; i<ssmesh->numen; i++)/* loop nodes of sub-submesh element */
{
/*----------------------------------------- sub-submesh element arrays */ 	        
  ssitope[i]   = ssmesh->ien.a.ia[iele][i];
  sslme[i]     = ssmesh->id.a.iv[ssitope[i]];
  ssxyze[0][i] = ele->e.f2->xyzssm.a.da[0][ssitope[i]];
  ssxyze[1][i] = ele->e.f2->xyzssm.a.da[1][ssitope[i]];
  
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
} /* end of f2_ssset */

/*!--------------------------------------------------------------------- 
\brief copy submesh solution to element array for fluid2

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
void f2_smcopy(DOUBLE  **smrhs,
               ELEMENT  *ele,
               INT       numeq,
               INT       numrhs)      
{
INT  i,j;

#ifdef DEBUG 
dstrc_enter("f2_smcopy");
#endif

for (i=0;i<numrhs;i++)
{
  for (j=0;j<numeq;j++)
  {
    ele->e.f2->solsm.a.da[j][i] = smrhs[j][i];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_smcopy */ 
   
/*!--------------------------------------------------------------------- 
\brief copy submesh element solution at (n+1) to (n) for fluid2

<pre>                                                       gravem 07/03

In this routine, the submesh element solution array at time step (n+1)
is copied to the respective array at time step (n). 

</pre>
\param   *ele      ELEMENT	   (i/o)  actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void                                                                       

------------------------------------------------------------------------*/
void f2_smcopy2(ELEMENT  *ele,
                INT       numeq,
                INT       numrhs)      
{
INT  i,j;

#ifdef DEBUG 
dstrc_enter("f2_smcopy2");
#endif

for (i=0;i<numrhs;i++)
{
  for (j=0;j<numeq;j++)
  {
    ele->e.f2->solsmn.a.da[j][i] = ele->e.f2->solsm.a.da[j][i];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_smcopy2 */ 
   
/*!--------------------------------------------------------------------- 
\brief routine to calculate small-scale pressure at int. p. for fluid2

<pre>                                                       gravem 07/03
				      
</pre>
\param   *smpreint  DOUBLE     (o)   small-scale pressure at int. point
\param  **pbubint   DOUBLE     (i)   pressure bubble functions at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    int        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f2_smprei(DOUBLE  *smpreint,     
               DOUBLE **pbubint,    
	       DOUBLE  *epre,     
	       INT      iel) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f2_smprei");
#endif

for (i=0;i<2;i++) /* loop spatial directions i */
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
} /* end of f2_smprei */

/*!--------------------------------------------------------------------- 
\brief routine to calculate s-s pressure deriv. at int. p. for fluid2

<pre>                                                       gravem 07/03
				      
In this routine, the derivatives of the small-scale pressure w.r.t x/y 
are calculated.
				      
</pre>
\param  **smpderxy  DOUBLE     (o)   s-s pressure deriv. at int. point
\param ***pbubderxy DOUBLE     (i)   pre. bubble fun. deriv. at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    int        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f2_smpder(DOUBLE  **smpderxy,     
               DOUBLE ***pbubderxy,    
	       DOUBLE   *epre,     
	       INT       iel) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f2_smpder");
#endif

for (i=0;i<2;i++) /* loop spatial directions i */
{
   smpderxy[0][i]=ZERO;
   smpderxy[1][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpderxy[0][i] += pbubderxy[i][0][j]*epre[j];
      smpderxy[1][i] += pbubderxy[i][1][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_smpder */

/*!--------------------------------------------------------------------- 
\brief routine to calc. 2nd s-s pressure deriv. at int. p. for fluid2

<pre>                                                       gravem 07/03
				      
In this routine, the 2nd derivatives of the small-scale pressure w.r.t 
x/y are calculated.
				      
</pre>
\param  **smpderxy2  DOUBLE     (o)   2nd s-s pre. deriv. at int. point
\param ***pbubderxy2 DOUBLE     (i)   2nd pre. bub. fun. deriv. at i. p.
\param   *epre       DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	     INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f2_smpder2(DOUBLE  **smpderxy2,    
                DOUBLE ***pbubderxy2,    
	        DOUBLE   *epre,      
	        INT       iel) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f2_smpder2");
#endif

for (i=0;i<3;i++) /* loop spatial directions i */
{
   smpderxy2[0][i]=ZERO;
   smpderxy2[1][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the l-s element */
   {
      smpderxy2[0][i] += pbubderxy2[i][0][j]*epre[j];
      smpderxy2[1][i] += pbubderxy2[i][1][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_smpder2 */

/*!--------------------------------------------------------------------- 
\brief print out bubble functions for display with CGS

<pre>                                                       gravem 07/03
				      
In this routine, the necessary information for displaying bubble 
functions with CGS is written out in several files 'mesh.*".
				      
</pre>
\param   *submesh    FLUID_ML_SMESH (i)   
\param  **smrhs      DOUBLE         (i)   bubble function solution array
\param    iel	     numrhs         (i)   number of bubble functions 
\return void                                                                       

------------------------------------------------------------------------*/
void f2_cgsbub(FLUID_ML_SMESH  *submesh,
	       DOUBLE	      **smrhs,
               INT              numrhs)      
{
INT           i,j,k,count,top;
DOUBLE        val,null;
FILE         *meshcoo,*meshbou,*meshtop; /* file identifications        */
FILE         *meshtem,*meshval;          /* file identifications        */

#ifdef DEBUG 
dstrc_enter("f2_cgsbub");
#endif

/*------------------------ write submesh coordinates in file 'mesh.coo' */
meshcoo = fopen("cgs_f/mesh.coo","w");
fprintf(meshcoo,"%2d \n",submesh->numnp);
for (i=0;i<submesh->numnp;i++)
{
  k=i+1;
  null=ZERO;
  fprintf(meshcoo,"%2d %4.3f %4.3f %4.3f \n",k,submesh->xyzpd.a.da[0][i],
                                               submesh->xyzpd.a.da[1][i],null);
}  
fclose(meshcoo);

/*-------------------------------------- count number of boundary nodes */
count=0;
for (i=0;i<submesh->numnp;i++)
{
  if(submesh->id.a.iv[i]==-1) count++;
}
  
/*--------------------- write submesh boundary nodes in file 'mesh.bou' */
meshbou = fopen("cgs_f/mesh.bou","w");
fprintf(meshbou,"%2d \n",count);
for (i=0;i<submesh->numnp;i++)
{
  k=i+1;
  if(submesh->id.a.iv[i]==-1) fprintf(meshbou,"%2d 1 1 0 \n",k);
}  
fclose(meshbou);

/*--------------------------- write submesh topology in file 'mesh.top' */
meshtop = fopen("cgs_f/mesh.top","w");
fprintf(meshtop,"%2d %2d \n",submesh->numele,submesh->numen);
for (i=0;i<submesh->numele;i++)
{
  k=i+1;
  fprintf(meshtop,"%2d",k);
  for (j=0;j<submesh->numen;j++)
  {
    top=submesh->ien.a.ia[i][j]+1;
    fprintf(meshtop,"%3d",top);
  }				      
  fprintf(meshtop,"\n");
}  
fclose(meshcoo);

/*---------------------- write solution information in file 'mesh.tem' */
meshtem = fopen("cgs_f/mesh.tem","w");
/*     number of nodes per element */
fprintf(meshtem,"%1d \n",submesh->numen);
/*     flag: 1=values per element, 2=? ,3=values per node  */      
fprintf(meshtem,"3 \n");
/*     number of different value groups */
fprintf(meshtem,"%2d \n",numrhs);
/*     value groups */
fprintf(meshtem,"VELBUB1 \n");
fprintf(meshtem,"VELBUB2 \n");
fprintf(meshtem,"VELBUB3 \n");
fprintf(meshtem,"VELBUB4 \n");
fprintf(meshtem,"PREBUB11 \n");
fprintf(meshtem,"PREBUB12 \n");
fprintf(meshtem,"PREBUB21 \n");
fprintf(meshtem,"PREBUB22 \n");
fprintf(meshtem,"PREBUB31 \n");
fprintf(meshtem,"PREBUB32 \n");
fprintf(meshtem,"PREBUB41 \n");
fprintf(meshtem,"PREBUB42 \n");
fprintf(meshtem,"FBUB1 \n");
fprintf(meshtem,"FBUB2 \n");
/*     number of timesteps */
fprintf(meshtem,"1 \n");
/*     number of nodes and elements, respectively */
fprintf(meshtem,"%2d \n",submesh->numnp);
fclose(meshtem);

/*--------------------------- write solution values in file 'mesh.tem' */
meshval = fopen("cgs_f/mesh.val","w");
for (i=0;i<numrhs;i++)
{
  count=0;
  for (j=0;j<submesh->numnp;j++)
  {
    if(submesh->id.a.iv[j]==-1) val=ZERO;
    else                        
    {
      val=smrhs[count][i];;
      count++;
    }  
    fprintf(meshbou,"%25.20f \n",val);
  }  
}  
fclose(meshval);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_cgsbub */ 

void f2_mlpermestif(DOUBLE         **estif,   
		    DOUBLE	 **emass, 
		    DOUBLE	 **tmp,   
		    INT		   iel,   
		    FLUID_DYN_CALC  *dynvar) 
{
INT    i,j,icol,irow;     /* simply some counters                       */
INT    nvdof;             /* number of vel dofs                         */
INT    npdof;             /* number of pre dofs                         */
INT    totdof;            /* total number of dofs                       */
DOUBLE thsl,thpl;         /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG 
dstrc_enter("f2_mlpermestif");
#endif

/*----------------------------------------------------- set some values */
nvdof  = NUM_F2_VELDOF*iel;
npdof  = iel;
totdof = (NUM_F2_VELDOF+1)*iel;
thsl   = dynvar->thsl;
thpl   = dynvar->thpl;

/*--------------------------------------------- copy estif to tmp-array *
                and multitply stiffness matrix with respective THETA*DT */
for (i=0;i<totdof;i++)
{
/*--------------------------------------------------------- Kvv and Kpv */
   for (j=0;j<nvdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl;
   } /* end of loop over j */
/*--------------------------------------------------------- Kvp and Kpp */
   for (j=nvdof;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thpl;
   } /* end of loop over j */
} /* end of loop over i */
/*------------------------------- add mass matrix for instationary case */
if (dynvar->nis==0)
{
   for (i=0;i<totdof;i++)
   {
      for (j=0;j<nvdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (dynvar->nis==0) */

/*------------------------------------------------------- rearrange Kvv */
irow = 0;
for (i=0;i<nvdof;i+=2)
{   
   icol = 0;
   for (j=0;j<nvdof;j+=2) 
   {
      estif[irow][icol]     = tmp[i][j];
      estif[irow+1][icol]   = tmp[i+1][j];
      estif[irow][icol+1]   = tmp[i][j+1];
      estif[irow+1][icol+1] = tmp[i+1][j+1];
      icol += 3;
   } /* end of loop over j */
   irow += 3;   
} /* end of loop over i */

/*------------------------------------------------------- rearrange Kvp */
irow = 0;
for (i=0;i<nvdof;i+=2)
{
   icol = 2;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow+1][icol] = tmp[i+1][j];
      icol += 3;     
   } /* end of loop over j */
   irow += 3;   
} /* end of loop over i */

/*------------------------------------------------------- rearrange Kpv */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   icol = 0;
   for (j=0;j<nvdof;j+=2)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow][icol+1] = tmp[i][j+1];
      icol += 3;     
   } /* end of loop over j */
   irow += 3;   
} /* end of loop over i */

/*------------------------------------------------------- rearrange Kpp */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   icol = 2;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol] = tmp[i][j];
      icol += 3;
   } /* end of loop over j */
   irow += 3;
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_mlpermestif */

void f2_mljaco(DOUBLE     *funct,    
               DOUBLE    **deriv,   
               DOUBLE    **xjm,     
               DOUBLE     *det,     
               ELEMENT    *ele,     
               INT         iel)

{
INT k;

#ifdef DEBUG 
dstrc_enter("f2_mljaco");
#endif	 

/*---------------------------------- determine jacobian at point r,s ---*/       
xjm[0][0] = 0.0 ;
xjm[0][1] = 0.0 ;
xjm[1][0] = 0.0 ;
xjm[1][1] = 0.0 ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     xjm[0][0] += deriv[0][k] * ele->node[k]->x[0] ;
     xjm[0][1] += deriv[0][k] * ele->node[k]->x[1] ;
     xjm[1][0] += deriv[1][k] * ele->node[k]->x[0] ;
     xjm[1][1] += deriv[1][k] * ele->node[k]->x[1] ;
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/        
*det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];
    
if(*det<0.0)
{   
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   dserror("NEGATIVE JACOBIAN DETERMINANT\n");
}
   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_mljaco */

void f2_mljaco3(DOUBLE    **xyze,
                DOUBLE	 *funct,    
                DOUBLE	**deriv,   
                DOUBLE	**xjm,     
                DOUBLE	 *det,  	
                INT	  iel,	
                ELEMENT	 *ele)
{
INT k;        /* just a counter                                         */

#ifdef DEBUG 
dstrc_enter("f2_mljaco3");
#endif	 

/*---------------------------------- determine jacobian at point r,s ---*/       
xjm[0][0] = ZERO ;
xjm[0][1] = ZERO ;
xjm[1][0] = ZERO ;
xjm[1][1] = ZERO ;

for (k=0; k<iel; k++) /* loop all nodes of the submesh element */
{
     xjm[0][0] += deriv[0][k] * xyze[0][k] ;
     xjm[0][1] += deriv[0][k] * xyze[1][k] ;
     xjm[1][0] += deriv[1][k] * xyze[0][k] ;
     xjm[1][1] += deriv[1][k] * xyze[1][k] ;
} /* end loop over all nodes of the submesh element */

/*------------------------------------------ determinant of jacobian ---*/        
*det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];
      
/*------------------------------------- check if determinant is zero ---*/        
if(*det<0.0)
{   
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   dserror("NEGATIVE JACOBIAN DETERMINANT\n");
}
   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_mljaco3 */

void f2_mlgcoor(DOUBLE     *funct,      
                ELEMENT    *ele,      
	        INT         iel,      
	        DOUBLE     *gcoor)
{
INT i;

#ifdef DEBUG 
dstrc_enter("f2_mlgcoor");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;

for(i=0;i<iel;i++) /* loop all nodes of the element */
{
   gcoor[0] += funct[i] * ele->node[i]->x[0];
   gcoor[1] += funct[i] * ele->node[i]->x[1];
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_mlgcoor */

void f2_mlgcoor2(DOUBLE     *funct,      
                 DOUBLE    **xyze,      
	         INT         iel,      
	         DOUBLE     *gcoor)
{
INT i;

#ifdef DEBUG 
dstrc_enter("f2_mlgcoor2");
#endif

gcoor[0]=ZERO;
gcoor[1]=ZERO;

for(i=0;i<iel;i++) /* loop all nodes of the element */
{
   gcoor[0] += funct[i] * xyze[0][i];
   gcoor[1] += funct[i] * xyze[1][i];
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_mlgcoor2 */

void f2_mlgder2(ELEMENT     *ele,     
	        DOUBLE	 **xjm,      
                DOUBLE	 **bi,     
	        DOUBLE	 **xder2,  
	        DOUBLE	 **derxy,  
	        DOUBLE	 **derxy2, 
                DOUBLE	 **deriv2, 
	        INT	   iel)
{
INT i;
DOUBLE x00,x01,x02,x10,x11,x12,x20,x21,x22;
DOUBLE det,dum;
DOUBLE r0,r1,r2;

#ifdef DEBUG 
dstrc_enter("f2_mlgder2");
#endif

/*--------------------------- calculate elements of jacobian_bar matrix */
x00 = xjm[0][0]*xjm[0][0];
x01 = xjm[0][1]*xjm[0][1];
x02 = TWO*xjm[0][0]*xjm[0][1];
x10 = xjm[1][0]*xjm[1][0];
x11 = xjm[1][1]*xjm[1][1];
x12 = TWO*xjm[1][0]*xjm[1][1];
x20 = xjm[0][0]*xjm[1][0];
x21 = xjm[0][1]*xjm[1][1];
x22 = xjm[0][0]*xjm[1][1] + xjm[0][1]*xjm[1][0];

/*-------------------------------------- inverse of jacobian_bar matrix */
det =   x00*x11*x22 + x01*x12*x20 + x10*x21*x02 \
      - x20*x11*x02 - x00*x12*x21 - x01*x10*x22 ;
dum = ONE/det;         
bi[0][0] =   dum*(x11*x22 - x12*x21);
bi[1][0] =  -dum*(x10*x22 - x20*x12);
bi[2][0] =   dum*(x10*x21 - x20*x11);
bi[0][1] =  -dum*(x01*x22 - x21*x02);
bi[1][1] =   dum*(x00*x22 - x02*x20);
bi[2][1] =  -dum*(x00*x21 - x20*x01);
bi[0][2] =   dum*(x01*x12 - x11*x02);
bi[1][2] =  -dum*(x00*x12 - x10*x02);
bi[2][2] =   dum*(x00*x11 - x01*x10);

/*---------------------------------------------------------- initialise*/
for (i=0;i<3;i++)
{
   xder2[i][0]=ZERO;
   xder2[i][1]=ZERO;
} /* end loop over i */
for (i=0;i<iel;i++)
{
   derxy2[0][i]=ZERO;
   derxy2[1][i]=ZERO;
   derxy2[2][i]=ZERO;
} /* end loop over iel */

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++) /* loop all nodes of the element */
{
   xder2[0][0] += deriv2[0][i] * ele->node[i]->x[0];
   xder2[1][0] += deriv2[1][i] * ele->node[i]->x[0];
   xder2[2][0] += deriv2[2][i] * ele->node[i]->x[0];
   xder2[0][1] += deriv2[0][i] * ele->node[i]->x[1];
   xder2[1][1] += deriv2[1][i] * ele->node[i]->x[1];
   xder2[2][1] += deriv2[2][i] * ele->node[i]->x[1];
} /* end loop over iel */

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++) /* loop all nodes of the element */
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i];
   
   derxy2[0][i] += bi[0][0]*r0 + bi[0][1]*r1 + bi[0][2]*r2;
   derxy2[1][i] += bi[1][0]*r0 + bi[1][1]*r1 + bi[1][2]*r2;
   derxy2[2][i] += bi[2][0]*r0 + bi[2][1]*r1 + bi[2][2]*r2;
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_mlgder2 */

void f2_mlcogder2(DOUBLE     **xyze,     
	          DOUBLE     **xjm,      
                  DOUBLE     **bi,     
	          DOUBLE     **xder2,  
	          DOUBLE     **derxy,  
	          DOUBLE     **derxy2, 
                  DOUBLE     **deriv2, 
	          INT          iel)
{
INT i;
DOUBLE x00,x01,x02,x10,x11,x12,x20,x21,x22;
DOUBLE det,dum;
DOUBLE r0,r1,r2;

#ifdef DEBUG 
dstrc_enter("f2_mlcogder2");
#endif

/*--------------------------- calculate elements of jacobian_bar matrix */
x00 = xjm[0][0]*xjm[0][0];
x01 = xjm[0][1]*xjm[0][1];
x02 = TWO*xjm[0][0]*xjm[0][1];
x10 = xjm[1][0]*xjm[1][0];
x11 = xjm[1][1]*xjm[1][1];
x12 = TWO*xjm[1][0]*xjm[1][1];
x20 = xjm[0][0]*xjm[1][0];
x21 = xjm[0][1]*xjm[1][1];
x22 = xjm[0][0]*xjm[1][1] + xjm[0][1]*xjm[1][0];

/*-------------------------------------- inverse of jacobian_bar matrix */
det =   x00*x11*x22 + x01*x12*x20 + x10*x21*x02 \
      - x20*x11*x02 - x00*x12*x21 - x01*x10*x22 ;
dum = ONE/det;         
bi[0][0] =   dum*(x11*x22 - x12*x21);
bi[1][0] =  -dum*(x10*x22 - x20*x12);
bi[2][0] =   dum*(x10*x21 - x20*x11);
bi[0][1] =  -dum*(x01*x22 - x21*x02);
bi[1][1] =   dum*(x00*x22 - x02*x20);
bi[2][1] =  -dum*(x00*x21 - x20*x01);
bi[0][2] =   dum*(x01*x12 - x11*x02);
bi[1][2] =  -dum*(x00*x12 - x10*x02);
bi[2][2] =   dum*(x00*x11 - x01*x10);

/*---------------------------------------------------------- initialise*/
for (i=0;i<3;i++)
{
   xder2[i][0]=ZERO;
   xder2[i][1]=ZERO;
} /* end loop over i */
for (i=0;i<iel;i++)
{
   derxy2[0][i]=ZERO;
   derxy2[1][i]=ZERO;
   derxy2[2][i]=ZERO;
} /* end loop over iel */

/*----------------------- determine 2nd derivatives of coord.-functions */
for (i=0;i<iel;i++) /* loop all nodes of the element */
{
   xder2[0][0] += deriv2[0][i] * xyze[0][i];
   xder2[1][0] += deriv2[1][i] * xyze[0][i];
   xder2[2][0] += deriv2[2][i] * xyze[0][i];
   xder2[0][1] += deriv2[0][i] * xyze[1][i];
   xder2[1][1] += deriv2[1][i] * xyze[1][i];
   xder2[2][1] += deriv2[2][i] * xyze[1][i];
} /* end loop over iel */

/*--------------------------------- calculate second global derivatives */
for (i=0;i<iel;i++) /* loop all nodes of the element */
{
   r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i];
   r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i];
   r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i];
   
   derxy2[0][i] += bi[0][0]*r0 + bi[0][1]*r1 + bi[0][2]*r2;
   derxy2[1][i] += bi[1][0]*r0 + bi[1][1]*r1 + bi[1][2]*r2;
   derxy2[2][i] += bi[2][0]*r0 + bi[2][1]*r1 + bi[2][2]*r2;
} /* end of loop over iel */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_mlcogder2 */

#endif
