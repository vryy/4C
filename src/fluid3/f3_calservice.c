/*!----------------------------------------------------------------------
\file
\brief service routines for fluid3 element 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>


/*----------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

static INT PREDOF = 3;

/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         genk 05/02

get the element velocities and the pressure at different times 

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the	 
      transformation in every time step 			 
				      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param  **eveln    DOUBLE	   (o)    ele vels at time n
\param  **evelng   DOUBLE	   (o)    ele vels at time n+g
\param   *epren    DOUBLE	   (o)    ele pres at time n
\param   *edeadn   DOUBLE          (o)    ele dead load at n (selfweight)
\param   *edeadng  DOUBLE          (o)    ele dead load at n+g (selfweight)
\param   *hasext   INT             (o)    flag for external loads
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,
                DOUBLE         **eveln,
	        DOUBLE         **evelng,
	        DOUBLE          *epren,
		DOUBLE          *edeadn,
		DOUBLE          *edeadng,
		INT             *hasext		
	      )
{
INT i,j,irow;
INT    actmat  ;    /* material number of the element                   */
DOUBLE dens;        /* density                                          */
NODE  *actnode;     /* actual node                                      */
GVOL  *actgvol;

#ifdef DEBUG 
dstrc_enter("f3_calset");
#endif


/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	 sol_increment[1][i]: solution at (n)                          |
 |	 sol_increment[2][i]: solution at (n+g)                        |
 |	 sol_increment[3][i]: solution at (n+1)                        |
 *---------------------------------------------------------------------*/


for(i=0;i<ele->numnp;i++) /* loop nodes */
{
   actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */      
   evelng[0][i]=actnode->sol_increment.a.da[3][0];
   evelng[1][i]=actnode->sol_increment.a.da[3][1]; 
   evelng[2][i]=actnode->sol_increment.a.da[3][2];
} /* end of loop over nodes */
   

if(dynvar->nif!=0) /* -> computation if time forces "on" --------------
                      -> velocities and pressure at (n) are needed ----*/
{
      for(i=0;i<ele->numnp;i++) /* loop nodes */
      {
         actnode=ele->node[i];
/*------------------------------------- set element velocities at (n) */	 
         eveln[0][i]=actnode->sol_increment.a.da[1][0];
	 eveln[1][i]=actnode->sol_increment.a.da[1][1];    
	 eveln[2][i]=actnode->sol_increment.a.da[1][2]; 
/*------------------------------------------------- set pressures (n) */   
         epren[i]   =actnode->sol_increment.a.da[1][PREDOF];      
      } /* end of loop over nodes */
} /* endif (dynvar->nif!=0) */

if(dynvar->nim!=0) /* -> computation of mass rhs "on" ------------------
                      -> vel(n)+a*acc(n) are needed --------------------*/
/* NOTE: if there is no classic time rhs (as described in WAW) the array
         eveln is misused and does NOT contain the velocity at time (n)
	 but rather a linear combination of old velocities and 
	 accelerations depending upon the time integration scheme!!!!!*/
{
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
         eveln[0][i] = actnode->sol_increment.a.da[2][0];
	 eveln[1][i] = actnode->sol_increment.a.da[2][1];
	 eveln[2][i] = actnode->sol_increment.a.da[2][2];
      } /* end of loop over nodes of element */  
} /* endif (dynvar->nim!=0) */		       
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
         edeadn[i]  = ZERO;
	 edeadng[i] = ZERO;
      }
      if (actgvol->neum->neum_type==neum_dead  &&
          actgvol->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadn[i]  = actgvol->neum->neum_val.a.dv[i]*dens;
	 edeadng[i] = actgvol->neum->neum_val.a.dv[i]*dens;
	 (*hasext)++;
      }
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_calset */

/*!--------------------------------------------------------------------- 
\brief routine to calculate velocities at integration point

<pre>                                                         genk 05/02
				      
</pre>
\param   *velint   DOUBLE        (o)   velocities at integration point
\param   *funct    DOUBLE        (i)   shape functions
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_veli(
             DOUBLE  *velint,    
             DOUBLE  *funct,    
	     DOUBLE **evel,     
	     INT      iel       
	    ) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_veli");
#endif

for (i=0;i<3;i++)
{
   velint[i]=ZERO;
   for (j=0;j<iel;j++)
   {
      velint[i] += funct[j]*evel[i][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_veli */

/*!--------------------------------------------------------------------- 
\brief routine to calculate pressure at integration point

<pre>                                                         genk 05/02
				      
</pre>
\param  *preint    DOUBLE        (o)   pressure at integration point
\param  *funct     DOUBLE        (i)   shape functions
\param  *epre      DOUBLE        (i)   pressure at element nodes
\param   iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_prei(
             DOUBLE  *preint,    
             DOUBLE  *funct,    
	     DOUBLE  *epre,     
	     INT      iel       
	    ) 
{
INT     j;

#ifdef DEBUG 
dstrc_enter("f3_prei");
#endif

*preint = ZERO;
for (j=0;j<iel;j++)
{
   *preint += funct[j] * epre[j];
} /* end of loop over j */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_prei */

/*!--------------------------------------------------------------------- 
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 05/02

In this routine the derivatives of the velocity w.r.t x/y are calculated
vderxy[0][2] = Ux,z  
				      
</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_vder(
             DOUBLE **vderxy,     
             DOUBLE **derxy,    
	     DOUBLE **evel,     
	     INT      iel       
	    ) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_vder");
#endif

for (i=0;i<3;i++)
{
   vderxy[0][i]=ZERO;
   vderxy[1][i]=ZERO;
   vderxy[2][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy[0][i] += derxy[i][j]*evel[0][j];
      vderxy[1][i] += derxy[i][j]*evel[1][j];
      vderxy[2][i] += derxy[i][j]*evel[2][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_vder */

/*!---------------------------------------------------------------------  
\brief routine to calculate 2nd velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the 2nd derivatives of the velocity
w.r.t x/y/z are calculated
   vderxy2[0][0] = Ux,xx 
   vderxy2[0][3] = Ux,xy 
   vderxy2[1][4] = Ux,xz 
   vderxy2[2][5] = Ux,yz 
				      
</pre>
\param  **vderxy2  DOUBLE        (o)   2nd velocity derivativs
\param  **derxy2   DOUBLE        (i)   2nd global derivatives
\param  **evel     DOUBLE        (i)   velocites at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_vder2(
             DOUBLE **vderxy2,   
             DOUBLE **derxy2,   
	     DOUBLE **evel,     
	     INT      iel       
	    ) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_vder2");
#endif

for (i=0;i<6;i++)
{
   vderxy2[0][i]=ZERO;
   vderxy2[1][i]=ZERO;
   vderxy2[2][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy2[0][i] += derxy2[i][j]*evel[0][j];
      vderxy2[1][i] += derxy2[i][j]*evel[1][j];
      vderxy2[2][i] += derxy2[i][j]*evel[2][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_vder2 */

/*!--------------------------------------------------------------------- 
\brief routine to calculate pressure derivatives at integration point

<pre>                                                         genk 05/02

In this routine derivatives of the pressure w.r.t x/y/z are calculated
				      
</pre>
\param   *pderxy   DOUBLE        (o)   pressure derivativs
\param  **derxy    DOUBLE        (i)   globael derivatives
\param   *epre     DOUBLE        (i)   pressure at element nodes
\param    iel	   INT           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_pder(
             DOUBLE  *pderxy,    
             DOUBLE **derxy,    
	     DOUBLE  *epre,     
	     INT      iel       
	    ) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("f3_pder");
#endif

for (i=0;i<3;i++)
{
   pderxy[i] =  ZERO;
   for (j=0;j<iel;j++)
   {
      pderxy[i] += derxy[i][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_pder */

/*!--------------------------------------------------------------------- 
\brief convective velocities 

<pre>                                                         genk 05/02		 

in this routine the convective velocity is calculated at the 
integration point:
 u * grad(u)
 e.g. 3D: COVx = Ux*Ux,x + Uy*Ux,y + Uz*Ux,z

</pre>
\param  **vderxy   DOUBLE        (o)   velocity derivativs
\param   *velint   DOUBLE        (i)   velocity at integration point
\param   *covint   DOUBLE        (i)   convective velocity at INT point
\return void                                                                       

------------------------------------------------------------------------*/
void f3_covi(
             DOUBLE **vderxy,    
             DOUBLE  *velint,   
	     DOUBLE  *covint    
	    ) 
{
INT     i,j;      
#ifdef DEBUG 
dstrc_enter("f3_covi");
#endif

for (i=0;i<3;i++)
{
   covint[i]=ZERO;
   for (j=0;j<3;j++)
   {
      covint[i] +=velint[j]*vderxy[i][j];
   } /* end of loop over j */
} /* end of loop over i */


/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_covi */

/*!--------------------------------------------------------------------- 
\brief permutation of element force vector 

<pre>                                                         genk 05/02		 

routine to rearrange the entries of the element force vector	  
this is necessary since we would like to use the existing assembly
routines for the RHS						  
hence a splitting of vel- and pre dof in the element force vector 
is not possible any more!!!!					  


</pre>
\param   *eforce   DOUBLE        (i/o) element force vector
\param  **tmp      DOUBLE        (i)   working array
\param    iel	   DOUBLE        (i)   number of nodes in this ele
\return void                                                                       

------------------------------------------------------------------------*/
void f3_permeforce(
		   DOUBLE    *eforce,
		   DOUBLE   **tmp,
		   INT        iel		   		   
	          ) 
{
INT i,irow;
INT nvdof;      /* number of vel dofs                                   */
INT totdof;     /* total number of dofs                                 */

#ifdef DEBUG 
dstrc_enter("f3_permeforce");
#endif

nvdof  = NUM_F3_VELDOF*iel;
totdof = (NUM_F3_VELDOF+1)*iel;

/*---------------------------------------------------- compute vel-dofs */
irow = 0;
for (i=0;i<nvdof;i+=3)
{   
   tmp[irow][0]   = eforce[i];
   tmp[irow+1][0] = eforce[i+1];
   tmp[irow+2][0] = eforce[i+2];
   irow += 4;   
} /* end of loop over i */

/*---------------------------------------------------- compute pre-dofs */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   tmp[irow][0] = eforce[i];
   irow += 4;
} /* end of loop over i */

/*------------------------------------------------- copy back to eforce */
for (i=0;i<totdof;i++)
{
   eforce[i] = tmp[i][0];
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_permeforce */

/*!--------------------------------------------------------------------- 
\brief permutation of element stiffness matrix

<pre>                                                         genk 05/02		 

routine to add galerkin and stabilisation parts of the elment	   
stiffness matrix and to rearrange its entries!  		   
this is necessary since we would like to use the existing assembly 
routines for the stiffness matrix				   
hence a splitting of vel- and pre dofs is not possible any more!!!!				  

</pre>
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array		
\param	  iel	  INT		 (i)   number of nodes in ele
\param	 *dynvar  FLUID_DYN_CALC
\return void                                                                       

------------------------------------------------------------------------*/
void f3_permestif(                  
		   DOUBLE         **estif,
		   DOUBLE         **emass,
		   DOUBLE         **tmp,
		   INT              iel,
		   FLUID_DYN_CALC  *dynvar		   		   
	          ) 
{
INT i,j,icol,irow;          /* simply some counters  	        	*/
INT nvdof;                  /* number of vel dofs 			*/
INT npdof;                  /* number of pre dofs 			*/
INT totdof;                 /* total number of dofs			*/
DOUBLE thsl;	            /* factor for LHS (THETA*DT)		*/

#ifdef DEBUG 
dstrc_enter("f3_permestif");
#endif

nvdof  = NUM_F3_VELDOF*iel;
npdof  = iel;
totdof = (NUM_F3_VELDOF+1)*iel;
thsl   = dynvar->thsl;

/*--------------------------------------------- copy estif to tmp-array *
                           and mutlitply stiffniss matrix with THETA*DT */
for (i=0;i<totdof;i++)
{
   for (j=0;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl;
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
} /* endif (dynvar->nis==0) *

/*--------------------------------------------------------- compute Kvv */
irow = 0;
for (i=0;i<nvdof;i+=3)
{   
   icol = 0;
   for (j=0;j<nvdof;j+=3)
   {
      estif[irow][icol]     = tmp[i][j];
      estif[irow+1][icol]   = tmp[i+1][j];
      estif[irow+2][icol]   = tmp[i+2][j];
      estif[irow][icol+1]   = tmp[i][j+1];
      estif[irow+1][icol+1] = tmp[i+1][j+1];
      estif[irow+2][icol+1] = tmp[i+2][j+1];
      estif[irow][icol+2]   = tmp[i][j+2];
      estif[irow+1][icol+2] = tmp[i+1][j+2];
      estif[irow+2][icol+2] = tmp[i+2][j+2];      
      icol += 4;
   } /* end of loop over j */
   irow += 4;   
} /* end of loop over i */

/*--------------------------------------------------------- compute Kvp */
irow = 0;
for (i=0;i<nvdof;i+=3)
{
   icol = 3;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow+1][icol] = tmp[i+1][j];
      estif[irow+2][icol] = tmp[i+2][j];
      icol += 4;     
   } /* end of loop over j */
   irow += 4;   
} /* end of loop over i */

/*--------------------------------------------------------- compute Kpv */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   icol = 0;
   for (j=0;j<nvdof;j+=3)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow][icol+1] = tmp[i][j+1];
      estif[irow][icol+2] = tmp[i][j+2];
      icol += 4;     
   } /* end of loop over j */
   irow += 4;   
} /* end of loop over i */

/*--------------------------------------------------------- compute Kpp */
irow = 3;
for (i=nvdof;i<totdof;i++)
{
   icol = 3;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol] = tmp[i][j];
      icol += 4;
   } /* end of loop over j */
   irow += 4;
} /* end of loop over i */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_permestif */



#endif
