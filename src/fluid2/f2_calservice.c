/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element 

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
static int PREDOF = 2;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*!--------------------------------------------------------------------- 
\brief set all arrays for element calculation

<pre>                                                         genk 04/02

get the element velocities and the pressure at different times 

NOTE: in contradiction to the old programm the kinematic pressure
      is stored in the solution history; so one can avoid the	 
      transformation in every time step 			 
				      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param  **eveln    double	   (o)    ele vels at time n
\param  **evelng   double	   (o)    ele vels at time n+g
\param   *epren    double	   (o)    ele pres at time n
\param   *edeadn   double          (o)    ele dead load at n (selfweight)
\param   *edeadng  double          (o)    ele dead load at n+g (selfweight)
\param   *hasext   int             (o)    flag for external loads
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,     
                double         **eveln,    
	        double         **evelng, 
	        double          *epren,
		double          *edeadn,
		double          *edeadng,
		int             *hasext
	      )
{
int i,j,irow;       /* simply some counters                             */
int    actmat  ;    /* material number of the element                   */
double dens;        /* density                                          */
NODE  *actnode;     /* actual node                                      */
GSURF *actgsurf;

#ifdef DEBUG 
dstrc_enter("f2_calset");
#endif


/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	 sol_increment[1][i]: solution at (n)                          |
 |	 sol_increment[2][i]: solution at (n+g)                        |
 |	 sol_increment[3][i]: solution at (n+1)                        |
 *---------------------------------------------------------------------*/


if(dynvar->isemim==0)  /* -> implicit time integration method ---------*/
{
   for(i=0;i<ele->numnp;i++) /* loop nodes of element */
   {
      actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */      
      evelng[0][i]=actnode->sol_increment.a.da[3][0];
      evelng[1][i]=actnode->sol_increment.a.da[3][1]; 
/*-------------------------------------- set supported pressures (n+1) */   
/*      if(actnode->dof[PREDOF]>numeq)
         epres[i]=actnode->sol_increment.a.da[3][PREDOF]; */
   } /* end of loop over nodes of element */
   
} /* endif (dynvar->isemim==0) */
else  /* -> semi-implicit time integration method 
       | for semi-impl. methods one needs extra suup. velocities-------*/
{   
   if(dynvar->itwost==0)   /* -> semi-implicit one-step ---------------*/
   {
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
/*--------------------- set element velocities (n) and supported (n+1) */ 	 
         evelng[0][i]=actnode->sol_increment.a.da[1][0];
	 evelng[1][i]=actnode->sol_increment.a.da[1][1];
/*	 if (actnode->dof[0]>numeq)
            evels[0][i]=actnode->sol_increment.a.da[3][0];
	 if (actnode->dof[1]>numeq)
            evels[1][i]=actnode->sol_increment.a.da[3][1];      
/*-------------------------------------- set supported pressures (n+1) */   
/*         if(actnode->dof[PREDOF]>numeq)
            epres[i]=actnode->sol_increment.a.da[3][PREDOF];           */
      } /* end of loop over nodes of element */
   } /* endif (dynvar->itwost==0) */
   else  /* -> semi-implicit two-step ---------------------------------*/
   {  
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
/*--------------- set element velocities (n+gamma) and supported (n+1) */	 
         evelng[0][i]=actnode->sol_increment.a.da[2][0];
	 evelng[1][i]=actnode->sol_increment.a.da[2][1];
/*	 if (actnode->dof[0]>numeq)
	    evels[0][i]=actnode->sol_increment.a.da[3][0];
	 if (actnode->dof[1]>numeq)
	    evels[1][i]=actnode->sol_increment.a.da[3][1];      
/*-------------------------------------- set supported pressures (n+1) */   
/*         if(actnode->dof[PREDOF]>numeq)
            epres[i]=actnode->sol_increment.a.da[3][PREDOF];      */
      } /* end of loop over nodes of element */
   } /* endif else */  
} /* endif else */

if(dynvar->nif!=0) /* -> computation if time forces "on" --------------
                      -> velocities and pressure at (n) are needed ----*/
{
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
/*------------------------------------- set element velocities at (n) */	 
         eveln[0][i]=actnode->sol_increment.a.da[1][0];
	 eveln[1][i]=actnode->sol_increment.a.da[1][1];    
/*------------------------------------------------- set pressures (n) */   
         epren[i]   =actnode->sol_increment.a.da[1][PREDOF];      
      } /* end of loop over nodes of element */  
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
	 edeadng[i] = ZERO;
      }
      if (actgsurf->neum->neum_type==neum_dead  &&
          actgsurf->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadn[i]  = actgsurf->neum->neum_val.a.dv[i]*dens;
	 edeadng[i] = actgsurf->neum->neum_val.a.dv[i]*dens;
	 (*hasext)++;
      }
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calset */

/*!--------------------------------------------------------------------- 
\brief routine to calculate velocities at integration point

<pre>                                                         genk 04/02
				      
</pre>
\param   *velint   double        (o)   velocities at integration point
\param   *funct    double        (i)   shape functions
\param  **evel     double        (i)   velocites at element nodes
\param    iel	   int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_veli(
             double  *velint,     
             double  *funct,    
	     double **evel,     
	     int      iel       
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_veli");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   velint[i]=ZERO; 
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      velint[i] += funct[j]*evel[i][j];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_veli */

/*!--------------------------------------------------------------------- 
\brief routine to calculate pressure at integration point

<pre>                                                         genk 04/02
				      
</pre>
\param  *preint    double        (o)   pressure at integration point
\param  *funct     double        (i)   shape functions
\param  *epre      double        (i)   pressure at element nodes
\param   iel	   int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_prei(
             double  *preint,     
             double  *funct,    
	     double  *epre,     
	     int      iel       
	    ) 
{
int     j;

#ifdef DEBUG 
dstrc_enter("f2_prei");
#endif

*preint = ZERO;
for (j=0;j<iel;j++)
{
   *preint += funct[j] * epre[j];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_prei */

/*!--------------------------------------------------------------------- 
\brief routine to calculate velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the derivatives of the velocity w.r.t x/y are calculated
				      
</pre>
\param  **vderxy   double        (o)   velocity derivativs
\param  **derxy    double        (i)   globael derivatives
\param  **evel     double        (i)   velocites at element nodes
\param    iel	   int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_vder(
             double **vderxy,     
             double **derxy,    
	     double **evel,     
	     int      iel       
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_vder");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   vderxy[0][i]=ZERO;
   vderxy[1][i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      vderxy[0][i] += derxy[i][j]*evel[0][j];
      vderxy[1][i] += derxy[i][j]*evel[1][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_vder */

/*!---------------------------------------------------------------------  
\brief routine to calculate 2nd velocity derivatives at integration point

<pre>                                                         genk 04/02

In this routine the 2nd derivatives of the velocity
w.r.t x/y are calculated
				      
</pre>
\param  **vderxy2  double        (o)   2nd velocity derivativs
\param  **derxy2   double        (i)   2nd global derivatives
\param  **evel     double        (i)   velocites at element nodes
\param    iel	   int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_vder2(
             double **vderxy2,    
             double **derxy2,    
	     double **evel,      
	     int      iel        
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_vder2");
#endif

for (i=0;i<3;i++)
{
   vderxy2[0][i]=ZERO;
   vderxy2[1][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy2[0][i] += derxy2[i][j]*evel[0][j];
      vderxy2[1][i] += derxy2[i][j]*evel[1][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_vder2 */

/*!--------------------------------------------------------------------- 
\brief routine to calculate pressure derivatives at integration point

<pre>                                                         genk 04/02

In this routine derivatives of the pressure w.r.t x/y are calculated
				      
</pre>
\param   *pderxy   double        (o)   pressure derivativs
\param  **derxy    double        (i)   globael derivatives
\param   *epre     double        (i)   pressure at element nodes
\param    iel	   int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_pder(
             double  *pderxy,     
             double **derxy,    
	     double  *epre,     
	     int      iel        
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_pder");
#endif

for (i=0;i<2;i++) /* loop over directions i */
{
   pderxy[i] =  ZERO;
   for (j=0;j<iel;j++) /* loop over nodes j of the element */
   {
      pderxy[i] += derxy[i][j]*epre[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_pder */

/*!--------------------------------------------------------------------- 
\brief convective velocities 

<pre>                                                         genk 04/02		 

in this routine the convective velocity is calculated at the 
integration point:
 u * grad(u)
 e.g. 2D: COVx = Ux*Ux,x + Uy*Ux,y

</pre>
\param  **vderxy   double        (o)   velocity derivativs
\param   *velint   double        (i)   velocity at integration point
\param   *covint   double        (i)   convective velocity at int point
\return void                                                                       

------------------------------------------------------------------------*/
void f2_covi(
             double **vderxy,    
             double  *velint,   
	     double  *covint    
	    ) 
{
int     i,j;      

#ifdef DEBUG 
dstrc_enter("f2_covi");
#endif

for (i=0;i<2;i++)
{
   covint[i]=ZERO;
   for (j=0;j<2;j++)
   {
      covint[i] +=velint[j]*vderxy[i][j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_covi */

/*!--------------------------------------------------------------------- 
\brief permutation of element force vector 

<pre>                                                         genk 04/02		 

routine to rearrange the entries of the element force vector	  
this is necessary since we would like to use the existing assembly
routines for the RHS						  
hence a splitting of vel- and pre dof in the element force vector 
is not possible any more!!!!					  


</pre>
\param   *eforce   double        (i/o) element force vector
\param  **tmp      double        (i)   working array
\param    iel	   double        (i)   number of nodes in this ele
\return void                                                                       

------------------------------------------------------------------------*/
void f2_permeforce( 
		   double   *eforce,  
		   double  **tmp,    
		   int       iel     
	          ) 
{
int i,irow;
int nvdof;      /* number of vel dofs                                   */
int totdof;     /* total number of dofs                                 */

#ifdef DEBUG 
dstrc_enter("f2_permeforce");
#endif

/*----------------------------------------------------- set some values */
nvdof  = NUM_F2_VELDOF*iel;
totdof = (NUM_F2_VELDOF+1)*iel;

/*-------------------------------------------------- rearrange vel-dofs */
irow = 0;
for (i=0;i<nvdof;i+=2)
{   
   tmp[irow][0]   = eforce[i];
   tmp[irow+1][0] = eforce[i+1];
   irow += 3;   
}

/*-------------------------------------------------- rearrange pre-dofs */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   tmp[irow][0] = eforce[i];
   irow += 3;
}

/*------------------------------------------------- copy back to eforce */
for (i=0;i<totdof;i++)
{
   eforce[i] = tmp[i][0];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_permeforce */

/*!--------------------------------------------------------------------- 
\brief permutation of element stiffness matrix

<pre>                                                         genk 04/02		 

routine to add galerkin and stabilisation parts of the elment	   
stiffness matrix and to rearrange its entries!  		   
this is necessary since we would like to use the existing assembly 
routines for the stiffness matrix				   
hence a splitting of vel- and pre dofs is not possible any more!!!!				  

</pre>
\param  **estif   double	 (i/o) ele stiffnes matrix
\param  **emass   double	 (i)   ele mass matrix
\param  **tmp     double	 (-)   working array		
\param	  iel	  int		 (i)   number of nodes in ele
\param	 *dynvar  FLUID_DYN_CALC
\return void                                                                       

------------------------------------------------------------------------*/
void f2_permestif(                  
		   double         **estif,   
		   double         **emass, 
		   double         **tmp,   
		   int              iel,   
		   FLUID_DYN_CALC  *dynvar		   		    
	          ) 
{
int    i,j,icol,irow;     /* simply some counters                       */
int    nvdof;             /* number of vel dofs                         */
int    npdof;             /* number of pre dofs                         */
int    totdof;            /* total number of dofs                       */
double thsl;              /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG 
dstrc_enter("f2_permestif");
#endif

/*----------------------------------------------------- set some values */
nvdof  = NUM_F2_VELDOF*iel;
npdof  = iel;
totdof = (NUM_F2_VELDOF+1)*iel;
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
} /* end of f2_permestif */



#endif
