#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#define PREDOF (2)
#define ZERO (0.0)
/*----------------------------------------------------------------------*
 |  set all arrays for element calculation                              |
 |  get all the element velocities / pressure at different times        |
 |  NOTE: in contradiction to the old programm the kinematic pressure   |
 |        is stored in the solution history; so one can avoid the       |
 |        transformation in every time step                             |
 |                                                           genk 04/02 |
 *----------------------------------------------------------------------*/
void f2_calset( 
                FLUID_DYN_CALC  *dynvar, 
	        ELEMENT         *ele,
                double         **eveln,
	        double         **evelng,
	        double          *epren
	      )
{
int i,j,irow;
NODE *actnode;

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
   for(i=0;i<ele->numnp;i++)
   {
      actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */      
      evelng[0][i]=actnode->sol_increment.a.da[3][0];
      evelng[1][i]=actnode->sol_increment.a.da[3][1]; 
/*-------------------------------------- set supported pressures (n+1) */   
/*      if(actnode->dof[PREDOF]>numeq)
         epres[i]=actnode->sol_increment.a.da[3][PREDOF]; */
   }
   
}
else  /* -> semi-implicit time integration method 
       | for semi-impl. methods one needs extra suup. velocities-------*/
{   
   if(dynvar->itwost==0)   /* -> semi-implicit one-step ---------------*/
   {
      for(i=0;i<ele->numnp;i++)
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
      } 
   }
   else  /* -> semi-implicit two-step ---------------------------------*/
   {  
      for(i=0;i<ele->numnp;i++)
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
      }
   }      
}

if(dynvar->nif!=0) /* -> computation if time forces "on" --------------
                      -> velocities and pressure at (n) are needed ----*/
{
      for(i=0;i<ele->numnp;i++)
      {
         actnode=ele->node[i];
/*------------------------------------- set element velocities at (n) */	 
         eveln[0][i]=actnode->sol_increment.a.da[1][0];
	 eveln[1][i]=actnode->sol_increment.a.da[1][1];    
/*------------------------------------------------- set pressures (n) */   
         epren[i]   =actnode->sol_increment.a.da[1][PREDOF];      
      }   
}		       

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_calset */


/*---------------------------------------------------------------------*
 | routine to calculate velocities at integration point                |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_veli(
             double  *velint,   /* velocities at integration point */ 
             double  *funct,    /* shape functions                 */
	     double **evel,     /* velocites at element nodes      */
	     int      iel       /* number of nodes in this element */
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_veli");
#endif

for (i=0;i<2;i++)
{
   velint[i]=ZERO;
   for (j=0;j<iel;j++)
   {
      velint[i] += funct[j]*evel[i][j];
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_veli */

/*---------------------------------------------------------------------*
 | routine to calculate pressure at integration point                |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_prei(
             double  *preint,   /* pressure at integration point */ 
             double  *funct,    /* shape functions                 */
	     double  *epre,     /* pressure at element nodes      */
	     int      iel       /* number of nodes in this element */
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

/*---------------------------------------------------------------------*
 | routine to calculate velocity derivatives at integration point      |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_vder(
             double **vderxy,   /* velocity derivativs             */ 
             double **derxy,    /* globael derivatives             */
	     double **evel,     /* velocites at element nodes      */
	     int      iel       /* number of nodes in this element */
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_vder");
#endif

for (i=0;i<2;i++)
{
   vderxy[0][i]=ZERO;
   vderxy[1][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy[0][i] += derxy[i][j]*evel[0][j];
      vderxy[1][i] += derxy[i][j]*evel[1][j];
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_vder */

/*---------------------------------------------------------------------*
 | routine to 2nd calculate velocity derivatives at integration point  |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_vder2(
             double **vderxy2,   /* velocity derivativs             */ 
             double **derxy2,    /* globael derivatives             */
	     double **evel,     /* velocites at element nodes      */
	     int      iel       /* number of nodes in this element */
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_vder2");
#endif

for (i=0;i<2;i++)
{
   vderxy2[0][i]=ZERO;
   vderxy2[1][i]=ZERO;
   for (j=0;j<iel;j++)
   {
      vderxy2[0][i] += derxy2[i][j]*evel[0][j];
      vderxy2[1][i] += derxy2[i][j]*evel[1][j];
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_vder2 */


/*---------------------------------------------------------------------*
 | routine to calculate pressure derivatives at integration point      |  
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_pder(
             double  *pderxy,   /* pressure derivativs             */ 
             double **derxy,    /* globael derivatives             */
	     double  *epre,     /* pressure at element nodes      */
	     int      iel       /* number of nodes in this element */
	    ) 
{
int     i,j;

#ifdef DEBUG 
dstrc_enter("f2_pder");
#endif

for (i=0;i<2;i++)
{
   pderxy[i] =  ZERO;
   for (j=0;j<iel;j++)
   {
      pderxy[i] += derxy[i][j]*epre[j];
   }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_pder */

/*---------------------------------------------------------------------*
 | routine to calculate convective velocities at integration point     |
 | u * grad(u)                                                         |
 | e.g. 2D: COVx = Ux*Ux,x + Uy*Ux,y                                   |
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_covi(
             double **vderxy,   /* velocity derivativs                 */ 
             double  *velint,   /* velocity at integration point       */
	     double  *covint    /* convective velocity at int point    */ 
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
   }
}


/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_covi */

/*---------------------------------------------------------------------*
 | routine to rearrange the entries of the element force vector        |
 | this is necessary since we would like to use the existing assembly  |
 | routines for the RHS                                                |
 | hence a splitting of vel- and pre dof in the element force vector   |
 | is not possible any more!!!!                                        |
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_permeforce(
		   double    *eforce,
		   double   **tmp,
		   int        iel		   		   
	          ) 
{
int i,irow;
int nvdof;      /* number of vel dofs */
int totdof;     /* total number of dofs */

#ifdef DEBUG 
dstrc_enter("f2_permeforce");
#endif

nvdof  = NUM_F2_VELDOF*iel;
totdof = (NUM_F2_VELDOF+1)*iel;

/*------------------------ compute vel-dofs */
irow = 0;
for (i=0;i<nvdof;i+=2)
{   
   tmp[irow][0]   = eforce[i];
   tmp[irow+1][0] = eforce[i+1];
   irow += 3;   
}

/*----------------------- compute pre-dofs */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   tmp[irow][0] = eforce[i];
   irow += 3;
}

/*----------------------- copy back to eforce */
for (i=0;i<totdof;i++)
{
   eforce[i] = tmp[i][0];
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_permeforce */


/*---------------------------------------------------------------------*
 | routine to add galerkin and stabilisation parts of the elment       |
 | stiffness matrix and to rearrange its entries!		       |
 | this is necessary since we would like to use the existing assembly  |
 | routines for the stiffness matrix                                   |
 | hence a splitting of vel- and pre dofs is not possible any more!!!! |  				       |
 |                                                         genk 04/02  |
 *---------------------------------------------------------------------*/
void f2_permestif(                  
		   double         **estif,
		   double         **emass,
		   double         **tmp,
		   int              iel,
		   FLUID_DYN_CALC  *dynvar		   		   
	          ) 
{
int i, j, icol, irow; 
int nvdof;      /* number of vel dofs */
int npdof;      /* number of pre dofs */
int totdof;     /* total number of dofs */
double thsl;

#ifdef DEBUG 
dstrc_enter("f2_permestif");
#endif

nvdof  = NUM_F2_VELDOF*iel;
npdof  = iel;
totdof = (NUM_F2_VELDOF+1)*iel;
thsl   = dynvar->thsl;

/*----------------------------------- copy estif to tmp-array *
                 and mutlitply stiffniss matrix with THETA*DT */
for (i=0;i<totdof;i++)
{
   for (j=0;j<totdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl;
   }
}
/*--------------------- add mass matrix for instationary case */
if (dynvar->nis==0)
{
   for (i=0;i<totdof;i++)
   {
      for (j=0;j<nvdof;j++)
      {
         tmp[i][j] += emass[i][j];
      }
   }
}

/*------------------------ compute Kvv */
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
   }
   irow += 3;   
}

/*----------------------- compute Kvp */
irow = 0;
for (i=0;i<nvdof;i+=2)
{
   icol = 2;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow+1][icol] = tmp[i+1][j];
      icol += 3;     
   }
   irow += 3;   
}

/*----------------------- compute Kpv */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   icol = 0;
   for (j=0;j<nvdof;j+=2)
   {
      estif[irow][icol]   = tmp[i][j];
      estif[irow][icol+1] = tmp[i][j+1];
      icol += 3;     
   }
   irow += 3;   
}

/*----------------------- compute Kpp */
irow = 2;
for (i=nvdof;i<totdof;i++)
{
   icol = 2;
   for (j=nvdof;j<totdof;j++)
   {
      estif[irow][icol] = tmp[i][j];
      icol += 3;
   }
   irow += 3;
}


/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of f2_permestif */




