/*!----------------------------------------------------------------------
\file
\brief service routines for fluid time algorithms

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par; 
/*!--------------------------------------------------------------------- 
\brief routine to check starting algorithm

<pre>                                                         genk 04/02

this routine conrols the starting algorithms schemes. For the first 'nums'
iteration steps a different time integration scheme may be used then 
during the rest of the simulation.
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC (i)  
\param *nfrastep 	INT           (o)  number of fract. steps

\return void 
\warning this routine is not completely tested yet!

------------------------------------------------------------------------*/
void fluid_startproc(
                          FLUID_DYNAMIC     *fdyn,
		          INT               *nfrastep 
		    )
{
static INT     iop_s;
static INT     iops_s;
static INT     itemax_s;
static DOUBLE  theta_s ;
static DOUBLE  thetas_s;

#ifdef DEBUG 
dstrc_enter("fluid_startproc");
#endif 

if (fdyn->step==1)  /* save parameter from input */
{
   iop_s    = fdyn->iop;
   iops_s   = fdyn->iops;
   itemax_s = fdyn->itemax;
   theta_s  = fdyn->theta;
   thetas_s = fdyn->thetas;
}
if (fdyn->step<=fdyn->nums) /* set parameter for starting phase */
{
   fdyn->iop    = iops_s;
   fdyn->theta  = thetas_s;
   *nfrastep    = 1;
   if (iops_s==5)
      *nfrastep = 3;
   if (iops_s==2 || iops_s==3)
      fdyn->itemax = 1;
   if (fdyn->step==fdyn->nums && iop_s==3)
       dserror ("starting algo for semi-impl. two-step method not implemented yet\n");
       /* set U(n-1) in last step of start-algo */
}
else if (fdyn->step==(fdyn->nums+1)) /* set original parameter */
{
   fdyn->iop    = iop_s;
   fdyn->theta  = theta_s;
   *nfrastep    = 1;
   if (iop_s==5)
      *nfrastep = 3;
   if (iop_s==2 || iop_s==3)
      fdyn->itemax = 1;
   else
      fdyn->itemax = itemax_s;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_startproc*/


/*!---------------------------------------------------------------------  
\brief calculating time independent time integration constants

<pre>                                                       chfoe 09/03

in this routine the constants for the time integration algorithms are 
calculated as far as they are independent of the time
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC   (i)  
\param *dynvar	        FLUID_DYN_CALC  (i/o)  

\return void 
\warning only ONE-STEP-THETA implemented up to now!

------------------------------------------------------------------------*/
void fluid_cons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar
		)
{

#ifdef DEBUG 
dstrc_enter("fluid_cons");
#endif
/*----------------------------------------------------------------------*/
dynvar->dtp = fdyn->dt;
/*----------------------------------------------------- check algorithm */
switch(fdyn->iop)
{
case 1:		/* gen alpha implementation 1 */
   if (fabs(fdyn->theta*TWO*fdyn->alpha_f/fdyn->alpha_m - ONE) > EPS10)
   {
       fdyn->theta = fdyn->alpha_m / fdyn->alpha_f * 0.5;
       printf("\nWarning: Theta, Alpha_m and Alpha_f do not satisfy 2nd order condition.\n");
       printf("         Theta is recalculated.\n");
       printf("\n Theta = Alpha_m / (2 Alpha_f) = %6.4f \n\n", fdyn->theta);
   }
   dynvar->dta = 0.0;
   dynvar->gen_alpha = 1;
   dynvar->omt = ONE-fdyn->theta;
break;
case 4:		/* one step theta */
   dynvar->dta = 0.0;
   dynvar->omt = ONE-fdyn->theta;
   dynvar->gen_alpha = 0;
break;
case 7:		/* 2nd order backward differencing BDF2 */
   dynvar->dta = 0.0;
   dynvar->omt = ONE-fdyn->theta;	
   dynvar->gen_alpha = 0;
   if(fdyn->adaptive)
      if(FABS(fdyn->thetas - 0.5) > EPS11)
         dswarning(1,5);   
break;
default:
   dserror ("constants for time algorithm not implemented yet!\n");
} /* end switch */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_cons*/ 


/*!---------------------------------------------------------------------  
\brief calculating time integration constants

<pre>                                                         genk 03/02

in this routine the constants for the time integration algorithms are 
calculated, here time dependent values only are calculated
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC   (i)  
\param *dynvar	        FLUID_DYN_CALC  (i/o)  

\return void 
\warning only ONE-STEP-THETA implemented up to now!

------------------------------------------------------------------------*/
void fluid_tcons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar
		)
{

#ifdef DEBUG 
dstrc_enter("fluid_tcons");
#endif

/*----------------------------------------------------- check algorithm */
switch (fdyn->iop)
{
case 1:		/* generalised alpha */
   if (fdyn->adaptive)
   {
      if(dynvar->dta == 0.0) dynvar->dta = fdyn->dt;   
   }
   else if (fdyn->adaptive==0)
   {
      dynvar->dta  = fdyn->dt;    
   }
   dynvar->thsl = dynvar->dta * fdyn->theta * fdyn->alpha_f / fdyn->alpha_m;
   dynvar->thpl = dynvar->thsl;
   dynvar->thsr = ZERO;
   dynvar->thpr = dynvar->thsr;
   dynvar->thnr = 0.0;	/* a factor of an addend needed for gen alpha 2 */
   dynvar->alpha = 1.0;	/* a factor needed for gen alpha 2 */
break;
case 4:		/* one step theta */
   if (fdyn->adaptive)
   {
      if(dynvar->dta == 0.0) dynvar->dta = fdyn->dt;   
   }
   else if (fdyn->adaptive==0)
   {
      dynvar->dta  = fdyn->dt;    
   }
   dynvar->thsl = dynvar->dta*fdyn->theta;
   dynvar->thpl = dynvar->thsl;
   dynvar->thsr = (ONE - fdyn->theta)*dynvar->dta;
   dynvar->thpr = dynvar->thsr;
   dynvar->thnr = 0.0;	/* a factor of an addend needed for gen alpha 2 */
   dynvar->alpha = 1.0;	/* a factor needed for gen alpha 2 */
   dynvar->theta = fdyn->theta;
break;
case 7:		/* 2nd order backward differencing BDF2 */
   fdyn->time_rhs = 0;	/* use mass rhs */
   if (fdyn->adaptive)
   {
      if(dynvar->dta == 0.0) dynvar->dta = fdyn->dt;
      dynvar->thsl = (DSQR(dynvar->dta) + dynvar->dta*dynvar->dtp) 
                     / (2.0*dynvar->dta + dynvar->dtp);
      dynvar->thpl = dynvar->thsl;
      dynvar->thsr = 0.0;
      dynvar->thpr = dynvar->thsr; 
      dynvar->thnr = 0.0;	/* a factor of an addend needed for gen alpha 2 */
      dynvar->alpha = 1.0;	/* a factor needed for gen alpha 2 	*/

   }
   else if (fdyn->adaptive==0)
   {
      dynvar->dta  = fdyn->dt;    
      dynvar->thsl = dynvar->dta*2.0/3.0;
      dynvar->thpl = dynvar->thsl;
      dynvar->thsr = 0.0;
      dynvar->thpr = dynvar->thsr;
      dynvar->thnr = 0.0;	/* a factor of an addend needed for gen alpha 2 */
      dynvar->alpha = 1.0;	/* a factor needed for gen alpha 2 	*/
   }
break;
default:
   dserror ("constants for time algorithm not implemented yet!\n");
}

/*----------------------------------------------- treatment of pressure */ 
if (fdyn->iprerhs!=1)
    dserror ("treatment of pressure not implemented yet!\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_tcons*/ 


/*!---------------------------------------------------------------------
\brief setting flags for nonlinear iteration schemes

<pre>                                                         genk 03/02

in this routine the flags for the different nonlinear iteration schemes
are set. Depending on the iteration schemes the evaluation of some 
termes in the LHS or the RHS have to turned on or off:

nik <->  EVALUATION OF LHS-MATRICES (w/o NONLINEAR TERM)    
nic <->  EVALUATION OF NONLINEAR LHS N-CONVECTIVE	    
nir <->  EVALUATION OF NONLINEAR LHS N-REACTION 	    
nie <->  EVALUATE ONLY LHS-TERMS FOR EXPLICIT VELOCITY      
nil <->  EVALUATION OF LUMPED MASS MATRIX (Mvv-lumped)      
nif <->  EVALUATION OF "TIME - RHS" (F-hat)		    
nii <->  EVALUATION OF "ITERATION - RHS"		    
nis <->  STATIONARY CASE (NO TIMEDEPENDENT TERMS)	   
			     
</pre>   
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param *dynvar	   FLUID_DYN_CALC  (i/o)  
\param  itnum      INT  	   (i)     actual number of iterations
\return void 
\warning up to now, only fixed-point like iteration checked!!!

------------------------------------------------------------------------*/
void fluid_icons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar,
		          INT                itnum           
		)
{

#ifdef DEBUG 
dstrc_enter("fluid_icons");
#endif

/*----------------------------------------------------- initialisation */
dynvar->nik=0;
dynvar->nic=0;
dynvar->nir=0;
dynvar->nie=0;
dynvar->nil=0;
dynvar->nif=0;
dynvar->nii=0;
dynvar->nis=0;
dynvar->nim=0;
dynvar->totarea=ZERO;

switch (fdyn->ite)
{
case 0:		/* no iteration */
   dynvar->sigma=ZERO;
   if(fdyn->time_rhs)	/* 'classic' time rhs as in W.A. Wall */
   {
      dynvar->nik=1;       
      dynvar->nic=2; 
      dynvar->nif=3;      /* KCF */
   }
   else if (fdyn->time_rhs == 0)	/* mass formulation of time rhs */
   {
      dynvar->nik=1;       
      dynvar->nic=2; 
      dynvar->nim=1;      /* KC(F) */
   }
   dserror("results with no nonlin. iteration not checked yet!\n");
break;
case 1:		/* fixed point like iteration */
   dynvar->sigma=ZERO;
   if(fdyn->time_rhs)	/* 'classic' time rhs as in W.A. Wall */
   {
      if (itnum>1)
      {
         dynvar->nik=1;  
         dynvar->nic=2;   /* KC */  
      }
      else
      {
         dynvar->nik=1;
         dynvar->nic=2;
         dynvar->nif=3;   /* KCF */
      }
   }
   else if (fdyn->time_rhs == 0)	/* mass formulation of time rhs */
   {
      dynvar->nik=1;  
      dynvar->nic=2;
      dynvar->nim=1;	/* KC(F) */
   }
break;
case 2:		/* Newton iteration */
   dynvar->sigma=ONE;
   if(fdyn->time_rhs)	/* 'classic' time rhs as in W.A. Wall */
   {
      if (itnum>1 || fdyn->iop==7) 	/* no time rhs for BDF2 */
      {
         dynvar->nik=1;
         dynvar->nic=2;
         dynvar->nir=3;
         dynvar->nii=4;  /* KCRI */
      }
      else
      {
         dynvar->nik=1;
         dynvar->nic=2;
         dynvar->nir=3;
         dynvar->nif=4;
         dynvar->nii=5;  /* KCRFI */
      }
   }
   else if (fdyn->time_rhs == 0)	/* mass formulation of time rhs */
   {
      dynvar->nik=1;
      dynvar->nic=2;
      dynvar->nir=3;
      dynvar->nii=4;
      dynvar->nim=1;	/* KCR(F)I */
   }
break;
case 3:		/* fixed point iteration */
   dserror("fixed point iteration not checked yet!!!\n");
   dynvar->sigma=-ONE;
   if(fdyn->time_rhs)	/* 'classic' time rhs as in W.A. Wall */
   {
      if (itnum>1)
      {
         dynvar->nik=1;
         dynvar->nii=2;  /* KI */     
      }
      else
      {
         dynvar->nik=1;
         dynvar->nif=2;      
         dynvar->nii=3;  /* KFI */      
      }
   }
   else if (fdyn->time_rhs == 0)	/* mass formulation of time rhs */
   {
      dynvar->nik=1;
      dynvar->nim=1;      
      dynvar->nii=3;  /* K(F)I */   
   }
break;
default:
   dserror("Unknown nonlinear iteration");
}

/*------------------------------ flags for free surface tension effects */
switch (fdyn->freesurf)
{
case 0: /* do nothing */
break;
case 1: /* explicit: include only for first iteration step in "time rhs" */
   if (itnum>1)
   {
      dynvar->surftens= 0;
      dynvar->fsstnif = 0; 
      dynvar->fsstnii = 0;
   }
   else
   {
      dynvar->surftens=fdyn->surftens;
      dynvar->fsstnif = fdyn->surftens; 
      dynvar->fsstnii = fdyn->surftens;
   }
break;
case 2:
   if (itnum>1)
   {
      dynvar->surftens=fdyn->surftens;
      dynvar->fsstnif = 0; 
      dynvar->fsstnii = fdyn->surftens;
   }
   else
   {
      dynvar->surftens=fdyn->surftens;
      dynvar->fsstnif = fdyn->surftens;
      dynvar->fsstnii = fdyn->surftens;
   }
break;
default:
      dserror("parameter fdyn->freesurf out of range!\n");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_icons*/

/*!--------------------------------------------------------------------- 
\brief initialisation of solution history

<pre>                                                         genk 04/02

in this routine the solution history of the fluid field is initialised.
The time-integration schemes require the solutions at different time-
steps. They are stored in 	   
node->sol_incement: solution history used for calculations
      sol_increment.a.da[0][i]: solution at (n-1)    
      sol_increment.a.da[1][i]: solution at (n)      
      sol_increment.a.da[2][i]: solution at (n+g)    
      sol_increment.a.da[3][i]: solution at (n+1)    
      sol_increment.a.da[4][i]: grid velocity	     
      sol_increment.a.da[5][i]: convective velocity at (n)  
      sol_increment.a.da[6][i]: convective velocity at (n+1)		    
In node->sol one findes the converged solutions of the time-steps, which
are written to the output-file and the pss-file.

If the initial fluid fuild comes from the input-file 'fluid-start-data',
these data are copied to sol_increment.
			     
</pre>   
\param *actfield   FIELD           (i)
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param  numr	   INT             (i/o)   number of rows in sol_incr
\param  str        FLUID_STRESS    (i)     flag for stress calculation
\return void 

------------------------------------------------------------------------*/
void fluid_init(
                          PARTITION	    *actpart,
                          INTRA	            *actintra,
			  FIELD             *actfield,  
                          FLUID_DYNAMIC     *fdyn,
                          CALC_ACTION       *action,
			  CONTAINER         *container,
		          INT                numr,
		          FLUID_STRESS       str	
	       )      
{
INT    i,j,k;           /* simply counters                              */
INT    actmat;          /* number of actual material                    */
INT    numdf;           /* number of dofs in this discretisation        */
INT    numnp_total;     /* total number of nodes in this discretisation */
INT    numele_total;    /* total number of elements in this discr.      */  
INT    numnp=0;         /* number of nodes at actual element            */
INT    found;
int    ldflag;
DOUBLE dens;            /* density                                      */
NODE  *actnode;         /* the actual node                              */
GNODE *actgnode;        /* the actual gnode                             */
ELEMENT *actele;        /* the actual element                           */
GSURF   *actgsurf;
GLINE   *actgline;
DLINE   *actdline;

/*----------------------------------------- variables for solitary wave */
DOUBLE u1,u2,u3,eta,p,c,g,H,d,x,y,t,fac,fac1,sech;

/* variables for beltrami */
DOUBLE    visc,a,x1,x2,x3;

#ifdef DEBUG 
dstrc_enter("fluid_init");
#endif

/*----------------------- set control variables for element evaluation */
fdyn->dynvar.ishape = 1;
fdyn->dynvar.iprerhs= fdyn->iprerhs;
  
numdf = fdyn->numdf;
numele_total = actfield->dis[0].numele;

/*-------------------------------- allocate space for solution history */
for (k=0;k<actfield->ndis;k++)
{
   numnp_total=actfield->dis[k].numnp;
   for (i=0;i<numnp_total;i++)
   {
      actnode=&(actfield->dis[k].node[i]);
      amredef(&(actnode->sol_increment),numr,actnode->numdf,"DA");
      amzero(&(actnode->sol_increment));
   }
}

/*---------------------- redefine solution history for projecton method */
if (fdyn->dyntyp==1)
for (i=0;i<actfield->dis[0].numnp;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   amredef(&(actnode->sol),actnode->sol.fdim,actnode->sol.sdim+1,"DA");
}

/*-------------------- check for stress calculation and allocate arrays */
switch (str)
{
case str_none: /* do nothing */
break;
#ifdef D_FSI
case str_fsicoupling: /* allocate stress field for elements with fsi-coupled nodes */
   for (i=0;i<numele_total;i++)
   {
      actele = &(actfield->dis[0].element[i]);
      numnp=actele->numnp;
      found=0;
      for (j=0;j<numnp;j++)
      {
         actnode=actele->node[j];
	 actgnode = actnode->gnode;
	 if (actgnode->mfcpnode[0]!=NULL) found=1;
      }
      if (found==1)
      {
#ifdef D_FLUID2
         if (numdf==3)    
            amdef("stress_ND",&(actele->e.f2->stress_ND),numnp,3,"DA");
#endif
#ifdef D_FLUID3
         if (numdf==4)
            amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
#endif	 
      }
   }
break;   
#endif
case str_all: /* allocate stress field for all elements */
   for (i=0;i<numele_total;i++)
   {
      actele = &(actfield->dis[0].element[i]);   
#ifdef D_FLUID2
      if (numdf==3)    
         amdef("stress_ND",&(actele->e.f2->stress_ND),numnp,3,"DA");
#endif
#ifdef D_FLUID3
      if (numdf==4)
         amdef("stress_ND",&(actele->e.f3->stress_ND),numnp,6,"DA");
#endif	 
   }
break;   
case str_liftdrag:  
#ifdef D_FLUID2      
   for (i=0;i<numele_total;i++)
   {
      actele = &(actfield->dis[0].element[i]); 
      if (numdf==3)
      {
         actgsurf=actele->g.gsurf;         
         ldflag=0;
	 for (j=0;j<actgsurf->ngline;j++)
         {
            actgline=actgsurf->gline[j];
            actdline=actgline->dline;
            if (actdline==NULL) continue;
            if (actdline->liftdrag==0) continue;
	    ldflag++;
	    break;
         }  
         if (ldflag>0)
            amdef("stress_ND",&(actele->e.f2->stress_ND),actele->numnp,3,"DA");
      }
   }
#endif	
break;   
default:
   dserror("option for 'str' not implemented yet!\n");
}

/*--------------- allocate curvature field for elements at free surface */
if (fdyn->surftens>0)
{
   for (i=0;i<numele_total;i++)
   {
      actele = &(actfield->dis[0].element[i]);      
#ifdef D_FLUID2
      if (numdf==3)
      {             
	 if (actele->e.f2->fs_on==0) continue;
	 {
	    amdef("kappa_ND",&(actele->e.f2->kappa_ND),numnp,2,"DA");
            amzero(&(actele->e.f2->kappa_ND));
         }
      }
#endif
#ifdef D_FLUID3
      if (numdf==4)
         dserror("curvature not implemented for 3D fluid elements!\n");
#endif
   }
}

/*--------- inherit the neuman conditions from design to discretization */
for (k=0; k<actfield->ndis; k++) inherit_design_dis_neum(&(actfield->dis[k]));

/*-------------------------------------------- create the initial field */
if (fdyn->init>=1)
{
   if (fdyn->init==1) /*------------ initial data from fluid_start.data */
      inp_fluid_start_data(actfield,fdyn);

   if (fdyn->init==2) /*-------------------- initial data from pss file */
      restart_read_fluiddyn(fdyn->resstep,fdyn,actfield,actpart,actintra,
                            action,container);
			    
   if (fdyn->init==6) /*--------------------------------- solitary wave */
   {
      actele = &(actfield->dis[0].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[0].numnp;
      /*-------------------------------------------- set some constants */
      dens   = mat[actmat].m.fluid->density;
      g      = -actele->g.gsurf->neum->neum_val.a.dv[1];
      d      = TEN;
      H      = ONE/FIVE*d;
      fac1   = sqrt((THREE*H)/(FOUR*d*d*d));
      c      = sqrt(g*d*(1+H/d));
      t      = ZERO;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
         actnode=&(actfield->dis[0].node[i]);
	 x   = actnode->x[0];
	 y   = actnode->x[1];
	 fac = fac1*(x-c*t);
	 /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 */ 
	 sech = cosh(fac);
	 sech = ONE/sech/sech;
	 eta  = d + H*sech;
	 /*---------------------------------------- modify initial values */ 	 
	 p    = dens*g*(eta-y);
	 u1   = sqrt(g*d)*H/d*sech;
	 u2   = sqrt(THREE*g*d)*(H/d)*sqrt(H/d)*(y/d)*sech*tanh(fac);
	 /*---------------------------------------- write values to nodes */
	 actnode->sol.a.da[0][0] = u1;
	 actnode->sol.a.da[0][1] = u2;
	 actnode->sol.a.da[0][2] = p ;
	 actnode->sol_increment.a.da[1][0] = u1;  
	 actnode->sol_increment.a.da[1][1] = u2;  
	 actnode->sol_increment.a.da[1][2] = p ;  
	 actnode->sol_increment.a.da[3][0] = u1;  
	 actnode->sol_increment.a.da[3][1] = u2;  
	 actnode->sol_increment.a.da[3][2] = p ;

      }    
   }
   if (fdyn->init==7) /*-------------------------- wavebreaking problem */
   {
      actele = &(actfield->dis[0].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[0].numnp;
      /*-------------------------------------------- set some constants */
      dens   = mat[actmat].m.fluid->density;
      g      = -actele->g.gsurf->neum->neum_val.a.dv[1];
      d      = TEN;
      H      = FIVE;
      fac1   = sqrt((THREE*H)/(FOUR*d*d*d));
      c      = sqrt(g*d*(1+H/d));
      t      = ZERO;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
         actnode=&(actfield->dis[0].node[i]);
	 x   = actnode->x[0];
	 y   = actnode->x[1];
	 fac = fac1*(x-c*t);
	 /* sech = sech(fac*xct)**2 = 1/cosh(fac*xct)**2 */ 
	 sech = cosh(fac);
	 sech = ONE/sech/sech;
	 eta  = d + H*sech;
	 /*---------------------------------------- modify initial values */ 	 
	 p    = dens*g*(eta-y);
	 u1   = sqrt(g*d)*H/d*sech;
	 u2   = sqrt(THREE*g*d)*(H/d)*sqrt(H/d)*(y/d)*sech*tanh(fac);
	 /*---------------------------------------- write values to nodes */
	 actnode->sol.a.da[0][0] = u1;
	 actnode->sol.a.da[0][1] = u2;
	 actnode->sol.a.da[0][2] = p ;
	 actnode->sol_increment.a.da[1][0] = u1;  
	 actnode->sol_increment.a.da[1][1] = u2;  
	 actnode->sol_increment.a.da[1][2] = p ;  
	 actnode->sol_increment.a.da[3][0] = u1;  
	 actnode->sol_increment.a.da[3][1] = u2;  
	 actnode->sol_increment.a.da[3][2] = p ;
      }    
   }

   if (fdyn->init==8) /* Beltrami flow */
   {
      actele = &(actfield->dis[0].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[0].numnp;
      /* set some constants */
      visc   = mat[actmat].m.fluid->viscosity;
      a      = PI/4.0;
      d      = PI/2.0;
      t      = 0.0;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
         actnode=&(actfield->dis[0].node[i]);
	 x1   = actnode->x[0];
	 x2   = actnode->x[1];
	 x3   = actnode->x[2];
	 /* calculate initial values */ 	 
	 p    = -a*a/2 * ( exp(2*a*x1) + exp(2*a*x2) + exp(2*a*x3) 
                + 2 * sin(a*x1 + d*x2) * cos(a*x3 + d*x1) * exp(a*(x2+x3))
                + 2 * sin(a*x2 + d*x3) * cos(a*x1 + d*x2) * exp(a*(x3+x1))
                + 2 * sin(a*x3 + d*x1) * cos(a*x2 + d*x3) * exp(a*(x1+x2)) )
                * exp(-2*visc*d*d*t);
	 u1   = -a * ( exp(a*x1) * sin(a*x2 + d*x3) +
                     exp(a*x3) * cos(a*x1 + d*x2) ) * exp(-visc*d*d*t);
	 u2   = -a * ( exp(a*x2) * sin(a*x3 + d*x1) +
                     exp(a*x1) * cos(a*x2 + d*x3) ) * exp(-visc*d*d*t);
	 u3   = -a * ( exp(a*x3) * sin(a*x1 + d*x2) +
                     exp(a*x2) * cos(a*x3 + d*x1) ) * exp(-visc*d*d*t);
	 /* write values to nodes */
	 actnode->sol.a.da[0][0] = u1;
	 actnode->sol.a.da[0][1] = u2;
	 actnode->sol.a.da[0][2] = u3;
	 actnode->sol.a.da[0][3] = p ;
	 actnode->sol_increment.a.da[0][0] = u1;  
	 actnode->sol_increment.a.da[0][1] = u2;  
	 actnode->sol_increment.a.da[0][2] = u3;  
	 actnode->sol_increment.a.da[0][3] = p ;  
	 actnode->sol_increment.a.da[1][0] = u1;  
	 actnode->sol_increment.a.da[1][1] = u2;  
	 actnode->sol_increment.a.da[1][2] = u3;  
	 actnode->sol_increment.a.da[1][3] = p ;  
	 actnode->sol_increment.a.da[3][0] = u1;  
	 actnode->sol_increment.a.da[3][1] = u2;  
	 actnode->sol_increment.a.da[3][2] = u3;  
	 actnode->sol_increment.a.da[3][3] = p ;
      }    
   }

   if (fdyn->init==9) /* Kim-Moin flow */
   {
      actele = &(actfield->dis[0].element[0]);
      actmat = actele->mat-1;
      numnp_total=actfield->dis[0].numnp;
      /* set some constants */
      visc   = mat[actmat].m.fluid->viscosity;
      a      = 2.0;
      t      = 0.0;
      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
         actnode=&(actfield->dis[0].node[i]);
	 x1   = actnode->x[0];
	 x2   = actnode->x[1];
	 /* calculate initial values */ 	 
	 p    = -1.0/4.0 * ( cos(2.0*a*PI*x1) + cos(2.0*a*PI*x2) ) * exp(-4.0*a*a*PI*PI*t*visc);
	 u1   = - cos(a*PI*x1) * sin(a*PI*x2) * exp(-2.0*a*a*PI*PI*t*visc); 
	 u2   = + sin(a*PI*x1) * cos(a*PI*x2) * exp(-2.0*a*a*PI*PI*t*visc);
	 /* write values to nodes */
	 actnode->sol.a.da[0][0] = u1;
	 actnode->sol.a.da[0][1] = u2;
	 actnode->sol.a.da[0][2] = p ;
	 actnode->sol_increment.a.da[0][0] = u1;  
	 actnode->sol_increment.a.da[0][1] = u2;  
	 actnode->sol_increment.a.da[0][2] = p ;  
	 actnode->sol_increment.a.da[1][0] = u1;  
	 actnode->sol_increment.a.da[1][1] = u2;  
	 actnode->sol_increment.a.da[1][2] = p ;  
	 actnode->sol_increment.a.da[3][0] = u1;  
	 actnode->sol_increment.a.da[3][1] = u2;  
	 actnode->sol_increment.a.da[3][2] = p ;
      }    
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_init */ 
   
/*!---------------------------------------------------------------------                                         
\brief storing results in solution history

<pre>                                                         genk 05/02

in this routine the results in the DIST_VECTOR are put to the nodes in
a certain place in ARRAY sol_increment.
Result has to be allreduced and is put to the whole field on each proc.
If necassary the norms for the iteration convergence check of the 
nonlinear iteration scheme are calculated.	   
			     
</pre>   
\param **actfield      FIELD	      (i)    actual field       
\param  *actintra      INTRA	      (i)    actual intra comm.
\param	*sol 	       DIST_VECTOR    (i)    solution vector
\param	 place         INT	      (i)    place in sol_incr.
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*vrat          DOUBLE	      (o)    vel.  conv. ratio
\param	*prat          DOUBLE	      (o)    pre.  conv. ratio
\param  *grat          DOUBLE         (o)    grid  conv. ratio
\param	*fdyn	       FLUID_DYNAMIC	     	
\return void 

------------------------------------------------------------------------*/
void fluid_result_incre(  
                          FIELD             *actfield,    
                          INTRA             *actintra,   
			  DIST_VECTOR       *sol,        
                          INT                place,      
			  SPARSE_ARRAY      *sysarray,      
			  SPARSE_TYP        *sysarray_typ,
			  DOUBLE            *vrat,        
			  DOUBLE            *prat,
                          DOUBLE            *grat,
			  FLUID_DYNAMIC     *fdyn           
		       )
{
INT      i,j;          /* simply some counters                         */
INT      max;       
INT      diff;
INT      dof;          /* actual dof number                            */
INT      predof;       /* number of pressure dof (2=2D; 3=3D)          */
INT      numeq_total;  /* total number of equations                    */
DOUBLE   dvnorm=ZERO;
DOUBLE   dpnorm=ZERO;
DOUBLE   dgnorm=ZERO;
DOUBLE    vnorm=ZERO;
DOUBLE    pnorm=ZERO;  /* values for norm calculation                  */
DOUBLE    gnorm=ZERO;
NODE    *actnode;      /* actual node                                  */
ARRAY    result_a;
DOUBLE  *result;       /* redundandent result vector                   */

#ifdef DEBUG 
dstrc_enter("fluid_result_incre");
#endif

/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
predof      = fdyn->numdf-1;

/*------------------------- allocate space to allreduce the DIST_VECTOR */
result = amdef("result",&result_a,numeq_total,1,"DV");
         amzero(&result_a);

/*------------------ copy distributed result to redundant result vector */
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );

if (fdyn->itchk==0) /* no iteration convergence check */
{
   /*-----------  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[0].numnp; i++)
   {
      actnode = &(actfield->dis[0].node[i]);
      /*------------------------------ enlarge sol_increment, if necessary */
      if (place >= actnode->sol_increment.fdim)
      {
         diff = place - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
	           actnode->sol_increment.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<actnode->numdf; j++) /* loop dofs */
      {
         dof = actnode->dof[j];
         if (dof>=numeq_total) continue;
         actnode->sol_increment.a.da[place][j] = result[dof];
      } /* end of loop over dofs */
     
   } /* end of loop over nodes */
} /* end if no iteration check */
/*----------------------------------------------------------------------------*/   
else if (fdyn->itchk==1) /* convergence check  */
{
   switch (fdyn->itnorm) /* switch to norm */
   {
   case 0: /* L_infinity norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /*------------------------------ enlarge sol_increment, if necessary */
         if (place >= actnode->sol_increment.fdim)
         {
            diff = place - actnode->sol_increment.fdim;
            max  = IMAX(diff,5);
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
	              actnode->sol_increment.sdim,"DA");
         } /* endif enlargement */
         for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
            if (j==predof) /* pressure dof */
	    {
               dpnorm = DMAX(dpnorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	       pnorm  = DMAX(pnorm, FABS(result[dof]));
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    } /* endif pressure dof */
	    else if (j>predof) /* grid velocity dof */
	    {
	       dgnorm = DMAX(dgnorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	       gnorm  = DMAX(gnorm,FABS(result[dof]));
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    }
	    else /* vel - dof */
	    {	       
	       dvnorm = DMAX(dvnorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	       vnorm  = DMAX(vnorm, FABS(result[dof]));
               actnode->sol_increment.a.da[place][j] = result[dof];
	    } /* endif vel dof */
         } /* end of loop over dofs */        
      } /* end of loop over nodes */      
   break; 
   /*-------------------------------------------------------------------------*/   
   case 1: /* L_1 norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /*------------------------------ enlarge sol_increment, if necessary */
         if (place >= actnode->sol_increment.fdim)
         {
            diff = place - actnode->sol_increment.fdim;
            max  = IMAX(diff,5);
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
	              actnode->sol_increment.sdim,"DA");
         } /* endif enlargement */
         for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
            if (j==predof) /* pressure dof */
	    {
               dpnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	       pnorm  += FABS(result[dof]);
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    } /* endif pressure dof */
            else if (j>predof) /* grid velocity dofs */
	    {
               dgnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	       gnorm  += FABS(result[dof]);
	       actnode->sol_increment.a.da[place][j] = result[dof];	       
	    } /* endif grid velocity dofs */
	    else /* vel - dof */
	    {	       
	       dvnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	       vnorm  += FABS(result[dof]);
               actnode->sol_increment.a.da[place][j] = result[dof];
	    } /* endif vel dof */
         } /* end of loop over dofs */        
      } /* end of loop over nodes */ 
   break; 
   /*-------------------------------------------------------------------------*/   
   case 2: /* L_2 norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         /*------------------------------ enlarge sol_increment, if necessary */
         if (place >= actnode->sol_increment.fdim)
         {
            diff = place - actnode->sol_increment.fdim;
            max  = IMAX(diff,5);
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,
	              actnode->sol_increment.sdim,"DA");
         }
         for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
            if (j==predof) /* pressure dof */
	    {
               dpnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	       pnorm  += pow(result[dof],2);
	       actnode->sol_increment.a.da[place][j] = result[dof];
	    } /* endif pressure dof */
	    else if (j>predof) /* grid velocity dofs */
	    {
               dgnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	       gnorm  += pow(result[dof],2);
	       actnode->sol_increment.a.da[place][j] = result[dof];	    
	    } /* endif grid velocity dofs */
	    else /* vel - dof */
	    {	       
	       dvnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	       vnorm  += pow(result[dof],2);
               actnode->sol_increment.a.da[place][j] = result[dof];
	    } /* endif vel dof */
         } /* end of loop over dofs */        
      } /* end of loop over nodes */
      dvnorm = sqrt(dvnorm);
       vnorm = sqrt( vnorm);
      dpnorm = sqrt(dpnorm);
       pnorm = sqrt( pnorm);  
      dgnorm = sqrt(dgnorm);
       gnorm = sqrt( gnorm);
   break; 
   /*-------------------------------------------------------------------*/   
   default:
      dserror("unknown norm for convergence check!\n");
   } /* end of switch(fdyn->itnorm) */  
   /*------------------------------------------- check for "ZERO-field" */
   if (vnorm<EPS5)
   {
      vnorm = ONE;
#ifdef DEBUT
      printf("ATTENTION: zero vel field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
   }
   if (pnorm<EPS5)
   {
      pnorm = ONE;
#ifdef DEBUG
      printf("ATTENTION: zero pre field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
   }
   if (gnorm<EPS5)
   {
      gnorm = ONE;
   }
   /*------------------------------------- set final convergence ratios */
   *vrat = dvnorm/vnorm;
   *prat = dpnorm/pnorm;
   if (fdyn->freesurf==2) *grat = dgnorm/gnorm;   
} /* endif convergence check */   
/*----------------------------------------------------------------------*/
else
{
   dserror("parameter itchk out of range!\n");
}

/*----------------------------------------------------------------------*/
amdel(&result_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_result_incre */

/*!--------------------------------------------------------------------- 
\brief calculating norms for steady state check

<pre>                                                         genk 05/02

in this routine the velocity and pressure norms for the steady state 
check are calculated:
   norm = ||U(n+1) - U(n)|| / ||U(n)||  		
      solution at (n+1): node->sol_increment[3][j]			
      solution at (n)  : node->sol_increment[1][j]			
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)   actual field
\param  numeq_total   INT	     (i)   total number of equations
\param *vrat	      DOUBLE	     (o)   vel.  conv. ratio
\param *prat	      DOUBLE	     (o)   pres. conv. ratio
\return void 

------------------------------------------------------------------------*/
void fluid_norm(          
                          FLUID_DYNAMIC     *fdyn, 	     
                          FIELD             *actfield,    
		          INT                numeq_total, 
                          DOUBLE            *vrat,        
		          DOUBLE            *prat          
	       )
{
INT         i,j;           /* simply some counters                      */
INT         numdf;         /* number of fluid dofs                      */
INT         numvel;        /* total number of vel-dofs                  */
INT         predof;        /* actual number of pres dof                 */
INT         numnp_total;   /* total number of fluid nodes               */
INT         actdof;        /* actual dof number                         */
DOUBLE      dvnorm=ZERO;   /* norms					*/
DOUBLE       vnorm=ZERO;   /* norms 					*/
DOUBLE      dpnorm=ZERO;   /* norms					*/
DOUBLE       pnorm=ZERO;   /* norms                                     */
NODE       *actnode;       /* actual node                               */

#ifdef DEBUG 
dstrc_enter("fluid_norm");
#endif

/*---------------------------------------------------- set some values */
numdf        = fdyn->numdf;
numnp_total  = actfield->dis[0].numnp;
predof       = numdf-1;
numvel       = numdf-1;

switch (fdyn->stnorm)
{
case 0: /* L_infinity norm */
   /*-------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numvel;j++) /* loop vel-dofs */
      {
	 actdof = actnode->dof[j];
         if (actdof>=numeq_total) continue;
         dvnorm = DMAX(dvnorm,FABS(actnode->sol_increment.a.da[3][j]  \
	                          -actnode->sol_increment.a.da[1][j]));
          vnorm = DMAX( vnorm,FABS(actnode->sol_increment.a.da[3][j]));
      } /* end of loop over vel-dofs */
      actdof = actnode->dof[predof];
      if (actdof>=numeq_total) continue;
      dpnorm = DMAX(dpnorm,FABS(actnode->sol_increment.a.da[3][predof]  \
	                       -actnode->sol_increment.a.da[1][predof]));
       pnorm = DMAX( pnorm,FABS(actnode->sol_increment.a.da[3][predof]));                
   } /* end of loop over nodes */
break; 
/*----------------------------------------------------------------------*/
case 1: /* L_1 norm */   
   /*--------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numvel;j++) /* loop vel-dofs */
      {
	 actdof = actnode->dof[j];
         if (actdof>=numeq_total) continue;
         dvnorm += FABS(actnode->sol_increment.a.da[3][j]  \
	               -actnode->sol_increment.a.da[1][j]);
          vnorm += FABS(actnode->sol_increment.a.da[3][j]);
      } /* end of loop over vel-dofs */
      actdof = actnode->dof[predof];
      if (actdof>=numeq_total) continue;
      dpnorm += FABS(actnode->sol_increment.a.da[3][predof]  \
	            -actnode->sol_increment.a.da[1][predof]);
       pnorm += FABS(actnode->sol_increment.a.da[3][predof]);                
   } /* end of loop over nodes */
break; 
/*----------------------------------------------------------------------*/
case 2: /* L_2 norm */   
   /*--------------------------------------------------- loop all nodes */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numvel;j++) /* loop vel-dofs */
      { 
	 actdof = actnode->dof[j];
         if (actdof>=numeq_total) continue;
         dvnorm += pow(actnode->sol_increment.a.da[3][j]  \
	               -actnode->sol_increment.a.da[1][j],2);
          vnorm += pow(actnode->sol_increment.a.da[3][j],2);
      } /* end of loop over vel-dofs */
      actdof = actnode->dof[predof];
      if (actdof>=numeq_total) continue;
      dpnorm += pow(actnode->sol_increment.a.da[3][predof]  \
	            -actnode->sol_increment.a.da[1][predof],2);
       pnorm += pow(actnode->sol_increment.a.da[3][predof],2);                
   } /* end of loop over nodes */
   dvnorm = sqrt(dvnorm);
    vnorm = sqrt( vnorm);
   dpnorm = sqrt(dpnorm);
    pnorm = sqrt( pnorm);   
break;
/*----------------------------------------------------------------------*/
default:
   dserror("unknown norm for steady state check!\n");
} /* end of switch(fdyn->stnorm) */  

/*---------------------------------------------- check for "ZERO-field" */
if (vnorm<EPS5)
{
   vnorm = ONE;
   printf("ATTENTION: zero vel field - norm <= 1.0e-5 set to 1.0!! \n");
}
if (pnorm<EPS5)
{
   pnorm = ONE;
   printf("ATTENTION: zero pre field - norm <= 1.0e-5 set to 1.0!! \n");
}

/*---------------------------------------- set final convergence ratios */
*vrat = dvnorm/vnorm;
*prat = dpnorm/pnorm;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_norm*/

/*!---------------------------------------------------------------------
\brief copy solution history

<pre>                                                         genk 05/02

in this routine the solution at postion 'from' in the nodal solution 
history is copied to the positon 'to'.	
The pressure is transformed from kinematic to real pressure
			     
</pre>     
\param *actfield      FIELD	     (i)  actual field
\param  disnum        INT	     (i)  number of the discr.
\param  arrayfrom     INT	     (i)  index of the "from-array"
\param  arrayto       INT	     (i)  index of the "to-array"
\param  from          INT            (i)  place in "from-array"
\param  to            INT            (i)  place in "to-array"
\param  numdf         INT            (i)  number of dofs -> 2D or 3D probl.
\return void 
\sa fluid_init()

------------------------------------------------------------------------*/
void fluid_sol_copy(       
                          FIELD             *actfield,
			  INT                disnum,
			  INT                arrayfrom,
			  INT                arrayto,  
                          INT                from,     
		          INT                to,       
		          INT                numdf      
		  )
{
INT         i,j;	/* simply some counters				*/
INT         diff,max;	/* integers for amredef				*/
INT         predof;	/* pressure dof					*/
DOUBLE      dens;	/* density					*/
NODE       *actnode;	/* actual node					*/
ELEMENT    *actele;	/* actual element				*/
DISCRET    *actdis;
ARRAY      *arrayf,*arrayt;

#ifdef DEBUG 
dstrc_enter("fluid_sol_copy");
#endif

predof=numdf-1;  

/* since different materials are not allowed  one can work with the 
   material parameters of any element ---------------------------------*/
actele = &(actfield->dis[disnum].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;

actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayfrom)
   {
   case 0:
      arrayf = &(actnode->sol);
   break;
   case 1:
      arrayf = &(actnode->sol_increment);
   break;
   case 2:
      arrayf = &(actnode->sol_residual);
   break;
   case 3:
      arrayf = &(actnode->sol_mf);
   break;   
   default:
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /*----------------------------------------- select correct arrayfrom */
   switch(arrayto)
   {
   case 0:
      arrayt = &(actnode->sol);
   break;
   case 1:
      arrayt = &(actnode->sol_increment);
   break;
   case 2:
      arrayt = &(actnode->sol_residual);
   break;
   case 3:
      arrayt = &(actnode->sol_mf);
   break;      
   default:
      dserror("Only 0,1,2,3 allowed for arrayfrom to select sol, sol_increment, sol_residual, sol_mf");
   }
   /* check the size of arrayf */
   if (from >= arrayf->fdim)
      dserror("Cannot copy from array, because place doesn't exist");
   /* check the size of arrayt */
   if (to >= arrayt->fdim)
   {
      diff = to - arrayt->fdim;
      /*max  = IMAX(diff,5);*/
      max = diff;
      amredef(arrayt,arrayt->fdim+max+1,arrayt->sdim,"DA");
   }
   for (j=0;j<actnode->numdf;j++)
      arrayt->a.da[to][j] = arrayf->a.da[from][j];  
   /*----------------------------------------------- transform pressure */
   arrayt->a.da[to][predof]*=dens;    
} /* end of loop over nodes */



/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_sol_copy*/

/*!---------------------------------------------------------------------
\brief transform pressure

<pre>                                                         genk 05/02

in this routine the pressure is transformed in the solution history.
  
  index = 0: actnode->sol
  index = 1: actnode->sol_increment
  index = 2: actnode->sol_residual
  index = 3: actnode->sol_mf
  
  option = 0: kinematic pressure -> real pressure
  option = 1: real pressure -> kinematic pressure
  
  (kintematic pressure) = (real pressure)/density
			     
</pre>     
\param *actfield      FIELD	     (i)  actual field
\param  disnum        INT	     (i)  number of the discr.
\param  index         INT	     (i)  index of the array
\param  actpos        INT            (i)  position where to transform
\param  predof        INT            (i)  postion of pressure dof
\return void 
\sa

------------------------------------------------------------------------*/
void fluid_transpres(       
                          FIELD             *actfield,
			  INT                disnum,
			  INT                index,
			  INT                actpos,
                          INT                predof,
			  INT                option     
		    )
{
INT         i;		/* simply some counters				*/
DOUBLE      dens;	/* density					*/
NODE       *actnode;	/* actual node					*/
ELEMENT    *actele;	/* actual element				*/
DISCRET    *actdis;	/* actual discretisation			*/
ARRAY      *array;	/* pointer to solution array			*/

#ifdef DEBUG 
dstrc_enter("fluid_transpres");
#endif


/* since different materials are not allowed  one can work with the 
   material parameters of any element ---------------------------------*/
actele = &(actfield->dis[disnum].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;

actdis = &(actfield->dis[disnum]);
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   /*----------------------------------------- select correct arrayfrom */
   switch(index)
   {
   case 0:
      array = &(actnode->sol);
   break;
   case 1:
      array = &(actnode->sol_increment);
   break;
   case 2:
      array = &(actnode->sol_residual);
   break;
   case 3:
      array = &(actnode->sol_mf);
   break;   
   default:
      dserror("Only 0,1,2,3 allowed for index to select sol, sol_increment, sol_residual, sol_mf");
   }
   dsassert(actpos<array->fdim,"cannot transform pressure\n");
   dsassert(predof<array->sdim,"cannot transform pressure\n");
   /*------------------------------------------ pressure transformation */
   switch (option)
   {
   case 0: /*kinematic pressure -> real pressure  */
      array->a.da[actpos][predof]*=dens;
   break;
   case 1: /*real pressure -> kinematic pressure  */
      array->a.da[actpos][predof]/=dens;
   break;
   default:
      dserror("option out of range: don't know what to do!\n");
   }
} /* end of loop over nodes */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_transpres*/


/*!---------------------------------------------------------------------
\brief steady state check

<pre>                                                         genk 05/02

in this routine the convergence ratios for the steady state check are 
calculated and the result is printed to the screen.
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)  actual field
\param  numeq_total   INT	     (i)  total number of equations
\return INT steady  

------------------------------------------------------------------------*/
INT fluid_steadycheck(    
                          FLUID_DYNAMIC     *fdyn, 	  
                          FIELD             *actfield,   
		          INT                numeq_total 
		     )
{
INT         steady=0;   /* flag for steady state                        */
DOUBLE      vrat,prat;  /* vel. & pres. ratios                          */

#ifdef DEBUG 
dstrc_enter("fluid_steadycheck");
#endif

/*------------------------------------------ determine the conv. ratios */
fluid_norm(fdyn,actfield,numeq_total,&vrat,&prat);

/*------------------------------------------------ output to the screen */
if (par.myrank==0)
{
   switch (fdyn->stnorm) 
   {  
   case 0:
      printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_in] \n",  
	        fdyn->sttol);
   break;
   case 1:
      printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_1 ] \n",
	        fdyn->sttol);
   break;
   case 2:
      printf("   --> steady state check   (tolerance[norm]):  %10.3E [L_2 ] \n",
	        fdyn->sttol);
   break;
   default:
      dserror("Norm for steady state check unknwon!\n");
   } /* end switch (fdyn->stnorm)  */
   printf("         velocities: %10.3E	   pressures:   %10.3E  \n", 
          vrat,prat);
} /* endif (par.myrank) */
/* check if the ratios are smaller than the given tolerance and set flag */	  
if (vrat<fdyn->sttol && prat<fdyn->sttol)    
{
   steady=1;
   if (par.myrank==0)
   {
      printf("\n");
      printf("    >>>>>> STEADY STATE REACHED <<<<<< \n");
      printf("\n");      
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return (steady);
} /* end of fluid_steadycheck*/

/*!---------------------------------------------------------------------
\brief iteration convergence check

<pre>                                                         genk 05/02

in this routine the iteration convergence ratios are compared with
the given tolerance. The result is printed out to the screen.
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param  vrat          DOUBLE  	     (i)  vel. conv. ratio
\param  prat          DOUBLE         (i)  pres. conv. ratio
\param  grat          DOUBLE         (i)  grid. conv. ratio
\param	itnum 	      INT	     (i)  actual numb. of iter steps
\param	te 	      DOUBLE	     (i)  time for element calcul.
\param	ts	      DOUBLE	     (i)  solver time
\return INT converged  

------------------------------------------------------------------------*/
INT fluid_convcheck(      
                          FLUID_DYNAMIC     *fdyn,   
                          DOUBLE             vrat,  
		          DOUBLE             prat,
			  DOUBLE             grat,  
                          INT                itnum, 
		          DOUBLE             te,    
		          DOUBLE             ts     
		   )
{
INT         converged=0;  /* flag for convergence check                  */

#ifdef DEBUG 
dstrc_enter("fluid_convcheck");
#endif

if (fdyn->itchk!=0)
{
   if (par.myrank==0) /* output to the screen */
   { 
      switch(fdyn->itnorm)
      {
      case 0: /* infinity norm */
         printf("|  %3d/%3d   | %10.3E[L_in]  | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n", 
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
      break;
      case 1: /* L_1 norm */
         printf("|  %3d/%3d   | %10.3E[L_1 ]  | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n", 
              itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
      break;
      case 2: /* L_2 norm */
         printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   | {te: %10.3E} {ts:%10.3E} \n", 
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
      break;
      default:
         dserror("Norm for nonlin. convergence check unknown!!\n");
      } /*end of switch(fdyn->itnorm) */
   } /* endif (par.myrank) */
   /*------------------------------------------------ convergence check */
   if (vrat<fdyn->ittol && prat<fdyn->ittol && grat<fdyn->ittol)
      converged=2;
   if (itnum==fdyn->itemax)
      converged++;
   if (converged==1 && par.myrank==0)
   {
      printf("---------------------------------------------------------------- \n");
      printf("|          >>>>>> not converged in itemax steps!               | \n");         
   }
   if (converged>0 && par.myrank==0)
   {
      printf("---------------------------------------------------------------- \n"); 
      printf("\n");
   }
} /* endif (fdyn->itchk) */
else 
{
   if (par.myrank==0) 
   {      
      printf("      iteration step: %3d / 3d \n",  
                              itnum, fdyn->itemax) ;          
      if(itnum==fdyn->itemax)
      {
	 printf("----------------------------------------------------------------- \n"); 
         printf("\n");   
      }
   }
} /* endif (par.myrank) */   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return (converged);
} /* end of fluid_convcheck*/

/*!---------------------------------------------------------------------
\brief print out to the screen

<pre>                                                         genk 05/02

time-integration parameters are printed out to the screen
			     
</pre>   
\param *fdyn 	        FLUID_DYNAMIC   (i)   
\param *dynvar	        FLUID_DYN_CALC  (i/o) 
\return void 

------------------------------------------------------------------------*/
void fluid_algoout(       
                          FLUID_DYNAMIC     *fdyn, 
                          FLUID_DYN_CALC    *dynvar
		  )
{

#ifdef DEBUG 
dstrc_enter("fluid_algoout");
#endif

printf("\n");

switch(fdyn->iop)
{
case 1:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Generalised-Alpha1  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;                  
case 4:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 7:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E     BDF2         STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;         
default:
   dserror("parameter out of range: IOP\n");
} /* end of swtich(fdyn->iop) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_algoout*/
/*!---------------------------------------------------------------------
\brief reduce shearstresses in parallel case

<pre>                                                       he 03/03

shearstress is calcuted within a parellel loop over the elements.
So the result is not available on all procs. Here they are reduced by
MPI to all procs!
			     
</pre>   
\param *actintra         INTRA           (i)   
\param *actfield	       FIELD           (i)    the actual field
\param  *dynvar          FLUID_DYN_CALC  (i)
\return void 

------------------------------------------------------------------------*/
void fluid_reduceshstr(INTRA             *actintra,
                         FIELD             *actfield,
                         FLUID_DYN_CALC    *dynvar)
{
#ifdef PARALLEL
INT      numnp_total;
INT      i;
NODE    *actnode;


#ifdef DEBUG 
dstrc_enter("fluid_reduceshstr");
#endif

/*------------------------------------------- get total number of nodes */
numnp_total = actfield->dis[0].numnp;  
for (i=0; i<numnp_total; i++) /* loop nodes */
{
   actnode = &(actfield->dis[0].node[i]);
   MPI_Bcast(&(actnode->fluid_varia->c_f_shear),1,MPI_DOUBLE,actnode->proc,
             actintra->MPI_INTRA_COMM);
/*------------------- compute shearvelocity for the scaned coordinates */
   if (FABS(actnode->x[0]-dynvar->coord_scale[0])<EPS7 && FABS(actnode->x[1]-dynvar->coord_scale[1])<EPS15)
   MPI_Bcast(&(dynvar->washvel),1,MPI_DOUBLE,actnode->proc,
             actintra->MPI_INTRA_COMM);
} 

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit(); 
#endif

#endif
return;
} /* end of fluid_reduceshstr*/
/*!---------------------------------------------------------------------
\brief reduce shearstresses in parallel case

<pre>                                                       he 03/03

shearstress is calcuted within a parellel loop over the elements.
So the result is not available on all procs. Here they are reduced by
MPI to all procs!
			     
</pre>   
\param *actintra         INTRA           (i)   
\param *actfield	       FIELD           (i)    the actual field
\return void 

------------------------------------------------------------------------*/
void fluid_nullshstr(INTRA             *actintra,
                        PARTITION         *actpart,
                        FIELD             *actfield)
{
#ifdef PARALLEL
INT      numnp_total;
INT      i;
INT      myrank,imyrank;
NODE    *actnode;


#ifdef DEBUG 
dstrc_enter("fluid_nullshstr");
#endif

/*------------------------------------------- get total number of nodes */
myrank = par.myrank;
imyrank= actintra->intra_rank;
numnp_total = actfield->dis[0].numnp;  

if (imyrank!=0 || myrank!=0)
for (i=0; i<actfield->dis[0].numnp; i++)
{
 actnode = &(actfield->dis[0].node[i]);
 actnode->fluid_varia->c_f_shear = ZERO;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif
return;
} /* end of fluid_nullshstr*/
/*!---------------------------------------------------------------------
\brief reduce fluid stresses

<pre>                                                         genk 05/02

fluid stresses are calcuted within a parellel loop over the elements.
So the results are not available on all procs. Here they are reduced by
MPI to all procs!
			     
</pre>   
\param *actintra         INTRA           (i)   
\param *actfield	 FIELD           (i)    the actual field
\param *actsolv          SOLVAR          (i)    the actual solver variables 
\param  numdf            INT             (i)    number of fluid dofs
\param  flag             INT             (i)    a flag
\return void 

------------------------------------------------------------------------*/
void fluid_reducestress(  
                          INTRA             *actintra,
                          FIELD             *actfield,
			  INT                numdf, 
			  FLUID_STRESS       str
		       )
{
#ifdef PARALLEL
INT      numele_total,numnp;
INT      i,j, coupled;
INT      dim;
DOUBLE  *buffer;
ELEMENT *actele;
GNODE   *actgnode;

#ifdef DEBUG 
dstrc_enter("fluid_reducestress");
#endif

numele_total = actfield->dis[0].numele;

switch(str)
{
#ifdef D_FSI
case str_fsicoupling: /* sent only if element has fsi-coupled nodes */
   for(i=0;i<numele_total;i++) /* loop elements */
   {
      actele=&(actfield->dis[0].element[i]);
      numnp = actele->numnp;
      /*--------------- loop nodes and if element has a fsi-coupled node */
      coupled=0;
      for (j=0;j<numnp;j++)
      {
         actgnode = actele->node[j]->gnode;
         /* check if there is a coupled struct node */
         if (actgnode->mfcpnode[genprob.numsf]==NULL) continue;
         coupled=1;
         break;    
      } /* end of loop over nodes */
      if (coupled==1)
      { 
	 /*------------------------------------------- set buffer and dim */
#ifdef D_FLUID2
         if (numdf==3)    
	 {
            buffer=&(actele->e.f2->stress_ND.a.da[0][0]);
	    dim = 3*numnp;
	 }
#endif
#ifdef D_FLUID3
         if (numdf==4)
         {
            buffer=&(actele->e.f2->stress_ND.a.da[0][0]);
	    dim = 6*numnp;
         }
#endif         
	 MPI_Bcast(buffer,dim,MPI_DOUBLE,actele->proc,actintra->MPI_INTRA_COMM);
      } /* endif (coupled==1) */
   } /* end of loop over all elements */
break;
#endif
default:
   dserror("invalid value for 'flag'\n");
} /* end switch (flag) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of fluid_reducestress*/


/*!---------------------------------------------------------------------
\brief calculate fluid acceleration field

<pre>                                                         chfoe 10/03

the fluid acceleration is calculated depending on the actual time 
stepping scheme. It is written to sol_increment[5][i]

One-step-Theta and Generalised Alpha:

acc(n+1) = (vel(n+1)-vel(n)) / (Theta * dt(n)) - (1/Theta -1) * acc(n)


BDF2:

	      2*dt(n)+dt(n-1)		     dt(n)+dt(n-1)
acc(n+1) = --------------------- vel(n+1) - --------------- vel(n) 
	   dt(n)*[dt(n)+dt(n-1)]	     dt(n)*dt(n-1)
	   
		    dt(n)  
	 + ----------------------- vel(n+1)
 	   dt(n-1)*[dt(n)+dt(n-1)]
 
 NOTE: For the One step Theta scheme the first 10 steps are performed
       by means of the BDF2 acceleration. This is necessary to get 
       zero initial field calculations running properly.
 
 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1) 
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)
			     
</pre>   
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\param *dynvar		FLUID_DYN_CALC	(i)
\param *fdyn		FLUID_DYNAMIC	(i)
\return void 

------------------------------------------------------------------------*/
void fluid_acceleration(	FIELD 		*actfield, 
				INT 	 	iop, 
				FLUID_DYN_CALC 	*dynvar,
				FLUID_DYNAMIC	*fdyn
			)
{
DOUBLE fact1, fact2, fact3, dta, dtp, sum;
#ifdef DEBUG 
dstrc_enter("fluid_acceleration");
#endif

switch (iop)
{
case 1:	/* Generalised Alpha time integration 				*/
case 4:	/* One step Theta time integration 				*/
   fact1 = 1.0 / (fdyn->theta*dynvar->dta);
   fact2 =-1.0 / fdyn->theta + 1.0;	/* = -1/Theta + 1		*/
   solserv_sol_zero(actfield,0,1,5);
   solserv_sol_add(actfield,0,1,1,3,5, fact1);
   solserv_sol_add(actfield,0,1,1,1,5,-fact1);
   solserv_sol_add(actfield,0,1,1,4,5, fact2); 
/*   dta = dynvar->dta;
   dtp = dynvar->dtp;
   if (dta*dtp < EPS15)
      dserror("Zero time step size!!!!!");
   sum = dta + dtp;
   fact1 = (2.0 * dta + dtp) / (dta*sum);
   fact2 =-sum / (dta*dtp);
   fact3 = dta / (dtp*sum);
   solserv_sol_zero(actfield,0,1,5);
   solserv_sol_add(actfield,0,1,1,3,5,fact1);
   solserv_sol_add(actfield,0,1,1,1,5,fact2);
   solserv_sol_add(actfield,0,1,1,0,5,fact3); */
break;
case 7:	/* 2nd order backward differencing BDF2				*/
   dta = dynvar->dta;
   dtp = dynvar->dtp;
   if (dta*dtp < EPS15)
      dserror("Zero time step size!!!!!");
   sum = dta + dtp;
   fact1 = (2.0 * dta + dtp) / (dta*sum);
   fact2 =-sum / (dta*dtp);
   fact3 = dta / (dtp*sum);
   solserv_sol_zero(actfield,0,1,5);
   solserv_sol_add(actfield,0,1,1,3,5,fact1);
   solserv_sol_add(actfield,0,1,1,1,5,fact2);
   solserv_sol_add(actfield,0,1,1,0,5,fact3);
break;
default:
   dserror("Time integration scheme unknown for calculation of acceleration!");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_acceleration */


/*!---------------------------------------------------------------------
\brief calculate fluid vector to be multiplied with the mass for time rhs

<pre>                                                         chfoe 10/03

This routine prepares the time rhs in mass form. Depending on the time 
integration scheme the linear combination of vel(n) and acc(n) is 
evaluated which is later multiplied by the mass matrix in order to give
the time rhs.
The linear combination is written to sol_increment[2]


Generalised Alpha:
		         Theta - alpha_m
vel(n) - dt * alpha_f * ---------------- * acc(n)
		            alpha_m

One-step-Theta:

vel(n) + dt*(1-Theta)*acc(n)


BDF2:

for constant time step:

4/3 vel(n) - 1/3 vel(n-1)

for adaptive time step:

/		   2	       \			  2
|	    (dt(n))	       |    		   (dt(n))	      
|1 + ------------------------- | vel(n) - ------------------------- vel(n-1)
\    dt(n-1)*[2*dt(n)+dt(n-1)] /	  dt(n-1)*[2*dt(n)+dt(n-1)]
 
 
 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1) 
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)
 
 
</pre>   
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\return void 

------------------------------------------------------------------------*/
void fluid_prep_rhs(FIELD 		*actfield, 
		    FLUID_DYN_CALC 	*dynvar,
		    FLUID_DYNAMIC	*fdyn)
{
DOUBLE 	fact;

#ifdef DEBUG 
dstrc_enter("fluid_prep_rhs");
#endif

switch (fdyn->iop)
{
case 1:	/* Generalised Alpha time integration 				*/
   fact = -dynvar->dta * fdyn->alpha_f * 
          (fdyn->theta - fdyn->alpha_m) / fdyn->alpha_m; 
   solserv_sol_copy(actfield,0,1,1,1,2);
   solserv_sol_add(actfield,0,1,1,5,2,fact);
break;
case 4:	/* One step Theta time integration 				*/
   fact = dynvar->dta * dynvar->omt;	/* = dt*(1-Theta)		*/
   solserv_sol_copy(actfield,0,1,1,1,2);
   solserv_sol_add(actfield,0,1,1,5,2,fact);
break;
case 7:	/* 2nd order backward differencing BDF2				*/
/*   solserv_sol_zero(actfield,0,1,2);
   solserv_sol_add(actfield,0,1,1,1,2,1.33333333333333333333);
   solserv_sol_add(actfield,0,1,1,0,2,-0.3333333333333333333); */
   fact = DSQR(dynvar->dta)/(dynvar->dtp*(2.0*dynvar->dta+dynvar->dtp));
   solserv_sol_zero(actfield,0,1,2);
   solserv_sol_add(actfield,0,1,1,1,2,1+fact);
   solserv_sol_add(actfield,0,1,1,0,2,-fact);
break;
default:
   dserror("Time integration scheme unknown for mass rhs!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_prep_rhs */



/*!---------------------------------------------------------------------
\brief calculate predictor for adaptive time stepping

<pre>                                                         chfoe 10/03

This routine evaluates the predicted velocity at the end of the next time
step. The prediction is performed via an explicit partner of the implicit
time stepping scheme used. 

One-step-Theta with Theta = 1/2 (Crack Nicolson)
and Generalised Alpha:
with 2nd order Adams-Bashford (AB2):
			     /					           \
		      dt(n) | /	  dt(n)    \	         dt(n)	           |
vel_P(n+1) = vel(n) + ----- || 2 + ------- | * acc(n) - ------- * acc(n-1) |
	 	        2   |\     dt(n-1) /            dt(n-1)            |
		            \					 	  /
			  
One step Theta with Forward Euler: (Theta != 1/2)

vel_P(n+1) = vel(n) + dt(n) * acc(n)

BDF2 with generalised leapfrog:

			     /	 dt(n) 	  \	    /dt(n) \2  
vel_P(n+1) = vel(n) + dt(n) | 1 + ------- |*acc(n)-|-------| * [u(n)-u(n-1)]
			    \ 	dt(n-1)   /	   \dt(n-1)/
 
 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1) 
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)
 sol_increment[6][i] .. predicted vel(n+1)
 
</pre>   
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\param *dynvar		FLUID_DYN_CALC	(i)	
\return void 

------------------------------------------------------------------------*/
void fluid_predictor(FIELD *actfield, INT iop, FLUID_DYN_CALC *dynvar)
{
DOUBLE 	fact1, fact2;

#ifdef DEBUG 
dstrc_enter("fluid_predictor");
#endif

if(dynvar->dtp < EPS15)
   dserror("'zero' previous time step size!");

switch (iop)
{
case 1:
case 4:	/* One step Theta time integration (including TR)		*/
   if (dynvar->theta == 0.5)	/* TR */
   {
      fact1 = dynvar->dta*0.5 * (2.0 + dynvar->dta/dynvar->dtp);
      fact2 =-dynvar->dta*0.5 * dynvar->dta/dynvar->dtp;
      solserv_sol_copy(actfield,0,1,1,1,6);
      solserv_sol_add(actfield,0,1,1,5,6,fact1);
      solserv_sol_add(actfield,0,1,1,4,6,fact2);
   }
   else
   {
      solserv_sol_copy(actfield,0,1,1,1,6);
      solserv_sol_add(actfield,0,1,1,5,6,dynvar->dta);
   }
break;
case 7:	/* 2nd order backward differencing BDF2				*/
   fact1 = dynvar->dta*(1.0+dynvar->dta/dynvar->dtp);
   fact2 = DSQR(dynvar->dta/dynvar->dtp);
   solserv_sol_copy(actfield,0,1,1,1,6);
   solserv_sol_add(actfield,0,1,1,5,6, fact1);
   solserv_sol_add(actfield,0,1,1,1,6,-fact2);
   solserv_sol_add(actfield,0,1,1,0,6, fact2);
break;
default:
   dserror("Time integration scheme unknown for adaptive time stepping!");
}
/*---- copy predicted velocities (no pressure) at (n+1) to sol_field ---*/
solserv_sol_copy(actfield,0,1,1,6,3);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_predictor */


/*!---------------------------------------------------------------------
\brief calculate local truncation error

<pre>                                                         chfoe 10/03

This routine evaluates the local truncation error of the actual time step
depending on the actual solution (sol_increment[3][i]) and the predicted
solution (sol_increment[6][i]). It acts according to the time integration
method iop


2nd Order Generalised Alpha:
with 2nd order Adams-Bashford:

		6*alpha_f - 2
LTE = ------------------------------- * [vel(n+1) - vel_P(n+1)]
      3*( dt(n-1)/dt(n) + 2*alpha_f )
      
One-step-Theta with Theta = 1/2 (Crack Nicolson):
with 2nd order Adams-Bashford:

       vel(n+1) - vel_P(n+1)
LTE = -----------------------
      3*( 1 + dt(n-1)/dt(n) )

One-step_Theta general (Theta != 1/2):
with 1st order Forward Euler:

       /      1    \
LTE = |1 - ------- | * [ vel(n+1) - vel_P(n+1) ]
      \    2*Theta /

BDF2 with generalised leapfrog:

	       /    dt(n-1) \2
	      | 1 + ------- |
	      \     dt(n)  /
LTE = -------------------------------------  * [ vel(n+1) - vel_P(n+1) ]
        dt(n-1)     /dt(n-1)\2    /dt(n-1)\3
     1+3------- + 4| -------| + 2| -------| 
         dt(n)     \  dt(n) /    \  dt(n) /
  
 the values at the nodes are
 sol_increment[0][i] .. vel(n-1)
 sol_increment[1][i] .. vel(n)
 sol_increment[2][i] .. lin. combin. of the other values
 sol_increment[3][i] .. vel(n+1) 
 sol_increment[4][i] .. acc(n-1)
 sol_increment[5][i] .. acc(n)
 sol_increment[6][i] .. predicted vel(n+1)
 sol_increment[7][i] .. local truncation error (lte)
 
</pre>   
\param *actfield	FIELD		(i)	the actual field
\param  iop		INT		(i)	flag, which scheme
\param *dynvar		FLUID_DYN_CALC	(i)	
\return void 

------------------------------------------------------------------------*/
void fluid_lte(	FIELD	 	*actfield, 
		INT 		 iop, 
		FLUID_DYN_CALC 	*dynvar,
		FLUID_DYNAMIC	*fdyn )
{
DOUBLE 	fact, ratio;

#ifdef DEBUG 
dstrc_enter("fluid_lte");
#endif

switch (iop)
{
case 1:	/* Generalised Alpha time integration 				*/
   fact = (6.0 * fdyn->alpha_f - 2.0) / ( 3.0*( dynvar->dtp/dynvar->dta 
           + 2.0*fdyn->alpha_f ) );
   solserv_sol_zero(actfield,0,1,7);
   solserv_sol_add(actfield,0,1,1,3,7, fact);
   solserv_sol_add(actfield,0,1,1,6,7,-fact); 
break;
case 4:	/* One step Theta time integration 				*/
   if (dynvar->theta == 0.5) 	/* TR! */
   {
      fact = 1.0/3.0*(1.0+dynvar->dtp/dynvar->dta);
      solserv_sol_zero(actfield,0,1,7);
      solserv_sol_add(actfield,0,1,1,3,7, fact);
      solserv_sol_add(actfield,0,1,1,6,7,-fact);     
   }
   else
   {
      fact = 1.0 - 1.0/( 2.0*dynvar->theta );
      solserv_sol_zero(actfield,0,1,7);
      solserv_sol_add(actfield,0,1,1,3,7, fact);
      solserv_sol_add(actfield,0,1,1,6,7,-fact);   
   }
break;
case 7:	/* 2nd order backward differencing BDF2				*/
   ratio = dynvar->dtp/dynvar->dta;
   fact  = DSQR(1.0+ratio)/
           ( 1.0 + 3.0*ratio + 4.0*DSQR(ratio) + 2.0*DSQR(ratio)*ratio );
   solserv_sol_zero(actfield,0,1,7);
   solserv_sol_add(actfield,0,1,1,3,7, fact);
   solserv_sol_add(actfield,0,1,1,6,7,-fact); 
break;
default:
   dserror("Time integration scheme unknown for adaptive time stepping!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_lte */



/*!---------------------------------------------------------------------
\brief calculate norm of local truncation error

<pre>                                                         chfoe 10/03

This routine evaluates the norm of the local truncation error and sets 
the new time step. It also decides whether or not a time step should be
repeated in order to secure accuracy.

NOTE: The proposed new time step is actually the proper time step size 
      for the recent step. Hence this is repeated in case the new 
      proposal is much smaller than the step size which was used.

NOTE: The repeat-timestep/don't-repeat-timestep-decision follows the 
      'heuristics' in GRESHO/SANI "Incompressible Flow and the Finite
      Element Method" p. 711 or W.A.Wall (dissertation) p. 98.

</pre>   
\param *actpart		PARTITION	(i) 	actual partition
\param *actfield	FIELD		(i)	the actual field
\param *actintra	INTRA		(i)
\param *dynvar		FLUID_DYN_CALC	(i/o)	
\param *fdyn		FLUID_DYNAMIC	(i/o)
\param *iststep		INT		(i/o)
\param *repeat		INT		(o)	flag, if to repeat
\param *repeated	INT		(i/o)	flag, if was repeated
\return void 

------------------------------------------------------------------------*/
void fluid_lte_norm(	
			PARTITION 	*actpart,
			INTRA		*actintra,
			FLUID_DYN_CALC	*dynvar,
                        FLUID_DYNAMIC	*fdyn,
			INT		*iststep,
			INT		*repeat,
			INT		*repeated
			)
{

INT	 i,j;		/* counters					*/
INT	 nvel;		/* number of free velocity dofs in global dir.	*/
INT	 get;		/* receive value (int) for MPI-process		*/
INT	 numveldof;

DOUBLE	 d_norm;	/* norm of LTE					*/
DOUBLE	 sum;
DOUBLE	 vel0[3];	/* 'characteristic' velocity in global dir. 	*/
DOUBLE	 recv;		/* receive value (double) for MPI-process	*/
DOUBLE   getvec[3];	/* vector to receive values 			*/
DOUBLE	 proposed_dt;	/* the adaptively proposed new delta t		*/
DOUBLE	 ratio;		/* ratio betw. recent step size and new proposal*/

NODE 	*actnode;	/* the actual node				*/
GNODE	*actgnode;	/* the corresponding gnode			*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("fluid_lte_norm");
#endif
/*----------------------------------------------------------------------*/
numveldof = fdyn->numdf - 1;
nvel = 0;
sum = 0.0;
vel0[0] = vel0[1] = vel0[2] = 0.0;

/*---------------------- get maximum velocity in x, y, (z)-direction ---*/
/* loop nodal points of this proc. */
for (i=0; i<actpart->pdis[0].numnp; i++)	
{
   actnode  = (actpart->pdis[0].node[i]);
   if(actnode->numdf > 4)
      dserror("Too many degrees of freedom!!");
   for (j=0;j<actnode->numdf-1;j++) /* loop all velocity dofs */    
   {
      vel0[j] = DMAX(actnode->sol_increment.a.da[3][j],vel0[j]);
   }
}

#ifdef PARALLEL
MPI_Reduce(vel0,getvec,3,MPI_DOUBLE,MPI_MAX,0,actintra->MPI_INTRA_COMM);
#endif

/*------------------------ truncate maximum velocity at lower bound! ---*/
for (i=0; i<3; i++) vel0[i]=DMAX(getvec[i],0.1);	

/*------------------------- determine norm of local truncation error ---*/
/* loop nodal points of this proc. */
for (i=0; i<actpart->pdis[0].numnp; i++)	
{
   actnode  = (actpart->pdis[0].node[i]);
   actgnode = actnode->gnode;
   for (j=0;j<numveldof;j++) /* loop all velocity dofs */    
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]!=0)
      continue; /* do nothing for dbc dofs*/
      sum += DSQR( actnode->sol_increment.a.da[7][j]
	   / ( FABS(actnode->sol_increment.a.da[3][j]) + vel0[j] ) );
      nvel++;
   }
}

#ifdef PARALLEL
MPI_Reduce(&sum,&recv,1,MPI_DOUBLE,MPI_SUM,0,actintra->MPI_INTRA_COMM);
sum = recv;
MPI_Reduce(&nvel,&get,1,MPI_INT,MPI_SUM,0,actintra->MPI_INTRA_COMM);
nvel = get;
#endif

/*-------------------- calculate time step size ratio on PROC 0 only ---*/
if (par.myrank == 0)
{
   d_norm = sqrt( sum / nvel );

   /*----------------------------------------- propose new time step ---*/
   if (d_norm < EPS15)	/* there was almost no error (-> no change)	*/
      ratio = 1.0;
   else
   {
      switch (fdyn->iop)
      {
      case 4: 	/* One step Theta					*/
         {
            if (dynvar->theta == 0.5) /* TR */
   	       ratio = pow((fdyn->lte/d_norm),0.33333333333333333333333);
	    else
               ratio = pow((fdyn->lte/d_norm),0.49);
         }
      break;
      case 1:	/* Generalised Alpha					*/
      case 7:	/* BDF2							*/
         ratio = pow((fdyn->lte/d_norm),0.33333333333333333333333);
      break;
      }	/* end switch (fdyn->iop) */
   }
}	/* end of par.myrank */

#ifdef PARALLEL
MPI_Bcast(&ratio,1,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif

if (*repeated)       /* step has already been repeated  	     */
{
   if (ratio < 0.5)   ratio = 0.5;
   else if (ratio > 1.5)   ratio = 1.5;
   *repeat = 0; 		     /* don't repeat time step       */
   proposed_dt = dynvar->dta * ratio;
   proposed_dt = DMIN(proposed_dt,fdyn->max_dt);
   dynvar->dt_prop = proposed_dt;
   *repeated = 0;
   goto end;
}

proposed_dt = dynvar->dta * ratio;

/* upper step size limit reached     */
if (proposed_dt > fdyn->max_dt && ratio < 5.0)
{
   *repeat = 0;
   dynvar->dt_prop = fdyn->max_dt;
   if (par.myrank == 0) printf("Time step size is cut to MAX_DT!\n");
   goto end;
}

/*----------------------------------------- act with the proposal ---*/
if (ratio <= 0.1)    /* massive reduction of time step -> warning    */
{
   *repeat = 1; 		     /* repeat time step again       */
   if (par.myrank == 0)
      printf("Warning: massive reduction in time step size, something's wrong!!\n");
   fdyn->time -= dynvar->dta;	     /* reset old time values	     */
   dynvar->acttime=fdyn->time;
   fdyn->step--;
   (*iststep--);
   dynvar->dta = proposed_dt;	     /* set smaller step size	     */
}
else if (ratio > 0.1 && ratio <= 0.8)/* significant reduction	     */
{
   *repeat = 1; 		     /* repeat time step again       */
   if (par.myrank == 0)
      printf("Significant time step size reduction. Step is repeated\n");
   fdyn->time -= dynvar->dta;	     /* reset old time values	     */
   dynvar->acttime=fdyn->time;
   fdyn->step--;
   (*iststep)--;
   dynvar->dta = proposed_dt;	     /* set smaller step size	     */
}
else if (ratio > 0.8 && ratio <= 1.5)/* about the same  	     */
{
   *repeat = 0; 		     /* don't repeat time step       */
   proposed_dt = DMIN(proposed_dt,fdyn->max_dt);
   dynvar->dt_prop = proposed_dt;
}
else if (ratio > 1.5 && ratio <= 5.0)/* significant enlargement      */
{
   *repeat = 1; 		     /* repeat time step again       */
   if (par.myrank == 0)
      printf("Significant time step size enlargement. Step is repeated\n");
   fdyn->time -= dynvar->dta;	     /* reset old time values	     */
   dynvar->acttime=fdyn->time;
   fdyn->step--;
   (*iststep)--;
   dynvar->dta = proposed_dt;
}
else if (ratio > 5.0)		     /* huge enlargement	     */
{
   if (dynvar->dta == fdyn->max_dt)  /* max timestep already reached */
   {
   *repeat = 0; 		     /* don't repeat time step       */
   proposed_dt = fdyn->max_dt;
   dynvar->dt_prop = proposed_dt;
   }
   else
   {
      *repeat = 1;		     /* repeat time step again       */
   if (par.myrank == 0)
      printf("Huge time step size enlargement. Step is repeated\n");
      fdyn->time -= dynvar->dta;     /* reset old time values	     */
      dynvar->acttime=fdyn->time;
      fdyn->step--;
      (*iststep)--;
      dynvar->dta = DMIN(5.0 * dynvar->dta,fdyn->max_dt);
   }
}
else
   dserror("Something's wrong in adaptive time stepping.");

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_lte_norm */



/*!--------------------------------------------------------------------- 
\brief calculating errors for beltrami and kim-moin

<pre>                                                         mn 02/04

In this routine different velocity and pressure error norms for the 
beltramii and kim-moin flow are calculated.
			     
</pre>   
\param *actfield      FIELD          (i)   actual field
\param  index         INT            (i)   index for type of flow
                                           index = 1: beltrami
                                           index = 2: kim-moin
\return void 

------------------------------------------------------------------------*/
void fluid_cal_error(
    FIELD             *actfield,
    INT                index
    )
{
  INT         i,j;            /* simply some counters */
  INT         numdf;          /* number of fluid dofs */
  INT         numvel;         /* total number of vel-dofs */
  INT         predof;         /* actual number of pres dof */
  INT         numnp_total;    /* total number of fluid nodes */
  INT         actdof;         /* actual dof number */

  DOUBLE      dvinorm=ZERO;   /* infinity norms */
  DOUBLE       vinorm=ZERO;   /* infinity norms */
  DOUBLE      dpinorm=ZERO;   /* infinity norms */
  DOUBLE       pinorm=ZERO;   /* infinity norms */


  DOUBLE      dv1norm=ZERO;   /* L1 norms */
  DOUBLE       v1norm=ZERO;   /* L1 norms */
  DOUBLE      dp1norm=ZERO;   /* L1 norms */
  DOUBLE       p1norm=ZERO;   /* L1 norms */

  DOUBLE      dv2norm=ZERO;   /* L2 norms */
  DOUBLE       v2norm=ZERO;   /* L2 norms */
  DOUBLE      dp2norm=ZERO;   /* L2 norms */
  DOUBLE       p2norm=ZERO;   /* L2 norms */

  NODE        *actnode;       /* actual node */
  ELEMENT     *actele;
  INT          actmat;
  DOUBLE       visc,a,d,t;
  DOUBLE       x1,x2,x3;
  DOUBLE       u[3],p;

  INT            numeq_total;
  FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */

#ifdef DEBUG 
  dstrc_enter("fluid_cal_error");
#endif

  /*---------------------------------------------------- set some values */
  fdyn = alldyn[genprob.numff].fdyn;
  numdf        = fdyn->numdf;
  numnp_total  = actfield->dis[0].numnp;
  numeq_total  = actfield->dis[0].numeq;
  predof       = numdf-1;
  numvel       = numdf-1;

  switch (index)
  {
    case 1: /* BELTRAMI */

      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
        actnode=&(actfield->dis[0].node[i]);

        /* calculate analytical solution */
        actele = actnode->element[0];
        actmat = actele->mat-1;
        /* set some constants */
        visc   = mat[actmat].m.fluid->viscosity;
        a      = PI/4.0;
        d      = PI/2.0;
        t      = fdyn->time;
        x1   = actnode->x[0];
        x2   = actnode->x[1];
        x3   = actnode->x[2];
        /* calculate analytical values */ 	 
        p    = -a*a/2 * ( exp(2*a*x1) + exp(2*a*x2) + exp(2*a*x3) 
            + 2 * sin(a*x1 + d*x2) * cos(a*x3 + d*x1) * exp(a*(x2+x3))
            + 2 * sin(a*x2 + d*x3) * cos(a*x1 + d*x2) * exp(a*(x3+x1))
            + 2 * sin(a*x3 + d*x1) * cos(a*x2 + d*x3) * exp(a*(x1+x2)) )
          * exp(-2*visc*d*d*t);
        u[0]   = -a * ( exp(a*x1) * sin(a*x2 + d*x3) +
            exp(a*x3) * cos(a*x1 + d*x2) ) * exp(-visc*d*d*t);
        u[1]   = -a * ( exp(a*x2) * sin(a*x3 + d*x1) +
            exp(a*x1) * cos(a*x2 + d*x3) ) * exp(-visc*d*d*t);
        u[2]   = -a * ( exp(a*x3) * sin(a*x1 + d*x2) +
            exp(a*x2) * cos(a*x3 + d*x1) ) * exp(-visc*d*d*t);

        /* infinity norm */
        for (j=0;j<numvel;j++) /* loop vel-dofs */
        {
          actdof = actnode->dof[j];
          if (actdof>=numeq_total) continue;
          dvinorm = DMAX(dvinorm,FABS(actnode->sol_increment.a.da[3][j] - u[j]));
          vinorm = DMAX( vinorm,FABS(u[j]));
        } /* end of loop over vel-dofs */

        actdof = actnode->dof[predof];
        if (actdof>=numeq_total) continue;
        dpinorm = DMAX(dpinorm,FABS(actnode->sol_increment.a.da[3][predof] - p));
        pinorm = DMAX( pinorm,FABS(p));                


        /* L1 norm */
        for (j=0;j<numvel;j++) /* loop vel-dofs */
        {
          actdof = actnode->dof[j];
          if (actdof>=numeq_total) continue;
          dv1norm += FABS(actnode->sol_increment.a.da[3][j] - u[j]);
          v1norm += FABS(u[j]);
        } /* end of loop over vel-dofs */

        actdof = actnode->dof[predof];
        if (actdof>=numeq_total) continue;
        dp1norm += FABS(actnode->sol_increment.a.da[3][predof] - p);
        p1norm += FABS(p);                


        /* L_2 norm */   
        actnode=&(actfield->dis[0].node[i]);
        for (j=0;j<numvel;j++) /* loop vel-dofs */
        { 
          actdof = actnode->dof[j];
          if (actdof>=numeq_total) continue;
          dv2norm += pow(actnode->sol_increment.a.da[3][j] - u[j],2);
          v2norm += pow(u[j],2);
        } /* end of loop over vel-dofs */

        actdof = actnode->dof[predof];
        if (actdof>=numeq_total) continue;
        dp2norm += pow(actnode->sol_increment.a.da[3][predof] - p,2);
        p2norm += pow(p,2);                

      } /* end ol LOOP all nodes */


      dv2norm = sqrt(dv2norm);
      v2norm = sqrt( v2norm);
      dp2norm = sqrt(dp2norm);
      p2norm = sqrt( p2norm);   


      /* check for "ZERO-field" */
      if (vinorm<EPS5)
      {
        vinorm = ONE;
        printf("ATTENTION: zero vel field - inf norm <= 1.0e-5 set to 1.0!! \n");
      }
      if (pinorm<EPS5)
      {
        pinorm = ONE;
        printf("ATTENTION: zero pre field - inf norm <= 1.0e-5 set to 1.0!! \n");
      }

      if (v1norm<EPS5)
      {
        v1norm = ONE;
        printf("ATTENTION: zero vel field - L1 norm <= 1.0e-5 set to 1.0!! \n");
      }
      if (p1norm<EPS5)
      {
        p1norm = ONE;
        printf("ATTENTION: zero pre field - L1 norm <= 1.0e-5 set to 1.0!! \n");
      }

      if (v2norm<EPS5)
      {
        v2norm = ONE;
        printf("ATTENTION: zero vel field - L2 norm <= 1.0e-5 set to 1.0!! \n");
      }
      if (p2norm<EPS5)
      {
        p2norm = ONE;
        printf("ATTENTION: zero pre field - L2 norm <= 1.0e-5 set to 1.0!! \n");
      }

      /*
         printf("\nabsolute vel error for Beltrami in inf norm:  %11.4E \n",  dvinorm);
         printf("absolute pre error for Beltrami in inf norm:  %11.4E \n\n",dpinorm);

         printf("absolute vel error for Beltrami in L1  norm:  %11.4E \n",  dv1norm);
         printf("absolute pre error for Beltrami in L1  norm:  %11.4E \n\n",dp1norm);

         printf("absolute vel error for Beltrami in L2  norm:  %11.4E \n",  dv2norm);
         printf("absolute pre error for Beltrami in L2  norm:  %11.4E \n\n",dp2norm);
         */


      printf("\nErrors for Beltrami flow:\n");
      printf("norm | abs. vel.  | rel. vel.  | abs. pre.  | rel. pre.  |\n");
      printf("----------------------------------------------------------\n");
      printf("inf  |%11.4E |%11.4E |%11.4E |%11.4E |\n",dvinorm,dvinorm/vinorm,dpinorm,dpinorm/pinorm);
      printf("L1   |%11.4E |%11.4E |%11.4E |%11.4E |\n",dv1norm,dv1norm/v1norm,dp1norm,dp1norm/p1norm);
      printf("L2   |%11.4E |%11.4E |%11.4E |%11.4E |\n\n",dv2norm,dv2norm/v2norm,dp2norm,dp2norm/p2norm);

      break;
      /* end of BELTRAMI */

    case 2: /* KIM-MOIN */

      for (i=0;i<numnp_total;i++) /* loop nodes */
      {
        actnode=&(actfield->dis[0].node[i]);

        /* calculate analytical solution */
        actele = actnode->element[0];
        actmat = actele->mat-1;
        /* set some constants */
        visc   = mat[actmat].m.fluid->viscosity;
        a      = 2.0;
        t      = fdyn->time;
        x1   = actnode->x[0];
        x2   = actnode->x[1];
        /* calculate analytical values */ 	 
        p      = -1.0/4.0 * ( cos(2.0*a*PI*x1) + cos(2.0*a*PI*x2) ) * exp(-4.0*a*a*PI*PI*t*visc);
        u[0]   = - cos(a*PI*x1) * sin(a*PI*x2) * exp(-2.0*a*a*PI*PI*t*visc);
        u[1]   = + sin(a*PI*x1) * cos(a*PI*x2) * exp(-2.0*a*a*PI*PI*t*visc);

        /* infinity norm */
        for (j=0;j<numvel;j++) /* loop vel-dofs */
        {
          actdof = actnode->dof[j];
          if (actdof>=numeq_total) continue;
          dvinorm = DMAX(dvinorm,FABS(actnode->sol_increment.a.da[3][j] - u[j]));
          vinorm = DMAX( vinorm,FABS(u[j]));
        } /* end of loop over vel-dofs */

        actdof = actnode->dof[predof];
        if (actdof>=numeq_total) continue;
        dpinorm = DMAX(dpinorm,FABS(actnode->sol_increment.a.da[3][predof] - p));
        pinorm = DMAX( pinorm,FABS(p));                


        /* L1 norm */
        for (j=0;j<numvel;j++) /* loop vel-dofs */
        {
          actdof = actnode->dof[j];
          if (actdof>=numeq_total) continue;
          dv1norm += FABS(actnode->sol_increment.a.da[3][j] - u[j]);
          v1norm += FABS(u[j]);
        } /* end of loop over vel-dofs */

        actdof = actnode->dof[predof];
        if (actdof>=numeq_total) continue;
        dp1norm += FABS(actnode->sol_increment.a.da[3][predof] - p);
        p1norm += FABS(p);                


        /* L_2 norm */   
        actnode=&(actfield->dis[0].node[i]);
        for (j=0;j<numvel;j++) /* loop vel-dofs */
        { 
          actdof = actnode->dof[j];
          if (actdof>=numeq_total) continue;
          dv2norm += pow(actnode->sol_increment.a.da[3][j] - u[j],2);
          v2norm += pow(u[j],2);
        } /* end of loop over vel-dofs */

        actdof = actnode->dof[predof];
        if (actdof>=numeq_total) continue;
        dp2norm += pow(actnode->sol_increment.a.da[3][predof] - p,2);
        p2norm += pow(p,2);                

      } /* end ol LOOP all nodes */


      dv2norm = sqrt(dv2norm);
      v2norm = sqrt( v2norm);
      dp2norm = sqrt(dp2norm);
      p2norm = sqrt( p2norm);   


      /* check for "ZERO-field" */
      if (vinorm<EPS5)
      {
        vinorm = ONE;
        printf("ATTENTION: zero vel field - inf norm <= 1.0e-5 set to 1.0!! \n");
      }
      if (pinorm<EPS5)
      {
        pinorm = ONE;
        printf("ATTENTION: zero pre field - inf norm <= 1.0e-5 set to 1.0!! \n");
      }

      if (v1norm<EPS5)
      {
        v1norm = ONE;
        printf("ATTENTION: zero vel field - L1 norm <= 1.0e-5 set to 1.0!! \n");
      }
      if (p1norm<EPS5)
      {
        p1norm = ONE;
        printf("ATTENTION: zero pre field - L1 norm <= 1.0e-5 set to 1.0!! \n");
      }

      if (v2norm<EPS5)
      {
        v2norm = ONE;
        printf("ATTENTION: zero vel field - L2 norm <= 1.0e-5 set to 1.0!! \n");
      }
      if (p2norm<EPS5)
      {
        p2norm = ONE;
        printf("ATTENTION: zero pre field - L2 norm <= 1.0e-5 set to 1.0!! \n");
      }

      /*
         printf("\nabsolute vel error for Beltrami in inf norm:  %11.4E \n",  dvinorm);
         printf("absolute pre error for Beltrami in inf norm:  %11.4E \n\n",dpinorm);

         printf("absolute vel error for Beltrami in L1  norm:  %11.4E \n",  dv1norm);
         printf("absolute pre error for Beltrami in L1  norm:  %11.4E \n\n",dp1norm);

         printf("absolute vel error for Beltrami in L2  norm:  %11.4E \n",  dv2norm);
         printf("absolute pre error for Beltrami in L2  norm:  %11.4E \n\n",dp2norm);
         */


      printf("\nErrors for Kim-Moin flow:\n");
      printf("norm | abs. vel.  | rel. vel.  | abs. pre.  | rel. pre.  |\n");
      printf("----------------------------------------------------------\n");
      printf("inf  |%11.4E |%11.4E |%11.4E |%11.4E |\n",dvinorm,dvinorm/vinorm,dpinorm,dpinorm/pinorm);
      printf("L1   |%11.4E |%11.4E |%11.4E |%11.4E |\n",dv1norm,dv1norm/v1norm,dp1norm,dp1norm/p1norm);
      printf("L2   |%11.4E |%11.4E |%11.4E |%11.4E |\n\n",dv2norm,dv2norm/v2norm,dp2norm,dp2norm/p2norm);

      break;
  /* end of KIM-MOIN */

    default:
      dserror("Unknown type of flow for error calculation");
      break;
  }
  
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of fluid_cal_error */


#endif
/*! @} (documentation module close)*/

