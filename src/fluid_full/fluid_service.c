/*!----------------------------------------------------------------------
\file
\brief service routines for fluid time algorithms

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
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
\param *nfrastep 	int           (o)  number of fract. steps

\return void 
\warning this routine is not completely tested yet!

------------------------------------------------------------------------*/
void fluid_startproc(
                     FLUID_DYNAMIC  *fdyn,
		     int            *nfrastep 
		    )
{
static int     iop_s;
static int     iops_s;
static int     itemax_s;
static double  theta_s ;
static double  thetas_s;

#ifdef DEBUG 
dstrc_enter("fluid_startproc");
#endif 

if (fdyn->step==1)  /* save parameter from input */
{
   iop_s = fdyn->iop;
   iops_s = fdyn->iops;
   itemax_s = fdyn->itemax;
   theta_s  = fdyn->theta;
   thetas_s = fdyn->thetas;
}
if (fdyn->step<=fdyn->nums) /* set parameter for starting phase */
{
   fdyn->iop = iops_s;
   fdyn->theta  = theta_s;
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
\brief calculating time integration constants

<pre>                                                         genk 03/02

in this routine the constants for the time integration algorithms are 
calculated
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC   (i)  
\param *dynvar	        FLUID_DYN_CALC  (i/o)  

\return void 
\warning only ONE-STEP-THETA implemented up to now!

------------------------------------------------------------------------*/
void fluid_tcons(FLUID_DYNAMIC *fdyn,
                 FLUID_DYN_CALC *dynvar)
{

#ifdef DEBUG 
dstrc_enter("fluid_tcons");
#endif

/*----------------------------------------------------- check algorithm */
if (fdyn->iop==4)   /* one step theta */
{
    dynvar->dta  = fdyn->dt;
    dynvar->thsl = fdyn->dt*fdyn->theta;
    dynvar->thpl = dynvar->thsl;
    dynvar->thsr = (ONE - fdyn->theta)*fdyn->dt;
    dynvar->thpr = dynvar->thsr;
}
else
   dserror ("constants for time algorithm not implemented yet!\n");

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
\param  itnum      int  	   (i)     actual number of iterations
\return void 
\warning up to now, only fixed-point like iteration checked!!!

------------------------------------------------------------------------*/
void fluid_icons(FLUID_DYNAMIC *fdyn,
                 FLUID_DYN_CALC *dynvar,
		 int itnum           
		)
{
int i;

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

if(fdyn->ite==0)           /* no iteration */
{
   dynvar->sigma=ZERO;
   dynvar->nik=1;       
   dynvar->nic=2; 
   dynvar->nif=3;      /* KCF */
   dserror("results with no nonlin. iteration not checked yet!\n");
}
else if (fdyn->ite==1)    /* fixed point like iteration */
{
   dynvar->sigma=ZERO;
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
else if (fdyn->ite==2)    /* Newton iteration */
{
   dserror("Newton iteration not checked yet!!!\n");
   dynvar->sigma=ONE;
   if (itnum>1)
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
else if (fdyn->ite==3)    /* fixed point iteration */
{
   dserror("fixed point iteration not checked yet!!!\n");
   dynvar->sigma=-ONE;
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
      sol_increment[0][i]: solution at (n-1)			    
      sol_increment[1][i]: solution at (n)			    
      sol_increment[2][i]: solution at (n+g)			    
      sol_increment[3][i]: solution at (n+1)			    
In node->sol one findes the converged solutions of the time-steps, which
are written to the output-file and the pss-file.

If the initial fluid fuild comes from the input-file 'fluid-start-data',
these data are copied to sol_increment.
			     
</pre>   
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param *dynvar	   FLUID_DYN_CALC  (i/o)  
\param  itnum      int  	   (i)     actual number of iterations
\return void 

------------------------------------------------------------------------*/
void fluid_init(
		FIELD  *actfield,  
                FLUID_DYNAMIC *fdyn	
	       )      
{
int    i,j;          /* simply counters                                 */
int    actmat;       /* number of actual material                       */
int    numdf;        /* number of dofs in this discretisation           */
int    numnp_total;  /* total number of nodes in this discretisation    */
int    numele_total; /* total number of elements in this discr.         */  
int    predof;       /* pressure dof number                             */
double dens;         /* density                                         */
NODE  *actnode;      /* the actual node                                 */
ELEMENT * actele;    /* the actual element                              */

#ifdef DEBUG 
dstrc_enter("fluid_init");
#endif

/*----------------------- set control variables for element evaluation */
fdyn->dynvar.itwost = 0;
fdyn->dynvar.isemim = 0;
fdyn->dynvar.ishape = 1;
fdyn->dynvar.iprerhs= fdyn->iprerhs;

if(fdyn->iop==3) 
   fdyn->dynvar.itwost = 1;
if(fdyn->iop==2 || fdyn->iop==3)   
  fdyn->dynvar.isemim = 1;  
  
/*---------------------------------------------------- set some values */
numdf        = fdyn->numdf;
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
predof       = numdf-1;

/*-------------------------------- allocate space for solution history */
for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   amredef(&(actnode->sol_increment),4,numdf,"DA");
   amzero(&(actnode->sol_increment));
}

/*---------------------- check if there are data from fluid_start_data */
if (fdyn->init>=1)
{
/*------------------------------ copy initial data to solution history */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
      actnode=&(actfield->dis[0].node[i]);
      for (j=0;j<numdf;j++) /* loop dofs */
      {      
         actnode->sol.a.da[0][j]          =fdyn->start.a.da[actnode->Id_loc][j];
	 actnode->sol_increment.a.da[1][j]=fdyn->start.a.da[actnode->Id_loc][j];
	 actnode->sol_increment.a.da[3][j]=fdyn->start.a.da[actnode->Id_loc][j];
      } /* end of loop over dofs */
   } /* end of loop over nodes */
/*--------------------------- transform pressure values of initial data 
             from real pressure to kinematic pressure ------------------*/
    for (i=0;i<numele_total;i++) /* loop elements */
    {
       actele = &(actfield->dis[0].element[i]);
       actmat = actele->mat-1;
       dens   = mat[actmat].m.fluid->density;
       for(j=0;j<actele->numnp;j++) /* loop nodes */
       {
          actnode=actele->node[j];
	  actnode->sol_increment.a.da[1][predof] /= dens;
	  actnode->sol_increment.a.da[3][predof] /= dens;
       } /* end of loop over nodes */
    } /* end of loop over elements */                                   
/*--------------------- free starting field and initial solution vector */
   amdel(&(fdyn->start));
}

/*--------- inherit the neuman conditions from design to discretization */
for (i=0; i<actfield->ndis; i++) inherit_design_dis_neum(&(actfield->dis[i]));

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
\param	 place         int	      (i)    place in sol_incr.
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*vrat          double	      (o)    vel.  conv. ratio
\param	*prat          double	      (o)    pre.  conv. ratio
\param	*fdyn	       FLUID_DYNAMIC	     	
\return void 

------------------------------------------------------------------------*/
void fluid_result_incre(FIELD         *actfield,    
                        INTRA         *actintra,   
			DIST_VECTOR   *sol,        
                        int            place,      
			SPARSE_ARRAY  *sysarray,      
			SPARSE_TYP    *sysarray_typ,
			double        *vrat,        
			double        *prat,       
			FLUID_DYNAMIC *fdyn           
		       )
{
int      i,j;
int      max;
int      diff;
int      dof;
int      predof;
int      numeq_total;
double  *result;
double   dvnorm=ZERO;
double   dpnorm=ZERO;
double    vnorm=ZERO;
double    pnorm=ZERO;
NODE    *actnode;
ARRAY    result_a;

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
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
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
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
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
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
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
            amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
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
   break; 
   /*-------------------------------------------------------------------*/   
   default:
      dserror("unknown norm for convergence check!\n");
   } /* end of switch(fdyn->itnorm) */  
   /*------------------------------------------- check for "ZERO-field" */
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
   /*------------------------------------- set final convergence ratios */
   *vrat = dvnorm/vnorm;
   *prat = dpnorm/pnorm;
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
\param  numeq_total   int	     (i)   total number of equations
\param *vrat	      double	     (o)   vel.  conv. ratio
\param *prat	      double	     (o)   pres. conv. ratio
\return void 

------------------------------------------------------------------------*/
void fluid_norm(FLUID_DYNAMIC *fdyn, 	     
                FIELD         *actfield,    
		int            numeq_total, 
                double        *vrat,        
		double        *prat         
	       )
{
int         i,j;           /* simply some counters                      */
int         numdf;         /* number of fluid dofs                      */
int         numvel;        /* total number of vel-dofs                  */
int         predof;        /* actual number of pres dof                 */
int         numnp_total;   /* total number of fluid nodes               */
int         actdof;        /* actual dof number                         */
double      dvnorm=ZERO;   /* norms					*/
double       vnorm=ZERO;   /* norms 					*/
double      dpnorm=ZERO;   /* norms					*/
double       pnorm=ZERO;   /* norms                                     */
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
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)  actual field
\param  from          int	     (i)  pos. in sol(_increment)
\param  to   	      int	     (i)  pos. in sol(_increment)
\param  flag 	      int	     (i)  simply a flag
\return void 
\sa fluid_init()

------------------------------------------------------------------------*/
void fluid_copysol(FLUID_DYNAMIC *fdyn, 
                   FIELD         *actfield,  
                   int            from,     
		   int            to,       
		   int            flag      
		  )
{
int         i,j;          /* simply some counters                       */
int         numnp_total;  /* total number of nodes                      */
int         diff,max;     /* integers for amredef                       */
int         numdf;        /* number of fluid dofs                       */
NODE       *actnode;      /* actual node                                */

#ifdef DEBUG 
dstrc_enter("fluid_copysol");
#endif

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[0].numnp;
numdf        = fdyn->numdf;

switch(flag)
{
case 0: /* copy from sol_increment[from] to sol_increment[to] */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
    /*------------------------------ enlarge sol_increment, if necessary */
      actnode=&(actfield->dis[0].node[i]);
      if (to >= actnode->sol_increment.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<numdf; j++) /* loop dofs */
      {
        actnode->sol_increment.a.da[to][j]=actnode->sol_increment.a.da[from][j];
      } /* end of loop over dofs */
   } /* end of loop over nodes */
break;
/*----------------------------------------------------------------------*/   
case 1: /* copy from sol_increment[from] to sol[to]*/
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
    /*------------------------------ enlarge sol, if necessary */
      actnode=&(actfield->dis[0].node[i]);
      if (to >= actnode->sol.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol),actnode->sol.fdim+max,actnode->sol.sdim,"DA");
      } /* endif enlargement */
      for (j=0;j<numdf;j++) /* loop dofs */
      {
         actnode->sol.a.da[to][j]=actnode->sol_increment.a.da[from][j];
      } /* end of loop over dofs */
   } /* end of loop over nodes */
break;
/*----------------------------------------------------------------------*/
default:
   dserror("don't know what to do, so I better stop! *gg* \n");
} /* end of switch(flag)

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_copysol*/

/*!---------------------------------------------------------------------                                         
\brief steady state check

<pre>                                                         genk 05/02

in this routine the convergence ratios for the steady state check are 
calculated and the result is printed to the screen.
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param *actfield      FIELD	     (i)  actual field
\param  numeq_total   int	     (i)  total number of equations
\return int steady  

------------------------------------------------------------------------*/
int fluid_steadycheck(FLUID_DYNAMIC *fdyn, 	  
                      FIELD         *actfield,   
		      int            numeq_total 
		     )
{
int         steady=0;   /* flag for steady state                        */
double      vrat,prat;  /* vel. & pres. ratios                          */

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
      printf("   --> steady state check   (tolerance[norm]):  %10.3#E [L_in] \n",  
	        fdyn->sttol);
   break;
   case 1:
      printf("   --> steady state check   (tolerance[norm]):  %10.3#E [L_1 ] \n",
	        fdyn->sttol);
   break;
   case 2:
      printf("   --> steady state check   (tolerance[norm]):  %10.3#E [L_2 ] \n",
	        fdyn->sttol);
   break;
   default:
      dserror("Norm for steady state check unknwon!\n");
   } /* end switch (fdyn->stnorm)  */
   printf("         velocities: %10.3#E	   pressures:   %10.3#E  \n", 
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

return ((int)(steady));
} /* end of fluid_steadycheck*/

/*!---------------------------------------------------------------------                                         
\brief iteration convergence check

<pre>                                                         genk 05/02

in this routine the iteration convergence ratios are compared with
the given tolerance. The result is printed out to the screen.
			     
</pre>   
\param *fdyn 	      FLUID_DYNAMIC  (i)   
\param  vrat          double  	     (i)  vel. conv. ratio
\param  prat          double         (i)  pres. conv. ratio
\param	itnum 	      int	     (i)  actual numb. of iter steps
\param	te 	      double	     (i)  time for element calcul.
\param	ts	      double	     (i)  solver time
\return int converged  

------------------------------------------------------------------------*/
int fluid_convcheck(FLUID_DYNAMIC *fdyn,   
                    double         vrat,  
		    double         prat,  
                    int            itnum, 
		    double         te,    
		    double         ts     
		   )
{
int         converged=0;  /* flag for convergence check                  */

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
         printf("|  %3d/%3d   | %10.3#E[L_in]  | %10.3#E   | %10.3#E  | {te: %10.3#E} {ts:%10.3#E} \n", 
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         printf("|            |                   |              |             | \n");
      break;
      case 1: /* L_1 norm */
         printf("|  %3d/%3d   | %10.3#E[L_1 ]  | %10.3#E   | %10.3#E  | {te: %10.3#E} {ts:%10.3#E} \n", 
              itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         printf("|            |                   |              |             | \n");
      break;
      case 2: /* L_2 norm */
         printf("|  %3d/%3d   | %10.3#E[L_2 ]  | %10.3#E   | %10.3#E  | {te: %10.3#E} {ts:%10.3#E} \n", 
                 itnum,fdyn->itemax,fdyn->ittol,vrat,prat,te,ts);
         printf("|            |                   |              |             | \n");
      break;
      default:
         dserror("Norm for nonlin. convergence check unknown!!\n");
      } /*end of switch(fdyn->itnorm) */
   } /* endif (par.myrank)
   /*------------------------------------------------ convergence check */
   if (vrat<fdyn->ittol && prat<fdyn->ittol)
      converged=2;
   if (itnum==fdyn->itemax)
      converged+=1;
   if (converged==1 && par.myrank==0)
   {
      printf("|            |                   |              |             | \n");
      printf("|          >>>>>> not converged in itemax steps!              | \n");        
      printf("|            |                   |              |             | \n");  
   }
} /* endif (fdyn->itchk) */
else if (par.myrank==0)
{
    printf("      iteration step: %3d / 3d \n",  
            itnum, fdyn->itemax) ;          
} /* endif (par.myrank) */   

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return ((int)(converged));
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
void fluid_algoout(FLUID_DYNAMIC  *fdyn, 
                   FLUID_DYN_CALC *dynvar
		  )
{

#ifdef DEBUG 
dstrc_enter("fluid_algoout");
#endif

printf("\n");

switch(fdyn->iop)
{
case 2:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  Semi-Impl-One-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 3:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  Semi-Impl-Two-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 4:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  One-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 5:
   printf("TIME: %11.4#E/%11.4#E  DT = %11.4#E  Fract-Step-Theta  STEP = %4d/%4d \n",
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
\brief extract digits from integer number

<pre>                                                         genk 04/02		     
</pre>   
\param  num	 int   (i)    integer number
\param *it	 int   (o)    integer on position "thousand"
\param *ih       int   (o)    integer on position "hundred"
\param *id       int   (o)    integer on position "ten"
\param *id       int   (o)    integer on position "one"
\return void 

------------------------------------------------------------------------*/
void intextract(
                int num,    
                int *it,    
		int *ih,    
		int *id,    
		int *io     
	       )
{
int nit, nih, nid, nio;

#ifdef DEBUG 
dstrc_enter("intextract");
#endif

nit = num/1000;
nih = (num-nit*1000)/100;
nid = (num-nit*1000-nih*100)/10;
nio = num -nit*1000-nih*100-nid*10;

*it=nit;
*ih=nih;
*id=nid;
*io=nio;

 /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of intextract*/ 

#endif
