/*!----------------------------------------------------------------------
\file
\brief service routines for fluid time algorithms

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
void fluid_tcons(         
                          FLUID_DYNAMIC     *fdyn,
                          FLUID_DYN_CALC    *dynvar
		)
{

#ifdef DEBUG 
dstrc_enter("fluid_tcons");
#endif

dynvar->omt = ONE-fdyn->theta;

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
dynvar->totarea=ZERO;

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
INT    numveldof;       /* number of veld dofs                          */
INT    found;
INT    num;
DOUBLE dens;            /* density                                      */
NODE  *actnode;         /* the actual node                              */
GNODE *actgnode;        /* the actual gnode                             */
ELEMENT *actele;        /* the actual element                           */

/*----------------------------------------- variables for solitary wave */
DOUBLE u1,u2,eta,p,c,g,H,d,x,y,t,fac,fac1,sech;

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
INT         i,j,k;        /* simply some counters                       */
INT         numnp_total;  /* total number of nodes                      */
INT         diff,max;     /* integers for amredef                       */
INT         predof;       /* pressure dof                               */
DOUBLE      dens;         /* density                                    */
NODE       *actnode;      /* actual node                                */
ELEMENT    *actele;       /* actual element                             */
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
INT         i,j,k;        /* simply some counters                       */
INT         numnp_total;  /* total number of nodes                      */
DOUBLE      dens;         /* density                                    */
NODE       *actnode;      /* actual node                                */
ELEMENT    *actele;       /* actual element                             */
DISCRET    *actdis;       /* actual discretisation                      */
ARRAY      *array;        /* pointer to solution array                  */

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
   } /* endif (par.myrank)
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
case 2:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Semi-Impl-One-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 3:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Semi-Impl-Two-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 4:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->step,fdyn->nstep);
break;
case 5:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Fract-Step-Theta  STEP = %4d/%4d \n",
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

#endif
/*! @} (documentation module close)*/

