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
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par; 
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!---------------------------------------------------------------------                                         
\brief initialisation of solution history

<pre>                                                        he  12/02

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
\return void 

------------------------------------------------------------------------*/
void fluid_init_tu(
		        FIELD  *actfield,  
                    FLUID_DYNAMIC *fdyn	
	            )      
{
int    i;            /* simply counter                                  */
int    numdf;        /* number of dofs in this discretisation           */
int    numnp_total;  /* total number of nodes in this discretisation    */
int    numele_total; /* total number of elements in this discr.         */  
NODE  *actnode;      /* the actual node                                 */

#ifdef DEBUG  
dstrc_enter("fluid_init_tu");
#endif

/*---------------------------------------------------- set some values */
numdf        = 1;
numnp_total  = actfield->dis[1].numnp;
numele_total = actfield->dis[1].numele;

/*-------------------------------- allocate space for solution history */
for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[1].node[i]);
   amredef(&(actnode->sol_increment),4,4,"DA");
   amzero(&(actnode->sol_increment));
}

/*---------------------------------------------------- set some values */
numdf        = fdyn->numdf;
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;

/*-------------------------------- allocate space for solution history */
for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   amredef(&(actnode->sol_increment),4,numdf*2,"DA");
   amzero(&(actnode->sol_increment));
   actnode->c_f_shear=ZERO;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
  
return;
} /* end of fluid_init_tu */ 
   
/*!---------------------------------------------------------------------                                         
\brief storing results in solution history

<pre>                                                        he  12/02

in this routine the results in the DIST_VECTOR are put to the nodes in
a certain place in ARRAY sol_increment.
Result has to be allreduced and is put to the whole field on each proc.
If necassary the norms for the iteration convergence check of the 
nonlinear iteration scheme are calculated.	   
			     
</pre>   
\param **actfield        FIELD            (i)    actual field       
\param  *actintra        INTRA	      (i)    actual intra comm.
\param	*sol 	       DIST_VECTOR      (i)    solution vector
\param	 place         int	      (i)    place in sol_incr.
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*kapepsrat    double	      (o)    kapeps  conv. ratio
\param	*fdyn	       FLUID_DYNAMIC	     	
\param lower_limit_kappa  double	      (i)    lower limit for kappa
\param lower_limit_eps    double	      (i)    lower limit for epsilon
\return void 

------------------------------------------------------------------------*/
void fluid_result_incre_tu(FIELD         *actfield,    
                           INTRA         *actintra,   
			         DIST_VECTOR   *sol,        
                           int            place,      
			         SPARSE_ARRAY  *sysarray,      
			         SPARSE_TYP    *sysarray_typ,
			         double        *kapepsrat,        
		               FLUID_DYNAMIC *fdyn,
                           double         lower_limit_kappa,
                           double         lower_limit_eps          
		              )
{
int      i,j;
int      max;
int      diff;
int      dof;
int      numeq_total;
double  *result;
double   dkapepsnorm=ZERO;
double    kapepsnorm=ZERO;
NODE    *actnode;
ARRAY    result_a;
FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */

#ifdef DEBUG 
dstrc_enter("fluid_result_incre_tu");
#endif

dynvar      = &(fdyn->dynvar);
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;

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

   switch (fdyn->itnorm) /* switch to norm */
   {
   case 0: /* L_infinity norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
         actnode = &(actfield->dis[1].node[i]);
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
	       if(dynvar->kapeps_flag==0)
             {
              if (result[dof]<lower_limit_kappa) result[dof]=lower_limit_kappa;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j];
              dkapepsnorm = DMAX(dkapepsnorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	        kapepsnorm  = DMAX(kapepsnorm, FABS(result[dof]));
              actnode->sol_increment.a.da[place][j] = result[dof];
             }
             if(dynvar->kapeps_flag==1)
             {
              if (result[dof]<lower_limit_eps) result[dof]=lower_limit_eps;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j+2];
              dkapepsnorm = DMAX(dkapepsnorm,FABS(result[dof]-actnode->sol_increment.a.da[place][j+2]));
	        kapepsnorm  = DMAX(kapepsnorm, FABS(result[dof]));
              actnode->sol_increment.a.da[place][j+2] = result[dof];
             }
          } /* end of loop over dofs */        
      } /* end of loop over nodes */      
   break; 
   /*-------------------------------------------------------------------------*/   
   case 1: /* L_1 norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
         actnode = &(actfield->dis[1].node[i]);
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
 	       if(dynvar->kapeps_flag==0)
             {
              if (result[dof]<lower_limit_kappa) result[dof]=lower_limit_kappa;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j];
              dkapepsnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	        kapepsnorm  += FABS(result[dof]);
              actnode->sol_increment.a.da[place][j] = result[dof];
             }
             if(dynvar->kapeps_flag==1)
             {
              if (result[dof]<lower_limit_eps) result[dof]=lower_limit_eps;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j+2];
              dkapepsnorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j+2]);
	        kapepsnorm  += FABS(result[dof]);
              actnode->sol_increment.a.da[place][j+2] = result[dof];   
             }

         } /* end of loop over dofs */        
      } /* end of loop over nodes */ 
   break; 
   /*-------------------------------------------------------------------------*/   
   case 2: /* L_2 norm */
      /*-----------  loop nodes and put the result back to the node structure */
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
         actnode = &(actfield->dis[1].node[i]);
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
             if(dynvar->kapeps_flag==0)
             { 
              if (result[dof]<lower_limit_kappa) result[dof]=lower_limit_kappa;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j];
              dkapepsnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	        kapepsnorm  += pow(result[dof],2);
              actnode->sol_increment.a.da[place][j] = result[dof];   
             }
             if(dynvar->kapeps_flag==1)
             { 
              if (result[dof]<lower_limit_eps) result[dof]=lower_limit_eps;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j+2];
              dkapepsnorm += pow(result[dof]-actnode->sol_increment.a.da[place][j+2],2);
	        kapepsnorm  += pow(result[dof],2);
              actnode->sol_increment.a.da[place][j+2] = result[dof];   
             }
          } /* end of loop over dofs */        
      } /* end of loop over nodes */
      dkapepsnorm = sqrt(dkapepsnorm);
       kapepsnorm = sqrt( kapepsnorm);
   break; 
   /*-------------------------------------------------------------------*/   
   default:
      dserror("unknown norm for convergence check!\n");
   } /* end of switch(fdyn->itnorm) */  
   /*------------------------------------------- check for "ZERO-field" */
   if (kapepsnorm<EPS10)
   {
      kapepsnorm = ONE;
      printf("ATTENTION: zero kapeps field - norm <= 1.0e-10 set to 1.0!! \n");
   }
/*------------------------------------- set final convergence ratios */
   *kapepsrat = dkapepsnorm/kapepsnorm;
/*----------------------------------------------------------------------*/
amdel(&result_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_result_incre_tu */

/*!---------------------------------------------------------------------                                         
\brief routine to set dirichlet boundary conditions at time <time>

<pre>                                                         he  12/02

in this routine the dirichlet boundary conditions for fluid2 and fluid3
elements are set at time <T=fdyn->time>.
the actual dirichlet values are written to the solution history of the
nodes:
    'actnode->sol_increment.a.da[3][j] 
                                 |        
                            time (n+1)                
                                     
</pre>
\param *actfield    FIELD         (i)  actual field (fluid)   
\param *fdyn	  FLUID_DYNAMIC (i)  
\param *lower_limit_kappa  double (o) lower limit for kappa  

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_set_check_tu(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       double lower_limit_kappa)
{
int        i,j;
int        numnp_total;              /* total number of fluid nodes     */
int        numele_total;             /* total number of fluid elements  */
int        numdf;	                   /* number of fluid dofs       	*/
double     k_2,int_lenght;
GNODE     *actgnode;	             /* actual GNODE		            */
NODE      *actnode;	             /* actual NODE		            */

#ifdef DEBUG 
dstrc_enter("fluid_set_check_tu");
#endif 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[1].numnp;
numele_total = actfield->dis[1].numele;
numdf        = 1;
int_lenght   = fdyn->lenght;

/*-------------------- loop all nodes and set actual dirichlet condition */
for (i=0;i<numnp_total;i++) 
{
   actnode  = &(actfield->dis[1].node[i]); 
   actgnode = actnode->gnode;      
   for (j=0;j<numdf;j++) /* loop dofs */
   {
      if (actgnode->dirich->dirich_onoff.a.iv[j]==0 || actgnode->dirich->dirich_onoff.a.iv[j+1]==0)
      {
      actnode->sol_increment.a.da[3][j]   = lower_limit_kappa*10000;
      actnode->sol_increment.a.da[1][j]   = actnode->sol_increment.a.da[3][j];
      actnode->sol_increment.a.da[3][j+1] = 0.005*int_lenght*sqrt(actnode->sol_increment.a.da[3][j]);
      actnode->sol_increment.a.da[2][j+1] = actnode->sol_increment.a.da[3][j+1];
      actnode->sol_increment.a.da[3][j+2] = 0.09 * pow(actnode->sol_increment.a.da[3][j],1.5) / (0.005*int_lenght);
      actnode->sol_increment.a.da[1][j+2] = actnode->sol_increment.a.da[3][j+2];
     }
   } /*end loop over nodes */
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_set_check_tu*/
   
/*!---------------------------------------------------------------------                                         
\brief storing results for eddy viscosity in solution history

<pre>                                                        he  12/02

</pre>   
\param **actfield      FIELD	            (i)    actual field       
\param	*sol 	     DIST_VECTOR        (i)    solution vector
\return void 

------------------------------------------------------------------------*/
void fluid_eddy_update(FIELD         *actfield, 
                       DIST_VECTOR   *sol   
                       )
{
int  i,j;
int  dof,actmat;
int  numeq_total;
double  visc,R_t,C_u;
ELEMENT *actele;
NODE    *actnode;

#ifdef DEBUG 
dstrc_enter("fluid_eddy_update");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;  
  
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
        actnode = &(actfield->dis[1].node[i]);
        actele  = actnode->element[0];
        actmat  = actele->mat-1;
        visc    = mat[actmat].m.fluid->viscosity;

        for (j=0; j<actnode->numdf; j++) /* loop dofs  */
        {
	   dof = actnode->dof[j];
         if (dof>=numeq_total) continue;
        
         if(actnode->sol_increment.a.da[3][2]!=0)
         {  
          R_t = pow(actnode->sol_increment.a.da[3][0],2)/
                (actnode->sol_increment.a.da[3][2]*visc);
          
          C_u = 0.09 * exp(-3.4/pow(1+R_t/50,2)); 
         
          actnode->sol_increment.a.da[3][1] = C_u*pow(actnode->sol_increment.a.da[3][0],2)/
                                              actnode->sol_increment.a.da[3][2];
        }
       }
      }   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_eddy_update */

/*!---------------------------------------------------------------------                                         
\brief storing results for eddy viscosity in solution history.
Write eddy-viscosity from step n+1 to n for production term. 

<pre>                                                        he  12/02

</pre>   
\param **actfield      FIELD	            (i)    actual field       
\return void 

------------------------------------------------------------------------*/
void fluid_eddy_pro(FIELD         *actfield 
                  )
{
int  i;
NODE    *actnode;

#ifdef DEBUG 
dstrc_enter("fluid_eddy_pro");
#endif
/*----------------------------------------------------------------------*/
  
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
        actnode = &(actfield->dis[1].node[i]);
         
        actnode->sol_increment.a.da[2][1] = actnode->sol_increment.a.da[3][1];
      }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_eddy_pro */


/*!---------------------------------------------------------------------                                         
\brief storing results for characteristic lenght in solution history
and calculate norms for the iteration convergence check

<pre>                                                        he  12/02

</pre>   
\param **actfield      FIELD	            (i)    actual field       
\param	*sol 	     DIST_VECTOR        (i)    solution vector
\param  lenghtrat	     double             (o)    conv. ratio
\param     *fdyn 	     FLUID_DYNAMIC      (i)   
\return void 

------------------------------------------------------------------------*/
void fluid_lenght_update(FIELD         *actfield, 
                         DIST_VECTOR   *sol,   
		             double        *lenghtrat, 
                         FLUID_DYNAMIC *fdyn      
                        )
{
NODE    *actnode;
int  i,j;
int  dof,actmat;
double   dlenghtnorm=ZERO;
double    lenghtnorm=ZERO;
int  numeq_total;
double  visc,R_t,C_u;
ELEMENT *actele;

#ifdef DEBUG 
dstrc_enter("fluid_lenght_update");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;  
  
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
         actnode = &(actfield->dis[1].node[i]);
         actele  = actnode->element[0];
         actmat  = actele->mat-1;
         visc    = mat[actmat].m.fluid->viscosity;
         for (j=0; j<actnode->numdf; j++) /* loop dofs  */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
             if(actnode->sol_increment.a.da[3][2]!=0)
             {
              R_t = pow(actnode->sol_increment.a.da[3][0],2)/
                    (actnode->sol_increment.a.da[3][2]*visc);
              C_u = 0.09 * exp(-3.4/pow(1+R_t/50,2)); 

              actnode->sol_increment.a.da[3][3] = 0.5*C_u*pow(actnode->sol_increment.a.da[3][0],1.5)/ actnode->sol_increment.a.da[3][2]
                                                 +0.5*actnode->sol_increment.a.da[1][3];
              switch (fdyn->itnorm) /* switch to norm */
              {
               case 0: /* L_infinity norm */
               dlenghtnorm  = DMAX(dlenghtnorm,FABS(actnode->sol_increment.a.da[3][3]-actnode->sol_increment.a.da[1][3]));
	         lenghtnorm   = DMAX(lenghtnorm, FABS(actnode->sol_increment.a.da[3][3]));
               actnode->sol_increment.a.da[1][3] = actnode->sol_increment.a.da[3][3];
               break; 
/*-------------------------------------------------------------------------*/   
               case 1: /* L_1 norm */
               dlenghtnorm  += FABS(actnode->sol_increment.a.da[3][3]-actnode->sol_increment.a.da[1][3]); 
	         lenghtnorm   += FABS(actnode->sol_increment.a.da[3][3]);
               actnode->sol_increment.a.da[1][3] = actnode->sol_increment.a.da[3][3];
               break; 
/*-------------------------------------------------------------------------*/   
               case 2: /* L_2 norm */
               dlenghtnorm  += pow(actnode->sol_increment.a.da[3][3]-actnode->sol_increment.a.da[1][3],2); 
	         lenghtnorm   += pow(actnode->sol_increment.a.da[3][3],2);
               actnode->sol_increment.a.da[1][3] = actnode->sol_increment.a.da[3][3];
               break; 
/*-------------------------------------------------------------------*/   
              default:
              dserror("unknown norm for convergence check!\n");
             } /* end of switch(fdyn->itnorm) */  
            } /* end if */ 
           } /* for j */
        }  /* for i */ 
      
       if (fdyn->itnorm==2)
       {
        dlenghtnorm = sqrt(dlenghtnorm);
        lenghtnorm  = sqrt(lenghtnorm);
       }

/*------------------------------------------- check for "ZERO-field" */
   if (lenghtnorm<EPS10)
   {
      lenghtnorm = ONE;
      printf("ATTENTION: zero lenght field - norm <= 1.0e-10 set to 1.0!! \n");
   }

/*------------------------------------- set final convergence ratios */
   *lenghtrat = dlenghtnorm/lenghtnorm;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_lenght_update */

/*!---------------------------------------------------------------------                                         
\brief iteration convergence check

<pre>                                                         he  12/02

in this routine the iteration convergence ratios are compared with
the given tolerance. The result is printed out to the screen.
			     
</pre>   
\param *fdyn 	        FLUID_DYNAMIC  (i)   
\param  kapepsrat       double  	     (i)  kapeps or lenght conv. ratio 
\param	itnum1      int	           (i)  actual numb. of iter steps
\param	te 	      double	     (i)  time for element calcul.
\param	ts	      double	     (i)  solver time
\return int converged  

------------------------------------------------------------------------*/
int fluid_convcheck_tu(FLUID_DYNAMIC *fdyn,   
                       double         kapepsrat,  
                       int            itnum1, 
		           double         te,    
		           double         ts     
		          )
{
int         converged=0;  /* flag for convergence check                  */

#ifdef DEBUG 
dstrc_enter("fluid_convcheck_tu");
#endif

   if (par.myrank==0) /* output to the screen */
   { 
      switch(fdyn->itnorm)
      {
      case 0: /* infinity norm */
         printf("|  %3d/%3d   | %10.3E[L_in]  | %10.3E     | {te: %10.3E} {ts:%10.3E} \n", 
                 itnum1,fdyn->itemax_ke,fdyn->ittol,kapepsrat,te,ts);
         printf("|            |                   |                |  \n");
      break;
      case 1: /* L_1 norm */
         printf("|  %3d/%3d   | %10.3E[L_1 ]  | %10.3E     | {te: %10.3E} {ts:%10.3E} \n", 
              itnum1,fdyn->itemax_ke,fdyn->ittol,kapepsrat,te,ts);
         printf("|            |                   |                |  \n");
      break;
      case 2: /* L_2 norm */
         printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E     | {te: %10.3E} {ts:%10.3E} \n", 
                 itnum1,fdyn->itemax_ke,fdyn->ittol,kapepsrat,te,ts);
         printf("|            |                   |                |   \n");
      break;
      default:
         dserror("Norm for nonlin. convergence check unknown!!\n");
      } /*end of switch(fdyn->itnorm) */
   } /* endif (par.myrank)
   /*------------------------------------------------ convergence check */
   if (kapepsrat<fdyn->ittol)
      converged=2;
   if (itnum1==fdyn->itemax_ke)
      converged+=1;
   if (converged==1 && par.myrank==0)
   {
      printf("|            |                   |                |    \n");
      printf("|      >>>>>> not converged in itemax steps!      |    \n");        
      printf("|            |                   |                |    \n");  
   }

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return ((int)(converged));
} /* end of fluid_convcheck_tu*/

/*!---------------------------------------------------------------------                                         
\brief iteration convergence check

<pre>                                                         he  12/02

in this routine the iteration convergence ratios are compared with
the given tolerance for rans.
			     
</pre>   
\param *fdyn 	        FLUID_DYNAMIC  (i)   
\param **actfield         FIELD	     (i)    actual field       
\para   itnum_check       int	           (i)    actual numb. of iter steps
\return int converged  

------------------------------------------------------------------------*/
int fluid_convcheck_test(
                     FLUID_DYNAMIC *fdyn, 
                     FIELD         *actfield, 
                     int            itnum_check 
		         )
{
NODE      *actnode;
int       i,j;
int      predof;
double   dvnorm=ZERO;
double   dpnorm=ZERO;
double    vnorm=ZERO;
double    pnorm=ZERO;
double    vrat,prat;
int       converged=0;  /* flag for convergence check                  */
int       numdf;

#ifdef DEBUG 
dstrc_enter("fluid_convcheck_test");
#endif
/*----------------------------------------------------------------------*/
predof      = fdyn->numdf-1;  
numdf       = fdyn->numdf;

      for (i=0; i<actfield->dis[0].numnp; i++)
      {
         actnode = &(actfield->dis[0].node[i]);
         for (j=0; j<actnode->numdf; j++) /* loop dofs  */
         {
          switch (fdyn->itnorm) /* switch to norm */
          {
           case 0: /* L_infinity norm */
           if (j==predof) /* pressure dof */
	     {
            dpnorm = DMAX(dpnorm,FABS(actnode->sol_increment.a.da[3][j]-actnode->sol_increment.a.da[3][j+numdf])); 
	      pnorm  = DMAX(pnorm, FABS(actnode->sol_increment.a.da[3][j]));
           } /* endif pressure dof */
	     else /* vel - dof */
	     {	       
	      dvnorm = DMAX(dvnorm,FABS(actnode->sol_increment.a.da[3][j]-actnode->sol_increment.a.da[3][j+numdf]));
	      vnorm  = DMAX(vnorm, FABS(actnode->sol_increment.a.da[3][j]));
	     } /* endif vel dof */
           break; 
            case 1: /* L_1 norm */
            if (j==predof) /* pressure dof */
	      {
             dpnorm += FABS(actnode->sol_increment.a.da[3][j]-actnode->sol_increment.a.da[3][j+numdf]);
	       pnorm  += FABS(actnode->sol_increment.a.da[3][j]);
	      } /* endif pressure dof */
	      else /* vel - dof */
	      {	       
	       dvnorm += FABS(actnode->sol_increment.a.da[3][j]-actnode->sol_increment.a.da[3][j+numdf]);
	       vnorm  += FABS(actnode->sol_increment.a.da[3][j]);
	      } /* endif vel dof */
           break; 
            case 2: /* L_2 norm */
            if (j==predof) /* pressure dof */
	      {
             dpnorm += pow(actnode->sol_increment.a.da[3][j]-actnode->sol_increment.a.da[3][j+numdf],2);
	       pnorm  += pow(actnode->sol_increment.a.da[3][j],2);
	      } /* endif pressure dof */
	      else /* vel - dof */
	      {	       
	       dvnorm += pow(actnode->sol_increment.a.da[3][j]-actnode->sol_increment.a.da[3][j+numdf],2);
	       vnorm  += pow(actnode->sol_increment.a.da[3][j],2);
	      } /* endif vel dof */
           } /* end of switch */  
         } /* end of loop over dofs */  
      } /* end of loop over nodes */      

       if (fdyn->itnorm==2)
       {
        dvnorm = sqrt(dvnorm);
         vnorm = sqrt(vnorm);
        dpnorm = sqrt(dpnorm);
         pnorm = sqrt( pnorm);   
       }

/*------------------------------------- set final convergence ratios */
   vrat = dvnorm/vnorm;
   prat = dpnorm/pnorm;

   if (vrat<fdyn->ittol && prat<fdyn->ittol)
      converged=2;
   if (itnum_check==fdyn->itemax_ke)
      converged+=1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return ((int)(converged));
} /* end of fluid_convcheck_test*/

/*!---------------------------------------------------------------------                                         
\brief setting flags for TIME-RHS (kappa)

<pre>                                                         he 12/02


nikapeps <->  EVALUATION OF "TIME - RHS" (F-hat)		    

</pre>   
\param *fdyn	   FLUID_DYNAMIC   (i)  
\param *dynvar	   FLUID_DYN_CALC  (i/o)  
\param  itnum1       int  	       (i)     actual number of iterations
\param  itnumke      int  	       (i)     actual number of iterations
\param  itnum_n      int  	       (i)     actual number of iterations
\return void 

------------------------------------------------------------------------*/
void fluid_icons_tu(FLUID_DYNAMIC *fdyn,
                    FLUID_DYN_CALC *dynvar,
		        int itnum1,
		        int itnumke,
                    int itnum_n
		        )
{

#ifdef DEBUG 
dstrc_enter("fluid_icons_tu");
#endif

/*----------------------------------------------------- initialisation */
dynvar->nis     = 0;
dynvar->kappan  = itnum1 + itnumke;  

dynvar->niturbu_n = 1;
if ((itnum1+itnum_n)>2)
{
 dynvar->niturbu_n=0;  
}

dynvar->niturbu_pro = 1;
if ((itnum1+itnumke)>2)
{
 dynvar->niturbu_pro=0;  
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_icons_tu*/


/*!---------------------------------------------------------------------                                         
\brief copy solution history

<pre>                                                         he 01/03

in this routine the solution at postion 'from' in the nodal solution 
history is copied to the positon 'to'.	
			     
</pre>   
\param *fdyn 	   FLUID_DYNAMIC  (i)   
\param *actfield     FIELD	      (i)  actual field
\param  from         int	      (i)  pos. in sol(_increment)
\param  to   	   int	      (i)  pos. in sol(_increment)
\param  flag 	   int	      (i)  simply a flag
\return void 

------------------------------------------------------------------------*/
void fluid_copysol_tu(FLUID_DYNAMIC *fdyn, 
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
double     visc;
NODE       *actnode;      /* actual node                                */
FLUID_DYN_CALC *dynvar;   /* pointer to fluid_dyn_calc                   */

#ifdef DEBUG 
dstrc_enter("fluid_copysol_tu");
#endif

dynvar       = &(fdyn->dynvar);

for(i=0; i<genprob.nmat; i++)
{
 if (mat[i].mattyp == m_fluid) visc=mat->m.fluid->viscosity;
} 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[1].numnp;
numdf        = 1;

switch(flag)
{
case 0: /* copy from sol_increment[from] to sol_increment[to] */
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
    /*------------------------------ enlarge sol_increment, if necessary */
      actnode=&(actfield->dis[1].node[i]);
      if (to >= actnode->sol_increment.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol_increment),actnode->sol_increment.fdim+max,actnode->sol_increment.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<numdf; j++) /* loop dofs */
      {
        actnode->sol_increment.a.da[to][j]  =actnode->sol_increment.a.da[from][j];
        actnode->sol_increment.a.da[to][j+2]=actnode->sol_increment.a.da[from][j+2];
      } /* end of loop over dofs */
   } /* end of loop over nodes */
break;
/*----------------------------------------------------------------------*/   
case 1: /* copy from sol_increment[from] to sol[to]*/
  if (dynvar->washvel == 0.0) dynvar->washvel=1.0;

  if(fdyn->turbu==2)
  {
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
    /*------------------------------ enlarge sol, if necessary */
      actnode=&(actfield->dis[1].node[i]);
      if (to >= actnode->sol.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol),actnode->sol.fdim+max,actnode->sol.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<numdf; j++) /* loop dofs */
      {
        actnode->sol.a.da[to][j]  =actnode->sol_increment.a.da[from][j]/pow(dynvar->washvel,2);
        actnode->sol.a.da[to][j+1]=actnode->sol_increment.a.da[from][j+1];
        actnode->sol.a.da[to][j+2]=actnode->sol_increment.a.da[from][j+2]*visc/pow(dynvar->washvel,4);
      } /* end of loop over dofs */
   } /* end of loop over nodes */
  }
  if(fdyn->turbu==3)
  {
   for (i=0;i<numnp_total;i++) /* loop nodes */
   {
    /*------------------------------ enlarge sol, if necessary */
      actnode=&(actfield->dis[1].node[i]);
      if (to >= actnode->sol.fdim)
      {
         diff = to - actnode->sol_increment.fdim;
         max  = IMAX(diff,5);
         amredef(&(actnode->sol),actnode->sol.fdim+max,actnode->sol.sdim,"DA");
      } /* endif enlargement */
      for (j=0; j<numdf; j++) /* loop dofs */
      {
        actnode->sol.a.da[to][j]  =actnode->sol_increment.a.da[from][j]/pow(dynvar->washvel,2);
        actnode->sol.a.da[to][j+1]=actnode->sol_increment.a.da[from][j+1];
        actnode->sol.a.da[to][j+2]=actnode->sol_increment.a.da[from][j+2]*visc/pow(dynvar->washvel,2);
      } /* end of loop over dofs */
   } /* end of loop over nodes */
  }
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
} /* end of fluid_copysol_tu*/

/*!---------------------------------------------------------------------                                         
\brief copy solution history

<pre>                                                        he  01/03

in this routine the solution at postion 'from' in the nodal solution 
history is copied to the positon 'to' and the result from is copied 
due to convergence check for RANS.	
			     
</pre>   
\param *fdyn 	    FLUID_DYNAMIC (i)   
\param *actfield      FIELD	     (i)  actual field
\param  from          int	     (i)  pos. in sol(_increment)
\param  to   	    int	     (i)  pos. in sol(_increment)
\return void 

------------------------------------------------------------------------*/

void fluid_copysol_test(FLUID_DYNAMIC *fdyn, 
                   FIELD         *actfield,  
                   int            from,     
                   int            to     
		       )
{
int         i,j;          /* simply some counters                       */
int         numnp_total;  /* total number of nodes                      */
int         diff,max;     /* integers for amredef                       */
int         numdf;        /* number of fluid dofs                       */
NODE       *actnode;      /* actual node                                */

#ifdef DEBUG 
dstrc_enter("fluid_copysol_test");
#endif

/*----------------------------------------------------- set some values */
numdf        = fdyn->numdf;
numnp_total  = actfield->dis[0].numnp;

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
       actnode->sol_increment.a.da[from][j+numdf]=actnode->sol_increment.a.da[from][j];
      } /* end of loop over dofs */
   } /* end of loop over nodes */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_copysol_test*/


/*!---------------------------------------------------------------------                                         
\brief print out to the screen

<pre>                                                         he 01/03

time-integration parameters are printed out to the screen
			     
</pre>   
\param *fdyn 	        FLUID_DYNAMIC   (i)   
\param *dynvar	        FLUID_DYN_CALC  (i/o) 
\return void 

------------------------------------------------------------------------*/
void fluid_algoout_tu(FLUID_DYNAMIC  *fdyn, 
                      FLUID_DYN_CALC *dynvar
		         )
{

#ifdef DEBUG 
dstrc_enter("fluid_algoout_tu");
#endif

printf("\n");

switch(fdyn->iop)
{
case 2:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Semi-Impl-One-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->stepke,fdyn->nstep);
break;
case 3:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Semi-Impl-Two-Step  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->stepke,fdyn->nstep);
break;
case 4:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->stepke,fdyn->nstep);
break;
case 5:
   printf("TIME: %11.4E/%11.4E  DT = %11.4E  Fract-Step-Theta  STEP = %4d/%4d \n",
          fdyn->time,fdyn->maxtime,dynvar->dta,fdyn->stepke,fdyn->nstep);
break;         
default:
   dserror("parameter out of range: IOP\n");
} /* end of swtich(fdyn->iop) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_algoout_tu*/

/*!---------------------------------------------------------------------                                         
\brief calculating time integration constants

<pre>                                                         he  01/03

in this routine the constants for the time integration algorithms are 
calculated
			     
</pre>   
\param *fdyn		FLUID_DYNAMIC   (i)  
\param *dynvar	        FLUID_DYN_CALC  (i/o)  

\return void 
\warning only ONE-STEP-THETA implemented up to now!

------------------------------------------------------------------------*/
void fluid_tcons_tu(FLUID_DYNAMIC *fdyn,
                    FLUID_DYN_CALC *dynvar)
{

#ifdef DEBUG 
dstrc_enter("fluid_tcons_tu");
#endif

/*----------------------------------------------------- check algorithm */
if (fdyn->iop==4)   /* one step theta */
{
    dynvar->dta  = fdyn->dt*0.05;
    dynvar->thsl = fdyn->dt*0.05*fdyn->theta;
    dynvar->thpl = dynvar->thsl;
    dynvar->thsr = (ONE - fdyn->theta)*fdyn->dt*0.05;
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
} /* end of fluid_tcons_tu*/ 


#endif
