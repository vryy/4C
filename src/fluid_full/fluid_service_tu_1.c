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
\param	 place         INT	      (i)    place in sol_incr.
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*kapomegarat   DOUBLE	      (o)    kapomega  conv. ratio
\param	*fdyn	       FLUID_DYNAMIC	     	
\param lower_limit_kappa  DOUBLE	      (i)    lower limit for kappa
\param lower_limit_omega  DOUBLE	      (i)    lower limit for omega
\return void 

------------------------------------------------------------------------*/
void fluid_result_incre_tu_1(FIELD       *actfield,    
                           INTRA         *actintra,   
			         DIST_VECTOR   *sol,        
                           INT            place,      
			         SPARSE_ARRAY  *sysarray,      
			         SPARSE_TYP    *sysarray_typ,
			         DOUBLE        *kapomegarat,        
		               FLUID_DYNAMIC *fdyn,
                           DOUBLE         lower_limit_kappa,
                           DOUBLE         lower_limit_omega         
		              )
{
INT      i,j;
INT      max;
INT      diff;
INT      dof;
INT      numeq_total;
DOUBLE  *result;
DOUBLE   dkapomenorm=ZERO;
DOUBLE    kapomenorm=ZERO;
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
	       if(dynvar->kapomega_flag==0)
             {
              if (result[dof]<lower_limit_kappa) result[dof]=lower_limit_kappa;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j];
              dkapomenorm = DMAX(dkapomenorm, FABS(result[dof]-actnode->sol_increment.a.da[place][j]));
	        kapomenorm  = DMAX(kapomenorm, FABS(result[dof]));
              actnode->sol_increment.a.da[place][j] = result[dof];
             }
             if(dynvar->kapomega_flag==1)
             {
              if (result[dof]<lower_limit_omega) result[dof]=lower_limit_omega;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j+2];
              dkapomenorm = DMAX(dkapomenorm,FABS(result[dof]-actnode->sol_increment.a.da[place][j+2]));
	        kapomenorm  = DMAX(kapomenorm, FABS(result[dof]));
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
 	       if(dynvar->kapomega_flag==0)
             {
              if (result[dof]<lower_limit_kappa) result[dof]=lower_limit_kappa;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j];
              dkapomenorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j]);
	        kapomenorm  += FABS(result[dof]);
              actnode->sol_increment.a.da[place][j] = result[dof];
             }
             if(dynvar->kapomega_flag==1)
             {
              if (result[dof]<lower_limit_omega) result[dof]=lower_limit_omega;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j+2];
              dkapomenorm += FABS(result[dof]-actnode->sol_increment.a.da[place][j+2]);
	        kapomenorm  += FABS(result[dof]);
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
             if(dynvar->kapomega_flag==0)
             { 
              if (result[dof]<lower_limit_kappa) result[dof]=lower_limit_kappa;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j];
              dkapomenorm += pow(result[dof]-actnode->sol_increment.a.da[place][j],2);
	        kapomenorm  += pow(result[dof],2);
              actnode->sol_increment.a.da[place][j] = result[dof];   
             }
             if(dynvar->kapomega_flag==1)
             { 
              if (result[dof]<lower_limit_omega) result[dof]=lower_limit_omega;
              result[dof]=0.5*result[dof]+0.5*actnode->sol_increment.a.da[place][j+2];
              dkapomenorm += pow(result[dof]-actnode->sol_increment.a.da[place][j+2],2);
	        kapomenorm  += pow(result[dof],2);
              actnode->sol_increment.a.da[place][j+2] = result[dof];   
             }
          } /* end of loop over dofs */        
      } /* end of loop over nodes */
      dkapomenorm = sqrt(dkapomenorm);
       kapomenorm = sqrt( kapomenorm);
   break; 
   /*-------------------------------------------------------------------*/   
   default:
      dserror("unknown norm for convergence check!\n");
   } /* end of switch(fdyn->itnorm) */  
   /*------------------------------------------- check for "ZERO-field" */
   if (kapomenorm<EPS10)
   {
      kapomenorm = ONE;
      printf("ATTENTION: zero kapeps field - norm <= 1.0e-10 set to 1.0!! \n");
   }
/*------------------------------------- set final convergence ratios */
   *kapomegarat = dkapomenorm/kapomenorm;
/*----------------------------------------------------------------------*/
amdel(&result_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_result_incre_tu_1 */

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
\param *lower_limit_kappa  DOUBLE (o) lower limit for kappa  
\param *lower_limit_omega  DOUBLE (o) lower limit for omega  

\return void                                                                             

------------------------------------------------------------------------*/
void fluid_set_check_tu_1(
                       FIELD  *actfield, 
                       FLUID_DYNAMIC *fdyn, 
                       DOUBLE lower_limit_kappa,
                       DOUBLE lower_limit_omega
                       )
{
INT        i,j;
INT        numnp_total;              /* total number of fluid nodes     */
INT        numele_total;             /* total number of fluid elements  */
INT        numdf;	                   /* number of fluid dofs       	*/
DOUBLE     k_2;
GNODE     *actgnode;	             /* actual GNODE		            */
NODE      *actnode;	             /* actual NODE		            */

#ifdef DEBUG 
dstrc_enter("fluid_set_check_tu");
#endif 

/*----------------------------------------------------- set some values */
numnp_total  = actfield->dis[1].numnp;
numele_total = actfield->dis[1].numele;
numdf        = 1;

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
      actnode->sol_increment.a.da[3][j+2] = lower_limit_omega*10000;
      actnode->sol_increment.a.da[1][j+2] = actnode->sol_increment.a.da[3][j+2];
      actnode->sol_increment.a.da[3][j+1] = actnode->sol_increment.a.da[3][j]/actnode->sol_increment.a.da[3][j+2];
      actnode->sol_increment.a.da[2][j+1] = actnode->sol_increment.a.da[3][j+1];
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
void fluid_eddy_update_1(FIELD         *actfield, 
                         DIST_VECTOR   *sol   
                        )
{
INT  i,j;
INT  dof,actmat;
INT  numeq_total;
DOUBLE  visc,Re_t;
NODE    *actnode;
ELEMENT *actele;

#ifdef DEBUG 
dstrc_enter("fluid_eddy_update_1");
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
          
          actnode->sol_increment.a.da[3][1] = actnode->sol_increment.a.da[3][0]/
                                              actnode->sol_increment.a.da[3][2];
         }
        }
      }   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of  fluid_eddy_update_1 */


/*!---------------------------------------------------------------------                                         
\brief storing results for characteristic lenght in solution history
and calculate norms for the iteration convergence check

<pre>                                                        he  12/02

</pre>   
\param **actfield      FIELD	            (i)    actual field       
\param	*sol 	     DIST_VECTOR        (i)    solution vector
\param  lenghtrat	     DOUBLE             (o)    conv. ratio
\param     *fdyn 	     FLUID_DYNAMIC      (i)   
\return void 

------------------------------------------------------------------------*/
void fluid_lenght_update_1(FIELD         *actfield, 
                          DIST_VECTOR   *sol,   
		              DOUBLE        *lenghtrat, 
                          FLUID_DYNAMIC *fdyn      
                         )
{
NODE    *actnode;
INT  i,j;
INT  dof;
DOUBLE   dlenghtnorm=ZERO;
DOUBLE    lenghtnorm=ZERO;
INT  numeq_total;

#ifdef DEBUG 
dstrc_enter("fluid_lenght_update_1");
#endif
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;  
  
      for (i=0; i<actfield->dis[1].numnp; i++)
      {
         actnode = &(actfield->dis[1].node[i]);

         for (j=0; j<actnode->numdf; j++) /* loop dofs  */
         {
	    dof = actnode->dof[j];
            if (dof>=numeq_total) continue;
             if(actnode->sol_increment.a.da[3][2]!=0)
             {
              actnode->sol_increment.a.da[3][3] = 0.5*sqrt(actnode->sol_increment.a.da[3][0])/ actnode->sol_increment.a.da[3][2]
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
} /* end of  fluid_lenght_update_1 */

#endif
