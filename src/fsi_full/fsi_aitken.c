/*!----------------------------------------------------------------------
\file
\brief AITKEN iteration

*----------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fsi_prototypes.h"    
#ifdef D_FSI
static DOUBLE done = 1.0E20;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;  
/*!---------------------------------------------------------------------*
\brief compute relaxation parameter via AITKEN iteration

<pre>                                                         genk 01/03

RELAX(i) = 1 - NU(i) (NU(i) at time (n+1))
NU(i) = N(i-1) + [NU(i-1) - 1)] * TOP/DEN
TOP = [DELTAd(i) - DETLAd(i+1)] * DELTAd(i+1)
DEN = [DELTAd(i) - DETLAd(i+1)]^2
DELTAd(i) = d(i-1) - d~(i)
d(i-1) ... relaxad interface displacements of last iteration step
d~(i)  ... unrelaxed interface displacements of this iteration step

see also Dissertation of D.P.MOK, chapter 6.4 
		     
</pre>

\param *structfield   FIELD	     (i)   structural field
\param *fsidyn 	      FSI_DYNAMIC    (i)   
\param  itnum         INT            (i)   actual iteration step
\return void     

------------------------------------------------------------------------*/
void fsi_aitken(  
                 FIELD          *structfield, 
                 FSI_DYNAMIC    *fsidyn, 
		 INT             itnum
	       )
{
INT            i,j;           /* some counters                          */
INT            numdf;         /* number of nodal dofs                   */
INT            dof;           /* actual dof                             */
static INT     numnp_total;   /* total number of nodes                  */
static INT     numdf_total;   /* total number of coupled dofs           */
static INT    *sid;           /* flags for structural interface dofs    */
static DOUBLE  nu;            /* actual AITKEN factor                   */
DOUBLE         top=ZERO;      
DOUBLE         den=ZERO;
DOUBLE         del2;				
DOUBLE       **sol_mf;        /* multifield nodal solution array        */
static ARRAY   del_a;      
static DOUBLE *del;           /* DELTAd(i)                              */
NODE          *actsnode;      /* actual structural node                 */

#ifdef DEBUG 
dstrc_enter("fsi_aitken");
#endif

if (itnum==0)
{
   if (fsidyn->step==1)   
   { 
      nu          = ZERO; 
      numnp_total = structfield->dis[0].numnp;
      sid         = fsidyn->sid.a.iv;
      numdf_total = fsidyn->sid.fdim;
      del         = amdef("del",&del_a,numdf_total,1,"DV");
      aminit(&del_a,&done);
   }
   else 
      aminit(&del_a,&done);  
}

/*----------------------- the solution history sol_mf of the structfield:
- sol_mf[0][j] holds the latest struct-displacements
- sol_mf[1][j] holds the (relaxed) displacements of the last iteration step

/*---------------------------------------------------------- loop nodes */
for (i=0;i<numnp_total;i++)
{
   actsnode  = &(structfield->dis[0].node[i]);
   numdf = actsnode->numdf; 
   sol_mf = actsnode->sol_mf.a.da;
   /*--------------------------------- loop dofs and check for coupling */
   for (j=0;j<numdf;j++)
   {
      dof = actsnode->dof[j];      
      dsassert(dof<numdf_total,"dofnumber not valid!\n");
      if (sid[dof]==0) continue;
      del2     = del[dof]; 
      del[dof] = sol_mf[1][j] - sol_mf[0][j];
      del2     = del2 - del[dof];
      top     += del2*del[dof];
      den     += del2*del2;        
   } /* end of loop over dofs */   
} /* end of loop over nodes */

nu = nu + (nu - ONE)*top/den;

/*---------------------------------------------------------- finalising */
fsidyn->relax = ONE - nu;
/*------------------------------------------------- output to the screen */
if (par.myrank==0)
{
printf("\n");
printf("AITKEN ITERATION: RELAX = %.5lf\n",fsidyn->relax);
printf("\n");
}

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fsi_aitken */

#endif

/*! @} (documentation module close)*/

