/*!----------------------------------------------------------------------
\file
\brief relaxation of structural interface displacements

*----------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fsi_prototypes.h"    
/*!---------------------------------------------------------------------
\brief relaxation of structural interface displacements

<pre>                                                         genk 01/03

Relaxation of structural interface displacements:

   d(i+1) = RELAX(i) * d~(i+1) + (1 - RELAX(i)) * d(i) 

   d(i+1)  ... new relaxed interface displacements
   d~(i+1) ... unrelaxed interface displacements
   d(i)    ... old relaxed interface displacements 
               (they are stored in sol_increment of the ALE-nodes)

   see dissertation of D.P. MOK chapter 6.2
   
   d~(i+1) = actsnode->sol_mf.a.da[0][j]
   d(i)    = actsnode->sol_mf.a.da[1][j]
      
   result is written to:	       
             actsnode->sol_mf.a.da[0][j]	       
		     
</pre>
\param *structfield   FIELD	     (i)   structural field
\param *fsidyn 	      FSI_DYNAMIC    (i)   
\return void                                                                             

------------------------------------------------------------------------*/
void fsi_relax_intdisp(
                            FIELD          *structfield, 
			    FSI_DYNAMIC    *fsidyn
		      )
{
int     i,j;              /* simply some counters		        */
int     numdf;            /* actual number of dofs		        */
int     numnp_total;      /* number of struct nodes		        */
int     numdf_total;      /* total number of struct dofs                */
int     dof;              /* actual dof                                 */
int    *sid;              /* structural interface dofs                  */
double  relax;            /* actual relaxation parameter omega	        */
double  fac;             
NODE   *actsnode;	  /* the actual struct	node 			*/

#ifdef DEBUG 
dstrc_enter("fsi_relax_intdisp");
#endif

/*----------------------------------------------------- set some values */
relax       = fsidyn->relax;
fac         = ONE-relax;
numnp_total = structfield->dis[0].numnp;
sid         = fsidyn->sid.a.iv;
numdf_total = fsidyn->sid.fdim;

/*---------------------------------------------------------- loop nodes */
for (i=0;i<numnp_total;i++)
{
   actsnode  = &(structfield->dis[0].node[i]);
   numdf = actsnode->numdf; 
   /*--------------------------------- loop dofs and check for coupling */
   for (j=0;j<numdf;j++)
   {
      dof = actsnode->dof[j];
      dsassert(dof<numdf_total,"dofnumber not valid!\n");
      if (sid[dof]==0) continue;
      actsnode->sol_mf.a.da[0][j] = relax*actsnode->sol_mf.a.da[0][j] 
                                  + fac*actsnode->sol_mf.a.da[1][j];
   } /* end of loop over dofs */   
} /* end of loop over nodes */

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fsi_relax_intdisp */

#endif

/*! @} (documentation module close)*/
