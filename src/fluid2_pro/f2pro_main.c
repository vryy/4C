/*!----------------------------------------------------------------------
\file
\brief main routine fluid2_pro element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2_PRO 
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2pro_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             genk 04/02    |
 *----------------------------------------------------------------------*/

/*!---------------------------------------------------------------------
\brief main fluid2_pro control routine

                                                        basol 10/02
</pre>
\param	*action	         CALC_ACTION    (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *lmass_global    ARRAY	        (o)   lumped mass matrix
\param  *gradopr_global  ARRAY	        (o)   gradient operator
\param  *etforce_global  ARRAY	        (o)   element time force 
\param  *etforce1_global ARRAY	        (o)   element time force part1
\param  *etforce2_global ARRAY	        (o)   element time force part2
\param  *etforce3_global ARRAY	        (o)   element time force part3
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param	 double	         dt             (i)   incremental time step

\return void

------------------------------------------------------------------------*/
void fluid2_pro(     PARTITION     *actpart,
                     INTRA         *actintra,
		     ELEMENT       *elev,
		     ELEMENT       *elep,
		     ARRAY         *estif_global,
		     ARRAY         *emass_global,
		     ARRAY         *lmass_global,
		     ARRAY         *gradopr_global,
		     ARRAY         *etforce_global,
   	             ARRAY         *eiforce_global,
		     ARRAY         *edforce_global,
		     ARRAY         *gforce_global,
		     CALC_ACTION   *action,
                     INT           *hasdirich
	       )
{
/*----------------------------------------------------------------------*/
#ifdef D_FLUID2_PRO 
/*----------------------------------------------------------------------*/
MATERIAL              *actmat;     /* actual material                   */
static FLUID_DATA      *data;      
FLUID_DYNAMIC         *fdyn;
static FLUID_DYN_CALC  *dynvar;
FIELD                 *actfield;   /* actual field                      */

#ifdef DEBUG 
dstrc_enter("fluid2_pro");
#endif
/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------------------ initialisation */
case calc_fluid_init:
   dynvar = &(alldyn[genprob.numff].fdyn->dynvar);
   data   = &(alldyn[genprob.numff].fdyn->dynvar.data);
   /*---------------------------------------- init the element routines */   
   f2_intg(data,0); 
   /*----------- f2_pro element is initialized, some arrays are defined */
   f2pro_calele(data,dynvar,NULL,NULL,estif_global,emass_global,
               lmass_global,gradopr_global,etforce_global,eiforce_global,
	       edforce_global,NULL,1);                                        
break;                    
/*---------------------------------------- calculate the A=(Ct)(Ml-1)(C)*/
case calc_fluid_amatrix:
   f2pro_calele(data,dynvar,elev,elep,estif_global,emass_global, 
                lmass_global,gradopr_global,etforce_global,eiforce_global,
		edforce_global,hasdirich,0);
break;                        
/*------------------------------------------- call the element routines */
case calc_fluid_f2pro:
   f2pro_calele(data,dynvar,elev,elep,estif_global,emass_global,
                lmass_global,gradopr_global,etforce_global,eiforce_global,
		edforce_global,hasdirich,0);   
break;                        
/*------------------------------------------- call the element routines */
case calc_fluid_f2pro_rhs_both:
   f2pro_calele(data,dynvar,elev,elep,estif_global,emass_global,
                lmass_global,gradopr_global,etforce_global,eiforce_global,
		edforce_global,hasdirich,0); 
break; 
/*----------------------------------------------------------------------*/
default:
   dserror("action unknown\n");
break;
} /* end swtich (*action) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
return; 
} /* end of fluid2_pro */
/*! @} (documentation module close)*/
