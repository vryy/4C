/*!----------------------------------------------------------------------
\file
\brief main routine fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
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
\brief main fluid2 control routine

<pre>                                                         genk 03/02
</pre>
\param  *actpart	 PARTITION    (i)	    
\param	*actintra	 INTRA        (i)
\param	*ele		 ELEMENT      (i)    actual element
\param	*estif_global	 ARRAY        (o)    element stiffness matrix
\param  *emass_global	 ARRAY        (o)    element mass matrix
\param	*etforce_global  ARRAY        (o)    element time force vector
\param  *eiforce_global  ARRAY        (o)    element iter force vecotr
\param	*edforce_global  ARRAY        (o)    ele dirichl. force vector
\param	*action	         CALC_ACTION  (i)
\param	*hasdirich	 int          (o)    flag
\param  *hasext          int          (o)    flag
\param  *container       CONTAINER    (i)    container
\return void

------------------------------------------------------------------------*/
void fluid2(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,             
            ELEMENT     *eleke,             
            ARRAY       *estif_global,   
            ARRAY       *emass_global,   
	    ARRAY       *etforce_global, 
	    ARRAY       *eiforce_global, 
	    ARRAY       *edforce_global, 
            CALC_ACTION *action,
	    int         *hasdirich,
	    int         *hasext,
            CONTAINER   *container       
	   )
{
/*----------------------------------------------------------------------*/
#ifdef D_FLUID2 
/*----------------------------------------------------------------------*/

static int              numff;      /* actual number of fluid field     */
static int              viscstr;
static FLUID_DATA      *data;      
static FLUID_DYN_CALC  *dynvar;
FIELD                  *actfield;   /* actual field                     */
#ifdef DEBUG 
dstrc_enter("fluid2");
#endif

/*------------------------------------------------- switch to do option */
switch (*action)
{
/*------------------------------------------------------ initialisation */
case calc_fluid_init:
/* ----------------------------------------- find number of fluid field */
   for (numff=0;numff<genprob.numfld;numff++)
   {
      actfield=&(field[numff]);
      if (actfield->fieldtyp==fluid)
      break;
   } /* end loop over numff */
   dynvar = &(alldyn[numff].fdyn->dynvar);
   data   = &(alldyn[numff].fdyn->dynvar.data);
   viscstr= alldyn[numff].fdyn->viscstr;
/*------------------------------------------- init the element routines */   
   f2_intg(data,0);
   f2_calele(data,dynvar,NULL,NULL,
             estif_global,emass_global,
	     etforce_global,eiforce_global,edforce_global,
	     NULL,NULL,0,1);
   f2_iedg(NULL,ele,-1,1);
break;

case calc_fluid_initvort:
/* ----------------------------------------- find number of fluid field */
   for (numff=0;numff<genprob.numfld;numff++)
   {
      actfield=&(field[numff]);
      if (actfield->fieldtyp==fluid)
      break;
   } /* end loop over numff */
   dynvar = &(alldyn[numff].fdyn->dynvar);
   data   = &(alldyn[numff].fdyn->dynvar.data);
/*------------------------------------------- init the element routines */
   f2_intg(data,0); 
   f2_calvort(data,dynvar,ele,1);  
break;

/*------------------------------------------- call the element routines */
case calc_fluid:
   f2_calele(data,dynvar,ele,eleke,
             estif_global,emass_global,
	     etforce_global,eiforce_global,edforce_global,
	     hasdirich,hasext,actintra->intra_rank,0);
break;

/*------------------------------------------- calculate fluid vorticity */
case calc_fluid_vort:
   f2_calvort(data,dynvar,ele,0);
break;

/*-------------------------------------------- calculate fluid stresses */
case calc_fluid_stress:
    f2_stress(container->str,viscstr,data,ele);
break;

/*--------------------------------- calculate curvature at free surface */
case calc_fluid_curvature:
   f2_curvature(data,dynvar,ele,actintra->intra_rank);
break;

/*---------------------------------------- calculate the shear stresses */
case calc_fluid_shearvelo:
   f2_shearstress(ele,dynvar);
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
} /* end of fluid2 */
/*! @} (documentation module close)*/
