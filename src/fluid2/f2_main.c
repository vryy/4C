/*!----------------------------------------------------------------------
\file
\brief main routine fluid2 element

------------------------------------------------------------------------*/
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
\return void

------------------------------------------------------------------------*/
void fluid2(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *ele,             
            ARRAY       *estif_global,   
            ARRAY       *emass_global,   
	    ARRAY       *etforce_global, 
	    ARRAY       *eiforce_global, 
	    ARRAY       *edforce_global, 
            CALC_ACTION *action,
	    int         *hasdirich       
	   )
{
/*----------------------------------------------------------------------*/
#ifdef D_FLUID2 
/*----------------------------------------------------------------------*/

static int              numff;      /* actual number of fluid field     */
MATERIAL               *actmat;     /* actual material                  */
static FLUID_DATA      *data;      
FLUID_DYNAMIC          *fdyn;
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
/*------------------------------------------- init the element routines */   
   f2_intg(data,0);
   f2_calele(data,dynvar,NULL,
             estif_global,emass_global,
	     etforce_global,eiforce_global,edforce_global,
	     NULL,1);
break;

/*------------------------------------------- call the element routines */
case calc_fluid:
   f2_calele(data,dynvar,ele,
             estif_global,emass_global,
	     etforce_global,eiforce_global,edforce_global,
	     hasdirich,0);
break;

/*------------------------------------------- calculate fluid vorticity */
case calc_fluid_vort:
   dserror("calculation of vorticity not implemented yet!\n");
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
