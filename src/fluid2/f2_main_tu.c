/*!----------------------------------------------------------------------
\file
\brief main routine fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
/*----------------------------------------------------------------------*/
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

<pre>                                                        he    11/02
</pre>
\param      *actpart	        PARTITION    (i)	    
\param	*actintra	        INTRA        (i)
\param	*eleke		  ELEMENT      (i)    actual element
\param	*elev		        ELEMENT      (i)    element with velocities
\param	*estif_global    ARRAY           (o)    element stiffness matrix
\param      *emass_global    ARRAY           (o)    element mass matrix
\param	*etforce_global  ARRAY           (o)    element time force vector
\param      *eiforce_global  ARRAY           (o)    element iter force vector
\param      *edforce_global  ARRAY           (o)    element dirich force vector
\param      *eproforce_global   ARRAY        (o)    element production force vector
\param	*action	        CALC_ACTION  (i)
\param	*hasdirich	        INT          (o)    flag
\param      *hasext             INT          (o)    flag
\param      *container	        CONTAINER    (i)	    
\return void

------------------------------------------------------------------------*/
void fluid2_tu(
            PARTITION   *actpart,
            INTRA       *actintra,
            ELEMENT     *eleke,  
            ELEMENT     *elev,            
            ARRAY       *estif_global,   
            ARRAY       *emass_global,   
	      ARRAY       *etforce_global, 
	      ARRAY       *eiforce_global, 
            ARRAY       *edforce_global, 
            ARRAY       *eproforce_global, 
            CALC_ACTION *action,
	      INT         *hasdirich,
	      INT         *hasext,
            CONTAINER   *container
	   )
{
/*----------------------------------------------------------------------*/

static INT              numff;      /* actual number of fluid field     */
MATERIAL               *actmat;     /* actual material                  */
static FLUID_DATA      *data;      
FLUID_DYNAMIC          *fdyn;
FIELD                  *actfield;   /* actual field                     */
INT                    start;
#ifdef DEBUG 
dstrc_enter("fluid2_tu");
#endif

/*--------------------------------------------------------------------- */
/*------------------------------------------------- KAPPA-EPSILON MODEL */
if(container->turbu == 2)
{
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
   data   = alldyn[numff].fdyn->data;
/*------------------------------------------- init the element routines */   
   f2_intg(data,0);
   f2_calele_tu(data,NULL,NULL,
                estif_global,emass_global,
	        etforce_global,eiforce_global,edforce_global,eproforce_global,
	        NULL,NULL,1);
break;
/*------------------------------------------- call the element routines */
case calc_fluid:
   f2_calele_tu(data,eleke,elev,
                estif_global,emass_global,
	          etforce_global,eiforce_global,edforce_global,eproforce_global,
	          hasdirich,hasext,0);
break;

/*----------------------------------------------------------------------*/
default:
   dserror("action unknown\n");
break;
} /* end swtich (*action) */

}

/*--------------------------------------------------------------------- */
/*--------------------------------------------------- KAPPA-OMEGA MODEL */
if(container->turbu == 3)
{
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
   data   = alldyn[numff].fdyn->data; 
/*------------------------------------------- init the element routines */   
   f2_intg(data,0);
   f2_calele_tu_1(data,NULL,NULL,
                  estif_global,emass_global,
	            etforce_global,eiforce_global,edforce_global,eproforce_global,
	            NULL,NULL,1);
break;
/*------------------------------------------- call the element routines */
case calc_fluid:
   f2_calele_tu_1(data,eleke,elev,
                  estif_global,emass_global,
	            etforce_global,eiforce_global,edforce_global,eproforce_global,
	            hasdirich,hasext,0);
break;

/*----------------------------------------------------------------------*/
default:
   dserror("action unknown\n");
break;
} /* end swtich (*action) */

}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return; 
} /* end of fluid2 */


#endif
