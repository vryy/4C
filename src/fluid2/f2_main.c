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

/*----------------------------------------------------------------------*
 | main fluid2  control routine                              genk 03/02 |
 *----------------------------------------------------------------------*/
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
#ifdef D_FLUID2 
static int              numff;      /* number of fluid field */
double                 *intforce;
MATERIAL               *actmat;
static F2_DATA         *data;
FLUID_DYNAMIC          *fdyn;
static FLUID_DYN_CALC  *dynvar;
FIELD                  *actfield;

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
   }
   dynvar = &(alldyn[numff].fdyn->dynvar);
   data   = &(alldyn[numff].fdyn->dynvar.data.f2data);
/*------------------------------------------- init the element routines */   
   f2_intg(NULL,data,0);
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

/*----------------------------------------------------------------------*/
default:
   dserror("action unknown");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return; 
} /* end of fluid2 */
