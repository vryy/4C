#include "../headers/standardtypes.h"
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
 | input of conditions                                    m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_conditions()
{
int  ierr;
int  i;
#ifdef DEBUG 
dstrc_enter("inp_conditions");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------ input of time curves */
inp_cond_curve();
/*----- input all element and nodal conditions as Dirich / Neumann .... */
for (i=0; i<genprob.numfld; i++)
{
   if (field[i].fieldtyp==structure)
   {
      inp_cond_nodal_struct(&(field[i]));
      inp_couple(&(field[i]));
      inp_cond_ele_struct(&(field[i]));
   }
   if (field[i].fieldtyp==fluid)
   {
      inp_cond_nodal_fluid(&(field[i]));
      inp_couple(&(field[i]));
      /*inp_cond_ele_fluid(&(field[i])); not yet implemented */
   }
   if (field[i].fieldtyp==ale)
   {
      inp_cond_nodal_ale(&(field[i]));
      inp_couple(&(field[i]));
      /*inp_cond_ele_ale(&(field[i])); not yet implemented */
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_conditions */



