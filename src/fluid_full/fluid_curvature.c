/*!----------------------------------------------------------------------
\file
\brief curvature at free surface 

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID 
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];
/*!---------------------------------------------------------------------
\brief curvature at free surface for fluid elements

<pre>                                                         genk 02/03

\param  *actfield   FIELD       (i)    active field 
\param  *actpart    PARTITION   (i)    my partition to this field
\param  *actintra   INTRA       (i)    my intra-communicator
\param  *action     CALC_ACTION (i)    calculation option
		     
</pre>
\return void                                                                       

------------------------------------------------------------------------*/
void fluid_curvature(FIELD        *actfield,    
                     PARTITION    *actpart,     
                     INTRA        *actintra,    
                     CALC_ACTION  *action)    
{
INT i;
ELEMENT *actele;

#ifdef DEBUG 
dstrc_enter("fluid_curvature");
#endif
/* =======================================================call elements */
/*---------------------------------------------- loop over all elements */
for (i=0; i<actpart->pdis[0].numele; i++)
{   
  /*------------------------------------ set pointer to active element */
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)/*======================= call element routines */
   {
   case el_fluid2: 
#ifdef D_FLUID2   
      if (actele->e.f2->fs_on==0) break;
      fluid2(actpart,actintra,actele,NULL,
             NULL,NULL,
             NULL,NULL,NULL,
	     action,NULL,NULL,NULL);
#endif
   break;
   default:
      dserror("curvature calculation only possible for fluid2!\n");
   }      
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_curvature */
#endif
/*! @} (documentation module close)*/
