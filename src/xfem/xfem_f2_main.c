#ifdef D_XFEM 
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "xfem_prototypes.h"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA            *alldyn;   
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





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void xfem_fluid2(
  PARTITION   *actpart,
  INTRA       *actintra,
  ELEMENT     *ele,             
  ARRAY       *estif_global,   
  ARRAY       *emass_global,   
  ARRAY       *etforce_global, 
  ARRAY       *eiforce_global, 
  ARRAY       *edforce_global, 
  CALC_ACTION *action,
  INT         *hasdirich,
  INT         *hasext,
  CONTAINER   *container       
  )
{
  static INT                 viscstr;
  static FLUID_DATA         *data;      
  static FLUID_DYNAMIC      *fdyn;
  
#ifdef DEBUG
  dstrc_enter("xfem_fluid2");
#endif
/*----------------------------------------------------------------------*/

  switch (*action)
  {
       /* initialize */
      case calc_fluid_init:
        /* access to the number of fluid field */
        fdyn  = alldyn[genprob.numff].fdyn;
        data  = alldyn[genprob.numff].fdyn->data;
        viscstr = alldyn[genprob.numff].fdyn->viscstr;
        xfem_f2_intg(data);
        /* init the element routines */   
        xfem_f2_calele(
          data,NULL,estif_global,emass_global,
          etforce_global,eiforce_global,edforce_global,
          NULL,NULL,0,0,1
          );
        f2_iedg(NULL,ele,-1,1);
        break;
      /* call the element routines */
      case calc_fluid:
        xfem_f2_calele(
          data,ele,estif_global,emass_global,
          etforce_global,eiforce_global,edforce_global,hasdirich,
          hasext,actintra->intra_rank,container->is_relax,0
          );
        break;
      default:
        dserror("action unknown\n");
        break;
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return; 
}  /* end of xfem_f2_main */
#endif
