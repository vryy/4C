/*!----------------------------------------------------------------------
\file
\brief main routine for fluid2_xfem

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM 
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "xfem_prototypes.h"
/*! 
\addtogroup XFEM 
*//*! @{ (documentation module open)*/



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



/*!---------------------------------------------------------------------
\brief main fluid2_xfem control routine

<pre>                                                            irhan 05/04
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
\param	*hasdirich	 INT          (o)    flag
\param  *hasext          INT          (o)    flag
\param  *container       CONTAINER    (i)    container
\return void

------------------------------------------------------------------------*/
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
/*! @} (documentation module close)*/
#endif
