/*-----------------------------------------------------------------------*/
/*!
\file
\brief Brief description.

  Very detailed description.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3/fluid3.h"

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

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;




/*-----------------------------------------------------------------------*/
/*!
  \brief main fluid3 control routine


  \param actpart        *PARTITION (i) current partition
  \param actintra       *INTRA    (i) current communicator
  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param estif          *DOUBLE   (i) element stiffness matrix
  \param emass          *DOUBLE   (i) element mass matrix
  \param eforce         *DOUBLE   (i) element force
  \param edforce        *DOUBLE   (i) element dirichlet force
  \param action         *CALC_ACTION (i) what should I do
  \param hasdirich      *INT      (i) flag if s.th. was written to edforce
  \param hasext         *INT      (i) flag if there are ext forces
  \param container      *CONTAINER (i) container
  \param loop            INT      (i) num of eles in ele[]

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void fluid3_fast(
    PARTITION   *actpart,
    INTRA       *actintra,
    ELEMENT     *ele[LOOPL],
    ARRAY       *estif_fast,
    ARRAY       *emass_fast,
    ARRAY       *eforce_fast,
    ARRAY       *edforce_fast,
    CALC_ACTION *action,
    INT         *hasdirich,
    INT         *hasext,
    CONTAINER   *container,
    INT          loop
    )
{
  static INT              viscstr;
  FLUID_DYNAMIC          *fdyn;
  INT                     i,l,ldflag;
  GVOL     *actgvol;
  GSURF    *actgsurf;
  DSURF    *actdsurf;
  ARRAY_POSITION* ipos;

#ifdef DEBUG
  dstrc_enter("fluid3_fast");
#endif

  ipos = &(field[genprob.numff].dis[container->disnum].ipos);

  /* switch to do option */
  switch (*action)
  {

    /* initialisation */
    case calc_fluid_init:
      /* find number of fluid field */
      fdyn   = alldyn[genprob.numff].fdyn;
      viscstr= alldyn[genprob.numff].fdyn->viscstr;

      /* init the element routines */
      f3_intg(0);

      f3fcalele(NULL, estif_fast,
          emass_fast,eforce_fast,edforce_fast,ipos,
          NULL,NULL,1,loop);
      break;


      /* call the element routines */
    case calc_fluid:
      f3fcalele(ele, estif_fast,
          emass_fast,eforce_fast,edforce_fast,ipos,
          hasdirich,hasext,0,loop);
      break;


      /* calculate fluid stresses */
    case calc_fluid_stress:

      /* calculate stresses only for elements belonging to this proc */
      /*if (par.myrank==ele[0]->proc)*/
        f3fstress(container->str,viscstr,ele,ipos,loop);
      break;


      /* calculate fluid stresses for lift&drag calculation */
    case calc_fluid_liftdrag:
      if (ele[0]->proc == actintra->intra_rank)
      {
        for (l=0;l<loop;l++)
        {
          ldflag=0;
          /* check if element is on liftdrag-dline */
          actgvol=ele[l]->g.gvol;
          for (i=0;i<actgvol->ngsurf;i++)
          {
            actgsurf=actgvol->gsurf[i];
            actdsurf=actgsurf->dsurf;
            if (actdsurf==NULL) continue;
            if (actdsurf->liftdrag==NULL) continue;
            ldflag++;
            break;
          }
          if (ldflag>0)
          {
            f3fstress(container->str,viscstr,ele,ipos,loop);
            f3fliftdrag(ele[l],container);
          }
        }
      }
      break;


    case calc_fluid_error:
      f3f_calerror(ele,hasext,container,ipos,loop);
      break;



    case write_restart:
      f3f_write_restart(ele[0],container->handsize,container->handles);
      break;


    case read_restart:
      f3f_read_restart(ele[0],container->handsize,container->handles);
      break;


    default:
      dserror("action unknown\n");
      break;
  } /* end switch (*action) */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of fluid3_fast */


#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/


