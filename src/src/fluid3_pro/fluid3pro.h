/*!---------------------------------------------------------------------
\file
\brief The 3D fluid element used for projection methods

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/

#ifndef CCADISCRET
#ifdef D_FLUID3_PRO

/*!---------------------------------------------------------------------
\brief fluid3

<pre>                                                         genk 03/02

In this structure all variables used for element evaluation by the 2D
fluid element fluid3_pro are stored.

</pre>

--------------------------------------------------------------------------*/
typedef struct _FLUID3_PRO
{

  INT                ntyp;     /*!< flag for element type: 1=quad; 2=tri    */
  INT                nGP[3];   /*!< number of gaussian points in rs direct. */
  INT                is_ale;   /*!< flag whether there is ale to me or not  */
#if 0
  struct _ELEMENT   *my_ale;   /*!< pointer to my ale ele, otherwise NULL   */
#endif
  enum   _DISMODE    dm;

  /* discontinuous pressure values stored at the element */
  DOUBLE* press;                /* pressure values */
  DOUBLE* pressm;
  DOUBLE* phi;                  /* pressure increments */
  INT*    dof;                  /* global dof numbers */
  INT*    ldof;                 /* processor local dof numbers */

  /*----- flag, if there is a lif&drag or fsi coupling line to this element */
  INT                force_on;

  ELEMENT* other;
#ifdef QUASI_NEWTON
  struct _ARRAY  estif;
#endif
} FLUID3_PRO;

#endif
#endif
