/*!---------------------------------------------------------------------
\file
\brief The 2D fluid element used for projection methods

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/

#ifndef CCADISCRET
#ifdef D_FLUID2_PRO

/*!---------------------------------------------------------------------
\brief fluid2

<pre>                                                         genk 03/02

In this structure all variables used for element evaluation by the 2D
fluid element fluid2_pro are stored.

</pre>

--------------------------------------------------------------------------*/
typedef struct _FLUID2_PRO
{

  INT                ntyp;     /*!< flag for element type: 1=quad; 2=tri    */
  INT                nGP[2];   /*!< number of gaussian points in rs direct. */
  INT                is_ale;   /*!< flag whether there is ale to me or not  */
  enum   _DISMODE    dm;

  /* discontinuous pressure values stored at the element */
  DOUBLE* press;                /* pressure values at (n+1) */
  DOUBLE* pressm;               /* pressure values at (n) */
  DOUBLE* phi;                  /* pressure increments */
  INT*    dof;                  /* global dof numbers */
  INT*    ldof;                 /* processor local dof numbers */

  /*----- flag, if there is a lif&drag or fsi coupling line to this element */
  INT                force_on;

  ELEMENT* other;
} FLUID2_PRO;

#endif
#endif
