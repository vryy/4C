/*======================================================================*/
/*!
\file
\brief Select proper material law

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>

\author mf
\date 10/06
*/
#ifndef CCADISCRET
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"
#ifdef D_MAT
#include "../materials/mat_prototypes.h"
#endif

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/07
*/
extern GENPROB genprob;


/*======================================================================*/
/*!
\brief Allocate internal variables due to non-elastic material laws
\author bborn
\date 03/07
*/
void so3_mat_init(PARTITION* part,  /*!< partition */
                  const MATERIAL* mat)  /*!< material */
{
  INT jdis;  /* discretisation loop jndex */
  INT iele;  /* element loop index */
  ELEMENT* actele;  /* pointer to current element */
  const MATERIAL* actmat;  /* material of element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_init");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate material internal variables stored in each element */
  /* loop over all discretisations of partition thermal field */
  for (jdis=0; jdis<part->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<part->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = part->pdis[jdis].element[iele];
      /* check if SOLID3 element */
      if (actele->eltyp == el_solid3)
      {
        /* get element material */
        actmat = &(mat[actele->mat-1]);
        /* switch according to element material */
        switch (actmat->mattyp)
        {
          /* St. Venant-Kirchhoff material */
          case m_stvenant:
            /* step over */
            break;
          /* Robison's material */
#ifdef D_TSI
          case m_vp_robinson:
            so3_mat_robinson_init(actele);
            break;
#endif
          /* catch the remainders */
          default:
            /* do nothing */
            break;
        }
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end so3_mat_init() */


/*======================================================================*/
/*!
\brief Deallocate internal variables due to non-elastic material laws
\author bborn
\date 03/07
*/
void so3_mat_final(PARTITION* part,  /*!< partition */
                   const MATERIAL* mat)  /*!< material */
{
  INT jdis;  /* discretisation loop jndex */
  INT iele;  /* element loop index */
  ELEMENT* actele;  /* pointer to current element */
  const MATERIAL* actmat;  /* material of element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_init");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocate material internal variables stored in each element */
  /* loop over all discretisations of partition thermal field */
  for (jdis=0; jdis<part->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<part->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = part->pdis[jdis].element[iele];
      /* check if SOLID3 element */
      if (actele->eltyp == el_solid3)
      {
        /* get element material */
        actmat = &(mat[actele->mat-1]);
        /* switch according to element material */
        switch (actmat->mattyp)
        {
          /* St. Venant-Kirchhoff material */
          case m_stvenant:
            /* step over */
            break;
          /* Robison's material */
#ifdef D_TSI
          case m_vp_robinson:
            so3_mat_robinson_final(actele);
            break;
#endif
          /* catch the remainders */
          default:
            /* do nothing */
            break;
        }
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end so3_mat_final() */



/*======================================================================*/
/*!
\brief Select proper material law

\param *container CONTAINER      (i)   container
\param *ele       ELEMENT        (i)   pointer to current element
\param *mat       MATERIAL       (i)   pointer to current material
\param ip         INT            (i)   current Gauss point index
\param *gds       SO3_GEODEFSTR  (i)   geom. & def. data at Gauss point
\param stress[]   DOUBLE         (o)   linear(Biot)/2.Piola-Kirchhoff stress
\param cmat[][]   DOUBLE         (o)   constitutive matrix
\return void

\author bborn
\date 03/07
*/
void so3_mat_sel(CONTAINER *container,
                 ELEMENT *ele,
                 MATERIAL *mat,
                 SO3_GPSHAPEDERIV *gpshade,
                 INT ip,
                 SO3_GEODEFSTR *gds,
                 DOUBLE stress[NUMSTR_SOLID3],
                 DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_sel");
#endif

  /*====================================================================*/
  /* ===> central material routines
   *
   *      These materials are supposed to be connected to the
   *      existant (or new?) central material routines.
   *      Right now, only the simple St.Venant-Kirchhoff material
   *      is included to test the element.
   *
   * ===> Do it, have fun. Or don't.
   */

  /*====================================================================*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    /*------------------------------------------------------------------*/
    /* St. Venant-Kirchhoff material */
    case m_stvenant:
      so3_mat_stvenant_sel(container, ele, mat, ip, gds, stress, cmat);
      break;
    /*------------------------------------------------------------------*/
    /* Robinson's visco-plastic temperature-dependent material */
#ifdef D_TSI
    case m_vp_robinson:
      so3_mat_robinson_sel(container, ele, mat->m.vp_robinson, gpshade,
                           ip, gds, stress, cmat);
      break;
#endif
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_mat_sel */


/*======================================================================*/
/*!
\brief Select material law to return stress (called for stress output)
\param *container CONTAINER      (i)   container
\param *ele       ELEMENT        (i)   pointer to current element
\param *mat       MATERIAL       (i)   pointer to current material
\param ip         INT            (i)   current Gauss point index
\param *gds       SO3_GEODEFSTR  (i)   geom. & def. data at Gauss point
\param stress[]   DOUBLE         (o)   linear(Biot)/2.Piola-Kirchhoff stress
\param cmat[][]   DOUBLE         (o)   constitutive matrix
\return void
\author bborn
\date 04/07
*/
void so3_mat_stress(CONTAINER *container,
                    ELEMENT *ele,
                    MATERIAL *mat,
                    INT ip,
                    SO3_GEODEFSTR *gds,
                    DOUBLE stress[NUMSTR_SOLID3],
                    DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_stress");
#endif

  /*====================================================================*/
  /* ===> central material routines
   *
   *      These materials are supposed to be connected to the
   *      existant (or new?) central material routines.
   *      Right now, only the simple St.Venant-Kirchhoff material
   *      is included to test the element.
   *
   * ===> Do it, have fun. Or don't.
   */

  /*====================================================================*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    /*------------------------------------------------------------------*/
    /* St. Venant-Kirchhoff material */
    case m_stvenant:
      so3_mat_stvenant_sel(container, ele, mat, ip, gds, stress, cmat);
      break;
    /*------------------------------------------------------------------*/
    /* Robinson's visco-plastic temperature-dependent material */
#ifdef D_TSI
    case m_vp_robinson:
      so3_mat_robinson_stress(container, ele, mat->m.vp_robinson, ip, gds,
                              stress, cmat);
      break;
#endif
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_mat_sel */

/*======================================================================*/
/*!
\brief Get density of material law

\param  mat       MATERIAL*   (i)   material data
\param  density   DOUBLE*     (o)   density value
\return void

\author bborn
\date 01/07
*/
void so3_mat_density(MATERIAL *mat,
                     DOUBLE *density)
{

#ifdef DEBUG
  dstrc_enter("so3_mat_density");
#endif

  /* switch material type */
  switch(mat->mattyp)
  {
    /* de St.Venant-Kirchhoff material */
    case m_stvenant:
      *density = mat->m.stvenant->density;
      break;
    /* kompressible neo-hooke */
    case m_neohooke:
      *density = mat->m.neohooke->density;
      break;
    /* porous linear elastic */
    case m_stvenpor:
      *density = mat->m.stvenpor->density;
      break;
    /* hyperelastic polyconvex material */
    case m_hyper_polyconvex:
      *density = mat ->m.hyper_polyconvex->density;
      break;
    /* Robinson's visco-plastic material */
    case m_vp_robinson:
      *density = mat ->m.vp_robinson->density;
      break;
    /* */
    default:
      dserror("Density of chosen material is not defined!");
      break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Test whether an update of material internal variables is required
\author bborn
\date 05/07
*/
void so3_mat_mivupdreq(ELEMENT* ele,
                       MATERIAL* mat,
                       INT* updreq)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_mivupdreq");
#endif

  /*--------------------------------------------------------------------*/
  /* check materials */
  switch (mat->mattyp)
  {
    /*------------------------------------------------------------------*/
    /* St. Venant-Kirchhoff material */
    case m_stvenant:
      /* do nothing */
      break;
    /*------------------------------------------------------------------*/
    /* Robinson's visco-plastic temperature-dependent material */
#ifdef D_TSI
    case m_vp_robinson:
      so3_mat_robinson_mivupdreq(ele, mat->m.vp_robinson, updreq);
      break;
#endif
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Iterative update material internal variables

\param  container CONTAINER*      (i)   container
\param  ele       ELEMENT*        (io)  pointer to current element
\param  mat       MATERIAL*       (i)   pointer to current material
\param  ip        INT             (i)   current Gauss point index
\param  gds       SO3_GEODEFSTR*  (i)   geom. & def. data at Gauss point
\param  epsii     DOUBLE[]        (i)   iterative strain increment
\return void

\author bborn
\date 05/07
*/
void so3_mat_mivupditer(const CONTAINER* container,
                        ELEMENT* ele,
                        const MATERIAL* mat,
                        const INT ip,
                        const SO3_GEODEFSTR* gds,
                        const DOUBLE epsii[NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_mivupditer");
#endif

  /*====================================================================*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    /*------------------------------------------------------------------*/
    /* St. Venant-Kirchhoff material */
    case m_stvenant:
      /* do nothing */
      break;
    /*------------------------------------------------------------------*/
    /* Robinson's visco-plastic temperature-dependent material */
#ifdef D_TSI
    case m_vp_robinson:
      so3_mat_robinson_mivupditer(container, ele, mat->m.vp_robinson,
                                  ip, epsii);
      break;
#endif
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Update material internal variables

\param *container CONTAINER      (i)   container
\param *ele       ELEMENT        (i)   pointer to current element
\param *mat       MATERIAL       (i)   pointer to current material
\param ip         INT            (i)   current Gauss point index
\param *gds       SO3_GEODEFSTR  (i)   geom. & def. data at Gauss point
\param stress[]   DOUBLE         (o)   linear(Biot)/2.Piola-Kirchhoff stress
\param cmat[][]   DOUBLE         (o)   constitutive matrix
\return void

\author bborn
\date 03/07
*/
void so3_mat_mivupdincr(const CONTAINER* container,
                        ELEMENT* ele,
                        const MATERIAL* mat)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_mivupdincr");
#endif

  /*====================================================================*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    /*------------------------------------------------------------------*/
    /* St. Venant-Kirchhoff material */
    case m_stvenant:
      /* do nothing */
      break;
    /*------------------------------------------------------------------*/
    /* Robinson's visco-plastic temperature-dependent material */
#ifdef D_TSI
    case m_vp_robinson:
      so3_mat_robinson_mivupdincr(container, ele, mat->m.vp_robinson);
      break;
#endif
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
#endif
