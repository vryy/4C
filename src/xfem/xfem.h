/*!----------------------------------------------------------------------
\file
\brief xfem.h

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM

#ifndef XFEM_H
#define XFEM_H

/*!
\addtogroup XFEM
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief XFEM execution flags for polygonization

<pre>                                                            irhan 05/04
This enumerator contains execution flags used in level set formulation
</pre>

*----------------------------------------------------------------------*/
typedef enum _XFEMPOLYFLAG
{
      xfem_poly_initialize,
      xfem_poly_material,
      xfem_poly_construct,
      xfem_poly_computeGP,
      xfem_poly_write,
      xfem_poly_open,
      xfem_poly_close
} XFEMPOLYFLAG;



/*!----------------------------------------------------------------------
\brief XFEM data structure

<pre>                                                            irhan 08/04
This structure contains options related to xfem formulation
</pre>

*----------------------------------------------------------------------*/
typedef struct _XFEM_DATA
{
  INT     xfem_on_off;
  INT     xfem_optimize;
  INT     xfem_area_check;
} XFEM_DATA;
/*! @} (documentation module close)*/

#endif
#endif
