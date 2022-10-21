/*----------------------------------------------------------------------------*/
/*! \file
\brief shell8

\level 1


*/
/*---------------------------------------------------------------------------*/
#include "../../headers/standardtypes.h"
#include "shell8.h"

#include "../../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 | get density out of material law                        m.gee 2/02    |
 *----------------------------------------------------------------------*/
void s8_getdensity(MATERIAL *mat, DOUBLE *density)
{
  /*----------------------------------------------------------------------*/
  /*------------------------------------------------ switch material type */
  switch (mat->mattyp)
  {
    case m_stvenant: /*-------------------------- ST.VENANT-KIRCHHOFF-MATERIAL */
      *density = mat->m.stvenant->density;
      break;
    case m_neohooke: /*------------------------------ kompressible neo-hooke */
      *density = mat->m.neohooke->density;
      break;
    case m_compogden: /*-------------------------------- kompressible ogden */
      *density = mat->m.compogden->density;
      break;
    case m_viscohyper: /*---------------------------viscos compressible ogden */
      *density = mat->m.viscohyper->density;
      break;
    default:
      dserror("Ilegal typ of material for this element");
      break;
  }
  /*----------------------------------------------------------------------*/
  return;
} /* end of s8_getdensity */
