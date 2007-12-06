/*!----------------------------------------------------------------------
 * \file
 * \brief experimental anisotropic formulation for cells
 *
 * <pre>
 * Maintainer: Robert Metzke
 * metzke@lnm.mw.tum.de
 * http://www.lnm.mw.tum.de/Members/metzke/
 * 089 289 15244
 * </pre>
 *----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
/*!
 * \addtogroup BRICK1
 *//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
 * \brief establish local material
 *
 * <pre>                                                              rm 11/07
 * </pre>
 *
 * \warning There is nothing special to this routine
 * \return void
 * \sa calling: ---; called by: c1_call_mat()
 *
 *----------------------------------------------------------------------*/
void c1_mat_quadratic_aniso (	
                                ELEMENT *ele,
				INT ip,
				MATERIAL *mat,
				DOUBLE *disd,
				DOUBLE *stress,
				DOUBLE **d) {
    
    printf("Hallo, ich bin's");
    exit(1);
    
}
#endif
/*! @} (documentation module close)*/
#endif /*D_BRICK1 */
