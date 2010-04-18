/*!
 \brief inf-sup fluid3 element definition.

 The elements is based on the normal fluid2 element. May functions are
 in common. The other ones have been adapted.

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
 */

#ifndef CCADISCRET
#ifdef D_FLUID3_IS

#ifndef FLUID3_IS_H
#define FLUID3_IS_H

#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"

/* caution! We use the fluid3 element structure here! This way we can
 * use all the fluid3 element functions, even is these access the
 * element data. But we cannot add any specific data... */
typedef FLUID3 FLUID3_IS;


void f3is_inp(ELEMENT* ele);
void fluid3_is(PARTITION   *actpart,
	       INTRA       *actintra,
	       ELEMENT     *ele,
	       ARRAY       *estif_global,
	       ARRAY       *emass_global,
	       ARRAY       *eforce_global,
	       ARRAY       *edforce_global,
	       CALC_ACTION *action,
	       INT         *hasdirich,
	       INT         *hasext,
	       CONTAINER   *container
  );

void f3is_calele(
  ELEMENT        *ele,
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext,
  INT             is_relax,
  INT             init
  );

void f3is_caleleres(
  ELEMENT        *ele,
  ARRAY          *eforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  );

void f3is_caleleres_relax(ELEMENT        *ele,
			  ARRAY          *estif_global,
			  ARRAY          *eforce_global,
			  ARRAY_POSITION *ipos,
			  INT            *hasdirich,
			  INT            *hasext);

void f3is_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **ehist,
  DOUBLE         **evelng,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos,
  INT             *hasext
  );

void f3is_calseta(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **ehist,
  DOUBLE         **evelng,
  DOUBLE         **ealecovng,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  ARRAY_POSITION  *ipos,
  INT             *hasext,
  INT              is_relax
  );

void f3is_int_usfem(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE         **estif,
  DOUBLE          *force,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE          *pfunct,
  DOUBLE         **pderiv,
  DOUBLE         **pderiv2,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **pderxy,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  DOUBLE         **vderxy,
  DOUBLE         **vderxy2,
  DOUBLE           visc,
  DOUBLE         **wa1,
  DOUBLE         **wa2,
  INT              is_relax
  );

void f3is_calmat( DOUBLE **estif,
		  DOUBLE  *eforce,
		  DOUBLE  *velint,
		  DOUBLE   histvec[3],
		  DOUBLE   gridvint[3],
		  DOUBLE   press,
		  DOUBLE **vderxy,
		  DOUBLE **vderxy2,
		  DOUBLE   gradp[3],
		  DOUBLE  *funct,
		  DOUBLE **derxy,
		  DOUBLE **derxy2,
		  DOUBLE  *pfunct,
		  DOUBLE **pderxy,
		  DOUBLE  *edeadng,
		  DOUBLE   fac,
		  DOUBLE   visc,
		  INT      iel,
		  INT     *hasext,
		  INT      isale,
		  INT      is_relax
  );

void f3is_int_res(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE          *force,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE          *pfunct,
  DOUBLE         **pderiv,
  DOUBLE         **pderiv2,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **pderxy,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  DOUBLE         **vderxy,
  DOUBLE         **vderxy2,
  DOUBLE           visc,
  DOUBLE         **wa1,
  DOUBLE         **wa2
  );

void f3is_calresvec(
  DOUBLE  *eforce,
  DOUBLE  *velint,
  DOUBLE   histvec[3],
  DOUBLE **vderxy,
  DOUBLE **vderxy2,
  DOUBLE  *funct,
  DOUBLE  *pfunct,
  DOUBLE **derxy,
  DOUBLE **derxy2,
  DOUBLE  *edeadng,
  DOUBLE   gridvint[3],
  DOUBLE   press,
  DOUBLE   gradp[3],
  DOUBLE   fac,
  DOUBLE   visc,
  INT      iel,
  INT     *hasext,
  INT      is_ale
  );

#endif
#endif
#endif
