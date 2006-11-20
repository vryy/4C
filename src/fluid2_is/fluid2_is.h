/*!
 \brief inf-sup fluid2 element definition.

 The elements is based on the normal fluid2 element. May functions are
 in common. The other ones have been adapted.

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
 */

#ifndef FLUID2_IS_H
#define FLUID2_IS_H

#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"

/* caution! We use the fluid2 element structure here! This way we can
 * use all the fluid2 element functions, even is these access the
 * element data. But we cannot add any specific data... */
typedef FLUID2 FLUID2_IS;


void f2is_inp(ELEMENT* ele);
void fluid2_is(PARTITION   *actpart,
	       INTRA       *actintra,
	       ELEMENT     *ele,
	       ELEMENT     *eleke,
	       ARRAY       *estif_global,
	       ARRAY       *emass_global,
	       ARRAY       *eforce_global,
	       ARRAY       *edforce_global,
	       CALC_ACTION *action,
	       INT         *hasdirich,
	       INT         *hasext,
	       CONTAINER   *container
  );

void f2is_calele(
	        ELEMENT        *ele,
                ELEMENT        *eleke,
                ARRAY          *estif_global,
                ARRAY          *emass_global,
	        ARRAY          *eforce_global,
		ARRAY          *edforce_global,
		INT            *hasdirich,
                INT            *hasext,
                INT             imyrank,
                ARRAY_POSITION *ipos,
		INT             is_relax,
		INT             init
  );

void f2is_calset(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE          *epren,
  DOUBLE          *edeadn,
  DOUBLE          *edeadng,
  ARRAY_POSITION *ipos,
  INT             *hasext
  );

void f2is_calseta(
  ELEMENT         *ele,
  DOUBLE         **xyze,
  DOUBLE         **eveln,
  DOUBLE         **evelng,
  DOUBLE         **evhist,
  DOUBLE         **ealecovn,
  DOUBLE         **ealecovng,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadn,
  DOUBLE          *edeadng,
  DOUBLE          *ekappan,
  DOUBLE          *ekappang,
  DOUBLE          *ephin,
  DOUBLE          *ephing,
  DOUBLE         **evnng,
  DOUBLE         **evnn,
  ARRAY_POSITION *ipos,
  INT             *hasext,
  INT              is_relax
  );


void f2is_int_usfem(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE         **estif,
  DOUBLE          *eforce,
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
  DOUBLE         **eveln,
  DOUBLE         **evhist,
  DOUBLE         **egridv,
  DOUBLE          *epren,
  DOUBLE          *edeadng,
  DOUBLE         **vderxy,
  DOUBLE         **vderxy2,
  DOUBLE           visc,
  DOUBLE         **wa1,
  DOUBLE         **wa2,
  DOUBLE           estress[3][MAXNOD_F2],
  INT              is_relax
  );


void f2is_calmat( DOUBLE **estif,
		DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE   histvec[2],
		DOUBLE   gridvint[2],
		DOUBLE   press,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
                DOUBLE   gradp[2],
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		  DOUBLE          *pfunct,
                DOUBLE  *edeadng,
		DOUBLE   fac,
		DOUBLE   visc,
		INT      iel,
                INT     *hasext,
                INT      isale,
                INT      is_relax
  );

void f2is_caleleres(
  ELEMENT        *ele,
  ARRAY          *eforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  );

void f2is_int_res(
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

void f2is_calresvec(
  DOUBLE  *eforce,
  DOUBLE  *velint,
  DOUBLE   histvec[2],
  DOUBLE **vderxy,
  DOUBLE **vderxy2,
  DOUBLE  *funct,
  DOUBLE  *pfunct,
  DOUBLE **derxy,
  DOUBLE **derxy2,
  DOUBLE  *edeadng,
  DOUBLE   aleconv[2],
  DOUBLE   press,
  DOUBLE   gradp[2],
  DOUBLE   fac,
  DOUBLE   visc,
  INT      iel,
  INT     *hasext,
  INT      is_ale
  );


#endif
