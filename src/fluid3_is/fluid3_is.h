#ifndef FLUID3_IS_H
#define FLUID3_IS_H

#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"

/* caution! We use the fluid2 element structure here! This way we can
 * use all the fluid2 element functions, even is these access the
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

#endif
