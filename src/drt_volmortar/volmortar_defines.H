/*!----------------------------------------------------------------------
\file volmortar_defines.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef VOLMORTAR_DEFINES_H_
#define VOLMORTAR_DEFINES_H_

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR INTEGRATION
#define VOLMORTARINTTOL   0.0       /* tolerance for assembling gp-values*/

// GEOMETRIC TOLERANCES
#define VOLMORTARELETOL   1.0e-12
#define VOLMORTARCUTTOL   1.0e-12
#define VOLMORTARCUT2TOL -1.0e-12

#endif /* VOLMORTAR_DEFINES_H_ */
