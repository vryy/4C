
/*----------------------------------------------------------------------*/
/*!
\file fsi_dyn.H

\brief Entry routines for FSI problems and some other problem types as well

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_DYN_H
#define FSI_DYN_H


void fluid_ale_drt();
void fluid_xfem_drt();
void fluid_fluid_ale_drt();
void fluid_fluid_fsi_drt();
void fluid_freesurf_drt();
void fsi_ale_drt();
void xfsi_drt();


#endif

