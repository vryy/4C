/*!-----------------------------------------------------------------------------------------------*
\file two_phase_defines.H

  \brief preprocessor define flags for two-phase flow problems

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

/*
 * two-phase flow flags:
 * to be defined only for special two-phase flow problems
 */


//-----------------------
// bubbly channel flow
//-----------------------
//#define BCF_SURFTENS_PRES_KINK //> use kink enrichment for pressure and apply surface tension (PROBLEMTYPE: Two_Phase_Flow_Surf)


//--------------------------
// multiscale two-phase flow
//--------------------------
//gmsh output of subgrid velocity
//#define SUGRVEL_OUTPUT //> don't forget to transfer the time step
