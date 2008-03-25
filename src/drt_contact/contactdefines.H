/*!----------------------------------------------------------------------
\file contactdefines.H

\brief A set of preprocessor defines for contact methods

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONTACTDEFINES_H
#define CONTACTDEFINES_H

#define CONTACTWNORMAL           /* flag for weighted nodal normal calc. */
//#define CONTACTDIMOUTPUT       /* flag for output of global matrix dim. */
//#define CONTACTCHECKHUEEBER    /* flag for S.Hüeber's quad. conv. check */ 

#define CONTACTPROJTOL 0.01      /* projection tol. for overlap decision */
#define CONTACTPROJLIM 1.0e-12   /* no problems with exact node projections */
#define CONTACTMAXITER 10        /* max. no. iterations for local Newton   */
#define CONTACTCONVTOL 1.0e-12   /* convergence tolerance for local Newton */
#define CONTACTNGP 5             /* no. Gauss points for integration */
#define CONTACTCRITDIST 0.5      /* critical distance for potential contact */

#define CONTACTINITIALSET 0      /* flag for initialization of active set */
                                 /* 0=empty, 1=full, 2=user-defined      */

// static solution routine
void contact_stru_static_drt();

#endif  // #ifndef CONTACTDEFINES_H
#endif  // #ifdef CCADISCRET
