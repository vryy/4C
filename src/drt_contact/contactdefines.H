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

//#define CONTACTDIMOUTPUT       /* flag for output of global matrix dim. */

#define CONTACTPROJTOL 0.05      /* projection tol. for overlap decision */
#define CONTACTPROJLIM 1.0e-6    /* no problems with exact node projections */
#define CONTACTMAXITER 10        /* max. no. iterations for local Newton   */
#define CONTACTCONVTOL 1.0e-10   /* convergence tolerance for local Newton */
#define CONTACTNGP 5             /* no. Gauss points for integration */
#define CONTACTCRITDIST 0.7      /* critical distance for potential contact */

// FIXME: Might be better to have critdist = f(element length)


#endif  // #ifndef CONTACTDEFINES_H
#endif  // #ifdef CCADISCRET
