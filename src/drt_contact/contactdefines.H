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


#define CONTACT_PROJTOL 0.05  		/* projection tol. for overlap decision */
#define CONTACT_MAXITER 10				/* max. no. iterations for local Newton   */
#define CONTACT_CONVTOL 1.0e-10		/* convergence tolerance for local Newton */
#define CONTACT_NGP 5							/* no. Gauss points for integration */
#define CONTACT_CRITDIST 0.5			/* critical distance for potential contact */


#endif  // #ifndef CONTACTDEFINES_H
#endif  // #ifdef CCADISCRET
