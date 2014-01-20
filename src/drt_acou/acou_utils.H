/*----------------------------------------------------------------------*/
/*!
\file acou_utils.H

\brief utility functions for scalar transport problems

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_UTILS_H
#define ACOU_UTILS_H

#include "../drt_inpar/inpar_acou.H"

namespace ACOU
{

/*!
\brief Fills in Butcher tableau coefficients for DIRK schemes

The coeffcients are stored for four differenkt diagonally implicit Runge-Kutta schemes.
It is important to note, that the coeffcients a[i][i] are equal for all i. This allows
for once-only assembly of the system matrix and reuse of the preconditioners. If new DIRK
schemes are added, make sure, that a[i][i] are all the same!
*/
void FillDIRKValues(INPAR::ACOU::DynamicType scheme,
                    double (&a)[6][6],
                    double (&b)[6],
                    unsigned int &q);

std::string DIRKTypeToString(INPAR::ACOU::DynamicType scheme);

} // namespace ACOU

#endif /* ACOU_UTILS_H */
