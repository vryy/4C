/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_validate.H

\brief preprocessor for exodusII format 

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bauer
            089 - 289-15252
</pre>

Validate a given BACI input file after all preprocessing steps

*/

#ifdef D_EXODUS
#ifdef CCADISCRET

#ifndef PRE_EXODUS_VALIDATE_H_
#define PRE_EXODUS_VALIDATE_H_

#include <string>

using namespace std;

namespace EXODUS{
void ValidateInputFile(string datfile);
}

#endif // PRE_EXODUS_VALIDATE_H_

#endif
#endif
