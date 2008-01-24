/*----------------------------------------------------------------------*/
/*!
\file pre_exodus.H

\brief preprocessor for exodusII format 

<pre>
Maintainer: Moritz & Georg
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Pre_exodus contains classes to open and preprocess exodusII files into the
drt of Baci. It uses the "valid-parameters"-list defined in Baci for preparing
a up-to-date Baci header and another file specifying element and boundary
specifications. As result either a preliminary input file set is suggestioned,
or the well-known .dat file is created.

*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS

#ifndef PRE_EXODUS_H_
#define PRE_EXODUS_H_

#ifdef CCADISCRET

//#include <stdio.h>
//#include <stdlib.h>
#include <string>
#include <iostream>
//#include <exodusII.h>
extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}




#endif /* of CCADISCRET */
#endif /*PRE_EXODUS_H_*/
#endif 
