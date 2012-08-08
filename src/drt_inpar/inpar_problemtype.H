/*----------------------------------------------------------------------*/
/*!
\file inpar_problemtype.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_PROBLEMTYPE_H
#define INPAR_PROBLEMTYPE_H

#include <string>
#include <map>
#include "../drt_lib/drt_globalproblem_enums.H"


namespace DRT
{

/// create map of problem name and problem type enum
std::map<std::string,PROBLEM_TYP> StringToProblemTypeMap();

/// return problem type enum for a given problem name
PROBLEM_TYP StringToProblemType(std::string name);


} // namespace DRT

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_PROBLEMTYPE_H
