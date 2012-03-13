/*!----------------------------------------------------------------------
\file global_inp_control.H
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifndef GLOBAL_INP_CONTROL_H
#define GLOBAL_INP_CONTROL_H

#include <string>

//> general input routine for baci
void ntainp_ccadiscret(
    std::string& inputfile_name,      //!< input file name
    std::string& outputfile_kenner,   //!< output file kenner
    std::string& restartfile_kenner   //!< restart file kenner
  );


#endif  // #ifndef GLOBAL_INP_CONTROL2_H

