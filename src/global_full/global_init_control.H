/*----------------------------------------------------------------------*/
/*!
\file global_init_control.H

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef GLOBAL_INIT_CONTROL_H
#define GLOBAL_INIT_CONTROL_H


#include <string>

//! general initialization and routine
void ntaini_ccadiscret(
    int argc,                         //!< number of command line arguments
    char** argv,                      //!< actual command line arguments
    std::string& inputfile_name,      //!< input file name
    std::string& outputfile_kenner,   //!< output file kenner
    std::string& restartfile_kenner); //!< restart file kenner


#endif
