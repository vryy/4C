/*----------------------------------------------------------------------*/
/*! \file

\brief general initialization routine of BACI

\maintainer Martin Kronbichler

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef GLOBAL_INIT_CONTROL_H
#define GLOBAL_INIT_CONTROL_H


#include <string>

//! general initialization and routine
void ntaini_ccadiscret(int argc,       //!< number of command line arguments
    char** argv,                       //!< actual command line arguments
    std::string& inputfile_name,       //!< input file name
    std::string& outputfile_kenner,    //!< output file kenner
    std::string& restartfile_kenner);  //!< restart file kenner


#endif
