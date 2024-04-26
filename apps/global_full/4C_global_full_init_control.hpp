/*----------------------------------------------------------------------*/
/*! \file

\brief general initialization routine of 4C


\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_GLOBAL_FULL_INIT_CONTROL_HPP
#define FOUR_C_GLOBAL_FULL_INIT_CONTROL_HPP


#include "4C_config.hpp"

#include <string>

//! general initialization and routine
void ntaini_ccadiscret(int argc,       //!< number of command line arguments
    char** argv,                       //!< actual command line arguments
    std::string& inputfile_name,       //!< input file name
    std::string& outputfile_kenner,    //!< output file kenner
    std::string& restartfile_kenner);  //!< restart file kenner

#endif
