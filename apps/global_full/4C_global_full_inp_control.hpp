/*----------------------------------------------------------------------*/
/*! \file

\brief Global control routine of baci

\level 0


*----------------------------------------------------------------------*/

#ifndef FOUR_C_GLOBAL_FULL_INP_CONTROL_HPP
#define FOUR_C_GLOBAL_FULL_INP_CONTROL_HPP

#include "4C_config.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

#include <string>

//> general input routine for baci
void ntainp_ccadiscret(std::string& inputfile_name,  //!< input file name
    std::string& outputfile_kenner,                  //!< output file kenner
    std::string& restartfile_kenner                  //!< restart file kenner
);

//> setup of parallel output as early as possible
void SetupParallelOutput(
    std::string& outputfile_kenner, Teuchos::RCP<Epetra_Comm> lcomm, int group);

#endif
