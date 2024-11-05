// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GLOBAL_FULL_INP_CONTROL_HPP
#define FOUR_C_GLOBAL_FULL_INP_CONTROL_HPP

#include "4C_config.hpp"

#include <Epetra_Comm.h>

#include <memory>
#include <string>

//> general input routine for 4C
void ntainp_ccadiscret(std::string& inputfile_name,  //!< input file name
    std::string& outputfile_kenner,                  //!< output file kenner
    std::string& restartfile_kenner                  //!< restart file kenner
);

//> setup of parallel output as early as possible
void setup_parallel_output(
    std::string& outputfile_kenner, std::shared_ptr<Epetra_Comm> lcomm, int group);

#endif
