// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_mpi_utils.hpp"

FOUR_C_NAMESPACE_OPEN

int Core::Communication::my_mpi_rank(const Epetra_Comm &comm) { return comm.MyPID(); }

int Core::Communication::num_mpi_ranks(const Epetra_Comm &comm) { return comm.NumProc(); }

void Core::Communication::barrier(const Epetra_Comm &comm) { comm.Barrier(); }

FOUR_C_NAMESPACE_CLOSE