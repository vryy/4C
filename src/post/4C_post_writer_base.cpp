// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_post_writer_base.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_post_common.hpp"

FOUR_C_NAMESPACE_OPEN

PostWriterBase::PostWriterBase(PostField* field, const std::string& filename)
    : field_(field),
      filename_(filename),
      myrank_(Core::Communication::my_mpi_rank(*field->problem()->get_comm())),
      numproc_(Core::Communication::num_mpi_ranks(*field->problem()->get_comm()))
{
}

FOUR_C_NAMESPACE_CLOSE
