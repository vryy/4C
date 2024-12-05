// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROMULTIPHASE_PARTITIONED_HPP
#define FOUR_C_POROMULTIPHASE_PARTITIONED_HPP

#include "4C_config.hpp"

#include "4C_poromultiphase_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitioned : public PoroMultiPhaseBase
  {
   public:
    PoroMultiPhasePartitioned(MPI_Comm comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder


  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASE


FOUR_C_NAMESPACE_CLOSE

#endif
