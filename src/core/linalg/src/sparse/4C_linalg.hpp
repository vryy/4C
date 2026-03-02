// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_HPP
#define FOUR_C_LINALG_HPP

#include "4C_config.hpp"

#include <Epetra_CombineMode.h>

#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*! \enum Core::LinAlg::DataAccess
   *  \brief Handling of data access (Copy or View)
   *
   *  If set to Core::LinAlg::DataAccess::Copy, user data will be copied at construction.
   *  If set to Core::LinAlg::DataAccess::Share, user data will be shared.
   *
   *  \note A separate Core::LinAlg::DataAccess is necessary in order to resolve
   *  possible ambiguity conflicts with the Epetra_DataAccess.
   *
   *  Use Core::LinAlg::DataAccess for construction of any Core::LINALG matrix object.
   *  Use plain 'Copy' or 'View' for construction of any Epetra matrix object.
   *
   */
  enum class DataAccess
  {
    Copy,  ///< deep copy
    Share  ///< Shared ownership to original data
  };

  /**
   * \brief Specifies the strategy used when combining distributed data.
   *
   * The CombineMode enum defines how values are combined when assembling or updating
   * distributed linear algebra objects (such as vectors or matrices) in parallel computations.
   * It determines how overlapping entries from different processes are merged.
   */
  enum class CombineMode
  {
    zero,
    insert,
    add,
    max,
    abs_max
  };

  /**
   * @brief Describes how an object is distributed across MPI ranks.
   *
   * This enumeration specifies whether an object is stored redundantly on every process or
   * partitioned across processes in a distributed-memory parallel environment.
   */
  enum class LocalGlobal
  {
    locally_replicated,
    globally_distributed
  };

  /**
   * \brief Mapping between internal CombineMode values and corresponding Epetra_CombineMode
   * constants.
   *
   * This lookup table allows translation between the internal linear algebra abstraction
   * and the Epetra library's specific combine mode definitions.
   */
  const std::unordered_map<CombineMode, Epetra_CombineMode> linalg_to_epetra_combine_mode = {
      {CombineMode::zero, Epetra_CombineMode::Zero},
      {CombineMode::insert, Epetra_CombineMode::Insert},
      {CombineMode::add, Epetra_CombineMode::Add},
      {CombineMode::max, Epetra_CombineMode::Epetra_Max},
      {CombineMode::abs_max, Epetra_CombineMode::AbsMax}};

  /**
   * \brief Converts a CombineMode enum value to its Epetra_CombineMode equivalent.
   *
   * @param combine_mode The CombineMode value to be converted.
   * @return The corresponding Epetra_CombineMode constant.
   */
  inline Epetra_CombineMode to_epetra_combine_mode(CombineMode combine_mode)
  {
    return linalg_to_epetra_combine_mode.at(combine_mode);
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
