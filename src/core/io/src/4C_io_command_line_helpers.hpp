// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_COMMAND_LINE_HELPERS_HPP
#define FOUR_C_IO_COMMAND_LINE_HELPERS_HPP

#include "4C_config.hpp"

#include "4C_linalg_multi_vector.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 * \brief Structure to hold legacy CLI option names.
 */
struct LegacyCliOptions
{
  std::vector<std::string> single_dash_legacy_names;
  std::vector<std::string> nodash_legacy_names;
};

/**
 * \brief Adapt legacy command line arguments.
 * This function converts legacy single-dash options (e.g. "-ngroup=2") into
 * their new form ("--ngroup=2") and combines legacy dashless positional
 * options (e.g. two positional arguments "restart=1" and "restart=2") into
 * comma-separated lists ("--restart=1,2").
 * \param args Input arguments (no program-name expected).
 * \param legacy_options Structure containing names of legacy options:
 *        - single_dash_legacy_names: Options that used a single dash and
 *          should be converted to the double-dash form (e.g. {"ngroup", "glayout"}).
 *        - nodash_legacy_names: Legacy dashless options that should be
 *          collected/combined (e.g. {"restart", "restartfrom"}).
 * \return Sanitized vector of arguments.
 */
FOUR_C_API(FOUR_C_CORE)
std::vector<std::string> adapt_legacy_cli_arguments(
    const std::vector<std::string>& args, LegacyCliOptions& legacy_options);

FOUR_C_NAMESPACE_CLOSE

#endif