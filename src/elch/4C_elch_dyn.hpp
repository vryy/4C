// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ELCH_DYN_HPP
#define FOUR_C_ELCH_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of electrochemistry problems */
void elch_dyn(int restart /* do we have to perform a restart?  */
);

/*! prints the 4C electrochemistry-module logo on the screen */
void printlogo();


FOUR_C_NAMESPACE_CLOSE

#endif
