// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_epetra_exceptions.hpp"

#include <Epetra_Object.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void set_trace_back_mode(int tracebackmode)
{
  static const int global_mode = Epetra_Object::GetTracebackMode();
  Epetra_Object::SetTracebackMode(std::max(global_mode, tracebackmode));
}
