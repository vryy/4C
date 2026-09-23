// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_FENV_HPP
#define FOUR_C_UTILS_FENV_HPP

#include "4C_config.hpp"

#include <cfenv>

FOUR_C_NAMESPACE_OPEN

#if defined(__APPLE__)

#if defined(__aarch64__)

static inline int fegetexcept(void)
{
  // Apple Silicon (ARM64) does not support hardware FP exception traps
  return 0;
}

static inline int fedisableexcept(unsigned int excepts)
{
  // Hardware traps for FP exceptions are unsupported on ARM64 macOS
  return 0;
}
static inline int feenableexcept(unsigned int excepts)
{
  return -1;  // Indicate unsupported
}

#elif defined(__x86_64__)

static inline int fegetexcept(void)
{
  // Intel Mac implementation reading MXCSR mask bits
  std::fenv_t fenv;
  std::fegetenv(&fenv);
  // In MXCSR, bits 7-12 are mask bits (0 = enabled/unmasked, 1 = disabled/masked)
  return (~(fenv.__mxcsr >> 7)) & FE_ALL_EXCEPT;
}

static inline int fedisableexcept(unsigned int excepts)
{
  fenv_t fenv;
  fegetenv(&fenv);
  // On x86, setting the exception mask bits disables traps
  unsigned int new_excepts = (excepts & FE_ALL_EXCEPT) << 7;
  fenv.__mxcsr |= new_excepts;
  fenv.__control |= (excepts & FE_ALL_EXCEPT);
  fesetenv(&fenv);
  return 0;
}

static inline int feenableexcept(unsigned int excepts)
{
  fenv_t fenv;
  fegetenv(&fenv);
  // Unmasking enables traps
  unsigned int new_excepts = (excepts & FE_ALL_EXCEPT) << 7;
  fenv.__mxcsr &= ~new_excepts;
  fenv.__control &= ~(excepts & FE_ALL_EXCEPT);
  fesetenv(&fenv);
  return 0;
}

#endif  // __aarch64__

#endif  // __APPLE__

FOUR_C_NAMESPACE_CLOSE

#endif
