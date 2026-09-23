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

#elif defined(_MSC_VER)

#include <cfloat>

// Helper: Translate standard FE_* flags to MSVC _EM_* mask bits
static inline unsigned int fenv_to_msvc_mask(int excepts)
{
  unsigned int mask = 0;
  if (excepts & FE_INVALID) mask |= _EM_INVALID;
  if (excepts & FE_DIVBYZERO) mask |= _EM_ZERODIVIDE;
  if (excepts & FE_OVERFLOW) mask |= _EM_OVERFLOW;
  if (excepts & FE_UNDERFLOW) mask |= _EM_UNDERFLOW;
  if (excepts & FE_INEXACT) mask |= _EM_INEXACT;
  return mask;
}

// Helper: Translate MSVC _EM_* mask bits back to FE_* flags
// Note: In MSVC control word, a bit is 1 if the exception is MASKED (disabled),
// and 0 if it is UNMASKED (enabled/trapping).
static inline int msvc_mask_to_fenv(unsigned int cw)
{
  int excepts = 0;
  if (!(cw & _EM_INVALID)) excepts |= FE_INVALID;
  if (!(cw & _EM_ZERODIVIDE)) excepts |= FE_DIVBYZERO;
  if (!(cw & _EM_OVERFLOW)) excepts |= FE_OVERFLOW;
  if (!(cw & _EM_UNDERFLOW)) excepts |= FE_UNDERFLOW;
  if (!(cw & _EM_INEXACT)) excepts |= FE_INEXACT;
  return excepts;
}

// Returns the previously enabled (unmasked) exceptions, or -1 on error.
inline int fegetexcept(void)
{
  unsigned int cw = 0;
  if (_controlfp_s(&cw, 0, 0) != 0)
  {
    return -1;
  }
  return msvc_mask_to_fenv(cw);
}

// Unmasks the given exceptions so hardware traps are generated.
// Returns the previously enabled (unmasked) exceptions, or -1 on error.
inline int feenableexcept(int excepts)
{
  unsigned int old_cw = 0;
  if (_controlfp_s(&old_cw, 0, 0) != 0)
  {
    return -1;
  }

  unsigned int mask = fenv_to_msvc_mask(excepts);
  unsigned int new_cw = 0;
  // Clearing the bits (passing 0 for value) UNMASKS them (enables trapping)
  if (_controlfp_s(&new_cw, 0, mask) != 0)
  {
    return -1;
  }

  return msvc_mask_to_fenv(old_cw);
}

// Masks the given exceptions so they are ignored / silent.
// Returns the previously enabled (unmasked) exceptions, or -1 on error.
inline int fedisableexcept(int excepts)
{
  unsigned int old_cw = 0;
  if (_controlfp_s(&old_cw, 0, 0) != 0)
  {
    return -1;
  }

  unsigned int mask = fenv_to_msvc_mask(excepts);
  unsigned int new_cw = 0;
  // Setting the bits (passing mask for value) MASKS them (disables trapping)
  if (_controlfp_s(&new_cw, mask, mask) != 0)
  {
    return -1;
  }

  return msvc_mask_to_fenv(old_cw);
}

#endif  // _MSC_VER

FOUR_C_NAMESPACE_CLOSE

#endif
