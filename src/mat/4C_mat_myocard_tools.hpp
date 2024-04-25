/*----------------------------------------------------------------------*/
/*! \file
\brief myocard tools

\level 3

*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_MYOCARD_TOOLS_HPP
#define FOUR_C_MAT_MYOCARD_TOOLS_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/// Myocard math tools general to all materials
///
/// \author cbert


/// \date 08/13

class MyocardTools
{
 public:
  /// construct empty material object
  MyocardTools();

  /// destructor
  virtual ~MyocardTools() = default;
  /// compute Heaviside step function
  double GatingFunction(const double Gate1, const double Gate2, const double p, const double var,
      const double thresh) const;

  /// compute gating variable 'y' from dy/dt = (y_inf-y)/y_tau
  double GatingVarCalc(const double dt, double y_0, const double y_inf, const double y_tau) const;

};  // Myocard_Tools

FOUR_C_NAMESPACE_CLOSE

#endif
