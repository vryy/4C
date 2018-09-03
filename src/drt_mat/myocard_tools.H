/*----------------------------------------------------------------------*/
/*!
\file myocard_tools.H
\brief myocard tools

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              cbert 08/13 |
 *----------------------------------------------------------------------*/
#ifndef MYOCARD_TOOLS_H
#define MYOCARD_TOOLS_H


/*----------------------------------------------------------------------*/
/// Myocard math tools general to all materials
///
/// \author cbert


/// \date 08/13

class Myocard_Tools
{
 public:
  /// construct empty material object
  Myocard_Tools();

  /// destructor
  virtual ~Myocard_Tools() { ; }

  /// compute Heaviside step function
  double GatingFunction(const double Gate1, const double Gate2, const double p, const double var,
      const double thresh) const;

  /// compute gating variable 'y' from dy/dt = (y_inf-y)/y_tau
  double GatingVarCalc(const double dt, double y_0, const double y_inf, const double y_tau) const;

};  // Myocard_Tools

#endif
