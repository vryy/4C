/*!----------------------------------------------------------------------
\file fluid_timint_topopt_genalpha.H
\brief

\level 3

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_TOPOPT_GENALPHA_H
#define FLD_TIMINT_TOPOPT_GENALPHA_H


#include "fluid_timint_topopt.H"
#include "fluid_timint_genalpha.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntTopOptGenAlpha : public TimIntGenAlpha, public TimIntTopOpt
{
public:
  /// Standard Constructor
  TimIntTopOptGenAlpha(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntTopOptGenAlpha();

  /*!
  \brief initialization

  */
  virtual void Init();


protected:

private:


}; // class TimIntTopOptGenAlpha

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_TOPOPT_GENALPHA_H
