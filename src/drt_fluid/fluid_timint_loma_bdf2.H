/*!----------------------------------------------------------------------
\file fluid_timint_loma_bdf2.H
\brief

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_LOMA_BDF2_H
#define FLD_TIMINT_LOMA_BDF2_H


#include "fluid_timint_loma.H"
#include "fluid_timint_bdf2.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntLomaBDF2 : public TimIntBDF2, public TimIntLoma
{
public:
  /// Standard Constructor
  TimIntLomaBDF2(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntLomaBDF2();

  /*!
  \brief initialization

  */
  virtual void Init();


protected:

private:


}; // class TimIntLomaBDF2

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_LOMA_BDF2_H
