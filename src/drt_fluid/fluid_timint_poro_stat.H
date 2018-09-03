/*!----------------------------------------------------------------------
\file fluid_timint_poro_stat.H
\brief

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_PORO_STAT_H
#define FLD_TIMINT_PORO_STAT_H


#include "fluid_timint_poro.H"
#include "fluid_timint_stat.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{
  class TimIntPoroStat : public TimIntStationary, public TimIntPoro
  {
   public:
    /// Standard Constructor
    TimIntPoroStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntPoroStat();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief return scheme-specific time integration parameter
    */
    virtual double TimIntParam() const { return 0.0; };

    /*!
    \brief read restart data
    */
    virtual void ReadRestart(int step);


   protected:
   private:
  };  // class TimIntPoroStat

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_PORO_STAT_H
