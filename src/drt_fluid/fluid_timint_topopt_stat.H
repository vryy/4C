/*-----------------------------------------------------------*/
/*! \file

\brief Stationary problem setup for topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_TOPOPT_STAT_H
#define FLD_TIMINT_TOPOPT_STAT_H


#include "fluid_timint_topopt.H"
#include "fluid_timint_stat.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntTopOptStat : public TimIntStationary, public TimIntTopOpt
  {
   public:
    /// Standard Constructor
    TimIntTopOptStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTopOptStat();

    /*!
    \brief initialization

    */
    virtual void Init();


   protected:
   private:
  };  // class TimIntTopOptStat

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_TOPOPT_STAT_H
