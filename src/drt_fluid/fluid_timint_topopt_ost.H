/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_TOPOPT_OST_H
#define FLD_TIMINT_TOPOPT_OST_H


#include "fluid_timint_topopt.H"
#include "fluid_timint_ost.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntTopOptOst : public TimIntOneStepTheta, public TimIntTopOpt
  {
   public:
    /// Standard Constructor
    TimIntTopOptOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTopOptOst();

    /*!
    \brief initialization

    */
    virtual void Init();


   protected:
   private:
  };  // class TimIntTopOptOst

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_TOPOPT_OST_H
