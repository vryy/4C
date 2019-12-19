/*-----------------------------------------------------------*/
/*! \file

\brief BDF-2 time integration for topology optimization

\maintainer Martin Kronbichler

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_TOPOPT_BDF2_H
#define FLD_TIMINT_TOPOPT_BDF2_H


#include "fluid_timint_topopt.H"
#include "fluid_timint_bdf2.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntTopOptBDF2 : public TimIntBDF2, public TimIntTopOpt
  {
   public:
    /// Standard Constructor
    TimIntTopOptBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTopOptBDF2();

    /*!
    \brief initialization

    */
    virtual void Init();


   protected:
   private:
  };  // class TimIntTopOptBDF2

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_TOPOPT_BDF2_H
