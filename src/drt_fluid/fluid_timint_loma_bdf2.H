/*-----------------------------------------------------------*/
/*! \file

\brief BDF2 time integrator for Low Mach number problems

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_LOMA_BDF2_H
#define FLD_TIMINT_LOMA_BDF2_H


#include "fluid_timint_loma.H"
#include "fluid_timint_bdf2.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntLomaBDF2 : public TimIntBDF2, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntLomaBDF2();

    /*!
    \brief initialization

    */
    virtual void Init();


   protected:
   private:
  };  // class TimIntLomaBDF2

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_LOMA_BDF2_H
