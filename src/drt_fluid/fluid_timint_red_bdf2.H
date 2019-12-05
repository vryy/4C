/*-----------------------------------------------------------*/
/*! \file

\brief BDF-2 time integration for reduced models

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_RED_BDF2_H
#define FLD_TIMINT_RED_BDF2_H


#include "fluid_timint_red.H"
#include "fluid_timint_bdf2.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntRedModelsBDF2 : public TimIntBDF2, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntRedModelsBDF2();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief read restart data

    */
    virtual void ReadRestart(int step);


   protected:
   private:
  };  // class TimIntRedModelsBDF2

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_RED_BDF2_H
