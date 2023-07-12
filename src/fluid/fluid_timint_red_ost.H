/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_TIMINT_RED_OST_H
#define FLUID_TIMINT_RED_OST_H


#include "fluid_timint_red.H"
#include "fluid_timint_ost.H"

#include "linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntRedModelsOst : public TimIntOneStepTheta, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntRedModelsOst();

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
  };  // class TimIntRedModelsOst

}  // namespace FLD


#endif  // FLUID_TIMINT_RED_OST_H
