/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef BACI_FLUID_TIMINT_RED_OST_H
#define BACI_FLUID_TIMINT_RED_OST_H


#include "baci_fluid_timint_red.H"
#include "baci_fluid_timint_ost.H"

#include "baci_linalg_utils_sparse_algebra_math.H"


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
