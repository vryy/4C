/*-----------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_RED_GENALPHA_H
#define BACI_FLUID_TIMINT_RED_GENALPHA_H


#include "baci_fluid_timint_genalpha.H"
#include "baci_fluid_timint_red.H"
#include "baci_linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntRedModelsGenAlpha : public TimIntGenAlpha, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    ~TimIntRedModelsGenAlpha() override;

    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief read restart data

    */
    void ReadRestart(int step) override;


   protected:
   private:
  };  // class TimIntRedModelsGenAlpha

}  // namespace FLD


#endif  // FLUID_TIMINT_RED_GENALPHA_H
