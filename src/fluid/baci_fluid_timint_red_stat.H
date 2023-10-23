/*-----------------------------------------------------------*/
/*! \file

\brief Stationary driver for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_RED_STAT_H
#define BACI_FLUID_TIMINT_RED_STAT_H


#include "baci_fluid_timint_red.H"
#include "baci_fluid_timint_stat.H"
#include "baci_linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntRedModelsStat : public TimIntStationary, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


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
  };  // class TimIntRedModelsStat

}  // namespace FLD


#endif  // FLUID_TIMINT_RED_STAT_H
