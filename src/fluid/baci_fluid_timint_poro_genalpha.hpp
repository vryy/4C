/*-----------------------------------------------------------*/
/*! \file

\brief generalized alpha time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_TIMINT_PORO_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_PORO_GENALPHA_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_genalpha.hpp"
#include "baci_fluid_timint_poro.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntPoroGenAlpha : public TimIntGenAlpha, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
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
    /*!
    \brief compute values at intermediate time steps for gen.-alpha

    */
    void GenAlphaIntermediateValues() override;

   private:
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
