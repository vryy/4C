/*-----------------------------------------------------------*/
/*! \file

\brief generalized alpha time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_TIMINT_PORO_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_PORO_GENALPHA_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_genalpha.hpp"
#include "4C_fluid_timint_poro.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntPoroGenAlpha : public TimIntGenAlpha, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroGenAlpha(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);

    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief read restart data

    */
    void read_restart(int step) override;

   protected:
    /*!
    \brief compute values at intermediate time steps for gen.-alpha

    */
    void gen_alpha_intermediate_values() override;

   private:
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
