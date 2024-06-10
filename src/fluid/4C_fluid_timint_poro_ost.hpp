/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_PORO_OST_HPP
#define FOUR_C_FLUID_TIMINT_PORO_OST_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_ost.hpp"
#include "4C_fluid_timint_poro.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntPoroOst : public TimIntOneStepTheta, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroOst(const Teuchos::RCP<Core::FE::Discretization>& actdis,
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
   private:
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
