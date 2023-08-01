/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_PORO_OST_H
#define BACI_FLUID_TIMINT_PORO_OST_H


#include "baci_fluid_timint_ost.H"
#include "baci_fluid_timint_poro.H"
#include "baci_linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntPoroOst : public TimIntOneStepTheta, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroOst(const Teuchos::RCP<DRT::Discretization>& actdis,
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
  };

}  // namespace FLD


#endif  // FLUID_TIMINT_PORO_OST_H
