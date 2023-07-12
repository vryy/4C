/*-----------------------------------------------------------*/
/*! \file

\brief generalized alpha time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_TIMINT_PORO_GENALPHA_H
#define FLUID_TIMINT_PORO_GENALPHA_H


#include "fluid_timint_poro.H"
#include "fluid_timint_genalpha.H"

#include "linalg_utils_sparse_algebra_math.H"


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


#endif  // FLUID_TIMINT_PORO_GENALPHA_H
