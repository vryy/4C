/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration scheme for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_PORO_OST_H
#define FLD_TIMINT_PORO_OST_H


#include "fluid_timint_poro.H"
#include "fluid_timint_ost.H"

#include "linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntPoroOst : public TimIntOneStepTheta, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
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


#endif  // #ifndef FLD_TIMINT_PORO_OST_H
