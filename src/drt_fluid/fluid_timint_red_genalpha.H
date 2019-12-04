/*-----------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time integration for reduced models

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_RED_GENALPHA_H
#define FLD_TIMINT_RED_GENALPHA_H


#include "fluid_timint_red.H"
#include "fluid_timint_genalpha.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntRedModelsGenAlpha : public TimIntGenAlpha, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntRedModelsGenAlpha();

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
  };  // class TimIntRedModelsGenAlpha

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_RED_GENALPHA_H
