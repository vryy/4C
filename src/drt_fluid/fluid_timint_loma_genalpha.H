/*-----------------------------------------------------------*/
/*! \file

\brief Gen-alpha time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_LOMA_GENALPHA_H
#define FLD_TIMINT_LOMA_GENALPHA_H


#include "fluid_timint_loma.H"
#include "fluid_timint_genalpha.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntLomaGenAlpha : public TimIntGenAlpha, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntLomaGenAlpha();

    /*!
    \brief initialization

    */
    virtual void Init();


   protected:
   private:
  };  // class TimIntLomaGenAlpha

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_LOMA_GENALPHA_H
