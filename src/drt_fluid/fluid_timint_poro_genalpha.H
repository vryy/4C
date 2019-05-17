/*-----------------------------------------------------------*/
/*!

\brief generalized alpha time integration scheme for porous fluid

\maintainer Christoph Ager

\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_TIMINT_PORO_GENALPHA_H_
#define FLUID_TIMINT_PORO_GENALPHA_H_


#include "fluid_timint_poro.H"
#include "fluid_timint_genalpha.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{
  class TimIntPoroGenAlpha : public TimIntGenAlpha, public TimIntPoro
  {
   public:
    /// Standard Constructor
    TimIntPoroGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntPoroGenAlpha();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief read restart data

    */
    virtual void ReadRestart(int step);

   protected:
    /*!
    \brief compute values at intermediate time steps for gen.-alpha

    */
    virtual void GenAlphaIntermediateValues();

   private:
  };  // class TimIntPoroGenAlpha

}  // namespace FLD



#endif /* FLUID_TIMINT_PORO_GENALPHA_H_ */
