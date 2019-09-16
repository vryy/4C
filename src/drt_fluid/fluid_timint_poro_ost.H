/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration scheme for porous fluid

\maintainer Christoph Ager

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_PORO_OST_H
#define FLD_TIMINT_PORO_OST_H


#include "fluid_timint_poro.H"
#include "fluid_timint_ost.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{
  class TimIntPoroOst : public TimIntOneStepTheta, public TimIntPoro
  {
   public:
    /// Standard Constructor
    TimIntPoroOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntPoroOst();

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
  };  // class TimIntPoroOst

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_PORO_OST_H
