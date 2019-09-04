/*-----------------------------------------------------------*/
/*! \file

\brief Fluid time integrator for FS3I-AC problems

\maintainer Christoph Ager

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_AC_H
#define FLD_TIMINT_AC_H


#include "fluidimplicitintegration.H"

namespace FLD
{
  class TimIntAC : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntAC(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntAC();

    /// read restart from step
    void ReadRestart(int step);

    /// write output
    virtual void Output();

   protected:
   private:
  };  // class TimIntAC

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_AC_H
