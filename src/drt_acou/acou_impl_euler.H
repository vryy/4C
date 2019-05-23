/*!----------------------------------------------------------------------
\file acou_impl_euler.H
\brief

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*----------------------------------------------------------------------*/

#ifndef ACOU_IMPL_EULER_H
#define ACOU_IMPL_EULER_H

#include "acou_impl.H"

namespace ACOU
{
  class TimIntImplEuler : public AcouImplicitTimeInt
  {
   public:
    /// constructor
    TimIntImplEuler(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// return the name of the time integration scheme
    virtual std::string Name() { return "Euler"; }
  };


}  // namespace ACOU

#endif /* ACOU_IMPL_EULER_H */
