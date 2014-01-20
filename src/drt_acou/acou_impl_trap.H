/*!----------------------------------------------------------------------
\file acou_impl_trap.H
\brief

<pre>
Maintainers: Svenja Schoeder
             schoeder@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#ifndef ACOU_IMPL_TRAP_H
#define ACOU_IMPL_TRAP_H

#include "acou_impl.H"

namespace ACOU
{

class TimIntImplTrap : public AcouImplicitTimeInt
{

public:
  /// constructor
  TimIntImplTrap(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      );

  /// return the name of the time integration scheme
  virtual std::string Name() { return "trapez" ; }

};


} // namespace ACOU

#endif /* ACOU_IMPL_TRAP_H */
