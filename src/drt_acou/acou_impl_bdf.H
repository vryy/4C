/*!----------------------------------------------------------------------
\file acou_impl_bdf.H
\brief

<pre>
Maintainers: Svenja Schoeder
             schoeder@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#ifndef ACOU_IMPL_BDF_H
#define ACOU_IMPL_BDF_H

#include "acou_impl.H"

namespace ACOU
{

class TimIntImplBDF : public AcouImplicitTimeInt
{

public:
  /// constructor
  TimIntImplBDF(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      );

  /// read restart information for given time step
  virtual void ReadRestart(int step);

  /// write restart
  virtual void WriteRestart();

  /// return the name of the time integration scheme
  virtual std::string Name();

private:
  /// order of the BDF scheme
  int order_;
};

} // namespace ACOU

#endif /* ACOU_IMPL_BDF_H */
