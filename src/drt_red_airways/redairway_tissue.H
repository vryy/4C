
/*!----------------------------------------------------------------------
\file redairway_tissue.H
\brief Control routine for coupled reduced airways and continuum tissue models

<pre>
Maintainer: Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*----------------------------------------------------------------------*/

#ifndef REDAIRWAYTISSUE_H
#define REDAIRWAYTISSUE_H

#include "../drt_adapter/adapter_algorithmbase.H"
#include <Teuchos_RefCountPtr.hpp>
#include "Epetra_Vector.h"


// forward declarations

namespace ADAPTER
{
  class StructureRedAirway;
}
namespace DRT
{
  class Condition;
}

/// airway algorithms
namespace AIRWAY
{
  class RedAirwayImplicitTimeInt;

  class RedAirwayTissue : public ADAPTER::AlgorithmBase
{
public:


  /// Standard Constructor
  RedAirwayTissue(const Epetra_Comm& comm,
                  const Teuchos::ParameterList& timeparams);

  /// Destructor
  virtual ~RedAirwayTissue() { return; };

  virtual void ReadRestart(int step);

  void SetupRedAirways();

  /// time integration of coupled problem
  void Integrate();

  void DoStructureStep();

  void DoRedAirwayStep();

  /// flag whether iteration between fields should be finished
  bool NotConverged();

  void OutputIteration(double pres_inc_norm, double flux_inc_norm);

  void UpdateAndOutput();

private:

  /// underlying structure
  Teuchos::RCP< ::ADAPTER::StructureRedAirway> structure_;

  Teuchos::RCP<RedAirwayImplicitTimeInt> redairways_;

  /// redundant vector of outlet pressures (new iteration step)
  Teuchos::RCP<Epetra_Vector> couppres_ip_;

  /// redundant vector of outlet pressures (old iteration step)
  Teuchos::RCP<Epetra_Vector> couppres_im_;

  /// redundant vector of outlet fluxes (new iteration step)
  Teuchos::RCP<Epetra_Vector> coupflux_ip_;

  /// redundant vector of outlet fluxes (old iteration step)
  Teuchos::RCP<Epetra_Vector> coupflux_im_;

};
}

#endif
