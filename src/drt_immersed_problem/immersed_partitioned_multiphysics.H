/*!----------------------------------------------------------------------
\file immersed_partitioned_multiphysics.H

\brief base class for all multifield partitioned immersed algorithms

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
 *----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_MULTIPHYSICS_H
#define IMMERSED_PARTITIONED_MULTIPHYSICS_H

#include "immersed_base.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>


namespace IMMERSED {

class ImmersedPartitionedMultiphysics : public NOX::Epetra::Interface::Required

{
public:

  explicit ImmersedPartitionedMultiphysics(Teuchos::ParameterList& params, const Epetra_Comm& comm);

  /// loop over time steps of partitioned scheme
  void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

  /// compute residual of partitioned scheme S^{-1}(F(d)) - d
  bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

};// class ImmersedPartitionedMultiphysics
}// namespace IMMERSED

#endif
