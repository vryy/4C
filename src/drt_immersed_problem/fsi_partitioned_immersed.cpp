  /*!----------------------------------------------------------------------
\file fsi_partitioned_immersed.cpp

\brief partitioned immersed fsi subclass

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "fsi_partitioned_immersed.H"
#include "../drt_lib/drt_globalproblem.H"

// relaxation related
#include "../drt_fsi/fsi_nox_aitken.H"
#include "../drt_fsi/fsi_nox_fixpoint.H"

#include "../drt_inpar/inpar_fsi.H"

FSI::PartitionedImmersed::PartitionedImmersed(const Epetra_Comm& comm)
  : Partitioned(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  SetupCoupling(fsidyn,comm);
}


void FSI::PartitionedImmersed::SetupCoupling(const Teuchos::ParameterList& fsidyn ,const Epetra_Comm& comm)
{
  coupsfm_ = Teuchos::null;
  matchingnodes_ = false;

}


void FSI::PartitionedImmersed::ExtractPreviousInterfaceSolution()
{
  // not necessary in immersed fsi.
  // overrides version in fsi_paritioned with "do nothing".
}

