  /*!----------------------------------------------------------------------
\file fsi_partitioned_immersed.cpp

\brief partitioned immersed fsi subclass

\level 1

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#include "fsi_partitioned_immersed.H"
#include "../drt_lib/drt_globalproblem.H"


FSI::PartitionedImmersed::PartitionedImmersed(const Epetra_Comm& comm)
  : Partitioned(comm)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
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

