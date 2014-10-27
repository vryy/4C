/*----------------------------------------------------------------------*/
/*!
\file ad_str_invana.cpp

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_str_invana.H"
#include "../drt_lib/drt_utils_timintmstep.H"
#include "../drt_inpar/inpar_structure.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureInvana::StructureInvana(Teuchos::RCP<Structure> structure):
  StructureTimeLoop(structure),
stime_(Teuchos::null),
sdis_(Teuchos::null),
svel_(Teuchos::null),
sizestorage_(0)
{
  stime_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<double>(0, 0, 0.0));
  sdis_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), false));
  svel_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), false));
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void ADAPTER::StructureInvana::ResizeStorage()
{
  // resize storage
  stime_->Resize(-(sizestorage_-1), 0, 0.0);
  // -> time is already stored for the next time step
  sdis_->Resize(-(sizestorage_-1), 0, DofRowMapView(), false);
  svel_->Resize(-(sizestorage_-1), 0, DofRowMapView(), false);

  return;
}

void ADAPTER::StructureInvana::PostPredict()
{
  sizestorage_+=1;

  ResizeStorage();

  return;
}

void ADAPTER::StructureInvana::PostUpdate()
{
  // state and time was already updated so Dispn() and Veln()
  // (pointing to "(*state_)(0)") give the results of
  // this time step
  sdis_->UpdateSteps(*Dispn());
  svel_->UpdateSteps(*Veln());
  stime_->UpdateSteps(TimeOld());

  return;
}

