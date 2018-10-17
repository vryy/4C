/*----------------------------------------------------------------------*/
/*!
 * \file particle_comm.cpp
 * \brief Communication management for particle data
 *
<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/
#include "particle_comm.H"

#include "Epetra_MpiComm.h"


/*----------------------------------------------------------------------*/
INVANA::ParticleComm::ParticleComm() : lnumparticles_(0), group_(-1) {}

/*----------------------------------------------------------------------*/
void INVANA::ParticleComm::Init(
    Teuchos::RCP<Epetra_Comm> gcomm, Teuchos::RCP<Epetra_Comm> lcomm, int lnumparticles, int group)
{
  gcomm_ = gcomm;
  lcomm_ = lcomm;
  lnumparticles_ = lnumparticles;
  group_ = group;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleComm::Setup()
{
  // color all the i-th procs
  int color = MPI_UNDEFINED;
  gcomm_->Barrier();

  color = lcomm_->MyPID() + 1;

  // Split Comm for inter group communication
  MPI_Comm intercomm;
  Epetra_MpiComm* mpicomm = dynamic_cast<Epetra_MpiComm*>(gcomm_.get());
  if (!mpicomm) dserror("dyncast failed");
  MPI_Comm_split(mpicomm->Comm(), color, gcomm_->MyPID(), &intercomm);
  gcomm_->Barrier();

  // With this we can comunicate between all proc i
  icomm_ = Teuchos::rcp(new Epetra_MpiComm(intercomm));

  // at the moment the global ranks are contiguous within each group
  // and ascending accros groups from group 0 to group N. Given this
  // the rank in icomm_ is the groupID. Check it since I rely on it.
  if (icomm_->MyPID() != group_) dserror("Mess up in the intergroup communicator!");

  return;
}
