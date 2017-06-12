/*!----------------------------------------------------------------------
\file springdashpot_manager.cpp

\brief Methods for spring and dashpot constraints / boundary conditions:

<pre>
\level 2

\maintainer Martin Pfaller
</pre>

*----------------------------------------------------------------------*/

#include <iostream>
#include "../drt_io/io_pstream.H" // has to go before io.H
#include "../drt_io/io.H"

#include "springdashpot_manager.H"
#include "springdashpot.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_globalproblem.H"

UTILS::SpringDashpotManager::SpringDashpotManager(Teuchos::RCP<DRT::Discretization> dis):
actdisc_(dis),
havespringdashpot_(false)
{
  // get all spring dashpot conditions
  std::vector<Teuchos::RCP<DRT::Condition> > springdashpots;
  actdisc_->GetCondition("RobinSpringDashpot",springdashpots);

  // number of spring dashpot conditions
  n_conds_ = (int)springdashpots.size();

  // check if spring dashpots are present
  if (n_conds_)
  {
    // yes! we have spring dashpot conditions!
    havespringdashpot_ = true;

    // new instance of spring dashpot BC with current condition for every spring dashpot condition
    for (int i=0; i<n_conds_; ++i)
      springs_.push_back(Teuchos::rcp(new SpringDashpot(actdisc_, springdashpots[i])));
  }

  return;
}

void UTILS::SpringDashpotManager::StiffnessAndInternalForces(
    Teuchos::RCP<LINALG::SparseOperator> stiff,
    Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln,
    Teuchos::ParameterList parlist)
{
  // evaluate all spring dashpot conditions
  for (int i=0; i<n_conds_; ++i)
  {
    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if(stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
      springs_[i]->EvaluateRobin(stiff, fint, disn, veln, parlist);
    if(stype == UTILS::SpringDashpot::cursurfnormal)
      springs_[i]->Evaluate(stiff, fint, disn, veln, parlist);
  }

  return;
}

void UTILS::SpringDashpotManager::Update()
{
  // update all spring dashpot conditions for each new time step
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->Update();

  return;
}

void UTILS::SpringDashpotManager::ResetPrestress(
    Teuchos::RCP<Epetra_Vector> dis)
{
  // loop over all spring dashpot conditions and reset them
  for (int i=0; i<n_conds_; ++i)
    springs_[i]->Reset(dis);

  return;
}

void UTILS::SpringDashpotManager::Output(
    Teuchos::RCP<IO::DiscretizationWriter> output,
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<Epetra_Vector> disp)
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> gap = Teuchos::rcp(new Epetra_Vector(*(actdisc_->NodeRowMap()),true));
  Teuchos::RCP<Epetra_MultiVector> normals = Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->NodeRowMap()),3,true));
  Teuchos::RCP<Epetra_MultiVector> springstress = Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->NodeRowMap()),3,true));

  // collect outputs from all spring dashpot conditions
  bool found_cursurfnormal = false;
  for (int i=0; i<n_conds_; ++i)
  {
    springs_[i]->OutputGapNormal(gap, normals, springstress);

    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();
    if(stype == UTILS::SpringDashpot::cursurfnormal)
      found_cursurfnormal = true;
  }

  // write vectors to output
  if (found_cursurfnormal)
  {
    output->WriteVector("gap", gap);
    output->WriteVector("curnormals", normals);
  }

  // write spring stress if defined in io-flag
  if (DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(),"OUTPUT_SPRING") == true)
    output->WriteVector("springstress", springstress);

  return;
}

void UTILS::SpringDashpotManager::OutputRestart(
    Teuchos::RCP<IO::DiscretizationWriter> output,
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<Epetra_Vector> disp)
{

  // row maps for export
  Teuchos::RCP<Epetra_Vector> springoffsetprestr = Teuchos::rcp(new Epetra_Vector(*actdisc_->DofRowMap()));
  Teuchos::RCP<Epetra_MultiVector> springoffsetprestr_old = Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->NodeRowMap()),3,true));

  // collect outputs from all spring dashpot conditions
  for (int i=0; i<n_conds_; ++i)
  {
    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if(stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
      springs_[i]->OutputPrestrOffset(springoffsetprestr);
    if(stype == UTILS::SpringDashpot::cursurfnormal)
      springs_[i]->OutputPrestrOffsetOld(springoffsetprestr_old);
  }

  // write vector to output for restart
  output->WriteVector("springoffsetprestr", springoffsetprestr);
  // write vector to output for restart
  output->WriteVector("springoffsetprestr_old", springoffsetprestr_old);

  // normal output as well
  Output(output, discret, disp);

  return;
}


/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpotManager::ReadRestart(IO::DiscretizationReader& reader,const double& time)
{

  Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(*actdisc_->DofRowMap()));
  Teuchos::RCP<Epetra_MultiVector> tempvecold = Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->NodeRowMap()),3,true));

  reader.ReadVector(tempvec, "springoffsetprestr");
  reader.ReadMultiVector(tempvecold, "springoffsetprestr_old");

  // loop over all spring dashpot conditions and set restart
  for (int i=0; i<n_conds_; ++i)
  {
    // get spring type from current condition
    const UTILS::SpringDashpot::SpringType stype = springs_[i]->GetSpringType();

    if(stype == UTILS::SpringDashpot::xyz or stype == UTILS::SpringDashpot::refsurfnormal)
      springs_[i]->SetRestart(tempvec);
    if(stype == UTILS::SpringDashpot::cursurfnormal)
      springs_[i]->SetRestartOld(tempvecold);
  }

  return;

}


