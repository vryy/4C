/*----------------------------------------------------------------------*/
/*! \file
\brief structural model evaluator for monolithic scalar-structure interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#include "ssi_str_model_evaluator_monolithic.H"

#include <Epetra_IntVector.h>

#include "../drt_adapter/adapter_coupling.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_io/io.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"

#include "../drt_ssi/ssi_monolithic.H"

#include "../drt_structure_new/str_model_evaluator_data.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"


/*----------------------------------------------------------------------*
 | constructor                                               fang 12/17 |
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::MonolithicSSI::MonolithicSSI(const Teuchos::RCP<const SSI::SSI_Mono>
        ssi_mono  //!< monolithic algorithm for scalar-structure interaction
    )
    : stresses_(Teuchos::null), ssi_mono_(ssi_mono)
{
  return;
}


/*----------------------------------------------------------------------*
 | calculate stresses and strains                            fang 12/17 |
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::MonolithicSSI::DetermineStressStrain()
{
  // extract raw data for element-wise stresses
  const std::vector<char>& stressdata = EvalData().StressData();

  // initialize map for element-wise stresses
  const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> stresses =
      Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>);

  // initialize position pointer
  std::vector<char>::size_type position(0);

  // loop over all row elements
  for (int i = 0; i < Discret().ElementRowMap()->NumMyElements(); ++i)
  {
    // initialize matrix for stresses associated with current element
    const Teuchos::RCP<Epetra_SerialDenseMatrix> stresses_ele =
        Teuchos::rcp(new Epetra_SerialDenseMatrix);

    // extract stresses
    DRT::ParObject::ExtractfromPack(position, stressdata, *stresses_ele);

    // store stresses
    (*stresses)[Discret().ElementRowMap()->GID(i)] = stresses_ele;
  }

  // export map to column format
  DRT::Exporter exporter(*Discret().ElementRowMap(),
      *Teuchos::rcp_dynamic_cast<DRT::Discretization>(DiscretPtr())->ElementColMap(),
      Discret().Comm());
  exporter.Export(*stresses);

  // evaluate nodal stresses
  stresses_->PutScalar(0.);
  Teuchos::ParameterList parameters;
  parameters.set("action", "postprocess_stress");
  parameters.set("gpstressmap", stresses);
  parameters.set("stresstype", "ndxyz");
  parameters.set("poststress", stresses_);
  Discret().Evaluate(
      parameters, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // initialize vectors for stresses and element numbers on slave and master sides of scatra-scatra
  // coupling interfaces
  const Teuchos::RCP<Epetra_MultiVector> stresses_slave(Teuchos::rcp(
      new Epetra_MultiVector(*ssi_mono_->CouplingAdapterStructure()->SlaveDofMap(), 2))),
      stresses_master(Teuchos::rcp(
          new Epetra_MultiVector(*ssi_mono_->CouplingAdapterStructure()->MasterDofMap(), 2)));
  Epetra_IntVector numelement_slave(*ssi_mono_->CouplingAdapterStructure()->SlaveDofMap()),
      numelement_master(*ssi_mono_->CouplingAdapterStructure()->MasterDofMap());

  // extract scatra-scatra interface coupling conditions
  std::vector<DRT::Condition*> conditions;
  Discret().GetCondition("S2ICoupling", conditions);

  // loop over all scatra-scatra interface coupling conditions
  for (unsigned icond = 0; icond < conditions.size(); ++icond)
  {
    // extract interface side
    const int side = conditions[icond]->GetInt("interface side");

    // set references depending on interface side
    Epetra_MultiVector& stresses_vector =
        side == INPAR::S2I::side_slave ? *stresses_slave : *stresses_master;
    Epetra_IntVector& numelement_vector =
        side == INPAR::S2I::side_slave ? numelement_slave : numelement_master;

    // extract nodal cloud of current condition
    const std::vector<int>* const nodegids = conditions[icond]->Nodes();
    if (!nodegids or !nodegids->size())
      dserror("Scatra-scatra interface coupling condition does not have a nodal cloud!");

    // loop over all nodes
    for (unsigned inode = 0; inode < nodegids->size(); ++inode)
    {
      // extract global ID of current node
      const int nodegid = (*nodegids)[inode];

      // process only nodes stored on calling processor
      if (Discret().HaveGlobalNode(nodegid))
      {
        // extract node
        const DRT::Node* const node = Discret().gNode(nodegid);
        if (node == NULL) dserror("Couldn't find node!");

        // process only nodes owned by calling processor
        if (node->Owner() == Discret().Comm().MyPID())
        {
          // extract local ID of current node
          const int nodelid = Discret().NodeRowMap()->LID(nodegid);
          if (nodelid < 0) dserror("Local ID not found!");

          // extract global ID of first degree of freedom associated with current node
          const int dofgid = Discret().Dof(0, node, 0);

          // extract corresponding local ID
          const int doflid = stresses_vector.Map().LID(dofgid);
          if (doflid < 0) dserror("Local ID not found!");

          // extract number of elements adjacent to current node
          const int numelement = node->NumElement();

          // scale and store stresses associated with current node
          (*stresses_vector(0))[doflid] = (*(*stresses_)(0))[nodelid] * numelement;
          (*stresses_vector(0))[doflid + 1] = (*(*stresses_)(1))[nodelid] * numelement;
          (*stresses_vector(0))[doflid + 2] = (*(*stresses_)(2))[nodelid] * numelement;
          (*stresses_vector(1))[doflid] = (*(*stresses_)(3))[nodelid] * numelement;
          (*stresses_vector(1))[doflid + 1] = (*(*stresses_)(4))[nodelid] * numelement;
          (*stresses_vector(1))[doflid + 2] = (*(*stresses_)(5))[nodelid] * numelement;

          // store number of elements adjacent to current node
          numelement_vector[doflid] = numelement;
          numelement_vector[doflid + 1] = numelement;
          numelement_vector[doflid + 2] = numelement;
        }
      }
    }
  }

  // communicate vectors from master to slave side
  const Teuchos::RCP<Epetra_MultiVector> stresses_temp(Teuchos::rcp(
      new Epetra_MultiVector(*ssi_mono_->CouplingAdapterStructure()->SlaveDofMap(), 2)));
  Epetra_IntVector numelement_temp(*ssi_mono_->CouplingAdapterStructure()->SlaveDofMap());
  ssi_mono_->CouplingAdapterStructure()->MasterToSlave(stresses_master, stresses_temp);
  ssi_mono_->CouplingAdapterStructure()->MasterToSlave(numelement_master, numelement_temp);

  // add stresses together
  stresses_slave->Update(1., *stresses_temp, 1.);

  // unscale stresses
  for (int lid = 0; lid < stresses_slave->MyLength(); ++lid)
    for (int ivec = 0; ivec < 2; ++ivec)
      (*(*stresses_slave)(ivec))[lid] /= numelement_slave[lid] + numelement_temp[lid];

  // communicate back
  ssi_mono_->CouplingAdapterStructure()->SlaveToMaster(stresses_slave, stresses_master);

  // loop over all scatra-scatra interface coupling conditions
  for (unsigned icond = 0; icond < conditions.size(); ++icond)
  {
    // extract interface side
    const int side = conditions[icond]->GetInt("interface side");

    // set reference depending on interface side
    Epetra_MultiVector& stresses_vector =
        side == INPAR::S2I::side_slave ? *stresses_slave : *stresses_master;

    // extract nodal cloud of current condition
    const std::vector<int>* const nodegids = conditions[icond]->Nodes();
    if (!nodegids or !nodegids->size())
      dserror("Scatra-scatra interface coupling condition does not have a nodal cloud!");

    // loop over all nodes
    for (unsigned inode = 0; inode < nodegids->size(); ++inode)
    {
      // extract global ID of current node
      const int nodegid = (*nodegids)[inode];

      // process only nodes stored on calling processor
      if (Discret().HaveGlobalNode(nodegid))
      {
        // extract node
        const DRT::Node* const node = Discret().gNode(nodegid);
        if (node == NULL) dserror("Couldn't find node!");

        // process only nodes owned by calling processor
        if (node->Owner() == Discret().Comm().MyPID())
        {
          // extract local ID of current node
          const int nodelid = Discret().NodeRowMap()->LID(nodegid);
          if (nodelid < 0) dserror("Local ID not found!");

          // extract global ID of first degree of freedom associated with current node
          const int dofgid = Discret().Dof(0, node, 0);

          // extract corresponding local ID
          const int doflid = stresses_vector.Map().LID(dofgid);
          if (doflid < 0) dserror("Local ID not found!");

          // write back stresses associated with current node
          (*(*stresses_)(0))[nodelid] = (*stresses_vector(0))[doflid];
          (*(*stresses_)(1))[nodelid] = (*stresses_vector(0))[doflid + 1];
          (*(*stresses_)(2))[nodelid] = (*stresses_vector(0))[doflid + 2];
          (*(*stresses_)(3))[nodelid] = (*stresses_vector(1))[doflid];
          (*(*stresses_)(4))[nodelid] = (*stresses_vector(1))[doflid + 1];
          (*(*stresses_)(5))[nodelid] = (*stresses_vector(1))[doflid + 2];
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | return dofrowmap                                          fang 12/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::MonolithicSSI::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}


/*----------------------------------------------------------------------*
 | output state                                              fang 12/17 |
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::MonolithicSSI::OutputStepState(IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

  // write nodal stresses
  iowriter.WriteVector("nodal_stresses_xyz", stresses_, IO::nodevector);

  return;
}


/*----------------------------------------------------------------------*
 | set up model evaluator                                    fang 12/17 |
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::MonolithicSSI::Setup()
{
  // check initialization
  CheckInit();

  // create vector for nodal stresses
  stresses_ = Teuchos::rcp(new Epetra_MultiVector(*Discret().NodeRowMap(), 6));

  // set flag
  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------*
 | write model-specific restart                              fang 12/17 |
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::MonolithicSSI::WriteRestart(
    IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
{
  // write nodal stresses
  iowriter.WriteVector("nodal_stresses_xyz", stresses_, IO::nodevector);

  return;
}
