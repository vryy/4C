/*----------------------------------------------------------------------*/
/*! \file
 \brief time integration schemes for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/


#include "scatra_timint_poromulti.H"

#include "../drt_io/io.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_utils.H"

/*----------------------------------------------------------------------*
 | constructor                                             vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMulti::ScaTraTimIntPoroMulti(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output), L2_projection_(false)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
  return;
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::Init() { return; }

/*----------------------------------------------------------------------*
 | set solution fields on given dof sets                    vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::SetL2FluxOfMultiFluid(
    Teuchos::RCP<const Epetra_MultiVector> multiflux, const int nds_flux)
{
  // set L2-projection to true
  L2_projection_ = true;

  // safety check
  if (nds_flux >= discret_->NumDofSets()) dserror("Too few dofsets on scatra discretization!");

  // store number of dof-set associated with velocity related dofs
  nds_vel_ = nds_flux;

  if (multiflux->NumVectors() % nsd_ != 0)
    dserror("Unexpected length of flux vector: %i", multiflux->NumVectors());

  const int totalnumdof = multiflux->NumVectors() / nsd_;

  std::string stateprefix = "flux";

  for (int curphase = 0; curphase < totalnumdof; ++curphase)
  {
    // initialize velocity vectors
    Teuchos::RCP<Epetra_Vector> phaseflux =
        LINALG::CreateVector(*discret_->DofRowMap(nds_flux), true);

    std::stringstream statename;
    statename << stateprefix << curphase;

    // loop all nodes on the processor
    for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
    {
      // get the processor local node
      DRT::Node* lnode = discret_->lRowNode(lnodeid);

      // get dofs associated with current node
      std::vector<int> nodedofs = discret_->Dof(nds_flux, lnode);

      if ((int)nodedofs.size() != nsd_)
        dserror(
            "Expected number of DOFs to be equal to the number of space dimensions for flux "
            "state!");

      for (int index = 0; index < nsd_; ++index)
      {
        // get global and local dof IDs
        const int gid = nodedofs[index];
        const int lid = phaseflux->Map().LID(gid);
        if (lid < 0) dserror("Local ID not found in map for given global ID!");

        const double value = (*(*multiflux)(curphase * nsd_ + index))[lnodeid];

        int err = phaseflux->ReplaceMyValue(lid, 0, value);
        if (err != 0) dserror("error while inserting a value into convel");
      }
    }

    // provide scatra discretization with convective velocity
    discret_->SetState(nds_flux, statename.str(), phaseflux);
  }

  return;

}  // ScaTraTimIntImpl::SetSolutionFields

/*----------------------------------------------------------------------*
 | set solution fields on given dof sets              kremheller  07/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::SetSolutionFieldOfMultiFluid(
    Teuchos::RCP<const Epetra_Vector> phinp_fluid, Teuchos::RCP<const Epetra_Vector> phin_fluid,
    const int nds_phi_fluid)
{
  if (nds_phi_fluid >= discret_->NumDofSets()) dserror("Too few dofsets on scatra discretization!");

  // TODO: this is a hack to allow evaluation of initial time derivative
  //      with check nds_vel_ != -1 since in the case without L2- projection
  //      the velocity field is directly calculated at GPs with the Darcy eqn.
  nds_vel_ = 1;

  // store number of dof-set
  nds_pres_ = nds_phi_fluid;
  // provide scatra discretization with fluid primary variable field
  discret_->SetState(nds_pres_, "phinp_fluid", phinp_fluid);
  discret_->SetState(nds_pres_, "phin_fluid", phin_fluid);
}

/*----------------------------------------------------------------------*
 | add parameters depending on the problem                  vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  // set dof set numbers
  // note: the velocity dof set is set by the standard time integrator

  // provide pressure field
  params.set<int>("ndspres", nds_pres_);
  // provide pressure field
  params.set<bool>("L2-projection", L2_projection_);

  return;
}

/*----------------------------------------------------------------------*
 |  write current state to BINIO                           vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::OutputState()
{
  // solution
  output_->WriteVector("phinp", phinp_);

  //  // convective velocity (written in case of coupled simulations since volmortar is now
  //  possible) if ( cdvel_ == INPAR::SCATRA::velocity_function or cdvel_ ==
  //  INPAR::SCATRA::velocity_Navier_Stokes)
  //  {
  //    Teuchos::RCP<const Epetra_Vector> convel = discret_->GetState(nds_vel_, "convective velocity
  //    field"); if(convel == Teuchos::null)
  //      dserror("Cannot get state vector convective velocity");
  //
  //    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
  //    Teuchos::RCP<Epetra_MultiVector> convel_multi = Teuchos::rcp(new
  //    Epetra_MultiVector(*discret_->NodeRowMap(),nsd_,true)); for (int inode=0;
  //    inode<discret_->NumMyRowNodes(); ++inode)
  //    {
  //      DRT::Node* node = discret_->lRowNode(inode);
  //      for (int idim=0; idim<nsd_; ++idim)
  //        (*convel_multi)[idim][inode] =
  //        (*convel)[convel->Map().LID(discret_->Dof(nds_vel_,node,idim))];
  //    }
  //
  //    output_->WriteVector("convec_velocity", convel_multi, IO::nodevector);
  //  }

  // displacement field
  if (isale_)
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = discret_->GetState(nds_disp_, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot extract displacement field from discretization");

    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
    Teuchos::RCP<Epetra_MultiVector> dispnp_multi =
        Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), nsd_, true));
    for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
    {
      DRT::Node* node = discret_->lRowNode(inode);
      for (int idim = 0; idim < nsd_; ++idim)
        (*dispnp_multi)[idim][inode] =
            (*dispnp)[dispnp->Map().LID(discret_->Dof(nds_disp_, node, idim))];
    }

    output_->WriteVector("dispnp", dispnp_multi, IO::nodevector);
  }

  return;
}  // ScaTraTimIntImpl::OutputState

/*----------------------------------------------------------------------*
 | problem specific output                             kremheller 10/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::OutputProblemSpecific()
{
  // oxygen partial pressure (if desired)
  OutputOxygenPartialPressure();

  return;
}

/*----------------------------------------------------------------------*
 | output of oxygen partial pressure                   kremheller 10/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMulti::OutputOxygenPartialPressure()
{
  // extract conditions for oxygen partial pressure
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition("PoroMultiphaseScatraOxyPartPressCalcCond", conditions);

  // perform all following operations only if there is at least one condition for oxygen partial
  // pressure
  if (conditions.size() > 0)
  {
    const Teuchos::RCP<Epetra_Vector> oxypartpress =
        Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(), true));

    // this condition is supposed to be for output of oxygen partial pressure over whole domain
    // it does not make sense to have more than one condition
    if (conditions.size() != 1)
      dserror("Should have only one PoroMultiphaseScatraOxyPartPressCalcCond per discretization");

    // extract nodal cloud from condition
    const std::vector<int>* nodegids = conditions[0]->Nodes();

    // output
    double Pb = 0.0;

    // read input from condition
    const int oxyscalar = conditions[0]->GetInt("SCALARID") - 1;
    const double CaO2_max = conditions[0]->GetDouble("CaO2_max");
    const double Pb50 = conditions[0]->GetDouble("Pb50");
    const double n = conditions[0]->GetDouble("n");
    const double alpha_eff = conditions[0]->GetDouble("alpha_bl_eff");
    const double rho_oxy = conditions[0]->GetDouble("rho_oxy");
    const double rho_bl = conditions[0]->GetDouble("rho_bl");

    // loop over all nodes
    for (unsigned inode = 0; inode < nodegids->size(); ++inode)
    {
      // extract global ID of current node
      const int nodegid((*nodegids)[inode]);
      // process only nodes stored by current processor
      if (discret_->HaveGlobalNode(nodegid))
      {
        // extract current node
        const DRT::Node* const node = discret_->gNode(nodegid);

        // process only nodes owned by current processor
        if (node->Owner() == discret_->Comm().MyPID())
        {
          // get dof
          int myoxydof = discret_->Dof(0, node, oxyscalar);
          const int lidoxydof = discret_->DofRowMap()->LID(myoxydof);
          if (lidoxydof < 0) dserror("Couldn't extract local ID of oxygen dof!");
          // compute CaO2
          const double CaO2 = (*phinp_)[lidoxydof] * rho_bl / rho_oxy;
          // compute Pb
          POROMULTIPHASESCATRA::UTILS::GetOxyPartialPressureFromConcentration<double>(
              Pb, CaO2, CaO2_max, Pb50, n, alpha_eff);
          // replace value
          oxypartpress->ReplaceGlobalValue(node->Id(), 0, Pb);
        }
      }
    }
    output_->WriteVector("oxypartpress", oxypartpress, IO::nodevector);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiOST::ScaTraTimIntPoroMultiOST(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                  vuong  08/16 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiOST::~ScaTraTimIntPoroMultiOST() { return; }

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiOST::Update(const int num)
{
  TimIntOneStepTheta::Update(num);
  ScaTraTimIntPoroMulti::Update(num);

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiBDF2::ScaTraTimIntPoroMultiBDF2(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntBDF2(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiBDF2::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntBDF2::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                  vuong  08/16 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiBDF2::~ScaTraTimIntPoroMultiBDF2() { return; }

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiBDF2::Update(const int num)
{
  TimIntBDF2::Update(num);
  ScaTraTimIntPoroMulti::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiGenAlpha::ScaTraTimIntPoroMultiGenAlpha(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiGenAlpha::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                  vuong  08/16 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiGenAlpha::~ScaTraTimIntPoroMultiGenAlpha() { return; }

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiGenAlpha::Update(const int num)
{
  TimIntGenAlpha::Update(num);
  ScaTraTimIntPoroMulti::Update(num);

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    vuong  08/16 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiStationary::ScaTraTimIntPoroMultiStationary(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      ScaTraTimIntPoroMulti(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntStationary(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiStationary::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntStationary::Init();
  ScaTraTimIntPoroMulti::Init();

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                  vuong  08/16 |
*-----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntPoroMultiStationary::~ScaTraTimIntPoroMultiStationary() { return; }

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                         vuong  08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntPoroMultiStationary::Update(const int num)
{
  TimIntStationary::Update(num);
  ScaTraTimIntPoroMulti::Update(num);

  return;
}
