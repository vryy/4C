/*-----------------------------------------------------------*/
/*! \file

\brief Stationary fluid problem with HDG discretization

\maintainer Andrea La Spina

\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_stat_hdg.H"
#include "fluid_volumetric_surfaceFlow_condition.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_turbulence/dyn_smag.H"
#include "../drt_fluid_turbulence/dyn_vreman.H"
#include "../drt_fluid_turbulence/boxfilter.H"
#include "../drt_fluid/fluid_utils.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_ele/fluid_ele_hdg.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                      als 01/18 |    // TODO als fix
 fluid_timint_stat_hdg because it is not working
 *----------------------------------------------------------------------*/
FLD::TimIntStationaryHDG::TimIntStationaryHDG(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntStationary(actdis, solver, params, output, alefluid),
      firstAssembly_(false)
{
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                      als 01/18 |
 *----------------------------------------------------------------------*/
void FLD::TimIntStationaryHDG::Init()
{
  DRT::DiscretizationHDG* hdgdis = dynamic_cast<DRT::DiscretizationHDG*>(discret_.get());
  if (hdgdis == NULL) dserror("Did not receive an HDG discretization");

  int elementndof = hdgdis->NumMyRowElements() > 0
                        ? dynamic_cast<DRT::ELEMENTS::FluidHDG*>(hdgdis->lRowElement(0))
                              ->NumDofPerElementAuxiliary()
                        : 0;

  // set degrees of freedom in the discretization
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0, elementndof, 0, false));
  discret_->AddDofSet(dofsetaux);
  discret_->FillComplete();

  // build velocity/pressure splitting
  std::set<int> conddofset;
  std::set<int> otherdofset;

  for (int j = 0; j < hdgdis->NumMyRowElements(); ++j)
  {
    std::vector<int> dof = hdgdis->Dof(0, hdgdis->lRowElement(j));
    dsassert(dof.size() >= 1, "Internal error: could not find HDG pressure dof");
    for (unsigned int i = 0; i < dof.size(); ++i) conddofset.insert(dof[i]);
  }
  for (int i = 0; i < hdgdis->NumMyRowFaces(); ++i)
  {
    std::vector<int> dof = hdgdis->Dof(0, hdgdis->lRowFace(i));
    for (unsigned int j = 0; j < dof.size(); ++j) otherdofset.insert(dof[j]);
  }

  std::vector<int> conddofmapvec;
  conddofmapvec.reserve(conddofset.size());
  conddofmapvec.assign(conddofset.begin(), conddofset.end());
  conddofset.clear();
  Teuchos::RCP<Epetra_Map> conddofmap =
      Teuchos::rcp(new Epetra_Map(-1, conddofmapvec.size(), &conddofmapvec[0], 0, hdgdis->Comm()));
  std::vector<int> otherdofmapvec;
  otherdofmapvec.reserve(otherdofset.size());
  otherdofmapvec.assign(otherdofset.begin(), otherdofset.end());
  otherdofset.clear();
  Teuchos::RCP<Epetra_Map> otherdofmap = Teuchos::rcp(
      new Epetra_Map(-1, otherdofmapvec.size(), &otherdofmapvec[0], 0, hdgdis->Comm()));
  velpressplitter_->Setup(*hdgdis->DofRowMap(), conddofmap, otherdofmap);

  // call Init()-functions of base classes
  // note: this order is important
  FLD::TimIntStationary::Init();
}


void FLD::TimIntStationaryHDG::Reset(bool completeReset, int numsteps, int iter)
{
  FluidImplicitTimeInt::Reset(completeReset, numsteps, iter);
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  intvelnp_ = LINALG::CreateVector(*intdofrowmap, true);
  if (discret_->Comm().MyPID() == 0)
    std::cout << "Number of degrees of freedom in HDG system: "
              << discret_->DofRowMap(0)->NumGlobalElements() << std::endl;
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                       als 01/18 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationaryHDG::SetCustomEleParamsAssembleMatAndRHS(
    Teuchos::ParameterList& eleparams)
{
  eleparams.set<bool>("needslocalupdate", !firstAssembly_);
}


/*----------------------------------------------------------------------*
| set old part of right hand side                             als 01/18 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationaryHDG::SetOldPartOfRighthandside()
{
  /*
     Stationary:

                   mom: hist_ = 0.0
                  (con: hist_ = 0.0)
  */

  hist_->PutScalar(0.0);

  // This code is entered at the beginning of the nonlinear iteration, so
  // store that the assembly to be done next is going to be the first one
  // (without combined vector update) for HDG.
  firstAssembly_ = true;
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                       als 01/18 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationaryHDG::SetStateTimInt()
{
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  Teuchos::RCP<Epetra_Vector> zerovec = LINALG::CreateVector(*intdofrowmap, true);

  discret_->SetState(0, "velaf", velnp_);
  discret_->SetState(1, "intvelaf", intvelnp_);  // TODO als fill in intvelnp_!
  discret_->SetState(1, "intaccam", zerovec);
  discret_->SetState(1, "intvelnp", intvelnp_);
}

/*----------------------------------------------------------------------*
| set integration-scheme-specific state                       als 01/18 |
*-----------------------------------------------------------------------*/
void FLD::TimIntStationaryHDG::ClearStateAssembleMatAndRHS()
{
  if (!firstAssembly_)
  {
    // Wrote into the state vector during element calls, need to transfer the
    // data back before it disappears when clearing the state (at least for nproc>1)
    const Epetra_Vector& intvelnpGhosted = *discret_->GetState(1, "intvelnp");
    for (int i = 0; i < intvelnp_->MyLength(); ++i)
      (*intvelnp_)[i] = intvelnpGhosted[intvelnpGhosted.Map().LID(intvelnp_->Map().GID(i))];
  }
  firstAssembly_ = false;
  FluidImplicitTimeInt::ClearStateAssembleMatAndRHS();
}

/*----------------------------------------------------------------------*
 |  set initial flow field for test cases              kronbichler 05/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntStationaryHDG::SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield, const int startfuncno)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;
  Teuchos::ParameterList initParams;
  initParams.set<int>("action", FLD::project_fluid_field);
  initParams.set("startfuncno", startfuncno);
  initParams.set<int>("initfield", initfield);
  // loop over all elements on the processor
  DRT::Element::LocationArray la(2);
  double error = 0;
  for (int el = 0; el < discret_->NumMyColElements(); ++el)
  {
    DRT::Element* ele = discret_->lColElement(el);

    ele->LocationVector(*discret_, la, false);
    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1, ele)) elevec2.Shape(discret_->NumDof(1, ele), 1);

    ele->Evaluate(initParams, *discret_, la[0].lm_, elemat1, elemat2, elevec1, elevec2, elevec3);

    // now fill the ele vector into the discretization
    for (unsigned int i = 0; i < la[0].lm_.size(); ++i)
    {
      const int lid = dofrowmap->LID(la[0].lm_[i]);
      if (lid >= 0)
      {
        if ((*velnp_)[lid] != 0) error += std::abs((*velnp_)[lid] - elevec1(i));
        (*velnp_)[lid] = elevec1(i);
        (*veln_)[lid] = elevec1(i);
        (*velnm_)[lid] = elevec1(i);
      }
    }

    if (ele->Owner() == discret_->Comm().MyPID())
    {
      std::vector<int> localDofs = discret_->Dof(1, ele);
      dsassert(localDofs.size() == static_cast<std::size_t>(elevec2.M()), "Internal error");
      for (unsigned int i = 0; i < localDofs.size(); ++i)
        localDofs[i] = intdofrowmap->LID(localDofs[i]);
      intvelnp_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
    }
  }
  double globerror = 0;
  discret_->Comm().SumAll(&error, &globerror, 1);
  if (discret_->Comm().MyPID() == 0)
    std::cout << "Error project when setting face twice: " << globerror << std::endl;
}


// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::TimIntStationaryHDG::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_time_parameter);
  eleparams.set<int>("Physical Type", physicaltype_);

  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  // set general element parameters
  eleparams.set("dt", dta_);
  eleparams.set("theta", theta_);
  eleparams.set("omtheta", 0.0);

  // set scheme-specific element parameters and vector values
  eleparams.set("total time", time_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}
