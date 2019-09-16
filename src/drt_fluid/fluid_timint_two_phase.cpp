/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration class for two-phase flow

\maintainer Christoph Ager

\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_two_phase.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       mw 05/14 |
 *----------------------------------------------------------------------*/
FLD::TimIntTwoPhase::TimIntTwoPhase(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      tpf_gradphi_curvaf_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                       mw 07/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhase::Init()
{
  tpf_gradphi_curvaf_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));

  if (convform_ == "conservative")
    dserror(
        "conservative formulation currently not supported for two-phase-flow flow within "
        "generalized-alpha time-integration scheme");

  // set some Two phase-specific parameters
  SetElementCustomParameter();
  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                     mw 07/14 |
*----------------------------------------------------------------------*/
FLD::TimIntTwoPhase::~TimIntTwoPhase() { return; }

/*----------------------------------------------------------------------*
 | set fields for two phase flow within iteration loop         mw 07/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhase::SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
    Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
    Teuchos::RCP<const Epetra_Vector> fsscalaraf, Teuchos::RCP<const Epetra_Vector> curvatureaf,
    Teuchos::RCP<const Epetra_MultiVector> gradphiaf, Teuchos::RCP<DRT::Discretization> scatradis)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  //--------------------------------------------------------------------------
  // Filling the scaaf-vector and scaam-vector at time n+alpha_F/n+1 and
  // n+alpha_M/n, respectively, with scalar at pressure dofs
  // Additionally, filling the scaam-vector at time n+alpha_M/n with
  // velocity at time n at velocity dofs for OST/BDF2
  // Filling the accam-vector at time n+alpha_M/n+1, respectively, with
  // scalar time derivative values at pressure dofs
  //--------------------------------------------------------------------------
  // get velocity values at time n in scaam-vector as copy from veln-vector
  scaam_->Update(1.0, *veln_, 0.0);

  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0, lscatranode);
    if (numscatradof > 1)
      dserror(
          "Number of scalars for level set problems is expected to be 1. Too many dofs in scalar "
          "variable.");
    const int globalscatradofid = scatradis->Dof(0, lscatranode, 0);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0) dserror("localdofid not found in map for given globaldofid");


    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(0, lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int numdof = discret_->NumDof(0, lnode);
    const int globaldofid = discret_->Dof(0, lnode, numdof - 1);
    const int localdofid = tpf_gradphi_curvaf_->Map().LID(globaldofid);
    if (localdofid < 0) dserror("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into scaaf_");

    value = (*scalaram)[localscatradofid];
    err = scaam_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into scaam_");

    if (scalardtam != Teuchos::null)
    {
      value = (*scalardtam)[localscatradofid];
    }
    else
    {
      value = 0.0;  // for safety reasons: set zeros in accam_
    }
    err = accam_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into accam_");

    //=====================================================================

    value = (*curvatureaf)[localscatradofid];
    err = tpf_gradphi_curvaf_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0)
      dserror("error while inserting value into tpf_gradphi_curvaf_, for curvature entry");

    for (int ndimcounter = 0; ndimcounter < 3; ndimcounter++)
    {
      value = (*gradphiaf)[ndimcounter][localscatradofid];
      err = tpf_gradphi_curvaf_->ReplaceMyValue(localdofid - 3 + ndimcounter, 0, value);
      if (err != 0)
        dserror("error while inserting value into tpf_gradphi_curvaf_, for gradphi entry");
    }

    //=====================================================================

    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (fsscalaraf != Teuchos::null)
        value = (*fsscalaraf)[localscatradofid];
      else
        dserror("Expected fine-scale scalar!");

      err = fsscaaf_->ReplaceMyValue(localdofid, 0, value);
      if (err != 0) dserror("error while inserting value into fsscaaf_");
    }
  }
  return;

}  // TimIntTwoPhase::SetIterScalarFields

/*----------------------------------------------------------------------*
 | set scalar fields                                       winter 05/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhase::SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
    Teuchos::RCP<const Epetra_Vector> curvaturenp, Teuchos::RCP<const Epetra_MultiVector> gradphinp,
    Teuchos::RCP<const Epetra_Vector> scatraresidual, Teuchos::RCP<DRT::Discretization> scatradis,
    const int whichscalar)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  //--------------------------------------------------------------------------
  // Filling the scaaf-vector with scalar at time n+1 at pressure dofs
  //--------------------------------------------------------------------------
  // loop all nodes on the processor
  for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0, lscatranode);
    int globalscatradofid(-1);
    if (whichscalar == (-1))
    {
      // default: always take the LAST scatra dof at each node
      globalscatradofid = scatradis->Dof(0, lscatranode, numscatradof - 1);
    }
    else
    {
      // respect the explicit wish of the user
      globalscatradofid = scatradis->Dof(0, lscatranode, whichscalar);
    }
    const int localscatradofid = scalarnp->Map().LID(globalscatradofid);
    if (localscatradofid < 0) dserror("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(0, lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int globaldofid = nodedofs[numdim_];
    const int localdofid = scaam_->Map().LID(globaldofid);
    if (localdofid < 0) dserror("localdofid not found in map for given globaldofid");

    value = (*scalarnp)[localscatradofid];
    err = scaaf_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0) dserror("error while inserting value into scaaf_");

    //=====================================================================

    value = (*curvaturenp)[localscatradofid];
    err = tpf_gradphi_curvaf_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0)
      dserror("error while inserting value into tpf_gradphi_curvaf_, for curvature entry");

    for (int ndimcounter = 0; ndimcounter < 3; ndimcounter++)
    {
      value = (*gradphinp)[ndimcounter][localscatradofid];
      err = tpf_gradphi_curvaf_->ReplaceMyValue(localdofid - 3 + ndimcounter, 0, value);
      if (err != 0)
        dserror("error while inserting value into tpf_gradphi_curvaf_, for gradphi entry");
    }

    //=====================================================================

    //--------------------------------------------------------------------------
    // Filling the trueresidual vector with scatraresidual at pre-dofs
    //--------------------------------------------------------------------------
    if (scatraresidual != Teuchos::null)
    {
      value = (*scatraresidual)[localscatradofid];
      trueresidual_->ReplaceMyValue(localdofid, 0, value);
    }
  }

  return;
}  // TimIntTwoPhase::SetScalarFields

// -------------------------------------------------------------------
// set two phase parameters                               winter 07/2014
// -------------------------------------------------------------------
void FLD::TimIntTwoPhase::SetElementCustomParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", FLD::set_two_phase_parameter);
  eleparams.sublist("SMEARED") = params_->sublist("SMEARED");
  eleparams.sublist("SURFACE TENSION") = params_->sublist("SURFACE TENSION");

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  return;
}  // TimIntTwoPhase::SetElementCustomParameter

/*----------------------------------------------------------------------*
 | Give problem dependent vectors to                           mw 05/14 |
 | discretization                                                       |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhase::AddProblemDependentVectors()
{
  discret_->SetState("tpf_gradphi_curvaf", tpf_gradphi_curvaf_);
  return;
}  // TimIntTwoPhase::AddProblemDependentVectors
