/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for two-phase flow


\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_two_phase_ost.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntTwoPhaseOst::TimIntTwoPhaseOst(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      TimIntOneStepTheta(actdis, solver, params, output, alefluid),
      TimIntTwoPhase(actdis, solver, params, output, alefluid),
      tpf_gradphi_curvn_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhaseOst::Init()
{
  tpf_gradphi_curvn_ = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();
  TimIntTwoPhase::Init();

  return;
}

/*----------------------------------------------------------------------*
 | set fields for two phase flow within iteration loop         mw 03/15 |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhaseOst::SetIterScalarFieldsn(Teuchos::RCP<const Epetra_Vector> curvaturen,
    Teuchos::RCP<const Epetra_MultiVector> gradphin, Teuchos::RCP<DRT::Discretization> scatradis)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;


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
    const int localscatradofid = curvaturen->Map().LID(globalscatradofid);
    if (localscatradofid < 0) dserror("localdofid not found in map for given globaldofid");

    // get the processor's local fluid node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);
    // get the global ids of degrees of freedom associated with this node
    nodedofs = discret_->Dof(0, lnode);
    // get global and processor's local pressure dof id (using the map!)
    const int numdof = discret_->NumDof(0, lnode);
    const int globaldofid = discret_->Dof(0, lnode, numdof - 1);
    const int localdofid = tpf_gradphi_curvn_->Map().LID(globaldofid);
    if (localdofid < 0) dserror("localdofid not found in map for given globaldofid");

    //=====================================================================

    value = (*curvaturen)[localscatradofid];
    err = tpf_gradphi_curvn_->ReplaceMyValue(localdofid, 0, value);
    if (err != 0)
      dserror("error while inserting value into tpf_gradphi_curvn_, for curvature entry");


    for (int ndimcounter = 0; ndimcounter < 3; ndimcounter++)
    {
      value = (*gradphin)[ndimcounter][localscatradofid];
      err = tpf_gradphi_curvn_->ReplaceMyValue(localdofid - 3 + ndimcounter, 0, value);
      if (err != 0)
        dserror("error while inserting value into tpf_gradphi_curvn_, for gradphi entry");
    }

    //=====================================================================
  }
  return;

}  // TimIntTwoPhase::SetIterScalarFieldsn

/*----------------------------------------------------------------------*
 | Give problem dependent vectors to                           mw 03/15 |
 | discretization                                                       |
 *----------------------------------------------------------------------*/
void FLD::TimIntTwoPhaseOst::AddProblemDependentVectors()
{
  discret_->SetState("tpf_gradphi_curvaf", tpf_gradphi_curvaf_);
  discret_->SetState("tpf_gradphi_curvn", tpf_gradphi_curvn_);
  return;
}  // TimIntTwoPhase::AddProblemDependentVectors


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntTwoPhaseOst::~TimIntTwoPhaseOst() { return; }
