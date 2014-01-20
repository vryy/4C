/*!----------------------------------------------------------------------
\file acou_impl_bdf.cpp
\brief

<pre>
Maintainers: Svenja Schoeder
             schoeder@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "acou_impl_bdf.H"
#include "acou_ele.H"
#include "acou_ele_action.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::TimIntImplBDF::TimIntImplBDF(
      const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output
      )
:AcouImplicitTimeInt(actdis,solver,params,output)
{
  order_ = 0;
  switch(dyna_)
  {
  case INPAR::ACOU::acou_bdf2:
  {
    order_ = 2;
    break;
  }
  case INPAR::ACOU::acou_bdf3:
  {
    order_ = 3;
    velnmm_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    intvelnmm_ = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    break;
  }
  case INPAR::ACOU::acou_bdf4:
  {
    order_ = 4;
    velnmm_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    intvelnmm_ = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    velnmmm_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    intvelnmmm_ = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    break;
  }
  default:
    dserror("Unknown time integration scheme for acoustical BDF time integrator");
    break;
  }
} // TimIntImplBDF

/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnp");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(intvelnp_,"intvelnp");
  reader.ReadVector(intveln_ ,"intveln");
  reader.ReadVector(intvelnm_ ,"intvelnm");
  if ( order_ > 2 )
  {
    reader.ReadVector(intvelnmm_ ,"intvelnmm");
    reader.ReadVector(velnmm_,"velnmm");
  }
  if ( order_ > 3 )
  {
    reader.ReadVector(intvelnmmm_ ,"intvelnmmm");
    reader.ReadVector(velnmmm_,"velnmmm");
  }
  return;
} // ReadRestart

/*----------------------------------------------------------------------*
 |  Initialization of algorithm to zero (public)         schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::SetInitialZeroField()
{
  // call base class function first
  ACOU::AcouImplicitTimeInt::SetInitialZeroField();

  // and additionally set BDF specific vectors to zero
  if ( order_ > 2 )
  {
    velnmm_->PutScalar(0.0);
    intvelnmm_->PutScalar(0.0);
  }
  if ( order_ > 3 )
  {
    velnmmm_->PutScalar(0.0);
    intvelnmmm_->PutScalar(0.0);
  }
  return;
} // SetInitialZeroField

/*----------------------------------------------------------------------*
 |  Initialization of algorithm by given function (pub)  schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::SetInitialField(int startfuncno, double pulse)
{
  // call base class function first
  ACOU::AcouImplicitTimeInt::SetInitialField(startfuncno, pulse);

  // and additionally initialize BDF specific vectors
  if ( order_ > 2 )
  {
    velnmm_->Update(1.0,*velnp_,0.0);
    intvelnmm_->Update(1.0,*intvelnp_,0.0);
  }
  if ( order_ > 3 )
  {
    velnmmm_->Update(1.0,*velnp_,0.0);
    intvelnmmm_->Update(1.0,*intvelnp_,0.0);
  }

  return;
} // SetInitialField

/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::SetInitialPhotoAcousticField(double pulse, Teuchos::RCP<Epetra_Vector> light, Teuchos::RCP<DRT::Discretization> scatradis, bool meshconform)
{
  // call base class function first
  ACOU::AcouImplicitTimeInt::SetInitialPhotoAcousticField(pulse, light, scatradis, meshconform);

  // and additionally initialize BDF specific vectors
  if ( order_ > 2 )
  {
    velnmm_->Update(1.0,*velnp_,0.0);
    intvelnmm_->Update(1.0,*intvelnp_,0.0);
  }
  if ( order_ > 3 )
  {
    velnmmm_->Update(1.0,*velnp_,0.0);
    intvelnmmm_->Update(1.0,*intvelnp_,0.0);
  }

  return;
} // SetInitialPhotoAcousticField

/*----------------------------------------------------------------------*
 |  Update Vectors (public)                              schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::TimeUpdate()
{
  // first update BDF specific vectors
  if ( order_ > 3)
  {
    intvelnmmm_->Update(1.0,*intvelnmm_,0.0);
    velnmmm_->Update(1.0,*velnmm_,0.0);
    intvelnmm_->Update(1.0,*intvelnm_,0.0);
    velnmm_->Update(1.0,*velnm_,0.0);
  }
  else if ( order_ > 2 )
  {
    intvelnmm_->Update(1.0,*intvelnm_,0.0);
    velnmm_->Update(1.0,*velnm_,0.0);
  }

  // call base class function to update remaining vectors
  ACOU::AcouImplicitTimeInt::TimeUpdate();

  return;
} // TimeUpdate

/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::WriteRestart()
{
  // call base class function first
  ACOU::AcouImplicitTimeInt::WriteRestart();

  // and additionally write BDF specific vectors
  if ( order_ > 2 )
  {
    output_->WriteVector("intvelnmm",intvelnmm_);
    output_->WriteVector("velnmm",velnmm_);
  }
  if ( order_ > 3 )
  {
    output_->WriteVector("intvelnmmm",intvelnmmm_);
    output_->WriteVector("velnmmm",velnmmm_);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Calculate system matrix (public)                     schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::AssembleMatAndRHS()
{
  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  residual_->Scale(0.0);

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------

  // set general vector values needed by elements
  discret_->ClearState();

  discret_->SetState("trace",veln_);
  discret_->SetState("trace_m",velnm_);

  Teuchos::RCP<Epetra_Vector> hist;
  switch(order_)
  {
  case 2:
  {
    hist = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    hist->Update(4.0/3.0,*intveln_,0.0);
    hist->Update(-1.0/3.0,*intvelnm_,1.0);
    discret_->SetState(1,"intvel",hist);
    eleparams.set<double>("dt",dtp_*2.0/3.0);
    break;
  }
  case 3:
  {
    hist = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    hist->Update(18.0/11.0,*intveln_,0.0);
    hist->Update(-9.0/11.0,*intvelnm_,1.0);
    hist->Update(2.0/11.0,*intvelnmm_,1.0);
    discret_->SetState(1,"intvel",hist);
    eleparams.set<double>("dt",dtp_*6.0/11.0);
    break;
  }
  case 4:
  {
    hist = LINALG::CreateVector(*(discret_->DofRowMap(1)),true);
    hist->Update(48.0/25.0,*intveln_,0.0);
    hist->Update(-36.0/25.0,*intvelnm_,1.0);
    hist->Update(16./25.0,*intvelnmm_,1.0);
    hist->Update(-3.0/25.0,*intvelnmmm_,1.0);
    discret_->SetState(1,"intvel",hist);
    eleparams.set<double>("dt",dtp_*12.0/25.0);
    break;
  }
  default:
    dserror("unknown order for AssembleMatAndRHS in TimIntImplBDF");
    break;
  }

  // call standard loop over elements
  bool resonly = false;// !(!bool(step_-1) || !bool(step_-restart_-1));

  eleparams.set<bool>("resonly",resonly);
  eleparams.set<int>("action",ACOU::calc_systemmat_and_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);

  discret_->Evaluate(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  discret_->ClearState();

  if(!resonly || adjoint_)
  {
    // absorbing boundary conditions
    std::string condname = "Absorbing";
    std::vector<DRT::Condition*> absorbingBC;
    discret_->GetCondition(condname,absorbingBC);
    if(absorbingBC.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_abc);
      discret_->EvaluateCondition(eleparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }
  sysmat_->Complete();

  return;
} // AssembleMatAndRHS

/*----------------------------------------------------------------------*
 | Update interior field and calculate residual (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::UpdateInteriorVariablesAndAssemebleRHS()
{
  // time measurement
  dtele_ = 0.0;
  const double tcpu=Teuchos::Time::wallTime();

  Teuchos::ParameterList eleparams;

  discret_->SetState(1,"intvel",intvelnp_); // intveln_ == intvelnp_ at this point
  discret_->SetState(1,"intvelm",intvelnm_);

  // set time step size and state vectors depending on the order
  switch(order_)
  {
  case 2:
  {
    eleparams.set<double>("dt",dtp_*2.0/3.0);
    break;
  }
  case 3:
  {
    discret_->SetState(1,"intvelmm",intvelnmm_);
    eleparams.set<double>("dt",dtp_*6.0/11.0);
    break;
  }
  case 4:
  {
    discret_->SetState(1,"intvelmm",intvelnmm_);
    discret_->SetState(1,"intvelmmm",intvelnmmm_);
    eleparams.set<double>("dt",dtp_*12.0/25.0);
    break;
  }
  default:
    dserror("unknown order for AssembleMatAndRHS in TimIntImplBDF");
    break;
  }

  // set remaining parameters
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);

  Teuchos::RCP<std::vector<double> > elevals = Teuchos::rcp(new std::vector<double>(discret_->NumGlobalElements(),0.0));
  eleparams.set<Teuchos::RCP<std::vector<double> > >("elevals",elevals);

  eleparams.set<int>("action",ACOU::update_secondary_solution_and_calc_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);

  discret_->SetState("trace",velnp_);
  discret_->SetState("trace_m",veln_);

  residual_->Scale(0.0);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);
  bool resonly = true;
  eleparams.set<bool>("resonly",resonly);

  // evaluate elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  // fill in error vector if required
  if(errormaps_)
  {
    std::vector<double> localvals = *(elevals.get());
    for(int el=0; el<discret_->NumMyRowElements(); ++el)
      error_->ReplaceGlobalValue(el,0,localvals[error_->Map().GID(el)]);
  }

  // update internal field
  const Epetra_Vector& intvelnpGhosted = *discret_->GetState(1,"intvel");
  for (int i=0; i<intvelnp_->MyLength(); ++i)
    (*intvelnp_)[i] = intvelnpGhosted[intvelnpGhosted.Map().LID(intvelnp_->Map().GID(i))];

  discret_->ClearState();

  // calculate boundary source term for inverse adjoint runs
  if(adjoint_)
  {
    // absorbing boundary conditions
    std::string condname = "Absorbing";
    std::vector<DRT::Condition*> absorbingBC;
    discret_->GetCondition(condname,absorbingBC);
    if(absorbingBC.size())
    {
      eleparams.remove("action",false);
      eleparams.set<int>("action",ACOU::calc_abc);
      discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condname);
    }
  }

  // end time measurement for element
  dtele_=Teuchos::Time::wallTime()-tcpu;

  return;
} // UpdateInteriorVariablesAndAssemebleRHS

/*----------------------------------------------------------------------*
 | Return the name of the time integrator       (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
std::string ACOU::TimIntImplBDF::Name()
{
  std::ostringstream s;
  s<<"BDF"<<order_;
  return s.str();
} // Name

