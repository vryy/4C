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
    break;
  }
  case INPAR::ACOU::acou_bdf4:
  {
    order_ = 4;
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

  reader.ReadVector(velnp_,"velnps");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");

  Teuchos::RCP<Epetra_Vector> intvelnmm;
  Teuchos::RCP<Epetra_Vector> intvelnmmm;
  if ( order_ > 2 )
  {
    intvelnmm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
    reader.ReadVector(intvelnmm,"intvelnmm");
    discret_->SetState(1,"intvelnmm",intvelnmm);
  }
  if ( order_ > 3 )
  {
    intvelnmmm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
    reader.ReadVector(intvelnmmm,"intvelnmmm");
    discret_->SetState(1,"intvelnmmm",intvelnmmm);
  }

  Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  Teuchos::RCP<Epetra_Vector> intveln = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  Teuchos::RCP<Epetra_Vector> intvelnm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  reader.ReadVector(intvelnp,"intvelnp");
  reader.ReadVector(intveln ,"intveln");
  reader.ReadVector(intvelnm ,"intvelnm");

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",ACOU::ele_init_from_restart);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  discret_->SetState(1,"intvelnp",intvelnp);
  discret_->SetState(1,"intveln",intveln);
  discret_->SetState(1,"intvelnm",intvelnm);
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  discret_->ClearState();

  return;
} // ReadRestart

/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::TimIntImplBDF::WriteRestart()
{
  Teuchos::RCP<Epetra_Vector> intvelnmm;
  Teuchos::RCP<Epetra_Vector> intvelnmmm;

  if ( order_ > 2 )
  {
    intvelnmm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
    discret_->SetState(1,"intvelnmm",intvelnmm);
  }
  if ( order_ > 3 )
  {
    intvelnmmm = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
    discret_->SetState(1,"intvelnmmm",intvelnmmm);
  }

  // call base class function first
  ACOU::AcouImplicitTimeInt::WriteRestart();

  // and additionally write BDF specific vectors
  if ( order_ > 2 )
    output_->WriteVector("intvelnmm",intvelnmm);
  if ( order_ > 3 )
    output_->WriteVector("intvelnmmm",intvelnmmm);

  return;
}


/*----------------------------------------------------------------------*
 | Return the name of the time integrator       (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
std::string ACOU::TimIntImplBDF::Name()
{
  std::ostringstream s;
  s<<"BDF"<<order_;
  return s.str();
} // Name

