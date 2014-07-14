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
    break;
  }
  case INPAR::ACOU::acou_bdf4:
  {
    order_ = 4;
    velnmm_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
    velnmmm_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
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

  dserror("TODO: adapt read restart function for bdf");

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
    velnmm_->PutScalar(0.0);

  if ( order_ > 3 )
    velnmmm_->PutScalar(0.0);

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
    velnmm_->Update(1.0,*velnp_,0.0);

  if ( order_ > 3 )
    velnmmm_->Update(1.0,*velnp_,0.0);


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
    velnmm_->Update(1.0,*velnp_,0.0);

  if ( order_ > 3 )
    velnmmm_->Update(1.0,*velnp_,0.0);


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
    velnmmm_->Update(1.0,*velnmm_,0.0);
    velnmm_->Update(1.0,*velnm_,0.0);
  }
  else if ( order_ > 2 )
  {
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
    output_->WriteVector("velnmm",velnmm_);

  if ( order_ > 3 )
    output_->WriteVector("velnmmm",velnmmm_);

  dserror("TODO: adapt write restart function for bdf");

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

