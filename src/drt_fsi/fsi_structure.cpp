/*----------------------------------------------------------------------*/
/*!
\file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "fsi_structure.H"

#include <vector>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::Structure::Structure(Teuchos::RCP<ParameterList> params,
                          Teuchos::RCP<DRT::Discretization> dis,
                          Teuchos::RCP<LINALG::Solver> solver,
                          Teuchos::RCP<IO::DiscretizationWriter> output)

  : StruGenAlpha(*params, *dis, *solver, *output),
    params_(params),
    solver_(solver),
    output_(output)
{
  FSI::UTILS::SetupInterfaceExtractor(*dis,"FSICoupling",interface_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Structure::PrepareTimeStep()
{
  string pred = params_->get<string>("predictor","constant");
  if (pred=="constant")
  {
    ConstantPredictor();
  }
  else if (pred=="consistent")
  {
    ConsistentPredictor();
  }
  else
    dserror("predictor %s unknown", pred.c_str());

  fextncopy_ = rcp(new Epetra_Vector(*fextn_));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Structure::SetInterfaceMap(Teuchos::RCP<Epetra_Map> im)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Structure::ExtractInterfaceDispn()
{
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(dis_);
  return idis;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Structure::ExtractInterfaceDisplacement()
{
  Teuchos::RCP<Epetra_Vector> idism = interface_.ExtractCondVector(dism_);
  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(dis_);

  double alphaf = params_->get<double>("alpha f", 0.459);
  idis->Update(1./(1.-alphaf),*idism,-alphaf/(1.-alphaf));

  return idis;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Structure::PredictInterfaceDisplacement()
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  Teuchos::RCP<Epetra_Vector> idis  = interface_.ExtractCondVector(dis_);

  switch (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR"))
  {
  case 1:
    // d(n)
    // nothing to do
    break;
  case 2:
    // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
    dserror("interface velocity v(n-1) not available");
    break;
  case 3:
  {
    // d(n)+dt*v(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(vel_);

    idis->Update(dt,*ivel,1.0);
    break;
  }
  case 4:
  {
    // d(n)+dt*v(n)+0.5*dt^2*a(n)
    double dt            = params_->get<double>("delta time"             ,0.01);

    Teuchos::RCP<Epetra_Vector> ivel  = interface_.ExtractCondVector(vel_);
    Teuchos::RCP<Epetra_Vector> iacc  = interface_.ExtractCondVector(acc_);

    idis->Update(dt,*ivel,0.5*dt*dt,*iacc,1.0);
    break;
  }
  default:
    dserror("unknown interface displacement predictor '%s'",
            fsidyn.get<string>("PREDICTOR").c_str());
  }

  return idis;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Structure::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
{
  // Play it save. In the first iteration everything is already set up
  // properly. However, all following iterations need to calculate the
  // stiffness matrix here. Furthermore we are bound to reset fextm_
  // before we add our special contribution.
  // So we calculate the stiffness anyway (and waste the available
  // stiffness in the first iteration).

  // iforce gets changed. You cannot use it after this call.

  double alphaf  = params_->get<double>("alpha f", 0.459);
  double alpham  = params_->get<double>("alpha m", 0.378);

  // restort initial state
  // Do we really want that? We could just as well start from the last
  // solution. This could be much closer...
  string pred = params_->get<string>("predictor","constant");
  if (pred=="constant")
  {
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
    velm_->Update(1.0,*vel_,0.0);
    accm_->Update(1.0,*acc_,0.0);
  }
  else if (pred=="consistent")
  {
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
    velm_->Update(1.-alphaf,*veln_,alphaf,*vel_,0.0);
    accm_->Update(1.-alpham,*accn_,alpham,*acc_,0.0);
  }
  else
    dserror("predictor %s unknown", pred.c_str());

  // reset of external forces is needed
  fextn_->Update(1.0, *fextncopy_, 0.0);
  interface_.InsertCondVector(iforce,fextn_);

  fextm_->Update(1.-alphaf,*fextn_,alphaf,*fext_,0.0);

  // and rebuild the stiffness matrix
  CalculateStiffness();
}


/*----------------------------------------------------------------------*
 | element call and effective stiffness calculation          u.kue 06/07|
 *----------------------------------------------------------------------*/
void FSI::Structure::CalculateStiffness()
{
  double time    = params_->get<double>("total time",0.0);
  double dt      = params_->get<double>("delta time",0.01);
  double timen   = time + dt;  // t_{n+1}
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
  bool   damping = params_->get<bool>  ("damping"   ,false);

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",timen);
    p.set("delta time",dt);
    // set vector values needed by elements
    discret_.ClearState();
    discret_.SetState("residual displacement",disi_);
    discret_.SetState("displacement",dism_);
    //discret_.SetState("velocity",velm_); // not used at the moment
    fint_->PutScalar(0.0);  // initialise internal force vector
    discret_.Evaluate(p,stiff_,fint_);
    discret_.ClearState();
    // do NOT finalize the stiffness matrix, add mass and damping to it later
  }

  //-------------------------------------------- compute residual forces
  // Res = M . A_{n+1-alpha_m}
  //     + C . V_{n+1-alpha_f}
  //     + F_int(D_{n+1-alpha_f})
  //     - F_{ext;n+1-alpha_f}
  // add mid-inertial force
  mass_->Multiply(false,*accm_,*fresm_);
  // add mid-viscous damping force
  if (damping)
  {
      RCP<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
      damp_->Multiply(false,*velm_,*fviscm);
      fresm_->Update(1.0,*fviscm,1.0);
  }

  // add static mid-balance
  fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);

  // blank residual at DOFs on Dirichlet BC
  Epetra_Vector fresmcopy(*fresm_);
  fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::Structure::Solve()
{
  std::string equil = params_->get<string>("equilibrium iteration","full newton");
  if (equil=="full newton")
  {
    FullNewton();
  }
  else if (equil=="modified newton")
  {
    ModifiedNewton();
  }
  else if (equil=="matrixfree newton")
  {
    MatrixFreeNewton();
  }
  else if (equil=="nonlinear cg")
  {
    NonlinearCG();
  }
  else if (equil=="ptc")
  {
    PTC();
  }
  else
    dserror("Unknown type of equilibrium iteration '%s'", equil.c_str());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Structure::RelaxationSolve(Teuchos::RCP<Epetra_Vector> iforce)
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double dt        = params_->get<double>("delta time"             ,0.01);
  bool   damping   = params_->get<bool>  ("damping"                ,false);
  double beta      = params_->get<double>("beta"                   ,0.292);
  double gamma     = params_->get<double>("gamma"                  ,0.581);
  double alpham    = params_->get<double>("alpha m"                ,0.378);
  double alphaf    = params_->get<double>("alpha f"                ,0.459);

  // set external forces to just the forces at the interface
  fextn_->PutScalar(0.0);
  interface_.InsertCondVector(iforce,fextn_);

  // we start from zero
  fextm_->Update(1.-alphaf,*fextn_,0.0);

  // This (re)creates the stiffness matrix at the current
  // configuration.
  CalculateStiffness();

  //------------------------------------------- effective rhs is fresm
  //---------------------------------------------- build effective lhs
  // (using matrix stiff_ as effective matrix)
  stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
  if (damping)
    stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
  stiff_->Complete();

  //----------------------- apply dirichlet BCs to system of equations
  disi_->PutScalar(0.0);  // Useful? depends on solver and more
  LINALG::ApplyDirichlettoSystem(stiff_,disi_,fextm_,zeros_,dirichtoggle_);

  //--------------------------------------------------- solve for disi
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  solver_->Solve(stiff_->Matrix(),disi_,fextm_,true,true);

  // we are just interested in the incremental interface displacements
  Teuchos::RCP<Epetra_Vector> idisi = interface_.ExtractCondVector(disi_);

//   double norm;
//   disi_->Norm2(&norm);
//   if (disi_->Map().Comm().MyPID()==0)
//     cout << "==> disi norm = " << norm << " <==\n";

  // just to make sure...
  disi_->PutScalar(0.0);

  return idisi;
}


#endif
