/*!----------------------------------------------------------------------
\file microstatic.cpp
\brief Static control for  microstructural problems in case of multiscale
analyses

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <Epetra_LinearProblem.h>
#include <Amesos_Klu.h>

#include "microstatic.H"

#include <vector>

#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_inpar/inpar_structure.H"


using namespace IO;

/*----------------------------------------------------------------------*
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |  ctor (public)|
 *----------------------------------------------------------------------*/
STRUMULTI::MicroStatic::MicroStatic(const int microdisnum,
                                    const double V0):
microdisnum_(microdisnum),
V0_(V0)
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  discret_ = DRT::Problem::Instance(microdisnum_)->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!discret_->Filled()) discret_->FillComplete();

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // time step size etc. need to be consistent in both input files, we
  // choose to use the ones defined in the macroscale input file
  // while other parameters (like output options, convergence checks)
  // can be used individually from the microscale input file
  const Teuchos::ParameterList& sdyn_micro  = DRT::Problem::Instance(microdisnum_)->StructuralDynamicParams();
  const Teuchos::ParameterList& sdyn_macro  = DRT::Problem::Instance()->StructuralDynamicParams();

  // we can use here the parameters of the macroscale input file
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();

  // i/o options should be read from the corresponding micro-file
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance(microdisnum_)->IOParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  solver_ = rcp (new LINALG::Solver(DRT::Problem::Instance(microdisnum_)->StructSolverParams(),
                                    discret_->Comm(),
                                    DRT::Problem::Instance()->ErrorFile()->Handle()));
  discret_->ComputeNullSpaceIfNecessary(solver_->Params());

  // -------------------------------------------------------------------
  // access dynamic / io / etc. parameters
  // -------------------------------------------------------------------
  // new time intgration implementation -> generalized alpha
  // parameters are located in a sublist
  if (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn_macro,"DYNAMICTYP") == INPAR::STR::dyna_genalpha)
  {
    const Teuchos::ParameterList& genalpha  = DRT::Problem::Instance()->StructuralDynamicParams().sublist("GENALPHA");
    beta_ = genalpha.get<double>("BETA");
    gamma_ = genalpha.get<double>("GAMMA");
    alpham_ = genalpha.get<double>("ALPHA_M");
    alphaf_ = genalpha.get<double>("ALPHA_F");
  }
  // old time integration implementation
  else if (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn_macro,"DYNAMICTYP") == INPAR::STR::dyna_gen_alfa)
  {
    beta_ = sdyn_macro.get<double>("BETA");
    gamma_ = sdyn_macro.get<double>("GAMMA");
    alpham_ = sdyn_macro.get<double>("ALPHA_M");
    alphaf_ = sdyn_macro.get<double>("ALPHA_F");
  }
  else
    dserror("multi-scale problems are only implemented for imr-like generalized alpha time integration schemes");

  INPAR::STR::PredEnum pred = Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn_micro, "PREDICT");
  pred_ = pred;

  INPAR::STR::ConvCheck convcheck = Teuchos::getIntegralValue<INPAR::STR::ConvCheck>(sdyn_micro, "CONV_CHECK");
  convcheck_ = convcheck;
  INPAR::STR::VectorNorm iternorm = Teuchos::getIntegralValue<INPAR::STR::VectorNorm>(sdyn_micro,"ITERNORM");
  iternorm_ = iternorm;

  time_ = 0.0;
  dt_ = sdyn_macro.get<double>("TIMESTEP");
  step_ = 0;
  numstep_ = sdyn_macro.get<int>("NUMSTEP");
  maxiter_ = sdyn_micro.get<int>("MAXITER");
  numiter_ = -1;

  tolres_ = sdyn_micro.get<double>("TOLRES");
  toldis_ = sdyn_micro.get<double>("TOLDISP");
  printscreen_ = false;

  restart_ = probtype.get<int>("RESTART");
  restartevry_ = sdyn_macro.get<int>("RESTARTEVRY");
  iodisp_ = Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP");
  resevrydisp_ = sdyn_micro.get<int>("RESEVRYDISP");
  INPAR::STR::StressType iostress = Teuchos::getIntegralValue<INPAR::STR::StressType>(ioflags,"STRUCT_STRESS");
  iostress_ = iostress;
  resevrystrs_ = sdyn_micro.get<int>("RESEVRYSTRS");
  INPAR::STR::StrainType iostrain = Teuchos::getIntegralValue<INPAR::STR::StrainType>(ioflags,"STRUCT_STRAIN");
  iostrain_ = iostrain;
  iosurfactant_ = Teuchos::getIntegralValue<int>(ioflags,"STRUCT_SURFACTANT");

  isadapttol_ = (getIntegralValue<int>(sdyn_micro,"ADAPTCONV")==1);
  adaptolbetter_ = sdyn_micro.get<double>("ADAPTCONV_BETTER");

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  // -------------------------------------------------------------------
  if (!discret_->Filled()) discret_->FillComplete();
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  myrank_ = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // create empty matrices
  // -------------------------------------------------------------------
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,true,true));

  // -------------------------------------------------------------------
  // create empty vectors
  // -------------------------------------------------------------------
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*dofrowmap,true);
  // vector of full length; for each component
  //                /  1   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  0   i-th DOF is free
  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  // opposite of dirichtoggle vector, ie for each component
  //                /  0   i-th DOF is supported, ie Dirichlet BC
  //    vector_i =  <
  //                \  1   i-th DOF is free
  invtoggle_ = LINALG::CreateVector(*dofrowmap,false);

  // displacements D_{n} at last time
  dis_ = LINALG::CreateVector(*dofrowmap,true);

  // displacements D_{n+1} at new time
  disn_ = LINALG::CreateVector(*dofrowmap,true);

  // mid-displacements D_{n+1-alpha_f}
  dism_ = LINALG::CreateVector(*dofrowmap,true);

  // iterative displacement increments IncD_{n+1}
  // also known as residual displacements
  disi_ = LINALG::CreateVector(*dofrowmap,true);

  // internal force vector F_int at different times
  fintm_ = LINALG::CreateVector(*dofrowmap,true);

  // dynamic force residual at mid-time R_{n+1-alpha}
  // also known as out-of-balance-force
  fresm_ = LINALG::CreateVector(*dofrowmap,false);

  // -------------------------------------------------------------------
  // create "empty" EAS history map
  //
  // -------------------------------------------------------------------
  {
    lastalpha_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
    oldalpha_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
    oldfeas_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
    oldKaainv_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
    oldKda_ = Teuchos::rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
  }

  // -------------------------------------------------------------------
  // call elements to calculate stiffness and mass
  // -------------------------------------------------------------------
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time_);
    p.set("delta time",dt_);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement",zeros_);
    discret_->SetState("displacement",dis_);

    // provide EAS history of the last step
    p.set("oldalpha", oldalpha_);
    p.set("oldfeas", oldfeas_);
    p.set("oldKaainv", oldKaainv_);
    p.set("oldKda", oldKda_);

    discret_->Evaluate(p,stiff_,null,fintm_,null,null);
    discret_->ClearState();
  }

  // Determine dirichtoggle_ and its inverse since boundary conditions for
  // microscale simulations are due to the MicroBoundary condition
  // (and not Dirichlet BC)

  STRUMULTI::MicroStatic::DetermineToggle();
  STRUMULTI::MicroStatic::SetUpHomogenization();

  // reaction force vector at different times
  freactm_ = LINALG::CreateVector(*pdof_,true);

  //----------------------- compute an inverse of the dirichtoggle vector
  invtoggle_->PutScalar(1.0);
  invtoggle_->Update(-1.0,*dirichtoggle_,1.0);

  if (V0 > 0.0)
    V0_ = V0;
  else
  {
    cout << "You have not specified the initial volume of the RVE with number "
         << microdisnum_ << ", therefore it will now be calculated.\n\n"
         << "CAUTION: This calculation works only for RVEs without holes penetrating the surface!\n"
         << endl;

    // -------------------------- Calculate initial volume of microstructure
    ParameterList p;
    // action for elements
    p.set("action","calc_init_vol");
    p.set("V0", 0.0);
    discret_->EvaluateCondition(p, null, null, null, null, null, "MicroBoundary");
    V0_ = p.get<double>("V0", -1.0);
    if (V0_ == -1.0)
      dserror("Calculation of initial volume failed");
  }

  // ------------------------- Calculate initial density of microstructure
  // the macroscopic density has to be averaged over the entire
  // microstructural reference volume

  // create the parameters for the discretization
  ParameterList par;
  // action for elements
  par.set("action","calc_homog_dens");
  // set density to zero
  par.set("homogdens", 0.0);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->Evaluate(par,null,null,null,null,null);
  discret_->ClearState();

  density_ = 1/V0_*par.get<double>("homogdens", 0.0);
  if (density_ == 0.0)
    dserror("Density determined from homogenization procedure equals zero!");

  return;
} // STRUMULTI::MicroStatic::MicroStatic


void STRUMULTI::MicroStatic::Predictor(LINALG::Matrix<3,3>* defgrd)
{
  if (pred_ == INPAR::STR::pred_constdis)
    PredictConstDis(defgrd);
  else if (pred_ == INPAR::STR::pred_tangdis)
    PredictTangDis(defgrd);
  else
    dserror("requested predictor not implemented on the micro-scale");
  return;
}


/*----------------------------------------------------------------------*
 |  do predictor step (public)                               mwgee 03/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::PredictConstDis(LINALG::Matrix<3,3>* defgrd)
{
  // apply new displacements at DBCs -> this has to be done with the
  // mid-displacements since the given macroscopic deformation
  // gradient is evaluated at the mid-point!
  {
    // dism then also holds prescribed new dirichlet displacements
    EvaluateMicroBC(defgrd, dism_);
    disn_->Update(1.0, *dism_, -alphaf_, *dis_, 0.);
    disn_->Scale(1.0/(1.0-alphaf_));
    discret_->ClearState();
  }

  //--------------------------------- set EAS internal data if necessary

  // this has to be done only once since the elements will remember
  // their EAS data until the end of the microscale simulation
  // (end of macroscopic iteration step)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","eas_set_multi");

    p.set("oldalpha", oldalpha_);
    p.set("oldfeas", oldfeas_);
    p.set("oldKaainv", oldKaainv_);
    p.set("oldKda", oldKda_);

    discret_->Evaluate(p,null,null,null,null,null);
  }

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time_);
    p.set("delta time",dt_);
    p.set("alpha f",alphaf_);
    // set vector values needed by elements
    discret_->ClearState();
    disi_->Scale(0.0);
    discret_->SetState("residual displacement",disi_);
    discret_->SetState("displacement",dism_);
    fintm_->PutScalar(0.0);  // initialise internal force vector

    discret_->Evaluate(p,stiff_,null,fintm_,null,null);
    discret_->ClearState();

    if (surf_stress_man_->HaveSurfStress())
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fintm_,stiff_);
    }

    // complete stiffness matrix
    stiff_->Complete();

    // set norm of displacement increments
    disinorm_ = 1.0e6;
  }

  //-------------------------------------------- compute residual forces
  // add static mid-balance
  fresm_->Update(-1.0,*fintm_,0.0);

  // extract reaction forces
  int err = freactm_->Import(*fresm_, *importp_, Insert);
  if (err)
    dserror("Importing reaction forces of prescribed dofs using importer returned err=%d",err);

  // blank residual at DOFs on Dirichlet BC
  Epetra_Vector fresmcopy(*fresm_);
  fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

  // store norm of residual
  resnorm_ = STR::AUX::CalculateVectorNorm(iternorm_, fresm_);

  return;
} // STRUMULTI::MicroStatic::Predictor()


/*----------------------------------------------------------------------*
 |  do predictor step (public)                                  lw 01/09|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::PredictTangDis(LINALG::Matrix<3,3>* defgrd)
{
  // for displacement increments on Dirichlet boundary
  Teuchos::RCP<Epetra_Vector> dbcinc
    = LINALG::CreateVector(*(discret_->DofRowMap()), true);

  // copy last converged displacements
  dbcinc->Update(1.0, *dism_, 0.0);

  // apply new displacements at DBCs -> this has to be done with the
  // mid-displacements since the given macroscopic deformation
  // gradient is evaluated at the mid-point!
  {
    // dbcinc then also holds prescribed new dirichlet displacements
    EvaluateMicroBC(defgrd, dbcinc);
    discret_->ClearState();
  }

  // subtract the displacements of the last converged step
  // DBC-DOFs hold increments of current step
  // free-DOFs hold zeros
  dbcinc->Update(-1.0, *dism_, 1.0);

  //--------------------------------- set EAS internal data if necessary

  // this has to be done only once since the elements will remember
  // their EAS data until the end of the microscale simulation
  // (end of macroscopic iteration step)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","eas_set_multi");

    p.set("oldalpha", oldalpha_);
    p.set("oldfeas", oldfeas_);
    p.set("oldKaainv", oldKaainv_);
    p.set("oldKda", oldKda_);

    discret_->Evaluate(p,null,null,null,null,null);
  }

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time_);
    p.set("delta time",dt_);
    p.set("alpha f",alphaf_);
    // set vector values needed by elements
    discret_->ClearState();
    disi_->PutScalar(0.0);
    discret_->SetState("residual displacement",disi_);
    discret_->SetState("displacement",dism_);
    fintm_->PutScalar(0.0);  // initialise internal force vector

    discret_->Evaluate(p,stiff_,null,fintm_,null,null);
    discret_->ClearState();

    if (surf_stress_man_->HaveSurfStress())
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fintm_,stiff_);
    }
  }

  stiff_->Complete();

  //-------------------------------------------- compute residual forces
  // add static mid-balance
  fresm_->Update(-1.0,*fintm_,0.0);

  // add linear reaction forces to residual
  {
    // linear reactions
    Teuchos::RCP<Epetra_Vector> freact
      = LINALG::CreateVector(*(discret_->DofRowMap()), true);
    stiff_->Multiply(false, *dbcinc, *freact);

    // add linear reaction forces due to prescribed Dirichlet BCs
    fresm_->Update(-1.0, *freact, 1.0);
  }

  // blank residual at DOFs on Dirichlet BC
  {
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
  }

  // apply Dirichlet BCs to system of equations
  disi_->PutScalar(0.0);
  stiff_->Complete();
  LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

  // solve for disi_
  // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
  solver_->Reset();
  solver_->Solve(stiff_->EpetraMatrix(), disi_, fresm_, true, true);
  solver_->Reset();

  // store norm of displacement increments
  disinorm_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);

  //---------------------------------- update mid configuration values
  // set Dirichlet increments in displacement increments
  disi_->Update(1.0, *dbcinc, 1.0);

  // displacements
  // note that disi is not Inc_D{n+1} but Inc_D{n+1-alphaf} since everything
  // on the microscale "lives" exclusively at the pseudo generalized
  // mid point! This is just a quasi-static problem!
  dism_->Update(1.0,*disi_,1.0);
  disn_->Update(1.0/(1.0-alphaf_),*disi_,1.0);

  // reset anything that needs to be reset at the element level

  // strictly speaking, this (as well as the resetting of disi) is not
  // mandatory here, we do it just to be in line with the classical
  // time intgrator sti. there tangdis is assumed to be a predictor only, no
  // update of EAS parameters etc is desired. perhaps this might be
  // changed when speed should be optimized later on.
  {
    // create the parameters for the discretization
    ParameterList p;
    p.set("action", "calc_struct_reset_istep");
    // go to elements
    discret_->Evaluate(p, Teuchos::null, Teuchos::null,
                       Teuchos::null, Teuchos::null, Teuchos::null);
    discret_->ClearState();
  }

  //------------- eval fint at interpolated state, eval stiffness matrix
  {
    // zero out stiffness
    stiff_->Zero();
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_nlnstiff");
    // other parameters that might be needed by the elements
    p.set("total time",time_);
    p.set("delta time",dt_);
    p.set("alpha f",alphaf_);
    // set vector values needed by elements
    discret_->ClearState();
    disi_->PutScalar(0.0);
    discret_->SetState("residual displacement",disi_);
    discret_->SetState("displacement",dism_);
    fintm_->PutScalar(0.0);  // initialise internal force vector

    discret_->Evaluate(p,stiff_,null,fintm_,null,null);
    discret_->ClearState();

    if (surf_stress_man_->HaveSurfStress())
    {
      p.set("surfstr_man", surf_stress_man_);
      surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fintm_,stiff_);
    }
  }

  //-------------------------------------------- compute residual forces
  // add static mid-balance
  fresm_->Update(-1.0,*fintm_,0.0);

  // extract reaction forces
  int err = freactm_->Import(*fresm_, *importp_, Insert);
  if (err)
    dserror("Importing reaction forces of prescribed dofs using importer returned err=%d",err);

  // blank residual at DOFs on Dirichlet BC
  Epetra_Vector fresmcopy(*fresm_);
  fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

  // store norm of residual
  resnorm_ = STR::AUX::CalculateVectorNorm(iternorm_, fresm_);

  return;
}

/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             mwgee 03/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::FullNewton()
{
  //=================================================== equilibrium loop
  numiter_ = 0;

  // if TangDis-Predictor is employed, the number of iterations needs
  // to be increased by one, since it involves already one solution of
  // the non-linear system!
  if (pred_ == INPAR::STR::pred_tangdis)
    numiter_++;

  // store norms of old displacements and maximum of norms of
  // internal, external and inertial forces (needed for relative convergence
  // check)
  CalcRefNorms();

  Epetra_Time timer(discret_->Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged() && numiter_<=maxiter_)
  {

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more

    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol_ && numiter_)
    {
      double worst = resnorm_;
      double wanted = tolres_;
      solver_->AdaptTolerance(wanted,worst,adaptolbetter_);
    }
    solver_->Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,numiter_==0);
    solver_->ResetTolerance();

    //---------------------------------- update mid configuration values
    // displacements
    // note that disi is not Inc_D{n+1} but Inc_D{n+1-alphaf} since everything
    // on the microscale "lives" exclusively at the pseudo generalized
    // mid point! This is just a quasi-static problem!
    dism_->Update(1.0,*disi_,1.0);
    disn_->Update(1.0/(1.0-alphaf_),*disi_,1.0);

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // other parameters that might be needed by the elements
      p.set("total time",time_);
      p.set("delta time",dt_);
      p.set("alpha f",alphaf_);
      // set vector values needed by elements
      discret_->ClearState();
      // we do not need to scale disi_ here with 1-alphaf (cf. strugenalpha), since
      // everything on the microscale "lives" at the pseudo generalized midpoint
      // -> we solve our quasi-static problem there and only update data to the "end"
      // of the time step after having finished a macroscopic dt
      discret_->SetState("residual displacement",disi_);
      discret_->SetState("displacement",dism_);
      fintm_->PutScalar(0.0);  // initialise internal force vector

      // provide EAS history of the last step (and a place to store
      // new EAS related stuff)
      p.set("oldalpha", oldalpha_);
      p.set("oldfeas", oldfeas_);
      p.set("oldKaainv", oldKaainv_);
      p.set("oldKda", oldKda_);

      discret_->Evaluate(p,stiff_,null,fintm_,null,null);
      discret_->ClearState();

      if (surf_stress_man_->HaveSurfStress())
      {
        p.set("surfstr_man", surf_stress_man_);
        surf_stress_man_->EvaluateSurfStress(p,dism_,disn_,fintm_,stiff_);
      }
    }

    // complete stiffness matrix
    stiff_->Complete();

    //------------------------------------------ compute residual forces
    // add static mid-balance
    fresm_->Update(-1.0,*fintm_,0.0);

    // extract reaction forces
    int err = freactm_->Import(*fresm_, *importp_, Insert);
    if (err)
      dserror("Importing reaction forces of prescribed dofs using importer returned err=%d",err);

    // blank residual DOFs which are on Dirichlet BC
    Epetra_Vector fresmcopy(*fresm_);
    fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);

    //---------------------------------------------- build residual norm
    disinorm_ = STR::AUX::CalculateVectorNorm(iternorm_, disi_);

    resnorm_ = STR::AUX::CalculateVectorNorm(iternorm_, fresm_);

  //--------------------------------- increment equilibrium loop index
    ++numiter_;
  }
  //============================================= end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter_>=maxiter_)
  {
     dserror("Newton unconverged in %d iterations",numiter_);
  }

  return;
} // STRUMULTI::MicroStatic::FullNewton()


/*----------------------------------------------------------------------*
 |  write output (public)                                       lw 02/08|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::Output(RefCountPtr<DiscretizationWriter> output,
                                    const double time,
                                    const int istep,
                                    const double dt)
{
  bool isdatawritten = false;

  //------------------------------------------------- write restart step
  if (restartevry_ and step_%restartevry_==0)
  {
    output->WriteMesh(istep,time);
    output->NewStep(istep, time);
    output->WriteVector("displacement",dis_);
    isdatawritten = true;

    if (surf_stress_man_->HaveSurfStress())
      surf_stress_man_->WriteRestart(istep, time);

    RCP<std::vector<char> > lastalphadata = rcp(new std::vector<char>());

    // note that the microstructure is (currently) serial only i.e. we
    // can use the GLOBAL number of elements!
    for (int i=0;i<discret_->NumGlobalElements();++i)
    {
      RCP<Epetra_SerialDenseMatrix> lastalpha;

      if ((*lastalpha_)[i]!=null)
      {
        lastalpha = (*lastalpha_)[i];
      }
      else
      {
        lastalpha = rcp(new Epetra_SerialDenseMatrix(1, 1));
      }
      DRT::ParObject::AddtoPack(*lastalphadata, *lastalpha);
    }
    output->WriteVector("alpha", *lastalphadata, *discret_->ElementColMap());
  }

  //----------------------------------------------------- output results
  if (iodisp_ && resevrydisp_ && step_%resevrydisp_==0 && !isdatawritten)
  {
    output->NewStep(istep, time);
    output->WriteVector("displacement",dis_);
    isdatawritten = true;

    if (surf_stress_man_->HaveSurfStress() && iosurfactant_)
      surf_stress_man_->WriteResults(istep, time);
  }

  //------------------------------------- do stress calculation and output
  if (resevrystrs_ and !(istep%resevrystrs_) and iostress_!=INPAR::STR::stress_none)
  {
    // create the parameters for the discretization
    ParameterList p;
    // action for elements
    p.set("action","calc_struct_stress");
    // other parameters that might be needed by the elements
    p.set("total time",time);
    p.set("delta time",dt);
    Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
    Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
    p.set("stress", stress);
    p.set("strain", strain);
    p.set("iostress", iostress_);
    p.set("iostrain", iostrain_);
    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("residual displacement",zeros_);
    discret_->SetState("displacement",dis_);
    discret_->Evaluate(p,null,null,null,null,null);
    discret_->ClearState();
    if (!isdatawritten) output->NewStep(istep, time);
    isdatawritten = true;

    switch (iostress_)
    {
    case INPAR::STR::stress_cauchy:
      output->WriteVector("gauss_cauchy_stresses_xyz",*stress,*discret_->ElementColMap());
      break;
    case INPAR::STR::stress_2pk:
      output->WriteVector("gauss_2PK_stresses_xyz",*stress,*discret_->ElementColMap());
      break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror ("requested stress type not supported");
    }

    switch (iostrain_)
    {
    case INPAR::STR::strain_ea:
      output->WriteVector("gauss_EA_strains_xyz",*strain,*discret_->ElementColMap());
      break;
    case INPAR::STR::strain_gl:
      output->WriteVector("gauss_GL_strains_xyz",*strain,*discret_->ElementColMap());
      break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain type not supported");;
    }
  }

  return;
} // STRUMULTI::MicroStatic::Output()


/*----------------------------------------------------------------------*
 |  read restart (public)                                       lw 03/08|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::ReadRestart(int step,
                                         RCP<Epetra_Vector> dis,
                                         RCP<std::map<int, RCP<Epetra_SerialDenseMatrix> > > lastalpha,
                                         RefCountPtr<UTILS::SurfStressManager> surf_stress_man,
                                         string name)
{
  RCP<IO::InputControl> inputcontrol = rcp(new IO::InputControl(name, true));
  IO::DiscretizationReader reader(discret_, inputcontrol, step);
  double time  = reader.ReadDouble("time");
  int    rstep = reader.ReadInt("step");
  if (rstep != step) dserror("Time step on file not equal to given step");

  reader.ReadVector(dis, "displacement");
  // It does not make any sense to read the mesh and corresponding
  // element based data because we surely have different element based
  // data at every Gauss point
  // reader.ReadMesh(step);

  // Override current time and step with values from file
  time_ = time;
  step_ = rstep;

  if (surf_stress_man->HaveSurfStress())
  {
    surf_stress_man->ReadRestart(rstep, name, true);
  }

  reader.ReadSerialDenseMatrix(lastalpha, "alpha");

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 03/07|
 *----------------------------------------------------------------------*/
STRUMULTI::MicroStatic::~MicroStatic()
{
  return;
}


void STRUMULTI::MicroStatic::DetermineToggle()
{
  int np = 0;   // number of prescribed (=boundary) dofs needed for the
                // creation of vectors and matrices for homogenization
                // procedure

  vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      vector<int> dofs = discret_->Dof(actnode);
      const unsigned numdf = dofs.size();

      for (unsigned j=0; j<numdf; ++j)
      {
        const int gid = dofs[j];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);

        if ((*dirichtoggle_)[lid] != 1.0)  // be careful not to count dofs more
                                           // than once since nodes belong to
                                           // several surfaces simultaneously
          ++np;

        (*dirichtoggle_)[lid] = 1.0;
      }
    }
  }

  np_ = np;
}

void STRUMULTI::MicroStatic::EvaluateMicroBC(LINALG::Matrix<3,3>* defgrd,
                                             RefCountPtr<Epetra_Vector> disp)
{
  vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("MicroBoundary condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int j=0; j<nnode; ++j)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[j])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[j]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[j]);

      // nodal coordinates
      const double* x = actnode->X();

      // boundary displacements are prescribed via the macroscopic
      // deformation gradient
      double disp_prescribed[3];
      LINALG::Matrix<3,3> Du(defgrd->A(),false);
      LINALG::Matrix<3,3> I(true);
      I(0,0)=-1.0;
      I(1,1)=-1.0;
      I(2,2)=-1.0;
      Du+=I;

      for (int k=0; k<3;k++)
      {
        double dis = 0.;

        for (int l=0;l<3;l++)
        {
          dis += Du(k, l) * x[l];
        }

        disp_prescribed[k] = dis;
      }

      vector<int> dofs = discret_->Dof(actnode);
      //cout << "dofs:\n" << dofs[0] << "\n" << dofs[1] << "\n" << dofs[2] << endl;

      for (int l=0; l<3; ++l)
      {
        const int gid = dofs[l];

        const int lid = disp->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        (*disp)[lid] = disp_prescribed[l];
      }
    }
  }
}

void STRUMULTI::MicroStatic::SetOldState(RefCountPtr<Epetra_Vector> dis,
                                         RefCountPtr<Epetra_Vector> dism,
                                         RefCountPtr<Epetra_Vector> disn,
                                         RefCountPtr<UTILS::SurfStressManager> surfman,
                                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > lastalpha,
                                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha,
                                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldfeas,
                                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKaainv,
                                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKda)
{
  dis_ = dis;
  dism_ = dism;
  disn_ = disn;
  surf_stress_man_ = surfman;

  // using RCP's here means we do not need to return EAS data explicitly
  lastalpha_ = lastalpha;
  oldalpha_  = oldalpha;
  oldfeas_   = oldfeas;
  oldKaainv_ = oldKaainv;
  oldKda_    = oldKda;
}

void STRUMULTI::MicroStatic::UpdateNewTimeStep(RefCountPtr<Epetra_Vector> dis,
                                               RefCountPtr<Epetra_Vector> dism,
                                               RefCountPtr<Epetra_Vector> disn,
                                               RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > alpha,
                                               RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha,
                                               RefCountPtr<UTILS::SurfStressManager> surf_stress_man)
{
  // these updates hold for an imr-like generalized alpha time integration
  // -> if another time integration scheme should be used, this needs
  // to be changed accordingly

  dis->Update(1.0, *disn, 0.0);
  dism->Update(1.0, *dis, 0.0);
  disn->Update(1.0, *dis, 0.0);

  if (surf_stress_man->HaveSurfStress())
  {
    surf_stress_man->Update();
  }

  const Epetra_Map* elemap = discret_->ElementRowMap();
  for (int i=0;i<elemap->NumMyElements();++i)
  {
    RCP<Epetra_SerialDenseMatrix> alphai  = (*alpha)[i];
    RCP<Epetra_SerialDenseMatrix> alphao = (*oldalpha)[i];

    if (alphai!=null && alphao!=null) // update only those elements with EAS
    {
      Epetra_BLAS::Epetra_BLAS blas;
      blas.SCAL(alphao->M() * alphao->N(), -alphaf_/(1.0-alphaf_), alphao->A());  // alphao *= -alphaf/(1.0-alphaf)
      blas.AXPY(alphao->M() * alphao->N(), 1.0/(1.0-alphaf_), alphai->A(), alphao->A());  // alphao += 1.0/(1.0-alphaf) * alpha
      blas.COPY(alphai->M() * alphai->N(), alphao->A(), alphai->A());  // alpha := alphao
    }
  }
}

void STRUMULTI::MicroStatic::SetTime(const double timen, const double dt, const int istep)
{
  time_ = timen;
  dt_ = dt;
  step_ = istep;
}

//RefCountPtr<Epetra_Vector> STRUMULTI::MicroStatic::ReturnNewDism() { return rcp(new Epetra_Vector(*dism_)); }

void STRUMULTI::MicroStatic::ClearState()
{
  dis_ = null;
  dism_ = null;
  disn_ = null;
}

void STRUMULTI::MicroStatic::SetUpHomogenization()
{
  int indp = 0;
  int indf = 0;

  ndof_ = discret_->DofRowMap()->NumMyElements();

  std::vector <int>   pdof(np_);
  std::vector <int>   fdof(ndof_-np_);        // changed this, previously this
                                              // has been just fdof(np_),
                                              // but how should that
                                              // work for ndof_-np_>np_???

  for (int it=0; it<ndof_; ++it)
  {
    if ((*dirichtoggle_)[it] == 1.0)
    {
      pdof[indp]=discret_->DofRowMap()->GID(it);
      ++indp;
    }
    else
    {
      fdof[indf]=discret_->DofRowMap()->GID(it);
      ++indf;
    }
  }

  // create map based on the determined dofs of prescribed and free nodes
  pdof_ = rcp(new Epetra_Map(-1, np_, &pdof[0], 0, discret_->Comm()));
  fdof_ = rcp(new Epetra_Map(-1, ndof_-np_, &fdof[0], 0, discret_->Comm()));

  // create importer
  importp_ = rcp(new Epetra_Import(*pdof_, *(discret_->DofRowMap())));
  importf_ = rcp(new Epetra_Import(*fdof_, *(discret_->DofRowMap())));

  // create vector containing material coordinates of prescribed nodes
  Epetra_Vector Xp_temp(*pdof_);

  vector<DRT::Condition*> conds;
  discret_->GetCondition("MicroBoundary", conds);
  for (unsigned i=0; i<conds.size(); ++i)
  {
    const vector<int>* nodeids = conds[i]->Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("MicroBoundary condition does not have nodal cloud");
    const int nnode = (*nodeids).size();

    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret_->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret_->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);

      // nodal coordinates
      const double* x = actnode->X();

      vector<int> dofs = discret_->Dof(actnode);

      for (int k=0; k<3; ++k)
      {
        const int gid = dofs[k];

        const int lid = disn_->Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);

        for (int l=0;l<np_;++l)
        {
          if (pdof[l]==gid)
            Xp_temp[l]=x[k];
        }
      }
    }
  }

  Xp_ = LINALG::CreateVector(*pdof_,true);
  *Xp_ = Xp_temp;

  // now create D and its transpose DT (following Miehe et al., 2002)

  Epetra_Map Dmap(9, 0, Epetra_SerialComm());
  D_ = rcp(new Epetra_MultiVector(Dmap, np_));

  for (int n=0;n<np_/3;++n)
  {
    Epetra_Vector* temp1 = (*D_)(3*n);
    (*temp1)[0] = (*Xp_)[3*n];
    (*temp1)[3] = (*Xp_)[3*n+1];
    (*temp1)[6] = (*Xp_)[3*n+2];
    Epetra_Vector* temp2 = (*D_)(3*n+1);
    (*temp2)[1] = (*Xp_)[3*n+1];
    (*temp2)[4] = (*Xp_)[3*n+2];
    (*temp2)[7] = (*Xp_)[3*n];
    Epetra_Vector* temp3 = (*D_)(3*n+2);
    (*temp3)[2] = (*Xp_)[3*n+2];
    (*temp3)[5] = (*Xp_)[3*n];
    (*temp3)[8] = (*Xp_)[3*n+1];
  }

  Epetra_MultiVector DT(*pdof_, 9);

  for (int n=0;n<np_/3;++n)
  {
    (*(DT(0)))[3*n]   = (*Xp_)[3*n];
    (*(DT(1)))[3*n+1] = (*Xp_)[3*n+1];
    (*(DT(2)))[3*n+2] = (*Xp_)[3*n+2];
    (*(DT(3)))[3*n]   = (*Xp_)[3*n+1];
    (*(DT(4)))[3*n+1] = (*Xp_)[3*n+2];
    (*(DT(5)))[3*n+2] = (*Xp_)[3*n];
    (*(DT(6)))[3*n]   = (*Xp_)[3*n+2];
    (*(DT(7)))[3*n+1] = (*Xp_)[3*n];
    (*(DT(8)))[3*n+2] = (*Xp_)[3*n+1];
  }

  rhs_ = rcp(new Epetra_MultiVector(*(discret_->DofRowMap()), 9));

  for (int i=0;i<9;++i)
  {
    ((*rhs_)(i))->Export(*(DT(i)), *importp_, Insert);
  }
}


/*----------------------------------------------------------------------*
 |  check convergence of Newton iteration (public)              lw 12/07|
 *----------------------------------------------------------------------*/
bool STRUMULTI::MicroStatic::Converged()
{
  if (convcheck_ == INPAR::STR::convcheck_absres_or_absdis)
  {
    return (disinorm_ < toldis_ or resnorm_ < tolres_);
  }
  else if (convcheck_ == INPAR::STR::convcheck_absres_and_absdis)
  {
    return (disinorm_ < toldis_ and resnorm_ < tolres_);
  }
  else if (convcheck_ == INPAR::STR::convcheck_relres_or_absdis)
  {
    return (disinorm_ < toldis_ or (resnorm_/ref_resnorm_) < tolres_);
  }
  else if (convcheck_ == INPAR::STR::convcheck_relres_and_absdis)
  {
    return (disinorm_ < toldis_ and (resnorm_/ref_resnorm_) < tolres_);
  }
  else if (convcheck_ == INPAR::STR::convcheck_relres_or_reldis)
  {
    return ((disinorm_/ref_disinorm_) < toldis_  or (resnorm_/ref_resnorm_) < tolres_);
  }
  else if (convcheck_ == INPAR::STR::convcheck_relres_and_reldis)
  {
    return ((disinorm_/ref_disinorm_) < toldis_ and (resnorm_/ref_resnorm_) < tolres_);
  }
  else if (convcheck_ == INPAR::STR::convcheck_mixres_or_mixdis)
  {
    return (((disinorm_/ref_disinorm_) < toldis_ or disinorm_ < toldis_) or
            ((resnorm_/ref_resnorm_) < tolres_ or resnorm_ < tolres_));
  }
  else if (convcheck_ == INPAR::STR::convcheck_mixres_and_mixdis)
  {
    return (((disinorm_/ref_disinorm_) < toldis_ or disinorm_ < toldis_) and
            ((resnorm_/ref_resnorm_) < tolres_ or resnorm_ < tolres_));
  }
  else
  {
    dserror("Requested convergence check not (yet) implemented");
    return true;
  }
}

/*----------------------------------------------------------------------*
 |  calculate reference norms for relative convergence checks   lw 12/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::CalcRefNorms()
{
  // The reference norms are used to scale the calculated iterative
  // displacement norm and/or the residual force norm. For this
  // purpose we only need the right order of magnitude, so we don't
  // mind evaluating the corresponding norms at possibly different
  // points within the timestep (end point, generalized midpoint).
  // In the beginning (construction of macroscale time integrator and
  // first macroscale predictor), macro displacements are generally
  // 0 leading to no load on the micro problem. Consequently, the
  // microscale reference norms are 0 in case of displacements, and
  // near 0 in case of the residual (sum of ndof numerical near zero
  // values). To enable convergence in these cases, the reference norm
  // is automatically set to 1 if the calculated values are below
  // the chosen tolerances. Simply testing against 0 only works for
  // the displacements, but not for the residual!

  ref_disinorm_ = STR::AUX::CalculateVectorNorm(iternorm_, dis_);
  if (ref_disinorm_ < toldis_) ref_disinorm_ = 1.0;

  double fintnorm = STR::AUX::CalculateVectorNorm(iternorm_, fintm_);
  double freactnorm = STR::AUX::CalculateVectorNorm(iternorm_, freactm_);
  ref_resnorm_ = max(fintnorm, freactnorm);
  if (ref_resnorm_ < tolres_) ref_resnorm_ = 1.0;
}

/*----------------------------------------------------------------------*
 |  print to screen and/or error file                           lw 12/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::PrintNewton(bool print_unconv, Epetra_Time timer)
{
  bool relres        = (convcheck_ == INPAR::STR::convcheck_relres_and_absdis ||
                        convcheck_ == INPAR::STR::convcheck_relres_or_absdis);
  bool relres_reldis = (convcheck_ == INPAR::STR::convcheck_relres_and_reldis ||
                        convcheck_ == INPAR::STR::convcheck_relres_or_reldis);

  if (relres)
  {
    resnorm_ /= ref_resnorm_;
  }
  if (relres_reldis)
  {
    resnorm_ /= ref_resnorm_;
    disinorm_  /= ref_disinorm_;
  }

  if (print_unconv)
  {
    if (printscreen_)
    {
      if (relres)
      {
        printf("      MICROSCALE numiter %2d scaled res-norm %10.5e absolute dis-norm %20.15E\n",numiter_+1, resnorm_, disinorm_);
        fflush(stdout);
      }
      else if (relres_reldis)
      {
        printf("      MICROSCALE numiter %2d scaled res-norm %10.5e scaled dis-norm %20.15E\n",numiter_+1, resnorm_, disinorm_);
        fflush(stdout);
      }
      else
        {
        printf("      MICROSCALE numiter %2d absolute res-norm %10.5e absolute dis-norm %20.15E\n",numiter_+1, resnorm_, disinorm_);
        fflush(stdout);
      }
    }
  }
  else
  {
    double timepernlnsolve = timer.ElapsedTime();

    if (relres)
    {
      printf("      MICROSCALE Newton iteration converged: numiter %d scaled res-norm %e absolute dis-norm %e time %10.5f\n\n",
             numiter_,resnorm_,disinorm_,timepernlnsolve);
      fflush(stdout);
    }
    else if (relres_reldis)
    {
      printf("      MICROSCALE Newton iteration converged: numiter %d scaled res-norm %e scaled dis-norm %e time %10.5f\n\n",
             numiter_,resnorm_,disinorm_,timepernlnsolve);
      fflush(stdout);
    }
    else
    {
      printf("      MICROSCALE Newton iteration converged: numiter %d absolute res-norm %e absolute dis-norm %e time %10.5f\n\n",
             numiter_,resnorm_,disinorm_,timepernlnsolve);
      fflush(stdout);
    }
  }
}

/*----------------------------------------------------------------------*
 |  print to screen                                             lw 12/07|
 *----------------------------------------------------------------------*/
void STRUMULTI::MicroStatic::PrintPredictor()
{
  if (convcheck_ == INPAR::STR::convcheck_absres_or_absdis && convcheck_ != INPAR::STR::convcheck_absres_and_absdis)
  {
    resnorm_ /= ref_resnorm_;
    cout << "      MICROSCALE Predictor scaled res-norm " << resnorm_ << endl;
  }
  else
  {
    cout << "      MICROSCALE Predictor absolute res-norm " << resnorm_ << endl;
  }
  fflush(stdout);
}




void STRUMULTI::MicroStatic::StaticHomogenization(LINALG::Matrix<6,1>* stress,
                                                  LINALG::Matrix<6,6>* cmat,
                                                  double *density,
                                                  LINALG::Matrix<3,3>* defgrd,
                                                  const bool mod_newton,
                                                  bool& build_stiff)
{
  // determine macroscopic parameters via averaging (homogenization) of
  // microscopic features accoring to Kouznetsova, Miehe etc.
  // this was implemented against the background of serial usage
  // -> if a parallel version of microscale simulations is EVER wanted,
  // carefully check if/what/where things have to change

  // split microscale stiffness into parts corresponding to prescribed
  // and free dofs -> see thesis of Kouznetsova (Computational
  // homogenization for the multi-scale analysis of multi-phase
  // materials, Eindhoven, 2002)

  // for calculating the stresses, we need to choose the
  // right three components of freactm_ corresponding to a single node and
  // take the inner product with the material coordinates of this
  // node. The sum over all boundary nodes delivers the first
  // Piola-Kirchhoff macroscopic stress which has to be transformed
  // into the second Piola-Kirchhoff counterpart.
  // All these complicated conversions are necessary since only for
  // the energy-conjugated pair of first Piola-Kirchhoff and
  // deformation gradient the averaging integrals can be transformed
  // into integrals over the boundaries only in case of negligible
  // inertial forces (which simplifies matters significantly) whereas
  // the calling macroscopic material routine demands a second
  // Piola-Kirchhoff stress tensor.

  // IMPORTANT: the RVE has to be centered around (0,0,0), otherwise
  // modifications of this approach are necessary.

  freactm_->Scale(-1.0);

  LINALG::Matrix<3,3> P(true);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      for (int n=0; n<np_/3; ++n)
      {
        P(i,j) += (*freactm_)[n*3+i]*(*Xp_)[n*3+j];
      }
      P(i,j) /= V0_;
    }
  }

  // determine inverse of deformation gradient

  LINALG::Matrix<3,3> F_inv(defgrd->A(),false);
  F_inv.Invert();

  // convert to second Piola-Kirchhoff stresses and store them in
  // vector format
  // assembly of stresses (cf Solid3 Hex8): S11,S22,S33,S12,S23,S13

  stress->Scale(0.);

  for (int i=0; i<3; ++i)
  {
    (*stress)(0) += F_inv(0, i)*P(i,0);                     // S11
    (*stress)(1) += F_inv(1, i)*P(i,1);                     // S22
    (*stress)(2) += F_inv(2, i)*P(i,2);                     // S33
    (*stress)(3) += F_inv(0, i)*P(i,1);                     // S12
    (*stress)(4) += F_inv(1, i)*P(i,2);                     // S23
    (*stress)(5) += F_inv(0, i)*P(i,2);                     // S13
  }

  if (build_stiff)
  {
    // The calculation of the consistent macroscopic constitutive tensor
    // follows
    //
    // C. Miehe, Computational micro-to-macro transitions for
    // discretized micro-structures of heterogeneous materials at finite
    // strains based on a minimization of averaged incremental energy.
    // Computer Methods in Applied Mechanics and Engineering 192: 559-591, 2003.

    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    Epetra_MultiVector cmatpf(D_->Map(), 9);

    Epetra_Vector x(*dofrowmap);
    Epetra_Vector y(*dofrowmap);

    // make a copy
    stiff_dirich_ = Teuchos::rcp(new LINALG::SparseMatrix(*stiff_));

    stiff_->ApplyDirichlet(dirichtoggle_);

    Epetra_LinearProblem linprob(&(*stiff_->EpetraMatrix()), &x, &y);
    int error=linprob.CheckInput();
    if (error)
      dserror("Input for linear problem inconsistent");
#ifndef HAVENOT_UMFPACK
    Amesos_Umfpack solver(linprob);
    int err = solver.NumericFactorization();   // LU decomposition of stiff_ only once
    if (err)
      dserror("Numeric factorization of stiff_ for homogenization failed");

    for (int i=0;i<9;++i)
    {
      x.PutScalar(0.0);
      y.Update(1.0, *((*rhs_)(i)), 0.0);
      solver.Solve();

      Epetra_Vector f(*dofrowmap);
      stiff_dirich_->Multiply(false, x, f);
      Epetra_Vector fexp(*pdof_);
      int err = fexp.Import(f, *importp_, Insert);
      if (err)
        dserror("Export of boundary 'forces' failed with err=%d", err);

      (cmatpf(i))->Multiply('N', 'N', 1.0/V0_, *D_, fexp, 0.0);
    }
#endif

    // We now have to transform the calculated constitutive tensor
    // relating first Piola-Kirchhoff stresses to the deformation
    // gradient into a constitutive tensor relating second
    // Piola-Kirchhoff stresses to Green-Lagrange strains.

    ConvertMat(cmatpf, F_inv, *stress, *cmat);

    // after having constructed the stiffness matrix, this need not be
    // done in case of modified Newton as nonlinear solver of the
    // macroscale until the next update of macroscopic time step, when
    // build_stiff is set to true in the micromaterialgp again!

    if (mod_newton == true)
      build_stiff = false;
  }
  // homogenized density was already determined in the constructor

  *density = density_;

  return;
}


#endif
