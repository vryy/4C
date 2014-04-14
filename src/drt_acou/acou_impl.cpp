/*!----------------------------------------------------------------------
\file acou_impl.cpp
\brief Main control routine for acoustic simulations

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "acou_impl.H"
#include "acou_ele.H"
#include "acou_visc_ele.H"
#include "acou_ele_action.H"
#include "acou_resulttest.H"
#include "acou_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"
#include "../drt_mat/scatra_mat.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::AcouImplicitTimeInt::AcouImplicitTimeInt(
  const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
  ):
  discret_        (actdis),
  solver_         (solver),
  params_         (params),
  output_         (output),
  adjoint_        (params_->get<bool>("adjoint")),
  myrank_         (actdis->Comm().MyPID()),
  time_           (0.0),
  step_           (0),
  restart_        (params_->get<int>("restart")),
  maxtime_        (params_->get<double>("MAXTIME")),
  stepmax_        (params_->get<int>   ("NUMSTEP")),
  uprestart_      (params_->get<int>("RESTARTEVRY", -1)),
  upres_          (params_->get<int>("UPRES", -1)),
  dyna_           (DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*params_,"TIMEINT")),
  phys_           (DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*params,"PHYSICAL_TYPE")),
  numdim_         (DRT::Problem::Instance()->NDim()),
  dtp_            (params_->get<double>("TIMESTEP")),
  dtele_          (0.0),
  dtsolve_        (0.0),
  invana_         (params_->get<bool>("invana")),
  errormaps_      (DRT::INPUT::IntegralValue<bool>(*params_,"ERRORMAPS")),
  calcerr_      (DRT::INPUT::IntegralValue<INPAR::ACOU::CalcError>(*params,"CALCERROR")),
  adjoint_rhs_    (Teuchos::null)
{

  if ( dtp_ == 0.0)  dserror("Can't work with time step size == 0.0");

  // create all vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  velnp_ = LINALG::CreateVector(*dofrowmap,true);
  veln_  = LINALG::CreateVector(*dofrowmap,true);
  velnm_ = LINALG::CreateVector(*dofrowmap,true);

  // create hdg related vectors for internal fields
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  intvelnp_  = LINALG::CreateVector(*intdofrowmap,true);
  intvelnm_  = LINALG::CreateVector(*intdofrowmap,true);
  intveln_   = LINALG::CreateVector(*intdofrowmap,true);

  // create vector containing element based error values
  if(errormaps_)
    error_ = LINALG::CreateVector(*(discret_->ElementRowMap()),true);


  if(adjoint_)
    adjoint_rhs_ = params_->get<Teuchos::RCP<Epetra_MultiVector> >("rhsvec");

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*dofrowmap,true);
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise

  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    //eleparams.set("total time",time_);
    discret_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);

  }

  // print user information which might not be known by everyone
  if (errormaps_ && !myrank_ )
    std::cout<<"Local postprocessing is only effective when temporal accuracy is of order k+2, here k = "<<DRT::ELEMENTS::Acou::degree<<". Did you choose your time integrator accordingly?"<<std::endl;

  // create system matrix
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  sysmat_->Zero();

  // Vector used for solution process
  residual_      = LINALG::CreateVector(*dofrowmap,true);

  output_->WriteMesh(0,0.0);
} // AcouImplicitTimeInt

/*----------------------------------------------------------------------*
 |  Desctructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::AcouImplicitTimeInt::~AcouImplicitTimeInt()
{

}

/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(velnp_,"velnps");
  reader.ReadVector(veln_, "veln");
  reader.ReadVector(velnm_,"velnm");
  reader.ReadVector(intvelnp_,"intvelnp");
  reader.ReadVector(intveln_ ,"intveln");
  reader.ReadVector(intvelnm_ ,"intvelnm");

  return;
} // ReadRestart

/*----------------------------------------------------------------------*
 |  Initialization of algorithm to zero (public)         schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::SetInitialZeroField()
{
  velnp_->PutScalar(0.0);
  veln_->PutScalar(0.0);
  velnm_->PutScalar(0.0);
  intvelnp_->PutScalar(0.0);
  intveln_->PutScalar(0.0);
  intvelnm_->PutScalar(0.0);
} // SetInitialZeroField()

/*----------------------------------------------------------------------*
 |  Initialization of algorithm by given function (pub)  schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::SetInitialField(int startfuncno, double pulse)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;

  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::project_field);
  initParams.set<int>("startfuncno",startfuncno);
  initParams.set<double>("pulse",pulse);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

  DRT::Element::LocationArray la(2);
  int err = 0;
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    elevec1.Scale(0.0);elevec2.Scale(0.0);
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);
    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1,ele))
      elevec2.Shape(discret_->NumDof(1,ele), 1);

    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);

    // now fill the ele vector into the discretization
    for (unsigned int i=0; i<la[0].lm_.size(); ++i)
      la[0].lm_[i] = dofrowmap->LID(la[0].lm_[i]);
    err += velnp_->ReplaceMyValues(la[0].lm_.size(),elevec1.A(),&la[0].lm_[0]);
    // err += veln_ ->ReplaceMyValues(la[0].lm_.size(),elevec1.A(),&la[0].lm_[0]);
    // err += velnm_->ReplaceMyValues(la[0].lm_.size(),elevec1.A(),&la[0].lm_[0]);

    std::vector<int> localDofs = discret_->Dof(1, ele);
    dsassert(localDofs.size() == static_cast<std::size_t>(elevec2.M()), "Internal error");
    for (unsigned int i=0; i<localDofs.size(); ++i)
      localDofs[i] = intdofrowmap->LID(localDofs[i]);
    intvelnp_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
    // intveln_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
    // intvelnm_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
  }

  veln_->Update(1.0,*velnp_,0.0);
  velnm_->Update(1.0,*velnp_,0.0);
  intveln_->Update(1.0,*intvelnp_,0.0);
  intvelnm_->Update(1.0,*intvelnp_,0.0);

  return;
} // SetInitialField

/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::SetInitialPhotoAcousticField(double pulse,
                                                       Teuchos::RCP<Epetra_Vector> light,
                                                       Teuchos::RCP<DRT::Discretization> scatradis,
                                                       bool meshconform)
{
  if (meshconform == false)
    dserror("non conforming meshes not yet implemented");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  const Epetra_Map* intdofrowmap = discret_->DofRowMap(1);

  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;

  DRT::Element::LocationArray la(2);

  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::project_optical_field);
  initParams.set<double>("pulse",pulse);
  initParams.set<bool>("mesh conform",meshconform);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

  int numoptele = scatradis->NumGlobalElements();

  int minoptelegid = scatradis->ElementRowMap()->MinAllGID();

  // export light vector to column map, this is necessary
  Teuchos::RCP<Epetra_Vector> lightcol = Teuchos::rcp(new Epetra_Vector(*(scatradis->DofColMap())));
  LINALG::Export(*light,*lightcol);

  // loop all optical elements (usually there are less optical elements than acoustical elements
  for(int optel=0; optel<numoptele; ++optel)
  {
    int myopteleowner = -1;
    int opteleowner = -1;
    int myacoueleowner = -1;
    int acoueleowner = -1;

    // determine owner of the optical element
    DRT::Element* optele = NULL;
    if ( scatradis->HaveGlobalElement(optel+minoptelegid) )
    {
      optele = scatradis->gElement(optel+minoptelegid);
      myopteleowner = optele->Owner();
      if ( myopteleowner != myrank_ ) myopteleowner = -1; // do not want to consider ghosted elements
    }
    // now, every proc knows the owner of this optical element
    scatradis->Comm().MaxAll(&myopteleowner,&opteleowner,1);

    // who is the owner of the associated acoustical element?
    DRT::Element* acouele = NULL;
    if( discret_->HaveGlobalElement(optel) )
    {
      acouele = discret_->gElement(optel);
      myacoueleowner = acouele->Owner();
      if ( myacoueleowner != myrank_ ) myacoueleowner = -1; // do not want to consider ghosted elements
    }
    // now, every proc knows the owner of this acoustical element
    discret_->Comm().MaxAll(&myacoueleowner,&acoueleowner,1);

    // easiest case: both elements are on the same processor
    if ( acoueleowner == opteleowner )
    {
      // the owning proc can do all his business
      if ( opteleowner == myrank_ )
      {
        elevec1.Scale(0.0);elevec2.Scale(0.0);
        acouele->LocationVector(*discret_,la,false);
        if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
          elevec1.Shape(la[0].lm_.size(), 1);
        if (elevec2.M() != discret_->NumDof(1,acouele))
          elevec2.Shape(discret_->NumDof(1,acouele), 1);

        // we need the absorption coefficient of the light element TODO: colmap or rowmap???
        double absorptioncoeff = static_cast <MAT::ScatraMat*>((optele->Material()).get())->ReaCoeff(scatradis->ElementColMap()->LID(optele->Id()));
        initParams.set<double>("absorption",absorptioncoeff);
        Teuchos::RCP<std::vector<double> > nodevals = Teuchos::rcp(new std::vector<double>);
        int numlightnode = optele->NumNode();
        (*nodevals).resize((numdim_+1)*numlightnode);

        // fill nodevals with node coords and nodebased solution values
        DRT::Node** lightnodes = optele->Nodes();
        for(int i=0; i<numlightnode; ++i)
        {
          for(int j=0; j<numdim_; ++j)
          {
            (*nodevals)[i*(numdim_+1)+j] = lightnodes[i]->X()[j];
          }

          int dof = scatradis->Dof(lightnodes[i],0);
          int lid = lightcol->Map().LID(dof);
          if ( lid < 0 )
            dserror("given dof is not stored on proc %d although map is colmap",myrank_);
          else
            (*nodevals)[i*(numdim_+1)+numdim_] = (*(lightcol.get()))[lid];
        }

        initParams.set<Teuchos::RCP<std::vector<double> > >("nodevals",nodevals);

        // evaluate the element
        acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);

        // fill evelvec1 and elevec2 into the global vectors
        int err = 0;
        std::vector<int> localDofs = discret_->Dof(1, acouele);
        dsassert(localDofs.size() == static_cast<std::size_t>(elevec2.M()), "Internal error");
        for (unsigned int i=0; i<localDofs.size(); ++i)
          localDofs[i] = intdofrowmap->LID(localDofs[i]);
        err += intvelnp_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
        if(err) dserror("could not replace my values");

        for (unsigned int i=0; i<la[0].lm_.size(); ++i)
        {
          const int lid = dofrowmap->LID(la[0].lm_[i]);
          if ( lid >= 0 )
            (*velnp_)[lid] = elevec1(i);
        }

      }
      discret_->Comm().Barrier(); // other procs please wait for the one, who did all the work
    } // if ( acoueleowner == opteleowner )
    else // the other case: the acoustical element and the optical element are owner by different procs
    {
      // now, we perform several steps:
      // 1. get the information from the optical element we need
      // 2. communicate this information to the proc owning the acoustical element
      // 3. the acoustical element owning proc shall do his business

      // this is the 1. step:
      Teuchos::RCP<std::vector<double> > nodevals = Teuchos::rcp(new std::vector<double>);
      int size = 0;
      double absorptioncoeff = 0.0;
      if ( myrank_ == opteleowner )
      {
        // we need the absorption coefficient of the light element
        absorptioncoeff = static_cast <MAT::ScatraMat*>((optele->Material()).get())->ReaCoeff(scatradis->ElementRowMap()->LID(optele->Id()));

        int numlightnode = optele->NumNode();
        size = (numdim_+1)*numlightnode;
        (*nodevals).resize(size);

        DRT::Node** lightnodes = optele->Nodes();
        for(int i=0; i<numlightnode; ++i)
        {
          for(int j=0; j<numdim_; ++j)
          {
            (*nodevals)[i*(numdim_+1)+j] = lightnodes[i]->X()[j];
          }

          int dof = scatradis->Dof(lightnodes[i],0);
          int lid = lightcol->Map().LID(dof);
          if ( lid < 0 )
            dserror("given dof is not stored on proc %d although map is colmap",myrank_);
          else
            (*nodevals)[i*(numdim_+1)+numdim_] = (*(lightcol.get()))[lid];
        }
      }

      // this is the 2. step:
      discret_->Comm().Broadcast(&size,1,opteleowner);
      discret_->Comm().Broadcast(&absorptioncoeff,1,opteleowner);
      if ( myrank_ != opteleowner)
        (*nodevals).resize(size);
      discret_->Comm().Broadcast(&(*nodevals)[0],size,opteleowner);

      // this is the 3. step:
      if ( myrank_ == acoueleowner )
      {
        elevec1.Scale(0.0);elevec2.Scale(0.0);
        acouele->LocationVector(*discret_,la,false);
        if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
          elevec1.Shape(la[0].lm_.size(), 1);
        if (elevec2.M() != discret_->NumDof(1,acouele))
          elevec2.Shape(discret_->NumDof(1,acouele), 1);

        initParams.set<double>("absorption",absorptioncoeff);
        initParams.set<Teuchos::RCP<std::vector<double> > >("nodevals",nodevals);

        acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);

        // fill evelvec1 and elevec2 into the global vectors
        int err = 0;
        std::vector<int> localDofs = discret_->Dof(1, acouele);
        dsassert(localDofs.size() == static_cast<std::size_t>(elevec2.M()), "Internal error");
        for (unsigned int i=0; i<localDofs.size(); ++i)
          localDofs[i] = intdofrowmap->LID(localDofs[i]);
        err += intvelnp_->ReplaceMyValues(localDofs.size(), elevec2.A(), &localDofs[0]);
        if(err) dserror("could not replace my values");

        for (unsigned int i=0; i<la[0].lm_.size(); ++i)
        {
          const int lid = dofrowmap->LID(la[0].lm_[i]);
          if ( lid >= 0 )
            (*velnp_)[lid] = elevec1(i);
        }
      }

    } // else ** if ( acoueleowner == opteleowner )
  } // for(int optel=0; optel<numoptele; ++optel)

  intveln_->Update(1.0,*intvelnp_,0.0);
  intvelnm_->Update(1.0,*intvelnp_,0.0);
  veln_->Update(1.0,*velnp_,0.0);
  velnm_->Update(1.0,*velnp_,0.0);

  return;
} // SetInitialPhotoAcousticField

/*----------------------------------------------------------------------*
 |  Print some information (public)                      schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::PrintInformationToScreen()
{
  if (!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    if(!adjoint_)
      std::cout << "INTEGRATION OF AN ACOUSTIC PROBLEM USING HDG" << std::endl;
    else
      std::cout << "INTEGRATION OF AN ADJOINT ACOUSTIC PROBLEM USING HDG" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "number of nodes             " << discret_->NumGlobalNodes() << std::endl;
    std::cout << "number of elements          " << discret_->NumGlobalElements() << std::endl;
    std::cout << "number of faces             " << discret_->NumGlobalFaces() << std::endl;
    std::cout << "number of trace unknowns    " << discret_->DofRowMap(0)->NumGlobalElements() << std::endl;
    std::cout << "number of interior unknowns " << discret_->DofRowMap(1)->NumGlobalElements() << std::endl;
    if(phys_ == INPAR::ACOU::acou_lossless)
      std::cout << "polynomial order            " << DRT::ELEMENTS::Acou::degree << std::endl;
    else
      std::cout << "polynomial order            " << DRT::ELEMENTS::AcouVisc::degree << std::endl;
    std::cout << "time step size              " << dtp_ << std::endl;
    std::cout << "time integration scheme     " << Name() << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    if ( (dyna_ == INPAR::ACOU::acou_dirk23 || dyna_ == INPAR::ACOU::acou_dirk34 ) && !myrank_)
      std::cout<<"warning: you're using "<<DIRKTypeToString(dyna_)<<", this method is only a-stable, be careful"<<std::endl;
    std::cout << std::endl;
  }
  return;
} // PrintInformationToScreen

/*----------------------------------------------------------------------*
 |  Time loop (public)                                   schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::Integrate(Teuchos::RCP<Epetra_MultiVector> history, Teuchos::RCP<LINALG::MapExtractor> splitter)
{
  // time measurement: integration
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::Integrate");

  // write some information for the curious user
  PrintInformationToScreen();

  // output of initial field (given by function for purely acoustic simulation or given by optics for PAT simulation)
  Output(history,splitter);

  // call elements to calculate system matrix/rhs and assemble
  AssembleMatAndRHS();

  // apply Dirichlet boundary conditions to system of equations
  ApplyDirichletToSystem();

  // time loop
  while (step_<stepmax_ and time_<maxtime_)
  {
    // increment time and step
    IncrementTimeAndStep();

    // output to screen
    OutputToScreen();

    // solve
    Solve();

    // update solution, current solution becomes old solution of next timestep
    TimeUpdate();

    // output of solution
    Output(history,splitter);

    // evaluate error
    EvaluateErrorComparedToAnalyticalSol();
  } // while (step_<stepmax_ and time_<maxtime_)

  if (!myrank_) printf("\n");

  return;
} // Integrate

/*----------------------------------------------------------------------*
 |  Solve the system for trace and then interior field   schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::Solve()
{
  // solve linear equation and timing
  const double tcpusolve=Teuchos::Time::wallTime();

  solver_->Solve(sysmat_->EpetraOperator(),velnp_,residual_,true,false,Teuchos::null);

  // update interior variables
  UpdateInteriorVariablesAndAssemebleRHS();
  ApplyDirichletToSystem();

  dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
  return;
} // Solve

/*----------------------------------------------------------------------*
 |  Dirichlet function (public)                    schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::ApplyDirichletToSystem()
{
  TEUCHOS_FUNC_TIME_MONITOR("      + apply DBC");
  LINALG::ApplyDirichlettoSystem(sysmat_,velnp_,residual_,Teuchos::null,zeros_,*(dbcmaps_->CondMap()));
  return;
} // ApplyDirichletToSystem

/*----------------------------------------------------------------------*
 |  Calculate system matrix (public)                     schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::AssembleMatAndRHS()
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::AssembleMatAndRHS");

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // reset residual and sysmat
  residual_->Scale(0.0);
  sysmat_->Zero();

  //----------------------------------------------------------------------
  // evaluate elements
  //----------------------------------------------------------------------

  // set general vector values needed by elements
  discret_->ClearState();

  discret_->SetState("trace",velnp_);
  discret_->SetState("trace_m",veln_);

  // set history variable
  discret_->SetState(1,"intvel",intvelnp_);
  discret_->SetState(1,"intvelm",intveln_);
  eleparams.set<double>("dt",dtp_);

  // call standard loop over elements
  bool resonly = false;// !(!bool(step_-1) || !bool(step_-restart_-1));

  eleparams.set<bool>("resonly",resonly);
  eleparams.set<int>("action",ACOU::calc_systemmat_and_residual);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs",adjoint_rhs_);
  eleparams.set<int>("step",step_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

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
 |  Update Vectors (public)                              schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::TimeUpdate()
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::TimeUpdate");

  velnm_->Update(1.0,*veln_ ,0.0);
  intvelnm_->Update(1.0,*intveln_ ,0.0);
  veln_ ->Update(1.0,*velnp_,0.0);
  intveln_ ->Update(1.0,*intvelnp_,0.0);

  return;
} // TimeUpdate

/*----------------------------------------------------------------------*
 | Update interior field and calculate residual (public) schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::UpdateInteriorVariablesAndAssemebleRHS()
{
  dtele_ = 0.0;

  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::UpdateInteriorVariablesAndAssemebleRHS");

  // get cpu time
  const double tcpu=Teuchos::Time::wallTime();

  // create parameterlist
  Teuchos::ParameterList eleparams;

  // fill in parameters and set states needed by elements
  discret_->SetState(1,"intvel",intvelnp_);
  eleparams.set<double>("dt",dtp_);
  eleparams.set<bool>("adjoint",adjoint_);
  eleparams.set<bool>("errormaps",errormaps_);
  eleparams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);

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

  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,residual_,Teuchos::null,Teuchos::null);

  // update the error vector
  if(errormaps_)
  {
    std::vector<double> localvals = *(elevals.get());
    for(int el=0; el<discret_->NumMyRowElements(); ++el)
      error_->ReplaceGlobalValue(el,0,localvals[error_->Map().GID(el)]);
  }

  // update internal field for parallel usage
  const Epetra_Vector& intvelnpGhosted = *discret_->GetState(1,"intvel");
  for (int i=0; i<intvelnp_->MyLength(); ++i)
    (*intvelnp_)[i] = intvelnpGhosted[intvelnpGhosted.Map().LID(intvelnp_->Map().GID(i))];

  discret_->ClearState();

  // calculate source term for adjoint simulation
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


namespace
{
  void getNodeVectorsHDG (DRT::Discretization               &dis,
                          const Teuchos::RCP<Epetra_Vector> &interiorValues,
                          const Teuchos::RCP<Epetra_Vector> &traceValues,
                          const int                          ndim,
                          Teuchos::RCP<Epetra_MultiVector>  &velocity,
                          Teuchos::RCP<Epetra_Vector>       &pressure,
                          Teuchos::RCP<Epetra_Vector>       &tracevel,
                          Teuchos::RCP<Epetra_Vector>       &cellPres,
                          INPAR::ACOU::PhysicalType         phys)
  {
    //if (pressure.get() == NULL || pressure->MyLength() != dis.NumMyRowNodes())
    {
      const Epetra_Map* nodemap = dis.NodeRowMap();
      pressure.reset(new Epetra_Vector(*nodemap));
      velocity.reset(new Epetra_MultiVector(*nodemap,3));
      tracevel.reset(new Epetra_Vector(pressure->Map()));
      cellPres.reset(new Epetra_Vector(*dis.ElementRowMap()));
    }

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::interpolate_hdg_to_node);
    dis.SetState(1,"intvel",interiorValues);
    dis.SetState(0,"trace",traceValues);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);

    std::vector<int> dummy;
    DRT::Element::LocationArray la(2);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(pressure->MyLength());
    velocity->PutScalar(0.);
    pressure->PutScalar(0.);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      if (interpolVec.M() == 0)
        interpolVec.Resize(ele->NumNode()*(ndim+2)+1);

      ele->Evaluate(params,dis,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i=0; i<ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = pressure->Map().LID(node->Id());

        if (localIndex < 0)
          continue;

        touchCount[localIndex]++;
        for (int d=0; d<ndim; ++d)
        {
          velocity->SumIntoMyValue(localIndex,d,interpolVec(i+d*ele->NumNode()));
        }
        (*pressure)[localIndex] += interpolVec(i+ndim*ele->NumNode());
        (*tracevel)[localIndex] += interpolVec(i+(ndim+1)*ele->NumNode());
      }

      const int eleIndex = dis.ElementRowMap()->LID(ele->Id());
      if (eleIndex >= 0)
        (*cellPres)[eleIndex] += interpolVec((ndim+2)*ele->NumNode());
    }

    for (int i=0; i<pressure->MyLength(); ++i)
    {
      (*pressure)[i] /= touchCount[i];
      for (int d=0; d<ndim; ++d)
        (*velocity)[d][i] /= touchCount[i];
      (*tracevel)[i] /= touchCount[i];
    }
    dis.ClearState();

    return;
  } // getNodeVectorsHDG
  void getNodalPsiHDG (DRT::Discretization               &dis,
                       const Teuchos::RCP<Epetra_Vector> &interiorValues,
                       const Teuchos::RCP<Epetra_Vector> &traceValues,
                       Teuchos::RCP<Epetra_Vector>       &psi,
                       INPAR::ACOU::PhysicalType         phys,
                       double                            dt)
  {

    const Epetra_Map* nodemap = dis.NodeRowMap();
    psi.reset(new Epetra_Vector(*nodemap));

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::interpolate_psi_to_node);
    dis.SetState(1,"intvel",interiorValues);
    dis.SetState(0,"trace",traceValues);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);
    params.set<double>("dt",dt);

    std::vector<int> dummy;
    DRT::Element::LocationArray la(2);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(psi->MyLength());

    psi->PutScalar(0.);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      if (interpolVec.M() == 0)
        interpolVec.Resize(ele->NumNode());

      ele->Evaluate(params,dis,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i=0; i<ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = psi->Map().LID(node->Id());

        if (localIndex < 0)
          continue;

        touchCount[localIndex]++;
        (*psi)[localIndex] += interpolVec(i);
      }

    }

    for (int i=0; i<psi->MyLength(); ++i)
      (*psi)[i] /= touchCount[i];

    dis.ClearState();

    return;
  } // getNodeVectorsHDG
  void getNodeVectorsHDGVisc(DRT::Discretization              &dis,
                            const Teuchos::RCP<Epetra_Vector> &interiorValues,
                            const Teuchos::RCP<Epetra_Vector> &traceValues,
                            const int                          ndim,
                            Teuchos::RCP<Epetra_MultiVector>  &velocitygradient,
                            Teuchos::RCP<Epetra_MultiVector>  &velocity,
                            Teuchos::RCP<Epetra_Vector>       &pressure,
                            Teuchos::RCP<Epetra_Vector>       &density,
                            Teuchos::RCP<Epetra_MultiVector>  &tracevelocity,
                            Teuchos::RCP<Epetra_Vector>       &cellPres,
                            Teuchos::RCP<Epetra_Vector>       &cellDensity,
                            INPAR::ACOU::PhysicalType         phys)
  {
    {
      const Epetra_Map* nodemap = dis.NodeRowMap();
      velocity.reset(new Epetra_MultiVector(*nodemap,3));
      velocitygradient.reset(new Epetra_MultiVector(*nodemap,6));
      pressure.reset(new Epetra_Vector(*nodemap));
      density.reset(new Epetra_Vector(*nodemap));
      tracevelocity.reset(new Epetra_MultiVector(*nodemap,3));
      cellPres.reset(new Epetra_Vector(*dis.ElementRowMap()));
      cellDensity.reset(new Epetra_Vector(*dis.ElementRowMap()));
    }

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::interpolate_hdg_to_node);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);
    dis.SetState(1,"intvel",interiorValues);
    dis.SetState(0,"trace",traceValues);

    std::vector<int> dummy;
    DRT::Element::LocationArray la(2);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(pressure->MyLength());

    velocity->PutScalar(0.0);
    pressure->PutScalar(0.0);
    density->PutScalar(0.0);
    tracevelocity->PutScalar(0.0);
    cellPres->PutScalar(0.0);
    cellDensity->PutScalar(0.0);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      if (interpolVec.M() == 0)
        interpolVec.Resize(ele->NumNode()*(2*ndim+2+6)+2);

      ele->Evaluate(params,dis,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i=0; i<ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = pressure->Map().LID(node->Id());

        if (localIndex < 0)
          continue;

        touchCount[localIndex]++;
        for (int d=0; d<ndim; ++d)
        {
          velocity->SumIntoMyValue(localIndex,d,interpolVec(d*ele->NumNode()+i));
          tracevelocity->SumIntoMyValue(localIndex,d,interpolVec((d+ndim)*ele->NumNode()+i));
        }
        for (int d=0; d<6; ++d)
          velocitygradient->SumIntoMyValue(localIndex,d,interpolVec(ele->NumNode()*(2*ndim+2+d)+i+2));
        (*pressure)[localIndex] += interpolVec(ele->NumNode()*(2*ndim)+i);
        (*density)[localIndex] += interpolVec(ele->NumNode()*(2*ndim+1)+i);
      }
      const int eleIndex = dis.ElementRowMap()->LID(ele->Id());
      if (eleIndex >= 0)
      {
        (*cellPres)[eleIndex] += interpolVec(ele->NumNode()*(2*ndim+2));
        (*cellDensity)[eleIndex] += interpolVec(ele->NumNode()*(2*ndim+2)+1);
      }
    } // for (int el=0; el<dis.NumMyColElements();++el)

    for (int i=0; i<pressure->MyLength(); ++i)
    {
      (*pressure)[i] /= touchCount[i];
      (*density)[i] /= touchCount[i];
      for (int d=0; d<ndim; ++d)
      {
        (*velocity)[d][i] /= touchCount[i];
        (*tracevelocity)[d][i] /= touchCount[i];
      }
      for (int d=0; d<6; ++d)
        (*velocitygradient)[d][i] /= touchCount[i];
    }
    dis.ClearState();

    return;
  } // getNodeVectorsHDGVisc

} // namespace


/*----------------------------------------------------------------------*
 |  Output (public)                                      schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::Output(Teuchos::RCP<Epetra_MultiVector> history, Teuchos::RCP<LINALG::MapExtractor> splitter)
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::Output");

  // output of solution
  if(phys_ == INPAR::ACOU::acou_lossless)
  {
    Teuchos::RCP<Epetra_Vector> interpolatedPressure, traceVel, cellPres;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity;

    getNodeVectorsHDG(*discret_, intvelnp_, velnp_, numdim_,
                      interpolatedVelocity, interpolatedPressure, traceVel, cellPres, phys_);

    if( history != Teuchos::null )
    {
      // std::cout<<"======= History written in step "<<step_<<std::endl;
      Teuchos::RCP<Epetra_Vector> temp = (splitter->ExtractCondVector(traceVel));

      for(int i=0; i<temp->MyLength(); ++i)
        history->ReplaceMyValue(i,step_,temp->operator [](i));
    }

    if (step_%upres_ == 0)
    {
      if (myrank_ == 0 && !invana_)
        std::cout<<"======= Output written in step "<<step_<<std::endl;
      // step number and time
      output_->NewStep(step_,time_);
      // write element data only once
      if (step_==0) output_->WriteElementData(true);

      output_->WriteVector("velnp",interpolatedVelocity);
      output_->WriteVector("pressure",interpolatedPressure);
      output_->WriteVector("par_vel",traceVel);
      output_->WriteVector("pressure_avg",cellPres);
      if(errormaps_) output_->WriteVector("error",error_);

      // add restart data
      if (uprestart_ != 0 && step_%uprestart_ == 0)
      {
        WriteRestart();
      }
    }
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> interpolatedPressure, cellPres, interpolatedDensity, cellDensity;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity, traceVelocity;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocityGradient;

    getNodeVectorsHDGVisc(*discret_, intvelnp_, velnp_, numdim_,
        interpolatedVelocityGradient,interpolatedVelocity,interpolatedPressure,interpolatedDensity,
        traceVelocity,cellPres,cellDensity,phys_);

    if (step_%upres_ == 0)
    {
      if (myrank_ == 0 && !invana_)
        std::cout<<"======= Output written in step "<<step_<<std::endl;
      // step number and time
      output_->NewStep(step_,time_);

      output_->WriteVector("velnp",interpolatedVelocity);
      output_->WriteVector("pressure",interpolatedPressure);
      output_->WriteVector("density",interpolatedDensity);
      output_->WriteVector("pressure_avg",cellPres);
      output_->WriteVector("density_avg",cellDensity);
      output_->WriteVector("trace_velocity",traceVelocity);
      output_->WriteVector("velocity_gradient",interpolatedVelocityGradient,output_->nodevector);

      // add restart data
      if (uprestart_ != 0 && step_%uprestart_ == 0)
      {
        WriteRestart();
      }
    }

  }
  return;
} // Output

/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::WriteRestart()
{
  if (myrank_ == 0 && !invana_)
    std::cout<<"======= Restart written in step "<<step_<<std::endl;

  output_->WriteVector("velnps", velnp_);
  output_->WriteVector("veln", veln_);
  output_->WriteVector("velnm",velnm_);
  output_->WriteVector("intvelnp", intvelnp_);
  output_->WriteVector("intveln", intveln_);
  output_->WriteVector("intvelnm",intvelnm_);

  return;
} // WriteRestart

/*----------------------------------------------------------------------*
 |  Output time step information (public)                schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::OutputToScreen()
{
  // output to screen
  if (!myrank_)
  {
    if(invana_)
      printf(".");
    else
    printf("TIME: %11.4E/%11.4E  DT = %11.4E %s STEP = %4d/%4d, ts=%10.3E, te=%10.3E \n",time_,maxtime_,dtp_,Name().c_str(),step_,stepmax_,dtsolve_,dtele_);
  }

  return;
} // OutputToScreen

/*----------------------------------------------------------------------*
 |  Calculate node based values (public)                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec)
{
  if(phys_ == INPAR::ACOU::acou_lossless)
  {
    Teuchos::RCP<Epetra_Vector> interpolatedPressure;
    getNodalPsiHDG(*discret_, intvelnp_, velnp_,interpolatedPressure, phys_, dtp_);

    for(int i=0; i<interpolatedPressure->MyLength(); ++i)
      outvec->ReplaceMyValue(i,0,interpolatedPressure->operator [](i));
  }
  else
  {
    dserror("NodalPsiField not yet implemented for lossy fluid!");
    // TODO: implement calculation of psi for lossy fluid
  }
  return;
} // NodalPsiField

/*----------------------------------------------------------------------*
 |  Calculate node based values (public)                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec)
{
  if(phys_ == INPAR::ACOU::acou_lossless)
  {
    Teuchos::RCP<Epetra_Vector> interpolatedPressure, traceVel, cellPres;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity;

    getNodeVectorsHDG(*discret_, intvelnp_, velnp_, numdim_,
                      interpolatedVelocity, interpolatedPressure, traceVel, cellPres, phys_);

    for(int i=0; i<traceVel->MyLength(); ++i)
      outvec->ReplaceMyValue(i,0,interpolatedPressure->operator [](i));
  }
  else
  {
    Teuchos::RCP<Epetra_Vector> interpolatedPressure, cellPres, interpolatedDensity, cellDensity;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity, traceVelocity;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocityGradient;

    getNodeVectorsHDGVisc(*discret_, intvelnp_, velnp_, numdim_,
        interpolatedVelocityGradient,interpolatedVelocity,interpolatedPressure,interpolatedDensity,
        traceVelocity,cellPres,cellDensity,phys_);

    for(int i=0; i<traceVelocity->MyLength(); ++i)
      outvec->ReplaceMyValue(i,0,interpolatedPressure->operator [](i));
  }
  return;
} // NodalPressurField

/*----------------------------------------------------------------------*
 |  Return discretization (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouImplicitTimeInt::EvaluateErrorComparedToAnalyticalSol()
{
  if(calcerr_ != INPAR::ACOU::calcerror_no)
  {
    // call element routine
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::calc_acou_error);
    params.set<double>("time",time_);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
    params.set<INPAR::ACOU::CalcError>("error calculation",calcerr_);
    params.set<int>("startfuncno",params_->get<int>("STARTFUNCNO"));

    discret_->SetState(1,"intvel",intvelnp_);
    discret_->SetState(0,"trace",velnp_);

    Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(3+3));

    // call loop over elements (assemble nothing)
    discret_->EvaluateScalars(params, errors);
    discret_->ClearState();

    // std::vector containing
    // [0]: relative L2 velocity error
    // [1]: relative L2 pressure error
    // [2]: relative H1 velocity error
    Teuchos::RCP<std::vector<double> > relerror = Teuchos::rcp(new std::vector<double>(3));

    if ( (*errors)[3] != 0.0 )
      (*relerror)[1] = sqrt((*errors)[1])/sqrt((*errors)[3]);
    else if ((*errors)[1] != 0.0)
      (*relerror)[1] = 1.0;
    else
      (*relerror)[1] = 0.0;

    std::cout<<"time "<<time_<<" relative L2 pressure error "<<(*relerror)[1]<<std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Return discretization (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ACOU::AcouImplicitTimeInt::Discretization()
{
  return Teuchos::rcp_dynamic_cast<DRT::Discretization>(discret_);
} // Discretization

/*----------------------------------------------------------------------*
 |  Create test field (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ACOU::AcouImplicitTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new AcouResultTest(*this));
} // CreateFieldTest
