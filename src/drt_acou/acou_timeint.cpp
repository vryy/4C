/*!----------------------------------------------------------------------
\file acou_timeint.cpp
\brief Base class functions for implicit and explicit time integration

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "acou_timeint.H"
#include "acou_ele_action.H"
#include "acou_ele.H"
#include "acou_resulttest.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/scatra_mat.H"
#include "../linalg/linalg_utils.H"
#include "../drt_scatra/scatra_timint_stat.H"


#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 03/15 |
 *----------------------------------------------------------------------*/
ACOU::AcouTimeInt::AcouTimeInt(
  const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
  ):
  discret_        (actdis),
  solver_         (solver),
  params_         (params),
  output_         (output),
  dyna_           (DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*params_,"TIMEINT")),
  phys_           (DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*params,"PHYSICAL_TYPE")),
  invana_         (params_->get<bool>("invana")),
  padaptivity_    (DRT::INPUT::IntegralValue<bool>(*params_,"P_ADAPTIVITY")),
  adjoint_        (params_->get<bool>("adjoint")),
  myrank_         (actdis->Comm().MyPID()),
  writemonitor_   (DRT::INPUT::IntegralValue<bool>(*params_,"WRITEMONITOR")),
  writestress_    (DRT::INPUT::IntegralValue<bool>(*params_,"WRITESTRESS")),
  time_           (0.0),
  step_           (0),
  restart_        (params_->get<int>("restart")),
  maxtime_        (params_->get<double>("MAXTIME")),
  stepmax_        (params_->get<int>("NUMSTEP")),
  uprestart_      (params_->get<int>("RESTARTEVRY", -1)),
  upres_          (params_->get<int>("RESULTSEVRY", -1)),
  numdim_         (DRT::Problem::Instance()->NDim()),
  dtp_            (params_->get<double>("TIMESTEP"))
{
  // create the global trace vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  velnp_ = LINALG::CreateVector(*dofrowmap,true);

  if(invana_ && (DRT::INPUT::IntegralValue<INPAR::ACOU::InvAnalysisType>(*params_,"INV_ANALYSIS")==INPAR::ACOU::pat_reduction))
    reduction_ = true;
  else
    reduction_ = false;


  if(invana_) // has to be provided for inverse analysis (forward and adjoint)
    monitor_manager_ = params_->get<Teuchos::RCP<PATMonitorManager> >("monitormanager");

  if(!invana_)
    params_->set<bool>("timereversal",false);

  if(params_->get<std::string>("PML_DEFINITION_FILE")=="none.txt")
    withpmls_ = false;
  else
    withpmls_ = true;

} // AcouTimeInt

/*----------------------------------------------------------------------*
 |  Desctructor (public)                                 schoeder 03/15 |
 *----------------------------------------------------------------------*/
ACOU::AcouTimeInt::~AcouTimeInt()
{}

/*----------------------------------------------------------------------*
 |  Print some information (public)                      schoeder 03/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::PrintInformationToScreen()
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
    std::cout << "time step size              " << dtp_ << std::endl;
    std::cout << "time integration scheme     " << Name() << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }
  return;
} // PrintInformationToScreen

/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::SetInitialZeroField()
{
  // we have to call an init for the elements (for inverse analysis, otherwise adjoint run starts from values of previous forward run)
  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::ele_init);
  initParams.set<bool>("withPML",withpmls_);
  initParams.set<bool>("padaptivity",false);
  initParams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  // discret_->Evaluate(initParams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  Epetra_SerialDenseVector elevec;
  Epetra_SerialDenseMatrix elemat;
  DRT::Element::LocationArray la(2);
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat,elemat,elevec,elevec,elevec);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::SetInitialPhotoAcousticField(Teuchos::RCP<SCATRA::TimIntStationary> scatraalgo)
{
  // we have to call an init for the elements first!
  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::ele_init);
  initParams.set<bool>("withPML",withpmls_);
  initParams.set<bool>("padaptivity",false);
  initParams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  // discret_->Evaluate(initParams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  Epetra_SerialDenseVector elevec;
  Epetra_SerialDenseMatrix elemat1;
  Epetra_SerialDenseMatrix elemat2;
  DRT::Element::LocationArray la(2);
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec,elevec,elevec);
  }

  // some quantities we need in the following
  int minacouelegid = discret_->ElementRowMap()->MinAllGID();
  int maxacouelegid = discret_->ElementRowMap()->MaxAllGID();

  for(int i = minacouelegid; i<=maxacouelegid; ++i)
  {
    std::vector<double> gausspointsinrealcoordinates;
    int lroot = -1;
    int lsize = 0;
    if(discret_->HaveGlobalElement(i))
    {
      if(discret_->gElement(i)->Owner()==myrank_) // no ghosts please
      {
        lroot = myrank_;
        // get the element
        DRT::Element* acouele = discret_->gElement(i);
        acouele->LocationVector(*discret_,la,false);

        // set action
        initParams.set<int>("action",ACOU::get_gauss_points);

        // evaluate element
        acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec,elevec,elevec);
        lsize = elemat1.M()*elemat1.N();
        gausspointsinrealcoordinates.resize(lsize);

        for(int r=0; r<elemat1.M(); ++r) // rows is dimension
          for(int c=0; c<elemat1.N(); ++c) // columns is gausspoints
            gausspointsinrealcoordinates[c+r*elemat1.N()] = elemat1(r,c);
      }
    }

    // communitcate the gausspoints to all processors
    int gsize = 0;
    discret_->Comm().MaxAll(&lsize,&gsize,1);
    gausspointsinrealcoordinates.resize(gsize);
    int groot = -1;
    discret_->Comm().MaxAll(&lroot,&groot,1);
    discret_->Comm().Broadcast(&gausspointsinrealcoordinates[0],gsize,groot);

    // bring the gauss points in the format scatra expects
    std::vector<std::vector<double> > gausspointsforscatra(gsize/numdim_,std::vector<double>(numdim_));
    for(int gp=0; gp<gsize/numdim_; ++gp)
      for(int d=0; d<numdim_; ++d)
      {
        gausspointsforscatra[gp][d] = gausspointsinrealcoordinates[gp+d*int(gsize/numdim_)];
      }

    // evaluate phi in the gausspoints (multiplied with reaction coefficient)
    std::vector<double> values(gausspointsforscatra.size());
    scatraalgo->GetPointsPhiValues(gausspointsforscatra,values,true,0);

    // go into the acoustic element and do the L2 projection
    if(discret_->HaveGlobalElement(i))
    {
      if(discret_->gElement(i)->Owner()==myrank_) // no ghosts please
      {
        initParams.set<int>("action",ACOU::project_optical_field);
        initParams.set<double*>("gpvalues",&values[0]);
        DRT::Element* acouele = discret_->gElement(i);
        acouele->LocationVector(*discret_,la,false);
        acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec,elevec,elevec);
      }
    }
  }

  return;
} // SetInitialPhotoAcousticField

/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                       schoeder 06/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::WriteRestart()
{
  if (myrank_ == 0 && !invana_)
    std::cout<<"======= Restart written in step "<<step_<<std::endl;

  output_->WriteVector("velnps", velnp_);

  // write internal field for which we need to create and fill the corresponding vectors
  // since this requires some effort, the WriteRestart method should not be used excessively!
  Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",ACOU::fill_restart_vecs);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<int>("useacouoptvecs",-6);

  discret_->SetState(1,"intvelnp",intvelnp);
  discret_->Evaluate(eleparams);

  Teuchos::RCP<const Epetra_Vector> matrix_state = discret_->GetState(1,"intvelnp");

  intvelnp->PutScalar(0.0);
  LINALG::Export(*matrix_state,*intvelnp);
  //output_->WriteVector("intvelnp",matrix_state);
  output_->WriteVector("intvelnp",intvelnp);

  discret_->ClearState(true);
  return;
} // WriteRestart

/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  reader.ReadVector(intvelnp,"intvelnp");
  reader.ReadVector(velnp_,"velnps");

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",ACOU::ele_init_from_restart);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  discret_->SetState(1,"intvelnp",intvelnp);
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  discret_->ClearState(true);

  return;
} // ReadRestart

/*----------------------------------------------------------------------*
 |  OutputDensityAndSpeedOfSound()                       schoeder 06/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::OutputDensityAndSpeedOfSound()
{
  // build the two vectors
  Teuchos::RCP<Epetra_Vector> densvec = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap()),false));
  Teuchos::RCP<Epetra_Vector> cvec = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap()),false));

  for (int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    // map in GetParameter calculates LID, so we need GID here       05/2017 birzle
    double dens = actele->Material()->Parameter()->GetParameter(0,actele->Id());
    double c = actele->Material()->Parameter()->GetParameter(1,actele->Id());
    densvec->operator [](i) = dens;
    cvec->operator [](i) = c;
  }
  output_->WriteVector("density",densvec);
  output_->WriteVector("speedofsound",cvec);

  return;
}

namespace
{
  void getNodeVectorsHDG (DRT::Discretization               &dis,
                          const Teuchos::RCP<Epetra_Vector> &traceValues,
                          const int                          ndim,
                          Teuchos::RCP<Epetra_MultiVector>  &velocity,
                          Teuchos::RCP<Epetra_Vector>       &pressure,
                          Teuchos::RCP<Epetra_Vector>       &tracevel,
                          Teuchos::RCP<Epetra_Vector>       &cellPres,
                          INPAR::ACOU::PhysicalType         phys,
                          bool                              padapt)
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
    params.set<int>("useacouoptvecs",-1);
    dis.SetState(0,"trace",traceValues);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);
    params.set<bool>("padaptivity",padapt);

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
    dis.ClearState(true);

    return;
  } // getNodeVectorsHDG

  void getNodeVectorsHDGSolid(DRT::Discretization              &dis,
                            const Teuchos::RCP<Epetra_Vector> &traceValues,
                            const int                          ndim,
                            Teuchos::RCP<Epetra_MultiVector>  &velocitygradient,
                            Teuchos::RCP<Epetra_MultiVector>  &velocity,
                            Teuchos::RCP<Epetra_Vector>       &pressure,
                            Teuchos::RCP<Epetra_MultiVector>  &tracevelocity,
                            Teuchos::RCP<Epetra_Vector>       &cellPres,
                            INPAR::ACOU::PhysicalType         phys,
                            bool                              writestress)
  {
    {
      const Epetra_Map* nodemap = dis.NodeRowMap();
      velocity.reset(new Epetra_MultiVector(*nodemap,3));
      velocitygradient.reset(new Epetra_MultiVector(*nodemap,ndim*ndim));
      pressure.reset(new Epetra_Vector(*nodemap));
      tracevelocity.reset(new Epetra_MultiVector(*nodemap,3));
      cellPres.reset(new Epetra_Vector(*dis.ElementRowMap()));
    }

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::interpolate_hdg_to_node);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);
    params.set<bool>("writestress",writestress);
    dis.SetState(0,"trace",traceValues);

    std::vector<int> dummy;
    DRT::Element::LocationArray la(2);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(pressure->MyLength());

    velocity->PutScalar(0.0);
    pressure->PutScalar(0.0);
    tracevelocity->PutScalar(0.0);
    cellPres->PutScalar(0.0);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      if (interpolVec.M() == 0)
        interpolVec.Resize(ele->NumNode()*(2*ndim+2+ndim*ndim)+2);

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
        for (int d=0; d<ndim*ndim; ++d)
          velocitygradient->SumIntoMyValue(localIndex,d,interpolVec(ele->NumNode()*(2*ndim+2+d)+i+2));
        (*pressure)[localIndex] += interpolVec(ele->NumNode()*(2*ndim)+i);
      }
      const int eleIndex = dis.ElementRowMap()->LID(ele->Id());
      if (eleIndex >= 0)
      {
        (*cellPres)[eleIndex] += interpolVec(ele->NumNode()*(2*ndim+2));
      }
    } // for (int el=0; el<dis.NumMyColElements();++el)

    for (int i=0; i<pressure->MyLength(); ++i)
    {
      (*pressure)[i] /= touchCount[i];
      for (int d=0; d<ndim; ++d)
      {
        (*velocity)[d][i] /= touchCount[i];
        (*tracevelocity)[d][i] /= touchCount[i];
      }
      for (int d=0; d<ndim*ndim; ++d)
        (*velocitygradient)[d][i] /= touchCount[i];
    }
    dis.ClearState(true);

    return;
  } // getNodeVectorsHDGSolid

} // namespace


/*----------------------------------------------------------------------*
 |  Output (public)                                      schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::Output()
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::Output");
  if (step_%upres_ == 0)
  {
    // output of solution
    Teuchos::RCP<Epetra_Vector> interpolatedPressure, traceVel, cellPres;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity;
    Teuchos::RCP<Epetra_MultiVector> traceVelocity;
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocityGradient;
    if(phys_ == INPAR::ACOU::acou_lossless)
    {
      getNodeVectorsHDG(*discret_, velnp_, numdim_,
          interpolatedVelocity, interpolatedPressure, traceVel, cellPres, phys_,padaptivity_);
    }
    else // if(phys_ == INPAR::ACOU::acou_solid)
    {
      getNodeVectorsHDGSolid(*discret_, velnp_, numdim_,
          interpolatedVelocityGradient,interpolatedVelocity,interpolatedPressure,
          traceVelocity,cellPres,phys_,writestress_);
    }

    //fill in pressure values into monitor file, if required
    if(!invana_  || reduction_)
      FillMonitorFile(interpolatedPressure);

    Teuchos::RCP<Epetra_Vector> dmap;
    if(padaptivity_)
    {
      dmap.reset(new Epetra_Vector(*discret_->ElementRowMap()));
      for(int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        dmap->operator [](i) = double(discret_->lRowElement(i)->Degree());
      }
    }

    if (myrank_ == 0 && !invana_)
      std::cout<<"======= Output written in step "<<step_<<std::endl;
    // step number and time
    output_->NewStep(step_,time_);
    // write element data only once
    if (step_==0)
    {
      output_->WriteElementData(true);
      if(phys_!=INPAR::ACOU::acou_solid)
        OutputDensityAndSpeedOfSound();
    }

    output_->WriteVector("velnp",interpolatedVelocity);
    output_->WriteVector("pressure",interpolatedPressure);
    //output_->WriteVector("pressure_avg",cellPres);
    if(phys_ == INPAR::ACOU::acou_lossless)
    {
      //output_->WriteVector("par_vel",traceVel);
    }
    else // (phys_ == INPAR::ACOU::acou_solid)
    {
      output_->WriteVector("trace_velocity",traceVelocity);
      //output_->WriteVector("stress",interpolatedVelocityGradient,IO::nodevector);
      if(numdim_==2)
      {
        Teuchos::RCP<Epetra_Vector> stress = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap()));
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(0),0.0);
        output_->WriteVector("stress_xx",stress);
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(1),0.0);
        output_->WriteVector("stress_xy",stress);
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(2),0.0);
        output_->WriteVector("stress_yx",stress);
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(3),0.0);
        output_->WriteVector("stress_yy",stress);
      }
      if(numdim_==3)
      {
        output_->WriteVector("stress_xx",Teuchos::rcp(interpolatedVelocityGradient->operator ()(0)));
        output_->WriteVector("stress_xy",Teuchos::rcp(interpolatedVelocityGradient->operator ()(1)));
        output_->WriteVector("stress_xz",Teuchos::rcp(interpolatedVelocityGradient->operator ()(2)));
        output_->WriteVector("stress_yx",Teuchos::rcp(interpolatedVelocityGradient->operator ()(3)));
        output_->WriteVector("stress_yy",Teuchos::rcp(interpolatedVelocityGradient->operator ()(4)));
        output_->WriteVector("stress_yz",Teuchos::rcp(interpolatedVelocityGradient->operator ()(5)));
        output_->WriteVector("stress_zx",Teuchos::rcp(interpolatedVelocityGradient->operator ()(6)));
        output_->WriteVector("stress_zy",Teuchos::rcp(interpolatedVelocityGradient->operator ()(7)));
        output_->WriteVector("stress_zz",Teuchos::rcp(interpolatedVelocityGradient->operator ()(8)));
      }

    }

    //if(errormaps_) output_->WriteVector("error",error_);
    if(padaptivity_) output_->WriteVector("degree",dmap);

    // add restart data
    if (uprestart_ != 0 && step_%uprestart_ == 0)
    {
      WriteRestart();
    }
  }

  return;
} // Output

/*----------------------------------------------------------------------*
 |  Calculate node based values (public)                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec)
{
  // call element routine for interpolate HDG to elements
  Teuchos::ParameterList params;
  params.set<int>("action",ACOU::interpolate_psi_to_node);
  discret_->SetState(0,"trace",velnp_);
  params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  params.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  params.set<double>("dt",dtp_);
  params.set<bool>("padaptivity",false);

  std::vector<int> dummy;
  DRT::Element::LocationArray la(2);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;
  Epetra_SerialDenseVector interpolVec;
  std::vector<unsigned char> touchCount(outvec->MyLength());

  outvec->PutScalar(0.);

  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);
    if (interpolVec.M() == 0)
      interpolVec.Resize(ele->NumNode());

    ele->Evaluate(params,*discret_,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

    // sum values on nodes into vectors and record the touch count (build average of values)
    for (int i=0; i<ele->NumNode(); ++i)
    {
      DRT::Node* node = ele->Nodes()[i];
      const int localIndex = outvec->Map().LID(node->Id());

      if (localIndex < 0)
        continue;

      touchCount[localIndex]++;
      (*outvec)[localIndex] += interpolVec(i);
    }
  }

  for (int i=0; i<outvec->MyLength(); ++i)
    (*outvec)[i] /= touchCount[i];

  discret_->ClearState(true);

  return;
} // NodalPsiField

/*----------------------------------------------------------------------*
 |  InitMonitorFile                                      schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::InitMonitorFile()
{
  if(writemonitor_)
  {
    FILE *fp = NULL;
    if(myrank_==0)
    {
      std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
      name.append(".monitor");
      fp = fopen(name.c_str(), "w");
      if(fp == NULL)
        dserror("Couldn't open file.");
    }

    //get condition
    std::string condname="PressureMonitor";
    std::vector<DRT::Condition*>pressuremon;
    discret_->GetCondition(condname,pressuremon);

    //if(pressuremon.size()>1) dserror("write of monitor file only implemented for one pressure monitor condition");
    unsigned int last = pressuremon.size()-1;
    const std::vector<int> pressuremonmics = *(pressuremon[last]->Nodes());

    int mics=pressuremonmics.size();
    int steps=0;
    if(dtp_*stepmax_<maxtime_)
      steps=stepmax_/upres_;
    else
      steps=maxtime_/dtp_/upres_+1;//+3; // first, last and int cut off

    if(myrank_ == 0)
    {
      fprintf(fp,"steps %d ",steps);
      fprintf(fp,"mics %d\n",mics);
    }

    int speakingproc=-1;
    int helptospeak=-1;
    const double* nod_coords;
    double coords[3];

    for(unsigned int n=0;n<pressuremonmics.size();++n)
    {

      if(discret_->HaveGlobalNode(pressuremonmics[n]))
      {
        helptospeak = myrank_;
        nod_coords = discret_->gNode(pressuremonmics[n])->X();
        coords[0]=nod_coords[0];
        coords[1]=nod_coords[1];
        coords[2]=nod_coords[2];
      }
      else
        helptospeak = 0;
      discret_->Comm().MaxAll(&helptospeak,&speakingproc,1);
      discret_->Comm().Broadcast(coords,3,speakingproc);

      if(myrank_==0)
        fprintf(fp,"%e %e %e\n",coords[0],coords[1],coords[2]);
    }
    if(myrank_ == 0)
    {
      fprintf(fp,"#\n#\n#\n");
      fprintf(fp,"%e",0.0); // first time
      for(int m=0; m<mics; ++m)
        fprintf(fp," %e",0.0);
      fprintf(fp,"\n");
      fclose(fp);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  FillMonitorFile                                      schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::FillMonitorFile(Teuchos::RCP<Epetra_Vector> ip)
{
  if(writemonitor_)
  {
    FILE *fp = NULL;
    if(myrank_==0)
    {
     std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
     name.append(".monitor");
     fp = fopen(name.c_str(), "a");
    }

    // get condition
    std::string condname="PressureMonitor";
    std::vector<DRT::Condition*>pressuremon;
    discret_->GetCondition(condname,pressuremon);
    unsigned int last = pressuremon.size()-1;
    const std::vector<int> pressuremonmics = *(pressuremon[last]->Nodes());
    int mics=pressuremonmics.size();

    if(myrank_==0) fprintf(fp,"%e ",time_);
    int helptospeak=-1;
    int speakingproc=-1;
    double pressure = 0.0;

    for(int n=0;n<mics;n++)
    {
      if(discret_->HaveGlobalNode(pressuremonmics[n]))
      {
        if(ip->Map().LID(pressuremonmics[n])>=0)
        {
          helptospeak=myrank_;
          pressure=ip->operator [](ip->Map().LID(pressuremonmics[n]));
        }
      }
      else
        helptospeak = -1;
      discret_->Comm().MaxAll(&helptospeak,&speakingproc,1);
      discret_->Comm().Broadcast(&pressure,1,speakingproc);

      if(myrank_==0) fprintf(fp,"%e ", pressure);
    }
    if(myrank_==0)
    {
      fprintf(fp,"\n");
      fclose(fp);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Return discretization (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ACOU::AcouTimeInt::Discretization()
{
  return Teuchos::rcp_dynamic_cast<DRT::Discretization>(discret_);
} // Discretization


/*----------------------------------------------------------------------*
 |  Create test field (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ACOU::AcouTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new AcouResultTest(*this));
} // CreateFieldTest
