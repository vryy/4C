/*!----------------------------------------------------------------------
\file pat_imagereconstruction.cpp


\brief image reconstruction

\level 3

<pre>
\level 3

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "pat_imagereconstruction.H"
#include "pat_utils.H"
#include "acou_expl.H"
#include "acou_impl_euler.H"
#include "acou_inv_resulttest.H"
#include "acou_ele.H"
#include "acou_ele_action.H"
#include "acou_sol_ele.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstruction::PatImageReconstruction(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout)
:
scatra_discret_(scatradis),
acou_discret_(acoudis),
scatraparams_(scatrapara),
acouparams_(acoupara),
scatrasolver_(scatrasolv),
acousolver_(acousolv),
scatraoutput_(scatraout),
acououtput_(acouout),
dyna_(DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*acouparams_,"TIMEINT")),
phys_(DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*acouparams_,"PHYSICAL_TYPE")),
name_(DRT::Problem::Instance()->OutputControlFile()->FileName()),
tstart_(Teuchos::Time::wallTime()),
tol_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("INV_TOL")),
iter_(0),
maxiter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_MAX_RUN")),
myrank_(acoudis->Comm().MyPID()),
output_count_(0),
last_acou_fw_output_count_(0),
meshconform_(DRT::INPUT::IntegralValue<bool>(*acouparams_,"MESHCONFORM")),
timereversal_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"TIMEREVERSAL")),
patchtype_(DRT::INPUT::IntegralValue<INPAR::ACOU::PatchType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"PATCHTYPE")),
fdcheck_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("FDCHECK"))),
J_(0.0),
J_start_(0.0),
error_(0.0),
error_start_(0.0),
overwrite_output_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),("OVERWRITEOUTPUT")))
{
  // set time reversal to false
  acouparams_->set<bool>("timereversal",false);
  acouparams_->set<bool>("reduction",false);

  // create necessary extra parameter list for scatra
  //{
  scatraextraparams_ = Teuchos::rcp(new Teuchos::ParameterList());
  scatraextraparams_->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());
  scatraextraparams_->set<bool>("isale",false);
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  scatraextraparams_->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
  scatraextraparams_->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
  scatraextraparams_->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
  scatraextraparams_->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");
  //}

  // initialize the needed vectors
  //{
  adjoint_psi_ = LINALG::CreateVector(*(acou_discret_->NodeRowMap()),true);
  phi_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  adjoint_phi_ = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  node_reac_ = LINALG::CreateVector(*(scatra_discret_->NodeRowMap()),true);
  //}

  // setup the line search
  linesearch_ = Teuchos::rcp(new PATLineSearch(Teuchos::rcp(this,false)));

  // setup the search direction handler
  reac_searchdirection_ = Teuchos::rcp(new PATSearchDirection(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")));
  reac_searchdirection_->Setup(scatra_discret_->ElementRowMap(),scatra_discret_->ElementRowMap());

  // create a values vector
  reac_vals_ = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  // fill values vector
  for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
  {
    DRT::Element* opti_ele = scatra_discret_->lRowElement(e);
    opti_ele->Material()->Parameter()->GetParameter(1,-1);
    reac_vals_->ReplaceMyValue(e,0,opti_ele->Material()->Parameter()->GetParameter(1,-1));
  }

  // create a gradient vector
  reac_objgrad_ = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  // read the material ids
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIPARAMLIST"));
  char* pEnd;
  while (pstream >> word2)
  {
    int id = std::strtol(&word2[0],&pEnd,10);
    if (*pEnd=='\0')
      opti_matids_.push_back(id);
  }

  // create value vector
  opti_opt_ind_ = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),true));
  for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
  {
    DRT::Element* sca_ele = scatra_discret_->lRowElement(e);
    int elematid = sca_ele->Material()->Parameter()->Id();    // check if element has material that is optimized
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      if(opti_matids_[i]==elematid)
      {
        opti_opt_ind_->ReplaceMyValue(e,0,1.0);
        break;
      }
    }
  }

  // create regularization
  if(DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE")!=INPAR::ACOU::pat_regula_none)
    reac_regula_ = Teuchos::rcp(new PATRegula(
        DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TIKHWEIGHT_MUA"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDWEIGHT_MUA"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDEPS_MUA"),
        scatra_discret_));


  // read monitor file, create multivector and map for measured values
  ReadMonitor(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("MONITORFILE"),acouparams_->get<double>("TIMESTEP"));

  // compute node based reaction vector
  ComputeNodeBasedReactionCoefficient();

  // check if we need the impulse response and if so, transform it to the used time step
  double dtimpresp = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("IMPULSERESPONSE_DT");
  if(dtimpresp!=0.0)
    conv_imp_resp_ = true;
  else
    conv_imp_resp_ = false;
  if(conv_imp_resp_)
  {
    // get file name in which the impulse response is held
    std::string impulseresponsefilename = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("IMPULSERESPONSE");

    // check if file is given
    if(impulseresponsefilename=="none.impresp")
      dserror("if you set IMPULSERESPONSE_DT != 0.0 you have to provide an impulse response in a file");

    // insert path to file if necessary
    if (impulseresponsefilename[0]!='/')
    {
      std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos!=std::string::npos)
      {
        std::string path = filename.substr(0,pos+1);
        impulseresponsefilename.insert(impulseresponsefilename.begin(), path.begin(), path.end());
      }
    }

    // open file
    FILE* file = fopen(impulseresponsefilename.c_str(),"rb");
    if (file==NULL) dserror("Could not open impulse response file %s",impulseresponsefilename.c_str());

    // read file
    char buffer[150000];
    fgets(buffer,150000,file);
    char* foundit = NULL;
    char* test = NULL;
    foundit = buffer;

    std::vector<double> impulseresponse;
    int num_imprespvals = 0;
    double norm = 0.0;
    do{
      impulseresponse.push_back(strtod(foundit,&foundit));
      test = fgets(buffer,150000,file);
      foundit = buffer;
      num_imprespvals++;
      norm += impulseresponse[num_imprespvals-1];
    } while(test!=NULL);

    double dtacou = acouparams_->get<double>("TIMESTEP");

    // in case time steps are the same just copy the impulse response
    if(dtimpresp==dtacou)
    {
      imp_resp_.Resize(num_imprespvals);
      for(int i=0; i<num_imprespvals; ++i)
        imp_resp_(i)=impulseresponse[i];
      if(norm!=0.0)
        imp_resp_.Scale(1.0/std::abs(norm));
    }
    else // otherwise interpolate the given impulse response to the dtacou_ timestep
    {
      double maxtime = dtimpresp*num_imprespvals;
      int num_baciimprespvals = maxtime/dtacou;
      imp_resp_.Resize(num_baciimprespvals);

      for(int i=0; i<num_baciimprespvals; ++i)
      {
        double actualt = i * dtacou; // we need values for this time point
        int impresindex = actualt / dtimpresp; // corresponds to this index
        imp_resp_(i) = impulseresponse[impresindex] + (impulseresponse[impresindex+1]-impulseresponse[impresindex])*(actualt - impresindex*dtimpresp)/(dtimpresp);
      }
      imp_resp_.Scale(double(num_imprespvals)/std::abs(norm)/double(num_baciimprespvals));
    }
  } // read impulse response end

  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",false);
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::ReplaceParams(Teuchos::RCP<Epetra_Vector> params)
{
  /*for(int i=0; i<params->MyLength(); ++i)
  {
    if(params->operator [](i)<0.0)
      params->ReplaceMyValue(i,0,0.0);
  }*/

  Teuchos::RCP<Epetra_Vector> paramscol = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementColMap()),false));
  LINALG::Export(*params,*paramscol);
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  reac_vals_->Update(1.0,*params,0.0);
  //int elematid = opti_ele->Material()->Parameter()->Id();
  for(unsigned int i=0; i<opti_matids_.size(); ++i)
  {
    Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
    actmat->Parameter()->SetParameter(1,paramscol);
  }

  // update node based vector
  ComputeNodeBasedReactionCoefficient();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::FDCheck()
{
  std::cout<<"FDCHECK"<<std::endl;
  if(scatra_discret_->Comm().NumProc()>1)
    dserror("FDCHECK only implemented for one processor");

  // reaction part
  {
    double J_before = J_;
    std::cout<<"reaction gradient according to adjoint analysis"<<std::endl;
    reac_objgrad_->Print(std::cout);

    Epetra_Vector fd_reac_grad(*scatra_discret_->ElementRowMap(),false);
    Teuchos::RCP<Epetra_Vector> perturb_reac_vals = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> reac_vals_before = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    reac_vals_before->Update(1.0,*reac_vals_,0.0);

    for(int i=0; i<reac_vals_->MyLength(); ++i)
    {
      double perturba = 1.0e-3;
      double perturbb = 1.0e-4;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = reac_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_reac_vals->Update(1.0,*reac_vals_before,0.0);
      perturb_reac_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_reac_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_reac_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"reaction gradient according to FD analysis"<<std::endl;
    fd_reac_grad.Print(std::cout);

    ReplaceParams(reac_vals_before);
    J_ = J_before;
  }

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstruction::EvalulateObjectiveFunction()
{
  // evaluate error contribution
  EvaluateError();
  J_ = error_;

  if(reac_regula_!=Teuchos::null)
    reac_regula_->Evaluate(reac_vals_,&J_);

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<J_-error_<<std::endl;

  return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::EvaluateGradient()
{
  // zero out gradient vector initially
  reac_objgrad_->Scale(0.0);

  // set quantities needed by the elements
  scatra_discret_->SetState("adjoint phi",adjoint_phi_);

  // fill and set psi vector
  Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  psi = CalculateAdjointOptiRhsvec(adjoint_psi_);
  scatra_discret_->SetState("psi",psi);

  // do the actual evaluation (including regularization)
  EvaluateReacGrad();

  // check gradient if required
  if(fdcheck_)
    FDCheck();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::CalculateGradDirNorm(const Epetra_Vector& bvector, const Epetra_Map& uniquemap, double* result)
{
  reac_objgrad_->Dot(bvector,result);
  return;
}


/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstruction::PerformIteration()
{
  linesearch_->Init(J_,reac_objgrad_,reac_searchdirection_->ComputeDirection(reac_objgrad_,reac_vals_,iter_),reac_vals_,scatra_discret_->ElementRowMap());
  bool reacsucc = linesearch_->Run();

  return reacsucc;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::Optimize()
{
  // initial guess with time reversal
  if(timereversal_)
    TimeReversalEstimate();

  // initial evaluation of everything
  InitialRun();

  // init
  bool success = true;

  // fitting loop
  do
  {
    if(!myrank_)
    {
      std::cout<<std::endl;
      std::cout<<"*********************************************************************************"<<std::endl;
      std::cout<<"iteration "<<iter_+1<<" of maximal "<<maxiter_<<" iterations "<<std::endl;
      std::cout<<"*********************************************************************************"<<std::endl;
      std::cout<<std::endl;
    }

    // update the sought parameters
    success = PerformIteration();

    // output some useful user information, like time consume, solution advance, ...
    OutputStats();

    // iteration count
    iter_++;

  } while ( J_>tol_ && iter_<maxiter_ && success );

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::InitialRun()
{
  // determine if we have to do the forward run or everything is zero
  double maxval = 0.0;
  reac_vals_->MaxValue(&maxval);
  if(maxval>1.0e-9)
  {
    // solve the standard problem
    SolveStandardScatra();
    SolveStandardAcou();
  }

  // calculate the error and the value of the objective function
  EvalulateObjectiveFunction();

  // set start values
  J_start_ = J_;
  error_start_ = error_;

  // solve the adjoint problem
  SolveAdjointAcou();
  SolveAdjointScatra();

  // calculate the gradient
  EvaluateGradient();

  ComputeParameterError();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::EvaluateReacGrad()
{
  // export solution vector to column map
  Teuchos::RCP<Epetra_Vector> phicol = LINALG::CreateVector(*scatra_discret_->DofColMap(),false);
  LINALG::Export(*phi_,*phicol);

  // loop elements
  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_reac);

    //initialize element vectors
    int ndof = actele->NumNode();
    Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
    Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
    Epetra_SerialDenseVector elevector1(ndof);
    Epetra_SerialDenseVector elevector2(ndof);
    Epetra_SerialDenseVector elevector3(ndof);

    DRT::Element::LocationArray la(scatra_discret_->NumDofSets());
    actele->LocationVector(*scatra_discret_,la,false);
    actele->Evaluate(p,*scatra_discret_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    //reuse elevector2
    for (int l=0; l<(int)la[0].lm_.size(); l++)
    {
      int lid = phicol->Map().LID(la[0].lm_.at(l));
      if (lid==-1) dserror("not found on this processor");
      elevector2[l] = (*phicol)[lid];
    }
    double val2 = elevector2.Dot(elevector1);
    reac_objgrad_->ReplaceMyValue(i,0,val2);

  }//loop elements

  // just to be safe
  reac_objgrad_->Multiply(1.0,*opti_opt_ind_,*reac_objgrad_,0.0);

  // evaluate the regularization gradients
  if(reac_regula_!=Teuchos::null)
    reac_regula_->EvaluateGradient(reac_vals_,reac_objgrad_);

  ConvertGradient(scatra_discret_,reac_objgrad_);

//  reac_objgrad_->Print(std::cout);
//
//    std::string soutname = name_;
//    soutname.append("_reac_grad");
//    scatraoutput_->NewResultFile(soutname,0);
//    scatraoutput_->WriteMesh(0,0.0);
//    scatraoutput_->NewStep(1,1.0);
//    scatraoutput_->WriteElementData(true);
//    scatraoutput_->WriteVector("rea_coeff",reac_objgrad_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::CheckNeighborsReacGrad(DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double reacval, double interval, Teuchos::RCP<Epetra_Vector> auxvals)
{
  // parallel version
  int lactelenodeids[4]={0,0,0,0};
  int gactelenodeids[4]={0,0,0,0};
  if(owner==myrank_)
  {
    if(actele->Shape()!=DRT::Element::quad4)
      dserror("distypes other than quad4 not yet implemented");

    for(int n=0; n<4; ++n)
      lactelenodeids[n]=actele->NodeIds()[n];
  }
  scatra_discret_->Comm().MaxAll(&lactelenodeids[0],&gactelenodeids[0],4);

  for(int n=0; n<4; ++n)
  {
    std::vector<int> toevaluate;
    if(scatra_discret_->HaveGlobalNode(gactelenodeids[n]))
    {
      DRT::Node* node = scatra_discret_->gNode(gactelenodeids[n]);
      for(int e=0; e<node->NumElement(); ++e)
      {
        DRT::Element* neighborele = node->Elements()[e];

        // is it real neighbor (only if they share 2 nodes)
        int share = 0;
        for(int a=0; a<4; ++a)
          for(int b=0; b<4; ++b)
          {
            if(gactelenodeids[a]==neighborele->NodeIds()[b])
              share++;
          }

        if(share == 4) // same element -> skip
          continue;
        else if(share == 1) // not really connected
            continue;
        else if(share == 2) // neighbor element
        {
          // if already evaluated, skip
          if(scatra_discret_->ElementRowMap()->LID(neighborele->Id())<0)
            continue;
          if(setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) <= set && setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) >= 0)
            continue;

          // determine reaction coefficient
          double neighborreac = reac_objgrad_->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id()));
          if(abs(neighborreac-reacval) <= interval)
          {
            setsids->ReplaceMyValue(scatra_discret_->ElementRowMap()->LID(neighborele->Id()),0,set);
            auxvals->ReplaceMyValue(scatra_discret_->ElementRowMap()->LID(neighborele->Id()),0,std::numeric_limits<double>::max());

            // this has to be checked and its neighbors too
            toevaluate.push_back(neighborele->Id());
          }
        }
        else
          dserror("this is strange");
      }
    }
    int lsize = toevaluate.size();
    int size = -1;
    scatra_discret_->Comm().MaxAll(&lsize,&size,1);
    if(toevaluate.size()!=unsigned(size))
      toevaluate.resize(size,0);
    std::vector<int> gtoeva(size);
    scatra_discret_->Comm().MaxAll(&toevaluate[0],&gtoeva[0],size);

    for(int s=0; s<size; ++s)
    {
      int llid = scatra_discret_->ElementRowMap()->LID(gtoeva[s]);
      int lid = -1;
      scatra_discret_->Comm().MaxAll(&llid,&lid,1);
      int lnbowner = -1;
      if(lid==llid) // owner
        lnbowner = myrank_;
      int nbowner = -1;
      scatra_discret_->Comm().MaxAll(&lnbowner,&nbowner,1);
      DRT::Element* neighborele = scatra_discret_->gElement(gtoeva[s]);
      CheckNeighborsReacGrad(neighborele,nbowner,setsids,set,reacval,interval,auxvals);
    }
  }

  /* serial version
  if(actele->Shape()==DRT::Element::quad4)
  {
    // ask each node for its elements
    for(int n=0; n<actele->NumNode(); ++n) // should be 4
    {
      DRT::Node* node = actele->Nodes()[n];
      for(int e=0; e<node->NumElement(); ++e)
      {
        DRT::Element* neighborele = node->Elements()[e];

        // is it real neighbor (only if they share 2 nodes)
        int share = 0;
        for(int a=0; a<actele->NumNode(); ++a)
          for(int b=0; b<neighborele->NumNode(); ++b)
          {
            if(actele->NodeIds()[a]==neighborele->NodeIds()[b])
              share++;
          }
        if(share == 4) // same element -> skip
          continue;
        else if(share == 2) // neighbor element
        {
          // if already evaluated, skip
          if(setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) <= set && setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) >= 0.0)
            continue;

          double neighborreac = reac_objgrad_->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id()));
          if(abs(neighborreac-reacval) <= interval)
          {
            setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) = set;
            auxvals->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) = 123456.789;
            CheckNeighbors(neighborele,owner,setsids,set,reacval,interval,auxvals);
          }
        }
        else if(share == 1) // not really connected
          continue;
        else
          dserror("to be implemented for parallel usage");
      }
    }
  }
  else
    dserror("distype not yet implemented");
  */

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::ConvertGradient(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> gradient)
{
  double L=8.0;

  // harmonic basis
  if(0)
  {
    // this is an implementation for circular harmonics

    // basis size
    int basis_size = iter_+5;

    // to calculate the gradient with respect to the values scaling the basis function, we need to scale each element gradient value with the value of the
    // basis function at that point

    // new gradient
    std::vector<double> newgrad(3*basis_size*basis_size);

    // loop basis functions to evaluate gradient contribution
    for(int i=0; i<basis_size; ++i)
    {
      for(int j=0; j<basis_size; ++j)
      {
        double loc_grad_contrib[3] = {0.0};
        for(int e=0; e<discret->NumMyRowElements(); ++e)
        {
          // get center coordinates
          std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(discret->lRowElement(e));

          // evaluate i-th shape function at these coordinates
          double N_ij_at_xy_cc = std::cos(double(i)*PI*xyz[0]/L)*std::cos(double(j)*PI*xyz[1]/L);
          double N_ij_at_xy_cs = std::cos(double(i)*PI*xyz[0]/L)*std::sin(double(j+1)*PI*xyz[1]/L);
          double N_ij_at_xy_ss = std::sin(double(i+1)*PI*xyz[0]/L)*std::sin(double(j+1)*PI*xyz[1]/L);

          // compute product of evaluated shape function and gradient value
          loc_grad_contrib[0] += N_ij_at_xy_cc * gradient->operator [](e);
          loc_grad_contrib[1] += N_ij_at_xy_cs * gradient->operator [](e);
          loc_grad_contrib[2] += N_ij_at_xy_ss * gradient->operator [](e);
        }
        double glo_grad_contrib[3] = {0.0};
        discret->Comm().SumAll(&loc_grad_contrib[0],&glo_grad_contrib[0],3);

        newgrad[i*basis_size+j] = glo_grad_contrib[0];
        newgrad[i*basis_size+j+1*basis_size*basis_size] = glo_grad_contrib[1];
        newgrad[i*basis_size+j+2*basis_size*basis_size] = glo_grad_contrib[2];

      }
    }

    // build new gradient
    for(int e=0; e<discret->NumMyRowElements(); ++e)
    {
      // get center coordinates
      std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(discret->lRowElement(e));

      double new_grad_ele = 0.0;
      for(int i=0; i<basis_size; ++i)
      {
        for(int j=0; j<basis_size; ++j)
        {
          // evaluate i-th shape function at these coordinates
          double N_ij_at_xy_cc = std::cos(double(i)*PI*xyz[0]/L)*std::cos(double(j)*PI*xyz[1]/L);
          double N_ij_at_xy_cs = std::cos(double(i)*PI*xyz[0]/L)*std::sin(double(j+1)*PI*xyz[1]/L);
          double N_ij_at_xy_ss = std::sin(double(i+1)*PI*xyz[0]/L)*std::sin(double(j+1)*PI*xyz[1]/L);

          new_grad_ele += N_ij_at_xy_cc * newgrad[i*basis_size+j]
                        + N_ij_at_xy_cs * newgrad[i*basis_size+j+1*basis_size*basis_size]
                        + N_ij_at_xy_ss * newgrad[i*basis_size+j+2*basis_size*basis_size];
        }
      }
      gradient->ReplaceMyValue(e,0,new_grad_ele);
    }
  }

  // lagrange basis
  if(0)
  {
    // basis size
    int basis_size = 5;//(iter_+2);

    // new gradient
    std::vector<double> newgrad(basis_size*basis_size);

    // loop basis functions to evaluate gradient contribution
    for(int i=0; i<basis_size; ++i)
    {
      for(int j=0; j<basis_size; ++j)
      {
        double loc_grad_contrib = 0.0;
        for(int e=0; e<discret->NumMyRowElements(); ++e)
        {
          // get center coordinates
          std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(discret->lRowElement(e));

          // evaluate i-th shape function at these coordinates
          double N_i_at_x = 1.0;
          double N_j_at_y = 1.0;

          for(int k=0; k<basis_size; ++k)
          {
            if(k!=i)
              N_i_at_x *= (xyz[0]-(-L/2.0+k*L/(basis_size-1)))/(-L/2.0+i*L/(basis_size-1)-(-L/2.0+k*L/(basis_size-1)));
            if(k!=j)
              N_j_at_y *= (xyz[1]-(-L/2.0+k*L/(basis_size-1)))/(L/2.0+j*L/(basis_size-1)-(L/2.0+k*L/(basis_size-1)));
          }

          // compute product of evaluated shape function and gradient value
          loc_grad_contrib += N_i_at_x * N_j_at_y * gradient->operator [](e);

        }
        double glo_grad_contrib = 0.0;
        discret->Comm().SumAll(&loc_grad_contrib,&glo_grad_contrib,1);

        newgrad[i*basis_size+j] = glo_grad_contrib;
      }
    }

    // build new gradient
    for(int e=0; e<discret->NumMyRowElements(); ++e)
    {
      // get center coordinates
      std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(discret->lRowElement(e));

      double new_grad_ele = 0.0;
      for(int i=0; i<basis_size; ++i)
      {
        for(int j=0; j<basis_size; ++j)
        {
          // evaluate i-th shape function at these coordinates
          double N_i_at_x = 1.0;
          double N_j_at_y = 1.0;

          for(int k=0; k<basis_size; ++k)
          {
            if(k!=i)
              N_i_at_x *= (xyz[0]-(-L/2.0+k*L/(basis_size-1)))/(-L/2.0+i*L/(basis_size-1)-(-L/2.0+k*L/(basis_size-1)));
            if(k!=j)
              N_j_at_y *= (xyz[1]-(-L/2.0+k*L/(basis_size-1)))/(-L/2.0+j*L/(basis_size-1)-(-L/2.0+k*L/(basis_size-1)));
          }

          new_grad_ele += N_i_at_x * N_j_at_y * newgrad[i*basis_size+j];
        }
      }
      gradient->ReplaceMyValue(e,0,new_grad_ele);
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionOptiSplit::PatImageReconstructionOptiSplit(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstruction(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout),
  sequenzeiter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("SEQUENZE"))
{

  // setup the search direction handler
  diff_searchdirection_ = Teuchos::rcp(new PATSearchDirection(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")));
  diff_searchdirection_->Setup(scatra_discret_->ElementRowMap(),scatra_discret_->ElementRowMap());

  // create a values vector
  diff_vals_ = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  // fill values vector
  for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
  {
    DRT::Element* opti_ele = scatra_discret_->lRowElement(e);
    opti_ele->Material()->Parameter()->GetParameter(1,-1);
    diff_vals_->ReplaceMyValue(e,0,opti_ele->Material()->Parameter()->GetParameter(0,-1));
  }

  // create a gradient vector
  diff_objgrad_ = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  // create regularization
  if(DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE")!=INPAR::ACOU::pat_regula_none)
    diff_regula_ = Teuchos::rcp(new PATRegula(
        DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TIKHWEIGHT_D"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDWEIGHT_D"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDEPS_D"),
        scatra_discret_));

  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",false);
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::ReplaceParams(Teuchos::RCP<Epetra_Vector> params)
{
  // check for negative entries
  if(reacordifforcorrho_==0)
    for(int i=0; i<params->MyLength(); ++i)
    {
      ;
      //if(params->operator [](i)<0.0)
        //params->ReplaceMyValue(i,0,0.0);
    }
  else
    for(int i=0; i<params->MyLength(); ++i)
    {
      if(params->operator [](i)<0.1)
        params->ReplaceMyValue(i,0,0.1);
      else if(params->operator [](i)>1.0)
        params->ReplaceMyValue(i,0,1.0);
    }

  Teuchos::RCP<Epetra_Vector> paramscol = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementColMap()),false));
  LINALG::Export(*params,*paramscol);
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if(reacordifforcorrho_==0)
  {
    reac_vals_->Update(1.0,*params,0.0);
    //int elematid = opti_ele->Material()->Parameter()->Id();
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
      actmat->Parameter()->SetParameter(1,paramscol);
    }
  }
  else if(reacordifforcorrho_==1)
  {
    diff_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
      actmat->Parameter()->SetParameter(0,paramscol);
    }
  }

  // update node based vector
  if(reacordifforcorrho_==0)
    ComputeNodeBasedReactionCoefficient();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos)
{
  reacordifforcorrho_ = 0;
  ReplaceParams(reacs);

  reacordifforcorrho_ = 1;
  ReplaceParams(diffs);

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstructionOptiSplit::EvalulateObjectiveFunction()
{
  // evaluate error contribution
  EvaluateError();
  J_ = error_;

  if(reac_regula_!=Teuchos::null)
    reac_regula_->Evaluate(reac_vals_,&J_);
  if(diff_regula_!=Teuchos::null)
    diff_regula_->Evaluate(diff_vals_,&J_);

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<J_-error_<<std::endl;

  return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::EvaluateGradient()
{
  // zero out gradient vector initially
  reac_objgrad_->Scale(0.0);
  diff_objgrad_->Scale(0.0);

  // set quantities needed by the elements
  scatra_discret_->SetState("adjoint phi",adjoint_phi_);

  // fill and set psi vector
  Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  psi = CalculateAdjointOptiRhsvec(adjoint_psi_);
  scatra_discret_->SetState("psi",psi);

  // do the actual evaluation (including regularization)
  EvaluateReacGrad();
  EvaluateDiffGrad();

  // check gradient if required
  if(fdcheck_)
    FDCheck();

  return;
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::EvaluateDiffGrad()
{
  // export solution vector to column map
  Teuchos::RCP<Epetra_Vector> phicol = LINALG::CreateVector(*scatra_discret_->DofColMap(),false);
  LINALG::Export(*phi_,*phicol);

  // loop elements
  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_diff);

    //initialize element vectors
    int ndof = actele->NumNode();
    Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
    Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
    Epetra_SerialDenseVector elevector1(ndof);
    Epetra_SerialDenseVector elevector2(ndof);
    Epetra_SerialDenseVector elevector3(ndof);

    DRT::Element::LocationArray la(scatra_discret_->NumDofSets());
    actele->LocationVector(*scatra_discret_,la,false);
    actele->Evaluate(p,*scatra_discret_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    //reuse elevector2
    for (int l=0; l<(int)la[0].lm_.size(); l++)
    {
      int lid=phicol->Map().LID(la[0].lm_.at(l));
      if (lid==-1) dserror("not found on this processor");
      elevector2[l] = (*phicol)[lid];
    }
    double val2 = elevector2.Dot(elevector1);
    diff_objgrad_->ReplaceMyValue(i,0,val2);

  }//loop elements

  // just to be safe
  diff_objgrad_->Multiply(1.0,*opti_opt_ind_,*diff_objgrad_,0.0);

  if(diff_regula_!=Teuchos::null)
    diff_regula_->EvaluateGradient(diff_vals_,diff_objgrad_);

  ConvertGradient(scatra_discret_,diff_objgrad_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::CheckNeighborsDiffGrad(DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double diffval, double interval, Teuchos::RCP<Epetra_Vector> auxvals)
{
  // parallel version
  int lactelenodeids[4]={0,0,0,0};
  int gactelenodeids[4]={0,0,0,0};
  if(owner==myrank_)
  {
    if(actele->Shape()!=DRT::Element::quad4)
      dserror("distypes other than quad4 not yet implemented");

    for(int n=0; n<4; ++n)
      lactelenodeids[n]=actele->NodeIds()[n];
  }
  scatra_discret_->Comm().MaxAll(&lactelenodeids[0],&gactelenodeids[0],4);

  for(int n=0; n<4; ++n)
  {
    std::vector<int> toevaluate;
    if(scatra_discret_->HaveGlobalNode(gactelenodeids[n]))
    {
      DRT::Node* node = scatra_discret_->gNode(gactelenodeids[n]);
      for(int e=0; e<node->NumElement(); ++e)
      {
        DRT::Element* neighborele = node->Elements()[e];

        // is it real neighbor (only if they share 2 nodes)
        int share = 0;
        for(int a=0; a<4; ++a)
          for(int b=0; b<4; ++b)
          {
            if(gactelenodeids[a]==neighborele->NodeIds()[b])
              share++;
          }

        if(share == 4) // same element -> skip
          continue;
        else if(share == 1) // not really connected
            continue;
        else if(share == 2) // neighbor element
        {
          // if already evaluated, skip
          if(scatra_discret_->ElementRowMap()->LID(neighborele->Id())<0)
            continue;
          if(setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) <= set && setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) >= 0.0)
            continue;

          // determine reaction coefficient
          double neighbordiff = diff_objgrad_->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id()));
          if(abs(neighbordiff-diffval) <= interval)
          {
            setsids->ReplaceMyValue(scatra_discret_->ElementRowMap()->LID(neighborele->Id()),0,set);
            auxvals->ReplaceMyValue(scatra_discret_->ElementRowMap()->LID(neighborele->Id()),0,std::numeric_limits<double>::max());

            // this has to be checked and its neighbors too
            toevaluate.push_back(neighborele->Id());
          }
        }
        else
          dserror("how can two quad4 elements share exactly 3 nodes??");
      }
    }
    int lsize = toevaluate.size();
    int size = -1;
    scatra_discret_->Comm().MaxAll(&lsize,&size,1);
    if(toevaluate.size()!=unsigned(size))
      toevaluate.resize(size,0);
    std::vector<int> gtoeva(size);
    scatra_discret_->Comm().MaxAll(&toevaluate[0],&gtoeva[0],size);

    for(int s=0; s<size; ++s)
    {
      int llid = scatra_discret_->ElementRowMap()->LID(gtoeva[s]);
      int lid = -1;
      scatra_discret_->Comm().MaxAll(&llid,&lid,1);
      int lnbowner = -1;
      if(lid==llid) // owner
        lnbowner = myrank_;
      int nbowner = -1;
      scatra_discret_->Comm().MaxAll(&lnbowner,&nbowner,1);
      DRT::Element* neighborele = scatra_discret_->gElement(gtoeva[s]);
      CheckNeighborsDiffGrad(neighborele,nbowner,setsids,set,diffval,interval,auxvals);
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionOptiSplit::PerformIteration()
{
  if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"REACTION LINE SEARCH"<<std::endl;
    std::cout<<std::endl;
  }

  reacordifforcorrho_ = 0;
  bool reacsucc = false;
  int reacseq = sequenzeiter_;
  if(iter_==0)
    reacseq *= 5;
  for(int i=0; i<reacseq; ++i)
  {
    if(!myrank_)
      std::cout<<"ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,reac_objgrad_,reac_searchdirection_->ComputeDirection(reac_objgrad_,reac_vals_,iter_),reac_vals_,scatra_discret_->ElementRowMap());
    reacsucc = linesearch_->Run();

    if(!myrank_)
    {
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** objective function value          "<<J_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** error value                       "<<error_<<std::endl;
      std::cout<<"*** output count                      "<<output_count_<<std::endl;
    }
    if(reacsucc==false)
      break;
  }

  std::cout<<std::endl;
  std::cout<<"DIFFUSION LINE SEARCH"<<std::endl;
  std::cout<<std::endl;

  reacordifforcorrho_ = 1;
  bool diffsucc = false;
  for(int i=0; i<sequenzeiter_; ++i)
  {
    if(!myrank_)
      std::cout<<"ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,diff_objgrad_,diff_searchdirection_->ComputeDirection(diff_objgrad_,diff_vals_,iter_),diff_vals_,scatra_discret_->ElementRowMap());
    diffsucc = linesearch_->Run();

    if(!myrank_)
    {
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** objective function value          "<<J_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** error value                       "<<error_<<std::endl;
      std::cout<<"*** output count                      "<<output_count_<<std::endl;
    }
    if(diffsucc==false)
      break;
  }

  return (reacsucc || diffsucc);
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::CalculateGradDirNorm(const Epetra_Vector& bvector, const Epetra_Map& uniquemap, double* result)
{
  if(reacordifforcorrho_==0)
    reac_objgrad_->Dot(bvector,result);
  else if(reacordifforcorrho_==1)
    diff_objgrad_->Dot(bvector,result);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplit::FDCheck()
{
  std::cout<<"FDCHECK"<<std::endl;
  if(scatra_discret_->Comm().NumProc()>1)
    dserror("FDCHECK only implemented for one processor");

  // reaction part
  {
    reacordifforcorrho_ = 0;
    double J_before = J_;
    std::cout<<"reaction gradient according to adjoint analysis"<<std::endl;
    reac_objgrad_->Print(std::cout);

    Epetra_Vector fd_reac_grad(*scatra_discret_->ElementRowMap(),false);
    Teuchos::RCP<Epetra_Vector> perturb_reac_vals = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> reac_vals_before = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    reac_vals_before->Update(1.0,*reac_vals_,0.0);

    for(int i=0; i<reac_vals_->MyLength(); ++i)
    {
      double perturba = 1.0e-3;
      double perturbb = 1.0e-4;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = reac_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_reac_vals->Update(1.0,*reac_vals_before,0.0);
      perturb_reac_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_reac_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_reac_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"reaction gradient according to FD analysis"<<std::endl;
    fd_reac_grad.Print(std::cout);

    ReplaceParams(reac_vals_before);
    J_ = J_before;
  }

  // diffusion part
  {
    reacordifforcorrho_ = 1;
    double J_before = J_;
    std::cout<<"diffusion gradient according to adjoint analysis"<<std::endl;
    diff_objgrad_->Print(std::cout);

    Epetra_Vector fd_diff_grad(*scatra_discret_->ElementRowMap(),false);
    Teuchos::RCP<Epetra_Vector> perturb_diff_vals = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> diff_vals_before = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    diff_vals_before->Update(1.0,*diff_vals_,0.0);

    for(int i=0; i<diff_vals_->MyLength(); ++i)
    {
      double perturba = 1.0e-3;
      double perturbb = 1.0e-4;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = diff_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_diff_vals->Update(1.0,*diff_vals_before,0.0);
      perturb_diff_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_diff_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_diff_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"diffusion gradient according to FD analysis"<<std::endl;
    fd_diff_grad.Print(std::cout);

    ReplaceParams(diff_vals_before);
    J_ = J_before;
  }

  return;
}


/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionOptiSplitAcouSplit::PatImageReconstructionOptiSplitAcouSplit(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstructionOptiSplit(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout)
{

  // read the material ids
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"ACOUPARAMLIST"));
  char* pEnd;
  while (pstream >> word2)
  {
    int id = std::strtol(&word2[0],&pEnd,10);
    if (*pEnd=='\0')
      acou_matids_.push_back(id);
  }

  // create value vector
  acou_opt_ind_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));
  for(int e=0; e<acou_discret_->NumMyRowElements(); ++e)
  {
    DRT::Element* acou_ele = acou_discret_->lRowElement(e);
    int elematid = acou_ele->Material()->Parameter()->Id();    // check if element has material that is optimized
    for(unsigned int i=0; i<acou_matids_.size(); ++i)
    {
      if(acou_matids_[i]==elematid)
      {
        acou_opt_ind_->ReplaceMyValue(e,0,1.0);
        break;
      }
    }
  }

  // init int for identification of update
  reacordifforcorrho_ = 0;

  // create value vector
  c_vals_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));

  // create value vector
  rho_vals_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));

  // fill values vector
  for(int e=0; e<acou_discret_->NumMyRowElements(); ++e)
  {
    DRT::Element* acou_ele = acou_discret_->lRowElement(e);
    c_vals_->ReplaceMyValue(e,0,acou_ele->Material()->Parameter()->GetParameter(1,-1));
    rho_vals_->ReplaceMyValue(e,0,acou_ele->Material()->Parameter()->GetParameter(0,-1));
  }

  // allocate gradients
  c_objgrad_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));
  rho_objgrad_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));

  // setup directions
  c_searchdirection_ = Teuchos::rcp(new PATSearchDirection(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")));
  c_searchdirection_->Setup(acou_discret_->ElementColMap(),acou_discret_->ElementRowMap());
  rho_searchdirection_ = Teuchos::rcp(new PATSearchDirection(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")));
  rho_searchdirection_->Setup(acou_discret_->ElementColMap(),acou_discret_->ElementRowMap());

  // create regularization
  if(DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE")!=INPAR::ACOU::pat_regula_none)
  {
    c_regula_ = Teuchos::rcp(new PATRegula(
          DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE"),
          acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TIKHWEIGHT_C"),
          acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDWEIGHT_C"),
          acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDEPS_C"),
          acou_discret_));
    rho_regula_ = Teuchos::rcp(new PATRegula(
        DRT::INPUT::IntegralValue<INPAR::ACOU::RegulaType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"REGULATYPE"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TIKHWEIGHT_RHO"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDWEIGHT_RHO"),
        acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("TVDEPS_RHO"),
        acou_discret_));
  }

  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",true);

  /* only for optimization with "correctly set" acoustical properties */
  if(0)
  {
    double mu_1=0.01;
    double D_1=0.5;
    double c_1=1.54;
    double rho_1=1.1;

    double mu_2=0.1;
    double D_2=0.3;
    double c_2=1.8;
    double rho_2=1.2;

    c_vals_->PutScalar(1.48);
    rho_vals_->PutScalar(1.00);
    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {
      // local id
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get element center coordinates
      std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(scatra_discret_->lRowElement(lid));

      // determine value to set
      double f = std::exp(-0.5*xyz[0]*xyz[0]-0.25*xyz[1]*xyz[1]);
      double mu_tat = (1.0-f)*mu_1+f*mu_2;
      double D_tat = (1.0-f)*D_1+f*D_2;
      double c_tat = (1.0-f)*c_1+f*c_2;
      double rho_tat = (1.0-f)*rho_1+f*rho_2;

      reac_vals_->ReplaceMyValue(lid,0,mu_tat);
      diff_vals_->ReplaceMyValue(lid,0,D_tat);
      c_vals_->ReplaceMyValue(lid,0,c_tat);
      rho_vals_->ReplaceMyValue(lid,0,rho_tat);
    }
    reacordifforcorrho_ = 0;
    ReplaceParams(reac_vals_);
    reacordifforcorrho_ = 1;
    ReplaceParams(diff_vals_);
    reacordifforcorrho_ = 2;
    ReplaceParams(c_vals_);
    reacordifforcorrho_ = 3;
    ReplaceParams(rho_vals_);

    acouparams_->set<bool>("acouopt",false);
  }

  if(0)
  {
    /*
    //acouknown strhet and lowhet
    double c_firstcircle = 1.6;
    //double c_firstcircle =  1.8;
    double c_secondcircle = 1.6;
    //double c_secondcircle = 1.8;
    double c_rect = 1.35;
    //double c_rect = 1.0;
    double c_soft = 1.54;
    //double c_soft = 1.6;
    double c_def = 1.48;
    double rho_firstcircle = 1.2;
    //double rho_firstcircle = 1.5;
    double rho_secondcircle = 1.2;
    //double rho_secondcircle = 1.5;
    double rho_rect = 1.0;
    //double rho_rect = 1.0;
    double rho_soft = 1.1;
    //double rho_soft = 1.2;
    double rho_def = 1.0;
    //double rho_def = 1.0;
    */
    /*
    // strhet good acouwrong
    double c_firstcircle = 1.6;
    double c_secondcircle = 1.6;
    double c_rect = 1.6;
    double c_soft =  1.6;
    double c_def = 1.48;
    double rho_firstcircle = 1.2;
    double rho_secondcircle = 1.2;
    double rho_rect = 1.2;
    double rho_soft = 1.2;
    double rho_def = 1.0;
    */

    // strhet good acouknown
    double c_firstcircle = 1.9;
    double c_secondcircle = 1.9;
    double c_rect = 1.8;
    double c_soft =  1.6;
    double c_def = 1.48;
    double rho_firstcircle = 1.5;
    double rho_secondcircle = 1.5;
    double rho_rect = 1.0;
    double rho_soft = 1.2;
    double rho_def = 1.0;

    /* strhet acouwrong
    double c_firstcircle = 1.6;
    double c_secondcircle = 1.6;
    double c_rect = 1.6;
    double c_soft =  1.6;
    double c_def = 1.48;
    double rho_firstcircle = 1.2;
    double rho_secondcircle = 1.2;
    double rho_rect = 1.2;
    double rho_soft = 1.2;
    double rho_def = 1.0;
    */
    /* lowhet acouwrong
    double c_firstcircle = 1.54;
    double c_secondcircle = 1.54;
    double c_rect = 1.54;
    double c_soft =  1.54;
    double c_def = 1.48;
    double rho_firstcircle = 1.1;
    double rho_secondcircle = 1.1;
    double rho_rect = 1.1;
    double rho_soft = 1.1;
    double rho_def = 1.0;
    */

    Teuchos::RCP<Epetra_Vector> c_tatparams = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));
    Teuchos::RCP<Epetra_Vector> rho_tatparams = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));

    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {

      double c_val = 0.0; // to be determined
      double rho_val = 0.0; // to be determined

      // local id
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get element center coordinates
      std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(scatra_discret_->lRowElement(lid));

      // check first circle:
      double p = sqrt((xyz[0]-2.)*(xyz[0]-2.)+(xyz[1]-1.5)*(xyz[1]-1.5));
      if(p<1.0)
      {
        c_val = c_firstcircle;
        rho_val = rho_firstcircle;
      }
      else
      {
        // check second circle
        p = sqrt((xyz[0]+2.)*(xyz[0]+2.)+(xyz[1]-1.5)*(xyz[1]-1.5));
        if(p<1.25)
        {
          c_val = c_secondcircle;
          rho_val = rho_secondcircle;
        }
        else
        {
          // check rectangle
          double g1 = 0.176327*(xyz[0]+1.260151405)-1.49148196;
          double g2 = 0.176327*(xyz[0]+0.99967914)-2.968693585;
          double g3 = -5.671281835*(xyz[0]+0.99967914)-2.968693585;
          double g4 = -5.671281835*(xyz[0]-2.93955187)-2.274100875;

          if(xyz[1]<g1 && xyz[1]>g2 && xyz[1]>g3 && xyz[1]<g4)
          {
            c_val = c_rect;
            rho_val = rho_rect;
          }
          else if( sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])<5. ) // check soft tissue circle
          {
            c_val = c_soft;
            rho_val = rho_soft;
          }
          else
          {
            c_val = c_def;
            rho_val = rho_def;
          }
        }
      }

      c_vals_->ReplaceMyValue(lid,0,c_val);
      rho_vals_->ReplaceMyValue(lid,0,rho_val);
    }

    reacordifforcorrho_ = 2;
    ReplaceParams(c_vals_);
    reacordifforcorrho_ = 3;
    ReplaceParams(rho_vals_);

    acouparams_->set<bool>("acouopt",false);
  }


  if(0) // forward calculation for mixed material
  {
    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {
      // local id
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get element center coordinates
      std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(scatra_discret_->lRowElement(lid));

      if(sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])>1.86)
      {
        // set default parameters
        reac_vals_->ReplaceMyValue(lid,0,1.);
        diff_vals_->ReplaceMyValue(lid,0,0.6);
        c_vals_->ReplaceMyValue(lid,0,1.);
        rho_vals_->ReplaceMyValue(lid,0,1.);
      }
      else
      {
        double expval=exp(-xyz[0]*xyz[0]-xyz[1]*xyz[1]);
        reac_vals_->ReplaceMyValue(lid,0,(1-expval)*1.+expval*2.);
        diff_vals_->ReplaceMyValue(lid,0,(1-expval)*0.6+expval*0.8);
        c_vals_->ReplaceMyValue(lid,0,(1-expval)*1.+expval*1.2);
        rho_vals_->ReplaceMyValue(lid,0,(1-expval)*1.+expval*1.5);
      }
    }
    reacordifforcorrho_ = 0;
    ReplaceParams(reac_vals_);
    reacordifforcorrho_ = 1;
    ReplaceParams(diff_vals_);

    reacordifforcorrho_ = 2;
    ReplaceParams(c_vals_);
    reacordifforcorrho_ = 3;
    ReplaceParams(rho_vals_);

  }
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::ReplaceParams(Teuchos::RCP<Epetra_Vector> params)
{
  if(reacordifforcorrho_==0)
    for(int i=0; i<params->MyLength(); ++i)
    {
      ;
      //if(params->operator [](i)<0.0)
        //params->ReplaceMyValue(i,0,0.0);
    }
  else
    for(int i=0; i<params->MyLength(); ++i)
    {
      if(params->operator [](i)<0.1)
        params->ReplaceMyValue(i,0,0.1);
      else if(params->operator [](i)>3.0)
        params->ReplaceMyValue(i,0,3.0);
    }

  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if(reacordifforcorrho_==0)
  {
    // export params to col map
    Teuchos::RCP<Epetra_Vector> paramscol = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementColMap()),false));
    LINALG::Export(*params,*paramscol);

    reac_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
      actmat->Parameter()->SetParameter(1,paramscol);
    }
  }
  else if(reacordifforcorrho_==1)
  {
    Teuchos::RCP<Epetra_Vector> paramscol = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementColMap()),false));
    LINALG::Export(*params,*paramscol);

    diff_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
      actmat->Parameter()->SetParameter(0,paramscol);
    }
  }
  else if(reacordifforcorrho_==2)
  {
    Teuchos::RCP<Epetra_Vector> paramscol = Teuchos::rcp(new Epetra_Vector(*(acou_discret_->ElementColMap()),false));
    LINALG::Export(*params,*paramscol);

    c_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<acou_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(acou_matids_[i]);
      actmat->Parameter()->SetParameter(1,paramscol);
    }
  }
  else if(reacordifforcorrho_==3)
  {
    Teuchos::RCP<Epetra_Vector> paramscol = Teuchos::rcp(new Epetra_Vector(*(acou_discret_->ElementColMap()),false));
    LINALG::Export(*params,*paramscol);

    rho_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<acou_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(acou_matids_[i]);
      actmat->Parameter()->SetParameter(0,paramscol);
    }
  }

  // update node based vector
  if(reacordifforcorrho_==0)
    ComputeNodeBasedReactionCoefficient();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos)
{
  reacordifforcorrho_ = 0;
  ReplaceParams(reacs);

  reacordifforcorrho_ = 1;
  ReplaceParams(diffs);

  reacordifforcorrho_ = 2;
  ReplaceParams(cs);

  reacordifforcorrho_ = 3;
  ReplaceParams(rhos);

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstructionOptiSplitAcouSplit::EvalulateObjectiveFunction()
{
  // evaluate error contribution
  EvaluateError();
  J_ = error_;

  if(reac_regula_!=Teuchos::null)
    reac_regula_->Evaluate(reac_vals_,&J_);
  if(diff_regula_!=Teuchos::null)
    diff_regula_->Evaluate(diff_vals_,&J_);

  if(c_regula_!=Teuchos::null)
    c_regula_->Evaluate(c_vals_,&J_);
  if(rho_regula_!=Teuchos::null)
    rho_regula_->Evaluate(rho_vals_,&J_);

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<J_-error_<<std::endl;

  return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::EvaluateGradient()
{
  // zero out gradient vector initially
  reac_objgrad_->Scale(0.0);
  diff_objgrad_->Scale(0.0);
  c_objgrad_->Scale(0.0);
  rho_objgrad_->Scale(0.0);

  // set quantities needed by the elements
  scatra_discret_->SetState("adjoint phi",adjoint_phi_);

  // fill and set psi vector
  Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  psi = CalculateAdjointOptiRhsvec(adjoint_psi_);
  scatra_discret_->SetState("psi",psi);

  // do the actual evaluation
  EvaluateReacGrad();
  EvaluateDiffGrad();
  EvaluateCGrad();
  EvaluateRhoGrad();

  // check gradient if required
  if(fdcheck_)
    FDCheck();

  // reduce basis
  ReduceBasis();

//  double norm1 = 0.0, norm2 = 0.0, norminf=0.0;
//  reac_objgrad_->Norm1(&norm1); reac_objgrad_->Norm2(&norm2); reac_objgrad_->NormInf(&norminf);
//  std::cout<<"reac "<<norm1<<" "<<norm2<<" "<<norminf<<std::endl;
//  diff_objgrad_->Norm1(&norm1); diff_objgrad_->Norm2(&norm2); diff_objgrad_->NormInf(&norminf);
//  std::cout<<"diff "<<norm1<<" "<<norm2<<" "<<norminf<<std::endl;
//  c_objgrad_->Norm1(&norm1); c_objgrad_->Norm2(&norm2); c_objgrad_->NormInf(&norminf);
//  std::cout<<"sos  "<<norm1<<" "<<norm2<<" "<<norminf<<std::endl;
//  rho_objgrad_->Norm1(&norm1); rho_objgrad_->Norm2(&norm2); rho_objgrad_->NormInf(&norminf);
//  std::cout<<"rho  "<<norm1<<" "<<norm2<<" "<<norminf<<std::endl;
//
  // output the solution
  std::string outname = name_;
  outname.append("_c_and_rho_grad");
  acououtput_->NewResultFile(outname,0);
  output_count_++;
  acououtput_->WriteMesh(0,0.0);
  acououtput_->NewStep(1,1.0);
  acououtput_->WriteElementData(true);
  acououtput_->WriteVector("density",rho_objgrad_);
  acououtput_->WriteVector("speedofsound",c_objgrad_);
//
//  std::string soutname = name_;
//  soutname.append("_reac_and_diff_grad");
//  scatraoutput_->NewResultFile(soutname,0);
//  output_count_++;
//  scatraoutput_->WriteMesh(0,0.0);
//  scatraoutput_->NewStep(1,1.0);
//  scatraoutput_->WriteElementData(true);
//  scatraoutput_->WriteVector("rea_coeff",reac_objgrad_);
//  scatraoutput_->WriteVector("diff_coeff",diff_objgrad_);
//
//  if(reducedbasis_ && setidsdiff != Teuchos::null) dserror("und stopp");

  return;
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::ReduceBasis()
{
  // first part: reduce the basis of the absorption coefficient, it always uses its own gradient
  if(patchtype_==INPAR::ACOU::pat_patch_none)
    return;

  // security check
  if(meshconform_==false) dserror("reduced basis not implemented for nonconforming meshes");

  // reduce the basis for the reaction coefficient first
  // create a vector to store the set ids
  Teuchos::RCP<Epetra_Vector> setids = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  if(patchtype_!=INPAR::ACOU::pat_patch_mixed)
  {
    double maxval = 0.0;
    double minval = 0.0;
    reac_objgrad_->MaxValue(&maxval);
    reac_objgrad_->MinValue(&minval);
    if(maxval == minval) dserror("tried to reduce basis, but gradient has only one value");
    double rangeval = std::abs(maxval - minval);

    // create a helper vector
    Teuchos::RCP<Epetra_Vector> auxvals = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
    auxvals->Update(1.0,*reac_objgrad_,0.0);

    // find minid
    int minid = -1;
    for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
    {
      if(reac_objgrad_->operator [](e)<=minval+1.0e-10 && opti_opt_ind_->operator [](e)!= 0.0)
        minid = e;
    }
    int global_minid = -1;
    scatra_discret_->Comm().MaxAll(&minid,&global_minid,1);
    int loc_owner = -1;
    if(minid==global_minid)
      loc_owner = myrank_;
    int owner = -1;
    scatra_discret_->Comm().MaxAll(&loc_owner,&owner,1);

    setids->PutScalar(-1.0);
    double minvalsetids = -1.0;

    // set all the set ids
    int i=0;
    while(minvalsetids<0.0)
    {
      double set = double(i);
      DRT::Element* actele = NULL;
      if(myrank_==owner) // only for owning processor
      {
        setids->ReplaceMyValue(minid,0,set);
        auxvals->ReplaceMyValue(minid,0,std::numeric_limits<double>::max());
        actele = scatra_discret_->lRowElement(minid);
      }
      CheckNeighborsReacGrad(actele,owner,setids,set,minval,0.1*rangeval,auxvals);

      // find next minimum value
      auxvals->MinValue(&minval);

      // find minid
      minid = -1;
      for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
      {
        if(auxvals->operator [](e)<=minval+1.0e-10 && opti_opt_ind_->operator [](e)!= 0.0)
          minid = e;
      }
      global_minid = -1;
      scatra_discret_->Comm().MaxAll(&minid,&global_minid,1);
      loc_owner = -1;
      if(minid==global_minid)
        loc_owner = myrank_;
      owner = -1;
      scatra_discret_->Comm().MaxAll(&loc_owner,&owner,1);

      setids->MinValue(&minvalsetids);
      i++;
    }

    if(!myrank_)
      std::cout<<"identified "<<i<<" sets using the reaction gradient for the reaction basis"<<std::endl;

    // now recalculate the entries in the gradients according to the sets
    for(int j=0; j<i; ++j) // i holds the number of sets right now
    {
      double lsetvalreac = 0.0;
      int lnumsetval = 0;
      for(int g=0; g<reac_objgrad_->MyLength(); ++g)
      {
        if(opti_opt_ind_->operator [](g)!= 0.0)
        {
          double reacgradval = reac_objgrad_->operator [](g);
          int set = setids->operator [](g);
          if(set==j)
          {
            lsetvalreac += reacgradval;
            lnumsetval++;
          }
        }
      }

      double gsetvalreac = 0.0;
      scatra_discret_->Comm().SumAll(&lsetvalreac,&gsetvalreac,1);
      int gnumsetval = 0;
      scatra_discret_->Comm().SumAll(&lnumsetval,&gnumsetval,1);

      if(gnumsetval!=0.0)
        gsetvalreac/=gnumsetval;

      for(int g=0; g<reac_objgrad_->MyLength(); ++g)
      {
        if(opti_opt_ind_->operator [](g)!= 0.0)
        {
          int set = setids->operator [](g);
          if(set==j)
            reac_objgrad_->ReplaceMyValue(g,0,gsetvalreac);
        }
      }
      scatra_discret_->Comm().Barrier();
    }
  } // reduction of reaction basis

  if(patchtype_==INPAR::ACOU::pat_patch_reacgrad)
  {
    // this is the easy case, because all other gradient just use the same patch as the absorption coefficient
    // but this makes only sense for sequenze_=1
    if(sequenzeiter_!=1 && !myrank_)
      std::cout<<"warning: patchreac makes only sense for sequenze_=1"<<std::endl;

    ReducedGradientsCalculation(setids,false);

  } // ** if(patchtype_==INPAR::ACOU::pat_patch_reacgrad)
  else if(patchtype_==INPAR::ACOU::pat_patch_self) // patch self
  {
    ReducePatchSelf(scatra_discret_,diff_objgrad_,opti_opt_ind_);
    ReducePatchSelf(acou_discret_,c_objgrad_,acou_opt_ind_);
    ReducePatchSelf(acou_discret_,rho_objgrad_,acou_opt_ind_);
  } // ** else if(patchtype_==INPAR::ACOU::pat_patch_self)
  else if(patchtype_==INPAR::ACOU::pat_patch_reacvals)
  {
    ReducePatchReacVals();
  }
  else if(patchtype_==INPAR::ACOU::pat_patch_mixed)
  {
    // diffusion coefficient with reaction gradient basis and speed of sound and mass density with reaction values

    // 1.) diffusion gradient
    /*double doublemaxset = 0.0;
    setids->MaxValue(&doublemaxset);
    for(int j=0; j<=int(doublemaxset); ++j) // i holds the number of sets right now
    {
      double lsetvaldiff = 0.0;
      int lnumsetvaldiff = 0;
      for(int g=0; g<diff_objgrad_->MyLength(); ++g)
        if(opti_opt_ind_->operator [](g)!=0.0)
        {
          double diffgradval = diff_objgrad_->operator [](g);
          int set = setids->operator [](g);
          if(set==j)
          {
            lsetvaldiff += diffgradval;
            lnumsetvaldiff++;
          }
        }

      double gsetvaldiff = 0.0;
      scatra_discret_->Comm().SumAll(&lsetvaldiff,&gsetvaldiff,1);
      int gnumsetvaldiff = 0;
      scatra_discret_->Comm().SumAll(&lnumsetvaldiff,&gnumsetvaldiff,1);

      if(gnumsetvaldiff!=0.0)
        gsetvaldiff/=gnumsetvaldiff;

      for(int g=0; g<diff_objgrad_->MyLength(); ++g)
        if(opti_opt_ind_->operator [](g)!=0.0)
        {
          int set = setids->operator [](g);
          if(set==j)
            diff_objgrad_->ReplaceMyValue(g,0,gsetvaldiff);
        }

      scatra_discret_->Comm().Barrier();
    }*/

    // 2.) acoustic gradients
    ReducePatchReacVals();
  }
  else
    dserror("invalid patch build type");

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::ReducedGradientsCalculation(Teuchos::RCP<Epetra_Vector> setids, bool acouonly)
{
  // create a vector to store the set ids
  Teuchos::RCP<Epetra_Vector> setidsacou = Teuchos::rcp(new Epetra_Vector(*(acou_discret_->ElementRowMap()),true));

  // calculate the acou setids from the reaction setids!
  for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
  {
    double lsetval = -1.0;
    if(scatra_discret_->ElementRowMap()->LID(i)>=0)
    {
      lsetval = setids->operator [](scatra_discret_->ElementRowMap()->LID(i));
    }
    double setval = 0.0;
    scatra_discret_->Comm().MaxAll(&lsetval,&setval,1);

    int agid = i - scatra_discret_->ElementRowMap()->MinAllGID() + acou_discret_->ElementRowMap()->MinAllGID();
    int alid = acou_discret_->ElementRowMap()->LID(agid);
    if(alid>=0)
      setidsacou->ReplaceMyValue(alid,0,setval);
  }

  // now recalculate the entries in the gradients according to the sets
  double doublemaxset = 0.0;
  setidsacou->MaxValue(&doublemaxset);
  for(int j=0; j<=int(doublemaxset); ++j) // i holds the number of sets right now
  {
    double lsetvaldiff = 0.0;
    double lsetvalc = 0.0;
    double lsetvalrho = 0.0;
    int lnumsetvalacou = 0;
    int lnumsetvaldiff = 0;
    for(int g=0; g<c_objgrad_->MyLength(); ++g)
      if(acou_opt_ind_->operator [](g)!=0.0) // only for those who play a role (not the water region)
      {
        double cgradval = c_objgrad_->operator [](g);
        double rhogradval = rho_objgrad_->operator [](g);
        int set = setidsacou->operator [](g);
        if(set==j)
        {
          lsetvalc += cgradval;
          lsetvalrho += rhogradval;
          lnumsetvalacou++;
        }
      }
    if(acouonly == false)
      for(int g=0; g<diff_objgrad_->MyLength(); ++g)
        if(opti_opt_ind_->operator [](g)!=0.0)
        {
          double diffgradval = diff_objgrad_->operator [](g);
          int set = setids->operator [](g);
          if(set==j)
          {
            lsetvaldiff += diffgradval;
            lnumsetvaldiff++;
          }
        }

    double gsetvalc = 0.0;
    double gsetvalrho = 0.0;
    double gsetvaldiff = 0.0;
    scatra_discret_->Comm().SumAll(&lsetvalc,&gsetvalc,1);
    scatra_discret_->Comm().SumAll(&lsetvalrho,&gsetvalrho,1);
    scatra_discret_->Comm().SumAll(&lsetvaldiff,&gsetvaldiff,1);
    int gnumsetvalacou = 0;
    int gnumsetvaldiff = 0;
    scatra_discret_->Comm().SumAll(&lnumsetvalacou,&gnumsetvalacou,1);
    scatra_discret_->Comm().SumAll(&lnumsetvaldiff,&gnumsetvaldiff,1);

    if(gnumsetvalacou!=0.0)
    {
      gsetvalc/=gnumsetvalacou;
      gsetvalrho/=gnumsetvalacou;
    }
    if(gnumsetvaldiff!=0.0)
      gsetvaldiff/=gnumsetvaldiff;

    for(int g=0; g<c_objgrad_->MyLength(); ++g)
      if(acou_opt_ind_->operator [](g)!=0.0) // only for those who play a role (not the water region)
      {
        int set = setidsacou->operator [](g);
        if(set==j)
        {
          c_objgrad_->ReplaceMyValue(g,0,gsetvalc);
          rho_objgrad_->ReplaceMyValue(g,0,gsetvalrho);
        }
      }
    if(acouonly == false)
      for(int g=0; g<diff_objgrad_->MyLength(); ++g)
        if(opti_opt_ind_->operator [](g)!=0.0)
        {
          int set = setids->operator [](g);
          if(set==j)
            diff_objgrad_->ReplaceMyValue(g,0,gsetvaldiff);
        }

    scatra_discret_->Comm().Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::ReducePatchSelf(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> patchvec, Teuchos::RCP<Epetra_Vector> opt_ind)
{

  // for reduced basis, the diffusion coefficient has to build patches according to the absorption coefficient distribution
  double maxval = 0.0;
  double minval = 0.0;
  patchvec->MaxValue(&maxval);
  patchvec->MinValue(&minval);

  double rangeval = std::abs(maxval - minval);
  if(rangeval == 0.0)
    if(!myrank_)
    dserror("rangeval is zero");

  // create a helper vector
  Teuchos::RCP<Epetra_Vector> auxvals = Teuchos::rcp(new Epetra_Vector(*(discret->ElementRowMap()),false));
  auxvals->Update(1.0,*patchvec,0.0);

  // find maxid
  int minid = -1;
  for(int e=0; e<discret->NumMyRowElements(); ++e)
  {
    if(patchvec->operator [](e)<=minval+1.0e-10 && opt_ind->operator [](e)!= 0.0)
      minid = e;
  }
  int global_minid = -1;
  discret->Comm().MaxAll(&minid,&global_minid,1);
  int loc_owner = -1;
  if(minid==global_minid)
    loc_owner = myrank_;
  int owner = -1;
  discret->Comm().MaxAll(&loc_owner,&owner,1);

  // create a vector to store the set ids
  Teuchos::RCP<Epetra_Vector> setids = Teuchos::rcp(new Epetra_Vector(*(discret->ElementRowMap()),false));
  setids->PutScalar(-1.0);
  double minvalsetids = -1.0;

  // set all the set ids
  int i=0;
  while(minvalsetids<0.0)
  {
    double set = double(i);
    DRT::Element* actele = NULL;
    if(myrank_==owner) // only for owning processor
    {
      setids->ReplaceMyValue(minid,0,set);
      auxvals->ReplaceMyValue(minid,0,std::numeric_limits<double>::max());
      actele = discret->lRowElement(minid);
    }
    CheckNeighborsPatchSelf(discret,patchvec,actele,owner,setids,set,minval,0.1*rangeval,auxvals);

    // find next minimum value
    auxvals->MinValue(&minval);

    // find maxid
    minid = -1;
    for(int e=0; e<discret->NumMyRowElements(); ++e)
    {
      if(auxvals->operator [](e)<=minval+1.0e-10 && opt_ind->operator [](e)!= 0.0)
        minid = e;
    }
    global_minid = -1;
    discret->Comm().MaxAll(&minid,&global_minid,1);
    loc_owner = -1;
    if(minid==global_minid)
      loc_owner = myrank_;
    owner = -1;
    discret->Comm().MaxAll(&loc_owner,&owner,1);

    setids->MinValue(&minvalsetids);
    i++;
  }

  if(!myrank_)
    std::cout<<"identified "<<i<<" sets using the gradient for the basis"<<std::endl;

  // now recalculate the entries in the gradients according to the sets
  for(int j=0; j<i; ++j) // i holds the number of sets right now
  {
    double lsetval = 0.0;
    int lnumsetval = 0;
    for(int g=0; g<patchvec->MyLength(); ++g)
    {
      double gradval = patchvec->operator [](g);
      int set = setids->operator [](g);
      if(set==j)
      {
        lsetval += gradval;
        lnumsetval++;
      }
    }

    double gsetval = 0.0;
    discret->Comm().SumAll(&lsetval,&gsetval,1);
    int gnumsetval = 0;
    discret->Comm().SumAll(&lnumsetval,&gnumsetval,1);

    if(gnumsetval!=0.0)
      gsetval/=gnumsetval;

    for(int g=0; g<patchvec->MyLength(); ++g)
    {
      if(opt_ind->operator [](g)!= 0.0)
      {
        int set = setids->operator [](g);
        if(set==j)
          patchvec->ReplaceMyValue(g,0,gsetval);
      }
    }
    discret->Comm().Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::CheckNeighborsPatchSelf(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> patchvec, DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double val, double interval, Teuchos::RCP<Epetra_Vector> auxvals)
{
  // parallel version
  int lactelenodeids[4]={0,0,0,0};
  int gactelenodeids[4]={0,0,0,0};
  if(owner==myrank_)
  {
    if(actele->Shape()!=DRT::Element::quad4)
      dserror("distypes other than quad4 not yet implemented");

    for(int n=0; n<4; ++n)
      lactelenodeids[n]=actele->NodeIds()[n];
  }
  discret->Comm().MaxAll(&lactelenodeids[0],&gactelenodeids[0],4);

  for(int n=0; n<4; ++n)
  {
    std::vector<int> toevaluate;
    if(discret->HaveGlobalNode(gactelenodeids[n]))
    {
      DRT::Node* node = discret->gNode(gactelenodeids[n]);
      for(int e=0; e<node->NumElement(); ++e)
      {
        DRT::Element* neighborele = node->Elements()[e];

        // is it real neighbor (only if they share 2 nodes)
        int share = 0;
        for(int a=0; a<4; ++a)
          for(int b=0; b<4; ++b)
          {
            if(gactelenodeids[a]==neighborele->NodeIds()[b])
              share++;
          }

        if(share == 4) // same element -> skip
          continue;
        else if(share == 1) // not really connected
            continue;
        else if(share == 2) // neighbor element
        {
          // if already evaluated, skip
          if(discret->ElementRowMap()->LID(neighborele->Id())<0)
            continue;
          if(setsids->operator [](discret->ElementRowMap()->LID(neighborele->Id())) <= set && setsids->operator [](discret->ElementRowMap()->LID(neighborele->Id())) >= 0.0)
            continue;

          // determine reaction coefficient
          double neighborval = patchvec->operator [](discret->ElementRowMap()->LID(neighborele->Id()));
          if(abs(neighborval-val) <= interval)
          {
            setsids->ReplaceMyValue(discret->ElementRowMap()->LID(neighborele->Id()),0,set);
            auxvals->ReplaceMyValue(discret->ElementRowMap()->LID(neighborele->Id()),0,std::numeric_limits<double>::max());

            // this has to be checked and its neighbors too
            toevaluate.push_back(neighborele->Id());
          }
        }
        else
          dserror("how can two quad4 elements share exactly 3 nodes??");
      }
    }
    int lsize = toevaluate.size();
    int size = -1;
    discret->Comm().MaxAll(&lsize,&size,1);
    if(toevaluate.size()!=unsigned(size))
      toevaluate.resize(size,0);
    std::vector<int> gtoeva(size);
    discret->Comm().MaxAll(&toevaluate[0],&gtoeva[0],size);

    for(int s=0; s<size; ++s)
    {
      int llid = discret->ElementRowMap()->LID(gtoeva[s]);
      int lid = -1;
      discret->Comm().MaxAll(&llid,&lid,1);
      int lnbowner = -1;
      if(lid==llid) // owner
        lnbowner = myrank_;
      int nbowner = -1;
      discret->Comm().MaxAll(&lnbowner,&nbowner,1);
      DRT::Element* neighborele = discret->gElement(gtoeva[s]);
      CheckNeighborsPatchSelf(discret,patchvec,neighborele,nbowner,setsids,set,val,interval,auxvals);
    }
  }
  return;

}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::ReducePatchReacVals()
{
  Teuchos::RCP<Epetra_Vector> setids = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  double maxval = 0.0;
  double minval = 0.0;
  reac_vals_->MaxValue(&maxval);
  reac_vals_->MinValue(&minval);
  if(maxval == minval)
  {
    if(!myrank_)
      std::cout<<"warning: tried to reduce basis, but gradient has only one value"<<std::endl;
    return;
  }
  double rangeval = std::abs(maxval - minval);

  // create a helper vector
  Teuchos::RCP<Epetra_Vector> auxvals = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  auxvals->Update(1.0,*reac_vals_,0.0);

  // find maxid
  int maxid = -1;
  for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
  {
    if(reac_vals_->operator [](e)>=maxval-1.0e-10 && opti_opt_ind_->operator [](e)!= 0.0)
      maxid = e;
  }
  int global_maxid = -1;
  scatra_discret_->Comm().MaxAll(&maxid,&global_maxid,1);
  int loc_owner = -1;
  if(maxid==global_maxid)
    loc_owner = myrank_;
  int owner = -1;
  scatra_discret_->Comm().MaxAll(&loc_owner,&owner,1);

  setids->PutScalar(-1.0);
  double minvalsetids = -1.0;

  // set all the set ids
  int i=0;
  while(minvalsetids<0.0)
  {
    double set = double(i);
    DRT::Element* actele = NULL;
    if(myrank_==owner) // only for owning processor
    {
      setids->ReplaceMyValue(maxid,0,set);
      auxvals->ReplaceMyValue(maxid,0,-std::numeric_limits<double>::max());
      actele = scatra_discret_->lRowElement(maxid);
    }
    CheckNeighborsPatchReacVals(actele,owner,setids,set,maxval,0.1*rangeval,auxvals);

    // find next maximum value
    auxvals->MaxValue(&maxval);

    // find maxid
    maxid = -1;
    for(int e=0; e<scatra_discret_->NumMyRowElements(); ++e)
    {
      if(auxvals->operator [](e)>=maxval-1.0e-10 && opti_opt_ind_->operator [](e)!= 0.0)
        maxid = e;
    }
    global_maxid = -1;
    scatra_discret_->Comm().MaxAll(&maxid,&global_maxid,1);
    loc_owner = -1;
    if(maxid==global_maxid)
      loc_owner = myrank_;
    owner = -1;
    scatra_discret_->Comm().MaxAll(&loc_owner,&owner,1);

    setids->MinValue(&minvalsetids);
    i++;
  }

  if(!myrank_)
    std::cout<<"identified "<<i<<" sets using the reaction values"<<std::endl;

  if(patchtype_==INPAR::ACOU::pat_patch_mixed)
    ReducedGradientsCalculation(setids,true);
  else
    ReducedGradientsCalculation(setids,false);

  return;
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::CheckNeighborsPatchReacVals(DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double reacval, double interval, Teuchos::RCP<Epetra_Vector> auxvals)
{
  // parallel version
  int lactelenodeids[4]={0,0,0,0};
  int gactelenodeids[4]={0,0,0,0};
  if(owner==myrank_)
  {
    if(actele->Shape()!=DRT::Element::quad4)
      dserror("distypes other than quad4 not yet implemented");

    for(int n=0; n<4; ++n)
      lactelenodeids[n]=actele->NodeIds()[n];
  }
  scatra_discret_->Comm().MaxAll(&lactelenodeids[0],&gactelenodeids[0],4);

  for(int n=0; n<4; ++n)
  {
    std::vector<int> toevaluate;
    if(scatra_discret_->HaveGlobalNode(gactelenodeids[n]))
    {
      DRT::Node* node = scatra_discret_->gNode(gactelenodeids[n]);
      for(int e=0; e<node->NumElement(); ++e)
      {
        DRT::Element* neighborele = node->Elements()[e];

        // is it real neighbor (only if they share 2 nodes)
        int share = 0;
        for(int a=0; a<4; ++a)
          for(int b=0; b<4; ++b)
          {
            if(gactelenodeids[a]==neighborele->NodeIds()[b])
              share++;
          }

        if(share == 4) // same element -> skip
          continue;
        else if(share == 1) // not really connected
            continue;
        else if(share == 2) // neighbor element
        {
          // if already evaluated, skip
          if(scatra_discret_->ElementRowMap()->LID(neighborele->Id())<0)
            continue;
          if(setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) <= set && setsids->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id())) >= 0)
            continue;

          // determine reaction coefficient
          double neighborreac = reac_vals_->operator [](scatra_discret_->ElementRowMap()->LID(neighborele->Id()));
          if(abs(neighborreac-reacval) <= interval)
          {
            setsids->ReplaceMyValue(scatra_discret_->ElementRowMap()->LID(neighborele->Id()),0,set);
            auxvals->ReplaceMyValue(scatra_discret_->ElementRowMap()->LID(neighborele->Id()),0,-std::numeric_limits<double>::max());

            // this has to be checked and its neighbors too
            toevaluate.push_back(neighborele->Id());
          }
        }
        else
          dserror("this is strange");
      }
    }
    int lsize = toevaluate.size();
    int size = -1;
    scatra_discret_->Comm().MaxAll(&lsize,&size,1);
    if(toevaluate.size()!=unsigned(size))
      toevaluate.resize(size,0);
    std::vector<int> gtoeva(size);
    scatra_discret_->Comm().MaxAll(&toevaluate[0],&gtoeva[0],size);

    for(int s=0; s<size; ++s)
    {
      int llid = scatra_discret_->ElementRowMap()->LID(gtoeva[s]);
      int lid = -1;
      scatra_discret_->Comm().MaxAll(&llid,&lid,1);
      int lnbowner = -1;
      if(lid==llid) // owner
        lnbowner = myrank_;
      int nbowner = -1;
      scatra_discret_->Comm().MaxAll(&lnbowner,&nbowner,1);
      DRT::Element* neighborele = scatra_discret_->gElement(gtoeva[s]);
      CheckNeighborsPatchReacVals(neighborele,nbowner,setsids,set,reacval,interval,auxvals);
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::EvaluateCGrad()
{
  // loop the row elements
  for (int i=0; i<acou_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele = acou_discret_->lRowElement(i);

    // we do not have to do the calculation for elements which are not optimized
    if (acou_opt_ind_->operator [](i)==0.0)
      continue;

    const DRT::ELEMENTS::Acou * hdgele = dynamic_cast<const DRT::ELEMENTS::Acou*>(actele);
    double val = hdgele->GetSoSGradient();

    // write it to the gradient
    c_objgrad_->ReplaceMyValue(i,0,val);
  }

  // just to be safe
  c_objgrad_->Multiply(1.0,*acou_opt_ind_,*c_objgrad_,0.0);

  // regularization
  if(c_regula_!=Teuchos::null)
    c_regula_->EvaluateGradient(c_vals_,c_objgrad_);

  ConvertGradient(acou_discret_,c_objgrad_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::EvaluateRhoGrad()
{
  // loop the row elements
  for (int i=0; i<acou_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele = acou_discret_->lRowElement(i);

    // we do not have to do the calculation for elements which are not optimized
    if (acou_opt_ind_->operator [](i)==0.0)
      continue;

    const DRT::ELEMENTS::Acou * hdgele = dynamic_cast<const DRT::ELEMENTS::Acou*>(actele);
    double val = hdgele->GetDensityGradient();

    // write it to the gradient
    rho_objgrad_->ReplaceMyValue(i,0,val);
  }

  // just to be safe
  rho_objgrad_->Multiply(1.0,*acou_opt_ind_,*rho_objgrad_,0.0);

  // regularization
  if(rho_regula_!=Teuchos::null)
    rho_regula_->EvaluateGradient(rho_vals_,rho_objgrad_);

  ConvertGradient(acou_discret_,rho_objgrad_);

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionOptiSplitAcouSplit::PerformIteration()
{
  if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"REACTION LINE SEARCH"<<std::endl;
    std::cout<<std::endl;
  }

  int reacseq = sequenzeiter_;
  int diffseq = sequenzeiter_;
  int cseq    = sequenzeiter_;
  int rhoseq  = sequenzeiter_;


  if(patchtype_==INPAR::ACOU::pat_patch_reacvals&&iter_==0)
    reacseq = 1;
  else if(patchtype_==INPAR::ACOU::pat_patch_mixed&&iter_<10)
  {
    reacseq = 1;
    diffseq = 1;
    cseq = 0;
    rhoseq = 0;
  }

  reacordifforcorrho_ = 0;
  bool reacsucc = false;
  for(int i=0; i<reacseq; ++i)
  {
    if(!myrank_)
      std::cout<<"REACTION ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,reac_objgrad_,reac_searchdirection_->ComputeDirection(reac_objgrad_,reac_vals_,iter_),reac_vals_,scatra_discret_->ElementRowMap());
    reacsucc = linesearch_->Run();

    if(!myrank_)
    {
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** objective function value          "<<J_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** error value                       "<<error_<<std::endl;
      std::cout<<"*** output count                      "<<output_count_<<std::endl;
    }
    ComputeParameterError();
    if(reacsucc==false)
      break;
  }

  if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"SOUND SPEED LINE SEARCH"<<std::endl;
    std::cout<<std::endl;
  }
  reacordifforcorrho_ = 2;
  bool csucc = false;
  for(int i=0; i<cseq; ++i)
  {
    if(!myrank_)
      std::cout<<"SOUND SPEED ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,c_objgrad_,c_searchdirection_->ComputeDirection(c_objgrad_,c_vals_,iter_),c_vals_,acou_discret_->ElementRowMap());
    csucc = linesearch_->Run();

    if(!myrank_)
    {
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** objective function value          "<<J_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** error value                       "<<error_<<std::endl;
      std::cout<<"*** output count                      "<<output_count_<<std::endl;
    }
    ComputeParameterError();
    if(csucc==false)
      break;
  }

  if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"DENSITY LINE SEARCH"<<std::endl;
    std::cout<<std::endl;
  }
  reacordifforcorrho_ = 3;
  bool rhosucc = false;
  for(int i=0; i<rhoseq; ++i)
  {
    if(!myrank_)
      std::cout<<"DENSITY ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,rho_objgrad_,rho_searchdirection_->ComputeDirection(rho_objgrad_,rho_vals_,iter_),rho_vals_,acou_discret_->ElementRowMap());
    rhosucc = linesearch_->Run();

    if(!myrank_)
    {
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** objective function value          "<<J_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** error value                       "<<error_<<std::endl;
      std::cout<<"*** output count                      "<<output_count_<<std::endl;
    }
    ComputeParameterError();
    if(rhosucc==false)
      break;
  }

  if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"DIFFUSION LINE SEARCH"<<std::endl;
    std::cout<<std::endl;
  }
  reacordifforcorrho_ = 1;
  bool diffsucc = false;
  for(int i=0; i<diffseq; ++i)
  {
    if(!myrank_)
      std::cout<<"DIFFUSION ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,diff_objgrad_,diff_searchdirection_->ComputeDirection(diff_objgrad_,diff_vals_,iter_),diff_vals_,scatra_discret_->ElementRowMap());
    diffsucc = linesearch_->Run();

    if(!myrank_)
    {
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** objective function value          "<<J_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** error value                       "<<error_<<std::endl;
      std::cout<<"*** output count                      "<<output_count_<<std::endl;
    }
    ComputeParameterError();
    if(diffsucc==false)
      break;
  }

  return (reacsucc || diffsucc || csucc || rhosucc);
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::CalculateGradDirNorm(const Epetra_Vector& bvector, const Epetra_Map& uniquemap, double* result)
{
  if(reacordifforcorrho_==0)
    reac_objgrad_->Dot(bvector,result);
  else if(reacordifforcorrho_==1)
    diff_objgrad_->Dot(bvector,result);
  else if(reacordifforcorrho_==2)
    c_objgrad_->Dot(bvector,result);
  else if(reacordifforcorrho_==3)
    rho_objgrad_->Dot(bvector,result);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::SampleObjectiveFunction()
{
  double reac_firstcircle = 0.25;
  double reac_secondcircle = 0.25;
  double reac_rect = 0.1;
  double reac_soft = 0.01;
  double reac_def = 0.01; // 0.004;
  double D_firstcircle = 0.3;
  double D_secondcircle = 0.3;
  double D_rect = 0.6;
  double D_soft = 0.5;
  double D_def = 0.5; //0.1;
  /*
  //double c_firstcircle = 1.6;
  double c_firstcircle =  1.8;
  //double c_secondcircle = 1.6;
  double c_secondcircle = 1.8;
  //double c_rect = 1.35;
  double c_rect = 1.0;
  //double c_soft = 1.54;
  double c_soft = 1.6;
  double c_def = 1.48;
  //double rho_firstcircle = 1.2;
  double rho_firstcircle = 1.5;
  //double rho_secondcircle = 1.2;
  double rho_secondcircle = 1.5;
  //double rho_rect = 1.0;
  double rho_rect = 1.0;
  //double rho_soft = 1.1;
  double rho_soft = 1.2;
  //double rho_def = 1.0;
  double rho_def = 1.0;
  */
  double c_firstcircle = 1.9;
  double c_secondcircle = 1.9;
  double c_rect = 1.8;
  double c_soft =  1.6;
  double c_def = 1.48;
  double rho_firstcircle = 1.5;
  double rho_secondcircle = 1.5;
  double rho_rect = 1.0;
  double rho_soft = 1.2;
  double rho_def = 1.0;

  int rmax = 10;
  int smax = 10;

  for(int s=0; s<=smax+2; ++s)
  for(int r=0; r<=rmax+2; ++r)
  //for(int r=0; r<=rmax; ++r)
  {
    double ratiocorrect_r = double(r)/double(rmax);
    double ratiocorrect_s = double(s)/double(smax);

    std::cout<<"run r "<<r<<" run s "<<s<<" ratiocorrect r "<<ratiocorrect_r<<" ratiocorrect s "<<ratiocorrect_s<<std::endl;
    //std::cout<<"run r "<<r<<" ratiocorrect r "<<ratiocorrect_r<<std::endl;

    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {
      double reac_val = 0.0; // to be determined
      double D_val = 0.0; // to be determined
      double c_val = 0.0; // to be determined
      double rho_val = 0.0; // to be determined

      // local id
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get element center coordinates
      std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(scatra_discret_->lRowElement(lid));

      // check first circle:
      double p = sqrt((xyz[0]-2.)*(xyz[0]-2.)+(xyz[1]-1.5)*(xyz[1]-1.5));
      if(p<1.0)
      {
        reac_val = reac_firstcircle;
        D_val = D_firstcircle;
        c_val = c_firstcircle;
        rho_val = rho_firstcircle;
      }
      else
      {
        // check second circle
        p = sqrt((xyz[0]+2.)*(xyz[0]+2.)+(xyz[1]-1.5)*(xyz[1]-1.5));
        if(p<1.25)
        {
          reac_val = reac_secondcircle;
          D_val = D_secondcircle;
          c_val = c_secondcircle;
          rho_val = rho_secondcircle;
        }
        else
        {
          // check rectangle
          double g1 = 0.176327*(xyz[0]+1.260151405)-1.49148196;
          double g2 = 0.176327*(xyz[0]+0.99967914)-2.968693585;
          double g3 = -5.671281835*(xyz[0]+0.99967914)-2.968693585;
          double g4 = -5.671281835*(xyz[0]-2.93955187)-2.274100875;

          if(xyz[1]<g1 && xyz[1]>g2 && xyz[1]>g3 && xyz[1]<g4)
          {
            reac_val = reac_rect;
            D_val = D_rect;
            c_val = c_rect;
            rho_val = rho_rect;
          }
          else if( sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])<5. ) // check soft tissue circle
          {
            reac_val = reac_soft;
            D_val = D_soft;
            c_val = c_soft;
            rho_val = rho_soft;
          }
          else
          {
            reac_val = reac_def;
            D_val = D_def;
            c_val = c_def;
            rho_val = rho_def;
          }
        }
      }
      //reac_val = reac_def;
      reac_val = ratiocorrect_r * reac_val + (1.0-ratiocorrect_r)*reac_soft;
      //D_val = D_def;
      //D_val = ratiocorrect_r * D_val + (1.0-ratiocorrect_r)*D_soft;
      //c_val = c_soft;
      c_val = ratiocorrect_s * c_val + (1.0-ratiocorrect_s)*c_soft;
      //rho_val = rho_soft ;
      //rho_val = ratiocorrect_r * rho_val + (1.0-ratiocorrect_r)*rho_soft;

      reac_vals_->ReplaceMyValue(lid,0,reac_val);
      diff_vals_->ReplaceMyValue(lid,0,D_val);
      c_vals_->ReplaceMyValue(lid,0,c_val);
      rho_vals_->ReplaceMyValue(lid,0,rho_val);
    }
    reacordifforcorrho_ = 0;
    ReplaceParams(reac_vals_);
    reacordifforcorrho_ = 1;
    ReplaceParams(diff_vals_);
    reacordifforcorrho_ = 2;
    ReplaceParams(c_vals_);
    reacordifforcorrho_ = 3;
    ReplaceParams(rho_vals_);

    SolveStandardScatra();
    SolveStandardAcou();
    EvalulateObjectiveFunction();

  }
  dserror("that is it");

  return;
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::FDCheck()
{
  // reaction part
  {
    reacordifforcorrho_ = 0;
    double J_before = J_;
    std::cout<<"reaction gradient according to adjoint analysis"<<std::endl;
    reac_objgrad_->Print(std::cout);

    Epetra_Vector fd_reac_grad(*scatra_discret_->ElementRowMap(),false);
    Teuchos::RCP<Epetra_Vector> perturb_reac_vals = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> reac_vals_before = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    reac_vals_before->Update(1.0,*reac_vals_,0.0);

    for(int i=0; i<reac_vals_->MyLength(); ++i)
    {
      double perturba = 0.0; //1.0e-3;
      double perturbb = 1.0e-6;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = reac_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_reac_vals->Update(1.0,*reac_vals_before,0.0);
      perturb_reac_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_reac_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_reac_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"reaction gradient according to FD analysis"<<std::endl;
    fd_reac_grad.Print(std::cout);

    ReplaceParams(reac_vals_before);
    J_ = J_before;
  }

  // diffusion part
  {
    reacordifforcorrho_ = 1;
    double J_before = J_;
    std::cout<<"diffusion gradient according to adjoint analysis"<<std::endl;
    diff_objgrad_->Print(std::cout);

    Epetra_Vector fd_diff_grad(*scatra_discret_->ElementRowMap(),false);
    Teuchos::RCP<Epetra_Vector> perturb_diff_vals = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> diff_vals_before = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),false));
    diff_vals_before->Update(1.0,*diff_vals_,0.0);

    for(int i=0; i<diff_vals_->MyLength(); ++i)
    {
      double perturba = 0.0;//1.0e-3;
      double perturbb = 1.0e-6;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = diff_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_diff_vals->Update(1.0,*diff_vals_before,0.0);
      perturb_diff_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_diff_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_diff_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"diffusion gradient according to FD analysis"<<std::endl;
    fd_diff_grad.Print(std::cout);

    ReplaceParams(diff_vals_before);
    J_ = J_before;
  }

  // sos part
  {
    double J_before = J_;
    std::cout<<"sos gradient according to adjoint analysis"<<std::endl;
    c_objgrad_->Print(std::cout);

    Epetra_Vector fd_c_grad(*acou_discret_->ElementRowMap(),true);
    Teuchos::RCP<Epetra_Vector> perturb_c_vals = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> c_vals_before = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));
    c_vals_before->Update(1.0,*c_vals_,0.0);

    reacordifforcorrho_ = 2;

    for(int i=0; i<c_vals_->MyLength(); ++i)
    {
      if(acou_opt_ind_->operator [](i) == 0.0) continue;

      double perturba =0.0;// 1.0e-3;
      double perturbb = 1.0e-6;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = c_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_c_vals->Update(1.0,*c_vals_before,0.0);
      perturb_c_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_c_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_c_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"sos gradient according to FD analysis"<<std::endl;
    fd_c_grad.Print(std::cout);

    ReplaceParams(c_vals_before);
    J_ = J_before;
  }

  // density part
  {
    double J_before = J_;
    std::cout<<"density gradient according to adjoint analysis"<<std::endl;
    rho_objgrad_->Print(std::cout);

    Epetra_Vector fd_rho_grad(*acou_discret_->ElementRowMap(),true);
    Teuchos::RCP<Epetra_Vector> perturb_rho_vals = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));
    Teuchos::RCP<Epetra_Vector> rho_vals_before = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));
    rho_vals_before->Update(1.0,*rho_vals_,0.0);

    reacordifforcorrho_ = 3;

    for(int i=0; i<rho_vals_->MyLength(); ++i)
    {
      if(acou_opt_ind_->operator [](i) == 0.0) continue;

      double perturba = 0.0;//1.0e-3;
      double perturbb = 1.0e-6;

      double pn=0.0;
      double p=0.0;
      double dp=0.0;

      p = rho_vals_->operator [](i);
      pn = p+p*perturba+perturbb;
      std::cout<<"i "<<i<<" p "<<p<<" disturbed "<<pn<<std::endl;
      perturb_rho_vals->Update(1.0,*rho_vals_before,0.0);
      perturb_rho_vals->ReplaceMyValue(i,0,pn);

      ReplaceParams(perturb_rho_vals);

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();

      dp=(J_before-J_)/(p-pn);
      std::cout<<"J_before - J_ "<<J_before-J_<<" p-pn "<<p-pn<<" val "<<dp<<std::endl;
      fd_rho_grad.ReplaceMyValue(i,0,dp);
    }
    std::cout<<"density gradient according to FD analysis"<<std::endl;
    fd_rho_grad.Print(std::cout);

    ReplaceParams(rho_vals_before);
    J_ = J_before;
  }

  dserror("that is it");

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouSplit::ComputeParameterError()
{
  return;

  if(acou_discret_->Comm().NumProc()>1)
  {
    if(!myrank_)
      std::cout<<"Function ComputeParameterError() skipped! Not yet implemented for parallel usage";
    return;
  }

  // this is implemented problem specific, here for test_recon.dat
  // for a different geometry, you have to implement the correct values here!
  double reac_error = 0.0;
  double diff_error = 0.0;
  double c_error = 0.0;
  double rho_error = 0.0;
  int skip = 0;

  for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
  {
    double reac_val = 0.0; // to be determined
    double D_val = 0.0; // to be determined
    double c_val = 0.0; // to be determined
    double rho_val = 0.0; // to be determined

    // local id
    int lid = scatra_discret_->ElementRowMap()->LID(i);

    // get element center coordinates
    std::vector<double> xyz = DRT::UTILS::ElementCenterRefeCoords(scatra_discret_->lRowElement(lid));

    // check if the center of this element is in the rectangular region
    double g1 = (1.90184747127+0.49231421291)/(3.032389891147+3.54545535887)*(xyz[0]+3.54545535887)-0.49231421291;
    double g2 = (1.90184747127+0.49231421291)/(3.032389891147+3.54545535887)*(xyz[0]+3.0322954880291)-1.90173750140464;
    double g3 = (-1.90173750140464+0.49231421291)/(-3.0322954880291+3.54545535887)*(xyz[0]+3.54545535887)-0.49231421291;
    double g4 = (-1.90173750140464+0.49231421291)/(-3.0322954880291+3.54545535887)*(xyz[0]-3.03238981147)+1.90184747127;

    if(xyz[1]<g1 && xyz[1]>g2 && xyz[1]>g3 && xyz[1]<g4)
    {
      reac_val = 0.1;
      D_val = 0.5;
      c_val = 2.0;
      rho_val = 2.0;
    }
    else
    {
      reac_val = 0.01;
      D_val = 0.2;
      c_val = 1.5;
      rho_val = 1.1;
    }

    // check if the element is close to the boundary of the inclusion

    double h_halbe = 0.23/2.0;
    if( ((xyz[1]<g1+h_halbe)&&(xyz[1]>g1-h_halbe)) ||
        ((xyz[1]<g2+h_halbe)&&(xyz[1]>g2-h_halbe)) ||
        ((xyz[1]<g3+h_halbe)&&(xyz[1]>g3-h_halbe)) ||
        ((xyz[1]<g4+h_halbe)&&(xyz[1]>g4-h_halbe)) )
    {
      skip++;
    }
    else
    {
      reac_error += (reac_vals_->operator [](lid)-reac_val)*(reac_vals_->operator [](lid)-reac_val);
      diff_error += (diff_vals_->operator [](lid)-D_val)*(diff_vals_->operator [](lid)-D_val);
      c_error += (c_vals_->operator [](lid)-c_val)*(c_vals_->operator [](lid)-c_val);
      rho_error += (rho_vals_->operator [](lid)-rho_val)*(rho_vals_->operator [](lid)-rho_val);
    }
    // this is to check if the musterloesung is correctly set
    /*reac_vals_->ReplaceMyValue(lid,0,reac_val);
    diff_vals_->ReplaceMyValue(lid,0,D_val);
    c_vals_->ReplaceMyValue(lid,0,c_val);
    rho_vals_->ReplaceMyValue(lid,0,rho_val);*/
  }

  // this is to check if the musterloesung is correctly set
  /*reacordifforcorrho_ = 0;
  ReplaceParams(reac_vals_);
  reacordifforcorrho_ = 1;
  ReplaceParams(diff_vals_);
  reacordifforcorrho_ = 2;
  ReplaceParams(c_vals_);
  reacordifforcorrho_ = 3;
  ReplaceParams(rho_vals_);

  std::string soutname = name_;
  soutname.append("_reac_and_diff");
  scatraoutput_->NewResultFile(soutname,0);
  output_count_++;
  scatraoutput_->WriteMesh(0,0.0);
  scatraoutput_->NewStep(1,1.0);
  scatraoutput_->WriteElementData(true);
  scatraoutput_->WriteVector("rea_coeff",reac_vals_);
  scatraoutput_->WriteVector("diff_coeff",diff_vals_);dserror("stop");
*/
  std::cout<<std::endl<<"AFTER "<<skip<<" SKIPS: ERRORS IN THE PARAMETER FIELDS (order: mu, D, c, rho): "<<reac_error<<" "<<diff_error<<" "<<c_error<<" "<<rho_error<<std::endl<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionOptiSplitAcouIdent::PatImageReconstructionOptiSplitAcouIdent(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstructionOptiSplit(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout),
  sequenzeiter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("SEQUENZE")),
  acouident_avg_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"ACOUIDENT_AVG"))
{
  c_vals_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));
  rho_vals_ = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),false));

  // fill values vector
  for(int e=0; e<acou_discret_->NumMyRowElements(); ++e)
  {
    DRT::Element* acou_ele = acou_discret_->lRowElement(e);
    c_vals_->ReplaceMyValue(e,0,acou_ele->Material()->Parameter()->GetParameter(1,-1));
    rho_vals_->ReplaceMyValue(e,0,acou_ele->Material()->Parameter()->GetParameter(0,-1));
  }

  // read the material ids
  std::string word2;
  std::istringstream pstream(Teuchos::getNumericStringParameter(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"ACOUPARAMLIST"));
  char* pEnd;
  while (pstream >> word2)
  {
    int id = std::strtol(&word2[0],&pEnd,10);
    if (*pEnd=='\0')
      acou_matids_.push_back(id);
  }

  // read materials
  ReadMaterials(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("SEGMENTATIONMATS"));

}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouIdent::ReplaceParams(Teuchos::RCP<Epetra_Vector> params)
{
  if(reacordifforcorrho_==0)
    for(int i=0; i<params->MyLength(); ++i)
    {
      if(params->operator [](i)<0.0)
        params->ReplaceMyValue(i,0,0.0);
    }
  else if(reacordifforcorrho_==1)
    for(int i=0; i<params->MyLength(); ++i)
    {
      if(params->operator [](i)<0.4)
        params->ReplaceMyValue(i,0,0.4);
      else if(params->operator [](i)>1.0)
        params->ReplaceMyValue(i,0,1.0);
    }
  else
    for(int i=0; i<params->MyLength(); ++i)
    {
      if(params->operator [](i)<0.6)
        params->ReplaceMyValue(i,0,0.6);
      else if(params->operator [](i)>2.0)
        params->ReplaceMyValue(i,0,2.0);
    }

  Teuchos::RCP<Epetra_Vector> paramscol;
  if(reacordifforcorrho_==0||reacordifforcorrho_==1)
    paramscol = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementColMap(),false));
  else
    paramscol = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementColMap(),false));
  LINALG::Export(*params,*paramscol);

  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if(reacordifforcorrho_==0)
  {
    reac_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
      actmat->Parameter()->SetParameter(1,paramscol);
    }
  }
  else if(reacordifforcorrho_==1)
  {
    diff_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<opti_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(opti_matids_[i]);
      actmat->Parameter()->SetParameter(0,paramscol);
    }
  }
  else if(reacordifforcorrho_==2)
  {
    c_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<acou_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(acou_matids_[i]);
      actmat->Parameter()->SetParameter(1,paramscol);
    }
  }
  else if(reacordifforcorrho_==3)
  {
    rho_vals_->Update(1.0,*params,0.0);
    for(unsigned int i=0; i<acou_matids_.size(); ++i)
    {
      Teuchos::RCP<MAT::PAR::Material> actmat = mats.at(acou_matids_[i]);
      actmat->Parameter()->SetParameter(0,paramscol);
    }
  }

  // update node based vector
  if(reacordifforcorrho_==0)
    ComputeNodeBasedReactionCoefficient();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouIdent::SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos)
{
  reacordifforcorrho_ = 0;
  ReplaceParams(reacs);

  reacordifforcorrho_ = 1;
  ReplaceParams(diffs);

  reacordifforcorrho_ = 2;
  ReplaceParams(cs);

  reacordifforcorrho_ = 3;
  ReplaceParams(rhos);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouIdent::ReadMaterials(std::string materialfilename)
{
  // read from the given file
  if (materialfilename=="none.material") dserror("No material file provided");

  // insert path to monitor file if necessary
  if (materialfilename[0]!='/')
  {
    std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    std::string::size_type pos = filename.rfind('/');
    if (pos!=std::string::npos)
    {
      std::string path = filename.substr(0,pos+1);
      materialfilename.insert(materialfilename.begin(), path.begin(), path.end());
    }
  }

  // open the file
  FILE* file = fopen(materialfilename.c_str(),"rb");
  if (file==NULL) dserror("Could not open material file %s",materialfilename.c_str());

  // prepare read in quantities
  char buffer[150000];
  fgets(buffer,150000,file);
  char* foundit = NULL;

  // read number of materials
  foundit = strstr(buffer,"nummats");
  foundit += strlen("nummats");
  nummats_ = strtol(foundit,&foundit,10);

  // prepare the materials
  materialtable_.resize(nummats_);
  for(unsigned i=0; i<nummats_; ++i)
    materialtable_[i].resize(4); // mu_a, D, c, rho

  // read the materials
  foundit = buffer;
  fgets(buffer,150000,file);
  for(unsigned int i=0; i<nummats_; ++i)
  {
    for(int j=0; j<4; ++j)
    {
      materialtable_[i][j] = strtod(foundit,&foundit);
    }
    fgets(buffer,150000,file);
    foundit = buffer;
  }

  std::cout<<"materialtable: "<<std::endl;
  for(unsigned int m=0; m<nummats_; ++m)
    std::cout<<materialtable_[m][0]<<" "<<materialtable_[m][1]<<" "<<materialtable_[m][2]<<" "<<materialtable_[m][3]<<" "<<std::endl;


  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionOptiSplitAcouIdent::PerformIteration()
{
  bool succ = PatImageReconstructionOptiSplit::PerformIteration();
  if(succ==false)
    return succ;

  // update acoustical parameters
  UpdateAcousticalParameters();

  // evaluate everything with the new acoustical parameters
  SolveStandardScatra();
  SolveStandardAcou();
  EvalulateObjectiveFunction();
  SolveAdjointAcou();
  SolveAdjointScatra();
  EvaluateGradient();

  return succ;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouIdent::UpdateAcousticalParameters()
{
  if(!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"ACOUSTICAL UPDATE"<<std::endl;
    std::cout<<std::endl;
  }

  Teuchos::RCP<Epetra_Vector> c_p = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));
  Teuchos::RCP<Epetra_Vector> rho_p = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap(),true));

  // loop the global scatra elements
  for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
  {
    // lid of the element
    int lid = scatra_discret_->ElementRowMap()->LID(i);

    // get the material parameter value
    double loc_D=0.0, loc_reac=0.0;
    if(lid>=0)
    {
      DRT::Element* actele = scatra_discret_->gElement(i);
      // map in GetParameter calculates LID, so we need GID here       05/2017 birzle
      loc_reac = actele->Material()->Parameter()->GetParameter(1,actele->Id());
      loc_D = actele->Material()->Parameter()->GetParameter(0,actele->Id());
    }
    double D = 0.0, reac = 0.0;
    scatra_discret_->Comm().SumAll(&loc_D,&D,1);
    scatra_discret_->Comm().SumAll(&loc_reac,&reac,1);

    // calculate the acoustical values which are required
    // first possibility: closest
    double c = 0.0, rho = 0.0;
    if(!acouident_avg_)
    {
      double abst = 1.0e6;
      int mat = -1;
      for(unsigned m=0; m<nummats_; ++m)
      {
        double abstm = sqrt(0.85*(reac-materialtable_[m][0])*(reac-materialtable_[m][0])+0.15*(D-materialtable_[m][1])*(D-materialtable_[m][1]));
        if(abstm<abst)
        {
          abst = abstm;
          mat = m;
        }
      }
      c = materialtable_[mat][2];
      rho = materialtable_[mat][3];
    }
    else // second possibility: average from two closest
    {
      double mindistance = 1000000.0;
      double mindistancesecond = mindistance;
      int index_mindistance = 0;
      int index_mindistancesecond = 0;

      for(unsigned int m=0; m<nummats_; ++m)
      {
        double currentdistance = sqrt(0.85*(reac-materialtable_[m][0])*(reac-materialtable_[m][0])+0.15*(D-materialtable_[m][1])*(D-materialtable_[m][1]));
        if(currentdistance<mindistance)
        {
          index_mindistancesecond = index_mindistance;
          index_mindistance = m;
          mindistancesecond = mindistance;
          mindistance = currentdistance;
        }
        else if(currentdistance<mindistancesecond)
        {
          mindistancesecond = currentdistance;
          index_mindistancesecond = m;
        }
      }

      double entiredistance = mindistance + mindistancesecond;
      c =   mindistancesecond/entiredistance * materialtable_[index_mindistance][2]
          + mindistance/entiredistance       * materialtable_[index_mindistancesecond][2];
      rho = mindistancesecond/entiredistance * materialtable_[index_mindistance][3]
          + mindistance/entiredistance       * materialtable_[index_mindistancesecond][3];
    }

    // write values to acoustical vector
    int agid = i - scatra_discret_->ElementRowMap()->MinAllGID() + acou_discret_->ElementRowMap()->MinAllGID();
    int alid = acou_discret_->ElementRowMap()->LID(agid);
    c_p->ReplaceMyValue(alid,0,c);
    rho_p->ReplaceMyValue(alid,0,rho);
  }
  reacordifforcorrho_=2;
  ReplaceParams(c_p);
  reacordifforcorrho_=3;
  ReplaceParams(rho_p);
  reacordifforcorrho_=0;

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiSplitAcouIdent::ComputeParameterError()
{
  std::cout<<"here could be an error evaluation "<< std::endl; // TODO
  return;
}

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionReduction::PatImageReconstructionReduction(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstruction(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout),
  reductioncuttime_(acoupara->sublist("PA IMAGE RECONSTRUCTION").get<double>("REDUCTIONCUTTIME"))
{
  if(reductioncuttime_==0.0)
    dserror("if you want to reduce domain and choose patreduction, set REDUCTIONCUTTIME to a reasonable non-zero value, thanks");
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionReduction::Optimize()
{
  // set parameter indicating that not the adjoint problem is solved
  acouparams_->set<bool>("adjoint",false);
  acouparams_->set<bool>("timereversal",true);
  acouparams_->set<bool>("reduction",true);

  // set parameter for acoustic time integration
  acouparams_->set<bool>("acouopt",false);

  Teuchos::RCP<Epetra_MultiVector> tempvec = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors(),true));
  acouparams_->set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",acou_rhsm_);

  // create time integrator
  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_expleuler:
  case INPAR::ACOU::acou_classrk4:
  case INPAR::ACOU::acou_lsrk45reg2:
  case INPAR::ACOU::acou_lsrk33reg2:
  case INPAR::ACOU::acou_lsrk45reg3:
  case INPAR::ACOU::acou_ssprk:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }
  // initialize all quantities to zero
  acoualgo_->SetInitialZeroField();

  // do the time integration
  acoualgo_->Integrate(acou_rhs_);

  // now read the values from the newly written monitor file, return them in time and overwrite the monitor file!
  if(myrank_==0)
  {
    std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
    name.append(".monitor");
    FILE* file = fopen(name.c_str(), "r");
    if (file==NULL) dserror("Could not open monitor file %s",name.c_str());

    char buffer[150000];
    fgets(buffer,150000,file);
    char* foundit = NULL;

    // read steps
    unsigned int nsteps = 0;
    foundit = strstr(buffer,"steps"); foundit += strlen("steps");
    nsteps = strtol(foundit,&foundit,10);
    std::vector<double> timesteps(nsteps);

    // read mics
    unsigned int nmics = 0;
    foundit = strstr(buffer,"mics"); foundit += strlen("mics");
    nmics = strtol(foundit,&foundit,10);

    // read measurement coordinates for every microphone
    std::vector<std::vector<double> > meascoords(nmics);
    for (unsigned int i=0; i<nmics; ++i)
    {
      meascoords[i].resize(3);
      fgets(buffer,150000,file);
      foundit = buffer;
      for(int j=0; j<3; ++j)
        meascoords[i][j] = strtod(foundit,&foundit);
    }

    // read in measured curve
    //{
      Epetra_SerialDenseVector mcurve(nmics*nsteps);

      // read comment lines
      foundit = buffer;
      fgets(buffer,150000,file);
      while(strstr(buffer,"#"))
        fgets(buffer,150000,file);

      // read in the values for each node
      unsigned int count = 0;
      for (unsigned int i=0; i<nsteps; ++i)
      {
        // read the time step
        timesteps[i] = strtod(foundit,&foundit);
        for (unsigned int j=0; j<nmics; ++j)
          mcurve[count++] = strtod(foundit,&foundit);
        fgets(buffer,150000,file);
        foundit = buffer;
      }
      if (count != nmics*nsteps) dserror("Number of measured pressure values wrong on input");
    //}
    fclose(file);

    // open the file to overwrite
    file = fopen(name.c_str(), "w");
    fprintf(file,"steps %d ",nsteps);
    fprintf(file,"mics %d\n",nmics);

    for(unsigned int n=0; n<nmics; ++n)
      fprintf(file,"%e %e %e \n",meascoords[n][0],meascoords[n][1],meascoords[n][2]);
    fprintf(file,"#\n#\n#\n");

    double finaltime = acoualgo_->Time();
    double dt = acoualgo_->TimeStep();
    for(unsigned int i=0; i<nsteps; ++i)
    {
      double time = finaltime-(nsteps-1-i)*dt;
      if(time<1e-3*dt)
        time = 0.0;
      fprintf(file,"%e ",time);
      for(unsigned int m=0; m<nmics; ++m)
        fprintf(file,"%e ",mcurve(m+(nsteps-1-i)*nmics));
      fprintf(file,"\n");
      if(time>reductioncuttime_)
        break;
    }
    fclose(file);

  } // if(myrank_==0)
  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ACOU::PatImageReconstruction::CreateFieldTest()
{
  return Teuchos::rcp(new AcouInvResultTest(*this));
}

/*----------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_MultiVector> ACOU::PatImageReconstruction::ElementMatVec()
{
  return reac_vals_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::ReadMonitor(std::string monitorfilename, double dtacou)
{
  // initialize acou_rhs_: we need a vector with the nodes of the boundary where
  // the pressure is monitored-> Read the monitor file and create a vector with
  // corresponding nodes OR take the boundary where absorbing bcs are prescribed!
  // we need this map extractor thing!
  // we deal with NODES here, not with DOFS

  std::string condname = "PressureMonitor";
  std::vector<DRT::Condition*> pressuremon;
  acou_discret_->GetCondition(condname,pressuremon);
  if(pressuremon.size()==0)
    dserror("you have to use pressure monitor conditions for inverse analysis!");
  const std::vector<int> pressuremonmics = *(pressuremon[0]->Nodes());
  std::vector<int> pressuremonmicsunique;
  std::vector<int> pressuremonmicscolumn;

  nodes_.resize(pressuremonmics.size());
  for(unsigned int i=0; i<pressuremonmics.size(); ++i)
    nodes_[i] = pressuremonmics[i];

  // create unique map
  acou_discret_->Comm().Barrier();
  for(int i=0; i<acou_discret_->Comm().NumProc(); ++i)
  {
    if(acou_discret_->Comm().MyPID() == i)
    {
      for(unsigned int j=0; j<pressuremonmics.size(); ++j)
      {
        if(acou_discret_->HaveGlobalNode(pressuremonmics[j]))
        {
          pressuremonmicscolumn.push_back(pressuremonmics[j]);
          if(acou_discret_->gNode(pressuremonmics[j])->Owner()==int(i))
            pressuremonmicsunique.push_back(pressuremonmics[j]);
        }
      }
    }
    acou_discret_->Comm().Barrier();
  }
  acou_discret_->Comm().Barrier();

  // create map
  abcnodes_map_ = Teuchos::rcp(new Epetra_Map(-1, pressuremonmicsunique.size(), &pressuremonmicsunique[0], 0, acou_discret_->Comm()));
  abcnodes_colmap_ = Teuchos::rcp(new Epetra_Map(-1, pressuremonmicscolumn.size(), &pressuremonmicscolumn[0], 0, acou_discret_->Comm()));

  // determine the number of vectors for monitoring
  int numvec = acouparams_->get<int>("NUMSTEP");
  int oderso = acouparams_->get<double>("MAXTIME")/dtacou;
  if ( oderso < numvec)
    numvec = oderso;

  acou_rhs_ = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,numvec+1,true));
  acou_rhsm_ = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,numvec+1,true));

  unsigned int nsteps = 0;
  unsigned int nmics = 0;

  std::vector<std::vector<double> > meascoords;
  Epetra_SerialDenseVector mcurve;

  // check for monitor file
  if (monitorfilename=="none.monitor")
    dserror("No monitor file provided");

  // insert path to monitor file if necessary
  if (monitorfilename[0]!='/')
  {
    std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
    std::string::size_type pos = filename.rfind('/');
    if (pos!=std::string::npos)
    {
      std::string path = filename.substr(0,pos+1);
      monitorfilename.insert(monitorfilename.begin(), path.begin(), path.end());
    }
  }

  // open monitor file and read it
  FILE* file = fopen(monitorfilename.c_str(),"rb");
  if (file==NULL) dserror("Could not open monitor file %s",monitorfilename.c_str());

  char buffer[150000];
  fgets(buffer,150000,file);
  char* foundit = NULL;

  // read steps
  foundit = strstr(buffer,"steps"); foundit += strlen("steps");
  nsteps = strtol(foundit,&foundit,10);
  std::vector<double> timesteps(nsteps);

  // read mics
  foundit = strstr(buffer,"mics"); foundit += strlen("mics");
  nmics = strtol(foundit,&foundit,10);

  // read measurement coordinates for every microphone
  meascoords.resize(nmics);
  for (unsigned int i=0; i<nmics; ++i)
  {
    meascoords[i].resize(3);
    fgets(buffer,150000,file);
    foundit = buffer;
    for(int j=0; j<3; ++j)
      meascoords[i][j] = strtod(foundit,&foundit);
  }

  // read in measured curve
  {
    mcurve = Epetra_SerialDenseVector(nmics*nsteps);

    // read comment lines
    foundit = buffer;
    fgets(buffer,150000,file);
    while(strstr(buffer,"#"))
      fgets(buffer,150000,file);

    // read in the values for each node
    unsigned int count = 0;
    for (unsigned int i=0; i<nsteps; ++i)
    {
      // read the time step
      timesteps[i] = strtod(foundit,&foundit);
      for (unsigned int j=0; j<nmics; ++j)
        mcurve[count++] = strtod(foundit,&foundit);
      fgets(buffer,150000,file);
      foundit = buffer;
    }
    if (count != nmics*nsteps) dserror("Number of measured pressure values wrong on input");
  }
  fclose(file);

  // interpolation
  // vector for the interpolated data
  Epetra_SerialDenseVector nodcurvinterpol(pressuremonmicsunique.size()*nsteps);

  if((pressuremonmicsunique.size())!=0)
  {
    unsigned int i=0, j=0, l=0;
    unsigned int must_set_curve=1;
    double help;
    double distance[nmics];
    double epsilon = acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("EPSILON");

    //if the user doesn't want to give an espilon as input, we'll calculate it individually
    if(epsilon==-1.0)
      epsilon=ReadMonitorGetEpsilon(pressuremonmicsunique.size());

    //interpolation-loop for every single nod
    for (i=0; i<pressuremonmicsunique.size(); ++i)
    {
      const double* nod_coords = acou_discret_->gNode(pressuremonmicsunique[i])->X();
      unsigned int M_1=0, M_2=0;

      for(j=0;j<nmics;j++)
      {
        distance[j] = ReadMonitorDelta(meascoords[j][0],meascoords[j][1],meascoords[j][2],nod_coords[0],nod_coords[1],nod_coords[2]);

        // if the nod is in an epsilon bubble of any of the microphones, the measured curve of this microphone and the nod's curve should be equal
        if(distance[j]<0.9*epsilon)
        {
          for(l=0;l<nsteps;l++)
            nodcurvinterpol[i*nsteps+l]=mcurve[j+l*nmics];
          must_set_curve=0;
        }
        else
          must_set_curve=1;
      }

      // finds those two microphones, that are the nearest ones to the actual nod
      if(must_set_curve)
      {
        help=distance[0];
        for(j=0;j<nmics;j++)
        {
          if(distance[j]<help)
          {
            help=distance[j];
            M_1=j;
          }
        }
        if((M_1+1)==nmics)
        {
          help=distance[M_1-1];
          M_2=M_1-1;
        }
        else
        {
          help=distance[M_1+1];
          M_2=M_1+1;
        }
        for(j=0;j<nmics;j++)
        {
          if(j==M_1)
            ++j;
          if(distance[j]<help&&j<nmics)
          {
            help=distance[j];
            M_2=j;
          }
        }
        ReadMonitorInterpol(nod_coords,meascoords, M_1,M_2, nmics,i, nsteps, mcurve,nodcurvinterpol);
      }
    }
  }

  double eps = dtacou/1000.0;
  if((numvec-1)*dtacou>timesteps[nsteps-1]+eps) dserror("You want to simulate till %.15f but your monitor file only provides values till %.15f! Fix it!",(numvec-1)*dtacou,timesteps[nsteps-1]);

  // every proc knows mcurve, now, we want to write mcurve to a Epetra_MultiVector in the same form as acou_rhs_
  // with the same parallel distribution!
  // and we want to interpolate measured values in case the monitored time step size is not the same as the one for the simulation
  acou_rhsm_->PutScalar(0.0);

  if( timesteps[0] != 0.0 )
    dserror("your measured values have to start at time 0.0");
  else if( timesteps[0] == 0.0 && timesteps[1] == dtacou ) // the standard case
  {
    for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
      if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
        for(unsigned int j=0; j<nsteps; ++j)
          acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nodcurvinterpol(i*nsteps+j)); // the proc who has this row, writes the value
  }
  else // we have to interpolate!
  {
    if( numvec < int(nsteps) )
    {
      if(  (dtacou/(timesteps[1]-timesteps[0]) - std::ceil(dtacou/(timesteps[1]-timesteps[0]))< 1e-16) && (dtacou/(timesteps[1]-timesteps[0]) - std::ceil(dtacou/(timesteps[1]-timesteps[0]))> -1e-16) ) // dtacou is a multiple of the monitor time step
      {
        int mult = std::ceil(dtacou/(timesteps[1]-timesteps[0]));
        for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
          if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
            for(int j=0; j<numvec; j++)
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nodcurvinterpol(i*nsteps+j*mult)); // the proc who has this row, writes the value
      }
      else
      {
        for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
          if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
          {
            for(int j=0; j<numvec; j++)
            {
              double actualt = j * dtacou; // we need values for this time
              int timeval = 0;
              // find next higher and next lower value
              while(actualt>timesteps[timeval]-eps)
              {
                timeval++;
              }

              // timesteps[timeval] has the next higher point in time
              // now interpolate from this and the value before
              if(timeval == 0)
              {
                acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,0.0);
              }
              else if(actualt<timesteps[timeval]+eps && actualt>timesteps[timeval]-eps) // then this is more or less it
              {
                acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nodcurvinterpol(i*nsteps+timeval));
              }
              else
              {
                double value = nodcurvinterpol(i*nsteps+(timeval-1)) + (nodcurvinterpol(i*nsteps+(timeval))-nodcurvinterpol(i*nsteps+(timeval-1))) * (actualt - timesteps[timeval-1]) / (timesteps[timeval]-timesteps[timeval-1]);
                acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,value);
              }
            }
          }
      }
      //else
      //  dserror("time step bigger than monitor time step but no multiple -> implement here!");
    }
    else
    {
      for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
      {
        if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
        {
          for(int j=0; j<numvec; ++j)
          {
            double actualt = j * dtacou; // we need values for this time
            int timeval = 0;

            // find next higher and next lower value
            while(actualt>timesteps[timeval]-eps)
            {
              timeval++;
            }

            // timesteps[timeval] has the next higher point in time
            // now interpolate from this and the value before
            if(timeval == 0)
            {
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,0.0);
            }
            else if(actualt<timesteps[timeval]+eps && actualt>timesteps[timeval]-eps) // then this is more or less it
            {
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nodcurvinterpol(i*nsteps+timeval));
            }
            else
            {
              double value = nodcurvinterpol(i*nsteps+(timeval-1)) + (nodcurvinterpol(i*nsteps+(timeval))-nodcurvinterpol(i*nsteps+(timeval-1))) * (actualt - timesteps[timeval-1]) / (timesteps[timeval]-timesteps[timeval-1]);
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,value);
            }
          } // for(int j=0; j<numvec; ++j)
        } // if( acou_discret_->HaveGlobalNode(nodes_[i]) )
      } // for(unsigned int i=0; i<nnodes; ++i)
      acou_discret_->Comm().Barrier();
    } // else ** if( numvec < nsteps_ )
  } // else ** if( timesteps[0] == dtacou || (timesteps_[0]==0.0 && timesteps_[1] = dtacou) )

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstruction::ReadMonitorDelta(double coord_M_x,double coord_M_y, double coord_M_z, double coord_N_x,double coord_N_y, double coord_N_z)
{
  double distance = sqrt((coord_M_x-coord_N_x)*(coord_M_x-coord_N_x)+(coord_M_y-coord_N_y)*(coord_M_y-coord_N_y)+(coord_M_z-coord_N_z)*(coord_M_z-coord_N_z));
  return distance;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::ReadMonitorInterpol(const double nod_coords[3],std::vector<std::vector<double> > mic_coords, unsigned int mic_1, unsigned int mic_2, int nmic, unsigned int nod, int timesteps, Epetra_SerialDenseVector& curve, Epetra_SerialDenseVector& inter_curve)
{
  double d1 = ReadMonitorDelta(mic_coords[mic_1][0],mic_coords[mic_1][1],mic_coords[mic_1][2],nod_coords[0],nod_coords[1],nod_coords[2]);
  double d2 = ReadMonitorDelta(mic_coords[mic_2][0],mic_coords[mic_2][1],mic_coords[mic_2][2],nod_coords[0],nod_coords[1],nod_coords[2]);
  double D2 = d2/(d1+d2);
  double D1 = d1/(d1+d2);

  for(int i=0;i<timesteps;i++)
    inter_curve[nod*timesteps+i]=D2*curve[mic_1+i*nmic]+D1*curve[mic_2+i*nmic];

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstruction::ReadMonitorGetEpsilon(int nnodes)
{
  double min_dis[nnodes];
  double dist;
  double min_abs;
  double eps;

  //creates a vector which contains the distance of every nod to its nearest neighbor
  for(int i=0;i<abcnodes_map_->NumMyElements();i++)
  {
    const double* nc= acou_discret_->gNode(abcnodes_map_->GID(i))->X();//acou_discret_->gNode(nodes_[i])->X();
    int iter=0;
    for(int j=0; j<nnodes; j++)
    {
      if(j==i)
        j++;
      if(j==nnodes)
        break;
      const double* ncc=acou_discret_->gNode(abcnodes_map_->GID(j))->X();
      dist = sqrt((nc[0]-ncc[0])*(nc[0]-ncc[0])+(nc[1]-ncc[1])*(nc[1]-ncc[1])+(nc[2]-ncc[2])*(nc[2]-ncc[2]));
      if(iter==0)
        min_dis[i]=dist;
      else if(dist<min_dis[i])
        min_dis[i]=dist;
      ++iter;
    }
  }

  //searches for the (absolute) smallest distance
  min_abs = min_dis[0];
  for(int i=0; i<nnodes; i++)
    if(min_abs>min_dis[i])
      min_abs=min_dis[i];

  eps = min_abs;

  return eps;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::SolveStandardScatra()
{
  // output for user
  scatra_discret_->Comm().Barrier();
  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SCALAR TRANSPORT PROBLEM - OPTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  // create and run scatra algorithm
  const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(*scatraparams_,"VELOCITYFIELD");
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatra_discret_->NumGlobalNodes()==0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      std::string outname = name_;
      outname.append("_invforward_opti");
      if(overwrite_output_)
        scatraoutput_->NewResultFile(outname,0);
      else
        scatraoutput_->NewResultFile(outname,output_count_);
      scatraoutput_->OverwriteResultFile();
      output_count_++;
      scatraoutput_->WriteMesh(0,0.0);

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      scatraalgo_ = Teuchos::rcp(new SCATRA::TimIntStationary(scatra_discret_, scatrasolver_, scatraparams_, scatraextraparams_, scatraoutput_));

      scatraalgo_->Init();
      scatraalgo_->Setup();
      scatraalgo_->SetVelocityField(1);

      scatraalgo_->TimeLoop();

      // output of elemental reaction coefficient
      //OutputReactionAndDiffusion();

      // store the solution vector
      phi_ = scatraalgo_->Phinp();

      break;
    }
    default:
      dserror("unknown velocity field type for transport of passive scalar in problem type Acoustics");
      break;
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::SolveStandardAcou()
{
  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SOUND TRANSPORT PROBLEM - ACOUSTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  // set parameter indicating that the forward problem is solved
  acouparams_->set<bool>("adjoint",false);

  std::string outname = name_;
  outname.append("_invforward_acou");
  if(overwrite_output_)
  {
    acououtput_->NewResultFile(outname,0);
    last_acou_fw_output_count_ = 0;
  }
  else
  {
    acououtput_->NewResultFile(outname,output_count_);
    last_acou_fw_output_count_ = output_count_;
  }
  output_count_++;

  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_expleuler:
  case INPAR::ACOU::acou_classrk4:
  case INPAR::ACOU::acou_lsrk45reg2:
  case INPAR::ACOU::acou_lsrk33reg2:
  case INPAR::ACOU::acou_lsrk45reg3:
  case INPAR::ACOU::acou_ssprk:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }
  acoualgo_->SetInitialPhotoAcousticField(phi_,scatra_discret_,meshconform_);

  // we have to call a slightly changed routine, which fills our history vector which we need for the adjoint problem
  acou_rhs_->Scale(0.0);

  // do the time integration
  acoualgo_->Integrate(acou_rhs_);

  // in case a impulse response of the detectors is given, convolve the simulated pressure curves with this thing
  if(conv_imp_resp_)
  {
    // create temporal vector to hold the convolved values
    Epetra_MultiVector acou_rhs_conv(acou_rhs_->Map(),acou_rhs_->NumVectors());
    for(int n=0; n<acou_rhs_->MyLength(); ++n) // for each node ("detector")
    {
      for(int j=0; j<acou_rhs_->NumVectors(); ++j) // for each time step
      {
        for(int i=0; i<imp_resp_.Length() && i<=j; ++i) // for the length of the impulse response (careful with the first values -> i<=j)
        {
          acou_rhs_conv.SumIntoMyValue(n,j,(acou_rhs_->operator ()(j-i)->operator [](n))*imp_resp_(i));
        }
      }
    }
    acou_rhs_->Update(1.0,acou_rhs_conv,0.0);
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::SolveAdjointAcou()
{
  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SOUND TRANSPORT PROBLEM - ADJOINT ACOUSTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  // set parameter indicating that the adjoint problem is solved
  acouparams_->set<bool>("adjoint",true);

  // set list of monitored nodes
  Teuchos::RCP<std::vector<int> > nodes_rcp= Teuchos::rcp(new std::vector<int> (nodes_.size()));
  for(unsigned int i=0; i<nodes_.size(); ++i)
    (*nodes_rcp)[i] = nodes_[i];
  acouparams_->set<Teuchos::RCP<std::vector<int> > >("monitorednodes",nodes_rcp);
  acouparams_->set<int>("outputcount",last_acou_fw_output_count_);
  acouparams_->set<std::string>("name",name_);

  // build difference vector for adjoint source term
  Teuchos::RCP<Epetra_MultiVector> tempvec = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors(),true));
  tempvec->Update(1.0,*acou_rhs_,0.0);
  tempvec->Update(-1.0,*acou_rhsm_,1.0);

  // in case a impulse response of the detectors is given, do ADJOINT convolution (do this before adjoint mapping with touchcount)
  if(conv_imp_resp_)
  {
    Epetra_MultiVector acou_rhs_conv(acou_rhs_->Map(),acou_rhs_->NumVectors());
    for(int n=0; n<acou_rhs_->MyLength(); ++n) // for each node ("detector")
    {
      for(int j=0; j<acou_rhs_->NumVectors(); ++j) // for each time step
      {
        for(int i=0; i<imp_resp_.Length() && (j+i)<acou_rhs_->NumVectors(); ++i) // for the length of the impulse response (careful with the first values -> i<=j)
        {
          acou_rhs_conv.SumIntoMyValue(n,j,(tempvec->operator ()(j+i)->operator [](n))*imp_resp_(i));
        }
      }
    }
    tempvec->Update(1.0,acou_rhs_conv,0.0);
  }

  // acou_rhs_ has to be scaled with weighting (adjoint of the mapping)
  Teuchos::RCP<Epetra_Vector> touchcountvec = LINALG::CreateVector(*abcnodes_map_);
  acoualgo_->FillTouchCountVec(touchcountvec);
  tempvec->Multiply(1.0,*touchcountvec,*tempvec,0.0);

  // set the difference between measured and simulated values
  Teuchos::RCP<Epetra_MultiVector> tosetvec = Teuchos::rcp(new Epetra_MultiVector(*abcnodes_colmap_,acou_rhsm_->NumVectors(),true));
  LINALG::Export(*tempvec,*tosetvec);
  acouparams_->set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",tosetvec);

  // prepare the output
  std::string outname = name_;
  outname.append("_invadjoint_acou");
  if(overwrite_output_)
    acououtput_->NewResultFile(outname,0);
  else
    acououtput_->NewResultFile(outname,output_count_);
  output_count_++;

  // create the acoustic algorithm
  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_expleuler:
  case INPAR::ACOU::acou_classrk4:
  case INPAR::ACOU::acou_lsrk45reg2:
  case INPAR::ACOU::acou_lsrk33reg2:
  case INPAR::ACOU::acou_lsrk45reg3:
  case INPAR::ACOU::acou_ssprk:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }

  // here the initial field is zero everywhere
  acoualgo_->SetInitialZeroField();

  // integrate the adjoint problem
  acoualgo_->Integrate();

  // give me psi which is needed for the source term of the adjoint optical problem
  adjoint_psi_->PutScalar(0.0);
  acoualgo_->NodalPsiField(adjoint_psi_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::SolveAdjointScatra()
{
  // output for the user
  if(!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "SCALAR TRANSPORT PROBLEM - ADJOINT OPTICAL SYSTEM " << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }

  // get a pointer to the system matrix
  Teuchos::RCP<LINALG::SparseMatrix> sysmatscatra = scatraalgo_->SystemMatrix(); // this matrix and the algorithm should still exist

  // create the right hand side vector for the adjoint optical problem
  Teuchos::RCP<Epetra_Vector> rhsvec = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  rhsvec = CalculateAdjointOptiRhsvec(adjoint_psi_);

  for(int i=0; i<node_reac_->MyLength(); ++i)
  {
    double mu_a = node_reac_->operator [](i);
    int dofgid = scatra_discret_->Dof(0,scatra_discret_->lRowNode(i),0);
    int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
    rhsvec->operator [](doflid) *= -mu_a;
  }

  // perform the element integration
  Teuchos::ParameterList eleparams;
  scatra_discret_->SetState("rhsnodebasedvals",rhsvec);
  eleparams.set<int>("action",SCATRA::calc_integr_pat_rhsvec);
  Teuchos::RCP<Epetra_Vector> b = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  scatra_discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,b,Teuchos::null,Teuchos::null);

  // consider Dirichlet boundaries in the right hand side vector
  std::string condname = "Dirichlet";
  std::vector<DRT::Condition*> dirichlets;
  scatra_discret_->GetCondition(condname,dirichlets);
  for(int nd=0; nd<scatra_discret_->NumMyRowNodes(); ++nd)
  {
    DRT::Node* opti_node = scatra_discret_->lRowNode(nd);
    int nodegid = opti_node->Id();
    for(unsigned int i=0; i<dirichlets.size(); ++i)
    {
      if (dirichlets[i]->ContainsNode(nodegid))
      {
        int dofgid = scatra_discret_->Dof(0,opti_node,0);
        int err = b->ReplaceGlobalValue(dofgid,0,0.0);
        if (err) dserror("could not replace global vector entry");
      }
    }
  }

  // solve the system
  scatrasolver_->Solve(sysmatscatra->EpetraOperator(),adjoint_phi_,b,true,true);

  // output the solution

  std::string outname = name_;
  outname.append("_invadjoint_opti");
  scatraoutput_->NewResultFile(outname,output_count_);
  output_count_++;
  scatraoutput_->WriteMesh(0,0.0);
  scatraoutput_->NewStep(1,1.0);
  scatraoutput_->WriteElementData(true);
  scatraoutput_->WriteVector("phinp",adjoint_phi_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::TimeReversalEstimate()
{
  //{ run the time reversal

  // set parameter indicating that not the adjoint problem is solved
  acouparams_->set<bool>("adjoint",false);
  acouparams_->set<bool>("timereversal",true);

  // initialize output
  std::string outname = name_;
  outname.append("_invforward_acou");
  if(overwrite_output_)
  {
    acououtput_->NewResultFile(outname,0);
    last_acou_fw_output_count_ = 0;
  }
  else
  {
    acououtput_->NewResultFile(outname,output_count_);
    last_acou_fw_output_count_ = output_count_;
  }
  output_count_++;

  // set parameter for acoustic time integration
  acouparams_->set<bool>("acouopt",false);
  acouparams_->set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",acou_rhsm_);

  // create time integrator
  switch(dyna_)
  {
  case INPAR::ACOU::acou_impleuler:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::TimIntImplEuler(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  case INPAR::ACOU::acou_expleuler:
  case INPAR::ACOU::acou_classrk4:
  case INPAR::ACOU::acou_lsrk45reg2:
  case INPAR::ACOU::acou_lsrk33reg2:
  case INPAR::ACOU::acou_lsrk45reg3:
  case INPAR::ACOU::acou_ssprk:
  {
    acoualgo_ = Teuchos::rcp(new ACOU::AcouExplicitTimeInt(acou_discret_,acousolver_,acouparams_,acououtput_));
    break;
  }
  default:
    dserror("Unknown time integration scheme for problem type Acoustics");
    break;
  }
  // initialize all quantities to zero
  acoualgo_->SetInitialZeroField();

  // do the time integration
  acoualgo_->Integrate(acou_rhs_);

  // reset parameter
  acouparams_->set<bool>("timereversal",false);

  //} time reversal run finished

  // now update the optical parameters
  // 1.) solve optical problem with initial guess for absorption coefficient
  // 2.) calculate mu_a as -p_0/Gamma/phi
  // 3.) bring these values to the parameter vector
  for(int i=0;i<10; ++i)
  {
    std::cout<<"TR ITERATION "<<i<<std::endl;

    // do step 1.
    SolveStandardScatra(); // phi now holds the optical solution values

    // do step 2. and 3.
    UpdateAbsorptionCoefficientFromTimeReversal();
  }
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::UpdateAbsorptionCoefficientFromTimeReversal()
{
  // we need a parameter list for the acoustical element evaluation
  Teuchos::ParameterList para;
  para.set<int>("action",ACOU::calc_average_pressure);
  para.set<bool>("padaptivity",false);
  para.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  para.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  para.set<bool>("mesh conform",meshconform_);
  para.set<int>("useacouoptvecs",-1);
  DRT::Element::LocationArray la(2);
  Epetra_SerialDenseVector elevec(1);
  Epetra_SerialDenseMatrix elemat;

  // the vector which we fill with absorption values
  Teuchos::RCP<Epetra_Vector> trparams = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap(),true));

  // do the business
  if(meshconform_)
  {
    for(int e=scatra_discret_->ElementRowMap()->MinAllGID(); e<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++e)
    {
      // find the owner of the optical element
      int myopteleowner = -1;
      int opteleowner = -1;
      DRT::Element* opti_ele = NULL;
      if(scatra_discret_->HaveGlobalElement(e))
      {
        opti_ele = scatra_discret_->gElement(e);
        myopteleowner = opti_ele->Owner();
        if(myopteleowner!=scatra_discret_->Comm().MyPID())
          myopteleowner = -1;
      }
      scatra_discret_->Comm().MaxAll(&myopteleowner,&opteleowner,1);

      // find the owner of the acoustical element
      int myacoueleowner = -1;
      int acoueleowner = -1;
      DRT::Element* acou_ele = NULL;
      if(acou_discret_->HaveGlobalElement(e-scatra_discret_->ElementRowMap()->MinAllGID()+acou_discret_->ElementRowMap()->MinAllGID()))
      {
        acou_ele = acou_discret_->gElement(e-scatra_discret_->ElementRowMap()->MinAllGID()+acou_discret_->ElementRowMap()->MinAllGID());
        myacoueleowner = acou_ele->Owner();
        if(myacoueleowner!=myrank_)
          myacoueleowner = -1;
      }
      acou_discret_->Comm().MaxAll(&myacoueleowner,&acoueleowner,1);

      if(acoueleowner == opteleowner)
      {
        // the owning processor can do all his business
        if(opteleowner==myrank_)
        {
          // get grueneisen
          double gamma = 1.0;

          // get average light flux from solution vector phi_
          double phi = 0.0;
          for(int i=0; i<opti_ele->NumNode(); ++i)
            phi += phi_->operator [](scatra_discret_->DofRowMap()->LID(scatra_discret_->Dof((opti_ele->Nodes()[i]),0)));
          phi /= double(opti_ele->NumNode());

          // get average pressure value from the acoustical element
          acou_ele->LocationVector(*acou_discret_,la,false);
          acou_ele->Evaluate(para,*acou_discret_,la[0].lm_,elemat,elemat,elevec,elevec,elevec);
          double pressure = elevec[0];

          // compute absorption coefficient
          double reac = -pressure/gamma/phi;
          if(reac<0.0)
            reac=0.0;

          // write absorption coefficient to parameter vector
          trparams->ReplaceMyValue(scatra_discret_->ElementColMap()->LID(e),0,reac);

        }
        // the other processors do not have to do anything
      }
      else
      {
        // optical and acoustical element are not owned by the same processor -> communicate acoustical values
        // optical owner does most of the business

        // get average pressure
        double locpress = 0.0;
        double pressure = 0.0;
        if(acoueleowner==myrank_)
        {
          acou_ele->Evaluate(para,*acou_discret_,la[0].lm_,elemat,elemat,elevec,elevec,elevec);
          pressure = elevec[0];
        }
        acou_discret_->Comm().SumAll(&locpress,&pressure,1);

        // compute absorption coefficient
        if(opteleowner==myrank_)
        {
          // get grueneisen
          double gamma = 1.0;

          // get average light flux from solution vector phi_
          double phi = 0.0;
          for(int i=0; i<opti_ele->NumNode(); ++i)
            phi += phi_->operator [](scatra_discret_->DofRowMap()->LID(scatra_discret_->Dof((opti_ele->Nodes()[i]),0)));
          phi /= double(opti_ele->NumNode());

          double reac = -pressure/gamma/phi;
          if(reac<0.0)
            reac=0.0;

          // write absorption coefficient to parameter vector
          trparams->ReplaceMyValue(scatra_discret_->ElementColMap()->LID(e),0,reac);
        }
      }
    }
  }
  else
    dserror("update of absorption coefficient not yet implemented for nonconforming mesh");

  // bring values to the elements
  ReplaceParams(trparams);

  // update the node based vector
  ComputeNodeBasedReactionCoefficient();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::EvaluateError()
{
  // contribution from difference between measured and simulated values

  // build difference vector
  Epetra_MultiVector tempvec = Epetra_MultiVector(*abcnodes_map_,acou_rhsm_->NumVectors());
  tempvec.Update(1.0,*acou_rhsm_,0.0);
  tempvec.Update(1.0,*acou_rhs_,-1.0);

  // take the square
  tempvec.Multiply(1.0,tempvec,tempvec,0.0);

  // build the norm of each vector
  Epetra_SerialDenseVector normvec(acou_rhsm_->NumVectors());
  tempvec.Norm1(normvec.Values());

  // sum all norms and do not forget factor 0.5
  error_ = 0.5 * normvec.Norm1();

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::OutputReactionAndDiffusion()
{

  // build the two vectors
  Teuchos::RCP<Epetra_Vector> reacvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  Teuchos::RCP<Epetra_Vector> diffvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);
    // map in GetParameter calculates LID, so we need GID here       05/2017 birzle
    double reac = actele->Material()->Parameter()->GetParameter(1,actele->Id());
    double diff = actele->Material()->Parameter()->GetParameter(0,actele->Id());
    reacvec->operator [](i) = reac;
    diffvec->operator [](i) = diff;
  }
  scatraoutput_->WriteVector("rea_coeff",reacvec);
  scatraoutput_->WriteVector("diff_coeff",diffvec);
  scatraoutput_->WriteInt("iteration",iter_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::ComputeNodeBasedReactionCoefficient()
{
  // fill the vectors
  int minnodeidscatra = scatra_discret_->NodeRowMap()->MinAllGID();

  for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd) // cannot loop scatra nodes, because they don't necessarily start with gid 0
  {
    // get node and owner
    int myoptnodeowner = -1;
    int optnodeowner = -1;
    DRT::Node* opti_node = NULL;
    if( scatra_discret_->HaveGlobalNode(nd+minnodeidscatra) )
    {
      opti_node = scatra_discret_->gNode(nd+minnodeidscatra);
      myoptnodeowner = opti_node->Owner();
      if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
    }
    scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
    if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
      continue;
    // here, every proc knows the owner and the gid of the optical node

    // now, every proc knows the speed of sound and density of this node -> write them to vector
    int nodelid = scatra_discret_->NodeRowMap()->LID(nd+minnodeidscatra);

    // we have to do the same procedure for the absorption coefficient with the scatra discretization
    int loc_numoptiele = 0;
    double loc_mu_a = 0.0;
    for(int roel = 0; roel<scatra_discret_->NumMyRowElements(); ++roel)
    {
      DRT::Element* roptele = scatra_discret_->lRowElement(roel);
      const int* nodeids = roptele->NodeIds();
      int numnode = roptele->NumNode();
      for(int i=0; i<numnode; ++i)
        if( nodeids[i] == nd+minnodeidscatra )
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(roptele->Material().get());
          loc_mu_a += actmat->ReaCoeff(roptele->Id());
          loc_numoptiele++;
        }
    }
    int glo_numoptiele = 0;
    double glo_mu_a = 0.0;
    scatra_discret_->Comm().SumAll(&loc_numoptiele,&glo_numoptiele,1);
    scatra_discret_->Comm().SumAll(&loc_mu_a,&glo_mu_a,1);
    glo_mu_a /= double(glo_numoptiele);

    // write this value to vect
    if( nodelid >= 0 ) // only on owning proc
      node_reac_->ReplaceMyValue(nodelid,0,glo_mu_a);
  } // for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ACOU::PatImageReconstruction::CalculateAdjointOptiRhsvec(Teuchos::RCP<Epetra_Vector> acounodevec)
{
  // this function is similar to the mapping in void ACOU::AcouImplicitTimeInt::SetInitialPhotoAcousticField
  // just the other way round
  Teuchos::RCP<Epetra_Vector> rhsvec = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);

  if(meshconform_)
  {
    int minscatranodegid = scatra_discret_->NodeRowMap()->MinAllGID();
    int minacounodegid = acou_discret_->NodeRowMap()->MinAllGID();
    for(int nd=0; nd<acou_discret_->NumGlobalNodes(); ++nd)
    {
      // get node and owner
      int myoptnodeowner = -1;
      int optnodeowner = -1;
      DRT::Node* opti_node = NULL;
      if( scatra_discret_->HaveGlobalNode(nd+minscatranodegid) )
      {
        opti_node = scatra_discret_->gNode(nd+minscatranodegid);
        myoptnodeowner = opti_node->Owner();
        if( myoptnodeowner != scatra_discret_->Comm().MyPID() ) myoptnodeowner = -1; // cannot use myrank_ because that is acou_discret_->Comm().MyPID()
      }
      scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);
      if( optnodeowner == -1 ) // in this case, this node does not exist in the scatra discretization
        continue;

      double loc_value = 0.0;
      if(acou_discret_->NodeRowMap()->LID(nd+minacounodegid)>-1)
      {
        loc_value = adjoint_psi_->operator [](acou_discret_->NodeRowMap()->LID(nd+minacounodegid));
      }
      double glo_value = 0.0;
      acou_discret_->Comm().SumAll(&loc_value,&glo_value,1);

      // ok, we got the value, we still need c, rho and mu_a, but they are stored on the nodemap of the scatra dis, so this should not be a problem
      if(scatra_discret_->Comm().MyPID() == optnodeowner)
      {
        int dofgid = scatra_discret_->Dof(0,opti_node,0);
        int doflid = scatra_discret_->DofRowMap()->LID(dofgid);
        int err = rhsvec->ReplaceMyValue(doflid,0,glo_value);
        if (err) dserror("could not replace local vector entry");
      }
    }
  }
  else
  {
    // export input vector to column map
    Teuchos::RCP<Epetra_Vector> acounodeveccol = Teuchos::rcp(new Epetra_Vector(*(acou_discret_->NodeColMap())),true);
    LINALG::Export(*acounodevec,*acounodeveccol);

    int numdim = DRT::Problem::Instance()->NDim();
    int minoptnodegid = scatra_discret_->NodeRowMap()->MinAllGID();
    for(int optnd=0; optnd<scatra_discret_->NumGlobalNodes(); ++optnd)
    {
      DRT::Node* optnode = NULL;
      int myoptnodeowner = -1;
      if(scatra_discret_->HaveGlobalNode(optnd+minoptnodegid))
      {
        optnode = scatra_discret_->gNode(optnd+minoptnodegid);
        myoptnodeowner = optnode->Owner();
        if(myoptnodeowner != myrank_) myoptnodeowner = -1;
      }
      int optnodeowner = -1;
      scatra_discret_->Comm().MaxAll(&myoptnodeowner,&optnodeowner,1);

      double optnodecoords[numdim];
      if(myrank_==optnodeowner)
      {
        for(int d=0; d<numdim; ++d)
          optnodecoords[d] = optnode->X()[d];
      }
      scatra_discret_->Comm().Broadcast(&optnodecoords[0],numdim,optnodeowner);

      double r = 0.0;
      for(int acouel=0; acouel<acou_discret_->NumMyRowElements(); ++acouel)
      {
        DRT::Element* ele = acou_discret_->lRowElement(acouel);
        // get the nodes of this element, and then check if acoustical node is inside
        if(ele->Shape()==DRT::Element::quad4)
        {
          double acounodecoords[4][numdim];
          double minmaxvals[2][numdim];
          for(int j=0; j<numdim; ++j)
          {
            minmaxvals[0][j] = 1.0e6; // minvals
            minmaxvals[1][j] = -1.0e6; // maxvals
          }
          for(int nd=0;nd<4;++nd) // quad4 has 4 nodes
            for(int d=0;d<numdim;++d)
            {
              acounodecoords[nd][d] = ele->Nodes()[nd]->X()[d];
              if(acounodecoords[nd][d] < minmaxvals[0][d]) minmaxvals[0][d]=acounodecoords[nd][d];
              if(acounodecoords[nd][d] > minmaxvals[1][d]) minmaxvals[1][d]=acounodecoords[nd][d];
            }
          // check, if acoustical node is in bounding box
          bool inside = true;
          for(int d=0;d<numdim;++d)
            if(optnodecoords[d]>minmaxvals[1][d]+5.0e-5 || optnodecoords[d]<minmaxvals[0][d]-5.0e-5)
              inside=false;
          if(inside)
          {
            // solve for xi by local Newton
            LINALG::Matrix<2,1> F(true);
            LINALG::Matrix<2,2> dFdxi(true);
            LINALG::Matrix<2,1> xi(true);
            LINALG::Matrix<2,1> deltaxi(true);
            double deltaxinorm = 0.0;
            int count = 0;
            do{
              count++;
              F(0) = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * acounodecoords[0][0]
                   + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * acounodecoords[1][0]
                   + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * acounodecoords[2][0]
                   + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * acounodecoords[3][0]  - optnodecoords[0];
              F(1) = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * acounodecoords[0][1]
                   + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * acounodecoords[1][1]
                   + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * acounodecoords[2][1]
                   + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * acounodecoords[3][1]  - optnodecoords[1] ;

              dFdxi(0,0) = - 0.25 * (1. - xi(1)) * acounodecoords[0][0]
                           + 0.25 * (1. - xi(1)) * acounodecoords[1][0]
                           + 0.25 * (1. + xi(1)) * acounodecoords[2][0]
                           - 0.25 * (1. + xi(1)) * acounodecoords[3][0] ;
              dFdxi(0,1) = - 0.25 * (1. - xi(0)) * acounodecoords[0][0]
                           - 0.25 * (1. + xi(0)) * acounodecoords[1][0]
                           + 0.25 * (1. + xi(0)) * acounodecoords[2][0]
                           + 0.25 * (1. - xi(0)) * acounodecoords[3][0] ;
              dFdxi(1,0) = - 0.25 * (1. - xi(1)) * acounodecoords[0][1]
                           + 0.25 * (1. - xi(1)) * acounodecoords[1][1]
                           + 0.25 * (1. + xi(1)) * acounodecoords[2][1]
                           - 0.25 * (1. + xi(1)) * acounodecoords[3][1] ;
              dFdxi(1,1) = - 0.25 * (1. - xi(1)) * acounodecoords[0][1]
                           - 0.25 * (1. + xi(1)) * acounodecoords[1][1]
                           + 0.25 * (1. + xi(1)) * acounodecoords[2][1]
                           + 0.25 * (1. - xi(1)) * acounodecoords[3][1] ;

              LINALG::FixedSizeSerialDenseSolver<2,2,1> inverser;
              inverser.SetMatrix(dFdxi);
              inverser.SetVectors(deltaxi,F);
              inverser.Solve();

              deltaxinorm = deltaxi.Norm2();
              xi.Update(-1.0,deltaxi,1.0);
            } while ( deltaxinorm > 1.0e-8 && count < 10 );
            if(!(count == 10 || xi.NormInf()>1.0+0.15))
            {
              // get the values!
              double values[4] = {0};
              for(int nd=0;nd<4;++nd)
              {
                int lid = acou_discret_->NodeColMap()->LID(ele->Nodes()[nd]->Id());
                if(lid<0)
                  dserror("node of element not on this processor");
                else
                  values[nd] = acounodeveccol->operator [](lid);
              }
              r = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * values[0]
                + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * values[1]
                + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * values[2]
                + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * values[3];
            }
          } // if(inside)
        }
        else dserror("up to now only implemented for quad4");
      } // for(int acouel=0; acouel<acou_discret_->NumMyRowElements(); ++acouel)

      // one processor might provide a value
      double glob_p_min = 0.0;
      scatra_discret_->Comm().MinAll(&r,&glob_p_min,1);
      double glob_p_max = 0.0;
      scatra_discret_->Comm().MaxAll(&r,&glob_p_max,1);
      // take higher absolute values
      double glob_p = 0.0;
      if(std::abs(glob_p_min)>std::abs(glob_p_max))
        glob_p = glob_p_min;
      else
        glob_p = glob_p_max;

      // set p value in node based vector
      if(myrank_==optnodeowner && glob_p != 0.0)
      {
        int dof = scatra_discret_->Dof(0,optnode,0);
        int lid = scatra_discret_->DofRowMap()->LID(dof);
        if(lid<0) dserror("cannot find dof for node %d ",optnd);

        int err = rhsvec->ReplaceMyValue(lid,0,glob_p);
        if (err) dserror("could not replace local vector entry");
      }
    } // for(int optnd=0; optnd<scatra_discret_->NumGlobalNodes(); ++optnd)
  }
  return rhsvec;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::OutputStats()
{
  if (!myrank_)
  {
    std::cout<<std::endl;
    std::cout<<"*** objective function value:             "<<J_<<std::endl;
    std::cout<<"*** relative objective function value:    "<<J_/J_start_<<std::endl;
    std::cout<<"*** error value:                          "<<error_<<std::endl;
    std::cout<<"*** relative error value:                 "<<error_/error_start_<<std::endl;
    std::cout<<"*** output count:                         "<<output_count_<<std::endl;
    std::cout<<"*** simulation time since start [h]:      "<<(Teuchos::Time::wallTime()-tstart_)/(60.0*60.0)<<std::endl;
    std::cout<<"*** parameters:                           "<<std::endl;
  }
  return;
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::ReadRestart(int restartoutputcount)
{
  output_count_ = restartoutputcount;

  // first step is to get rid of the "-1" which one automatically gets from restart
  name_.erase(name_.size()-2,2);

  // create file name
  std::string scainputfilename;
  {
    std::ostringstream temp;
    if(overwrite_output_)
      temp<<name_<<"_invforward_opti_run_"<<0;
    else
      temp<<name_<<"_invforward_opti_run_"<<restartoutputcount;
    scainputfilename = temp.str();
  }
  std::string acouinputfilename;
  {
    std::ostringstream temp;
    if(overwrite_output_)
      temp<<name_<<"_invforward_acou_run_"<<0;
    else
      temp<<name_<<"_invforward_acou_run_"<<restartoutputcount+1;
    acouinputfilename = temp.str();
  }

  if(!myrank_)
  {
    std::cout<<"READING RESTART FROM OUTPUTCOUNT "<<restartoutputcount<<" FROM FILES:"<<std::endl;
    std::cout<<scainputfilename<<std::endl;
    std::cout<<acouinputfilename<<std::endl;
    std::cout<<std::endl;
  }

  Teuchos::RCP<IO::InputControl> scainputfile = Teuchos::rcp(new IO::InputControl(scainputfilename,scatra_discret_->Comm()));
  Teuchos::RCP<IO::InputControl> acouinputfile = Teuchos::rcp(new IO::InputControl(acouinputfilename,acou_discret_->Comm()));

  IO::DiscretizationReader scareader(scatra_discret_,scainputfile,1);
  IO::DiscretizationReader acoureader(acou_discret_,acouinputfile,0);

  iter_ = scareader.ReadInt("iteration");

  // vectors
  Teuchos::RCP<Epetra_Vector> reacs = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap()));
  Teuchos::RCP<Epetra_Vector> diffs = Teuchos::rcp(new Epetra_Vector(*scatra_discret_->ElementRowMap()));
  Teuchos::RCP<Epetra_Vector> cs = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap()));
  Teuchos::RCP<Epetra_Vector> rhos = Teuchos::rcp(new Epetra_Vector(*acou_discret_->ElementRowMap()));

  // read all vectors
  scareader.ReadVector(reacs,"rea_coeff");
  scareader.ReadVector(diffs,"diff_coeff");
  acoureader.ReadVector(rhos,"density");
  acoureader.ReadVector(cs,"speedofsound");

  SetRestartParameters(reacs,diffs,cs,rhos);

  return;
}


/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos)
{
  //reac_vals_->Update(1.0,*reacs,0.0);
  ReplaceParams(reacs);
  return;
}

