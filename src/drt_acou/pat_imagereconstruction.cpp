/*!----------------------------------------------------------------------
\file pat_imagereconstruction.cpp

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "pat_imagereconstruction.H"
#include "pat_matpar_manager.H"
#include "pat_utils.H"
#include "acou_expl.H"
#include "acou_impl_euler.H"
#include "acou_inv_resulttest.H"
#include "acou_ele.H"
#include "acou_ele_action.H"
#include "acou_sol_ele.H"

#include "../drt_inv_analysis/invana_utils.H"
#include "../drt_inv_analysis/regularization_base.H"
#include "../drt_inv_analysis/regularization_tikhonov.H"
#include "../drt_inv_analysis/regularization_totalvariation.H"
#include "../drt_inv_analysis/regularization_tvdtikh.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/scatra_mat.H"
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
J_(0.0),
J_start_(0.0),
error_(0.0),
error_start_(0.0)
{
  // set time reversal to false
  acouparams_->set<bool>("timereversal",false);

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

  // setup the search direction handler
  searchdirection_ = Teuchos::rcp(new PATSearchDirection(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")));

  // setup the line search
  linesearch_ = Teuchos::rcp(new PATLineSearch(Teuchos::rcp(this,false)));

  // create material manager
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(invp,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
      scatra_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerPerElement(scatra_discret_));
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
      scatra_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerUniform(scatra_discret_));
    break;
    default:
      dserror("choose a valid method of parametrization for the material parameter field");
    break;
  }
  scatra_matman_->Init(invp);
  scatra_matman_->Setup();

  // read monitor file, create multivector and map for measured values
  ReadMonitor(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("MONITORFILE"),acouparams_->get<double>("TIMESTEP"));

  // compute node based reaction vector
  ComputeNodeBasedReactionCoefficient();
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
void ACOU::PatImageReconstruction::SampleObjectiveFunction()
{
  // reaction coefficient
  double firstcircle = 5.0;
  double secondcircle = 0.12;
  double rect = 0.1;
  double soft = 0.3;
  double def = 0.1;

  // diffusion coefficient to default
  //double D = 0.1;
  double reac = 0.004;

  int rmax = 20;

  for(int r=0; r<rmax; ++r)
  {
    double ratiocorrect = double(r+1)/double(rmax);
    std::cout<<"run "<<r<<" ratiocorrect "<<ratiocorrect<<std::endl;

    Teuchos::RCP<Epetra_MultiVector> tatparams = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));

    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {
      double val = 0.0; // to be determined

      // local id
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get element center coordinates
      std::vector<double> xyz = scatra_discret_->lRowElement(lid)->ElementCenterRefeCoords();

      // check first circle:
      double p = sqrt((xyz[0]-4.)*(xyz[0]-4.)+(xyz[1]-3.)*(xyz[1]-3.));
      if(p<2.0)
        val = firstcircle;
      else
      {
        // check second circle
        p = sqrt((xyz[0]+4.)*(xyz[0]+4.)+(xyz[1]-3.)*(xyz[1]-3.));
        if(p<2.5)
          val = secondcircle;
        else
        {
          // check rectangle
          double g1 = 0.176327*(xyz[0]+2.52030281)-2.98296392;
          double g2 = 0.176327*(xyz[0]+1.99935828)-5.93738717;
          double g3 = -5.671281835*(xyz[0]+1.99935828)-5.93738717;
          double g4 = -5.671281835*(xyz[0]-5.87910374)-4.54820175;

          if(xyz[1]<g1 && xyz[1]>g2 && xyz[1]>g3 && xyz[1]<g4)
            val = rect;
          else if( sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])<10. ) // check soft tissue circle
            val = soft;
          else
            val = def;
        }
      }
      val = ratiocorrect * val + (1.0-ratiocorrect)*def;

      int elematid = scatra_discret_->lRowElement(lid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::OptMatParManagerPerElement>(scatra_matman_)->WriteValuesToVector(elematid,i,val,reac,tatparams);
    }
    // bring the parameters to the materials
    scatra_matman_->ReplaceParams(*tatparams);

    SolveStandardScatra();
    SolveStandardAcou();
    EvalulateObjectiveFunction();

  }
  dserror("that is it");
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::InitialRun()
{
  // determine if we have to do the forward run or everything is zero
  double maxval = 0.0;
  scatra_matman_->GetParams()->MaxValue(&maxval);
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

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::InitialRun()
{
  InitConcentrations();
  ACOU::PatImageReconstruction::InitialRun();

  return;
}
/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionOpti::PatImageReconstructionOpti(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstruction(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout)
{
  // get inverse analysis parameters
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // create regularization manager
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp,"REGULARIZATION"))
  {
  case INPAR::INVANA::stat_inv_reg_none:
    break;
  case INPAR::INVANA::stat_inv_reg_tikhonov:
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
  break;
  case INPAR::INVANA::stat_inv_reg_totalvariation:
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
  break;
  case INPAR::INVANA::stat_inv_reg_tvdtikh:
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
  break;
  default:
    dserror("no valid regularization type provided");
  break;
  }
  if (scatra_regman_!=Teuchos::null)
  {
    scatra_regman_->Init(scatra_discret_,scatra_matman_->GetConnectivityData());
    scatra_regman_->Setup(invp);
  }

  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",false);

  // allocate gradient
  scatra_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));

  // setup direction
  searchdirection_->Setup(scatra_matman_->ParamLayoutMap(),scatra_matman_->ParamLayoutMapUnique(),1);
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOpti::ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params)
{
  // bring the parameters to the materials
  scatra_matman_->ReplaceParams(*params);

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstructionOpti::EvalulateObjectiveFunction()
{
  // evaluate error contribution
  EvaluateError();
  J_ = error_;

  // evaluate optical regularization
  if(scatra_regman_ != Teuchos::null)
    scatra_regman_->Evaluate(*(scatra_matman_->GetParams()),&J_);

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<J_-error_<<std::endl;

  return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOpti::EvaluateGradient()
{
  // zero out gradient vector initially
  scatra_objgrad_->Scale(0.0);

  // contribution form regularization
  if(scatra_regman_ != Teuchos::null)
    scatra_regman_->EvaluateGradient(*(scatra_matman_->GetParams()),scatra_objgrad_);

  // set quantities needed by the elements
  scatra_discret_->SetState("adjoint phi",adjoint_phi_);
  scatra_discret_->SetState("phi",phi_);

  // fill and set psi vector
  Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
  psi = CalculateAdjointOptiRhsvec(adjoint_psi_);
  scatra_discret_->SetState("psi",psi);

  // do the actual evaluation
  scatra_matman_->AddEvaluate(0.0,scatra_objgrad_);
  scatra_matman_->Finalize(scatra_objgrad_);

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionOpti::PerformIteration()
{
  linesearch_->Init(J_,scatra_objgrad_,searchdirection_->ComputeDirection(scatra_objgrad_,scatra_matman_->GetParams(),iter_),scatra_matman_->GetParams(),scatra_matman_->ParamLayoutMapUnique());
  return linesearch_->Run();
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOpti::CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result)
{
  INVANA::MVDotProduct(*scatra_objgrad_,bvector,uniquemap,result);
  return;
}

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionOptiAcou::PatImageReconstructionOptiAcou(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstructionOpti(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout)
{
  // create material manager
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  Teuchos::ParameterList list(invp);
  list.set<std::string>("PARAMLIST",acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("ACOUPARAMLIST"));
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(list,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
      acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerPerElement(acou_discret_));
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
      acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerUniform(acou_discret_));
    break;
    default:
      dserror("choose a valid method of parametrization for the material parameter field");
    break;
  }
  acou_matman_->Init(list);
  acou_matman_->Setup();

  // create regularization manager
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp,"REGULARIZATION"))
  {
  case INPAR::INVANA::stat_inv_reg_none:
    break;
  case INPAR::INVANA::stat_inv_reg_tikhonov:
    acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
  break;
  case INPAR::INVANA::stat_inv_reg_totalvariation:
    acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
  break;
  case INPAR::INVANA::stat_inv_reg_tvdtikh:
    acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
  break;
  default:
    dserror("no valid regularization type provided");
  break;
  }
  if (acou_regman_!=Teuchos::null)
  {
    acou_regman_->Init(acou_discret_,acou_matman_->GetConnectivityData());
    acou_regman_->Setup(invp);
  }

  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",true);

  // allocate gradient
  acou_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));

  // setup direction
  acou_searchdirection_ = Teuchos::rcp(new PATSearchDirection(DRT::INPUT::IntegralValue<INPAR::ACOU::OptimizationType>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"OPTIMIZATION")));
  acou_searchdirection_->Setup(acou_matman_->ParamLayoutMap(),acou_matman_->ParamLayoutMapUnique(),1);
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiAcou::SampleObjectiveFunction()
{
  // rho coefficient
  double reac_firstcircle = 0.32;
  double reac_secondcircle = 0.32;
  double reac_rect = 0.2;
  double reac_soft = 0.02;
  double reac_def = 0.02; // 0.004;
  double D_firstcircle = 0.12;
  double D_secondcircle = 0.12;
  double D_rect = 0.2;
  double D_soft = 0.3;
  double D_def = 0.3; //0.1;
  double c_firstcircle = 1.8;
  double c_secondcircle = 1.8;
  double c_rect = 0.7;
  double c_soft = 1.6;
  double c_def = 1.6;//1.48;
  double rho_firstcircle = 1.5;
  double rho_secondcircle = 1.5;
  double rho_rect = 1.0;
  double rho_soft = 1.2;
  double rho_def = 1.2; //1.0;

  int rmax = 20;

  for(int r=0; r<=rmax; ++r)
  {
    double ratiocorrect = double(r)/double(rmax);
    std::cout<<"run "<<r<<" ratiocorrect "<<ratiocorrect<<std::endl;

    Teuchos::RCP<Epetra_MultiVector> sca_tatparams = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
    Teuchos::RCP<Epetra_MultiVector> acou_tatparams = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));

    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {
      double reac_val = 0.0; // to be determined
      double D_val = 0.0; // to be determined
      double c_val = 0.0; // to be determined
      double rho_val = 0.0; // to be determined

      // local id
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get element center coordinates
      std::vector<double> xyz = scatra_discret_->lRowElement(lid)->ElementCenterRefeCoords();

      // check first circle:
      double p = sqrt((xyz[0]-4.)*(xyz[0]-4.)+(xyz[1]-3.)*(xyz[1]-3.));
      if(p<2.0)
      {
        reac_val = reac_firstcircle;
        D_val = D_firstcircle;
        c_val = c_firstcircle;
        rho_val = rho_firstcircle;
      }
      else
      {
        // check second circle
        p = sqrt((xyz[0]+4.)*(xyz[0]+4.)+(xyz[1]-3.)*(xyz[1]-3.));
        if(p<2.5)
        {
          reac_val = reac_secondcircle;
          D_val = D_secondcircle;
          c_val = c_secondcircle;
          rho_val = rho_secondcircle;
        }
        else
        {
          // check rectangle
          double g1 = 0.176327*(xyz[0]+2.52030281)-2.98296392;
          double g2 = 0.176327*(xyz[0]+1.99935828)-5.93738717;
          double g3 = -5.671281835*(xyz[0]+1.99935828)-5.93738717;
          double g4 = -5.671281835*(xyz[0]-5.87910374)-4.54820175;

          if(xyz[1]<g1 && xyz[1]>g2 && xyz[1]>g3 && xyz[1]<g4)
          {
            reac_val = reac_rect;
            D_val = D_rect;
            c_val = c_rect;
            rho_val = rho_rect;
          }
          else if( sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])<10. ) // check soft tissue circle
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
      reac_val = ratiocorrect * reac_val + (1.0-ratiocorrect)*reac_def;
      //D_val = D_def;
      D_val = ratiocorrect * D_val + (1.0-ratiocorrect)*D_def;
      //c_val = c_def;
      c_val = ratiocorrect * c_val + (1.0-ratiocorrect)*c_def;
      //rho_val = rho_def ;
      rho_val = ratiocorrect * rho_val + (1.0-ratiocorrect)*rho_def;


      int elematid = scatra_discret_->lRowElement(lid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::OptMatParManagerPerElement>(scatra_matman_)->WriteValuesToVector(elematid,i,D_val,reac_val,sca_tatparams);

      int agid = i - scatra_discret_->ElementRowMap()->MinAllGID() + acou_discret_->ElementRowMap()->MinAllGID();
      int alid = acou_discret_->ElementRowMap()->LID(agid);
      if(alid >=0)
      {
        int elematid = acou_discret_->lRowElement(alid)->Material()->Parameter()->Id();
        Teuchos::rcp_dynamic_cast<ACOU::AcouMatParManagerPerElement>(acou_matman_)->WriteValuesToVector(elematid,agid,rho_val,c_val,acou_tatparams);
      }
    }
    // bring the parameters to the materials
    acou_matman_->ReplaceParams(*acou_tatparams);
    scatra_matman_->ReplaceParams(*sca_tatparams);

    SolveStandardScatra();
    SolveStandardAcou();
    EvalulateObjectiveFunction();

  }
  dserror("that is it");
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiAcou::ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params)
{
  if(optoracou_==false)
    scatra_matman_->ReplaceParams(*params);
  else
    acou_matman_->ReplaceParams(*params);

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstructionOptiAcou::EvalulateObjectiveFunction()
{
  // evaluate error contribution
  EvaluateError();
  J_ = error_;

  // evaluate optical and acoustical regularization
  if(scatra_regman_ != Teuchos::null)
    scatra_regman_->Evaluate(*(scatra_matman_->GetParams()),&J_);
  if(acou_regman_ != Teuchos::null)
    acou_regman_->Evaluate(*(acou_matman_->GetParams()),&J_);

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<J_-error_<<std::endl;

  return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiAcou::EvaluateGradient()
{
  // evaluate scatra gradient
  ACOU::PatImageReconstructionOpti::EvaluateGradient();

  // evaluate acoustic gradient
  acou_objgrad_->Scale(0.0);

  // contribution from regularization
  if(acou_regman_ != Teuchos::null)
    acou_regman_->EvaluateGradient(*(acou_matman_->GetParams()),acou_objgrad_);

  // standard contribution
  acou_matman_->AddEvaluate(0.0,acou_objgrad_);
  acou_matman_->Finalize(acou_objgrad_);

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionOptiAcou::PerformIteration()
{
  std::cout<<std::endl;
  std::cout<<"OPTICAL LINE SEARCH"<<std::endl;
  std::cout<<std::endl;

  optoracou_ = false;
  linesearch_->Init(J_,scatra_objgrad_,searchdirection_->ComputeDirection(scatra_objgrad_,scatra_matman_->GetParams(),iter_),scatra_matman_->GetParams(),scatra_matman_->ParamLayoutMapUnique());
  bool optisucc = linesearch_->Run();

  std::cout<<std::endl;
  std::cout<<"ACOUSTICAL LINE SEARCH"<<std::endl;
  std::cout<<std::endl;

  optoracou_ = true;
  linesearch_->Init(J_,acou_objgrad_,acou_searchdirection_->ComputeDirection(acou_objgrad_,acou_matman_->GetParams(),iter_),acou_matman_->GetParams(),acou_matman_->ParamLayoutMapUnique());
  bool acousucc = linesearch_->Run();

  return (optisucc || acousucc);
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiAcou::CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result)
{
  if(optoracou_==false)
  {
    INVANA::MVDotProduct(*scatra_objgrad_,bvector,uniquemap,result);
  }
  else
  {
    INVANA::MVDotProduct(*acou_objgrad_,bvector,uniquemap,result);
  }

  return;
}

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionOptiAcouIdent::PatImageReconstructionOptiAcouIdent(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstructionOpti(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout),
  sequenzeiter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("SEQUENZE"))
{
  // create material manager
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  Teuchos::ParameterList list(invp);
  list.set<std::string>("PARAMLIST",acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("ACOUPARAMLIST"));
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(list,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
      acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerPerElement(acou_discret_));
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
      acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerUniform(acou_discret_));
    break;
    default:
      dserror("choose a valid method of parametrization for the material parameter field");
    break;
  }
  acou_matman_->Init(list);
  acou_matman_->Setup();

  // read materials
  ReadMaterials(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("SEGMENTATIONMATS"));

}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionOptiAcouIdent::ReadMaterials(std::string materialfilename)
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

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionOptiAcouIdent::PerformIteration()
{
  bool succ = false;

  for(int i=0; i<sequenzeiter_; ++i)
  {
    std::cout<<"ITERATION "<<i<<std::endl;
    linesearch_->Init(J_,scatra_objgrad_,searchdirection_->ComputeDirection(scatra_objgrad_,scatra_matman_->GetParams(),iter_),scatra_matman_->GetParams(),scatra_matman_->ParamLayoutMapUnique());
    succ = linesearch_->Run();

    std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
    std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;

    if(succ==false) return succ;
  }
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
void ACOU::PatImageReconstructionOptiAcouIdent::UpdateAcousticalParameters()
{
  Teuchos::RCP<Epetra_MultiVector> acou_p = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(), true));

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
      loc_reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
      loc_D = actele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(actele->Id()));
    }
    double D=0.0, reac = 0.0;
    scatra_discret_->Comm().SumAll(&loc_D,&D,1);
    scatra_discret_->Comm().SumAll(&loc_reac,&reac,1);

    // calculate the acoustical values which are required
    // first possibility: closest
    double c = 0.0, rho = 0.0;
    if(1)
    {
      double abst = 1.0e6;
      int mat = -1;
      for(unsigned m=0; m<nummats_; ++m)
      {
        double abstm = sqrt((reac-materialtable_[m][0])*(reac-materialtable_[m][0])+0.01*(D-materialtable_[m][1])*(D-materialtable_[m][1]));
        if(abstm<abst)
        {
          abst = abstm;
          mat = m;
        }
      }
      c = materialtable_[mat][2];
      rho = materialtable_[mat][3];
    }
    else // second possibility: average from all
    {

    }

    // write values to acoustical vector
    int agid = i - scatra_discret_->ElementRowMap()->MinAllGID() + acou_discret_->ElementRowMap()->MinAllGID();
    int alid = acou_discret_->ElementRowMap()->LID(agid);
    if(alid >=0)
    {
      int elematid = acou_discret_->lRowElement(alid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::AcouMatParManagerPerElement>(acou_matman_)->WriteValuesToVector(elematid,agid,rho,c,acou_p);
    }
  }
  acou_matman_->ReplaceParams(*acou_p);

  return;
}

/*----------------------------------------------------------------------*/
ACOU::PatImageReconstructionSegmentation::PatImageReconstructionSegmentation(
  Teuchos::RCP<DRT::Discretization>      scatradis,
  Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
  Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
  Teuchos::RCP<Teuchos::ParameterList>   acoupara,
  Teuchos::RCP<LINALG::Solver>           scatrasolv,
  Teuchos::RCP<LINALG::Solver>           acousolv,
  Teuchos::RCP<IO::DiscretizationWriter> scatraout,
  Teuchos::RCP<IO::DiscretizationWriter> acouout)
: PatImageReconstruction(scatradis,acoudis,scatrapara,acoupara,scatrasolv,acousolv,scatraout,acouout),
  penaltyparam_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("EQUALITYPENALTY")),
  penalty_(0.0),
  regularization_(0.0),
  gradweight_(DRT::INPUT::IntegralValue<bool>(acouparams_->sublist("PA IMAGE RECONSTRUCTION"),"GRADWEIGHTING")),
  sequenzeiter_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("SEQUENZE")),
  sequenze_(-1)
{
  // create material manager
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  Teuchos::ParameterList list(invp);
  list.set<std::string>("PARAMLIST",acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("ACOUPARAMLIST"));
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(list,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
      acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerPerElement(acou_discret_));
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
      acou_matman_ = Teuchos::rcp(new ACOU::AcouMatParManagerUniform(acou_discret_));
    break;
    default:
      dserror("choose a valid method of parametrization for the material parameter field");
    break;
  }
  acou_matman_->Init(list);
  acou_matman_->Setup();

  list.set<std::string>("PARAMLIST","3 REAC");
  k_dummy_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerPerElement(scatra_discret_));
  k_dummy_matman_->Init(list);
  k_dummy_matman_->Setup();

  // create regularization manager
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvRegularization>(invp,"REGULARIZATION"))
  {
  case INPAR::INVANA::stat_inv_reg_none:
    break;
  case INPAR::INVANA::stat_inv_reg_tikhonov:
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
    acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTikhonov());
  break;
  case INPAR::INVANA::stat_inv_reg_totalvariation:
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
    acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
  break;
  case INPAR::INVANA::stat_inv_reg_tvdtikh:
    scatra_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
    acou_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariationTikhonov());
  break;
  default:
    dserror("no valid regularization type provided");
  break;
  }
  if (scatra_regman_!=Teuchos::null)
  {
    scatra_regman_->Init(scatra_discret_,scatra_matman_->GetConnectivityData());
    scatra_regman_->Setup(invp);
    acou_regman_->Init(acou_discret_,acou_matman_->GetConnectivityData());
    acou_regman_->Setup(invp);
    k_regman_ = Teuchos::rcp(new INVANA::RegularizationTotalVariation());
    k_regman_->Init(scatra_discret_,k_dummy_matman_->GetConnectivityData());
    k_regman_->Setup(invp);
  }

  // read materials
  ReadMaterials(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("SEGMENTATIONMATS"));

  // create concentration vector and gradient
  //k_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),nummats_,true));
  //k_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),nummats_,true));
  k_ = Teuchos::rcp(new Epetra_MultiVector(*(k_dummy_matman_->ParamLayoutMap()),nummats_,true));
  k_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(k_dummy_matman_->ParamLayoutMap()),nummats_,true));


  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",true);

  // setup direction
  {
    Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(*scatra_discret_->ElementRowMap()));
    searchdirection_->Setup(map,map,nummats_);
  }
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::ReadMaterials(std::string materialfilename)
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

  return;
}

/*----------------------------------------------------------------------*/
double ACOU::PatImageReconstructionSegmentation::EvalulateObjectiveFunction()
{
  // evaluate error contribution
  EvaluateError();
  J_ = error_;

  // evaluate optical and acoustical regularization

  regularization_ = 0.0;
  if(scatra_regman_ != Teuchos::null)
    scatra_regman_->Evaluate(*(scatra_matman_->GetParams()),&regularization_);
  if(acou_regman_ != Teuchos::null)
    acou_regman_->Evaluate(*(acou_matman_->GetParams()),&regularization_);
 /*
  regularization_ = 0.0;
  if(scatra_regman_ != Teuchos::null)
    for(unsigned int i=0; i<nummats_; ++i)
      k_regman_->Evaluate(*k_->operator ()(i),&regularization_);
 */

  J_ += regularization_;

  // evaluate penalty term
  Teuchos::RCP<Epetra_Vector> kvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  kvec->PutScalar(-1.0);
  for(unsigned int m=0; m<nummats_; ++m)
    for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
      kvec->operator [](i) += 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;
  double val = 0.0;
  kvec->Norm2(&val);
  penalty_ = 0.5 * penaltyparam_ * val * val;
  J_ += penalty_;

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<regularization_<<" penalty "<<penalty_<<std::endl;

  if(sequenze_==0)
    return error_+regularization_;
  else if(sequenze_==1)
    return penalty_;
  else
    return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::EvaluateGradient()
{
  Teuchos::RCP<Epetra_MultiVector> scatra_objgrad = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
  Teuchos::RCP<Epetra_MultiVector> acou_objgrad = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));

  if(sequenze_==0 || sequenze_==-1)
  {
    // calculate scatra gradient
    {
      // set quantities needed by the elements
      scatra_discret_->SetState("adjoint phi",adjoint_phi_);
      scatra_discret_->SetState("phi",phi_);

      // fill and set psi vector
      Teuchos::RCP<Epetra_Vector> psi = LINALG::CreateVector(*(scatra_discret_->DofRowMap()),true);
      psi = CalculateAdjointOptiRhsvec(adjoint_psi_);
      scatra_discret_->SetState("psi",psi);

      // do the actual evaluation
      scatra_matman_->AddEvaluate(0.0,scatra_objgrad);
      scatra_matman_->Finalize(scatra_objgrad);
    }

    // calculate acoustic gradient
    {
      acou_matman_->AddEvaluate(0.0,acou_objgrad);
      acou_matman_->Finalize(acou_objgrad);
    }

    // contribution form regularization

    if(scatra_regman_ != Teuchos::null)
      scatra_regman_->EvaluateGradient(*(scatra_matman_->GetParams()),scatra_objgrad);

    // contribution from regularization
    if(acou_regman_ != Teuchos::null)
      acou_regman_->EvaluateGradient(*(acou_matman_->GetParams()),acou_objgrad);

/*
    if(scatra_regman_ != Teuchos::null)
      for(unsigned int i=0; i<nummats_; ++i)
      {
        Teuchos::RCP<Epetra_Vector> tempreggrad = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),true));
        k_regman_->EvaluateGradient(*k_->operator ()(i),tempreggrad);
        k_objgrad_->operator ()(i)->Update(1.0,*tempreggrad,1.0);
      }
*/
  }

  // build the concentration gradient
  int numeleentries = k_->MyLength();
  if(sequenze_ == 0)
  {
    double relerror = error_/error_start_;
    double weightreac = 20.0*relerror;
    double weightdiff = 1.0*relerror;
    double weightc    = 1.0*relerror;
    double weightrho  = 1.0*relerror;
    double sumweight = (weightreac + weightdiff + weightc + weightrho)/4.0;

    for(int i=0; i<numeleentries; ++i)
    {
      for(unsigned int m=0; m<nummats_; ++m)
      {
        double k_grad_val;
        if(gradweight_==false)
          k_grad_val = scatra_objgrad->operator ()(0)->operator [](i) * materialtable_[m][0]
                     + scatra_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][1]
                     + acou_objgrad->operator ()(0)->operator [](i) * materialtable_[m][2]
                     + acou_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][3];
        else
          k_grad_val = (weightreac * scatra_objgrad->operator ()(0)->operator [](i) * materialtable_[m][0]
                     + weightdiff * scatra_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][1]
                     + weightc * acou_objgrad->operator ()(0)->operator [](i) * materialtable_[m][2]
                     + weightrho * acou_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][3])/sumweight;

        // chain rule for meta parametrization
        double k = k_->operator ()(m)->operator [](i);
        k_grad_val *= 1./PI/(k*k+1.);

        // write it in the k gradient vector
        k_objgrad_->operator ()(m)->operator [](i) = k_grad_val;
      }
    }
  }
  else if(sequenze_ == 1)
  {
    for(int i=0; i<numeleentries; ++i)
    {
      // evaluate constraint violation
      double equ = -1.0;
      for(unsigned int m=0; m<nummats_; ++m)
        equ += 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;

      for(unsigned int m=0; m<nummats_; ++m)
      {
        // contribution from penalty term
        double k_grad_val = penaltyparam_ * equ;

        // chain rule for meta parametrization
        double k = k_->operator ()(m)->operator [](i);
        k_grad_val *= 1./PI/(k*k+1.);

        // write it in the k gradient vector
        k_objgrad_->operator ()(m)->operator [](i) = k_grad_val;
      }
    }
  }
  else if(sequenze_ == -1)
  {
    double relerror = error_/error_start_;
    double weightreac = 100.0*relerror;
    double weightdiff = 1.0*relerror;
    double weightc    = 20.0*relerror;
    double weightrho  = 10.0*relerror;
    double sumweight = (weightreac + weightdiff + weightc + weightrho)/4.0;

    for(int i=0; i<numeleentries; ++i)
    {
      // evaluate constraint violation
      double equ = -1.0;
      for(unsigned int m=0; m<nummats_; ++m)
        equ += 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;

      for(unsigned int m=0; m<nummats_; ++m)
      {
        double k_grad_val;
        if(gradweight_==false)
          k_grad_val = scatra_objgrad->operator ()(0)->operator [](i) * materialtable_[m][0]
                     + scatra_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][1]
                     + acou_objgrad->operator ()(0)->operator [](i) * materialtable_[m][2]
                     + acou_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][3];
        else
          k_grad_val = (weightreac * scatra_objgrad->operator ()(0)->operator [](i) * materialtable_[m][0]
                     + weightdiff * scatra_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][1]
                     + weightc * acou_objgrad->operator ()(0)->operator [](i) * materialtable_[m][2]
                     + weightrho * acou_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][3])/sumweight;

        // contribution from penalty term
        k_grad_val += penaltyparam_ * equ;

        // chain rule for meta parametrization
        double k = k_->operator ()(m)->operator [](i);
        k_grad_val *= 1./PI/(k*k+1.);

        // write it in the k gradient vector
        k_objgrad_->operator ()(m)->operator [](i) = k_grad_val;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::InitConcentrations()
{
  if(!timereversal_)
  {
    // the first material is assumed to be the default material and thus has highest concentration (default concentration defcon 0.9)
    double defcon = 0.9;
    double c = tan(PI*(defcon-0.5));
    k_->operator ()(0)->PutScalar(c);

    // the others make up the rest with equal parts
    c = tan(PI*((1.0-defcon)/(double(nummats_)-1.)-0.5));
    for(unsigned int m=1; m<nummats_; ++m)
      k_->operator ()(m)->PutScalar(c);
  }
  else
  {
    // in case of time reversal, diffusivity, speed of sound and density are on default, the reaction coefficient is set
    Teuchos::RCP<Epetra_MultiVector> sca_p = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(), true));
    Teuchos::RCP<Epetra_MultiVector> acou_p = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(), true));

    // loop the global scatra elements
    for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
    {
      // lid of the element
      int lid = scatra_discret_->ElementRowMap()->LID(i);

      // get reaction coefficient from this element
      double reac;
      if(lid>=0)
      {
        reac = scatra_discret_->gElement(i)->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(i));

        // find the two closest materials from the material table
        unsigned int closest = 0, secondclosest = 0;
        double closestdist = 1.0e6, secondclosestdist = 1.0e6;
        for(unsigned int m=0; m<nummats_; ++m)
          if(std::abs(reac-materialtable_[m][0])<closestdist)
          {
            closestdist = std::abs(reac-materialtable_[m][0]);;
            closest = m;
          }
          else if(std::abs(reac-materialtable_[m][0])<secondclosestdist)
          {
            secondclosestdist = std::abs(reac-materialtable_[m][0]);;
            secondclosest = m;
          }

        for(unsigned int m=0; m<nummats_; ++m)
          if(m==closest)
          {
            double k_closest = secondclosestdist/(secondclosestdist+closestdist)-0.005;
            if(k_closest <0.0) k_closest = (0.01)/(double(nummats_)-2.);
            k_->operator ()(closest)->operator [](lid) = tan(PI*(k_closest-0.5));
          }
          else if(m==secondclosest)
          {
            double k_secondclosest = closestdist/(secondclosestdist+closestdist)-0.005;
            if(k_secondclosest <0.0) k_secondclosest = (0.01)/(double(nummats_)-2.);
            k_->operator ()(secondclosest)->operator [](lid) = tan(PI*(k_secondclosest-0.5));
          }
          else
          {
            k_->operator ()(m)->operator [](lid) = tan(PI*((0.01)/(double(nummats_)-2.)-0.5));
          }

      }
    }
  }

  ReplaceParams(k_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params)
{
  k_->Update(1.0,*params,0.0);

  Teuchos::RCP<Epetra_MultiVector> sca_p = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(), true));
  Teuchos::RCP<Epetra_MultiVector> acou_p = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(), true));

  // loop the global scatra elements
  for(int i=scatra_discret_->ElementRowMap()->MinAllGID(); i<=scatra_discret_->ElementRowMap()->MaxAllGID(); ++i)
  {
    // lid of the element
    int lid = scatra_discret_->ElementRowMap()->LID(i);

    // determine local material values
    double loc_mat[4] = {0.0};
    if(lid>=0)
    {
      // determine material value as sum of concentration and correspondent value
      for(unsigned int m=0; m<nummats_; ++m)
      {
        double meta_k = params->operator ()(m)->operator [](lid);
        // calculate from meta to normal
        double k = 1./PI*atan(meta_k)+0.5;
        loc_mat[0] += k * materialtable_[m][0];
        loc_mat[1] += k * materialtable_[m][1];
        loc_mat[2] += k * materialtable_[m][2];
        loc_mat[3] += k * materialtable_[m][3];
      }
    }

    // communicate material values to all processors
    double glo_mat[4] = {0.0};
    scatra_discret_->Comm().SumAll(&loc_mat[0],&glo_mat[0],4);

    // write scatra values to scatra vector
    if(lid>=0)
    {
      int elematid = scatra_discret_->lRowElement(lid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::OptMatParManagerPerElement>(scatra_matman_)->WriteValuesToVector(elematid,i,glo_mat[1],glo_mat[0],sca_p);
    }

    // write acou values to acou vector
    int agid = i - scatra_discret_->ElementRowMap()->MinAllGID() + acou_discret_->ElementRowMap()->MinAllGID();
    int alid = acou_discret_->ElementRowMap()->LID(agid);
    if(alid >=0)
    {
      int elematid = acou_discret_->lRowElement(alid)->Material()->Parameter()->Id();
      Teuchos::rcp_dynamic_cast<ACOU::AcouMatParManagerPerElement>(acou_matman_)->WriteValuesToVector(elematid,agid,glo_mat[3],glo_mat[2],acou_p);
    }
  }

  // bring the parameters to the materials
  scatra_matman_->ReplaceParams(*sca_p);
  acou_matman_->ReplaceParams(*acou_p);

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PatImageReconstructionSegmentation::PerformIteration()
{
  bool succ = false;
  if(0)
  {
    if(sequenzeiter_>0) // sequential optimization
    {
      // optimize the parameters without penalty
      sequenze_ = 0;
      EvalulateObjectiveFunction();
      EvaluateGradient();
      bool parasucc;
      for(int i=0; i<sequenzeiter_; ++i)
      {
        std::cout<<"SEQUENZE "<<sequenze_<<" iteration "<<i<<std::endl;
        linesearch_->Init(error_+regularization_,k_objgrad_,searchdirection_->ComputeDirection(k_objgrad_,k_,i),k_,k_dummy_matman_->ParamLayoutMapUnique());
        parasucc = linesearch_->Run();
        std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
        std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
        std::cout<<"*** penalty value                     "<<penalty_<<std::endl;
        if(parasucc==false)
          dserror("sequenze 0 failed line search");
      }

      // optimize the penalty only
      sequenze_ = 1;
      EvalulateObjectiveFunction();
      EvaluateGradient();
      bool penalsucc;
      for(int i=0; i<sequenzeiter_; ++i)
      {
        std::cout<<"SEQUENZE "<<sequenze_<<" iteration "<<i<<std::endl;
        linesearch_->Init(penalty_,k_objgrad_,searchdirection_->ComputeDirection(k_objgrad_,k_,i),k_,k_dummy_matman_->ParamLayoutMapUnique());
        penalsucc = linesearch_->Run();
        std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
        std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
        std::cout<<"*** penalty value                     "<<penalty_<<std::endl;
        if(penalty_<1e-3)
        {
          penalsucc = true;
          break;
        }
        if(penalsucc==false || penalty_<1e-3)
          break;
      }

      SolveStandardScatra();
      SolveStandardAcou();

      succ = (parasucc&&penalsucc);
    }
    else
    {
      sequenze_ = -1;
      EvaluateGradient();
      linesearch_->Init(J_,k_objgrad_,searchdirection_->ComputeDirection(k_objgrad_,k_,iter_),k_,k_dummy_matman_->ParamLayoutMapUnique());
      succ = linesearch_->Run();
    }
  }
  else
  {
    for(int i=0; i<5; ++i)
    {
      sequenze_ = 0; // error and regularization
      std::cout<<"SEQUENZE "<<sequenze_<<" iteration "<<i<<std::endl;
      EvaluateGradient();
      EvalulateObjectiveFunction();
      linesearch_->Init(J_,k_objgrad_,searchdirection_->ComputeDirection(k_objgrad_,k_,iter_),k_,k_dummy_matman_->ParamLayoutMapUnique());
      succ = linesearch_->Run();
      std::cout<<"*** relative objective function value "<<J_/J_start_<<std::endl;
      std::cout<<"*** relative error value              "<<error_/error_start_<<std::endl;
      std::cout<<"*** penalty value                     "<<penalty_<<std::endl;
      if(succ==false) break;
    }
    if(succ==true)
    {
      EqualityCompensation();

      SolveStandardScatra();
      SolveStandardAcou();
      EvalulateObjectiveFunction();
      SolveAdjointAcou();
      SolveAdjointScatra();
      EvaluateGradient();
    }
  }

  return succ;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::EqualityCompensation()
{
  for(int i=0; i<k_->MyLength(); ++i)
  {
    double equ = 0.0;
    for(int m=0; m<k_->NumVectors(); ++m)
      equ += 1./PI *atan(k_->operator ()(m)->operator [](i))+0.5;
    for(int m=0; m<k_->NumVectors(); ++m)
    {
      double val = (1./PI *atan(k_->operator ()(m)->operator [](i)) +0.5) / equ;
      k_->operator ()(m)->operator [](i) = tan((val-0.5)*PI);
    }
  }
  ReplaceParams(k_);
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result)
{
  INVANA::MVDotProduct(*k_objgrad_,bvector,uniquemap,result);
  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::OutputReactionAndDiffusion()
{
  /*
  // build the two vectors
  Teuchos::RCP<Epetra_Vector> reacvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  Teuchos::RCP<Epetra_Vector> diffvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);
    double reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
    double diff = actele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(actele->Id()));
    reacvec->operator [](i) = reac;
    diffvec->operator [](i) = diff;
  }
  //scatraoutput_->WriteVector("rea_coeff",reacvec);
  //scatraoutput_->WriteVector("diff_coeff",diffvec);
  */

  for(unsigned int m=0; m<nummats_; ++m)
  {
    std::ostringstream nameos;
    nameos<< "k_" << m+1;
    Teuchos::RCP<Epetra_Vector> kvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
    for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
      kvec->operator [](i) = 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;
    scatraoutput_->WriteVector(nameos.str(),kvec);
  }

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
  Teuchos::RCP<Epetra_Vector> reacvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);
    double reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
    reacvec->operator [](i) = reac;
  }
  return reacvec;
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

  // determine the number of vectors for monitoring
  // this is naive: later on, we won't be able to store everything at once and we have to implement
  // a smarter approach to reduce storage requirements
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
  if (monitorfilename=="none.monitor") dserror("No monitor file provided");
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
        if(distance[j]<=epsilon)
        {
          for(l=0;l<nsteps;l++)
            nodcurvinterpol[i*nsteps+l]=mcurve[j+l*nmics];
          must_set_curve=0;
        }
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
      if(dtacou/(timesteps[1]-timesteps[0]) - std::ceil(dtacou/(timesteps[1]-timesteps[0]))< 1e-16 ) // dtacou is a multiple of the monitor time step
      {
        int mult = std::ceil(dtacou/(timesteps[1]-timesteps[0]));
        for(unsigned int i=0; i<pressuremonmicsunique.size(); ++i)
          if( acou_discret_->HaveGlobalNode(pressuremonmicsunique[i]) )
            for(int j=0; j<numvec; j++)
              acou_rhsm_->ReplaceGlobalValue(pressuremonmicsunique[i],j,nodcurvinterpol(i*nsteps+j*mult)); // the proc who has this row, writes the value
      }
      else
        dserror("time step bigger than monitor time step but no multiple -> implement here!");
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
              acou_rhsm_->ReplaceGlobalValue(nodes_[i],j,value);
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

  eps = min_abs*10.0;

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
      scatraoutput_->NewResultFile(outname,output_count_);
      output_count_++;
      scatraoutput_->WriteMesh(0,0.0);

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      scatraalgo_ = Teuchos::rcp(new SCATRA::TimIntStationary(scatra_discret_, scatrasolver_, scatraparams_, scatraextraparams_, scatraoutput_));

      scatraalgo_->Init();
      scatraalgo_->SetVelocityField(1);

      scatraalgo_->TimeLoop();

      // output of elemental reaction coefficient
      OutputReactionAndDiffusion();

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
  // set parameter indicating that the forward problem is solved
  acouparams_->set<bool>("adjoint",false);

  std::string outname = name_;
  outname.append("_invforward_acou");
  acououtput_->NewResultFile(outname,output_count_);
  last_acou_fw_output_count_ = output_count_;
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
  acoualgo_->Integrate(acou_rhs_,abcnodes_map_);

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstruction::SolveAdjointAcou()
{
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

  // acou_rhs_ has to be scaled with weighting (adjoint of the mapping)
  Teuchos::RCP<Epetra_Vector> touchcountvec = LINALG::CreateVector(*abcnodes_map_);
  acoualgo_->FillTouchCountVec(touchcountvec);
  tempvec->Multiply(1.0,*touchcountvec,*tempvec,0.0);

  // set the difference between measured and simulated values
  acouparams_->set<Teuchos::RCP<Epetra_MultiVector> >("rhsvec",tempvec);

  // prepare the output
  std::string outname = name_;
  outname.append("_invadjoint_acou");
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
  acououtput_->NewResultFile(outname,output_count_);
  last_acou_fw_output_count_ = output_count_;
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
  acoualgo_->Integrate(acou_rhs_,abcnodes_map_);

  // reset parameter
  acouparams_->set<bool>("timereversal",false);

  //} time reversal run finished

  // now update the optical parameters
  // 1.) solve optical problem with initial guess for absorption coefficient
  // 2.) calculate mu_a as -p_0/Gamma/phi
  // 3.) bring these values to the parameter vector

  // do step 1.
  SolveStandardScatra(); // phi now holds the optical solution values

  // do step 2. and 3.
  UpdateAbsorptionCoefficientFromTimeReversal();

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
  Teuchos::RCP<Epetra_MultiVector> trparams = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));

  // do the business
  if(meshconform_)
  {
    for(int e=scatra_discret_->ElementRowMap()->MinAllGID(); e<scatra_discret_->ElementRowMap()->MaxAllGID(); ++e)
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
          int elematid = opti_ele->Material()->Parameter()->Id();
          double diff = opti_ele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(e));
          Teuchos::rcp_dynamic_cast<ACOU::OptMatParManagerPerElement>(scatra_matman_)->WriteValuesToVector(elematid,e,diff,reac,trparams);
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
          int elematid = opti_ele->Material()->Parameter()->Id();
          double diff = opti_ele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(e));
          Teuchos::rcp_dynamic_cast<ACOU::OptMatParManagerPerElement>(scatra_matman_)->WriteValuesToVector(elematid,e,diff,reac,trparams);
        }
      }
    }
  }
  else
    dserror("update of absorption coefficient not yet implemented for nonconforming mesh");

  // bring values to the elements
  scatra_matman_->ReplaceParams(*trparams);

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
  /*
  // build the two vectors
  Teuchos::RCP<Epetra_Vector> reacvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  Teuchos::RCP<Epetra_Vector> diffvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));

  for (int i=0; i<scatra_discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = scatra_discret_->lRowElement(i);
    double reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
    double diff = actele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(actele->Id()));
    reacvec->operator [](i) = reac;
    diffvec->operator [](i) = diff;
  }
  scatraoutput_->WriteVector("rea_coeff",reacvec);
  scatraoutput_->WriteVector("diff_coeff",diffvec);
   */
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
          loc_mu_a += actmat->ReaCoeff(scatra_discret_->ElementColMap()->LID(roptele->Id()));
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
      if(acou_discret_->NodeRowMap()->LID(nd)>-1)
      {
        loc_value = adjoint_psi_->operator [](acou_discret_->NodeRowMap()->LID(nd));
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
