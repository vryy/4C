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
J_(0.0),
J_start_(0.0),
error_(0.0),
error_start_(0.0)
{
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
  switch(DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvMatParametrization>(invp,"PARAMETRIZATION"))
  {
    case INPAR::INVANA::stat_inv_mp_elementwise:
      acou_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerPerElement(acou_discret_));
    break;
    case INPAR::INVANA::stat_inv_mp_uniform:
      acou_matman_ = Teuchos::rcp(new ACOU::OptMatParManagerUniform(acou_discret_));
    break;
    default:
      dserror("choose a valid method of parametrization for the material parameter field");
    break;
  }
  acou_matman_->Init(invp);
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
  acou_searchdirection_->Setup(acou_matman_->ParamLayoutMap(),acou_matman_->ParamLayoutMapUnique(),1);
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
  optoracou_ = false;
  linesearch_->Init(J_,scatra_objgrad_,searchdirection_->ComputeDirection(scatra_objgrad_,scatra_matman_->GetParams(),iter_),scatra_matman_->GetParams(),scatra_matman_->ParamLayoutMapUnique());
  bool optisucc = linesearch_->Run();

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
  penaltyparam_(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<double>("EQUALITYPENALTY"))
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
  }

  // read materials
  ReadMaterials(acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<std::string>("SEGMENTATIONMATS"));

  // create concentration vector and gradient
  k_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),nummats_,true));
  k_objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*(scatra_discret_->ElementRowMap()),nummats_,true));

  // set parameter for aocustic time integration
  acouparams_->set<bool>("acouopt",true);

  // initialize the material values
  InitConcentrations();

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
  if(scatra_regman_ != Teuchos::null)
    scatra_regman_->Evaluate(*(scatra_matman_->GetParams()),&J_);
  if(acou_regman_ != Teuchos::null)
    acou_regman_->Evaluate(*(acou_matman_->GetParams()),&J_);

  // evaluate penalty term
  Teuchos::RCP<Epetra_Vector> kvec = Teuchos::rcp(new Epetra_Vector(*(scatra_discret_->ElementRowMap()),false));
  kvec->PutScalar(-1.0);
  for(unsigned int m=0; m<nummats_; ++m)
    for(int i=0; i<scatra_discret_->NumMyRowElements(); ++i)
      kvec->operator [](i) += 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;
  double val = 0.0;
  kvec->Norm2(&val);
  double penaltycontrib = 0.5 * penaltyparam_ * val * val;
  J_ += penaltycontrib;

  // output
  if(!myrank_)
    std::cout<<"objective function value "<<J_<<" error value "<<error_<<" regularization "<<J_-error_-penaltycontrib<<" penalty "<<penaltycontrib<<std::endl;

  return J_;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::EvaluateGradient()
{
  // calculate scatra gradient
  Teuchos::RCP<Epetra_MultiVector> scatra_objgrad = Teuchos::rcp(new Epetra_MultiVector(*(scatra_matman_->ParamLayoutMap()), scatra_matman_->NumVectors(),true));
  {
    // contribution form regularization
    if(scatra_regman_ != Teuchos::null)
      scatra_regman_->EvaluateGradient(*(scatra_matman_->GetParams()),scatra_objgrad);

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
  Teuchos::RCP<Epetra_MultiVector> acou_objgrad = Teuchos::rcp(new Epetra_MultiVector(*(acou_matman_->ParamLayoutMap()), acou_matman_->NumVectors(),true));
  {
    // contribution from regularization
    if(acou_regman_ != Teuchos::null)
      acou_regman_->EvaluateGradient(*(acou_matman_->GetParams()),acou_objgrad);

    // standard contribution
    acou_matman_->AddEvaluate(0.0,acou_objgrad);
    acou_matman_->Finalize(acou_objgrad);
  }

  // build the concentration gradient
  int numeleentries = k_->MyLength();
  for(int i=0; i<numeleentries; ++i)
  {
    // evaluate constraint violation
    double equ = -1.0;
    for(unsigned int m=0; m<nummats_; ++m)
      equ += 1./PI*atan(k_->operator ()(m)->operator [](i))+0.5;

    for(unsigned int m=0; m<nummats_; ++m)
    {
      double k_grad_val = scatra_objgrad->operator ()(0)->operator [](i) * materialtable_[m][0]
                        + scatra_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][1]
                        + acou_objgrad->operator ()(0)->operator [](i) * materialtable_[m][2]
                        + acou_objgrad->operator ()(0)->operator [](i+numeleentries) * materialtable_[m][3];

      // contribution from penalty term
      k_grad_val += penaltyparam_ * equ;

      // chain rule for meta parametrization
      double k = k_->operator ()(m)->operator [](i);
      k_grad_val *= 1./PI/(k*k+1.);

      // write it in the k gradient vector
      k_objgrad_->operator ()(m)->operator [](i) = k_grad_val;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::InitConcentrations()
{
  // the first material is assumed to be the default material and thus has highest concentration (default concentration defcon 0.9)
  double defcon = 0.9;
  double c = tan(PI*(defcon-0.5));
  k_->operator ()(0)->PutScalar(c);

  // the others make up the rest with equal parts
  c = tan(PI*((1.0-defcon)/(double(nummats_)-1.)-0.5));
  for(unsigned int m=1; m<nummats_; ++m)
    k_->operator ()(m)->PutScalar(c);

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
  linesearch_->Init(J_,k_objgrad_,searchdirection_->ComputeDirection(k_objgrad_,k_,iter_),k_,Teuchos::rcp(new Epetra_Map(*scatra_discret_->ElementRowMap())));
  return linesearch_->Run();
}

/*----------------------------------------------------------------------*/
void ACOU::PatImageReconstructionSegmentation::CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result)
{
  INVANA::MVDotProduct(*k_objgrad_,bvector,uniquemap,result);
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
    int dofgid = scatra_discret_->Dof(scatra_discret_->lRowNode(i),0);
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
        int dofgid = scatra_discret_->Dof(opti_node,0);
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
    double reac = actele->Material()->Parameter()->GetParameter(1,scatra_discret_->ElementColMap()->LID(actele->Id()));
    double diff = actele->Material()->Parameter()->GetParameter(0,scatra_discret_->ElementColMap()->LID(actele->Id()));
    reacvec->operator [](i) = reac;
    diffvec->operator [](i) = diff;
  }
  scatraoutput_->WriteVector("rea_coeff",reacvec);
  scatraoutput_->WriteVector("diff_coeff",diffvec);

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
        int dofgid = scatra_discret_->Dof(opti_node,0);
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
        int dof = scatra_discret_->Dof(optnode,0);
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
}
