/*----------------------------------------------------------------------*/
/*!
 * \file gen_inv_analysis.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gee
            089 - 289-15239
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "gen_inv_analysis.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_io/io_hdf.H"
#include "../drt_mat/material.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_matelast/elast_coupanisoexpotwo.H"
#include "../drt_matelast/elast_coupanisoneohooketwo.H"
#include "../drt_matelast/elast_coupblatzko.H"
#include "../drt_matelast/elast_couplogneohooke.H"
#include "../drt_matelast/elast_isoexpo.H"
#include "../drt_matelast/elast_isomooneyrivlin.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_isoyeoh.H"
#include "../drt_matelast/elast_volpenalty.H"
#include "../drt_matelast/elast_vologden.H"
#include "../drt_matelast/elast_volsussmanbathe.H"

using namespace std;
using namespace DRT;
using namespace MAT;


#include "../drt_structure/stru_resulttest.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::GenInvAnalysis::GenInvAnalysis(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    solver_(solver),
    output_(output),
    sti_(Teuchos::null)
{
  int myrank = dis->Comm().MyPID();
  
  reset_out_count_=0;

  // input parameters structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // input parameters inverse analysis
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  //  tolerance for the curve fitting
  tol_ = iap.get<double>("INV_ANA_TOL");
  
  /* this is how a monitor file should look like:
  
steps 25 nnodes 5 
11 2 0 1
19 2 0 1
27 2 0 1 
35 2 0 1
42 2 0 1
# any number of comment lines, but only at this position
# x and y disps at 5 nodes
# time node 11 x y 19 x y 27 x y 35 x y 42 x y
#
2.000000e-02   -9.914936e-02    7.907108e-02   -2.421293e-01    3.475674e-01   -2.018457e-01    6.164171e-01   2.622612e-02    5.351817e-01  2.623198e-03    6.524309e-02
4.000000e-02   -2.060081e-01    1.561979e-01   -4.495438e-01    6.425663e-01   -3.561049e-01    1.099916e+00   7.625965e-02    9.675498e-01  4.770331e-02    1.287581e-01
6.000000e-02   -3.063373e-01    2.286715e-01   -6.147396e-01    8.904266e-01   -4.628068e-01    1.481550e+00   1.405767e-01    1.313955e+00  1.144058e-01    1.871826e-01
8.000000e-02   -3.995311e-01    2.965706e-01   -7.500371e-01    1.103765e+00   -5.399646e-01    1.795528e+00   2.080302e-01    1.600827e+00  1.902863e-01    2.404236e-01
1.000000e-01   -4.871174e-01    3.603335e-01   -8.656673e-01    1.291507e+00   -6.003143e-01    2.063556e+00   2.732486e-01    1.846735e+00  2.684362e-01    2.890503e-01
1.200000e-01   -5.705636e-01    4.204134e-01   -9.684391e-01    1.459720e+00   -6.514304e-01    2.298985e+00   3.338117e-01    2.063424e+00  3.449541e-01    3.336880e-01
1.400000e-01   -6.508979e-01    4.772118e-01   -1.062618e+00    1.612613e+00   -6.975891e-01    2.510286e+00   3.888783e-01    2.258437e+00  4.178564e-01    3.749320e-01
1.600000e-01   -7.287384e-01    5.310691e-01   -1.150807e+00    1.753191e+00   -7.411234e-01    2.703069e+00   4.384501e-01    2.436805e+00  4.863525e-01    4.133178e-01
1.800000e-01   -8.044079e-01    5.822715e-01   -1.234572e+00    1.883657e+00   -7.832539e-01    2.881204e+00   4.829452e-01    2.602019e+00  5.503314e-01    4.493058e-01
2.000000e-01   -8.780469e-01    6.310615e-01   -1.314853e+00    2.005671e+00   -8.245805e-01    3.047458e+00   5.229514e-01    2.756584e+00  6.100265e-01    4.832787e-01
2.200000e-01   -9.496984e-01    6.776474e-01   -1.392223e+00    2.120513e+00   -8.653673e-01    3.203878e+00   5.590896e-01    2.902354e+00  6.658185e-01    5.155483e-01
2.400000e-01   -1.019365e+00    7.222113e-01   -1.467040e+00    2.229192e+00   -9.057055e-01    3.352017e+00   5.919446e-01    3.040735e+00  7.181279e-01    5.463664e-01
2.600000e-01   -1.087042e+00    7.649147e-01   -1.539546e+00    2.332515e+00   -9.456044e-01    3.493080e+00   6.220354e-01    3.172817e+00  7.673616e-01    5.759360e-01
2.800000e-01   -1.152735e+00    8.059027e-01   -1.609916e+00    2.431140e+00   -9.850398e-01    3.628021e+00   6.498057e-01    3.299461e+00  8.138890e-01    6.044213e-01
3.000000e-01   -1.216466e+00    8.453065e-01   -1.678289e+00    2.525611e+00   -1.023980e+00    3.757606e+00   6.756268e-01    3.421355e+00  8.580353e-01    6.319557e-01
3.200000e-01   -1.278278e+00    8.832457e-01   -1.744784e+00    2.616378e+00   -1.062399e+00    3.882462e+00   6.998053e-01    3.539062e+00  9.000811e-01    6.586487e-01
3.400000e-01   -1.338225e+00    9.198289e-01   -1.809512e+00    2.703823e+00   -1.100276e+00    4.003104e+00   7.225922e-01    3.653044e+00  9.402670e-01    6.845912e-01
3.600000e-01   -1.396378e+00    9.551554e-01   -1.872574e+00    2.788270e+00   -1.137603e+00    4.119965e+00   7.441921e-01    3.763688e+00  9.787981e-01    7.098591e-01
3.800000e-01   -1.452811e+00    9.893160e-01   -1.934068e+00    2.869995e+00   -1.174379e+00    4.233408e+00   7.647719e-01    3.871317e+00  1.015850e+00    7.345163e-01
4.000000e-01   -1.507606e+00    1.022393e+00   -1.994085e+00    2.949237e+00   -1.210611e+00    4.343745e+00   7.844680e-01    3.976208e+00  1.051572e+00    7.586177e-01
4.200000e-01   -1.560847e+00    1.054463e+00   -2.052713e+00    3.026205e+00   -1.246311e+00    4.451242e+00   8.033921e-01    4.078601e+00  1.086093e+00    7.822100e-01
4.400000e-01   -1.612617e+00    1.085593e+00   -2.110034e+00    3.101081e+00   -1.281495e+00    4.556132e+00   8.216362e-01    4.178700e+00  1.119525e+00    8.053342e-01
4.600000e-01   -1.662996e+00    1.115847e+00   -2.166128e+00    3.174023e+00   -1.316181e+00    4.658615e+00   8.392767e-01    4.276685e+00  1.151964e+00    8.280258e-01
4.800000e-01   -1.712064e+00    1.145282e+00   -2.221065e+00    3.245174e+00   -1.350390e+00    4.758870e+00   8.563772e-01    4.372714e+00  1.183494e+00    8.503160e-01
5.000000e-01   -1.759895e+00    1.173949e+00   -2.274916e+00    3.314659e+00   -1.384141e+00    4.857055e+00   8.729911e-01    4.466927e+00  1.214188e+00    8.722325e-01


  */
  
  // open monitor file and read it
  vector<double> timesteps;
  ndofs_ = 0;
  {
    char* foundit = NULL;
    string filename = iap.get<string>("MONITORFILE");
    if (filename=="none.monitor") dserror("No monitor file provided");
    FILE* file = fopen(filename.c_str(),"rb");
    char buffer[150000];
    fgets(buffer,150000,file);
    // read steps
    foundit = strstr(buffer,"steps"); foundit += strlen("steps");
    nsteps_ = strtol(foundit,&foundit,10);
    timesteps.resize(nsteps_);
    // read nnodes
    foundit = strstr(buffer,"nnodes"); foundit += strlen("nnodes");
    nnodes_ = strtol(foundit,&foundit,10);
    // read nodes
    nodes_.resize(nnodes_);
    dofs_.resize(nnodes_);
    for (int i=0; i<nnodes_; ++i)
    {
      fgets(buffer,150000,file);
      foundit = buffer;
      nodes_[i] = strtol(foundit,&foundit,10);
      int ndofs = strtol(foundit,&foundit,10);
      ndofs_ += ndofs;
      dofs_[i].resize(ndofs,-1);
      if (!myrank) printf("Monitored node %d ndofs %d dofs ",nodes_[i],(int)dofs_[i].size());
      for (int j=0; j<ndofs; ++j)
      {
        dofs_[i][j] = strtol(foundit,&foundit,10);
        if (!myrank) printf("%d ",dofs_[i][j]);
      }
      if (!myrank) printf("\n");
    }

    // total number of measured dofs
    nmp_    = ndofs_*nsteps_;
    mcurve_ = Epetra_SerialDenseVector(nmp_);

    if (!myrank) printf("nsteps %d ndofs %d\n",nsteps_,ndofs_);

    // time step from input file (do I need this?)
    tstep_  = sdyn.get<double>("TIMESTEP");
    
    // read comment lines
    foundit = buffer;
    fgets(buffer,150000,file);
    while(strstr(buffer,"#")) 
      fgets(buffer,150000,file);
    
    // read in the values for each node in dirs directions
    int count=0;
    for (int i=0; i<nsteps_; ++i)
    {
      // read the time step
      timesteps[i] = strtod(foundit,&foundit);
      for (int j=0; j<ndofs_; ++j)
        mcurve_[count++] = strtod(foundit,&foundit);
      fgets(buffer,150000,file);
      foundit = buffer;
    }
    if (count != nmp_) dserror("Number of measured disps wrong on input");
  }
  
  // error: diference of the measured to the calculated curve
  error_  = 1.0E6;
  error_o_= 1.0E6;

  // trainings parameter
  mu_ = iap.get<double>("INV_INITREG");

  kappa_multi_=1.0;

  // read material parameters from input file
  ReadInParameters();

  // Number of moterial parameters
  np_ = p_.Length();

  // controlling parameter
  numb_run_ =  0;     // counter of how many runs were made in the inverse analysis
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::GenInvAnalysis::Integrate()
{
  int myrank = discret_->Comm().MyPID();
  int max_itter = 100;
  output_->NewResultFile((numb_run_-1));
  // fitting loop
  do
  {
    if (!myrank)
    {
      cout << "#################################################################" << endl;
      cout << "########################### making Jacobian matrix ##############" <<endl;
      cout << "#################################################################" << endl;
      printf("Measured parameters nmp_ %d # parameters to fit np_ %d\n",nmp_,np_);
    }

    // pertubation of material parameter (should be relativ to the value that is perturbed)
    vector<double> perturb(np_,0.0);
    const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
    double alpha = iap.get<double>("INV_ALPHA");
    double beta  = iap.get<double>("INV_BETA");
    for (int i=0; i<np_; ++i) 
    {
      perturb[i] = alpha + beta * p_[i];
      if (!myrank) printf("perturbation[%d] %15.10e\n",i,perturb[i]);
    }
    
    // as the actual inv analysis is on proc 0 only, do only provide storage on proc 0
    Epetra_SerialDenseMatrix cmatrix;
    if (!myrank) cmatrix.Shape(nmp_,np_+1);
    
    // loop over parameters to fit and build cmatrix
    for (int i=0; i<np_+1;i++)
    {
      if (!myrank)
        cout << "--------------------------- run "<< i+1 << " of: " << np_+1 <<" -------------------------" <<endl;
      // make current set of material parameters
      Epetra_SerialDenseVector p_cur = p_;
      // perturb parameter i
      if (i!= np_) p_cur[i] = p_[i] + perturb[i];
      
      // put perturbed material parameters to material laws
      SetParameters(p_cur);
      
      // compute nonlinear problem and obtain computed displacements
      // output only for last run
      Epetra_SerialDenseVector cvector;
      if ( i != np_) cvector = CalcCvector(false);
      else           cvector = CalcCvector(true);
      
      // copy displacements to sensitivity matrix
      if (!myrank)
        for (int j=0; j<nmp_;j++)
          cmatrix(j,i) = cvector[j];
    }
    
    discret_->Comm().Barrier();
    if (!myrank)
      CalcNewParameters(cmatrix,perturb);
    discret_->Comm().Barrier();
      
    // set new material parameters
    SetParameters(p_);
    numb_run_++;
    discret_->Comm().Broadcast(&error_o_,1,0);
    discret_->Comm().Broadcast(&error_,1,0);
    discret_->Comm().Broadcast(&numb_run_,1,0);
    if (error_>tol_ && numb_run_<max_itter) 
      output_->NewResultFile(numb_run_);
  } while (error_>tol_ && numb_run_<max_itter);

  return;
}

//---------------------------------------------------------------------------------------------
void STR::GenInvAnalysis::CalcNewParameters(Epetra_SerialDenseMatrix& cmatrix, vector<double>& perturb)
{
  // initalization of the Jacobian and other storage
  Epetra_SerialDenseMatrix sto(np_,  np_);
  Epetra_SerialDenseVector delta_p(np_);
  Epetra_SerialDenseVector tmp(np_);
  Epetra_SerialDenseVector rcurve(nmp_);
  Epetra_SerialDenseVector ccurve(nmp_);

  // copy column with unperturbed values to extra vector
  for (int i=0; i<nmp_; i++) ccurve[i] = cmatrix(i,np_);

  // remove the extra column np_+1 with unperturbed values
  cmatrix.Reshape(nmp_,np_);
  
  // reuse the cmatrix array as J to save storage
  Epetra_SerialDenseMatrix& J = cmatrix;
  
  for (int i=0; i<nmp_; i++)
    for (int j=0; j<np_; j++)
    {
      J(i,j) -= ccurve[i]; 
      J(i,j) /= perturb[j];
    }

  sto.Multiply('T','N',1.0,J,J,0.0);

  // do regularization by adding artifical lumped mass
  for (int i=0; i<np_; i++)
    sto(i,i) += mu_*sto(i,i);

  // compute residual displacement (measured vs. computed)
  for (int i=0; i<nmp_; i++)
    rcurve[i] = mcurve_[i] - ccurve[i];

  // delta_p = (J.T*J+mu*diag(J^T * J))^{-1} * J^T * r
  tmp.Multiply('T','N',1.0,J,rcurve,0.0);
  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(sto);
  solver.SetVectors(delta_p,tmp);
  solver.SolveToRefinedSolution(true);
  solver.Solve();
  
  for (int i=0;i<np_;i++)
    p_[i] += delta_p[i];

  // dependent on the # of steps
  error_o_ = error_;
  error_   = rcurve.Norm2()/sqrt(nmp_);
  if (!numb_run_) error_o_ = error_;


  //Adjust training parameter
  mu_ *= (error_/error_o_);

  PrintStorage(delta_p);

  // return cmatrix to previous size and zero out
  cmatrix.Shape(nmp_,np_+1);

  return;
}


//--------------------------------------------------------------------------------------------
Epetra_SerialDenseVector STR::GenInvAnalysis::CalcCvector(bool outputtofile)
{
  int myrank = discret_->Comm().MyPID();
  
  // create a StruGenAlpha solver
  ParameterList genalphaparams;
  StruGenAlpha::SetDefaults(genalphaparams);
  const Teuchos::ParameterList& sdyn    = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& ioflags = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();

  // get input parameter lists
  genalphaparams.set<string>("DYNAMICTYP",sdyn.get<string>("DYNAMICTYP"));
  genalphaparams.set<bool>  ("damping",(! (sdyn.get<std::string>("DAMPING") == "no"
                                        || sdyn.get<std::string>("DAMPING") == "No"
                                        || sdyn.get<std::string>("DAMPING") == "NO")));
  genalphaparams.set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
  genalphaparams.set<double>("damping factor M",sdyn.get<double>("M_DAMP"));
#ifdef STRUGENALPHA_BE
  genalphaparams.set<double>("delta",sdyn.get<double>("DELTA"));
#endif
  genalphaparams.set<double>("gamma",sdyn.get<double>("GAMMA"));
  genalphaparams.set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
  genalphaparams.set<double>("alpha f",sdyn.get<double>("ALPHA_F"));
  genalphaparams.set<double>("total time",0.0);
  genalphaparams.set<double>("delta time",sdyn.get<double>("TIMESTEP"));
  genalphaparams.set<double>("max time",sdyn.get<double>("MAXTIME"));
  genalphaparams.set<int>   ("step",0);
  genalphaparams.set<int>   ("nstep",sdyn.get<int>("NUMSTEP"));
  genalphaparams.set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
  genalphaparams.set<int>   ("num iterations",-1);
  genalphaparams.set<string>("convcheck", sdyn.get<string>("CONV_CHECK"));
  genalphaparams.set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));
  genalphaparams.set<double>("tolerance residual",sdyn.get<double>("TOLRES"));
  genalphaparams.set<double>("tolerance constraint",sdyn.get<double>("TOLCONSTR"));
  genalphaparams.set<double>("UZAWAPARAM",sdyn.get<double>("UZAWAPARAM"));
  genalphaparams.set<double>("UZAWATOL",sdyn.get<double>("UZAWATOL"));
  genalphaparams.set<int>   ("UZAWAMAXITER",sdyn.get<int>("UZAWAMAXITER"));
  genalphaparams.set<INPAR::STR::ConSolveAlgo>("UZAWAALGO",getIntegralValue<INPAR::STR::ConSolveAlgo>(sdyn,"UZAWAALGO"));
  genalphaparams.set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
  genalphaparams.set<int>   ("io disp every nstep",sdyn.get<int>("RESEVRYDISP"));
  genalphaparams.set<bool>  ("ADAPTCONV",getIntegralValue<int>(sdyn,"ADAPTCONV")==1);
  genalphaparams.set<double>("ADAPTCONV_BETTER",sdyn.get<double>("ADAPTCONV_BETTER"));
  INPAR::STR::StressType iostress = Teuchos::getIntegralValue<INPAR::STR::StressType>(ioflags,"STRUCT_STRESS");
  genalphaparams.set<INPAR::STR::StressType>("io structural stress", iostress);
  genalphaparams.set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));
  INPAR::STR::StrainType iostrain = Teuchos::getIntegralValue<INPAR::STR::StrainType>(ioflags,"STRUCT_STRAIN");
  genalphaparams.set<INPAR::STR::StrainType>("io structural strain", iostrain);
  genalphaparams.set<bool>  ("io surfactant",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_SURFACTANT"));
  genalphaparams.set<int>   ("restart",probtype.get<int>("RESTART"));
  genalphaparams.set<int>   ("write restart every",sdyn.get<int>("RESTARTEVRY"));
  genalphaparams.set<bool>  ("print to screen",false);
  genalphaparams.set<bool>  ("print to err",true);
  genalphaparams.set<FILE*> ("err file",DRT::Problem::Instance()->ErrorFile()->Handle());
  genalphaparams.set<bool>  ("LOADLIN",false);
  INPAR::STR::ControlType controltype = Teuchos::getIntegralValue<INPAR::STR::ControlType>(sdyn,"CONTROLTYPE");
  genalphaparams.set<INPAR::STR::ControlType>("CONTROLTYPE",controltype);
  {
    vector<int> controlnode;
    std::istringstream contnode(Teuchos::getNumericStringParameter(sdyn,"CONTROLNODE"));
    std::string word;
    while (contnode >> word)
      controlnode.push_back(std::atoi(word.c_str()));
    if ((int)controlnode.size() != 3) dserror("Give proper values for CONTROLNODE in input file");
    genalphaparams.set("CONTROLNODE",controlnode[0]);
    genalphaparams.set("CONTROLDOF",controlnode[1]);
    genalphaparams.set("CONTROLCURVE",controlnode[2]);
  }
  genalphaparams.set<string>("equilibrium iteration","full newton");
  switch (Teuchos::getIntegralValue<INPAR::STR::PredEnum>(sdyn,"PREDICT"))
  {
    case INPAR::STR::pred_vague:
      dserror("You have to define the predictor");
      break;
    case INPAR::STR::pred_constdis:
      genalphaparams.set<string>("predictor","consistent");
      break;
    case INPAR::STR::pred_constdisvelacc:
      genalphaparams.set<string>("predictor","constant");
      break;
    case INPAR::STR::pred_tangdis:
      genalphaparams.set<string>("predictor","tangdis");
      break;
    default:
      dserror("Cannot cope with choice of predictor");
      break;
  }
  
  RCP<StruGenAlpha> tintegrator = rcp(new StruGenAlpha(genalphaparams,*discret_,*solver_,*output_));
  int    step    = genalphaparams.get<int>("step",0);
  int    nstep   = genalphaparams.get<int>("nstep",-1);
  //double time    = genalphaparams.get<double>("total time",0.0);
  double maxtime = genalphaparams.get<double>("max time",0.0);
  string pred    = genalphaparams.get<string>("predictor","constant");
  string equil   = genalphaparams.get<string>("equilibrium iteration","full newton");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  Epetra_SerialDenseVector cvector;
  if (!myrank) cvector.Size(nmp_);

  // load controled Newton only
  for (int i=step; i<nstep; ++i)
  {
    if      (predictor==1) tintegrator->ConstantPredictor();
    else if (predictor==2) tintegrator->ConsistentPredictor();
    tintegrator->FullNewton();
    tintegrator->Update();
    if (outputtofile) tintegrator->Output();
    if (!myrank) printf("Step %d\n",i);
    Epetra_SerialDenseVector cvector_arg = GetCalculatedCurve(*(tintegrator->Disp()));
    if (!myrank)
      for (int j=0; j<ndofs_; ++j)
        cvector[i*ndofs_+j] = cvector_arg[j];

    double time = genalphaparams.get<double>("total time",0.0);
    if (time>=maxtime) break;
  }
  return cvector;
}


/*----------------------------------------------------------------------*/
/* */
Epetra_SerialDenseVector STR::GenInvAnalysis::GetCalculatedCurve(Epetra_Vector& disp)
{
  int myrank = discret_->Comm().MyPID();
  
  // values observed at current time step
  Epetra_SerialDenseVector lcvector_arg(ndofs_);
  Epetra_SerialDenseVector gcvector_arg(ndofs_);
  int count=0;
  for (int i=0; i<nnodes_; ++i)
  {
    int gnode = nodes_[i];
    if (!discret_->NodeRowMap()->MyGID(gnode)) 
    {
      count += (int)dofs_[i].size();
      continue;
    }
    DRT::Node* node = discret_->gNode(gnode);
    if (!node) dserror("Cannot find node on this proc");
    if (myrank != node->Owner()) 
    {
      count += (int)dofs_[i].size();
      continue;
    }
    for (int j=0; j<(int)dofs_[i].size(); ++j)
    {
      int ldof = dofs_[i][j];
      int gdof = discret_->Dof(node,ldof);
      if (!disp.Map().MyGID(gdof)) dserror("Cannot find dof on this proc");
      lcvector_arg[count+j] = disp[disp.Map().LID(gdof)];
    }
    count += (int)dofs_[i].size();
  }
  discret_->Comm().SumAll(&lcvector_arg[0],&gcvector_arg[0],ndofs_);
  
  return gcvector_arg;
}

/*----------------------------------------------------------------------*/
/* */
void STR::GenInvAnalysis::PrintStorage(Epetra_SerialDenseVector delta_p)
{
  int myrank = discret_->Comm().MyPID();
  // store the error and mu_

  // storing current material parameters
  p_s_.Reshape(numb_run_+1,  np_);
  for (int i=0; i<np_; i++)
    p_s_(numb_run_, i)=p_(i);

  // storing current material parameter increments
  delta_p_s_.Reshape(numb_run_+1,  np_);
  for (int i=0; i<np_; i++)
    delta_p_s_(numb_run_, i)=delta_p(i);

  mu_s_.Resize(numb_run_+1);
  mu_s_(numb_run_)=mu_;

  error_s_.Resize(numb_run_+1);
  error_s_(numb_run_) = error_;
  // print error and parameter

  if (!myrank)
  {
    cout << endl;
    printf("#################################################################\n");
    printf("############################ Inverse Analysis run ## %3i ########\n",  numb_run_);
    printf("#################################################################\n");

    for (int i=0; i < numb_run_+1; i++)
    {
      printf("Error    : ");
      printf("%10.5e\n", error_s_(i));

      printf("Parameter: ");
      for (int j=0; j < delta_p.Length(); j++)
        printf("%10.5e ", p_s_(i, j));
      printf("\n");

      printf("mu       : ");
      printf("%10.5e", mu_s_(i));
      printf("\n");
      printf("#################################################################\n");
    }
  }
  return;
}


//-----------------------------------------------------------------------------------
void STR::GenInvAnalysis::ReadInParameters()
{
  int myrank = discret_->Comm().MyPID();
  
  // loop all materials in problem
  map<int,RCP<MAT::PAR::Material> >& mats = const_cast<map<int,RCP<MAT::PAR::Material> >& >
                                        (*(DRT::Problem::Instance()->Materials()->Map()));
  
  if (!myrank) printf("No. material laws considered : %d\n",mats.size());
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        int j = p_.Length();
        p_.Resize(j+2);
        p_(j)   = params->youngs_;
        p_(j+1) = params->beta_;
        //p_(j+2) = params->nue_; // need also change resize above to invoke nue
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i=0; i<nummat; ++i)
        {
          int id = (*matids)[i];
          RCP<MAT::PAR::Material> actelastmat = mats[id];
          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              int j = p_.Length();
              p_.Resize(j+3);
              p_(j)   = params->c1_;
              p_(j+1) = params->c2_;
              p_(j+2) = params->c3_;
            }
            break;
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              int j = p_.Length();
              p_.Resize(j+1);
              p_(j)   = params->kappa_;
              // p_(j+1) = params->beta_; // need also change resize above to invoke beta_
            }
            break;
            default:
              dserror("Unknown type of elasthyper material");
            break;
          }
        }
      }
      case INPAR::MAT::mes_isoyeoh: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      case INPAR::MAT::mes_vologden: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      default:
        // ignore unknown materials ?
        dserror("Unknown type of material");
      break;
    }
  }
  return;
}
//--------------------------------------------------------------------------------------
void STR::GenInvAnalysis::SetParameters(Epetra_SerialDenseVector p_cur)
{
  int myrank = discret_->Comm().MyPID();

  // parameters are evaluated on proc 0 only? This should not be neccessary....
  discret_->Comm().Broadcast(&p_cur[0],np_,0);

  // loop all materials in problem
  map<int,RCP<MAT::PAR::Material> >& mats = const_cast<map<int,RCP<MAT::PAR::Material> >& >
                                        (*(DRT::Problem::Instance()->Materials()->Map()));
  int count=0;
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        // This is a tiny little bit brutal!!!
        const_cast<double&>(params->youngs_) = p_cur[count];
        const_cast<double&>(params->beta_)   = p_cur[count+1];
        //const_cast<double&>(params->nue_)    = p_cur[count+2];
        if (!myrank) printf("MAT::PAR::AAAneohooke %20.15e %20.15e\n",p_cur[count],p_cur[count+1]);
        count += 2;
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i=0; i<nummat; ++i)
        {
          int id = (*matids)[i];
          RCP<MAT::PAR::Material> actelastmat = mats[id];
          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const_cast<double&>(params->c1_) = p_cur[count];
              const_cast<double&>(params->c2_) = p_cur[count+1];
              const_cast<double&>(params->c3_) = p_cur[count+2];
              if (!myrank) printf("MAT::ELASTIC::PAR::IsoYeoh %20.15e %20.15e %20.15e\n",params->c1_,params->c2_,params->c3_);
              count += 3;
            }
            break;
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const_cast<double&>(params->kappa_) = p_cur[count];
              //const_cast<double&>(params->beta_) = p_cur[count+1];
              if (!myrank) printf("MAT::ELASTIC::PAR::VolOgden %20.15e %20.15e\n",params->kappa_,params->beta_);
              count += 1;
            }
            break;
            default:
              dserror("Unknown type of elasthyper material");
            break;
          }
        }
      }
      break;
      case INPAR::MAT::mes_isoyeoh: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      case INPAR::MAT::mes_vologden: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      default:
        // ignore unknown materials ?
        dserror("Unknown type of material");
      break;
    }
  }
  

  return;
}



#endif
