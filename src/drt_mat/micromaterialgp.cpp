#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "micromaterialgp.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"

#include "../drt_stru_multi/microstrugenalpha.H"
#include "../io/io_drt_micro.H"

using namespace std;
using namespace Teuchos;
using namespace IO;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | defined in global_control.c                                          |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA  *alldyn;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
//extern struct _IO_FLAGS     ioflags;

extern struct _MATERIAL    *mat;


RefCountPtr<MicroStruGenAlpha> MAT::MicroMaterialGP::microgenalpha_;


/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterialGP::MicroMaterialGP(const int gp, const int ele_ID)
  : gp_(gp),
    ele_ID_(ele_ID)
{
  RefCountPtr<DRT::Problem> microproblem = DRT::Problem::Instance(1);
  RefCountPtr<DRT::Discretization> microdis = microproblem->Dis(0, 0);
  disp_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  vel_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
  acc_ = LINALG::CreateVector(*microdis->DofRowMap(),true);
}

/// destructor

MAT::MicroMaterialGP::~MicroMaterialGP()
{ }


/// Set up microscale generalized alpha

void MAT::MicroMaterialGP::SetUpMicroGenAlpha()
{
  //cout << "SetUpMicroGenAlpha\n";

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[0];
  STRUCT_DYNAMIC* sdyn     = alldyn[genprob.numsf].sdyn;

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> genalphaparams = rcp(new ParameterList());
  MicroStruGenAlpha::SetDefaults(*genalphaparams);

  genalphaparams->set<bool>  ("damping",sdyn->damp);
  genalphaparams->set<double>("damping factor K",sdyn->k_damp);
  genalphaparams->set<double>("damping factor M",sdyn->m_damp);

  genalphaparams->set<double>("beta",sdyn->beta);
  genalphaparams->set<double>("gamma",sdyn->gamma);
  genalphaparams->set<double>("alpha m",sdyn->alpha_m);
  genalphaparams->set<double>("alpha f",sdyn->alpha_f);

  genalphaparams->set<double>("total time",0.0);
  genalphaparams->set<double>("delta time",sdyn->dt);
  genalphaparams->set<int>   ("step",0);
  genalphaparams->set<int>   ("nstep",sdyn->nstep);
  genalphaparams->set<int>   ("max iterations",sdyn->maxiter);
  genalphaparams->set<int>   ("num iterations",-1);
  genalphaparams->set<double>("tolerance displacements",sdyn->toldisp);

  // takes values "full newton" , "modified newton" , "nonlinear cg"
  genalphaparams->set<string>("equilibrium iteration","full newton");

  // takes values "constant" "consistent"
  genalphaparams->set<string>("predictor","constant");

  microgenalpha_ = rcp(new MicroStruGenAlpha(genalphaparams,actdis,solver));
}


/// perform microscale simulation

void MAT::MicroMaterialGP::PerformMicroSimulation(const Epetra_SerialDenseMatrix* defgrd,
                                                  Epetra_SerialDenseVector* stress,
                                                  Epetra_SerialDenseMatrix* cmat,
                                                  double* density,
                                                  const double time)
{
  //cout << "PerformMicroSimulation\n";

  // if MicroDiscretizationWriter is not yet initialized -> set up
  // and write mesh

  if (micro_output_ == null)
  {
    RefCountPtr<DRT::Discretization> actdis = DRT::Problem::Instance(1)->Dis(genprob.numsf,0);

    micro_output_ = rcp(new MicroDiscretizationWriter(actdis, 1, ele_ID_, gp_));

    micro_output_->WriteMesh(0, 0.);

    // initialize total time, time step number and set time step size

    STRUCT_DYNAMIC* sdyn     = alldyn[genprob.numsf].sdyn;
    timen_ = 0.;
    istep_ = 0;
    dt_    = sdyn->dt;
  }

  // if derived generalized alpha class for microscale simulations is
  // not yet initialized -> set up

  if (microgenalpha_ == null)
  {
    MAT::MicroMaterialGP::SetUpMicroGenAlpha();
  }

  // set displacements, velocities and accelerations from last time step
  microgenalpha_->SetOldState(disp_, vel_, acc_);

  // check if we have to update absolute time and step number
  if (time != timen_)
  {
    if (timen_ != 0.)              // we don't want to write results
                                   // after just having constructed
                                   // the StruGenAlpha class
      microgenalpha_->Output(micro_output_, timen_, istep_);
    timen_ = time;
    istep_++;
  }

  // set current absolute time and step number
  microgenalpha_->SetTime(timen_, istep_);

  microgenalpha_->ConstantPredictor();

  // set boundary conditions derived from macroscale
  microgenalpha_->EvaluateMicroBC(defgrd);

  microgenalpha_->FullNewton();
  microgenalpha_->Update();

  // save calculated displacements, velocities and accelerations
  disp_ = microgenalpha_->ReturnNewDisp();
  vel_  = microgenalpha_->ReturnNewVel();
  acc_  = microgenalpha_->ReturnNewAcc();

  //cout << "current displacements in microscale of gp " << gp_ << ":\n" << *disp_ << "\n";

  // clear displacements in MicroStruGenAlpha for next usage
  microgenalpha_->ClearState();

  MAT::MicroMaterialGP::Homogenization(stress, cmat, density, defgrd);
}


/// determine macroscopic parameters via averaging (homogenization) of
/// microscopic features

void MAT::MicroMaterialGP::Homogenization(Epetra_SerialDenseVector* stress,
                                          Epetra_SerialDenseMatrix* cmat,
                                          double *density,
                                          const Epetra_SerialDenseMatrix* defgrd)
{
  // deformation gradient is only needed for test purposes and can be
  // eliminated later on for real calculations
  //cout << "In homogenization\n";

  double Emod = 100.0;    // Young's modulus (modulus of elasticity)
  double nu = 0.0;      // Poisson's ratio (Querdehnzahl)
  (*density) = 1.0;     // density, returned to evaluate mass matrix
  double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */

  Epetra_SerialDenseMatrix cauchygreen(3,3);
  cauchygreen.Multiply('T','N',1.0,*defgrd,*defgrd,1.0);

  Epetra_SerialDenseVector glstrain(6);
  glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
  glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
  glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
  glstrain(3) = cauchygreen(0,1);
  glstrain(4) = cauchygreen(1,2);
  glstrain(5) = cauchygreen(2,0);


  /* write non-zero components */
  (*cmat)(0,0) = mfac*(1.0-nu);
  (*cmat)(0,1) = mfac*nu;
  (*cmat)(0,2) = mfac*nu;
  (*cmat)(1,0) = mfac*nu;
  (*cmat)(1,1) = mfac*(1.0-nu);
  (*cmat)(1,2) = mfac*nu;
  (*cmat)(2,0) = mfac*nu;
  (*cmat)(2,1) = mfac*nu;
  (*cmat)(2,2) = mfac*(1.0-nu);
  /* ~~~ */
  (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);

  // evaluate stresses
  (*cmat).Multiply('N',glstrain,(*stress));

  //cout << "stresses: \n" << *stress << "\n";
}


#if 0
/// test routine for calculating stresses, constitutive matrix and density in
/// case of St Venant Kirchhoff material

void MAT::MicroMaterialGP::CalcStressStiffDens (Epetra_SerialDenseVector* stress,
                                                Epetra_SerialDenseMatrix* cmat,
                                                double* density)
{
  double Emod = 1.0;    // Young's modulus (modulus of elasticity)
  double nu = 0.0;      // Poisson's ratio (Querdehnzahl)
  (*density) = 0.0;     // density, returned to evaluate mass matrix
  double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
  /* write non-zero components */
  (*cmat)(0,0) = mfac*(1.0-nu);
  (*cmat)(0,1) = mfac*nu;
  (*cmat)(0,2) = mfac*nu;
  (*cmat)(1,0) = mfac*nu;
  (*cmat)(1,1) = mfac*(1.0-nu);
  (*cmat)(1,2) = mfac*nu;
  (*cmat)(2,0) = mfac*nu;
  (*cmat)(2,1) = mfac*nu;
  (*cmat)(2,2) = mfac*(1.0-nu);
  /* ~~~ */
  (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);

  // evaluate stresses
  (*cmat).Multiply('N',(*glstrain),(*stress));   // sigma = C . epsilon
}
#endif

#endif
#endif
#endif
