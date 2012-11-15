/*!-------------------------------------------------------------------
\file drt_surfstress_manager.cpp

\brief Class controlling surface stresses due to interfacial phenomena
and containing all necessary history data

<pre>
Maintainer: Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*--------------------------------------------------------------------*/


#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "drt_surfstress_manager.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"

/*-------------------------------------------------------------------*
 |  ctor (public)                                            lw 12/07|
 *-------------------------------------------------------------------*/
UTILS::SurfStressManager::SurfStressManager(Teuchos::RCP<DRT::Discretization> discret,
                                            ParameterList sdynparams,
                                            const std::string file_prefix):
  discret_(discret)
{
  std::vector<DRT::Condition*> surfstresscond(0);
  discret_->GetCondition("SurfaceStress", surfstresscond);
  if (surfstresscond.size())
  {
    switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams,"DYNAMICTYP"))
    {
    case INPAR::STR::dyna_genalpha:
      if (DRT::INPUT::IntegralValue<INPAR::STR::MidAverageEnum>(sdynparams.sublist("GENALPHA"), "GENAVG")
             != INPAR::STR::midavg_imrlike)
        dserror("Surface stresses are only implemented for imr-like generalized alpha");
      alphaf_ = sdynparams.sublist("GENALPHA").get<double>("ALPHA_F");
      break;
    default:
      dserror("Surface stresses are only implemented for imr-like generalized alpha");
      break;
    }

    havesurfstress_ = true;
    timen_ = 0.;

    surfrowmap_ = DRT::UTILS::ConditionElementMap(*discret, "SurfaceStress", false);
    Teuchos::RCP<Epetra_Map> surfcolmap = DRT::UTILS::ConditionElementMap(*discret, "SurfaceStress", true);

    // We start with interfacial area and concentration = 0.
    // This is wrong but does not make a difference
    // since we apply the equilibrium concentration gradually, thus we
    // do not need these history variables needed for the dynamic model
    // in the beginning.
    A_current_       = Teuchos::rcp(new Epetra_Vector(*surfcolmap,true));
    A_last_          = Teuchos::rcp(new Epetra_Vector(*surfcolmap,true));
    con_current_     = Teuchos::rcp(new Epetra_Vector(*surfcolmap,true));
    con_last_        = Teuchos::rcp(new Epetra_Vector(*surfcolmap,true));
    gamma_current_   = Teuchos::rcp(new Epetra_Vector(*surfcolmap,true));
    gamma_last_      = Teuchos::rcp(new Epetra_Vector(*surfcolmap,true));


    std::vector<string> conditions_to_copy;
    std::string condname = "SurfaceStress";
    conditions_to_copy.push_back(condname);

    surfdiscret_ = DRT::UTILS::CreateDiscretizationFromCondition(discret_, condname,
                                                                 "boundary", "BELE3",
                                                                 conditions_to_copy);
    surfdiscret_->FillComplete();
    std::string outfile = file_prefix + "_" + condname;

    const int ndim = DRT::Problem::Instance()->NDim();
    Teuchos::RCP<IO::OutputControl> condcontrol = Teuchos::rcp(new IO::OutputControl(surfdiscret_->Comm(),
                                                                                     "none",                   // we do not have a problem type
                                                                                     "Polynomial",
                                                                                     "debug-output",           // no input file either
                                                                                     outfile,                  // an output file name is needed
                                                                                     ndim,
                                                                                     0,                        // restart is meaningless here
                                                                                     1000));                   // we never expect to get 1000 iterations

    surfoutput_=Teuchos::rcp(new IO::DiscretizationWriter(surfdiscret_,condcontrol));
    surfoutput_->WriteMesh(0,0.0);

  }
  else
    havesurfstress_ = false;
}


/*-------------------------------------------------------------------*
 |  Write Results to file                                    lw 12/08|
 *-------------------------------------------------------------------*/
void UTILS::SurfStressManager::WriteResults(const int istep, const double timen)
{
  surfoutput_->NewStep(istep, timen);

  // The column map based vectors used for calculations are exported
  // to row map based ones needed for writing

  RCP<Epetra_Vector> A_last_row = LINALG::CreateVector(*surfrowmap_,true);
  RCP<Epetra_Vector> con_last_row = LINALG::CreateVector(*surfrowmap_,true);
  RCP<Epetra_Vector> gamma_last_row = LINALG::CreateVector(*surfrowmap_,true);

  LINALG::Export(*A_last_, *A_last_row);
  LINALG::Export(*con_last_, *con_last_row);
  LINALG::Export(*gamma_last_, *gamma_last_row);

  surfoutput_->WriteVector("gamma", gamma_last_row, IO::DiscretizationWriter::elementvector);
  surfoutput_->WriteVector("Gamma", con_last_row, IO::DiscretizationWriter::elementvector);
  surfoutput_->WriteVector("Area", A_last_row, IO::DiscretizationWriter::elementvector);
}


/*-------------------------------------------------------------------*
 |  Read restart from binary file                            lw 12/08|
 *-------------------------------------------------------------------*/
void UTILS::SurfStressManager::ReadRestart(const int step,
                                           const std::string file_prefix,
                                           const bool serial)
{
  // figure out how the file we restart from is called
  std::string condname = "SurfaceStress";
  std::string restartname = file_prefix + "_" + condname;

  Teuchos::RCP<IO::InputControl> surfinpcontrol = Teuchos::rcp(new IO::InputControl(restartname, serial));

  IO::DiscretizationReader reader(surfdiscret_, surfinpcontrol, step);
  int    rstep = reader.ReadInt("step");
  if (rstep != step) dserror("Time step on file not equal to given step");

  // The row map based vectors written in case of restart are exported
  // to column based ones needed for calculations again

  RCP<Epetra_Vector> A = LINALG::CreateVector(*surfrowmap_,true);
  RCP<Epetra_Vector> con = LINALG::CreateVector(*surfrowmap_,true);
  RCP<Epetra_Vector> gamma = LINALG::CreateVector(*surfrowmap_,true);

  reader.ReadVector(A, "Area");
  reader.ReadVector(con, "Gamma");
  reader.ReadVector(gamma, "gamma");

  LINALG::Export(*A, *A_last_);
  LINALG::Export(*con, *con_last_);
  LINALG::Export(*gamma, *gamma_last_);
}


/*-------------------------------------------------------------------*
| (public)						     lw 12/07|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| interfacial phenomena                                              |
*--------------------------------------------------------------------*/

void UTILS::SurfStressManager::EvaluateSurfStress(ParameterList& p,
                                                  const RCP<Epetra_Vector> dism,
                                                  const RCP<Epetra_Vector> disn,
                                                  RCP<Epetra_Vector> fint,
                                                  RCP<LINALG::SparseOperator> stiff)
{
  // action for elements
  p.set("action","calc_surfstress_stiff");

  discret_->ClearState();
  discret_->SetState("displacement",dism);
  discret_->SetState("new displacement",disn);

  double currtime = p.get<double>("total time", 0.);
  double dt = p.get<double>("delta time", 0.);
  if (fabs(currtime-timen_)>0.5*dt)  // This is true for every new time step as well as for the first
                                     // time step after reading restart, where step_ still equals 0.
                                     // Note that currtime and timen_ actually _have_ to be equal unless one
                                     // enters a new time step (because timen_ is always just a copy).
  {
    timen_ = currtime;
    p.set("newstep", true);
  }

  discret_->EvaluateCondition(p,stiff,null,fint,null,null,"SurfaceStress");

  return;
}

/*-------------------------------------------------------------------*
| (public)						     lw 03/08|
|                                                                    |
| update surface area and concentration                              |
*--------------------------------------------------------------------*/

void UTILS::SurfStressManager::Update()
{
  A_last_->Update(1.0, *A_current_, 0.0);
  con_last_->Update(1.0, *con_current_, 0.0);
  gamma_last_->Update(1.0, *gamma_current_, 0.0);
}

/*-------------------------------------------------------------------*
| (public)						     lw 12/07|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level                                                   |
*--------------------------------------------------------------------*/

void UTILS::SurfStressManager::StiffnessAndInternalForces(const int curvenum,
                                                          const double& A,
                                                          const RCP<Epetra_SerialDenseVector> Adiff,
                                                          const RCP<Epetra_SerialDenseMatrix> Adiff2,
                                                          const double& Anew,
                                                          const RCP<Epetra_SerialDenseVector> Adiffnew,
                                                          Epetra_SerialDenseVector& fint,
                                                          Epetra_SerialDenseMatrix& K_surf,
                                                          const int ID,
                                                          const double time,
                                                          const double dt,
                                                          const int surface_flag,
                                                          const double const_gamma,
                                                          const double k1xC,
                                                          const double k2,
                                                          const double m1,
                                                          const double m2,
                                                          const double gamma_0,
                                                          const double gamma_min,
                                                          const double gamma_min_eq,
                                                          const double con_quot_max,
                                                          const double con_quot_eq,
                                                          const bool newstep)
{
  double gamma, dgamma;
  int LID = A_last_->Map().LID(ID);
  double t_end = DRT::Problem::Instance()->Curve(curvenum).end();

  if (surface_flag==0)                                   // SURFACTANT
  {
    (*A_current_)[LID] = Anew;

    /*-----------calculation of current surface stress and its partial
     *-----------------derivative with respect to the interfacial area */

    if (time <= t_end)         /* gradual application of surface stress */
    {
      gamma = gamma_0-m1*con_quot_eq;
      dgamma = 0.;
      (*con_current_)[LID] = con_quot_eq;
      (*gamma_current_)[LID] = gamma;
    }
    else
    {
      SurfactantModel(ID, gamma, dgamma, dt, k1xC, k2, m1, m2,
                      gamma_0, gamma_min, gamma_min_eq, con_quot_max, newstep);
    }
  }

  else if (surface_flag==1)                             // SURFACE TENSION
  {
    gamma = const_gamma;
    dgamma = 0.;
  }

  else
    dserror("Surface flag not implemented");

  double curvefac = 1.;

  /*------------gradual application of surface stresses via time curve*/
  if (time <= t_end)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  double ndof = Adiff->Length();

  for (int i=0;i<ndof;++i)
  {
    fint[i] = gamma*(*Adiff)[i]*curvefac;

    for (int j=0;j<ndof;++j)
    {
      K_surf(i,j) = (dgamma*(*Adiff)[i]*(*Adiffnew)[j]+gamma*(*Adiff2)(i,j))*curvefac;
    }
  }

  return;
}


/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Calculate current interfacial concentration of surfactant          |
| molecules and corresponding surface stresses                       |
*--------------------------------------------------------------------*/
void UTILS::SurfStressManager::SurfactantModel(
                const int ID,                // (i) ID of surface condition
                double& gamma,               // (o) surface stress
                double& dgamma,              // (o) derivative of surface stress
                const double dt,             // (i) timestep size
                double k1xC,                 // (i) adsorption coefficent k1 times bulk concentration C
                double k2,                   // (i) desorption coefficient
                const double m1,             // (i) 1st isothermal slope
                const double m2,             // (i) 2nd isothermal slope
                const double gamma_0,        // (i) surface tension of water
                const double gamma_min,      // (i) mimimum surface stress
                const double gamma_min_eq,   // (i) mimimum equilibrium surface stress
                const double con_max,        // (i) max. surfactant concentration
                const bool newstep)          // (i) flag for new step (predictor)
{
  const int LID = A_last_->Map().LID(ID);
  const double A_last = (*A_last_)[LID];
  const double A_new = (*A_current_)[LID];
  const double con_last = (*con_last_)[LID];
  const double M_last = con_last*A_last;     // mass of surfactant molecules (n)

  /*-------------------------------------------------------------------*
   |     DETERMINATION OF CURRENT NON-DIMENSIONALIZED INTERFACIAL      |
   |       SURFACTANT CONCENTRATION CON_QUOT_NEW  (n+1-alphaf)         |
   *-------------------------------------------------------------------*/

  /* use of non-dimensionalized interfacial surfactant concentration:
   * con_quot=Gamma/Gamma_min_eq
   * mass of surfactant molecules M_old is also divided
   * by minimum equilibrium interfacial surfactant concentration!*/
  double con_new = 0.;
  double dcon_new = 0.;

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_last < 1.)
  {

    /* assumed continuous drop of adsorption/desorption coefficients
     * when approaching insoluble regime */

    if (con_last > 0.98)
    {
      k1xC *= (1.-con_last)/0.02;
      k2 *= (1.-con_last)/0.02;
    }

    con_new = (1.0/dt*M_last+k1xC*A_new)/(A_new*(1.0/dt+k1xC+k2));
    dcon_new = -M_last/(dt*A_new*A_new*(1/dt+k1xC+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_last < con_max)
    {
      con_new = M_last/A_new;
      dcon_new = -con_new/A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      // we need to check whether we just entered a new time step
      // because in this case, A_old equals A_new and we don't
      // want to base the decision on entering a new regime on this comparison!
      if (A_last < A_new && newstep == false)
      {
        con_new = M_last/A_new;
        dcon_new = -con_new/A_new;
      }
      else
      {
        con_new = con_max;
        dcon_new = 0.;
      }
    }
  }

  /* interfacial surfactant concentration must not be larger than the
   * maximum interfacial concentration */
  if (con_new > con_max)
  {
    con_new = con_max;
  }

  /*-------------------------------------------------------------------*
   |        DETERMINATION OF CURRENT GENERALIZED SURFACE ENERGY        |
   *-------------------------------------------------------------------*/

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_new < 1.)
  {
    gamma  = gamma_0-m1*con_new;
    dgamma = -m1*dcon_new;
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_new < con_max)
    {
      gamma  = gamma_min_eq-m2*(con_new-1.);
      dgamma = -m2*dcon_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      gamma  = gamma_min;
      dgamma = 0.;
    }
  }

  /* save history variables */
  (*con_current_)[LID] = con_new;
  (*gamma_current_)[LID] = gamma;

  gamma = (1.-alphaf_)*(*gamma_current_)[LID]+alphaf_*(*gamma_last_)[LID];

  return;
}


/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Write restart (current gamma, Gamma, A)                            |
*--------------------------------------------------------------------*/
void UTILS::SurfStressManager::WriteRestart(const int istep, const double timen)
{
  surfoutput_->WriteMesh(istep,timen);
  WriteResults(istep, timen);
}


