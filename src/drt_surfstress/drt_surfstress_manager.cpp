/*---------------------------------------------------------------------*/
/*! \file

\brief surface stresses due to interfacial phenomena

\maintainer Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/


#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "drt_surfstress_manager.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"

/*-------------------------------------------------------------------*
 |  ctor (public)                                            lw 12/07|
 *-------------------------------------------------------------------*/
UTILS::SurfStressManager::SurfStressManager(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::ParameterList sdynparams, const std::string file_prefix)
    : discret_(discret)
{
  std::vector<DRT::Condition*> surfstresscond(0);
  discret_->GetCondition("SurfaceStress", surfstresscond);
  if (surfstresscond.size())
  {
    switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP"))
    {
      case INPAR::STR::dyna_genalpha:
        if (DRT::INPUT::IntegralValue<INPAR::STR::MidAverageEnum>(
                sdynparams.sublist("GENALPHA"), "GENAVG") != INPAR::STR::midavg_trlike)
          dserror("Surface stresses are only implemented for tr-like generalized alpha");
        break;
      default:
        dserror("Surface stresses are only implemented for tr-like generalized alpha");
        break;
    }

    havesurfstress_ = true;
    timen_ = 0.;

    surfrowmap_ = DRT::UTILS::ConditionElementMap(*discret, "SurfaceStress", false);
    Teuchos::RCP<Epetra_Map> surfcolmap =
        DRT::UTILS::ConditionElementMap(*discret, "SurfaceStress", true);

    // We start with interfacial area and concentration = 0.
    // This is wrong but does not make a difference
    // since we apply the equilibrium concentration gradually, thus we
    // do not need these history variables needed for the dynamic model
    // in the beginning.
    A_current_ = Teuchos::rcp(new Epetra_Vector(*surfcolmap, true));
    A_last_ = Teuchos::rcp(new Epetra_Vector(*surfcolmap, true));
    con_current_ = Teuchos::rcp(new Epetra_Vector(*surfcolmap, true));
    con_last_ = Teuchos::rcp(new Epetra_Vector(*surfcolmap, true));
    gamma_current_ = Teuchos::rcp(new Epetra_Vector(*surfcolmap, true));

    std::vector<std::string> conditions_to_copy;
    std::string condname = "SurfaceStress";
    conditions_to_copy.push_back(condname);

    Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase> discreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
    surfdiscret_ = discreator->CreateMatchingDiscretizationFromCondition(
        *discret_, condname, "boundary", "BELE3_3", conditions_to_copy);

    surfdiscret_->FillComplete(true, true, true);
    std::string outfile = file_prefix + "_" + condname;

    const int ndim = DRT::Problem::Instance()->NDim();
    Teuchos::RCP<IO::OutputControl> condcontrol =
        Teuchos::rcp(new IO::OutputControl(surfdiscret_->Comm(),
            "none",  // we do not have a problem type
            ShapeFunctionType::shapefunction_polynomial,
            "debug-output",  // no input file either
            outfile,         // an output file name is needed
            ndim,
            0,     // restart is meaningless here
            1000,  // we never expect to get 1000 iterations
            DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(), "OUTPUT_BIN")));
    surfdiscret_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(surfdiscret_)));
    surfoutput_ = surfdiscret_->Writer();
    surfoutput_->SetOutput(condcontrol);
    surfoutput_->WriteMesh(0, 0.0);
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

  Teuchos::RCP<Epetra_Vector> A_row = LINALG::CreateVector(*surfrowmap_, true);
  Teuchos::RCP<Epetra_Vector> con_row = LINALG::CreateVector(*surfrowmap_, true);
  Teuchos::RCP<Epetra_Vector> gamma_row = LINALG::CreateVector(*surfrowmap_, true);

  LINALG::Export(*A_current_, *A_row);
  LINALG::Export(*con_current_, *con_row);
  LINALG::Export(*gamma_current_, *gamma_row);

  surfoutput_->WriteVector("gamma", gamma_row, IO::elementvector);
  surfoutput_->WriteVector("Gamma", con_row, IO::elementvector);
  surfoutput_->WriteVector("Area", A_row, IO::elementvector);
}


/*-------------------------------------------------------------------*
 |  Read restart from binary file                            lw 12/08|
 *-------------------------------------------------------------------*/
void UTILS::SurfStressManager::ReadRestart(
    const int step, const std::string file_prefix, const bool serial)
{
  // figure out how the file we restart from is called
  std::string condname = "SurfaceStress";
  std::string restartname = file_prefix + "_" + condname;

  Teuchos::RCP<IO::InputControl> surfinpcontrol =
      Teuchos::rcp(new IO::InputControl(restartname, serial));

  IO::DiscretizationReader reader(surfdiscret_, surfinpcontrol, step);
  int rstep = reader.ReadInt("step");
  if (rstep != step) dserror("Time step on file not equal to given step");

  // The row map based vectors written in case of restart are exported
  // to column based ones needed for calculations again

  Teuchos::RCP<Epetra_Vector> A = LINALG::CreateVector(*surfrowmap_, true);
  Teuchos::RCP<Epetra_Vector> con = LINALG::CreateVector(*surfrowmap_, true);

  reader.ReadVector(A, "Area");
  reader.ReadVector(con, "Gamma");

  LINALG::Export(*A, *A_last_);
  LINALG::Export(*con, *con_last_);
}


/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| interfacial phenomena                                              |
*--------------------------------------------------------------------*/

void UTILS::SurfStressManager::EvaluateSurfStress(Teuchos::ParameterList& p,
    const Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<LINALG::SparseOperator> stiff)
{
  // action for elements
  p.set("action", "calc_surfstress_stiff");

  discret_->ClearState();
  discret_->SetState("new displacement", disn);

  double currtime = p.get<double>("total time", 0.);
  double dt = p.get<double>("delta time", 0.);
  if (fabs(currtime - timen_) >
      0.5 * dt)  // This is true for every new time step as well as for the first
                 // time step after reading restart, where step_ still equals 0.
                 // Note that currtime and timen_ actually _have_ to be equal unless one
                 // enters a new time step (because timen_ is always just a copy).
  {
    timen_ = currtime;
    p.set("newstep", true);
  }

  discret_->EvaluateCondition(
      p, stiff, Teuchos::null, fint, Teuchos::null, Teuchos::null, "SurfaceStress");

  return;
}

/*-------------------------------------------------------------------*
| (public)                                                   lw 03/08|
|                                                                    |
| update surface area and concentration                              |
*--------------------------------------------------------------------*/

void UTILS::SurfStressManager::Update()
{
  A_last_->Update(1.0, *A_current_, 0.0);
  con_last_->Update(1.0, *con_current_, 0.0);
}

/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level                                                   |
*--------------------------------------------------------------------*/

void UTILS::SurfStressManager::StiffnessAndInternalForces(const int curvenum, const double& A,
    const Teuchos::RCP<Epetra_SerialDenseVector> Adiff,
    const Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2, Epetra_SerialDenseVector& fint,
    Epetra_SerialDenseMatrix& K_surf, const int ID, const double time, const double dt,
    const int surface_flag, const double const_gamma, const double k1xC, const double k2,
    const double m1, const double m2, const double gamma_0, const double gamma_min,
    const double gamma_min_eq, const double con_quot_max, const double con_quot_eq,
    const bool newstep)
{
  double gamma, dgamma;
  int LID = A_last_->Map().LID(ID);
  dserror("Here some proper solution has to implemented!");
  double t_end = 23.0;  // TODO: fix that, previously asked for Curve().end() to be replaced by
                        // problem specific parameter -> already discussed with Lena

  if (surface_flag == 0)  // SURFACTANT
  {
    (*A_current_)[LID] = A;

    /*-----------calculation of current surface stress and its partial
     *-----------------derivative with respect to the interfacial area */

    if (time <= t_end) /* gradual application of surface stress */
    {
      // start with regime 3 (start with maximal concentration) Birzle 12/15
      // If one is interested just in the stationary process, use this version.
      // So no build up time to get in the stationary process is necassary.
      // For more information see Bachelor Thesis by Elias Lochner Page 19 and 20.
      gamma = gamma_min;
      dgamma = 0.;
      (*con_current_)[LID] = con_quot_max;
      (*gamma_current_)[LID] = gamma;

      // old version: start with regime 1 (start with equilibrium concentration)
      // this version is equal to the experiments of Otis
      // gamma = gamma_0-m1*con_quot_eq;
      // dgamma = 0.;
      // (*con_current_)[LID] = con_quot_eq;
      // (*gamma_current_)[LID] = gamma;
    }
    else
    {
      SurfactantModel(ID, gamma, dgamma, dt, k1xC, k2, m1, m2, gamma_0, gamma_min, gamma_min_eq,
          con_quot_max, newstep);
    }
  }

  else if (surface_flag == 1)  // SURFACE TENSION
  {
    gamma = const_gamma;
    dgamma = 0.;
  }

  else
    dserror("Surface flag not implemented");

  double curvefac = 1.;

  /*------------gradual application of surface stresses via time curve*/
  if (time <= t_end) curvefac = DRT::Problem::Instance()->Funct(curvenum).EvaluateTime(time);

  double ndof = Adiff->Length();

  for (int i = 0; i < ndof; ++i)
  {
    fint[i] = gamma * (*Adiff)[i] * curvefac;

    for (int j = 0; j < ndof; ++j)
    {
      K_surf(i, j) = (dgamma * (*Adiff)[i] * (*Adiff)[j] + gamma * (*Adiff2)(i, j)) * curvefac;
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
void UTILS::SurfStressManager::SurfactantModel(const int ID,  // (i) ID of surface condition
    double& gamma,                                            // (o) surface stress
    double& dgamma,                                           // (o) derivative of surface stress
    const double dt,                                          // (i) timestep size
    double k1xC,                // (i) adsorption coefficent k1 times bulk concentration C
    double k2,                  // (i) desorption coefficient
    const double m1,            // (i) 1st isothermal slope
    const double m2,            // (i) 2nd isothermal slope
    const double gamma_0,       // (i) surface tension of water
    const double gamma_min,     // (i) mimimum surface stress
    const double gamma_min_eq,  // (i) mimimum equilibrium surface stress
    const double con_max,       // (i) max. surfactant concentration
    const bool newstep)         // (i) flag for new step (predictor)
{
  const int LID = A_last_->Map().LID(ID);
  const double A_last = (*A_last_)[LID];
  const double A_new = (*A_current_)[LID];
  const double con_last = (*con_last_)[LID];
  const double M_last = con_last * A_last;  // mass of surfactant molecules (n)

  /*-------------------------------------------------------------------*
   |     DETERMINATION OF CURRENT NON-DIMENSIONALIZED INTERFACIAL      |
   |       SURFACTANT CONCENTRATION CON_QUOT_NEW  (n+1)                |
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
      k1xC *= (1. - con_last) / 0.02;
      k2 *= (1. - con_last) / 0.02;
    }

    con_new = (1.0 / dt * M_last + k1xC * A_new) / (A_new * (1.0 / dt + k1xC + k2));
    dcon_new = -M_last / (dt * A_new * A_new * (1 / dt + k1xC + k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_last < con_max)
    {
      con_new = M_last / A_new;
      dcon_new = -con_new / A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      // we need to check whether we just entered a new time step
      // because in this case, A_old equals A_new and we don't
      // want to base the decision on entering a new regime on this comparison!
      if (A_last < A_new && newstep == false)
      {
        con_new = M_last / A_new;
        dcon_new = -con_new / A_new;
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
    gamma = gamma_0 - m1 * con_new;
    dgamma = -m1 * dcon_new;
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_new < con_max)
    {
      gamma = gamma_min_eq - m2 * (con_new - 1.);
      dgamma = -m2 * dcon_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      gamma = gamma_min;
      dgamma = 0.;
    }
  }

  /* save history variables */
  (*con_current_)[LID] = con_new;
  (*gamma_current_)[LID] = gamma;

  return;
}


/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Write restart (current gamma, Gamma, A)                            |
*--------------------------------------------------------------------*/
void UTILS::SurfStressManager::WriteRestart(const int istep, const double timen)
{
  surfoutput_->WriteMesh(istep, timen);
  WriteResults(istep, timen);
}
