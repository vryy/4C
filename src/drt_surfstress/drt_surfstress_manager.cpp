/*!-------------------------------------------------------------------
\file drt_surfstress_manager.cpp

\brief Class controlling surface stresses due to interfacial phenomena
and containing all necessary history data

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "drt_surfstress_manager.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include <cstdlib>

/*-------------------------------------------------------------------*
 |  ctor (public)                                            lw 12/07|
 *-------------------------------------------------------------------*/
DRT::SurfStressManager::SurfStressManager(DRT::Discretization& discret):
discret_(discret)
{
  surfrowmap_ = DRT::UTILS::GeometryElementMap(discret, "SurfaceStress", false);
  RCP<Epetra_Map> surfcolmap = DRT::UTILS::GeometryElementMap(discret, "SurfaceStress", true);

  // we start with interfacial area (A_old_) and concentration
  // (con_quot_) = 0. this is wrong but does not make a difference
  // since we apply the equilibrium concentration gradually, thus we
  // do not need these history variables needed for the dynamic model
  // in the beginning.
  A_old_temp_    = rcp(new Epetra_Vector(*surfcolmap,true));
  A_old_         = rcp(new Epetra_Vector(*surfcolmap,true));
  con_quot_temp_ = rcp(new Epetra_Vector(*surfcolmap,true));
  con_quot_      = rcp(new Epetra_Vector(*surfcolmap,true));
}

/*-------------------------------------------------------------------*
| (public)						     lw 12/07|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| interfacial phenomena                                              |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::EvaluateSurfStress(ParameterList& p,
                                                RefCountPtr<Epetra_Vector> disp,
                                                RefCountPtr<Epetra_Vector> fint,
                                                RefCountPtr<LINALG::SparseMatrix> stiff)
{
  // action for elements
  p.set("action","calc_surfstress_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  discret_.EvaluateCondition(p,stiff,null,fint,null,null,"SurfaceStress");

  return;
}

/*-------------------------------------------------------------------*
| (public)						     lw 03/08|
|                                                                    |
| update surface area and concentration                              |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::Update()
{
  A_old_->Update(1.0, *A_old_temp_, 0.0);
  con_quot_->Update(1.0, *con_quot_temp_, 0.0);
}

/*-------------------------------------------------------------------*
| (public)						     lw 03/08|
|                                                                    |
| write restart                                                      |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::GetHistory(RCP<Epetra_Vector> A_old_row,
                                        RCP<Epetra_Vector> con_quot_row)
{
  LINALG::Export(*A_old_, *A_old_row);
  LINALG::Export(*con_quot_, *con_quot_row);
}

/*-------------------------------------------------------------------*
| (public)						     lw 12/07|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level                                                   |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::StiffnessAndInternalForces(const int curvenum,
                                                        const double& A,
                                                        const Epetra_SerialDenseVector& Adiff,
                                                        const Epetra_SerialDenseMatrix& Adiff2,
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
                                                        const double alphaf,
                                                        const bool newstep)
{
  double gamma, dgamma;
  int LID = A_old_->Map().LID(ID);
  double t_end = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).end();

  if (surface_flag==0)                                   // SURFACTANT
  {
    if ((*A_old_)[LID] == 0.)
      (*A_old_)[LID] = A;

    /* update for imr-like gen-alpha*/
    (*A_old_temp_)[LID] = (A - alphaf*(*A_old_)[LID])/(1.-alphaf);

    /*-----------calculation of current surface stress and its partial
     *-----------------derivative with respect to the interfacial area */

    if (time <= t_end)         /* gradual application of surface stress */
    {
      gamma = gamma_0-m1*con_quot_eq;
      dgamma = 0.;
      (*con_quot_temp_)[LID] = con_quot_eq;
    }
    else
    {
      SurfactantModel(ID, gamma, dgamma, dt, k1xC, k2, m1, m2,
                      gamma_0, gamma_min, gamma_min_eq, con_quot_max, alphaf, newstep);
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
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  double ndof = Adiff.Length();

  for (int i=0;i<ndof;++i)
  {
    for (int j=0;j<ndof;++j)
    {
      K_surf(i,j) = (dgamma*Adiff[i]*Adiff[j]+gamma*Adiff2(i,j))*curvefac;
    }
  }

  /*------calculation of current internal force due to surface energy*/
  for (int i=0;i<ndof;++i)
  {
    fint[i] = gamma*Adiff[i]*curvefac;
  }

  return;
}


/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Calculate current interfacial concentration of surfactant          |
| molecules and corresponding surface stresses                       |
*--------------------------------------------------------------------*/
void DRT::SurfStressManager::SurfactantModel(
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
                const double con_quot_max,   // (i) max. surfactant concentration
                const double alphaf,         // (i) generalized-alpha parameter alphaf
                const bool newstep)          // (i) flag for new step (predictor)
{
  const int LID = A_old_->Map().LID(ID);
  const double A_old = (*A_old_)[LID];
  const double A_new = (*A_old_temp_)[LID];

  /*-------------------------------------------------------------------*
   |     DETERMINATION OF CURRENT NON-DIMENSIONALIZED INTERFACIAL      |
   |       SURFACTANT CONCENTRATION CON_QUOT_NEW  (n+1-alphaf)         |
   *-------------------------------------------------------------------*/

  /* use of non-dimensionalized interfacial surfactant concentration:
   * con_quot=Gamma/Gamma_min_eq
   * mass of surfactant molecules M_old is also divided
   * by minimum equilibrium interfacial surfactant concentration!*/

  const double con_quot_old = (*con_quot_)[LID];
  const double M_old = con_quot_old*A_old;       // mass of surfactant molecules of last time step
  double con_quot_alphaf = 0.;
  double con_quot_new = 0.;

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_quot_old < 1.)
  {

    /* assumed continuous drop of adsorption/desorption coefficients
     * when approaching insoluble regime */

    if (con_quot_old > 0.98)
    {
      k1xC *= (1.-con_quot_old)/0.02;
      k2 *= (1.-con_quot_old)/0.02;
    }

    con_quot_new = (1.0/dt*M_old+k1xC*A_new)/(A_new*(1.0/dt+k1xC+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_quot_old < con_quot_max)
    {
      con_quot_new = M_old/A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      // we need to check whether we just entered a new time step
      // because in this case, A_old nearly equals A_new and we don't
      // want to base the decision on entering a new regime on this comparison!
      if (A_old < A_new && newstep == false)
      {
        con_quot_new = M_old/A_new;
      }
      else
      {
        con_quot_new = con_quot_max;
      }
    }
  }

  // for imr-like gen-alpha time discretization
  con_quot_alphaf = (1.0-alphaf)*con_quot_new + alphaf*con_quot_old;

  /* interfacial surfactant concentration must not be larger than the
   * maximum interfacial concentration */
  if (con_quot_alphaf > con_quot_max)
  {
    con_quot_alphaf = con_quot_max;
  }

  /* check of con_quot_new is done after update of con_quot_alphaf in
   * order to improve estimate of con_quot_alphaf (there is a kink in
   * the transition regime) */
  if (con_quot_new > con_quot_max)
  {
    con_quot_new = con_quot_max;
  }
  (*con_quot_temp_)[LID] = con_quot_new;


  /*-------------------------------------------------------------------*
   |  DETERMINATION OF CURRENT GENERALIZED SURFACE STRESS (n+1-alphaf) |
   *-------------------------------------------------------------------*/

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_quot_alphaf < 1.)
  {
    if (con_quot_old > 0.98)
    {
      k1xC *= (1.-con_quot_old)/0.02;
      k2 *= (1.-con_quot_old)/0.02;
    }

    gamma  = gamma_0-m1*con_quot_alphaf;
    dgamma = m1/dt*M_old/(A_new*A_new*(1/dt+k1xC+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if ((con_quot_max - con_quot_alphaf) > 1.0e-06)
    {
      gamma  = gamma_min_eq-m2*(con_quot_alphaf-1.);
      dgamma = m2*con_quot_new/A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      gamma  = gamma_min;
      dgamma = 0.;
    }
  }

  return;
}

#endif

