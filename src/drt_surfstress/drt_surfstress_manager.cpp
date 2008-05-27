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

void DRT::SurfStressManager::GetHistory(RCP<Epetra_Vector> A_old_temp_row,
                                        RCP<Epetra_Vector> con_quot_temp_row)
{
  // Note that the temporal vectors need to be written since in the
  // final ones we still have the data of the former step. The column
  // map based vector used for calculations is exported to a row map
  // based one needed for writing.

  LINALG::Export(*A_old_temp_, *A_old_temp_row);
  LINALG::Export(*con_quot_temp_, *con_quot_temp_row);
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
                                                        const double con_quot_eq)
{
  double gamma, dgamma;
  int LID = A_old_->Map().LID(ID);
  double t_end = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).end();

  if (surface_flag==0)                                   // SURFACTANT
  {
    /*------------------------------------------------- initialization */
    (*con_quot_temp_)[LID] = 0.0;
    (*A_old_temp_)[LID] = A;

    /*-----------calculation of current surface stress and its partial
     *-----------------derivative with respect to the interfacial area */

    if (time <= t_end)         /* gradual application of surface stress */
    {
      gamma = gamma_0-m1*con_quot_eq;
      (*con_quot_temp_)[LID] = con_quot_eq;
      dgamma = 0.;
    }
    else
    {
      (*con_quot_temp_)[LID] = (*con_quot_)[LID];
      gamma_calc(gamma, dgamma, (*con_quot_temp_)[LID], (*A_old_)[LID],
                 (*A_old_temp_)[LID], dt, k1xC, k2, m1, m2, gamma_0,
                 gamma_min, gamma_min_eq, con_quot_max);
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
void DRT::SurfStressManager::gamma_calc(
                double& gamma,                // (o) surface stress
                double& dgamma,               // (o) derivative of surface stress
                double& con_quot,             // (i/o) surfactant concentration
                const double A_old,          // (i) surface area of last time step
                const double A_new,          // (i) current surface area
                const double dt,             // (i) timestep size
                const double k1xC,           // (i) adsorption coefficent k1
                                              //     times bulk concentration C
                const double k2,             // (i) desorption coefficient
                const double m1,             // (i) 1st isothermal slope
                const double m2,             // (i) 2nd isothermal slope
                const double gamma_0,        // (i) surface tension of water
                const double gamma_min,      // (i) mimimum surface stress
                const double gamma_min_eq,   // (i) mimimum equilibrium surface stress
                const double con_quot_max)   // (i) max. surfactant concentration
{
  /* use of non-dimensionalized interfacial surfactant concentration:
   * con_quot=Gamma/Gamma_min_eq
   * mass of surfactant molecules M_old is also divided
   * by minimum equilibrium interfacial surfactant concentration!*/

  double M_old = con_quot*A_old;     // mass of surfactant molecules of
                                     // last time step

  /*-------------------------------------------------------------------*
   |     DETERMINATION OF CURRENT NON-DIMENSIONALIZED INTERFACIAL      |
   |             SURFACTANT CONCENTRATION CON_QUOT_NEW                 |
   *-------------------------------------------------------------------*/

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_quot < 1.)
  {
    /* assumed continuous drop of adsorption/desorption coefficients
     * when approaching insoluble regime (cf. Morris et al. 2001,
     * p108) is only needed for high bulk concentrations C! To
     * simplify matters, drop is left out */

    con_quot = (1.0/dt*M_old+k1xC*A_new)/(A_new*(1.0/dt+k1xC+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_quot < con_quot_max)
    {
      con_quot = M_old/A_new;
    }
    /*----------------------------Regime 3: "Squeeze out"/Film collapse*/
    else
    {
      if (A_new < A_old)
      {
        con_quot = con_quot_max;
      }
      else
      {
        con_quot = M_old/A_new;
      }
    }
  }

  /* interfacial surfactant concentration must not be larger than the
   * maximum interfacial concentration */

  if (con_quot > con_quot_max)
  {
    con_quot = con_quot_max;
  }

  /*-------------------------------------------------------------------*
   |        DETERMINATION OF CURRENT GENERALIZED SURFACE ENERGY        |
   *-------------------------------------------------------------------*/

  /*----------------Regime 1: Langmuir kinetics (adsorption/desorption)*/

  if (con_quot < 1.)
  {
    gamma  = gamma_0-m1*con_quot;
    dgamma = m1/dt*M_old/(A_new*A_new*(1/dt+k1xC+k2));
  }
  else
  {
    /*------------------------------------Regime 2: Insoluble monolayer*/
    if (con_quot < con_quot_max)
    {
      gamma  = gamma_min_eq-m2*(con_quot-1.);
      dgamma = m2*con_quot/A_new;
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

