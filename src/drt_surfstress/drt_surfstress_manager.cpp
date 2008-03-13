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
    // choose what to assemble
    p.set("assemble matrix 1",true);
    p.set("assemble matrix 2",false);
    p.set("assemble vector 1",true);
    p.set("assemble vector 2",false);
    p.set("assemble vector 3",false);

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
| (public)                                                   lw 03/08|
|                                                                    |
| Calculate interfacial area and its derivatives with respect        |
| to the nodal displacements                                         |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::SurfaceCalc(
  const vector<Epetra_SerialDenseMatrix>& dx,       // (i) current Jacobian
  const vector<Epetra_SerialDenseMatrix>& der,      // (i) derivatives of shape functions
  const Epetra_SerialDenseMatrix& gpcoord,          // (i) coordinates of Gauss points
  const Epetra_SerialDenseVector& gpweight,         // (i) Gaussian weights
  const int ngp,                                    // (i) number of Gauss points
  const int nnode,                                  // (i) number of surface nodes
  double& A,                                        // (o) surface area
  Epetra_SerialDenseVector& Adiff,                  // (o) first partial derivatives of area
  Epetra_SerialDenseMatrix& Adiff2)                 // (o) second partial derivatives of area

{
  int ndof = nnode*3;

  blitz::Array<double,1> det(3);
  blitz::Array<double,2> ddet(3,ndof);
  blitz::Array<double,3> ddet2(3,ndof,ndof);
  blitz::Array<double,1> jacobi_deriv(ndof);
  double Jac;

  /*-----------------------------------------------------initialization*/
  A = 0.;

  for (int i=0;i<ndof;++i)
  {
    Adiff[i] = 0.;
    for (int j=0;j<ndof;++j)
    {
      Adiff2(i,j) = 0.;
    }
  }

  ddet2 = 0.0;

  for (int gpid = 0; gpid < ngp; ++gpid)      // loop over integration points
  {
    Epetra_SerialDenseMatrix deriv = der[gpid];
    Epetra_SerialDenseMatrix dxyzdrs = dx[gpid];

    det(0) = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    det(1) = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    det(2) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

    Jac = sqrt( det(0)*det(0) + det(1)*det(1) + det(2)*det(2) );

    A += Jac*gpweight[gpid];

    /*--------------- derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int i=0;i<nnode;++i)
    {
      ddet(0,3*i)   = 0.;
      ddet(0,3*i+1) = deriv(i,0)*dxyzdrs(1,2)-deriv(i,1)*dxyzdrs(0,2);
      ddet(0,3*i+2) = deriv(i,1)*dxyzdrs(0,1)-deriv(i,0)*dxyzdrs(1,1);

      ddet(1,3*i)   = deriv(i,1)*dxyzdrs(0,2)-deriv(i,0)*dxyzdrs(1,2);
      ddet(1,3*i+1) = 0.;
      ddet(1,3*i+2) = deriv(i,0)*dxyzdrs(1,0)-deriv(i,1)*dxyzdrs(0,0);

      ddet(2,3*i)   = deriv(i,0)*dxyzdrs(1,1)-deriv(i,1)*dxyzdrs(0,1);
      ddet(2,3*i+1) = deriv(i,1)*dxyzdrs(0,0)-deriv(i,0)*dxyzdrs(1,0);
      ddet(2,3*i+2) = 0.;

      jacobi_deriv(i*3)   = 1/Jac*(det(2)*ddet(2,3*i  )+det(1)*ddet(1,3*i  ));
      jacobi_deriv(i*3+1) = 1/Jac*(det(2)*ddet(2,3*i+1)+det(0)*ddet(0,3*i+1));
      jacobi_deriv(i*3+2) = 1/Jac*(det(0)*ddet(0,3*i+2)+det(1)*ddet(1,3*i+2));
    }

    /*--- calculation of first derivative of current interfacial area
     *----------------------------- with respect to the displacements */
    for (int i=0;i<ndof;++i)
    {
      Adiff[i] += jacobi_deriv(i)*gpweight[gpid];
    }


    /*-------- second derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int n=0;n<nnode;++n)
    {
      for (int o=0;o<nnode;++o)
      {
        ddet2(0,n*3+1,o*3+2) = deriv(n,0)*deriv(o,1)-deriv(n,1)*deriv(o,0);
        ddet2(0,n*3+2,o*3+1) = - ddet2(0,n*3+1,o*3+2);

        ddet2(1,n*3  ,o*3+2) = deriv(n,1)*deriv(o,0)-deriv(n,0)*deriv(o,1);
        ddet2(1,n*3+2,o*3  ) = - ddet2(1,n*3,o*3+2);

        ddet2(2,n*3  ,o*3+1) = deriv(n,0)*deriv(o,1)-deriv(n,1)*deriv(o,0);
        ddet2(2,n*3+1,o*3  ) = - ddet2(2,n*3,o*3+1);
      }
    }

    /*- calculation of second derivatives of current interfacial areas
     *----------------------------- with respect to the displacements */
    for (int i=0;i<ndof;++i)
    {
      int var1, var2;

      if (i%3==0)           // displacement in x-direction
      {
        var1 = 1;
        var2 = 2;
      }
      else if ((i-1)%3==0)  // displacement in y-direction
      {
        var1 = 0;
        var2 = 2;
      }
      else if ((i-2)%3==0)  // displacement in z-direction
      {
        var1 = 0;
        var2 = 1;
      }
      else
      {
          dserror("calculation of second derivatives of interfacial area failed");
          exit(1);
      }

      for (int j=0;j<ndof;++j)
      {
        Adiff2(i,j) += (-1/Jac*jacobi_deriv(j)*jacobi_deriv(i)+1/Jac*
                       (ddet(var1,i)*ddet(var1,j)+det(var1)*ddet2(var1,i,j)+
                        ddet(var2,i)*ddet(var2,j)+det(var2)*ddet2(var2,i,j)))*gpweight[gpid];
      }
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

