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

#include "drt_surfstress_manager.H"
#include "../drt_lib/linalg_utils.H"

/*-------------------------------------------------------------------*
 |  ctor (public)                                            lw 12/07|
 *-------------------------------------------------------------------*/
DRT::SurfStressManager::SurfStressManager(DRT::Discretization& discret,
                                          const int numsurf):
discret_(discret)
{
  Epetra_Map surfmap(numsurf, 0, discret_.Comm());
  time_ = 0.0;
  A_old_temp_    = rcp(new Epetra_Vector(surfmap,true));
  A_old_         = rcp(new Epetra_Vector(surfmap,true));
  con_quot_temp_ = rcp(new Epetra_Vector(surfmap,true));
  con_quot_      = rcp(new Epetra_Vector(surfmap,true));
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
                                                RefCountPtr<Epetra_CrsMatrix> stiff)
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

    discret_.EvaluateCondition(p,stiff,fint,null,"SurfaceStress");

    return;
}

/*-------------------------------------------------------------------*
| (public)						     lw 12/07|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level                                                   |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::StiffnessAndInternalForces(const Epetra_SerialDenseMatrix& xs,
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
  //cout << "ID: " << ID << "\n";

  Epetra_SerialDenseVector Adiff(12);
  Epetra_SerialDenseMatrix Adiff2(12,12);
  double gamma, dgamma;

  /*---------------------------------------------- update if necessary */
  if (time != time_)
  {
    (*A_old_)[ID] = (*A_old_temp_)[ID];
    (*con_quot_)[ID] = (*con_quot_temp_)[ID];
    time_ = time;
  }

  /*--------------------------------------------------- initialization */
  (*A_old_temp_)[ID] = 0.0;
  (*con_quot_temp_)[ID] = 0.0;

  soh8_surface_calc(xs, (*A_old_temp_)[ID], Adiff, Adiff2);

  if (surface_flag==0)                                   // SURFACTANT
  {
    /*-----------calculation of current surface stress and its partial
     *-----------------derivative with respect to the interfacial area */
    if (time<100.*dt)         /* gradual application of surface stress */
    {
      gamma = gamma_0-m1*con_quot_eq;
      (*con_quot_temp_)[ID] = con_quot_eq;
      dgamma = 0.;
    }
    else
    {
      (*con_quot_temp_)[ID] = (*con_quot_)[ID];
      gamma_calc(gamma, dgamma, (*con_quot_temp_)[ID], (*A_old_)[ID], (*A_old_temp_)[ID], dt, k1xC, k2, m1,
               m2, gamma_0, gamma_min, gamma_min_eq, con_quot_max);
    }
  }

  else if (surface_flag==1)                             // SURFACE TENSION
  {
    gamma = const_gamma;
    dgamma = 0.;
  }

  else
    dserror("Surface flag not implemented");

  for (int i=0;i<12;i++)
  {
    for (int j=0;j<12;j++)
    {
      K_surf(i,j) = dgamma*Adiff[i]*Adiff[j]+gamma*Adiff2(i,j);
    }
  }
  if (time<100.*dt)     /* gradual application of surface stress */
  {
    K_surf.Scale(-2*pow(time/(99.*dt), 3)+3*pow(time/(99.*dt), 2));
  }

//   cout << "Ksurf:\n" << K_surf << "\n";

  /*------calculation of current internal force due to surface energy*/
  for (int i=0;i<12;i++)
  {
    fint[i] = gamma*Adiff[i];
  }
  if (time<100.*dt)       /* gradual application of surface stress */
  {
    fint.Scale(-2*pow(time/(99.*dt), 3)+3*pow(time/(99.*dt), 2));
  }

//   cout << "fsurf:\n" << fint << "\n";

  return;
}

/*-------------------------------------------------------------------*
| (public)                                                   lw 12/07|
|                                                                    |
| Calculate interfacial area and its derivatives with respect        |
| to the nodal displacements                                         |
*--------------------------------------------------------------------*/

void DRT::SurfStressManager::soh8_surface_calc(
      const Epetra_SerialDenseMatrix& xs,    // (i) element coords
      double& A,                              // (o) surface area
      Epetra_SerialDenseVector& Adiff,       // (o) first partial derivatives of area
      Epetra_SerialDenseMatrix& Adiff2)      // (o) second partial derivatives of area

{
  double det[3], ddet[3][12], ddet2[3][12][12], jacobi_deriv[12];
  double Jac;

  /*-----------------------------------------------------initialization*/
  A = 0.;

  for (int i=0;i<12;i++)
  {
    Adiff[i] = 0.;
    for (int j=0;j<12;j++)
    {
      Adiff2(i,j) = 0.;
    }
  }

  for (int n=0;n<3;n++)
    {
      for (int o=0;o<12;o++)
      {
        for (int p=0;p<12;p++)
          ddet2[n][o][p]=0.;
      }
    }

  /*
  ** Here we integrate a 4-node surface with 2x2 Gauss Points
  */
  const int ngp = 4;

  // gauss parameters
  const double gpweight = 1.0;
  const double gploc    = 1.0/sqrt(3.0);
  Epetra_SerialDenseMatrix gpcoord (ngp,2);
  gpcoord(0,0) = - gploc;
  gpcoord(0,1) = - gploc;
  gpcoord(1,0) =   gploc;
  gpcoord(1,1) = - gploc;
  gpcoord(2,0) = - gploc;
  gpcoord(2,1) =   gploc;
  gpcoord(3,0) =   gploc;
  gpcoord(3,1) =   gploc;

  for (int gpid = 0; gpid < 4; ++gpid)      // loop over integration points
  {
    double r = gpcoord(gpid,0);
    double s = gpcoord(gpid,1);

    // get shape functions and derivatives of element surface
    vector<double> funct(ngp);                // 4 shape function values

    // shape functions for 4 nodes
    funct[0] = 0.25 * (1.0-r) * (1.0-s);
    funct[1] = 0.25 * (1.0+r) * (1.0-s);
    funct[2] = 0.25 * (1.0+r) * (1.0+s);
    funct[3] = 0.25 * (1.0-r) * (1.0+s);
    // derivatives of 4 shape functions wrt 2 directions
    Epetra_SerialDenseMatrix deriv(4,2);
    deriv(0,0) = -0.25 * (1.0-s);
    deriv(0,1) = -0.25 * (1.0-r);
    deriv(1,0) =  0.25 * (1.0-s);
    deriv(1,1) = -0.25 * (1.0+r);
    deriv(2,0) =  0.25 * (1.0+s);
    deriv(2,1) =  0.25 * (1.0+r);
    deriv(3,0) = -0.25 * (1.0+s);
    deriv(3,1) =  0.25 * (1.0-r);

    // compute dXYZ / drs
    Epetra_SerialDenseMatrix dxyzdrs(2,3);
    dxyzdrs.Multiply('T','N',1.0,deriv,xs,1.0);

    det[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    det[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    det[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

    Jac = sqrt( det[0]*det[0] + det[1]*det[1] + det[2]*det[2] );

    A += Jac*gpweight;

    /*--------------- derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int i=0;i<4;i++)
    {
      ddet[0][3*i]   = 0.;
      ddet[0][3*i+1] = deriv(i,0)*dxyzdrs(1,2)-deriv(i,1)*dxyzdrs(0,2);
      ddet[0][3*i+2] = deriv(i,1)*dxyzdrs(0,1)-deriv(i,0)*dxyzdrs(1,1);

      ddet[1][3*i]   = deriv(i,1)*dxyzdrs(0,2)-deriv(i,0)*dxyzdrs(1,2);
      ddet[1][3*i+1] = 0.;
      ddet[1][3*i+2] = deriv(i,0)*dxyzdrs(1,0)-deriv(i,1)*dxyzdrs(0,0);

      ddet[2][3*i]   = deriv(i,0)*dxyzdrs(1,1)-deriv(i,1)*dxyzdrs(0,1);
      ddet[2][3*i+1] = deriv(i,1)*dxyzdrs(0,0)-deriv(i,0)*dxyzdrs(1,0);
      ddet[2][3*i+2] = 0.;

      jacobi_deriv[i*3]   = 1/Jac*(det[2]*ddet[2][3*i]+det[1]*ddet[1][3*i]);
      jacobi_deriv[i*3+1] = 1/Jac*(det[2]*ddet[2][3*i+1]+det[0]*ddet[0][3*i+1]);
      jacobi_deriv[i*3+2] = 1/Jac*(det[0]*ddet[0][3*i+2]+det[1]*ddet[1][3*i+2]);
    }
//     cout << "ddet:\n";
//     for (int i=0;i<3;++i)
//     {
//       for (int j=0;j<12;++j)
//         cout << ddet[i][j] << "    ";
//       cout << "\n";
//     }
//     cout << "jacobi_deriv:\n";
//     for (int i=0;i<12;++i)
//     {
//       cout << jacobi_deriv[i] << "\n";
//     }

    /*--- calculation of first derivative of current interfacial area
     *----------------------------- with respect to the displacements */
    for (int i=0;i<12;i++)
    {
      Adiff[i] += jacobi_deriv[i]*gpweight;
    }


    /*-------- second derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int n=0;n<4;n++)
    {
      for (int o=0;o<4;o++)
      {
        ddet2[0][n*3+1][o*3+2] = deriv(0,n)*deriv(1,o)-deriv(1,n)*deriv(0,o);
        ddet2[0][n*3+2][o*3+1] = -ddet2[0][n*3+1][o*3+2];

        ddet2[1][n*3][o*3+2]   = deriv(1,n)*deriv(0,o)-deriv(0,n)*deriv(1,o);
        ddet2[1][n*3+2][o*3]   = -ddet2[1][n*3][o*3+2];

        ddet2[2][n*3][o*3+1]   = deriv(0,n)*deriv(1,o)-deriv(1,n)*deriv(0,o);
        ddet2[2][n*3+1][o*3]   = -ddet2[2][n*3][o*3+1];
      }
    }

    /*- calculation of second derivatives of current interfacial areas
     *----------------------------- with respect to the displacements */
    int var1, var2;

    for (int i=0;i<12;i++)
    {
      if (i==0 || i==3 || i==6 || i==9)
      {
        var1 = 1;
        var2 = 2;
      }
      if (i==1 || i==4 || i==7 || i==10)
      {
        var1 = 0;
        var2 = 2;
      }
      if (i==2 || i==5 || i==8 || i==11)
      {
        var1 = 0;
        var2 = 1;
      }

      for (int j=0;j<12;j++)
      {
        Adiff2(i,j) += (-1/Jac*jacobi_deriv[j]*jacobi_deriv[i]+1/Jac*
                       (ddet[var1][i]*ddet[var1][j]+det[var1]*ddet2[var1][i][j]+
                        ddet[var2][i]*ddet[var2][j]+det[var2]*ddet2[var2][i][j]))*gpweight;
      }
    }
  }
//   cout << "A: " << A << "\n";
//   cout << "Adiff:\n" << Adiff << "\nAdiff2:\n" << Adiff2 << "\n";
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

