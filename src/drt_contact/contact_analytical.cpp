/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of analytical solutions for convergence analysis of contact / meshtying methods

\level 2

\maintainer Matthias Mayr
*/
/*-----------------------------------------------------------------------*/

#include <math.h>
#include "contact_analytical.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  Analytical solutions for 2D problems                      popp 06/11|
 *----------------------------------------------------------------------*/
void CONTACT::AnalyticalSolutions2D(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<2, 1>& uanalyt,
    LINALG::Matrix<4, 1>& epsanalyt, LINALG::Matrix<2, 2>& derivanalyt)
{
  // get corresponding input parameter
  const Teuchos::ParameterList& listcmt = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::CONTACT::ErrorNorms entype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ErrorNorms>(listcmt, "ERROR_NORMS");

  //------------------------------------------------------------------------
  // available analytical solutions (enum ErrorNorms)
  //------------------------------------------------------------------------
  // errornorms_none        = no output of error norms
  // errornorms_zero        = 2D/3D zero analytical solution
  // errornorms_bending     = 2D/3D beam bending example
  // errornorms_sphere      = 3D pressurized hollow sphere example
  // errornorms_thicksphere = 3D thick pressurized hollow sphere example
  //------------------------------------------------------------------------

  //----------------------------------------------------------------------
  // no error norm computation
  //----------------------------------------------------------------------
  if (entype == INPAR::CONTACT::errornorms_none)
  {
    dserror("ERROR: Error norm computation switched off. You should not be here.");
  }

  //----------------------------------------------------------------------
  // 2D default (=zero) solution
  //----------------------------------------------------------------------
  else if (entype == INPAR::CONTACT::errornorms_zero)
  {
    // displacements
    uanalyt(0, 0) = 0.0;
    uanalyt(1, 0) = 0.0;

    // strains
    epsanalyt(0, 0) = 0.0;
    epsanalyt(1, 0) = 0.0;
    epsanalyt(2, 0) = 0.0;
    epsanalyt(3, 0) = 0.0;

    // displacement derivatives
    derivanalyt(0, 0) = 0.0;
    derivanalyt(0, 1) = 0.0;
    derivanalyt(1, 0) = 0.0;
    derivanalyt(1, 1) = 0.0;
  }

  //----------------------------------------------------------------------
  // 2D beam bending solution (plane stress)
  // (see Timoshenko and Goodier, Theory of Elasticity, 1970, p. 284)
  //----------------------------------------------------------------------
  else if (entype == INPAR::CONTACT::errornorms_bending)
  {
    // model parameters
    const double h = 2.0;
    const double p = 100.0;
    const double nu = 0.3;
    const double E = 1000.0;

    // displacements
    uanalyt(0, 0) = ((2 * p) / (E * h)) * (pos(0, 0) * pos(1, 0));
    uanalyt(1, 0) = (p / (E * h)) * (-pos(0, 0) * pos(0, 0) - nu * pos(1, 0) * pos(1, 0));

    // strains
    epsanalyt(0, 0) = ((2 * p) / (E * h)) * pos(1, 0);
    epsanalyt(1, 0) = -((2 * p * nu) / (E * h)) * pos(1, 0);
    epsanalyt(2, 0) = 0.0;
    epsanalyt(3, 0) = 0.0;

    // displacement derivatives
    derivanalyt(0, 0) = ((2 * p) / (E * h)) * pos(1, 0);
    derivanalyt(0, 1) = ((2 * p) / (E * h)) * pos(0, 0);
    derivanalyt(1, 0) = -((2 * p) / (E * h)) * pos(0, 0);
    derivanalyt(1, 1) = -((2 * p * nu) / (E * h)) * pos(1, 0);
  }

  else if (entype == INPAR::CONTACT::errornorms_infiniteplate)
  {
    double E = 1.0e5;
    double nue = 0.3;
    const double s_inf = 10.0;  // tensile stress at infinity
    const double a = 1.0;       // hole radius
    const double x = -pos(0, 0);
    const double y = -pos(1, 0);
    LINALG::Matrix<3, 1> stress_analytical;

    // plane strain modification
    E /= 1. - nue * nue;
    nue /= 1. - nue;

    stress_analytical(0) =
        0.5 * s_inf / ((x * x + y * y) * (x * x + y * y) * (x * x + y * y) * (x * x + y * y)) *
        (2.0 * x * x * x * x * x * x * x * x + 8.0 * x * x * x * x * x * x * y * y +
            12.0 * x * x * x * x * y * y * y * y - 5.0 * x * x * x * x * x * x * a * a +
            7.0 * x * x * x * x * a * a * y * y + 13.0 * x * x * a * a * y * y * y * y +
            3.0 * x * x * x * x * a * a * a * a - 18.0 * x * x * a * a * a * a * y * y +
            8.0 * y * y * y * y * y * y * x * x + 2.0 * y * y * y * y * y * y * y * y +
            1.0 * y * y * y * y * y * y * a * a + 3.0 * y * y * y * y * a * a * a * a);
    stress_analytical(1) =
        -1.0 / 2.0 * s_inf * a * a /
        ((x * x + y * y) * (x * x + y * y) * (x * x + y * y) * (x * x + y * y)) *
        (11.0 * x * x * x * x * y * y + 9.0 * x * x * y * y * y * y - 3.0 * y * y * y * y * y * y -
            18.0 * a * a * x * x * y * y + 3.0 * a * a * y * y * y * y -
            1.0 * x * x * x * x * x * x + 3.0 * a * a * x * x * x * x);
    stress_analytical(2) =
        -s_inf * x * a * a /
        ((x * x + y * y) * (x * x + y * y) * (x * x + y * y) * sqrt(x * x + y * y)) *
        (-5.0 * x * x * x * x - 2.0 * x * x * y * y + 3.0 * y * y * y * y + 6.0 * a * a * x * x -
            6.0 * a * a * y * y) *
        sqrt(y * y / (x * x + y * y));

    LINALG::Matrix<3, 3> CmatInv(true);

    // plane stress
    CmatInv(0, 0) = 1.;
    CmatInv(1, 1) = 1.;
    CmatInv(0, 1) = -nue;
    CmatInv(1, 0) = -nue;
    CmatInv(2, 2) = 2. * (1. + nue);
    CmatInv.Scale(1. / E);

    LINALG::Matrix<3, 1> strain_analytical;
    strain_analytical.Multiply(CmatInv, stress_analytical);

    // strains
    epsanalyt(0, 0) = strain_analytical(0);
    epsanalyt(1, 0) = strain_analytical(1);
    epsanalyt(2, 0) = .5 * strain_analytical(2);
    epsanalyt(3, 0) = .5 * strain_analytical(2);

    double kappa = (3. - nue) / (1. + nue);
    double mu = E / (2. * (1. + nue));
    uanalyt.Clear();
    uanalyt(0) = +s_inf * a *
                 (-((kappa + 1.) * x / a) +
                     2. * a *
                         (-(kappa + 1.) * x * pow((x * x + y * y), -1. / 2.) -
                             cos(3. * acos(x * pow((x * x + y * y), -1. / 2.)))) *
                         pow((x * x + y * y), -1. / 2.) +
                     2. * (int)pow(a, 3) * cos(3. * acos(x * pow((x * x + y * y), -1. / 2.))) *
                         pow((x * x + y * y), -3. / 2.)) /
                 mu / 8.;
    uanalyt(1) = s_inf * a / mu *
                 ((kappa - 3) / a * sqrt(x * x + y * y) * sqrt(1. - x * x / (x * x + y * y)) +
                     2. * a *
                         ((1 - kappa) * sqrt(1. - x * x / (x * x + y * y)) +
                             sin(3. * acos(x * pow(x * x + y * y, -1. / 2.)))) *
                         pow(x * x + y * y, -1. / 2.) -
                     2. * pow(a, 3.) * sin(3. * acos(x * pow(x * x + y * y, -1. / 2.))) *
                         pow(x * x + y * y, -3. / 2.)) /
                 8.;
  }

  //----------------------------------------------------------------------
  // Other cases
  //----------------------------------------------------------------------
  else
  {
    dserror("ERROR: Other 2D analytical solutions not yet implemented");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Analytical solutions for 3D problems                      popp 06/11|
 *----------------------------------------------------------------------*/
void CONTACT::AnalyticalSolutions3D(const LINALG::Matrix<3, 1>& pos, LINALG::Matrix<3, 1>& uanalyt,
    LINALG::Matrix<6, 1>& epsanalyt, LINALG::Matrix<3, 3>& derivanalyt)
{
  // get corresponding input parameter
  const Teuchos::ParameterList& listcmt = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::CONTACT::ErrorNorms entype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ErrorNorms>(listcmt, "ERROR_NORMS");

  //------------------------------------------------------------------------
  // available analytical solutions (enum ErrorNorms)
  //------------------------------------------------------------------------
  // errornorms_none        = no output of error norms
  // errornorms_zero        = 2D/3D zero analytical solution
  // errornorms_bending     = 2D/3D beam bending example
  // errornorms_sphere      = 3D pressurized hollow sphere example
  // errornorms_thicksphere = 3D thick pressurized hollow sphere example
  //------------------------------------------------------------------------

  //----------------------------------------------------------------------
  // no error norm computation
  //----------------------------------------------------------------------
  if (entype == INPAR::CONTACT::errornorms_none)
  {
    dserror("ERROR: Error norm computation switched off. You should not be here.");
  }

  //----------------------------------------------------------------------
  // 3D default (=zero) solution
  //----------------------------------------------------------------------
  if (entype == INPAR::CONTACT::errornorms_zero)
  {
    // displacements
    uanalyt(0, 0) = 0.0;
    uanalyt(1, 0) = 0.0;
    uanalyt(2, 0) = 0.0;

    // strains
    epsanalyt(0, 0) = 0.0;
    epsanalyt(1, 0) = 0.0;
    epsanalyt(2, 0) = 0.0;
    epsanalyt(3, 0) = 0.0;
    epsanalyt(4, 0) = 0.0;
    epsanalyt(5, 0) = 0.0;

    // displacement derivatives
    derivanalyt(0, 0) = 0.0;
    derivanalyt(0, 1) = 0.0;
    derivanalyt(0, 2) = 0.0;
    derivanalyt(1, 0) = 0.0;
    derivanalyt(1, 1) = 0.0;
    derivanalyt(1, 2) = 0.0;
    derivanalyt(2, 0) = 0.0;
    derivanalyt(2, 1) = 0.0;
    derivanalyt(2, 2) = 0.0;
  }

  //----------------------------------------------------------------------
  // 3D beam bending solution
  // (see Timoshenko and Goodier, Theory of Elasticity, 1970, p. 284)
  //----------------------------------------------------------------------
  else if (entype == INPAR::CONTACT::errornorms_bending)
  {
    // model parameters
    const double h = 2.0;
    const double p = 100.0;
    const double nu = 0.3;
    const double E = 1000.0;

    // displacements
    uanalyt(0, 0) = ((2 * p) / (E * h)) * (pos(0, 0) * pos(1, 0));
    uanalyt(1, 0) = (p / (E * h)) * (-pos(0, 0) * pos(0, 0) - nu * pos(1, 0) * pos(1, 0) +
                                        nu * pos(2, 0) * pos(2, 0));
    uanalyt(2, 0) = -((2 * p * nu) / (E * h)) * (pos(1, 0) * pos(2, 0));

    // strains
    epsanalyt(0, 0) = ((2 * p) / (E * h)) * pos(1, 0);
    epsanalyt(1, 0) = -((2 * p * nu) / (E * h)) * pos(1, 0);
    epsanalyt(2, 0) = -((2 * p * nu) / (E * h)) * pos(1, 0);
    epsanalyt(3, 0) = 0.0;
    epsanalyt(4, 0) = 0.0;
    epsanalyt(5, 0) = 0.0;

    // displacement derivatives
    derivanalyt(0, 0) = ((2 * p) / (E * h)) * pos(1, 0);
    derivanalyt(0, 1) = ((2 * p) / (E * h)) * pos(0, 0);
    derivanalyt(0, 2) = 0.0;
    derivanalyt(1, 0) = -((2 * p) / (E * h)) * pos(0, 0);
    derivanalyt(1, 1) = -((2 * p * nu) / (E * h)) * pos(1, 0);
    derivanalyt(1, 2) = ((2 * p * nu) / (E * h)) * pos(2, 0);
    derivanalyt(2, 0) = 0.0;
    derivanalyt(2, 1) = -((2 * p * nu) / (E * h)) * pos(2, 0);
    derivanalyt(2, 2) = -((2 * p * nu) / (E * h)) * pos(1, 0);
  }

  //----------------------------------------------------------------------
  // 3D pressurized sphere solution
  // (see Bower, Applied Mechanics of Solids, 2009, p. 197)
  //----------------------------------------------------------------------
  else if (entype == INPAR::CONTACT::errornorms_sphere ||
           entype == INPAR::CONTACT::errornorms_thicksphere)
  {
    // model parameters
    double a = 0.9;
    double b = 1.1;
    const double pi = 1.0;
    const double nu = 0.3;
    const double E = 1.0;

    // change geometry for thick version
    if (entype == INPAR::CONTACT::errornorms_thicksphere)
    {
      a = 0.5;
      b = 2.0;
    }

    // radial coordinate r = sqrt(x^2+y^2+z^2)
    const double r = sqrt(pos(0, 0) * pos(0, 0) + pos(1, 0) * pos(1, 0) + pos(2, 0) * pos(2, 0));

    // azimuthal angle phi \in [-PI,PI]
    const double phi = atan2(pos(1, 0), pos(0, 0));

    // polar angle theta \in [0,PI]
    const double theta = acos(pos(2, 0) / r);

    // transformation matrix S
    LINALG::Matrix<3, 3> trafo;
    trafo(0, 0) = sin(theta) * cos(phi);
    trafo(0, 1) = cos(theta) * cos(phi);
    trafo(0, 2) = -sin(phi);
    trafo(1, 0) = sin(theta) * sin(phi);
    trafo(1, 1) = cos(theta) * sin(phi);
    trafo(1, 2) = cos(phi);
    trafo(2, 0) = cos(theta);
    trafo(2, 1) = -sin(theta);
    trafo(2, 2) = 0.0;

    // transformation matrix J^(-1)
    LINALG::Matrix<3, 3> Jinv;
    Jinv(0, 0) = sin(theta) * cos(phi);
    Jinv(0, 1) = sin(theta) * sin(phi);
    Jinv(0, 2) = cos(theta);
    Jinv(1, 0) = cos(theta) * cos(phi) / r;
    Jinv(1, 1) = cos(theta) * sin(phi) / r;
    Jinv(1, 2) = -sin(theta) / r;
    Jinv(2, 0) = -sin(phi) / (r * sin(theta));
    Jinv(2, 1) = cos(phi) / (r * sin(theta));
    Jinv(2, 2) = 0.0;

    // solution coefficients
    const double B = (pi * (1 + nu) * (1 - 2 * nu) * a * a * a * b * b * b) /
                     (E * (1 + nu) * a * a * a + 2 * E * (1 - 2 * nu) * b * b * b);
    const double A = -B / (b * b * b);

    // displacements
    LINALG::Matrix<3, 1> usphere;
    usphere(0, 0) = (A * r) + (B / (r * r));
    usphere(1, 0) = 0.0;
    usphere(2, 0) = 0.0;
    uanalyt.MultiplyNN(trafo, usphere);

    // strains
    LINALG::Matrix<3, 3> epssphere;
    epssphere(0, 0) = A - 2 * B / (r * r * r);
    epssphere(0, 1) = 0.0;
    epssphere(0, 2) = 0.0;
    epssphere(1, 0) = 0.0;
    epssphere(1, 1) = usphere(0, 0) / r;
    epssphere(1, 2) = 0.0;
    epssphere(2, 0) = 0.0;
    epssphere(2, 1) = 0.0;
    epssphere(2, 2) = usphere(0, 0) / r;

    LINALG::Matrix<3, 3> temp1, temp2;
    temp1.MultiplyNT(epssphere, trafo);
    temp2.MultiplyNN(trafo, temp1);

    epsanalyt(0, 0) = temp2(0, 0);
    epsanalyt(1, 0) = temp2(1, 1);
    epsanalyt(2, 0) = temp2(2, 2);
    epsanalyt(3, 0) = 2 * temp2(0, 1);
    epsanalyt(4, 0) = 2 * temp2(1, 2);
    epsanalyt(5, 0) = 2 * temp2(0, 2);

    // displacement derivatives
    LINALG::Matrix<3, 3> derivsphere;
    derivsphere(0, 0) = sin(theta) * cos(phi) * (A - 2 * B / (r * r * r));
    derivsphere(0, 1) = cos(theta) * cos(phi) * (A * r + B / (r * r));
    derivsphere(0, 2) = -sin(theta) * sin(phi) * (A * r + B / (r * r));
    derivsphere(1, 0) = sin(theta) * sin(phi) * (A - 2 * B / (r * r * r));
    derivsphere(1, 1) = cos(theta) * sin(phi) * (A * r + B / (r * r));
    derivsphere(1, 2) = sin(theta) * cos(phi) * (A * r + B / (r * r));
    derivsphere(2, 0) = cos(theta) * (A - 2 * B / (r * r * r));
    derivsphere(2, 1) = -sin(theta) * (A * r + B / (r * r));
    derivsphere(2, 2) = 0.0;
    derivanalyt.MultiplyNN(derivsphere, Jinv);

    // std::cout << "\nA: " << A << " B: " << B << " u: " << usphere(0,0) << std::endl;
    // std::cout << "x: " << pos(0,0) << " y: " << pos(1,0) << " z: " << pos(2,0) << std::endl;
    // std::cout << "r: " << r << " phi: " << phi*180/M_PI << "° theta: " << theta*180/M_PI << "°"
    // << std::endl; std::cout << "u: " << usphere(0,0) << " v: " << usphere(1,0) << "w: " <<
    // usphere(2,0) << std::endl; std::cout << trafo << std::endl << std::endl;
  }

  else if (entype == INPAR::CONTACT::errornorms_infiniteplate)
  {
    const double E = 100.0;
    const double nue = 0.0;
    const double s_inf = 1.0;  // tensile stress at infinity
    const double a = 2.0;      // hole radius
    const double x = -pos(0, 0);
    const double y = -pos(1, 0);
    LINALG::Matrix<3, 1> stress_analytical;


    stress_analytical(0) =
        0.5 * s_inf / ((x * x + y * y) * (x * x + y * y) * (x * x + y * y) * (x * x + y * y)) *
        (2.0 * x * x * x * x * x * x * x * x + 8.0 * x * x * x * x * x * x * y * y +
            12.0 * x * x * x * x * y * y * y * y - 5.0 * x * x * x * x * x * x * a * a +
            7.0 * x * x * x * x * a * a * y * y + 13.0 * x * x * a * a * y * y * y * y +
            3.0 * x * x * x * x * a * a * a * a - 18.0 * x * x * a * a * a * a * y * y +
            8.0 * y * y * y * y * y * y * x * x + 2.0 * y * y * y * y * y * y * y * y +
            1.0 * y * y * y * y * y * y * a * a + 3.0 * y * y * y * y * a * a * a * a);
    stress_analytical(1) =
        -1.0 / 2.0 * s_inf * a * a /
        ((x * x + y * y) * (x * x + y * y) * (x * x + y * y) * (x * x + y * y)) *
        (11.0 * x * x * x * x * y * y + 9.0 * x * x * y * y * y * y - 3.0 * y * y * y * y * y * y -
            18.0 * a * a * x * x * y * y + 3.0 * a * a * y * y * y * y -
            1.0 * x * x * x * x * x * x + 3.0 * a * a * x * x * x * x);
    stress_analytical(2) =
        -s_inf * x * a * a /
        ((x * x + y * y) * (x * x + y * y) * (x * x + y * y) * sqrt(x * x + y * y)) *
        (-5.0 * x * x * x * x - 2.0 * x * x * y * y + 3.0 * y * y * y * y + 6.0 * a * a * x * x -
            6.0 * a * a * y * y) *
        sqrt(y * y / (x * x + y * y));

    LINALG::Matrix<3, 3> CmatInv(true);

    // plane stress
    CmatInv(0, 0) = 1.;
    CmatInv(1, 1) = 1.;
    CmatInv(0, 1) = -nue;
    CmatInv(1, 0) = -nue;
    CmatInv(2, 2) = 2. * (1. + nue);
    CmatInv.Scale(1. / E);

    LINALG::Matrix<3, 1> strain_analytical;
    strain_analytical.Multiply(CmatInv, stress_analytical);

    // strains
    epsanalyt(0, 0) = strain_analytical(0);
    epsanalyt(1, 0) = strain_analytical(1);
    epsanalyt(2, 0) = 0.;
    epsanalyt(4, 0) = 0.;
    epsanalyt(5, 0) = 0.;
    epsanalyt(3, 0) = strain_analytical(2);
  }


  //----------------------------------------------------------------------
  // Other cases
  //----------------------------------------------------------------------
  else
  {
    dserror("ERROR: Other 3D analytical solutions not yet implemented");
  }

  return;
}
