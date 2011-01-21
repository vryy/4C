/*!----------------------------------------------------------------------
\file mortar_analytical.H

\brief A set of analytical solutions for convergence analysis of mortar methods

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_ANALYTICAL_H
#define MORTAR_ANALYTICAL_H

#include <math.h>
#include "../linalg/linalg_fixedsizematrix.H"

//------------------------------------------------------------------------
// available analytical solutions
//------------------------------------------------------------------------
// 0 = no output of error norms
// 1 = 3D beam bending example
// 2 = 3D pressurized hollow sphere example
// 3 = 3D pressurized hollow sphere example (thick version)
//------------------------------------------------------------------------
#define MORTARANALYTICALSOL 0

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms
*/
namespace MORTAR
{
  /*!
  \brief Analytical solutions for 3D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  inline void AnalyticalSolutions3D(const LINALG::Matrix<3,1>& pos,
                             LINALG::Matrix<3,1>& uanalyt,
                             LINALG::Matrix<6,1>& epsanalyt,
                             LINALG::Matrix<3,3>& derivanalyt)
  {
    //----------------------------------------------------------------------
    // 3D beam bending solution
    // (see Timoshenko and Goodier, Theory of Elasticity, 1970, p. 284)
    //----------------------------------------------------------------------
    if (MORTARANALYTICALSOL==1)
    {
      // model parameters
      double h= 2.0;
      double p = 100.0;
      double nu = 0.3;
      double E = 1000.0;

      // displacements
      uanalyt(0,0) =  ((2*p)/(E*h)) * (pos(0,0)*pos(1,0));
      uanalyt(1,0) =  (p/(E*h)) * (-pos(0,0)*pos(0,0) - nu*pos(1,0)*pos(1,0) + nu*pos(2,0)*pos(2,0));
      uanalyt(2,0) = -((2*p*nu)/(E*h)) * (pos(1,0)*pos(2,0));

      // strains
      epsanalyt(0,0) =  ((2*p)/(E*h)) * pos(1,0);
      epsanalyt(1,0) = -((2*p*nu)/(E*h)) * pos(1,0);
      epsanalyt(2,0) = -((2*p*nu)/(E*h)) * pos(1,0);
      epsanalyt(3,0) =  0.0;
      epsanalyt(4,0) =  0.0;
      epsanalyt(5,0) =  0.0;

      // displacement derivatives
      derivanalyt(0,0) =  ((2*p)/(E*h)) * pos(1,0);
      derivanalyt(0,1) =  ((2*p)/(E*h)) * pos(0,0);
      derivanalyt(0,2) =  0.0;
      derivanalyt(1,0) = -((2*p)/(E*h)) * pos(0,0);
      derivanalyt(1,1) = -((2*p*nu)/(E*h)) * pos(1,0);
      derivanalyt(1,2) =  ((2*p*nu)/(E*h)) * pos(2,0);
      derivanalyt(2,0) =  0.0;
      derivanalyt(2,1) = -((2*p*nu)/(E*h)) * pos(2,0);
      derivanalyt(2,2) = -((2*p*nu)/(E*h)) * pos(1,0);
    }

    //----------------------------------------------------------------------
    // 3D pressurized sphere solution
    // (see Bower, Applied Mechanics of Solids, 2009, p. 197)
    //----------------------------------------------------------------------
    else if (MORTARANALYTICALSOL==2 || MORTARANALYTICALSOL==3)
    {
      // model parameters
      double a = 0.9;
      double b = 1.1;
      double pi = 1.0;
      double nu = 0.3;
      double E = 1.0;

      // change geometry for thick version
      if (MORTARANALYTICALSOL==3) { a = 0.5; b = 2.0; }

      // radial coordinate r = sqrt(x^2+y^2+z^2)
      double r = sqrt(pos(0,0)*pos(0,0)+pos(1,0)*pos(1,0)+pos(2,0)*pos(2,0));

      // azimuthal angle phi \in [-PI,PI]
      double phi = atan2(pos(1,0),pos(0,0));

      // polar angle theta \in [0,PI]
      double theta = acos(pos(2,0)/r);

      // transformation matrix S
      LINALG::Matrix<3,3> trafo;
      trafo(0,0) =  sin(theta)*cos(phi);
      trafo(0,1) =  cos(theta)*cos(phi);
      trafo(0,2) = -sin(phi);
      trafo(1,0) =  sin(theta)*sin(phi);
      trafo(1,1) =  cos(theta)*sin(phi);
      trafo(1,2) =  cos(phi);
      trafo(2,0) =  cos(theta);
      trafo(2,1) = -sin(theta);
      trafo(2,2) =  0.0;

      // transformation matrix J^(-1)
      LINALG::Matrix<3,3> Jinv;
      Jinv(0,0) =  sin(theta)*cos(phi);
      Jinv(0,1) =  sin(theta)*sin(phi);
      Jinv(0,2) =  cos(theta);
      Jinv(1,0) =  cos(theta)*cos(phi)/r;
      Jinv(1,1) =  cos(theta)*sin(phi)/r;
      Jinv(1,2) = -sin(theta)/r;
      Jinv(2,0) = -sin(phi)/(r*sin(theta));
      Jinv(2,1) =  cos(phi)/(r*sin(theta));
      Jinv(2,2) =  0.0;

      // solution coefficients
      double B = (pi*(1+nu)*(1-2*nu)*a*a*a*b*b*b) / (E*(1+nu)*a*a*a + 2*E*(1-2*nu)*b*b*b);
      double A = -B / (b*b*b);

      // displacements
      LINALG::Matrix<3,1> usphere;
      usphere(0,0) = (A*r) + (B/(r*r));
      usphere(1,0) = 0.0;
      usphere(2,0) = 0.0;
      uanalyt.MultiplyNN(trafo,usphere);

      // strains
      LINALG::Matrix<3,3> epssphere;
      epssphere(0,0) = A - 2*B/(r*r*r);
      epssphere(0,1) = 0.0;
      epssphere(0,2) = 0.0;
      epssphere(1,0) = 0.0;
      epssphere(1,1) = usphere(0,0)/r;
      epssphere(1,2) = 0.0;
      epssphere(2,0) = 0.0;
      epssphere(2,1) = 0.0;
      epssphere(2,2) = usphere(0,0)/r;

      LINALG::Matrix<3,3> temp1, temp2;
      temp1.MultiplyNT(epssphere,trafo);
      temp2.MultiplyNN(trafo,temp1);

      epsanalyt(0,0) = temp2(0,0);
      epsanalyt(1,0) = temp2(1,1);
      epsanalyt(2,0) = temp2(2,2);
      epsanalyt(3,0) = 2*temp2(0,1);
      epsanalyt(4,0) = 2*temp2(1,2);
      epsanalyt(5,0) = 2*temp2(0,2);

      // displacement derivatives
      LINALG::Matrix<3,3> derivsphere;
      derivsphere(0,0) =  sin(theta)*cos(phi)*(A - 2*B/(r*r*r));
      derivsphere(0,1) =  cos(theta)*cos(phi)*(A*r + B/(r*r));
      derivsphere(0,2) = -sin(theta)*sin(phi)*(A*r + B/(r*r));
      derivsphere(1,0) =  sin(theta)*sin(phi)*(A - 2*B/(r*r*r));
      derivsphere(1,1) =  cos(theta)*sin(phi)*(A*r + B/(r*r));
      derivsphere(1,2) =  sin(theta)*cos(phi)*(A*r + B/(r*r));
      derivsphere(2,0) =  cos(theta)*(A - 2*B/(r*r*r));
      derivsphere(2,1) = -sin(theta)*(A*r + B/(r*r));
      derivsphere(2,2) =  0.0;
      derivanalyt.MultiplyNN(derivsphere,Jinv);

      //cout << "\nA: " << A << " B: " << B << " u: " << usphere(0,0) << endl;
      //cout << "x: " << pos(0,0) << " y: " << pos(1,0) << " z: " << pos(2,0) << endl;
      //cout << "r: " << r << " phi: " << phi*180/M_PI << "° theta: " << theta*180/M_PI << "°" << endl;
      //cout << "u: " << usphere(0,0) << " v: " << usphere(1,0) << "w: " << usphere(2,0) << endl;
      //cout << trafo << endl << endl;
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

} // namespace MORTAR

#endif  // #ifndef MORTAR_ANALYTICAL_H
#endif  // #ifdef CCADISCRET
