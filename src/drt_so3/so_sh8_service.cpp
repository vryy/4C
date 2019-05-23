/*----------------------------------------------------------------------*/
/**
\brief Service functionality for the solid shell8 elements

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
\level 1
*/
/*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "so_sh8.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_io/io_gmsh.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/material.H"


/*----------------------------------------------------------------------*
 |  find shell-thickness direction via Jacobian                maf 07/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::ThicknessDirection DRT::ELEMENTS::So_sh8::sosh8_findthickdir()
{
  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(false);  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }
  // vector of df(origin), ie parametric derivatives of shape functions
  // evaluated at the origin (r,s,t)=(0,0,0)
  const double df0_vector[] = {-0.125, -0.125, -0.125, +0.125, -0.125, -0.125, +0.125, +0.125,
      -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125, -0.125, +0.125, +0.125,
      +0.125, +0.125, -0.125, +0.125, +0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> df0(df0_vector);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  // (J0_i^A) = (X^A_{,i})^T
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac0;
  jac0.MultiplyNN(df0, xrefe);
  // compute inverse of Jacobian at element origin
  // (Jinv0_A^i) = (X^A_{,i})^{-T}
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> iJ0(jac0);
  iJ0.Invert();

  // separate "stretch"-part of J-mapping between parameter and global space
  // (G0^ji) = (Jinv0^j_B) (krondelta^BA) (Jinv0_A^i)
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac0stretch;
  jac0stretch.MultiplyTN(iJ0, iJ0);
  const double r_stretch = sqrt(jac0stretch(0, 0));
  const double s_stretch = sqrt(jac0stretch(1, 1));
  const double t_stretch = sqrt(jac0stretch(2, 2));

  // minimal stretch equivalents with "thinnest" direction
  // const double max_stretch = max(r_stretch, max(s_stretch, t_stretch));
  double max_stretch = -1.0;

  ThicknessDirection thickdir = none;  // of actual element
  int thick_index = -1;

  if (r_stretch >= s_stretch and r_stretch >= t_stretch)
  {
    max_stretch = r_stretch;
    if ((max_stretch / s_stretch <= 1.5) || (max_stretch / t_stretch <= 1.5))
    {
      // std::cout << "ID: " << this->Id() << ", has aspect ratio of: ";
      // std::cout << max_stretch / s_stretch << " , " << max_stretch / t_stretch << std::endl;
      // dserror("Solid-Shell element geometry has not a shell aspect ratio");
      return undefined;
    }
    thickdir = autor;
    thick_index = 0;
  }
  else if (s_stretch > r_stretch and s_stretch >= t_stretch)
  {
    max_stretch = s_stretch;
    if ((max_stretch / r_stretch <= 1.5) || (max_stretch / t_stretch <= 1.5))
    {
      // std::cout << "ID: " << this->Id() << ", has aspect ratio of: ";
      // std::cout << max_stretch / r_stretch << " , " << max_stretch / t_stretch << std::endl;
      // dserror("Solid-Shell element geometry has not a shell aspect ratio");
      return undefined;
    }
    thickdir = autos;
    thick_index = 1;
  }
  else if (t_stretch > r_stretch and t_stretch > s_stretch)
  {
    max_stretch = t_stretch;
    if ((max_stretch / r_stretch <= 1.5) || (max_stretch / s_stretch <= 1.5))
    {
      // std::cout << "ID: " << this->Id() << ", has aspect ratio of: ";
      // std::cout << max_stretch / r_stretch << " , " << max_stretch / s_stretch << std::endl;
      // dserror("Solid-Shell element geometry has not a shell aspect ratio");
      return undefined;
    }
    thickdir = autot;
    thick_index = 2;
  }

  if (thick_index == -1)
    dserror("Trouble with thick_index=%d %g,%g,%g,%g", thick_index, r_stretch, s_stretch, t_stretch,
        max_stretch);

  // thickness-vector in parameter-space, has 1.0 in thickness-coord
  LINALG::Matrix<NUMDIM_SOH8, 1> loc_thickvec(true);
  loc_thickvec(thick_index) = 1.0;
  // thickness-vector in global coord is J times local thickness-vector
  // (X^A) = (J0_i^A)^T . (xi_i)
  LINALG::Matrix<NUMDIM_SOH8, 1> glo_thickvec;
  glo_thickvec.MultiplyTN(jac0, loc_thickvec);
  // return doubles of thickness-vector
  thickvec_.resize(3);
  thickvec_[0] = glo_thickvec(0);
  thickvec_[1] = glo_thickvec(1);
  thickvec_[2] = glo_thickvec(2);

  return thickdir;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::So_sh8::sosh8_calcaspectratio()
{
  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(false);  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }
  // vector of df(origin), ie parametric derivatives of shape functions
  // evaluated at the origin (r,s,t)=(0,0,0)
  const double df0_vector[] = {-0.125, -0.125, -0.125, +0.125, -0.125, -0.125, +0.125, +0.125,
      -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125, -0.125, +0.125, +0.125,
      +0.125, +0.125, -0.125, +0.125, +0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> df0(df0_vector);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  // (J0_i^A) = (X^A_{,i})^T
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac0;
  jac0.MultiplyNN(df0, xrefe);
  // compute inverse of Jacobian at element origin
  // (Jinv0_A^i) = (X^A_{,i})^{-T}
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> iJ0(jac0);
  iJ0.Invert();

  // separate "stretch"-part of J-mapping between parameter and global space
  // (G0^ji) = (Jinv0^j_B) (krondelta^BA) (Jinv0_A^i)
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac0stretch;
  jac0stretch.MultiplyTN(iJ0, iJ0);
  const double r_stretch = sqrt(jac0stretch(0, 0));
  const double s_stretch = sqrt(jac0stretch(1, 1));
  const double t_stretch = sqrt(jac0stretch(2, 2));

  // return an averaged aspect ratio
  if (r_stretch >= s_stretch and r_stretch >= t_stretch)
  {
    //    return std::min(r_stretch/s_stretch,r_stretch/t_stretch);
    return 0.5 * (r_stretch / s_stretch + r_stretch / t_stretch);
  }
  else if (s_stretch > r_stretch and s_stretch >= t_stretch)
  {
    //    return std::min(s_stretch/r_stretch,s_stretch/t_stretch);
    return 0.5 * (s_stretch / r_stretch + s_stretch / t_stretch);
  }
  else if (t_stretch > r_stretch and t_stretch > s_stretch)
  {
    //    return std::min(t_stretch/s_stretch,t_stretch/r_stretch);
    return 0.5 * (t_stretch / s_stretch + t_stretch / r_stretch);
  }

  return 0.0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::ThicknessDirection DRT::ELEMENTS::So_sh8::sosh8_enfthickdir(
    LINALG::Matrix<NUMDIM_SOH8, 1>& thickdirglo)
{
  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(false);  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }
  // vector of df(origin), ie parametric derivatives of shape functions
  // evaluated at the origin (r,s,t)=(0,0,0)
  const double df0_vector[NUMDOF_SOH8 * NUMNOD_SOH8] = {-0.125, -0.125, -0.125, +0.125, -0.125,
      -0.125, +0.125, +0.125, -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125,
      -0.125, +0.125, +0.125, +0.125, +0.125, -0.125, +0.125, +0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> df0(df0_vector);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  // (J0_i^A) = (X^A_{,i})^T
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac0(false);
  jac0.MultiplyNN(df0, xrefe);

  // compute inverse of Jacobian at element origin
  // (Jinv0_A^i) = (X^A_{,i})^{-T}
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> iJ0(jac0);
  iJ0.Invert();

  // make enforced global thickness direction a unit vector
  const double thickdirglolength = thickdirglo.Norm2();
  thickdirglo.Scale(1.0 / thickdirglolength);

  // pull thickness direction from global to contra-variant local
  // (dxi^i) = (Jinv0_A^i)^T . (dX^A)
  LINALG::Matrix<NUMDIM_SOH8, 1> thickdirlocsharp(false);
  thickdirlocsharp.MultiplyTN(iJ0, thickdirglo);

#if 0
  // metric tensor
  // (G0_ji) = (J0_j^B) (delta_BA) (J0^A_i)
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> metrflat(false);
  metrflat.MultiplyNT(jac0,jac0);

  // co-variant local enforced thickness direction
  // (dxi_j) = (G0_ji) (dxi^i)
  LINALG::Matrix<NUMDIM_SOH8,1> thickdirlocflat(false);
  thickdirlocflat.MultiplyNN(metrflat,thickdirlocsharp);

  // thickdirloclength
  const double thickdirloclength = thickdirlocsharp.Dot(thickdirlocflat);

  // check if transformation was successful
  if (fabs(thickdirglolength-thickdirloclength)>EPS6)
    dserror("Transformation erroneous: Vector length is not in-variant: %g!=%g",
            thickdirglolength, thickdirloclength);
#endif

  // identify parametric co-ordinate closest to enforced thickness direction
  int thick_index = -1;
  double thickdirlocmax = 0.0;
  for (int i = 0; i < NUMDIM_SOH8; ++i)
  {
    if (fabs(thickdirlocsharp(i)) > thickdirlocmax)
    {
      thickdirlocmax = fabs(thickdirlocsharp(i));
      thick_index = i;
    }
  }
  const double tol = 0.9;  // should be larger than 1/sqrt(2)=0.707
  // check if parametric co-ordinate is clear
  if (thickdirlocmax < tol * thickdirlocsharp.Norm2())
    dserror(
        "could not clearly identify a parametric direction pointing along enforced thickness "
        "direction");

  ThicknessDirection thickdir = none;  // of actual element
  if (thick_index == 0)
  {
    thickdir = autor;
  }
  else if (thick_index == 1)
  {
    thickdir = autos;
  }
  else if (thick_index == 2)
  {
    thickdir = autot;
  }
  else
  {
    dserror("Trouble with thick_index=%g", thick_index);
  }

  // thickness-vector in parameter-space, has 1.0 in thickness-coord
  LINALG::Matrix<NUMDIM_SOH8, 1> loc_thickvec(true);
  loc_thickvec(thick_index) = 1.0;
  // thickness-vector in global coord is J times local thickness-vector
  // (X^A) = (J0_i^A)^T . (xi_i)
  LINALG::Matrix<NUMDIM_SOH8, 1> glo_thickvec;
  glo_thickvec.MultiplyTN(jac0, loc_thickvec);
  // return doubles of thickness-vector
  thickvec_.resize(3);
  thickvec_[0] = glo_thickvec(0);
  thickvec_[1] = glo_thickvec(1);
  thickvec_[2] = glo_thickvec(2);

  return thickdir;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::sosh8_gmshplotlabeledelement(const int LabelIds[NUMNOD_SOH8])
{
  std::stringstream filename;
  filename << "solidelement" << this->Id() << ".gmsh";
  std::ofstream f_system("solidelement.gmsh");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" One Solid Element \" {" << std::endl;
  gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(this->thickdir_, this) << std::endl;
  // plot vector from 1st node to 5th node which is parametric t-dir
  std::vector<double> X15(3);
  X15[0] = this->Nodes()[4]->X()[0] - this->Nodes()[0]->X()[0];
  X15[1] = this->Nodes()[4]->X()[1] - this->Nodes()[0]->X()[1];
  X15[2] = this->Nodes()[4]->X()[2] - this->Nodes()[0]->X()[2];
  gmshfilecontent << "VP(" << std::scientific << this->Nodes()[0]->X()[0] << ",";
  gmshfilecontent << std::scientific << this->Nodes()[0]->X()[1] << ",";
  gmshfilecontent << std::scientific << this->Nodes()[0]->X()[2] << ")";
  gmshfilecontent << "{" << std::scientific << X15[0] << "," << X15[1] << "," << X15[2] << "};"
                  << std::endl;
  gmshfilecontent << "};" << std::endl;
  gmshfilecontent << "View \" LabelIds \" {" << std::endl;
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    gmshfilecontent << "SP(" << std::scientific << this->Nodes()[i]->X()[0] << ",";
    gmshfilecontent << std::scientific << this->Nodes()[i]->X()[1] << ",";
    gmshfilecontent << std::scientific << this->Nodes()[i]->X()[2] << ")";
    gmshfilecontent << "{" << LabelIds[i] << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  gmshfilecontent << "View \" I order \" {" << std::endl;
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    gmshfilecontent << "SP(" << std::scientific << this->Nodes()[i]->X()[0] << ",";
    gmshfilecontent << std::scientific << this->Nodes()[i]->X()[1] << ",";
    gmshfilecontent << std::scientific << this->Nodes()[i]->X()[2] << ")";
    gmshfilecontent << "{" << i << "};" << std::endl;
  }
  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>>
DRT::ELEMENTS::So_sh8::sosh8_derivs_sdc()
{
  std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs(NUMGPT_SOH8);
  // (r,s,t) gp-locations of fully integrated linear 8-node Hex
  const double r[NUMNOD_SOH8] = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
  const double s[NUMNOD_SOH8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
  const double t[NUMNOD_SOH8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};
  //  const double t[NUMNOD_SOH8] = {0.0,0.0,0.0,0.0, 0.0, 0.0, 0.0, 0.0};
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    // df wrt to r for each node(0..7) at each gp [i]
    (derivs[i])(0, 0) = -(1.0 - s[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(0, 1) = (1.0 - s[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(0, 2) = (1.0 + s[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(0, 3) = -(1.0 + s[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(0, 4) = -(1.0 - s[i]) * (1.0 + t[i]) * 0.125;
    (derivs[i])(0, 5) = (1.0 - s[i]) * (1.0 + t[i]) * 0.125;
    (derivs[i])(0, 6) = (1.0 + s[i]) * (1.0 + t[i]) * 0.125;
    (derivs[i])(0, 7) = -(1.0 + s[i]) * (1.0 + t[i]) * 0.125;

    // df wrt to s for each node(0..7) at each gp [i]
    (derivs[i])(1, 0) = -(1.0 - r[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(1, 1) = -(1.0 + r[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(1, 2) = (1.0 + r[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(1, 3) = (1.0 - r[i]) * (1.0 - t[i]) * 0.125;
    (derivs[i])(1, 4) = -(1.0 - r[i]) * (1.0 + t[i]) * 0.125;
    (derivs[i])(1, 5) = -(1.0 + r[i]) * (1.0 + t[i]) * 0.125;
    (derivs[i])(1, 6) = (1.0 + r[i]) * (1.0 + t[i]) * 0.125;
    (derivs[i])(1, 7) = (1.0 - r[i]) * (1.0 + t[i]) * 0.125;

    // df wrt to t for each node(0..7) at each gp [i]
    (derivs[i])(2, 0) = -(1.0 - r[i]) * (1.0 - s[i]) * 0.125;
    (derivs[i])(2, 1) = -(1.0 + r[i]) * (1.0 - s[i]) * 0.125;
    (derivs[i])(2, 2) = -(1.0 + r[i]) * (1.0 + s[i]) * 0.125;
    (derivs[i])(2, 3) = -(1.0 - r[i]) * (1.0 + s[i]) * 0.125;
    (derivs[i])(2, 4) = (1.0 - r[i]) * (1.0 - s[i]) * 0.125;
    (derivs[i])(2, 5) = (1.0 + r[i]) * (1.0 - s[i]) * 0.125;
    (derivs[i])(2, 6) = (1.0 + r[i]) * (1.0 + s[i]) * 0.125;
    (derivs[i])(2, 7) = (1.0 - r[i]) * (1.0 + s[i]) * 0.125;
  }
  return derivs;
}
