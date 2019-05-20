/*----------------------------------------------------------------------------*/
/*!
\brief wall EAS

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
// macros


/*----------------------------------------------------------------------*/
// headers
#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/stvenantkirchhoff.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  setup of constant EAS data (private)                       mgit 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_eassetup(Epetra_SerialDenseMatrix& boplin0,
    Epetra_SerialDenseVector& F0,           // deformation gradient at origin
    Epetra_SerialDenseMatrix& xjm0,         // jacobian matrix at origin
    double& detJ0,                          // det of Jacobian at origin
    const Epetra_SerialDenseMatrix& xrefe,  // material element coords
    const Epetra_SerialDenseMatrix& xcure,  // current element coords
    const DRT::Element::DiscretizationType& distype)

{
  // derivatives at origin
  Epetra_SerialDenseMatrix deriv0;
  deriv0.Shape(2, NumNode());

  DRT::UTILS::shape_function_2D_deriv1(deriv0, 0.0, 0.0, distype);

  // compute jacobian matrix at origin
  memset(xjm0.A(), 0, xjm0.N() * xjm0.M() * sizeof(double));
  for (int k = 0; k < NumNode(); k++)
  {
    xjm0(0, 0) += deriv0(0, k) * xrefe(0, k);  // X,r += (X,r)^k
    xjm0(0, 1) += deriv0(0, k) * xrefe(1, k);  // Y,r += (Y,r)^k
    xjm0(1, 0) += deriv0(1, k) * xrefe(0, k);  // X,s += (X,s)^k
    xjm0(1, 1) += deriv0(1, k) * xrefe(1, k);  // Y,s += (Y,s)^k
  }

  /*------------------------------------------ determinant of jacobian ---*/
  detJ0 = xjm0[0][0] * xjm0[1][1] - xjm0[1][0] * xjm0[0][1];

  if (detJ0 < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");


  // compute boplin at origin (boplin0)

  int inode0;
  int dnode0;
  double xji0[2][2];
  /*---------------------------------------------- inverse of jacobian ---*/
  xji0[0][0] = xjm0(1, 1) / detJ0;
  xji0[0][1] = -xjm0(0, 1) / detJ0;
  xji0[1][0] = -xjm0(1, 0) / detJ0;
  xji0[1][1] = xjm0(0, 0) / detJ0;
  /*----------------------------- get operator boplin of global derivatives -*/
  /*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,X    0   |
       |   0    Nk,Y |
       | Nk,Y    0   |
       |  0     Nk,X |
  */
  for (inode0 = 0; inode0 < NumNode(); inode0++)
  {
    dnode0 = inode0 * 2;

    boplin0(0, dnode0 + 0) = deriv0(0, inode0) * xji0[0][0] + deriv0(1, inode0) * xji0[0][1];
    boplin0(1, dnode0 + 1) = deriv0(0, inode0) * xji0[1][0] + deriv0(1, inode0) * xji0[1][1];
    boplin0(2, dnode0 + 0) = boplin0(1, dnode0 + 1);
    boplin0(3, dnode0 + 1) = boplin0(0, dnode0 + 0);
  }


  // compute displacement-based deformation gradient at origin (F0)

  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,X  |
        |  1 + Uy,Y  |
        |      Ux,Y  |
        |      Uy,X  |
  */
  F0[0] = 1;
  F0[1] = 1;
  F0[2] = 0;
  F0[3] = 0;
  for (int inode = 0; inode < NumNode(); inode++)
  {
    F0[0] += boplin0(0, 2 * inode) * (xcure(0, inode) - xrefe(0, inode));
    F0[1] += boplin0(1, 2 * inode + 1) * (xcure(1, inode) - xrefe(1, inode));
    F0[2] += boplin0(2, 2 * inode) * (xcure(0, inode) - xrefe(0, inode));
    F0[3] += boplin0(3, 2 * inode + 1) * (xcure(1, inode) - xrefe(1, inode));
  } /* end of loop over nodes */


  return;
}  // end of w1_eassetup


/*----------------------------------------------------------------------*
 | get the enhanced deformation gradient and                            |
 | also the operators G, W0 and Z                   (private) mgit 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_call_defgrad_enh(Epetra_SerialDenseMatrix& F_enh,
    const Epetra_SerialDenseMatrix xjm0, const Epetra_SerialDenseMatrix xjm, const double detJ0,
    const double det, const Epetra_SerialDenseVector F0, const Epetra_SerialDenseMatrix alpha,
    const double e1, const double e2, Epetra_SerialDenseMatrix& G, Epetra_SerialDenseMatrix& W0,
    const Epetra_SerialDenseMatrix boplin0, Epetra_SerialDenseMatrix& Z)
{
  // EAS

  // The interpolation functions for the 4 EAS modes are based on
  //
  // Q1E4:
  //     M1 = r 0    M2 = 0 s    M3 = 0 0    M4 = 0 0
  //          0 0         0 0         r 0         0 s
  //
  // Q1ET4:
  //     M1 = r 0    M2 = 0 r    M3 = 0 0    M4 = 0 0
  //          0 0         0 0         s 0         0 s
  //
  //     M = M1*alpha1 + M2*alpha2 + M3*alpha3 + M4*alpha4
  //

  Epetra_SerialDenseMatrix M;
  M.Shape(2, 2);

  Epetra_SerialDenseMatrix M_temp;
  M_temp.Shape(2, 2);

  // fill up 4 EAS matrices at each GP
  if (eastype_ == eas_q1e4)
  {
    M(0, 0) = e1 * alpha(0, 0);  // r
    M(0, 1) = e2 * alpha(1, 0);  // s
    M(1, 0) = e1 * alpha(2, 0);  // r
    M(1, 1) = e2 * alpha(3, 0);  // s
  }
  else if (eastype_ == eas_q1et4)
  {
    M(0, 0) = e1 * alpha(0, 0);  // r
    M(0, 1) = e1 * alpha(1, 0);  // r
    M(1, 0) = e2 * alpha(2, 0);  // s
    M(1, 1) = e2 * alpha(3, 0);  // s
  }
  else
  {
    dserror("Cannot handle EAS type=%d", eastype_);
  }

  // inverse of jacobian matrix at element origin
  Epetra_SerialDenseMatrix xjm_inv0;
  xjm_inv0.Shape(2, 2);

  xjm_inv0(0, 0) = xjm0(1, 1) / detJ0;
  xjm_inv0(0, 1) = -xjm0(0, 1) / detJ0;
  xjm_inv0(1, 0) = -xjm0(1, 0) / detJ0;
  xjm_inv0(1, 1) = xjm0(0, 0) / detJ0;

  // build A-operator
  //
  // Q1E4
  //    A = det(J_o)/det(J) * J_o^T . sum_i^neas ( M_i alpha_i ) . J_o^{-T}
  // Q1ET4
  //    A = det(J_o)/det(J) * J_o . sum_i^neas ( M_i alpha_i ) . J_o^{-T}
  Epetra_SerialDenseMatrix A(2, 2);  // A operator
  M_temp.Multiply('N', 'T', 1.0, M, xjm_inv0, 0.0);
  if (eastype_ == eas_q1e4)
    A.Multiply('T', 'N', detJ0 / det, xjm0, M_temp, 0.0);
  else if (eastype_ == eas_q1et4)
    A.Multiply('N', 'N', detJ0 / det, xjm0, M_temp, 0.0);
  else
    dserror("Cannot handle EAS type=%d", eastype_);

  // enhanced deformation gradient at origin (four rows, one column)

  F_enh(0, 0) = A(0, 0) * F0(0) + A(1, 0) * F0(2);
  F_enh(1, 0) = A(1, 1) * F0(1) + A(0, 1) * F0(3);
  F_enh(2, 0) = A(0, 1) * F0(0) + A(1, 1) * F0(2);
  F_enh(3, 0) = A(1, 0) * F0(1) + A(0, 0) * F0(3);

  // get operator W0

  // write matrix A in a different way (matrix 4x4)

  Epetra_SerialDenseMatrix A_big;
  A_big.Shape(4, 4);  // entries are zero

  A_big(0, 0) = A(0, 0);
  A_big(0, 2) = A(1, 0);
  A_big(1, 1) = A(1, 1);
  A_big(1, 3) = A(0, 1);
  A_big(2, 0) = A(0, 1);
  A_big(2, 2) = A(1, 1);
  A_big(3, 1) = A(1, 0);
  A_big(3, 3) = A(0, 0);

  // multiplication A_big x boplin0

  W0.Multiply('N', 'N', 1.0, A_big, boplin0, 0.0);

  // calculation operators G and Z, therfore matrices A are needed
  // without alphas

  // vector M_ges, includes the matrices M1 to M4
  std::vector<Epetra_SerialDenseMatrix> M_ges(Wall1::neas_);

  // vector A_ges, includes the matrices A1 to A4
  std::vector<Epetra_SerialDenseMatrix> A_ges(4);

  for (int ieas = 0; ieas < Wall1::neas_; ieas++)
  {
    (M_ges[ieas]).Shape(2, 2);

    memset((M_ges[ieas]).A(), 0, (M_ges[ieas]).N() * (M_ges[ieas]).M() * sizeof(double));

    (A_ges[ieas]).Shape(2, 2);
  }

  // fill M-Matrixes, not including eas-parameters alpha

  if (eastype_ == eas_q1e4)
  {
    (M_ges[0])(0, 0) = e1;
    (M_ges[1])(0, 1) = e2;
    (M_ges[2])(1, 0) = e1;
    (M_ges[3])(1, 1) = e2;
  }
  else if (eastype_ == eas_q1et4)
  {
    (M_ges[0])(0, 0) = e1;
    (M_ges[1])(0, 1) = e1;
    (M_ges[2])(1, 0) = e2;
    (M_ges[3])(1, 1) = e2;
  }
  else
  {
    dserror("Cannot handle EAS type=%d", eastype_);
  }

  // declaration of matrix (4x4) without eas-parameter alpha

  Epetra_SerialDenseMatrix Awa_big;
  Awa_big.Shape(4, 4);

  // declaration of matrix WO without eas-parameter alpha

  Epetra_SerialDenseMatrix W0wa;
  W0wa.Shape(4, 2 * NumNode());


  for (int i = 0; i < Wall1::neas_; i++)  // loop over eas-parameter
  {
    M_temp.Multiply('N', 'T', 1.0, M_ges[i], xjm_inv0, 0.0);
    if (eastype_ == eas_q1e4)
      A_ges[i].Multiply('T', 'N', detJ0 / det, xjm0, M_temp, 0.0);
    else if (eastype_ == eas_q1et4)
      A_ges[i].Multiply('N', 'N', detJ0 / det, xjm0, M_temp, 0.0);
    else
      dserror("Cannot handle EAS type=%d", eastype_);

    // fill G-operator

    G(0, i) = A_ges[i](0, 0) * F0(0) + A_ges[i](1, 0) * F0(2);
    G(1, i) = A_ges[i](1, 1) * F0(1) + A_ges[i](0, 1) * F0(3);
    G(2, i) = A_ges[i](0, 1) * F0(0) + A_ges[i](1, 1) * F0(2);
    G(3, i) = A_ges[i](1, 0) * F0(1) + A_ges[i](0, 0) * F0(3);


    memset(Awa_big.A(), 0, Awa_big.N() * Awa_big.M() * sizeof(double));

    Awa_big(0, 0) = A_ges[i](0, 0);
    Awa_big(0, 2) = A_ges[i](1, 0);
    Awa_big(1, 1) = A_ges[i](1, 1);
    Awa_big(1, 3) = A_ges[i](0, 1);
    Awa_big(2, 0) = A_ges[i](0, 1);
    Awa_big(2, 2) = A_ges[i](1, 1);
    Awa_big(3, 1) = A_ges[i](1, 0);
    Awa_big(3, 3) = A_ges[i](0, 0);

    // calculate operator W0wa without eas-parameters alpha

    W0wa.Multiply('N', 'N', 1.0, Awa_big, boplin0, 0.0);

    // fill Z-operator

    for (int indof = 0; indof < NumNode(); indof = indof + 2)
    {
      Z(indof, i) = W0wa(0, indof);
      Z(indof + 1, i) = W0wa(1, indof + 1);
    }


  }  // end loop over eas parameter

  return;
}  // end of w1_call_fenh


/*----------------------------------------------------------------------*
 |total deformation gradient and green lagrange strain (private)mgit 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_call_defgrad_tot(const Epetra_SerialDenseMatrix& F_enh,
    Epetra_SerialDenseMatrix& F_tot, const Epetra_SerialDenseVector& F,
    Epetra_SerialDenseVector& strain)
{
  // total deformation gradient in matrix notation
  F_tot(0, 0) = F(0) + F_enh(0, 0);
  F_tot(0, 2) = F(2) + F_enh(2, 0);
  F_tot(1, 1) = F(1) + F_enh(1, 0);
  F_tot(1, 2) = F(3) + F_enh(3, 0);
  F_tot(2, 1) = F(2) + F_enh(2, 0);
  F_tot(2, 2) = F(0) + F_enh(0, 0);
  F_tot(3, 0) = F(3) + F_enh(3, 0);
  F_tot(3, 2) = F(1) + F_enh(1, 0);


  /*-----------------------calculate Green-Lagrange strain ------------*/
  strain[0] = 0.5 * (F_tot(0, 0) * F_tot(0, 0) + F_tot(3, 0) * F_tot(3, 0) - 1.0);
  strain[1] = 0.5 * (F_tot(0, 2) * F_tot(0, 2) + F_tot(3, 2) * F_tot(3, 2) - 1.0);
  strain[2] = 0.5 * (F_tot(0, 0) * F_tot(0, 2) + F_tot(3, 0) * F_tot(3, 2));
  strain[3] = strain[2];

  return;
}  // end of w1_call_defgrad_tot

/*-----------------------------------------------------------------------------*
 |first piola-kirchhoff stress vector                       (private)mgit 02/08|
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_stress_eas(const Epetra_SerialDenseMatrix& stress,
    const Epetra_SerialDenseMatrix& F_tot, Epetra_SerialDenseMatrix& p_stress)
{
  /*-------------reduce stress matrix-----------------------------------------*/

  LINALG::SerialDenseMatrix stress_red(3, 1, false);  // 2. piola-krichhoff (vector (3-dim))

  stress_red(0, 0) = stress(0, 0);  // S_11
  stress_red(1, 0) = stress(1, 1);  // S_22
  stress_red(2, 0) = stress(0, 2);  // S_12 (=S_21)

  /*-first piola-kirchhoff stress vector--------------------------------------*/
  p_stress.Multiply('N', 'N', 1.0, F_tot, stress_red, 0.0);

  return;
}  // end of w1_p_stress


/*-----------------------------------------------------------------------------*
| calculate stiffness matrix kdd                                     mgit 03/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_kdd(const Epetra_SerialDenseMatrix& boplin,
    const Epetra_SerialDenseMatrix& W0, const Epetra_SerialDenseMatrix& F_tot,
    const Epetra_SerialDenseMatrix& C, const Epetra_SerialDenseMatrix& stress,
    Epetra_SerialDenseMatrix& FCF, Epetra_SerialDenseMatrix& estif, const double fac)
{
  // contitutive matrix (3x3)
  LINALG::SerialDenseMatrix C_red(3, 3, false);
  C_red(0, 0) = C(0, 0);
  C_red(0, 1) = C(0, 1);
  C_red(0, 2) = C(0, 2);
  C_red(1, 0) = C(1, 0);
  C_red(1, 1) = C(1, 1);
  C_red(1, 2) = C(1, 2);
  C_red(2, 0) = C(2, 0);
  C_red(2, 1) = C(2, 1);
  C_red(2, 2) = C(2, 2);

  // BplusW = B+W0
  LINALG::SerialDenseMatrix BplusW(4, 2 * NumNode(), false);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 2 * NumNode(); j++) BplusW(i, j) = boplin(i, j) + W0(i, j);

  // Temp (4x3) = F*C
  LINALG::SerialDenseMatrix Temp(4, 3, true);
  Temp.Multiply('N', 'N', 1.0, F_tot, C_red, 0.0);

  // FCF^T (4x4) = Temp*F^T
  memset(FCF.A(), 0, FCF.N() * FCF.M() * sizeof(double));
  FCF.Multiply('N', 'T', 1.0, Temp, F_tot, 0.0);

  // Temp1 (4x8) = FCF^T * (B+W0)
  LINALG::SerialDenseMatrix Temp1(4, 2 * NumNode(), true);
  Temp1.Multiply('N', 'N', 1.0, FCF, BplusW, 0.0);

  // Temp3 (4x8) = S*(B+W0)
  LINALG::SerialDenseMatrix Temp3(4, 2 * NumNode(), true);
  Temp3.Multiply('N', 'N', 1.0, stress, BplusW, 0.0);

  // Kdd = (B+W0)^T*FCF^T*(B+W0) + (B+W0)^T*S*(B+W0)
  estif.Multiply('T', 'N', fac, BplusW, Temp1, 1.0);
  estif.Multiply('T', 'N', fac, BplusW, Temp3, 1.0);

  return;
}  // DRT::ELEMENTS::Wall1::w1_kdd


/*-----------------------------------------------------------------------------*
| calculate matrix kda                                               mgit 03/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_kda(const Epetra_SerialDenseMatrix& FCF,
    const Epetra_SerialDenseMatrix& W0, const Epetra_SerialDenseMatrix& boplin,
    const Epetra_SerialDenseMatrix& stress, const Epetra_SerialDenseMatrix& G,
    const Epetra_SerialDenseMatrix& Z, Epetra_SerialDenseMatrix& Kda,
    const Epetra_SerialDenseMatrix& p_stress, const double fac)
{
  // BplusW = B+W0
  LINALG::SerialDenseMatrix BplusW(4, 2 * NumNode(), false);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 2 * NumNode(); j++) BplusW(i, j) = boplin(i, j) + W0(i, j);

  // Temp1 = FCF^T*G
  LINALG::SerialDenseMatrix Temp1(4, Wall1::neas_, true);
  Temp1.Multiply('N', 'N', 1.0, FCF, G, 0.0);

  // Temp3 = S*(G)
  LINALG::SerialDenseMatrix Temp3(4, Wall1::neas_, true);
  Temp3.Multiply('N', 'N', 1.0, stress, G, 0.0);

  // Kda (8x4) = (B+W0)^T*FCF^T*G) + (B+W0)^T*S*G + PZ
  Kda.Multiply('T', 'N', fac, BplusW, Temp1, 1.0);
  Kda.Multiply('T', 'N', fac, BplusW, Temp3, 1.0);
  // Temp5 = fac * P*Z
  for (int i = 0; i < NumNode(); i++)
  {
    for (int ieas = 0; ieas < Wall1::neas_; ieas++)
    {
      Kda(i * 2, ieas) +=
          (p_stress(0, 0) * Z(i * 2, ieas) + p_stress(2, 0) * Z(i * 2 + 1, ieas)) * fac;
      Kda(i * 2 + 1, ieas) +=
          (p_stress(3, 0) * Z(i * 2, ieas) + p_stress(1, 0) * Z(i * 2 + 1, ieas)) * fac;
    }
  }

  return;
}  // DRT::ELEMENTS::Wall1::w1_kda


/*-----------------------------------------------------------------------------*
| calculate matrix kaa                                               mgit 03/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_kaa(const Epetra_SerialDenseMatrix& FCF,
    const Epetra_SerialDenseMatrix& stress, const Epetra_SerialDenseMatrix& G,
    Epetra_SerialDenseMatrix& Kaa, const double fac)
{
  // Temp1 = FCF*G
  LINALG::SerialDenseMatrix Temp1(4, Wall1::neas_, true);
  Temp1.Multiply('N', 'N', 1.0, FCF, G, 0.0);

  // Temp3 = S*G
  LINALG::SerialDenseMatrix Temp3(4, Wall1::neas_, true);
  Temp3.Multiply('N', 'N', 1.0, stress, G, 0.0);

  // Kaa = G^T*FCF^T*G + G^T*S*G
  Kaa.Multiply('T', 'N', fac, G, Temp1, 1.0);
  Kaa.Multiply('T', 'N', fac, G, Temp3, 1.0);

  return;
}  // DRT::ELEMENTS::Wall1::w1_kaa


/*-----------------------------------------------------------------------------*
| calculate internal forces fint(displacements u) and feas           mgit 03/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_fint_eas(const Epetra_SerialDenseMatrix& W0,
    const Epetra_SerialDenseMatrix& boplin, const Epetra_SerialDenseMatrix& G,
    const Epetra_SerialDenseMatrix& p_stress, Epetra_SerialDenseVector& intforce,
    Epetra_SerialDenseVector& feas, const double fac)

{
  // BplusW = B+W0
  LINALG::SerialDenseMatrix BplusW(4, 2 * NumNode(), false);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 2 * NumNode(); j++) BplusW(i, j) = boplin(i, j) + W0(i, j);

  // Temp1 (8x1) = (BL+W0)^T*p_stress
  LINALG::SerialDenseMatrix Temp1(2 * NumNode(), 1, true);
  Temp1.Multiply('T', 'N', 1.0, BplusW, p_stress, 0.0);

  for (int i = 0; i < 2 * NumNode(); i++) intforce(i) += fac * Temp1(i, 0);

  // Temp2 = G^T*p_stress
  LINALG::SerialDenseMatrix Temp2(Wall1::neas_, 1, true);
  Temp2.Multiply('T', 'N', 1.0, G, p_stress, 0.0);

  for (int i = 0; i < Wall1::neas_; i++) feas(i) += fac * Temp2(i, 0);

  return;
}  // DRT::ELEMENTS::Wall1::w1_fint_eas


/*----------------------------------------------------------------------*/
