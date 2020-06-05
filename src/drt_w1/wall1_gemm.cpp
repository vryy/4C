/*----------------------------------------------------------------------------*/
/*! \file
\brief Routines for generalised energy-momentum method
       FixMe This file is currently unsupported in the new structural
       time integration, since the corresponding time integration is
       not yet (re-)implemented.                   04/17 hiermeier

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include "wall1.H"



/*======================================================================*/
/* evaluate the element forces and stiffness and mass for GEMM */
void DRT::ELEMENTS::Wall1::FintStiffMassGEMM(Teuchos::ParameterList& params,
    const std::vector<int>& lm, const std::vector<double>& dispo, const std::vector<double>& disp,
    const std::vector<double>& residual, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force,
    Epetra_SerialDenseMatrix* elestress, Epetra_SerialDenseMatrix* elestrain,
    Teuchos::RCP<const MAT::Material> material, const INPAR::STR::StressType iostress,
    const INPAR::STR::StrainType iostrain)
{
  // constants
  // element porperties
  const int numnode = NumNode();
  const int edof = numnode * Wall1::noddof_;
  const DiscretizationType distype = Shape();
  // Gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);
  // GEMM coefficients
  const double gemmalphaf = params.get<double>("alpha f");
  const double gemmxi = params.get<double>("xi");
  // density if mass is calculated
  const double density = (massmatrix) ? (material->Density()) : 0.0;

  // general arrays
  Epetra_SerialDenseVector shpfct(numnode);  // shape functions at Gauss point
  Epetra_SerialDenseMatrix shpdrv(
      Wall1::numdim_, numnode);  // parametric derivatives of shape funct. at Gauss point
  Epetra_SerialDenseMatrix Xjm(
      Wall1::numdim_, Wall1::numdim_);  // material-to-parameter-space Jacobian
  double Xjdet;                         // determinant of #Xjm
  Epetra_SerialDenseMatrix boplin(4, edof);

  Epetra_SerialDenseVector Fuvo(4);  // disp-based def.grad. vector at t_{n}
  Epetra_SerialDenseVector Fuv(4);   // disp-based def.grad. vector at t_{n+1}

  LINALG::SerialDenseVector Evo(4, false);  // Green-Lagrange strain vector at t_{n}
  LINALG::SerialDenseVector Ev(4, false);   // Green-Lagrange strain vector at t_{n+1}
  LINALG::SerialDenseVector& Evm = Evo;     // Green-Lagrange mid-strain vector

  Epetra_SerialDenseMatrix Xe(Wall1::numdim_, numnode);  // material/initial element co-ordinates
  Epetra_SerialDenseMatrix xeo(
      Wall1::numdim_, numnode);  // spatial/current element co-ordinates at t_{n}
  Epetra_SerialDenseMatrix xe(
      Wall1::numdim_, numnode);  // spatial/current element co-ordinates at t_{n+1}
  LINALG::SerialDenseMatrix bopo(Wall1::numstr_, edof, false);  // non-linear B-op at t_{n}
  LINALG::SerialDenseMatrix bop(Wall1::numstr_, edof, false);   // non-linear B-op at t_{n+1}
  LINALG::SerialDenseMatrix& bopm = bopo;                       // non-linear mid-B-op
  Epetra_SerialDenseMatrix Smm(4, 4);  // 2nd Piola-Kirchhoff mid-stress matrix  // CHECK THIS:
                                       // STRESS MATRIX SHOULD NOT EXIST IN EFFICIENT CODE
  Epetra_SerialDenseMatrix C(4, 4);

  // for EAS, in any case declare variables, sizes etc. only allocated in EAS version
  Epetra_SerialDenseMatrix* alphao = NULL;     // EAS alphas at t_{n}
  Epetra_SerialDenseMatrix* alpha = NULL;      // EAS alphas at t_{n+1}
  Epetra_SerialDenseMatrix* oldfeas = NULL;    // EAS history
  Epetra_SerialDenseMatrix* oldKaainv = NULL;  // EAS history
  Epetra_SerialDenseMatrix* oldKda = NULL;     // EAS history
  Epetra_SerialDenseMatrix* oldKad = NULL;     // EAS history
  Epetra_SerialDenseMatrix Fenhvo;             // EAS matrix Fenhv
  Epetra_SerialDenseMatrix Fenhv;              // EAS matrix Fenhv
  LINALG::SerialDenseMatrix Fmo;               // total def.grad. matrix at t_{n}
  LINALG::SerialDenseMatrix Fm;                // total def.grad. matrix at t_{n+1}
  LINALG::SerialDenseMatrix& Fmm = Fmo;        // total mid-def.grad. matrix
  Epetra_SerialDenseMatrix Pvmm;               // first Piola-Kirchhoff stress vector
  Epetra_SerialDenseMatrix Xjm0;               // Jacobian Matrix (origin)
  double Xjdet0;                               // determinant of #Xjm0
  Epetra_SerialDenseVector Fuv0o;              // deformation gradient at origin at t_{n}
  Epetra_SerialDenseVector Fuv0;               // deformation gradient at origin at t_{n+1}
  Epetra_SerialDenseMatrix boplin0;            // B-operator (origin)
  LINALG::SerialDenseMatrix W0o;               // W-operator (origin) at t_{n}
  LINALG::SerialDenseMatrix W0;                // W-operator (origin) at t_{n+1}
  LINALG::SerialDenseMatrix& W0m = W0o;        // mid-W-operator (origin)
  LINALG::SerialDenseMatrix Go;                // G-operator at t_{n}
  LINALG::SerialDenseMatrix G;                 // G-operator at t_{n+1}
  LINALG::SerialDenseMatrix& Gm = Go;          // mid-G-operator
  Epetra_SerialDenseMatrix Z;                  // Z-operator
  Epetra_SerialDenseMatrix FmCF;               // FCF^T
  Epetra_SerialDenseMatrix Kda;                // EAS matrix Kda
  Epetra_SerialDenseMatrix Kad;                // EAS matrix Kad
  Epetra_SerialDenseMatrix Kaa;                // EAS matrix Kaa
  Epetra_SerialDenseVector feas;               // EAS portion of internal forces

  // element co-ordinates
  for (int k = 0; k < numnode; ++k)
  {
    Xe(0, k) = Nodes()[k]->X()[0];
    Xe(1, k) = Nodes()[k]->X()[1];
    xeo(0, k) = Xe(0, k) + dispo[k * Wall1::noddof_ + 0];
    xeo(1, k) = Xe(1, k) + dispo[k * Wall1::noddof_ + 1];
    xe(0, k) = Xe(0, k) + disp[k * Wall1::noddof_ + 0];
    xe(1, k) = Xe(1, k) + disp[k * Wall1::noddof_ + 1];
  }

  // set-up EAS parameters
  if (iseas_)
  {
    // allocate EAS quantities
    Fenhvo.Shape(4, 1);
    Fenhv.Shape(4, 1);
    Fmo.Shape(4, 3);
    Fm.Shape(4, 3);
    Pvmm.Shape(4, 1);
    Xjm0.Shape(2, 2);
    Fuv0o.Size(4);
    Fuv0.Size(4);
    boplin0.Shape(4, edof);
    W0o.Shape(4, edof);
    W0.Shape(4, edof);
    Go.Shape(4, Wall1::neas_);
    G.Shape(4, Wall1::neas_);
    Z.Shape(edof, Wall1::neas_);
    FmCF.Shape(4, 4);
    Kda.Shape(edof, Wall1::neas_);
    Kad.Shape(Wall1::neas_, edof);
    Kaa.Shape(Wall1::neas_, Wall1::neas_);
    feas.Size(Wall1::neas_);

    // EAS Update of alphas:
    // the current alphas are (re-)evaluated out of
    // Kaa and Kda of previous step to avoid additional element call.
    // This corresponds to the (innermost) element update loop
    // in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    alphao =
        data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // get alpha of last converged state
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // get alpha of previous iteration

    // get stored EAS history
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    oldKad = data_.GetMutable<Epetra_SerialDenseMatrix>("Kad");
    if ((not alpha) or (not oldKaainv) or (not oldKda) or (not oldKad) or (not oldfeas))
      dserror("Missing EAS history-data");

    // we need the (residual) displacement at the previous step
    Epetra_SerialDenseVector res_d(edof);
    for (int i = 0; i < edof; ++i)
    {
      res_d(i) = residual[i];
    }

    // update enhanced strain scales by condensation
    // add Kda . res_d to feas
    (*oldfeas).Multiply('N', 'N', 1.0, (*oldKad), res_d, 1.0);
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    (*alpha).Multiply('N', 'N', -1.0, (*oldKaainv), (*oldfeas), 1.0);

    // derivatives at origin
    DRT::UTILS::shape_function_2D_deriv1(shpdrv, 0.0, 0.0, distype);
    // material-to-parameter space Jacobian at origin
    w1_jacobianmatrix(Xe, shpdrv, Xjm0, &Xjdet0, numnode);
    // calculate linear B-operator at origin
    w1_boplin(boplin0, shpdrv, Xjm0, Xjdet0, numnode);
    // displ.-based def.grad. at origin
    w1_defgrad(Fuv0o, Ev, Xe, xeo, boplin0, numnode);  // at t_{n}
    w1_defgrad(Fuv0, Ev, Xe, xe, boplin0, numnode);    // at t_{n+1}
  }

  // integration loops over element domain
  for (int ip = 0; ip < intpoints.nquad; ++ip)
  {
    // Gaussian point and weight at it
    const double xi1 = intpoints.qxg[ip][0];
    const double xi2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(shpfct, xi1, xi2, distype);
    DRT::UTILS::shape_function_2D_deriv1(shpdrv, xi1, xi2, distype);

    // compute Jacobian matrix
    w1_jacobianmatrix(Xe, shpdrv, Xjm, &Xjdet, numnode);

    // integration factor
    double fac = wgt * Xjdet * thickness_;

    // compute mass matrix
    if (massmatrix)
    {
      double facm = fac * density;
      for (int a = 0; a < numnode; a++)
      {
        for (int b = 0; b < numnode; b++)
        {
          (*massmatrix)(2 * a, 2 * b) += facm * shpfct(a) * shpfct(b);         /* a,b even */
          (*massmatrix)(2 * a + 1, 2 * b + 1) += facm * shpfct(a) * shpfct(b); /* a,b odd  */
        }
      }
    }

    // calculate linear B-operator
    w1_boplin(boplin, shpdrv, Xjm, Xjdet, numnode);

    // calculate defgrad F^u, Green-Lagrange-strain E^u
    w1_defgrad(Fuvo, Evo, Xe, xeo, boplin, numnode);  // at t_{n}
    w1_defgrad(Fuv, Ev, Xe, xe, boplin, numnode);     // at t_{n+1}

    // calculate non-linear B-operator in current configuration
    w1_boplin_cure(bopo, boplin, Fuvo, Wall1::numstr_,
        edof);  // at t_{n} // CHECK THIS: NOT SURE IF bopo NEEDED
    w1_boplin_cure(bop, boplin, Fuv, Wall1::numstr_, edof);  // at t_{n+1}

    // EAS: The deformation gradient is enhanced
    if (iseas_)
    {
      // calculate the enhanced deformation gradient and
      // also the operators G, W0 and Z
      w1_call_defgrad_enh(Fenhvo, Xjm0, Xjm, Xjdet0, Xjdet, Fuv0o, *alphao, xi1, xi2, Go, W0o,
          boplin0, Z);  // at t_{n}
      w1_call_defgrad_enh(Fenhv, Xjm0, Xjm, Xjdet0, Xjdet, Fuv0, *alpha, xi1, xi2, G, W0, boplin0,
          Z);  // at t_{n+1}

      // total deformation gradient F, and total Green-Lagrange-strain E
      w1_call_defgrad_tot(Fenhvo, Fmo, Fuvo, Evo);  // at t_{n}
      w1_call_defgrad_tot(Fenhv, Fm, Fuv, Ev);      // at t_{n+1}
    }

    // mid-def.grad.
    // F_m = (1.0-gemmalphaf)*F_{n+1} + gemmalphaf*F_{n}
    Fmm.Update((1.0 - gemmalphaf), Fm, gemmalphaf);  // remember same pointer: F_m = F_{n}

    // non-linear mid-B-operator
    // B_m = (1.0-gemmalphaf)*B_{n+1} + gemmalphaf*B_{n}
    bopm.Update((1.0 - gemmalphaf), bop, gemmalphaf);  // remember same pointer B_m = B_{n}

    // mid-strain GL vector
    // E_m = (1.0-gemmalphaf+gemmxi)*E_{n+1} + (gemmalphaf-gemmxi)*E_n
    Evm.Update((1.0 - gemmalphaf + gemmxi), Ev,
        (gemmalphaf - gemmxi));  // remember same pointer: E_m = E_{n}

    // extra mid-quantities for case of EAS
    if (iseas_)
    {
      // mid-G-operator : G_m = 0.5*G_{n+1} + 0.5*G_{n}
      Gm.Update(0.5, G, 0.5);  // remember same pointer: G_m = G_{n}

      // mid-W0-operator : W_{0;m} = 0.5*W_{0;n+1} + 0.5*W_{0;n}
      W0m.Update(0.5, W0, 0.5);  // remember same pointer: W0_m = W0_{n}
    }

    // call material law
    if (material->MaterialType() == INPAR::MAT::m_stvenant)
      w1_call_matgeononl(Evm, Smm, C, Wall1::numstr_, material, params, ip);
    else
      dserror("It must be St.Venant-Kirchhoff material.");

    // return Gauss point strains (only in case of stress/strain output)
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("no strain data available");
        for (int i = 0; i < Wall1::numstr_; ++i) (*elestrain)(ip, i) = Ev(i);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      case INPAR::STR::strain_ea:
      default:
        dserror("requested strain type not supported");
        break;
    }

    // return stresses at Gauss points (only in case of stress/strain output)
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("no stress data available");
        (*elestress)(ip, 0) = Smm(0, 0);  // 2nd Piola-Kirchhoff stress S_{11}
        (*elestress)(ip, 1) = Smm(1, 1);  // 2nd Piola-Kirchhoff stress S_{22}
        (*elestress)(ip, 2) = Smm(0, 2);  // 2nd Piola-Kirchhoff stress S_{12}
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("no stress data available");
        if (iseas_)
          StressCauchy(ip, Fm(0, 0), Fm(1, 1), Fm(1, 1), Fm(1, 2), Smm, elestress);
        else
          StressCauchy(ip, Fuv[0], Fuv[1], Fuv[2], Fuv[3], Smm, elestress);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not supported");
        break;
    }

    // stiffness and internal force
    if (iseas_)
    {
      // first mid-mid Piola-Kirchhoff stress vector P_{mm} = F_m . S_m
      w1_stress_eas(Smm, (Fmm), Pvmm);

      // stiffness matrix kdd
      if (stiffmatrix)
        TangFintByDispGEMM(
            gemmalphaf, gemmxi, fac, boplin, W0m, W0, Fmm, Fm, C, Smm, FmCF, *stiffmatrix);
      // matrix kda
      TangFintByEnhGEMM(gemmalphaf, gemmxi, fac, boplin, W0m, FmCF, Smm, G, Z, Pvmm, Kda);
      // matrix kad (this is NOT kda, because GEMM produces non-symmetric tangent!)
      TangEconByDispGEMM(gemmalphaf, gemmxi, fac, boplin, W0, FmCF, Smm, Gm, Z, Pvmm, Kad);
      // matrix kaa
      TangEconByEnhGEMM(gemmalphaf, gemmxi, fac, FmCF, Smm, G, Gm, Kaa);
      // nodal forces
      if (force) w1_fint_eas(W0m, boplin, Gm, Pvmm, *force, feas, fac);
    }
    else
    {
      // element stiffness matrix constribution at current Gauss point
      if (stiffmatrix)
        TangFintByDispGEMM(gemmalphaf, gemmxi, fac, bopm, bop, C, boplin, Smm, *stiffmatrix);

      // nodal forces fi from integration of stresses
      if (force) w1_fint(Smm, bopm, *force, fac, edof);
    }

  }  // for (int ip=0; ip<totngp; ++ip)


  // EAS technology: static condensation
  // subtract EAS matrices from disp-based Kdd to "soften" element
  if ((iseas_) and (force) and (stiffmatrix))
  {
    // we need the inverse of Kaa
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix(Kaa);
    solve_for_inverseKaa.Invert();

    Epetra_SerialDenseMatrix KdaKaa(edof, Wall1::neas_);  // temporary Kda.Kaa^{-1}
    KdaKaa.Multiply('N', 'N', 1.0, Kda, Kaa, 1.0);

    // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
    if (stiffmatrix) (*stiffmatrix).Multiply('N', 'N', -1.0, KdaKaa, Kad, 1.0);

    // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
    if (force) (*force).Multiply('N', 'N', -1.0, KdaKaa, feas, 1.0);

    // store current EAS data in history
    for (int i = 0; i < Wall1::neas_; ++i)
      for (int j = 0; j < Wall1::neas_; ++j) (*oldKaainv)(i, j) = Kaa(i, j);

    for (int i = 0; i < edof; ++i)
      for (int j = 0; j < Wall1::neas_; ++j)
      {
        (*oldKda)(i, j) = Kda(i, j);
        (*oldfeas)(j, 0) = feas(j);
      }
  }

  // good Bye
  return;
}

/*======================================================================*/
/* elastic and initial displacement stiffness */
void DRT::ELEMENTS::Wall1::TangFintByDispGEMM(const double& alphafgemm, const double& xigemm,
    const double& fac, const Epetra_SerialDenseMatrix& bopm, const Epetra_SerialDenseMatrix& bopn,
    const Epetra_SerialDenseMatrix& C, const Epetra_SerialDenseMatrix& boplin,
    const Epetra_SerialDenseMatrix& Smm, Epetra_SerialDenseMatrix& estif)
{
  // constants
  const int nd = Wall1::noddof_ * NumNode();  // number of element DOFs
  const int numeps = Wall1::numnstr_;

  // elastic and initial displacement stiffness
  // perform B_m^T . C . B_{n+1}, whereas B_{n+1} = F_{n+1}^T . B_L */
  {
    const double faceu = (1.0 - alphafgemm + xigemm) * fac;
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < nd; j++)
        for (int k = 0; k < numeps; k++)
          for (int m = 0; m < numeps; m++) estif(i, j) += bopm(k, i) * C(k, m) * bopn(m, j) * faceu;
  }

  // geometric stiffness part
  // perform B_L^T * S_m * B_L
  {
    const double fackg = (1.0 - alphafgemm) * fac;
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < nd; j++)
        for (int r = 0; r < numeps; r++)
          for (int m = 0; m < numeps; m++)
            estif(i, j) += boplin(r, i) * Smm(r, m) * boplin(m, j) * fackg;
  }

  // see you
  return;
}  // TangFintByDispGEMM

/*======================================================================*/
/* calcuate tangent (f_{int;m}),d */
void DRT::ELEMENTS::Wall1::TangFintByDispGEMM(const double& alphafgemm, const double& xigemm,
    const double& fac, const Epetra_SerialDenseMatrix& boplin, const Epetra_SerialDenseMatrix& W0m,
    const Epetra_SerialDenseMatrix& W0, const Epetra_SerialDenseMatrix& Fmm,
    const Epetra_SerialDenseMatrix& Fm, const Epetra_SerialDenseMatrix& C,
    const Epetra_SerialDenseMatrix& Smm, Epetra_SerialDenseMatrix& FmCF,
    Epetra_SerialDenseMatrix& estif)
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

  // FdotC (4 x 3) : F_m . C
  LINALG::SerialDenseMatrix FmC(4, 3, true);
  FmC.Multiply('N', 'N', 1.0, Fmm, C_red, 0.0);

  // FmCF (4 x 4) : ( F_m . C ) . F_{n+1}^T
  memset(FmCF.A(), 0, FmCF.N() * FmCF.M() * sizeof(double));
  FmCF.Multiply('N', 'T', 1.0, FmC, Fm, 0.0);

  // BplusW (4 x edof) :  B_L + W0_{n+1}
  LINALG::SerialDenseMatrix BplusW(4, 2 * NumNode(), true);
  BplusW += boplin;  // += B_L
  BplusW += W0;      // += W0_{n+1}

  // FmCFBW (4 x 8) : (Fm . C . F_{n+1}^T) . (B_L + W0_{n+1})
  LINALG::SerialDenseMatrix FmCFBW(4, 2 * NumNode(), true);
  FmCFBW.Multiply('N', 'N', 1.0, FmCF, BplusW, 0.0);

  // SmBW (4 x 8) : S_m . (B_L + W0_{n+1})
  LINALG::SerialDenseMatrix SmBW(4, 2 * NumNode(), true);
  SmBW.Multiply('N', 'N', 1.0, Smm, BplusW, 0.0);

  // BplusW (4 x 8) :  B_L + W0_{m}
  BplusW.Update(1.0, boplin, 0.0);
  BplusW.Update(1.0, W0m, 1.0);

  // k_{dd} (8 x 8) :
  // k_{dd} += fac * (B_L + W0_m)^T . (Fm . C . F_{n+1}^T) . (B_L + W0_m)
  estif.Multiply('T', 'N', (1.0 - alphafgemm + xigemm) * fac, BplusW, FmCFBW, 1.0);
  // k_{dd} += fac * (B_L+W0_m)^T . S_m . (B_L+W0_{n+1})
  estif.Multiply('T', 'N', (1.0 - alphafgemm) * fac, BplusW, SmBW, 1.0);

  // that's it
  return;
}  // DRT::ELEMENTS::Wall1::w1_kdd

/*======================================================================*/
/* calculate tangent (f_{int;m}),alpha */
void DRT::ELEMENTS::Wall1::TangFintByEnhGEMM(const double& alphafgemm, const double& xigemm,
    const double& fac, const Epetra_SerialDenseMatrix& boplin, const Epetra_SerialDenseMatrix& W0m,
    const Epetra_SerialDenseMatrix& FmCF, const Epetra_SerialDenseMatrix& Smm,
    const Epetra_SerialDenseMatrix& G, const Epetra_SerialDenseMatrix& Z,
    const Epetra_SerialDenseMatrix& Pvmm, Epetra_SerialDenseMatrix& kda)
{
  // FmCFG (4 x 4) : (F_m . C . F_{n+1}^T) . G_{n+1}
  LINALG::SerialDenseMatrix FmCFG(4, Wall1::neas_, true);
  FmCFG.Multiply('N', 'N', 1.0, FmCF, G, 0.0);

  // SmG (4 x 4) : S_m . G_{n+1}
  LINALG::SerialDenseMatrix SmG(4, Wall1::neas_, true);
  SmG.Multiply('N', 'N', 1.0, Smm, G, 0.0);

  // BplusW (4 x 8) :  B_L + W0_{m}
  LINALG::SerialDenseMatrix BplusW(4, 2 * NumNode(), true);
  BplusW += boplin;
  BplusW += W0m;

  // k_{da} (8 x 4) :
  // k_{da} += fac * (B_l+W0_m)^T . (F_m . C . F_{n+1}^T) . G_{n+1}
  kda.Multiply('T', 'N', (1.0 - alphafgemm + xigemm) * fac, BplusW, FmCFG, 1.0);
  // k_{da} += fac * (B_l+W0_m)^T . S_m . G_{n+1}
  kda.Multiply('T', 'N', (1.0 - alphafgemm) * fac, BplusW, SmG, 1.0);
  // k_{da} += fac * \bar{\bar{P}}_{mm} . Z_{n+1}
  for (int i = 0; i < NumNode(); i++)
  {
    for (int ieas = 0; ieas < Wall1::neas_; ieas++)
    {
      kda(i * 2, ieas) +=
          0.5 * fac * (Pvmm(0, 0) * Z(i * 2, ieas) + Pvmm(2, 0) * Z(i * 2 + 1, ieas));
      kda(i * 2 + 1, ieas) +=
          0.5 * fac * (Pvmm(3, 0) * Z(i * 2, ieas) + Pvmm(1, 0) * Z(i * 2 + 1, ieas));
    }
  }

  // ciao
  return;
}


/*======================================================================*/
/* calculate tangent (s_m),d */
void DRT::ELEMENTS::Wall1::TangEconByDispGEMM(const double& alphafgemm, const double& xigemm,
    const double& fac, const Epetra_SerialDenseMatrix& boplin, const Epetra_SerialDenseMatrix& W0,
    const Epetra_SerialDenseMatrix& FmCF, const Epetra_SerialDenseMatrix& Smm,
    const Epetra_SerialDenseMatrix& Gm, const Epetra_SerialDenseMatrix& Z,
    const Epetra_SerialDenseMatrix& Pvmm, Epetra_SerialDenseMatrix& kad)
{
  // BplusW (4 x 8) :  B_L + W0_{n+1}
  LINALG::SerialDenseMatrix BplusW(4, 2 * NumNode(), true);
  BplusW += boplin;
  BplusW += W0;

  // FmCFBW (4 x 8) : (F_m . C . F_{n+1}^T) . (B_L + W0_{n+1})
  LINALG::SerialDenseMatrix FmCFBW(4, 2 * NumNode(), true);
  FmCFBW.Multiply('N', 'N', 1.0, FmCF, BplusW, 0.0);

  // SmGBW (4 x 8) : S_m . G_{n+1} . (B_L + W0_{n+1})
  LINALG::SerialDenseMatrix SmGBW(4, 2 * NumNode(), true);
  SmGBW.Multiply('N', 'N', 1.0, Smm, BplusW, 0.0);

  // k_{ad} (4 x 8) :
  // k_{ad} += fac * G_{m}^T . (F_m . C . F_{n+1}^T) . (B_lin+W0_{n+1})
  kad.Multiply('T', 'N', (1.0 - alphafgemm + xigemm) * fac, Gm, FmCFBW, 1.0);
  // k_{ad} += fac *  G_{m}^T . S_m . (B_l+W0_{n+1})^T
  kad.Multiply('T', 'N', (1.0 - alphafgemm) * fac, Gm, SmGBW, 1.0);
  // k_{ad} += fac * (\bar{\bar{P}}_{mm} . Z_{n+1})^T
  for (int i = 0; i < NumNode(); i++)
  {
    for (int ieas = 0; ieas < Wall1::neas_; ieas++)
    {
      kad(ieas, i * 2) +=
          0.5 * fac * (Pvmm(0, 0) * Z(i * 2, ieas) + Pvmm(2, 0) * Z(i * 2 + 1, ieas));
      kad(ieas, i * 2 + 1) +=
          0.5 * fac * (Pvmm(3, 0) * Z(i * 2, ieas) + Pvmm(1, 0) * Z(i * 2 + 1, ieas));
    }
  }

  // ciao
  return;
}

/*======================================================================*/
/* calculate tangent (s_m),alpha */
void DRT::ELEMENTS::Wall1::TangEconByEnhGEMM(const double& alphafgemm, const double& xigemm,
    const double& fac, const Epetra_SerialDenseMatrix& FmCF, const Epetra_SerialDenseMatrix& Smm,
    const Epetra_SerialDenseMatrix& G, const Epetra_SerialDenseMatrix& Gm,
    Epetra_SerialDenseMatrix& kaa)
{
  // FmCFG : (F_m . C . F_{n+1}^T) . G_{n+1}
  LINALG::SerialDenseMatrix FmCFG(4, Wall1::neas_, true);
  FmCFG.Multiply('N', 'N', 1.0, FmCF, G, 0.0);

  // SmG : S_m  . G_{n+1}
  LINALG::SerialDenseMatrix SmG(4, Wall1::neas_, true);
  SmG.Multiply('N', 'N', 1.0, Smm, G, 0.0);

  // k_{aa} (4 x 4) :
  // k_{aa} += fac * G_m^T . (F_m . C . F_{n+1}^T) . G_{n+1}
  kaa.Multiply('T', 'N', (1.0 - alphafgemm + xigemm) * fac, Gm, FmCFG, 1.0);
  // k_{aa} += fac * G_m^T . S_m . G_{n+1}
  kaa.Multiply('T', 'N', (1.0 - alphafgemm) * fac, Gm, SmG, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
