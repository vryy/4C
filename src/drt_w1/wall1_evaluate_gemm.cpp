/*----------------------------------------------------------------------*/
/*!
\file wall1_evaluate_gemm.cpp
\brief Routines for generalised energy-momentum method

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET
#ifdef D_WALL1

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include "wall1.H"



/*======================================================================*/
/* evaluate the element forces and stiffness and mass for GEMM */
void DRT::ELEMENTS::Wall1::FintStiffMassGEMM(
  const ParameterList& params,
  const std::vector<int>& lm,
  const std::vector<double>& dispo,
  const std::vector<double>& disp,
  const std::vector<double>& residual,
  Epetra_SerialDenseMatrix* stiffmatrix,
  Epetra_SerialDenseMatrix* massmatrix,
  Epetra_SerialDenseVector* force,
  Epetra_SerialDenseMatrix* elestress,
  Epetra_SerialDenseMatrix* elestrain,
  struct _MATERIAL* material,
  const bool cauchy
)
{
  // constants
  // element porperties
  const int numnode = NumNode();
  const int edof = numnode * Wall1::noddof_;
  const DiscretizationType distype = Shape();
  // Gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints = getIntegrationPoints2D(gaussrule_);
  // GEMM coefficients
  const double gemmalphaf = params.get<double>("alpha f");
  const double gemmxi = params.get<double>("xi");

  // BLAS dummy
  Epetra_BLAS::Epetra_BLAS blas;

  // general arrays
  Epetra_SerialDenseVector shpfct(numnode);  // shape functions at Gauss point
  Epetra_SerialDenseMatrix shpdrv(Wall1::numdim_,numnode);  // parametric derivatives of shape funct. at Gauss point
  Epetra_SerialDenseMatrix Xjm(Wall1::numdim_,Wall1::numdim_);  // material-to-parameter-space Jacobian
  double Xjdet;  // determinant of #Xjm
  Epetra_SerialDenseMatrix boplin(4,edof);

  Epetra_SerialDenseVector Fuvo(4);  // disp-based def.grad. vector at t_{n}
  Epetra_SerialDenseVector Fuv(4);  // disp-based def.grad. vector at t_{n+1}

  Epetra_SerialDenseVector Evo(4);  // Green-Lagrange strain vector at t_{n}
  Epetra_SerialDenseVector Ev(4);  // Green-Lagrange strain vector at t_{n+1}
  Epetra_SerialDenseVector& Evm = Evo;  // Green-Lagrange mid-strain vector

  Epetra_SerialDenseMatrix Xe(Wall1::numdim_,numnode);  // material/initial element co-ordinates
  Epetra_SerialDenseMatrix xeo(Wall1::numdim_,numnode);  // spatial/current element co-ordinates at t_{n}
  Epetra_SerialDenseMatrix xe(Wall1::numdim_,numnode);  // spatial/current element co-ordinates at t_{n+1}
  Epetra_SerialDenseMatrix bopo(Wall1::numstr_,edof);  // non-linear B-op at t_{n}
  Epetra_SerialDenseMatrix bop(Wall1::numstr_,edof);  // non-linear B-op at t_{n+1}
  Epetra_SerialDenseMatrix& bopm = bopo;  // non-linear mid-B-op
  Epetra_SerialDenseMatrix Smm(4,4);  // 2nd Piola-Kirchhoff mid-stress matrix  // CHECK THIS: STRESS MATRIX SHOULD NOT EXIST IN EFFICIENT CODE
  Epetra_SerialDenseMatrix C(4,4);

  // for EAS, in any case declare variables, sizes etc. only allocated in EAS version
  Epetra_SerialDenseMatrix* alphao;  // EAS alphas at t_{n}
  Epetra_SerialDenseMatrix* alpha;  // EAS alphas at t_{n+1}
  Epetra_SerialDenseMatrix* Fenhvo;  // EAS matrix Fenhv
  Epetra_SerialDenseMatrix* Fenhv;  // EAS matrix Fenhv
  Epetra_SerialDenseMatrix* Fmo;  // total def.grad. matrix at t_{n}
  Epetra_SerialDenseMatrix* Fm;  // total def.grad. matrix at t_{n+1}
  Epetra_SerialDenseMatrix* Fmm;  // total mid-def.grad. matrix
  Epetra_SerialDenseMatrix* Pvmm;  // first Piola-Kirchhoff stress vector
  Epetra_SerialDenseMatrix* Xjm0;  // Jacobian Matrix (origin)
  double Xjdet0;  // determinant of #Xjm0
  Epetra_SerialDenseVector* Fuv0o;  // deformation gradient at origin at t_{n}
  Epetra_SerialDenseVector* Fuv0;  // deformation gradient at origin at t_{n+1}
  Epetra_SerialDenseMatrix* boplin0; // B-operator (origin)
  Epetra_SerialDenseMatrix* W0o;  // W-operator (origin) at t_{n}
  Epetra_SerialDenseMatrix* W0;  // W-operator (origin) at t_{n+1}
  Epetra_SerialDenseMatrix* W0m;  // mid-W-operator (origin)
  Epetra_SerialDenseMatrix* Go;  // G-operator at t_{n}
  Epetra_SerialDenseMatrix* G;  // G-operator at t_{n+1}
  Epetra_SerialDenseMatrix* Gm;  // mid-G-operator
  Epetra_SerialDenseMatrix* Z;  // Z-operator
  Epetra_SerialDenseMatrix* FmCF;  // FCF^T
  Epetra_SerialDenseMatrix* Kda;  // EAS matrix Kda
  Epetra_SerialDenseMatrix* Kad;  // EAS matrix Kad
  Epetra_SerialDenseMatrix* Kaa;  // EAS matrix Kaa
  Epetra_SerialDenseVector* feas; // EAS portion of internal forces
  Epetra_SerialDenseMatrix* oldfeas;  // EAS history
  Epetra_SerialDenseMatrix* oldKaainv;  // EAS history
  Epetra_SerialDenseMatrix* oldKda;  // EAS history
  Epetra_SerialDenseMatrix* oldKad;  // EAS history

  // ------------------------------------ check calculation of mass matrix
  double density = (massmatrix) ? Density(material) : 0.0;

  // element co-ordinates
  for (int k=0; k<numnode; ++k)
  {
    Xe(0,k) = Nodes()[k]->X()[0];
    Xe(1,k) = Nodes()[k]->X()[1];
    xeo(0,k) = Xe(0,k) + dispo[k*Wall1::noddof_+0];
    xeo(1,k) = Xe(1,k) + dispo[k*Wall1::noddof_+1];
    xe(0,k) = Xe(0,k) + disp[k*Wall1::noddof_+0];
    xe(1,k) = Xe(1,k) + disp[k*Wall1::noddof_+1];
  }

  // set-up EAS parameters
  if (iseas_)
  {
    // allocate EAS quantities
    Fenhvo = new Epetra_SerialDenseMatrix(4,1);
    Fenhv = new Epetra_SerialDenseMatrix(4,1);
    Fmo = new Epetra_SerialDenseMatrix(4,3);
    Fm = new Epetra_SerialDenseMatrix(4,3);
    Fmm = Fmo;  // convenience
    Pvmm = new Epetra_SerialDenseMatrix(4,1);
    Xjm0 = new Epetra_SerialDenseMatrix(2,2);
    Fuv0o = new Epetra_SerialDenseVector(4);
    Fuv0 = new Epetra_SerialDenseVector(4);
    boplin0 = new Epetra_SerialDenseMatrix(4,edof);
    W0o = new Epetra_SerialDenseMatrix(4,edof);
    W0 = new Epetra_SerialDenseMatrix(4,edof);
    W0m = W0o;  // convenience pointer
    Go = new Epetra_SerialDenseMatrix(4,Wall1::neas_); 
    G = new Epetra_SerialDenseMatrix(4,Wall1::neas_); 
    Gm = Go;  // convenience pointer
    Z = new Epetra_SerialDenseMatrix(edof,Wall1::neas_);
    FmCF = new Epetra_SerialDenseMatrix(4,4);
    Kda = new Epetra_SerialDenseMatrix(edof,Wall1::neas_);
    Kad = new Epetra_SerialDenseMatrix(Wall1::neas_,edof);
    Kaa = new Epetra_SerialDenseMatrix(Wall1::neas_,Wall1::neas_);
    feas = new Epetra_SerialDenseVector(Wall1::neas_);

    // EAS Update of alphas:
    // the current alphas are (re-)evaluated out of
    // Kaa and Kda of previous step to avoid additional element call.
    // This corresponds to the (innermost) element update loop
    // in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");   // get alpha of last converged state
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get alpha of previous iteration

    // get stored EAS history
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    oldKad = data_.GetMutable<Epetra_SerialDenseMatrix>("Kad");
    if ( (not alpha) or (not oldKaainv) or (not oldKda) or (not oldKad) or (not oldfeas) )
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
    // calculate linear B-operator at origin
    w1_boplin(*boplin0, shpdrv, *Xjm0, Xjdet0, numnode);
    // displ.-based def.grad. at origin
    w1_defgrad(*Fuv0o, Ev, Xe, xeo, *boplin0, numnode);  // at t_{n}
    w1_defgrad(*Fuv0, Ev, Xe, xe, *boplin0, numnode);  // at t_{n+1}

    // evaluation of EAS variables (which are constant for the following):
    // -> M defining interpoxlation of enhanced strains alpha, evaluated at GPs
    // -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    // -> T0^{-T}
    //w1_eassetup(*boplin0, *Fuv0, *Xjm0, Xjdet0, Xe, xe, distype);
  }


  // integration loops over element domain
  for (int ip=0; ip<intpoints.nquad; ++ip)
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
      for (int a=0; a<numnode; a++)
      {
        for (int b=0; b<numnode; b++)
        {
          (*massmatrix)(2*a,2*b) += facm * shpfct(a) * shpfct(b); /* a,b even */
          (*massmatrix)(2*a+1,2*b+1) += facm * shpfct(a) * shpfct(b); /* a,b odd  */
        }
      }
    }

    // calculate linear B-operator
    w1_boplin(boplin, shpdrv, Xjm, Xjdet, numnode);

    // calculate defgrad F^u, Green-Lagrange-strain E^u
    w1_defgrad(Fuvo, Evo, Xe, xeo, boplin, numnode);  // at t_{n+1}
    w1_defgrad(Fuv, Ev, Xe, xe, boplin, numnode);  // at t_{n+1}

    // calculate non-linear B-operator in current configuration
    w1_boplin_cure(bopo, boplin, Fuvo, Wall1::numstr_, edof);  // at t_{n} // CHECK THIS: NOT SURE IF bopo NEEDED
    w1_boplin_cure(bop, boplin, Fuv, Wall1::numstr_, edof);  // at t_{n+1}

    // EAS: The deformation gradient is enhanced
    if (iseas_)
    {
      // calculate the enhanced deformation gradient and
      // also the operators G, W0 and Z
      w1_call_defgrad_enh(*Fenhvo, *Xjm0, Xjm, Xjdet0, Xjdet, *Fuv0, *alphao, xi1, xi2, *Go, *W0o, *boplin0, *Z); // at t_{n}
      w1_call_defgrad_enh(*Fenhv, *Xjm0, Xjm, Xjdet0, Xjdet, *Fuv0, *alpha, xi1, xi2, *G, *W0, *boplin0, *Z); // at t_{n+1}

      // total deformation gradient F, and total Green-Lagrange-strain E
      w1_call_defgrad_tot(*Fenhvo, *Fmo, Fuvo, Evo);  // at t_{n}
      w1_call_defgrad_tot(*Fenhv, *Fm, Fuv, Ev);  // at t_{n+1}
    }

    // mid-def.grad.
    // F_m = (1.0-gemmalphaf)*F_{n+1} + gemmalphaf*F_{n}
    {
      const int totdim = (*Fmm).M() * (*Fmm).N();
      // remember same pointer: F_m = F_{n}
      blas.SCAL(totdim, gemmalphaf, (*Fmm).A());  // F_m *= gemmalphaf = gemmalphaf*F_{n}
      blas.AXPY(totdim, (1.0-gemmalphaf), (*Fm).A(), (*Fmm).A());  // F_m += (1.0-gemmalphaf)*F_{n+1}
    }

    // non-linear mid-B-operator
    // B_m = (1.0-gemmalphaf)*B_{n+1} + gemmalphaf*B_{n}
    {
      const int totdim = bopm.M() * bopm.N();
      // remember same pointer B_m = B_{n}
      blas.SCAL(totdim, gemmalphaf, bopm.A());  // B_m *= gemmalphaf = gemmalphaf*B_{n}
      blas.AXPY(totdim, (1.0-gemmalphaf), bop.A(), bopm.A());  // B_m += (1.0-gemmalphaf)*B_{n+1}
    }

    // mid-strain GL vector
    // E_m = (1.0-gemmalphaf+gemmxi)*E_{n+1} + (gemmalphaf-gemmxi)*E_n
    {
      const int totdim = Evm.M() * Evm.N();
      // remember same pointer: E_m = E_{n}
      blas.SCAL(totdim, (gemmalphaf-gemmxi), Evm.A());  // E_m *= (gemmalphaf-gemmxi) = (gemmalphaf-gemmxi)*E_n
      blas.AXPY(totdim, (1.0-gemmalphaf+gemmxi), Ev.A(), Evm.A());  // E_m += (1.0-gemmalphaf+gemmxi)*E_{n+1}
    }

    // extra mid-quantities for case of EAS
    if (iseas_)
    {
      // mid-G-operator
      // G_m = 0.5*G_{n+1} + 0.5*G_{n}
      {
        const int totdim = (*Gm).M() * (*Gm).N();
        // remember same pointer: G_m = G_{n}
        blas.SCAL(totdim, 0.5, (*Gm).A());  // G_m *= 0.5 = 0.5*G_{n}
        blas.AXPY(totdim, 0.5, (*G).A(), (*Gm).A());  // G_m += 0.5*G_{n+1}
      }

      // mid-W0-operator
      // W_{0;m} = 0.5*W_{0;n+1} + 0.5*W_{0;n}
      {
        const int totdim = (*W0m).M() * (*W0m).N();
        // remember same pointer: W0_m = W0_{n}
        blas.SCAL(totdim, 0.5, (*W0m).A());   // W0_m *= 0.5 = 0.5*W0_{n}
        blas.AXPY(totdim, 0.5, (*W0).A(), (*W0m).A());  // W0_m += 0.5*W0_{n+1}
      }
    }

    // call material law
    if (material->mattyp == m_stvenant)
      w1_call_matgeononl(Evm, Smm, C, Wall1::numstr_, material);
    else
      dserror("It must be St.Venant-Kirchhoff material.");

    // return Gauss point strains (only in case of stress/strain output)
    if (elestrain)
    {
      for (int i = 0; i < Wall1::numstr_; ++i)
        (*elestrain)(ip,i) = Ev(i);
    }

    // return stresses at Gauss points (only in case of stress/strain output)
    if (elestress)
    {
      if (cauchy)
      {
        if (iseas_)
          StressCauchy(ip, 
                       (*Fm)(0,0), (*Fm)(1,1), (*Fm)(1,1), (*Fm)(1,2), 
                       Smm, elestress);
        else
          StressCauchy(ip, Fuv[0], Fuv[1], Fuv[2], Fuv[3], Smm, elestress);
      }
      else
      {
        (*elestress)(ip,0) = Smm(0,0);  // 2nd Piola-Kirchhoff stress S_{11}
        (*elestress)(ip,1) = Smm(1,1);  // 2nd Piola-Kirchhoff stress S_{22}
        (*elestress)(ip,2) = Smm(0,2);  // 2nd Piola-Kirchhoff stress S_{12}
      }
    }

    // stiffness and internal force
    if (iseas_)
    {
      // first mid-mid Piola-Kirchhoff stress vector P_{mm} = F_m . S_m
      w1_stress_eas(Smm, (*Fmm), (*Pvmm));

      // stiffness matrix kdd
      if (stiffmatrix)
        TangFintByDispGEMM(gemmalphaf, gemmxi, fac, 
                           boplin, (*W0m), (*W0),
                           (*Fmm), (*Fm), C, Smm,
                           (*FmCF), *stiffmatrix);
      // matrix kda
      TangFintByEnhGEMM(gemmalphaf, gemmxi, fac,
                        boplin, (*W0m), 
                        (*FmCF), Smm, (*G),
                        (*Z), (*Pvmm), 
                        (*Kda));
      // matrix kad
      TangEconByDispGEMM(gemmalphaf, gemmxi, fac,
                         boplin, (*W0), 
                         (*FmCF), Smm, (*Gm),
                         (*Z), (*Pvmm), 
                         (*Kad));
      // matrix kaa
      TangEconByEnhGEMM(gemmalphaf, gemmxi, fac, (*FmCF), Smm, (*G), (*Gm), (*Kaa));
      // nodal forces
      if (force) w1_fint_eas((*W0m), boplin, (*Gm), (*Pvmm), *force, (*feas), fac);
    }
    else
    {
      // geometric part of stiffness matrix kg
      if (stiffmatrix)
      {
        const double gemmfac = fac * (1.0-gemmalphaf);
        w1_kg(*stiffmatrix, boplin, Smm, gemmfac, edof, Wall1::numstr_);
      }
      // elastic+displacement stiffness matrix keu
      if (stiffmatrix)
      {
        const double gemmfac = fac * (1.0-gemmalphaf+gemmxi);
        w1_keu(*stiffmatrix, bopm, C, gemmfac, edof, Wall1::numstr_);
      }
      // nodal forces fi from integration of stresses
      if (force) w1_fint(Smm, bop, *force, fac, edof);
    }

  } // for (int ip=0; ip<totngp; ++ip)


  // EAS technology: static condensation
  // subtract EAS matrices from disp-based Kdd to "soften" element
  if ( (iseas_) and (force) and (stiffmatrix) )
  {
    // we need the inverse of Kaa
    Epetra_SerialDenseSolver solve_for_inverseKaa;
    solve_for_inverseKaa.SetMatrix((*Kaa));
    solve_for_inverseKaa.Invert();

    Epetra_SerialDenseMatrix KdaKaa(edof,Wall1::neas_); // temporary Kda.Kaa^{-1}
    KdaKaa.Multiply('N', 'N', 1.0, (*Kda), (*Kaa), 1.0);

    // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
    if (stiffmatrix) (*stiffmatrix).Multiply('N', 'N', -1.0, KdaKaa, (*Kad), 1.0);

    // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
    if (force) (*force).Multiply('N', 'N', -1.0, KdaKaa, (*feas), 1.0);

    // store current EAS data in history
    for (int i=0; i<Wall1::neas_; ++i)
      for (int j=0; j<Wall1::neas_; ++j)
        (*oldKaainv)(i,j) = (*Kaa)(i,j);

    for (int i=0; i<edof; ++i)
      for (int j=0; j<Wall1::neas_; ++j)
      {
        (*oldKda)(i,j) = (*Kda)(i,j);
        (*oldfeas)(j,0) = (*feas)(j);
      }
  }

  // clean EAS data
  if (iseas_)
  {
    delete Fenhvo;
    delete Fenhv;
    delete Fmo;
    delete Fm;
    delete Pvmm;
    delete Xjm0;
    delete Fuv0o;
    delete Fuv0;
    delete boplin0;
    delete W0o;
    delete W0;
    delete Go;
    delete G;
    delete Z;
    delete FmCF;
    delete Kda;
    delete Kad;
    delete Kaa;
    delete feas;
  }

  // good Bye
  return;
}

/*======================================================================*/
/* calcuate tangent (f_{int;m}),d */
void DRT::ELEMENTS::Wall1::TangFintByDispGEMM(
  const double& alphafgemm,
  const double& xigemm,
  const double& fac,
  const Epetra_SerialDenseMatrix& boplin,
  const Epetra_SerialDenseMatrix& W0m,
  const Epetra_SerialDenseMatrix& W0,
  const Epetra_SerialDenseMatrix& Fmm,
  const Epetra_SerialDenseMatrix& Fm,
  const Epetra_SerialDenseMatrix& C,
  const Epetra_SerialDenseMatrix& Smm,
  Epetra_SerialDenseMatrix& FmCF,
  Epetra_SerialDenseMatrix& estif
)
{
  // BLAS dummy
  Epetra_BLAS::Epetra_BLAS blas;

  // contitutive matrix (3x3)
  Epetra_SerialDenseMatrix C_red(3,3);
  C_red(0,0) = C(0,0);  C_red(0,1) = C(0,1);  C_red(0,2) = C(0,2);
  C_red(1,0) = C(1,0);  C_red(1,1) = C(1,1);  C_red(1,2) = C(1,2);
  C_red(2,0) = C(2,0);  C_red(2,1) = C(2,1);  C_red(2,2) = C(2,2);

  // FdotC (4 x 3) : F_m . C 
  Epetra_SerialDenseMatrix FmC(4,3);
  FmC.Multiply('N', 'N', 1.0, Fmm, C_red, 0.0);

  // FmCF (4 x 4) : ( F_m . C ) . F_{n+1}^T
  FmCF.Multiply('N', 'T', 1.0, FmC, Fm, 0.0);

  // BplusW (4 x edof) :  B_L + W0_{n+1}
  Epetra_SerialDenseMatrix BplusW(4,2*NumNode());
  {
    const int totdim = BplusW.M() * BplusW.N();
    blas.AXPY(totdim, 1.0, boplin.A(), BplusW.A());  // += B_L
    blas.AXPY(totdim, 1.0, W0.A(), BplusW.A());  // += W0_{n+1}
  }

  // Temp1 (4 x 8) : (Fm . C . F_{n+1}^T) . (B_L + W0_{n+1})
  Epetra_SerialDenseMatrix temp1(4,2*NumNode());
  temp1.Multiply('N', 'N', 1.0, FmCF, BplusW, 0.0);

  // Temp3 (4 x 8) : S_m . (B_L + W0_{n+1})
  Epetra_SerialDenseMatrix temp3(4,2*NumNode());
  temp3.Multiply('N', 'N', 1.0, Smm, BplusW, 0.0);

  // BplusW (4 x 8) :  B_L + W0_{m}
  {
    //const int totdim = BplusW.M() * BplusW.N();
    //blas.SCAL(totdim, 0.0, BplusW.A());  // BplusW *= 0.0
    //blas.AXPY(totdim, 1.0, boplin.A(), BplusW.A());  // BplusW += B_L
    //blas.AXPY(totdim, 1.0, W0m.A(), BplusW.A());  // BplusW += W0_{m}
    BplusW.Scale(0.0);
    BplusW += boplin;
    BplusW += W0m;
  } 

  // k_{dd} (8 x 8) :
  // k_{dd} += fac * (B_L + W0_m)^T . (Fm . C . F_{n+1}^T) . (B_L + W0_m)
  estif.Multiply('T', 'N', (1.0-alphafgemm+xigemm)*fac, BplusW, temp1, 1.0);
  // k_{dd} += fac * (B_L+W0_m)^T . S_m . (B_L+W0_{n+1})
  estif.Multiply('T', 'N', (1.0-alphafgemm)*fac, BplusW, temp3, 1.0);

  // that's it
  return;
}  // DRT::ELEMENTS::Wall1::w1_kdd

/*======================================================================*/
/* calculate tangent (f_{int;m}),alpha */
void DRT::ELEMENTS::Wall1::TangFintByEnhGEMM(
  const double& alphafgemm,
  const double& xigemm,
  const double& fac,
  const Epetra_SerialDenseMatrix& boplin,
  const Epetra_SerialDenseMatrix& W0m,
  const Epetra_SerialDenseMatrix& FmCF,
  const Epetra_SerialDenseMatrix& Smm,
  const Epetra_SerialDenseMatrix& G,
  const Epetra_SerialDenseMatrix& Z,
  const Epetra_SerialDenseMatrix& Pvmm,
  Epetra_SerialDenseMatrix& kda
)
{
  // BLAS dummy
  Epetra_BLAS::Epetra_BLAS blas;

  // FmCFG (4 x 4) : (F_m . C . F_{n+1}^T) . G_{n+1}
  Epetra_SerialDenseMatrix FmCFG(4,Wall1::neas_);
  FmCFG.Multiply('N', 'N', 1.0, FmCF, G, 0.0);

  // SmG (4 x 4) : S_m . G_{n+1}
  Epetra_SerialDenseMatrix SmG(4,Wall1::neas_);
  SmG.Multiply('N', 'N', 1.0, Smm, G, 0.0);

  // BplusW (4 x 8) :  B_L + W0_{m}
  Epetra_SerialDenseMatrix BplusW(4,2*NumNode());
  BplusW += boplin;
  BplusW += W0m;

  // PZ (8 x 4) : \bar{\bar{P}}_{mm} . Z_{n+1}
  Epetra_SerialDenseMatrix PZ(2*NumNode(),Wall1::neas_);
  for (int i=0; i<NumNode(); i++)
  {
    for (int ieas=0; ieas<Wall1::neas_; ieas++)
    {
      PZ(i*2,ieas) = Pvmm(0,0)*Z(i*2,ieas) + Pvmm(2,0)*Z(i*2+1,ieas);
      PZ(i*2+1,ieas) = Pvmm(3,0)*Z(i*2,ieas) + Pvmm(1,0)*Z(i*2+1,ieas);
    }
  }

  // k_{da} (8 x 4) :
  // k_{da} += fac * (B_l+W0_m)^T . (F_m . C . F_{n+1}^T) . G_{n+1}
  kda.Multiply('T', 'N', (1.0-alphafgemm+xigemm)*fac, BplusW, FmCFG, 1.0);
  // k_{da} += fac * (B_l+W0_m)^T . S_m . G_{n+1}
  kda.Multiply('T', 'N', (1.0-alphafgemm)*fac, BplusW, SmG, 1.0);
  // k_{da} += fac * \bar{\bar{P}}_{mm} . Z_{n+1}
  blas.AXPY(kda.N()*kda.M(), 0.5*fac, PZ.A(), kda.A());

  // ciao
  return;
}


/*======================================================================*/
/* calculate tangent (s_m),d */
void DRT::ELEMENTS::Wall1::TangEconByDispGEMM(
  const double& alphafgemm,
  const double& xigemm,
  const double& fac,
  const Epetra_SerialDenseMatrix& boplin,
  const Epetra_SerialDenseMatrix& W0,
  const Epetra_SerialDenseMatrix& FmCF,
  const Epetra_SerialDenseMatrix& Smm,
  const Epetra_SerialDenseMatrix& Gm,
  const Epetra_SerialDenseMatrix& Z,
  const Epetra_SerialDenseMatrix& Pvmm,
  Epetra_SerialDenseMatrix& kad
)
{
  // BLAS dummy
  Epetra_BLAS::Epetra_BLAS blas;

  // BplusW (4 x 8) :  B_L + W0_{n+1}
  Epetra_SerialDenseMatrix BplusW(4,2*NumNode());
  BplusW += boplin;
  BplusW += W0;

  // FmCFBW (4 x 8) : (F_m . C . F_{n+1}^T) . (B_L + W0_{n+1})
  Epetra_SerialDenseMatrix FmCFBW(4,Wall1::neas_);
  FmCFBW.Multiply('N', 'N', 1.0, FmCF, BplusW, 0.0);

  // SmGBW (4 x 8) : S_m . G_{n+1} . (B_L + W0_{n+1})
  Epetra_SerialDenseMatrix SmGBW(4,Wall1::neas_);
  SmGBW.Multiply('N', 'N', 1.0, Smm, BplusW, 0.0);

  // ZP (4 x 8) : (\bar{\bar{P}}_{mm} . Z_{n+1})^T
  Epetra_SerialDenseMatrix ZP(Wall1::neas_,2*NumNode());
  for (int i=0; i<NumNode(); i++)
  {
    for (int ieas=0; ieas<Wall1::neas_; ieas++)
    {
      ZP(ieas,i*2) = Pvmm(0,0)*Z(i*2,ieas) + Pvmm(2,0)*Z(i*2+1,ieas);
      ZP(ieas,i*2+1) = Pvmm(3,0)*Z(i*2,ieas) + Pvmm(1,0)*Z(i*2+1,ieas);
    }
  }

  // k_{ad} (4 x 8) :
  // k_{ad} += fac * G_{m}^T . (F_m . C . F_{n+1}^T) . (B_lin+W0_{n+1})
  kad.Multiply('T', 'N', (1.0-alphafgemm+xigemm)*fac, Gm, FmCFBW, 1.0);
  // k_{ad} += fac *  G_{m}^T . S_m . (B_l+W0_{n+1})^T
  kad.Multiply('T', 'N', (1.0-alphafgemm)*fac, Gm, SmGBW, 1.0);
  // k_{ad} += fac * (\bar{\bar{P}}_{mm} . Z_{n+1})^T
  blas.AXPY(kad.N()*kad.M(), 0.5*fac, ZP.A(), kad.A());

  // ciao
  return;
}

/*======================================================================*/
/* calculate tangent (s_m),alpha */
void DRT::ELEMENTS::Wall1::TangEconByEnhGEMM(
  const double& alphafgemm,
  const double& xigemm,
  const double& fac,
  const Epetra_SerialDenseMatrix& FmCF,
  const Epetra_SerialDenseMatrix& Smm,
  const Epetra_SerialDenseMatrix& G,
  const Epetra_SerialDenseMatrix& Gm,
  Epetra_SerialDenseMatrix& kaa
)
{
  // FmCFG : (F_m . C . F_{n+1}^T) . G_{n+1}
  Epetra_SerialDenseMatrix FmCFG(4,Wall1::neas_);
  FmCFG.Multiply('N', 'N', 1.0, FmCF, G, 0.0);

  // SmG : S_m  . G_{n+1}
  Epetra_SerialDenseMatrix SmG(4,Wall1::neas_);
  SmG.Multiply('N', 'N', 1.0, Smm, G, 0.0);

  // k_{aa} (4 x 4) :
  // k_{aa} += fac * G_m^T . (F_m . C . F_{n+1}^T) . G_{n+1}
  kaa.Multiply('T', 'N', (1.0-alphafgemm+xigemm)*fac, Gm, FmCFG, 1.0);
  // k_{aa} += fac * G_m^T . S_m . G_{n+1}
  kaa.Multiply('T', 'N', (1.0-alphafgemm)*fac, Gm, SmG, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
#endif  // D_WALL1
#endif  // CCADISCRET
