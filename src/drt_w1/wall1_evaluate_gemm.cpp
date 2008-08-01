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

#include "wall1.H"


/*======================================================================*/
/* evaluate the element forces and stiffness and mass for GEMM */
void DRT::ELEMENTS::Wall1::FintStiffMassGEMM(
  const std::vector<int>& lm,
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
  const int numnode = NumNode();
  const int nd = numnode * Wall1::noddof_;
  const DiscretizationType distype = this->Shape();
  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints = getIntegrationPoints2D(gaussrule_);

  // general arrays
  Epetra_SerialDenseVector funct(numnode);
  Epetra_SerialDenseMatrix deriv(Wall1::numdim_,numnode);
  Epetra_SerialDenseMatrix xjm(Wall1::numdim_,Wall1::numdim_);
  Epetra_SerialDenseMatrix boplin(4,2*numnode);
  Epetra_SerialDenseVector F(4);
  Epetra_SerialDenseVector strain(4);
  double det;
  Epetra_SerialDenseMatrix xrefe(2,numnode);
  Epetra_SerialDenseMatrix xcure(2,numnode);
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure(numeps,nd);
  Epetra_SerialDenseMatrix stress(4,4);
  Epetra_SerialDenseMatrix C(4,4);

  // for eas, in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha;  // EAS alphas
  Epetra_SerialDenseMatrix* F_enh;  // EAS matrix F_enh
  Epetra_SerialDenseMatrix* F_tot;  // EAS vector F_tot
  Epetra_SerialDenseMatrix* p_stress;  // first piola-kirchhoff stress vector
  Epetra_SerialDenseMatrix* xjm0;  // Jacobian Matrix (origin)
  Epetra_SerialDenseVector* F0;  // Deformation Gradient (origin)
  Epetra_SerialDenseMatrix* boplin0; // B operator (origin)
  Epetra_SerialDenseMatrix* W0;  // W operator (origin)
  Epetra_SerialDenseMatrix* G;  // G operator
  Epetra_SerialDenseMatrix* Z;  // Z operator
  Epetra_SerialDenseMatrix* FCF;  // FCF^T
  Epetra_SerialDenseMatrix* Kda;  // EAS matrix Kda
  Epetra_SerialDenseMatrix* Kaa;  // EAS matrix Kaa
  Epetra_SerialDenseVector* feas; // EAS portion of internal forces
  double detJ0;  // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv; // EAS history
  Epetra_SerialDenseMatrix* oldKda;    // EAS history

  // ------------------------------------ check calculation of mass matrix
  double density = (massmatrix) ? Density(material) : 0.0;

  // element co-ordinates
  for (int k=0; k<numnode; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];
    xcure(0,k) = xrefe(0,k) + disp[k*Wall1::noddof_+0];
    xcure(1,k) = xrefe(1,k) + disp[k*Wall1::noddof_+1];
  }

  // set-up EAS parameters
  if (iseas_ == true)
  {
    // allocate EAS quantities
    F_enh = new Epetra_SerialDenseMatrix(4,1);
    F_tot = new Epetra_SerialDenseMatrix(4,3);
    p_stress = new Epetra_SerialDenseMatrix(4,1);
    xjm0 = new Epetra_SerialDenseMatrix(2,2);
    F0 = new Epetra_SerialDenseVector(4);
    boplin0 = new Epetra_SerialDenseMatrix(4,2*numnode);
    W0 = new Epetra_SerialDenseMatrix(4,2*numnode);
    G = new Epetra_SerialDenseMatrix(4,Wall1::neas_);
    Z = new Epetra_SerialDenseMatrix(2*numnode,Wall1::neas_);
    FCF = new Epetra_SerialDenseMatrix(4,4);
    Kda = new Epetra_SerialDenseMatrix(2*numnode,Wall1::neas_);
    Kaa = new Epetra_SerialDenseMatrix(Wall1::neas_,Wall1::neas_);
    feas = new Epetra_SerialDenseVector(Wall1::neas_);

    // EAS Update of alphas:
    // the current alphas are (re-)evaluated out of
    // Kaa and Kda of previous step to avoid additional element call.
    // This corresponds to the (innermost) element update loop
    // in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get alpha of previous iteration

    // get stored EAS history
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    if (!alpha || !oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history-data");

    // we need the (residual) displacement at the previous step
    Epetra_SerialDenseVector res_d(2*numnode);
    for (int i = 0; i < (2*numnode); ++i) {
      res_d(i) = residual[i];
    }

    // add Kda . res_d to feas
    (*oldfeas).Multiply('T','N',1.0,(*oldKda),res_d,1.0);
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);

    // evaluation of EAS variables (which are constant for the following):
    // -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    // -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    // -> T0^{-T}
    w1_eassetup((*boplin0), (*F0), (*xjm0), detJ0, xrefe, xcure, distype);
  }


  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    // Gaussian point and weight at it
    const double e1 = intpoints.qxg[ip][0];
    const double e2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

    // compute jacobian Matrix
    w1_jacobianmatrix(xrefe, deriv, xjm, &det, numnode);

    // integration factor
    double fac = wgt * det * thickness_;

    // compute mass matrix
    if (massmatrix)
    {
      double facm = fac * density;
      for (int a=0; a<numnode; a++)
      {
        for (int b=0; b<numnode; b++)
        {
          const double mab = facm * funct(a) * funct(b);
          (*massmatrix)(2*a,2*b) += mab; /* a,b even */
          (*massmatrix)(2*a+1,2*b+1) += mab; /* a,b odd  */
        }
      }
    }

    // calculate operator Blin
    w1_boplin(boplin, deriv, xjm, det, numnode);

    // calculate defgrad F^u, Green-Lagrange-strain E^u
    w1_defgrad(F, strain, xrefe, xcure, boplin, numnode);

    // calculate defgrad F in matrix notation and Blin in current conf.
    w1_boplin_cure(b_cure, boplin, F, numeps, nd);

    // EAS technology: "enhance the deformation gradient"
    if (iseas_)
    {
      // calculate the enhanced deformation gradient and
      // also the operators G, W0 and Z
      w1_call_defgrad_enh((*F_enh), (*xjm0), xjm, detJ0, det, (*F0), *alpha, e1, e2, (*G), (*W0), (*boplin0), (*Z));

      // total deformation gradient, Green-Lagrange-strain E^F
      w1_call_defgrad_tot((*F_enh), (*F_tot), F, strain);
    }

    // call material law
    w1_call_matgeononl(strain, stress, C, numeps, material);

    // return gp strains (only in case of stress/strain output)
    if (elestrain)
    {
      for (int i = 0; i < Wall1::numstr_; ++i)
        (*elestrain)(ip,i) = strain(i);
    }

    // return gp stresses (only in case of stress/strain output)
    if (elestress)
    {
      if (cauchy)
      {
        if (iseas_)
          StressCauchy(ip, (*F_tot)(0,0), (*F_tot)(1,1), (*F_tot)(1,1), (*F_tot)(1,2), stress, elestress);
        else
          StressCauchy(ip, F[0], F[1], F[2], F[3], stress, elestress);
      }
      else
      {
        (*elestress)(ip,0) = stress(0,0);
        (*elestress)(ip,1) = stress(1,1);
        (*elestress)(ip,2) = stress(0,2);
      }
    }

    // stiffness and internal force
    if (iseas_)
    {
      // first Piola-Kirchhoff stress vector
      w1_stress_eas(stress, (*F_tot), (*p_stress));

      // stiffness matrix kdd
      if (stiffmatrix) w1_kdd(boplin, (*W0), (*F_tot), C, stress, (*FCF), *stiffmatrix, fac);
      // matrix kda
      w1_kda((*FCF), (*W0), boplin, stress, (*G), (*Z), (*Kda), (*p_stress), fac);
      // matrix kaa
      w1_kaa((*FCF), stress, (*G), (*Kaa), fac);
      // nodal forces
      if (force) w1_fint_eas((*W0), boplin, (*G), (*p_stress), *force, (*feas), fac);
    }
    else
    {
      // geometric part of stiffness matrix kg
      if (stiffmatrix) w1_kg(*stiffmatrix, boplin, stress, fac, nd, numeps);
      // elastic+displacement stiffness matrix keu
      if (stiffmatrix) w1_keu(*stiffmatrix, b_cure, C, fac, nd, numeps);
      // nodal forces fi from integration of stresses
      if (force) w1_fint(stress, b_cure, *force, fac, nd);
    }

  } // for (int ip=0; ip<totngp; ++ip)


  // EAS technology: static condensation
  // subtract EAS matrices from disp-based Kdd to "soften" element
  if ( (force) and (stiffmatrix) )
  {
    if (iseas_)
    {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix((*Kaa));
      solve_for_inverseKaa.Invert();

      Epetra_SerialDenseMatrix KdaKaa(nd,Wall1::neas_); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('N', 'N', 1.0, (*Kda), (*Kaa), 1.0);


      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
      if (stiffmatrix) (*stiffmatrix).Multiply('N', 'T', -1.0, KdaKaa, (*Kda), 1.0);

      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      if (force) (*force).Multiply('N', 'N', -1.0, KdaKaa, (*feas), 1.0);

      // store current EAS data in history
      for (int i=0; i<Wall1::neas_; ++i)
        for (int j=0; j<Wall1::neas_; ++j)
          (*oldKaainv)(i,j) = (*Kaa)(i,j);

      for (int i=0; i<nd; ++i)
        for (int j=0; j<Wall1::neas_; ++j)
        {
          (*oldKda)(i,j) = (*Kda)(i,j);
          (*oldfeas)(j,0) = (*feas)(j);
        }
    }
  }

  // clean EAS data
  if (iseas_)
  {
    delete F_enh;
    delete F_tot;
    delete p_stress;
    delete xjm0;
    delete F0;
    delete boplin0;
    delete W0;
    delete G;
    delete Z;
    delete FCF;
    delete Kda;
    delete Kaa;
    delete feas;
  }

  // good Bye
  return;
}


/*----------------------------------------------------------------------*/
#endif  // D_WALL1
#endif  // CCADISCRET
