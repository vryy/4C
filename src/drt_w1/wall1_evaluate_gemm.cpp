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
  Epetra_SerialDenseVector Evm(4);  // Green-Lagrange mid-strain vector

  Epetra_SerialDenseMatrix Xe(Wall1::numdim_,numnode);  // material/initial element co-ordinates
  Epetra_SerialDenseMatrix xeo(Wall1::numdim_,numnode);  // spatial/current element co-ordinates at t_{n}
  Epetra_SerialDenseMatrix xe(Wall1::numdim_,numnode);  // spatial/current element co-ordinates at t_{n+1}
  Epetra_SerialDenseMatrix bopo(Wall1::numstr_,edof);  // non-linear B-op at t_{n}
  Epetra_SerialDenseMatrix bop(Wall1::numstr_,edof);  // non-linear B-op at t_{n+1}
  Epetra_SerialDenseMatrix bopm(Wall1::numstr_,edof);  // non-linear mid-B-op
  Epetra_SerialDenseMatrix pk2mm(4,4);  // 2nd Piola-Kirchhoff mid-stress matrix  // CHECK THIS: STRESS MATRIX SHOULD NOT EXIST IN EFFICIENT CODE
  Epetra_SerialDenseMatrix C(4,4);

  // for EAS, in any case declare variables, sizes etc. only allocated in EAS version
  Epetra_SerialDenseMatrix* alphao;  // EAS alphas at t_{n}
  Epetra_SerialDenseMatrix* alpha;  // EAS alphas at t_{n+1}
  Epetra_SerialDenseMatrix* Fenhvo;  // EAS matrix Fenhv
  Epetra_SerialDenseMatrix* Fenhv;  // EAS matrix Fenhv
  Epetra_SerialDenseMatrix* Fmo;  // total def.grad. matrix at t_{n}
  Epetra_SerialDenseMatrix* Fm;  // total def.grad. matrix at t_{n+1}
  Epetra_SerialDenseMatrix* Ftoto;  // EAS vector Ftot at t_{n}

  Epetra_SerialDenseMatrix* pk1sts;  // first piola-kirchhoff pk2mm vector
  Epetra_SerialDenseMatrix* Xjm0;  // Jacobian Matrix (origin)
  double Xjdet0;  // determinant of #Xjm0
  Epetra_SerialDenseVector* Fuv0o;  // deformation gradient at origin at t_{n}
  Epetra_SerialDenseVector* Fuv0;  // deformation gradient at origin at t_{n+1}
  Epetra_SerialDenseMatrix* boplin0; // B-operator (origin)
  Epetra_SerialDenseMatrix* W0;  // W-operator (origin)
  Epetra_SerialDenseMatrix* G;  // G-operator
  Epetra_SerialDenseMatrix* Z;  // Z-operator
  Epetra_SerialDenseMatrix* FCF;  // FCF^T
  Epetra_SerialDenseMatrix* Kda;  // EAS matrix Kda
  Epetra_SerialDenseMatrix* Kaa;  // EAS matrix Kaa
  Epetra_SerialDenseVector* feas; // EAS portion of internal forces
  Epetra_SerialDenseMatrix* oldfeas;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv; // EAS history
  Epetra_SerialDenseMatrix* oldKda;    // EAS history

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

    pk1sts = new Epetra_SerialDenseMatrix(4,1);
    Xjm0 = new Epetra_SerialDenseMatrix(2,2);
    Fuv0o = new Epetra_SerialDenseVector(4);
    Fuv0 = new Epetra_SerialDenseVector(4);
    boplin0 = new Epetra_SerialDenseMatrix(4,edof);
    W0 = new Epetra_SerialDenseMatrix(4,edof);
    G = new Epetra_SerialDenseMatrix(4,Wall1::neas_);
    Z = new Epetra_SerialDenseMatrix(edof,Wall1::neas_);
    FCF = new Epetra_SerialDenseMatrix(4,4);
    Kda = new Epetra_SerialDenseMatrix(edof,Wall1::neas_);
    Kaa = new Epetra_SerialDenseMatrix(Wall1::neas_,Wall1::neas_);
    feas = new Epetra_SerialDenseVector(Wall1::neas_);

    // Get quantities of last converged step
    Ftoto = new Epetra_SerialDenseMatrix(4,3);

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
    if ( (not alpha) or (not oldKaainv) or (not oldKda) or (not oldfeas) )
      dserror("Missing EAS history-data");

    // we need the (residual) displacement at the previous step
    Epetra_SerialDenseVector res_d(edof);
    for (int i = 0; i < edof; ++i) 
    {
      res_d(i) = residual[i];
    }

    // add Kda . res_d to feas
    (*oldfeas).Multiply('T','N',1.0,(*oldKda),res_d,1.0);
    // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);

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
      w1_call_defgrad_enh(*Fenhvo, *Xjm0, Xjm, Xjdet0, Xjdet, *Fuv0, *alpha, xi1, xi2, *G, *W0, *boplin0, *Z); // at t_{n}
      w1_call_defgrad_enh(*Fenhv, *Xjm0, Xjm, Xjdet0, Xjdet, *Fuv0, *alpha, xi1, xi2, *G, *W0, *boplin0, *Z); // at t_{n+1}

      // total deformation gradient F, and total Green-Lagrange-strain E
      w1_call_defgrad_tot(*Fenhvo, *Fmo, Fuvo, Evo);  // at t_{n}
      w1_call_defgrad_tot(*Fenhv, *Fm, Fuv, Ev);  // at t_{n+1}
    }

    // create algorithmic mid-strain ... after all it's GEMM
    if (iseas_)
    {
      dserror("Do something here");
    }
    else
    {
      // non-linear mid-B-operator
      // B_m = (1.0-gemmalphaf)*B_{n+1} + gemmalphaf*B_{n}
      {
        const int totdim = bopm.M() * bopm.N();
        blas.SCAL(totdim, 0.0, bopm.A());  // B_m *= 0.0;
        blas.AXPY(totdim, (1.0-gemmalphaf), bop.A(), bopm.A());  // B_m += (1.0-gemmalphaf)*B_{n+1}
        blas.AXPY(totdim, gemmalphaf, bopo.A(), bopm.A());  // B_m += gemmalphaf*B_{n}
      }

      // mid-strain GL vector
      // E_m = (1.0-gemmalphaf+gemmxi)*E_{n+1} + (gemmalphaf-gemmxi)*E_n
      {
        const int totdim = Evm.M() * Evm.N();
        blas.SCAL(totdim, 0.0, Evm.A());  // E_m *= 0.0
        blas.AXPY(totdim, (1.0-gemmalphaf+gemmxi), Ev.A(), Evm.A());  // E_m += (1.0-gemmalphaf+gemmxi)*E_{n+1}
        blas.AXPY(totdim, (gemmalphaf-gemmxi), Evo.A(), Evm.A());  // E_m += (gemmalphaf-gemmxi)*E_n
      }
    }
    

    // call material law
    if (material->mattyp == m_stvenant)
      w1_call_matgeononl(Evm, pk2mm, C, Wall1::numstr_, material);
    else
      dserror("It must be St.Venant-Kirchhoff material.");

    // return Gauss point strains (only in case of pk2mm/strain output)
    if (elestrain)
    {
      for (int i = 0; i < Wall1::numstr_; ++i)
        (*elestrain)(ip,i) = Ev(i);
    }

    // return gp stresses (only in case of pk2mm/strain output)
    if (elestress)
    {
      if (cauchy)
      {
        if (iseas_)
          StressCauchy(ip, 
                       (*Fm)(0,0), (*Fm)(1,1), (*Fm)(1,1), (*Fm)(1,2), 
                       pk2mm, elestress);
        else
          StressCauchy(ip, Fuv[0], Fuv[1], Fuv[2], Fuv[3], pk2mm, elestress);
      }
      else
      {
        (*elestress)(ip,0) = pk2mm(0,0);
        (*elestress)(ip,1) = pk2mm(1,1);
        (*elestress)(ip,2) = pk2mm(0,2);
      }
    }

    // stiffness and internal force
    if (iseas_)
    {
      // first Piola-Kirchhoff pk2mm vector
      w1_stress_eas(pk2mm, (*Fm), (*pk1sts));

      // stiffness matrix kdd
      if (stiffmatrix) w1_kdd(boplin, (*W0), (*Fm), C, pk2mm, (*FCF), *stiffmatrix, fac);
      // matrix kda
      w1_kda((*FCF), (*W0), boplin, pk2mm, (*G), (*Z), (*Kda), (*pk1sts), fac);
      // matrix kaa
      w1_kaa((*FCF), pk2mm, (*G), (*Kaa), fac);
      // nodal forces
      if (force) w1_fint_eas((*W0), boplin, (*G), (*pk1sts), *force, (*feas), fac);
    }
    else
    {
      // geometric part of stiffness matrix kg
      if (stiffmatrix) w1_kg(*stiffmatrix, boplin, pk2mm, fac, edof, Wall1::numstr_);
      // elastic+displacement stiffness matrix keu
      if (stiffmatrix) w1_keu(*stiffmatrix, bop, C, fac, edof, Wall1::numstr_);
      // nodal forces fi from integration of stresses
      if (force) w1_fint(pk2mm, bop, *force, fac, edof);
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

      Epetra_SerialDenseMatrix KdaKaa(edof,Wall1::neas_); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('N', 'N', 1.0, (*Kda), (*Kaa), 1.0);

      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
      if (stiffmatrix) (*stiffmatrix).Multiply('N', 'T', -1.0, KdaKaa, (*Kda), 1.0);

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
  }

  // clean EAS data
  if (iseas_)
  {
    delete Fenhvo;
    delete Fenhv;
    delete Fmo;
    delete Fm;
    delete pk1sts;
    delete Xjm0;
    delete Fuv0o;
    delete Fuv0;
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
