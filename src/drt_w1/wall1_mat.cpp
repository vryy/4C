/*======================================================================*/
/*!
\file wall1_mat.cpp
\brief material laws

\level 1

<pre>
\maintainer Michael Hiermeier
            hiermeier@lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
// macros


/*----------------------------------------------------------------------*/
// headers
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/structporo.H"

#include "../drt_poroelast/poroelast_utils.H"

#include "wall1.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Constitutive matrix C and stresses (private)                mgit 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_call_matgeononl(
    const Epetra_SerialDenseVector& strain,      ///< Green-Lagrange strain vector
    Epetra_SerialDenseMatrix& stress,            ///< stress vector
    Epetra_SerialDenseMatrix& C,                 ///< elasticity matrix
    const int numeps,                            ///< number of strains
    Teuchos::RCP<const MAT::Material> material,  ///< the material data
    Teuchos::ParameterList& params               ///< element parameter list
)
{
  if (material->MaterialType() == INPAR::MAT::m_structporo or
      material->MaterialType() == INPAR::MAT::m_structpororeaction or
      material->MaterialType() == INPAR::MAT::m_structpororeactionECM)
  {
    Teuchos::RCP<const MAT::StructPoro> actmat =
        Teuchos::rcp_static_cast<const MAT::StructPoro>(material);
    // setup is done in so3_poro
    // actmat->Setup(NUMGPT_SOH8);
    material = actmat->GetMaterial();
  }

  /*--------------------------- call material law -> get tangent modulus--*/
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*----------------------- linear elastic ---*/
    {
      const MAT::StVenantKirchhoff* actmat =
          static_cast<const MAT::StVenantKirchhoff*>(material.get());
      double ym = actmat->Youngs();
      double pv = actmat->PoissonRatio();


      /*-------------- some comments, so that even fluid people are able to
         understand this quickly :-)
         the "strain" vector looks like:

             | EPS_xx |
             | EPS_yy |
             | EPS_xy |
             | EPS_yx | */

      switch (wtype_)
      {
          /*---------------------------------material-tangente-- plane stress ---*/
        case plane_stress:
        {
          const double e1 = ym / (1. - pv * pv);
          const double e2 = pv * e1;
          const double e3 = e1 * (1. - pv) / 2.;

          C(0, 0) = e1;
          C(0, 1) = e2;
          C(0, 2) = 0.;
          C(0, 3) = 0.;

          C(1, 0) = e2;
          C(1, 1) = e1;
          C(1, 2) = 0.;
          C(1, 3) = 0.;

          C(2, 0) = 0.;
          C(2, 1) = 0.;
          C(2, 2) = e3;
          C(2, 3) = e3;

          C(3, 0) = 0.;
          C(3, 1) = 0.;
          C(3, 2) = e3;
          C(3, 3) = e3;

          break;
        }

          /*----------- material-tangente - plane strain, rotational symmetry ---*/
        case plane_strain:
        {
          const double c1 = ym / (1.0 + pv);
          const double b1 = c1 * pv / (1.0 - 2.0 * pv);
          const double a1 = b1 + c1;

          C(0, 0) = a1;
          C(0, 1) = b1;
          C(0, 2) = 0.;
          C(0, 3) = 0.;

          C(1, 0) = b1;
          C(1, 1) = a1;
          C(1, 2) = 0.;
          C(1, 3) = 0.;

          C(2, 0) = 0.;
          C(2, 1) = 0.;
          C(2, 2) = c1 / 2.;
          C(2, 3) = c1 / 2.;

          C(3, 0) = 0.;
          C(3, 1) = 0.;
          C(3, 2) = c1 / 2;
          C(3, 3) = c1 / 2;

          break;
        }

        default:
        {
          dserror("Only plane strain and plane stress options exist for wtype_");
          break;
        }
      }  // switch(wtype_)

      /*-------------------------- evaluate 2.PK-stresses -------------------*/
      /*------------------ Summenschleife -> += (2.PK stored as vecor) ------*/

      Epetra_SerialDenseVector svector;
      svector.Size(3);

      for (int k = 0; k < 3; k++)
      {
        for (int i = 0; i < numeps; i++)
        {
          svector(k) += C(k, i) * strain(i);
        }
      }
      /*------------------ 2.PK stored as matrix -----------------------------*/
      stress(0, 0) = svector(0);
      stress(0, 2) = svector(2);
      stress(1, 1) = svector(1);
      stress(1, 3) = svector(2);
      stress(2, 0) = svector(2);
      stress(2, 2) = svector(1);
      stress(3, 1) = svector(2);
      stress(3, 3) = svector(0);

      break;
    }

    case INPAR::MAT::m_neohooke: /*----- neo-Hookean material (popp 07/08) ---*/
    {
      const MAT::NeoHooke* actmat = static_cast<const MAT::NeoHooke*>(material.get());
      // get material parameters
      const double ym = actmat->Youngs();        // Young's modulus
      const double nu = actmat->PoissonRatio();  // Poisson's ratio

      switch (wtype_)
      {
        case plane_stress:
        {
          // Green-Lagrange Strain Tensor
          LINALG::SerialDenseMatrix E(3, 3);
          E(0, 0) = strain[0];
          E(1, 1) = strain[1];
          E(2, 2) = 0.0;
          E(0, 1) = strain[3];
          E(1, 0) = strain[3];
          E(1, 2) = 0.0;
          E(2, 1) = 0.0;
          E(0, 2) = 0.0;
          E(2, 0) = 0.0;

          // Right Cauchy-Green Tensor  CG = 2 * E + I
          LINALG::SerialDenseMatrix CG(E);
          CG.Scale(2.0);
          CG(0, 0) += 1.0;
          CG(1, 1) += 1.0;

          // define material constants c1 and beta
          const double c1 = 0.5 * ym / (2 * (1 + nu));
          const double beta = nu / (1 - 2 * nu);

          // S(2,2)=0 -> condition for unknown CG(2,2)
          double temp = CG(0, 0) * CG(1, 1) - CG(0, 1) * CG(1, 0);
          CG(2, 2) = std::pow(temp, -beta / (1 + beta));

          // Principal Invariant I3 = det(CG)
          const double I3 = CG(0, 0) * CG(1, 1) * CG(2, 2) + CG(0, 1) * CG(1, 2) * CG(2, 0) +
                            CG(0, 2) * CG(1, 0) * CG(2, 1) -
                            (CG(0, 2) * CG(1, 1) * CG(2, 0) + CG(0, 1) * CG(1, 0) * CG(2, 2) +
                                CG(0, 0) * CG(1, 2) * CG(2, 1));

          // Calculation of CG^-1 (CGinv)
          LINALG::SerialDenseMatrix CGinv(CG);
          LINALG::SymmetricInverse(CGinv, 3);

          // PK2 Stresses
          LINALG::SerialDenseMatrix PK2(3, 3);
          int i, j;
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
              PK2(i, j) = 2.0 * c1 * (-pow(I3, -beta) * CGinv(i, j));
            }
          PK2(0, 0) += (2.0 * c1);
          PK2(1, 1) += (2.0 * c1);
          PK2(2, 2) += (2.0 * c1);

          // Transfer PK2 to our matrix form
          Epetra_SerialDenseVector svector(3);
          svector[0] = PK2(0, 0);
          svector[1] = PK2(1, 1);
          svector[2] = PK2(0, 1);

          stress(0, 0) = svector(0);
          stress(0, 2) = svector(2);
          stress(1, 1) = svector(1);
          stress(1, 3) = svector(2);
          stress(2, 0) = svector(2);
          stress(2, 2) = svector(1);
          stress(3, 1) = svector(2);
          stress(3, 3) = svector(0);

          dserror("Plane stress NeoHooke material not yet implemented!");

          break;
        }
        case plane_strain:
        {
          // Green-Lagrange Strain Tensor
          LINALG::SerialDenseMatrix E(3, 3);
          E(0, 0) = strain[0];
          E(1, 1) = strain[1];
          E(2, 2) = 0.0;
          E(0, 1) = strain[3];
          E(1, 0) = strain[3];
          E(1, 2) = 0.0;
          E(2, 1) = 0.0;
          E(0, 2) = 0.0;
          E(2, 0) = 0.0;

          // Right Cauchy-Green Tensor  CG = 2 * E + I
          LINALG::SerialDenseMatrix CG(E);
          CG.Scale(2.0);
          CG(0, 0) += 1.0;
          CG(1, 1) += 1.0;

          // E(2,2)=0 -> condition for unknown CG(2,2)
          CG(2, 2) = 1.0;

          // define material constants c1 and beta
          const double c1 = 0.5 * ym / (2 * (1 + nu));
          const double beta = nu / (1 - 2 * nu);

          // Principal Invariant I3 = det(CG)
          const double I3 = CG(0, 0) * CG(1, 1) * CG(2, 2) + CG(0, 1) * CG(1, 2) * CG(2, 0) +
                            CG(0, 2) * CG(1, 0) * CG(2, 1) -
                            (CG(0, 2) * CG(1, 1) * CG(2, 0) + CG(0, 1) * CG(1, 0) * CG(2, 2) +
                                CG(0, 0) * CG(1, 2) * CG(2, 1));

          // Calculation of CG^-1 (CGinv)
          LINALG::SerialDenseMatrix CGinv(CG);
          LINALG::SymmetricInverse(CGinv, 3);

          // PK2 Stresses
          LINALG::SerialDenseMatrix PK2(3, 3);
          int i, j;
          for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
              PK2(i, j) = 2.0 * c1 * (-pow(I3, -beta) * CGinv(i, j));
            }
          PK2(0, 0) += (2.0 * c1);
          PK2(1, 1) += (2.0 * c1);
          PK2(2, 2) += (2.0 * c1);

          // Transfer PK2 to our matrix form
          Epetra_SerialDenseVector svector(3);
          svector[0] = PK2(0, 0);
          svector[1] = PK2(1, 1);
          svector[2] = PK2(0, 1);

          stress(0, 0) = svector(0);
          stress(0, 2) = svector(2);
          stress(1, 1) = svector(1);
          stress(1, 3) = svector(2);
          stress(2, 0) = svector(2);
          stress(2, 2) = svector(1);
          stress(3, 1) = svector(2);
          stress(3, 3) = svector(0);

          // Elasticity Tensor
          const double delta6 = 4. * c1 * beta * pow(I3, -beta);
          const double delta7 = 4. * c1 * pow(I3, -beta);

          int k, l;
          LINALG::SerialDenseMatrix ET(9, 9);
          LINALG::SerialDenseMatrix cmat(6, 6);

          for (k = 0; k < 3; k++)
            for (l = 0; l < 3; l++)
            {
              ET(k, l) = delta6 * (CGinv(0, 0) * CGinv(k, l)) +
                         delta7 * 0.5 * (CGinv(0, k) * CGinv(0, l) + CGinv(0, l) * CGinv(0, k));
              ET(k + 3, l) = delta6 * (CGinv(1, 0) * CGinv(k, l)) +
                             delta7 * 0.5 * (CGinv(1, k) * CGinv(0, l) + CGinv(1, l) * CGinv(0, k));
              ET(k + 3, l + 3) =
                  delta6 * (CGinv(1, 1) * CGinv(k, l)) +
                  delta7 * 0.5 * (CGinv(1, k) * CGinv(1, l) + CGinv(1, l) * CGinv(1, k));
              ET(k + 6, l) = delta6 * (CGinv(2, 0) * CGinv(k, l)) +
                             delta7 * 0.5 * (CGinv(2, k) * CGinv(0, l) + CGinv(2, l) * CGinv(0, k));
              ET(k + 6, l + 3) =
                  delta6 * (CGinv(2, 1) * CGinv(k, l)) +
                  delta7 * 0.5 * (CGinv(2, k) * CGinv(1, l) + CGinv(2, l) * CGinv(1, k));
              ET(k + 6, l + 6) =
                  delta6 * (CGinv(2, 2) * CGinv(k, l)) +
                  delta7 * 0.5 * (CGinv(2, k) * CGinv(2, l) + CGinv(2, l) * CGinv(2, k));
            }

          cmat(0, 0) = ET(0, 0);
          cmat(0, 1) = ET(1, 1);
          cmat(0, 2) = ET(2, 2);
          cmat(0, 3) = ET(1, 0);
          cmat(0, 4) = ET(2, 1);
          cmat(0, 5) = ET(2, 0);

          cmat(1, 0) = ET(3, 3);
          cmat(1, 1) = ET(4, 4);
          cmat(1, 2) = ET(5, 5);
          cmat(1, 3) = ET(4, 3);
          cmat(1, 4) = ET(5, 4);
          cmat(1, 5) = ET(5, 3);

          cmat(2, 0) = ET(6, 6);
          cmat(2, 1) = ET(7, 7);
          cmat(2, 2) = ET(8, 8);
          cmat(2, 3) = ET(7, 6);
          cmat(2, 4) = ET(8, 7);
          cmat(2, 5) = ET(8, 6);

          cmat(3, 0) = ET(3, 0);
          cmat(3, 1) = ET(4, 1);
          cmat(3, 2) = ET(5, 2);
          cmat(3, 3) = ET(4, 0);
          cmat(3, 4) = ET(5, 1);
          cmat(3, 5) = ET(5, 0);

          cmat(4, 0) = ET(6, 3);
          cmat(4, 1) = ET(7, 4);
          cmat(4, 2) = ET(8, 5);
          cmat(4, 3) = ET(7, 3);
          cmat(4, 4) = ET(8, 4);
          cmat(4, 5) = ET(8, 3);

          cmat(5, 0) = ET(6, 0);
          cmat(5, 1) = ET(7, 1);
          cmat(5, 2) = ET(8, 2);
          cmat(5, 3) = ET(7, 0);
          cmat(5, 4) = ET(8, 1);
          cmat(5, 5) = ET(8, 0);

          // Transfer elasticity tensor to our matrix form
          C(0, 0) = cmat(0, 0);
          C(0, 1) = cmat(0, 1);
          C(0, 2) = cmat(0, 3);
          C(0, 3) = cmat(0, 3);

          C(1, 0) = cmat(1, 0);
          C(1, 1) = cmat(1, 1);
          C(1, 2) = cmat(1, 3);
          C(1, 3) = cmat(1, 3);

          C(2, 0) = cmat(3, 0);
          C(2, 1) = cmat(3, 1);
          C(2, 2) = cmat(3, 3);
          C(2, 3) = cmat(3, 3);

          C(3, 0) = cmat(3, 0);
          C(3, 1) = cmat(3, 1);
          C(3, 2) = cmat(3, 3);
          C(3, 3) = cmat(3, 3);

          break;
        }
        default:
        {
          dserror("Only plane strain and plane stress options exist for wtype_");
          break;
        }
      }


      // dserror("2D NeoHooke not yet implemented");
      break;
    }

    case INPAR::MAT::m_elasthyper:  // general hyperelastic matrial (bborn, 06/09)
      // case INPAR::MAT::m_stvenant:  // st.venant-kirchhoff-material
      {
        MaterialResponse3dPlane(stress, C, strain, params);
        break;
      }

    default:
    {
      dserror("Invalid type of material law for wall element");
      break;
    }
  }  // switch(material->MaterialType())

  return;
}  // DRT::ELEMENTS::Wall1::w1_call_matgeononl

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::MaterialResponse3dPlane(Epetra_SerialDenseMatrix& stress,
    Epetra_SerialDenseMatrix& C, const Epetra_SerialDenseVector& strain,
    Teuchos::ParameterList& params)
{
  // make 3d equivalent of Green-Lagrange strain
  LINALG::Matrix<6, 1> gl(false);
  GreenLagrangePlane3d(strain, gl);

  // call 3d stress response
  LINALG::Matrix<6, 1> pk2(true);   // must be zerofied!!!
  LINALG::Matrix<6, 6> cmat(true);  // must be zerofied!!!
  MaterialResponse3d(&pk2, &cmat, &gl, params);

  // dimension reduction type
  if (wtype_ == plane_strain)
  {
    ;  // go on
  }
  else if (wtype_ == plane_stress)
  {
    // strain vector above bears final values on: E_{11},E_{22},E_{12}
    // strain vector above bears initial guesses on: E_{33},E_{23},E_{31}

    // initial plane stress error
    double pserr = std::sqrt(pk2(2) * pk2(2) + pk2(4) * pk2(4) + pk2(5) * pk2(5));

    // make Newton-Raphson iteration to identify
    // E_{33},E_{23},E_{31} which satisfy S_{33}=S_{23}=S_{31}=0
    int i = 0;
    const double tol = EPS6;
    const int n = 10;
    // working arrays
    LINALG::Matrix<3, 3> crr(
        true);  // LHS // constitutive matrix of restraint compo
                // this matrix needs to be zeroed out for further usage
                // in case the following while loop is entirely skipped during runtime
    LINALG::Matrix<3, 1> rr(false);  // RHS // stress residual of restraint compo
    LINALG::Matrix<3, 1> ir(false);  // SOL  // restraint strain components
    // the Newton-Raphson loop
    while ((pserr > tol) and (i < n))
    {
      // build sub-system a.b=c to solve
      crr(0, 0) = cmat(2, 2);
      crr(0, 1) = cmat(2, 4);
      crr(0, 2) = cmat(2, 5);
      crr(1, 0) = cmat(4, 2);
      crr(1, 1) = cmat(4, 4);
      crr(1, 2) = cmat(4, 5);
      crr(2, 0) = cmat(5, 2);
      crr(2, 1) = cmat(5, 4);
      crr(2, 2) = cmat(5, 5);
      rr(0) = -pk2(2);
      rr(1) = -pk2(4);
      rr(2) = -pk2(5);
      // solution
      // an in-place inversion is used, 'coz the inverse is needed below
      crr.Invert();
      ir.Multiply(crr, rr);
      // update
      gl(2) += ir(0);
      gl(4) += 2.0 * ir(1);  // NOT SURE ABOUT 2.0, LACKED TESTING MATERIAL
      gl(5) += 2.0 * ir(2);  // NOT SURE ABOUT 2.0, LACKED TESTING MATERIAL

      // call for new 3d stress response
      pk2.Clear();   // must be blanked!!
      cmat.Clear();  // must be blanked!!
      MaterialResponse3d(&pk2, &cmat, &gl, params);

      // current plane stress error
      pserr = std::sqrt(pk2(2) * pk2(2) + pk2(4) * pk2(4) + pk2(5) * pk2(5));

      // increment loop index
      i += 1;
    }

    // check if convergence was reached
    if ((i >= n) and (pserr > tol))
    {
      dserror("Failed to identify plane stress solution");
    }
    else
    {
      // static condensation
      // The restraint strains E_{33},E_{23},E_{31} have been made
      // dependent on free strains E_{11},E_{22},E_{12}
      // --- with an implicit function.
      // Thus the effect of the linearisation with respect to the
      // dependent strains must be added onto the free strains.
      LINALG::Matrix<3, 3> cfr(false);
      cfr(0, 0) = cmat(0, 2);
      cfr(0, 1) = cmat(0, 4);
      cfr(0, 2) = cmat(0, 5);
      cfr(1, 0) = cmat(1, 2);
      cfr(1, 1) = cmat(1, 4);
      cfr(1, 2) = cmat(1, 5);
      cfr(2, 0) = cmat(3, 2);
      cfr(2, 1) = cmat(3, 4);
      cfr(2, 2) = cmat(3, 5);
      LINALG::Matrix<3, 3> crrrf(false);
      crrrf.MultiplyNT(crr, cfr);
      LINALG::Matrix<3, 3> cfrrrrf(false);
      cfrrrrf.MultiplyNN(cfr, crrrf);
      // update constitutive matrix of free components
      cmat(0, 0) -= cfrrrrf(0, 0);
      cmat(0, 1) -= cfrrrrf(0, 1);
      cmat(0, 3) -= cfrrrrf(0, 2);
      cmat(1, 0) -= cfrrrrf(1, 0);
      cmat(1, 1) -= cfrrrrf(1, 1);
      cmat(1, 3) -= cfrrrrf(1, 2);
      cmat(3, 0) -= cfrrrrf(2, 0);
      cmat(3, 1) -= cfrrrrf(2, 1);
      cmat(3, 3) -= cfrrrrf(2, 2);
    }
  }
  else
  {
    dserror("Dimension reduction type wtype_=%d is not available.", wtype_);
  }

  // transform 2nd Piola--Kirchhoff stress back to 2d stress matrix
  memset(stress.A(), 0, stress.M() * stress.N() * sizeof(double));     // zerofy
  stress(0, 0) = stress(3, 3) = pk2(0);                                // S_{11}
  stress(1, 1) = stress(2, 2) = pk2(1);                                // S_{22}
  stress(0, 2) = stress(1, 3) = stress(3, 1) = stress(2, 0) = pk2(3);  // S_{12}

  // transform elasticity matrix  back to 2d matrix
  C(0, 0) = cmat(0, 0);  // C_{1111}
  C(0, 1) = cmat(0, 1);  // C_{1122}
  C(0, 2) = cmat(0, 3);  // C_{1112}
  C(0, 3) = cmat(0, 3);  // C_{1112} = C_{1121}

  C(1, 0) = cmat(1, 0);  // C_{2211}
  C(1, 1) = cmat(1, 1);  // C_{2222}
  C(1, 2) = cmat(1, 3);  // C_{2212}
  C(1, 3) = cmat(1, 3);  // C_{2221} = C_{2212}

  C(2, 0) = cmat(3, 0);  // C_{1211}
  C(2, 1) = cmat(3, 1);  // C_{1222}
  C(2, 2) = cmat(3, 3);  // C_{1212}
  C(2, 3) = cmat(3, 3);  // C_{1221} = C_{1212}

  C(3, 0) = cmat(3, 0);  // C_{2111} = C_{1211}
  C(3, 1) = cmat(3, 1);  // C_{2122} = C_{1222}
  C(3, 2) = cmat(3, 3);  // C_{2112} = C_{1212}
  C(3, 3) = cmat(3, 3);  // C_{2121} = C_{1212}

  // leave this dump
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::MaterialResponse3d(LINALG::Matrix<6, 1>* stress,
    LINALG::Matrix<6, 6>* cmat, const LINALG::Matrix<6, 1>* glstrain,
    Teuchos::ParameterList& params)
{
  SolidMaterial()->Evaluate(NULL, glstrain, params, stress, cmat, Id());

  return;
}

/*-----------------------------------------------------------------------------*
| deliver internal/strain energy                                    bborn 08/08|
*-----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Wall1::EnergyInternal(Teuchos::RCP<const MAT::Material> material,
    Teuchos::ParameterList& params, const Epetra_SerialDenseVector& Ev)
{
  // switch material type
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // linear elastic
    {
      Epetra_SerialDenseMatrix Cm(Wall1::numnstr_, Wall1::numnstr_);  // elasticity matrix
      Epetra_SerialDenseMatrix Sm(Wall1::numnstr_, Wall1::numnstr_);  // 2nd PK stress matrix
      w1_call_matgeononl(Ev, Sm, Cm, Wall1::numnstr_, material, params);
      Epetra_SerialDenseVector Sv(Wall1::numnstr_);  // 2nd PK stress vector
      Sv(0) = Sm(0, 0);
      Sv(1) = Sm(1, 1);
      Sv(2) = Sv(3) = Sm(0, 2);
      return 0.5 * Sv.Dot(Ev);
    }
    break;
    case INPAR::MAT::m_neohooke:
    {
      dserror(
          "Not implemented for material type INPAR::MAT::m_neohooke. "
          "This material type should go away, anyway.");
    }
    break;
    case INPAR::MAT::m_elasthyper:
    {
      // transform the 2d Green-Lagrange strains into 3d notation
      LINALG::Matrix<6, 1> glstrain(true);
      GreenLagrangePlane3d(Ev, glstrain);

      // strain energy
      double psi = 0.0;

      // call material for evaluation of strain energy function
      SolidMaterial()->StrainEnergy(glstrain, psi, Id());

      return psi;
    }
    break;
    case INPAR::MAT::m_structporo:
    case INPAR::MAT::m_structpororeaction:
    case INPAR::MAT::m_structpororeactionECM:
    {
      // transform the 2d Green-Lagrange strains into 3d notation
      LINALG::Matrix<6, 1> glstrain(true);
      GreenLagrangePlane3d(Ev, glstrain);

      // strain energy
      double psi = 0.0;

      // call material for evaluation of strain energy function
      SolidMaterial()->StrainEnergy(glstrain, psi, Id());

      return psi;
    }
    break;
    default:
    {
      dserror("Illegal type of material for this element");
      return 0.0;
    }
    break;
  }  // end of switch (material->MaterialType())

  dserror(
      "You should never end up here, since all possible cases should be "
      "covered by the material selection.");

  return 0.0;
}

/*-----------------------------------------------------------------------------*
| deliver kinetic energy                                            bborn 08/08|
*-----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Wall1::EnergyKinetic(
    const Epetra_SerialDenseMatrix& mass, const std::vector<double>& vel)
{
  double kin = 0.0;
  for (int i = 0; i < 2 * NumNode(); ++i)
    for (int j = 0; j < 2 * NumNode(); ++j) kin += 0.5 * vel[i] * mass(i, j) * vel[j];
  return kin;
}

/*----------------------------------------------------------------------*/
