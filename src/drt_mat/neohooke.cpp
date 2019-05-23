/*!----------------------------------------------------------------------
\brief contains the functions to establish local material law
       stress-strain law for isotropic material for a 3D hex element
       following compressible Neo Hookean material law
       see Holzapfel, Nonlinear Solid Mechanics, pp. 243
       example input line:
       MAT 1 MAT_Struct_NeoHooke  YOUNG 100.0 NUE 0.3 DENS 1.0
\level 3

\maintainer Fabian Braeu
\param  Epetra_SerialDenseVector* glstrain      (i) Green-Lagrange strains
\param  Epetra_SerialDenseVector* stress        (o) ele stress vector
\param  Epetra_SerialDenseMatrix* cmat          (o) constitutive matrix
*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "neohooke.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::NeoHooke::NeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->GetDouble("YOUNG")),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::NeoHooke::CreateMaterial()
{
  return Teuchos::rcp(new MAT::NeoHooke(this));
}

MAT::NeoHookeType MAT::NeoHookeType::instance_;


DRT::ParObject* MAT::NeoHookeType::Create(const std::vector<char>& data)
{
  MAT::NeoHooke* neo = new MAT::NeoHooke();
  neo->Unpack(data);
  return neo;
}


/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
MAT::NeoHooke::NeoHooke() : params_(NULL)
{
  dserror(
      "This material law - NEOHOOKE - is maintained only inside the Elasthyper Toolbox.\n"
      "If you want to use this law, the material input line should read :\n"
      "MAT 1   MAT_ElastHyper   NUMMAT 1 MATIDS 2 DENS 0.1\n"
      "MAT 2   ELAST_CoupNeoHooke YOUNG 1.0 NUE 0.3\n");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::NeoHooke::NeoHooke(MAT::PAR::NeoHooke* params) : params_(params)
{
  dserror(
      "This material law - NEOHOOKE - is maintained only inside the Elasthyper Toolbox.\n"
      "If you want to use this law, the material input line should read :\n"
      "MAT 1   MAT_ElastHyper   NUMMAT 1 MATIDS 2 DENS 0.1\n"
      "MAT 2   ELAST_CoupNeoHooke YOUNG 1.0 NUE 0.3\n");
}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::NeoHooke::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::NeoHooke::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::NeoHooke*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 |  Calculate stress and constitutive tensor (Neo Hookean law) gee 10/09|
 *----------------------------------------------------------------------*/
void MAT::NeoHooke::Evaluate(
    const LINALG::Matrix<6, 1>& glstrain, LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress)
{
  // get material parameters
  const double ym = params_->youngs_;        // Young's modulus
  const double nu = params_->poissonratio_;  // Poisson's ratio

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6, 1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg(0) += 1.0;
  rcg(1) += 1.0;
  rcg(2) += 1.0;

  // 3rd invariant, determinant
  const double I3 = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                    0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                    0.25 * rcg(0) * rcg(4) * rcg(4);

  //----------------------------------------------------------------------
  // invert C
  LINALG::Matrix<6, 1> invc(false);

  double invdet = 1.0 / I3;

  invc(0) = rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4);
  invc(1) = rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5);
  invc(2) = rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3);
  invc(3) = 0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2);
  invc(4) = 0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4);
  invc(5) = 0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1);

  invc.Scale(invdet);

  // Material Constants c1 and beta
  const double c1 = 0.5 * ym / (2.0 * (1.0 + nu));
  const double beta = nu / (1.0 - 2.0 * nu);

  // energy function
  // Psi = c1/beta (I3^{-beta} - 1) + c1 ( I1-3 )

  // Second Piola-Kirchhoff stress tensor
  // S = -2 c1 I3^{-beta} C^{-1} + 2 c1 Identity
  const double fac = std::pow(I3, -beta);
  stress = invc;
  stress.Scale(-2.0 * c1 * fac);  // volumetric part
  const double iso = 2.0 * c1;    // isochoric part
  stress(0) += iso;
  stress(1) += iso;
  stress(2) += iso;

  // material tensor
  // C = 4 c1 beta I3^{-beta} C^{-1} dyad C^{-1} + 4 c1 I3^{-beta} C^{-1} boeppel C^{-1}
  // where `boeppel' is called `Holzapfelproduct' below
  const double delta6 = 4.0 * c1 * beta * fac;
  const double delta7 = 4.0 * c1 * fac;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) cmat(i, j) = delta6 * invc(i) * invc(j);
  AddtoCmatHolzapfelProduct(cmat, invc, delta7);

  return;
}  // end of neohooke evaluate

/*----------------------------------------------------------------------*
 |  Calculate stress and constitutive tensor (Neo Hookean law)  rm 08/07|
 *----------------------------------------------------------------------*/
void MAT::NeoHooke::Evaluate(const Epetra_SerialDenseVector* glstrain_e,
    Epetra_SerialDenseMatrix* cmat_e, Epetra_SerialDenseVector* stress_e)
{
#ifdef DEBUG
  // are these always supplied?
  if (!glstrain_e || !cmat_e || !stress_e) dserror("Data missing upon input in material neo hooke");
#endif

  // this is temporary as long as the material does not have a
  // Matrix-type interface
  const LINALG::Matrix<6, 1> glstrain(glstrain_e->A(), true);
  LINALG::Matrix<6, 6> cmat(cmat_e->A(), true);
  LINALG::Matrix<6, 1> stress(stress_e->A(), true);

  // get material parameters
  const double ym = params_->youngs_;        // Young's modulus
  const double nu = params_->poissonratio_;  // Poisson's ratio

  // Green-Lagrange Strain Tensor
  LINALG::Matrix<3, 3> E(false);
  E(0, 0) = glstrain(0);
  E(1, 1) = glstrain(1);
  E(2, 2) = glstrain(2);
  E(0, 1) = 0.5 * glstrain(3);
  E(1, 0) = 0.5 * glstrain(3);
  E(1, 2) = 0.5 * glstrain(4);
  E(2, 1) = 0.5 * glstrain(4);
  E(0, 2) = 0.5 * glstrain(5);
  E(2, 0) = 0.5 * glstrain(5);

  // Right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<3, 3> C(E);
  C.Scale(2.0);
  C(0, 0) += 1.0;
  C(1, 1) += 1.0;
  C(2, 2) += 1.0;

  // Principal Invariants I1 = tr(C) and I3 = det(C)
  // double I1 = C(0,0)+C(1,1)+C(2,2); // Needed only for energy
  const double I3 =
      C(0, 0) * C(1, 1) * C(2, 2) + C(0, 1) * C(1, 2) * C(2, 0) + C(0, 2) * C(1, 0) * C(2, 1) -
      (C(0, 2) * C(1, 1) * C(2, 0) + C(0, 1) * C(1, 0) * C(2, 2) + C(0, 0) * C(1, 2) * C(2, 1));

  // Calculation of C^-1 (Cinv)
  // C no longer needed, so we destroy it by inplace inversion
  C.Invert();

  // Material Constants c1 and beta
  const double c1 = 0.5 * ym / (2 * (1 + nu));
  const double beta = nu / (1 - 2 * nu);

  // Energy
  // double W = c1/beta * (pow(J,-beta) - 1) + c1 (I1-3);

  // PK2 Stresses
  LINALG::Matrix<3, 3> PK2(false);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      PK2(i, j) = 2.0 * c1 * (-pow(I3, -beta) * C(i, j));
    }
  PK2(0, 0) += (2.0 * c1);
  PK2(1, 1) += (2.0 * c1);
  PK2(2, 2) += (2.0 * c1);


  // Transfer PK2 tensor to stress vector
  stress(0) = PK2(0, 0);
  stress(1) = PK2(1, 1);
  stress(2) = PK2(2, 2);
  stress(3) = PK2(0, 1);
  stress(4) = PK2(1, 2);
  stress(5) = PK2(0, 2);

  // Elasticity Tensor
  const double delta6 = 4. * c1 * beta * pow(I3, -beta);
  const double delta7 = 4. * c1 * pow(I3, -beta);

  LINALG::Matrix<9, 9> ET(false);


  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
    {
      ET(k, l) =
          delta6 * (C(0, 0) * C(k, l)) + delta7 * 0.5 * (C(0, k) * C(0, l) + C(0, l) * C(0, k));
      ET(k + 3, l) =
          delta6 * (C(1, 0) * C(k, l)) + delta7 * 0.5 * (C(1, k) * C(0, l) + C(1, l) * C(0, k));
      ET(k + 3, l + 3) =
          delta6 * (C(1, 1) * C(k, l)) + delta7 * 0.5 * (C(1, k) * C(1, l) + C(1, l) * C(1, k));
      ET(k + 6, l) =
          delta6 * (C(2, 0) * C(k, l)) + delta7 * 0.5 * (C(2, k) * C(0, l) + C(2, l) * C(0, k));
      ET(k + 6, l + 3) =
          delta6 * (C(2, 1) * C(k, l)) + delta7 * 0.5 * (C(2, k) * C(1, l) + C(2, l) * C(1, k));
      ET(k + 6, l + 6) =
          delta6 * (C(2, 2) * C(k, l)) + delta7 * 0.5 * (C(2, k) * C(2, l) + C(2, l) * C(2, k));
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

  return;
}  // end of neohooke evaluate


/*----------------------------------------------------------------------*
 |  Calculate the inverse of a 2nd order tensor                 rm 08/07|
 *----------------------------------------------------------------------*/
void MAT::NeoHooke::InverseTensor(
    const Epetra_SerialDenseMatrix& M, Epetra_SerialDenseMatrix& Minv, const double I3)
{
  if (I3 == 0.0)
  {
    dserror("Right Cauchy Green not invertable in Neo Hooke material law");
  }
  else
  {
    Minv(0, 0) = 1 / I3 * (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2));
    Minv(1, 0) = -1 / I3 * (M(0, 1) * M(2, 2) - M(2, 1) * M(0, 2));
    Minv(2, 0) = 1 / I3 * (M(0, 1) * M(1, 2) - M(1, 1) * M(0, 2));
    Minv(0, 1) = -1 / I3 * (M(1, 0) * M(2, 2) - M(2, 0) * M(1, 2));
    Minv(1, 1) = 1 / I3 * (M(0, 0) * M(2, 2) - M(2, 0) * M(0, 2));
    Minv(2, 1) = -1 / I3 * (M(0, 0) * M(1, 2) - M(1, 0) * M(0, 2));
    Minv(0, 2) = 1 / I3 * (M(1, 0) * M(2, 1) - M(2, 0) * M(1, 1));
    Minv(1, 2) = -1 / I3 * (M(0, 0) * M(2, 1) - M(2, 0) * M(0, 1));
    Minv(2, 2) = 1 / I3 * (M(0, 0) * M(1, 1) - M(1, 0) * M(0, 1));
  }
  return;
}
