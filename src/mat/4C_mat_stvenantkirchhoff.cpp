/*----------------------------------------------------------------------*/
/*! \file
\brief
St. Venant-Kirchhoff material

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_stvenantkirchhoff.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::StVenantKirchhoff::StVenantKirchhoff(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->Get<double>("YOUNG")),
      poissonratio_(matdata->Get<double>("NUE")),
      density_(matdata->Get<double>("DENS"))
{
  if (youngs_ <= 0.) FOUR_C_THROW("Young's modulus must be greater zero");
  if (poissonratio_ >= 0.5 || poissonratio_ < -1.)
    FOUR_C_THROW("Poisson's ratio must be in [-1;0.5)");
}

Teuchos::RCP<CORE::MAT::Material> MAT::PAR::StVenantKirchhoff::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StVenantKirchhoff(this));
}

MAT::StVenantKirchhoffType MAT::StVenantKirchhoffType::instance_;


CORE::COMM::ParObject* MAT::StVenantKirchhoffType::Create(const std::vector<char>& data)
{
  auto* stvenantk = new MAT::StVenantKirchhoff();
  stvenantk->Unpack(data);
  return stvenantk;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::StVenantKirchhoff::StVenantKirchhoff() : params_(nullptr) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::StVenantKirchhoff::StVenantKirchhoff(MAT::PAR::StVenantKirchhoff* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::StVenantKirchhoff*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
// computes isotropic elasticity tensor in matrix notion for 3d
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::SetupCmat(CORE::LINALG::Matrix<6, 6>& cmat)
{
  // get material parameters
  const double Emod = params_->youngs_;      // Young's modulus (modulus of elasticity)
  const double nu = params_->poissonratio_;  // Poisson's ratio (Querdehnzahl)

  FillCmat(cmat, Emod, nu);
}

void MAT::StVenantKirchhoff::FillCmat(
    CORE::LINALG::Matrix<6, 6>& cmat, const double Emod, const double nu)
{
  // isotropic elasticity tensor C in Voigt matrix notation
  //                     [ 1-nu     nu     nu |       0       0       0    ]
  //                     [        1-nu     nu |       0       0       0    ]
  //         E           [               1-nu |       0       0       0    ]
  // C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~~~~  ~~~~~~ ]
  //     (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0       0    ]
  //                     [                    |         (1-2*nu)/2    0    ]
  //                     [ symmetric          |                 (1-2*nu)/2 ]
  //
  const double mfac = Emod / ((1.0 + nu) * (1.0 - 2.0 * nu));  // factor

  // clear the material tangent
  cmat.Clear();
  // write non-zero components
  cmat(0, 0) = mfac * (1.0 - nu);
  cmat(0, 1) = mfac * nu;
  cmat(0, 2) = mfac * nu;
  cmat(1, 0) = mfac * nu;
  cmat(1, 1) = mfac * (1.0 - nu);
  cmat(1, 2) = mfac * nu;
  cmat(2, 0) = mfac * nu;
  cmat(2, 1) = mfac * nu;
  cmat(2, 2) = mfac * (1.0 - nu);
  // ~~~
  cmat(3, 3) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(4, 4) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(5, 5) = mfac * 0.5 * (1.0 - 2.0 * nu);
}


/*----------------------------------------------------------------------*
//calculates stresses using one of the above method to evaluate the elasticity tensor
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Evaluate(const CORE::LINALG::SerialDenseVector* glstrain_e,
    CORE::LINALG::SerialDenseMatrix* cmat_e, CORE::LINALG::SerialDenseVector* stress_e)
{
  // this is temporary as long as the material does not have a
  // Matrix-type interface
  const CORE::LINALG::Matrix<6, 1> glstrain(glstrain_e->values(), true);
  CORE::LINALG::Matrix<6, 6> cmat(cmat_e->values(), true);
  CORE::LINALG::Matrix<6, 1> stress(stress_e->values(), true);

  SetupCmat(cmat);
  // evaluate stresses
  stress.MultiplyNN(cmat, glstrain);  // sigma = C . epsilon
}


/*----------------------------------------------------------------------*
//calculates stresses using one of the above method to evaluate the elasticity tensor
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  SetupCmat(*cmat);
  // evaluate stresses
  stress->MultiplyNN(*cmat, *glstrain);  // sigma = C . epsilon
}


/*----------------------------------------------------------------------*
 |  Calculate strain energy                                    gee 10/09|
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::StrainEnergy(
    const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, const int gp, const int eleGID)
{
  CORE::LINALG::Matrix<6, 6> cmat(true);
  SetupCmat(cmat);

  CORE::LINALG::Matrix<6, 1> stress(true);
  stress.MultiplyNN(cmat, glstrain);

  for (int k = 0; k < 6; ++k) psi += glstrain(k) * stress(k);
  psi /= 2.0;
}

FOUR_C_NAMESPACE_CLOSE
