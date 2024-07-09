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
Mat::PAR::StVenantKirchhoff::StVenantKirchhoff(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      youngs_(matdata.parameters.get<double>("YOUNG")),
      poissonratio_(matdata.parameters.get<double>("NUE")),
      density_(matdata.parameters.get<double>("DENS"))
{
  if (youngs_ <= 0.) FOUR_C_THROW("Young's modulus must be greater zero");
  if (poissonratio_ >= 0.5 || poissonratio_ < -1.)
    FOUR_C_THROW("Poisson's ratio must be in [-1;0.5)");
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::StVenantKirchhoff::create_material()
{
  return Teuchos::rcp(new Mat::StVenantKirchhoff(this));
}

Mat::StVenantKirchhoffType Mat::StVenantKirchhoffType::instance_;


Core::Communication::ParObject* Mat::StVenantKirchhoffType::create(const std::vector<char>& data)
{
  auto* stvenantk = new Mat::StVenantKirchhoff();
  stvenantk->unpack(data);
  return stvenantk;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::StVenantKirchhoff::StVenantKirchhoff() : params_(nullptr) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::StVenantKirchhoff::StVenantKirchhoff(Mat::PAR::StVenantKirchhoff* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::StVenantKirchhoff::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::StVenantKirchhoff::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
  {
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::StVenantKirchhoff*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
// computes isotropic elasticity tensor in matrix notion for 3d
 *----------------------------------------------------------------------*/
void Mat::StVenantKirchhoff::setup_cmat(Core::LinAlg::Matrix<6, 6>& cmat)
{
  // get material parameters
  const double Emod = params_->youngs_;      // Young's modulus (modulus of elasticity)
  const double nu = params_->poissonratio_;  // Poisson's ratio (Querdehnzahl)

  fill_cmat(cmat, Emod, nu);
}

void Mat::StVenantKirchhoff::fill_cmat(
    Core::LinAlg::Matrix<6, 6>& cmat, const double Emod, const double nu)
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
  cmat.clear();
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
void Mat::StVenantKirchhoff::evaluate(const Core::LinAlg::SerialDenseVector* glstrain_e,
    Core::LinAlg::SerialDenseMatrix* cmat_e, Core::LinAlg::SerialDenseVector* stress_e)
{
  // this is temporary as long as the material does not have a
  // Matrix-type interface
  const Core::LinAlg::Matrix<6, 1> glstrain(glstrain_e->values(), true);
  Core::LinAlg::Matrix<6, 6> cmat(cmat_e->values(), true);
  Core::LinAlg::Matrix<6, 1> stress(stress_e->values(), true);

  setup_cmat(cmat);
  // evaluate stresses
  stress.multiply_nn(cmat, glstrain);  // sigma = C . epsilon
}


/*----------------------------------------------------------------------*
//calculates stresses using one of the above method to evaluate the elasticity tensor
 *----------------------------------------------------------------------*/
void Mat::StVenantKirchhoff::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  setup_cmat(*cmat);
  // evaluate stresses
  stress->multiply_nn(*cmat, *glstrain);  // sigma = C . epsilon
}


/*----------------------------------------------------------------------*
 |  Calculate strain energy                                    gee 10/09|
 *----------------------------------------------------------------------*/
void Mat::StVenantKirchhoff::strain_energy(
    const Core::LinAlg::Matrix<6, 1>& glstrain, double& psi, const int gp, const int eleGID)
{
  Core::LinAlg::Matrix<6, 6> cmat(true);
  setup_cmat(cmat);

  Core::LinAlg::Matrix<6, 1> stress(true);
  stress.multiply_nn(cmat, glstrain);

  for (int k = 0; k < 6; ++k) psi += glstrain(k) * stress(k);
  psi /= 2.0;
}

FOUR_C_NAMESPACE_CLOSE
