/*----------------------------------------------------------------------*/
/*! \file
\brief St.Venant Kirchhoff with an additional temperature dependent term
       describing heat expansion

       example input line:
       MAT 1   MAT_Struct_ThrStVenantK YOUNGNUM 2 YOUNG 1.48e8 1.48e5 NUE 0.3 DENS
       9.130e-6 THEXPANS 1.72e-5 INITTEMP 293.15

\level 2

*/
/*----------------------------------------------------------------------*
 | headers                                                   dano 02/10 |
 *----------------------------------------------------------------------*/
#include "4C_mat_thermostvenantkirchhoff.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |                                                           dano 02/10 |
 *----------------------------------------------------------------------*/
Mat::PAR::ThermoStVenantKirchhoff::ThermoStVenantKirchhoff(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_((matdata->Get<std::vector<double>>("YOUNG"))),
      poissonratio_(matdata->Get<double>("NUE")),
      density_(matdata->Get<double>("DENS")),
      thermexpans_(matdata->Get<double>("THEXPANS")),
      capa_(matdata->Get<double>("CAPA")),
      conduct_(matdata->Get<double>("CONDUCT")),
      thetainit_(matdata->Get<double>("INITTEMP")),
      thermomat_(matdata->Get<int>("THERMOMAT"))
{
  if (poissonratio_ >= 0.5 || poissonratio_ < -1.)
    FOUR_C_THROW("Poisson's ratio must be in [-1;0.5)");
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ThermoStVenantKirchhoff::create_material()
{
  return Teuchos::rcp(new Mat::ThermoStVenantKirchhoff(this));
}


Mat::ThermoStVenantKirchhoffType Mat::ThermoStVenantKirchhoffType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ThermoStVenantKirchhoffType::Create(
    const std::vector<char>& data)
{
  auto* thrstvenantk = new Mat::ThermoStVenantKirchhoff();
  thrstvenantk->Unpack(data);
  return thrstvenantk;
}


/*----------------------------------------------------------------------*
 |  constructor (public)                                     dano 02/10 |
 *----------------------------------------------------------------------*/
Mat::ThermoStVenantKirchhoff::ThermoStVenantKirchhoff() : params_(nullptr), thermo_(Teuchos::null)
{
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 02/10 |
 *----------------------------------------------------------------------*/
Mat::ThermoStVenantKirchhoff::ThermoStVenantKirchhoff(Mat::PAR::ThermoStVenantKirchhoff* params)
    : params_(params), thermo_(Teuchos::null)
{
  create_thermo_material_if_set();
}

void Mat::ThermoStVenantKirchhoff::create_thermo_material_if_set()
{
  const int thermoMatId = this->params_->thermomat_;
  if (thermoMatId != -1)
  {
    auto mat = Mat::Factory(thermoMatId);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to create thermo material, id=%d", thermoMatId);
    thermo_ = Teuchos::rcp_dynamic_cast<Mat::Trait::Thermo>(mat);
  }
}


/*----------------------------------------------------------------------*
 |  Pack (public)                                            dano 02/10 |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}  // Pack()


/*----------------------------------------------------------------------*
 |  Unpack (public)                                          dano 02/10 |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ThermoStVenantKirchhoff*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());

      create_thermo_material_if_set();
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}  // Unpack()


/*----------------------------------------------------------------------*
 | calculates stresses using one of the above method to      dano 02/10 |
 | evaluate the elasticity tensor                                       |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::Evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // fixme this backwards compatibility modification should be moved outside
  // use initial value as a default value
  double temperature = params.get<double>("scalartemp", params_->thetainit_);
  Reinit(defgrd, glstrain, temperature, gp);  // fixme call this before

  setup_cmat(*cmat);
  // purely mechanical part
  stress->MultiplyNN(*cmat, *current_glstrain_);

  // additive thermal part
  double Tref = params_->thetainit_;
  double m = st_modulus();

  // loop over the element nodes, non-zero entries only in main directions
  for (int i = 0; i < 3; ++i) (*stress)(i, 0) += m * (current_temperature_ - Tref);

}  // STR_Evaluate()

/*----------------------------------------------------------------------*
 | calculates strain energy                                 seitz 11/15 |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::StrainEnergy(
    const Core::LinAlg::Matrix<6, 1>& glstrain, double& psi, const int gp, const int eleGID)
{
  if (youngs_is_temp_dependent())
    FOUR_C_THROW("Calculation of strain energy only for constant Young's modulus");
  Core::LinAlg::Matrix<6, 6> cmat;
  setup_cmat(cmat);
  Core::LinAlg::Matrix<6, 1> s;
  s.Multiply(cmat, glstrain);
  psi += .5 * s.Dot(glstrain);
}

void Mat::ThermoStVenantKirchhoff::Evaluate(const Core::LinAlg::Matrix<3, 1>& gradtemp,
    Core::LinAlg::Matrix<3, 3>& cmat, Core::LinAlg::Matrix<3, 1>& heatflux) const
{
  thermo_->Evaluate(gradtemp, cmat, heatflux);
}

void Mat::ThermoStVenantKirchhoff::Evaluate(const Core::LinAlg::Matrix<2, 1>& gradtemp,
    Core::LinAlg::Matrix<2, 2>& cmat, Core::LinAlg::Matrix<2, 1>& heatflux) const
{
  thermo_->Evaluate(gradtemp, cmat, heatflux);
}

void Mat::ThermoStVenantKirchhoff::Evaluate(const Core::LinAlg::Matrix<1, 1>& gradtemp,
    Core::LinAlg::Matrix<1, 1>& cmat, Core::LinAlg::Matrix<1, 1>& heatflux) const
{
  thermo_->Evaluate(gradtemp, cmat, heatflux);
}

void Mat::ThermoStVenantKirchhoff::ConductivityDerivT(Core::LinAlg::Matrix<3, 3>& dCondDT) const
{
  thermo_->ConductivityDerivT(dCondDT);
}

void Mat::ThermoStVenantKirchhoff::ConductivityDerivT(Core::LinAlg::Matrix<2, 2>& dCondDT) const
{
  thermo_->ConductivityDerivT(dCondDT);
}

void Mat::ThermoStVenantKirchhoff::ConductivityDerivT(Core::LinAlg::Matrix<1, 1>& dCondDT) const
{
  thermo_->ConductivityDerivT(dCondDT);
}

double Mat::ThermoStVenantKirchhoff::CapacityDerivT() const { return thermo_->CapacityDerivT(); }

void Mat::ThermoStVenantKirchhoff::Reinit(double temperature, unsigned gp)
{
  current_temperature_ = temperature;
  if (thermo_ != Teuchos::null) thermo_->Reinit(temperature, gp);
}
void Mat::ThermoStVenantKirchhoff::ResetCurrentState()
{
  if (thermo_ != Teuchos::null) thermo_->ResetCurrentState();
}

void Mat::ThermoStVenantKirchhoff::CommitCurrentState()
{
  if (thermo_ != Teuchos::null) thermo_->CommitCurrentState();
}

void Mat::ThermoStVenantKirchhoff::Reinit(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, double temperature, unsigned gp)
{
  current_glstrain_ = glstrain;
  Reinit(temperature, gp);
}

void Mat::ThermoStVenantKirchhoff::GetdSdT(Core::LinAlg::Matrix<6, 1>* dS_dT)
{
  // total derivative of stress (mechanical + thermal part) wrt. temperature
  // calculate derivative of cmat w.r.t. T_{n+1}
  Core::LinAlg::Matrix<6, 6> cmat_T(false);
  get_cmat_at_tempnp_t(cmat_T);

  // evaluate meachnical stress part
  // \f \sigma = {\mathbf C}_{,T} \,\varepsilon_{\rm GL} \f
  dS_dT->MultiplyNN(cmat_T, *current_glstrain_);

  // calculate the temperature difference
  // Delta T = T - T_0
  const double deltaT = current_temperature_ - params_->thetainit_;

  // calculate derivative of ctemp w.r.t. T_{n+1}
  Core::LinAlg::Matrix<6, 1> ctemp_T(false);
  get_cthermo_at_tempnp_t(ctemp_T);

  // temperature dependent stress part
  // sigma = C_T . Delta T = m . I . Delta T
  dS_dT->Update(deltaT, ctemp_T, 1.0);

  setup_cthermo(ctemp_T);
  dS_dT->Update(1.0, ctemp_T, 1.0);
}

void Mat::ThermoStVenantKirchhoff::stress_temperature_modulus_and_deriv(
    Core::LinAlg::Matrix<6, 1>& stm, Core::LinAlg::Matrix<6, 1>& stm_dT)
{
  setup_cthermo(stm);
  get_cthermo_at_tempnp_t(stm_dT);
}

/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 02/10 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::setup_cmat(Core::LinAlg::Matrix<6, 6>& cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  double Emod = 0.0;

  if (youngs_is_temp_dependent())
  {
    Emod = get_mat_parameter_at_tempnp(&(params_->youngs_), current_temperature_);
  }
  // young's modulus is constant
  else
    Emod = params_->youngs_[0];

  // Poisson's ratio (Querdehnzahl)
  const double nu = params_->poissonratio_;

  StVenantKirchhoff::FillCmat(cmat, Emod, nu);
}

/*----------------------------------------------------------------------*
 | calculates stress-temperature modulus                     dano 04/10 |
 *----------------------------------------------------------------------*/
double Mat::ThermoStVenantKirchhoff::st_modulus() const
{
  const double Emod = youngs_is_temp_dependent()
                          ? get_mat_parameter_at_tempnp(&(params_->youngs_), current_temperature_)
                          : params_->youngs_[0];

  // initialise the parameters for the lame constants
  const double nu = params_->poissonratio_;

  // initialise the thermal expansion coefficient
  const double thermexpans = params_->thermexpans_;

  // plane strain, rotational symmetry
  // E / (1+nu)
  const double c1 = Emod / (1.0 + nu);
  // (E*nu) / ((1+nu)(1-2nu))
  const double b1 = c1 * nu / (1.0 - 2.0 * nu);

  // build the lame constants
  //            E
  //   mu = --------
  //        2*(1+nu)
  //                  E*nu
  //   lambda = ----------------
  //            (1+nu)*(1-2*nu)
  //
  //  \f \mu =  \frac{E}{2(1+\nu)} \f
  const double mu = 0.5 * c1;
  // lambda
  // \f \frac{E\,\nu}{(1-2\nu)(1+\nu)} \f
  const double lambda = b1;

  // stress-temperature modulus
  // \f m\, = \, -(2\,\cdot \mu \, +\, 3\cdot\lambda)\cdot\varalpha_T \f
  const double stmodulus = (-1.0) * (2.0 * mu + 3.0 * lambda) * thermexpans;

  return stmodulus;

}  // st_modulus()


/*----------------------------------------------------------------------*
 | computes temperature dependent isotropic                  dano 05/10 |
 | elasticity tensor in matrix notion for 3d, second(!) order tensor    |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::setup_cthermo(Core::LinAlg::Matrix<6, 1>& ctemp)
{
  double m = st_modulus();
  FillCthermo(ctemp, m);
}

void Mat::ThermoStVenantKirchhoff::FillCthermo(Core::LinAlg::Matrix<6, 1>& ctemp, double m)
{
  // isotropic elasticity tensor C_temp in Voigt matrix notation C_temp = m I
  //
  // Matrix-notation for 3D case
  //       [ m      0      0 ]
  // C_T = [ 0      m      0 ]
  //       [ 0      0      m ]
  //
  // in Vector notation
  // C_T = [m, m, m, 0, 0, 0]^T
  //
  // write non-zero components

  // clear the material tangent, equal to PutScalar(0.0), but faster
  ctemp.Clear();

  // loop over the element nodes, non-zero entries only in main directions
  for (int i = 0; i < 3; ++i) ctemp(i, 0) = m;
  // else zeros
}
// setup_cthermo()


/*----------------------------------------------------------------------*
 | return temperature-dependent material parameter           dano 01/13 |
 | at current temperature --> polynomial type, cf. robinson material    |
 *----------------------------------------------------------------------*/
double Mat::ThermoStVenantKirchhoff::get_mat_parameter_at_tempnp(
    const std::vector<double>* paramvector,  // (i) given parameter is a vector
    const double& tempnp                     // tmpr (i) current temperature
) const
{
  // polynomial type

  // initialise the temperature dependent material parameter
  double parambytempnp = 0.0;
  double tempnp_pow = 1.0;

  // Param = a + b . T + c . T^2 + d . T^3 + ...
  // with T: current temperature
  for (double i : (*paramvector))
  {
    // calculate coefficient of variable T^i
    parambytempnp += i * tempnp_pow;
    // for the higher polynom increase the exponent of the temperature
    tempnp_pow *= tempnp;
  }

  // return temperature-dependent material parameter
  return parambytempnp;

}  // get_mat_parameter_at_tempnp()


/*----------------------------------------------------------------------*
 | calculate derivative of material parameter with respect   dano 01/13 |
 | to the current temperature --> polynomial type                       |
 *----------------------------------------------------------------------*/
double Mat::ThermoStVenantKirchhoff::get_mat_parameter_at_tempnp_t(
    const std::vector<double>* paramvector,  // (i) given parameter is a vector
    const double& tempnp                     // tmpr (i) current temperature
) const
{
  // polynomial type

  // initialise the temperature dependent material parameter
  double parambytempnp = 0.0;
  double tempnp_pow = 1.0;

  // Param = a + b . T + c . T^2 + d . T^3 + ...
  //       = a + b . N_T . T_{n+1} + c . (N_T . T_{n+1})^2 + d . (N_T . T_{n+1})^3 + ...
  // with T: current scalar-valued temperature, T = N_T . T_{n+1}

  // calculate derivative of E(T_{n+1}) w.r.t. T_{n+1}
  // d(Param)/dT . Delta T = b . N_T + 2 . c . T . N_T + 3 . d . T^2 . N_T + ...
  //                       = ( b + 2 . c . T + 3 . d . T^2 + ...) . N_T
  //                       = parambytempnp . N_T

  // the first, constant term has no influence for the derivative
  // --> start with i=1(!): neglect first term of paramvector
  for (unsigned i = 1; i < (*paramvector).size(); ++i)
  {
    // calculate coefficient of variable T^i
    parambytempnp += i * (*paramvector)[i] * tempnp_pow;
    // for the higher polynom increase the exponent of the temperature
    tempnp_pow *= tempnp;
  }

  // return derivative of temperature-dependent material parameter w.r.t. T_{n+1}
  return parambytempnp;

}  // get_mat_parameter_at_tempnp_t()


/*----------------------------------------------------------------------*
 | calculate linearisation of stress-temperature modulus     dano 04/10 |
 | w.r.t. T_{n+1} for k_dT, k_TT                                        |
 *----------------------------------------------------------------------*/
double Mat::ThermoStVenantKirchhoff::get_st_modulus_t() const
{
  // build the derivative of the stress-temperature modulus w.r.t. T_{n+1}
  // m = - (2 . mu + 3 . lambda) . varalpha_T
  //   = - (2 . nu / ((1+nu)(1-2nu)) + 3 / (2 . (1+nu))) . varalpha_T . E(T)
  double stmodulus_T = 0.0;


  if (youngs_is_temp_dependent())
  {
    const double Ederiv = get_mat_parameter_at_tempnp_t(&(params_->youngs_), current_temperature_);

    // initialise the parameters for the lame constants
    const double nu = params_->poissonratio_;

    // initialise the thermal expansion coefficient
    const double thermexpans = params_->thermexpans_;

    // plane strain, rotational symmetry
    // E / (1+nu)
    const double c1 = Ederiv / (1.0 + nu);
    // (E . nu) / ((1+nu)(1-2nu))
    const double b1 = c1 * nu / (1.0 - 2.0 * nu);

    // build the lame constants
    //         E
    // mu = --------  --> \f \mu = \frac{E}{2(1+\nu)} \f
    //      2*(1+nu)
    const double mu = 0.5 * c1;
    //              E*nu
    // lambda = --------------- --> \f \frac{E\,\nu}{(1-2\nu)(1+\nu)} \f
    //          (1+nu)*(1-2*nu)
    const double lambda = b1;

    // build the derivative of the stress-temperature modulus w.r.t. T_{n+1}
    // m = -(2 . mu + 3 . lambda) . varalpha_T
    stmodulus_T = (-1.0) * (2.0 * mu + 3.0 * lambda) * thermexpans;
  }
  // else (young_temp == false)
  // constant young's modulus, i.e. independent of T, no linearisation, return 0;

  return stmodulus_T;

}  // get_st_modulus_t()

/*----------------------------------------------------------------------*
 | computes thermal derivative of the isotropic elasticity   dano 01/13 |
 | tensor in matrix notion for 3d for k_dT                              |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::get_cmat_at_tempnp_t(Core::LinAlg::Matrix<6, 6>& derivcmat)
{
  // clear the material tangent, identical to PutScalar(0.0)
  derivcmat.Clear();

  if (youngs_is_temp_dependent())
  {
    const double Ederiv = get_mat_parameter_at_tempnp_t(&(params_->youngs_), current_temperature_);
    // Poisson's ratio (Querdehnzahl)
    const double nu = params_->poissonratio_;

    StVenantKirchhoff::FillCmat(derivcmat, Ederiv, nu);
  }
}


/*----------------------------------------------------------------------*
 | computes linearisation of temperature-dependent isotropic dano 01/13 |
 | elasticity tensor in matrix notion for 3d, 2nd order tensor          |
 *----------------------------------------------------------------------*/
void Mat::ThermoStVenantKirchhoff::get_cthermo_at_tempnp_t(Core::LinAlg::Matrix<6, 1>& derivctemp)
{
  double m_T = get_st_modulus_t();

  FillCthermo(derivctemp, m_T);
}


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
