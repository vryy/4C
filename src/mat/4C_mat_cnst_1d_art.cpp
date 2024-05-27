/*----------------------------------------------------------------------*/
/*! \file
\brief Material for a 1D artery, contains its initial diameter, thickness, dynamic
       viscosity and density of the fluid flowing in it, Young's modulus and Poisson ratio and
       external constant tissue pressures for the nodes

\level 3

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_cnst_1d_art.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Cnst1dArt::Cnst1dArt(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->Get<double>("VISCOSITY")),
      density_(matdata->Get<double>("DENS")),
      young_(matdata->Get<double>("YOUNG")),
      nue_(matdata->Get<double>("NUE")),
      th_(matdata->Get<double>("TH")),
      pext1_(matdata->Get<double>("PEXT1")),
      pext2_(matdata->Get<double>("PEXT2")),
      viscositylaw_(viscositylaw_undefined),
      diameterlaw_(diameterlaw_undefined),
      blood_visc_scale_diam_to_microns_(matdata->Get<double>("BLOOD_VISC_SCALE_DIAM_TO_MICRONS")),
      diameter_law_funct_(matdata->Get<int>("VARYING_DIAMETER_FUNCTION")),
      collapse_threshold_(matdata->Get<double>("COLLAPSE_THRESHOLD"))
{
  const std::string& typestring_visc = matdata->Get<std::string>("VISCOSITYLAW");

  if (typestring_visc == "CONSTANT")
    viscositylaw_ = viscositylaw_constant;
  else if (typestring_visc == "BLOOD")
    viscositylaw_ = viscositylaw_blood;
  else
    FOUR_C_THROW(
        "wrong type of viscosity law for artery material, only CONSTANT and BLOOD are valid");

  const std::string& typestring_diam = matdata->Get<std::string>("VARYING_DIAMETERLAW");

  if (typestring_diam == "CONSTANT")
    diameterlaw_ = diameterlaw_constant;
  else if (typestring_diam == "BY_FUNCTION")
    diameterlaw_ = diameterlaw_by_function;
  else
    FOUR_C_THROW(
        "wrong type of diameter law for artery material, only CONSTANT and BY_FUNCTION are valid");
}

Teuchos::RCP<CORE::MAT::Material> MAT::PAR::Cnst1dArt::create_material()
{
  return Teuchos::rcp(new MAT::Cnst1dArt(this));
}


MAT::Cnst1dArtType MAT::Cnst1dArtType::instance_;


CORE::COMM::ParObject* MAT::Cnst1dArtType::Create(const std::vector<char>& data)
{
  MAT::Cnst1dArt* cnst_art = new MAT::Cnst1dArt();
  cnst_art->Unpack(data);
  return cnst_art;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Cnst1dArt::Cnst1dArt()
    : params_(nullptr), diam_init_(0.0), diam_(0.0), diam_previous_time_step_(0.0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Cnst1dArt::Cnst1dArt(MAT::PAR::Cnst1dArt* params)
    : params_(params), diam_init_(0.0), diam_(0.0), diam_previous_time_step_(0.0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Cnst1dArt::Pack(CORE::COMM::PackBuffer& data) const
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
  AddtoPack(data, diam_init_);
  AddtoPack(data, diam_);
  AddtoPack(data, diam_previous_time_step_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Cnst1dArt::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Cnst1dArt*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // diameter
  ExtractfromPack(position, data, diam_init_);
  ExtractfromPack(position, data, diam_);
  ExtractfromPack(position, data, diam_previous_time_step_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Cnst1dArt::Viscosity() const
{
  switch (params_->viscositylaw_)
  {
    case MAT::PAR::ArteryViscosityLaw::viscositylaw_constant:
      return params_->viscosity_;
    case MAT::PAR::ArteryViscosityLaw::viscositylaw_blood:
      return calculate_blood_viscosity(
          diam_ * params_->blood_visc_scale_diam_to_microns_, params_->viscosity_);
    default:
      FOUR_C_THROW("Unknown viscosity law for 1D artery element");
  }
  return 0.0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Cnst1dArt::calculate_blood_viscosity(const double diam, const double plasmavisc) const
{
  // parameters
  const double hd = 0.45;
  const double dOff = 2.4;
  const double dCrit = 10.5;
  const double d50 = 100;
  const double eAmp = 1.1;
  const double eWidth = 0.03;
  const double ePeak = 0.6;
  const double eHD = 1.18;
  const double wMax = 2.6;

  std::vector<double> viscpar(6);

  // In vitro visocity params
  viscpar[0] = 220;
  viscpar[1] = -1.3;
  viscpar[2] = 3.2;
  viscpar[3] = -2.44;
  viscpar[4] = -0.06;
  viscpar[5] = 0.645;

  double wAs = 0.;
  if (dOff < diam)
  {
    wAs = wMax * (diam - dOff) / (diam + d50 - 2 * dOff);
  }

  double wPeak = 0.;
  if (diam > dOff && diam <= dCrit)
  {
    wPeak = eAmp * (diam - dOff) / (dCrit - dOff);
  }
  else if (dCrit < diam)
  {
    wPeak = eAmp * exp(-eWidth * (diam - dCrit));
  }

  const double wPH = wAs + wPeak * ePeak;
  const double wEFF = wAs + wPeak * (1 + hd * eHD);
  const double dPH = diam - 2 * wPH;

  // relative apparent blood viscosity for a fixed hematocrit value of 0.45
  const double eta45 = viscpar[0] * exp(viscpar[1] * dPH) + viscpar[2] +
                       viscpar[3] * exp(viscpar[4] * pow(dPH, viscpar[5]));

  // effective viscosity \eta_vivo = \eta_45 *(D/D_eff)^4
  // finally, blood viscosity = \eta_vivo * visc_plasma
  const double visc = eta45 * pow(diam / (diam - 2 * wEFF), 4) * plasmavisc;

  return visc;
}

FOUR_C_NAMESPACE_CLOSE
