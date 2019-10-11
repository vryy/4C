/*----------------------------------------------------------------------*/
/*! \file
\brief Material for a 1D artery, contains its initial diameter, thickness, dynamic
       viscosity and density of the fluid flowing in it, Young's modulus and Poisson ratio and
       external constant tissue pressures for the nodes

\level 3

\maintainer Johannes Kremheller
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "cnst_1d_art.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Cnst_1d_art::Cnst_1d_art(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("VISCOSITY")),
      density_(matdata->GetDouble("DENS")),
      young_(matdata->GetDouble("YOUNG")),
      nue_(matdata->GetDouble("NUE")),
      th_(matdata->GetDouble("TH")),
      pext1_(matdata->GetDouble("PEXT1")),
      pext2_(matdata->GetDouble("PEXT2")),
      viscositylaw_(viscositylaw_undefined),
      blood_visc_scale_diam_to_microns_(matdata->GetDouble("BLOOD_VISC_SCALE_DIAM_TO_MICRONS"))
{
  const std::string* typestring = matdata->Get<std::string>("VISCOSITYLAW");

  if (*typestring == "CONSTANT")
    viscositylaw_ = constant;
  else if (*typestring == "BLOOD")
    viscositylaw_ = blood;
}

Teuchos::RCP<MAT::Material> MAT::PAR::Cnst_1d_art::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Cnst_1d_art(this));
}


MAT::Cnst_1d_artType MAT::Cnst_1d_artType::instance_;


DRT::ParObject* MAT::Cnst_1d_artType::Create(const std::vector<char>& data)
{
  MAT::Cnst_1d_art* cnst_art = new MAT::Cnst_1d_art();
  cnst_art->Unpack(data);
  return cnst_art;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Cnst_1d_art::Cnst_1d_art() : params_(NULL), diam_(0.0) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Cnst_1d_art::Cnst_1d_art(MAT::PAR::Cnst_1d_art* params) : params_(params), diam_(0.0) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Cnst_1d_art::Pack(DRT::PackBuffer& data) const
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
  AddtoPack(data, diam_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Cnst_1d_art::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Cnst_1d_art*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // diameter
  ExtractfromPack(position, data, diam_);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Cnst_1d_art::Viscosity() const
{
  switch (params_->viscositylaw_)
  {
    case MAT::PAR::ArteryViscosityLaw::constant:
      return params_->viscosity_;
    case MAT::PAR::ArteryViscosityLaw::blood:
      return CalculateBloodViscosity(
          diam_ * params_->blood_visc_scale_diam_to_microns_, params_->viscosity_);
    default:
      dserror("Unknown viscosity law for 1D artery element");
  }
  return 0.0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Cnst_1d_art::CalculateBloodViscosity(const double diam, const double plasmavisc) const
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
