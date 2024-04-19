/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material with simplified chemical kinetics

\level 2

*/
/*----------------------------------------------------------------------*/


#include "baci_mat_ferech_pv.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::FerEchPV::FerEchPV(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      refvisc_(*matdata->Get<double>("REFVISC")),
      reftemp_(*matdata->Get<double>("REFTEMP")),
      suthtemp_(*matdata->Get<double>("SUTHTEMP")),
      pranum_(*matdata->Get<double>("PRANUM")),
      reacratecon_(*matdata->Get<double>("REACRATECON")),
      pvcrit_(*matdata->Get<double>("PVCRIT")),
      unbshc_(*matdata->Get<double>("UNBSHC")),
      burshc_(*matdata->Get<double>("BURSHC")),
      unbtemp_(*matdata->Get<double>("UNBTEMP")),
      burtemp_(*matdata->Get<double>("BURTEMP")),
      unbdens_(*matdata->Get<double>("UNBDENS")),
      burdens_(*matdata->Get<double>("BURDENS")),
      mod_(*matdata->Get<double>("MOD"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::FerEchPV::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FerEchPV(this));
}


MAT::FerEchPVType MAT::FerEchPVType::instance_;


CORE::COMM::ParObject* MAT::FerEchPVType::Create(const std::vector<char>& data)
{
  MAT::FerEchPV* ferech_pv = new MAT::FerEchPV();
  ferech_pv->Unpack(data);
  return ferech_pv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::FerEchPV::FerEchPV() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::FerEchPV::FerEchPV(MAT::PAR::FerEchPV* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::FerEchPV::Pack(CORE::COMM::PackBuffer& data) const
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::FerEchPV::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::FerEchPV*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeTemperature(const double provar) const
{
  const double temperature = UnbTemp() + provar * (BurTemp() - UnbTemp());

  return temperature;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeDensity(const double provar) const
{
  // BML hypothesis
  const double density = UnbDens() + provar * (BurDens() - UnbDens());

  // equation of state
  // const double density = UnbDens()*BurDens()/(BurDens() + provar * (UnbDens() - BurDens()));

  return density;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeFactor(const double provar) const
{
  // BML hypothesis
  const double factor = (UnbDens() - BurDens()) / (UnbDens() + provar * (BurDens() - UnbDens()));

  // equation of state
  // const double factor = (UnbDens() - BurDens())/(BurDens() + provar * (UnbDens() - BurDens()));

  return factor;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeShc(const double provar) const
{
  const double shc = UnbShc() + provar * (BurShc() - UnbShc());

  return shc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeViscosity(const double temp) const
{
  // modified version by Poinsot and Veynante (2005): constant viscosity
  double visc = RefVisc();

  // original version by Ferziger and Echekki (1993): linear variation
  if (Mod() < 1e-15) visc *= temp / RefTemp();
  // modified version by Hartmann et al. (2010): Sutherland law
  else if (Mod() > (1.0 + 1e-15))
    visc *= std::pow((temp / RefTemp()), 1.5) * ((RefTemp() + SuthTemp()) / (temp + SuthTemp()));

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeDiffusivity(const double temp) const
{
  // modified version by Poinsot and Veynante (2005): constant diffusivity
  double diffus = RefVisc() / PraNum();

  // original version by Ferziger and Echekki (1993): linear variation
  if (Mod() < 1e-15) diffus *= temp / RefTemp();
  // modified version by Hartmann et al. (2010): Sutherland law
  else if (Mod() > (1.0 + 1e-15))
    diffus *= std::pow((temp / RefTemp()), 1.5) * ((RefTemp() + SuthTemp()) / (temp + SuthTemp()));

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeReactionCoeff(const double provar) const
{
  // no reaction coefficient if progress variable lower than critical value
  double reacoeff = 0.0;

  // Heaviside function loop
  if (provar > PvCrit())
  {
    // original version by Ferziger and Echekki (1993)
    reacoeff = ReacRateCon();

    // modified version by Poinsot and Veynante (2005)
    if ((Mod() > (1.0 - 1e-15)) and (Mod() < (1.0 + 1e-15)))
    {
      // BML hypothesis
      reacoeff *= UnbDens() / (UnbDens() + provar * (BurDens() - UnbDens()));

      // equation of state
      // reacoeff *= (BurDens() + provar * (UnbDens() - BurDens()))/BurDens();
    }
    else if (Mod() > (1.0 + 1e-15))
    {
      // modified version by Hartmann et al. (2010)
      reacoeff *= UnbDens() * UnbDens() / (BurDens() + provar * (UnbDens() - BurDens()));
    }
  }

  return reacoeff;
}

FOUR_C_NAMESPACE_CLOSE
