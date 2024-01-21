/*----------------------------------------------------------------------*/
/*! \file
\brief
Nonlinear viscosity according to a modified power law


\level 3
*/
/*----------------------------------------------------------------------*/


#include "baci_mat_modpowerlaw.H"

#include "baci_global_data.H"
#include "baci_mat_par_bundle.H"

#include <vector>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ModPowerLaw::ModPowerLaw(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      m_cons_(matdata->GetDouble("MCONS")),
      delta_(matdata->GetDouble("DELTA")),
      a_exp_(matdata->GetDouble("AEXP")),
      density_(matdata->GetDouble("DENSITY"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ModPowerLaw::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ModPowerLaw(this));
}

MAT::ModPowerLawType MAT::ModPowerLawType::instance_;


CORE::COMM::ParObject* MAT::ModPowerLawType::Create(const std::vector<char>& data)
{
  MAT::ModPowerLaw* powLaw = new MAT::ModPowerLaw();
  powLaw->Unpack(data);
  return powLaw;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ModPowerLaw::ModPowerLaw() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ModPowerLaw::ModPowerLaw(MAT::PAR::ModPowerLaw* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ModPowerLaw::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ModPowerLaw::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ModPowerLaw*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

BACI_NAMESPACE_CLOSE
