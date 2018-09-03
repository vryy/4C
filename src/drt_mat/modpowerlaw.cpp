/*----------------------------------------------------------------------*/
/*!
\file modpowerlaw.cpp

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235

\level 3
*/
/*----------------------------------------------------------------------*/


#include <vector>

#include "modpowerlaw.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

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


DRT::ParObject* MAT::ModPowerLawType::Create(const std::vector<char>& data)
{
  MAT::ModPowerLaw* powLaw = new MAT::ModPowerLaw();
  powLaw->Unpack(data);
  return powLaw;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ModPowerLaw::ModPowerLaw() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ModPowerLaw::ModPowerLaw(MAT::PAR::ModPowerLaw* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ModPowerLaw::Pack(DRT::PackBuffer& data) const
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
/*----------------------------------------------------------------------*/
void MAT::ModPowerLaw::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
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
        params_ = static_cast<MAT::PAR::ModPowerLaw*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
