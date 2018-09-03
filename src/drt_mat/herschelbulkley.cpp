/*----------------------------------------------------------------------*/
/*!
\file herschelbulkley.cpp

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/


#include <vector>

#include "herschelbulkley.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::HerschelBulkley::HerschelBulkley(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      tau0_(matdata->GetDouble("TAU_0")),
      kfac_(matdata->GetDouble("KFAC")),
      nexp_(matdata->GetDouble("NEXP")),
      mexp_(matdata->GetDouble("MEXP")),
      lolimshearrate_(matdata->GetDouble("LOLIMSHEARRATE")),
      uplimshearrate_(matdata->GetDouble("UPLIMSHEARRATE")),
      density_(matdata->GetDouble("DENSITY"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::HerschelBulkley::CreateMaterial()
{
  return Teuchos::rcp(new MAT::HerschelBulkley(this));
}


MAT::HerschelBulkleyType MAT::HerschelBulkleyType::instance_;


DRT::ParObject* MAT::HerschelBulkleyType::Create(const std::vector<char>& data)
{
  MAT::HerschelBulkley* herbul = new MAT::HerschelBulkley();
  herbul->Unpack(data);
  return herbul;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::HerschelBulkley::HerschelBulkley() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::HerschelBulkley::HerschelBulkley(MAT::PAR::HerschelBulkley* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::HerschelBulkley::Pack(DRT::PackBuffer& data) const
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
void MAT::HerschelBulkley::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::HerschelBulkley*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
